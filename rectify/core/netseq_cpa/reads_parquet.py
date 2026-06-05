"""Per-read parquet sidecar for ``rectify netseq-cpa`` — the unified second arm.

The Tier-1 pileup writes one row per *position* (``pileup.tsv.gz``). This module
adds the orthogonal *per-read* view (``reads.parquet``): one row per NET-seq read
that contributed a CPA call, carrying the fields a downstream poly(A)/CPA analysis
needs without re-walking the BAM:

  read_id        QNAME (SRA description + ``_pa<N>`` rescue suffix stripped)
  chrom          standardized (chrom-map applied) reference name
  cpa_pos        0-based corrected CPA coordinate (walkback output)
  gene_strand    RNA/gene strand ('+'/'-'), NOT the BAM is_reverse strand
  oaNT_tail_len  decontaminated non-templated poly(A) length (capped at maxlen)
  at_cpa         within +-reg of a DRS CPA modal position; **null** when no
                 ``--drs-clusters`` map was supplied (unknowable != not-at-CPA)
  tier           'mapped' (Tier-1) | 'rescued' (poly-T-trim re-map)
  mapq           alignment MAPQ. NB tier=='rescued' rows FLOOR at the rescue
                 ``min_mapq`` (default 3) — quantify_rescue drops lower — so the
                 mapq column is not a uniform-quality axis across tiers.

Design — **never let the sidecar cost the science output**:
  * ``pyarrow`` missing at import -> writer constructs as a no-op (logs once).
  * a write/flush exception mid-run -> the writer DISABLES itself, logs loudly,
    and returns; it never propagates. The caller emits per read INSIDE the
    Tier-1 loop, but ``pileup.tsv.gz`` is only written after that loop finishes,
    so a propagating parquet error would throw away the entire (expensive) pass.
    A missing/partial parquet is recoverable; a lost 50M-read walk is not.

Streaming: rows buffer in column lists and flush as a parquet RecordBatch every
``batch_size`` reads (bounded memory on 50M+ read BAMs). chrom/gene_strand/tier
are Parquet dictionary-encoded; snappy compression.
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable, List, Optional

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
    HAS_PYARROW = True
except Exception:  # pragma: no cover - exercised only on envs without pyarrow
    pa = None
    pq = None
    HAS_PYARROW = False

COLUMNS = [
    "read_id", "chrom", "cpa_pos", "gene_strand",
    "oaNT_tail_len", "at_cpa", "tier", "mapq",
]

BATCH_SIZE_DEFAULT = 1_000_000


def _schema():
    return pa.schema([
        ("read_id", pa.string()),
        ("chrom", pa.string()),
        ("cpa_pos", pa.int32()),
        ("gene_strand", pa.string()),
        ("oaNT_tail_len", pa.int16()),
        ("at_cpa", pa.bool_()),          # nullable: null = no DRS map supplied
        ("tier", pa.string()),
        ("mapq", pa.int16()),
    ])


class ParquetReadWriter:
    """Streaming per-read parquet writer. Failure-tolerant by construction.

    Use as: ``w = ParquetReadWriter(path, log=...)`` then ``w.add(...)`` per read,
    then ``w.close()``. ``w.total`` is the count of rows actually persisted.
    ``w.disabled`` is True if pyarrow was absent or a write failed (the run
    continues regardless).
    """

    def __init__(
        self,
        path: str | Path,
        *,
        batch_size: int = BATCH_SIZE_DEFAULT,
        log: Optional[Callable[[str], None]] = None,
    ) -> None:
        self.path = str(path)
        self.batch_size = max(1, int(batch_size))
        self._log = log or (lambda _m: None)
        self._buf: dict[str, List] = {c: [] for c in COLUMNS}
        self._n = 0
        self.total = 0
        self.disabled = not HAS_PYARROW
        self._writer = None
        self._schema = _schema() if HAS_PYARROW else None
        if not HAS_PYARROW:
            self._log("reads.parquet disabled: pyarrow not available "
                      "(pip install --only-binary=:all: pyarrow)")

    def add(
        self,
        *,
        read_id: str,
        chrom: str,
        cpa_pos: int,
        gene_strand: str,
        oaNT_tail_len: int,
        at_cpa: Optional[bool],
        tier: str,
        mapq: int,
    ) -> None:
        """Buffer one per-read record. Never raises; disables on first failure."""
        if self.disabled:
            return
        b = self._buf
        b["read_id"].append(read_id)
        b["chrom"].append(chrom)
        b["cpa_pos"].append(cpa_pos)
        b["gene_strand"].append(gene_strand)
        b["oaNT_tail_len"].append(oaNT_tail_len)
        b["at_cpa"].append(at_cpa)
        b["tier"].append(tier)
        b["mapq"].append(mapq)
        self._n += 1
        if self._n >= self.batch_size:
            self._flush()

    def _flush(self) -> None:
        if self.disabled or self._n == 0:
            self._reset()
            return
        try:
            arrays = [pa.array(self._buf[c], type=self._schema.field(c).type)
                      for c in COLUMNS]
            batch = pa.record_batch(arrays, schema=self._schema)
            if self._writer is None:
                self._writer = pq.ParquetWriter(
                    self.path, self._schema, compression="snappy",
                    use_dictionary=["chrom", "gene_strand", "tier"],
                )
            self._writer.write_batch(batch)
            self.total += self._n
        except Exception as e:  # disable, keep the run alive (see module docstring)
            self.disabled = True
            self._log(f"reads.parquet WRITE FAILED after {self.total} rows: {e}; "
                      "disabling parquet sidecar, run continues")
            self._close_writer_quietly()
        finally:
            self._reset()

    def _reset(self) -> None:
        for v in self._buf.values():
            v.clear()
        self._n = 0

    def _close_writer_quietly(self) -> None:
        if self._writer is not None:
            try:
                self._writer.close()
            except Exception:
                pass
            self._writer = None

    def close(self) -> int:
        """Flush the tail, close the file, return rows persisted (``total``)."""
        if not self.disabled:
            self._flush()
        self._close_writer_quietly()
        return self.total
