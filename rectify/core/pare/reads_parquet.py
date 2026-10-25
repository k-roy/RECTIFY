"""Per-read parquet sidecar for ``rectify pare`` — the netseq-cpa schema + 5'P.

Same failure-tolerant design as
:mod:`rectify.core.netseq_cpa.reads_parquet` (a missing/partial parquet is
recoverable; a lost 50M-read walk is not), with two PARE-specific columns
appended:

  five_p_pos     0-based 5'P cut-site coordinate (strand-aware read start).
                 tier=='mapped': the aligned RNA-5' terminus (see five_p_clip).
                 tier=='rescued': the re-mapped anchor's non-poly-A end.
  five_p_clip    soft/hard clip length at the RNA 5' terminus. > 0 means the
                 5'P placement is an alignment artifact (workspace invariant:
                 gate 5'-end analyses on five_p_clip == 0). Rescued rows are 0
                 by construction (the anchor is a clean re-map).

One row per PARE read that contributed to Tier-1 (ALL used mapped reads, tail
or no tail — so the (5'P, pA-junction) pairing is complete) or to the rescue
tier. Kept as a separate writer class rather than a parameterization of the
netseq one so the validated netseq arm stays untouched.
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable, List, Optional

from ..netseq_cpa.reads_parquet import BATCH_SIZE_DEFAULT, HAS_PYARROW, pa, pq

COLUMNS = [
    "read_id", "chrom", "cpa_pos", "gene_strand",
    "oaNT_tail_len", "at_cpa", "tier", "mapq",
    "five_p_pos", "five_p_clip",
]


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
        ("five_p_pos", pa.int32()),
        ("five_p_clip", pa.int16()),
    ])


class PareReadWriter:
    """Streaming per-read parquet writer for the PARE arm. Failure-tolerant.

    Use as: ``w = PareReadWriter(path, log=...)`` then ``w.add(...)`` per read,
    then ``w.close()``. ``w.total`` is the count of rows actually persisted;
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
        five_p_pos: int,
        five_p_clip: int,
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
        b["five_p_pos"].append(five_p_pos)
        b["five_p_clip"].append(five_p_clip)
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
