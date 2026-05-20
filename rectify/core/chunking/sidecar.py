"""Read-number sidecar storage for cross-aligner read identity.

The sidecar maps a per-sample globally unique integer ``read_num`` to the
source FASTQ identity and comment fields. ``rectify split`` assigns read
numbers with a simple global counter in input FASTQ order, starting at zero.
That counter is stable across the round-robin chunk layout because the chunk
index is derived from ``read_num % n_chunks``; no per-chunk bit packing is
needed.

Schema:
  read_num: int64 primary key
  original_qname: string, FASTQ header between ``@`` and first whitespace
  fastq_comment: string, FASTQ header after first whitespace, or ``""``
  chunk_id: string, e.g. ``chunk_003_of_016``
  seq_md5: binary(16), MD5 digest of the sequence line
  qual_md5: binary(16), MD5 digest of the quality line

When PyArrow is available, this module writes true Parquet with zstd
compression and 100k-row groups. The local development environment used for
some tests may lack PyArrow; in that case a compact pickle fallback is used
behind the same API. The file extension and public contract stay unchanged.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import pickle
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional

logger = logging.getLogger(__name__)

ROW_GROUP_SIZE = 100_000
SCHEMA_VERSION = "1.0"


@dataclass(frozen=True)
class SidecarRow:
    read_num: int
    original_qname: str
    fastq_comment: str
    chunk_id: str
    seq_md5: bytes
    qual_md5: bytes


@dataclass(frozen=True)
class VerificationResult:
    read_num: int
    seq_matches: bool
    qual_matches: bool
    row: SidecarRow

    @property
    def ok(self) -> bool:
        return self.seq_matches and self.qual_matches


def split_fastq_header(header: str) -> tuple:
    """Return ``(qname, comment)`` from a FASTQ header line."""
    h = header.rstrip('\n\r')
    if h.startswith('@'):
        h = h[1:]
    cut = len(h)
    for ws in (' ', '\t'):
        i = h.find(ws)
        if 0 <= i < cut:
            cut = i
    qname = h[:cut]
    comment = h[cut:].lstrip(' \t') if cut < len(h) else ''
    return qname, comment


def format_fastq_header_with_rn(original_qname: str, fastq_comment: str, read_num: int) -> str:
    """Build a FASTQ header carrying the RN aux tag in the comment field."""
    rn = f"RN:i:{int(read_num)}"
    if fastq_comment:
        return f"@{original_qname}\t{rn}\t{fastq_comment}\n"
    return f"@{original_qname}\t{rn}\n"


def parse_rn_from_comment(comment: str) -> Optional[int]:
    """Extract ``RN:i:<int>`` from a FASTQ comment string, if present."""
    for tok in comment.split():
        if tok.startswith('RN:i:'):
            try:
                return int(tok[5:])
            except ValueError:
                return None
    return None


def parse_rn_from_fastq_header(header: str) -> tuple:
    """Return ``(qname, rn)`` from a FASTQ header line."""
    qname, comment = split_fastq_header(header)
    return qname, parse_rn_from_comment(comment)


def _md5_bytes(text: str) -> bytes:
    return hashlib.md5(text.rstrip('\n\r').encode()).digest()


def sha256_file(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        while True:
            chunk = fh.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def provenance_path_for(sidecar_path: str | Path) -> Path:
    path = Path(sidecar_path)
    return path.with_suffix('.PROVENANCE.json')


class ReadNumSidecarWriter:
    def __init__(
        self,
        path: str | Path,
        sample_id: str,
        *,
        row_group_size: int = ROW_GROUP_SIZE,
        source_fastq: str | Path | None = None,
        source_fastq_sha256: str | None = None,
        chunk_count: int | None = None,
    ):
        self.path = Path(path)
        self.sample_id = sample_id
        self.row_group_size = row_group_size
        self.source_fastq = str(source_fastq) if source_fastq is not None else None
        self.source_fastq_sha256 = source_fastq_sha256
        self.chunk_count = chunk_count
        self._buffer: List[SidecarRow] = []
        self._n_rows = 0
        self._closed = False
        self._pa = None
        self._writer = None
        self._fallback_rows: List[SidecarRow] = []
        self.path.parent.mkdir(parents=True, exist_ok=True)

        try:
            import pyarrow as pa  # type: ignore
            import pyarrow.parquet as pq  # type: ignore
            self._pa = pa
            self._schema = pa.schema([
                ('read_num', pa.int64()),
                ('original_qname', pa.string()),
                ('fastq_comment', pa.string()),
                ('chunk_id', pa.string()),
                ('seq_md5', pa.binary(16)),
                ('qual_md5', pa.binary(16)),
            ])
            self._writer = pq.ParquetWriter(
                self.path,
                self._schema,
                compression='zstd',
                compression_level=3,
                use_dictionary=['original_qname', 'chunk_id'],
            )
        except ImportError:
            self._schema = None
            logger.warning(
                "pyarrow is not available; using non-parquet fallback for %s",
                self.path,
            )

    def add(
        self,
        read_num: int,
        original_qname: str,
        fastq_comment: str,
        chunk_id: str,
        seq: str,
        qual: str,
    ) -> None:
        if self._closed:
            raise RuntimeError("ReadNumSidecarWriter is already closed")
        row = SidecarRow(
            read_num=int(read_num),
            original_qname=str(original_qname),
            fastq_comment=str(fastq_comment or ''),
            chunk_id=str(chunk_id),
            seq_md5=_md5_bytes(seq),
            qual_md5=_md5_bytes(qual),
        )
        self._buffer.append(row)
        if len(self._buffer) >= self.row_group_size:
            self._flush()

    def _flush(self) -> None:
        if not self._buffer:
            return
        rows = self._buffer
        self._buffer = []
        self._n_rows += len(rows)
        if self._writer is None:
            self._fallback_rows.extend(rows)
            return
        pa = self._pa
        table = pa.Table.from_pydict(
            {
                'read_num': [r.read_num for r in rows],
                'original_qname': [r.original_qname for r in rows],
                'fastq_comment': [r.fastq_comment for r in rows],
                'chunk_id': [r.chunk_id for r in rows],
                'seq_md5': [r.seq_md5 for r in rows],
                'qual_md5': [r.qual_md5 for r in rows],
            },
            schema=self._schema,
        )
        self._writer.write_table(table, row_group_size=self.row_group_size)

    def close(self) -> None:
        if self._closed:
            return
        self._flush()
        if self._writer is not None:
            self._writer.close()
        else:
            tmp = self.path.with_suffix(self.path.suffix + '.tmp')
            with open(tmp, 'wb') as fh:
                pickle.dump(
                    {'schema_version': SCHEMA_VERSION, 'rows': self._fallback_rows},
                    fh,
                    protocol=pickle.HIGHEST_PROTOCOL,
                )
            os.replace(tmp, self.path)
        self._write_provenance()
        self._closed = True

    def _write_provenance(self) -> None:
        prov_path = provenance_path_for(self.path)
        parquet_sha = sha256_file(self.path) if self.path.exists() else None
        prov = {
            'schema_version': SCHEMA_VERSION,
            'sample_id': self.sample_id,
            'chunker_version': 'rectify.read_num_sidecar.v1',
            'total_reads': self._n_rows,
            'chunk_count': self.chunk_count,
            'source_fastq': self.source_fastq,
            'source_fastq_sha256': self.source_fastq_sha256,
            'sidecar_path': str(self.path),
            'sidecar_sha256': parquet_sha,
        }
        tmp = prov_path.with_suffix(prov_path.suffix + '.tmp')
        with open(tmp, 'w') as fh:
            json.dump(prov, fh, indent=2, sort_keys=True)
            fh.write('\n')
        os.replace(tmp, prov_path)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
        return False


class ReadNumSidecar:
    def __init__(self, rows: Iterable[SidecarRow]):
        self._by_read_num: Dict[int, SidecarRow] = {}
        self._by_qname: Dict[str, SidecarRow] = {}
        for row in rows:
            self._by_read_num[row.read_num] = row
            self._by_qname.setdefault(row.original_qname, row)

    @classmethod
    def open(cls, path: str | Path) -> "ReadNumSidecar":
        path = Path(path)
        try:
            import pyarrow.parquet as pq  # type: ignore
            table = pq.read_table(path)
            cols = table.to_pydict()
            rows = [
                SidecarRow(
                    read_num=int(rn),
                    original_qname=str(qn),
                    fastq_comment=str(comment or ''),
                    chunk_id=str(chunk_id),
                    seq_md5=bytes(seq_md5),
                    qual_md5=bytes(qual_md5),
                )
                for rn, qn, comment, chunk_id, seq_md5, qual_md5 in zip(
                    cols['read_num'],
                    cols['original_qname'],
                    cols['fastq_comment'],
                    cols['chunk_id'],
                    cols['seq_md5'],
                    cols['qual_md5'],
                )
            ]
            return cls(rows)
        except ImportError:
            with open(path, 'rb') as fh:
                payload = pickle.load(fh)
            return cls(payload['rows'])

    def lookup(self, read_num: int) -> SidecarRow:
        return self._by_read_num[int(read_num)]

    def lookup_many(self, read_nums: Iterable[int]) -> dict:
        return {int(rn): self.lookup(int(rn)) for rn in read_nums}

    def lookup_by_qname(self, qname: str) -> SidecarRow:
        return self._by_qname[qname]

    def verify(self, read_num: int, seq: str, qual: str) -> VerificationResult:
        row = self.lookup(read_num)
        result = VerificationResult(
            read_num=int(read_num),
            seq_matches=(row.seq_md5 == _md5_bytes(seq)),
            qual_matches=(row.qual_md5 == _md5_bytes(qual)),
            row=row,
        )
        if not result.ok:
            logger.warning(
                "Read-num sidecar fingerprint mismatch for RN=%s "
                "(seq=%s, qual=%s)",
                read_num, result.seq_matches, result.qual_matches,
            )
        return result

    def __len__(self) -> int:
        return len(self._by_read_num)
