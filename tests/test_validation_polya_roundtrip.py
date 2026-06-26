"""Poly-A trim/restore round-trip integrity for the DRS validation bundle.

The bundle's artifacts must come from ONE coherent build so the poly-A
trim → align → restore loop is exact:

    dorado_source (untrimmed)  --trim-polya-->  validation_reads (trimmed) + parquet
    original_seq_len == len(dorado_source read)
    trimmed_seq_len  == len(validation_reads read)
    original - trimmed == polya_len == len(trimmed_3prime_seq)

If the parquet is from a different basecalling/build than the aligned reads
(as was the case 2026-06-25 — parquet orig 412 vs actual 399; "trimmed"
validation_reads identical in length to the untrimmed source), the restored
poly-A soft-clip is wrong and every validation figure's green tail is bogus.

This test guards that invariant. See dev handoff: coherent-bundle-rebuild.
"""
from pathlib import Path

import pysam
import pytest

import rectify

DATA = Path(rectify.__file__).parent / "data" / "validation"
DORADO_SRC = DATA / "validation_reads_dorado_source.bam"
VAL_BAM = DATA / "validation_reads.bam"
PARQUET_TSV = (
    Path(rectify.__file__).parent.parent
    / "scripts/validation_data/rebuild_2026_05/trimmed"
    / "validation_reads_polya_trim_metadata.tsv"
)


def _bam_lengths(bam_path: Path) -> dict:
    """{bare_uuid8: primary query_length}."""
    out = {}
    with pysam.AlignmentFile(str(bam_path), "rb") as bf:
        for r in bf:
            if r.is_secondary or r.is_supplementary:
                continue
            out[r.query_name[:8]] = r.query_length
    return out


def _load_parquet_tsv() -> dict:
    """{bare_uuid8: (original_seq_len, trimmed_seq_len, polya_len, tail_seq)}."""
    rows = {}
    with PARQUET_TSV.open() as f:
        hdr = f.readline().rstrip("\n").split("\t")
        oi, ti = hdr.index("original_seq_len"), hdr.index("trimmed_seq_len")
        pi, si = hdr.index("polya_len"), hdr.index("trimmed_3prime_seq")
        for line in f:
            c = line.rstrip("\n").split("\t")
            rows[c[0][:8]] = (int(c[oi]), int(c[ti]), int(c[pi]), c[si])
    return rows


@pytest.mark.skipif(not PARQUET_TSV.exists(), reason="trim parquet not present")
@pytest.mark.skipif(not DORADO_SRC.exists(), reason="dorado source BAM not present")
class TestPolyARoundTrip:
    """Each invariant is a separate assertion so failures are granular.

    These were XFAIL while the bundle carried build skew; the 2026-06-26
    coherent-source fix (dorado_source ← build-X combined source; uLTRA
    no-trim fix restoring cat9_minus_2 to its full 640 bp) makes them pass,
    so the markers were removed — they are now LIVE guards against
    re-introducing build skew or aligner read-trimming.
    """

    def test_parquet_original_matches_dorado_source(self):
        pq = _load_parquet_tsv()
        src = _bam_lengths(DORADO_SRC)
        bad = {k: (pq[k][0], src.get(k)) for k in pq
               if k in src and pq[k][0] != src[k]}
        assert not bad, f"parquet original_seq_len != dorado_source length: {bad}"

    def test_aligned_reads_are_the_trimmed_reads(self):
        pq = _load_parquet_tsv()
        val = _bam_lengths(VAL_BAM)
        bad = {k: (pq[k][1], val.get(k)) for k in pq
               if k in val and pq[k][1] != val[k]}
        assert not bad, f"validation_reads length != parquet trimmed_seq_len: {bad}"

    def test_parquet_internally_consistent(self):
        """orig - trimmed == len(trimmed_3prime_seq), and polya_len <= len(tail).

        ``polya_len`` is the poly-A base count; ``trimmed_3prime_seq`` is the full
        trimmed 3' segment (poly-A plus any trailing non-A adapter base), so the
        trimmed *length* equals len(tail), while polya_len may be a hair smaller.
        Parquet-internal — holds even though the cross-artifact checks fail.
        """
        pq = _load_parquet_tsv()
        bad = {}
        for k, (orig, trim, pa, tail) in pq.items():
            if (orig - trim) != len(tail) or pa > len(tail):
                bad[k] = (orig, trim, pa, len(tail))
        assert not bad, f"parquet internal round-trip broken (orig-trim vs len(tail); polya_len>len): {bad}"
