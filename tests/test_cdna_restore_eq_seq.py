"""Regression test for restore_eq_seq: undo samtools ``calmd -e`` (SEQ='=').

``rectify align`` runs ``samtools calmd -e`` which rewrites reference-matching M
bases to ``=`` in SEQ. The cDNA Stage-1 consensus writer used to copy ``=``
verbatim into the consensus FASTQ, producing unmappable reads that silently
dropped ~99% of cDNA molecules at Stage-C align2. ``restore_eq_seq`` recovers the
real bases from the reference (they are NOT in the BAM: MD stores match run
lengths, not bases).
"""
import pysam

from rectify.core.cdna.consensus import restore_eq_seq


class _MockFasta:
    """Duck-typed stand-in for pysam.FastaFile (restore_eq_seq only calls fetch)."""
    def __init__(self, seqs):
        self._seqs = seqs

    def fetch(self, name, start, end):
        return self._seqs[name][start:end]


def _read(cigar, seq, ref_start=2, is_reverse=False):
    h = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.0"}, "SQ": [{"SN": "chr1", "LN": 12}]}
    )
    a = pysam.AlignedSegment(h)
    a.query_name = "r"
    a.flag = 16 if is_reverse else 0
    a.reference_id = 0
    a.reference_start = ref_start
    a.mapping_quality = 60
    a.cigarstring = cigar
    a.query_sequence = seq
    return a


REF = {"chr1": "ACGTACGTACGT"}


def test_restore_eq_softclip_match_mismatch_indel():
    # aligned at pos2; span fetched = ref[2:9] = "GTACGTA"
    #   2S "NN" | 4M "====" (all match GTAC) | 1I "G" | 3M "=A=" (mismatch idx1) | 2S "TT"
    a = _read("2S4M1I3M2S", "NN" + "====" + "G" + "=A=" + "TT")
    assert a.reference_end == 9
    out = restore_eq_seq(a, _MockFasta(REF))
    # 2S kept | GTAC | I kept 'G' | 3M: G,(A keep),A | 2S kept
    assert out == "NN" + "GTAC" + "G" + "GAA" + "TT"


def test_restore_eq_with_deletion():
    # 3M "===" over ref[2:5]="GTA", 2D (ref only), 3M "===" over ref[7:10]="TAC"
    a = _read("3M2D3M", "===" + "===")
    assert a.reference_end == 2 + 3 + 2 + 3  # 10
    out = restore_eq_seq(a, _MockFasta(REF))
    assert out == "GTA" + "TAC"


def test_no_eq_passthrough():
    a = _read("6M", "GTACGT")
    assert restore_eq_seq(a, _MockFasta(REF)) == "GTACGT"


def test_no_fasta_returns_raw():
    a = _read("6M", "==A===")
    assert restore_eq_seq(a, None) == "==A==="


def test_none_seq_returns_none():
    a = _read("6M", "GTACGT")
    a.query_sequence = None
    assert restore_eq_seq(a, _MockFasta(REF)) is None
