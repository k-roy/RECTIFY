"""Unit tests for rectify.core.bam.read_edits.realign_exon_blocks.

Exercises the targeted X-at-HP pre-check added in the lazy-consensus WT
(Stage 1 of lazy_corrected_consensus_plan.md).  Five cases:

1. Two exon blocks separated by N — only the one with X-at-HP is realigned.
2. X op exists but the ref position is NOT a homopolymer (run < min_run=3).
3. Clean read (no X ops at all) — returns False immediately.
4. Block with X-at-HP but larger than max_block_bp=120 — skipped, no realign.
5. Leading soft-clip followed by exon block with X-at-HP — block_start
   tracking survives the S op; the exon is realigned, S is preserved.

All tests use small synthetic genomes (dict {chr: str}) and hand-constructed
pysam.AlignedSegment records — no real BAMs required.

The synthetic genome for tests 1, 2, 4, 5 is::

    positions 0-1   : GG  (non-HP)
    positions 2-4   : TTT (homopolymer run of 3)
    positions 5-6   : CG
    positions 7-10  : TTTT
    positions 11-60 : A*50 (intron filler — not directly aligned)
    positions 61-69 : CGCGTACGT (clean block B — no HP run ≥ 3)
    positions 70+   : C*200

Run with:
    pytest tests/test_realign_exon_blocks.py -v
"""

import pysam
import pytest

from rectify.core.bam.read_edits import realign_exon_blocks


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_CHROM = "chr1"

# Synthetic genome string.
# T run at positions 2–4 (length 3), used by tests 1, 4, 5.
# Block-B region at 61–69 is CGCGTACGT — no run ≥ 3 anywhere in that span.
_GENOME_STR = (
    "GG"          # 0-1   non-HP
    + "TTT"       # 2-4   homopolymer run of 3
    + "CG"        # 5-6
    + "TTTT"      # 7-10
    + "A" * 50    # 11-60 intron filler
    + "CGCGTACGT" # 61-69 block B, all distinct bases
    + "C" * 200   # 70+
)
_GENOME = {_CHROM: _GENOME_STR}
_GENOME_LEN = len(_GENOME_STR)

_HEADER = pysam.AlignmentHeader.from_dict({
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"LN": _GENOME_LEN, "SN": _CHROM}],
})


def _make_read(name, ref_start, query_seq, cigartuples):
    """Return a minimal mapped AlignedSegment.  Sanity-checked for common pitfalls."""
    read = pysam.AlignedSegment(_HEADER)
    read.query_name = name
    read.flag = 0
    read.reference_id = 0  # index into header SQ list -> chr1
    read.reference_start = ref_start
    read.mapping_quality = 60
    read.query_sequence = query_seq
    read.cigartuples = cigartuples
    # Sanity guard: misconfigurations cause silent no-ops rather than exceptions.
    assert not read.is_unmapped, "read is_unmapped — construction error"
    assert read.reference_name == _CHROM, f"got {read.reference_name!r}"
    return read


def _find_n_idx(cigartuples):
    """Return the index of the first N (intron skip) op, or None."""
    for i, (op, _) in enumerate(cigartuples):
        if op == 3:
            return i
    return None


# ---------------------------------------------------------------------------
# Test 1 — only the block with X-at-HP is realigned; the clean block is not
# ---------------------------------------------------------------------------

class TestOnlyBlockWithXAtHpIsRealigned:
    """Two exon blocks separated by an N intron.  Block A has X-at-HP (position
    inside the TTT homopolymer run); block B has only = ops.  After
    realign_exon_blocks the function must return True, block A's CIGAR ops must
    change, and block B's CIGAR ops must be byte-identical to the original.
    """

    def _build_read(self):
        # Block A spans ref positions 0–10 (11 ref bases).
        # Cigar for block A: 2= 1X 1D 7=
        #   - 2=   : ref[0:2]  = GG  matches query[0:2] = GG
        #   - 1X   : ref[2]    = T   mismatches query[2] = A  (T run len=3 >= 3: HP)
        #   - 1D   : ref[3]    = T   skipped (deletion in query)
        #   - 7=   : ref[4:11] = TCGTTTT matches query[3:10] = TCGTTTT
        # Q_block = 2+1+7 = 10, R_block = 2+1+1+7 = 11 ✓
        #
        # Block B spans ref positions 61–69 (9 ref bases).
        # Cigar: 9=  (perfect match to _GENOME_STR[61:70])
        # N intron: 50 ref bases (positions 11–60)
        qry_block_a = "GG" + "A" + "TCGTTTT"   # 10 query bases
        qry_block_b = _GENOME_STR[61:70]         # 9 query bases; exact match
        return _make_read(
            name="two_block_x_at_hp",
            ref_start=0,
            query_seq=qry_block_a + qry_block_b,
            cigartuples=[
                (7, 2), (8, 1), (2, 1), (7, 7),  # block A: 2=, 1X, 1D, 7=
                (3, 50),                           # 50-bp intron N
                (7, 9),                            # block B: 9=
            ],
        )

    def test_returns_true(self):
        read = self._build_read()
        assert realign_exon_blocks(read, _GENOME) is True

    def test_block_a_cigar_ops_changed(self):
        read = self._build_read()
        orig_n_idx = _find_n_idx(read.cigartuples)
        orig_block_a = list(read.cigartuples[:orig_n_idx])
        realign_exon_blocks(read, _GENOME)
        new_n_idx = _find_n_idx(read.cigartuples)
        new_block_a = list(read.cigartuples[:new_n_idx])
        assert new_block_a != orig_block_a, (
            f"Block A should have changed; orig={orig_block_a} new={new_block_a}"
        )

    def test_block_b_cigar_ops_unchanged(self):
        """Block B has only = ops and no X; its CIGAR must be byte-identical after the call."""
        read = self._build_read()
        orig_n_idx = _find_n_idx(read.cigartuples)
        orig_block_b = list(read.cigartuples[orig_n_idx + 1:])
        realign_exon_blocks(read, _GENOME)
        new_n_idx = _find_n_idx(read.cigartuples)
        new_block_b = list(read.cigartuples[new_n_idx + 1:])
        assert new_block_b == orig_block_b, (
            f"Block B should be unchanged; orig={orig_block_b} new={new_block_b}"
        )


# ---------------------------------------------------------------------------
# Test 2 — X op not at a homopolymer ref position → no realign
# ---------------------------------------------------------------------------

class TestXNotAtHpIsNotRealigned:
    """The block contains an X op but the reference position beneath the X has a
    homopolymer run length < min_run (default 3).  realign_exon_blocks must
    return False and leave the CIGAR unchanged.
    """

    def _build_read(self):
        # Genome has no HP run of length >= 3 around position 5 (ref='C').
        # 'GGTTTCG...' -> pos 5 = 'C', run of C = 1 (no adjacent C in _GENOME_STR).
        # Put X at pos 5 so _homo_run(5) = 1 < 3.
        # Cigar: 5= 1X 4= (10 Q, 10 R), ref positions 0–9.
        # ref[0:10] = 'GGTTTCGTTT'
        # q[5] must differ from ref[5]='C' -> use 'A'
        ref_slice = _GENOME_STR[0:10]   # GGTTTCGTTT
        qry = ref_slice[:5] + "A" + ref_slice[6:]  # change pos 5 C->A
        return _make_read(
            name="x_not_at_hp",
            ref_start=0,
            query_seq=qry,
            cigartuples=[(7, 5), (8, 1), (7, 4)],  # 5=, 1X, 4=
        )

    def test_returns_false(self):
        read = self._build_read()
        assert realign_exon_blocks(read, _GENOME) is False

    def test_cigar_unchanged(self):
        read = self._build_read()
        orig = list(read.cigartuples)
        realign_exon_blocks(read, _GENOME)
        assert list(read.cigartuples) == orig


# ---------------------------------------------------------------------------
# Test 3 — clean read with no X ops → returns False fast
# ---------------------------------------------------------------------------

class TestCleanReadReturnsFalseFast:
    """A read with only = (match) CIGAR ops and no X has nothing to realign.
    The pre-check must exit early with False, leaving the CIGAR untouched.
    """

    def _build_read(self):
        # 12= perfectly matching the first 12 bases of the genome.
        qry = _GENOME_STR[:12]
        return _make_read(
            name="clean_read_no_x",
            ref_start=0,
            query_seq=qry,
            cigartuples=[(7, 12)],  # 12=
        )

    def test_returns_false(self):
        read = self._build_read()
        assert realign_exon_blocks(read, _GENOME) is False

    def test_cigar_unchanged(self):
        read = self._build_read()
        orig = list(read.cigartuples)
        realign_exon_blocks(read, _GENOME)
        assert list(read.cigartuples) == orig


# ---------------------------------------------------------------------------
# Test 4 — block > max_block_bp → DP skipped, no realign
# ---------------------------------------------------------------------------

class TestBlockLargerThanMaxBlockBpSkipped:
    """The block's R_block exceeds max_block_bp=120.  The pre-check WILL add
    the block to target_block_starts (it sees X-at-HP), but the rebuild loop's
    'if R_block > max_block_bp: continue' skips the DP.  Function returns False
    and the CIGAR is unchanged.
    """

    def _build_read(self):
        # 126-base block: 2= 1X 123=  (Q=126, R=126 > max_block_bp=120)
        # X at ref pos 2 (ref='T', run=3 — the TTT at positions 2-4).
        # q[2] must differ from ref[2]='T' -> use 'A'.
        # ref[0:126] = 'GG' + 'TTT' + 'CG' + 'TTTT' + 'A'*50 + 'CGCGTACGT' + 'C'*52
        # That's 126 chars of the genome (positions 0-125).
        ref_slice = _GENOME_STR[0:126]
        qry = ref_slice[:2] + "A" + ref_slice[3:]  # change pos 2 T->A (X there)
        assert len(qry) == 126
        return _make_read(
            name="big_block_over_max",
            ref_start=0,
            query_seq=qry,
            cigartuples=[(7, 2), (8, 1), (7, 123)],  # 2=, 1X, 123=  (R=126 > 120)
        )

    def test_returns_false(self):
        read = self._build_read()
        assert realign_exon_blocks(read, _GENOME) is False

    def test_cigar_unchanged(self):
        read = self._build_read()
        orig = list(read.cigartuples)
        realign_exon_blocks(read, _GENOME)
        assert list(read.cigartuples) == orig

    def test_r_block_exceeds_limit(self):
        """Sanity: verify the block is actually > 120 so the skip is real."""
        # 2=, 1X, 123= -> R=2+1+123=126
        assert 2 + 1 + 123 > 120


# ---------------------------------------------------------------------------
# Test 5 — leading soft-clip; block_start tracked correctly across S op
# ---------------------------------------------------------------------------

class TestBlockStartTrackedAcrossSoftClip:
    """A 5-base leading soft-clip precedes an exon block with X-at-HP.
    The pre-check's block_start index (set to idx+1 after S) and the rebuild
    loop's block_start (captured after processing S) must agree.

    After the call:
    - The function returns True (the exon block was realigned and CIGAR changed).
    - The leading 5S op is preserved byte-for-byte in the new CIGAR.
    """

    def _build_read(self):
        # Same block-A geometry as test 1 (2=1X1D7= at ref start=0), but
        # prepend a 5-base soft-clip.  Soft-clip query bases can be anything
        # (they don't consume the reference).
        # cigar: 5S 2= 1X 1D 7=  (Q=5+10=15, R=0+11=11)
        qry_soft = "TTTTT"
        qry_block_a = "GG" + "A" + "TCGTTTT"  # 10 query bases (same as test 1)
        return _make_read(
            name="softclip_then_x_at_hp",
            ref_start=0,
            query_seq=qry_soft + qry_block_a,
            cigartuples=[
                (4, 5),                           # 5S  soft-clip
                (7, 2), (8, 1), (2, 1), (7, 7),  # block A: 2=, 1X, 1D, 7=
            ],
        )

    def test_returns_true(self):
        read = self._build_read()
        assert realign_exon_blocks(read, _GENOME) is True

    def test_exon_block_realigned(self):
        read = self._build_read()
        orig_exon = list(read.cigartuples[1:])  # everything after the 5S
        realign_exon_blocks(read, _GENOME)
        new_exon = list(read.cigartuples[1:])
        assert new_exon != orig_exon, (
            f"Exon block should have changed; orig={orig_exon} new={new_exon}"
        )

    def test_soft_clip_preserved(self):
        read = self._build_read()
        realign_exon_blocks(read, _GENOME)
        first_op = read.cigartuples[0]
        assert first_op == (4, 5), (
            f"Leading 5S op should be preserved; got {first_op}"
        )
