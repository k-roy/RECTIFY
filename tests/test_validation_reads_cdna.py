"""Validation test suite for the bundled ONT cDNA example reads.

Mirrors :mod:`tests.test_validation_reads` but for the cDNA arm. The
bundled BAM at ``rectify/data/validation/validation_reads_cdna.bam``
contains 12 primary reads from
``/Users/kevinroy/work/ont_cdna/test_data/wt_rep1.chrI.bam`` covering
three distinct cDNA correction scenarios:

  * ``cat1_cdna_polya``      — 3'-end poly-A indel artifact (4 reads,
                                strand-balanced)
  * ``cat4_cdna_false_junc`` — spurious terminal N-op near the 3' end
                                (4 reads, strand-balanced)
  * ``cat5_cdna_chimeric``   — multi-segment chimeric alignment (4 reads,
                                strand-balanced)

The tests call ``walk_back_anchor_and_tail`` directly on each read and
assert the (anchor_ref_pos, tail_len) tuple. Expected values were
locked in by ``scripts/validation_data/characterize_baseline.py`` on
the simple-core wiring (no artifact guards). When the cDNA wrapper is
later wired through ``walkback_3prime_guarded_with_qpos`` the expected
values may shift; the test file is the place to record those changes.

Build the BAM with::

    python scripts/validation_data/build_validation_bam.py \\
        --manifest scripts/validation_data/cdna_manifest.tsv \\
        --output rectify/data/validation/validation_reads_cdna.bam
"""
from __future__ import annotations

from pathlib import Path

import pysam
import pytest

from rectify.config import CHROM_TO_GENOME


# Per-XV expected (anchor_ref_pos, tail_len) tuples, locked from the
# simple-core baseline. The reads were classified into XG categories by
# manual inspection of the source BAM; the (anchor, tail) values reflect
# what ``walk_back_anchor_and_tail`` returns today with no artifact guards.
EXPECTED_ANCHOR_AND_TAIL = {
    # cat1_cdna_polya — terminal poly-A indel scenario
    "cdna_cat1_plus_1":  (6330, 32),
    "cdna_cat1_plus_2":  (8796, 23),
    "cdna_cat1_minus_1": (7040, 36),
    "cdna_cat1_minus_2": (7113, 18),
    # cat4_cdna_false_junc — terminal N-op near 3' end
    "cdna_cat4_plus_1":  (35763, 0),
    "cdna_cat4_plus_2":  (73445, 36),
    "cdna_cat4_minus_1": (45832, 32),
    "cdna_cat4_minus_2": (71758, 32),
    # cat5_cdna_chimeric — multi-segment SA alignment
    "cdna_cat5_plus_1":  (34332, 44),
    "cdna_cat5_plus_2":  (59577, 80),
    "cdna_cat5_minus_1": (35106, 76),
    "cdna_cat5_minus_2": (37146, 58),
}

EXPECTED_CATEGORIES = {
    "cdna_cat1_plus_1":  "cat1_cdna_polya",
    "cdna_cat1_plus_2":  "cat1_cdna_polya",
    "cdna_cat1_minus_1": "cat1_cdna_polya",
    "cdna_cat1_minus_2": "cat1_cdna_polya",
    "cdna_cat4_plus_1":  "cat4_cdna_false_junc",
    "cdna_cat4_plus_2":  "cat4_cdna_false_junc",
    "cdna_cat4_minus_1": "cat4_cdna_false_junc",
    "cdna_cat4_minus_2": "cat4_cdna_false_junc",
    "cdna_cat5_plus_1":  "cat5_cdna_chimeric",
    "cdna_cat5_plus_2":  "cat5_cdna_chimeric",
    "cdna_cat5_minus_1": "cat5_cdna_chimeric",
    "cdna_cat5_minus_2": "cat5_cdna_chimeric",
}


def _bam_path() -> Path:
    p = Path(__file__).resolve().parent.parent / "rectify" / "data" / "validation" / "validation_reads_cdna.bam"
    if not p.exists():
        pytest.skip(f"cDNA validation BAM not present: {p}")
    return p


def _genome_path() -> Path:
    p = (
        Path(__file__).resolve().parent.parent / "rectify" / "data" / "genomes"
        / "saccharomyces_cerevisiae"
        / "S288C_reference_sequence_R64-5-1_20240529.fsa.gz"
    )
    if not p.exists():
        pytest.skip(f"Bundled S288C genome not present: {p}")
    return p


@pytest.fixture(scope="module")
def genome():
    from rectify.utils.genome import load_genome
    return load_genome(str(_genome_path()))


@pytest.fixture(scope="module")
def reads_by_xv():
    """Load all primary reads from the bundled cDNA BAM, keyed by XV label."""
    with pysam.AlignmentFile(str(_bam_path()), "rb") as bam:
        out = {}
        for r in bam:
            if r.is_secondary or r.is_supplementary or r.is_unmapped:
                continue
            try:
                xv = r.get_tag("XV")
            except KeyError:
                continue
            out[xv] = r
    return out


# ---------------------------------------------------------------------------
# BAM integrity
# ---------------------------------------------------------------------------
class TestCdnaBamIntegrity:
    def test_twelve_primary_reads(self, reads_by_xv):
        assert len(reads_by_xv) == 12, (
            f"Expected 12 reads in cDNA validation BAM, got {len(reads_by_xv)}"
        )

    def test_all_labels_present(self, reads_by_xv):
        assert set(reads_by_xv.keys()) == set(EXPECTED_ANCHOR_AND_TAIL.keys())

    def test_category_tags(self, reads_by_xv):
        for xv, expected_xg in EXPECTED_CATEGORIES.items():
            assert reads_by_xv[xv].get_tag("XG") == expected_xg, (
                f"{xv}: XG mismatch (expected {expected_xg}, "
                f"got {reads_by_xv[xv].get_tag('XG')})"
            )

    def test_strand_balance(self, reads_by_xv):
        """Each category has 2 plus + 2 minus reads."""
        from collections import defaultdict
        per_cat_strand = defaultdict(lambda: {"+": 0, "-": 0})
        for xv, r in reads_by_xv.items():
            cat = r.get_tag("XG")
            strand = "-" if r.is_reverse else "+"
            per_cat_strand[cat][strand] += 1
        for cat, counts in per_cat_strand.items():
            assert counts["+"] == 2, f"{cat}: expected 2 plus, got {counts['+']}"
            assert counts["-"] == 2, f"{cat}: expected 2 minus, got {counts['-']}"


# ---------------------------------------------------------------------------
# Walkback characterization — baseline values from simple-core wiring
# ---------------------------------------------------------------------------
def _resolve_chrom_seq(genome, chrom):
    return genome.get(chrom) or genome.get(CHROM_TO_GENOME.get(chrom, ""), "")


class TestCdnaWalkbackBaseline:
    """Each XV-labelled read must produce the locked-in (anchor, tail_len)
    when run through the current cDNA walkback wrapper. The expected values
    were captured by
    ``scripts/validation_data/characterize_baseline.py --arm cdna``.

    The XV label encodes which strand pattern the wrapper should follow:
    ``*_plus_*`` reads → ``orient='fwd'``, ``*_minus_*`` reads →
    ``orient='rev'``. (In the live cDNA pipeline ``orient`` is derived from
    SSP/UMI bridge detection; here we use the strand of the test read as a
    deterministic proxy.)
    """

    @pytest.mark.parametrize("xv", list(EXPECTED_ANCHOR_AND_TAIL.keys()))
    def test_anchor_and_tail(self, xv, reads_by_xv, genome):
        from rectify.core.commands.cdna_correct_command import walk_back_anchor_and_tail

        read = reads_by_xv[xv]
        chrom_seq = _resolve_chrom_seq(genome, read.reference_name)
        assert chrom_seq, f"No genome sequence for {read.reference_name}"
        orient = "fwd" if not read.is_reverse else "rev"
        anchor, tail_len = walk_back_anchor_and_tail(read, chrom_seq, orient)
        expected_anchor, expected_tail = EXPECTED_ANCHOR_AND_TAIL[xv]
        assert (anchor, tail_len) == (expected_anchor, expected_tail), (
            f"{xv} ({orient}, is_reverse={int(read.is_reverse)}): "
            f"expected anchor={expected_anchor} tail={expected_tail}, "
            f"got anchor={anchor} tail={tail_len}"
        )
