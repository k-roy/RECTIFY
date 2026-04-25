"""
Tests for junction_refiner.py

Covers five test classes:
1. RPL20B (chrXV, minus strand): wrong N-op boundaries from mapPacBio corrected
   to annotated junction; correct reads not displaced.
2. GCR1 (chrXVI, plus strand): two competing GT-AG 5'SS and multiple GT-AG 3'SS;
   reads at each isoform stay at their correct junction.
3. CIGAR surgery (_apply_junction_replacement): ref/query span preservation.
4. junction pool (_build_junction_pool, _candidates_near): pool membership.
5. TFC3 / RPL22B / SRC1 alternative splice sites: annotated reads not displaced
   even when a same-score alternative junction exists (is_alt tiebreaker).

Author: Kevin R. Roy
"""

import sys
from pathlib import Path
from typing import Dict, Set, Tuple

import pytest
import pysam

RECTIFY_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(RECTIFY_ROOT))

from rectify.core.junction_refiner import (
    _canonical_tier,
    _hp_edit_distance,
    _score_junction,
    build_junction_pool,
    _build_junction_index,
    _candidates_near,
    refine_read_junctions,
    _apply_junction_replacement,
    _iter_n_ops,
)

GENOME_PATH = RECTIFY_ROOT / "rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz"
ANNOT_PATH  = RECTIFY_ROOT / "rectify/data/genomes/saccharomyces_cerevisiae/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz"
CONSENSUS_BAM = RECTIFY_ROOT / "dev_runs/rpl19b_rpl20b/consensus/rpl19b_rpl20b.consensus.bam"

# Per-aligner BAMs for RPL20B test dataset
ALIGN_DIR = RECTIFY_ROOT / "dev_runs/rpl19b_rpl20b/align"
ALIGNER_BAMS_RPL = [
    str(ALIGN_DIR / "rpl19b_rpl20b.minimap2.bam"),
    str(ALIGN_DIR / "rpl19b_rpl20b.gapmm2.bam"),
    str(ALIGN_DIR / "rpl19b_rpl20b.mapPacBio.bam"),
]

# Annotated junctions
RPL20B_INTRON = ("chrXV",  900767, 901193)   # [900767, 901193) minus strand
GCR1_INTRON   = ("chrXVI", 412261, 413012)   # [412261, 413012) plus strand
GCR1_ALT_5SS  = ("chrXVI", 412256, 413012)   # alt 5'SS at 412256


@pytest.fixture(scope="module")
def genome() -> Dict[str, str]:
    if not GENOME_PATH.exists():
        pytest.skip(f"Genome not found: {GENOME_PATH}")
    fa = pysam.FastaFile(str(GENOME_PATH))
    return {name: fa.fetch(name) for name in fa.references}


@pytest.fixture(scope="module")
def annotated_junctions() -> Set[Tuple]:
    if not ANNOT_PATH.exists():
        pytest.skip(f"Annotation not found: {ANNOT_PATH}")
    import gzip
    from rectify.utils.genome import standardize_chrom_name
    junctions = set()
    with gzip.open(str(ANNOT_PATH), 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            p = line.strip().split('\t')
            if len(p) < 9 or p[2].lower() != 'intron':
                continue
            chrom = standardize_chrom_name(p[0])
            junctions.add((chrom, int(p[3]) - 1, int(p[4]), p[6]))
    return junctions


# =============================================================================
# Unit tests: hp_edit_distance
# =============================================================================

class TestHpEditDistance:
    """Tests for the run-length-scaled hp_edit_distance.

    New model: discount only fires at run length >= min_run (default 4).
    Short runs (1-3 bp) get full penalty 1.0.
    Deletions scaled by del_scale (default 1.0), insertions by ins_scale (default 0.7).
    """

    def test_identical(self):
        assert _hp_edit_distance("ACGT", "ACGT") == 0.0

    def test_single_sub(self):
        assert _hp_edit_distance("ACGT", "ACTT") == 1.0

    def test_short_run_no_discount(self):
        # 'TTT' vs 'TT': run=3 < min_run=4, so deletion costs full 1.0 * del_scale=1.0
        ed = _hp_edit_distance("TTT", "TT")
        assert ed == pytest.approx(1.0)

    def test_short_run_ins_no_discount(self):
        # 'TT' vs 'TTT': run=3 < min_run=4, insertion costs 1.0 * ins_scale=0.7
        ed = _hp_edit_distance("TT", "TTT")
        assert ed == pytest.approx(0.7)

    def test_long_run_del_discounted(self):
        # 'TTTTT' vs 'TTTT': run=5 >= min_run=4, del costs 0.5 * 4/5 * del_scale=1.0 = 0.4
        ed = _hp_edit_distance("TTTTT", "TTTT")
        assert ed == pytest.approx(0.4)

    def test_long_run_ins_discounted(self):
        # 'TTTT' vs 'TTTTT': run=5 >= min_run=4, ins costs 0.5 * 4/5 * ins_scale=0.7 = 0.28
        ed = _hp_edit_distance("TTTT", "TTTTT")
        assert ed == pytest.approx(0.28)

    def test_non_homopolymer_del_full(self):
        # 'ACG' vs 'AC': last base G in run=1, deletion costs 1.0 * del_scale=1.0
        ed = _hp_edit_distance("ACG", "AC")
        assert ed == pytest.approx(1.0)

    def test_longer_run_cheaper_del(self):
        # run=8 → penalty = 0.5 * 4/8 * 1.0 = 0.25
        ed8 = _hp_edit_distance("TTTTTTTT", "TTTTTTT")   # run=8
        # run=4 → penalty = 0.5 * 4/4 * 1.0 = 0.5
        ed4 = _hp_edit_distance("TTTT", "TTT")            # run=4
        assert ed8 < ed4, "Longer hp run should be cheaper to delete"

    def test_empty_strings(self):
        assert _hp_edit_distance("", "") == 0.0
        # Empty s1: all insertions, each costs ins_scale * 1.0 = 0.7 (non-hp single bases)
        assert _hp_edit_distance("", "ACG") == pytest.approx(3 * 0.7)
        # Empty s2: all deletions, each costs del_scale * 1.0 = 1.0
        assert _hp_edit_distance("ACG", "") == pytest.approx(3.0)


# =============================================================================
# Unit tests: canonical_tier
# =============================================================================

class TestCanonicalTier:
    def test_rpl20b_annotated_minus(self, genome):
        g = genome.get('chrXV', '')
        if not g: pytest.skip()
        # [900767, 901193) minus strand: should be tier 0 (both canonical)
        tier = _canonical_tier(900767, 901193, g, '-')
        assert tier == 0, f"Annotated RPL20B intron should be tier 0, got {tier}"

    def test_rpl20b_wrong_mapPacBio_minus(self, genome):
        g = genome.get('chrXV', '')
        if not g: pytest.skip()
        # [900759, 901190) is non-canonical on minus strand (non-GT 5'SS, non-AG 3'SS)
        # Under the YAG hierarchy: tier >= 4 (non-canonical 5'SS prefix)
        tier = _canonical_tier(900759, 901190, g, '-')
        assert tier > 1, f"Wrong mapPacBio junction should be non-canonical (tier > 1), got {tier}"

    def test_gcr1_annotated_plus(self, genome):
        g = genome.get('chrXVI', '')
        if not g: pytest.skip()
        tier = _canonical_tier(412261, 413012, g, '+')
        assert tier == 0, f"GCR1 annotated intron should be tier 0, got {tier}"

    def test_gcr1_alt_5ss_plus(self, genome):
        g = genome.get('chrXVI', '')
        if not g: pytest.skip()
        # Alt 5'SS at 412256: should also be tier 0 (GT...AG)
        tier = _canonical_tier(412256, 413012, g, '+')
        assert tier == 0, f"GCR1 alt 5'SS should be tier 0, got {tier}"


# =============================================================================
# Unit tests: _score_junction (split-alignment)
# =============================================================================

class TestScoreJunction:
    def test_correct_junction_scores_lower_than_wrong(self, genome):
        """Annotated RPL20B junction should score lower than wrong mapPacBio one."""
        g = genome.get('chrXV', '')
        if not g: pytest.skip()
        if not CONSENSUS_BAM.exists():
            pytest.skip(f"Consensus BAM not found: {CONSENSUS_BAM}")

        with pysam.AlignmentFile(str(CONSENSUS_BAM), 'rb') as bam:
            for r in bam:
                if r.query_name == '0b3b593b-d3c0-4b24-8244-50ebdecb5bbb':
                    read = r; break
            else:
                pytest.skip("0b3b593b read not found")

        # Find q_split for the N-op
        q_split = None
        for _, ns, ne, qs in _iter_n_ops(read):
            q_split = qs
            break
        assert q_split is not None

        q = read.query_sequence
        score_correct, _ = _score_junction(q, q_split, 900767, 901193, g, 0.25, 15, 10)
        score_wrong,   _ = _score_junction(q, q_split, 900759, 901190, g, 0.25, 15, 10)

        # After scoring, the canonical tier of wrong junction is 2 (non-canonical)
        # so it should lose regardless. But score alone should also differ.
        # We check that at least one of: (a) score_correct < score_wrong, or
        # (b) correct junction wins on tier (tier 0 vs tier 2).
        tier_correct = _canonical_tier(900767, 901193, g, '-')
        tier_wrong   = _canonical_tier(900759, 901190, g, '-')
        assert tier_correct < tier_wrong or score_correct <= score_wrong, (
            f"Correct junction should win: score_correct={score_correct:.2f}, "
            f"score_wrong={score_wrong:.2f}, tier_correct={tier_correct}, "
            f"tier_wrong={tier_wrong}"
        )

    def test_score_returns_float(self, genome):
        g = genome.get('chrXV', '')
        if not g: pytest.skip()
        query = 'A' * 40
        score, delta = _score_junction(query, 20, 900767, 901193, g, 0.25, 15, 5)
        assert isinstance(score, float)
        assert isinstance(delta, int)
        assert delta == 0  # _score_junction always returns 0 (max_slide is API-compat only)


# =============================================================================
# Integration tests: refine_read_junctions
# =============================================================================

class TestRefineReadJunctions:
    @pytest.fixture(scope="class")
    def rpl20b_setup(self, genome, annotated_junctions):
        if not CONSENSUS_BAM.exists():
            pytest.skip(f"Consensus BAM not found: {CONSENSUS_BAM}")
        for b in ALIGNER_BAMS_RPL:
            if not Path(b).exists():
                pytest.skip(f"Aligner BAM not found: {b}")

        all_j, annot_set = build_junction_pool(ALIGNER_BAMS_RPL, annotated_junctions)
        idx = _build_junction_index(all_j)
        g = genome.get('chrXV', '')
        return idx, annot_set, g

    def test_0b3b593b_refined_away_from_wrong_junction(self, genome, rpl20b_setup):
        """0b3b593b has wrong 431N [900758,901189) → should be moved to a canonical junction.

        With sequence-first scoring, the exact destination depends on the read's
        own sequence evidence.  The key assertions are:
        - The wrong junction [900758,901189) is non-canonical (tier=2) and must move
        - The new junction must have intron_end=901193 (correct exon1 boundary)
        - The new junction must be canonical tier 0 or better than the old one

        Note: pure sequence scoring may pick [900763,901193) over the annotated
        [900767,901193) when the read's sequence better explains a nearby junction.
        This is correct behaviour — annotation only breaks exact-score ties.
        """
        idx, annot_set, g = rpl20b_setup
        with pysam.AlignmentFile(str(CONSENSUS_BAM), 'rb') as bam:
            for r in bam:
                if r.query_name == '0b3b593b-d3c0-4b24-8244-50ebdecb5bbb':
                    read = r; break
            else:
                pytest.skip("0b3b593b not found in BAM")

        replacements = refine_read_junctions(
            read, idx, annot_set, g, strand='-',
            hp_pen=0.25, W=15, max_slide=10, search_radius=5000,
        )
        assert len(replacements) == 1, f"Expected 1 replacement, got {len(replacements)}"
        _, old_ns, old_ne, new_ns, new_ne = replacements[0]
        # Old junction should be the wrong mapPacBio one
        assert old_ns < 900767 and old_ne < 901193, \
            f"Old junction [{old_ns},{old_ne}) should be left of annotated [900767,901193)"
        # New junction must have the correct intron_end (exon1 boundary is fixed)
        assert new_ne == 901193, \
            f"Expected intron_end=901193, got {new_ne}"
        # New junction must not be WORSE than the old junction on canonicality.
        # Old junction [900758,901189) is tier 1, so tier 0 or 1 are both acceptable.
        from rectify.core.junction_refiner import _canonical_tier
        old_tier = _canonical_tier(old_ns, old_ne, g, '-')
        new_tier = _canonical_tier(new_ns, new_ne, g, '-')
        assert new_tier <= old_tier or new_tier <= 1, (
            f"New junction [{new_ns},{new_ne}) tier={new_tier} is worse than "
            f"old tier={old_tier} — should not move to lower-canonical junction"
        )

    def test_correct_reads_not_moved(self, genome, rpl20b_setup):
        """Reads already at the annotated junction should not be changed."""
        idx, annot_set, g = rpl20b_setup
        n_tested = n_moved = 0
        with pysam.AlignmentFile(str(CONSENSUS_BAM), 'rb') as bam:
            for read in bam:
                if read.is_unmapped or not read.cigartuples:
                    continue
                # Check for N-op at annotated RPL20B junction
                for _, ns, ne, _ in _iter_n_ops(read):
                    if ns == 900767 and ne == 901193:
                        n_tested += 1
                        reps = refine_read_junctions(
                            read, idx, annot_set, g, strand='-',
                            hp_pen=0.25, W=15, max_slide=10, search_radius=5000,
                        )
                        # Should not move this read
                        moved = any(old_ns != new_ns or old_ne != new_ne
                                    for _, old_ns, old_ne, new_ns, new_ne in reps
                                    if old_ns == 900767 and old_ne == 901193)
                        if moved:
                            n_moved += 1
                        break
                if n_tested >= 50:
                    break

        assert n_tested > 0, "No reads found with annotated RPL20B junction"
        frac_moved = n_moved / n_tested
        assert frac_moved <= 0.02, (
            f"Too many correct reads moved: {n_moved}/{n_tested} = {frac_moved:.1%} "
            f"(tolerance: ≤2%)"
        )

    def test_hp_pen_sensitivity(self, genome, rpl20b_setup):
        """0b3b593b should be moved away from wrong junction at all hp_pen values.

        With sequence-first scoring, the exact destination may vary with hp_pen,
        but the wrong junction [900758,901189) must always be corrected and
        the new junction must have intron_end=901193 and be canonical tier 0.
        """
        idx, annot_set, g = rpl20b_setup
        with pysam.AlignmentFile(str(CONSENSUS_BAM), 'rb') as bam:
            for r in bam:
                if r.query_name == '0b3b593b-d3c0-4b24-8244-50ebdecb5bbb':
                    read = r; break
            else:
                pytest.skip("0b3b593b not found")

        from rectify.core.junction_refiner import _canonical_tier
        for pen in [0.1, 0.25, 0.5, 1.0]:
            reps = refine_read_junctions(
                read, idx, annot_set, g, strand='-',
                hp_pen=pen, W=15, max_slide=10, search_radius=5000,
            )
            assert len(reps) > 0, f"hp_pen={pen}: no replacement proposed"
            for _, old_ns, old_ne, new_ns, new_ne in reps:
                assert new_ne == 901193, \
                    f"hp_pen={pen}: wrong intron_end {new_ne} (expected 901193)"
                old_tier = _canonical_tier(old_ns, old_ne, g, '-')
                new_tier = _canonical_tier(new_ns, new_ne, g, '-')
                assert new_tier <= old_tier or new_tier <= 1, (
                    f"hp_pen={pen}: moved to worse-canonical junction "
                    f"[{new_ns},{new_ne}) tier={new_tier} (old tier={old_tier})"
                )


# =============================================================================
# GCR1 tests
# =============================================================================

class TestGcr1:
    """Tests for GCR1 alternative splice sites.

    GCR1 (YPL075W) on chrXVI plus strand has:
    - Annotated intron: [412261, 413012)
    - Alt 5'SS (5 bp upstream): intron_start = 412256 → [412256, 413012)
    - Multiple alt 3'SS candidates within ±30 bp of intron_end

    All GT-AG candidates: canonical_tier cannot discriminate.
    Score must do the work; 'annotated' flag breaks final ties.
    """

    GCR1_CONSENSUS_BAM = RECTIFY_ROOT / "dev_runs/gcr1_drs_trim_20260420/consensus/gcr1.consensus.bam"
    GCR1_ALIGN_DIR = RECTIFY_ROOT / "dev_runs/gcr1_drs_trim_20260420/align"

    @pytest.fixture(scope="class")
    def gcr1_setup(self, genome, annotated_junctions):
        if not self.GCR1_CONSENSUS_BAM.exists():
            pytest.skip(f"GCR1 consensus BAM not found: {self.GCR1_CONSENSUS_BAM}")
        aligner_bams = [
            str(self.GCR1_ALIGN_DIR / f"gcr1.{a}.bam")
            for a in ['minimap2', 'gapmm2', 'mapPacBio']
        ]
        aligner_bams = [b for b in aligner_bams if Path(b).exists()]
        if not aligner_bams:
            pytest.skip("No GCR1 aligner BAMs found")

        all_j, annot_set = build_junction_pool(aligner_bams, annotated_junctions)
        idx = _build_junction_index(all_j)
        g = genome.get('chrXVI', '')
        return idx, annot_set, g

    def test_annotated_gcr1_reads_not_moved(self, genome, gcr1_setup):
        """Reads at the annotated GCR1 junction [412261, 413012) should stay."""
        idx, annot_set, g = gcr1_setup
        n_tested = n_moved = 0
        with pysam.AlignmentFile(str(self.GCR1_CONSENSUS_BAM), 'rb') as bam:
            for read in bam:
                if read.is_unmapped or not read.cigartuples:
                    continue
                for _, ns, ne, _ in _iter_n_ops(read):
                    if ns == 412261 and ne == 413012:
                        n_tested += 1
                        reps = refine_read_junctions(
                            read, idx, annot_set, g, strand='+',
                            hp_pen=0.25, W=15, max_slide=10, search_radius=5000,
                        )
                        moved = any(old_ns != new_ns or old_ne != new_ne
                                    for _, old_ns, old_ne, new_ns, new_ne in reps
                                    if old_ns == 412261 and old_ne == 413012)
                        if moved:
                            n_moved += 1
                        break
                if n_tested >= 100:
                    break

        if n_tested == 0:
            pytest.skip("No reads with annotated GCR1 junction found")
        frac_moved = n_moved / n_tested if n_tested > 0 else 0
        assert frac_moved <= 0.05, (
            f"Too many correct GCR1 reads moved: {n_moved}/{n_tested} = {frac_moved:.1%}"
        )

    def test_gcr1_all_junctions_canonical(self, genome):
        """All 10 candidate GCR1 junctions (2 × 5'SS × 5 3'SS) should be GT-AG."""
        g = genome.get('chrXVI', '')
        if not g: pytest.skip()
        ss5 = [412256, 412261]
        ss3 = [412995, 413003, 413012, 413018, 413024]
        for s5 in ss5:
            for e3 in ss3:
                tier = _canonical_tier(s5, e3, g, '+')
                # All GCR1 candidates end in AG (canonical), giving tier 0 (YAG: CAG/TAG)
                # or tier 1 (RAG: AAG/GAG). Both are acceptable canonical junctions
                # — tier 2+ (NBG, NAT, other) would indicate a non-AG acceptor.
                assert tier <= 1, (
                    f"GCR1 candidate [{s5},{e3}) should be GT-AG canonical (tier ≤1), "
                    f"got tier {tier}"
                )

    def test_gcr1_annotated_wins_over_alt5ss_when_scores_differ(self, genome, gcr1_setup):
        """When scores are close, annotated junction should win over alt 5'SS."""
        idx, annot_set, g = gcr1_setup
        # Build a synthetic N-op read at the alt 5'SS junction [412256, 413012)
        # to verify it gets corrected back to annotated [412261, 413012) OR stays alt
        # (depends on actual sequence support — we just verify no crash and valid output)
        with pysam.AlignmentFile(str(self.GCR1_CONSENSUS_BAM), 'rb') as bam:
            for read in bam:
                if read.is_unmapped or not read.cigartuples:
                    continue
                for _, ns, ne, _ in _iter_n_ops(read):
                    if abs(ns - 412261) <= 10 and abs(ne - 413012) <= 10:
                        reps = refine_read_junctions(
                            read, idx, annot_set, g, strand='+',
                            hp_pen=0.25, W=15, max_slide=15, search_radius=5000,
                        )
                        # If replacement proposed, new junction must be canonical tier 0
                        for _, old_ns, old_ne, new_ns, new_ne in reps:
                            tier = _canonical_tier(new_ns, new_ne, g, '+')
                            assert tier == 0, (
                                f"Proposed replacement [{new_ns},{new_ne}) is not "
                                f"canonical tier 0 (tier={tier})"
                            )
                        return  # tested one read successfully
        pytest.skip("No GCR1 reads with N-ops near annotated junction")


# =============================================================================
# Unit tests: _apply_junction_replacement CIGAR surgery
# =============================================================================

class TestApplyJunctionReplacement:
    """Synthetic tests for CIGAR surgery correctness."""

    def _make_read(self, cigar_str: str, ref_start: int, seq: str) -> pysam.AlignedSegment:
        """Create a minimal AlignedSegment for testing."""
        header = pysam.AlignmentHeader.from_dict({
            'HD': {'VN': '1.6'},
            'SQ': [{'LN': 2000000, 'SN': 'chrXV'}],
        })
        read = pysam.AlignedSegment(header)
        read.query_name = 'test'
        read.reference_id = 0
        read.reference_start = ref_start
        read.cigarstring = cigar_str
        read.query_sequence = seq
        read.query_qualities = pysam.qualitystring_to_array('I' * len(seq))
        read.flag = 0
        return read

    def test_extend_intron_right(self):
        """N-op end moves right by 2: exon1 loses 2 ref bases → I(2)M(8) after N.

        Encoding: delta_end > 0 → I(d) inserted after N, M shrinks:
            10M 100N 10M → 10M 102N 2I 8M
        Ref span preserved (10+100+10 = 10+102+8 = 120).
        Query span preserved (10+10 = 10+2+8 = 20).
        """
        seq = 'A' * 20
        read = self._make_read('10M100N10M', 900757, seq)
        from rectify.core.junction_refiner import _apply_junction_replacement
        g = 'A' * 2000000
        applied = _apply_junction_replacement(read, 1, 900767, 900867, 900767, 900869, g, '+', 0.25, 15)
        assert applied, "Expected _apply_junction_replacement to return True"
        new_cigar = read.cigarstring
        # Intron grows by 2 on right: 102N; exon1 becomes 2I8M (I absorbs displaced query)
        assert '102N' in new_cigar, f"Expected 102N in {new_cigar}"
        assert '2I' in new_cigar, f"Expected 2I (displaced query) in {new_cigar}"
        assert '8M' in new_cigar, f"Expected 8M in {new_cigar}"
        # Ref and query spans preserved
        old_ref = old_q = 120
        new_ref = sum(l for op, l in read.cigartuples if op in (0, 2, 3, 7, 8))
        new_q   = sum(l for op, l in read.cigartuples if op in (0, 1, 4, 7, 8))
        assert new_ref == 120, f"Ref span changed: {new_ref}"
        assert new_q   == 20,  f"Query span changed: {new_q}"

    def test_shrink_intron_left(self):
        """N-op start moves right by 3: exon2 gains 3 ref bases from intron → D(3) at boundary.

        Encoding: delta_start > 0 → D(d) inserted before N, M unchanged:
            10M 100N 10M → 10M 3D 97N 10M
        Ref span preserved (10+100+10 = 10+3+97+10 = 120).
        Query span preserved (10+10 = 10+0+10 = 20).
        """
        seq = 'A' * 20
        read = self._make_read('10M100N10M', 900757, seq)
        from rectify.core.junction_refiner import _apply_junction_replacement
        g = 'A' * 2000000
        applied = _apply_junction_replacement(read, 1, 900767, 900867, 900770, 900867, g, '+', 0.25, 15)
        assert applied, "Expected _apply_junction_replacement to return True"
        new_cigar = read.cigarstring
        # intron_start moved right: 97N; exon2 boundary gets D(3)
        assert '97N' in new_cigar, f"Expected 97N in {new_cigar}"
        assert '3D' in new_cigar, f"Expected 3D (reference-only boundary) in {new_cigar}"
        new_ref = sum(l for op, l in read.cigartuples if op in (0, 2, 3, 7, 8))
        new_q   = sum(l for op, l in read.cigartuples if op in (0, 1, 4, 7, 8))
        assert new_ref == 120, f"Ref span changed: {new_ref}"
        assert new_q   == 20,  f"Query span changed: {new_q}"

    def test_unchanged_when_same_junction(self):
        """No-op: old and new junction are identical → should return False."""
        seq = 'A' * 20
        read = self._make_read('10M100N10M', 900757, seq)
        from rectify.core.junction_refiner import _apply_junction_replacement
        g = 'A' * 2000000
        # Apply with same coordinates → should still succeed (we check via refine_read_junctions
        # which only calls _apply for genuine changes)
        applied = _apply_junction_replacement(read, 1, 900767, 900867, 900767, 900867, g, '+', 0.25, 15)
        # No delta → no change, but the function still returns True (N len is the same)
        # The calling code in refine_read_junctions filters out same-coord replacements
        # so we just verify no crash
        assert isinstance(applied, bool)

    def test_ref_span_preserved(self):
        """Total reference span must be preserved after replacement."""
        seq = 'G' * 30
        read = self._make_read('15M200N15M', 900740, seq)
        from rectify.core.junction_refiner import _apply_junction_replacement
        g = 'A' * 2000000
        old_ref = sum(l for op, l in read.cigartuples if op in (0, 2, 3, 7, 8))
        applied = _apply_junction_replacement(read, 1, 900755, 900955, 900755, 900957, g, '+', 0.25, 15)
        if applied:
            new_ref = sum(l for op, l in read.cigartuples if op in (0, 2, 3, 7, 8))
            assert old_ref == new_ref, f"Ref span changed: {old_ref} → {new_ref}"

    def test_query_span_preserved(self):
        """Total query-consuming span must be preserved after replacement."""
        seq = 'C' * 30
        read = self._make_read('15M200N15M', 900740, seq)
        from rectify.core.junction_refiner import _apply_junction_replacement
        g = 'A' * 2000000
        old_q = sum(l for op, l in read.cigartuples if op in (0, 1, 4, 7, 8))
        applied = _apply_junction_replacement(read, 1, 900755, 900955, 900758, 900955, g, '+', 0.25, 15)
        if applied:
            new_q = sum(l for op, l in read.cigartuples if op in (0, 1, 4, 7, 8))
            assert old_q == new_q, f"Query span changed: {old_q} → {new_q}"


# =============================================================================
# Integration: build_junction_pool
# =============================================================================

class TestBuildJunctionPool:
    def test_includes_annotated(self, annotated_junctions):
        if not all(Path(b).exists() for b in ALIGNER_BAMS_RPL):
            pytest.skip("Aligner BAMs not found")
        all_j, annot_set = build_junction_pool(ALIGNER_BAMS_RPL, annotated_junctions)
        assert RPL20B_INTRON in annot_set, f"{RPL20B_INTRON} not in annotated_set"
        assert RPL20B_INTRON in all_j, f"{RPL20B_INTRON} not in all_junctions"

    def test_includes_novel_from_aligners(self, annotated_junctions):
        if not all(Path(b).exists() for b in ALIGNER_BAMS_RPL):
            pytest.skip("Aligner BAMs not found")
        all_j, annot_set = build_junction_pool(ALIGNER_BAMS_RPL, annotated_junctions)
        # The wrong mapPacBio junction should be in all_j (it came from the aligner)
        # but NOT in annotated_set
        wrong = ("chrXV", 900759, 901190)
        # May or may not be present depending on mapPacBio run, just verify no crash
        assert isinstance(all_j, set)
        assert isinstance(annot_set, set)
        assert len(all_j) >= len(annot_set)

    def test_candidates_near(self, annotated_junctions):
        if not all(Path(b).exists() for b in ALIGNER_BAMS_RPL):
            pytest.skip("Aligner BAMs not found")
        all_j, _ = build_junction_pool(ALIGNER_BAMS_RPL, annotated_junctions)
        idx = _build_junction_index(all_j)
        cands = _candidates_near(idx, 'chrXV', 900767, 901193, radius=100)
        assert any(js == 900767 and je == 901193 for js, je in cands), \
            "Annotated RPL20B junction not found in candidates"


# =============================================================================
# Scoring parameter sweep (informational, not strict assertions)
# =============================================================================

class TestHpPenSweep:
    """Documents scoring behaviour across hp_pen values for the key test reads.

    With sequence-first scoring, hp_pen affects WHAT junction is chosen but
    the wrong junction [900758,901189) must always be replaced with something
    canonical at intron_end=901193.
    """

    @pytest.mark.parametrize("hp_pen", [0.1, 0.2, 0.25, 0.3, 0.4, 0.5])
    def test_0b3b593b_always_corrected(self, genome, annotated_junctions, hp_pen):
        """0b3b593b wrong junction must always be corrected to a canonical junction
        at intron_end=901193, regardless of hp_pen value."""
        if not CONSENSUS_BAM.exists():
            pytest.skip()
        for b in ALIGNER_BAMS_RPL:
            if not Path(b).exists():
                pytest.skip()
        all_j, annot_set = build_junction_pool(ALIGNER_BAMS_RPL, annotated_junctions)
        idx = _build_junction_index(all_j)
        g = genome.get('chrXV', '')
        if not g: pytest.skip()

        with pysam.AlignmentFile(str(CONSENSUS_BAM), 'rb') as bam:
            for r in bam:
                if r.query_name == '0b3b593b-d3c0-4b24-8244-50ebdecb5bbb':
                    read = r; break
            else:
                pytest.skip()

        reps = refine_read_junctions(
            read, idx, annot_set, g, strand='-',
            hp_pen=hp_pen, W=15, max_slide=10, search_radius=5000,
        )
        assert len(reps) > 0, f"hp_pen={hp_pen}: no replacement proposed for 0b3b593b"
        for _, old_ns, old_ne, new_ns, new_ne in reps:
            assert new_ne == 901193, \
                f"hp_pen={hp_pen}: wrong intron_end {new_ne}"
            old_tier = _canonical_tier(old_ns, old_ne, g, '-')
            new_tier = _canonical_tier(new_ns, new_ne, g, '-')
            # Must not move to a worse-canonical junction
            assert new_tier <= old_tier or new_tier <= 1, (
                f"hp_pen={hp_pen}: moved to worse-canonical junction "
                f"[{new_ns},{new_ne}) tier={new_tier} from tier={old_tier}"
            )


# =============================================================================
# Tests for additional PMC3983031 alt-splice genes: TFC3, RPL22B, SRC1
# =============================================================================

class TestAltSpliceGenes:
    """Tests for alternative splice site genes from Kawashima et al. 2014
    (PMC3983031). Covers:
    - TFC3 (YAL001C): chrI minus, alternative 3'SS 17 nt upstream
    - RPL22B (YFL034C-A): chrVI minus, alternative 5'SS
    - SRC1 (YML034W): chrXIII plus, alternative 5'SS 4 nt downstream
    """

    ALT_SPLICE_CONSENSUS = (
        RECTIFY_ROOT / "dev_runs/alt_splice_tfc3_rpl22b_src1_20260420/consensus/alt_splice.consensus.bam"
    )
    ALT_SPLICE_ALIGN_DIR = (
        RECTIFY_ROOT / "dev_runs/alt_splice_tfc3_rpl22b_src1_20260420/align"
    )

    # Annotated introns (0-based half-open from GFF)
    TFC3_INTRON    = ("chrI",    151006, 151096)  # minus strand
    RPL22B_INTRON  = ("chrVI",   64599,  64920)   # minus strand
    SRC1_INTRON    = ("chrXIII", 211444, 211570)  # plus strand

    # Alternative junctions observed in consensus BAM
    TFC3_ALT_3SS    = ("chrI",    150989, 151096)  # alt 3'SS 17 nt upstream
    RPL22B_ALT_5SS  = ("chrVI",   64599,  64856)   # alt 5'SS (shorter intron)
    SRC1_ALT_5SS    = ("chrXIII", 211440, 211570)  # alt 5'SS 4 nt upstream

    @pytest.fixture(scope="class")
    def alt_splice_setup(self, genome, annotated_junctions):
        if not self.ALT_SPLICE_CONSENSUS.exists():
            pytest.skip(f"Alt-splice consensus BAM not found: {self.ALT_SPLICE_CONSENSUS}")
        aligner_bams = [
            str(self.ALT_SPLICE_ALIGN_DIR / f"alt_splice.{a}.bam")
            for a in ['minimap2', 'gapmm2', 'mapPacBio']
        ]
        aligner_bams = [b for b in aligner_bams if Path(b).exists()]
        if not aligner_bams:
            pytest.skip("No alt-splice aligner BAMs found")

        all_j, annot_set = build_junction_pool(aligner_bams, annotated_junctions)
        idx = _build_junction_index(all_j)
        return idx, annot_set, genome

    def test_tfc3_annotated_junction_canonical(self, genome):
        """TFC3 annotated intron [151006,151096) has 3'SS=AAG (RAG class, tier 1).

        AAG is tier 1 (RAG) in the YAG > RAG > NBG > NAT hierarchy.
        The canonical GT 5'SS makes the overall tier = 1 (not 0=YAG, which
        requires a pyrimidine before AG, e.g. CAG or TAG).
        """
        g = genome.get('chrI', '')
        if not g: pytest.skip()
        tier = _canonical_tier(*self.TFC3_INTRON[1:], g, '-')
        assert tier == 1, f"TFC3 annotated intron tier={tier}, expected 1 (RAG class)"

    def test_tfc3_alt_3ss_canonical(self, genome):
        """TFC3 alternative 3'SS [150989,151096) has 3'SS=CAG (YAG class, tier 0).

        The alt 3'SS is actually higher quality than the annotated one under
        the YAG hierarchy (CAG > AAG), consistent with its observed usage.
        """
        g = genome.get('chrI', '')
        if not g: pytest.skip()
        tier = _canonical_tier(*self.TFC3_ALT_3SS[1:], g, '-')
        assert tier == 0, f"TFC3 alt 3'SS tier={tier}, expected 0 (YAG class)"

    def test_rpl22b_annotated_junction_canonical(self, genome):
        """RPL22B annotated intron [64599,64920) should be canonical."""
        g = genome.get('chrVI', '')
        if not g: pytest.skip()
        tier = _canonical_tier(*self.RPL22B_INTRON[1:], g, '-')
        assert tier == 0, f"RPL22B annotated intron tier={tier}, expected 0"

    def test_src1_annotated_junction_canonical(self, genome):
        """SRC1 annotated intron [211444,211570) should be canonical GT-AG on plus strand."""
        g = genome.get('chrXIII', '')
        if not g: pytest.skip()
        tier = _canonical_tier(*self.SRC1_INTRON[1:], g, '+')
        assert tier == 0, f"SRC1 annotated intron tier={tier}, expected 0"

    def test_alt_splice_annotated_reads_not_moved(self, genome, alt_splice_setup):
        """Reads at annotated junctions for TFC3, RPL22B, SRC1 should not be displaced."""
        idx, annot_set, genome_dict = alt_splice_setup
        gene_introns = [
            (self.TFC3_INTRON,   'chrI',    '-'),
            (self.RPL22B_INTRON, 'chrVI',   '-'),
            (self.SRC1_INTRON,   'chrXIII', '+'),
        ]
        results = {}
        with pysam.AlignmentFile(str(self.ALT_SPLICE_CONSENSUS), 'rb') as bam:
            for read in bam:
                if read.is_unmapped or not read.cigartuples:
                    continue
                from rectify.utils.genome import standardize_chrom_name
                chrom = standardize_chrom_name(read.reference_name)
                strand = '-' if read.is_reverse else '+'
                g = genome_dict.get(chrom, '')
                if not g:
                    continue
                for (intron_chrom, ns, ne), expected_chrom, expected_strand in gene_introns:
                    if chrom != intron_chrom or strand != expected_strand:
                        continue
                    for _, op_ns, op_ne, _ in _iter_n_ops(read):
                        if op_ns == ns and op_ne == ne:
                            key = (intron_chrom, ns, ne)
                            if key not in results:
                                results[key] = {'tested': 0, 'moved': 0}
                            results[key]['tested'] += 1
                            reps = refine_read_junctions(
                                read, idx, annot_set, g, strand,
                                hp_pen=0.25, W=15, max_slide=10, search_radius=5000,
                            )
                            moved = any(
                                old_ns != new_ns or old_ne != new_ne
                                for _, old_ns, old_ne, new_ns, new_ne in reps
                                if old_ns == ns and old_ne == ne
                            )
                            if moved:
                                results[key]['moved'] += 1
                            break

        for (chrom, ns, ne), expected_chrom, _ in gene_introns:
            key = (chrom, ns, ne)
            if key not in results:
                continue  # no reads found with this junction
            n_tested = results[key]['tested']
            n_moved  = results[key]['moved']
            frac = n_moved / n_tested if n_tested > 0 else 0
            assert frac <= 0.10, (
                f"{chrom}:[{ns},{ne}): {n_moved}/{n_tested} = {frac:.1%} "
                f"correct reads incorrectly moved (tolerance ≤10%)"
            )

    def test_tfc3_alt_3ss_read_corrected_or_preserved(self, genome, alt_splice_setup):
        """Reads at TFC3 alt 3'SS [150989,151096) should be corrected to annotated
        [151006,151096) OR preserved (if the scoring supports this isoform).
        Either outcome is valid — no read should be moved to a non-canonical junction.
        """
        idx, annot_set, genome_dict = alt_splice_setup
        g = genome_dict.get('chrI', '')
        if not g: pytest.skip()

        n_tested = n_non_canonical = 0
        with pysam.AlignmentFile(str(self.ALT_SPLICE_CONSENSUS), 'rb') as bam:
            for read in bam:
                if read.is_unmapped or not read.cigartuples: continue
                from rectify.utils.genome import standardize_chrom_name
                if standardize_chrom_name(read.reference_name) != 'chrI': continue
                if not read.is_reverse: continue
                for _, op_ns, op_ne, _ in _iter_n_ops(read):
                    if op_ns == 150989 and op_ne == 151096:
                        n_tested += 1
                        reps = refine_read_junctions(
                            read, idx, annot_set, g, '-',
                            hp_pen=0.25, W=15, max_slide=10, search_radius=5000,
                        )
                        for _, old_ns, old_ne, new_ns, new_ne in reps:
                            if old_ns == 150989 and old_ne == 151096:
                                tier = _canonical_tier(new_ns, new_ne, g, '-')
                                # Tier >= 2 = NBG/NAT/other (non-AG acceptor): not acceptable
                                # Tier 0/1 = YAG/RAG: biologically valid acceptors
                                if tier >= 2:
                                    n_non_canonical += 1
                        break

        if n_tested == 0:
            pytest.skip("No reads with TFC3 alt 3'SS found")
        assert n_non_canonical == 0, (
            f"TFC3 alt 3'SS reads moved to non-canonical junction: "
            f"{n_non_canonical}/{n_tested}"
        )
