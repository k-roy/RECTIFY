"""ISSUE-002: the 5' rescue writer must draw its N-op where the rescue meant it.

`extend_read_5prime_for_junction_rescue` derives intron_len from the LIVE
alignment edge and never reads `five_prime_intron_clip_pos`. For a read whose 5'
end lies inside the rescued intron those differ, so extend fabricates an intron
running to the read's OLD 5' edge instead of to the annotated acceptor the TSV
names. The gate is "extend's N-op would land on icp", not "icp < 0" — the latter
regresses the reads that reach the right junction via `five_prime_upstream_trim`.

Run with:
    pytest tests/test_2f_rescue_writer_routing.py -v
"""

import os
from pathlib import Path
from typing import List, Tuple

import pysam
import pytest

from rectify.core.bam.bam_writer import (
    REFUSAL_EXTEND,
    REFUSAL_NONCANONICAL,
    REFUSAL_REROUTE,
    REFUSAL_SOFTCLIP_ONLY,
    _n_op_intervals,
    apply_5prime_rescue_surgery,
    apply_corrected_edits_to_read,
    predict_5prime_rescue_refusal,
)
from rectify.core.bam.read_edits import (
    extend_read_5prime_for_junction_rescue,
    projected_5prime_rescue_intron_edge,
)
from rectify.core.splice.splice_aware_5prime import (
    _get_intronic_query_bases,
    _intronic_bases_favour_intron,
)


# exon1 [0, 100)  |  intron [100, 200) GT..AG  |  exon2 [200, 300)
GENOME_SEQ = (
    'ACGTTGCATGCAGTCCATGACGTTGCATGCAGTCCATGACGTTGCATGCAGTCCATG'
    'ACGTTGCATGCAGTCCATGACGTTGCATGCAGTCCATGAC'          # pos 0-96
    + 'TAG'                                              # pos 97-99 (exon1 tail)
    + 'GT' + 'C' * 96 + 'AG'                             # pos 100-199 (intron)
    + 'TTGCACCGATTGCACCGAT' * 5 + 'TTGCA'                # pos 200-299 (exon2)
)
GENOME = {'chrR': GENOME_SEQ}
assert len(GENOME_SEQ) == 300
assert GENOME_SEQ[100:102] == 'GT' and GENOME_SEQ[198:200] == 'AG'


def _make_read(start: int, cigar: List[Tuple[int, int]], strand: str,
               seq: str, name: str = 'r') -> pysam.AlignedSegment:
    hdr = pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.6'},
        'SQ': [{'SN': 'chrR', 'LN': 300}],
    })
    read = pysam.AlignedSegment(hdr)
    read.query_name = name
    read.reference_name = 'chrR'
    read.reference_start = start
    read.cigartuples = cigar
    read.is_reverse = (strand == '-')
    read.is_unmapped = False
    read.is_secondary = False
    read.is_supplementary = False
    read.mapping_quality = 60
    read.query_sequence = seq
    return read


def _intronic_5prime_read():
    """Plus-strand read whose 5' edge is 50 bp INSIDE the intron, with a 12 nt
    5' soft clip: `12S 100M` at 150 spans ref [150, 250)."""
    seq = GENOME_SEQ[88:100] + GENOME_SEQ[150:250]
    return _make_read(150, [(4, 12), (0, 100)], '+', seq, name='intronic5p')


def _correction(**kw):
    base = {
        'five_prime_rescued': True,
        'five_prime_position': 99,        # last base of exon 1 (intron_start - 1)
        'five_prime_soft_clip': 12,
        'five_prime_exon_cigar': '',
        'five_prime_upstream_trim': 0,
        'five_prime_intron_clip_pos': -1,
        'reanchor_clip_len': 0,
        'strand': '+',
        'corrected_3prime': 249,
    }
    base.update(kw)
    return base


# ---------------------------------------------------------------------------
# The predicate
# ---------------------------------------------------------------------------

class TestProjectedIntronEdge:
    def test_plus_edge_is_reference_start(self):
        read = _intronic_5prime_read()
        assert projected_5prime_rescue_intron_edge(read, 12, '+', 0) == 150

    def test_plus_edge_includes_an_applicable_upstream_trim(self):
        read = _intronic_5prime_read()
        assert projected_5prime_rescue_intron_edge(read, 12, '+', 3) == 153

    def test_upstream_trim_is_ignored_when_extend_could_not_apply_it(self):
        """Mirrors extend's own guard: the body op must be M/=/X and long enough."""
        read = _make_read(150, [(4, 12), (2, 5), (0, 95)], '+',
                          GENOME_SEQ[88:100] + GENOME_SEQ[155:250])
        assert projected_5prime_rescue_intron_edge(read, 12, '+', 3) == 150

    def test_minus_edge_is_reference_end(self):
        read = _make_read(150, [(0, 100), (4, 12)], '-',
                          GENOME_SEQ[150:250] + GENOME_SEQ[88:100])
        assert projected_5prime_rescue_intron_edge(read, 12, '-', 0) == 250
        assert projected_5prime_rescue_intron_edge(read, 12, '-', 4) == 246

    def test_none_when_extend_would_not_run(self):
        read = _make_read(150, [(0, 100)], '+', GENOME_SEQ[150:250])
        assert projected_5prime_rescue_intron_edge(read, 0, '+', 0) is None
        assert projected_5prime_rescue_intron_edge(read, 12, '+', 0) is None


# ---------------------------------------------------------------------------
# Routing
# ---------------------------------------------------------------------------

class TestWriterRouting:
    def test_extend_alone_would_fabricate_an_intron_to_the_old_edge(self):
        """The defect, demonstrated on the helper in isolation."""
        read = _intronic_5prime_read()
        assert extend_read_5prime_for_junction_rescue(
            read, 99, 12, '+', exon_cigar_str='12M') is True
        # N runs exon-1 -> the read's ORIGINAL start (150), not the acceptor (200).
        assert _n_op_intervals(read) == ((100, 150),)

    def test_writer_routes_to_reroute_and_lands_on_the_annotated_acceptor(self):
        read = _intronic_5prime_read()
        # 12 clip bases + 50 body bases mapped inside [150, 200)
        assert len(_get_intronic_query_bases(read, 200, '+')) == 62
        corr = _correction(five_prime_intron_clip_pos=200,
                           five_prime_exon_cigar='62M')
        assert apply_corrected_edits_to_read(read, corr, GENOME) is True
        assert _n_op_intervals(read) == ((100, 200),)
        assert read.reference_start == 38          # 99 - 62 + 1
        assert read.cigarstring == '62M100N50M'

    def test_extend_still_runs_when_its_n_op_lands_on_icp(self):
        """The `icp < 0` gate agent2 tested would have regressed this shape."""
        read = _make_read(200, [(4, 12), (0, 100)], '+',
                          GENOME_SEQ[88:100] + GENOME_SEQ[200:300])
        corr = _correction(five_prime_intron_clip_pos=200,  # == reference_start
                           five_prime_exon_cigar='12M',
                           corrected_3prime=299)
        assert apply_corrected_edits_to_read(read, corr, GENOME) is True
        assert _n_op_intervals(read) == ((100, 200),)
        assert read.reference_start == 88

    def test_extend_still_runs_when_the_rescue_published_no_icp(self):
        read = _make_read(200, [(4, 12), (0, 100)], '+',
                          GENOME_SEQ[88:100] + GENOME_SEQ[200:300])
        corr = _correction(five_prime_exon_cigar='12M', corrected_3prime=299)
        assert apply_corrected_edits_to_read(read, corr, GENOME) is True
        assert _n_op_intervals(read) == ((100, 200),)

    def test_no_fabricated_intron_when_both_helpers_refuse(self):
        """Conservative fallback: an exon CIGAR whose query span cannot match the
        intronic run leaves the read alone rather than drawing a junction."""
        read = _intronic_5prime_read()
        before = read.cigartuples
        corr = _correction(five_prime_intron_clip_pos=200,
                           five_prime_exon_cigar='7M')   # wrong query span
        apply_corrected_edits_to_read(read, corr, GENOME)
        # softclip_intronic_tail_5prime may hide the intronic bases, but under no
        # circumstance may an intron to the OLD edge appear.
        assert (100, 150) not in _n_op_intervals(read)
        for start, end in _n_op_intervals(read):
            assert (start, end) not in ((100, 150),), before


# ---------------------------------------------------------------------------
# 006c — the Case-4 unspliced guard, shared with the sequence-rescue branch
# ---------------------------------------------------------------------------

class TestUnsplicedGuardShared:
    def _unspliced_read(self):
        """Read continuing straight into the intron: its intronic bases ARE the
        intron sequence, so no rescue should splice it by fiat."""
        seq = GENOME_SEQ[130:250]
        return _make_read(130, [(0, 120)], '+', seq, name='unspliced')

    def _spliced_read(self):
        """Read whose intron-mapped bases are really exon-1 sequence."""
        seq = GENOME_SEQ[80:100] + GENOME_SEQ[200:300]
        return _make_read(180, [(0, 120)], '+', seq, name='spliced')

    def test_unspliced_read_is_flagged(self):
        read = self._unspliced_read()
        assert _intronic_bases_favour_intron(
            read, GENOME_SEQ, 100, 200, 130, '+') is True

    def test_spliced_read_is_not_flagged(self):
        read = self._spliced_read()
        assert _intronic_bases_favour_intron(
            read, GENOME_SEQ, 100, 200, 180, '+') is False

    def test_strict_mode_does_not_veto_on_a_tie(self):
        """Cases 1/2 arrive with a scored, canonical-donor-preferred sequence
        match; a tie in the intronic bases decides nothing and must not discard
        it. Case 4 (no sequence evidence of its own) keeps the tie-favours-
        unspliced rule."""
        # Intron continuation and exon-1 tail identical -> the bases decide nothing.
        tie_genome = 'A' * 300
        read = _make_read(150, [(0, 100)], '+', 'A' * 100, name='tie')
        assert _intronic_bases_favour_intron(
            read, tie_genome, 100, 200, 150, '+', strict=False) is True
        assert _intronic_bases_favour_intron(
            read, tie_genome, 100, 200, 150, '+', strict=True) is False

    def test_none_when_nothing_to_compare(self):
        read = _make_read(200, [(0, 100)], '+', GENOME_SEQ[200:300])
        assert _intronic_bases_favour_intron(
            read, GENOME_SEQ, 100, 200, 200, '+') is None


# ---------------------------------------------------------------------------
# Human RNA004 panel (git-ignored inputs -> skipif)
# ---------------------------------------------------------------------------

_PANEL_DIR = Path(os.environ.get(
    'RECTIFY_SUMNER_PANEL_DIR',
    '/Users/kevinroy/work/rectify/dev/sumner_misplaced_panel_20260904'))
_PANEL_BAM = _PANEL_DIR / '175_panel_orig.bam'
_PANEL_FA = _PANEL_DIR / 'ref' / 'chr5.fa'
_PANEL_GTF = _PANEL_DIR / 'ref' / 'gencode.v48.basic.chr5.gtf'

_HAVE_PANEL = all(p.exists() for p in (_PANEL_BAM, _PANEL_FA, _PANEL_GTF))


@pytest.mark.skipif(not _HAVE_PANEL,
                    reason='Sumner human panel inputs are git-ignored; set '
                           'RECTIFY_SUMNER_PANEL_DIR to run')
class TestSumnerPanelTsvBamAgreement:
    """Baselines recorded 2026-09-05 on 4eefd1f: TSV<->BAM agreement 22/38
    overall and 2/17 on the softclip branch. After ISSUE-006 + ISSUE-002:
    23/27 and 2/2."""

    @pytest.fixture(scope='class')
    def replay(self):
        import copy

        import rectify.core.splice.splice_aware_5prime as s5
        from rectify.core.consensus.consensus import load_annotated_junctions

        saved = s5.standardize_chrom_name
        s5.standardize_chrom_name = lambda c: c
        try:
            fa = pysam.FastaFile(str(_PANEL_FA))
            chrom = fa.references[0]
            genome = {chrom: fa.fetch(chrom)}
            cand = {(chrom,) + tuple(j[1:])
                    for j in load_annotated_junctions(str(_PANEL_GTF))}
            out = []
            for read in pysam.AlignmentFile(str(_PANEL_BAM)):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                strand = '-' if read.is_reverse else '+'
                ct = read.cigartuples
                clip = ((ct[-1][1] if ct[-1][0] == 4 else 0) if strand == '-'
                        else (ct[0][1] if ct[0][0] == 4 else 0))
                res = s5.rescue_3ss_truncation(copy.deepcopy(read), genome,
                                               cand, strand=strand)
                icp = -1
                rj = res.get('rescued_junction')
                if res['rescued'] and rj:
                    _c, i_s, i_e = rj
                    a5 = ((read.reference_end - 1) if strand == '-'
                          else read.reference_start)
                    if i_s <= a5 < i_e:
                        icp = i_s if strand == '-' else i_e
                corr = {
                    'five_prime_rescued': bool(res['rescued']),
                    'five_prime_position': res['five_prime_corrected'],
                    'five_prime_soft_clip': clip,
                    'five_prime_exon_cigar': res.get('five_prime_exon_cigar', ''),
                    'five_prime_upstream_trim': res.get('five_prime_upstream_trim', 0),
                    'reanchor_clip_len': res.get('reanchor_clip_len', 0),
                    'five_prime_intron_clip_pos': icp,
                    'strand': strand,
                    'corrected_3prime': ((read.reference_end - 1) if strand == '+'
                                         else read.reference_start),
                }
                before = _n_op_intervals(read)
                edited = copy.deepcopy(read)
                apply_corrected_edits_to_read(edited, corr, genome)
                out.append({
                    'read_id': read.query_name,
                    'rescued': bool(res['rescued']),
                    'rescue_type': res['rescue_type'],
                    'tsv_j': (rj[1], rj[2]) if rj else None,
                    'strand': strand,
                    'icp': icp,
                    'before': before,
                    'after': _n_op_intervals(edited),
                })
            return out
        finally:
            s5.standardize_chrom_name = saved

    def test_a_rescue_either_draws_its_own_junction_or_draws_nothing(self, replay):
        """The ISSUE-002 invariant. A rescued read may end up unmodified (both
        surgery helpers refused — the conservative outcome), but it must never
        gain an N-op that is NOT the junction the rescue chose.

        Excluded here: an N-op the reroute produced by MERGING with a pre-existing
        N it had partially trimmed. That is ISSUE-007, whose own test asserts the
        merge is refused; before that fix the merged N inherits the aligner's
        boundary (panel row 14b2473a: chose 133993018-134004894, wrote
        133993015-...).
        """
        offenders = []
        for row in replay:
            if not row['rescued']:
                continue
            lost_edges = {e for n in row['before'] if n not in row['after'] for e in n}
            for n in row['after']:
                if n in row['before'] or n == row['tsv_j']:
                    continue
                if n[0] in lost_edges or n[1] in lost_edges:
                    continue   # merged with a trimmed pre-existing N -> ISSUE-007
                offenders.append((row['read_id'][:12], row['rescue_type'],
                                  'wrote', n, 'chose', row['tsv_j']))
        assert not offenders, (
            'writer drew a junction the rescue did not choose:\n  '
            + '\n  '.join(map(str, offenders)))

    def test_softclip_branch_agreement_is_total(self, replay):
        sc = [r for r in replay if r['rescued'] and r['rescue_type'] == 'softclip']
        agree = [r for r in sc if r['tsv_j'] in r['after']]
        assert len(agree) == len(sc), (
            f'softclip-branch TSV<->BAM agreement {len(agree)}/{len(sc)} '
            '(2/17 before ISSUE-002)')

    def test_every_disagreement_is_a_refusal_or_a_merge(self, replay):
        """A rescued row whose TSV junction is absent from the BAM may only be:

        * a REFUSAL — both surgery helpers declined and the read is untouched
          (the conservative outcome Fix 1 introduces: no junction is better than
          a fabricated one). This leaves a known residual, that the TSV still
          reports a rescue the BAM does not carry; the TSV is written by
          bam_processor before the writer runs, so closing it is a separate
          change.
        * a MERGE — the reroute fused the new N with a pre-existing one it had
          partially trimmed, so the written N inherits the aligner's far edge.
          That is ISSUE-007; its own test asserts the merge is refused.

        Anything else is a fabricated junction — the ISSUE-002 defect.
        """
        offenders = []
        for row in replay:
            if not row['rescued'] or row['tsv_j'] in row['after']:
                continue
            new = [n for n in row['after'] if n not in row['before']]
            if not new:
                continue                      # refusal
            lost_edges = {e for n in row['before'] if n not in row['after'] for e in n}
            if all(n[0] in lost_edges or n[1] in lost_edges for n in new):
                continue                      # merge (ISSUE-007)
            offenders.append((row['read_id'][:12], row['rescue_type'],
                              row['tsv_j'], new))
        assert not offenders, offenders

    def test_interior_n_ops_are_never_displaced(self, replay):
        """The plus-strand reroute recomputes reference_start; if exon_ref_span
        were off by one, EVERY interior N-op would shift and the canonical guard
        would revert whole reads for the wrong reason.

        An N-op inside the region the reroute legitimately trims away (above the
        clip boundary on `-`, below it on `+`) is expected to disappear; every
        other pre-existing N-op must survive at identical coordinates.
        """
        offenders = []
        for row in replay:
            if not row['rescued']:
                continue
            after = set(row['after'])
            for n in row['before']:
                if n in after:
                    continue
                if row['icp'] >= 0 and (
                        (row['strand'] == '-' and n[0] >= row['icp'])
                        or (row['strand'] == '+' and n[1] <= row['icp'])):
                    continue   # inside the trimmed intronic tail
                # A terminal N whose far edge was extended onto the rescued
                # acceptor keeps its near edge; that is the surgery, not drift.
                if any(m[0] == n[0] or m[1] == n[1] for m in after - set(row['before'])):
                    continue
                offenders.append((row['read_id'][:12], row['strand'], n,
                                  row['icp'], row['before'], row['after']))
        assert not offenders, offenders


# ---------------------------------------------------------------------------
# The writer verdict reaches the TSV (ISSUE-002, orchestrator follow-up)
# ---------------------------------------------------------------------------

class TestRefusalIsReportedNotSwallowed:
    """The corrected TSV is written a pipeline stage BEFORE any BAM writer runs,
    so a refused surgery would leave the row claiming a rescue the BAM does not
    carry — and a live TSV consumer would draw a junction that is not there."""

    def test_a_drawn_rescue_reports_no_refusal(self):
        read = _make_read(200, [(4, 12), (0, 100)], '+',
                          GENOME_SEQ[88:100] + GENOME_SEQ[200:300])
        corr = _correction(five_prime_intron_clip_pos=200,
                           five_prime_exon_cigar='12M', corrected_3prime=299)
        modified, refusal = apply_5prime_rescue_surgery(read, corr, GENOME)
        assert modified is True and refusal == ''
        assert _n_op_intervals(read) == ((100, 200),)

    def test_a_surgery_that_draws_nothing_reports_a_token(self):
        """Whichever guard stops it, the row must never be left silently
        claiming a rescue."""
        read = _intronic_5prime_read()
        corr = _correction(five_prime_intron_clip_pos=-1,
                           five_prime_exon_cigar='7M')   # wrong query span, no icp
        modified, refusal = apply_5prime_rescue_surgery(read, corr, GENOME)
        assert refusal in (REFUSAL_REROUTE, REFUSAL_EXTEND, REFUSAL_NONCANONICAL)
        assert modified is False

    def test_noncanonical_destination_is_reported(self):
        """A canonical-guard revert must surface as its own token, not as silence."""
        genome = {'chrR': 'A' * 100 + 'TT' + 'C' * 96 + 'TT' + 'G' * 100}
        read = _make_read(200, [(4, 12), (0, 100)], '+', 'A' * 112)
        corr = _correction(five_prime_intron_clip_pos=200,
                           five_prime_exon_cigar='12M', corrected_3prime=299)
        modified, refusal = apply_5prime_rescue_surgery(read, corr, genome)
        assert refusal == REFUSAL_NONCANONICAL
        assert modified is False
        assert _n_op_intervals(read) == ()

    def test_softclip_fallback_is_a_partial_refusal(self):
        """`softclip_intronic_tail_5prime` succeeds by HIDING the intronic bases;
        no junction is drawn, so the junction claim must still be retracted even
        though the read was modified."""
        read = _intronic_5prime_read()
        corr = _correction(five_prime_intron_clip_pos=200,
                           five_prime_exon_cigar='')     # no exon CIGAR -> softclip path
        modified, refusal = apply_5prime_rescue_surgery(read, corr, GENOME)
        assert modified is True
        assert refusal == REFUSAL_SOFTCLIP_ONLY
        assert not [n for n in _n_op_intervals(read) if n == (100, 200)]

    def test_prediction_matches_the_surgery_and_leaves_the_read_alone(self):
        for corr, expect_token in (
            (_correction(five_prime_intron_clip_pos=200,
                         five_prime_exon_cigar='62M'), ''),
            (_correction(five_prime_intron_clip_pos=-1,
                         five_prime_exon_cigar='7M'), True),
        ):
            probe = _intronic_5prime_read()
            before_cigar, before_start = probe.cigartuples, probe.reference_start
            predicted = predict_5prime_rescue_refusal(probe, corr, GENOME)
            # The prediction must not mutate the caller's read.
            assert probe.cigartuples == before_cigar
            assert probe.reference_start == before_start
            # ...and must equal what the surgery actually produces.
            actual_read = _intronic_5prime_read()
            _, actual = apply_5prime_rescue_surgery(actual_read, corr, GENOME)
            assert predicted == actual
            if expect_token is True:
                assert predicted != ''
            else:
                assert predicted == expect_token

    def test_prediction_is_empty_for_an_unrescued_row(self):
        read = _intronic_5prime_read()
        assert predict_5prime_rescue_refusal(
            read, _correction(five_prime_rescued=False), GENOME) == ''

    def test_tsv_header_carries_the_column_last(self):
        from rectify.core.bam.output import (
            CORRECTION_TSV_HEADER, correction_result_to_tsv_row,
        )
        # ISSUE-017 (2026-09-05) appended the two provenance/evidence columns
        # after it, ISSUE-026 invariant D the exon-2 prefix after those; every
        # earlier column keeps its absolute index.
        assert CORRECTION_TSV_HEADER[-4:] == [
            'five_prime_rescue_refused', 'five_prime_landing_annotated',
            'five_prime_novel_evidence', 'five_prime_exon2_prefix']
        row = correction_result_to_tsv_row({
            'read_id': 'r', 'chrom': 'chrR', 'strand': '+',
            'original_3prime': 1, 'corrected_3prime': 1,
            'ambiguity_min': 0, 'ambiguity_max': 0, 'ambiguity_range': 0,
            'correction_applied': [], 'confidence': 'high', 'qc_flags': [],
            'five_prime_rescue_refused': REFUSAL_REROUTE,
        })
        assert len(row) == len(CORRECTION_TSV_HEADER)
        assert row[-4] == REFUSAL_REROUTE
        assert row[-3:] == ['', '', '0']

    def test_all_three_writers_share_one_implementation(self):
        """write_corrected_bam / write_softclipped_bam / write_dual_bam used to
        carry copy-pasted 5' blocks, and the ISSUE-002 fix first landed in only
        one of them. Assert the duplication is gone."""
        import inspect

        from rectify.core.bam import bam_writer as bw

        src = inspect.getsource(bw)
        # The icp gate must be CALLED exactly once, inside the shared helper.
        assert src.count('projected_5prime_rescue_intron_edge(') == 1, (
            'the icp gate should have exactly one call site '
            '(apply_5prime_rescue_surgery); more means a writer grew its own copy')
        assert src.count('apply_5prime_rescue_surgery(read, correction, genome)') == 3, (
            'all three writers must delegate to the shared helper')
        for fn in (bw.write_corrected_bam, bw.write_softclipped_bam, bw.write_dual_bam):
            body = inspect.getsource(fn)
            assert 'extend_read_5prime_for_junction_rescue(' not in body, fn.__name__
