"""Module 2F acceptance guards: minimum informative clip length (ISSUE-006) and
the writer's canonical-destination guard.

Both landed together because the tester's clip-length table on the Sumner human
RNA004 chr5 slice shows a length floor ALONE is not enough: the <10 nt bins are
>=80% non-canonical and hold 92% of all rescues, but even 15-29 nt rescues are
40-58% non-canonical, so the floor is paired with a destination guard.

Run with:
    pytest tests/test_2f_min_informative_clip.py -v
"""

import os
from pathlib import Path
from typing import List, Tuple

import pysam
import pytest

from rectify.core.bam.bam_writer import (
    _n_op_intervals,
    _revert_selfinflicted_noncanonical_n,
)
from rectify.core.splice.splice_aware_5prime import (
    MAX_RESCUE_JUNCTIONS,
    MAX_SS_SHIFT,
    DEFAULT_JUNCTION_PROXIMITY_BP,
    _clip_is_periodic,
    _clip_search_refused,
    min_informative_clip_bp,
    rescue_3ss_truncation,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_read(chrom: str, start: int, cigar: List[Tuple[int, int]],
               strand: str, seq: str, name: str = 'r') -> pysam.AlignedSegment:
    hdr = pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.6'},
        'SQ': [{'SN': chrom, 'LN': 3_000_000}],
    })
    read = pysam.AlignedSegment(hdr)
    read.query_name = name
    read.reference_name = chrom
    read.reference_start = start
    read.cigartuples = cigar
    read.is_reverse = (strand == '-')
    read.is_unmapped = False
    read.is_secondary = False
    read.is_supplementary = False
    read.mapping_quality = 60
    read.query_sequence = seq
    return read


# A toy plus-strand locus with a real (non-homopolymer) exon-1 tail so the clip
# under test is genuinely informative:
#   exon1  [0, 40)    ends with the 20 nt EXON1_TAIL
#   intron [40, 140)  GT ... AG
#   exon2  [140, 240)
EXON1_TAIL = 'ACGTTGCATGCAGTCCATG'          # 19 nt, aperiodic
_GENOME_SEQ = (
    ('T' * (40 - len(EXON1_TAIL) - 1)) + 'A' + EXON1_TAIL   # pos 0-39
    + 'GT' + 'N' * 96 + 'AG'                                # pos 40-139 (intron)
    + 'C' * 100                                             # pos 140-239 (exon2)
)
GENOME = {'chrT': _GENOME_SEQ}
JUNCTION = {('chrT', 40, 140)}


def _plus_read_with_clip(clip_len: int) -> pysam.AlignedSegment:
    """Read starting at exon2 whose 5' soft clip is the last `clip_len` bases of
    exon 1 — i.e. a perfect (ed_exon == 0) exon-1 match at the annotated donor."""
    clip = _GENOME_SEQ[40 - clip_len:40]
    return _make_read('chrT', 140, [(4, clip_len), (0, 60)], '+',
                      clip + 'C' * 60, name=f'clip{clip_len}')


# ---------------------------------------------------------------------------
# The bound itself
# ---------------------------------------------------------------------------

class TestMinInformativeClipBp:
    def test_lands_at_ten_on_the_shipped_constants(self):
        """W = 25 x 31 x 11 = 8525 windows, alpha = 0.01 -> 19.70 bits -> 10 nt.

        10 nt is also the shortest 5' soft clip among the bundled yeast cat3
        fixtures (10 / 13 / 16 / 22 nt), so the bound and the FN guard agree.
        """
        assert min_informative_clip_bp() == 10

    def test_is_derived_from_the_search_space_not_hardcoded(self):
        """A smaller search space must lower the bar; the constant must track it."""
        import math

        for prox in (0, 5, 10, 50):
            w = MAX_RESCUE_JUNCTIONS * (2 * MAX_SS_SHIFT + 1) * (prox + 1)
            expected = max(1, int(math.ceil(math.log2(w / 0.01) / 2.0)))
            assert min_informative_clip_bp(prox) == expected
        assert (min_informative_clip_bp(0)
                <= min_informative_clip_bp(DEFAULT_JUNCTION_PROXIMITY_BP)
                <= min_informative_clip_bp(1000))


# ---------------------------------------------------------------------------
# ISSUE-006 (a): short clips are not evidence
# ---------------------------------------------------------------------------

class TestShortClipsRefused:
    @pytest.mark.parametrize('clip_len', [1, 2, 3, 5, 9])
    def test_sub_floor_clip_not_rescued(self, clip_len):
        """A 1-9 nt clip reaches ed_exon == 0 trivially; it must not rescue.

        Regression for the 42% fire rate: 17 of 26 rescued rows on the Sumner
        panel had a <=3 nt clip.
        """
        read = _plus_read_with_clip(clip_len)
        r = rescue_3ss_truncation(read, GENOME, JUNCTION, strand='+')
        assert r['rescue_type'] != 'softclip'
        assert not (r['rescued'] and r['rescue_type'] == 'softclip')

    @pytest.mark.parametrize('clip_len', [10, 12, 15, 19])
    def test_at_or_above_floor_still_rescues(self, clip_len):
        """The floor must not disable the rescue it is guarding — an informative
        clip at or above it still reaches exon 1 (FN side of the guard)."""
        read = _plus_read_with_clip(clip_len)
        r = rescue_3ss_truncation(read, GENOME, JUNCTION, strand='+')
        assert r['rescued'] is True, f'{clip_len} nt informative clip was refused'
        assert r['rescue_type'] == 'softclip'
        assert r['rescued_junction'] == ('chrT', 40, 140)
        assert r['five_prime_corrected'] == 39   # intron_start - 1

    def test_boundary_is_exactly_the_derived_value(self):
        floor = min_informative_clip_bp()
        below = rescue_3ss_truncation(_plus_read_with_clip(floor - 1),
                                      GENOME, JUNCTION, strand='+')
        at = rescue_3ss_truncation(_plus_read_with_clip(floor),
                                   GENOME, JUNCTION, strand='+')
        assert below['rescue_type'] != 'softclip'
        assert at['rescue_type'] == 'softclip' and at['rescued'] is True

    def test_structural_paths_stay_live_when_the_search_is_refused(self):
        """The floor refuses the SEQUENCE search, not the read.

        A minus-strand read whose 5' edge sits inside the intron must still reach
        the Case-4 intronic snap even though its clip is below the floor —
        clearing `rescue_seq` instead of flagging it silently suppressed 11 of the
        Sumner panel's intronic-snap rescues.
        """
        # 5' edge (reference_start, plus strand) 10 bases inside the intron.
        read = _make_read('chrT', 50, [(4, 2), (0, 90)], '+',
                          'CA' + _GENOME_SEQ[50:140], name='snap')
        r = rescue_3ss_truncation(read, GENOME, JUNCTION, strand='+')
        assert r['rescue_type'] != 'softclip'
        # The structural (non-sequence) evidence path is still reachable.
        assert r['rescue_type'] in ('intronic_snap', 'proximity', 'none')


# ---------------------------------------------------------------------------
# ISSUE-006 (b): the is_repeat_expansion refusal, re-keyed on informativeness
# ---------------------------------------------------------------------------

class TestClipComplexityCriterion:
    @pytest.mark.parametrize('seq', ['AAG' * 5, 'AAG' * 13, 'CTT' * 30, 'AAG' * 130])
    def test_repeat_expansions_refused_at_every_length(self, seq):
        """The old rule needed >= 30 bp; a 15 nt (AAG)5 slipped through."""
        assert _clip_search_refused(seq, '+') is True
        assert _clip_is_periodic(seq, '+') is True

    @pytest.mark.parametrize('seq', ['A' * 16, 'AG' * 8, 'A' * 40])
    def test_homopolymer_and_dinucleotide_refused(self, seq):
        assert _clip_search_refused(seq, '+') is True

    def test_informative_long_clip_is_assessed_not_refused(self):
        """The inversion ISSUE-006 names: a genuinely informative 33-388 nt clip
        was refused by a LENGTH threshold. 40 nt of aperiodic sequence must pass."""
        seq = (EXON1_TAIL + 'GGATCCTAGCTAGGCATTAGCACC')[:40]
        assert len(seq) == 40
        assert _clip_search_refused(seq, '+') is False
        assert _clip_is_periodic(seq, '+') is False

    def test_yeast_cat3_clips_are_not_refused(self):
        """The four bundled cat3 fixtures' actual 5' soft clips (per-aligner BAMs).
        A change that refuses any of these fails the FN guard."""
        for clip in ('GAGGAAAAAT',                 # cat3_plus_1,  10 nt
                     'TCATATGTAGACA',              # cat3_minus_2, 13 nt
                     'TTTTTCTTTGCTTAAA',           # cat3_minus_1, 16 nt
                     'GGCTGACAAGTCATCATTGAAG'):    # cat3_plus_2,  22 nt
            assert len(clip) >= min_informative_clip_bp()
            assert _clip_is_periodic(clip, '+') is False, clip

    def test_period_criterion_is_silent_below_its_own_resolution(self):
        """min_self_match_period needs a 12 bp overlap, so it cannot fire on a
        clip the length floor has already refused — no double judgement."""
        assert _clip_is_periodic('AGAGAGAGAG', '+') is False   # 10 nt, period 2
        assert _clip_is_periodic('AGAGAGAGAGAGAG', '+') is True  # 14 nt, period 2


# ---------------------------------------------------------------------------
# agent2 Fix 3 — the writer's canonical-destination guard
# ---------------------------------------------------------------------------

class TestWriterCanonicalGuard:
    # exon1 [0,20) | intron [20,120) | exon2 [120,220)
    GUARD_GENOME = {'chrG': 'A' * 20 + 'GT' + 'C' * 96 + 'AG' + 'T' * 100}

    def _read(self, cigar, start=0):
        return _make_read('chrG', start, cigar, '+', 'A' * 200, name='g')

    def test_reverts_a_writer_invented_noncanonical_intron(self):
        pre = self._read([(0, 20), (0, 100)])
        pre_nops = _n_op_intervals(pre)
        assert pre_nops == ()
        # "Surgery": the writer draws an N-op at 21-121 (T..CA -> non-canonical).
        post = self._read([(0, 20), (3, 100), (0, 80)], start=1)
        reverted = _revert_selfinflicted_noncanonical_n(
            post, self.GUARD_GENOME, pre_nops, pre.cigartuples, pre.reference_start)
        assert reverted is True
        assert post.cigartuples == pre.cigartuples
        assert post.reference_start == pre.reference_start

    def test_keeps_a_writer_invented_canonical_intron(self):
        pre = self._read([(0, 20), (0, 100)])
        pre_nops = _n_op_intervals(pre)
        post = self._read([(0, 20), (3, 100), (0, 80)], start=0)   # N = 20-120, GT..AG
        reverted = _revert_selfinflicted_noncanonical_n(
            post, self.GUARD_GENOME, pre_nops, pre.cigartuples, pre.reference_start)
        assert reverted is False
        assert post.cigartuples == [(0, 20), (3, 100), (0, 80)]

    def test_never_judges_a_preexisting_aligner_n_op(self):
        """Not a motif filter on aligner evidence: an N-op present BEFORE the
        surgery is untouched however non-canonical its motif."""
        pre = self._read([(0, 20), (3, 100), (0, 80)], start=1)    # non-canonical
        pre_nops = _n_op_intervals(pre)
        post = self._read([(0, 20), (3, 100), (0, 80)], start=1)
        reverted = _revert_selfinflicted_noncanonical_n(
            post, self.GUARD_GENOME, pre_nops, pre.cigartuples, pre.reference_start)
        assert reverted is False

    def test_atac_is_canonical(self):
        genome = {'chrG': 'A' * 20 + 'AT' + 'C' * 96 + 'AC' + 'T' * 100}
        pre = self._read([(0, 20), (0, 100)])
        post = self._read([(0, 20), (3, 100), (0, 80)], start=0)
        assert _revert_selfinflicted_noncanonical_n(
            post, genome, _n_op_intervals(pre), pre.cigartuples,
            pre.reference_start) is False

    def test_no_genome_is_a_no_op(self):
        pre = self._read([(0, 20), (0, 100)])
        post = self._read([(0, 20), (3, 100), (0, 80)], start=1)
        assert _revert_selfinflicted_noncanonical_n(
            post, None, _n_op_intervals(pre), pre.cigartuples,
            pre.reference_start) is False


# ---------------------------------------------------------------------------
# Human RNA004 panel fixture (git-ignored inputs -> skipif)
# ---------------------------------------------------------------------------

_PANEL_DIR = Path(os.environ.get(
    'RECTIFY_SUMNER_PANEL_DIR',
    '/Users/kevinroy/work/rectify/dev/sumner_misplaced_panel_20260904'))
_PANEL_BAM = _PANEL_DIR / '175_panel_orig.bam'
_PANEL_FA = _PANEL_DIR / 'ref' / 'chr5.fa'
_PANEL_GTF = _PANEL_DIR / 'ref' / 'gencode.v48.basic.chr5.gtf'
_PANEL_TSV = _PANEL_DIR / '175_panel.tsv'
_HOLDOUT_BAM = _PANEL_DIR / 'holdout' / 'chr5_holdout1k.bam'

_HAVE_PANEL = all(p.exists() for p in (_PANEL_BAM, _PANEL_FA, _PANEL_GTF, _PANEL_TSV))
_HAVE_HOLDOUT = _HOLDOUT_BAM.exists() and _PANEL_FA.exists() and _PANEL_GTF.exists()

_CANONICAL_PAIRS = {('GT', 'AG'), ('CT', 'AC'), ('GC', 'AG'),
                    ('CT', 'GC'), ('AT', 'AC'), ('GT', 'AT')}


def _replay_2f(bam_path, want_categories=False):
    """Run Module 2F + the bam_writer 5' surgery over *bam_path* and return one
    row per primary read. The 3' fields are chosen so the 3' clip is a no-op, so
    every difference is attributable to the modules under test."""
    import copy
    import csv

    import rectify.core.splice.splice_aware_5prime as s5
    from rectify.core.bam.bam_writer import apply_corrected_edits_to_read
    from rectify.core.consensus.consensus import load_annotated_junctions

    # These inputs were produced with the tester's RECTIFY_CHROM_VERBATIM patch;
    # the real chrom fix is ISSUE-001 and belongs to another change.
    saved = s5.standardize_chrom_name
    s5.standardize_chrom_name = lambda c: c
    try:
        fa = pysam.FastaFile(str(_PANEL_FA))
        chrom = fa.references[0]
        genome = {chrom: fa.fetch(chrom)}
        ann_raw = load_annotated_junctions(str(_PANEL_GTF))
        cand = {(chrom,) + tuple(j[1:]) for j in ann_raw}
        ann = {(j[1], j[2]) for j in ann_raw}
        cats = {}
        if want_categories:
            cats = {r['read_id']: r['category'] for r in
                    csv.DictReader(open(_PANEL_TSV), delimiter='\t')}
        out = []
        for read in pysam.AlignmentFile(str(bam_path)):
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
                a5 = (read.reference_end - 1) if strand == '-' else read.reference_start
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
                'category': cats.get(read.query_name, ''),
                # The reanchor pre-pass can collapse a mangled 5' edge into a
                # LONGER leading soft clip; that clip is what the rescue
                # actually searches, so it is the length the floor judges.
                'clip': max(clip, int(res.get('reanchor_clip_len', 0) or 0)),
                'raw_clip': clip,
                'rescued': bool(res['rescued']),
                'rescue_type': res['rescue_type'],
                'before': before,
                'after': _n_op_intervals(edited),
            })
        return {'rows': out, 'seq': genome[chrom], 'ann': ann}
    finally:
        s5.standardize_chrom_name = saved


@pytest.mark.skipif(not _HAVE_PANEL,
                    reason='Sumner human panel inputs are git-ignored; set '
                           'RECTIFY_SUMNER_PANEL_DIR to run')
class TestSumnerHumanPanel:
    """FP/FN guard on the 100-read human RNA004 panel.

    Baseline recorded 2026-09-05 on this branch's parent (4eefd1f): 12
    non-canonical N-ops created by the 5' rescue, 38/100 reads rescued, 0 of the
    5 C1 controls touched. Post-fix: 0 non-canonical created.
    """

    @pytest.fixture(scope='class')
    def replay(self):
        return _replay_2f(_PANEL_BAM, want_categories=True)

    def test_no_noncanonical_intron_is_created(self, replay):
        seq = replay['seq']
        offenders = []
        for row in replay['rows']:
            for start, end in row['after']:
                if (start, end) in row['before']:
                    continue
                pair = (seq[start:start + 2].upper(), seq[end - 2:end].upper())
                if pair not in _CANONICAL_PAIRS:
                    offenders.append((row['read_id'][:12], row['category'],
                                      start, end, pair))
        assert not offenders, (
            'Module 2F created non-canonical introns (baseline before ISSUE-006 '
            f'+ the writer guard was 12):\n  ' +
            '\n  '.join(map(str, offenders)))

    def test_c1_controls_stay_untouched(self, replay):
        c1 = [r for r in replay['rows'] if r['category'] == 'C1']
        assert len(c1) == 5
        for row in c1:
            assert row['rescued'] is False, row['read_id']
            assert row['after'] == row['before'], row['read_id']

    def test_short_clip_rows_no_longer_rescue_on_the_softclip_branch(self, replay):
        floor = min_informative_clip_bp()
        bad = [r['read_id'][:12] for r in replay['rows']
               if r['rescue_type'] == 'softclip' and r['clip'] < floor
               and r['rescued']]
        assert not bad, f'sub-floor softclip rescues survived: {bad}'

    def test_the_panel_still_rescues_something(self, replay):
        """The panel is ADVERSARIAL — its rows were selected because RECTIFY
        changed them, and 92% of its rescues are 1-3 nt clips by construction —
        so it cannot show that the module is still useful. That is the hold-out's
        job (TestChr5Holdout below). Here we only assert the module did not go
        completely silent."""
        assert any(r['rescued'] for r in replay['rows'])


@pytest.mark.skipif(not _HAVE_HOLDOUT,
                    reason='chr5 hold-out is git-ignored; set '
                           'RECTIFY_SUMNER_PANEL_DIR to run')
class TestChr5Holdout:
    """1,000 UNSELECTED primary chr5 reads (panel reads excluded) — the set that
    shows whether the guards generalize rather than fit the adversarial panel.

    Baseline recorded 2026-09-05 on 4eefd1f: 228/1000 rescued, 222 introns
    created of which 158 non-canonical and only 25 annotated (11%), 3 pre-existing
    N-ops lost; scorer FP_total 186 / TP_total 25.
    """

    @pytest.fixture(scope='class')
    def replay(self):
        return _replay_2f(_HOLDOUT_BAM)

    def _created(self, replay):
        out = []
        for row in replay['rows']:
            for n in row['after']:
                if n not in row['before']:
                    out.append(n)
        return out

    def test_creates_no_noncanonical_intron(self, replay):
        """FP guard. Baseline: 158 of 222 created introns were non-canonical."""
        seq = replay['seq']
        bad = [n for n in self._created(replay)
               if (seq[n[0]:n[0] + 2].upper(),
                   seq[n[1] - 2:n[1]].upper()) not in _CANONICAL_PAIRS]
        assert not bad, f'{len(bad)} non-canonical introns created: {bad[:5]}'

    def test_destroys_no_preexisting_junction(self, replay):
        """FP guard. Three mechanisms deleted aligner-called junctions: Case 4
        snapping a read across its OWN junctions, the writer's reroute trimming
        the intronic tail, and the reroute's N-merge fusing a partially trimmed N
        (ISSUE-007). Baseline lost 3 N-ops; after the informative-clip floor alone
        it was 18. Now 0 — the read's junction count never drops."""
        lost = []
        for row in replay['rows']:
            after = set(row['after'])
            for n in row['before']:
                if n not in after and any(m[0] == n[0] or m[1] == n[1] for m in after):
                    continue      # a terminal N extended onto the rescued acceptor
                if n not in after:
                    lost.append((row['read_id'][:12], n))
        assert not lost, f'{len(lost)} pre-existing N-ops destroyed: {lost[:6]}'

    def test_no_read_loses_junctions(self, replay):
        """The ISSUE-007 invariant stated on read counts rather than coordinates:
        a 5' rescue may add an intron, never subtract one."""
        shrank = [(r['read_id'][:12], len(r['before']), len(r['after']))
                  for r in replay['rows'] if len(r['after']) < len(r['before'])]
        assert not shrank, shrank

    def test_rescue_still_fires_and_lands_on_annotation(self, replay):
        """FN guard, and the one the panel cannot give. Baseline: 25 of 222
        created introns were annotated (11%)."""
        created = self._created(replay)
        rescued = [r for r in replay['rows'] if r['rescued']]
        # ISSUE-028 invariant E (2026-09-06) re-pinned this guard: at 4993253 the
        # hold-out rescued 47/1000 and created 27 introns (27 annotated); with the
        # evidence floor (identity >= 0.8, >= 18 bits) it rescues 29 and creates
        # 24 (24 annotated, 0 non-canonical) — three annotated creations whose
        # placed block carried < 18 bits are refused, the rest of the drop is
        # rescues that created no intron. Measured, not chosen; Kevin's call.
        assert len(rescued) >= 25, (
            f'only {len(rescued)}/1000 hold-out reads rescued — the guards have '
            'disabled the module rather than aiming it')
        assert created, 'no introns created at all'
        annotated = [n for n in created if n in replay['ann']]
        assert len(annotated) >= 22, (
            f'{len(annotated)} annotated introns created, baseline was 25 (24 under invariant E)')
        assert len(annotated) / len(created) >= 0.70, (
            f'only {len(annotated)}/{len(created)} created introns are annotated '
            '(baseline 25/222 = 11%)')


# ---------------------------------------------------------------------------
# A rescue may not displace a canonical junction the aligner already called
# ---------------------------------------------------------------------------

class TestNoDisplacementOfAlignerCanonicalJunction:
    """Shape of Sumner hold-out read 34625d8e.

    Its aligner CIGAR carried seven annotated CT-AC junctions; the correction
    replaced the 5'-most, 91378775-91382812, with 91378775-91380924 — one edge
    from that junction and the other from a DIFFERENT annotated junction
    (91380924-91418708), 1,888 nt inside the original intron. Not annotated, but
    CT-AC, so a motif check on the written coordinates cannot see that a real
    junction was destroyed to build it. Only the read's own CIGAR can.
    """

    # exon [0,1000) | intron A [1000,3000) CT..AC | exon [3000,3400)
    # candidate C [2500,3200) also CT..AC — overlaps A, shares neither edge.
    @staticmethod
    def _genome():
        seq = list('ACGTTGCAT' * 400)          # 3600 nt, aperiodic
        def put(at, s):
            seq[at:at + len(s)] = list(s)
        put(1000, 'CT'); put(2998, 'AC')       # intron A, canonical on '-'
        put(2500, 'CT'); put(3198, 'AC')       # candidate C, canonical on '-'
        return {'chrD': ''.join(seq)}

    @staticmethod
    def _read(with_n_op):
        """Minus-strand read whose 5' end (reference_end - 1) sits at ~3004."""
        hdr = pysam.AlignmentHeader.from_dict(
            {'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'chrD', 'LN': 3600}]})
        r = pysam.AlignedSegment(hdr)
        r.query_name = 'displace'
        r.reference_name = 'chrD'
        r.is_reverse = True
        r.is_unmapped = r.is_secondary = r.is_supplementary = False
        r.mapping_quality = 60
        if with_n_op:
            r.reference_start = 900
            r.cigartuples = [(0, 100), (3, 2000), (0, 5), (4, 12)]
            r.query_sequence = 'A' * 117
        else:
            r.reference_start = 2900
            r.cigartuples = [(0, 105), (4, 12)]
            r.query_sequence = 'A' * 117
        return r

    def test_candidate_overlapping_an_aligner_canonical_junction_is_refused(self):
        from rectify.core.splice.overhang_informativeness import COUNTERS
        import rectify.core.splice.splice_aware_5prime as s5

        genome = self._genome()
        cands = {('chrD', 1000, 3000), ('chrD', 2500, 3200)}
        COUNTERS['candidate_would_displace_canonical'] = 0
        res = s5.rescue_3ss_truncation(self._read(True), genome, cands, strand='-')
        assert COUNTERS.get('candidate_would_displace_canonical', 0) >= 1, (
            'the overlapping candidate was not filtered')
        assert res.get('displaced_canonical_refused') is not False or not res['rescued']
        # Whatever else happens, the read must not be rescued ONTO the displacing
        # junction.
        assert res.get('rescued_junction') != ('chrD', 2500, 3200)

    def test_the_same_junction_within_proximity_is_not_a_displacement(self):
        """A boundary nudge of the read's OWN junction is an adjustment, not a
        displacement — the `already_has_n` equivalence."""
        import rectify.core.splice.splice_aware_5prime as s5
        from rectify.core.splice.overhang_informativeness import COUNTERS

        genome = self._genome()
        COUNTERS['candidate_would_displace_canonical'] = 0
        s5.rescue_3ss_truncation(
            self._read(True), genome, {('chrD', 1002, 2998)}, strand='-')
        assert COUNTERS.get('candidate_would_displace_canonical', 0) == 0

    def test_a_read_with_no_junction_is_unaffected(self):
        """The legitimate case the rule must not touch: the aligner called no
        junction, so the rescue is free to create the annotated one."""
        import rectify.core.splice.splice_aware_5prime as s5
        from rectify.core.splice.overhang_informativeness import COUNTERS

        genome = self._genome()
        COUNTERS['candidate_would_displace_canonical'] = 0
        s5.rescue_3ss_truncation(
            self._read(False), genome,
            {('chrD', 1000, 3000), ('chrD', 2500, 3200)}, strand='-')
        assert COUNTERS.get('candidate_would_displace_canonical', 0) == 0

    def test_plus_strand_mirror(self):
        """Same displacement, opposite strand. 34625d8e is minus-strand and the
        first reproduction attempt reasoned plus-strand-first, which put the
        moved boundary ~1,900 nt from the 5' edge instead of 475 — so both
        orientations are pinned."""
        import rectify.core.splice.splice_aware_5prime as s5
        from rectify.core.splice.overhang_informativeness import COUNTERS

        genome = self._genome()
        hdr = pysam.AlignmentHeader.from_dict(
            {'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'chrD', 'LN': 3600}]})
        r = pysam.AlignedSegment(hdr)
        r.query_name = 'displace_plus'
        r.reference_name = 'chrD'
        r.is_reverse = False
        r.is_unmapped = r.is_secondary = r.is_supplementary = False
        r.mapping_quality = 60
        # 5' end (reference_start) just below intron A, which the read spans.
        r.reference_start = 995
        r.cigartuples = [(4, 12), (0, 5), (3, 2000), (0, 100)]
        r.query_sequence = 'A' * 117
        COUNTERS['candidate_would_displace_canonical'] = 0
        res = s5.rescue_3ss_truncation(
            r, genome, {('chrD', 1000, 3000), ('chrD', 900, 2500)}, strand='+')
        assert COUNTERS.get('candidate_would_displace_canonical', 0) >= 1
        assert res.get('rescued_junction') != ('chrD', 900, 2500)

    def test_the_pool_arm_supplies_strandless_candidates(self):
        """Why the first hand-back was wrong. `load_annotated_junctions` yields
        4-tuples carrying the junction's strand, and the per-candidate filter
        drops a '+'-annotated junction for a minus-strand read — which is why
        34625d8e looked unreachable. `bam_processor` adds POOL junctions as
        3-tuples with no strand, so the filter cannot exclude them and the
        junction IS a candidate on a real run. The rule must therefore hold for
        a strand-less candidate."""
        import rectify.core.splice.splice_aware_5prime as s5
        from rectify.core.splice.overhang_informativeness import COUNTERS

        genome = self._genome()
        COUNTERS['candidate_would_displace_canonical'] = 0
        s5.rescue_3ss_truncation(
            self._read(True), genome,
            {('chrD', 1000, 3000), ('chrD', 2500, 3200)},   # 3-tuples: no strand
            strand='-')
        assert COUNTERS.get('candidate_would_displace_canonical', 0) >= 1


def _pseudo_random_seq(n, seed=1729):
    """Deterministic aperiodic ACGT string (no numpy, no test-order coupling)."""
    out = []
    x = seed
    for _ in range(n):
        x = (1103515245 * x + 12345) & 0x7FFFFFFF
        out.append('ACGT'[(x >> 16) & 3])
    return ''.join(out)


class TestOffContigRescueRefused:
    """ISSUE-015: 2F rescued an 8,134-nt 5' soft clip on the 16,569-nt chrM into
    an 8,835 M "exon" at reference_start -5843. `pysam.index` then refused the
    BAM ("pos=-5843 cannot be indexed") and a 1.5 h `correct` exited 1 with the
    corrected BAM unindexed. One record in 145,000 killed the run.
    """

    CHROM_LEN = 16_569

    @staticmethod
    def _genome(n):
        # Keyed by the STANDARDIZED name: rescue_3ss_truncation looks the contig
        # up as standardize_chrom_name(read.reference_name), and 'chrM' maps to
        # yeast's 'chrMito' (ISSUE-001's collision, agent S1). Without this the
        # rescue returns at the "no genome sequence" guard and the contig-bounds
        # rule under test is never reached.
        from rectify.utils.genome import standardize_chrom_name
        return {standardize_chrom_name('chrM'): _pseudo_random_seq(n, seed=7)}

    def _read(self, clip, start, strand='+'):
        hdr = pysam.AlignmentHeader.from_dict(
            {'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'chrM', 'LN': self.CHROM_LEN}]})
        r = pysam.AlignedSegment(hdr)
        r.query_name = 'b06afeb3'
        r.reference_name = 'chrM'
        r.reference_start = start
        r.is_reverse = (strand == '-')
        r.is_unmapped = r.is_secondary = r.is_supplementary = False
        r.mapping_quality = 60
        r.cigartuples = ([(4, clip), (0, 296)] if strand == '+'
                         else [(0, 296), (4, clip)])
        # Genuinely aperiodic: a homopolymer clip is refused earlier by the
        # complexity short-circuit, and so is a tandem repeat of any period
        # <= 32 -- neither would exercise the contig-bounds rule at all.
        r.query_sequence = _pseudo_random_seq(clip + 296)
        return r

    def test_the_issue_015_record_shape_is_refused(self):
        """Stock `chrM:4561 8134S296M…`: the exon would start at 4561 - 8134 =
        -3573, before the contig. Refused at ADMISSION, with a token."""
        import rectify.core.splice.splice_aware_5prime as s5
        from rectify.core.splice.overhang_informativeness import COUNTERS

        genome = self._genome(self.CHROM_LEN)
        read = self._read(8134, 4561)
        before_cigar, before_start = read.cigartuples, read.reference_start
        COUNTERS['clip_exceeds_contig'] = 0
        res = s5.rescue_3ss_truncation(
            read, genome, {('chrM', 2990, 5270)}, strand='+')
        assert res['rescued'] is False
        assert res.get('clip_refused') == 'clip_exceeds_contig'
        assert COUNTERS.get('clip_exceeds_contig', 0) == 1
        # the BAM record itself is untouched
        assert read.cigartuples == before_cigar
        assert read.reference_start == before_start

    def test_minus_strand_mirror_off_the_high_edge(self):
        import rectify.core.splice.splice_aware_5prime as s5

        genome = self._genome(self.CHROM_LEN)
        read = self._read(8134, self.CHROM_LEN - 400, strand='-')
        res = s5.rescue_3ss_truncation(
            read, genome, {('chrM', 2990, 5270)}, strand='-')
        assert res.get('clip_refused') == 'clip_exceeds_contig'

    def test_a_long_but_in_bounds_clip_is_still_assessed(self):
        """The rule is geometric, not a length ban: GENCODE chr5 first exons run
        to 12,358 nt, so a multi-kb clip that FITS must still reach the rescue."""
        import rectify.core.splice.splice_aware_5prime as s5

        genome = self._genome(self.CHROM_LEN)
        read = self._read(8134, 12_000)          # 12000 - 8134 = 3866 >= 0
        res = s5.rescue_3ss_truncation(
            read, genome, {('chrM', 2990, 5270)}, strand='+')
        assert res.get('clip_refused') is None

    def test_absurd_clip_length_is_capped(self):
        import rectify.core.splice.splice_aware_5prime as s5

        genome = self._genome(200_000)
        read = self._read(s5.MAX_RESCUABLE_CLIP_BP + 1, 150_000)
        res = s5.rescue_3ss_truncation(
            read, genome, {('chrM', 100_000, 120_000)}, strand='+')
        assert res.get('clip_refused') == 'clip_too_long'

    def test_writer_never_emits_a_negative_position(self):
        """The invariant behind the admission rule: whichever helper ran, a
        surgery that placed the read off the contig is reverted. The ISSUE-015
        read reached `reroute`, not `extend`, which is why 17d1c35's off-edge
        refusal (inside extend) did not cover it."""
        from rectify.core.bam.bam_writer import (
            REFUSAL_OFF_CONTIG, apply_5prime_rescue_surgery,
        )

        genome = self._genome(self.CHROM_LEN)
        read = self._read(8134, 4561)
        before_cigar, before_start = read.cigartuples, read.reference_start
        corr = {
            'five_prime_rescued': True,
            'five_prime_position': 2990,
            'five_prime_soft_clip': 8134,
            'five_prime_exon_cigar': '8835M',
            'five_prime_upstream_trim': 0,
            'five_prime_intron_clip_pos': 5270,
            'reanchor_clip_len': 0,
            'strand': '+',
        }
        modified, refusal = apply_5prime_rescue_surgery(read, corr, genome)
        assert read.reference_start >= 0
        if refusal == REFUSAL_OFF_CONTIG:
            assert modified is False
            assert read.cigartuples == before_cigar
            assert read.reference_start == before_start
