"""Tests for the triage clip legs (landscape §4 step 2 wiring).

Covers: (a) default-off identity — with ``clip_legs_enable=False`` the triage
output is byte-identical and no clip-leg/gate code runs; (b) refusal as a
first-class outcome — poly(A) and repeat-expansion clips are counted and the
read passes through unchanged, never sequence-searched; (c) the terminal-clip
leg produces a resolver proposal that enters the strict hp_ed adjudication;
(d) the 5'-clip leg routes through the shared informativeness gate (assessed
exactly ONCE) into the Cat3 rescue machinery and its proposal is adjudicated.
"""

import random
from pathlib import Path

import pysam
import pytest

import rectify
from rectify.core.consensus.triage import (
    REASON_CLIP_3P,
    REASON_CLIP_5P,
    TriagePolicy,
    classify_read,
    triage_realign_bam,
)
from rectify.core.splice.overhang_informativeness import COUNTERS, reset_counters

DATA = Path(rectify.__file__).parent / 'data'
VALIDATION_BAM = DATA / 'validation' / 'aligners' / 'validation_reads.minimap2.bam'
_SCER = DATA / 'genomes' / 'saccharomyces_cerevisiae'
GENOME_FSA = _SCER / 'S288C_reference_sequence_R64-5-1_20240529.fsa.gz'
ANNOTATION = _SCER / 'saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz'

# ---------------------------------------------------------------------------
# Synthetic genome: one planted + intron [1200, 1500) (GT..AG) on 'chrI'
# (a real chrom name so standardize_chrom_name is a no-op end to end).
# ---------------------------------------------------------------------------

GLEN = 4200
DON, ACC = 1200, 1500          # intron [1200, 1500), length 300


def _make_genome():
    rng = random.Random(20260811)
    seq = [rng.choice('ACGT') for _ in range(GLEN)]
    seq[DON:DON + 2] = ['G', 'T']
    seq[ACC - 2:ACC] = ['A', 'G']
    return ''.join(seq)


GENOME_SEQ = _make_genome()
GENOME = {'chrI': GENOME_SEQ}
ANNOTATED = {('chrI', DON, ACC, '+')}
HEADER = pysam.AlignmentHeader.from_dict(
    {'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'chrI', 'LN': GLEN}]})


@pytest.fixture(autouse=True)
def _fresh_counters():
    reset_counters()
    yield
    reset_counters()


def _read(name, query, cigar, start, reverse=False):
    r = pysam.AlignedSegment(header=HEADER)
    r.query_name = name
    r.query_sequence = query
    r.flag = 16 if reverse else 0
    r.reference_id = 0
    r.reference_start = start
    r.mapping_quality = 60
    r.cigarstring = cigar
    return r


def _write_bam(path, reads):
    with pysam.AlignmentFile(str(path), 'wb', header=HEADER) as out:
        for r in reads:
            out.write(r)
    return str(path)


def _records(path):
    with pysam.AlignmentFile(str(path), 'rb', check_sq=False) as bam:
        return [r.to_string() for r in bam.fetch(until_eof=True)]


def _run(bam_path, out_path, **kw):
    kw.setdefault('genome', GENOME)
    kw.setdefault('annotated_junctions', ANNOTATED)
    kw.setdefault('sort_and_index', False)
    return triage_realign_bam(str(bam_path), str(out_path), **kw)


LEGS_ON = TriagePolicy(clip_legs_enable=True)


# ---------------------------------------------------------------------------
# (a) Default-off identity
# ---------------------------------------------------------------------------

class TestDefaultOffIdentity:
    def test_flag_defaults_false(self):
        assert TriagePolicy().clip_legs_enable is False

    @pytest.mark.skipif(not VALIDATION_BAM.exists(),
                        reason='bundled validation BAM absent')
    def test_bundled_output_identical_and_no_leg_code_runs(self, tmp_path):
        from rectify.utils.genome import load_genome
        from rectify.core.consensus.consensus import load_annotated_junctions
        genome = load_genome(GENOME_FSA)
        annotated = load_annotated_junctions(str(ANNOTATION))

        out_default = tmp_path / 'default.bam'
        rows_d, stats_d = triage_realign_bam(
            str(VALIDATION_BAM), str(out_default), genome=genome,
            annotated_junctions=annotated, sort_and_index=False)

        out_off = tmp_path / 'flag_off.bam'
        rows_o, stats_o = triage_realign_bam(
            str(VALIDATION_BAM), str(out_off), genome=genome,
            annotated_junctions=annotated,
            policy=TriagePolicy(clip_legs_enable=False), sort_and_index=False)

        # Byte-identical records + identical rows, and the informativeness
        # gate was never touched (no clip-leg code ran on the default path).
        assert _records(out_default) == _records(out_off)
        assert rows_d == rows_o
        assert COUNTERS['assessed'] == 0
        for key in ('clip5_leg', 'clip3_leg', 'clip5_refused', 'clip3_refused',
                    'clip5_proposed', 'clip3_proposed',
                    'clip5_accepted', 'clip3_accepted'):
            assert stats_d[key] == 0
            assert stats_o[key] == 0


# ---------------------------------------------------------------------------
# (b) Refusal is a first-class outcome (counted; read unchanged; no search)
# ---------------------------------------------------------------------------

class TestRefusals:
    def test_polya_3p_clip_refused_read_unchanged(self, tmp_path):
        # 60M over exon-1 end + a 40 bp poly(A) 3' clip (> max_clip_3p).
        query = GENOME_SEQ[DON - 60:DON] + 'A' * 40
        r = _read('polyA3', query, '60M40S', DON - 60)
        assert REASON_CLIP_3P in classify_read(r, GENOME).reasons
        bam = _write_bam(tmp_path / 'in.bam', [_read('polyA3', query, '60M40S', DON - 60)])
        out = tmp_path / 'out.bam'
        rows, stats = _run(bam, out, policy=LEGS_ON)
        assert stats['clip3_leg'] == 1
        assert stats['clip3_refused'] == 1
        assert stats['clip3_proposed'] == 0
        assert stats['clip3_accepted'] == 0
        assert COUNTERS['candidates_evaluated'] == 0   # never sequence-searched
        assert _records(out) == _records(bam)          # passthrough unchanged

    def test_repeat_expansion_5p_clip_refused_read_unchanged(self, tmp_path):
        # (AAG)14 = 42 bp 5' clip (> max_clip_5p) — basecaller artifact class.
        query = 'AAG' * 14 + GENOME_SEQ[ACC:ACC + 60]
        r = _read('rep5', query, '42S60M', ACC)
        assert REASON_CLIP_5P in classify_read(r, GENOME).reasons
        bam = _write_bam(tmp_path / 'in.bam', [_read('rep5', query, '42S60M', ACC)])
        out = tmp_path / 'out.bam'
        rows, stats = _run(bam, out, policy=LEGS_ON)
        assert stats['clip5_leg'] == 1
        assert stats['clip5_refused'] == 1
        assert stats['clip5_proposed'] == 0
        assert _records(out) == _records(bam)

    def test_polya_5p_clip_refused_via_assessment(self, tmp_path):
        # Poly(A) passes the repeat-expansion check but is refused by
        # assess_overhang (W_max < 1) — the counter proves the gate ran ONCE.
        query = 'A' * 40 + GENOME_SEQ[ACC:ACC + 60]
        bam = _write_bam(tmp_path / 'in.bam', [_read('polyA5', query, '40S60M', ACC)])
        out = tmp_path / 'out.bam'
        rows, stats = _run(bam, out, policy=LEGS_ON)
        assert stats['clip5_refused'] == 1
        assert COUNTERS['assessed'] == 1
        assert COUNTERS['refused'] == 1
        assert _records(out) == _records(bam)


# ---------------------------------------------------------------------------
# (c) Terminal-clip leg: resolver proposal enters adjudication
# ---------------------------------------------------------------------------

class TestTerminalClipLeg:
    def test_exon2_clip_resolved_and_accepted(self, tmp_path):
        # Aligned over exon-1 end; 40 bp right clip = exon-2 head. The
        # resolver proposes M60 N300 M40; hp_ed drops 40 (S) -> 0: accepted.
        query = GENOME_SEQ[DON - 60:DON] + GENOME_SEQ[ACC:ACC + 40]
        bam = _write_bam(tmp_path / 'in.bam',
                         [_read('clip3', query, '60M40S', DON - 60)])
        out = tmp_path / 'out.bam'
        rows, stats = _run(bam, out, policy=LEGS_ON)
        assert stats['clip3_leg'] == 1
        assert stats['clip3_refused'] == 0
        assert stats['clip3_proposed'] == 1
        assert stats['clip3_accepted'] == 1
        with pysam.AlignmentFile(str(out), 'rb', check_sq=False) as f:
            rec = next(f.fetch(until_eof=True))
        assert rec.cigartuples == [(0, 60), (3, ACC - DON), (0, 40)]
        assert rec.get_tag('XJ').startswith(f'{DON}-{ACC}:')
        row = {r['read_id']: r for r in rows}['clip3']
        assert row['realigned'] and row['accepted']

    def test_garbage_clip_no_proposal_read_unchanged(self, tmp_path):
        # Informative (passes the gate) but matches nothing: no proposal,
        # nothing accepted — the read passes through unchanged.
        rng = random.Random(99)
        garbage = ''.join(rng.choice('ACGT') for _ in range(40))
        query = GENOME_SEQ[DON - 60:DON] + garbage
        bam = _write_bam(tmp_path / 'in.bam',
                         [_read('garbage3', query, '60M40S', DON - 60)])
        out = tmp_path / 'out.bam'
        rows, stats = _run(bam, out, policy=LEGS_ON)
        assert stats['clip3_leg'] == 1
        assert stats['clip3_proposed'] == 0
        assert stats['clip3_accepted'] == 0
        assert _records(out) == _records(bam)


# ---------------------------------------------------------------------------
# (d) 5'-clip leg: gate (exactly once) -> Cat3 rescue -> adjudication
# ---------------------------------------------------------------------------

class TestFivePrimeClipLeg:
    def test_exon1_tail_clip_rescued_and_accepted(self, tmp_path):
        # + strand read starting at intron_end with a 32 bp 5' clip that IS
        # the exon-1 tail. The Cat3 rescue proposes 32M 300N 60M; accepted.
        query = GENOME_SEQ[DON - 32:DON] + GENOME_SEQ[ACC:ACC + 60]
        bam = _write_bam(tmp_path / 'in.bam',
                         [_read('cat3', query, '32S60M', ACC)])
        out = tmp_path / 'out.bam'
        rows, stats = _run(bam, out, policy=LEGS_ON)
        assert stats['clip5_leg'] == 1
        assert stats['clip5_refused'] == 0
        assert stats['clip5_proposed'] == 1
        assert stats['clip5_accepted'] == 1
        # ONE refusal discipline: the overhang was assessed exactly once.
        assert COUNTERS['assessed'] == 1
        with pysam.AlignmentFile(str(out), 'rb', check_sq=False) as f:
            rec = next(f.fetch(until_eof=True))
        assert rec.reference_start == DON - 32
        assert rec.cigartuples == [(0, 32), (3, ACC - DON), (0, 60)]
        row = {r['read_id']: r for r in rows}['cat3']
        assert row['realigned'] and row['accepted']

    def test_minus_strand_exon1_tail_clip_rescued(self, tmp_path):
        # Mirror geometry on a minus-strand read: the 5' clip is the RIGHT
        # clip (transcript orientation) and the transcript's first exon sits
        # at HIGHER coordinates. The read aligns [DON-60, DON) with its 5'
        # end at reference_end (= intron_start); the 32 bp clip is the
        # cross-intron exon head at [ACC, ACC+32). The rescue proposes
        # 60M 300N 32M in place; adjudication accepts.
        query = GENOME_SEQ[DON - 60:DON] + GENOME_SEQ[ACC:ACC + 32]
        bam = _write_bam(tmp_path / 'in.bam',
                         [_read('minus5', query, '60M32S', DON - 60, reverse=True)])
        out = tmp_path / 'out.bam'
        rows, stats = _run(bam, out, policy=LEGS_ON)
        assert stats['clip5_leg'] == 1          # routed (5' = right clip on '-')
        assert stats['clip5_proposed'] == 1
        assert stats['clip5_accepted'] == 1
        assert COUNTERS['assessed'] == 1        # gate ran exactly once
        with pysam.AlignmentFile(str(out), 'rb', check_sq=False) as f:
            rec = next(f.fetch(until_eof=True))
        assert rec.reference_start == DON - 60
        assert rec.cigartuples == [(0, 60), (3, ACC - DON), (0, 32)]
