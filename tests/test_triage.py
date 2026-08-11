"""Tests for the consensus triage layer (rectify.core.consensus.triage).

Covers the read-evidence classifier (high-confidence bypass vs triage), the
strict-improvement hp_ed re-entry gate, and an end-to-end smoke on the bundled
validation reads (classify + motif-blind junction re-align).
"""

import itertools
from pathlib import Path

import pysam
import pytest

import rectify
from rectify.core.consensus.triage import (
    LABEL_HIGH_CONFIDENCE,
    LABEL_TRIAGE,
    REASON_CLIP_5P,
    REASON_JUNCTION_ERRORS,
    REASON_UNANNOTATED_JUNCTION,
    TriagePolicy,
    classify_bam,
    classify_read,
    reentry_accept,
    triage_realign_bam,
)

DATA = Path(rectify.__file__).parent / 'data'
VALIDATION_BAM = DATA / 'validation' / 'aligners' / 'validation_reads.minimap2.bam'
_SCER = DATA / 'genomes' / 'saccharomyces_cerevisiae'
GENOME_FSA = _SCER / 'S288C_reference_sequence_R64-5-1_20240529.fsa.gz'
ANNOTATION = _SCER / 'saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz'


# ---------------------------------------------------------------------------
# Constructed-read fixtures (no HP runs: cycling ACGT — every error weight 1.0)
# ---------------------------------------------------------------------------

CHROM = 'chrT'
GENOME_SEQ = ''.join(itertools.islice(itertools.cycle('ACGT'), 400))
GENOME = {CHROM: GENOME_SEQ}
HEADER = pysam.AlignmentHeader.from_dict(
    {'HD': {'VN': '1.6'}, 'SQ': [{'SN': CHROM, 'LN': len(GENOME_SEQ)}]})

# Layout: exon1 = [50, 70), intron = [70, 170), exon2 = [170, 190).
JUNC = (CHROM, 70, 170, '+')
ANNOTATED = {JUNC}


def _spliced_read(query=None, cigar='20M100N20M', start=50, name='r1'):
    read = pysam.AlignedSegment(header=HEADER)
    read.query_name = name
    read.reference_id = 0
    read.reference_start = start
    read.mapping_quality = 60
    read.cigarstring = cigar
    if query is None:
        query = GENOME_SEQ[50:70] + GENOME_SEQ[170:190]
    read.query_sequence = query
    read.query_qualities = pysam.qualitystring_to_array('I' * len(query))
    read.flag = 0
    return read


class TestClassifier:
    def test_clean_read_is_high_confidence(self):
        res = classify_read(_spliced_read(), GENOME, {JUNC[:3]}, TriagePolicy())
        assert res.label == LABEL_HIGH_CONFIDENCE
        assert res.reasons == []
        assert res.n_junctions == 1
        assert res.n_unannotated == 0
        assert res.junction_proximal_errors == 0.0

    def test_junction_proximal_mismatches_triage(self):
        # Mutate two bases inside the 5bp exon-side windows: query idx 18
        # (ref 68, donor window) and query idx 22 (ref 172, acceptor window).
        q = list(GENOME_SEQ[50:70] + GENOME_SEQ[170:190])
        for i in (18, 22):
            q[i] = 'A' if q[i] != 'A' else 'C'
        res = classify_read(_spliced_read(query=''.join(q)), GENOME,
                            {JUNC[:3]}, TriagePolicy())
        assert res.label == LABEL_TRIAGE
        assert REASON_JUNCTION_ERRORS in res.reasons
        assert res.junction_proximal_errors == pytest.approx(2.0)
        assert res.junction_leg

    def test_distal_mismatch_does_not_triage(self):
        # A mismatch far from the junction (query idx 5 -> ref 55) is not a
        # junction-proximal signal.
        q = list(GENOME_SEQ[50:70] + GENOME_SEQ[170:190])
        q[5] = 'A' if q[5] != 'A' else 'C'
        res = classify_read(_spliced_read(query=''.join(q)), GENOME,
                            {JUNC[:3]}, TriagePolicy())
        assert res.junction_proximal_errors == 0.0
        assert res.label == LABEL_HIGH_CONFIDENCE

    def test_softclip_5p_triage(self):
        clip = 'A' * 40
        q = clip + GENOME_SEQ[50:70] + GENOME_SEQ[170:190]
        res = classify_read(
            _spliced_read(query=q, cigar='40S20M100N20M'), GENOME,
            {JUNC[:3]}, TriagePolicy())
        assert res.label == LABEL_TRIAGE
        assert REASON_CLIP_5P in res.reasons
        assert res.clip_5p == 40
        assert not res.junction_leg  # clip-only triage waits for the clip legs

    def test_unannotated_junction_triage_and_policy_off(self):
        res = classify_read(_spliced_read(), GENOME, set(), TriagePolicy())
        assert res.label == LABEL_TRIAGE
        assert res.reasons == [REASON_UNANNOTATED_JUNCTION]
        assert res.n_unannotated == 1

        off = TriagePolicy(triage_unannotated_junctions=False)
        res2 = classify_read(_spliced_read(), GENOME, set(), off)
        assert res2.label == LABEL_HIGH_CONFIDENCE

    def test_no_annotation_set_never_flags_unannotated(self):
        res = classify_read(_spliced_read(), GENOME, None, TriagePolicy())
        assert res.label == LABEL_HIGH_CONFIDENCE
        assert res.n_unannotated == 0


class TestReentryGate:
    def _plain(self, start, n_mut, name):
        q = list(GENOME_SEQ[start:start + 40])
        for i in range(n_mut):
            q[i] = 'A' if q[i] != 'A' else 'C'
        read = pysam.AlignedSegment(header=HEADER)
        read.query_name = name
        read.reference_id = 0
        read.reference_start = start
        read.mapping_quality = 60
        read.cigarstring = '40M'
        read.query_sequence = ''.join(q)
        read.flag = 0
        return read

    def test_identical_signature_rejected(self):
        a = self._plain(50, 2, 'x')
        b = self._plain(50, 0, 'x')  # same start + cigar => same signature
        assert reentry_accept(a, b, GENOME) is False

    def test_strictly_better_accepted(self):
        orig = self._plain(50, 2, 'x')       # hp_ed 2
        cand = self._plain(100, 0, 'x')      # hp_ed 0, different placement
        assert reentry_accept(orig, cand, GENOME) is True

    def test_worse_and_tied_rejected(self):
        orig = self._plain(50, 2, 'x')
        worse = self._plain(100, 5, 'x')
        tied = self._plain(100, 2, 'x')
        assert reentry_accept(orig, worse, GENOME) is False
        assert reentry_accept(orig, tied, GENOME) is False


@pytest.mark.skipif(not VALIDATION_BAM.exists(), reason='bundled validation BAM absent')
class TestBundledSmoke:
    @pytest.fixture(scope='class')
    def genome(self):
        from rectify.utils.genome import load_genome
        return load_genome(GENOME_FSA)

    @pytest.fixture(scope='class')
    def annotated(self):
        from rectify.core.consensus.consensus import load_annotated_junctions
        return load_annotated_junctions(str(ANNOTATION))

    def test_classify_bundled_bam(self, genome, annotated):
        results = classify_bam(str(VALIDATION_BAM), genome, annotated)
        assert len(results) > 0
        labels = {r.label for r in results}
        assert labels <= {LABEL_HIGH_CONFIDENCE, LABEL_TRIAGE}
        n_hc = sum(1 for r in results if r.label == LABEL_HIGH_CONFIDENCE)
        assert n_hc + sum(1 for r in results if r.label == LABEL_TRIAGE) == len(results)

    def test_triage_realign_end_to_end(self, genome, annotated, tmp_path):
        out_bam = tmp_path / 'triaged.bam'
        rows, stats = triage_realign_bam(
            str(VALIDATION_BAM), str(out_bam),
            genome=genome,
            annotated_junctions=annotated,
            realign=True,
            sort_and_index=True,
        )
        assert out_bam.exists()
        assert stats['classified'] == len(rows) > 0
        assert stats['high_confidence'] + stats['triaged'] == stats['classified']
        assert stats['accepted'] <= stats['realigned'] <= stats['junction_leg'] <= stats['triaged']

        # Read count preserved: triage swaps records, never adds or drops.
        def _count(path):
            with pysam.AlignmentFile(str(path), 'rb', check_sq=False) as bam:
                return sum(1 for _ in bam.fetch(until_eof=True))
        assert _count(out_bam) == _count(VALIDATION_BAM)

        # Every accepted read must be labelled + junction-leg in the TSV rows.
        by_id = {r['read_id']: r for r in rows}
        for r in rows:
            if r['accepted']:
                assert r['realigned']
                assert r['label'] == LABEL_TRIAGE
        assert set(by_id) == {r['read_id'] for r in rows}
