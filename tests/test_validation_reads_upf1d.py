"""Validation suite for the upf1Δ parallel example reads (validation_reads_upf1d.*).

Parallel to tests/test_validation_reads.py but sourced from upf1Δ rep1 DRS
(poly-A trimmed, 5-aligner panel built on Sherlock 2026-05-29). Purely additive:
the wt 36-read set and its tests are untouched. Same correct-first pipeline
(`rectify correct` per aligner → merge_corrected_tsvs), keyed by XV label.

Run with:
    pytest tests/test_validation_reads_upf1d.py -v
"""
import csv
import subprocess
import sys
from pathlib import Path

import pytest
import pysam

ALIGNERS = ['minimap2', 'gapmm2', 'mapPacBio', 'deSALT', 'uLTRA']


# ---------------------------------------------------------------------------
# Helpers (upf1d-specific bundle paths)
# ---------------------------------------------------------------------------

def get_validation_bam() -> Path:
    data_dir = Path(__file__).parent.parent / 'rectify' / 'data' / 'validation'
    bam = data_dir / 'validation_reads_upf1d.bam'
    if not bam.exists():
        pytest.skip(f'upf1d validation BAM not found: {bam}')
    return bam


def load_reads(bam_path: Path) -> dict:
    bam = pysam.AlignmentFile(str(bam_path), 'rb')
    reads = {}
    for r in bam:
        if r.is_secondary or r.is_supplementary:
            continue
        try:
            reads[r.get_tag('XV')] = r
        except KeyError:
            pass
    bam.close()
    return reads


def _discover_per_aligner_bams(bam_path: Path) -> dict:
    aligner_dir = bam_path.parent / 'aligners'
    out = {}
    for a in ALIGNERS:
        p = aligner_dir / f'validation_reads_upf1d.{a}.bam'
        if p.exists():
            out[a] = p
    return out


def _load_all_rows(tsv_path: Path) -> dict:
    rows: dict = {}
    with open(tsv_path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            rows.setdefault(row['read_id'], []).append(row)
    return rows


def _run_correct_first_pipeline(per_aligner_bams, genome_path, annotation_path,
                                tmp_dir, *, with_module_2h):
    """Mirror of the wt fixture: per-aligner correct (sequential, M1 memory
    rule) → merge_corrected_tsvs. Returns (merged_tsv, per_aligner_corrected)."""
    per_tsv, per_cbam, failed = {}, {}, []
    for aligner, bam in per_aligner_bams.items():
        out_tsv = tmp_dir / f'{aligner}.tsv'
        out_bam = tmp_dir / f'{aligner}.corrected.bam'
        cmd = [sys.executable, '-m', 'rectify.cli', 'correct', str(bam),
               '--genome', str(genome_path), '-o', str(out_tsv),
               '--write-corrected-bam', str(out_bam), '--emit-merged-tsv']
        if annotation_path is not None:
            cmd += ['--annotation', str(annotation_path)]
        if with_module_2h:
            for an, ap in per_aligner_bams.items():
                cmd += ['--aligner-bams', f'{an}:{ap}']
        res = subprocess.run(cmd, capture_output=True)
        if res.returncode != 0:
            failed.append((aligner, res.stderr.decode(errors='replace')))
            continue
        per_tsv[aligner] = out_tsv
        if out_bam.exists():
            per_cbam[aligner] = str(out_bam)
    if failed:
        msgs = '\n'.join(f'--- {a} ---\n{e}' for a, e in failed)
        pytest.fail(f'rectify correct failed for {len(failed)} aligner(s):\n{msgs}')

    from rectify.core.consensus.corrected_consensus import merge_corrected_tsvs
    merged = tmp_dir / 'merged_corrected.tsv'
    merge_corrected_tsvs(per_tsv, merged, per_aligner_corrected_bams=per_cbam)
    return merged, per_cbam


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope='module')
def bam_path():
    return get_validation_bam()


@pytest.fixture(scope='module')
def raw_reads(bam_path):
    return load_reads(bam_path)


@pytest.fixture(scope='module')
def genome_path():
    from rectify.data import get_bundled_genome_path
    p = get_bundled_genome_path('saccharomyces_cerevisiae')
    if p is None:
        pytest.skip('Bundled S. cerevisiae genome not available')
    return p


@pytest.fixture(scope='module')
def annotation_path():
    from rectify.data import get_bundled_annotation_path
    return get_bundled_annotation_path('saccharomyces_cerevisiae')


@pytest.fixture(scope='module')
def _correction_outputs(bam_path, genome_path, annotation_path, tmp_path_factory):
    """Run the correct-first pipeline twice (no Module 2H, and with it).
    Returns (no2h_rows, with2h_rows) as dict[read_id, list[row]]."""
    per_aligner_bams = _discover_per_aligner_bams(bam_path)
    if not per_aligner_bams:
        pytest.skip('No per-aligner BAMs found alongside validation_reads_upf1d.bam')
    tmp_no = tmp_path_factory.mktemp('upf1d_no2h')
    no2h, _ = _run_correct_first_pipeline(
        per_aligner_bams, genome_path, annotation_path, tmp_no, with_module_2h=False)
    with2h_rows = None
    if annotation_path is not None:
        tmp_w = tmp_path_factory.mktemp('upf1d_2h')
        w, _ = _run_correct_first_pipeline(
            per_aligner_bams, genome_path, annotation_path, tmp_w, with_module_2h=True)
        with2h_rows = _load_all_rows(w)
    return _load_all_rows(no2h), with2h_rows


@pytest.fixture(scope='module')
def corrected_all_rows(_correction_outputs):
    """All TSV rows per read_id (no Module 2H)."""
    return _correction_outputs[0]


@pytest.fixture(scope='module')
def corrected(_correction_outputs):
    """First (primary) row per read_id, no Module 2H."""
    return {rid: rs[0] for rid, rs in _correction_outputs[0].items()}


@pytest.fixture(scope='module')
def corrected_with_aligner_bams(_correction_outputs):
    """First row per read_id from the Module-2H run (cat9 junction refinement)."""
    w = _correction_outputs[1]
    if w is None:
        pytest.skip('annotation unavailable; Module 2H not runnable')
    return {rid: rs[0] for rid, rs in w.items()}


def _row_for(corrected, raw_reads, label):
    rid = raw_reads[label].query_name
    return corrected.get(rid)


# ---------------------------------------------------------------------------
# Structural
# ---------------------------------------------------------------------------

_CAT_XG = {
    'cat1': 'cat1_indel', 'cat2': 'cat2_softclip', 'cat3': 'cat3_junction',
    'cat5': 'cat5_multi_intron', 'cat6': 'cat6_chimeric', 'cat7': 'cat7_alt_splice',
    'cat8': 'cat8_atract', 'cat9': 'cat9_junction_refine',
}


class TestBundleIntegrity:
    def test_32_reads(self, raw_reads):
        assert len(raw_reads) == 32

    def test_eight_categories_four_each(self, raw_reads):
        from collections import Counter
        c = Counter(lbl.split('_')[0] for lbl in raw_reads)
        assert dict(c) == {k: 4 for k in _CAT_XG}, dict(c)

    def test_strand_balance(self, raw_reads):
        plus = sum(1 for lbl in raw_reads if '_plus' in lbl)
        minus = sum(1 for lbl in raw_reads if '_minus' in lbl)
        assert plus == 16 and minus == 16

    def test_category_tags(self, raw_reads):
        for lbl, r in raw_reads.items():
            assert r.get_tag('XG') == _CAT_XG[lbl.split('_')[0]], lbl


# ---------------------------------------------------------------------------
# Category 3 — 5' junction rescue (exercises the rescue DP cap + early-exit fix)
# raw_5prime = reference_start (+) / reference_end-1 (-) of the stored winner aln.
# ---------------------------------------------------------------------------

class TestCategory3JunctionRescue:
    # label -> (strand, raw_5prime, expected_five_prime_position)
    CASES = {
        'cat3_plus_1':  ('+', 282665, 282139),
        'cat3_plus_2':  ('+', 148284, 148193),
        'cat3_minus_1': ('-', 307332, 307765),
        'cat3_minus_2': ('-', 712725, 713155),
    }

    @pytest.mark.parametrize('label', list(CASES))
    def test_five_prime_rescued_applied(self, corrected, raw_reads, label):
        row = _row_for(corrected, raw_reads, label)
        assert row is not None, f'{label}: no corrected row'
        assert 'five_prime_rescued' in row.get('correction_applied', ''), \
            f'{label}: correction_applied={row.get("correction_applied")!r}'

    @pytest.mark.parametrize('label', list(CASES))
    def test_five_prime_position_set_and_moved(self, corrected, raw_reads, label):
        _, raw5, _ = self.CASES[label]
        row = _row_for(corrected, raw_reads, label)
        fp = row.get('five_prime_position', '')
        assert fp not in ('', 'NA', None), f'{label}: five_prime_position unset'
        assert int(fp) != raw5, f'{label}: 5p not moved from raw {raw5}'

    @pytest.mark.parametrize('label', list(CASES))
    def test_five_prime_exact_position(self, corrected, raw_reads, label):
        _, _, expected = self.CASES[label]
        row = _row_for(corrected, raw_reads, label)
        assert int(row['five_prime_position']) == expected, \
            f'{label}: five_prime_position {row.get("five_prime_position")} != {expected}'


# ---------------------------------------------------------------------------
# Category 1 — A-tract indel / polya_walkback (3' end)
# ---------------------------------------------------------------------------

class TestCategory1IndelCorrection:
    # label -> (strand, expected_corrected_3prime)
    CASES = {
        'cat1_plus_1':  ('+', 708336),
        'cat1_plus_2':  ('+', 1239967),
        'cat1_minus_1': ('-', 1162861),
        'cat1_minus_2': ('-', 381917),
    }

    @pytest.mark.parametrize('label', list(CASES))
    def test_indel_correction_applied(self, corrected, raw_reads, label):
        row = _row_for(corrected, raw_reads, label)
        assert row is not None, f'{label}: no corrected row'
        assert 'indel_correction' in row.get('correction_applied', ''), \
            f'{label}: correction_applied={row.get("correction_applied")!r}'

    @pytest.mark.parametrize('label', list(CASES))
    def test_corrected_3prime_exact(self, corrected, raw_reads, label):
        _, expected = self.CASES[label]
        row = _row_for(corrected, raw_reads, label)
        assert int(row['corrected_3prime']) == expected, \
            f'{label}: corrected_3prime {row.get("corrected_3prime")} != {expected}'


# ---------------------------------------------------------------------------
# Category 2 — soft-clip rescue at an under-called homopolymer (3' end)
# ---------------------------------------------------------------------------

class TestCategory2SoftClipRescue:
    CASES = {  # label -> expected corrected_3prime
        'cat2_plus_1': 723682, 'cat2_plus_2': 66060,
        'cat2_minus_1': 233917, 'cat2_minus_2': 236091,
    }

    @pytest.mark.parametrize('label', list(CASES))
    def test_softclip_rescue_applied(self, corrected, raw_reads, label):
        row = _row_for(corrected, raw_reads, label)
        assert row is not None, f'{label}: no corrected row'
        assert 'softclip_rescue' in row.get('correction_applied', ''), \
            f'{label}: correction_applied={row.get("correction_applied")!r}'

    @pytest.mark.parametrize('label', list(CASES))
    def test_corrected_3prime_exact(self, corrected, raw_reads, label):
        row = _row_for(corrected, raw_reads, label)
        assert int(row['corrected_3prime']) == self.CASES[label], \
            f'{label}: corrected_3prime {row.get("corrected_3prime")} != {self.CASES[label]}'


# ---------------------------------------------------------------------------
# Category 5 — multi-intron reads (2 annotated introns; structural)
# The chimeric "no single aligner reconstructs" variant is rare in this dataset;
# upf1Δ cat5 reads are genuine 2-intron transcripts (Xg>=3). Structural only.
# ---------------------------------------------------------------------------

class TestCategory5MultiIntron:
    LABELS = ['cat5_plus_1', 'cat5_plus_2', 'cat5_minus_1', 'cat5_minus_2']

    def test_all_present(self, raw_reads):
        for lbl in self.LABELS:
            assert lbl in raw_reads, f'{lbl} missing'

    @pytest.mark.parametrize('label', LABELS)
    def test_segment_count(self, raw_reads, label):
        """Xg = segment count (>=2 means >=1 intron). Stored alignment carries
        the multi-intron structure."""
        r = raw_reads[label]
        xg = r.get_tag('Xg') if r.has_tag('Xg') else 0
        assert xg >= 2, f'{label}: expected Xg>=2, got {xg}'

    @pytest.mark.parametrize('label', LABELS)
    def test_stored_has_intron(self, raw_reads, label):
        r = raw_reads[label]
        assert any(op == 3 for op, _ in (r.cigartuples or [])), \
            f'{label}: stored alignment should contain an N op'


# ---------------------------------------------------------------------------
# Category 6 — single-aligner intron span (mapPacBio spans the 5' intron)
# ---------------------------------------------------------------------------

class TestCategory6SimpleChimeric:
    EXPECTED_INTRON = {  # label -> (start, end)
        'cat6_plus_1': (45644, 45977), 'cat6_plus_2': (45644, 45977),
        'cat6_minus_1': (89132, 89440), 'cat6_minus_2': (592416, 592768),
    }

    @pytest.mark.parametrize('label', list(EXPECTED_INTRON))
    def test_xm_tag(self, raw_reads, label):
        r = raw_reads[label]
        xm = r.get_tag('Xm') if r.has_tag('Xm') else None
        assert xm is not None and xm >= 1, f'{label}: expected Xm>=1, got {xm}'

    @pytest.mark.parametrize('label', list(EXPECTED_INTRON))
    def test_stored_spans_intron_no_5p_clip(self, raw_reads, label):
        r = raw_reads[label]
        cig = r.cigartuples or []
        assert any(op == 3 for op, _ in cig), f'{label}: stored should have N op'
        clip5 = (cig[-1][1] if (r.is_reverse and cig[-1][0] == 4) else
                 cig[0][1] if (not r.is_reverse and cig[0][0] == 4) else 0)
        assert clip5 == 0, f"{label}: stored mapPacBio should have no 5' soft-clip, got {clip5}S"

    @pytest.mark.parametrize('label', list(EXPECTED_INTRON))
    def test_corrected_intron_present(self, corrected, raw_reads, label):
        row = _row_for(corrected, raw_reads, label)
        if row is None:
            pytest.skip(f'{label} not in correction output')
        s, e = self.EXPECTED_INTRON[label]
        assert f'{s}-{e}' in (row.get('junctions', '') or '').split(';'), \
            f'{label}: expected intron {s}-{e} in {row.get("junctions")!r}'


# ---------------------------------------------------------------------------
# Category 7 — non-canonical, unannotated splice junctions (mapPacBio)
# ---------------------------------------------------------------------------

class TestCategory7AltSplice:
    EXPECTED_JUNCTION = {
        'cat7_plus_1': '246695-246734', 'cat7_plus_2': '248644-248687',
        'cat7_minus_1': '1489311-1489390', 'cat7_minus_2': '476768-476896',
    }

    @pytest.mark.parametrize('label', list(EXPECTED_JUNCTION))
    def test_has_junction(self, corrected, raw_reads, label):
        row = _row_for(corrected, raw_reads, label)
        if row is None:
            pytest.skip(f'{label} not in correction output')
        assert int(row.get('n_junctions', 0) or 0) >= 1, \
            f'{label}: expected n_junctions>=1, got {row.get("n_junctions")}'

    @pytest.mark.parametrize('label', list(EXPECTED_JUNCTION))
    def test_junction_coordinates(self, corrected, raw_reads, label):
        row = _row_for(corrected, raw_reads, label)
        if row is None:
            pytest.skip(f'{label} not in correction output')
        assert self.EXPECTED_JUNCTION[label] in (row.get('junctions', '') or ''), \
            f'{label}: expected {self.EXPECTED_JUNCTION[label]} in {row.get("junctions")!r}'

    @pytest.mark.parametrize('label', list(EXPECTED_JUNCTION))
    def test_no_five_prime_rescue(self, corrected, raw_reads, label):
        row = _row_for(corrected, raw_reads, label)
        if row is None:
            pytest.skip(f'{label} not in correction output')
        assert 'five_prime_rescued' not in row.get('correction_applied', ''), \
            f'{label}: five_prime_rescued unexpectedly applied'


# ---------------------------------------------------------------------------
# Category 8 — A-tract 3' end (single-peak, fraction=1.0)
# ---------------------------------------------------------------------------

class TestCategory8Atract:
    LABELS = ['cat8_plus_1', 'cat8_plus_2', 'cat8_minus_1', 'cat8_minus_2']

    @pytest.mark.parametrize('label', LABELS)
    def test_single_peak_fraction(self, corrected_all_rows, raw_reads, label):
        rid = raw_reads[label].query_name
        rows = corrected_all_rows.get(rid, [])
        if not rows:
            pytest.skip(f'{label} not in correction output')
        assert len(rows) == 1, f'{label}: expected 1 output row, got {len(rows)}'
        frac = float(rows[0].get('fraction', 1.0) or 1.0)
        assert abs(frac - 1.0) < 0.01, f'{label}: fraction {frac} != 1.0'

    @pytest.mark.parametrize('label', LABELS)
    def test_atract_ambiguity_applied(self, corrected, raw_reads, label):
        row = _row_for(corrected, raw_reads, label)
        assert 'atract_ambiguity' in row.get('correction_applied', ''), \
            f'{label}: correction_applied={row.get("correction_applied")!r}'


# ---------------------------------------------------------------------------
# Category 9 — Module-2H N-op junction boundary refinement (--aligner-bams)
# ---------------------------------------------------------------------------

class TestCategory9JunctionRefinement:
    REFINED_INTRON = {  # label -> (start, end) annotated boundary after refinement
        'cat9_plus_1': '439093-439323', 'cat9_plus_2': '110879-110948',
        'cat9_minus_1': '414753-415259', 'cat9_minus_2': '407027-407122',
    }

    @pytest.mark.parametrize('label', list(REFINED_INTRON))
    def test_refined_to_annotated_boundary(self, corrected_with_aligner_bams, raw_reads, label):
        rid = raw_reads[label].query_name
        row = corrected_with_aligner_bams.get(rid)
        if row is None:
            pytest.skip(f'{label} not in Module-2H correction output')
        assert self.REFINED_INTRON[label] in (row.get('junctions', '') or ''), \
            f'{label}: refined intron {self.REFINED_INTRON[label]} not in {row.get("junctions")!r}'
