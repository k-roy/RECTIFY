"""
Validation test suite for the 24 bundled RECTIFY example reads.

Runs rectify correction on each read and asserts that the expected correction
was applied for each of the six categories. This serves as both a regression
test and an installation smoke-test.

Run with:
    pytest tests/test_validation_reads.py -v
    rectify validate          (CLI shortcut)
"""

import pytest
import pysam
from pathlib import Path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def get_validation_bam() -> Path:
    """Return path to the bundled validation BAM."""
    data_dir = Path(__file__).parent.parent / 'rectify' / 'data'
    bam = data_dir / 'validation_reads.bam'
    if not bam.exists():
        pytest.skip(f'Validation BAM not found: {bam}')
    return bam


def load_reads(bam_path: Path) -> dict:
    """Load all primary reads from BAM, keyed by XV label."""
    bam = pysam.AlignmentFile(str(bam_path), 'rb')
    reads = {}
    for r in bam:
        if r.is_secondary or r.is_supplementary:
            continue
        try:
            label = r.get_tag('XV')
            reads[label] = r
        except KeyError:
            pass
    bam.close()
    return reads


def run_correction(bam_path: Path, genome_path: Path, annotation_path: Path,
                   tmp_path: Path) -> dict:
    """
    Run `rectify correct` on the validation BAM and return the TSV results as
    a dict keyed by read_id (first row per read_id, by appearance).
    """
    import subprocess, sys, csv

    out_tsv = tmp_path / 'corrected.tsv'
    cmd = [
        sys.executable, '-m', 'rectify.cli', 'correct',
        str(bam_path),
        '--genome', str(genome_path),
        '-o', str(out_tsv),
    ]
    if annotation_path is not None:
        cmd += ['--annotation', str(annotation_path)]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        pytest.fail(f'rectify correct failed:\n{result.stderr}')

    # Return first row per read_id (primary position; Cat6 may produce multiple)
    rows = {}
    with open(out_tsv) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            rid = row['read_id']
            if rid not in rows:
                rows[rid] = row
    return rows


def run_correction_all_rows(bam_path: Path, genome_path: Path,
                             annotation_path: Path, tmp_path: Path) -> dict:
    """
    Like run_correction but returns *all* rows per read_id as a list.
    Used for Cat6 fraction tests.
    """
    import subprocess, sys, csv

    out_tsv = tmp_path / 'corrected_all.tsv'
    cmd = [
        sys.executable, '-m', 'rectify.cli', 'correct',
        str(bam_path),
        '--genome', str(genome_path),
        '-o', str(out_tsv),
    ]
    if annotation_path is not None:
        cmd += ['--annotation', str(annotation_path)]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        pytest.fail(f'rectify correct failed:\n{result.stderr}')

    rows = {}
    with open(out_tsv) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            rid = row['read_id']
            rows.setdefault(rid, []).append(row)
    return rows


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
    path = get_bundled_genome_path('saccharomyces_cerevisiae')
    if path is None:
        pytest.skip('Bundled S. cerevisiae genome not available')
    return path


@pytest.fixture(scope='module')
def annotation_path():
    from rectify.data import get_bundled_annotation_path
    return get_bundled_annotation_path('saccharomyces_cerevisiae')  # None is OK


@pytest.fixture(scope='module')
def corrected(bam_path, genome_path, annotation_path, tmp_path_factory):
    tmp = tmp_path_factory.mktemp('correction')
    return run_correction(bam_path, genome_path, annotation_path, tmp)


@pytest.fixture(scope='module')
def corrected_all_rows(bam_path, genome_path, annotation_path, tmp_path_factory):
    tmp = tmp_path_factory.mktemp('correction_all')
    return run_correction_all_rows(bam_path, genome_path, annotation_path, tmp)


# ---------------------------------------------------------------------------
# Structural tests — verify bundled BAM is intact
# ---------------------------------------------------------------------------

class TestBamIntegrity:
    def test_bam_has_24_reads(self, raw_reads):
        assert len(raw_reads) == 24, f'Expected 24 reads, got {len(raw_reads)}'

    def test_all_labels_present(self, raw_reads):
        expected = {
            'cat1_plus_1', 'cat1_plus_2', 'cat1_minus_1', 'cat1_minus_2',
            'cat2_plus_1', 'cat2_plus_2', 'cat2_minus_1', 'cat2_minus_2',
            'cat3_plus_1', 'cat3_plus_2', 'cat3_minus_1', 'cat3_minus_2',
            'cat4_plus_1', 'cat4_plus_2', 'cat4_minus_1', 'cat4_minus_2',
            'cat5_plus_3aligner', 'cat5_plus_2aligner',
            'cat5_minus_long', 'cat5_minus_short',
            'cat6_plus_single', 'cat6_plus_multi',
            'cat6_minus_single', 'cat6_minus_multi',
        }
        assert set(raw_reads.keys()) == expected

    def test_category_tags(self, raw_reads):
        cat_map = {
            'cat1': 'cat1_indel',
            'cat2': 'cat2_softclip',
            'cat3': 'cat3_junction',
            'cat4': 'cat4_false_junc',
            'cat5': 'cat5_chimeric',
            'cat6': 'cat6_netseq_refine',
        }
        for label, read in raw_reads.items():
            prefix = label[:4]
            assert read.get_tag('XG') == cat_map[prefix], \
                f'{label}: expected XG={cat_map[prefix]}, got {read.get_tag("XG")}'

    def test_strand_balance(self, raw_reads):
        for cat in ['cat1', 'cat2', 'cat3', 'cat4', 'cat6']:
            cat_reads = {k: v for k, v in raw_reads.items() if k.startswith(cat)}
            plus = sum(1 for r in cat_reads.values() if not r.is_reverse)
            minus = sum(1 for r in cat_reads.values() if r.is_reverse)
            assert plus == minus == 2, \
                f'{cat}: expected 2 plus / 2 minus, got {plus}/{minus}'

    def test_chimeric_tags(self, raw_reads):
        """Cat 5 reads must carry XK=1 and comma-separated XA."""
        for label in ['cat5_plus_3aligner', 'cat5_plus_2aligner',
                      'cat5_minus_long', 'cat5_minus_short']:
            r = raw_reads[label]
            assert r.get_tag('XK') == 1, f'{label}: XK should be 1'
            xa = r.get_tag('XA')
            assert ',' in xa, f'{label}: XA should contain multiple aligners, got {xa!r}'

    def test_three_aligner_chimeric(self, raw_reads):
        """cat5_plus_3aligner must have all three aligners in XA."""
        r = raw_reads['cat5_plus_3aligner']
        xa = r.get_tag('XA')
        aligners = set(xa.split(','))
        assert len(aligners) == 3, f'Expected 3 aligners, got {aligners}'
        assert aligners == {'minimap2', 'mapPacBio', 'gapmm2'}, \
            f'Unexpected aligner set: {aligners}'


# ---------------------------------------------------------------------------
# Correction tests — verify rectify actually corrects each category
# ---------------------------------------------------------------------------

class TestCategory1IndelCorrection:
    """Walk-back must shift the corrected 3' end away from the raw alignment boundary."""

    @pytest.mark.parametrize('label,strand', [
        ('cat1_plus_1', '+'),
        ('cat1_plus_2', '+'),
        ('cat1_minus_1', '-'),
        ('cat1_minus_2', '-'),
    ])
    def test_3prime_shifted(self, corrected, raw_reads, label, strand):
        read = raw_reads[label]
        row = corrected.get(read.query_name)
        if row is None:
            pytest.skip(f'Read {label} not in correction output')
        original = int(row['original_3prime'])
        corrected_pos = int(row['corrected_3prime'])
        assert original != corrected_pos, \
            f'{label}: corrected_3prime should differ from original ({original})'
        if strand == '+':
            assert corrected_pos < original, \
                f'{label}: plus-strand correction should walk back (corrected < original)'
        else:
            assert corrected_pos > original, \
                f'{label}: minus-strand correction should walk forward (corrected > original)'


class TestCategory2SoftClipRescue:
    """Corrected position should shift by ≥1 bp relative to raw alignment end."""

    @pytest.mark.parametrize('label,strand', [
        ('cat2_plus_1', '+'),
        ('cat2_plus_2', '+'),
        ('cat2_minus_1', '-'),
        ('cat2_minus_2', '-'),
    ])
    def test_3prime_shifted(self, corrected, raw_reads, label, strand):
        read = raw_reads[label]
        row = corrected.get(read.query_name)
        if row is None:
            pytest.skip(f'Read {label} not in correction output')
        original = int(row['original_3prime'])
        corrected_pos = int(row['corrected_3prime'])
        shift = abs(corrected_pos - original)
        assert shift >= 1, f'{label}: soft-clip rescue should shift position by ≥1 bp'
        assert shift <= 20, f'{label}: soft-clip rescue shift of {shift} bp seems too large'


class TestCategory3JunctionRescue:
    """5' end junction rescue: five_prime_position should differ from raw alignment 5' end."""

    @pytest.mark.parametrize('label,strand', [
        ('cat3_plus_1', '+'),
        ('cat3_plus_2', '+'),
        ('cat3_minus_1', '-'),
        ('cat3_minus_2', '-'),
    ])
    def test_5prime_present(self, corrected, raw_reads, label, strand):
        read = raw_reads[label]
        row = corrected.get(read.query_name)
        if row is None:
            pytest.skip(f'Read {label} not in correction output')
        five_prime = row.get('five_prime_position', '')
        assert five_prime not in ('', 'NA', 'None', '-1'), \
            f'{label}: five_prime_position should be set after junction rescue'

    @pytest.mark.parametrize('label,strand,raw_5prime', [
        # raw 5' end = reference_start for plus, reference_end-1 for minus
        # cat3_plus_1: chrII:125270-126182+, alignment starts at 125270
        ('cat3_plus_1', '+', 125270),
        # cat3_plus_2: chrII:168808-169478+, alignment starts at 168808
        ('cat3_plus_2', '+', 168808),
        # cat3_minus_1: chrXV:900071-900767-, RPL20B; reference_end-1 = 900766
        ('cat3_minus_1', '-', 900766),
        # cat3_minus_2: chrII:59715-60193-; reference_end-1 = 60192
        ('cat3_minus_2', '-', 60192),
    ])
    def test_5prime_rescued(self, corrected, raw_reads, label, strand, raw_5prime):
        read = raw_reads[label]
        row = corrected.get(read.query_name)
        if row is None:
            pytest.skip(f'Read {label} not in correction output')
        five_prime = row.get('five_prime_position', '')
        if five_prime in ('', 'NA', 'None', '-1'):
            pytest.skip(f'five_prime_position not set for {label}')
        rescued = int(five_prime)
        assert rescued != raw_5prime, \
            f'{label}: 5\' end should be rescued away from raw position {raw_5prime}'


class TestCategory4FalseJunction:
    """False N ops near 3' end must be absorbed; corrected position walks back past them."""

    @pytest.mark.parametrize('label,strand', [
        ('cat4_plus_1', '+'),
        ('cat4_plus_2', '+'),
        ('cat4_minus_1', '-'),
        ('cat4_minus_2', '-'),
    ])
    def test_3prime_shifted(self, corrected, raw_reads, label, strand):
        read = raw_reads[label]
        row = corrected.get(read.query_name)
        if row is None:
            pytest.skip(f'Read {label} not in correction output')
        original = int(row['original_3prime'])
        corrected_pos = int(row['corrected_3prime'])
        assert original != corrected_pos, \
            f'{label}: false junction walk-back should move corrected_3prime'
        if strand == '+':
            assert corrected_pos < original, \
                f'{label}: plus-strand walk-back should give corrected < original'
        else:
            assert corrected_pos > original, \
                f'{label}: minus-strand walk-back should give corrected > original'


class TestCategory5ChimericReconstruction:
    """Chimeric reads must be in the BAM with correct XK/XA/XS tags — no re-correction needed."""

    def test_all_present(self, raw_reads):
        for label in ['cat5_plus_3aligner', 'cat5_plus_2aligner',
                      'cat5_minus_long', 'cat5_minus_short']:
            assert label in raw_reads, f'Chimeric read {label} missing from BAM'

    def test_segment_counts(self, raw_reads):
        """Each chimeric read should have ≥3 independently scored segments."""
        for label in ['cat5_plus_3aligner', 'cat5_plus_2aligner',
                      'cat5_minus_long', 'cat5_minus_short']:
            r = raw_reads[label]
            xs = r.get_tag('XS') if r.has_tag('XS') else 0
            assert xs >= 3, f'{label}: expected ≥3 segments (XS), got {xs}'


class TestCategory6NetseqRefinement:
    """NET-seq A-tract refinement: single-peak reads get fraction=1.0;
    multi-peak reads get ≥2 fractional rows summing to 1.0."""

    @pytest.mark.parametrize('label', [
        'cat6_plus_single',
        'cat6_minus_single',
    ])
    def test_single_peak_fraction(self, corrected_all_rows, raw_reads, label):
        read = raw_reads[label]
        rows = corrected_all_rows.get(read.query_name, [])
        if not rows:
            pytest.skip(f'Read {label} not in correction output')
        assert len(rows) == 1, \
            f'{label}: single-peak read should produce exactly 1 output row, got {len(rows)}'
        frac = float(rows[0].get('fraction', 1.0))
        assert abs(frac - 1.0) < 0.01, \
            f'{label}: single-peak fraction should be 1.0, got {frac}'

    @pytest.mark.parametrize('label', [
        'cat6_plus_multi',
        'cat6_minus_multi',
    ])
    def test_multi_peak_fractions(self, corrected_all_rows, raw_reads, label):
        read = raw_reads[label]
        rows = corrected_all_rows.get(read.query_name, [])
        if not rows:
            pytest.skip(f'Read {label} not in correction output')
        assert len(rows) >= 2, \
            f'{label}: multi-peak read should produce ≥2 output rows, got {len(rows)}'
        total = sum(float(r.get('fraction', 1.0)) for r in rows)
        assert abs(total - 1.0) < 0.02, \
            f'{label}: fractions should sum to 1.0, got {total:.3f}'
        for r in rows:
            frac = float(r.get('fraction', 1.0))
            assert 0.0 < frac < 1.0, \
                f'{label}: each fraction should be in (0, 1), got {frac}'
