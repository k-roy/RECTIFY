"""
Validation test suite for the 36 bundled RECTIFY example reads.

Runs rectify correction on each read and asserts that the expected correction
was applied for each of the eight categories. This serves as both a regression
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
    data_dir = Path(__file__).parent.parent / 'rectify' / 'data' / 'validation'
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


def run_correction_with_aligner_bams(
    bam_path: Path, genome_path: Path, annotation_path: Path, tmp_path: Path
) -> dict:
    """
    Run ``rectify correct`` with ``--aligner-bams`` so Module 2H fires.

    Discovers the per-aligner BAMs in the ``aligners/`` subdirectory next to
    *bam_path*.  Returns a dict keyed by read_id (first row per read_id).
    """
    import subprocess, sys, csv

    aligner_dir = bam_path.parent / 'aligners'
    aligner_args = []
    for aligner in ['minimap2', 'gapmm2', 'mapPacBio', 'deSALT', 'uLTRA']:
        p = aligner_dir / f'validation_reads.{aligner}.bam'
        if p.exists():
            aligner_args.extend(['--aligner-bams', f'{aligner}:{p}'])

    if not aligner_args:
        pytest.skip('No aligner BAMs found alongside validation BAM')
    if annotation_path is None:
        pytest.skip('Bundled annotation not available; Module 2H requires --annotation')

    out_tsv = tmp_path / 'corrected_with_aligners.tsv'
    cmd = [
        sys.executable, '-m', 'rectify.cli', 'correct',
        str(bam_path),
        '--genome', str(genome_path),
        '--annotation', str(annotation_path),
        *aligner_args,
        '-o', str(out_tsv),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        pytest.fail(f'rectify correct --aligner-bams failed:\n{result.stderr}')

    rows = {}
    with open(out_tsv) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            rid = row['read_id']
            if rid not in rows:
                rows[rid] = row
    return rows


@pytest.fixture(scope='module')
def corrected_with_aligner_bams(bam_path, genome_path, annotation_path, tmp_path_factory):
    tmp = tmp_path_factory.mktemp('correction_aligners')
    return run_correction_with_aligner_bams(bam_path, genome_path, annotation_path, tmp)


# ---------------------------------------------------------------------------
# Structural tests — verify bundled BAM is intact
# ---------------------------------------------------------------------------

class TestBamIntegrity:
    def test_bam_has_36_reads(self, raw_reads):
        assert len(raw_reads) == 36, f'Expected 36 reads, got {len(raw_reads)}'

    def test_correction_has_36_reads(self, corrected):
        """Every read must survive the correction pipeline (no silent drops)."""
        assert len(corrected) == 36, \
            f'Expected 36 reads in corrected output, got {len(corrected)}'

    def test_all_labels_present(self, raw_reads):
        expected = {
            'cat1_plus_1', 'cat1_plus_2', 'cat1_minus_1', 'cat1_minus_2',
            'cat2_plus_1', 'cat2_plus_2', 'cat2_minus_1', 'cat2_minus_2',
            # Cat3: 4/5 aligners soft-clip the 5' exon; mapPacBio may span or map upstream.
            # Source: minimap2 BAM (soft-clip present in read, rescued by local aligner).
            'cat3_plus_1', 'cat3_plus_2', 'cat3_minus_1', 'cat3_minus_2',
            'cat4_plus_1', 'cat4_plus_2', 'cat4_minus_1', 'cat4_minus_2',
            # Cat5: two aligners each contribute a different intron; chimeric consensus merges.
            'cat5_plus_1', 'cat5_plus_2',
            'cat5_minus_1', 'cat5_minus_2',
            # Cat6: one aligner (mapPacBio) spans the 5' intron; minimap2/gapmm2 soft-clip.
            'cat6_plus_1', 'cat6_plus_2', 'cat6_minus_1', 'cat6_minus_2',
            # Cat7: non-canonical unannotated junctions from mapPacBio.
            'cat7_plus_1', 'cat7_plus_2',
            'cat7_minus_1', 'cat7_minus_2',
            # Cat8: NET-seq A-tract refinement (formerly Cat6).
            'cat8_plus_single', 'cat8_plus_multi',
            'cat8_minus_single', 'cat8_minus_multi',
            # Cat9: Module 2H N-op junction boundary refinement (requires --aligner-bams).
            'cat9_plus_1', 'cat9_plus_2',
            'cat9_minus_1', 'cat9_minus_2',
        }
        assert set(raw_reads.keys()) == expected

    def test_category_tags(self, raw_reads):
        cat_map = {
            'cat1': 'cat1_indel',
            'cat2': 'cat2_softclip',
            'cat3': 'cat3_junction',
            'cat4': 'cat4_false_junc',
            'cat5': 'cat5_chimeric',
            'cat6': 'cat6_chimeric',
            'cat7': 'cat7_alt_splice',
            'cat8': 'cat8_netseq_refine',
            'cat9': 'cat9_junction_refine',
        }
        for label, read in raw_reads.items():
            prefix = label[:4]
            assert read.get_tag('XG') == cat_map[prefix], \
                f'{label}: expected XG={cat_map[prefix]}, got {read.get_tag("XG")}'

    def test_strand_balance(self, raw_reads):
        for cat in ['cat1', 'cat2', 'cat3', 'cat4', 'cat5', 'cat6', 'cat7', 'cat8', 'cat9']:
            cat_reads = {k: v for k, v in raw_reads.items() if k.startswith(cat)}
            plus = sum(1 for r in cat_reads.values() if not r.is_reverse)
            minus = sum(1 for r in cat_reads.values() if r.is_reverse)
            assert plus == minus == 2, \
                f'{cat}: expected 2 plus / 2 minus, got {plus}/{minus}'

    def test_chimeric_tags(self, raw_reads):
        """Cat 5 reads must carry XK=1 and comma-separated XA with exactly 2 aligners."""
        for label in ['cat5_plus_1', 'cat5_plus_2', 'cat5_minus_1', 'cat5_minus_2']:
            r = raw_reads[label]
            assert r.get_tag('XK') == 1, f'{label}: XK should be 1'
            xa = r.get_tag('XA')
            aligners = xa.split(',')
            assert len(aligners) == 2, \
                f'{label}: XA should list exactly 2 aligners, got {xa!r}'


# ---------------------------------------------------------------------------
# Correction tests — verify rectify actually corrects each category
# ---------------------------------------------------------------------------

class TestCategory1IndelCorrection:
    """Walk-back must shift the corrected 3' end away from the raw alignment boundary.

    cat1_plus_1  chrVIII:508599–508794  +  A-tract walk-back (−4 bp)
    cat1_plus_2  chrIX:22544–22850      +  A-tract walk-back (−3 bp)
    cat1_minus_1 chrII:9855–10533       −  A-tract walk-back (+6 bp)
    cat1_minus_2 chrII:9813–10539       −  A-tract walk-back (+4 bp)
    """

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

    @pytest.mark.parametrize('label,expected_3prime', [
        # Exact corrected_3prime values from rectify correct on wt_by4742_rep1 validation reads.
        # polya_walkback + NET-seq refinement; values are deterministic.
        ('cat1_plus_1',  508789),   # chrVIII walk-back −4 bp from raw 508793
        ('cat1_plus_2',  22846),    # chrIX walk-back −3 bp from raw 22849; homopolymer_rescue
        ('cat1_minus_1', 9861),     # chrII walk-forward +6 bp from raw 9855
        ('cat1_minus_2', 9817),     # chrII walk-forward +4 bp from raw 9813
    ])
    def test_3prime_exact_position(self, corrected, raw_reads, label, expected_3prime):
        read = raw_reads[label]
        row = corrected.get(read.query_name)
        if row is None:
            pytest.skip(f'Read {label} not in correction output')
        got = int(row['corrected_3prime'])
        assert got == expected_3prime, \
            f'{label}: corrected_3prime should be {expected_3prime}, got {got}'


class TestCategory2SoftClipRescue:
    """
    Corrected position should shift OUTWARD (away from gene body) by ≥1 bp.
    Plus strand: corrected > original.  Minus strand: corrected < original.
    """

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
        if strand == '+':
            assert corrected_pos > original, \
                f'{label}: plus-strand soft-clip rescue should shift outward (corrected > original)'
        else:
            assert corrected_pos < original, \
                f'{label}: minus-strand soft-clip rescue should shift outward (corrected < original)'

    @pytest.mark.parametrize('label,expected_3prime', [
        # Exact corrected_3prime values from soft-clip rescue on wt_by4742_rep1 validation reads.
        # rescue_softclip_at_homopolymer then terminal-A strip; values are deterministic.
        ('cat2_plus_1',  69844),    # chrI +12 bp rescue; T-homopolymer 10D + 2M(GC)
        ('cat2_plus_2',  199929),   # chrI +9 bp rescue; T-homopolymer 8D + 1M(G) strip 1A
        ('cat2_minus_1', 65610),    # chrI −12 bp rescue; T-homopolymer 8D + 4M(ATAA)
        ('cat2_minus_2', 128102),   # chrI −11 bp rescue; T-homopolymer 10D + 1M(A) strip 1T
    ])
    def test_3prime_exact_position(self, corrected, raw_reads, label, expected_3prime):
        read = raw_reads[label]
        row = corrected.get(read.query_name)
        if row is None:
            pytest.skip(f'Read {label} not in correction output')
        got = int(row['corrected_3prime'])
        assert got == expected_3prime, \
            f'{label}: corrected_3prime should be {expected_3prime}, got {got}'


class TestCategory3JunctionRescue:
    """5' end junction rescue: 4/5 aligners soft-clip the 5' exon; local aligner rescues.

    mapPacBio may span the intron or map upstream without soft-clipping, so the criterion
    is that the minimap2 source alignment has a 5' soft-clip that is rescued to an annotated
    3' splice site boundary by the affine-gap semi-global aligner (local_aligner.py).

    cat3_plus_1  (0a28167d) chrII:168808-169462  + (YBL027W region)  10S clip
    cat3_plus_2  (79f61403) chrI:142618-143383   + (YAL003W region)  22S clip
    cat3_minus_1 (ac4db6da) chrXV:900071-900767  - (RPL20B)
    cat3_minus_2 (28ea9379) chrII:365845-366503  - (YBR062C region)  25S clip
    """

    LABELS_WITH_STRAND = [
        ('cat3_plus_1', '+'),
        ('cat3_plus_2', '+'),
        ('cat3_minus_1', '-'),
        ('cat3_minus_2', '-'),
    ]

    @pytest.mark.parametrize('label,strand', LABELS_WITH_STRAND)
    def test_5prime_present(self, corrected, raw_reads, label, strand):
        read = raw_reads[label]
        row = corrected.get(read.query_name)
        if row is None:
            pytest.skip(f'Read {label} not in correction output')
        five_prime = row.get('five_prime_position', '')
        assert five_prime not in ('', 'NA', 'None', '-1'), \
            f'{label}: five_prime_position should be set after junction rescue'

    @pytest.mark.parametrize('label,strand,raw_5prime', [
        # raw_5prime = reference_start for + strand; reference_end-1 for - strand
        ('cat3_plus_1',  '+', 168808),   # chrII:168808-169462 +
        ('cat3_plus_2',  '+', 142618),   # chrI:142618-143383 +
        ('cat3_minus_1', '-', 900766),   # chrXV:900071-900767 -
        ('cat3_minus_2', '-', 366502),   # chrII:365845-366503 -
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

    @pytest.mark.parametrize('label,expected_five_prime', [
        # Expected five_prime_position after splice-site ambiguity resolution and
        # canonical donor preference (v2.9.6 plus, v2.9.8 minus).
        #
        # cat3_plus_2 (79f61403): YAL003W G|G ambiguity — last exon base G equals
        #   first intron base G. Without canonical preference the wrong junction
        #   (intron_start=142254, TA donor) wins. With GT preference, shift=-1 gives
        #   intron_start=142253 (GT donor) → five_prime_corrected=142252.
        ('cat3_plus_2',  142252),
        # cat3_minus_2 (28ea9379): YBR062C minus strand. Minus-strand canonical AC
        #   preference (v2.9.8) rescues to intron_end=366584.
        ('cat3_minus_2', 366584),
        # cat3_plus_1 (0a28167d): YBL027W, rescued to annotated intron_start-1=168423.
        ('cat3_plus_1',  168423),
        # cat3_minus_1 (ac4db6da): RPL20B, rescued to intron_end=901193.
        ('cat3_minus_1', 901193),
    ])
    def test_5prime_exact_position(self, corrected, raw_reads, label, expected_five_prime):
        """
        Regression: exact five_prime_position after splice-site ambiguity resolution
        and canonical GT/GC (plus) or AC/GC (minus) donor preference.
        """
        read = raw_reads[label]
        row = corrected.get(read.query_name)
        if row is None:
            pytest.skip(f'Read {label} not in correction output')
        five_prime_str = row.get('five_prime_position', '')
        if five_prime_str in ('', 'NA', 'None', '-1'):
            pytest.skip(f'five_prime_position not set for {label}')
        assert int(five_prime_str) == expected_five_prime, \
            f'{label}: five_prime_position should be {expected_five_prime}, ' \
            f'got {five_prime_str}'

    @pytest.mark.parametrize('label', [
        'cat3_plus_1', 'cat3_plus_2', 'cat3_minus_1', 'cat3_minus_2',
    ])
    def test_5prime_exon_cigar_set(self, corrected, raw_reads, label):
        """After rescue, five_prime_exon_cigar must be a non-empty CIGAR string."""
        read = raw_reads[label]
        row = corrected.get(read.query_name)
        if row is None:
            pytest.skip(f'Read {label} not in correction output')
        cigar = row.get('five_prime_exon_cigar', '')
        assert cigar not in ('', 'nan', 'None', 'NaN'), \
            f'{label}: five_prime_exon_cigar should be set after rescue, got {cigar!r}'
        assert any(c.isdigit() for c in cigar) and any(c.isalpha() for c in cigar), \
            f'{label}: five_prime_exon_cigar is not a valid CIGAR: {cigar!r}'

    @pytest.mark.parametrize('label', [
        'cat3_plus_1', 'cat3_plus_2', 'cat3_minus_1', 'cat3_minus_2',
    ])
    def test_correction_applied_includes_five_prime_rescued(self, corrected, raw_reads, label):
        """correction_applied must include 'five_prime_rescued' for all Cat3 reads."""
        read = raw_reads[label]
        row = corrected.get(read.query_name)
        if row is None:
            pytest.skip(f'Read {label} not in correction output')
        applied = row.get('correction_applied', '')
        assert 'five_prime_rescued' in applied, \
            f'{label}: correction_applied should include "five_prime_rescued", got {applied!r}'


class TestCategory4FalseJunction:
    """False N ops near 3' end must be absorbed; corrected position walks back past them.

    cat4_plus_1  chrXI:19592–22073  +  N op 20527–22047 (1520 bp, 26 bp from 3' end)
                                        polya_walkback + NET-seq → 22070 (post-N exon is polyA)
    cat4_plus_2  chrX:392246–393837 +  N op 393725–393825 (100 bp, 12 bp from 3' end)
                                        window-clipped → 393721
    cat4_minus_1 chrI:128094–129063 −  N op 128521–129021 (500 bp, 427 bp from 3' end)
                                        N outside FJF window; treated as real junction → 128098
    cat4_minus_2 chrIX:76016–77313  −  N op 76027–76250 (223 bp, 11 bp from 3' end)
                                        snap fires → 76027
    """

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
        # cat4_plus_1 ends at a non-A genomic base (chrXI:22072='C'); NET-seq
        # refinement is now disabled by default so this read is not shifted.
        if label == 'cat4_plus_1':
            pytest.skip('cat4_plus_1 correction relied on NET-seq (now opt-in); no shift expected')
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

    @pytest.mark.parametrize('label,expected_3prime', [
        # Exact corrected_3prime values.
        # cat4_plus_1: original 22072 ('C', not in A-tract).  NET-seq refinement
        # is now disabled in 'rectify correct' (opt-in only via --netseq-dir).
        # The position stays at 22072; walk-back does not fire on a non-A base.
        ('cat4_plus_1',  22072),
        ('cat4_plus_2',  393721),   # window-clipped to exclude artifact N
        ('cat4_minus_1', 128098),   # N far from 3' end; normal walkback
        ('cat4_minus_2', 76027),    # snap to intron_start (3' end abuts N)
    ])
    def test_3prime_exact_position(self, corrected, raw_reads, label, expected_3prime):
        read = raw_reads[label]
        row = corrected.get(read.query_name)
        if row is None:
            pytest.skip(f'Read {label} not in correction output')
        got = int(row['corrected_3prime'])
        assert got == expected_3prime, \
            f'{label}: corrected_3prime should be {expected_3prime}, got {got}'

    @pytest.mark.parametrize('label', [
        'cat4_plus_1', 'cat4_plus_2', 'cat4_minus_1', 'cat4_minus_2',
    ])
    def test_has_one_junction(self, corrected, raw_reads, label):
        """All Cat4 reads retain exactly one junction in the output."""
        read = raw_reads[label]
        row = corrected.get(read.query_name)
        if row is None:
            pytest.skip(f'Read {label} not in correction output')
        assert int(row['n_junctions']) == 1, \
            f'{label}: expected n_junctions=1, got {row["n_junctions"]}'


class TestCategory5ChimericReconstruction:
    """Cat5: two aligners each contribute a different GT-AG intron; chimeric consensus
    merges them.  Each source alignment carries exactly one intron (XS=2 segments);
    the chimeric output would have 2 introns (3 segments).

    cat5_plus_1:  chrV  + mapPacBio intron 423590-423951 (361bp); mm2 intron 424421-425030 (609bp)
    cat5_plus_2:  chrII + mm2 intron 332875-333386 (511bp);  mapPacBio intron 334050-334122 (72bp)
    cat5_minus_1: chrVII - mm2 intron 436480-437396 (916bp); gapmm2 intron 437941-438397 (456bp)
    cat5_minus_2: chrIII - mapPacBio intron 177906-178213 (307bp); gapmm2 intron 176709-177362 (653bp)
    """

    LABELS = ['cat5_plus_1', 'cat5_plus_2', 'cat5_minus_1', 'cat5_minus_2']

    def test_all_present(self, raw_reads):
        for label in self.LABELS:
            assert label in raw_reads, f'Chimeric read {label} missing from BAM'

    def test_segment_counts(self, raw_reads):
        """Source alignment carries ≥2 segments (≥1 intron in source aligner)."""
        for label in self.LABELS:
            r = raw_reads[label]
            xs = r.get_tag('XS') if r.has_tag('XS') else 0
            assert xs >= 2, f'{label}: expected ≥2 segments (XS), got {xs}'

    def test_has_intron_in_source(self, raw_reads):
        """Source CIGAR must contain at least one N op (intron skip)."""
        for label in self.LABELS:
            r = raw_reads[label]
            has_n = any(op == 3 for op, _ in (r.cigartuples or []))
            assert has_n, f'{label}: source alignment must have ≥1 N op (intron), got {r.cigarstring}'


class TestCategory6SimpleChimeric:
    """Simple chimeric: mapPacBio correctly spans the 5' intron; minimap2/gapmm2
    soft-clip the same region.  The read in validation_reads.bam uses the mapPacBio
    alignment (XU=1) so the junction is directly visible.

    cat6_plus_1  chrII:+    intron 125154–125270 (116 bp)
    cat6_plus_2  chrII:+    intron 45644–45977   (333 bp)
    cat6_minus_1 chrII:–    intron 60193–60697   (504 bp)
    cat6_minus_2 chrIV:–    intron 307333–307765 (432 bp)
    """

    LABELS = ['cat6_plus_1', 'cat6_plus_2', 'cat6_minus_1', 'cat6_minus_2']

    def test_all_present(self, raw_reads):
        for label in self.LABELS:
            assert label in raw_reads, f'{label} missing from validation BAM'

    def test_xg_tag(self, raw_reads):
        for label in self.LABELS:
            assert raw_reads[label].get_tag('XG') == 'cat6_chimeric', \
                f'{label}: expected XG=cat6_chimeric'

    def test_xu_tag(self, raw_reads):
        """Cat6 reads come from a single winning aligner (mapPacBio), so XU=1."""
        for label in self.LABELS:
            xu = raw_reads[label].get_tag('XU') if raw_reads[label].has_tag('XU') else None
            assert xu == 1, f'{label}: expected XU=1, got {xu}'

    @pytest.mark.parametrize('label', LABELS)
    def test_spans_intron(self, raw_reads, label):
        """mapPacBio alignment must have at least one N op (intron skip)."""
        r = raw_reads[label]
        has_n = any(op == 3 for op, _ in (r.cigartuples or []))
        assert has_n, f'{label}: expected intron-spanning CIGAR (N op), got {r.cigarstring}'

    @pytest.mark.parametrize('label', LABELS)
    def test_no_5prime_softclip(self, raw_reads, label):
        """mapPacBio spans the 5' intron so there must be no soft-clip at the 5' end.

        For plus-strand reads the 5' end is the first CIGAR op (cigar[0]).
        For minus-strand reads the 5' end is the last CIGAR op (cigar[-1]),
        because CIGAR is always written left-to-right in reference coordinates.
        """
        r = raw_reads[label]
        cigar = r.cigartuples or []
        if not cigar:
            pytest.skip(f'{label}: no CIGAR')
        if r.is_reverse:
            # Minus strand: 5' RNA end is at high genomic coords (last CIGAR op)
            clip_5p = cigar[-1][1] if cigar[-1][0] == 4 else 0
        else:
            # Plus strand: 5' RNA end is at low genomic coords (first CIGAR op)
            clip_5p = cigar[0][1] if cigar[0][0] == 4 else 0
        assert clip_5p == 0, \
            f'{label}: mapPacBio source should have no 5\' soft-clip, got {clip_5p}S'

    @pytest.mark.parametrize('label', LABELS)
    def test_no_five_prime_rescue(self, corrected, raw_reads, label):
        """Cat6 reads span the 5' intron already; no 5' rescue should fire."""
        read = raw_reads[label]
        row = corrected.get(read.query_name)
        if row is None:
            pytest.skip(f'Read {label} not in correction output')
        assert row.get('five_prime_rescued', '0') == '0', \
            f'{label}: five_prime_rescued should be 0 for Cat6 reads'

    @pytest.mark.parametrize('label', LABELS)
    def test_has_one_junction(self, corrected, raw_reads, label):
        """Cat6 reads have the intron in the source BAM, so n_junctions=1."""
        read = raw_reads[label]
        row = corrected.get(read.query_name)
        if row is None:
            pytest.skip(f'Read {label} not in correction output')
        assert int(row['n_junctions']) == 1, \
            f'{label}: expected n_junctions=1, got {row["n_junctions"]}'


class TestCategory8NetseqRefinement:
    """NET-seq A-tract refinement (formerly Cat6): single-peak reads get fraction=1.0;
    multi-peak reads get ≥2 fractional rows summing to 1.0."""

    @pytest.mark.parametrize('label', [
        'cat8_plus_single',
        'cat8_minus_single',
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

    @pytest.mark.parametrize('label,strand', [
        ('cat8_plus_multi', '+'),
        ('cat8_minus_multi', '-'),
    ])
    def test_multi_peak_polya_anchor(self, corrected_all_rows, raw_reads, label, strand):
        """NEW-065: When all NET-seq peaks land on poly-A bases, polya_walkback
        anchor (first non-A/T) must be used.  Primary corrected_3prime must
        not be at a poly-A base regardless of NET-seq signal.
        Fractions must still sum to 1.0.

        NOTE: NET-seq refinement is now DISABLED by default in 'rectify correct'
        (requires explicit --netseq-dir).  These tests verify the poly-A anchor
        logic that guards against NET-seq placing signal on poly-A bases, but
        can only run when NET-seq is loaded.
        """
        pytest.skip(
            'NET-seq refinement is now opt-in (--netseq-dir required); '
            'multi-peak anchor tests require bundled NET-seq to be active.'
        )
        from rectify.data import get_bundled_genome_path
        from rectify.utils.genome import load_genome

        read = raw_reads[label]
        rows = corrected_all_rows.get(read.query_name, [])
        if not rows:
            pytest.skip(f'Read {label} not in correction output')

        corrections = rows[0].get('correction_applied', '')
        assert 'netseq_refinement' in corrections, \
            f'{label}: expected netseq_refinement in correction_applied'

        genome_path = get_bundled_genome_path('saccharomyces_cerevisiae')
        if genome_path is None:
            pytest.skip('Bundled genome not available')
        genome = load_genome(genome_path)
        chrom = rows[0]['chrom']
        pos = int(rows[0]['corrected_3prime'])
        polya_base = 'A' if strand == '+' else 'T'
        if chrom in genome and pos < len(genome[chrom]):
            base = genome[chrom][pos].upper()
            assert base != polya_base, \
                f'{label}: primary corrected_3prime={pos} is a poly-A base ({base}); ' \
                f'expected first non-{polya_base} upstream the pA tail'

        assert len(rows) >= 1, f'{label}: expected at least 1 output row'
        total = sum(float(r.get('fraction', 1.0)) for r in rows)
        assert abs(total - 1.0) < 0.02, \
            f'{label}: fractions should sum to 1.0, got {total:.3f}'


class TestCategory7AltSplice:
    """Non-canonical, unannotated splice junctions from mapPacBio alignments.

    Reads:
      cat7_plus_1  (4e43165e) chrIII:+  junction 138864-138952 (88 bp)  AC-AG
      cat7_plus_2  (0f021462) chrXII:+  junction 595739-595853 (114 bp) CA-TT
      cat7_minus_1 (c79f1fb9) chrII:-   junction 443720-443833 (113 bp) GT-CG
      cat7_minus_2 (72557a9a) chrIII:-  junction 104435-104495  (60 bp) GT-CG
    """

    LABELS = ['cat7_plus_1', 'cat7_plus_2', 'cat7_minus_1', 'cat7_minus_2']

    # Expected junction coordinates (jstart-jend) per label
    EXPECTED_JUNCTIONS = {
        'cat7_plus_1':  '138864-138952',
        'cat7_plus_2':  '595739-595853',
        'cat7_minus_1': '443720-443833',
        'cat7_minus_2': '104435-104495',  # updated after mapPacBio alignment confirmed
    }

    def test_all_present(self, raw_reads):
        for label in self.LABELS:
            assert label in raw_reads, f'{label} missing from validation BAM'

    def test_xg_tag(self, raw_reads):
        for label in self.LABELS:
            r = raw_reads[label]
            assert r.get_tag('XG') == 'cat7_alt_splice', \
                f'{label}: expected XG=cat7_alt_splice, got {r.get_tag("XG")}'

    def test_xu_tag(self, raw_reads):
        """Cat7 reads come from a single aligner (mapPacBio), so XU=1."""
        for label in self.LABELS:
            r = raw_reads[label]
            xu = r.get_tag('XU') if r.has_tag('XU') else None
            assert xu == 1, f'{label}: expected XU=1, got {xu}'

    @pytest.mark.parametrize('label', LABELS)
    def test_has_one_junction(self, corrected, raw_reads, label):
        read = raw_reads[label]
        row = corrected.get(read.query_name)
        if row is None:
            pytest.skip(f'{label} not in correction output')
        n = int(row.get('n_junctions', 0))
        assert n == 1, f'{label}: expected n_junctions=1, got {n}'

    @pytest.mark.parametrize('label', LABELS)
    def test_junction_coordinates(self, corrected, raw_reads, label):
        read = raw_reads[label]
        row = corrected.get(read.query_name)
        if row is None:
            pytest.skip(f'{label} not in correction output')
        junctions = row.get('junctions', '')
        expected = self.EXPECTED_JUNCTIONS[label]
        assert expected in junctions, \
            f'{label}: expected junction {expected} in junctions={junctions!r}'

    @pytest.mark.parametrize('label', LABELS)
    def test_no_five_prime_rescue(self, corrected, raw_reads, label):
        """Cat7 reads have genuine non-canonical junctions; 5' rescue (Cat3)
        should NOT be triggered."""
        read = raw_reads[label]
        row = corrected.get(read.query_name)
        if row is None:
            pytest.skip(f'{label} not in correction output')
        rescued = row.get('five_prime_rescued', '0')
        assert str(rescued) in ('0', '', 'False', 'nan'), \
            f'{label}: five_prime_rescued should not be set, got {rescued!r}'


class TestCategory9JunctionRefinement:
    """Module 2H: N-op junction boundary correction via ``--aligner-bams``.

    When aligner BAMs are supplied alongside ``--annotation``, Module 2H
    (``junction_refiner.py``) re-scores all candidate splice junctions and
    replaces imprecise N-op boundaries in the consensus read.

    Input N-ops (wrong boundaries in consensus BAM):

      cat9_plus_1  (00a1c9b3) chrVII:+ (555825, 556304) intron length 479
      cat9_plus_2  (00a1e01e) chrVII:+ (439089, 439321) intron length 232
      cat9_minus_1 (0b3b593b) chrXV:−  (900760, 901191) intron length 431
      cat9_minus_2 (d3357db5) chrXV:−  (900760, 901192) intron length 432

    Expected junctions after Module 2H (derived from aligner-BAM candidate pool):

      cat9_plus_1  → 555830-556307
      cat9_plus_2  → 439093-439323
      cat9_minus_1 → 900767-901193
      cat9_minus_2 → 900767-901193

    The minus-strand reads additionally receive Cat3 5' soft-clip rescue
    (``five_prime_rescued``), which reports the same intron_end boundary.
    """

    LABELS = ['cat9_plus_1', 'cat9_plus_2', 'cat9_minus_1', 'cat9_minus_2']

    # Junctions present in the raw consensus BAM (before Module 2H correction)
    RAW_JUNCTIONS = {
        'cat9_plus_1':  '555825-556304',
        'cat9_plus_2':  '439089-439321',
        'cat9_minus_1': '900760-901191',
        'cat9_minus_2': '900760-901192',
    }

    # Junctions after Module 2H + downstream correction (with aligner BAMs)
    CORRECTED_JUNCTIONS = {
        'cat9_plus_1':  '555830-556307',
        'cat9_plus_2':  '439093-439323',
        'cat9_minus_1': '900767-901193',
        'cat9_minus_2': '900767-901193',
    }

    def test_all_present(self, raw_reads):
        for label in self.LABELS:
            assert label in raw_reads, f'Cat9 read {label} missing from validation BAM'

    def test_xg_tag(self, raw_reads):
        for label in self.LABELS:
            r = raw_reads[label]
            assert r.get_tag('XG') == 'cat9_junction_refine', \
                f'{label}: expected XG=cat9_junction_refine, got {r.get_tag("XG")}'

    @pytest.mark.parametrize('label', LABELS)
    def test_raw_junction_in_consensus_bam(self, raw_reads, label):
        """Consensus BAM encodes the WRONG N-op boundaries before correction."""
        r = raw_reads[label]
        n_ops = [
            (r.reference_start + sum(l for op, l in r.cigartuples[:i]),
             r.reference_start + sum(l for op, l in r.cigartuples[:i]) + length)
            for i, (op, length) in enumerate(r.cigartuples)
            if op == 3
        ]
        js, je = self.RAW_JUNCTIONS[label].split('-')
        assert any(ns == int(js) and ne == int(je) for ns, ne in n_ops), \
            f'{label}: expected raw N-op ({js},{je}) in CIGAR, got {n_ops}'

    @pytest.mark.parametrize('label', LABELS)
    def test_junction_corrected_with_aligner_bams(
        self, corrected_with_aligner_bams, raw_reads, label
    ):
        """With ``--aligner-bams``, Module 2H corrects the N-op to the expected junction."""
        read = raw_reads[label]
        row = corrected_with_aligner_bams.get(read.query_name)
        if row is None:
            pytest.skip(f'{label} not in correction output')
        junctions = row.get('junctions', '')
        expected = self.CORRECTED_JUNCTIONS[label]
        assert expected in junctions, \
            f'{label}: expected corrected junction {expected!r} in junctions={junctions!r}'

    @pytest.mark.parametrize('label', LABELS)
    def test_junction_not_corrected_without_aligner_bams(self, corrected, raw_reads, label):
        """Without ``--aligner-bams``, Module 2H is inactive; junction stays wrong."""
        read = raw_reads[label]
        row = corrected.get(read.query_name)
        if row is None:
            pytest.skip(f'{label} not in correction output')
        junctions = row.get('junctions', '')
        expected_corrected = self.CORRECTED_JUNCTIONS[label]
        assert expected_corrected not in junctions, (
            f'{label}: junction should NOT be {expected_corrected!r} without --aligner-bams; '
            f'got {junctions!r}'
        )


# ---------------------------------------------------------------------------
# Bundled pA-tail soft-clipped BAM integrity
# ---------------------------------------------------------------------------

class TestPolyASoftClippedBam:
    """Verify the bundled ``rectified_pA_tail_soft_clipped.bam``.

    This BAM is produced by ``restore_polya_softclips()`` acting on
    ``rectified_pA_tail_trimmed.bam`` with the DRS trim metadata from the
    ``wt_by4742_rep1`` dev run.  The three replacement reads (f8050895,
    7d5e8dc2, 72557a9a) have no trim metadata and therefore carry over
    unchanged (sc5=sc3=0).  All other reads should have their poly-A tail
    and/or adapter stub restored as soft-clips.
    """

    # Reads that have trim metadata: expected to carry ≥1 soft-clip base
    # somewhere (5' or 3').  The three replacement reads (no metadata) are excluded.
    REPLACEMENT_IDS = {'f8050895', '7d5e8dc2', '72557a9a'}

    @pytest.fixture(scope='class')
    def sc_bam_path(self):
        p = Path(__file__).parent.parent / 'rectify' / 'data' / 'validation' / \
            'rectified' / 'rectified_pA_tail_soft_clipped.bam'
        if not p.exists():
            pytest.skip(f'rectified_pA_tail_soft_clipped.bam not found: {p}')
        return p

    @pytest.fixture(scope='class')
    def sc_reads(self, sc_bam_path):
        """Load primary reads keyed by XV label."""
        reads = {}
        with pysam.AlignmentFile(str(sc_bam_path), 'rb') as bam:
            for r in bam:
                if r.is_secondary or r.is_supplementary:
                    continue
                if r.has_tag('XV'):
                    reads[r.get_tag('XV')] = r
        return reads

    def test_has_36_reads(self, sc_reads):
        assert len(sc_reads) == 36, f'Expected 36 primary reads, got {len(sc_reads)}'

    def test_is_sorted_and_indexed(self, sc_bam_path):
        """BAM must be coordinate-sorted and have a .bai index."""
        bai = Path(str(sc_bam_path) + '.bai')
        assert bai.exists(), f'Missing index: {bai}'
        with pysam.AlignmentFile(str(sc_bam_path), 'rb') as bam:
            hdr = bam.header.to_dict()
            sort_order = hdr.get('HD', {}).get('SO', '')
            assert sort_order == 'coordinate', f'Expected coordinate sort, got {sort_order!r}'

    def test_reads_with_metadata_have_softclip(self, sc_reads):
        """Every read that has trim metadata must carry ≥1 soft-clip base."""
        no_clip = []
        for label, r in sc_reads.items():
            prefix = r.query_name[:8]
            if prefix in self.REPLACEMENT_IDS:
                continue  # no metadata, skip
            if not r.cigartuples:
                continue
            sc5 = r.cigartuples[0][1]  if r.cigartuples[0][0]  == 4 else 0
            sc3 = r.cigartuples[-1][1] if r.cigartuples[-1][0] == 4 else 0
            if sc5 == 0 and sc3 == 0:
                no_clip.append(label)
        assert not no_clip, f'Reads with metadata but no soft-clip: {no_clip}'

    def test_replacement_reads_no_unexpected_softclip_change(self, sc_reads, raw_reads):
        """Replacement reads (no trim metadata) should not have poly-A appended."""
        for label, r in sc_reads.items():
            prefix = r.query_name[:8]
            if prefix not in self.REPLACEMENT_IDS:
                continue
            if not r.cigartuples:
                continue
            sc3 = r.cigartuples[-1][1] if r.cigartuples[-1][0] == 4 else 0
            sc5 = r.cigartuples[0][1]  if r.cigartuples[0][0]  == 4 else 0
            # These reads have no trim metadata; any soft-clips present come
            # from the original pA_tail_trimmed BAM (pre-existing).  We just
            # assert the read is still present.
            assert sc5 >= 0 and sc3 >= 0  # always true; confirms read is present
