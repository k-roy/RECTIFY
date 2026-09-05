"""Station C's length pre-gate bound (ISSUE-013).

``multi_aligner.derive_max_intron`` takes ``2 x the LONGEST annotated intron``
and clamps to 500,000. That is right for an aligner's hard ``-G`` cap, which
must not amputate real biology, but it makes Station C's length pre-gate inert
on any organism with a long tail: GENCODE v48 chr5's longest annotated intron is
772,519 bp, so the rule saturates the clamp and NO junction can ever trip the
pre-gate — the one term meant to stop the physically impossible was dead on
every human run.

``derive_pool_gate_max_intron`` uses a high QUANTILE instead. Measured:

    annotation              n       max      p99.5 x 2   derive_max_intron
    R64 GFF3 (yeast)      378     2,483          5,000               5,000
    R64 GTF  (yeast)      378     2,483          5,000               5,000
    GENCODE v48 chr5   38,849   772,519        310,100   500,000 (clamped)

The yeast identity is not a coincidence to be trusted — it is asserted below.
"""

from pathlib import Path

import pytest

import rectify
from rectify.core.align.multi_aligner import derive_max_intron
from rectify.core.consensus.station_c import (
    PoolGateConfig, annotated_intron_lengths, derive_pool_gate_max_intron,
    pool_gate)

_SCER = Path(rectify.__file__).parent / 'data' / 'genomes' / 'saccharomyces_cerevisiae'
YEAST_GFF = _SCER / 'saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz'
YEAST_GTF = _SCER / 'saccharomyces_cerevisiae_R64-5-1_20240529.gtf'


# ---------------------------------------------------------------------------
# Yeast: byte-identical to the historical bound
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not YEAST_GFF.exists(), reason='bundled yeast GFF absent')
@pytest.mark.parametrize('path', [YEAST_GFF, YEAST_GTF])
def test_yeast_pre_gate_is_unchanged(path):
    """The whole point of a quantile is that it must not move yeast. With 378
    annotated introns the top 0.5% IS the 2,483 bp maximum, so p99.5 x 2
    reproduces the historical 5,000 exactly — same value the aligner derives."""
    if not path.exists():
        pytest.skip(f'{path.name} absent')
    assert derive_pool_gate_max_intron(path) == 5000
    assert derive_pool_gate_max_intron(path) == derive_max_intron(str(path))

    lens = annotated_intron_lengths(path)
    assert len(lens) == 378 and lens[-1] == 2483
    assert lens == sorted(lens)


# ---------------------------------------------------------------------------
# A human-shaped annotation: live bound, well under the ceiling
# ---------------------------------------------------------------------------

#: Tail shaped to the MEASURED GENCODE v48 chr5 distribution rather than to
#: whatever makes the test pass: there, 38,849 introns give p99 = 119,187,
#: p99.5 = 155,029, p99.9 = 261,984, max = 772,519. So ~0.5% of introns exceed
#: ~150 kb and only a handful reach the hundreds of kb. A fatter synthetic tail
#: (a decade of outliers in 1,000 introns) pushes p99.5 into the outliers and
#: clamps at 500,000 again — which is a fixture artifact, not the real shape.
_TAIL = [100_000, 120_000, 130_000, 140_000, 150_000,     # the p99.5 band
         300_000, 400_000, 500_000, 600_000, 772_519]     # the rare extremes


def _human_like_gtf(tmp_path):
    """A GENCODE-shaped GTF with a human long tail: 990 introns of 5 kb plus
    the ten in ``_TAIL``. 2x the LONGEST clamps at 500,000; 2x p99.5 must not."""
    p = tmp_path / 'humanlike.gtf'
    attrs = 'gene_id "G{0}"; transcript_id "T{0}";'
    lines = []
    pos = 1

    def add(i, gap):
        nonlocal pos
        a = attrs.format(i)
        lines.append(f'chr5\tT\texon\t{pos}\t{pos + 199}\t.\t+\t.\t{a}')
        lines.append(f'chr5\tT\texon\t{pos + 200 + gap}\t{pos + 400 + gap}\t.\t+\t.\t{a}')
        pos += gap + 20000

    for i in range(990):
        add(i, 5000)
    for j, gap in enumerate(_TAIL):
        add(2000 + j, gap)
    p.write_text('\n'.join(lines) + '\n')
    return p


def test_a_human_like_annotation_yields_a_live_bound(tmp_path):
    gtf = _human_like_gtf(tmp_path)
    lens = annotated_intron_lengths(gtf)
    assert len(lens) == 1000 and lens[-1] == 772_519

    # The aligner's rule saturates the ceiling: pre-gate dead.
    assert derive_max_intron(str(gtf)) == 500_000

    q = derive_pool_gate_max_intron(gtf)
    assert q < 500_000, (
        'the whole defect is a pre-gate pinned at the 500,000 ceiling, where no '
        'junction can ever trip it')
    assert 5000 < q < 400_000


def test_a_400kb_junction_demotes_under_the_quantile_gate(tmp_path):
    """End to end: the bound is not just a number, it changes the verdict."""
    gtf = _human_like_gtf(tmp_path)
    q = derive_pool_gate_max_intron(gtf)
    assert q < 400_000 <= 500_000, 'fixture must straddle the two bounds'

    genome = {'chr5': 'ACGT' * 250_000}        # 1 Mb, sequence irrelevant here
    cfg_q = PoolGateConfig(max_intron=q)
    cfg_old = PoolGateConfig(max_intron=derive_max_intron(str(gtf)))

    import pysam
    bam = tmp_path / 'big.bam'
    header = {'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'chr5', 'LN': 1_000_000}]}
    seq = genome['chr5']
    with pysam.AlignmentFile(str(bam), 'wb', header=header) as out:
        for i in range(2):
            a = pysam.AlignedSegment()
            a.query_name = f'r{i}'
            a.reference_id = 0
            a.reference_start = 1000
            a.mapping_quality = 60
            a.cigarstring = '50M400000N50M'
            a.query_sequence = seq[1000:1050] + seq[401050:401100]
            out.write(a)

    rows_q, _ = pool_gate(str(bam), genome, gtf, cfg=cfg_q)
    rows_old, _ = pool_gate(str(bam), genome, gtf, cfg=cfg_old)
    assert rows_q and rows_old, 'the 400 kb junction must be censused either way'
    assert rows_q[0]['over_max_intron'] == 1
    assert rows_q[0]['verdict'] == 'demote_orthogonal_evidence'
    # ... and under the aligner's saturated bound it is NOT flagged at all
    assert rows_old[0]['over_max_intron'] == 0


# ---------------------------------------------------------------------------
# Shape / degenerate inputs
# ---------------------------------------------------------------------------

def test_no_annotation_and_no_introns_fall_back(tmp_path):
    assert derive_pool_gate_max_intron(None) == 5000
    empty = tmp_path / 'empty.gtf'
    empty.write_text('# nothing here\n')
    assert derive_pool_gate_max_intron(empty) == 5000
    assert annotated_intron_lengths(empty) == []


def test_the_bound_is_clamped_and_rounded(tmp_path):
    """Clamped to the same [1000, 500000] window the aligner uses, so the
    pre-gate can never end up looser than the aligner's own cap."""
    tiny = tmp_path / 'tiny.gtf'
    tiny.write_text(
        'chr5\tT\texon\t1\t100\t.\t+\t.\tgene_id "G"; transcript_id "T";\n'
        'chr5\tT\texon\t131\t230\t.\t+\t.\tgene_id "G"; transcript_id "T";\n')
    assert annotated_intron_lengths(tiny) == [30]
    assert derive_pool_gate_max_intron(tiny) == 1000      # 60 -> clamped up

    assert derive_pool_gate_max_intron(YEAST_GFF, multiplier=1000) == 500_000


def test_explicit_intron_features_take_precedence(tmp_path):
    """Same precedence as derive_max_intron: an annotation carrying `intron`
    features is read from those, not from exon gaps."""
    p = tmp_path / 'both.gff3'
    p.write_text(
        '##gff-version 3\n'
        'chr5\tT\tintron\t101\t400\t.\t+\t.\tID=i1\n'
        'chr5\tT\texon\t1\t100\t.\t+\t.\tParent=t1\n'
        'chr5\tT\texon\t100001\t100100\t.\t+\t.\tParent=t1\n')
    assert annotated_intron_lengths(p) == [300]   # not the 99,900 exon gap
