"""Human RNA004 DRS regression guard — Module 2H on the Sumner chr5 panel.

The yeast validation set cannot see the defects human data exposed: ~95 % of
intron-bearing yeast genes have ONE intron, so a walker that never advanced past
an N op looked correct, and yeast chromosome names are exactly the ones the
arabic->roman fallback was written for.

This module runs the Module 2H stage the way `rectify correct` runs it — contigs
registered first, annotation, pool, refine — over TWO datasets, and asserts the
FP/FN properties that were violated before 2026-09-05:

  FP  no annotated canonical junction is moved off annotation;
      no non-canonical N-op is created;
      the observed junction pool holds no coordinate that exists in no read.
  FN  the beneficial corrections (novel/drifted -> annotated) still happen.

Every test runs once per dataset (see DATASETS): the 100-read adversarial panel
and the 1,000-read unselected hold-out.  The hold-out is the one that earns
trust — the panel's reads were chosen BECAUSE RECTIFY changed them, so passing
on the panel alone shows only that the specific failures were addressed.  Both
of the FP classes fixed in this series were caught by one dataset and invisible
in the other, in both directions.

The BAMs and the chr5 reference are git-ignored (a 184 MB FASTA and a 47 MB
GTF), so the module skips when the reference is absent and each dataset skips
individually when its BAM is.  Point RECTIFY_SUMNER_PANEL_DIR at the directory
to run it elsewhere.

Author: Kevin R. Roy (agent S1)
"""

import os
import sys
from pathlib import Path

import pytest

RECTIFY_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(RECTIFY_ROOT))

PANEL_DIR = Path(os.environ.get(
    'RECTIFY_SUMNER_PANEL_DIR',
    '/Users/kevinroy/work/rectify/dev/sumner_misplaced_panel_20260904',
))
PANEL_BAM = PANEL_DIR / '175_panel_orig.bam'
HOLDOUT_BAM = PANEL_DIR / 'holdout' / 'chr5_holdout1k.bam'
REF_FA = PANEL_DIR / 'ref' / 'chr5.fa'
REF_GTF = PANEL_DIR / 'ref' / 'gencode.v48.basic.chr5.gtf'

# Two datasets, and the difference between them matters. The PANEL is
# adversarial: 100 reads selected BECAUSE `rectify correct` changed their
# junctions, so it over-represents every failure mode. The HOLD-OUT is 1,000
# random primary chr5 reads with the panel reads excluded — nothing about it was
# chosen by looking at RECTIFY's output, so it is the number that generalizes.
# A fix that only moves the panel has not been shown to help.
DATASETS = {'panel': PANEL_BAM, 'holdout': HOLDOUT_BAM}

pytestmark = pytest.mark.skipif(
    not (REF_FA.exists() and REF_GTF.exists()),
    reason=f"chr5 reference not present under {PANEL_DIR} (git-ignored data)",
)

# Human canonical splice dinucleotides, read on the genomic PLUS strand: the
# minus-strand form of each is listed alongside (GT-AG / CT-AC, GC-AG / CT-GC,
# AT-AC / GT-AT — the last pair is the U12 class).
CANONICAL = {('GT', 'AG'), ('CT', 'AC'),
             ('GC', 'AG'), ('CT', 'GC'),
             ('AT', 'AC'), ('GT', 'AT')}

_N = 3
_SAM_REF_CONSUMING = {0, 2, 3, 7, 8}

# Baselines recorded 2026-09-05 at d8b27e8 (the full fix series). MAX_* are
# ceilings that must not be exceeded; MIN_* are floors that must not be lost.
# Raising a ceiling or lowering a floor is a policy change and needs a reason in
# the commit message — that is the whole point of writing them down.
# `max_phantom` is the MEASURED residual, not a round number with headroom: both
# datasets sit at the D-abutting-N count exactly (panel 4 of 13 observed-only
# pool junctions, hold-out 5 of 46), so any new phantom coordinate fails here.
BASELINE = {
    #            harmful       beneficial          phantom pool junctions
    'panel':   {'max_harmful': 0, 'min_beneficial': 7,  'max_phantom': 4},
    'holdout': {'max_harmful': 0, 'min_beneficial': 12, 'max_phantom': 5},
}


def _n_ops(read):
    pos = read.reference_start
    out = []
    for op, ln in read.cigartuples or []:
        if op == _N:
            out.append((pos, pos + ln))
        if op in _SAM_REF_CONSUMING:
            pos += ln
    return out


def _motif(seq, s, e):
    return (seq[s:s + 2].upper(), seq[e - 2:e].upper())


@pytest.fixture(scope="module", params=sorted(DATASETS))
def refined(request, tmp_path_factory):
    """Run Module 2H over one dataset; return everything the assertions need."""
    dataset = request.param
    bam_path = DATASETS[dataset]
    if not bam_path.exists():
        pytest.skip(f"{dataset} BAM not present: {bam_path} (git-ignored data)")
    import pysam
    from rectify.utils.genome import load_genome, register_genome_contigs_from_fasta
    from rectify.core.consensus.consensus import load_annotated_junctions
    from rectify.core.splice.junction_refiner import (
        build_junction_pool, refine_bam_junctions,
    )
    from rectify.core.commands.correct_command import select_penalty_tables

    # Ordering is part of what is under test (ISSUE-001): contigs first.
    register_genome_contigs_from_fasta(str(REF_FA))
    annot = load_annotated_junctions(str(REF_GTF))
    pool, annot_set = build_junction_pool([str(bam_path)], annot)
    genome = load_genome(str(REF_FA))
    # The same tables `rectify correct --organism homo_sapiens` now selects
    # (ISSUE-005) — the empirical HP costs are what put a 1 nt acceptor slide
    # inside the noise floor where the annotated-canonical gate can hold it.
    jpt, spt = select_penalty_tables({'organism': 'homo_sapiens'})
    assert jpt is not None, "bundled human DRS penalty table should be found"

    out_bam = str(tmp_path_factory.mktemp(dataset) / "refined.bam")
    refine_bam_junctions(
        input_bam=str(bam_path), output_bam=out_bam,
        aligner_bams=[str(bam_path)], annotated_junctions=annot, genome=genome,
        prebuilt_junction_pool=pool, prebuilt_annotated_set=annot_set,
        penalty_table_path=jpt, str_penalty_table_path=spt,
        sort_and_index=True, n_workers=1,
    )

    before, after, observed = {}, {}, set()
    with pysam.AlignmentFile(str(bam_path)) as bam:
        for r in bam:
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            before[r.query_name] = (r.reference_name, _n_ops(r))
            observed.update((r.reference_name, s, e) for s, e in _n_ops(r))
    with pysam.AlignmentFile(out_bam) as bam:
        for r in bam:
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            after[r.query_name] = _n_ops(r)

    annot3 = {(str(j[0]), int(j[1]), int(j[2])) for j in annot if len(j) >= 3}
    return {
        'dataset': dataset, 'baseline': BASELINE[dataset],
        'before': before, 'after': after, 'observed': observed,
        'annot3': annot3, 'pool': pool, 'seq': genome['chr5'],
    }


def _moves(refined):
    """(read, chrom, old, new) for every N-op whose coordinates changed."""
    out = []
    for name, (chrom, old) in refined['before'].items():
        new = refined['after'].get(name)
        if new is None or len(new) != len(old):
            continue
        for o, n in zip(old, new):
            if o != n:
                out.append((name, chrom, o, n))
    return out


# ---------------------------------------------------------------------------
# Chromosome names (ISSUE-001)
# ---------------------------------------------------------------------------

def test_pool_is_keyed_by_the_reads_own_chromosome_name(refined):
    """A pool keyed 'chrV' against 'chr5' reads silently refines nothing."""
    assert {c for c, _s, _e in refined['pool']} == {'chr5'}


# ---------------------------------------------------------------------------
# Pool coordinates (ISSUE-004)
# ---------------------------------------------------------------------------

def test_pool_holds_no_junction_that_exists_in_no_read(refined):
    """Every observed (non-annotated) pool junction must come from a real N-op.

    A small residual is expected and allowed: `_merge_del_into_intron`
    deliberately folds a D abutting an N into the intron, so such a junction is
    stored one D-length from the read's raw N-op.  Before the cursor fix this
    number was 262 of 308 (85 %).
    """
    observed_only = {j for j in refined['pool'] if j not in refined['annot3']}
    phantom = observed_only - refined['observed']
    assert len(phantom) <= refined['baseline']['max_phantom'], \
        (refined['dataset'], len(phantom), sorted(phantom)[:10])


# ---------------------------------------------------------------------------
# FP guards
# ---------------------------------------------------------------------------

def test_no_annotated_junction_is_moved_off_annotation(refined):
    harmful = [
        (name, old, new) for name, chrom, old, new in _moves(refined)
        if (chrom,) + old in refined['annot3']
        and (chrom,) + new not in refined['annot3']
    ]
    assert len(harmful) <= refined['baseline']['max_harmful'], \
        (refined['dataset'], harmful)


def test_module_2h_creates_no_non_canonical_junction(refined):
    seq = refined['seq']
    created = [
        (name, old, new, _motif(seq, *new)) for name, _c, old, new in _moves(refined)
        if _motif(seq, *new) not in CANONICAL and _motif(seq, *old) in CANONICAL
    ]
    assert created == [], refined['dataset']


def test_read_count_is_preserved(refined):
    assert set(refined['after']) == set(refined['before'])


# ---------------------------------------------------------------------------
# FN guard — the corrections must still happen
# ---------------------------------------------------------------------------

def test_drifted_junctions_are_still_pulled_onto_annotation(refined):
    beneficial = [
        (name, old, new) for name, chrom, old, new in _moves(refined)
        if (chrom,) + old not in refined['annot3']
        and (chrom,) + new in refined['annot3']
    ]
    assert len(beneficial) >= refined['baseline']['min_beneficial'], \
        (refined['dataset'], len(beneficial), beneficial)
