"""Human RNA004 DRS regression guard — Module 2H on the Sumner chr5 panel.

The yeast validation set cannot see the defects human data exposed: ~95 % of
intron-bearing yeast genes have ONE intron, so a walker that never advanced past
an N op looked correct, and yeast chromosome names are exactly the ones the
arabic->roman fallback was written for.

This module runs the Module 2H stage the way `rectify correct` runs it — contigs
registered first, annotation, pool, refine — over the 100-read panel and asserts
the FP/FN properties that were violated before 2026-09-05:

  FP  no annotated canonical junction is moved off annotation;
      no non-canonical N-op is created;
      the observed junction pool holds no coordinate that exists in no read.
  FN  the beneficial corrections (novel/drifted -> annotated) still happen.

The panel and its chr5 reference are git-ignored (a 184 MB FASTA and a 47 MB
GTF), so the whole module skips when they are absent.  Point
RECTIFY_SUMNER_PANEL_DIR at the directory to run it elsewhere.

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
REF_FA = PANEL_DIR / 'ref' / 'chr5.fa'
REF_GTF = PANEL_DIR / 'ref' / 'gencode.v48.basic.chr5.gtf'

pytestmark = pytest.mark.skipif(
    not (PANEL_BAM.exists() and REF_FA.exists() and REF_GTF.exists()),
    reason=f"Sumner human panel not present under {PANEL_DIR} (git-ignored data)",
)

# Human canonical splice dinucleotides, read on the genomic PLUS strand: the
# minus-strand form of each is listed alongside (GT-AG / CT-AC, GC-AG / CT-GC,
# AT-AC / GT-AT — the last pair is the U12 class).
CANONICAL = {('GT', 'AG'), ('CT', 'AC'),
             ('GC', 'AG'), ('CT', 'GC'),
             ('AT', 'AC'), ('GT', 'AT')}

_N = 3
_SAM_REF_CONSUMING = {0, 2, 3, 7, 8}

# Baselines recorded on 2026-09-05 with 7d571f5 + 250a608 + 21dbc57 + 26263d6.
MAX_HARMFUL = 0
MIN_BENEFICIAL = 5


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


@pytest.fixture(scope="module")
def refined(tmp_path_factory):
    """Run Module 2H over the panel; return everything the assertions need."""
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
    pool, annot_set = build_junction_pool([str(PANEL_BAM)], annot)
    genome = load_genome(str(REF_FA))
    # The same tables `rectify correct --organism homo_sapiens` now selects
    # (ISSUE-005) — the empirical HP costs are what put a 1 nt acceptor slide
    # inside the noise floor where the annotated-canonical gate can hold it.
    jpt, spt = select_penalty_tables({'organism': 'homo_sapiens'})
    assert jpt is not None, "bundled human DRS penalty table should be found"

    out_bam = str(tmp_path_factory.mktemp("panel") / "refined.bam")
    refine_bam_junctions(
        input_bam=str(PANEL_BAM), output_bam=out_bam,
        aligner_bams=[str(PANEL_BAM)], annotated_junctions=annot, genome=genome,
        prebuilt_junction_pool=pool, prebuilt_annotated_set=annot_set,
        penalty_table_path=jpt, str_penalty_table_path=spt,
        sort_and_index=True, n_workers=1,
    )

    before, after, observed = {}, {}, set()
    with pysam.AlignmentFile(str(PANEL_BAM)) as bam:
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
    assert len(phantom) <= 5, sorted(phantom)[:10]


# ---------------------------------------------------------------------------
# FP guards
# ---------------------------------------------------------------------------

def test_no_annotated_junction_is_moved_off_annotation(refined):
    harmful = [
        (name, old, new) for name, chrom, old, new in _moves(refined)
        if (chrom,) + old in refined['annot3']
        and (chrom,) + new not in refined['annot3']
    ]
    assert len(harmful) <= MAX_HARMFUL, harmful


def test_module_2h_creates_no_non_canonical_junction(refined):
    seq = refined['seq']
    created = [
        (name, old, new, _motif(seq, *new)) for name, _c, old, new in _moves(refined)
        if _motif(seq, *new) not in CANONICAL and _motif(seq, *old) in CANONICAL
    ]
    assert created == []


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
    assert len(beneficial) >= MIN_BENEFICIAL, beneficial
