"""Chromosome-name standardization must not romanize a non-yeast genome (ISSUE-001).

``standardize_chrom_name`` carries a yeast arabic->roman fallback: ``chr1`` ->
``chrI``.  Applied to a human genome it rewrites ``chr5`` to ``chrV`` and —
worse — ``chr10`` to ``chrX``, which does not merely mislabel but MERGES those
reads with the real human chrX.  The guard is the loaded-genome contig registry:
a name the genome already has is returned verbatim, and the roman fallback fires
only when the ROMANIZED name is itself a contig of that genome.

The registry is process-global module state, so every test here restores it.

Author: Kevin R. Roy (agent S1)
"""

import sys
from pathlib import Path

import pytest

RECTIFY_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(RECTIFY_ROOT))

from rectify.utils import genome as gu  # noqa: E402
from rectify.utils.genome import (  # noqa: E402
    known_contig_count,
    register_genome_contigs,
    standardize_chrom_name,
)

HUMAN = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
YEAST = ["chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII",
         "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI",
         "chrMito"]


@pytest.fixture
def registry():
    """Give each test a clean registry and put the old one back."""
    saved = set(gu._KNOWN_CONTIGS)
    gu._KNOWN_CONTIGS.clear()
    yield gu._KNOWN_CONTIGS
    gu._KNOWN_CONTIGS.clear()
    gu._KNOWN_CONTIGS.update(saved)


# ---------------------------------------------------------------------------
# The trap, stated
# ---------------------------------------------------------------------------

def test_empty_registry_is_the_legacy_yeast_mode(registry):
    """Documented, not endorsed: with nothing registered the fallback fires."""
    assert known_contig_count() == 0
    assert standardize_chrom_name("chr5") == "chrV"
    assert standardize_chrom_name("chr10") == "chrX"     # collides with real chrX


# ---------------------------------------------------------------------------
# Human
# ---------------------------------------------------------------------------

def test_human_contigs_round_trip_verbatim(registry):
    register_genome_contigs(HUMAN)
    for name in HUMAN:
        assert standardize_chrom_name(name) == name


def test_human_chr10_does_not_become_chrX(registry):
    """The dangerous one: chrX exists in human, so this silently MERGES reads."""
    register_genome_contigs(HUMAN)
    assert standardize_chrom_name("chr10") == "chr10"
    assert standardize_chrom_name("chrX") == "chrX"


def test_chr17_to_chr22_are_untouched_either_way(registry):
    """The roman map only covers 1..16, so these were never rewritten —
    which is exactly why a partial fix would look like it worked."""
    for name in [f"chr{i}" for i in range(17, 23)]:
        assert standardize_chrom_name(name) == name          # empty registry
    register_genome_contigs(HUMAN)
    for name in [f"chr{i}" for i in range(17, 23)]:
        assert standardize_chrom_name(name) == name


def test_a_name_absent_from_a_human_genome_is_still_not_romanized(registry):
    """An annotation contig the FASTA lacks must not fall back to the roman map."""
    register_genome_contigs(HUMAN)
    assert standardize_chrom_name("chr5_GL000000v2_random") == \
        "chr5_GL000000v2_random"


# ---------------------------------------------------------------------------
# Yeast — the behavior that must survive
# ---------------------------------------------------------------------------

def test_yeast_still_romanizes(registry):
    register_genome_contigs(YEAST)
    assert standardize_chrom_name("chr1") == "chrI"
    assert standardize_chrom_name("chr10") == "chrX"
    assert standardize_chrom_name("chr16") == "chrXVI"
    assert standardize_chrom_name("chrM") == "chrMito"


def test_yeast_roman_names_pass_through(registry):
    register_genome_contigs(YEAST)
    for name in YEAST:
        assert standardize_chrom_name(name) == name


def test_bare_roman_numerals_still_gain_the_chr_prefix(registry):
    register_genome_contigs(YEAST)
    assert standardize_chrom_name("IV") == "chrIV"


# ---------------------------------------------------------------------------
# The count helper commands log
# ---------------------------------------------------------------------------

def test_known_contig_count(registry):
    assert known_contig_count() == 0
    register_genome_contigs(HUMAN)
    assert known_contig_count() == len(HUMAN)
    register_genome_contigs(HUMAN)          # idempotent union
    assert known_contig_count() == len(HUMAN)


def test_correct_registers_contigs_before_module_2h():
    """The ordering IS the fix: registration must precede the 2H block.

    Module 2H loads the annotation, builds the pool and keys the chrom index; if
    contigs are registered after that, those three are keyed 'chrV' while the
    reads (standardized once load_genome has run) are keyed 'chr5', and every
    candidate lookup misses.
    """
    src = (RECTIFY_ROOT / 'rectify/core/commands/correct_command.py').read_text()
    reg = src.index('_register_contigs(str(config[')
    mod2h = src.index('Module 2H: Junction N-op boundary refinement')
    assert reg < mod2h, "contig registration must come before the Module 2H block"
