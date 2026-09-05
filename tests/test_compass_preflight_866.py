"""planning 873 / 861 §S3 / 862: the COMPASS + deSALT index pre-flight.

Three failures this replaces, all of which exit 0 or hang rather than saying what is wrong:

* ``--Scer`` bundles a genome and a GTF but **no index of any kind**, and nothing in rectify builds
  one. On Hoffman2 that made ``run-all --short-read --Scer`` deadlock for six hours (hisat2's Perl
  wrapper orphans itself holding the stdout pipe).
* A dropped arm still exits 0, so the panel silently shrinks with no record of which arms ran.
* 🔴 The bundled GTF has **no ``exon`` features** (``mRNA`` 11,125 · ``CDS`` 7,072 · ``gene`` 6,613 ·
  ``intron`` 378). STAR's ``--sjdbGTFfile`` and ``hisat2_extract_splice_sites.py`` both read
  ``exon`` by default, so an "annotated" index built from it builds cleanly, aligns at a normal
  rate, and contains zero annotated junctions.

Pure functions over text and paths: no aligner, no index, no genome larger than a few lines.
"""
import gzip
from pathlib import Path

import pytest

from rectify.core.align.compass_preflight import (
    compass_preflight,
    gtf_exon_status,
    gtf_feature_census,
    synthesize_exon_gtf,
)

BUNDLED = Path(__file__).resolve().parents[1] / 'rectify/data/genomes/saccharomyces_cerevisiae'
GTF = BUNDLED / 'saccharomyces_cerevisiae_R64-5-1_20240529.gtf'
FASTA = BUNDLED / 'S288C_reference_sequence_R64-5-1_20240529.fsa.gz'


def _gtf(tmp_path, rows, name='a.gtf'):
    p = tmp_path / name
    p.write_text(''.join(
        f'chrI\tSGD\t{feat}\t{lo}\t{hi}\t.\t+\t.\tgene_id "{g}"; transcript_id "{g}_mRNA";\n'
        for feat, lo, hi, g in rows))
    return p


# ── the census / status classifier ────────────────────────────────────────────
def test_census_counts_features_and_skips_comments(tmp_path):
    p = tmp_path / 'x.gtf'
    p.write_text('#!genome-build R64\n'
                 'chrI\tSGD\texon\t1\t100\t.\t+\t.\tgene_id "A";\n'
                 'chrI\tSGD\texon\t200\t300\t.\t+\t.\tgene_id "A";\n'
                 'chrI\tSGD\tCDS\t1\t100\t.\t+\t.\tgene_id "A";\n'
                 'short line\n')
    c = gtf_feature_census(p)
    assert c['exon'] == 2 and c['CDS'] == 1 and 'short line' not in c


def test_census_reads_gzip(tmp_path):
    p = tmp_path / 'x.gtf.gz'
    with gzip.open(p, 'wt') as fh:
        fh.write('chrI\tSGD\texon\t1\t100\t.\t+\t.\tgene_id "A";\n')
    assert gtf_feature_census(p)['exon'] == 1


def test_census_of_an_unreadable_file_is_empty_not_an_exception():
    """A pre-flight must never be the thing that stops a run."""
    assert gtf_feature_census('/no/such/file.gtf') == {}


@pytest.mark.parametrize('rows,expected', [
    ([('exon', 1, 100, 'A'), ('CDS', 1, 100, 'A')], 'exon'),
    ([('CDS', 1, 100, 'A'), ('mRNA', 1, 100, 'A')], 'cds_only'),
    ([('gene', 1, 100, 'A')], 'unusable'),
])
def test_exon_status_classifies(tmp_path, rows, expected):
    assert gtf_exon_status(_gtf(tmp_path, rows))[0] == expected


def test_exon_status_missing_for_none_and_absent_file():
    assert gtf_exon_status(None)[0] == 'missing'
    assert gtf_exon_status('/no/such/file.gtf')[0] == 'missing'


@pytest.mark.skipif(not GTF.exists(), reason='bundled GTF not present')
def test_the_bundled_scer_gtf_is_the_exon_less_case():
    """This is the whole point: --Scer's own annotation cannot build an annotated index."""
    status, census = gtf_exon_status(GTF)
    assert status == 'cds_only'
    assert census.get('exon', 0) == 0
    assert census['CDS'] > 5000 and census['mRNA'] > 10000


# ── the repair ────────────────────────────────────────────────────────────────
def test_synthesize_adds_one_exon_per_cds_keyed_by_transcript(tmp_path):
    src = _gtf(tmp_path, [('CDS', 1, 100, 'A'), ('CDS', 200, 300, 'A'), ('CDS', 1, 50, 'B')])
    dst = tmp_path / 'out.gtf'
    stats = synthesize_exon_gtf(src, dst)
    assert stats['exon_written'] == 3 and stats['transcripts'] == 2
    out = gtf_feature_census(dst)
    assert out['exon'] == 3 and out['CDS'] == 3          # CDS lines are kept, not replaced
    # coordinates and strand are copied verbatim
    exons = [l.split('\t') for l in dst.read_text().splitlines() if l.split('\t')[2] == 'exon']
    assert sorted((int(e[3]), int(e[4])) for e in exons) == [(1, 50), (1, 100), (200, 300)]
    assert all(e[6] == '+' for e in exons)


def test_synthesize_is_a_no_op_when_exons_already_exist(tmp_path):
    src = _gtf(tmp_path, [('exon', 1, 100, 'A'), ('CDS', 1, 100, 'A')])
    dst = tmp_path / 'out.gtf'
    stats = synthesize_exon_gtf(src, dst)
    assert stats['exon_written'] == 0 and stats['exon_already'] == 1
    assert gtf_feature_census(dst)['exon'] == 1          # not doubled
    assert dst.read_text() == src.read_text()            # byte copy


def test_synthesize_handles_gff3_parent_attributes(tmp_path):
    src = tmp_path / 'a.gff'
    src.write_text('chrI\tSGD\tCDS\t1\t100\t.\t-\t.\tID=cds1;Parent=tx1,tx2\n'
                   'chrI\tSGD\tCDS\t200\t300\t.\t-\t.\tID=cds2;Parent=tx1\n')
    stats = synthesize_exon_gtf(src, tmp_path / 'o.gff')
    assert stats['exon_written'] == 2 and stats['transcripts'] == 1   # first Parent only


@pytest.mark.skipif(not GTF.exists(), reason='bundled GTF not present')
def test_synthesizing_from_the_bundled_gtf_makes_it_usable(tmp_path):
    dst = tmp_path / 'scer_with_exons.gtf'
    stats = synthesize_exon_gtf(GTF, dst)
    assert stats['exon_written'] == stats['cds_in'] > 5000
    assert gtf_exon_status(dst)[0] == 'exon'


# ── the report ────────────────────────────────────────────────────────────────
@pytest.mark.skipif(not FASTA.exists(), reason='bundled genome not present')
def test_scer_reports_every_arm_missing_and_the_exon_warning():
    rep = compass_preflight(FASTA, GTF, desalt=True)
    assert not rep.ok
    assert {a.arm for a in rep.missing} >= {'STAR', 'HISAT2', 'deSALT'}
    text = rep.render()
    assert 'SILENTLY UNANNOTATED' in text
    assert 'would be DROPPED' in text
    # every missing arm prints a runnable build command
    for a in rep.missing:
        assert a.build_cmd and str(a.path) not in ('', '.')


def test_present_index_is_reported_ok(tmp_path):
    """A built index must not be reported missing — otherwise the check is noise."""
    gdir = tmp_path / 'ref'
    gdir.mkdir()
    genome = gdir / 'g.fa'
    genome.write_text('>chrI\nACGT\n')
    star = gdir / 'STAR_annotated_150_bp_SJDB_index'
    star.mkdir()
    (star / 'SAindex').write_text('x')
    rep = compass_preflight(genome, None, arms=['STAR'])
    assert [a.present for a in rep.arms] == [True]
    assert rep.missing == []


def test_arms_argument_limits_what_is_checked(tmp_path):
    genome = tmp_path / 'g.fa'
    genome.write_text('>chrI\nACGT\n')
    rep = compass_preflight(genome, None, arms=['STAR'])
    assert [a.arm for a in rep.arms] == ['STAR']
    rep2 = compass_preflight(genome, None, arms=[], desalt=True)
    assert [a.arm for a in rep2.arms] == ['deSALT']


def test_desalt_placeholder_dir_is_not_counted_as_built(tmp_path):
    """An empty desalt_index/ reads as present to a naive directory check."""
    genome = tmp_path / 'g.fa'
    genome.write_text('>chrI\nACGT\n')
    (tmp_path / 'desalt_index').mkdir()
    rep = compass_preflight(genome, None, arms=[], desalt=True)
    assert rep.arms[0].present is False
    (tmp_path / 'desalt_index' / 'ref.seq').write_text('x')
    assert compass_preflight(genome, None, arms=[], desalt=True).arms[0].present is True


def test_missing_annotation_alone_does_not_make_the_report_not_ok(tmp_path):
    """No annotation is a legitimate configuration; an exon-less one is not."""
    gdir = tmp_path / 'ref'
    gdir.mkdir()
    genome = gdir / 'g.fa'
    genome.write_text('>chrI\nACGT\n')
    star = gdir / 'STAR_annotated_150_bp_SJDB_index'
    star.mkdir()
    (star / 'SAindex').write_text('x')
    assert compass_preflight(genome, None, arms=['STAR']).ok is True
    bad = _gtf(tmp_path, [('CDS', 1, 100, 'A')])
    assert compass_preflight(genome, bad, arms=['STAR']).ok is False


# ── the wiring ────────────────────────────────────────────────────────────────
def test_run_all_exposes_require_compass_index():
    import argparse
    from rectify.core.commands.run_command import create_run_parser
    p = argparse.ArgumentParser()
    create_run_parser(p.add_subparsers(dest='command'))
    # run-all requires -o/--output-dir; omitting it makes argparse SystemExit(2).
    base = ['run-all', 'r.fq', '-o', 'out', '--short-read']
    assert p.parse_args(base).require_compass_index is False
    assert p.parse_args(base + ['--require-compass-index']).require_compass_index is True


def test_alignment_stage_runs_the_preflight_before_launching_the_panel():
    import inspect
    from rectify.core.commands.run import stages
    assert 'require_compass_index' in inspect.signature(stages._run_alignment).parameters
    src = inspect.getsource(stages._run_alignment)
    assert src.index('compass_preflight(') < src.index('-aligner consensus')
    from rectify.core.commands.run import single_sample
    assert "require_compass_index=getattr(args, 'require_compass_index', False)" in \
        inspect.getsource(single_sample)
