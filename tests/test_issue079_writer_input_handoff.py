"""ISSUE-079: the correction TSV travels with the alignment state it was decided against.

The per-aligner stage corrects each arm AFTER Module 2H has moved its junction placements, but
the deferred consumers (lazy merge scoring, the final rectified-BAM writer) used to replay the TSV
onto the ORIGINAL arm BAM. Every 2H placement was undone in the final BAM while the merged TSV kept
describing it. These tests drive the real stage -> merge -> final writer, nothing mocked, and
compare where the read's bases land.
"""
from argparse import Namespace
from pathlib import Path

import pysam
import pytest

from rectify.core.bam import writer_input
from rectify.core.commands import correct_command
from rectify.core.commands.run.stages import _run_correction_per_aligner
from rectify.core.consensus.corrected_consensus import (
    _read_corrected_tsv_or_manifest, _stage_raw_bams, merge_corrected_tsvs,
    write_corrected_consensus_bam,
)
from tests.test_ultracode_block_shift import block_input

ARMS = ('minimap2', 'uLTRA')


def _placements(path):
    """query_name -> (start, cigar, SEQ, QUAL, every (query, reference) cell)."""
    with pysam.AlignmentFile(str(path)) as bam:
        return {
            r.query_name: (r.reference_start, r.cigarstring, r.query_sequence,
                           tuple(r.query_qualities), tuple(r.get_aligned_pairs(matches_only=True)))
            for r in bam
        }


@pytest.fixture
def two_arm_library(tmp_path):
    return build_two_arm_library(tmp_path)


def build_two_arm_library(tmp_path):
    """Both strands x terminal/internal 16-base block; one read 2H must move, one it must not."""
    raw = tmp_path / 'raw'
    raw.mkdir()
    genome, reads, truth = {}, [], {}
    for reverse in (False, True):
        for internal in (False, True):
            original, expected, g, _ = block_input(reverse, internal, origin='shifted')
            genome[original.reference_name] = g
            for role, source in (('needs_2h', original), ('stable', expected)):
                r = source.__copy__()
                r.query_name = f'handoff_{int(reverse)}_{int(internal)}_{role}'
                for tag in ('MD', 'NM', 'AS', 'cs'):  # no MD: indel correction stays off
                    r.set_tag(tag, None)
                reads.append(r)
                truth[r.query_name] = (expected.reference_start, expected.cigarstring)
    header = pysam.AlignmentHeader.from_references(list(genome), [len(s) for s in genome.values()])
    bams = {}
    for arm in ARMS:
        bams[arm] = raw / f'{arm}.bam'
        with pysam.AlignmentFile(str(bams[arm]), 'wb', header=header) as out:
            for r in sorted(reads, key=lambda r: (r.reference_name, r.reference_start, r.query_name)):
                out.write(pysam.AlignedSegment.fromstring(r.to_string(), header))
        pysam.index(str(bams[arm]))
    fasta = tmp_path / 'genome.fa'
    fasta.write_text(''.join(f'>{c}\n{s}\n' for c, s in genome.items()))
    pysam.faidx(str(fasta))
    gtf = []
    for i, r in enumerate(reads):
        if not r.query_name.endswith('stable'):
            continue
        strand = '-' if r.is_reverse else '+'
        attrs = f'gene_id "g{i}"; transcript_id "t{i}"; gene_name "G{i}";'
        gtf.append(f'{r.reference_name}\tfx\tgene\t{r.reference_start + 1}\t{r.reference_end}'
                   f'\t.\t{strand}\t.\t{attrs}')
        pos = r.reference_start
        for op, n in r.cigartuples:
            if op == 0:
                gtf.append(f'{r.reference_name}\tfx\texon\t{pos + 1}\t{pos + n}\t.\t{strand}\t.\t{attrs}')
            if op in (0, 2, 3):
                pos += n
    annotation = tmp_path / 'genes.gtf'
    annotation.write_text('\n'.join(gtf) + '\n')
    return bams, fasta, annotation, genome, truth


def _stage_args(emit):
    return Namespace(threads=1, aligner_concurrency='1', streaming=False,
                     write_corrected_bam=emit, write_softclip_bam=emit,
                     organism='saccharomyces_cerevisiae', drs=True)


def _merge_and_write(out, tsvs, corrected, replay_bams, genome):
    merged = out / 'merged.tsv'
    with _stage_raw_bams({a: str(p) for a, p in replay_bams.items()}) as staged:
        merge_corrected_tsvs(
            tsvs, merged, summary_tsv=out / 'selection.tsv',
            per_aligner_corrected_bams={a: str(p) for a, p in corrected.items()} or None,
            per_aligner_raw_bams=staged if not corrected else None,
            genome=genome, lazy_scoring_workers=1)
        write_corrected_consensus_bam(staged, tsvs, merged, out / 'final.bam', genome,
                                      threads=1, strict=True)
    return merged, out / 'final.bam'


@pytest.mark.parametrize('emit', [False, True], ids=['no_arm_bams', 'arm_bams'])
def test_final_bam_keeps_every_2h_placement(two_arm_library, tmp_path, emit):
    bams, fasta, annotation, genome, truth = two_arm_library
    out = tmp_path / 'run'
    out.mkdir()
    tsvs, corrected = _run_correction_per_aligner(bams, out, fasta, annotation, _stage_args(emit))
    assert set(tsvs) == set(ARMS) and bool(corrected) == emit

    replay = writer_input.resolve_writer_inputs(tsvs, bams)
    for arm in ARMS:
        assert Path(replay[arm]).name == 'corrected_reads.writer_input.bam'
        assert writer_input.load(tsvs[arm])['refinement']['refined'] == 4
    merged, final = _merge_and_write(out, tsvs, corrected, replay, genome)

    raw, got = _placements(bams['minimap2']), _placements(final)
    assert set(got) == set(truth)
    for name, (start, cigar) in truth.items():
        assert got[name][:2] == (start, cigar), name
        assert got[name][2:4] == raw[name][2:4], f'{name}: SEQ/QUAL changed'
    moved = [n for n in truth if raw[n][:2] != truth[n]]
    assert len(moved) == 4, 'the fixture must contain reads 2H actually moves'

    # The merged TSV and the final BAM describe the same introns.
    rows = _read_corrected_tsv_or_manifest(merged).set_index('read_id')
    with pysam.AlignmentFile(str(final)) as bam:
        for r in bam:
            pos, n_ops = r.reference_start, []
            for op, n in r.cigartuples:
                if op == 3:
                    n_ops.append(f'{pos}-{pos + n}')
                if op in (0, 2, 3):
                    pos += n
            assert sorted(str(rows.loc[r.query_name, 'junctions']).split(';')) == sorted(n_ops)

    # Writable control: the old caller contract (raw arms) reproduces the defect on these reads,
    # so the assertions above are testing the handoff and not an inert fixture.
    legacy = tmp_path / 'legacy'
    legacy.mkdir()
    _, legacy_final = _merge_and_write(legacy, tsvs, corrected, bams, genome)
    undone = _placements(legacy_final)
    assert all(undone[n][:2] == raw[n][:2] for n in moved)
    assert all(len(set(undone[n][4]) - set(got[n][4])) == 16 for n in moved)


def test_resume_reuses_only_a_vouched_pair(two_arm_library, tmp_path, monkeypatch):
    bams, fasta, annotation, _, _ = two_arm_library
    out = tmp_path / 'run'
    out.mkdir()
    calls = []
    real_run = correct_command.run
    monkeypatch.setattr(correct_command, 'run',
                        lambda a: (calls.append(Path(str(a.output)).parent.name), real_run(a))[1])

    def stage():
        calls.clear()
        tsvs, _ = _run_correction_per_aligner(bams, out, fasta, annotation, _stage_args(False))
        return tsvs, sorted(calls)

    tsvs, ran = stage()
    assert ran == sorted(ARMS)
    assert stage()[1] == [], 'an intact pair is reused'

    arm_dir = tsvs['uLTRA'].parent
    # legacy output: a TSV with no receipt must be rebuilt, never replayed onto raw geometry
    writer_input.receipt_path(tsvs['uLTRA']).unlink()
    assert stage()[1] == ['uLTRA']
    # the retained BAM is gone (e.g. excluded from a scratch sync)
    (arm_dir / 'corrected_reads.writer_input.bam').unlink()
    assert stage()[1] == ['uLTRA']
    # same name, different bytes
    retained = arm_dir / 'corrected_reads.writer_input.bam'
    retained.write_bytes(retained.read_bytes() + b'\0')
    assert stage()[1] == ['uLTRA']
    # the TSV changed under the receipt
    region = arm_dir / 'corrected_reads.region_000.tsv'
    region.write_text(region.read_text() + '\n')
    assert stage()[1] == ['uLTRA']
    # an evicted pair cannot feed the writer again
    assert writer_input.evict(tsvs['uLTRA']) and not retained.exists()
    with pytest.raises(writer_input.WriterInputError, match='evicted'):
        writer_input.resolve_writer_input(tsvs['uLTRA'], bams['uLTRA'])
    assert stage()[1] == ['uLTRA']


def test_failed_correction_cannot_pass_off_the_previous_tsv(two_arm_library, tmp_path, monkeypatch):
    bams, fasta, annotation, _, _ = two_arm_library
    out = tmp_path / 'run'
    out.mkdir()
    tsvs, _ = _run_correction_per_aligner(bams, out, fasta, annotation, _stage_args(False))
    # The candidate pool changes (one arm realigned), so every arm must be corrected again ...
    with pysam.AlignmentFile(str(bams['uLTRA'])) as src:
        records, header = list(src)[:-1], src.header
    with pysam.AlignmentFile(str(bams['uLTRA']), 'wb', header=header) as dst:
        for r in records:
            dst.write(r)
    pysam.index(str(bams['uLTRA']))

    # ... and this time correction dies. The old TSVs are still on disk and still parse.
    def boom(_args):
        raise RuntimeError('node failure')
    monkeypatch.setattr(correct_command, 'run', boom)
    again, _ = _run_correction_per_aligner(bams, out, fasta, annotation, _stage_args(False))
    assert again == {}
    assert all(p.exists() for p in tsvs.values())


def test_source_receipt_when_2h_does_not_run(two_arm_library, tmp_path):
    """No annotation -> no 2H -> the arm BAM itself is the writer input; nothing is copied."""
    bams, fasta, _, _, _ = two_arm_library
    out = tmp_path / 'run'
    out.mkdir()
    tsvs, _ = _run_correction_per_aligner(bams, out, fasta, None, _stage_args(False))
    replay = writer_input.resolve_writer_inputs(tsvs, bams)
    assert replay == {a: str(p) for a, p in bams.items()}
    assert not list(out.rglob('*.writer_input.bam'))
    assert not list(out.rglob('*junction_refined_*'))


def test_staging_keeps_same_named_arms_apart(tmp_path):
    paths = {}
    for arm in ARMS:
        (tmp_path / arm).mkdir()
        paths[arm] = tmp_path / arm / 'corrected_reads.writer_input.bam'
        paths[arm].write_text(arm)
    with _stage_raw_bams({a: str(p) for a, p in paths.items()}, scratch_root=str(tmp_path)) as staged:
        assert len(set(staged.values())) == len(ARMS)
        assert {a: Path(p).read_text() for a, p in staged.items()} == {a: a for a in ARMS}


def test_an_arm_that_needed_sorting_keeps_its_corrected_bam(two_arm_library, tmp_path):
    """ISSUE-084: an arm with no .bai is coordinate-sorted first, and its corrected BAM is then
    named after the SORTED copy. The collector looked for the unsorted name, so the arm came back
    without a corrected BAM — and an arm without one scores inf in the merge and can never win."""
    bams, fasta, annotation, _, _ = two_arm_library
    Path(str(bams['uLTRA']) + '.bai').unlink()
    out = tmp_path / 'run'
    out.mkdir()
    tsvs, corrected = _run_correction_per_aligner(bams, out, fasta, annotation, _stage_args(True))
    assert set(tsvs) == set(corrected) == set(ARMS)
    # ... and again on the resume path, which collects without having sorted anything.
    _, corrected_again = _run_correction_per_aligner(bams, out, fasta, annotation, _stage_args(True))
    assert corrected_again == corrected
