"""Sequence-placement regressions for the conservative terminal-tail post-pass."""
import json

import pysam
import pytest

from rectify.core.splice.terminal_tail_placement import (
    conditional_exact_bound, refine_terminal_tail, TAG,
)


def fixture(strand='+', previous_intron=False, tail=15, prefix='GTATATTTTC'):
    genome = list('C' * 600)
    body = [(0, 20), (3, 30), (0, 20)] if previous_intron else [(0, 40)]
    ns = 100 + sum(n for op, n in body if op in (0, 3))
    ne = ns + 100
    genome[ns:ns + len(prefix)] = prefix
    genome[ne:ne + len(prefix)] = 'C' * len(prefix)
    seq = 'C' * 40 + prefix + 'A' * tail
    ops = body + [(3, 100), (0, len(prefix)), (4, tail)]
    start = 100
    if strand == '-':
        tr = str.maketrans('ACGT', 'TGCA')
        genome = list(''.join(genome).translate(tr)[::-1])
        seq = seq.translate(tr)[::-1]
        start = 600 - (ne + len(prefix))
        ops = ops[::-1]
        ns, ne = 600 - ne, 600 - ns
    genome = ''.join(genome)
    header = pysam.AlignmentHeader.from_references(['chrT'], [600])
    read = pysam.AlignedSegment(header)
    read.query_name = 'terminal_test'
    read.reference_id = 0
    read.reference_start = start
    read.flag = 16 if strand == '-' else 0
    read.mapping_quality = 60
    read.query_sequence = seq
    read.query_qualities = pysam.qualitystring_to_array('I' * len(seq))
    read.cigartuples = ops
    read.set_tag('NM', 9)
    read.set_tag('MD', '40C')
    read.set_tag('XA', 0)
    return read, genome, (ns, ne)


def apply(read, genome, strand, annotation=()):
    return refine_terminal_tail(read, lambda c, s, e: genome[s:e],
        rna_strand=strand, protocol='ont-cdna', annotated_junctions=set(annotation))


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('previous_intron', [False, True])
@pytest.mark.parametrize('equals', [False, True])
def test_exact_native_preserves_every_body_base(strand, previous_intron, equals):
    read, genome, junction = fixture(strand, previous_intron)
    literal = read.query_sequence
    qualities = list(read.query_qualities)
    before = dict((q, r) for q, r in read.get_aligned_pairs() if q is not None)
    if equals:
        encoded = list(literal)
        for q, r in before.items():
            if r is not None and literal[q] == genome[r]:
                encoded[q] = '='
        read.query_sequence = ''.join(encoded)
        read.query_qualities = qualities
    result = apply(read, genome, strand)
    assert result.applied, result
    assert read.query_sequence == literal
    assert list(read.query_qualities) == qualities
    assert not read.has_tag('NM') and not read.has_tag('MD')
    assert read.get_tag('XA') == 0
    assert sum(op == 3 for op, _ in read.cigartuples) == int(previous_intron)
    after = dict((q, r) for q, r in read.get_aligned_pairs() if q is not None)
    # Body is the first/last forty query bases. Every upstream N is retained.
    body = range(40) if strand == '+' else range(len(literal) - 40, len(literal))
    assert all(before[q] == after[q] for q in body)
    evidence = json.loads(read.get_tag(TAG))
    assert evidence['native_matches'] == 10
    assert evidence['native_matches'] > evidence['current_matches']
    assert evidence['tail_splits_searched'] == 1


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('prefix', ['A', 'G' * 9, 'GTATATTTTC'])
def test_annotated_short_terminal_exon_protected_before_tail_gate(strand, prefix):
    read, genome, (ns, ne) = fixture(strand, prefix=prefix, tail=30)
    before = read.to_string()
    result = apply(read, genome, strand, [('chrT', ns, ne, strand)])
    assert result.reason == 'annotated_junction'
    assert read.to_string() == before


@pytest.mark.parametrize('strand', ['+', '-'])
def test_true_novel_terminal_exon_and_native_tie_are_unchanged(strand):
    read, genome, (ns, ne) = fixture(strand)
    g = list(genome)
    qref = {q: r for q, r in read.get_aligned_pairs() if q is not None and r is not None}
    query_range = range(40, 50) if strand == '+' else range(15, 25)
    for q in query_range:
        g[qref[q]] = read.query_sequence[q]
    genome = ''.join(g)
    before = read.to_string()
    result = apply(read, genome, strand)
    assert result.reason == 'no_positional_gain'
    assert read.to_string() == before


@pytest.mark.parametrize('strand', ['+', '-'])
def test_insufficient_tail_and_nonexact_native_are_unchanged(strand):
    for tail in (10, 14):
        read, genome, _ = fixture(strand, tail=tail)
        before = read.to_string()
        assert apply(read, genome, strand).reason == 'insufficient_terminal_tail'
        assert read.to_string() == before
    read, genome, (ns, ne) = fixture(strand)
    p = ns if strand == '+' else ne - 1
    genome = genome[:p] + ('C' if strand == '+' else 'G') + genome[p + 1:]
    before = read.to_string()
    assert apply(read, genome, strand).reason == 'native_not_exact'
    assert read.to_string() == before


def test_missing_context_chimera_tag_and_low_information_refuse():
    read, genome, _ = fixture()
    for protocol, strand, annotations in [('quantseq-rev', '+', set()),
                                         ('ont-cdna', None, set()),
                                         ('drs', '+', None)]:
        before = read.to_string()
        result = refine_terminal_tail(read, lambda c,s,e: genome[s:e],
            rna_strand=strand, protocol=protocol, annotated_junctions=annotations)
        assert not result.applied and read.to_string() == before
    for tag, value in [('SA', 'chrT,1,+,10M,60,0;'), (TAG, 'foreign')]:
        read, genome, _ = fixture()
        read.set_tag(tag, value)
        before = read.to_string()
        assert not apply(read, genome, '+').applied
        assert read.to_string() == before
    read, genome, _ = fixture(prefix='T' * 12)
    before = read.to_string()
    assert apply(read, genome, '+').reason == 'conditional_null_not_supported'
    assert read.to_string() == before
    assert conditional_exact_bound('GTATATTTTC') == pytest.approx(2 / (2520 * 0.8))


def test_conditional_partition_changes_boundary_decision():
    # Unconditioned 2/C(10,5)=.00794 would pass .01. Fixing the non-A final
    # base halves the allowed arrangements and correctly makes this fail.
    assert conditional_exact_bound('AAAAATTTTT') > .01


def test_modern_metadata_contract_and_streaming_postpass(tmp_path):
    from rectify.core.splice.terminal_tail_placement import (
        modern_cdna_rna_strand, run_terminal_tail_postpass, terminal_tail_output_for,
        terminal_tail_context,
    )
    from rectify.core.commands.run.helpers import _collect_per_aligner_bams
    read, genome, _ = fixture('-')
    assert modern_cdna_rna_strand(read) is None
    read.set_tag('XN', 1)
    assert modern_cdna_rna_strand(read) is None
    read.set_tag('XT', 2)
    read.set_tag('XR', 'source-read-id')
    assert modern_cdna_rna_strand(read) == '-'
    read.set_tag('XN', '1')
    assert modern_cdna_rna_strand(read) is None
    read.set_tag('XN', 1)
    fasta = tmp_path / 'genome.fa'
    fasta.write_text('>chrT\n' + genome + '\n')
    pysam.faidx(str(fasta))
    source = tmp_path / 'sample.minimap2.bam'
    with pysam.AlignmentFile(str(source), 'wb', header=read.header) as sink:
        sink.write(read)
    output = tmp_path / 'sample.minimap2.terminal_tail.bam'
    selected, counts = run_terminal_tail_postpass(source, fasta, output, set())
    assert selected == str(output) and counts['applied'] == 1
    with pysam.AlignmentFile(str(output), 'rb') as bam:
        emitted = next(bam)
        assert emitted.has_tag(TAG)
        assert not any(op == 3 for op, _ in emitted.cigartuples)
    assert terminal_tail_output_for(source, output) == output
    assert terminal_tail_output_for(source, output, context={"changed": True}) is None
    assert _collect_per_aligner_bams('sample', tmp_path)['minimap2'] == source
    assert _collect_per_aligner_bams('sample', tmp_path, terminal_tail_context=terminal_tail_context())['minimap2'] == output
    # Replacing the input with a legacy record invalidates disk discovery and
    # incurs no rewritten output during a fresh post-pass.
    read.set_tag('XN', None)
    with pysam.AlignmentFile(str(source), 'wb', header=read.header) as sink:
        sink.write(read)
    assert terminal_tail_output_for(source, output) is None
    def unnecessary_annotation_load():
        raise AssertionError('legacy context must not load full annotation')
    selected, counts = run_terminal_tail_postpass(source, fasta, output, unnecessary_annotation_load)
    assert selected == str(source) and not counts.get('applied')
    assert _collect_per_aligner_bams('sample', tmp_path, terminal_tail_context=terminal_tail_context())['minimap2'] == source


@pytest.mark.parametrize('strand', ['+', '-'])
def test_hardclip_and_contig_edge_stand_down(strand):
    read, genome, _ = fixture(strand)
    if strand == '+':
        read.cigartuples = read.cigartuples + [(5, 3)]
    else:
        read.cigartuples = [(5, 3)] + read.cigartuples
    before = read.to_string()
    assert apply(read, genome, strand).reason == 'unanchored_or_incomplete_terminal_block'
    assert read.to_string() == before
    read, genome, (ns, ne) = fixture(strand)
    before = read.to_string()
    def short_reference(chrom, start, end):
        return genome[start:end - 1]
    result = refine_terminal_tail(read, short_reference, rna_strand=strand,
        protocol='ont-cdna', annotated_junctions=set())
    assert result.reason == 'contig_edge'
    assert read.to_string() == before


def test_large_terminal_intron_is_not_materialized():
    read, genome, _ = fixture()
    # Native fetch is fixed at the body edge; the old terminal block can be
    # megabases away without expanding every N base into aligned pairs.
    read.cigartuples = [(0, 40), (3, 2_000_000), (0, 10), (4, 15)]
    from rectify.core.splice.terminal_tail_placement import _query_placements
    assert len(list(_query_placements(read))) == len(read.query_sequence)
    def fetch(chrom, start, end):
        return genome[start:end] if end <= len(genome) else 'C' * (end - start)
    result = refine_terminal_tail(read, fetch, rna_strand='+', protocol='drs',
                                  annotated_junctions=set())
    assert result.applied and read.cigarstring == '50M15S'


def test_consensus_receipt_rejects_changed_context_and_underlying_input(tmp_path):
    from shutil import copyfile
    from rectify.core.splice.terminal_tail_placement import (
        run_terminal_tail_postpass, terminal_tail_context,
        write_terminal_selection_receipt, terminal_selection_matches,
    )
    read, genome, _ = fixture()
    for tag, value in [('XN', 1), ('XT', 2), ('XR', 'source')]:
        read.set_tag(tag, value)
    fasta = tmp_path / 'g.fa'
    fasta.write_text('>chrT\n' + genome + '\n')
    pysam.faidx(str(fasta))
    annotation = tmp_path / 'g.gtf'
    annotation.write_text('')
    raw = tmp_path / 'raw.bam'
    with pysam.AlignmentFile(str(raw), 'wb', header=read.header) as sink:
        sink.write(read)
    terminal = tmp_path / 'terminal.bam'
    context = terminal_tail_context(fasta, annotation)
    selected, _ = run_terminal_tail_postpass(raw, fasta, terminal, set(), context)
    consensus = tmp_path / 'multialigned.bam'
    copyfile(selected, consensus)
    assert not terminal_selection_matches(consensus, context)
    write_terminal_selection_receipt(consensus, {'minimap2': selected}, context)
    assert terminal_selection_matches(consensus, context)
    annotation.write_text('# new annotated junction context\n')
    assert not terminal_selection_matches(consensus, terminal_tail_context(fasta, annotation))
    assert not terminal_selection_matches(consensus, None)  # protocol disabled
    # Same selected terminal BAM cannot hide a rewritten underlying raw arm.
    with pysam.AlignmentFile(str(raw), 'wb', header=read.header) as sink:
        read.query_name = 'replacement'
        sink.write(read)
    assert not terminal_selection_matches(consensus, context)


def test_ordinary_gzip_postpass_uses_shared_indexed_reference(tmp_path):
    import gzip
    from rectify.core.align.reference import open_alignment_reference
    from rectify.core.splice.terminal_tail_placement import run_terminal_tail_postpass
    read, genome, _ = fixture()
    for tag, value in [('XN', 1), ('XT', 2), ('XR', 'source')]:
        read.set_tag(tag, value)
    gz = tmp_path / 'reference.fa.gz'
    with gzip.open(str(gz), 'wt') as f:
        f.write('>chrT\n' + genome + '\n')
    original = gz.read_bytes()
    source = tmp_path / 'raw.bam'
    with pysam.AlignmentFile(str(source), 'wb', header=read.header) as f:
        f.write(read)
    chosen, counts = run_terminal_tail_postpass(source, gz, tmp_path / 'fixed.bam', set())
    assert counts['applied'] == 1
    assert gz.read_bytes() == original
    with open_alignment_reference(gz, tmp_path) as reference:
        assert reference.fetch('chrT') == genome
    assert len(list(tmp_path.glob('*.reference.fa'))) == 1


def _modern_arm(tmp_path):
    """A one-record arm BAM the post-pass edits (modern RNA-sense cDNA), its genome and an annotation."""
    read, genome, _ = fixture()
    for tag, value in [('XN', 1), ('XT', 2), ('XR', 'source')]:
        read.set_tag(tag, value)
    fasta = tmp_path / 'g.fa'
    fasta.write_text('>chrT\n' + genome + '\n')
    pysam.faidx(str(fasta))
    annotation = tmp_path / 'g.gtf'
    annotation.write_text('')
    arm = tmp_path / 'sample.minimap2.bam'
    with pysam.AlignmentFile(str(arm), 'wb', header=read.header) as sink:
        sink.write(read)
    return arm, fasta, annotation


@pytest.mark.parametrize('value', [None, '', '0', 'off'])
def test_align_skips_the_postpass_unless_switched_on(tmp_path, monkeypatch, value):
    # A17 is off by default until its reads have had a read-level review (Kevin, queue card D1,
    # 2026-10-07): `rectify align` leaves every arm as aligned and run-all derives no context.
    import argparse
    from rectify.core.commands.align_command import _terminal_tail_postpass
    from rectify.core.splice.terminal_tail_placement import terminal_tail_run_context
    if value is None:
        monkeypatch.delenv('RECTIFY_TERMINAL_TAIL', raising=False)
    else:
        monkeypatch.setenv('RECTIFY_TERMINAL_TAIL', value)
    arm, fasta, annotation = _modern_arm(tmp_path)
    results = {'minimap2': str(arm), 'uLTRA': None}
    args = argparse.Namespace(genome=fasta, annotation=annotation, output_dir=tmp_path)
    assert _terminal_tail_postpass(args, results, 'sample') is None
    assert results == {'minimap2': str(arm), 'uLTRA': None}
    assert not list(tmp_path.glob('*terminal_tail*'))
    assert terminal_tail_run_context(fasta, annotation) is None


def test_align_runs_the_postpass_when_switched_on_and_run_all_agrees(tmp_path, monkeypatch):
    import argparse
    from rectify.core.commands.align_command import _terminal_tail_postpass
    from rectify.core.splice.terminal_tail_placement import terminal_tail_run_context
    monkeypatch.setenv('RECTIFY_TERMINAL_TAIL', '1')
    arm, fasta, annotation = _modern_arm(tmp_path)
    results = {'minimap2': str(arm)}
    args = argparse.Namespace(genome=fasta, annotation=annotation, output_dir=tmp_path)
    context = _terminal_tail_postpass(args, results, 'sample')
    assert context is not None and context == terminal_tail_run_context(fasta, annotation)
    assert results['minimap2'] == str(tmp_path / 'sample.minimap2.terminal_tail.bam')
    with pysam.AlignmentFile(results['minimap2'], 'rb') as bam:
        assert next(bam).has_tag(TAG)
    # The switch does not widen the protocol: short-read, dT-primed and unannotated runs stay out.
    assert terminal_tail_run_context(fasta, annotation, short_read=True) is None
    assert terminal_tail_run_context(fasta, annotation, dt_primed_cdna=True) is None
    assert terminal_tail_run_context(fasta, None) is None


def test_switch_state_binds_the_run_all_selection_receipt(tmp_path, monkeypatch):
    # A run-all resume with the switch off keeps an off-run selection (no rebuild loop) and ignores
    # post-pass files on disk; an on-run selection stands down once the switch is off.
    from shutil import copyfile
    from rectify.core.commands.run.helpers import _collect_per_aligner_bams
    from rectify.core.splice.terminal_tail_placement import (
        run_terminal_tail_postpass, terminal_selection_matches, terminal_tail_run_context,
        write_terminal_selection_receipt,
    )
    arm, fasta, annotation = _modern_arm(tmp_path)
    consensus = tmp_path / 'sample.multialigned.bam'
    copyfile(arm, consensus)
    monkeypatch.delenv('RECTIFY_TERMINAL_TAIL', raising=False)
    write_terminal_selection_receipt(consensus, {'minimap2': str(arm)}, terminal_tail_run_context(fasta, annotation))
    assert terminal_selection_matches(consensus, terminal_tail_run_context(fasta, annotation))
    monkeypatch.setenv('RECTIFY_TERMINAL_TAIL', '1')
    on = terminal_tail_run_context(fasta, annotation)
    assert not terminal_selection_matches(consensus, on)
    selected, counts = run_terminal_tail_postpass(arm, fasta, tmp_path / 'sample.minimap2.terminal_tail.bam', set(), on)
    assert counts['applied'] == 1
    write_terminal_selection_receipt(consensus, {'minimap2': selected}, on)
    assert terminal_selection_matches(consensus, on)
    assert str(_collect_per_aligner_bams('sample', tmp_path, terminal_tail_context=on)['minimap2']) == selected
    monkeypatch.delenv('RECTIFY_TERMINAL_TAIL')
    off = terminal_tail_run_context(fasta, annotation)
    assert off is None
    assert not terminal_selection_matches(consensus, off)
    assert str(_collect_per_aligner_bams('sample', tmp_path, terminal_tail_context=off)['minimap2']) == str(arm)
