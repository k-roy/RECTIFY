"""Selected query paths own boundary edges; no score or result mocks."""
import copy
import random

import pysam
import pytest

from rectify.core.consensus.chimeric_consensus import (
    ChimericSegment, build_chimeric_cigar, build_chimeric_read,
    build_query_ref_map, cigar_to_events, select_best_chimeric,
)
from rectify.core.consensus.consensus import run_consensus_selection
from rectify.core.consensus.sequence import decoded_alignment_copy
from rectify.core.splice.microexon_provenance import (
    load_microexon_calls, record_microexon_draws,
)
from rectify.utils.genome import register_genome_contigs


def _rc(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]

def make_family(reverse=False, micro=False, shared='I', soft=False):
    """Two genuine-source splice paths, each wrong at one different boundary."""
    length, start = 1800, 400
    middle = 20 if micro else 66
    first, last, n1, n2 = 36, 40, 75, 90
    split = 8 if micro else 29
    gap = 2 if shared == 'D' else 0
    rest = middle - split - gap
    j1 = (start + first, start + first + n1)
    j2 = (j1[1] + middle, j1[1] + middle + n2)
    end = j2[1] + last
    rng = random.Random(670920)
    reference = list(''.join(rng.choice('ACGT') for _ in range(length)))
    for s, e in (j1, j2):
        reference[s:s+2] = 'GT'
        reference[e-2:e] = 'AG'
        reference[s-5:s] = 'CCGTC'
        reference[e-7:e-2] = 'AACCA'
    reference = ''.join(reference)
    genome = {'c67P': reference, 'c67M': _rc(reference)}
    register_genome_contigs(genome)
    header = pysam.AlignmentHeader.from_references(list(genome), [length, length])
    op = 1 if shared == 'I' else 2
    truth_ops = [(0, first), (3, n1), (0, split), (op, 2),
                 (0, rest), (3, n2), (0, last)]
    source_ops = {
        'minimap2': [(0, first), (3, n1), (0, split), (op, 2),
                    (0, rest-5), (3, n2), (0, last+5)],
        'uLTRA': [(0, first-3), (3, n1), (0, split+3), (op, 2),
                  (0, rest), (3, n2), (0, last)],
    }
    query, p = '', start
    for code, n in truth_ops:
        if code == 0:
            query += reference[p:p+n]
            p += n
        elif code == 1:
            query += 'TC'
        elif code in (2, 3):
            p += n
    if soft:
        query = 'GCA' + query + 'CGCT'
    def record(ops):
        r = pysam.AlignedSegment(header)
        r.query_name = f'c67_r{int(reverse)}_m{int(micro)}_{shared}_s{int(soft)}'
        r.reference_id = int(reverse)
        r.reference_start = length-end if reverse else start
        r.flag = 16 if reverse else 0
        r.mapping_quality = 55
        r.query_sequence = _rc(query) if reverse else query
        r.query_qualities = [27 + i % 7 for i in range(len(query))]
        if soft:
            ops = [(4, 3)] + ops + [(4, 4)]
        r.cigartuples = ops[::-1] if reverse else ops
        r.set_tag('XN', 1)
        r.set_tag('XT', 2)
        return r
    arms = {arm: record(ops) for arm, ops in source_ops.items()}
    truth = record(truth_ops)
    annotated = {(c, s, e, strand)
                 for c, strand, pairs in [('c67P', '+', [j1, j2]),
                                          ('c67M', '-', [(length-e, length-s) for s,e in (j1,j2)])]
                 for s,e in pairs}
    if micro:
        # Prior validated provenance is fixture input; this is a selection
        # ownership test, not a mocked/live StationB detection claim.
        intron, exon = [j1[0], j2[1]], [j1[1], j2[0]]
        if reverse:
            intron, exon = [length-intron[1], length-intron[0]], [length-exon[1], length-exon[0]]
        for r in arms.values():
            record_microexon_draws(r, [{'intron': intron, 'exons': [exon],
                                       'alternatives': 'fixture_prior_draw'}])
    return arms, truth, genome, annotated


def _junctions(read):
    p = read.reference_start
    out = []
    for op, n in read.cigartuples:
        if op == 3:
            out.append((p, p+n))
        if op in (0, 2, 3, 7, 8):
            p += n
    return out


@pytest.mark.parametrize('reverse', [False, True])
@pytest.mark.parametrize('reorder', [False, True])
@pytest.mark.parametrize('micro', [False, True])
@pytest.mark.parametrize('shared', ['I', 'D'])
@pytest.mark.parametrize('soft', [False, True])
def test_selected_source_path_is_realized(reverse, reorder, micro, shared, soft):
    arms, truth, genome, annotation = make_family(reverse, micro, shared, soft)
    if reorder:
        arms = dict(reversed(list(arms.items())))
    before = {a:r.to_string() for a,r in arms.items()}
    result = select_best_chimeric(arms, genome, annotation)
    assert not result.is_fallback
    assert result.is_chimeric
    anchor = arms[result.anchor_aligner]
    out = build_chimeric_read(anchor, result.chimeric_ref_start, result.chimeric_cigar,
                              result, anchor.header, anchor_read=anchor, aligner_reads=arms)
    assert out.query_sequence == truth.query_sequence
    assert out.query_qualities == truth.query_qualities
    assert out.reference_start == truth.reference_start
    assert out.cigartuples == truth.cigartuples
    assert build_query_ref_map(out) == build_query_ref_map(truth)
    for q,p in build_query_ref_map(out).items():
        assert out.query_sequence[q] == genome[out.reference_name][p]
    assert _junctions(out) == _junctions(truth)
    if micro:
        calls = load_microexon_calls(out)
        assert calls
        selected_source_edges = {
            (event.r_start,event.r_end)
            for segment in result.all_segment_scores
            for event in segment.cigar_events[segment.winning_aligner]
            if event.op == 3
        }
        # Gap-created bridges intentionally have no prior StationB credit.
        assert {tuple(j) for call in calls for j in call['junctions']} == (
            set(_junctions(truth)) & selected_source_edges)
    assert before == {a:r.to_string() for a,r in arms.items()}


@pytest.mark.parametrize('reverse', [False,True])
def test_agreed_microexon_keeps_both_source_junction_credits(reverse):
    arms, truth, genome, annotation = make_family(reverse,micro=True)
    for read in arms.values():
        read.cigartuples = truth.cigartuples
    result = select_best_chimeric(arms,genome,annotation)
    anchor = arms[result.anchor_aligner]
    out = build_chimeric_read(anchor,result.chimeric_ref_start,result.chimeric_cigar,
                              result,anchor.header,anchor_read=anchor,aligner_reads=arms)
    calls = load_microexon_calls(out)
    assert {tuple(j) for call in calls for j in call['junctions']} == set(_junctions(truth))


def _segments(cuts):
    return [ChimericSegment(a,b,kind,winning_aligner=arm) for a,b,kind,arm in cuts]


@pytest.mark.parametrize('chain', [[(3, 23)], [(2, 23)], [(3,20),(2,3)], [(2,3),(3,20)]])
def test_required_reference_edges_are_retained(chain):
    events = {'a': cigar_to_events([(0,10)]+chain+[(0,10)],100)}
    segments = _segments([(0,10,'agreement','a'),(10,20,'agreement','a')])
    assert build_chimeric_cigar(segments,events,{}) == (100,[(0,10)]+chain+[(0,10)])


@pytest.mark.parametrize('chain', [[(3, 23)], [(2, 23)], [(3,20),(2,3)], [(2,3),(3,20)]])
def test_only_already_traversed_alternative_edge_is_omitted(chain):
    events = {'a':cigar_to_events([(0,10)]+chain+[(0,15)],100),
              'b':cigar_to_events([(0,15)]+chain+[(0,10)],100)}
    segments = _segments([(0,15,'interior','a'),(15,25,'agreement','b')])
    frozen = copy.deepcopy(events)
    assert build_chimeric_cigar(segments,events,{}) == (100,[(0,10)]+chain+[(0,15)])
    assert events == frozen


@pytest.mark.parametrize('bad', ['partial', 'nonagreement', 'query_gap', 'regression'])
def test_incompatible_paths_keep_refusal(bad):
    events = {'a':cigar_to_events([(0,10),(3,23),(0,15)],100),
              'b':cigar_to_events([(0,15),(3,23),(0,10)],100)}
    cuts = [(0,15,'interior','a'),(15,25,'agreement','b')]
    if bad == 'partial':
        events['b'] = cigar_to_events([(0,15),(3,25),(0,10)],100)
    elif bad == 'nonagreement':
        cuts[1] = (15,25,'interior','b')
    elif bad == 'query_gap':
        cuts = [(0,15,'interior','a'),(15,20,'agreement','b'),(21,25,'agreement','b')]
    elif bad == 'regression':
        events['b'] = cigar_to_events([(0,25)],90)
    assert build_chimeric_cigar(_segments(cuts),events,{}) == (None,[])


def test_internal_and_terminal_junctions_preserve_one_base_agreement():
    # An actual terminal N with a one-base final exon still needs its edge.
    ops = [(0,10),(3,23),(0,1),(3,29),(0,1)]
    events = {'a':cigar_to_events(ops,100)}
    cuts = [(0,10,'agreement','a'),(10,11,'agreement','a'),(11,12,'agreement','a')]
    assert build_chimeric_cigar(_segments(cuts),events,{}) == (100,ops)


@pytest.mark.parametrize('shared', ['I','D'])
def test_native_selection_encoding_and_order_parity(tmp_path, shared):
    arms, _, genome, annotation = make_family(False, shared=shared)
    header = next(iter(arms.values())).header
    fasta = tmp_path/'genome.fa'
    fasta.write_text(''.join(f'>{c}\n{s}\n' for c,s in genome.items()))
    pysam.faidx(str(fasta))
    expected = {}
    paths = {mode:{} for mode in ('literal','calmd')}
    for arm in arms:
        source = tmp_path/f'{arm}.original.bam'
        with pysam.AlignmentFile(str(source),'wb',header=header) as bam:
            for reverse in (False,True):
                inputs, truth, _, _ = make_family(reverse, shared=shared)
                expected[truth.query_name] = truth
                bam.write(inputs[arm])
        for mode in paths:
            path = tmp_path/f'{arm}.{mode}.bam'
            path.write_bytes(pysam.calmd('-b', *(['-e'] if mode=='calmd' else []),str(source),str(fasta)))
            paths[mode][arm] = str(path)
    baseline = None
    for mode in paths:
        for reorder in (False,True):
            inputs = dict(reversed(list(paths[mode].items()))) if reorder else paths[mode]
            output = tmp_path/f'{mode}.{reorder}.selected.bam'
            run_consensus_selection(inputs,genome,str(output),annotated_junctions=annotation,
                                    n_workers=1,use_chimeric=True)
            with pysam.AlignmentFile(str(output)) as bam:
                observed = {r.query_name:r for r in bam}
            assert observed.keys() == expected.keys()
            states = {}
            for name, r in observed.items():
                truth = expected[name]
                r = decoded_alignment_copy(r,genome)
                assert r.query_sequence == truth.query_sequence
                assert r.query_qualities == truth.query_qualities
                assert r.cigartuples == truth.cigartuples
                assert build_query_ref_map(r) == build_query_ref_map(truth)
                states[name] = (r.reference_name,r.reference_start,r.cigarstring,
                                r.query_sequence,r.qual)
            if baseline is None:
                baseline = states
            assert states == baseline
