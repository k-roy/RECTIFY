"""Supplied interior query placements must not lose to motif-only rewards."""
from array import array
import copy
import random

import pysam

from rectify.core.consensus import chimeric_consensus as cc


def rc(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def make_family(seed=80911, reverse=False, family='non_equivalent_window', hard=False, explicit=False):
    """Two complementary two-intron alignments of one physical query."""
    rng = random.Random(seed)
    g = list(''.join(rng.choice('ACGT') for _ in range(2200)))
    start, a, b, c, d, stop = 401, 444, 541, 612, 725, 764
    for left, right in ((a, b), (c, d)):
        g[left:left+2], g[right-2:right] = 'GT', 'AG'
    g[a-4:a], g[b-6:b] = 'GTTA', 'AGGGAG'
    if family == 'equal_product_control':
        g[a-3:a] = g[b-3:b]
    elif family == 'no_window_control':
        g[a-4] = 'C'
    g = ''.join(g)
    query = g[start:a] + g[b:b+23] + 'CT' + g[b+23:c] + g[d:stop]
    genome = {'placement': rc(g) if reverse else g}
    query = rc(query) if reverse else query
    header = pysam.AlignmentHeader.from_references(['placement'], [len(g)])
    cigar = {
        'minimap2': [(0,43),(3,97),(0,23),(1,2),(0,45),(3,113),(0,42)],
        'uLTRA': [(0,40),(3,97),(0,26),(1,2),(0,48),(3,113),(0,39)],
        'truth': [(0,43),(3,97),(0,23),(1,2),(0,48),(3,113),(0,39)],
    }
    reads = {}
    for arm, ops in cigar.items():
        read = pysam.AlignedSegment(header)
        read.query_name = f'placement_{seed}_{int(reverse)}_{family}'
        read.flag = 16 if reverse else 0
        read.reference_id = 0
        read.reference_start = len(g)-stop if reverse else start
        read.mapping_quality = 60
        read.cigartuples = list(reversed(ops)) if reverse else ops
        read.query_sequence = query
        read.query_qualities = [23+i%17 for i in range(len(query))]
        if explicit:
            replaced = []
            for ev in cc.cigar_to_events(read.cigartuples, read.reference_start):
                if ev.op != 0:
                    replaced.append((ev.op, ev.length))
                else:
                    for offset in range(ev.length):
                        op = 7 if query[ev.q_start+offset] == genome['placement'][ev.r_start+offset] else 8
                        if replaced and replaced[-1][0] == op:
                            replaced[-1] = (op, replaced[-1][1]+1)
                        else:
                            replaced.append((op,1))
            read.cigartuples = replaced
        if hard:
            read.cigartuples = [(5,7)] + read.cigartuples + [(5,11)]
        read.set_tag('ZZ', array('H', [0,65535]))
        reads[arm] = read
    annotations = {('placement',len(g)-r,len(g)-l,'-') if reverse else ('placement',l,r,'+')
                   for l,r in ((a,b),(c,d))}
    return genome, reads, annotations


def select(reads, genome, annotations=None):
    source = {k:v for k,v in reads.items() if k != 'truth'}
    result = cc.select_best_chimeric(source, genome, annotations)
    anchor = source[result.anchor_aligner]
    output = cc.build_chimeric_read(anchor, result.chimeric_ref_start,
                                   result.chimeric_cigar, result, anchor.header,
                                   anchor_read=anchor, aligner_reads=source)
    return result, output


def query_map(read):
    return {q:r for q,r in read.get_aligned_pairs() if q is not None}


def known_matches(read, genome):
    return sum(read.query_sequence[q] == genome[read.reference_name][r]
               for q,r in query_map(read).items() if r is not None)


def test_supplied_query_evidence_both_strands_orders_and_cigar_spellings():
    for seed in (80911,92743):
        for reverse in (False,True):
            for explicit in (False,True):
                genome, reads, annotations = make_family(seed,reverse,explicit=explicit)
                for backward in (False,True):
                    source = dict(reversed(list(reads.items()))) if backward else reads
                    saved = {k:v.to_string() for k,v in source.items()}
                    result, output = select(source,genome,annotations)
                    assert not result.is_fallback
                    assert query_map(output) == query_map(reads['truth'])
                    assert known_matches(output,genome) == 153
                    assert saved == {k:v.to_string() for k,v in source.items()}
                    assert output.query_sequence == reads['truth'].query_sequence
                    assert output.query_qualities == reads['truth'].query_qualities
                    if not reverse:
                        changed = [s for s in result.all_segment_scores if s.selection_reason == 'query_evidence']
                        assert changed
                        for segment in changed:
                            assert segment.scores['uLTRA'].dominated_by == ['minimap2']
                            assert segment.scores['uLTRA'].n_mismatches == 3
                            assert segment.scores['minimap2'].n_matches == 3
                            assert cc._classify_term(segment.scores['minimap2'],segment.scores['uLTRA']) == 'query_evidence'


def test_equal_product_ambiguity_is_not_dominated():
    for reverse in (False,True):
        genome, reads, annotations = make_family(reverse=reverse,family='equal_product_control')
        result, output = select(reads,genome,annotations)
        assert known_matches(output,genome) == 153
        start,end = (112,115) if reverse else (40,43)
        equivalent = [s for s in result.all_segment_scores if s.q_start < end and s.q_end > start]
        assert not any(score.dominated_by for s in equivalent for score in s.scores.values())
        assert not any(s.selection_reason == 'query_evidence' for s in equivalent)


def test_hardclip_and_array_metadata_stay_in_original_frame():
    for reverse in (False,True):
        genome, reads, annotations = make_family(reverse=reverse,hard=True)
        result, output = select(reads,genome,annotations)
        assert output.cigartuples[0] == (5,7) and output.cigartuples[-1] == (5,11)
        assert query_map(output) == query_map(reads['truth'])
        assert output.get_tag('ZZ') == array('H',[0,65535])
        assert output.get_tag('ZZ').typecode == 'H'


def test_receipt_reads_bases_for_M_equals_X_and_unknowns():
    for op in (0,7,8):
        events = cc.cigar_to_events([(op,4)],0)
        actual = cc.score_segment(events,'interior','c',{'c':'ACGT'},query_sequence='ATGN')
        assert (actual.n_matches,actual.n_mismatches,actual.n_unknown) == (2,1,1)
        absent = cc.score_segment(events,'interior','c',{'c':'ACGT'})
        assert (absent.n_matches,absent.n_mismatches,absent.n_unknown) == (0,0,4)
        clipped = cc.score_segment(events,'interior','c',{'c':'AC'},query_sequence='ACGT')
        assert (clipped.n_matches,clipped.n_mismatches,clipped.n_unknown) == (2,0,2)


def test_gap_unknown_and_terminal_comparisons_keep_existing_policy():
    # No candidate score is invented here: test the evidence-only precondition.
    clean = cc.cigar_to_events([(0,6)],0)
    assert cc._gap_free_match_profile(clean,1,4,'ACGTAC','ACGTAC') == (True,True,True)
    for cigar in ([(0,1),(1,1),(0,4)],[(0,1),(4,1),(0,4)],
                  [(0,1),(2,1),(0,5)],[(0,4),(2,1),(0,2)]):
        assert cc._gap_free_match_profile(cc.cigar_to_events(cigar,0),1,4,'ACGTAC','ACGTACGT') is None
    for query,reference in [(None,'ACGTAC'),('ACNTAC','ACGTAC'),('ACGTAC','ACNTAC'),('AC','ACGTAC'),('ACGTAC','')]:
        assert cc._gap_free_match_profile(clean,1,4,query,reference) is None
    for position in ('five_prime','three_prime'):
        events = cc.cigar_to_events([(4,2),(0,4)],0)
        score = cc.score_segment(events,position,'c',{'c':'ACGT'},query_sequence='TTTTTT')
        assert score.score == (-4 if position == 'five_prime' else 6)
        genome, reads, _ = make_family()
        segment = cc.ChimericSegment(40,43,position)
        cc._mark_dominated_placements(segment,{k:cc.cigar_to_events(v.cigartuples,v.reference_start) for k,v in reads.items()},reads,genome['placement'])
        assert segment.scores == {}


def test_conflicting_per_base_tradeoff_and_different_query_are_incomparable():
    header = pysam.AlignmentHeader.from_references(['c'],[20])
    reads={}
    for name,start in [('first',0),('second',5)]:
        r=pysam.AlignedSegment(header); r.query_name='same'; r.reference_id=0
        r.reference_start=start; r.cigartuples=[(0,3)]; r.query_sequence='AAA'
        reads[name]=r
    events={k:cc.cigar_to_events(v.cigartuples,v.reference_start) for k,v in reads.items()}
    # Match profiles TTF and FTT: neither dominates, despite shared query.
    segment=cc.ChimericSegment(0,3,'interior',scores={k:cc.SegmentScore(k,0) for k in reads})
    cc._mark_dominated_placements(segment,events,reads,'AACGGCAA')
    assert not any(s.dominated_by for s in segment.scores.values())
    reads['second'].query_sequence='CAA'
    cc._mark_dominated_placements(segment,events,reads,'AACGGCAA')
    assert not any(s.dominated_by for s in segment.scores.values())


def make_three_arm_family(reverse=False):
    """A dominated score winner plus an incomparable, higher-scored frontier arm."""
    rng=random.Random(740703)
    g=list(''.join(rng.choice('ACGT') for _ in range(500)))
    # q40:46 = GTATAC. Old profiles are TFFTFF (A), TTFTFF (B), FFTTFF (C).
    # A gets canonical+5; C gets canonical+5/noncanonical-3; B's N belongs
    # to the following agreement edge, hence its interval policy score is0.
    g[140:146]='GTCTCA'
    g[198:205]='AGAAACC'
    g[238:246]='AGGCCTCA'
    g=''.join(g)
    query=g[100:140]+'GTATAC'+g[246:286]
    genome={'three':rc(g) if reverse else g}
    query=rc(query) if reverse else query
    header=pysam.AlignmentHeader.from_references(['three'],[len(g)])
    source={}
    # This order is material only to ties, and the demonstrated scores differ.
    for name,cigar in [('minimap2',[(0,40),(3,100),(0,46)]),
                       ('uLTRA',[(0,40),(3,60),(0,3),(3,40),(0,43)]),
                       ('GMAP',[(0,46),(3,100),(0,40)])]:
        r=pysam.AlignedSegment(header);r.query_name='three_arm_'+str(int(reverse))
        r.flag=16 if reverse else 0;r.reference_id=0;r.reference_start=500-286 if reverse else 100
        r.cigartuples=list(reversed(cigar)) if reverse else cigar
        r.mapping_quality=60;r.query_sequence=query;r.query_qualities=[35]*len(query)
        source[name]=r
    return genome,source


def test_changed_winner_itself_dominates_old_policy_winner_with_three_arms():
    genome,reads=make_three_arm_family()
    result,output=select(reads,genome)
    segment=next(s for s in result.all_segment_scores if s.position=='interior')
    assert (segment.q_start,segment.q_end)==(40,46)
    assert {k:v.score for k,v in segment.scores.items()}=={'minimap2':5,'uLTRA':2,'GMAP':0}
    assert segment.scores['minimap2'].dominated_by==['GMAP']
    assert segment.scores['uLTRA'].dominated_by==[]
    assert segment.winning_aligner=='GMAP' and segment.selection_reason=='query_evidence'
    assert not result.is_fallback
    assert query_map(output)==query_map(reads['GMAP'])
    old_cells={q for q,r in query_map(reads['minimap2']).items()
               if r is not None and reads['minimap2'].query_sequence[q]==genome['three'][r]}
    new_cells={q for q,r in query_map(output).items()
               if r is not None and output.query_sequence[q]==genome['three'][r]}
    assert old_cells < new_cells


def make_fallback_frontier_family(reverse=False):
    """Independent refuter's actual three-arm seam/fallback counterexample."""
    rng=random.Random(749220)
    ref=list(''.join(rng.choice('ACGT') for _ in range(1200)))
    middle='ACGTCAGTGACTCGAGTACCGTGA'
    positions={'A':250,'B':450,'C':650}
    profiles={'A':[True]*16+[False]*8,'B':[True]*8+[False]*8+[True]*8,'C':[False]*8+[True]*8+[False]*8}
    ref[150:152]='GT';ref[898:900]='AG'
    for arm,pos in positions.items():
        ref[pos:pos+24]=[x if ok else {'A':'C','C':'G','G':'T','T':'A'}[x]
                         for x,ok in zip(middle,profiles[arm])]
        ref[pos-2:pos]='TT' if arm=='A' else 'AG'
        ref[pos+24:pos+26]='AA' if arm=='A' else 'GT'
    ref=''.join(ref);query=ref[110:150]+middle+ref[900:940]
    if reverse:ref,query=rc(ref),rc(query)
    genome={'fallback_frontier':ref}
    h=pysam.AlignmentHeader.from_references(['fallback_frontier'],[1200])
    arms={}
    for arm,pos in positions.items():
        r=pysam.AlignedSegment(h);r.query_name='frontier';r.reference_id=0
        r.reference_start=260 if reverse else 110;r.flag=16 if reverse else 0
        r.mapping_quality=60
        ops=[(0,40),(3,pos-150),(0,24),(3,900-pos-24),(0,40)]
        r.cigartuples=ops[::-1] if reverse else ops
        r.query_sequence=query;r.query_qualities=[29]*104
        r.set_tag('ZZ',array('H',[2,19]));r.set_tag('XO','rev' if reverse else 'fwd')
        arms[arm]=r
    return genome,arms


def test_final_fallback_and_new_assembly_cannot_lose_legacy_exact_cells():
    cases=[(False,('A','C','B'),'B',True),
           (True,('C','A','B'),'C',False),(True,('C','B','A'),'C',False)]
    for reverse,order,legacy_arm,legacy_fallback in cases:
        genome,reads=make_fallback_frontier_family(reverse)
        result,output=select({arm:reads[arm] for arm in order},genome)
        assert query_map(output)==query_map(reads[legacy_arm])
        assert result.is_fallback==legacy_fallback
        assert output.query_sequence==reads[legacy_arm].query_sequence
        assert output.query_qualities==reads[legacy_arm].query_qualities
        assert output.get_tag('ZZ')==reads[legacy_arm].get_tag('ZZ')
        assert output.get_tag('XO')==reads[legacy_arm].get_tag('XO')


def test_final_comparison_refuses_gap_clip_unknown_and_quality_tradeoffs():
    genome,reads,_=make_family()
    candidate=cc._single_aligner_result('truth',reads,reads['truth'].query_name)
    legacy=cc._single_aligner_result('uLTRA',reads,reads['uLTRA'].query_name)
    assert cc._result_dominates_legacy(candidate,legacy,reads,genome)
    for change in ('quality','query','deletion','clip'):
        inputs=copy.deepcopy(reads)
        if change=='quality':inputs['truth'].query_qualities=[10]*len(inputs['truth'].query_sequence)
        elif change=='query':inputs['truth'].query_sequence='N'+inputs['truth'].query_sequence[1:]
        elif change=='deletion':
            # No altered query placement may buy evidence by changing a D.
            ops=inputs['truth'].cigartuples;inputs['truth'].cigartuples=ops[:1]+[(2,1)]+ops[1:]
        else:
            ops=inputs['truth'].cigartuples;inputs['truth'].cigartuples=[(4,1),(0,42)]+ops[1:]
        altered=cc._single_aligner_result('truth',inputs,inputs['truth'].query_name)
        assert not cc._result_dominates_legacy(altered,legacy,inputs,genome)
