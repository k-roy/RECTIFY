"""Unmocked whole-block 2H placement and atomic consumer regressions."""
from array import array
import pickle

import pysam

from rectify.core.splice.block_shift import (
    AtomicReadPlan, apply_atomic_plan, decoded_copy, plan_block_shifts,
)
from rectify.core.splice.junction_refiner import (
    _apply_replacements_to_read, _iter_n_ops, refine_read_junctions,
)
from rectify.utils.genome import register_genome_contigs

RC = str.maketrans('ACGT', 'TGCA')


def block_input(reverse=False, internal=False, encoding='literal', origin='shifted'):
    bases = list(('CGTACGATGCTAGTCGACTG' * 30)[:500])
    bases[80:170] = list('GT' + 'C' * 86 + 'AG')
    bases[170:206] = list('ACGTAGACGTAC' + 'TTGGCCAATTAACCGGATCTGACT')
    bases[186:188] = 'GT'
    bases[192:194] = 'GT'
    bases[284:286] = 'AG'
    forward = ''.join(bases)
    g = forward.translate(RC)[::-1] if reverse else forward
    chrom = 'blockM' if reverse else 'blockP'
    register_genome_contigs({chrom: g})
    header = pysam.AlignmentHeader.from_references([chrom], [500])
    ops, expected_ops = [(0, 40), (3, 90), (0, 16)], [(0, 40), (3, 96), (0, 16)]
    q = forward[40:80] + forward[176:192] if origin == 'shifted' else forward[40:80] + forward[170:186]
    if internal:
        ops += [(3, 100), (0, 24)]
        expected_ops += [(3, 94), (0, 24)]
        q += forward[286:310]
    read = pysam.AlignedSegment(header)
    read.query_name = f'block_{reverse}_{internal}_{encoding}_{origin}'
    read.reference_id = 0
    read.reference_start = 500 - (310 if internal else 186) if reverse else 40
    read.flag = 16 if reverse else 0
    read.mapping_quality = 60
    read.cigartuples = list(reversed(ops)) if reverse else ops
    read.query_sequence = q.translate(RC)[::-1] if reverse else q
    read.query_qualities = array('B', [30 + i % 10 for i in range(len(q))])
    read.set_tag('ZA', array('h', [-2, 300]))
    read.set_tag('ZZ', 'preserve')
    read.set_tag('MD', '0A1')
    read.set_tag('NM', 9)
    read.set_tag('AS', 40)
    read.set_tag('cs', ':2*ac')
    expected = read.__copy__()
    expected.cigartuples = list(reversed(expected_ops)) if reverse else expected_ops
    if reverse and not internal:
        expected.reference_start -= 6
    pool = {chrom: sorted({(s, e) for r in (read, expected) for _, s, e, _ in _iter_n_ops(r)})}
    if encoding == 'equals':
        seq = list(read.query_sequence)
        for qi, ri in read.get_aligned_pairs(matches_only=True):
            if seq[qi] == g[ri]:
                seq[qi] = '='
        qual = read.query_qualities
        read.query_sequence = ''.join(seq)
        read.query_qualities = qual
    return read, expected, g, pool


def production_case(reverse=False, internal=False, encoding='literal', origin='shifted', **kwargs):
    read, expected, g, pool = block_input(reverse, internal, encoding, origin)
    changes = refine_read_junctions(read, pool, set(), g, '-' if reverse else '+', **kwargs)
    out, applied = _apply_replacements_to_read(read, changes, g, '-' if reverse else '+', .25, 15)
    return read, expected, out, changes, applied, g, pool


def test_whole_block_literal_equals_both_strands_and_internal_boundaries():
    for reverse in (False, True):
        for internal in (False, True):
            for encoding in ('literal', 'equals'):
                read, expected, out, plans, applied, g, pool = production_case(reverse, internal, encoding)
                assert applied and len(plans) == 1 and isinstance(plans[0], AtomicReadPlan)
                assert (out.reference_start, out.cigartuples) == (expected.reference_start, expected.cigartuples)
                assert out.query_sequence == expected.query_sequence
                assert out.query_qualities == expected.query_qualities
                assert out.get_aligned_pairs() == expected.get_aligned_pairs()
                assert out.get_tag('ZA') == array('h', [-2, 300]) and out.get_tag('ZZ') == 'preserve'
                assert all(not out.has_tag(t) for t in ('MD', 'NM', 'AS', 'cs'))
                assert read.has_tag('MD')  # original is not mutated
                evidence = plans[0].evidence[0]
                assert evidence.new_matches == evidence.left_anchor == evidence.right_anchor == 16
                assert evidence.old_matches < evidence.new_matches
                assert evidence.search_trials == 101 * (2 if internal else 1)
                assert len(evidence.changed_junction_ordinals) == (2 if internal else 1)
                assert len(plans[0].junction_map) == (2 if internal else 1)
                changed = set(range(evidence.query_start, evidence.query_end))
                before_map = dict(decoded_copy(read, g).get_aligned_pairs(matches_only=True))
                after_map = dict(out.get_aligned_pairs(matches_only=True))
                assert all(before_map[q] == after_map[q] for q in before_map if q not in changed)


def test_correct_stock_and_noisy_blocks_stay():
    for reverse in (False, True):
        for internal in (False, True):
            read, _, out, plans, applied, _, _ = production_case(reverse, internal, origin='stock')
            assert not plans and not applied and out.cigartuples == read.cigartuples
    read, _, g, pool = block_input()
    q = list(read.query_sequence)
    q[45] = next(b for b in 'ACGT' if b != q[45])
    qual = read.query_qualities
    read.query_sequence = ''.join(q)
    read.query_qualities = qual
    assert not refine_read_junctions(read, pool, set(), g, '+')


def test_unsupported_internal_partner_and_annotated_incumbents_hold():
    for reverse in (False, True):
        read, expected, g, pool = block_input(reverse, True)
        primary = -1 if reverse else 0
        target = list(_iter_n_ops(expected))[primary][1:3]
        original = {(s, e) for _, s, e, _ in _iter_n_ops(read)}
        incomplete = {read.reference_name: sorted(original | {target})}
        assert not refine_read_junctions(read, incomplete, set(), g, '-' if reverse else '+')
        annotated = {(read.reference_name, s, e) for s, e in original}
        assert not refine_read_junctions(read, pool, annotated, g, '-' if reverse else '+')
        annotated.update((read.reference_name, s, e) for s, e in pool[read.reference_name])
        assert refine_read_junctions(read, pool, annotated, g, '-' if reverse else '+')


def test_research_policy_cap_and_information_refusals():
    for kwargs in ({'motif_blind': True}, {'hold_margin': 1.0},
                   {'hp_drift_margin': 1.0}, {'microhom_drift_margin': 1.0},
                   {'drift_near_tie_cap': 1.0}, {'drift_positional_gate': 1.0},
                   {'max_candidates_per_nop': 1}, {'max_boundary_shift': 5},
                   {'max_junction_size': 95}):
        assert not production_case(**kwargs)[3], kwargs
    read, _, g, pool = block_input()
    # A six-base exact destination is not enough for the full 101-endpoint
    # search. Do not lower the gate to admit a length-based positive fixture.
    read.cigartuples = [(0, 40), (3, 90), (0, 6)]
    read.query_sequence = read.query_sequence[:46]
    assert not refine_read_junctions(read, pool, set(), g, '+')


def test_atomic_roundtrip_pickle_stale_state_and_no_mixed_plan():
    read, expected, out, plans, applied, g, _ = production_case(encoding='equals')
    assert applied
    recovered = pickle.loads(pickle.dumps(plans[0]))
    same, ok = apply_atomic_plan(read, recovered, g)
    assert ok and same.to_string() == out.to_string()
    for kind in ('position', 'cigar', 'sequence', 'quality', 'name'):
        stale = read.__copy__()
        if kind == 'position': stale.reference_start += 1
        elif kind == 'cigar': stale.cigartuples = [(0, 39), (3, 90), (0, 17)]
        elif kind == 'sequence': stale.query_sequence = 'N' + stale.query_sequence[1:]
        elif kind == 'quality': stale.query_qualities = array('B', [1] * len(stale.query_sequence))
        else: stale.query_name += '_different'
        before = stale.to_string()
        unchanged, ok = apply_atomic_plan(stale, recovered, g)
        assert not ok and unchanged.to_string() == before
    unchanged, ok = _apply_replacements_to_read(read, plans + [(1, 80, 170, 80, 176)], g, '+', .25, 15)
    assert not ok and unchanged.to_string() == read.to_string()


def test_legacy_changed_shared_boundary_is_reserved():
    read, _, g, pool = block_input(internal=True)
    evolved = read.__copy__()
    # Concrete geometry representing a prior successful edit. This tests only
    # reservation, not the biological merits or reachability of that edit.
    evolved.cigartuples = [(0, 40), (3, 90), (0, 16), (3, 101), (0, 24)]
    assert plan_block_shifts(read, evolved, pool, set(), g, '+') is None


def test_whole_window_repeat_refused_even_if_absent_from_pool():
    read, _, g, pool = block_input()
    # Add a second exact placement at +30, without a junction-pool entry.
    duplicate = list(g)
    duplicate[200:216] = read.query_sequence[40:56]
    assert not refine_read_junctions(read, pool, set(), ''.join(duplicate), '+')


def test_clipped_query_frame_and_gapped_blocks():
    for reverse in (False, True):
        read, expected, g, pool = block_input(reverse)
        for r in (read, expected):
            ops = r.cigartuples
            q = r.query_sequence
            r.cigartuples = [(5, 7), (4, 3)] + ops + [(4, 4), (5, 2)]
            r.query_sequence = 'GCA' + q + 'ACTG'
            r.query_qualities = array('B', [31] * len(r.query_sequence))
        plans = refine_read_junctions(read, pool, set(), g, '-' if reverse else '+')
        out, applied = _apply_replacements_to_read(read, plans, g, '-' if reverse else '+', .25, 15)
        assert applied and out.cigartuples == expected.cigartuples
        assert out.get_aligned_pairs() == expected.get_aligned_pairs()
        assert out.query_sequence == expected.query_sequence
    read, _, g, pool = block_input()
    q = read.query_sequence
    read.cigartuples = [(0, 40), (3, 90), (0, 8), (1, 1), (0, 8)]
    read.query_sequence = q[:48] + 'A' + q[48:]
    # Legacy may consider I; direct fallback itself must reject gapped blocks.
    assert plan_block_shifts(read, read, pool, set(), g, '+') is None


def test_contig_identity_reference_and_chimeric_state():
    read, _, out, plans, _, g, pool = production_case()
    plan = plans[0]
    other_header = pysam.AlignmentHeader.from_references(['other_contig'], [500])
    stale = pysam.AlignedSegment.fromstring(read.to_string().replace('blockP', 'other_contig'), other_header)
    assert stale.reference_id == read.reference_id
    assert not apply_atomic_plan(stale, plan, g)[1]
    corrupt = list(g)
    corrupt[176] = 'N'
    assert not apply_atomic_plan(read, plan, ''.join(corrupt))[1]
    for tag, value in (('ms', 3), ('de', .1), ('dv', .2), ('UQ', 9)):
        read.set_tag(tag, value)
    clean, ok = apply_atomic_plan(read, plan, g)
    assert ok and all(not clean.has_tag(t) for t in ('ms', 'de', 'dv', 'UQ'))
    read.set_tag('SA', 'elsewhere,1,+,56M,60,0;')
    assert not refine_read_junctions(read, pool, set(), g, '+')
    read.set_tag('SA', None)
    read.flag |= 2048
    assert not refine_read_junctions(read, pool, set(), g, '+')


def test_apply_time_late_sa_and_motif_reference_changes_fail_closed():
    for reverse in (False, True):
        read, _, _, plans, _, genome, _ = production_case(reverse, True)
        plan = plans[0]
        read.set_tag('SA', 'other,1,+,80M,60,0;')
        before = read.to_string()
        unchanged, applied = apply_atomic_plan(read, plan, genome)
        assert not applied and unchanged.to_string() == before
        read.set_tag('SA', None)
        for ordinal in plan.evidence[0].changed_junction_ordinals:
            _, _, start, end = plan.junction_map[ordinal]
            for position in (start, end - 1):
                corrupt = list(genome)
                corrupt[position] = 'N'
                before = read.to_string()
                unchanged, applied = apply_atomic_plan(read, plan, ''.join(corrupt))
                assert not applied and unchanged.to_string() == before


def test_existing_microexon_provenance_reserves_block_placement():
    for reverse in (False, True):
        for internal in (False, True):
            for encoding in ('literal', 'equals'):
                read, _, g, pool = block_input(reverse, internal, encoding)
                read.set_tag('Xb', 'unknown-or-older-format')
                before = read.to_string()
                assert not refine_read_junctions(read, pool, set(), g, '-' if reverse else '+')
                assert read.to_string() == before
    read, _, _, plans, _, g, _ = production_case()
    read.set_tag('Xb', 'newly-added-after-planning')
    before = read.to_string()
    unchanged, applied = apply_atomic_plan(read, plans[0], g)
    assert not applied and unchanged.to_string() == before
