"""ISSUE-016 — Station C's census must account for EVERY N-op.

On real corrected output (850388e, 145k reads) the tester matched the 797
junctions RECTIFY created against ``<prefix>.pool_gate.tsv`` and found 693
ABSENT at every op budget (0/2/4) and support gate (1/2). Neither knob was
binding, because the table cannot show two things: it lists NON-annotated
junctions only (2F/2H create mostly annotated ones), and it keys junctions at
the LEFTMOST ambiguity-equivalent coordinate while the BAM's N-op sits at the
motif coordinate. On the 100-read panel: 10 created, 0 in the table, 231 of
271 censused junctions annotated.

These tests pin the ledger: refusal reasons per N-op, the annotated table,
raw coordinates beside the leftmost keys, and ``--attribute``'s statuses.
The genome is the hermetic 139 bp locus from test_station_c_human_gtf.
"""
import json

import pysam
import pytest

from rectify.core.consensus.station_c import (
    ATTRIBUTION_STATUSES,
    CensusLedger,
    PoolGateConfig,
    _anchor_run,
    _anchor_walk,
    _canonicalize,
    _refusal_reason,
    census_bam,
    load_junction_list,
    pool_gate,
    write_pool_gate_outputs,
)

EX1 = 'ACGTTGCATTACGGCATTGCACGTTACGGTACCATGCATA'          # 40
INTRON = 'GT' + 'CATGACTGACTTGCACGATTGCAAGTACCTGATGCACGTTGCACGTACCTGACTG'[:56] + 'AG'  # 59
EX2 = 'TTGCACGTAACCGGTTACGATCGGATCCATTACGGCATGA'          # 40
CHRSEQ = EX1 + INTRON + EX2
GENOME = {'chr5': CHRSEQ}
J = (40, 99)                     # the annotated GT..AG intron, 0-based half-open
NOVEL = (20, 79)                 # an unannotated N-op planted by the reads below


@pytest.fixture
def gtf(tmp_path):
    p = tmp_path / 'mini.gtf'
    attrs = 'gene_id "G1"; transcript_id "T1"; exon_number {};'
    p.write_text(
        f'chr5\tX\texon\t1\t{J[0]}\t.\t+\t.\t{attrs.format(1)}\n'
        f'chr5\tX\texon\t{J[1] + 1}\t{len(CHRSEQ)}\t.\t+\t.\t{attrs.format(2)}\n'
    )
    return p


def _qlen(cigar):
    return sum(ln for op, ln in cigar if op in (0, 1, 4, 7, 8))


def _read(header, name, start, cigar):
    a = pysam.AlignedSegment(header)
    a.query_name = name
    a.reference_id = 0
    a.reference_start = start
    a.cigartuples = cigar
    a.query_sequence = (CHRSEQ * 3)[start:start + _qlen(cigar)]   # length is what matters
    a.mapping_quality = 60
    return a


@pytest.fixture
def bam(tmp_path):
    """Six N-ops with one outcome each (min_anchor 8, max_ops 2, max_bp 30)."""
    header = pysam.AlignmentHeader.from_dict(
        {'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'chr5', 'LN': len(CHRSEQ)}]})
    reads = [
        # A: clean annotated junction, two reads -> censused, annotated table
        _read(header, 'A0', 0, [(0, 40), (3, 59), (0, 40)]),
        _read(header, 'A1', 0, [(0, 40), (3, 59), (0, 40)]),
        # B: left run 5 then the read starts -> refused L=read_end
        _read(header, 'B', 35, [(0, 5), (3, 59), (0, 40)]),
        # C: left run 5 behind a soft clip -> refused L=softclip
        _read(header, 'C', 35, [(4, 30), (0, 5), (3, 59), (0, 40)]),
        # D: right walk 3M I 2M D 2M I -> third indel op exceeds max_ops=2 with
        #    only 7 aligned bases reached -> refused R=indel_ops
        _read(header, 'D', 0, [(0, 40), (3, 59), (0, 3), (1, 1), (0, 2), (2, 1),
                               (0, 2), (1, 1), (0, 30)]),
        # E: right walk 4M then a 31 bp deletion -> bp budget exceeded -> R=indel_bp
        _read(header, 'E', 0, [(0, 40), (3, 59), (0, 4), (2, 31), (0, 5)]),
        # F: the same annotated junction reached over a 1I compensating indel
        #    (2H's signature) -> censused, adj_indel_r = I1
        _read(header, 'F', 0, [(0, 40), (3, 59), (0, 10), (1, 1), (0, 29)]),
        # G: an unannotated junction at NOVEL -> reported in the table
        _read(header, 'G', 0, [(0, 20), (3, 59), (0, 40)]),
    ]
    p = tmp_path / 'ledger.bam'
    with pysam.AlignmentFile(p, 'wb', header=header) as out:
        for r in reads:
            out.write(r)
    return p


# ---------------------------------------------------------------------------
# The walk and its reasons
# ---------------------------------------------------------------------------

def test_anchor_run_is_the_two_tuple_view_of_anchor_walk():
    ops = [(0, 40), (3, 59), (0, 3), (1, 1), (0, 2), (2, 1), (0, 2), (1, 1), (0, 30)]
    for direction in (-1, +1):
        assert _anchor_run(ops, 1, direction) == _anchor_walk(ops, 1, direction)[:2]
    assert _anchor_walk(ops, 1, -1) == (40, '', 'read_end')
    assert _anchor_walk(ops, 1, +1) == (7, 'I1', 'indel_ops')
    assert _anchor_walk([(0, 40), (3, 59), (0, 4), (2, 31), (0, 5)], 1, +1) == (4, 'D31', 'indel_bp')
    assert _anchor_walk([(4, 30), (0, 5), (3, 59), (0, 40)], 2, -1) == (5, '', 'softclip')
    assert _anchor_walk([(0, 40), (3, 59), (0, 6), (3, 20), (0, 40)], 1, +1) == (6, '', 'n_op')


def test_refusal_reason_names_the_short_sides():
    assert _refusal_reason(5, 40, 'softclip', 'read_end', 8) == 'L=softclip'
    assert _refusal_reason(40, 7, 'read_end', 'indel_ops', 8) == 'R=indel_ops'
    assert _refusal_reason(3, 4, 'n_op', 'read_end', 8) == 'L=n_op+R=read_end'


# ---------------------------------------------------------------------------
# The ledger through census_bam
# ---------------------------------------------------------------------------

def test_census_ledger_accounts_for_every_n_op(bam):
    ledger = CensusLedger()
    J_ = census_bam(str(bam), GENOME, PoolGateConfig(), ledger=ledger)
    assert ledger.n_reads == 8
    assert ledger.n_ops_seen == 8
    assert ledger.n_ops_censused + ledger.n_ops_refused == ledger.n_ops_seen
    assert ledger.n_ops_refused == 4
    assert ledger.reason_counts() == {
        'L=read_end': 1, 'L=softclip': 1, 'R=indel_bp': 1, 'R=indel_ops': 1}
    # Every refusal is keyed on the RAW N-op coordinates of the annotated intron.
    assert set(ledger.refusals) == {('chr5', 40, 99)}
    assert ledger.best_anchor[('chr5', 40, 99)] == 7          # read D came closest
    # The censused key remembers which raw coordinates its reads carried.
    ann_key = ('chr5', *_canonicalize(CHRSEQ, 40, 99, 30))
    assert ann_key in J_ and J_[ann_key]['support'] == 3      # A0, A1, F
    assert ledger.majority_raw(ann_key) == (40, 99, 1)
    assert ledger.majority_raw(('chr5', 1, 2)) == (1, 2, 0)   # unknown key: itself


def test_census_bam_without_a_ledger_is_unchanged(bam):
    assert census_bam(str(bam), GENOME, PoolGateConfig()) == \
        census_bam(str(bam), GENOME, PoolGateConfig(), ledger=CensusLedger())


# ---------------------------------------------------------------------------
# pool_gate: annotated table, raw coordinates, refusal rows, attribution
# ---------------------------------------------------------------------------

def _keys(rows):
    return {(r['chrom'], r['start'], r['end']) for r in rows}


def test_annotated_junctions_are_listed_not_just_counted(bam, gtf):
    rows, summary = pool_gate(str(bam), GENOME, gtf)
    ann_key = ('chr5', *_canonicalize(CHRSEQ, 40, 99, 30))
    assert ann_key not in _keys(rows), "the main table stays non-annotated"
    assert summary['n_annotated'] == 1
    ann_rows = summary['tables']['annotated']
    assert [(r['chrom'], r['start'], r['end']) for r in ann_rows] == [ann_key]
    r = ann_rows[0]
    assert (r['start_raw'], r['end_raw']) == (40, 99)
    assert r['support'] == 3 and r['adj_indel_r'] == 'I1' and r['n_adj_indel'] == 1


def test_reported_rows_carry_raw_coordinates(bam, gtf):
    rows, _ = pool_gate(str(bam), GENOME, gtf)
    assert len(rows) == 1
    r = rows[0]
    assert (r['start_raw'], r['end_raw']) == NOVEL
    assert (r['start'], r['end']) == _canonicalize(CHRSEQ, *NOVEL, 30)
    assert r['n_raw_variants'] == 1


def test_refusal_rows_attribute_each_refused_raw_junction(bam, gtf):
    _, summary = pool_gate(str(bam), GENOME, gtf)
    ref = summary['tables']['census_refusals']
    assert len(ref) == 1
    r = ref[0]
    assert (r['chrom'], r['start'], r['end']) == ('chr5', 40, 99)
    assert r['n_ops'] == 4 and r['best_anchor'] == 7
    assert set(r['reasons'].split(';')) == {
        'L=read_end:1', 'L=softclip:1', 'R=indel_bp:1', 'R=indel_ops:1'}
    assert r['annotated'] == 1 and r['censused_elsewhere'] == 1 and r['in_table'] == 0
    c = summary['census']
    assert (c['n_ops_seen'], c['n_ops_censused'], c['n_ops_refused']) == (8, 4, 4)
    assert c['n_refused_junctions'] == 1 and c['n_refused_junctions_censused_elsewhere'] == 1


def test_attribution_statuses(bam, gtf):
    novel = ('chr5', *NOVEL)
    listed = [('chr5', 40, 99), novel, ('chr5', 5, 30), ('chr9', 1, 2)]
    rows, summary = pool_gate(str(bam), GENOME, gtf, attribute=listed)
    att = {(r['chrom'], r['start'], r['end']): r for r in summary['tables']['attribution']}
    assert att[('chr5', 40, 99)]['status'] == 'annotated'
    assert att[('chr5', 40, 99)]['support'] == 3
    assert 'L=softclip:1' in att[('chr5', 40, 99)]['reasons']   # refused elsewhere, still shown
    assert att[novel]['status'] == 'reported'
    assert att[novel]['verdict'] == rows[0]['verdict']
    assert att[('chr5', 5, 30)]['status'] == 'not_seen'
    assert att[('chr9', 1, 2)]['status'] == 'chrom_missing'
    assert summary['census']['attribution'] == {
        'reported': 1, 'annotated': 1, 'not_seen': 1, 'chrom_missing': 1}
    assert all(s in ATTRIBUTION_STATUSES for s in summary['census']['attribution'])


def test_attribution_refused_when_no_read_anchors_it(tmp_path, gtf):
    header = pysam.AlignmentHeader.from_dict(
        {'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'chr5', 'LN': len(CHRSEQ)}]})
    p = tmp_path / 'refused.bam'
    with pysam.AlignmentFile(p, 'wb', header=header) as out:
        out.write(_read(header, 'C', 35, [(4, 30), (0, 5), (3, 59), (0, 40)]))
    _, summary = pool_gate(str(p), GENOME, gtf, attribute=[('chr5', 40, 99)])
    a = summary['tables']['attribution'][0]
    assert a['status'] == 'refused' and a['reasons'] == 'L=softclip:1'


def test_attribution_is_keyed_like_the_census(bam, gtf):
    """A listed junction given at a non-leftmost equivalent coordinate still
    resolves: the attribution canonicalizes exactly as the census does."""
    s, e = NOVEL
    l_amb = s - _canonicalize(CHRSEQ, s, e, 30)[0]
    if l_amb == 0:
        pytest.skip('NOVEL has no left ambiguity on this sequence')
    _, summary = pool_gate(str(bam), GENOME, gtf, attribute=[('chr5', s - l_amb, e - l_amb)])
    assert summary['tables']['attribution'][0]['status'] == 'reported'


# ---------------------------------------------------------------------------
# Outputs + the list loader
# ---------------------------------------------------------------------------

def test_write_outputs_adds_the_side_tables_and_keeps_the_json_small(bam, gtf, tmp_path):
    listed = [('chr5', 40, 99), ('chr5', *NOVEL)]
    rows, summary = pool_gate(str(bam), GENOME, gtf, attribute=listed)
    tsv, js = write_pool_gate_outputs(rows, summary, tmp_path / 's1')
    head = tsv.read_text().splitlines()[0].split('\t')
    assert head[-3:] == ['start_raw', 'end_raw', 'n_raw_variants']
    assert head[:3] == ['chrom', 'start', 'end']
    for name in ('s1.pool_gate.annotated.tsv', 's1.census_refusals.tsv', 's1.attribution.tsv'):
        assert (tmp_path / name).exists(), name
    ann = (tmp_path / 's1.pool_gate.annotated.tsv').read_text().splitlines()
    assert ann[0].split('\t')[:5] == ['chrom', 'start', 'end', 'start_raw', 'end_raw']
    assert len(ann) == 2
    with open(js) as fh:
        blob = json.load(fh)
    assert blob['tables']['annotated'].endswith('s1.pool_gate.annotated.tsv')
    assert blob['census']['n_ops_refused'] == 4
    # The in-memory summary handed in is not mutated into paths.
    assert isinstance(summary['tables']['annotated'], list)


def test_write_outputs_without_attribution_writes_no_attribution_file(bam, gtf, tmp_path):
    rows, summary = pool_gate(str(bam), GENOME, gtf)
    _, js = write_pool_gate_outputs(rows, summary, tmp_path / 's2')
    assert not (tmp_path / 's2.attribution.tsv').exists()
    with open(js) as fh:
        assert json.load(fh)['tables']['attribution'] is None


def test_load_junction_list_formats(tmp_path):
    j = tmp_path / 'a.json'
    j.write_text(json.dumps([['chr5', 40, 99], {'chrom': 'chr5', 'start': 20, 'end': 79},
                             'chr5:1-2', ['chr5', 40, 99]]))
    assert load_junction_list(j) == [('chr5', 40, 99), ('chr5', 20, 79), ('chr5', 1, 2)]
    d = tmp_path / 'b.json'
    d.write_text(json.dumps({'created': ['chr5:40-99', ['chr5', 20, 79]]}))
    assert load_junction_list(d) == [('chr5', 40, 99), ('chr5', 20, 79)]
    k = tmp_path / 'c.json'
    k.write_text(json.dumps({'chr5:40-99': 3, 'chr5:20-79': 1}))
    assert load_junction_list(k) == [('chr5', 40, 99), ('chr5', 20, 79)]
    ev = tmp_path / 'fpfn_events.tsv'
    ev.write_text('category\tread\tchrom\torig_junction\tnew_junction\tmotif_orig\tmotif_new\n'
                  'TP_rescue_annot\tr1\tchr5\t\t40-99\t\tGT-AG\n'
                  'FN_drift\tr2\tchr5\t20-79\t\tGT-AG\t\n'
                  'FP_added_nov\tr3\tchr5\t\t20-79\t\tGT-AG\n')
    assert load_junction_list(ev) == [('chr5', 40, 99), ('chr5', 20, 79)]
    headed = tmp_path / 'h.tsv'
    headed.write_text('chrom\tstart\tend\tnote\nchr5\t40\t99\tx\n')
    assert load_junction_list(headed) == [('chr5', 40, 99)]
    bed = tmp_path / 'j.bed'
    bed.write_text('# comment\nchr5\t40\t99\nchr5\t20\t79\tname\n')
    assert load_junction_list(bed) == [('chr5', 40, 99), ('chr5', 20, 79)]
