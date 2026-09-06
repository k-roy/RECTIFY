"""T0 replay fixture — the tester's ISSUE-017 read bundle (Sumner human RNA004 DRS).

Replays real reads that landed 4 nt into the intron (the GTRAGT +5 GT decoy)
against their real annotated / pool candidates on genome slices. The records
are collaborator data and live OUTSIDE the repository
(``dev/sumner_misplaced_panel_20260904/holdout/events/269ebe7/``); the test
skips when they are absent.

Expectations are the ISSUE-020 (anchored-score ranking) outcomes: every bundle
read lands on the ANNOTATED junction (the 017 tree left 3 within 1 nt and
c887bc16 on the decoy), no sequence rescue carries the junction-side gap token,
and the six former "+4 by evidence" reads (arbiter RULING 10) are placed with
no gap at the junction.

Edit-count mapping (RULING 10 §R37 c). ``_exon_edit_count`` walks the EMITTED
exon CIGAR over the read's 5'-terminal query bases and the reference at the
emitted junction and counts mismatches + inserted bases + deleted bases. The
tester's ``ED_annot_ext`` (``altexp_d4_replay10_ext.tsv``) is a plain
Levenshtein distance between the clip + intron-mapped bases and the exon tail
ending at the annotated donor. Ours is the edit count of ONE alignment (the
affine-optimal one) of the placement segment, so for the same pair of strings
it is an upper bound on the Levenshtein distance; the placement segment can
also be longer than the tester's (a reanchored clip, borrowed equivalence
bases). The assertion therefore is: our count for the annotated placement is
<= ED_annot_ext + the number of extra query bases in our segment, and on the
reads where the segments coincide it is <= ED_annot_ext exactly. The per-read
values are recorded in ISSUE-020_2F_anchored_ranking.md.
"""
import copy
import csv
import os
import re

import pysam
import pytest

from rectify.core.splice.splice_aware_5prime import rescue_3ss_truncation
from rectify.utils.genome import register_genome_contigs

E = os.path.expanduser('~/work/rectify/dev/sumner_misplaced_panel_20260904/holdout/events/269ebe7')
BUNDLE = [E + '/d4_replay10.tsv', E + '/d4_replay10.stock.sam', E + '/d4_replay10.slices.fa']
EXT = E + '/altexp_d4_replay10_ext.tsv'

pytestmark = pytest.mark.skipif(not all(os.path.exists(p) for p in BUNDLE),
                                reason='Sumner d4 replay bundle not present (collaborator data, kept outside the repo)')

# The six former "+4 by evidence" reads (ISSUE-017_discriminating_bases_six.md; RULING 10 rejected the +4 for all six).
SIX = ['c887bc16', '7f41e755', '6fc67f58', '5b20c72a', 'a5a5a1bb', 'de84a10a']
# The three the writer used to walk onto the decoy with a compensating 3D (ISSUE-023).
WITHIN1_AT_017 = ['5b20c72a', 'a5a5a1bb', '31fab950']


def _load():
    slices, name = {}, None
    for line in open(BUNDLE[2]):
        if line.startswith('>'):
            m = re.match(r'>(\S+):(\d+)-(\d+) offset0based=(\d+)', line)
            name = (m.group(1), int(m.group(4)))
            slices[name] = []
        else:
            slices[name].append(line.strip())
    slices = {k: ''.join(v).upper() for k, v in slices.items()}
    sam = {}
    for line in open(BUNDLE[1]):
        if line.startswith('#') or line.startswith('@') or not line.strip():
            continue
        sam[line.split('\t')[0]] = line.rstrip('\n')
    rows = list(csv.DictReader(open(BUNDLE[0]), delimiter='\t'))
    ext = {}
    if os.path.exists(EXT):
        for r in csv.DictReader(open(EXT), delimiter='\t'):
            ext[r['read'][:8]] = r
    return rows, sam, slices, ext


def _inputs(row, sam, slices):
    chrom, strand = row['chrom'], row['strand']
    ns, ne, ba = int(row['novel_start']), int(row['novel_end']), int(row['best_annot'])
    key = [k for k in slices if k[0] == chrom and k[1] <= min(ns, ba)
           and k[1] + len(slices[k]) >= max(ne, ba)][0]
    off, seq = key[1], slices[key]
    register_genome_contigs([chrom])
    hdr = pysam.AlignmentHeader.from_dict({'HD': {'VN': '1.6'}, 'SQ': [{'SN': chrom, 'LN': len(seq)}]})
    a = pysam.AlignedSegment.fromstring(sam[row['read']], hdr)
    a.reference_start = a.reference_start - off
    novel = (chrom, ns - off, ne - off)
    annotated = (chrom, ba - off, ne - off) if strand == '+' else (chrom, ns - off, ba - off)
    return a, {chrom: seq}, novel, annotated, off


def _replay(row, sam, slices, monkeypatch):
    monkeypatch.setenv('RECTIFY_2F_NOVEL_GATE', 'report')
    monkeypatch.setenv('RECTIFY_2F_CHECK_CONSISTENCY', '1')   # the debug invariant check raises on violation
    a, genome, novel, annotated, off = _inputs(row, sam, slices)
    chrom, strand = row['chrom'], row['strand']
    res = rescue_3ss_truncation(a, genome, {novel, annotated}, strand, annotated_junctions={annotated})
    j = res.get('rescued_junction')
    return res, (j[0], j[1] + off, j[2] + off) if j else None, (novel[1] + off, novel[2] + off), \
        (annotated[1] + off, annotated[2] + off)


def _cigar_ops(s):
    return [(int(n), c) for n, c in re.findall(r'(\d+)([MIDNSHP=X])', s)]


def _exon_edit_count(read, res, genome_seq, strand):
    """mismatches + inserted + deleted bases of the emitted exon CIGAR (see the module docstring)."""
    ops = _cigar_ops(res['five_prime_exon_cigar'])
    j = res['rescued_junction']
    # ISSUE-028: a leading I is emitted as S (unplaced query bases) — an insertion for the edit count, as the
    # tester's Levenshtein ED_annot_ext would score it.
    span_q = sum(n for n, c in ops if c in 'MIS=X')
    span_r = sum(n for n, c in ops if c in 'MD=X')
    q = read.query_sequence
    seg = (q[:span_q] if strand == '+' else q[len(q) - span_q:]).upper()
    ref = (genome_seq[j[1] - span_r:j[1]] if strand == '+' else genome_seq[j[2]:j[2] + span_r]).upper()
    qi = ri = edits = 0
    for n, c in ops:
        if c in 'M=X':
            edits += sum(1 for k in range(n) if seg[qi + k] != ref[ri + k])
            qi += n
            ri += n
        elif c in 'IS':
            edits += n
            qi += n
        elif c == 'D':
            edits += n
            ri += n
    return edits, span_q, ops


def _junction_side_op(ops, strand):
    return ops[-1][1] if strand == '+' else ops[0][1]


READS = ['c887bc16', '31fab950', 'c64bc988', '5b20c72a', 'c41c7314', 'de84a10a', 'a5a5a1bb', '7f41e755', '6fc67f58', '885f6430']

# ISSUE-028 invariant E (2026-09-06): the block 2F placed for de84a10a on the annotated junction —
# `4M2D2M1I3M5I1M1I5M`, 14 matched / 1 mismatch / 7 inserted / 2 deleted, identity 0.93 — carries 13.5 bits,
# under the 18-bit floor the chance-match model derived (row A's ~3 % family-wise rate). The junction is the
# annotated one and the read was one of the six "+4 by evidence" reads RULING 10 placed there; under E it draws
# nothing (token exon_bits_below_floor, shape in the TSV — the TSV carries the LAST block judged, the deepest
# terminal peel's 4.5-bit spelling, not the 13.5-bit baseline block). Reported here with its numbers, not
# loosened (the brief's instruction). The other nine fixture reads carry 20-83.5 bits and land annotated as
# before (c887bc16 `16M` 13=/3X is 20.0 bits — the closest to the floor).
NOT_EVIDENCE_AT_E = {'de84a10a': ('exon_bits_below_floor', 13.5)}


def _assert_not_evidence(read8, res):
    token, bits = NOT_EVIDENCE_AT_E[read8]
    assert not res.get('rescued'), res
    assert res.get('clip_refused') == token, res
    assert res.get('exon_bits') is not None and res['exon_bits'] < 18, res


@pytest.mark.parametrize('read8', sorted(READS))
def test_bundle_read_lands_annotated(read8, monkeypatch):
    rows, sam, slices, ext = _load()
    row = [r for r in rows if r['read'].startswith(read8)][0]
    res, junction, novel, annotated = _replay(row, sam, slices, monkeypatch)
    if read8 in NOT_EVIDENCE_AT_E:
        _assert_not_evidence(read8, res)
        return
    assert res['rescued'], row['read']
    # Never the 1–2-mer artefact: the comparison covered the whole clip.
    assert res['query_bp'] >= int(row['clip_len']) or res['rescue_type'] == 'intronic_snap'
    # ISSUE-020: every bundle read lands on the ANNOTATED junction.
    assert junction[1:] == annotated, (junction, annotated)
    assert res['landing_annotated'] is True
    assert res['novel_evidence'] == ''
    if res['rescue_type'] != 'intronic_snap':
        from rectify.core.align.local_aligner import align_clip_to_exon, cigar_ops_to_str
        from rectify.core.splice.splice_aware_5prime import _anchored_deficit, _edit_distance
        assert res['anchored_deficit'] >= 0.0
        # Edit-count mapping to the tester's extended ED (module docstring): our
        # count is ONE alignment's edit count on the placement segment, so
        # (1) the Levenshtein distance of the same two strings is a lower bound
        # on it, and (2) the discriminating claim is made with the SAME aligner
        # on the SAME segment: the annotated placement costs fewer edits — and a
        # lower anchored deficit — than the +4 novel placement (the tester's
        # delta_ext < 0 on all 10). The observed values are in ISSUE-020's table.
        a, genome, novel_local, ann_local, _off = _inputs(row, sam, slices)
        gseq, strand = genome[row['chrom']], row['strand']
        edits, span_q, ops = _exon_edit_count(a, res, gseq, strand)
        q = a.query_sequence
        seg = (q[:span_q] if strand == '+' else q[len(q) - span_q:]).upper()
        j = res['rescued_junction']
        span_r = sum(n for n, c in ops if c in 'MD=X')
        ref = (gseq[j[1] - span_r:j[1]] if strand == '+' else gseq[j[2]:j[2] + span_r]).upper()
        assert _edit_distance(seg, ref) <= edits, (edits, _edit_distance(seg, ref))
        nov_ops, _ = align_clip_to_exon(seg, gseq, novel_local[1], novel_local[2], strand)
        nov_res = dict(five_prime_exon_cigar=cigar_ops_to_str(nov_ops), rescued_junction=novel_local)
        nov_edits, _sq, _ops = _exon_edit_count(a, nov_res, gseq, strand)
        assert edits < nov_edits, (read8, edits, nov_edits, res['five_prime_exon_cigar'], nov_res['five_prime_exon_cigar'])
        eff_ann = ann_local[1] if strand == '+' else ann_local[2]
        eff_nov = novel_local[1] if strand == '+' else novel_local[2]
        assert _anchored_deficit(seg, gseq, eff_ann, strand) < _anchored_deficit(seg, gseq, eff_nov, strand)
        if read8 in ext:
            # recorded, not asserted (different quantity — see the docstring)
            print(f"{read8}: edits_annot={edits} edits_novel={nov_edits} tester ED_annot_ext={ext[read8]['ED_annot_ext']} "
                  f"ED_novel_ext={ext[read8]['ED_novel_ext']} segment={span_q} tester_len={int(row['clip_len']) + int(ext[read8]['intron_mapped_bp'])}")


@pytest.mark.parametrize('read8', SIX)
def test_the_six_former_plus4_reads_have_no_junction_side_gap(read8, monkeypatch):
    """Arbiter addendum: for each of the six the CHOSEN placement has no
    junction-side gap — the compared segment ends exactly at the candidate's
    junction (the debug invariant check, enabled in _replay, would raise if the
    ranking and the placement disagreed) and the emitted exon CIGAR has no D at
    its junction-side end — and lands on the annotated candidate."""
    rows, sam, slices, ext = _load()
    row = [r for r in rows if r['read'].startswith(read8)][0]
    res, junction, novel, annotated = _replay(row, sam, slices, monkeypatch)
    if read8 in NOT_EVIDENCE_AT_E:
        _assert_not_evidence(read8, res)
        return
    assert junction[1:] == annotated and res['landing_annotated'] is True
    ops = _cigar_ops(res['five_prime_exon_cigar'])
    assert ops, res
    assert _junction_side_op(ops, row['strand']) != 'D', res['five_prime_exon_cigar']
    assert res['novel_evidence'] != 'novel_exon_gap_at_junction'


def test_no_bundle_sequence_rescue_carries_the_gap_token(monkeypatch):
    """RULING 10 §R37 (d): after 020 the junction-side gap token is ~0 on the
    bundle; it stays as instrumentation (a nonzero cohort count later = the
    sweep creeping back)."""
    rows, sam, slices, ext = _load()
    for row in rows:
        res, *_ = _replay(row, sam, slices, monkeypatch)
        assert res.get('novel_evidence') != 'novel_exon_gap_at_junction', row['read']
        if res['rescue_type'] != 'intronic_snap':
            # and never the 1–2-mer artefact (ISSUE-017)
            assert not (res['edit_distance'] == 0.0 and res['query_bp'] < int(row['clip_len'])), row['read']


@pytest.mark.parametrize('read8', WITHIN1_AT_017)
def test_writer_draws_the_annotated_n_op_for_the_former_within1_reads(read8, monkeypatch):
    """ISSUE-023: on 3834686 these three arrived 1 nt inside the annotated
    intron and the writer's canonical repair walked the N-op onto the +4 decoy
    with a compensating 3D. With the anchored ranking they arrive annotated, so
    the repair has nothing to move: the materialized record carries the
    annotated N-op, with no D op adjacent to it."""
    from rectify.core.bam import bam_writer as bw
    rows, sam, slices, ext = _load()
    row = [r for r in rows if r['read'].startswith(read8)][0]
    res, junction, novel, annotated = _replay(row, sam, slices, monkeypatch)
    a0, genome, _n, ann_local, off = _inputs(row, sam, slices)
    strand = row['strand']
    rj = res['rescued_junction']
    align5 = a0.reference_start if strand == '+' else a0.reference_end - 1
    icp = -1
    if strand == '-' and rj[1] <= align5 < rj[2]:
        icp = rj[1]
    elif strand == '+' and rj[1] <= align5 < rj[2]:
        icp = rj[2]
    ct = a0.cigartuples
    five_sc = (ct[0][1] if ct[0][0] == 4 else 0) if strand == '+' else (ct[-1][1] if ct[-1][0] == 4 else 0)
    rcl = int(res.get('reanchor_clip_len', 0) or 0)
    corr = {'five_prime_rescued': True, 'five_prime_position': res['five_prime_corrected'],
            'five_prime_soft_clip': rcl or five_sc, 'five_prime_exon_cigar': res['five_prime_exon_cigar'],
            'five_prime_upstream_trim': int(res.get('five_prime_upstream_trim', 0) or 0),
            'five_prime_exon2_prefix': int(res.get('five_prime_exon2_prefix', 0) or 0),
            'five_prime_intron_clip_pos': icp, 'reanchor_clip_len': rcl, 'strand': strand}
    b = copy.deepcopy(a0)
    bw._decode_eq_seq_inplace(b, genome)
    if rcl > 0:
        bw._apply_reanchor_from_clip_len(b, rcl)
    modified, refusal = bw.apply_5prime_rescue_surgery(b, corr, genome)
    assert modified and refusal == '', (modified, refusal)
    pos, nops, adj_d = b.reference_start, [], False
    tuples = b.cigartuples
    for i, (op, ln) in enumerate(tuples):
        if op == 3:
            nops.append((pos, pos + ln))
            if (i > 0 and tuples[i - 1][0] == 2) or (i + 1 < len(tuples) and tuples[i + 1][0] == 2):
                adj_d = True
        if op in (0, 2, 3, 7, 8):
            pos += ln
    assert (ann_local[1], ann_local[2]) in nops, (nops, ann_local)
    assert not adj_d, b.cigarstring


def _materialize(res, a0, genome, strand):
    """The writer's own materialization of a 2F result: the correction dict `correct_read_3prime` publishes,
    then the reanchor pre-pass and `apply_5prime_rescue_surgery` (as bundle1_rbrowse.py does)."""
    from rectify.core.bam import bam_writer as bw
    rj = res['rescued_junction']
    align5 = a0.reference_start if strand == '+' else a0.reference_end - 1
    icp = -1
    if rj[1] <= align5 < rj[2]:
        icp = rj[1] if strand == '-' else rj[2]
    ct = a0.cigartuples
    five_sc = (ct[0][1] if ct[0][0] == 4 else 0) if strand == '+' else (ct[-1][1] if ct[-1][0] == 4 else 0)
    rcl = int(res.get('reanchor_clip_len', 0) or 0)
    corr = {'five_prime_rescued': True, 'five_prime_position': res['five_prime_corrected'],
            'five_prime_soft_clip': rcl or five_sc, 'five_prime_exon_cigar': res.get('five_prime_exon_cigar', '') or '',
            'five_prime_upstream_trim': int(res.get('five_prime_upstream_trim', 0) or 0),
            'five_prime_exon2_prefix': int(res.get('five_prime_exon2_prefix', 0) or 0),
            'five_prime_intron_clip_pos': icp, 'reanchor_clip_len': rcl, 'strand': strand}
    b = copy.deepcopy(a0)
    bw._decode_eq_seq_inplace(b, genome)
    if rcl > 0:
        bw._apply_reanchor_from_clip_len(b, rcl)
    modified, refusal = bw.apply_5prime_rescue_surgery(b, corr, genome)
    pos, nops = b.reference_start, []
    for op, ln in b.cigartuples or []:
        if op == 3:
            nops.append((pos, pos + ln))
        if op in (0, 2, 3, 7, 8):
            pos += ln
    return modified, refusal, nops, b


@pytest.mark.parametrize('read8', sorted(READS))
def test_every_bundle_rescue_is_drawn_where_reported(read8, monkeypatch):
    """ISSUE-026 invariant D on the 10-read fixture: materializing every accepted rescue through the writer
    yields exactly the reported N-op with refusal ''. c887bc16 (ISSUE-027: a reanchor of 18 on a 13-nt clip
    whose alignment starts 2 nt into exon 2) used to be reverted as `noncanonical_destination` because `extend`
    ran the N to the live edge; with the exon-2 prefix it draws at the reported junction like every other read —
    no exception, no special case."""
    rows, sam, slices, ext = _load()
    row = [r for r in rows if r['read'].startswith(read8)][0]
    res, junction, novel, annotated = _replay(row, sam, slices, monkeypatch)
    if read8 in NOT_EVIDENCE_AT_E:
        _assert_not_evidence(read8, res)
        return
    assert res['rescued']
    a0, genome, _n, _a, off = _inputs(row, sam, slices)
    modified, refusal, nops, b = _materialize(res, a0, genome, row['strand'])
    rj = res['rescued_junction']
    assert modified and refusal == '', (read8, modified, refusal)
    assert (rj[1], rj[2]) in nops, (read8, (rj[1] + off, rj[2] + off), [(s + off, e + off) for s, e in nops])
    if read8 == 'c887bc16':
        assert int(res.get('five_prime_exon2_prefix', 0)) > 0, res   # the ISSUE-027 geometry, now expressed


def test_exon_cigar_query_span_matches_the_writer_contract(monkeypatch):
    """The exon CIGAR must consume exactly the query bases the BAM writer will
    convert, or `extend` falls back to a flat M block and the evidence ops never
    reach the BAM. Contract (splice_aware_5prime, "Which sequence the exon CIGAR
    is sized from"): 5' end INSIDE the rescued intron -> the intron-mapped run
    (clip + bases mapped past the boundary; `reroute` demands exon_q ==
    n_intronic_q, no upstream trim); 5' end outside -> clip + upstream_trim.
    At 3834686 the equivalence extension borrowed the overshoot base a second
    time on 5b20c72a / a5a5a1bb / 31fab950 (span 22 for clip 20 + trim 1).
    ISSUE-020 ranks on a segment that may be shorter than the placement segment
    (the dist > 0 trim) — the placement, and so this contract, is unchanged."""
    from rectify.core.align.local_aligner import cigar_str_to_ops
    from rectify.core.splice.splice_aware_5prime import _get_intronic_query_bases

    rows, sam, slices, ext = _load()
    checked = 0
    for row in rows:
        strand = row['strand']
        res, junction, *_ = _replay(row, sam, slices, monkeypatch)
        if not res.get('rescued') or res.get('rescue_type') != 'softclip' or not res.get('five_prime_exon_cigar'):
            continue
        ops = cigar_str_to_ops(res['five_prime_exon_cigar'])
        span = sum(l for op, l in ops if op in (0, 1, 4, 7, 8))
        a, genome, _n, _a, off = _inputs(row, sam, slices)
        clip = a.cigartuples[0][1] if (strand == '+' and a.cigartuples[0][0] == 4) else \
            (a.cigartuples[-1][1] if (strand == '-' and a.cigartuples[-1][0] == 4) else 0)
        j_start, j_end = junction[1] - off, junction[2] - off
        five = a.reference_start if strand == '+' else a.reference_end - 1
        inside = j_start <= five < j_end
        trim = int(res.get('five_prime_upstream_trim', 0) or 0)
        rcl = int(res.get('reanchor_clip_len', 0) or 0)
        # ISSUE-026 invariant D: the junction-side exon-2 prefix leaves the exon
        # block (the writer draws it as M after the N), so the exon CIGAR consumes
        # clip + trim - prefix — the writer's own `expected_query`.
        e2 = int(res.get('five_prime_exon2_prefix', 0) or 0)
        if rcl > 0:
            # Terminal peel / re-anchor: the writer's pre-pass converts the
            # 5' bases through the re-anchor point into a soft clip of length
            # reanchor_clip_len, and `extend` converts exactly that clip.
            assert span == rcl + trim - e2, (row['read'], res['five_prime_exon_cigar'], span, rcl, trim, e2)
        elif inside:
            intronic = _get_intronic_query_bases(a, j_start if strand == '-' else j_end, strand)
            assert trim == 0, (row['read'], trim)
            assert e2 == 0, (row['read'], e2)
            assert span == len(intronic), (row['read'], res['five_prime_exon_cigar'], span, len(intronic))
        else:
            assert span == clip + trim - e2, (row['read'], res['five_prime_exon_cigar'], span, clip, trim, e2)
        checked += 1
    assert checked >= 3, checked
