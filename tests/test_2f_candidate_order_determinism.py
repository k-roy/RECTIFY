"""ISSUE-019 — candidate order is a documented sort, never set order.

(a) PRIOR: among candidates that contain the read's 5' end, the ANNOTATED
intron is preferred by the Case-4 intronic snap; novel introns rank by the
smallest snap distance; and on an equal-ED tie in the sequence rescue the
annotated candidate wins before any geometry tiebreaker.
(b) DETERMINISM: candidates carry a chromosome string and Python randomizes
str hashing per process, so any set-order dependence could flip between runs
and the byte-identity tier would see phantom churn. The same Case-4 tie is
run under two PYTHONHASHSEED values in fresh interpreters and must give
byte-identical output.
"""
import json
import os
import subprocess
import sys
import textwrap

import pysam

from rectify.core.splice.splice_aware_5prime import rescue_3ss_truncation

EXON1_TAIL = 'ACGTTGCATGCAGTCCATG'
_GENOME_SEQ = (
    ('T' * (40 - len(EXON1_TAIL) - 1)) + 'A' + EXON1_TAIL   # exon1  [0, 40)
    + 'GTAAGT' + 'N' * 92 + 'AG'                            # intron [40, 140)
    + 'C' * 100                                             # exon2  [140, 240)
)
GENOME = {'chrT': _GENOME_SEQ}
ANNOTATED = ('chrT', 40, 140)
NOVEL_NEAR = ('chrT', 44, 140)      # a pool intron 4 nt in (the GTRAGT decoy), also contains the 5' end
NOVEL_FAR = ('chrT', 20, 140)       # another containing intron, farther donor


def _read(start, seq, name='r'):
    hdr = pysam.AlignmentHeader.from_dict({'HD': {'VN': '1.6'},
                                           'SQ': [{'SN': 'chrT', 'LN': 3_000_000}]})
    r = pysam.AlignedSegment(hdr)
    r.query_name = name
    r.reference_name = 'chrT'
    r.reference_start = start
    r.cigartuples = [(0, len(seq))]
    r.mapping_quality = 60
    r.query_sequence = seq
    return r


def _case4_tie():
    """Five exon-1 bases mapped inside every candidate intron at the shared
    acceptor (5' end at 135: 95 nt past the annotated donor, 91 past the
    decoy's, 115 past the far one): a 5-mer is below the informative floor, so
    the Case-4 snap decides among three containing introns. ISSUE-028
    (2026-09-06): the read used to start at 46 and carry the whole N-filled
    intron, so the snap re-placed a 94-base intron-mapped run onto exon 1 — a
    block with no identity to exon 1 that the evidence floor now refuses; five
    exon-1 bases keep the fixture's question (which containing intron wins)
    with a block below the informative floor, which rests on the structural
    prior as before."""
    read = _read(135, _GENOME_SEQ[35:40] + 'C' * 60, 'tie')
    return rescue_3ss_truncation(read, GENOME, {NOVEL_FAR, NOVEL_NEAR, ANNOTATED}, '+',
                                 annotated_junctions={ANNOTATED})


def test_case4_snap_prefers_the_annotated_intron():
    res = _case4_tie()
    assert res['rescued'] and res['rescue_type'] == 'intronic_snap'
    assert res['rescued_junction'] == ANNOTATED
    assert res['landing_annotated'] is True


def test_case4_novel_only_picks_the_smallest_snap_distance():
    # nine exon-1 bases (a 5-mer favors the far intron's own sequence here; nine favor exon 1 for both candidates)
    read = _read(131, _GENOME_SEQ[31:40] + 'C' * 60, 'novel_only')
    res = rescue_3ss_truncation(read, GENOME, {NOVEL_FAR, NOVEL_NEAR}, '+',
                                annotated_junctions=set())
    assert res['rescued'] and res['rescue_type'] == 'intronic_snap'
    # NOVEL_NEAR's donor is 4 nt closer to the 5' end than NOVEL_FAR's -> smaller snap distance
    assert res['rescued_junction'] == NOVEL_NEAR


_CHILD = textwrap.dedent('''
    import json, sys
    sys.path.insert(0, %r)
    import test_2f_candidate_order_determinism as t
    res = t._case4_tie()
    print(json.dumps({k: (list(v) if isinstance(v, tuple) else v) for k, v in res.items()
                      if k in ("rescued", "rescue_type", "rescued_junction", "five_prime_corrected",
                               "five_prime_exon_cigar", "landing_annotated", "novel_evidence")},
                     sort_keys=True))
''')


def test_case4_tie_is_identical_under_two_hash_seeds():
    here = os.path.dirname(os.path.abspath(__file__))
    outs = []
    for seed in ('1', '2'):
        env = dict(os.environ, PYTHONHASHSEED=seed)
        env.pop('RECTIFY_2F_NOVEL_GATE', None)
        p = subprocess.run([sys.executable, '-c', _CHILD % here], env=env,
                           capture_output=True, text=True, timeout=300)
        assert p.returncode == 0, p.stderr[-2000:]
        outs.append(p.stdout.strip().splitlines()[-1])
    assert outs[0] == outs[1], outs
    assert json.loads(outs[0])['rescued_junction'] == list(ANNOTATED)
