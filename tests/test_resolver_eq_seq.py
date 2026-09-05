"""calmd '='-compressed SEQ must not survive a resolver rewrite (item A9).

A '=' byte means "identical to the reference AT THIS POSITION", so it is only
interpretable under the alignment it was written for. The resolver already
decodes '=' for SCORING (``_decoded_query``, and the ``inside_seq`` the caller
builds) but used to write the record back with the original '=' bytes under a
NEW alignment — so every base whose reference position moved (the ``k_inside``
bases re-assigned across a junction, and everything an arbiter shift/snap/split
displaces) decoded against the wrong reference downstream.

'='-compressed input is the production case, not an exotic one: ``samtools
calmd -e`` writes it, and the 748 fixture's minimap2 BAM carries '=' on
42,409/42,409 reads.

Pinned here: a rewritten read comes out with real letters equal to the decoded
input, its qualities intact; a read the resolver does not rewrite keeps its
SEQ byte-identical; and a read with no '=' is untouched either way.
"""

import random

import pysam
import pytest

from rectify.core.align.overhang_resolver import (
    ResolverConfig,
    ResolverStats,
    resolve_read,
    run_overhang_resolver,
)
from rectify.core.splice.overhang_informativeness import reset_counters
from rectify.core.splice.splice_site_index import SpliceSiteIndex

# ---------------------------------------------------------------------------
# Genome: the planning/644 T4c mis-anchored-acceptor geometry, which is the
# only one that exercises k_inside > 0 (aligned bases re-assigned ACROSS the
# junction — the bases whose '=' interpretation actually changes), plus a decoy
# acceptor 12 bp inside the intron tail for the arbiter's Case A shift.
# ---------------------------------------------------------------------------

GLEN = 4000
D, E = 1200, 1500
E_DEC = 1488


def _make_genome():
    rng = random.Random(90901)
    seq = [rng.choice('ACGT') for _ in range(GLEN)]

    def put(pos, s):
        seq[pos:pos + len(s)] = list(s)

    put(D, 'GT')
    put(E - 2, 'AG')
    put(E_DEC - 2, 'AG')
    # exon-1 end whose last 2 nt are AG, duplicated as the intron tail — the
    # local homology that makes an aligner anchor 10 bp into the intron.
    put(D - 10, ''.join(rng.choice('ACGT') for _ in range(8)) + 'AG')
    seq[E - 10:E] = seq[D - 10:D]
    return ''.join(seq)


GENOME_SEQ = _make_genome()
GENOME = {'chrI': GENOME_SEQ}


@pytest.fixture(scope='module')
def index():
    return SpliceSiteIndex.build(GENOME)


@pytest.fixture(autouse=True)
def _fresh_counters():
    reset_counters()
    yield
    reset_counters()


def _header():
    return pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.6', 'SO': 'queryname'},
        'SQ': [{'SN': 'chrI', 'LN': GLEN}],
    })


def eq_encode(query, cigartuples, reference_start, genome_seq=GENOME_SEQ):
    """What ``samtools calmd -e`` writes: an M base identical to the reference
    becomes '='. Soft clips and insertions are never '='-encoded."""
    out = []
    qi, ref = 0, reference_start
    for op, ln in cigartuples:
        if op in (0, 7, 8):
            for i in range(ln):
                b = query[qi + i]
                out.append('=' if b.upper() == genome_seq[ref + i].upper() else b)
            qi += ln
            ref += ln
        elif op in (1, 4):
            out.append(query[qi:qi + ln])
            qi += ln
        elif op in (2, 3):
            ref += ln
    return ''.join(out)


def _read(name, query, cigar, ref_start, quals=True, header=None):
    r = pysam.AlignedSegment(header or _header())
    r.query_name = name
    r.query_sequence = query
    r.flag = 0
    r.reference_id = 0
    r.reference_start = ref_start
    r.mapping_quality = 60
    r.cigartuples = cigar
    if quals:
        # a recognisable, non-constant pattern so a silent reset is visible
        r.query_qualities = pysam.qualitystring_to_array(
            ''.join(chr(33 + 10 + (i % 30)) for i in range(len(query))))
    return r


def _resolve(read, index, **cfg_kw):
    stats = ResolverStats()
    changed = resolve_read(read, GENOME, index,
                           ResolverConfig(alpha=0.01, max_intron=5000, **cfg_kw),
                           stats)
    qlen = sum(ln for op, ln in read.cigartuples if op in (0, 1, 4, 7, 8))
    assert qlen == len(read.query_sequence), (
        f'query length broken by the rewrite: {read.cigartuples}')
    return changed, stats


# ---------------------------------------------------------------------------
# The k_inside case — bases genuinely re-assigned across the junction.
# ---------------------------------------------------------------------------

def _inside_edge_read(encoded, quals=True):
    """A clip whose true acceptor sits 10 bp INSIDE the aligned block."""
    true_query = GENOME_SEQ[D - 30:D - 10] + GENOME_SEQ[E - 10:E] + \
        GENOME_SEQ[E:E + 50]
    cigar = [(4, 20), (0, 60)]
    q = eq_encode(true_query, cigar, E - 10) if encoded else true_query
    return _read('inside_edge', q, cigar, E - 10, quals=quals), true_query


class TestKInsideRewrite:
    def test_eq_input_comes_out_decoded(self, index):
        r, true_query = _inside_edge_read(encoded=True)
        assert '=' in r.query_sequence, 'the fixture is not actually encoded'
        changed, stats = _resolve(r, index)
        assert changed, stats.as_dict()
        assert stats.extra.get('resolved_inside_edge') == 1, (
            'this fixture must exercise k_inside > 0 — otherwise it does not '
            'test the bases whose reference position actually moved')
        assert '=' not in r.query_sequence
        assert r.query_sequence == true_query
        assert stats.extra.get('seq_eq_decoded') == 1

    def test_qualities_survive_the_seq_assignment(self, index):
        # pysam clears query_qualities on any query_sequence assignment
        r, _ = _inside_edge_read(encoded=True)
        before = list(r.query_qualities)
        assert _resolve(r, index)[0]
        assert r.query_qualities is not None, 'qualities were dropped'
        assert list(r.query_qualities) == before

    def test_read_without_qualities_is_fine(self, index):
        r, true_query = _inside_edge_read(encoded=True, quals=False)
        assert _resolve(r, index)[0]
        assert r.query_sequence == true_query

    def test_plain_input_is_byte_identical_to_the_decoded_one(self, index):
        # the '=' path must not change the ANSWER, only the spelling
        enc, _ = _inside_edge_read(encoded=True)
        plain, _ = _inside_edge_read(encoded=False)
        assert _resolve(enc, SpliceSiteIndex.build(GENOME))[0]
        assert _resolve(plain, SpliceSiteIndex.build(GENOME))[0]
        assert enc.query_sequence == plain.query_sequence
        assert enc.cigartuples == plain.cigartuples
        assert enc.reference_start == plain.reference_start


# ---------------------------------------------------------------------------
# The arbiter's rewrites move bases too.
# ---------------------------------------------------------------------------

class TestArbiterRewrite:
    def test_case_a_shift_decodes_the_seq(self, index):
        # spliced at (D, E) but the CIGAR asserts the decoy acceptor: exon-2's
        # query is misaligned over the intron tail, so a boundary shift moves
        # query bases across the junction by 12 bp.
        true_query = GENOME_SEQ[D - 60:D] + GENOME_SEQ[E:E + 60]
        cigar = [(0, 60), (3, E_DEC - D), (0, 60)]
        q = eq_encode(true_query, cigar, D - 60)
        assert '=' in q
        r = _read('shift', q, cigar, D - 60)
        changed, stats = _resolve(r, index)
        assert changed and stats.extra.get('arb_shifted') == 1, stats.as_dict()
        assert '=' not in r.query_sequence
        assert r.query_sequence == true_query
        assert stats.extra.get('seq_eq_decoded') == 1


# ---------------------------------------------------------------------------
# Passthrough must not be touched.
# ---------------------------------------------------------------------------

class TestPassthrough:
    def test_unrewritten_eq_read_keeps_its_bytes(self, index):
        # a clean linear read: nothing to resolve, nothing to arbitrate
        true_query = GENOME_SEQ[100:220]
        cigar = [(0, 120)]
        q = eq_encode(true_query, cigar, 100)
        assert q == '=' * 120
        r = _read('clean', q, cigar, 100)
        quals = list(r.query_qualities)
        changed, stats = _resolve(r, index)
        assert not changed
        assert r.query_sequence == q, 'passthrough SEQ was rewritten'
        assert list(r.query_qualities) == quals
        assert 'seq_eq_decoded' not in stats.extra

    def test_refused_clip_keeps_its_bytes(self, index):
        # a poly(A) clip is refused by the information gate — still passthrough
        true_query = 'A' * 30 + GENOME_SEQ[E:E + 60]
        cigar = [(4, 30), (0, 60)]
        q = eq_encode(true_query, cigar, E)
        r = _read('polya', q, cigar, E)
        changed, stats = _resolve(r, index)
        assert not changed
        assert r.query_sequence == q
        assert 'seq_eq_decoded' not in stats.extra

    def test_read_without_eq_is_never_touched(self, index):
        r, true_query = _inside_edge_read(encoded=False)
        changed, stats = _resolve(r, index)
        assert changed
        assert r.query_sequence == true_query
        assert 'seq_eq_decoded' not in stats.extra


# ---------------------------------------------------------------------------
# Through the real driver, so htslib gets a say.
# ---------------------------------------------------------------------------

class TestDriverRoundTrip:
    def test_eq_bam_round_trips(self, tmp_path):
        fasta = tmp_path / 'genome.fa'
        with open(fasta, 'w') as fh:
            fh.write('>chrI\n' + GENOME_SEQ + '\n')
        header = _header()
        enc, true_query = _inside_edge_read(encoded=True)
        clean_q = eq_encode(GENOME_SEQ[100:220], [(0, 120)], 100)
        reads = [
            _read('r1_eq_resolved', enc.query_sequence, enc.cigartuples,
                  enc.reference_start, header=header),
            _read('r2_eq_passthrough', clean_q, [(0, 120)], 100, header=header),
        ]
        in_bam = tmp_path / 'in.bam'
        with pysam.AlignmentFile(in_bam, 'wb', header=header) as fh:
            for r in reads:
                fh.write(r)
        out_bam = tmp_path / 'out.bam'
        run_overhang_resolver(str(in_bam), str(fasta), str(out_bam))

        with pysam.AlignmentFile(out_bam, 'rb') as fh:
            got = {r.query_name: r for r in fh}
        assert got['r1_eq_resolved'].query_sequence == true_query
        assert '=' not in got['r1_eq_resolved'].query_sequence
        assert got['r1_eq_resolved'].query_qualities is not None
        # untouched record keeps the compression it arrived with
        assert got['r2_eq_passthrough'].query_sequence == clean_q
