"""The silent-discard guard: `correct` must not throw away CIGAR-level work quietly.

Background (planning/691). ``rectify correct`` makes two edits that live in the
CIGAR and nowhere else — the 5' rescue N-op and the 3' hard clip. They reach a
BAM only when an output BAM was requested. Without a BAM flag the command still
exits 0 and its stats block still reports the corrections, so "no corrections
were needed" and "the corrections were computed and discarded" are
indistinguishable. A 51-sample DRS cohort was corrected that way, discarding a
121,215-read 5' rescue per sample; the browser then drew raw soft clips under a
page labelled "rectify-corrected" and it read as a broken 5' rescue.

TEST-DESIGN NOTE — this guard is in the class that fails ONLY when it fires, so
a happy-path test proves nothing. Two rules are enforced below:

  * every silent-case assertion also asserts the guard's PRECONDITION was
    reached (corrections were actually applied), so the test cannot pass
    vacuously by correcting nothing; and
  * the integration test checks the REMEDIATION the warning tells the user to
    perform — re-running with ``--write-corrected-bam`` — actually produces the
    N-op, because a warning that recommends a no-op is worse than no warning.
"""

import subprocess
import sys
from pathlib import Path

import pysam
import pytest

from rectify.core.commands.correct_command import build_no_bam_output_warning
from rectify.core.bam.processing_stats import ProcessingStats

import rectify

_DATA = Path(rectify.__file__).parent / 'data' / 'genomes' / 'saccharomyces_cerevisiae'
_FSA = _DATA / 'S288C_reference_sequence_R64-5-1_20240529.fsa'

# RPL22B / YFL034C-A, chrVI, minus strand. GFF 1-based intron 64600-64920
# => 0-based half-open (64599, 64920), length 321.
_I_S, _I_E = 64599, 64920
_BODY = 250


# ---------------------------------------------------------------------------
# Unit: the decision itself
# ---------------------------------------------------------------------------

def _stats(five_prime=0, shifts=0):
    s = ProcessingStats()
    s.ends_five_prime_rescued = five_prime
    s.total_position_shifts = shifts
    return s


def test_warns_when_five_prime_rescues_are_discarded():
    lines = build_no_bam_output_warning({}, _stats(five_prime=121_215))
    assert lines is not None
    assert 'NO CORRECTED BAM WAS WRITTEN' in lines[0]
    # The count must be in the message — a bare "something was discarded" does
    # not tell the operator whether it mattered.
    assert '121,215' in ' '.join(lines)
    assert '--write-corrected-bam' in ' '.join(lines)


def test_warns_when_only_three_prime_shifts_are_discarded():
    lines = build_no_bam_output_warning({}, _stats(shifts=42))
    assert lines is not None
    assert '42' in ' '.join(lines)


@pytest.mark.parametrize('key', ['corrected_bam', 'softclipped_bam', 'output_bam'])
def test_silent_when_any_bam_output_requested(key):
    """Anti-vacuity: the same stats that DO warn must go silent, and only
    because a BAM was requested."""
    stats = _stats(five_prime=121_215)
    assert build_no_bam_output_warning({}, stats) is not None      # precondition reached
    assert build_no_bam_output_warning({key: '/tmp/out.bam'}, stats) is None


def test_silent_when_nothing_was_corrected():
    assert build_no_bam_output_warning({}, _stats()) is None


def test_silent_when_stats_missing():
    """A crashed/absent stats object must not manufacture a warning."""
    assert build_no_bam_output_warning({}, None) is None


# ---------------------------------------------------------------------------
# Integration: the wiring, and the remediation the warning recommends
# ---------------------------------------------------------------------------

def _load_chrom(path, want='chrVI'):
    import re
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith('>'):
            if name:
                seqs[name] = ''.join(buf)
            hdr = line[1:].strip()
            m = re.search(r'chromosome=([IVXA-Za-z]+)', hdr)
            name = 'chr' + m.group(1) if m else hdr.split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = ''.join(buf)
    return seqs


@pytest.fixture(scope='module')
def rpl22b_bam(tmp_path_factory):
    """Reads truncated at the RPL22B 3'SS, 5' clip = the upstream exon.

    Minus strand, so in reference-forward BAM coordinates the transcript's 5'
    end is the TRAILING soft clip.
    """
    genome = _load_chrom(_FSA)
    g = genome['chrVI']
    contigs = [(c, len(s)) for c, s in genome.items()]
    out = tmp_path_factory.mktemp('rpl22b') / 'rpl22b.bam'

    header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
              'SQ': [{'SN': c, 'LN': n} for c, n in contigs]}
    recs = []
    for clip in (1, 8, 57):          # 1 bp must rescue, not only comfortable clips
        a = pysam.AlignedSegment()
        a.query_name = f'rpl22b_clip{clip}'
        a.query_sequence = g[_I_S - _BODY:_I_S] + g[_I_E:_I_E + clip]
        a.flag = 16
        a.reference_id = [c for c, _ in contigs].index('chrVI')
        a.reference_start = _I_S - _BODY
        a.mapping_quality = 60
        a.cigartuples = [(0, _BODY), (4, clip)]
        a.query_qualities = pysam.qualitystring_to_array('I' * len(a.query_sequence))
        recs.append(a)
    with pysam.AlignmentFile(str(out), 'wb', header=header) as f:
        for a in recs:
            f.write(a)
    pysam.index(str(out))
    return out


def _run_correct(bam, outdir, *extra):
    cmd = [sys.executable, '-m', 'rectify', 'correct', str(bam), '--Scer',
           '-o', str(outdir / 'c.tsv'), *extra]
    return subprocess.run(cmd, capture_output=True, text=True, timeout=900)


def test_cli_warns_then_the_remediation_works(rpl22b_bam, tmp_path):
    # Deliberately NOT marked `slow` (~25 s). The defect this guards against
    # survived precisely because nothing exercised the path by default; a
    # regression test that only runs on request would reproduce that.
    # 1. No BAM flag: the guard must fire, and must fire because rescues really
    #    happened (not because the guard is trigger-happy).
    d1 = tmp_path / 'no_bam'
    d1.mkdir()
    r1 = _run_correct(rpl22b_bam, d1)
    assert r1.returncode == 0, r1.stderr[-3000:]
    combined = r1.stdout + r1.stderr
    assert "5' soft-clip rescued:" in combined
    assert 'NO CORRECTED BAM WAS WRITTEN' in r1.stderr
    # precondition: a non-zero rescue count reached the message
    assert "0 read(s) had a 5' end rescued" not in r1.stderr
    # and nothing BAM-shaped was produced
    assert not list(d1.glob('*.bam'))

    # 2. The remediation the warning recommends: same input, add the flag.
    #    Must go silent AND must actually materialise the N-op.
    d2 = tmp_path / 'with_bam'
    d2.mkdir()
    out_bam = d2 / 'c.bam'
    r2 = _run_correct(rpl22b_bam, d2, '--write-corrected-bam', str(out_bam))
    assert r2.returncode == 0, r2.stderr[-3000:]
    assert 'NO CORRECTED BAM WAS WRITTEN' not in r2.stderr
    assert out_bam.exists()

    with pysam.AlignmentFile(str(out_bam)) as f:
        cigars = {a.query_name: a.cigarstring for a in f}
    assert len(cigars) == 3
    for name, cig in cigars.items():
        clip = int(name.rsplit('clip', 1)[1])
        # 321 = the RPL22B intron; the clip becomes an exon block beyond it.
        assert f'321N{clip}M' in cig, f'{name}: {cig}'
