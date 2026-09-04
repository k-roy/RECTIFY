#!/usr/bin/env python3
"""833 G-8 — every COMPASS short-read arm must receive rectify's DERIVED
--max-intron, not the aligner's own default.

Before this fix `align_command` and `run_multi_aligner` called `run_star`,
`run_hisat2`, `run_magicblast` and `run_gsnap` with no intron cap, so each arm
used its own default — STAR ~589 kb (winBinNbits-derived), HISAT2 200,000,
Magic-BLAST 500,000, GSNAP 200,000 — while rectify's own derived cap on
S. cerevisiae is 5,000 bp.  The arms therefore emitted junctions rectify itself
classifies as impossible, and those candidates entered the consensus pool:
17 of the 25 contig-swap records forensically examined in planning/861 §S8e had
a HISAT2 candidate carrying the same recurring 5,584 nt N-op on the borrowed
contig.

The argv is the ONLY auditable surface for this: `multi_aligner` logs
`' '.join(cmd[:6])`, the first six tokens, so no run log can show these flags
(planning/861 §S10, G-22).  These tests therefore assert on the constructed
command, not on a log.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest

from rectify.core.align import multi_aligner as ma


MAXI = 5000


def _pairs(cmd):
    """argv as {flag: value} for '--flag value' style, plus bare tokens."""
    out = {}
    for i, tok in enumerate(cmd):
        if str(tok).startswith('-') and i + 1 < len(cmd) and not str(cmd[i + 1]).startswith('-'):
            out[str(tok)] = str(cmd[i + 1])
    return out


def test_star_se_gets_alignintronmax():
    cmd = ma._build_star_cmd('r1.fq.gz', None, '/idx', '/out.', 4,
                             read_length=75, max_intron=MAXI)
    assert _pairs(cmd)['--alignIntronMax'] == str(MAXI)
    # SE has no mate gap to cap.
    assert '--alignMatesGapMax' not in cmd


def test_star_pe_also_gets_alignmatesgapmax():
    cmd = ma._build_star_cmd('r1.fq.gz', 'r2.fq.gz', '/idx', '/out.', 4,
                             read_length=75, max_intron=MAXI)
    p = _pairs(cmd)
    assert p['--alignIntronMax'] == str(MAXI)
    assert p['--alignMatesGapMax'] == str(MAXI)


def test_star_omits_the_cap_when_unset():
    """0 = leave STAR's own default alone (back-compat for direct callers)."""
    cmd = ma._build_star_cmd('r1.fq.gz', 'r2.fq.gz', '/idx', '/out.', 4)
    assert '--alignIntronMax' not in cmd


def test_hisat2_gets_max_intronlen():
    cmd = ma._build_hisat2_cmd('r1.fq.gz', 'r2.fq.gz', '/idx', '/ss.txt',
                               '/out.sam', 4, 20, MAXI, '/novel.txt', '/sum.txt')
    p = _pairs(cmd)
    assert p['--max-intronlen'] == str(MAXI)
    assert p['--max-intronlen'] != '200000', "hisat2's own default leaked through"


def test_hisat2_wrapper_default_is_rectifys_not_hisat2s():
    """The wrapper default must not be hisat2's 200,000 — a caller that forgets
    to pass max_intron is exactly how G-8 happened."""
    import inspect
    sig = inspect.signature(ma.run_hisat2)
    assert sig.parameters['max_intron'].default == ma.DEFAULT_MAX_INTRON
    assert sig.parameters['max_intron'].default != 200000


def test_magicblast_gets_max_intron_length():
    cmd = ma._build_magicblast_cmd('r1.fq', 'r2.fq', '/db', '/out.sam',
                                   threads=12, max_intron=MAXI)
    assert _pairs(cmd)['-max_intron_length'] == str(MAXI)
    # -out must remain the last pair; the wrapper's callers rely on the SAM path.
    assert cmd[-2] == '-out'


def test_gsnap_gets_localsplicedist():
    cmd = ma._build_gsnap_cmd('r1.fq.gz', 'r2.fq.gz', '/gdir', 'gv', '/out.sam',
                              4, max_intron=MAXI)
    assert f'--localsplicedist={MAXI}' in cmd
    # --pairmax-rna bounds fragment length + intron, not the intron: leave it.
    assert not any(str(t).startswith('--pairmax-rna') for t in cmd)


@pytest.mark.parametrize('aligner,kwarg', [
    ('STAR_default', 'max_intron'),
    ('STAR_noncanonical', 'max_intron'),
    ('HISAT2_default', 'max_intron'),
    ('HISAT2_noncanonical', 'max_intron'),
    ('magicblast', 'max_intron'),
    ('gsnap', 'max_intron'),
])
def test_dispatch_forwards_max_intron_to_every_compass_arm(
        aligner, kwarg, tmp_path, monkeypatch):
    """The builders being right is not enough — G-8 was a CALL SITE omission."""
    seen = {}

    def _rec(name):
        def _f(**kw):
            seen[name] = kw
            return str(tmp_path / f'{name}.bam')
        return _f

    for fn in ('run_star', 'run_hisat2', 'run_magicblast', 'run_gsnap'):
        monkeypatch.setattr(ma, fn, _rec(fn), raising=True)

    ma.run_multi_aligner(
        reads_path=str(tmp_path / 'r1.fq.gz'),
        genome_path=str(tmp_path / 'g.fa'),
        output_dir=str(tmp_path / 'out'),
        sample_name='s',
        aligners=[aligner],
        reads2_path=str(tmp_path / 'r2.fq.gz'),
        max_intron=MAXI,
    )
    assert seen, f"{aligner} was skipped — the assertion below would be vacuous"
    kw = next(iter(seen.values()))
    assert kw.get(kwarg) == MAXI, (
        f"{aligner} was dispatched without {kwarg}; it would fall back to the "
        f"aligner's own default (got {kw.get(kwarg)!r})"
    )


def test_derived_cap_on_the_bundled_yeast_annotation_is_5000():
    """Guard the premise: the derived value must actually be the tight one."""
    from pathlib import Path
    import rectify
    gff = (Path(rectify.__file__).parent / 'data'
           / 'saccharomyces_cerevisiae_R64-5-1_20240529.gff')
    if not gff.exists():
        pytest.skip('bundled annotation not present')
    assert ma.derive_max_intron(str(gff)) == 5000
