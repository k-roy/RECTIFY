"""Provenance regression: the junction-overhang table header must report the
concordance threshold actually ENFORCED, not the requested default.

`extract_concordant_junction_overhangs` silently clamps `n_concordant` down to
the number of aligners loaded (``n_concordant = max(2, len(dfs))``) when fewer
aligners are present. The header writer previously printed the unclamped argument,
producing the internally-impossible ``n_aligners: 3  n_concordant: 4`` seen in the
shipped homo_sapiens table. This pins the clamped value into the header.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from rectify.core.splice import calibrate_junction_overhang as cjo


def _write_aligner_tsv(path: Path) -> None:
    """One spliced read crossing the same junction, with a generous overhang so
    it survives the >=0 overhang filter. Same read/junction across all aligners
    makes it concordant."""
    pd.DataFrame({
        'read_id': ['readA'],
        'alignment_start': [1000],
        'alignment_end': [3000],
        'junctions': ['1500-2000'],
    }).to_csv(path, sep='\t', index=False)


def test_header_reports_clamped_n_concordant(tmp_path):
    specs = []
    for aligner in ('minimap2', 'uLTRA', 'deSALT'):  # 3 aligners
        p = tmp_path / f'{aligner}.tsv'
        _write_aligner_tsv(p)
        specs.append(f'{aligner}:{p}')

    out = tmp_path / 'overhang_table.tsv'
    # Request 4-concordance with only 3 aligners → must clamp to 3 in the header.
    cjo.calibrate(tsv_spec=specs, output=out, n_concordant=4)

    header = next(
        line for line in out.read_text().splitlines()
        if line.startswith('# n_aligners:')
    )
    assert 'n_aligners: 3' in header
    assert 'n_concordant: 3' in header, (
        f'header must report the clamped (enforced) concordance, got: {header!r}'
    )
    # Guard against the original bug specifically.
    assert 'n_concordant: 4' not in header
