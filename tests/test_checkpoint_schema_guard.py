"""A resume against a checkpoint written by a different corrected-TSV schema must fail
loudly, not produce a current header over foreign rows (ledger D3 risk, 2026-09-05)."""
import json
from pathlib import Path

import pytest

from rectify.core.bam import parallel as p
from rectify.core.bam.output import CORRECTION_TSV_HEADER


def _write_region(tmp_path: Path, ncols: int) -> Path:
    ck = tmp_path / 'ckpt'
    ck.mkdir()
    row = ['x'] * ncols
    row[0] = 'read_1'
    (ck / 'region_0000.tsv').write_text('\t'.join(row) + '\n')
    (ck / 'region_0000.stats.json').write_text(json.dumps({}))
    (ck / 'region_0000.done').write_text('done\n')
    return ck


def test_old_schema_checkpoint_is_refused_loudly(tmp_path):
    ck = _write_region(tmp_path, len(CORRECTION_TSV_HEADER) - 5)
    out = tmp_path / 'out.tsv'
    with pytest.raises(RuntimeError) as exc:
        p._rebuild_output_from_region_files(ck, str(out), 1, p.ProcessingStats())
    msg = str(exc.value)
    assert 'different schema' in msg and str(ck) in msg
    assert f'{len(CORRECTION_TSV_HEADER)} columns' in msg


def test_blank_lines_are_ignored_not_fatal(tmp_path):
    ck = _write_region(tmp_path, len(CORRECTION_TSV_HEADER))
    with open(ck / 'region_0000.tsv', 'a') as fh:
        fh.write('\n')
    out = tmp_path / 'out.tsv'
    try:
        p._rebuild_output_from_region_files(ck, str(out), 1, p.ProcessingStats())
    except RuntimeError as e:  # a schema refusal here would be the bug
        pytest.fail(f'blank line treated as a schema mismatch: {e}')
    except Exception:
        # other parse failures on the dummy row are outside this test's contract
        pass
