# HANDOFF — Fix NaN-junctions bug in `merge_corrected_tsvs`

**Date:** 2026-05-20
**Branch:** `drs-validation-rebuild`
**Discovered by:** Opus 4.7 during Commit A regression sweep
**Estimated effort:** 5-10 minutes (2-line fix + verify)
**Owner:** unassigned — any agent or Kevin

---

## Symptom

`pytest -m "not slow"` on `drs-validation-rebuild` HEAD `cb2fe6c` reports:

```
FAILED tests/test_corrected_consensus_tiebreaker.py::test_paralog_tiebreaker_picks_multi_aligner_consensus
FAILED tests/test_corrected_consensus_tiebreaker.py::test_majority_consensus_picks_chrXIV_even_when_outlier_has_wider_span
= 2 failed, 900 passed, 25 skipped, 4 deselected, 1 xfailed, 1 warning in 35.75s =
```

Both fail with:

```
AttributeError: 'float' object has no attribute 'split'
rectify/core/consensus/corrected_consensus.py:830
```

## Root cause

`rectify/core/consensus/corrected_consensus.py:826-830`:

```python
def _eff_key(_row):
    _juncs = _row.get('junctions', '') or ''
    # Normalize: sort the donor-acceptor tuples so semicolon-order
    # differences across aligners don't shadow real equivalence.
    _parts = tuple(sorted(p for p in _juncs.split(';') if p))
```

`_row` is a pandas Series. When a per-aligner TSV doesn't include the `junctions` column (or includes it but with a missing value for some rows), `_row.get('junctions', '')` returns **`float('nan')`** — not `''`. Then `nan or ''` evaluates to `nan` because **NaN is truthy in Python's `or` semantics** (`bool(float('nan'))` is `True`). The subsequent `.split(';')` then fails on a float.

The two failing tests build minimal per-aligner TSVs via `_write_tsv(...)` + `_make_row(...)`; the helper produces rows that lack the `junctions` column entirely. Pandas fills the missing column with NaN when merging across aligners. Production runs hit this less often because real per-aligner TSVs from `rectify correct` always include the `junctions` column, but the path is still latently unsafe.

## Fix

Replace the `or ''` idiom with an explicit isna check.

`rectify/core/consensus/corrected_consensus.py:826-830` (current):

```python
def _eff_key(_row):
    _juncs = _row.get('junctions', '') or ''
    _parts = tuple(sorted(p for p in _juncs.split(';') if p))
```

Proposed:

```python
def _eff_key(_row):
    _juncs = _row.get('junctions', '')
    if not isinstance(_juncs, str):  # catches NaN (float) and None
        _juncs = ''
    _parts = tuple(sorted(p for p in _juncs.split(';') if p))
```

(Equivalent alternative: `if pd.isna(_juncs): _juncs = ''`. The `isinstance` form avoids importing pandas in this hot per-row callback.)

## Acceptance

- [ ] Both `test_corrected_consensus_tiebreaker` tests pass.
- [ ] No other tests regress (`pytest -m "not slow"` shows 902 passed in total, assuming Commit A is already landed; otherwise 900 + 2 newly-fixed).
- [ ] Commit message: `fix(consensus): handle NaN junctions column in _eff_key merge tiebreaker`.

## Why this is in its own commit

Out of scope for Commit A (which only adds infrastructure under `rectify/core/bam/` + `rectify/core/splice/junction_scoring.py`). `corrected_consensus.py` was untouched by Commit A; this fix is for a pre-existing latent bug surfaced by tests written separately. Keep commits surgical per `CLAUDE.md`.

## Dispatch options

- **Direct fix by Kevin/Opus:** 1 edit + verify, ~5 min.
- **Sonnet subagent:** trivial; a paste-the-fix-and-run-pytest task. Briefing would be ~30 lines.
- **Bundle with another small follow-up commit** if there are other 1-2 line bugs queued.
