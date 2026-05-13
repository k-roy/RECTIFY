# HANDOFF — QuantSeq REV walkback integration

**Status:** core implemented; QuantSeq wrapper + CLI wiring pending
**Filed:** 2026-05-13
**Prior context:** `.github/ISSUE_walkback_integration.md` (NEW-075)
**For:** the agent picking up the QuantSeq REV piece

---

## Background in 30 seconds

`rectify correct --short-read --dT-primed-cDNA` produces output where ~92% of QuantSeq REV reads emit `correction_applied="none"`. The cause: Module 2E consults only the reference genome and never compares read bases to genome bases, so it can't detect when QuantSeq's internal priming has dragged the alignment a few bases into a downstream poly-A stretch. Existing workaround: a standalone `11_polya_walkback_recompute.py` script that produces `_corrected_3ends_wb.tsv` files — but it emits **BAM strand instead of gene/RNA strand**, which silently corrupts antisense classification at convergent loci (see the ZWF1 antisense vs YNL242W misassignment that came out of the IGV validation pass on Han 2023).

---

## What was done in this pass (already on GitHub)

### New module: `rectify/core/correct/walkback.py`

Protocol-agnostic core. Public surface:

- `walkback_3prime(read, chrom_seq, three_prime_side, stop_base="A") -> (orig, corr, applied)`
  - Walks aligned (qp, rp) pairs from `three_prime_side` inward.
  - Three-case terminal gate up front: anchored non-stop → no-op; anchored stop-base genomic match → no-op (don't drift); mismatch → walk.
  - Returns reference positions only. Knows nothing about gene strand.
- `walkback_drs(read, chrom_seq) -> (orig, corr, applied, gene_strand)` — ONT direct-RNA wrapper. BAM strand passes through.
- Constants `THREE_PRIME_SIDE_RIGHT`, `THREE_PRIME_SIDE_LEFT`, `APPLIED_WALKBACK`, `APPLIED_NONE`.

### Test suite: `tests/test_walkback_readvsref.py`

11 passing tests covering: terminal-gate Case 1 (anchored non-stop), Case 2 (anchored genomic-A), Case 3 right-side over-extension, V-primer-tip artifact (right + left), DRS wrapper on both strands, input validation, unmapped-edge case.

### What did NOT change

- `rectify/core/correct_command.py` is untouched — the new walkback core is opt-in via wrapper modules. Existing 718 unit tests pass unmodified.
- `rectify/core/atract_detector.py` is also untouched — the legacy reference-only A-tract ambiguity calc still runs. The new walkback is a *correction* that runs alongside (or eventually replaces) it.

---

## Your task: add the QuantSeq REV wrapper + CLI wiring

### Goal

When the user runs

```bash
rectify correct --short-read --dT-primed-cDNA <bam> ...
```

each read's 3'-end should be **corrected by the new read-vs-reference walkback**, and the `corrected_3ends.tsv` output column `strand` should hold the **gene/RNA strand**, not the BAM strand. Acceptance criteria below.

### File to create: `rectify/core/correct/protocols/quantseq_rev.py`

```python
"""QuantSeq REV protocol wrapper for the read-vs-reference walkback.

QuantSeq-3' REV chemistry produces cDNA antisense to mRNA. After alignment,
BAM strand is the *opposite* of the gene strand:

    is_reverse = False  →  read aligns to the genome '+' strand  →  gene is '-'
    is_reverse = True   →  read aligns to the genome '-' strand  →  gene is '+'

Poly-A tail (after pysam reverse-complementation) appears as A's at:
    is_reverse = False  →  the LEFT side of the alignment (lowest ref coord)
    is_reverse = True   →  the RIGHT side (highest ref coord)
"""
from __future__ import annotations

from typing import Tuple

import pysam

from rectify.core.correct.walkback import (
    THREE_PRIME_SIDE_LEFT,
    THREE_PRIME_SIDE_RIGHT,
    walkback_3prime,
)


def walkback_quantseq_rev(
    read: pysam.AlignedSegment, chrom_seq: str
) -> Tuple[int, int, str, str]:
    """Walkback for QuantSeq-3' REV.

    Returns ``(original_3prime, corrected_3prime, applied, gene_strand)``.

    Gene strand is the inverse of BAM strand. The 3' poly-A side of the
    aligned read maps to whichever genomic end the antisense cDNA puts it.
    """
    if read.is_reverse:
        gene_strand = "+"
        side = THREE_PRIME_SIDE_RIGHT
    else:
        gene_strand = "-"
        side = THREE_PRIME_SIDE_LEFT
    orig, corr, applied = walkback_3prime(read, chrom_seq, side, stop_base="A")
    return orig, corr, applied, gene_strand
```

### File to modify: `rectify/core/correct_command.py`

Hook the new walkback into the existing `--short-read --dT-primed-cDNA` flow. The current short-read path is at roughly lines 354–369 (look for `if is_short_read: ...` and `is_dt_primed`). You need to:

1. Import `walkback_quantseq_rev` from the new module.
2. After the existing per-read processing, call `walkback_quantseq_rev(read, chrom_seq)` to get `(orig, corr, applied, gene_strand)`.
3. Write **`gene_strand`** (not BAM strand) to the `strand` column of `corrected_3ends.tsv`.
4. Set `correction_applied = applied` (will be `polya_walkback_readgenome` or `none`).
5. The original Module 2E ambiguity-window logic can stay (still useful diagnostic), but the *correction* is now the new walkback's `corr`.

Add a feature flag if you want to be conservative — e.g., `--use-readvsref-walkback` defaulting to `True` for `--dT-primed-cDNA` and warn the user if they disable it.

### File to create: `tests/test_quantseq_rev_walkback.py`

Mirror the structure of `tests/test_walkback_readvsref.py::TestDrsWrapper` but with the **inverted strand mapping**:

```python
def test_quantseq_rev_plus_strand_gene_uses_right_side():
    """is_reverse=True → gene is '+' strand, 3' end at the RIGHT side."""
    # genome + 3 ref bases + 6 A's (overextension); V-primer terminal G
    ref = "X" * 1000 + "ACGTC" + "AAAAAA" + "X" * 100
    read = _make_read(start=1000, seq="ACGTC" + "AAAAA" + "G", cigar=((0, 11),), is_reverse=True)
    orig, corr, applied, gene_strand = walkback_quantseq_rev(read, ref)
    assert gene_strand == "+"
    assert applied == APPLIED_WALKBACK
    assert corr == 1004

def test_quantseq_rev_minus_strand_gene_uses_left_side():
    """is_reverse=False → gene is '-' strand, 3' end at the LEFT side."""
    # genome + V-primer G + 5 A's + ref body
    ref = "X" * 1000 + "AAAAAC" + "GTAC" + "X" * 100
    read = _make_read(start=1000, seq="GAAAAA" + "GTAC", cigar=((0, 10),), is_reverse=False)
    orig, corr, applied, gene_strand = walkback_quantseq_rev(read, ref)
    assert gene_strand == "-"
    assert applied == APPLIED_WALKBACK
    assert corr == 1006
```

Plus an integration test on real Han 2023 data — see acceptance criteria below.

---

## Acceptance criteria (from the ISSUE, copied verbatim)

- [ ] `walkback.py` core implements the read-vs-reference algorithm; passes a synthetic-BAM unit test where a read aligned with 5 bp polyA over-extension into a genomic A-tract is corrected back to the true CPA site. — **DONE in this pass**, see `tests/test_walkback_readvsref.py::TestWalkbackFires`.
- [ ] `protocols/quantseq_rev.py` wrapper: strand column in `*_corrected_3ends.tsv` is gene/RNA strand (verified by checking convergent locus chrXIV:196,150–196,200 in Han 2023 ends up on `–` strand for the Ysh1AA-dominant cluster). — **YOUR JOB**
- [ ] On Han 2023 wt_R1 50k subsample, `correction_applied="A-tract walking"` count rises substantially above the current ~8% baseline. — **YOUR JOB**
- [ ] Existing `rectify correct` calls without `--dT-primed-cDNA` / `--direct-rna` / `--pcb114-cdna` produce byte-identical output (no regression on non-protocol-specific use). — **gate this with an opt-in flag if needed**
- [ ] Add ONT protocol wrappers (`ont_drs`, `ont_cdna`) as separate follow-up issues; this issue ships QuantSeq REV first. — **`ont_drs` wrapper is DONE in this pass as `walkback_drs`; you don't need to do it.**

---

## Reference implementations (read these before coding)

1. **`rectify/core/correct/walkback.py`** — the core you'll be calling. Especially the docstring on `walkback_3prime`, which lays out the three terminal-gate cases.
2. **`tests/test_walkback_readvsref.py`** — the test patterns you'll be mirroring.
3. **`rectify/core/cdna_correct_command.py::walk_back_anchor_and_tail`** (lines 290–357) — the original working implementation that the core was derived from. Useful for cross-checking edge cases.
4. **Standalone reference impl** on H2 at `/u/project/guillom/kevinroy/projects/TRT/scripts/rectify/han2023/11_polya_walkback_recompute.py` (lines 56–195). This is what the lab has been running as a post-hoc fix. Read the docstring; it explains the strand convention and the V-primer-tip artifact handling. **Do NOT copy line 208** — that's the strand bug we are explicitly fixing.

---

## Validation plan

### Unit-level

Run `pytest tests/test_quantseq_rev_walkback.py -v --no-cov`. Should be ≥4 tests, all passing.

### Integration: Han 2023 wt_R1 50k subsample

On H2, run:

```bash
# Get the existing consensus BAM (or subsample to 50k)
BAM=/u/project/guillom/kevinroy/projects/TRT/processed_data/external_data/han2023_3end_seq/han2023_chunks_20260430/wt_R1.bam

# Run the new short-read path
/u/home/k/kevinroy/.conda/envs/rectify/bin/python3 -m rectify correct \
    --short-read --dT-primed-cDNA \
    --reference $REF --annotation $GFF \
    --out-dir /u/scratch/k/kevinroy/quantseq_walkback_test \
    $BAM

# Check correction_applied count
awk -F'\t' 'NR>1 {print $NF}' /u/scratch/k/kevinroy/quantseq_walkback_test/*_corrected_3ends.tsv \
  | sort | uniq -c
```

Expected: `polya_walkback_readgenome` count rises from ~8% baseline to substantially higher (target ≥50% of reads, but the exact number depends on internal-priming rate in the sample).

### Integration: convergent locus check

Pick a read in `wt_R1` that aligns near chrXIV:196,150–196,200 (the convergent ZWF1 / YNL242W locus). Confirm:

- BAM strand: `+` (from a `-` strand gene)
- Old output `strand` column: `+` (BAM strand — wrong, would be assigned to YNL242W)
- New output `strand` column: `-` (gene strand — correctly assigned to ZWF1 antisense)

Use the same locus from the original ISSUE for a known reproduction.

### Regression: non-protocol path

```bash
# Without --dT-primed-cDNA or --short-read — pure DRS path
diff <(rectify correct $DRS_BAM --out-dir A) <(rectify correct $DRS_BAM --out-dir B --use-readvsref-walkback=false)
```

(Or whatever flag you wire.) Expected: byte-identical.

---

## Strand-mapping cheat-sheet

Stick this above your monitor while you work:

| chemistry | `is_reverse` | gene strand | 3'-end side | poly-A appears at |
|---|---|---|---|---|
| ONT direct-RNA | False | + | right | reference_end - 1 |
| ONT direct-RNA | True | - | left | reference_start |
| ONT PCR-cDNA | depends on adapter orient | depends | depends | (see `cdna_correct_command.py`) |
| **QuantSeq REV** | **False** | **-** | **LEFT** | **reference_start** |
| **QuantSeq REV** | **True** | **+** | **RIGHT** | **reference_end - 1** |

The QuantSeq REV rows are inverted from DRS — that's the whole point of this issue.

---

## How to ship

1. Make changes on a branch:
   ```bash
   cd /Users/kevinroy/work/rectify
   git checkout -b quantseq-rev-walkback
   ```
2. Implement + test locally (M1 has the full test infrastructure).
3. Run the full suite: `python3 -m pytest --no-cov` (~80 s) — should still pass the original 718 + 11 from this pass + your new ones.
4. Validate on Han 2023 (above).
5. Update `.github/ISSUE_walkback_integration.md`: check off the QuantSeq REV criteria.
6. Push to GitHub master.
7. The H2 auto-sync launchd cron (`com.kevinroy.h2_rectify_sync`) picks it up within 10 min. No manual H2 action required.

---

## Notes from the original implementation pass

- The cDNA `walk_back_anchor_and_tail` returns `(canonical_cleavage_anchor, polyA_tail_length)` — the polyA length is *separately* useful and is written to the `XA:i` BAM tag for ONT PCR-cDNA. For QuantSeq REV that tag doesn't get emitted (the pipeline is short-read TSV-only), but if you want to compute tail length too, the algorithm is: count basecalled A's between the corrected anchor and the soft-clipped 3' end.
- The "stop_base='A'" parameter is on the core for future-proofing. Right now nothing changes it; everyone walks past A's. Leave it alone.
- pyabpoa is NOT involved in this codepath (no consensus generation). Only pysam + the reference fasta. The H2 environment already has both working.
- The QuantSeq REV path uses `--short-read`, which also disables a bunch of other modules (poly-A trim, indel correction). That's by design — short reads don't need them.
