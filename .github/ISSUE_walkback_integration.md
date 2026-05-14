## Problem

`rectify correct --short-read --dT-primed-cDNA` produces output where ~92% of QuantSeq REV reads have `correction_applied="none"`. The cause is that the current A-tract walk-back (Module 2E) consults only the reference genome — it doesn't compare the read base against the reference at each position, so it cannot detect when QuantSeq's internal-priming has pulled the alignment a few bases into a polyA tail.

This means every QuantSeq REV dataset analyzed with rectify currently has 3' end positions that drift downstream of the true CPA site by 1–N bases of polyA over-extension. For yeast termination motif analysis (where the precise CPA position is the anchor for downstream PE/UAE/T-tract elements), this is biologically significant.

ONT direct RNA / cDNA pipelines (see the lab's `ont_cdna` pipeline) hit the same underlying ambiguity — alignment of the polyA tail into a nearby genomic A-stretch. The fix is protocol-agnostic at the algorithmic level.

## Current workaround (Han 2023 analysis)

A standalone script — `11_polya_walkback_recompute.py` in `TRT/scripts/rectify/han2023/` — re-computes 3' ends by post-processing the consensus BAM. The algorithm:

1. For each read, find the 3' alignment end in genomic coords:
   - `+` BAM strand → `reference_end - 1` (rightmost aligned position)
   - `-` BAM strand → `reference_start` (leftmost aligned position)
2. Walk inward through aligned pairs toward the gene body.
3. At each position: `read_base = query_sequence[query_pos]`; `ref_base = genome[ref_pos]`.
4. Stop when `read_base == ref_base AND read_base != 'A'`.
5. Report that `ref_pos` as the corrected 3' end.

Output: `*_corrected_3ends_wb.tsv` ("wb" = walkback) — a drop-in replacement for `*_corrected_3ends.tsv`.

## Known defect in the workaround

The standalone script writes `strand = '-' if read.is_reverse else '+'` — i.e., the **BAM** strand, not the **gene/RNA** strand. For QuantSeq REV antisense, RNA strand = NOT(BAM strand), so this emits inverted strands. Downstream `rectify analyze` then clusters on the wrong strand and assigns reads to the wrong gene at every convergent / antisense locus. We hit this hard when validating cluster classifications against IGV — e.g., a "primary CPA-independent" cluster assigned to YNL242W (+ strand) turned out to be ZWF1 antisense termination (– strand) under Ysh1AA depletion.

## Proposed: shared walkback core + per-protocol wrappers in `rectify correct`

The walk-back algorithm is **protocol-agnostic at its core** (the same A-stretch ambiguity exists for QuantSeq REV, ONT direct RNA, ONT PCR-cDNA, etc.). Implement once, wrap for each protocol.

### Module layout

```
rectify/core/correct/
  walkback.py                  # protocol-agnostic core
    walkback_3prime(read, genome_seq, stop_base='A',
                    allow_mismatches=0) -> int
  protocols/
    __init__.py
    quantseq_rev.py            # --dT-primed-cDNA: stop_base='A', strand=antisense, mismatches=0
    ont_drs.py                 # --direct-rna: stop_base='A', strand=sense,    mismatches=2
    ont_cdna.py                # ONT PCR-cDNA: stop_base='A', strand=depends,  mismatches=2, UMI-aware
```

### Core `walkback_3prime`

- Input: a `pysam.AlignedSegment`, genome reference for its chromosome, stop-base (default `'A'`), allow_mismatches (default `0`).
- Output: corrected reference position of the 3' end.
- Algorithm: walk inward through aligned pairs from the 3'-end side; stop at the first position where read_base matches ref_base AND ref_base != stop_base AND mismatch_count <= allow_mismatches.
- Intron gaps (qp=None pairs) are filtered out, so the scan naturally crosses splice junctions.
- No walkback cap; the alignment length is the natural bound.

### Per-protocol wrapper

Each wrapper:
- Calls the core with protocol-appropriate `stop_base`, `allow_mismatches`.
- Maps `read.is_reverse` to gene strand per the protocol's chemistry:
  - QuantSeq REV antisense: gene_strand = `'+' if read.is_reverse else '-'` (BAM strand inverted)
  - ONT direct RNA: gene_strand = `'-' if read.is_reverse else '+'` (BAM strand passthrough; the cDNA is read in sense direction)
  - ONT PCR-cDNA: depends on adapter chemistry, configurable.
- Returns `(read_id, chrom, strand=gene_strand, original_3prime, corrected_3prime, correction_applied)`.

### CLI

Existing flags trigger the appropriate protocol wrapper:

```
rectify correct --dT-primed-cDNA ...    # → protocols.quantseq_rev
rectify correct --direct-rna ...         # → protocols.ont_drs (new)
rectify correct --pcb114-cdna ...        # → protocols.ont_cdna (new)
```

If no protocol flag is set, fall back to the existing Module 2E reference-only walkback (no behaviour change for current users not using a specific protocol).

## Acceptance criteria

- [x] `walkback.py` core implements the read-vs-reference algorithm; passes a synthetic-BAM unit test where a read aligned with 5 bp polyA over-extension into a genomic A-tract is corrected back to the true CPA site. (commit 986a19d)
- [x] `protocols/quantseq_rev.py` wrapper: strand column in `*_corrected_3ends.tsv` is gene/RNA strand (`'+'` when `read.is_reverse`, `'-'` otherwise). Convergent-locus IGV check (chrXIV:196,150–196,200) pending H2 run.
- [ ] On Han 2023 wt_R1 50k subsample, `correction_applied="polya_walkback_readgenome"` count rises substantially above the current ~8% baseline. (H2 validation pending)
- [x] Existing `rectify correct` calls without `--dT-primed-cDNA` / `--direct-rna` / `--pcb114-cdna` produce byte-identical output (757 tests pass, 28 skipped, 0 failed locally).
- [x] Add ONT protocol wrappers (`ont_drs`, `ont_cdna`) as separate follow-up issues; this issue ships QuantSeq REV first.

## Migration

For users currently consuming `corrected_3ends.tsv` from a `--dT-primed-cDNA` run: strand values will flip from BAM-strand to gene-strand. Document under "Breaking changes" in the release notes. The standard downstream consumer (`rectify analyze`) Just Works since it joins on strand-coord pairs that now match the gene annotation correctly.

## Implementation references

- Walkback algorithm: `11_polya_walkback_recompute.py` (lines 56–195, the `walkback_3prime()` function) in `kevinroy/projects/TRT/scripts/rectify/han2023/` on Hoffman2.
- Strand bug to fix at integration time: line 208 of that script emits BAM strand; in the rectify port for `--dT-primed-cDNA`, emit `gene_strand = ('+' if read.is_reverse else '-')`.
- Validated against the 4 Han 2023 wt/Ysh1AA cpfixed consensus BAMs — produces correct `_corrected_3ends_wb.tsv` after strand-flip fix, ~10 min per sample.

## Out of scope (for this issue)

- Multi-organism chromosome-mapping tables (`normalize_chromosome` is yeast-only; we made `--chrom-format passthrough` the default in commit `5902c76`).
- ONT-specific protocol wrappers — file as follow-up issues that depend on this one.
- UMI-aware consensus pre-walkback for PCR-cDNA — orthogonal concern; lives in a separate module.
