# EXPLORED ROUTES — dead ends and why they closed

Routes recorded here were explored seriously, with a working apparatus and a real empirical
test. The verdict is substrate+commit-specific: if the root cause changes (different substrate,
different algorithm version), revisit the relevant section rather than treating it as permanently
closed. Each entry states the mechanism that killed it, so a future reader can judge whether
that mechanism still applies.

---

## 1. Round-2 cDNA "discovery → assignment" alignment

**Dates:** scoped 2026-06-14; A549 kill-gate 2026-06-14; brain RNA004 kill-gate 2026-06-15
**Spec:** `dev/specs/SPEC_round2_cdna_discovery_assignment_20260614.md`
**Apparatus:** `dev/round2_phase0/` (harness + cluster scripts)
**Cluster workdirs:**
- A549: `sherlock:/oak/.../projects/rectify_round2_phase0/round2/`
- Brain: `sherlock:/oak/.../projects/rectify_round2_phase0/brain_rna004/`

### What was proposed

After Round-1 genome alignment + correction, a subset of reads are weak (high HP-ED, failed
anchor gate, heavy soft-clips at junctions). Hypothesis: aligning those reads to a library of
pre-spliced cDNA isoforms with a **contiguous** aligner (no intron discovery required), then
lifting the alignment back to genome coordinates, would place them better than any genome
aligner — recovering exonic bases soft-clipped at junction seeds and placing micro-exons that
seed-aligners skip.

The scoring is symmetric: N-ops cost 0 in both the genome and lifted-cDNA CIGAR, so the cDNA
gets no free bonus for the missing intron. A win requires genuine exonic HP-ED reduction.

The library admits MANE isoforms ∪ Round-1 read-supported de-novo chains that pass the k≥10
perfect-match anchor gate — so novel junctions enter, GT-AG is not imposed, and the Round-2
aligner is not restricted to canonical isoforms.

### Phase-0 kill-gate design

The spec mandated an empirical kill-gate before any production build. GO required all four:
- `n_cdna_wins > 0`
- `net_hp_ed_reduction > 0`
- `trivial_win_leak ≤ 0.05` (≥95% of wins carry an attributable structural cause)
- `ultra_specific_win_frac > 0.5` (wins beat uLTRA specifically, not just weaker aligners)

Additional falsifiers: CPA regression (>12 bp 3′ shift), novelty suppression.

### Substrate 1 — SG-NEx A549 directRNA, RNA002 (2026-06-14)

**Dataset:** rep1+4+5+6 pooled, chr5+6. Ensembl GRCh38 primary_assembly + GTF.109.
**7 loci:** TGFBI, IK, HSD17B4, ANXA6, TRIO (56 junctions), XPO5 (31 junctions), DIAPH1.
**Round-1 panel:** minimap2 + uLTRA + deSALT, per-aligner `rectify correct`.
**F.0 lift validation:** 833 real reads round-tripped, 660+/173−, **0 mismatches**.

| Metric | Value | GO threshold |
|---|---|---|
| Eligible reads | 4,579 | — |
| cDNA wins (vs corrected genome) | 417 | >0 ✓ |
| Structural wins (attributable cause) | **7** | meaningful |
| Leak (no attributable cause) | **90%** | ≤5% ❌ |
| CPA-regression falsifier | 0 | =0 ✓ |
| Net HP-ED reduction | +1,060 | >0 ✓ |

Three independent controls all agree:
1. **Corrected-genome baseline:** only 7/4,579 reads (0.15%) get a genuine structural win.
2. **Single-exon null:** single-exon reads (should never win) won at 8.87%+ from alignment
   jitter — the false-win floor already exceeds the 5% leak bar before any spliced wins are
   counted.
3. **Genome-unmapped bias check:** 0/844,262 genome-unmapped reads get a confident spliced
   rescue (all 17 spliced cDNA hits are 68–83%-identity force-fits of low-quality reads).

**Caveat:** A549 is SRRM4-low (micro-exon-poor). Multi-junction use case was well-powered
and tested-and-weak. Micro-exon use case was under-tested.

### Substrate 2 — ENCODE DLPFC ctrl85 (brain), RNA004 (2026-06-15)

**Dataset:** ENCFF043YFO, chr5+6 (107k + 98k reads). GENCODE GRCh38, chr-naming.
**7 loci:** DIAPH1, TRIO, MYO10, ABLIM3, KIF3A, MPC1, COL23A1.
**Micro-exon inclusion (confirmed in brain data):**

| Gene | Micro-exon | Inclusion ratio |
|---|---|---|
| TRIO | 28 bp | 95.8% (23/24 reads) |
| DIAPH1 | 9 bp | 100% (7/7) |
| ABLIM3 | 8 bp | 100% (10/10) |
| KIF3A | 9 bp | 23.1% (3/13) |
| MPC1 | 4 bp | 0% — excluded |

**Round-1 panel:** minimap2 + uLTRA + deSALT, per-aligner `rectify correct` (1,081 reads,
100% high-confidence). **F.0 lift validation:** 98 real reads, 0+/88−, **0 mismatches**.

| Metric | Value | GO threshold |
|---|---|---|
| Eligible reads | 1,065 | — |
| cDNA wins | 297 | >0 ✓ |
| Structural wins | **7** | meaningful |
| Leak (no attributable cause) | **94.9%** | ≤5% ❌ |
| CPA-regression falsifier | **5** | =0 ❌ |
| Net HP-ED reduction | +787.6 | >0 ✓ |

**NO-GO.** Two falsifiers fire: leak at 94.9% and 5 CPA-regression wins.

### Why it closed

**Root cause:** on both substrates, vanilla minimap2 already handles essentially all reads at
these loci. There is no headroom-creating population — a pool of reads that the genome aligner
systematically fails to seed across a junction but a contiguous cDNA aligner would place
cleanly. The 297/417 "wins" are dominated by alignment jitter: the cDNA candidate scores 2–3
HP-ED points better by chance, not by recovering exonic sequence. The 7 structural wins are
identical in count across substrates of very different sizes (1,065 vs 4,579 reads), strongly
suggesting they are a noise-floor artifact.

The micro-exon sub-case was specifically targeted in the brain retest: TRIO (28 bp, 95.8%
included), DIAPH1 (9 bp, 100%), ABLIM3 (8 bp, 100%) are all expressed and included. The
result (7 structural wins, 94.9% leak) is indistinguishable from A549. The micro-exon use
case does not create a rescuable population on ONT DRS at these loci.

### What stays valid (do not discard with the feature)

- **The apparatus is correct and tested.** The lift-over (`liftover.py`), win-guard
  (`score_phase0.py`), F.0 round-trip gate (`liftover_roundtrip.py`), and projection with
  transcript-contiguity enforcement (`cdna_models.py`) all work correctly — F.0 returned 0
  mismatches on real reads on both substrates. The two lift bugs caught by F.0 (boundary-break
  N-insertion; semantic CIGAR comparison) were fixed and guarded by `test_projection_local.py`.
- **The corrected-genome baseline is the right baseline.** A raw genome baseline is a
  false-GO generator: `rectify correct`'s 5′/soft-clip rescue competes with the cDNA round's
  exact mechanism. Always score cDNA vs the corrected winner + corrected uLTRA.
- **The kill-gate design worked.** Phase 0 cost ~2h of cluster time and caught the dead end
  before any production library/lift-over infrastructure was built. The gate (leak ≤5% +
  structural-cause attribution + CPA-regression check) is the right bar.

### Conditions that would justify revisiting

The NO-GO is substrate- and algorithm-specific. It would be worth revisiting if:
- A substrate can be identified where a substantial fraction of reads span junctions that
  minimap2 consistently fails to seed (e.g., extremely compact micro-exons <6 bp, or reads
  from a highly divergent strain/species where k-mer seeding breaks down systematically).
- The aligner panel changes in a way that leaves a larger weakened-read pool (e.g., if
  uLTRA is replaced with a less annotation-sensitive aligner).
- Direct evidence from IGV shows reads that visually span a junction but have a soft-clip
  at one end that a cDNA aligner fixes — that would be an existence proof justifying a
  targeted re-run on those specific reads.

Without such evidence, do not re-open: the apparatus is ready and the prior verdicts are
the answer.
