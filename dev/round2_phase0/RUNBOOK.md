# Round-2 cDNA — Phase-0 Empirical Kill-Gate — RUNBOOK

**Spec:** `dev/specs/SPEC_round2_cdna_discovery_assignment_20260614.md` (MASTER) + the two designer drafts.
**Status (2026-06-14):** harness built + tested (**15/15**); substrate verified; Round-1/Round-2
empirical run scripted, **not yet executed → the GO/NO-GO verdict does NOT exist yet** (the tests
validate the machinery on synthetic data, not the hypothesis). **This is Phase 0 only** — everything in Phases 1–3
is conditional on the GO/NO-GO verdict below. *If Phase 0 returns NO-GO → stop; do not build the
production library/lift-over.*

---

## 0. What Phase 0 answers (and the GO/NO-GO bar)

The single riskiest assumption (spec §11): **does a non-trivial subset of Round-1-weak reads
get placed strictly BETTER by a contiguous aligner on a pre-spliced cDNA than by the genome
consensus — including uLTRA's own per-read record — with the anchor gate ON?**

**GO** requires ALL of (encoded in `score_phase0.Verdict.go()`):
- `n_cdna_wins > 0` (a rescuable subset exists),
- `net_hp_ed_reduction > 0` (measurable benefit),
- `trivial_win_leak ≤ 0.05` (>95% of wins carry an attributable cause; >5% no-cause ⇒ reverse-gaming ⇒ NO-GO, spec §4.2.1),
- `ultra_specific_win_frac > 0.5` (wins beat **uLTRA specifically**, not just the consensus winner — spec §9 must-fix; else redundant with annotation-guided uLTRA).

Falsifiers that also force NO-GO (spec §4.2): CPA regression; novelty suppression (a
library-absent novel-junction read force-fit onto a library isoform); no net HP-ED reduction.

---

## 1. Substrate decision (Kevin, 2026-06-14): **fresh human DRS run**

Rationale: the feature targets the compact human genome (micro-exons, multi-junction reads) —
yeast is intron-poor and a NO-GO there would be substrate-limited, not feature-limited. No
existing human Round-1 5-aligner run (with uLTRA) was found on either cluster, so Phase 0 builds
a fresh, locus-scoped human run.

### Verified resources (all confirmed reachable 2026-06-14)
| Resource | Path / locator | Notes |
|---|---|---|
| **Dataset** | SG-NEx **A549 directRNA replicate1_run1** (+rep4/5/6 for depth) | public AWS Open Data; **RNA002/R9.4.1** (not RNA004 — caveat) |
| Genome BAM (cheap-subset trick) | `https://sg-nex-data.s3.amazonaws.com/data/sequencing_data_ont/bam/genome/<SAMPLE>/<SAMPLE>.bam` | ~608 MB each; **fetch via `wget`** (Sherlock login has egress; system `curl` is broken) |
| Reference genome | `/oak/.../Projects/split_tap/newvolume/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz` | Ensembl, **chr naming `5`** — MATCHES the SG-NEx BAM |
| Annotation (Ensembl) | `/oak/.../Projects/split_tap/newvolume/genomes/Homo_sapiens.GRCh38.109.gtf.gz` | chr naming `5` — use this (not GENCODE `chr5`) to avoid a rename |
| Aligners | rectify conda env: `samtools minimap2 uLTRA deSALT` present; **gmap MISSING** (Round-2 tie-breaker only, not needed) | `source /home/groups/larsms/users/kevinroy/anaconda3/etc/profile.d/conda.sh; conda activate rectify` |
| Working dir | `/oak/.../Users/kevinroy/projects/rectify_round2_phase0/` | `data/ logs/ coverage/` |

**chr-naming:** BAM + genome fasta are Ensembl `5`/`6`. Keep everything Ensembl-style; the
Ensembl GTF.109 matches. (GENCODE `chr5` would need a `sed 's/^chr//'`.)

**env trap:** `conda activate rectify` trips `JAVA_HOME: unbound variable` under `set -u`
(CLAUDE.md `set -u` trap). All scripts here use `set -eo pipefail`, activate, THEN `set -u`.

---

## 2. The harness (built + tested on M1, 10/10 passing)

| File | Role |
|---|---|
| `liftover.py` | transcript→genome lift-over with N-op insertion; `Block`/`CdnaModel` block-map; minus-strand revcomp+CIGAR-reversal; the §2.5 invariants incl. the **orientation invariant** that catches a missed revcomp. Emits a duck-typed `LiftedRead`. |
| `score_phase0.py` | the win-guard + **BLOCKER-1 fix** + verdict. Imports the **production** scorers (`_cigar_hp_edit_distance`, `_cigar_aligned_bases`, `_cigar_min_junction_anchor`) — so the verdict uses the same metric production winner-selection uses. |
| `test_phase0.py` | 15 spec-mandated unit tests (lift both strands, **boundary-break N-insertion regression**, ≥3 N-ops, E5/E6 boundary I/D, minus-strand multi-junction E11, orientation-invariant catch, trivial-win defeat §1.3, micro-exon win, IR-lose, no-shrink veto, BLOCKER-1, gate-off safety, verdict aggregation). |

> **Lift-over bug fixed 2026-06-14 (advisor catch):** N-insertion originally fired only when a
> *single* op spanned a block boundary, so ops breaking *exactly* at an exon boundary
> (`[(M,20),(X,1),(M,19)]`) silently dropped the intron — fusing exons and corrupting HP-ED,
> worst on micro-exon reads (Phase 0's target). Now N-insertion is driven by the transcript walk
> crossing a block boundary, tracked globally across the whole CIGAR. Regression test:
> `test_boundary_break_emits_n`.

**Run the tests:** `cd ~/work/rectify && /Users/kevinroy/miniconda3/bin/python dev/round2_phase0/test_phase0.py`
(base miniconda python imports numpy + rectify on M1).

**Key design choices baked in:**
- **BLOCKER-1 (spec §7):** `anchor_gate_passes()` gates on the **RAW** `min_junction_anchor`,
  never on `_effective_chimera_ok`. A candidate with `five_prime_rescued=1` does **not** get a
  free pass — the test `test_blocker1_rescued_does_not_bypass_gate` proves a sub-K anchor +
  rescued candidate is REJECTED.
- **Gate-off safety (spec §3.2):** `anchor_gate_passes(..., min_junction_anchor_bp=0)` raises —
  Round-2 refuses to run with the gate off.
- **Lift-over compatibility:** `LiftedRead` exposes exactly the 4 attributes the production
  scorers read (`cigartuples`, `reference_name`, `reference_start`, `query_sequence`).

---

## 3. Empirical findings so far (the gate already redirected the work)

1. **rep1 is shallow** (~71k reads across chr5+chr6: 41,367 + 30,522). The 6 a-priori microexon
   genes (DIAPH1/MYO10/ABLIM3/KIF3A/TRIO/MPC1 — see `apriori_stress_loci_NOTES.md`) have **~0
   coverage** in rep1: they aren't expressed in A549 at usable depth. **The empirical gate caught
   this before any library was built** (exactly its purpose).
2. **But the substrate IS viable:** **27,357 spliced reads** (N in CIGAR) in rep1 chr5+6, with
   deep multi-junction peaks. Top spliced-read loci (50 kb bins) → genes:
   `6:73.5M EEF1A1 (2474)`, `5:82.25M ATG10 (1042)`, `6:33.25M RPS18 (1024)`, `5:181.2M TRIM41`,
   plus several deep peaks not yet gene-mapped.
3. ⇒ **Locus selection must be data-driven** (coverage ∩ stress-structure), not an a-priori list.
   Deepening pool (rep1+4+5+6) downloading to `data/a549_pooled_chr5_6.bam`.

---

## 4. Remaining procedure (resume here)

### Step A — finish depth pool (running)
`logs/02_deepen.log` → produces `data/a549_pooled_chr5_6.bam` (rep1+4+5+6, chr5+6).

### Step B — data-driven locus selection (coverage ∩ stress-structure)
The principled selection (A549 is SRRM4-low → don't bet on neuronal microexons; use short
internal exons in the isoform A549 actually transcribes — see `apriori_stress_loci_NOTES.md`):
1. From the Ensembl GTF, list internal exons ≤30 bp (and genes with ≥3 junctions in <3 kb) on chr5/6.
2. Intersect with **per-gene spliced-read depth** in the pooled BAM (`samtools view -c` per gene).
3. Keep genes with (depth ≥ ~50 spliced reads) AND (≥1 short internal exon OR dense junctions).
4. **Confirm short-exon INCLUSION**: `samtools depth` over the short exon vs a flanking
   constitutive exon (PSI proxy) — the exon must actually be included in A549 reads.
Output → `coverage/selected_loci.tsv` (gene, chrom, strand, exon chain, stress feature).

### Step C — build the (hand) cDNA library for the selected loci
For each selected gene: emit the MANE/canonical exon-chain cDNA **+ any de-novo gate-passed
chain observed in Round-1**, padded 1–2 kb each side (real genomic seq, intergenic-clipped),
with the `CdnaModel` block-map (`liftover.Block`). Validate per-block 1:1 identity
(`CdnaModel.validate(genome)`). **Spike test:** a read with a known non-canonical gate-passed
junction must survive into the library (spec §3 — don't let a collapse step re-impose GT-AG).

### Step D — Round-1 (genome) on the selected-locus reads — **CHUNKED, gate ON**
- Extract pooled reads at the selected loci → fastq (`samtools view -b <regions>` → `samtools fastq -F 0x900`).
- Run the genome aligner panel (minimap2 `-ax splice -uf`, uLTRA [annotation-guided], deSALT)
  → per-aligner `rectify correct` → per-aligner corrected BAMs/TSVs, **`min_junction_anchor_bp=10`**.
- **MANDATORY chunking** (CLAUDE.md): array task per locus/sample, idempotent skip-checks,
  `--partition=owners` + AVX-512 constraint (`sherlock-sbatch` skill emits this). Keep uLTRA's
  per-read record (needed for the uLTRA-specific-win metric).

### Step E — identify the Round-1-weak subset (post-hoc, no re-run)
Per read spanning ≥1 junction AND any of: effective soft-clip ≥ 20 bp · high HP-ED · failed
anchor gate · low MAPQ / multi-aligner disagreement (spec §6). Compute with the production
`_cigar_*` functions on the existing corrected BAMs.

### Step F.0 — lift-over validation on REAL reads (do BEFORE trusting any verdict)
Before scoring, round-trip-test the lift-over on real SG-NEx reads (advisor): take a read
minimap2 already spliced cleanly across a multi-exon gene, synthesize that gene's cDNA + block
map, derive the read's cDNA-local alignment, `lift` it back, and assert the lifted CIGAR +
`reference_start` reproduce the read's original genome alignment. This exercises lift-over on
real-error (multi-junction, I/D/X) reads on both strands — the synthetic unit tests cannot.

### Step F — Round-2 (cDNA) on the weak subset → lift → score
- Align weak reads to the padded cDNA library: `minimap2 -ax map-ont -N 100 -p 0.8 -Y` (NO splice,
  NO `-G`) + mapPacBio (`intronlen=0 maxindel=100`). `--for-only` only if pre-oriented (DRS: yes).
- `liftover.lift(cdna_model, cdna_cigar, cdna_ref_start, read_query)` → `LiftedRead`; assert the
  §2.5 invariants (esp. orientation invariant on a known-clean read).
- Build `Candidate.from_read(...)` for the lifted cDNA, the genome winner, and uLTRA's record;
  `score_phase0.run_verdict(triples, min_junction_anchor_bp=10)` → `Verdict.summary()`.

### Step G — verdict + provenance
- Print `Verdict.summary()`; GO/NO-GO per §0. Emit `PROVENANCE.json` (dataset replicates,
  refs+versions, commit, `min_junction_anchor_bp`, ε, selected loci, render date).
- **NO-GO ⇒ STOP** and record why (which falsifier). GO ⇒ proceed to Phase 1.

---

## 5. Open risks / decisions
- **RNA002 not RNA004** — only public A549 DRS option; intron/micro-exon placement is
  chemistry-independent, but flag if the verdict is marginal.
- **uLTRA needs the Ensembl GTF index** (annotation-guided) — build once for chr5/6.
- **Pad must be real genomic sequence** (never N) or the 3' poly-A walkback breaks (designer-1 §2.2).
- **mapPacBio HP-ED gaming** (memory `feedback_hped_overcall_gaming`): kept here as a *secondary*
  Round-2 aligner with the anchor gate ON; watch its sole-wins in the attribution table.
