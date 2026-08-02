# ONT PCR-cDNA: orientation, UMI collapse, and what is / isn't implemented

**Status:** design note, 2026-08-02. Written after the `[[541]]` strand fix, to capture nuances that
were repeatedly re-derived (and once mis-derived) by successive sessions.
Companion to `cdna_correct.md` (the Path A algorithm) — read that first for the chemistry.

---

## 1. There are TWO cDNA routes. Know which one you are debugging.

| | **Path A — UMI-aware** | **Path B — pre-UMI-collapse** |
|---|---|---|
| stages | `correct-cdna` → `align -y` → `cdna-analyze` | `trim-cdna-polya` → aligner → `correct --ONT-cDNA` |
| unit of output | one consensus per **molecule** | one row per **read** |
| orientation carried as | `XO:Z:fwd\|rev`, `XY:Z:umi_captured_fwd\|rev` | `ro:A:S\|A\|B\|U` (added 2026-08-01) |
| gene strand from | `cdna_analyze_command.py` `{fwd:'+', rev:'-'}` | `protocols/ont_cdna.resolve_rna_strand()` |
| used for | isoform / molecule-level work | internal-poly(A) work, where per-read tail and site
  *distributions* are wanted and UMI collapse is deliberately skipped (contract C9) |

Both paths present reads to the aligner in **both genomic orientations**. Neither normalises.

---

## 2. 🔴 `orient` is the GENE STRAND — not "which strand of the duplex was sequenced"

This is the single most misread thing in this subsystem, and the Stage-2 design in `cdna_correct.md`
was written on the opposite assumption.

`orient` is derived from **BAM SEQ**, which pysam always presents in reference (+) orientation. So both
strands of one duplex — which align to the same locus with *opposite* `is_reverse` — produce the **same**
BAM SEQ pattern and therefore the **same** `orient`.

**Measured** (`WT_BY4742_rep1`, `align1/*.multialigned.bam`, chrI, 60,000 reads through the real
`extract_read_info`, gene model from **CDS**):

```
orient x annotated GENE STRAND          orient x is_reverse (Type 1)
  orient=fwd  gene=+ : 57503 (96.1%)      orient=fwd  is_reverse=fwd : 34025
  orient=fwd  gene=- :   119 ( 0.2%)      orient=fwd  is_reverse=rev : 18490   <-- both, same orient
  orient=rev  gene=+ :   168 ( 0.3%)      orient=rev  is_reverse=fwd :   871
  orient=rev  gene=- :  2051 ( 3.4%)      orient=rev  is_reverse=rev :  1111
=> orient matches gene strand in 99.5% of reads
```

**Consequences**
- `{fwd:'+', rev:'-'}` in `cdna_analyze_command` is correct, and `protocols/ont_cdna.ORIENT_TO_STRAND`
  is locked to it by a test.
- The two strands of one molecule differ in **`is_reverse`**, *not* in `orient`.
- ⚠️ **Therefore an `orient=fwd` + `orient=rev` pair sharing a UMI is two loci on OPPOSITE genomic
  strands — a collision, not one molecule.** See §3.

---

## 3. The deferred "cross-orientation merge" — premise refuted; the useful half already ships

`cdna_correct.md` §Stage 2 documents a sense+antisense pair-overlap merge that pairs an `orient=fwd`
consensus with an `orient=rev` consensus sharing a UMI. Given §2 that pairs opposite-strand loci. Its
own empirical note is consistent with exactly that: of 1,086 cross-orient UMI matches on chrI at Lev≤3,
**only 51 spanned ≤3 kb; 1,035 were random collisions**. The span filter was described as mandatory
precisely because the pairing is mostly noise.

**What people actually want from it — "all reads sharing a UMI collapse into one molecule, including
both strands of the duplex" — is ALREADY IMPLEMENTED.** The Stage-1 bucket key is

```python
bucket = (r.chrom, r.anchor // anchor_window, r.orient, r.read_type)   # cluster.py
```

`is_reverse` is deliberately **absent**, so both duplex strands land in the same bucket and UMI-cluster
together; `consensus.poa_consensus_strand_aware` then builds a per-strand sub-consensus (each strand has
its own reliable end — §4) and POAs the two together. The `XB:Z:<n_top>/<n_bottom>` tag records the split.

**Measured on the real chrI Stage-1 output** (`cdna_stage1/chrI/stage1_consensus.fastq.gz`):

```
clusters = 73,035    multi-read = 13,636 (18.7%)
clusters containing BOTH duplex strands = 8,185  = 60.0% of all multi-read clusters
```

⇒ Cross-strand UMI collapse works today. **Do not re-implement it.** What is genuinely NOT collapsed is
`read_type` (§5).

---

## 4. Pore truncation is asymmetric between the two duplex strands

`_constants.py`: *"the basecalled 3' end goes through the pore first and is reliable."* Truncation
therefore eats the basecalled **5'** end, and the two strands lose **opposite ends of the mRNA**:

| sequenced strand | basecalled 5' (lost on truncation) | reliable end |
|---|---|---|
| sense (SSP…polyA) | SSP / cap side = **mRNA 5'** | **mRNA 3' (polyA)** |
| antisense (VNP…SSP) | VNP / polyT side = **mRNA 3'** | **mRNA 5'** |

This is why `poa_consensus_strand_aware` splits by strand before merging — each strand contributes the
half it sequenced accurately — and it is the real content of the Stage-2 "opposite ends" observation.

**Implications**
- A read that *has* a detected tail necessarily reached that end, so both tailed classes stay
  trustworthy. Truncation instead determines **which reads become tail-less**, and those are exactly
  the reads whose 3' end is a pore artifact rather than a cleavage site.
- Path B measurement: the tail-less (`gene_overlap`) class shows **19.0%** of 3' ends >250 nt inside the
  CDS versus **~1.7–2.2%** for both tailed classes. Gate on
  `strand_evidence ∈ {polyA_3p, polyT_5p}` for 3'-end work.
- 🔴 **If reads are ever reoriented (§6), the orientation label MUST be retained**, or this asymmetry
  becomes invisible: a reoriented antisense read's 3' end is the truncation-prone end while a sense
  read's is the reliable one, and nothing downstream could tell them apart.

---

## 5. What is genuinely still un-collapsed: Type-1 ↔ Type-2

`read_type` **is** in the Stage-1 bucket key, so Type-1 (SSP+UMI captured) and Type-2 (SSP truncated
away, `umi=""`) never merge — they are only *linked* post-hoc by `cdna-analyze`'s `XL` tag. Type-2 reads
are precisely the truncation products of §4, so this is the same phenomenon.

Merging them is the open work. It is **not** a UMI join (Type-2 has no UMI) — it has to be positional +
orientation, which carries the same collision risk that made the cross-orient span filter mandatory.
Any implementation needs a false-merge control measured up front, on the model of the chrI 1086-vs-51
analysis. Type-2 was ~9% of reads in the chrI sample above.

---

## 6. Reorientation (canonicalising all reads to mRNA sense): analysed, not adopted

- **Never planned.** No commit on any branch contains "reorient"; no TODO, no spec. The nearest ancestor
  is the §3 merge, whose goal was molecule reconstruction.
- **DRS-equivalence was an explicit goal and it was met differently** — commit `0c145ef` emits a
  DRS-identical `corrected_reads.tsv` by *projecting* `{fwd:'+', rev:'-'}`, not by transforming reads.
- **It does NOT improve mapping.** Reverse-complementing is a no-op for minimap2: same alignment, same
  score, mirrored CIGAR, only `is_reverse` differs. **Trimming** non-genomic 5'/3' elements is what
  improves mapping; reorientation is a *convention* change.
- **It would genuinely simplify the tail detectors.** Verified: `rc(_CRTA_RC[-12:]) == _CRTA[:12]`, so
  `find_cdna_3p_polya_and_adapter(rc(read))` reproduces `find_cdna_5p_polyt(read)` exactly
  (`polya_len=22, adapter=124, pass=1` on a real antisense read). The 5' detector becomes deletable.
- **Costs:** inverts `XO`'s meaning for five consumers (`cdna_analyze_command`, `read_info`'s two
  orient-branching walkers, `pretrim_consensus`, `protocols/ont_cdna` + its drift-lock test); and reads
  with no tail evidence (~22% in Path B) cannot be reoriented from sequence, so a partial reorientation
  leaves a mixed population needing *both* conventions — worse than either pure choice.

**If adopted:** cheapest at Path A Stage 1 (`io.write_stage1_fastq` already rewrites the sequence, so
emitting the RC for `rev` clusters is nearly free), as ONE coordinated change, retaining `XO` per §4,
with an explicit decision on the tail-less reads.

---

## 7. Trimming asymmetry in Path B (open, straightforward)

Path B's trim strips the **sense** read's 3' (poly-A + adapter) and, since 2026-08-01, the **antisense**
read's 5' (adapter + barcode + UMI + poly-T). It does **not** strip the sense read's 5' SSP/UMI — which
is why post-fix reads still show **53.9%** with ≥100 nt soft-clips (down from 99.5%). Path A's
`pretrim_consensus` strips both ends of both orientations. Closing this is squarely the
"trim everything non-genomic" principle; it does not affect 3'-end calls but does affect mapping quality
and any 5'/TSS analysis.

---

## 8. Reference trap that invalidated an earlier version of these numbers

In `saccharomyces_cerevisiae_R64-5-1_20240529.gff` a **`gene` feature is the union span of all mRNA
isoforms and includes UTRs** — it is *not* the ORF (TDH3: `gene` 2749 bp vs `CDS` 999 bp). Build gene
models from **`CDS`** for anything measuring distance from a terminus. Sanity check: on DRS, **≥90% of
3' ends must fall 0–250 nt downstream of the stop codon** (CDS → 93.1%; `gene` → 11.2%).
Full account: `planning/541_ont_cdna_strand_fix.md` §3b.
