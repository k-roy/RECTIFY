# RECTIFY Validation Reads

This directory contains `validation_reads.bam` and its index — a curated set of 20 nanopore
direct RNA-seq reads from *Saccharomyces cerevisiae* (strain BY4742, WT) that each exemplify
one of the four RECTIFY correction categories or chimeric alignment reconstruction.

All reads are drawn from a single DRS experiment aligned to the S288C R64-5-1 reference genome.
Each read is tagged with `XV` (validation label) and `XG` (correction category) for easy
programmatic access.

---

## Category 1 — 3' End Indel Correction

The aligner introduces deletions (and occasionally insertions) at the boundary between the
poly(A) tail and a genomic A-tract, shifting the apparent 3' end downstream of the true
cleavage and polyadenylation site (CPA). RECTIFY walks backward from the soft-clip boundary,
skipping A's and deletions, until it reaches the first non-A/T agreement — the true CPA site.

| Label | Read ID | Chrom | Strand | Coords | Key artifact |
|-------|---------|-------|--------|--------|--------------|
| `cat1_plus_1` | `9747f421-2a30-435b-84da-9925649998f9` | chrVIII | + | 508,599–508,794 | 3×D ops in last 85 bp; 19-A run at ref_end |
| `cat1_plus_2` | `0821dc9e-5f02-4308-a4c7-e818ce2c0328` | chrIX | + | 22,544–22,850 | 4×1D in last 100 bp; 5-A run at ref_end |
| `cat1_minus_1` | `fc492516-6e83-4420-b76c-4addf55e5a3a` | chrI | − | 5,068–6,331 | 9 bp poly-T clip (100% T); D12 at +14 bp inside poly-A tract |
| `cat1_minus_2` | `21a45f3a-7f72-4e83-b083-7c3858967e15` | chrXII | − | 451,530–451,719 | 9 bp poly-T clip (100% T); D12 deletes `AAAAAAAAAAAT` at +14 bp |

**Expected correction:** `three_prime_corrected` < `reference_end` (plus strand) or
`three_prime_corrected` > `reference_start` (minus strand) — i.e. the CPA site is walked
back from the raw alignment boundary.

---

## Category 2 — Soft-Clip Rescue at Homopolymer Boundary

Nanopore basecallers systematically under-call homopolymer runs. At poly(A)/poly(T) CPA
boundaries the aligner soft-clips the tail one base short, leaving the true 3' end one
position beyond the alignment. RECTIFY identifies this pattern and extends the alignment
through the remaining homopolymer bases.

| Label | Read ID | Chrom | Strand | Coords | Key artifact |
|-------|---------|-------|--------|--------|--------------|
| `cat2_plus_1` | `7e72d78c-a6ff-41e3-8603-5fa46a5211a6` | chrVIII | + | 69,703–71,916 | 17 bp 3' clip; 10-T run starts at ref_end; first clip base = C |
| `cat2_plus_2` | `2587252b-3734-49fa-973e-e6e519389348` | chrIX | + | 4,714–6,968 | 9 bp 3' clip (`AAAAAAAAA`); 8-T run starts at ref_end; A/T boundary |
| `cat2_minus_1` | `47a52aea-9413-40c3-861c-c97d9160c58d` | chrI | − | 23,726–24,735 | 15 bp 5' clip (93% T); 11-A run straddles ref_start |
| `cat2_minus_2` | `90300e53-3388-4d53-9629-301aa9c19ffa` | chrIX | − | 21,257–21,809 | 11 bp 5' clip (100% T); unbroken 11-A run crosses ref_start; up to 6 bp ambiguity |

**Expected correction:** `correction_type` includes `softclip_rescue`; corrected CPA site
shifts by 1–6 bp relative to raw alignment boundary.

---

## Category 3 — 5' End Junction Rescue

Nanopore reads that begin near a splice junction have their 5'-most bases soft-clipped rather
than placed in the upstream exon. RECTIFY identifies soft-clipped sequences at the 5' end,
locates the nearest canonical splice donor, and extends the alignment through the intron to
recover the true transcription start position.

| Label | Read ID | Chrom | Strand | Coords | Key artifact |
|-------|---------|-------|--------|--------|--------------|
| `cat3_plus_1` | `443f3f26-b8f8-47e3-9bf8-5f363409da28` | chrVIII | + | 104,292–105,214 | 620 bp 5' clip; canonical 398 bp GT-AG intron |
| `cat3_plus_2` | `2b9f301c-d755-4c25-9728-6f01a34c2372` | chrII | + | 45,625–46,429 | 49 bp 5' clip; canonical 333 bp GT-AG intron |
| `cat3_minus_1` | `a5cb1f45-8acd-4bb7-a1b0-f8b9a856187b` | chrVI | − | 53,083–54,797 | 12 bp rightmost clip (RNA 5' end); canonical 309 bp CT-AC intron |
| `cat3_minus_2` | `345ff641-573b-4c92-958e-fcb5f945a89f` | chrIX | − | 98,413–99,413 | 11 bp rightmost clip; canonical 290 bp CT-AC intron |

**Expected correction:** `five_prime_corrected` differs from raw alignment 5' end; clipped
bases re-placed in upstream exon via the rescue alignment.

---

## Category 4 — False Junction Walk-Back

Poly(A) tails cause aligners to introduce spurious N (intron-skip) operations to align tail
bases against downstream A-tracts, creating phantom splice junctions near the 3' end of the
read. The same walk-back algorithm that handles indel artifacts (Category 1) transparently
absorbs these N operations — they require no separate detection step.

| Label | Read ID | Chrom | Strand | Coords | Key artifact |
|-------|---------|-------|--------|--------|--------------|
| `cat4_plus_1` | `22e25c29-ed88-42f3-8bc8-4d0b8d3b3697` | chrXI | + | 19,592–22,073 | 1,520 bp N near 3' end; acceptor = `AAAAAA`; 24-A run downstream |
| `cat4_plus_2` | `09b04cdd-4573-40cf-95fb-6110595cfc89` | chrX | + | 392,246–393,837 | 100 bp N, 112 bp from 3' end; CT/AA flanks (non-canonical); 16-A run at acceptor |
| `cat4_minus_1` | `5b387eb1-b81c-4635-b30e-db4da94a6813` | chrVI | − | 54,492–56,010 | 745 bp N at 15 bp from ref_start; TC/AG flanks; 80% T in ref before N |
| `cat4_minus_2` | `a9706bbe-b2b1-485f-ad7c-bec458c3f448` | chrIX | − | 76,016–77,313 | 223 bp N at 11 bp from ref_start; GT/AG flanks (plus-strand canonical = wrong for −); 14-T run |

**Expected correction:** `three_prime_corrected` is walked back past the spurious N op to the
true CPA site; `correction_type` includes `false_junction`.

---

## Category 5 — Chimeric Alignment Reconstruction

These reads were processed with `rectify align --chimeric-consensus`. For each read, different
aligners produced better results for different segments. RECTIFY identified sync points where
the aligners agreed, independently scored each inter-sync segment, and assembled the final
alignment from the winning aligner per segment.

All reads are from the chrI validation dataset (SRR32518283, minimap2 + mapPacBio + gapmm2)
and carry the `XK:i:1` tag (chimeric) and `XA` tag (comma-separated winning aligners).

| Label | Read ID | Chrom | Strand | Coords | Aligners used | Segments |
|-------|---------|-------|--------|--------|---------------|----------|
| `cat5_plus_3aligner` | `SRR32518283.3077333` | chrI | + | 9,411–9,965 | mapPacBio + minimap2 + gapmm2 | 9 |
| `cat5_plus_2aligner` | `SRR32518283.2923588` | chrI | + | 7,109–9,030 | minimap2 + mapPacBio | 32 |
| `cat5_minus_long` | `SRR32518283.1221212` | chrI | − | 4,923–6,316 | minimap2 + mapPacBio | 28 |
| `cat5_minus_short` | `SRR32518283.3516917` | chrI | − | 5,044–5,683 | minimap2 + mapPacBio | 18 |

**Expected tag:** `XK:i:1`; `XA` contains multiple aligner names separated by commas;
`XS` (segment count) ≥ 3.

---

## BAM Tags

| Tag | Type | Description |
|-----|------|-------------|
| `XV` | String | Validation label (e.g. `cat1_plus_1`) |
| `XG` | String | Correction category (e.g. `cat1_indel`) |
| `XK` | Integer | 1 = chimeric alignment, 0 = single aligner (Cat 5 only) |
| `XA` | String | Winning aligner(s), comma-separated (Cat 5 only) |
| `XS` | Integer | Number of independently scored segments (Cat 5 only) |

---

## Reference Genome

Reads are aligned to *S. cerevisiae* S288C reference, R64-5-1 (2024-05-29):

```
S288C_reference_sequence_R64-5-1_20240529.chrnames.fsa
MD5: (see SGD release page)
Source: https://www.yeastgenome.org/
```

Chromosome names use the `chrI`…`chrXVI` convention (not NCBI accessions).

---

## Usage

```python
import importlib.resources
import pysam

# Access bundled BAM
with importlib.resources.path('rectify.data', 'validation_reads.bam') as bam_path:
    bam = pysam.AlignmentFile(str(bam_path), 'rb')
    for read in bam:
        label = read.get_tag('XV')
        category = read.get_tag('XG')
        print(f'{label}: {read.reference_name}:{read.reference_start}-{read.reference_end}')
```

Or via the CLI validation command:

```bash
rectify validate --show-reads
```
