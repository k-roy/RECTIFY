# QNAME Mutation Survey — RECTIFY 4-Aligner Panel

**Date:** 2026-05-20
**Scope:** Characterize how each aligner in the RECTIFY production panel mutates
read IDs (QNAMEs) between FASTQ input and BAM output.
**Method:** 53-read synthetic adversarial FASTQ + 3 real-data spot-checks.
**Platform:** M1 (minimap2, mapPacBio, deSALT, gapmm2 binaries); Sherlock for
production BAM verification. All synthetic runs executed locally — character-level
QNAME behaviour is deterministic across platforms.

---

## 1. Per-Aligner Mutation Summary

### minimap2 (raw, no -y flag)

Raw binary (minimap2 -ax splice -uf -k14) truncates the FASTQ header at the
**first whitespace** before placing the QNAME into SAM. Any content after the
first space or tab is silently dropped and does NOT propagate to the BAM. Tab-
separated FASTQ comment tags (@id\tXA:Z:foo) are also stripped because the
tab counts as whitespace. QNAMEs >254 chars cause `samtools view` to fail with
"query name too long". NUL bytes (\x00) in the QNAME body act as a C-string
terminator, silently truncating the QNAME. Duplicate QNAMEs pass through
unchanged and both copies appear in the BAM. All printable non-whitespace special
characters in the QNAME body are preserved.

**Gap:** The raw binary provides no mechanism to carry FASTQ auxiliary tags into
the BAM — they are always dropped.

### minimap2 (wrapped, production -y flag)

The rectify run_minimap2() wrapper uses -y (plus -ax splice -uf -k14).
With -y, minimap2 copies tab-separated SAM-format tags from the FASTQ comment
into the BAM aux records. QNAME handling is otherwise identical to no-y:
whitespace truncation, 254-char limit crash, NUL truncation all apply. FASTQ tags
are correctly promoted to BAM aux fields (XA:Z, XC:i, XF:Z all preserved).
Space-separated comments (Dorado runid/ch/start_time) are still dropped from
the QNAME because minimap2 treats the first token as the read name regardless of
-y; they do NOT become aux tags.

**Gap:** Dorado-style space-separated metadata is silently dropped and
cannot be recovered from the BAM. The -y flag only propagates tab-separated
SAM-format tags.

### mapPacBio (raw binary, no rectify sanitizer)

Raw mapPacBio.sh copies the **full FASTQ header line** verbatim into the SAM
QNAME field -- spaces, tabs, and all trailing content are included. This violates
the SAM spec (QNAME may not contain spaces; max 254 chars). samtools view
then fails with "query name too long" for Dorado-style headers. Tab characters
are converted to underscores by BBMap before emission (so @id\tXA:Z:foo
becomes id_XA:Z:foo in the QNAME). This COMMENT_LEAKED behaviour affects all
records with any comment: space-only trailing space, Dorado metadata, Illumina
pair tags, tab-separated aux tags.

**Gap (raw):** Every FASTQ record with a comment produces a malformed BAM QNAME.
No length cap is applied; records with 300+ char headers cause samtools failures.

### mapPacBio (wrapped, via rectify)

The rectify run_map_pacbio() wrapper includes a FASTQ pre-sanitizer (commits
838293c and e8c8070) that strips everything after the first space in each
header line and truncates to 254 chars before passing the FASTQ to mapPacBio.

**Critical bug confirmed:** The sanitizer is gated on a per-file first-record
check (_need_san = _first.startswith('@') and ' ' in _first). If the first
record in the FASTQ has a bare QNAME (no space), _need_san = False and the
ENTIRE file is processed unsanitized. In the synthetic test FASTQ, bare_id is
first -- the sanitizer fires for 0 records, and all subsequent reads with Dorado
or Illumina comments leak into the BAM QNAME field unchanged. This means the
wrapped mapPacBio still emits COMMENT_LEAKED QNAMEs for any FASTQ whose first
read has no comment.

In addition, tab-separated comments are not handled by the sanitizer at all -- the
' ' in _first check passes for space-delimited reads but the tab case is missed,
so @id\tXA:Z:foo becomes id_XA:Z:foo in the QNAME even with the sanitizer active.

A post-alignment _pt:i:N suffix strip (lines 696-716) correctly removes BBMap's
poly-A tail annotation from QNAMEs after sorting. This works on production DRS data
(confirmed: 0 occurrences of _pt:i: in first 1000 QNAMEs of set1 wt_rep1).

**Gap (wrapped):** (1) Sanitizer is per-file, not per-record -- mixed-comment
FASTQs are not fully sanitized. (2) Tab-delimited FASTQ comments bypass the
sanitizer entirely.

### gapmm2 (raw and wrapped)

gapmm2 requires a ts:A:+ or ts:A:- strand tag in the FASTQ comment to
determine read orientation. Reads without this tag produce an empty PAF output
(all 53 synthetic reads: 0 records). This is the reason the previous survey
recorded all-CRASHED for gapmm2 -- it is not a crash but a silent 0-record output
due to missing strand annotation.

When reads DO have the ts: tag (real DRS data), gapmm2 outputs the QNAME into
PAF field[0] as the first whitespace-delimited token. The _paf_to_bam converter
uses PAF field[0] directly as the BAM QNAME. The _clean_fastq() pre-processor
(called by the wrapped version) strips all FASTQ comments before gapmm2 sees
them. Both raw and wrapped invocations produce IDENTITY QNAMEs (bare UUID) for
real DRS data.

Confirmed from set1 wt_rep1 production: gapmm2 BAM has bare UUIDs, matching
minimap2 and deSALT.

**Gap:** gapmm2 produces 0 reads for any input FASTQ without ts:A:+/- in the
comment. For non-DRS inputs (QuantSeq REV, cDNA without ts tag), the aligner
silently contributes nothing to consensus.

### deSALT (raw and wrapped)

Both raw and wrapped deSALT pass QNAMEs through identically as IDENTITY. The
_clean_fastq() pre-processor strips all FASTQ comments before deSALT sees them,
so the BAM QNAME is always the first whitespace token of the original header.
Duplicate QNAMEs are deduplicated by _clean_fastq() (second occurrence is dropped).
All special characters in the QNAME body are preserved.

**Gap:** Duplicate UUID deduplication is silent -- no warning in the TSV. Second
occurrence is simply dropped.

### splice_aware_5prime

NOT a separate aligner binary. splice_aware_5prime.py (Module 2F) is a
post-consensus CIGAR surgery module. The production "splice-aware" aligner is
the standard run_minimap2() wrapper (-ax splice -uf -k14 -y), which produces
IDENTITY QNAMEs. See the minimap2 (wrapped) entry.

---

## 2. Edge Cases That Crash Any Aligner

- **QNAME >254 chars:** samtools view exits 1 with "query name too long" for
  minimap2 (raw_no_y and wrapped) and mapPacBio (raw and wrapped). deSALT and
  gapmm2 use _clean_fastq() and are insulated from over-length comments.
- **NUL byte in QNAME body (\x00):** All aligners using C string handling
  truncate the QNAME at the NUL byte silently. Observed for minimap2 and
  mapPacBio (QNAME id_nul\x00byte -> id_nul).

---

## 3. Gaps in Current Sanitizer / Compatibility Check

### _normalize_bam_read_name() gaps

1. **tab-delimited FASTQ comments in mapPacBio QNAME:** id_tab_XA_XA:Z:foo
   appears in the BAM QNAME. _normalize_bam_read_name() does NOT strip
   _XA:Z:foo because it is not in the regex. This is a live gap: tab-comment
   mapPacBio BAMs will fail the K-way merge unless the regex is extended.

2. **space-delimited Dorado metadata in mapPacBio QNAME:** After samtools sort,
   spaces become underscores. id_dorado_style_runid=abc_ch=42_... in BAM QNAME
   is not covered by the current regex. Only _pt:i:N and _N_length=N are stripped.

### _check_read_name_compatibility() gap

The 50%-overlap check fires at consensus time, not alignment time. Sanitizer bugs
surface late, requiring a full re-alignment to fix.

### mapPacBio pre-sanitizer per-file gating bug

_need_san is evaluated once on the first record only. Any FASTQ whose first
record lacks a space comment silently bypasses sanitization for the entire file.

---

## 4. Production Impact (Part E)

**Sample:** set1 wt_rep1 DRS (production, v3_20260429 pipeline, Sherlock).

| Aligner | BAM read count |
|---|---|
| minimap2 | 11,126,320 |
| mapPacBio | 11,120,709 |
| gapmm2 | 8,958,909 |
| deSALT | 10,389,816 |
| corrected_3ends.tsv (23 chunks) | 9,776,410 rows |

**Conclusion:** Bug is NOT currently firing in DRS production. The _pt:i:N strip
is working (0 occurrences of _pt:i: in first 1000 mapPacBio QNAMEs). The 9.8M
corrected output from 11.1M minimap2 reads (88%) represents normal consensus
attrition. No QNAME-mismatch dropout detected.

The mapPacBio sanitizer per-file gating bug is a **latent risk** for production
DRS FASTQs: if the first chunk read happens to have no Dorado comment (e.g., a
bare UUID from an older Dorado version), the sanitizer skips the entire chunk and
subsequent reads will have COMMENT_LEAKED QNAMEs. This has not been observed in
set1 but is structurally possible.

---

## 5. Real-Data Sanity Checks (Part D)

| Protocol | Source | QNAME form | Pattern |
|---|---|---|---|
| DRS (ONT Dorado) | set1 wt_rep1 chunk FASTQ | @<uuid>\tpt:i:<N> | All 4 aligners produce bare UUID (IDENTITY); _clean_fastq strips tab+pt:i; mapPacBio sanitizer strips space comment (when first read has space) |
| QuantSeq REV | BY4742 QSrev wt_R1 BBmap (pre-fix) | @SRR.accession <N> length=76 | BBmap emits full comment; BWA emits bare accession -- MISMATCH confirmed; _normalize_bam_read_name handles underscore-encoded form |
| PCR-cDNA (ONT) | validation_reads_cdna.bam | @<uuid> (bare) | IDENTITY across all aligners (no comment in these reads) |

---

## 6. Recommendations

### (a) Fix mapPacBio sanitizer to be per-record (highest priority)

Multi-aligner commit cb2fe6c added trd=t to BBmap (fixes QuantSeq) and the
_normalize_bam_read_name() extensions. The mapPacBio sanitizer (838293c) is
still per-file. Change the gating:

```python
# CURRENT (per-file -- wrong):
_need_san = _first.startswith('@') and ' ' in _first
if _need_san:
    # sanitize whole file
```

to always run the per-record sanitization pass unconditionally. The cost is one
extra FASTQ read pass (~1-2 seconds per chunk), which is negligible.

### (b) Extend _normalize_bam_read_name() to cover tab-comment mapPacBio artifacts

Add to _UNDERSCORE_COMMENT_RE coverage for:
- _XA:Z:<value> and similar cDNA-pipeline tag patterns
- _runid=<value>... (Dorado metadata, underscore-encoded by samtools sort)

Or apply a generic "strip everything after first _XX:T:" pattern that matches
any SAM aux tag encoded as underscore after samtools sort.

### (c) Add post-alignment FASTQ-to-BAM QNAME validator

After each per-aligner run, validate that:
1. QNAME in BAM = first whitespace token of FASTQ header (normalized).
2. No BAM QNAME contains a space (SAM spec violation).
3. BAM record count is within 10% of FASTQ read count (unmapped allowance).

This validator would catch the per-file gating bug at alignment time rather than
at consensus time (currently where _check_read_name_compatibility() fires).

---

## Data Provenance

- Synthetic FASTQ: make_synthetic_fastq.py -- 53-read safe + 3-read adversarial
- Reference: S288C_reference_sequence_R64-5-1_20240529.fsa (R64-5-1)
- Aligner versions: minimap2 2.28-r1209, mapPacBio 39.26 (BBMap), deSALT 1.5.6, gapmm2 25.4.5
- Production verification: Sherlock set1 wt_rep1 (v3_20260429 pipeline)
- Full characterization table: results.tsv (393 rows, 5 aligner keys)
