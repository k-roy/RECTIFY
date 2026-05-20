# AGENT_FIXES.md

Fast coordination log for active debugging sessions across M1 / H2 / Sherlock.
**Read this before touching pipeline code. Update it when you find a bug.**
Archive entries into CHANGELOG.md when the session wave is done.

---

## [2026-05-20] mapPacBio QNAME sanitizer — ONT Dorado FASTQs

**Status:** FIXED (`e8c8070`)

**Affects:** any dataset sequenced with ONT Dorado (mex67aa, wtaa, and any
future DRS/cDNA deposits from the Nanopore). Does NOT affect cDNA pipeline
FASTQs produced by `rectify trim-polya` (those already have bare UUID QNAMEs).

**Symptom 1 — before fix 1:**
```
samtools view exit code 1: [E::sam_parse1] query name too long
mapPacBio failed after ~14400s
```
Dorado embeds full run metadata in the FASTQ header
(`@<uuid> runid=... ch=... flow_cell_id=... basecall_gpu=...`).
mapPacBio.sh copies the full header verbatim into SAM QNAME (346 chars);
SAM spec limit is 254. `samtools view -bS` exits 1 on every read.

**Symptom 2 — introduced by fix 1 (`838293c`), fixed in `382fcc7`:**
```
BBMap AssertionError: Error in mpb_san.fastq, line N, sequence line is blank
BBMap terminated in an error state
```
`split(' ', 1)[0]` on a no-space header (e.g. `__mpbsplit_*` sub-reads)
retains the trailing `\n`, so appending `\n` produced a double-newline
(blank sequence line).

**Symptom 3 — introduced by fix 2 (`382fcc7`), fixed in `e8c8070`:**
```
BBMap validateQualityLength: quality string length != sequence length
This can be bypassed with the flag 'tossbrokenreads' or 'nullifybrokenquality'
```
Line-by-line `startswith('@')` mis-identified quality score lines beginning
with `@` (valid FASTQ — Phred Q31) as headers and truncated them.

**Root fix (`e8c8070`):** parse sanitized FASTQ in 4-line blocks; only the
header line (block line 1) is ever touched. Quality lines are written as-is
regardless of content.

**Safe to pull:** strict correctness improvement; any FASTQ that worked
before still works.

---

*To add an entry: symptom (exact error string), root cause (one sentence),
fix commit, safe-to-pull verdict.*
