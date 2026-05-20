# RECTIFY QNAME / Tag Edge-Case Punch List

Audit date: 2026-05-20  
Branch: `drs-validation-rebuild`  
Scope: cross-aligner consistency and edge cases beyond the mapPacBio QNAME
sanitizer chain already documented in AGENT_FIXES.md.  
Method: read-only code inspection (no code was modified).

---

## What I checked but found CLEAN

- **uLTRA / deSALT sort order**: both emit coordinate-sorted BAMs, but
  `_ensure_name_sorted` in `multi_aligner.py` inspects the BAM header SO tag
  and re-sorts to name order before the K-way merge. No reads are skipped.
- **deSALT deduplication**: `_dedup_desalt_bam` keyed on
  `(query_name, flag, chrom, pos, cigar)` correctly handles exact-duplicate
  primaries deSALT is known to emit.
- **gapmm2 PAF→BAM sequence injection**: `_paf_to_bam` in `multi_aligner.py`
  correctly reverse-complements sequences when `strand == '-'` (lines 1180–
  1185). The fallback `_restore_sequence_from_aligner_reads` is only reached
  if the FASTQ dict lookup misses, which should not occur in normal runs.
- **bam_writer.py lookup key**: `write_corrected_bam` uses raw `read.query_name`
  as the lookup key. Because per-aligner BAMs are auto-sanitized BEFORE
  correction runs, QNAMEs are already in normalized form when this code
  executes. Safe for fresh runs.
- **bam_processor.py pt-tag handling**: guarded with `try/except KeyError`
  (lines 496–499); not aligner-specific and will not crash on any aligner.
- **splice_aware_5prime.py**: produces only TSV corrections; all BAM writes
  go through `bam_writer.py`. No direct QNAME manipulation. Clean.
- **Empty-BAM validator behavior**: `is_clean()` returning True when
  `n_checked == 0` (`qname_validator.py:92–102`) is intentional — deSALT
  crash-recovery BAMs are empty and should not fail validation.
- **`_natural_sort_key` in K-way merge**: applied to all aligner names
  uniformly; the merge key is normalized before comparison at line 295.

---

## Edge Cases (sorted by severity)

### 1. [SILENT DATA CORRUPTION] Chimeric write path skips QNAME normalization

**File:** `rectify/core/consensus/consensus.py:459`  
**Also:** `rectify/core/consensus/chimeric_consensus.py:737`

**Issue:** The chimeric write path (`use_chimeric=True`) writes `out_read`
without normalizing its QNAME. The non-chimeric path at line 489 correctly
calls `_normalize_bam_read_name`, but line 459 omits this step.

**Evidence:**
```python
# consensus.py line 459 — chimeric path
out_bam.write(out_read)          # no normalization

# consensus.py line 489 — non-chimeric path (correct)
best_read.query_name = _normalize_bam_read_name(best_read.query_name or '')
```

The template read fed to `build_chimeric_read` is ANY aligner's record
whose sequence length matches the chimeric CIGAR (lines 401–413 in
`consensus.py`); `chimeric_consensus.py:737` copies the template's raw
`query_name` with no strip:
```python
out.query_name = template_read.query_name  # line 737, raw copy
```

**Impact:** If mapPacBio wins the chimeric template selection, the output
BAM contains `UUID_pt:i:N` QNAMEs. All downstream tools that key on
QNAME — `write_corrected_bam`, `merge_corrected_tsvs`, `samtools markdup` —
see a different ID than the non-chimeric consensus BAM for the same read.
The per-read correction is silently dropped.

**Suggested fix:** Add immediately after `build_chimeric_read` returns:
```python
out_read.query_name = _normalize_bam_read_name(out_read.query_name or '')
```
Mirror the pattern at non-chimeric line 489.

---

### 2. [SILENT DATA CORRUPTION, conditional] TSV-side normalization is narrower than BAM-side normalization

**File:** `rectify/core/consensus/corrected_consensus.py:48` (`_bare_uuid`)  
**File:** `rectify/core/consensus/corrected_consensus.py:518` (`_normalize_read_id`)  
**Compared against:** `rectify/core/consensus/consensus.py:171–183` (`_normalize_bam_read_name`)

**Issue:** `_bare_uuid` (used as the BAM-side join key in
`_read_hp_edit_distances`) and `_normalize_read_id` (applied to the TSV
`read_id` column in `_load_tsv`) both strip ONLY the `_pt:i:N` suffix.
`_normalize_bam_read_name` additionally strips BBmap `_<N>_length=<N>`,
generic SAM aux leaks (`_XA:Z:`, `_XC:i:`, etc.), and Dorado metadata keys.

**Evidence:**
```python
# corrected_consensus.py:48 — only pt:i:N
def _bare_uuid(name: str) -> str:
    if "_pt:i:" in name:
        return name.split("_pt:i:")[0]
    ...

# consensus.py:171 — comprehensive
_UNDERSCORE_COMMENT_RE = re.compile(
    r'(?:_pt:i:\d+|_\d+_length=\d+|_[A-Za-z][A-Za-z0-9]:[AZifHB]:'
    r'|_(?:runid|ch|start_time|sampleid|flow_cell_id|read|parent_read_id'
    r'|model_version_id|barcode)=).*$'
)
```

**Impact:** For runs with `RECTIFY_SKIP_POST_ALIGN_VALIDATION=1`, old TSVs
from pre-validator pipeline runs, or future aligners that emit BBmap-style or
Dorado-style comments, the HP-ED BAM→TSV join and `merge_corrected_tsvs`
cross-aligner join will silently drop rows for any read whose QNAME contains
one of the non-pt:i:N suffixes. The read is not corrected and not flagged.

**Suggested fix:** Import and reuse `_normalize_bam_read_name` from
`consensus.py` in `corrected_consensus.py` for both `_bare_uuid` and
`_normalize_read_id`. Remove the two narrow implementations and replace with
calls to the canonical normalizer.

---

### 3. [SILENT DATA LOSS, cDNA pipeline only] cDNA pipeline FASTQ tags lost for all non-minimap2 winners

**File:** `rectify/core/align/multi_aligner.py:274` (minimap2 `-y` flag)  
**File:** `rectify/core/align/multi_aligner.py` (`_sanitize_mpb_fastq`, `_clean_fastq`)

**Issue:** Only `run_minimap2` passes `-y` to propagate FASTQ comment tags
(XA, XC, XF, XU, XK, XI, XS, XR) into BAM aux records. `_sanitize_mpb_fastq`
strips all FASTQ comments before calling mapPacBio, and `_clean_fastq` strips
all comments before calling deSALT and gapmm2. uLTRA has no `-y` equivalent.

**Evidence:**
```python
# multi_aligner.py:274 — only minimap2
cmd = ['minimap2', '-y', '-ax', 'splice:hq', ...]
# all other wrappers: no -y; _clean_fastq / _sanitize_mpb_fastq strip tags
```

**Impact:** When mapPacBio, deSALT, gapmm2, or uLTRA wins consensus for a
cDNA read, the output BAM record has no XA/XC/XF/XU/XR aux tags. Any
downstream analysis that reads cDNA pipeline metadata from the consensus BAM
(tail length, cluster size, full-length tier, UMI, read IDs) silently gets
missing values. The data was never lost from the source FASTQ; it was simply
not propagated.

**Suggested fix (preferred):** After consensus selection, inject the cDNA tags
from the per-read FASTQ metadata into the winning BAM record, regardless of
which aligner won. This is a post-consensus tag-injection pass keyed on the
normalized QNAME → FASTQ comment map. Alternatively, add `-y` equivalents to
each aligner wrapper that supports tag passthrough, and strip only known-bad
tags rather than all FASTQ comments.

---

### 4. [COSMETIC / DETERMINISM] Duplicate-primary tiebreaker comment claims MAPQ; code uses only N-ops

**File:** `rectify/core/consensus/consensus.py:300–320`

**Issue:** The comment at lines 300–301 states the duplicate-primary
tiebreaker is "most N-ops (most complete splice chain), then highest MAPQ".
The code at lines 315–320 uses only `_n_ops` with `max()` and no secondary
MAPQ key.

**Evidence:**
```python
# comment (line 300–301):
# "prefer the alignment with the most N-ops (most complete splice chain),
#  then highest MAPQ."

# code (line 315–320):
group[aligner] = max(
    candidates,
    key=lambda r: _n_ops(r),   # MAPQ never consulted
)
```

**Impact:** When two mapPacBio primaries have identical N-op counts, Python's
`max()` returns the first one encountered (iterator order), which may vary
across Python versions or pysam sort stability. Reproducibility concern, not
a correctness failure. No silent data corruption.

**Suggested fix:** Either update the comment to remove the MAPQ claim, or
implement the tiebreaker: `key=lambda r: (_n_ops(r), r.mapping_quality)`.

---

### 5. [SILENT DATA LOSS, narrow] Sequence restore does not check orientation match

**File:** `rectify/core/consensus/consensus.py:362–373`

**Issue:** `_restore_sequence_from_aligner_reads` copies `query_sequence`
from the first donor whose sequence length matches the CIGAR query span
(line 366), without checking whether the donor's `is_reverse` flag matches
`best_read.is_reverse`.

**Evidence:**
```python
# consensus.py:362–373
for donor_read in aligner_reads.values():
    seq = donor_read.query_sequence
    if seq is None:
        continue
    if expected_len > 0 and len(seq) != expected_len:
        ...
        continue
    best_read.query_sequence = seq          # no strand check
    best_read.query_qualities = donor_read.query_qualities
    return
```

**Impact:** This function fires only when gapmm2 wins and its PAF→BAM record
has `query_sequence=None` AND the per-read FASTQ dict lookup misses. In that
narrow scenario, if the only same-length donor is mapped on the opposite
strand, the inserted sequence is the reverse-complement of the correct
sequence. The SEQ in the output BAM is wrong; all base-level correction
downstream (mismatch/indel calls) operates on a reversed string.

Under normal operation (gapmm2 FASTQ dict populated, other aligners unmapped
or length-matched on the same strand) this path is not reached. The risk is
elevated for reads with very short unique sequence lengths shared across
strands.

**Suggested fix:** Add a strand-match check before accepting a donor:
```python
if donor_read.is_reverse != best_read.is_reverse:
    continue
```

---

### 6. [DIAGNOSTIC GAP] Validator samples first N primary reads sequentially, not randomly

**File:** `rectify/core/align/qname_validator.py:131–164`

**Issue:** `_scan_bam` iterates the BAM from position 0 and breaks after
`sample_size` (default 200) primary reads. If aligner-specific QNAME mutations
cluster in the tail of the BAM (e.g., a secondary pass or a chromosome that
sorts last), the sample misses them and reports `is_clean() == True`.

**Evidence:**
```python
# qname_validator.py:162
if r.n_checked >= sample_size:
    break
```

The validator then skips auto-repair and the per-aligner BAM is passed to
consensus with mutated QNAMEs in the unseen tail.

**Impact:** Known mutations (mapPacBio `_pt:i:N`) are uniform across the BAM,
so the sequential sample catches them reliably. The risk is for hypothetical
future aligner behaviors where mutations are position-dependent. Low severity
for the current aligner panel; higher if new aligners are added.

**Suggested fix:** For robustness, either: (a) sample at evenly-spaced offsets
by pre-scanning record count and choosing systematic indices, or (b) use a
coordinate-sort-agnostic approach: hash the QNAME and sample every Nth
primary read (`if hash(qn) % stride == 0`). Alternatively, document the
known limitation in the module docstring.

---

## What I couldn't check

- **Runtime behavior of chimeric path in production BAMs**: the chimeric
  consensus path (`use_chimeric=True`) is an opt-in flag. I could not
  determine from the codebase alone how frequently it is invoked in the
  Chanfreau lab's standard `rectify correct` runs, or whether mapPacBio
  is the template winner often enough for issue #1 to have already
  corrupted existing output BAMs. A one-shot check (`samtools view
  consensus.bam | awk '{print $1}' | grep '_pt:i:'`) on any set1/2/3
  consensus BAM produced with `--chimeric` would confirm or rule out
  exposure.

- **BBmap wrapper in the active aligner panel**: `_normalize_bam_read_name`
  covers BBmap's `_<N>_length=<N>` suffix. I did not find a `run_bbmap`
  wrapper in `multi_aligner.py` for the active panel (minimap2, mapPacBio,
  gapmm2, deSALT, uLTRA). If BBmap is not currently wired in, issues #2
  and #3 have no BBmap exposure. Confirm with `grep -r 'bbmap\|bwa'
  rectify/core/align/`.

- **`_fallback_simple_selection` in chimeric_consensus.py**: line 1051 sets
  `read_id = best_read.query_name` (unnormalized). This fallback fires when
  no segment-level chimeric assembly succeeds. I did not trace all downstream
  consumers of this `read_id` field to confirm whether normalization is needed
  there too, or whether the field is only used for logging.

- **Provenance / SIDECAR audit**: this audit focused on QNAME and tag
  consistency. It did not inspect whether consensus BAM outputs are
  accompanied by PROVENANCE.json sidecars per the lab standard
  (`feedback_provenance_alongside_outputs.md`).
