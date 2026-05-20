# Sonnet Briefing — `read_num` + sidecar architecture (RN aux tag hybrid)

**You are a Sonnet 4.x subagent.** Goal: implement Option B from the
sidecar design discussion (RN auxiliary tag hybrid). One PR scope.
Read this entire document before writing code.

**Working dir:** `/Users/kevinroy/work/rectify`
**Branch:** `drs-validation-rebuild` (NOT `master`)
**Cluster context:** M1 is the source of truth. Heavy work belongs on
Sherlock (`sherlock-sbatch` skill) or Hoffman2 (`h2-qsub` skill); the
default partition is `larsms` for ~0 queue wait.

---

## 0. Hard prereqs — verify before any code

Run each command. If any fails, HALT and report to Opus.

```bash
# (a) The QNAME hardening commit is at HEAD or ancestor
git log --oneline -5 | grep -q '45cbc13.*qname' && echo OK

# (b) Cross-cluster heads match (the doc in HANDOFF.md may be stale —
#     trust git, not the doc).
cd /Users/kevinroy/work/rectify && git log --oneline -1
git ls-remote origin drs-validation-rebuild | awk '{print substr($1,1,7)}'
ssh hoffman2 'bash -lc "cd /u/home/k/kevinroy/software/rectify && git log --oneline -1"'
ssh sherlock  'bash --norc --noprofile -c "cd /oak/stanford/groups/larsms/Users/kevinroy/software/rectify && git log --oneline -1"'

# (c) The existing QNAME test suite passes (regression baseline)
python -m pytest tests/test_qname_sanitizer_and_validator.py \
                 tests/test_normalize_read_name.py \
                 tests/test_consensus_selection.py \
                 tests/test_corrected_consensus_tiebreaker.py \
                 tests/test_chimeric_consensus.py \
                 tests/test_parallel_aligner_schedule.py \
                 tests/test_gapmm2_seq_restore.py -q
# Expected: 113 passed, 4 skipped, ~14s
```

---

## 1. Read-first list

Read these in order. Do NOT start coding until you've read all six.

| Doc | Why |
|---|---|
| `HANDOFF.md` (repo root) | Design discussion, A-vs-B analysis, the chosen recommendation (B). This briefing is the *execution* of the HANDOFF. |
| `scripts/diagnostics/qname_mutation_survey/EDGE_CASES.md` issue #3 | The motivating symptom (cDNA FASTQ tags lost for non-minimap2 winners). Your work makes this issue go away. |
| `scripts/diagnostics/qname_mutation_survey/REPORT.md` | Per-aligner tag-stripping behavior. Tells you exactly which wrappers need a post-alignment RN injection pass and which carry RN natively. |
| `AGENT_FIXES.md` "[2026-05-20] QNAME pipeline hardening" | Context for what landed at `45cbc13` and what the existing defensive layer does (you'll reuse it as belt-and-suspenders). |
| `rectify/core/align/qname_validator.py` | The streaming-rewrite pattern at `_rewrite_bam_normalize_qnames` is the template for the RN injection pass — read it before writing yours. |
| `CLAUDE.md` (project root) | "Verify-don't-rewrite", "Surgical staging", "M1 memory discipline", cluster choice rules. Non-negotiable. |

---

## 2. Pre-recorded design decisions (do NOT re-litigate)

These were decided in the design session. The HANDOFF.md has the full
trade-off analysis. Your job is execution, not redesign.

### 2.1 Variant choice: B (RN aux tag hybrid), not A (full anonymization)

- Chunker keeps original QNAMEs in the FASTQ.
- Chunker adds an `RN:i:<num>` tag in the FASTQ-comment portion of
  the header (tab-separated, alongside any existing X-tags from the
  cDNA pipeline).
- Aligners that natively propagate FASTQ comments to BAM aux records
  (minimap2 with `-y`) carry RN automatically.
- Aligners whose wrappers strip FASTQ comments (mapPacBio, deSALT,
  gapmm2, uLTRA, bwa) get a small post-alignment RN injection pass.
- Sidecar parquet maps `read_num → original_qname + fastq_comment + chunk_id`.
- Consensus emission restores original X-tags onto the winning record
  via the sidecar (RN-keyed lookup, with QNAME fallback if RN is absent).
- Existing UUID-QNAME BAMs without RN remain valid — the K-way merge
  prefers RN when present, falls back to normalized QNAME otherwise.

### 2.2 Backwards compatibility — strict

- Every code path that currently consumes a per-aligner BAM must
  continue to work on RN-less BAMs. Test fixtures and production set1/
  set2/set3 BAMs must NOT need regeneration.
- When ANY aligner BAM in a K-way merge lacks RN, the merge falls
  through to the existing `_normalize_bam_read_name`-based path. No
  cliff edges; no required migration step.

### 2.3 RN integer space

- Per-sample globally unique. Per-chunk-prefix form: `<chunk_id>_<seq>`
  encoded as int64 by packing chunk into the high bits, sequence into
  the low bits. Or — simpler — use a global counter coordinated by the
  chunker. **Decide in §3.2** but document the choice in the sidecar
  schema docstring.

### 2.4 Sidecar location

- One parquet per sample, written by the chunker, named
  `<sample>.read_num_sidecar.parquet`. Lives alongside the per-aligner
  BAMs in the same output directory.
- A PROVENANCE.json sidecar (lab convention,
  `feedback_provenance_alongside_outputs.md` memory) accompanies the
  parquet with: chunker version, total reads, chunk count, sha256 of
  source FASTQ, sha256 of the parquet itself.

### 2.5 Scope of "tags carried in the sidecar"

- The full FASTQ comment string (everything after the first whitespace
  in the FASTQ header), stored verbatim. Restoration parses on demand.
- This covers the cDNA pipeline's `XA/XC/XF/XU/XR/XI/XK/XS` namespace
  without RECTIFY having to know each tag.
- Empty string for FASTQs with no comment (e.g., bare-UUID DRS).

---

## 3. Architecture overview (1 page)

```
                    ┌─────────────────────┐
   Original FASTQ ─►│   rectify split     │──► derived FASTQs (per chunk)
   (untouched)      │   + RN assignment   │     with @<qname>\tRN:i:<n>\t<orig-comment>
                    └──────────┬──────────┘
                               │
                               ▼
                    sample.read_num_sidecar.parquet
                    columns: read_num, original_qname,
                             fastq_comment, chunk_id,
                             seq_md5, qual_md5
                               │
   ┌───────────────────────────┴──────────────────────────────┐
   │                                                          │
   ▼                                                          ▼
minimap2 (-y carries RN natively)            mapPacBio/deSALT/gapmm2/uLTRA/bwa
   │                                         (wrappers strip comments;
   ▼                                          post-alignment RN-injection pass
per-aligner BAM with RN:i tag                 keyed on QNAME → RN dict
                                              built at chunk time)
                                              │
                                              ▼
                                              per-aligner BAM with RN:i tag
   │                                          │
   └─────────────┬────────────────────────────┘
                 │
                 ▼
        K-way merge in consensus._iter_name_grouped_bams
        - Key on RN:i when ALL inputs have it (fast int compare)
        - Fall back to _normalize_bam_read_name(QNAME) when any input lacks it

                 │
                 ▼
        Consensus.bam emission
        - For winning read, read RN (or QNAME)
        - Look up original_qname + fastq_comment in sidecar
        - Restore comment-derived tags onto output record
        - (Existing _normalize_bam_read_name on winner's QNAME still runs)
```

---

## 4. Sidecar parquet schema

Create `rectify/core/chunking/sidecar.py`. Schema:

| Column | Type | Notes |
|---|---|---|
| `read_num` | `int64` | Primary key. Globally unique within sample. |
| `original_qname` | `string` | Verbatim FASTQ QNAME (everything between `@` and first whitespace). |
| `fastq_comment` | `string` | Everything after the first whitespace in the FASTQ header, verbatim. Empty string if no comment. |
| `chunk_id` | `string` | E.g., `chunk_03_of_24`. For traceability + per-chunk parallelism. |
| `seq_md5` | `binary(16)` | MD5 of the FASTQ sequence line. Verification only. |
| `qual_md5` | `binary(16)` | MD5 of the FASTQ quality line. Verification only. |

Row-group size: 100k. Compression: zstd level 3. Index column:
`read_num`.

A reader API:

```python
class ReadNumSidecar:
    @classmethod
    def open(cls, path: str | Path) -> "ReadNumSidecar": ...
    def lookup(self, read_num: int) -> SidecarRow: ...  # raises KeyError
    def lookup_many(self, read_nums: Iterable[int]) -> dict[int, SidecarRow]: ...
    def lookup_by_qname(self, qname: str) -> SidecarRow: ...  # fallback path
    def __len__(self) -> int: ...

@dataclass(frozen=True)
class SidecarRow:
    read_num: int
    original_qname: str
    fastq_comment: str
    chunk_id: str
    seq_md5: bytes
    qual_md5: bytes
```

The writer API (called from the chunker):

```python
class ReadNumSidecarWriter:
    def __init__(self, path: str | Path, sample_id: str): ...
    def add(self, read_num: int, original_qname: str,
            fastq_comment: str, chunk_id: str,
            seq: str, qual: str) -> None: ...
    def close(self) -> None: ...
    def __enter__(self): ...
    def __exit__(self, *args): ...
```

---

## 5. Implementation order (phases with acceptance criteria)

Each phase ends with a green pytest run BEFORE moving to the next.

### Phase 1 — Sidecar module + tests (no integration yet)

**Files to create:**
- `rectify/core/chunking/__init__.py` (if not present)
- `rectify/core/chunking/sidecar.py` — schema, writer, reader
- `tests/test_read_num_sidecar.py`

**Tests:**
- Round-trip: writer adds 100k rows, reader retrieves them all by
  `read_num` and by `original_qname`.
- Empty comment handled.
- Fingerprint mismatch detection (provide a row with the wrong seq;
  `verify` should report it).
- Large-batch lookup: `lookup_many(range(0, 100_000))` returns a dict
  of length 100k.
- Streaming write does not require holding all rows in memory
  (write 1M rows, peak RSS < 200 MB).

**Acceptance:** `pytest tests/test_read_num_sidecar.py -q` passes.

### Phase 2 — Chunker integration

**Files to modify:**
- `rectify/core/commands/split_command.py` — invoke the sidecar writer
  during chunking; rewrite each FASTQ header to inject the RN tag.
- `rectify/core/align/multi_aligner.py: extract_fastq_chunk` (if it's
  the per-chunk FASTQ writer in the mapPacBio chunk path — verify with
  `grep -n extract_fastq_chunk rectify/`) — same treatment.
- Possibly `rectify/core/commands/run_all_command.py` if it has its own
  chunking — verify with `grep -n "def.*chunk\|split" rectify/core/commands/`.

**Tests:**
- New `tests/test_chunker_emits_sidecar.py`:
  - Chunk a synthetic 100-read FASTQ; assert the sidecar parquet has
    100 rows.
  - Assert the derived chunk FASTQs have `RN:i:<num>` in the header.
  - Assert reading the chunk FASTQ + the sidecar reconstructs the
    original input (round-trip via fingerprint check).
  - Assert chunking a FASTQ with existing X-tag comments
    (`@uuid\tXA:Z:foo\tXC:i:42`) preserves them AND adds the RN tag.

**Acceptance:** the new chunker test passes; existing chunker tests
still pass; the production cDNA pipeline chunker (if separate) still
runs end-to-end on the validation fixture under
`rectify/data/validation/`.

### Phase 3 — Aligner wrapper RN propagation

**Files to modify:** `rectify/core/align/multi_aligner.py`. For each
wrapper, decide which of two paths applies:

| Aligner | Path | Reason |
|---|---|---|
| minimap2 | Native: `-y` is already present; RN flows automatically. | No code change beyond ensuring `-y` is set. |
| mapPacBio | Inject post-alignment | `_sanitize_mpb_fastq` STRIPS comments. Capture `qname → rn` dict during sanitize; injection pass after BAM write. |
| deSALT, gapmm2 | Inject post-alignment | `_clean_fastq` STRIPS comments. Capture dict before; inject after. |
| uLTRA, bwa | Inject post-alignment | No comment-passthrough mechanism. Same as above. |

**Injection-pass pattern** (template — adapt for each wrapper):

```python
def _inject_rn_into_bam(bam_path: str, qname_to_rn: Mapping[str, int]) -> int:
    """Streaming pysam rewrite that adds RN:i:<num> to every record.

    Pattern mirrors qname_validator._rewrite_bam_normalize_qnames.
    Returns count of records that got an RN tag (records whose QNAME
    is not in the dict are left untagged — should be 0 in normal runs).
    """
    import pysam
    src = Path(bam_path)
    tmp = src.with_suffix('.rn_injected.tmp.bam')
    n_tagged = 0
    n_missing = 0
    with pysam.AlignmentFile(str(src), 'rb') as bin_, \
         pysam.AlignmentFile(str(tmp), 'wb', header=bin_.header) as bout:
        for read in bin_:
            rn = qname_to_rn.get(read.query_name or '')
            if rn is not None:
                read.set_tag('RN', rn, value_type='i')
                n_tagged += 1
            else:
                n_missing += 1
            bout.write(read)
    tmp.replace(src)
    if n_missing > 0:
        logger.warning(
            "[%s] RN injection: %d/%d records had no QNAME→RN mapping",
            bam_path, n_missing, n_tagged + n_missing,
        )
    return n_tagged
```

**Wiring:** for each comment-stripping wrapper, before stripping, build
the `qname_to_rn` dict by parsing the FASTQ comment for `RN:i:<num>`.
After the wrapper's existing post-alignment work completes (including
today's `validate_post_alignment_qnames` call), call
`_inject_rn_into_bam` with the dict.

**Tests:**
- New `tests/test_rn_propagation.py`:
  - Build a synthetic FASTQ with RN tags.
  - Mock or use the real minimap2 wrapper on M1; assert every BAM
    record has `RN:i`.
  - Mock the mapPacBio path's `_inject_rn_into_bam` call; assert
    the dict is built correctly and every record gets tagged.

**Acceptance:** `tests/test_rn_propagation.py` passes; existing
`tests/test_qname_sanitizer_and_validator.py` + `test_normalize_read_name.py`
still pass.

### Phase 4 — Consensus merge: prefer RN, fall back to QNAME

**Files to modify:** `rectify/core/consensus/consensus.py`
- `_iter_name_grouped_bams`: detect at iterator-construction time
  whether ALL input BAMs have RN on their primary records. If yes,
  switch the merge key to the integer RN. If no, use the existing
  `_normalize_bam_read_name(QNAME)` path.
- The "all have RN" detection: read the first non-secondary,
  non-supplementary, non-unmapped record from each input BAM, check for
  `RN` tag. Set a `use_rn_key: bool` and use a single code path with a
  `key_fn` injected.

**Tests:**
- Extend `tests/test_consensus_selection.py` with:
  - All-RN merge: 3 mock per-aligner BAMs with RN, K-way merge joins
    correctly by RN.
  - Mixed: 2 BAMs have RN, 1 doesn't → falls back to QNAME path.
  - All-QNAME (existing behavior): unchanged.

**Acceptance:** consensus tests pass; the bare-fallback path is
unchanged for existing BAMs.

### Phase 5 — Consensus emission: restore original tags

**Files to modify:** `rectify/core/consensus/consensus.py` (near the
existing winner normalize at line 489 and the chimeric mirror at
line 459).

**Behavior:** if a sidecar is present (passed via a new optional
parameter to the consensus orchestrator), then for each winning read:
1. Read its `RN` tag (or fall back to QNAME normalize → sidecar QNAME index).
2. Look up the original `fastq_comment` in the sidecar.
3. Parse the comment for SAM-format tags (`XA:Z:foo\tXC:i:42\t...`)
   and set them on the output record IF not already present (don't
   overwrite an existing tag the aligner produced — minimap2-winning
   records will already have them via `-y`).

**Surface area:** the consensus command (`rectify correct` or whatever
emits the consensus.bam) needs a `--read-num-sidecar <path>` flag.
Default = autodetect (look for `<sample>.read_num_sidecar.parquet`
alongside the input BAMs). Backwards compatible: no sidecar = no
restoration pass.

**Tests:**
- `tests/test_consensus_tag_restoration.py`:
  - Run a mock consensus with a sidecar; assert non-minimap2 winners
    have their XA/XC/XF/XU/XR tags restored.
  - Run without a sidecar; assert no restoration happens and output
    is identical to today's behavior.
  - Run with a sidecar that's missing the winning read's row; assert
    a warning is logged and the read is written without tags.

**Acceptance:** new tag-restoration tests pass; existing consensus
tests pass.

### Phase 6 — End-to-end smoke

Run the validation fixture (under `rectify/data/validation/`) through
the full pipeline with and without the sidecar. Compare outputs:
- With sidecar: every consensus.bam record should have RN + the cDNA
  X-tags from the source FASTQ.
- Without sidecar: behavior identical to HEAD pre-this-PR.

If the validation fixture's source FASTQ doesn't currently have X-tags
in its headers (it likely doesn't), construct a small synthetic
parallel — a 100-read FASTQ with `\tXA:Z:test\tXC:i:1` on each header —
and run the pipeline on it.

**Acceptance:** consensus.bam records match expected tag sets;
existing fixtures regenerate identically when the sidecar is absent.

### Phase 7 — Docs + AGENT_FIXES entry

- Update `AGENT_FIXES.md` with a new entry documenting the
  architecture, what it solves, the backwards-compat story, and any
  open issues.
- Update `HANDOFF.md` to note that the sidecar work is now implemented
  (or replace HANDOFF.md with a new handoff for whatever's next).
- Add a short user-facing section to docs:
  - When and why to use `--read-num-sidecar`.
  - The "I have an old BAM without RN, will it still work?" FAQ.

---

## 6. Operational guardrails

- **M1 memory discipline:** never >2 concurrent rectify/samtools/pandas
  subprocesses. Default heavy iteration to Sherlock via the
  `sherlock-sbatch` skill (per CLAUDE.md). Sidecar reads of 10M rows
  are pandas/parquet — keep an eye on RSS.
- **Surgical staging:** always `git add <explicit paths>`. Never
  `git add -A` or `git add .` — Kevin keeps WIP in the tree.
- **Verify-don't-rewrite:** every claim about existing behavior gets
  a tool call (grep, pytest, samtools view). Don't extrapolate from
  comments. The QNAME hardening commit `45cbc13` shipped after exactly
  this discipline; the validator's regression-test density is the
  template.
- **No new aligner binaries.** This work touches only Python; the
  aligners are unchanged.
- **No commits without explicit user instruction.** Land changes in
  the working tree, run tests, report. Kevin will commit when ready.
- **AGENT_FIXES.md check:** if you hit a crash that's not obviously
  yours, `grep -i <symptom> AGENT_FIXES.md` before spending >10 min
  diagnosing.

---

## 7. Out of scope (do NOT do these)

- Full QNAME anonymization (Option A). Decision was B.
- Migration of existing UUID-QNAME BAMs on disk. Backwards
  compatibility is a hard requirement; no migration needed.
- Touching today's defensive layer (`qname_validator.py`,
  `_sanitize_mpb_fastq`, `_normalize_bam_read_name`). It stays as
  belt-and-suspenders; the sidecar work goes ON TOP of it.
- Fixing the pre-existing Commit B manifest-split breakage in
  `test_quantseq_rev_integration.py` / `test_validation_reads.py`.
  Out of scope for the sidecar PR.
- Restructuring the cDNA pipeline's own internal sidecar mechanisms
  (`project_cdna_pipeline.md` memory). If the cDNA pipeline has its
  own chunker, integrate with it; don't replace it.
- Removing `_normalize_bam_read_name` from the merge path. It's the
  fallback for RN-less BAMs.

---

## 8. Definition of Done

- All seven phases complete with green tests.
- The full QNAME suite from §0 (113 tests) still passes.
- New tests: `test_read_num_sidecar.py`, `test_chunker_emits_sidecar.py`,
  `test_rn_propagation.py`, `test_consensus_tag_restoration.py` — all
  passing.
- End-to-end smoke (phase 6) confirms: with sidecar, every consensus
  record has the cDNA X-tags regardless of which aligner won;
  without sidecar, output is byte-identical to today.
- AGENT_FIXES.md entry added documenting the new architecture.
- Working-tree diff is surgical (touches only the files listed in
  the phases; no spurious whitespace changes, no `git add -A`
  contamination).
- A short return message to Opus summarizing: phase status, test
  counts, anything that surprised you, and any deferred items. Under
  400 words.

---

## 9. Known pitfalls (from this session's design discussion)

- **Sequence/quality mutation by aligners.** gapmm2 PAF→BAM emits
  no SEQ; deSALT hard-clips; `is_reverse=True` records carry RC'd
  `query_sequence`. The seq/qual fingerprint check in §4 is for
  *verification only*, not the primary key. Don't write code that
  fails on a fingerprint mismatch — log it and continue. The primary
  key is `read_num` (integer, immutable).
- **PCR duplicates.** Two reads with identical (seq, qual) are
  routine in PCR-cDNA libraries. They have distinct `read_num`s and
  are distinguishable. Tests should cover this case explicitly.
- **Tab-aux interaction with mapPacBio.** Today's `_sanitize_mpb_fastq`
  strips tab separators. Your RN-capture pass needs to read the tab
  comment BEFORE the sanitize step. Pattern: do the read in the same
  loop as the sanitize, building the dict and writing the cleaned
  FASTQ in one pass.
- **The validator runs AFTER the aligner wrapper.** Today's
  `validate_post_alignment_qnames` does its work post-alignment. Your
  RN injection happens AFTER the validator (or simultaneously with
  the streaming rewrite). Be careful about ordering: the validator
  might rewrite QNAMEs, which would invalidate your `qname → rn`
  dict if you built it from the *original* FASTQ. Solution: build the
  dict against the SANITIZED QNAMEs (the form they'll have in the
  BAM post-validator).
- **uLTRA/bwa wrappers are less tested.** Today's validator was
  wired into them but they aren't in the default 4-aligner panel for
  DRS. Don't over-invest in their RN injection until you have a real
  use case.

---

## 10. If you get stuck

Halt and report to Opus with:
- What you've done so far (commit hashes if any, files staged if not).
- The exact error or contradiction.
- What you were about to do next, and why you stopped.

Do NOT discard work in progress. Do NOT take destructive shortcuts.
The QNAME defensive layer at `45cbc13` is the safety net — even if the
sidecar work is incomplete, the pipeline still works correctly because
the existing validator + normalizer + auto-sanitize are in place.
