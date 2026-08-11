# MICROHOM AUDIT V5 — signal-correctness — Auditor B

**Task:** POSITIONAL-SIGNAL + `_semiglobal_ed` CORRECTNESS.
Adversarially test that `_semiglobal_ed` and `_positional_signal` are correct and do
what the docstring claims. Verdict CLEAR/HOLD with reproducing cases.

**Mode:** READ-ONLY. Repo at worktree `agent-a25a2c1e784ad37dc`. Python
`/Users/kevinroy/miniconda3/bin/python`.

**Auditor:** B (independent of A — no coordination).

---

## Checkpoints

### CP0 — record created (harness not yet written)
- Read `junction_refiner.py` — target functions located:
  - `_semiglobal_ed(query, ref)` — lines 580-597
  - `_positional_signal(genome_seq, q, q_split, ne, new_je, W=28)` — lines 600-629
  - Wired via `drift_positional_gate` in `refine_read_junctions` (lines 915-930): a
    drift-flagged would-be-veto is SPARED when `psig >= gate` (gate default 0.0 = off).

### CP1 — semantics traced (source read; harness design fixed)
- **`_semiglobal_ed(query, ref)`** (lines 580-597): DP edit distance. `prev = range(n+1)`
  (row 0 = deletion of ref prefix is FREE for j, i.e. NOT — wait: row 0 is `list(range(n+1))`
  so cost of matching empty query against ref[:j] is j... NO. Re-examine below in CP2).
  Returns `min(prev)` = best over all ref-prefix END positions → FREE REF SUFFIX. Query is
  fully consumed (rows 1..m all present). Left edge `cur = [i] + ...` charges i for consuming
  i query chars against empty ref-prefix → HARD-ANCHORED at ref[0] (no free ref prefix). This
  matches the docstring: "all query to a PREFIX of ref, free ref suffix, hard-anchored ref[0]".
- **`q`, `q_split`**: `q = read.query_sequence` (pysam = ALWAYS forward-genome orientation,
  reverse-complemented for is_reverse reads). `q_split` = # query bases consumed before the
  N-op, counting aligned regions + adjacent (trailing) soft-clips (from `_iter_n_ops`, lines
  366-401). So `q[q_split:]` = the forward-genome DOWNSTREAM exon segment (3'-genome side of
  the intron), REGARDLESS of strand.
- **`ne`, `new_je`**: genome coords (0-based) of incumbent / moved intron END (right boundary).
  `genome[ne:]` = forward-genome downstream exon body. So `rescue = q[q_split:]` vs `genome[ne:]`
  is a FORWARD-GENOME-consistent comparison on BOTH strands. The docstring's "acceptor" framing
  is biologically loose for minus strand (there `ne` = biological DONOR), but the edit-distance
  comparison is COORDINATE-CONSISTENT because q, CIGAR, ne all live in forward-genome frame.
  → MINUS-STRAND geometry appears CORRECT despite loose docstring wording. TO VERIFY empirically.
- **The "free-k soft-clip escape"** the CLOSE claims to remove: `_score_junction`
  (junction_scoring.py:983) loops `k in [0,L)`, aligns `rescue[k:]` to moved exon2 — the first
  k rescue bases are SKIPPED FOR FREE. So a wrong move whose exon2 matches a few bases in still
  scores ~0. `_semiglobal_ed` hard-anchors query[0] too (no free-k) → forces leading rescue
  bases to align → exposes discriminating exon2 bases. Docstring rationale is STRUCTURALLY sound.

**Harness plan (adversarial probes):**
- (a) `_semiglobal_ed` edge cases: empty q/ref, q longer than ref, all-mismatch, free-suffix,
  indel vs subst, hard-anchor at ref[0]. Compare vs an independent brute-force / known values.
- (b) `_positional_signal` sign correctness: build reads where it gives WRONG sign
  (false-positive spare of fab drift; false-negative veto of real cryptic). Window W truncation,
  rescue off genome end, ties (=0), None cases.
- (c) MINUS-STRAND: construct is_reverse read, verify exon2 geometry (forward q vs genome[ne:]).
- (d) acceptor-centric: both-boundary move where ignoring donor causes a wrong spare.

Next: write harness to scratch, run edge-case battery, persist numbers.

### CP2 — PROBE (a) DONE: `_semiglobal_ed` primitive is CORRECT
Harness: `scratchpad/audit_v5/signal-correctness_B/harness.py`; output `out_a.txt`.
- **A1 random cross-check: 4000/4000 agree** with an independent brute-force
  "min over all ref-prefixes of full Levenshtein". No mismatch across m∈[0,12], n∈[0,18].
- A2 explicit edges all OK: empty q→0, empty ref→len(q), exact-prefix→0 (free suffix),
  all-mismatch len4→4, q-longer-all-mismatch(6,3)→6, q-longer-exact-then-over→2,
  single subst→1, 1 indel→1.
- **A3 HARD-ANCHOR at ref[0] CONFIRMED**: prepend non-matching base → ed 0→1 (a free-prefix
  "fitting" distance would stay 0). No free query/ref prefix.
- **A4 FREE-SUFFIX CONFIRMED**: appending arbitrary bases to ref leaves ed unchanged (0=0=0).
- A5 m>n lower bound ed≥m−n: 2000/2000 hold.
- **VERDICT for the primitive: matches docstring exactly. No fault in `_semiglobal_ed`.**

### CP3 — PROBES (b)(c)(d) DONE: unit-level results
Outputs `out_b.txt`, `out_c.txt`, `out_d.txt`.

**(b) sign correctness — NORMAL case CORRECT, TRUNCATION cases WRONG at unit level:**
- B1 None cases OK (new_je==ne → None; empty rescue → None).
- B2 TRUE POSITIVE (read==moved exon2): sig=+18 → spare correct. OK.
- B3 TRUE NEGATIVE (read==incumbent exon2): sig=−16 → no spare. OK.
- B4 ambiguous read: sig=0 → not spared at gate≥1 (conservative). OK.
- **B5 GENOME-END TRUNCATION (advisor candidate 1): sig=+2 → WRONG SPARE.** Incumbent
  `ne=254` sits 4 bp from contig end (len 258) so `ref_inc` truncates to 4 bases; the read
  truly matches the incumbent (its available bases) but `ed_inc=20` (forced insertions from
  the missing genome) > `ed_mov=18` → positive → spurious spare. Deterministic unit fault.
- **B6 RESCUE OFF MOVED-END: sig=0 → FALSE VETO** of a real cryptic whose moved exon2 is
  truncated at the contig end (`new_je=70`, len 74). Deterministic unit fault.
- B7 W-truncation: sig=0 when the ONLY discriminating base is past rescue position 28 (W=28).
  Signal is blind past W → falls through to margin/cap (conservative, not a wrong sign).

**(c) MINUS STRAND — geometry CORRECT (verified with a real pysam is_reverse read):**
- q_split=10=len(exonA); q[q_split:]==exonB (forward-genome downstream exon); ne=160 correct.
- Read-matches-incumbent → sig=−6 (veto); read-matches-moved → sig=+6 (spare). Both signs
  CORRECT. The forward-genome frame (pysam q, CIGAR, ne all forward) makes the ed comparison
  strand-consistent. Docstring "acceptor" wording is biologically loose on minus strand (ne is
  the DONOR there) but produces no sign error. NO minus-strand signal bug.

**(d) both-boundary donor blind spot (advisor candidate 2): CONFIRMED at unit level.**
- Acceptor 160→155, donor 110→105 (both ≤50). Signal AS CODED (q_split=110): sig=+5 → SPARE.
  The signal has NO donor term and uses the incumbent q_split, so a positive acceptor-ed spares
  the WHOLE both-boundary move including an unsupported donor shift. Matches the docstring's
  own admission (it argues a donor term would hurt genuine cryptics) — so this is BY DESIGN, but
  the design leaves a real hole IF both-boundary moves reach the gate.

**Next (DECISIVE): reachability.** A unit-level wrong sign only merits HOLD if the pipeline can
reach it. Must determine:
  (R1) can `ne`/`new_je` land within W+6 of a contig end AND survive candidate filters? (B5/B6)
  (R2) can `refine_read_junctions` emit a BOTH-boundary move that reaches the positional gate?
       (docstring claims donor-only ties keep incumbent — but both-boundary may win on acceptor.)
  (R3) is B5's truncation asymmetry (incumbent truncated, moved not) decision-flipping in situ?

### CP4 — DECISIVE reachability of B5/B6 genome-end truncation: NOT REACHABLE on real data
Data: bundled yeast fai (`…/S288C_reference_sequence_R64-5-1_20240529.fsa.gz.fai`),
`annotation.junc.bed` (402 lines / 370 usable), GFF introns (`intron`+`five_prime_UTR_intron`,
370). Metric: distance from each intron_END (`ne`, the coord `ref_inc = genome[ne:ne+W+6]`
extends RIGHT from) to the right contig end. Truncation needs (len − ne) < W+6 = 34 (+50 for
a candidate shift → 84).
- **min(len − intron_end) = 3929 bp** (chrVIII intron 558714–562643), IDENTICAL for
  annotation.junc.bed and the canonical GFF.
- **0 / 370 junctions within 34 bp; 0 / 370 within 84 bp.** Closest junction is ~46× the
  truncation threshold from any contig end. Left boundary min-to-left-end = 1897 bp.
- Smallest yeast contig = chrMito 85,779 bp; smallest nuclear = chrI 230,218 bp.
- **⇒ B5/B6 (genome-end truncation asymmetry) is a LATENT / theoretical fault: the wrong sign
  is deterministic at unit level but UNREACHABLE for real chromosome-scale references.** It could
  only manifest on a fragmented assembly / scaffold < ~100 bp, or a junction placed within ~84 bp
  of a contig break — not present in the shipped Scer reference, and the guard is DEFAULT-OFF +
  Scer-validated. NON-BLOCKING (worth a defensive clamp: skip/return None when the ref window
  truncates below rescue length).

### CP5 — (R2) both-boundary + interior drift reachability (corrected construction)
First construction (reach2.py) FAILED to reach the gate — advisor caught the bug: I overwrote
genome[160:184] with poly-A AFTER writing the moved exon2 at [155:179], collapsing the two
windows to the same string → both candidates scored 10.000 → tie → `is_alt` keeps incumbent →
`moves=False` → gate never reached; `move_microhom=0.000` (not even drift-flagged). This collision
is INTRINSIC to small-δ drift: genome[new_je:] and genome[ne:] share everything past max(ne,new_je),
so for the small shifts the guard targets the two references are nearly identical and the signal is
near zero BY GEOMETRY. Rebuilding cleanly (no overwrite) below to confirm the advisor's prediction:
interior fabricated drift gets NO false spare from TWO independent mechanisms (scorer free-k ties →
move not selected; hard-anchored psig<0 → not spared even if selected).

