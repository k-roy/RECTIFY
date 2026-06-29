#!/usr/bin/env python3
"""Error-realism INJECTOR — the 3 correlation layers real ONT has that pbsim3 /
a marginal error table LACK (measured 2026-06-27; see SPEC §ERROR-REALISM).

This is NOT a read simulator. It takes an ALREADY-CONSTRUCTED read — a clean
absolute-truth read from ``controlled.py`` OR a pbsim3 read — and corrupts it
with a generative error process whose MARGINAL rate matches the empirical table
but whose CORRELATION STRUCTURE matches real ONT on three measured axes:

  1. per-read OVER-DISPERSION  (hot-read mixture)            -> layer 1
  2. within-read BURST / clustering                          -> layer 2 (2-state HMM)
  3. fat-tailed multi-base INDEL-RUN length                  -> layer 3

Measured real-yeast-DRS vs rate-matched pbsim3 gaps (SPEC §ERROR-REALISM), the
calibration TARGETS this injector is fit against on the M1:
  * length-adj per-read dispersion index   real 9.95   vs pbsim 3.55
  * within-read sub-5bp-gap EXCESS          real 2.83x  vs pbsim ~1.2x
  * fraction of indel events >= 2bp         real 0.39   vs pbsim-DRS 0.19

================================ TRUTH RULE ================================
The injected errors are a SEPARATE per-read ERROR TRACK (``List[ErrorEvent]``),
NOT ``IndelTruth``. Recording stochastic background errors as ``IndelTruth``
re-triggers the scorer's per-read ``has_unexplained`` gate (the SESSION-2 pbsim
lesson — scattered MAF errors zeroed exact-indel concordance) and would break
the (C)/(D) Tier-1 assertions. The error track is ALSO exactly PI-#2's per-read
"hotness" label (exonic error density) for the novel-junction-FDR test.

=========================== CALIBRATION STATUS ============================
The realistic params produced by ``calibrate_params`` are PLACEHOLDER. They are
fit by a numerical (coordinate-descent) loop to the MEASURED real-yeast targets
above — but those targets are an UPPER BOUND (read-vs-reference conflates RNA-mod
basecall error + minimap2 pile-up near indels/junctions, an ALIGNMENT artifact a
pre-alignment injector must NOT reproduce; SPEC). MAGNITUDE is therefore UNKNOWN
until calibrated against SIRV/LRGASP ABSOLUTE truth on Sherlock.

The M1 ``self_check`` is a COMPOSITION / INTERACTION check ONLY:
  * NULL params reproduce the marginal / Poisson baseline (dispersion ~1, gap
    excess ~1, indel runs ~single-base) — backward-compat;
  * realistic params move each statistic in the right direction AND the layers
    compose without one zeroing another.
It is NOT a realism validation. Realism is established only by SIRV magnitude
calibration AND a DISTRIBUTIONAL comparison (multiple quantiles / full gap-size
and run-length shapes, fit-some / hold-out-others) — both Sherlock-gated.

Author: Kevin R. Roy
"""
from __future__ import annotations

import math
import random
from dataclasses import dataclass, field, replace
from typing import Dict, List, Optional, Sequence, Tuple

BASES = "ACGT"


# ---------------------------------------------------------------------------
# The per-read ERROR TRACK (NOT IndelTruth — see module docstring TRUTH RULE)
# ---------------------------------------------------------------------------
@dataclass
class ErrorEvent:
    """One contiguous background error EVENT in CLEAN-sequence (ref) coordinates.

    A maximal contiguous error run = ONE event with a ``length`` — this is the
    canonical unit for BOTH the dispersion (count of events) and the gap
    (distance between event starts) statistics, so multi-base indel runs (layer
    3) do NOT inflate the sub-5bp-gap count (they are one event, and the next
    event starts >= length away). This is the run-grouping the empirical
    profiler's ``_group_pos_into_runs`` produces on a BAM, so the same statistic
    is computable from real alignments.
    """
    pos: int        # 0-based start in the CLEAN (pre-injection) sequence
    kind: str       # 'sub' | 'ins' | 'del'
    length: int     # bases substituted / inserted / deleted (run length, >= 1)
    bases: str      # sub/ins: introduced bases; del: the deleted reference bases


# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
@dataclass
class InjectorParams:
    """Generative-process parameters. ``NULL`` reproduces the marginal baseline.

    The marginal is the per-position event-INITIATION hazard ``base_rate``; the
    calibration loop tunes it so the MEASURED per-base error rate hits the target
    (it does not try to preserve the marginal analytically across layers — the
    loop measures the combined output, which is the honest way given the layers
    interact; advisor note #2).
    """
    base_rate: float = 0.019           # per-position error-event initiation hazard

    # error TYPE split among events (from the empirical/marginal table)
    frac_sub: float = 0.55
    frac_ins: float = 0.15
    frac_del: float = 0.30

    # --- layer 1: per-read rate multiplier (over-dispersion) ---
    # 'const' (NULL), 'gamma' (E[m]=1, Var=1/shape), or 'mixture' (cold + hot tail)
    read_rate_dist: str = "const"
    gamma_shape: float = 1.0
    mix_hot_frac: float = 0.0          # 'mixture': P(read is hot)
    mix_hot_mult: float = 1.0          # 'mixture': hot-read multiplier (cold derived)

    # --- layer 2: 2-state burst HMM (cold/hot), modulates the LOCAL hazard ---
    burst_on: bool = False
    p_cold_to_hot: float = 0.0         # burst FREQUENCY knob
    p_hot_to_cold: float = 0.5         # 1/p = mean burst length (geometry, FIXED on M1)
    hot_factor: float = 1.0            # hazard multiplier inside a burst (cold derived)

    # --- layer 3: indel run-length pmf over lengths >= 1 (fat tail) ---
    indel_run_pmf: Tuple[Tuple[int, float], ...] = ((1, 1.0),)

    @classmethod
    def null(cls) -> "InjectorParams":
        """Backward-compat baseline: constant rate, no burst, single-base indels
        => marginal Poisson errors with ~single-base indel runs."""
        return cls(read_rate_dist="const", burst_on=False, indel_run_pmf=((1, 1.0),))

    # -- derived quantities --
    def _stationary_hot(self) -> float:
        a, b = self.p_cold_to_hot, self.p_hot_to_cold
        return 0.0 if (a + b) == 0 else a / (a + b)

    def cold_factor(self) -> float:
        """Cold-state hazard multiplier chosen so the stationary-weighted mean
        factor == 1 (the burst layer redistributes hazard, it does not add it)."""
        if not self.burst_on:
            return 1.0
        ph = self._stationary_hot()
        pc = 1.0 - ph
        if pc <= 0:
            return 1.0
        cf = (1.0 - ph * self.hot_factor) / pc
        return max(0.0, cf)


# ---------------------------------------------------------------------------
# The generative process
# ---------------------------------------------------------------------------
def _draw_read_multiplier(p: InjectorParams, rng: random.Random) -> float:
    if p.read_rate_dist == "const":
        return 1.0
    if p.read_rate_dist == "gamma":
        # E[m] = shape * scale = 1 with scale = 1/shape; Var = 1/shape.
        return rng.gammavariate(p.gamma_shape, 1.0 / p.gamma_shape)
    if p.read_rate_dist == "mixture":
        if rng.random() < p.mix_hot_frac:
            return p.mix_hot_mult
        # cold multiplier so the mean is 1: f*hot + (1-f)*cold = 1
        f = p.mix_hot_frac
        cold = (1.0 - f * p.mix_hot_mult) / (1.0 - f) if f < 1.0 else 1.0
        return max(0.0, cold)
    raise ValueError(f"unknown read_rate_dist {p.read_rate_dist!r}")


def _draw_indel_run(p: InjectorParams, rng: random.Random) -> int:
    r, acc = rng.random(), 0.0
    for length, prob in p.indel_run_pmf:
        acc += prob
        if r <= acc:
            return length
    return p.indel_run_pmf[-1][0]


def _draw_kind(p: InjectorParams, rng: random.Random) -> str:
    r = rng.random()
    if r < p.frac_sub:
        return "sub"
    if r < p.frac_sub + p.frac_ins:
        return "ins"
    return "del"


def _walk_events(n: int, params: InjectorParams, rng: random.Random,
                 clean_seq: Optional[str] = None
                 ) -> Tuple[Optional[List[str]], List[ErrorEvent]]:
    """Core generative walk over ``n`` clean positions.

    If ``clean_seq`` is given, also builds the dirty-sequence fragments (returned
    as a list to ``"".join``); if it is None this is the EVENTS-ONLY fast path
    (calibration / measurement) — no base generation, no string building.
    """
    mult = _draw_read_multiplier(params, rng)
    cold_f = params.cold_factor()
    state_hot = False  # HMM starts cold
    build = clean_seq is not None
    events: List[ErrorEvent] = []
    out: Optional[List[str]] = [] if build else None
    rnd = rng.random
    burst = params.burst_on
    p_ch, p_hc = params.p_cold_to_hot, params.p_hot_to_cold
    hot_f, br = params.hot_factor, params.base_rate
    i = 0
    while i < n:
        if burst:
            if state_hot:
                if rnd() < p_hc:
                    state_hot = False
            elif rnd() < p_ch:
                state_hot = True
            factor = hot_f if state_hot else cold_f
        else:
            factor = 1.0

        hazard = br * mult * factor
        if hazard > 1.0:
            hazard = 1.0

        if rnd() >= hazard:
            if build:
                out.append(clean_seq[i])
            i += 1
            continue

        kind = _draw_kind(params, rng)
        if kind == "sub":
            if build:
                alt = rng.choice([b for b in BASES if b != clean_seq[i]])
                out.append(alt)
            else:
                alt = ""
            events.append(ErrorEvent(pos=i, kind="sub", length=1, bases=alt))
            i += 1
        elif kind == "ins":
            L = _draw_indel_run(params, rng)
            if build:
                ins = "".join(rng.choice(BASES) for _ in range(L))
                out.append(ins)
                out.append(clean_seq[i])  # insertion does not consume the ref base
            else:
                ins = ""
            events.append(ErrorEvent(pos=i, kind="ins", length=L, bases=ins))
            i += 1
        else:  # del
            L = min(_draw_indel_run(params, rng), n - i)
            if L <= 0:
                if build:
                    out.append(clean_seq[i])
                i += 1
                continue
            deleted = clean_seq[i:i + L] if build else ""
            events.append(ErrorEvent(pos=i, kind="del", length=L, bases=deleted))
            i += L                         # deletion consumes L ref bases
    return out, events


def inject(clean_seq: str, params: InjectorParams, rng: random.Random
           ) -> Tuple[str, List[ErrorEvent]]:
    """Corrupt ``clean_seq`` and return ``(dirty_seq, events)``.

    Source-agnostic: ``clean_seq`` may be a controlled-construction read OR a
    pbsim3 read (advisor note #5 — the before/after-pbsim comparison needs this).
    ``events`` is the separate ERROR TRACK in CLEAN-sequence coordinates.
    """
    out, events = _walk_events(len(clean_seq), params, rng, clean_seq=clean_seq)
    return "".join(out or []), events


# ---------------------------------------------------------------------------
# Measurement — the SAME three statistics, source-agnostic (advisor note #5)
# ---------------------------------------------------------------------------
def _geometric_p_lt(r: float, k: int) -> float:
    """P(gap < k) for inter-event gaps under an iid hazard ``r`` per position.
    Gaps are differences of consecutive event START positions (>= 1), so the
    null is Geometric on {1,2,...}: P(gap=g) = (1-r)^(g-1) r => P(gap<k) =
    1-(1-r)^(k-1)."""
    if r <= 0:
        return 0.0
    return 1.0 - (1.0 - r) ** (k - 1)


def measure_error_structure(per_read_events: Sequence[Sequence[ErrorEvent]],
                            read_lengths: Sequence[int],
                            gap_lt: int = 5) -> Dict[str, float]:
    """Compute the three validation statistics from an ERROR TRACK.

    Source-agnostic: ``per_read_events`` may come from ``inject`` (exact truth,
    no alignment) OR from a BAM adapter (CIGAR/MD error runs grouped by
    ``_group_pos_into_runs``). The unit is the error EVENT throughout.
    """
    counts = [len(ev) for ev in per_read_events]
    lengths = list(read_lengths)
    total_events = sum(counts)
    total_len = sum(lengths)
    mean_rate = (total_events / total_len) if total_len else 0.0

    # 1) length-adjusted per-read dispersion index (exposure-adjusted index of
    #    dispersion): mean over reads of (count - rate*len)^2 / (rate*len).
    #    -> 1.0 for Poisson, > 1 for over-dispersion. Reads with rate*len==0 skip.
    disp_terms = []
    rates = []
    for c, L in zip(counts, lengths):
        if L <= 0:
            continue
        rates.append(c / L)
        expected = mean_rate * L
        if expected > 0:
            disp_terms.append((c - expected) ** 2 / expected)
    dispersion_index = (sum(disp_terms) / len(disp_terms)) if disp_terms else 0.0

    # per-read rate tail (the layer-1 tail statistic)
    rates_sorted = sorted(rates)
    def _q(q: float) -> float:
        if not rates_sorted:
            return 0.0
        idx = min(len(rates_sorted) - 1, int(q * len(rates_sorted)))
        return rates_sorted[idx]
    median_rate = _q(0.5)
    p90_rate = _q(0.9)
    p90_over_median = (p90_rate / median_rate) if median_rate > 0 else 0.0

    # 2) within-read sub-gap_lt EXCESS over the geometric null at the marginal rate
    gaps: List[int] = []
    for ev in per_read_events:
        starts = sorted(e.pos for e in ev)
        for a, b in zip(starts, starts[1:]):
            gaps.append(b - a)
    if gaps:
        obs_frac = sum(1 for g in gaps if g < gap_lt) / len(gaps)
        null_frac = _geometric_p_lt(mean_rate, gap_lt)
        gap_excess = (obs_frac / null_frac) if null_frac > 0 else 0.0
    else:
        obs_frac = null_frac = gap_excess = 0.0

    # 3) fraction of INDEL events >= 2bp
    indel_lens = [e.length for ev in per_read_events for e in ev if e.kind in ("ins", "del")]
    frac_ge2 = (sum(1 for L in indel_lens if L >= 2) / len(indel_lens)) if indel_lens else 0.0

    # length-INVARIANT over-dispersion summary (guard #2): dispersion_index is
    # ~ 1 + rate*L*v with v = Var(per-read multiplier), so v ~ (disp-1)/mean_events
    # is comparable across read sources with different length distributions. (For
    # real data it conflates the multiplier + burst contributions, but is the
    # length-normalized number to compare real vs pbsim+injector on.)
    mean_events = (total_events / len(counts)) if counts else 0.0
    overdisp_v = ((dispersion_index - 1.0) / mean_events) if mean_events > 0 else 0.0
    return {
        "n_reads": float(len(counts)),
        "mean_ref_len": (total_len / len(counts)) if counts else 0.0,
        "marginal_event_rate": mean_rate,
        "mean_events_per_read": mean_events,
        "dispersion_index": dispersion_index,
        "overdisp_v": overdisp_v,
        "p90_over_median": p90_over_median,
        "sub5_gap_excess": gap_excess,
        "obs_gap_frac_lt": obs_frac,
        "null_gap_frac_lt": null_frac,
        "indel_run_frac_ge2": frac_ge2,
        "n_indel_events": float(len(indel_lens)),
    }


# ---------------------------------------------------------------------------
# BAM adapter — the SAME measurement on REAL / pbsim / pbsim+injector reads.
# The event UNIT matches the injector's generative unit: each mismatch = one
# 'sub' event (length 1); each CIGAR I/D op = one 'ins'/'del' event of its
# length. Positions are REFERENCE coordinates; the read-length denominator is
# the reference-aligned span (M/=/X/D), matching the injector's clean-ref length.
# Kept as a pure function (no pysam) so the grouping is M1-unit-testable.
# ---------------------------------------------------------------------------
def events_from_alignment(cigartuples: Sequence[Tuple[int, int]],
                          mismatch_ref_positions: Sequence[int],
                          ref_start: int = 0,
                          exonic_coords: bool = False) -> Tuple[List[ErrorEvent], int]:
    """Build the ERROR TRACK + reference-aligned span from a parsed alignment.

    ``cigartuples`` are pysam-style (op, length) with op codes
    0=M 1=I 2=D 3=N 4=S 5=H 7== 8=X. ``mismatch_ref_positions`` are the 0-based
    reference positions of substitutions (from MD / get_aligned_pairs(with_seq)).
    Returns (events sorted by ref pos, ref_span). Insertions are placed at the
    current reference position; N (intron) and clips are skipped (not errors).

    ``exonic_coords`` (the SPLICED-READ fix): when True, the intron (N) length is
    EXCLUDED from ``ref_span`` and every event position is mapped to EXON-LOCAL
    (intron-collapsed) coordinates. This matters because the rate (events/span)
    and the within-read gap/clustering statistics must be measured along the read
    MOLECULE, not along the genome — otherwise a spliced read's span is inflated
    by its introns (deflating the rate) and an error gap straddling an intron is
    counted as a huge spurious gap (corrupting sub5_gap_excess / dispersion). With
    no N op the result is identical to the default, so callers measuring rate /
    clustering can always pass ``exonic_coords=True`` safely. Default stays False
    for backward compatibility with callers that need absolute genome positions.
    """
    events: List[ErrorEvent] = []
    rpos = ref_start
    ref_span = 0
    introns: List[Tuple[int, int]] = []   # (abs_start, length)
    mm = set(mismatch_ref_positions)
    for op, length in cigartuples:
        if op in (0, 7, 8):          # M / = / X consume both
            ref_span += length
            rpos += length
        elif op == 2:                # D
            events.append(ErrorEvent(pos=rpos, kind="del", length=length, bases=""))
            ref_span += length
            rpos += length
        elif op == 1:                # I — at current ref position
            events.append(ErrorEvent(pos=rpos, kind="ins", length=length, bases=""))
        elif op == 3:                # N intron — skip (not an error), consumes ref
            introns.append((rpos, length))
            ref_span += length
            rpos += length
        # 4=S, 5=H, 6=P consume neither ref nor count
    for p in mm:                     # substitutions as individual len-1 events
        events.append(ErrorEvent(pos=p, kind="sub", length=1, bases=""))
    if exonic_coords and introns:
        total_intron = sum(L for _, L in introns)
        def _to_exon(p: int) -> int:   # subtract intron length lying before p
            return p - sum(L for s, L in introns if s < p)
        events = [ErrorEvent(pos=_to_exon(e.pos), kind=e.kind,
                             length=e.length, bases=e.bases) for e in events]
        ref_span -= total_intron
    events.sort(key=lambda e: e.pos)
    return events, ref_span


def events_from_bam_read(read, exonic_coords: bool = False) -> Tuple[List[ErrorEvent], int]:
    """pysam wrapper: derive the ERROR TRACK from an aligned read.

    Substitutions need the MD tag (``get_aligned_pairs(with_seq=True)`` returns a
    lowercase ref base at mismatches). Requires the BAM to carry MD (minimap2
    ``--MD`` or ``calmd``); raises if absent so a silent zero is impossible.

    ``exonic_coords`` (see ``events_from_alignment``): exclude introns from the
    span and map positions to exon-local coords — REQUIRED for correct rate /
    clustering on spliced (junction-spanning) reads."""
    mm = []
    for qpos, rpos, refbase in read.get_aligned_pairs(with_seq=True):
        if qpos is not None and rpos is not None and refbase is not None and refbase.islower():
            mm.append(rpos)
    return events_from_alignment(read.cigartuples or [], mm,
                                 ref_start=read.reference_start,
                                 exonic_coords=exonic_coords)


def head_tail_autocorr(params: InjectorParams, rng: random.Random,
                       n_reads: int = 3000, read_len: int = 800,
                       frac: float = 0.3) -> float:
    """Synthetic analogue of `error_realism_validate autocorr-bam`: head-vs-tail
    error-density Pearson r across reads (truth track, no alignment). This is the
    GLOBAL-vs-LOCAL discriminator the marginals (overdisp_v, gap5x) cannot give —
    the global multiplier (layer 1) scales both windows together (high r); local
    bursts (layer 2) do not (r~0). Pick PLACEHOLDER params so this matches the
    real-BAM r (~0.3, measured 2026-06-27), which pins the global fraction."""
    w = int(read_len * frac)
    h0, h1 = 0, w
    t0, t1 = read_len - w, read_len
    hd, td = [], []
    for _ in range(n_reads):
        _, ev = _walk_events(read_len, params, rng, clean_seq=None)
        hd.append(sum(1 for e in ev if h0 <= e.pos < h1) / w)
        td.append(sum(1 for e in ev if t0 <= e.pos < t1) / w)
    n = len(hd)
    mh, mt = sum(hd) / n, sum(td) / n
    shh = sum((x - mh) ** 2 for x in hd)
    stt = sum((x - mt) ** 2 for x in td)
    sht = sum((x - mh) * (y - mt) for x, y in zip(hd, td))
    return sht / math.sqrt(shh * stt) if shh > 0 and stt > 0 else 0.0


def simulate_and_measure(params: InjectorParams, rng: random.Random,
                         n_reads: int = 2000, read_len: int = 1000,
                         ) -> Dict[str, float]:
    """Inject into ``n_reads`` random clean reads and measure the structure.

    Uses the EVENTS-ONLY fast path (no dirty-sequence construction): the three
    statistics depend only on event positions/kinds/lengths, not on the actual
    substituted/inserted bases."""
    per_read: List[List[ErrorEvent]] = []
    lengths: List[int] = []
    for _ in range(n_reads):
        _, ev = _walk_events(read_len, params, rng, clean_seq=None)
        per_read.append(ev)
        lengths.append(read_len)
    return measure_error_structure(per_read, lengths)


# ---------------------------------------------------------------------------
# Calibration — coordinate descent against the COMBINED output (handles the
# cross-layer contamination the advisor flagged: we never solve a layer in
# isolation, we measure the joint result each round).
# ---------------------------------------------------------------------------
# Targets RE-MEASURED 2026-06-27 with THIS module's measure_error_structure on the
# real BY4742 DRS BAM (minimap2 read-vs-genome, window [600,1000], 50k reads) — the
# guard-#1 re-derivation. The new code REPRODUCED the lost err_corr.py verdict on
# dispersion (9.64 vs 9.95) and indel-run (0.39 vs 0.39); sub5_gap_excess came out
# 5.28 (vs the prior 2.83 — a definitional difference: this stat is gaps between
# error EVENTS in ref coords, grouped per the profiler). We use OUR numbers since
# real/pbsim/injector are all measured by this same code.
# CALIBRATION METRIC = overdisp_v (LENGTH-INVARIANT). The dispersion index is
# ~ 1 + rate*L*v (LINEAR in read length L), so it is NOT comparable across sources
# with different length distributions; v = Var(per-read multiplier) ~
# (disp-1)/mean_events IS comparable, and is the parameter (gamma_shape ~ 1/v).
# CAVEATS still in force: these are read-vs-reference (UPPER BOUND — RNA-mod +
# minimap2 pile-up inflate sub5_gap_excess especially); the SIRV absolute-truth
# re-fit sets the true magnitude. Error TYPE split below is still PLACEHOLDER.
REAL_TARGETS = {
    "overdisp_v": 0.70,            # length-invariant; gamma_shape ~ 1/v
    "dispersion_index": 9.64,      # at mean_ref_len ~807; length-DEPENDENT, informational
    "sub5_gap_excess": 5.28,       # incl. alignment-artifact inflation (read-vs-ref); SIRV-gated
    "indel_run_frac_ge2": 0.39,
    "marginal_event_rate": 0.0153,
    "head_tail_autocorr": 0.30,    # measured on the real BAM (frac 0.3, window [600,1000]);
                                   # the GLOBAL-fraction constraint PI #2 depends on (overdisp_v
                                   # + gap5x cannot identify it). PI-#2-critical + magnitude-robust.
}

# ===================== HAND-PICKED PLACEHOLDER PARAMS =====================
# Locked, NOT auto-fit (advisor: don't auto-fit a SIRV-gated magnitude). Chosen
# from the M1 probe to match the real-data SPLIT — head_tail_autocorr ~= 0.30 (the
# PI-#2-critical global fraction) — with gap5x ELEVATED and run>=2 ~= 0.39, at the
# real marginal rate. KNOWN MODEL LIMITATION (key input to the SIRV-fit cycle): the
# 2-state burst HMM CANNOT jointly reproduce real's (overdisp_v 0.70, gap5x ~5,
# autocorr 0.30) — longer bursts (needed for gap5x) raise per-read hot-fraction
# variance => raise autocorr, so strong clustering and moderate global correlation
# are COUPLED here, whereas real DECOUPLES them (real over-dispersion is mostly
# LOCAL). Matching real likely needs a 3-state / self-exciting (Hawkes) burst — do
# that in the SIRV cycle, NOT now. This placeholder prioritizes the autocorr split
# (so it does NOT overstate the PI-#2 covariate); achieved overdisp_v ~0.27 < real
# 0.70 is the documented model gap. PLACEHOLDER-PENDING-SIRV.
def placeholder_params(base_rate: float = 0.0153) -> "InjectorParams":
    return InjectorParams(
        base_rate=base_rate,
        read_rate_dist="gamma", gamma_shape=5.0,
        burst_on=True, hot_factor=8.0, p_hot_to_cold=0.2, p_cold_to_hot=0.025,
        indel_run_pmf=_indel_pmf_for_frac_ge2(REAL_TARGETS["indel_run_frac_ge2"]),
    )


def _indel_pmf_for_frac_ge2(target: float) -> Tuple[Tuple[int, float], ...]:
    """A fat-tailed run-length pmf (geometric on {1,2,3,...}) with P(len>=2) =
    target. Geometric P(len=k) = (1-q)^(k-1) q => P(>=2) = 1-q => q = 1-target."""
    q = max(1e-3, min(0.999, 1.0 - target))
    pmf = []
    cum = 0.0
    for k in range(1, 9):
        pk = (1.0 - q) ** (k - 1) * q
        pmf.append((k, pk))
        cum += pk
    # dump the residual tail mass onto the last bin so probs sum to 1
    pmf.append((9, max(0.0, 1.0 - cum)))
    return tuple(pmf)


def calibrate_params(rng: random.Random,
                     targets: Optional[Dict[str, float]] = None,
                     n_reads: int = 800, read_len: int = 600,
                     rounds: int = 3, base_rate: float = 0.0153,
                     verbose: bool = False) -> Tuple[InjectorParams, Dict[str, float]]:
    """Fit PLACEHOLDER realistic params by coordinate descent on the joint output.

    Knobs (one per target, to stay just-identified — the HMM burst GEOMETRY
    ``p_hot_to_cold`` is FIXED on the M1; the full fit of the burst to the
    gap-SIZE distribution is the Sherlock step, advisor note #2):
      * layer 1 -> gamma_shape         (smaller => higher overdisp_v)
      * layer 2 -> p_cold_to_hot       (higher  => higher sub5 gap excess)
      * layer 3 -> indel_run_pmf       (closed-form for P(>=2))
    The layer-1 target is the LENGTH-INVARIANT overdisp_v (NOT dispersion_index,
    which is ~linear in read length); calibrate at the source's real base_rate so
    mean_events matches. Returns (params, measured-structure). PLACEHOLDER.

    SCAFFOLD STATUS: this auto-fit is NOT used to lock the M1 placeholder (we
    hand-pick `placeholder_params` instead — advisor: don't auto-fit a SIRV-gated,
    alignment-contaminated magnitude). It is kept as the scaffold for the eventual
    SIRV absolute-truth fit, where the gap5x target is clean. NOTE the 2-state HMM
    cannot jointly hit (overdisp_v, gap5x, head_tail_autocorr) — see
    `placeholder_params`; the SIRV fit will likely need a 3-state/self-exciting burst
    and should fit the gap-SIZE distribution + the autocorr, not these 2 marginals.
    """
    t = dict(REAL_TARGETS)
    if targets:
        t.update(targets)

    p = InjectorParams(
        base_rate=base_rate,
        read_rate_dist="gamma", gamma_shape=0.3,
        burst_on=True, p_cold_to_hot=0.05, p_hot_to_cold=0.4, hot_factor=8.0,
        indel_run_pmf=_indel_pmf_for_frac_ge2(t["indel_run_frac_ge2"]),
    )

    def measure(pp: InjectorParams) -> Dict[str, float]:
        return simulate_and_measure(pp, random.Random(rng.random()),
                                    n_reads=n_reads, read_len=read_len)

    # layer 3 is closed-form; lock it once.
    p = replace(p, indel_run_pmf=_indel_pmf_for_frac_ge2(t["indel_run_frac_ge2"]))

    lo_shape, hi_shape = 0.02, 5.0
    lo_pch, hi_pch = 0.005, 0.5
    m = measure(p)
    for _ in range(rounds):
        # --- layer 1: bisect gamma_shape on overdisp_v (monotone decreasing) ---
        a, b = lo_shape, hi_shape
        for _ in range(9):
            mid = math.sqrt(a * b)
            mm = measure(replace(p, gamma_shape=mid))
            if mm["overdisp_v"] > t["overdisp_v"]:
                a = mid          # too dispersed => raise shape
            else:
                b = mid
        p = replace(p, gamma_shape=math.sqrt(a * b))

        # --- layer 2: bisect p_cold_to_hot on sub5 gap excess ---
        # MONOTONICITY (empirically, M1 probe 2026-06-27): LOWER p_cold_to_hot =>
        # MORE clustering (rare high-contrast bursts), the OPPOSITE of the naive
        # reading. So if measured < target (too little clustering), LOWER it.
        a, b = lo_pch, hi_pch
        for _ in range(9):
            mid = math.sqrt(a * b)
            mm = measure(replace(p, p_cold_to_hot=mid))
            if mm["sub5_gap_excess"] < t["sub5_gap_excess"]:
                b = mid          # too little clustering => LOWER p_cold_to_hot
            else:
                a = mid
        p = replace(p, p_cold_to_hot=math.sqrt(a * b))

        m = measure(p)
        if verbose:
            print(f"[calib] shape={p.gamma_shape:.3f} pch={p.p_cold_to_hot:.4f} "
                  f"disp={m['dispersion_index']:.2f} gap5x={m['sub5_gap_excess']:.2f} "
                  f"run>=2={m['indel_run_frac_ge2']:.2f}")
    return p, m


# ---------------------------------------------------------------------------
# Composition / interaction self-check (M1-local). NOT a realism validation.
# ---------------------------------------------------------------------------
def self_check(seed: int = 7, verbose: bool = True) -> bool:
    rng = random.Random(seed)
    ok = True

    def emit(msg: str):
        if verbose:
            print(msg)

    # (1) NULL params reproduce the marginal / Poisson baseline.
    null = InjectorParams.null()
    null = replace(null, base_rate=0.02)
    mn = simulate_and_measure(null, random.Random(11), n_reads=1500, read_len=600)
    null_ok = (0.8 <= mn["dispersion_index"] <= 1.25
               and 0.85 <= mn["sub5_gap_excess"] <= 1.15
               and mn["indel_run_frac_ge2"] <= 0.02)
    emit(f"[self-check] (1) NULL baseline: disp={mn['dispersion_index']:.2f} "
         f"gap5x={mn['sub5_gap_excess']:.2f} run>=2={mn['indel_run_frac_ge2']:.3f} "
         f"-> {'PASS' if null_ok else 'FAIL'} (expect ~1 / ~1 / ~0, Poisson single-base)")
    ok &= null_ok

    # (2) Each layer moves its OWN statistic in the right direction, in isolation.
    base = replace(InjectorParams.null(), base_rate=0.02)
    l1 = replace(base, read_rate_dist="gamma", gamma_shape=0.25)
    m1 = simulate_and_measure(l1, random.Random(12), n_reads=1500, read_len=600)
    l1_ok = m1["dispersion_index"] > 2.0 * mn["dispersion_index"]
    emit(f"[self-check] (2a) layer1 over-dispersion: disp {mn['dispersion_index']:.2f} "
         f"-> {m1['dispersion_index']:.2f} -> {'PASS' if l1_ok else 'FAIL'}")
    ok &= l1_ok

    l2 = replace(base, burst_on=True, p_cold_to_hot=0.05, p_hot_to_cold=0.4, hot_factor=8.0)
    m2 = simulate_and_measure(l2, random.Random(13), n_reads=1500, read_len=600)
    l2_ok = m2["sub5_gap_excess"] > 1.4 * mn["sub5_gap_excess"]
    emit(f"[self-check] (2b) layer2 burst clustering: gap5x {mn['sub5_gap_excess']:.2f} "
         f"-> {m2['sub5_gap_excess']:.2f} -> {'PASS' if l2_ok else 'FAIL'}")
    ok &= l2_ok

    l3 = replace(base, indel_run_pmf=_indel_pmf_for_frac_ge2(0.4))
    m3 = simulate_and_measure(l3, random.Random(14), n_reads=1500, read_len=600)
    l3_ok = m3["indel_run_frac_ge2"] > 0.25
    emit(f"[self-check] (2c) layer3 indel runs: run>=2 {mn['indel_run_frac_ge2']:.3f} "
         f"-> {m3['indel_run_frac_ge2']:.3f} -> {'PASS' if l3_ok else 'FAIL'}")
    ok &= l3_ok

    # (3) The hand-picked PLACEHOLDER composes — all layers ON produce REACHABLE
    #     composition targets AND no layer zeroes another. NOT a realism check and
    #     NOT pinned to the contaminated/SIRV-gated magnitude (advisor): we assert
    #     elevated/reachable bands + the autocorr SPLIT matches real (~0.30).
    pp = placeholder_params()
    m = simulate_and_measure(pp, random.Random(15), n_reads=2500, read_len=800)
    ac = head_tail_autocorr(pp, random.Random(16), n_reads=2500, read_len=800)
    v_ok = m["overdisp_v"] > 0.1                       # over-dispersion present
    gap_ok = m["sub5_gap_excess"] > 2.5                # clustering elevated
    run_ok = abs(m["indel_run_frac_ge2"] - REAL_TARGETS["indel_run_frac_ge2"]) <= 0.06
    ac_ok = abs(ac - REAL_TARGETS["head_tail_autocorr"]) <= 0.15   # global SPLIT ~ real
    nonzero = m["overdisp_v"] > 0.0 and m["n_indel_events"] > 0
    comp_ok = v_ok and gap_ok and run_ok and ac_ok and nonzero
    emit(f"[self-check] (3) hand-picked PLACEHOLDER composes: "
         f"overdisp_v={m['overdisp_v']:.3f}(>0.1) gap5x={m['sub5_gap_excess']:.2f}(>2.5) "
         f"run>=2={m['indel_run_frac_ge2']:.3f}(~0.39) autocorr={ac:.3f}(~0.30 real) "
         f"-> {'PASS' if comp_ok else 'FAIL'}")
    emit("[self-check] NOTE: (3) is COMPOSITION + the global SPLIT matching real — NOT a "
         "magnitude/realism check (gap5x magnitude is SIRV-gated; the 2-state HMM can't "
         "jointly match real overdisp_v+gap5x+autocorr, see placeholder_params).")
    ok &= comp_ok

    emit("\n" + ("SELF-CHECK PASSED — composition/interaction green (NOT realism)"
                 if ok else "SELF-CHECK FAILED"))
    return ok


if __name__ == "__main__":
    import sys
    sys.exit(0 if self_check() else 1)
