"""PCB114.24 chemistry + adapter-detection constants shared across cdna modules."""
from __future__ import annotations

import re as _re

# ---- PCB114.24 chemistry ----------------------------------------------------
SSP_FWD = "TTTCTGTTGGTGCTGATATTGCT"
SSP_RC  = "AGCAATATCAGCACCAACAGAAA"
UMI_LEN = 27
BRIDGE_LEN = 3   # the rGrGrG of the TSO that becomes the GGG bridge in BAM SEQ.

# ---- Adapter / polyA detection ---------------------------------------------
# End-concordance check (XF tier): the polyA/polyT homopolymer at the OTHER end
# of the read (opposite the SSP/UMI side). For PCB114 chemistry, the basecalled
# 3' end goes through the pore first and is reliable.
END_WINDOW_BP = 200            # unanchored window (last/first N bp of BAM SEQ)
MIN_HOMOPOLYMER_UNANCHORED = 10  # unanchored A/T-run threshold (XF=1 tier)
MIN_HOMOPOLYMER_ANCHORED = 6     # anchored threshold (XF=2 tier)
ANCHOR_UPSTREAM_WIN = 30         # bp window upstream/downstream of anchor to test for polyA/polyT

# Adapter-start anchor — found empirically at ~72 bp from read 3' end in fwd-orient
# reads, 0 occurrences in 230 kb of chrI yeast genome (essentially zero FP combined
# with polyA).
ANCHOR_FWD = "TTGCGGGCGGCGG"   # fwd-orient reads: in last ~300 bp of BAM SEQ, polyA in 30 bp UPSTREAM
ANCHOR_RC  = "CCGCCGCCCGCAA"   # rev-orient reads: in first ~300 bp of BAM SEQ, polyT in 30 bp DOWNSTREAM
ANCHOR_LEN = len(ANCHOR_FWD)
ANCHOR_SEARCH_WIN = 300        # window within which to search for the anchor
# v1.11: fuzzy anchor via edlib HW mode (handles ONT R10.4.1 sub/ins/del). Lev≤2
# on 13-bp anchor → ~95–97% sensitivity (vs ~73% exact) with effectively zero FP
# when combined with the polyA-window constraint. Lev=3 starts to admit incidental
# yeast-genome hits.
ANCHOR_MAX_EDIT = 2

# v0.9 (planning/681): fuzzy SSP search for the consensus pretrim. The exact 23-mer
# `find()` misses the ~29% of consensus molecules whose SSP carries an ONT error, but a
# fuzzy 23-mer search must be WINDOW-GATED to the end where the SSP belongs: an ungated
# edlib HW hit over a ~2 kb adapter-free consensus lands in genomic sequence and trims real
# mRNA. The consensus SSP sits within the first ~100 nt (fwd) / last ~60 nt (rev) — 300 is
# generous.
#
# k=3 is justified by MEASUREMENT, not by arithmetic against ANCHOR_MAX_EDIT. (For scale it
# is marginally the more conservative of the two — 3/23 = 13% of the pattern vs 2/13 = 15% —
# but the note above already warns that Lev=3 admits incidental yeast-genome hits on the
# 13-mer, so a length ratio is not the argument.) The argument is the A/B realignment in
# planning/681 CP5: with this setting, mean aligned reference span moved 1506 -> 1505 nt
# (-1 nt) while ~58 nt of soft clip per molecule disappeared and mapped fraction held at
# 99.9%. A fuzzy search eating real mRNA would have shown up there. Re-run that A/B before
# raising k.
SSP_SEARCH_WIN = 300
SSP_MAX_EDIT = 3

POLY_A_UNANCH_RE = _re.compile(r"A{%d,}" % MIN_HOMOPOLYMER_UNANCHORED)
POLY_T_UNANCH_RE = _re.compile(r"T{%d,}" % MIN_HOMOPOLYMER_UNANCHORED)
POLY_A_ANCH_RE   = _re.compile(r"A{%d,}" % MIN_HOMOPOLYMER_ANCHORED)
POLY_T_ANCH_RE   = _re.compile(r"T{%d,}" % MIN_HOMOPOLYMER_ANCHORED)

# Reverse-complement translation table — used by revcomp() and any module that
# needs to flip a BAM-orient sequence back to basecalled orient.
COMPLEMENT_TABLE = str.maketrans("ACGTacgtN", "TGCAtgcaN")
