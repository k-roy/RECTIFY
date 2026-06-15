#!/usr/bin/env python3
"""Tests for the terminal (AAG/GAA)n triplet-repeat basecaller-artifact strip in
the DRS pre-trim (``rectify.core.commands.drs_trim_command``).

Background: Dorado v5.2.0 / RNA004 mis-basecalls the poly-A homopolymer tail as a
low-period multi-base repeat (period-3 (AAG/GAA)n). The poly-A scan cannot see it
(no terminal A-run), so it is carried into alignment and force-aligned. The strip
peels the artifact block before the poly-A scan, recovering any genuine poly-A
behind it (the "blocked poly-A" case) and removing the junk otherwise.

The strip MUST be a strict no-op on clean poly-A tails (period-1 homopolymer).

NOTE: this is the correct file for ``find_polya_and_adapter`` / ``drs_trim_command``
coverage. ``tests/test_polya_trimming.py`` tests a different module
(``rectify.core.polya.polya_trimmer``).
"""

import pytest

from rectify.core.commands.drs_trim_command import (
    find_polya_and_adapter,
    _find_terminal_repeat_block,
)

# A non-poly-A, non-repeating "gene body". Ends in two non-A bases ("CG") so the
# phase walk breaks cleanly at the body/artifact junction (two consecutive
# off-phase non-A bases -> stop), making repeat_len deterministic in tests.
BODY = "GCTAGCTAGCTAGCTACGATCGATCGATCG"
assert BODY[-1] != "A" and BODY[-2] != "A"


# ---------------------------------------------------------------------------
# _find_terminal_repeat_block — direct unit tests
# ---------------------------------------------------------------------------

class TestFindTerminalRepeatBlock:

    def test_pure_gaa_block(self):
        block_len, motif = _find_terminal_repeat_block(BODY + "GAA" * 10)
        assert block_len == 30
        assert motif == "AAG"

    def test_preserves_polya_run_before_block(self):
        # [body][polyA20][(GAA)10]: block boundary anchored at first G, so the
        # 20-A run is NOT absorbed into the peeled block.
        block_len, motif = _find_terminal_repeat_block(BODY + "A" * 20 + "GAA" * 10)
        assert block_len == 30
        assert motif == "AAG"

    def test_clean_polya_is_noop(self):
        # Pure poly-A tail -> block always includes the trailing A-run -> dominant
        # motif is single-base -> rejected (strict no-op guarantee).
        block_len, motif = _find_terminal_repeat_block(BODY + "A" * 30)
        assert block_len == 0
        assert motif == ""

    def test_pure_homopolymer_is_noop(self):
        assert _find_terminal_repeat_block("A" * 50) == (0, "")
        assert _find_terminal_repeat_block("T" * 50) == (0, "")

    def test_incidental_short_codon_not_stripped(self):
        # Two incidental GAA codons (6 bp) at a genuine 3' end: below min_len.
        block_len, motif = _find_terminal_repeat_block(BODY + "GAAGAA")
        assert block_len == 0

    def test_min_len_gate(self):
        # Default min_len=15 (>= 5 tandem triplet copies), chosen for the SMA
        # repeat-disease context (see _REPEAT_MIN_LEN rationale).
        assert _find_terminal_repeat_block(BODY + "GAA" * 5)[0] == 15  # 15 bp: kept
        assert _find_terminal_repeat_block(BODY + "GAA" * 4)[0] == 0   # 12 bp: dropped
        # An explicit looser floor still works.
        assert _find_terminal_repeat_block(BODY + "GAA" * 3, min_len=9)[0] == 9

    def test_off_phase_basecall_error_tolerated(self):
        # One substitution inside the block (off-phase C where an A was expected)
        # must not truncate the strip — purity stays above min_frac.
        block = list("GAA" * 10)
        block[13] = "C"  # an 'A' position -> off-phase non-A
        seq = BODY + "".join(block)
        block_len, motif = _find_terminal_repeat_block(seq)
        assert block_len == 30
        assert motif == "AAG"

    def test_disabled_via_empty_returns_when_too_short(self):
        assert _find_terminal_repeat_block("GAA") == (0, "")


# ---------------------------------------------------------------------------
# find_polya_and_adapter — integration of the strip with the poly-A scan
# ---------------------------------------------------------------------------

class TestFindPolyaAndAdapterRepeatStrip:

    def test_recovers_polya_behind_artifact(self):
        # [body][polyA20][(GAA)10] -> strip the GAA, recover the 20-A poly-A.
        seq = BODY + "A" * 20 + "GAA" * 10
        polya_len, adapter_seq, last_base, adapter_pass, repeat_len, repeat_motif = \
            find_polya_and_adapter(seq, strip_repeat=True)
        assert repeat_len == 30
        assert repeat_motif == "AAG"
        assert polya_len == 20

    def test_strips_artifact_with_no_recoverable_polya(self):
        # [body][(GAA)10] -> strip the GAA; no poly-A behind it (polya_len 0).
        seq = BODY + "GAA" * 10
        polya_len, _, _, _, repeat_len, repeat_motif = find_polya_and_adapter(
            seq, strip_repeat=True
        )
        assert repeat_len == 30
        assert repeat_motif == "AAG"
        assert polya_len == 0

    def test_clean_polya_unchanged(self):
        # [body][polyA30] -> no strip, full poly-A detected (strict no-op).
        seq = BODY + "A" * 30
        polya_len, _, _, _, repeat_len, repeat_motif = find_polya_and_adapter(
            seq, strip_repeat=True
        )
        assert repeat_len == 0
        assert repeat_motif == ""
        assert polya_len == 30

    def test_incidental_codon_not_stripped(self):
        seq = BODY + "GAAGAA"
        _, _, _, _, repeat_len, _ = find_polya_and_adapter(seq, strip_repeat=True)
        assert repeat_len == 0

    def test_strip_disabled_is_passthrough(self):
        # With the strip off, the artifact is invisible to the poly-A scan exactly
        # as before this feature: repeat_len 0 and the GAA is not peeled.
        seq = BODY + "GAA" * 10
        polya_len, _, _, _, repeat_len, repeat_motif = find_polya_and_adapter(
            seq, strip_repeat=False
        )
        assert repeat_len == 0
        assert repeat_motif == ""
        # poly-A scan sees the trailing "AA" of the last GAA at most (pre-existing).
        assert polya_len <= 2

    def test_strip_off_matches_legacy_on_clean_polya(self):
        seq = BODY + "A" * 30
        on = find_polya_and_adapter(seq, strip_repeat=True)
        off = find_polya_and_adapter(seq, strip_repeat=False)
        assert on[0] == off[0] == 30
        assert on[4] == off[4] == 0

    def test_return_arity(self):
        # All callers unpack a 6-tuple; guard against accidental arity drift.
        result = find_polya_and_adapter("A" * 30)
        assert len(result) == 6
