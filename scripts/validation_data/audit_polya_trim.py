#!/usr/bin/env python3
"""Audit `rectify trim-polya` on every validation read.

For each record in `validation_reads_dorado_source.bam` we reconstruct the
RNA-oriented 3' end of the read (last 150 bp in 5'->3' order) and run the
canonical detector `find_polya_and_adapter` from `drs_trim_command`. The
result is a TSV with one row per read plus a summary of invariant flags.

Flags:
  A   adapter detected (pass > 0) but trim skipped because polya_len < 1
  B   trim applied AND boundary_base == 'A' -- would violate "trim to
      first non-A" (impossible by construction of `_scan_polya`; flagged
      as a sanity check on the audit itself, not the trimmer)
  C   Pass 2 (iterative peel for A-miscall stubs) triggered
  D   No stub (pass == 0) and the read does not end in A -- usually means
      the read is dorado-pre-trimmed or terminates inside the body

`source` column is `dorado` for reads coming from the raw Dorado BAM and
`minimap2_fallback` for the two reads (`cat1_plus_1`, `cat1_minus_1`)
that `build_dorado_source.py` pulled from the post-trim minimap2 merged
BAM. For those the audit reflects a sequence the trimmer already ran on
once and is informational only.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional

import pysam

# rectify package import (script lives outside the package)
_REPO = Path(__file__).resolve().parents[2]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from rectify.core.commands.drs_trim_command import find_polya_and_adapter  # noqa: E402
from rectify.utils.genome import reverse_complement  # noqa: E402

# Reads in the dorado-source BAM that originated from the post-trim minimap2
# merged BAM rather than the raw Dorado output (their sequences are already
# trimmed). See HANDOFF.md §1a / build_dorado_source.py Phase 2.
_MINIMAP2_FALLBACK_XV = {"cat1_plus_1", "cat1_minus_1"}


def _rna_three_prime_window(seq: str, is_reverse: bool, window: int = 150) -> str:
    """Return the last ``window`` bp of the read in RNA 5'->3' orientation.

    For a forward-strand BAM record the BAM seq is already RNA-oriented;
    the last ``window`` bp is `seq[-window:]`. For a reverse-strand record
    the BAM seq is the reverse-complement of the RNA, so the last ``window``
    bp of RNA = reverse_complement of the FIRST ``window`` bp of BAM seq.
    """
    if not seq:
        return ""
    if is_reverse:
        head = seq[:window]
        return reverse_complement(head)
    return seq[-window:]


def audit(bam_path: Path, output_tsv: Path, window: int = 150) -> int:
    rows: list[dict] = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam:
            if read.is_secondary or read.is_supplementary:
                continue
            seq = read.query_sequence or ""
            if not seq:
                continue
            xv = read.get_tag("XV") if read.has_tag("XV") else "?"
            strand = "-" if read.is_reverse else "+"
            rna_window = _rna_three_prime_window(seq, read.is_reverse, window)
            polya_len, adapter_seq, last_base, adapter_pass = find_polya_and_adapter(
                rna_window
            )
            trim_applied = polya_len >= 1
            # `last_base` returned by the function is window[-1] (= the
            # very last base of the RNA = the outermost A of the pA tail
            # when one exists). For the "did the trim walk all A's"
            # invariant we want the FIRST base UPSTREAM of the
            # trim cut — i.e. ``rna_window[: trim_cut][-1]`` where
            # ``trim_cut = polya_len + len(adapter_seq)`` counts from end.
            trim_cut = polya_len + len(adapter_seq)
            kept = rna_window[: len(rna_window) - trim_cut] if trim_cut <= len(rna_window) else ""
            post_trim_boundary = kept[-1] if kept else "-"
            source = ("minimap2_fallback"
                      if xv in _MINIMAP2_FALLBACK_XV else "dorado")
            rows.append({
                "XV": xv,
                "qname8": read.query_name[:8],
                "strand": strand,
                "last_50bp_RNA": rna_window[-50:],
                "polya_len": polya_len,
                "adapter_seq": adapter_seq or "-",
                "pass": adapter_pass,
                "trim_applied": "Y" if trim_applied else "N",
                "raw_last_base": last_base or "-",
                "post_trim_boundary": post_trim_boundary,
                "source": source,
            })

    rows.sort(key=lambda r: r["XV"])

    # Write TSV
    cols = ["XV", "qname8", "strand", "last_50bp_RNA", "polya_len",
            "adapter_seq", "pass", "trim_applied",
            "raw_last_base", "post_trim_boundary", "source"]
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(output_tsv, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    # ----- Summary -----
    n = len(rows)
    by_pass = {0: 0, 1: 0, 2: 0}
    for r in rows:
        by_pass[r["pass"]] += 1
    trim_yes = sum(1 for r in rows if r["trim_applied"] == "Y")
    flag_A = [r for r in rows if r["pass"] > 0 and r["polya_len"] == 0]
    flag_B = [r for r in rows
              if r["trim_applied"] == "Y" and r["post_trim_boundary"] == "A"]
    flag_C = [r for r in rows if r["pass"] == 2]
    flag_D = [r for r in rows
              if r["pass"] == 0 and r["polya_len"] == 0
              and r["raw_last_base"] not in ("A", "-")]
    fallback = [r for r in rows if r["source"] == "minimap2_fallback"]

    print(f"\nAudit: {n} reads from {bam_path.name}")
    print(f"  source = dorado:             {n - len(fallback)}")
    print(f"  source = minimap2_fallback:  {len(fallback)}  "
          f"({', '.join(r['XV'] for r in fallback)})")
    print(f"  trim_applied (polya_len>=1): {trim_yes}/{n}")
    print(f"  pass 0 / 1 / 2:              "
          f"{by_pass[0]} / {by_pass[1]} / {by_pass[2]}")
    print()
    print(f"  flag A (adapter detected, trim skipped): {len(flag_A)}")
    for r in flag_A:
        print(f"    {r['XV']:<14} {r['qname8']} pass={r['pass']} "
              f"adapter={r['adapter_seq']!r} last_50={r['last_50bp_RNA']!r}")
    print(f"  flag B (trim_applied AND post_trim_boundary==A): {len(flag_B)} "
          f"(expected 0)")
    for r in flag_B:
        print(f"    {r['XV']:<14} {r['qname8']} polya_len={r['polya_len']} "
              f"post_trim_boundary={r['post_trim_boundary']} "
              f"last_50={r['last_50bp_RNA']!r}")
    print(f"  flag C (Pass 2 triggered):                  {len(flag_C)}")
    for r in flag_C:
        print(f"    {r['XV']:<14} {r['qname8']} polya_len={r['polya_len']} "
              f"adapter={r['adapter_seq']!r}")
    print(f"  flag D (no-stub but raw_last_base non-A):   {len(flag_D)}")
    for r in flag_D:
        print(f"    {r['XV']:<14} {r['qname8']} raw_last={r['raw_last_base']} "
              f"last_50={r['last_50bp_RNA']!r}")
    print(f"\nWrote {output_tsv}")
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("bam", type=Path,
                   help="dorado-source BAM (validation_reads_dorado_source.bam)")
    p.add_argument("--output", "-o", type=Path,
                   default=Path("validation_read_review/polya_trim_audit.tsv"),
                   help="output TSV path (default: validation_read_review/polya_trim_audit.tsv)")
    p.add_argument("--window", type=int, default=150,
                   help="3' window length passed to find_polya_and_adapter (default 150)")
    args = p.parse_args(argv)
    if not args.bam.exists():
        print(f"ERROR: BAM not found: {args.bam}", file=sys.stderr)
        return 1
    return audit(args.bam, args.output, args.window)


if __name__ == "__main__":
    sys.exit(main())
