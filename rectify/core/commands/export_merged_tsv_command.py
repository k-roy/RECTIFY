"""rectify export-merged-tsv — concatenate a manifest into a single merged TSV.

For users with downstream scripts that expect a single concatenated
corrected_reads.tsv. Reads a manifest, validates its sha256s, and emits
the merged form. Idempotent: re-running on the same manifest produces
the same bytes.

This is the back-compat shim introduced in Commit B (2026-05-20) when the
default output of ``rectify correct`` switched from a monolithic
``corrected_reads.tsv`` to a ``corrected_reads.manifest.tsv`` pointing to
per-region TSVs.
"""

import argparse
import hashlib
import sys
from pathlib import Path
from typing import Optional


def _sha256_of_file(path: Path) -> str:
    """Return hex SHA-256 digest of file contents."""
    h = hashlib.sha256()
    with path.open('rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def run(args) -> int:
    """Concatenate a corrected_reads.manifest.tsv into a single merged TSV."""
    manifest_path = Path(args.manifest)
    if not manifest_path.exists():
        print(f"error: manifest not found: {manifest_path}", file=sys.stderr)
        return 1

    from rectify.core.bam.tsv_partition import load_manifest
    try:
        entries = load_manifest(manifest_path)
    except (ValueError, FileNotFoundError) as exc:
        print(f"error: cannot read manifest: {exc}", file=sys.stderr)
        return 1

    if not entries:
        print("error: manifest is empty — nothing to merge.", file=sys.stderr)
        return 2

    out_path = Path(args.output) if args.output else (
        manifest_path.parent / "corrected_reads.tsv"
    )

    # Validate sha256s upfront; fail fast on any mismatch.
    for entry in entries:
        region_tsv = entry['tsv_path']  # absolute Path from load_manifest
        if not region_tsv.exists():
            print(f"error: region TSV missing: {region_tsv}", file=sys.stderr)
            return 3
        observed = _sha256_of_file(region_tsv)
        expected = entry['sha256']
        if observed != expected:
            print(
                f"error: sha256 mismatch on {region_tsv.name}: "
                f"expected {expected[:16]}.., got {observed[:16]}..",
                file=sys.stderr,
            )
            return 4

    # Concatenate. Header from the first region TSV; skip headers on the rest.
    total_rows = 0
    with out_path.open('w') as out_fh:
        for i, entry in enumerate(entries):
            region_tsv = entry['tsv_path']
            with region_tsv.open() as in_fh:
                first_line = next(in_fh)
                if i == 0:
                    out_fh.write(first_line)  # keep header on first file only
                for line in in_fh:
                    out_fh.write(line)
                    total_rows += 1

    print(f"wrote merged TSV: {out_path} ({total_rows:,} rows)")
    return 0


def create_export_merged_tsv_parser(subparsers) -> argparse.ArgumentParser:
    """Wire the ``export-merged-tsv`` subcommand into the given subparsers."""
    parser = subparsers.add_parser(
        'export-merged-tsv',
        help='Concatenate a corrected_reads.manifest.tsv into a single merged TSV (back-compat shim)',
        description=(
            'Reads a corrected_reads.manifest.tsv (produced by rectify correct in default '
            'manifest-only mode) and emits the legacy concatenated corrected_reads.tsv form '
            'expected by older downstream scripts. Idempotent. Equivalent to `cat` over the '
            'per-region TSVs but with sha256 validation against the manifest.\n\n'
            'Example:\n'
            '    rectify export-merged-tsv results/sample/corrected_reads.manifest.tsv\n'
            '    rectify export-merged-tsv results/sample/corrected_reads.manifest.tsv '
            '-o /tmp/merged.tsv'
        ),
    )
    parser.add_argument(
        'manifest',
        type=str,
        help='Path to corrected_reads.manifest.tsv',
    )
    parser.add_argument(
        '-o', '--output',
        type=str,
        default=None,
        help='Output TSV path (default: <manifest_dir>/corrected_reads.tsv)',
    )
    parser.set_defaults(func=run)
    return parser
