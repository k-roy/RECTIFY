#!/usr/bin/env python3
"""Per-region TSV emission + manifest.

Each region worker writes its own corrected_reads.region_NNN.tsv.
The manifest is the canonical artifact; a single concatenated TSV is
optionally emitted for back-compat via --emit-merged-tsv (which lives in
Commit B's correct_command wiring, NOT here).

Author: Kevin R. Roy
"""

import hashlib
import os
from pathlib import Path
from typing import Dict, List, Optional

from .regions import RegionPlan


class RegionTsvWriter:
    """Context manager for writing per-region corrected TSV files.

    Writes rows to a per-region TSV file while incrementally computing
    a SHA-256 hash and counting rows. The hash and row count are available
    after __exit__ for inclusion in the manifest.

    Usage::

        with RegionTsvWriter(plan, header) as w:
            w.write_row(row)
        manifest_entry = (plan.region_id, w.path, w.n_rows, w.sha256_hex)
    """

    def __init__(self, region_plan: RegionPlan, header: List[str]):
        self.plan = region_plan
        self.header = header
        self.path = region_plan.tmp_dir / f"{region_plan.region_id}.tsv"
        self._fh = None
        self._n_rows = 0
        self._sha256 = hashlib.sha256()  # hash as we write, no separate pass

    def __enter__(self) -> "RegionTsvWriter":
        self._fh = open(self.path, 'w')
        line = '\t'.join(self.header) + '\n'
        self._fh.write(line)
        self._sha256.update(line.encode())
        return self

    def write_row(self, row: List[str]) -> None:
        """Write a single data row (list of string values)."""
        line = '\t'.join(row) + '\n'
        self._fh.write(line)
        self._sha256.update(line.encode())
        self._n_rows += 1

    def __exit__(self, *args) -> None:
        if self._fh is not None:
            self._fh.close()
            self._fh = None
        # Caller reads .n_rows and .sha256_hex for manifest entry

    @property
    def n_rows(self) -> int:
        return self._n_rows

    @property
    def sha256_hex(self) -> str:
        return self._sha256.hexdigest()


# Manifest header columns
_MANIFEST_HEADER = ['region_id', 'chrom', 'start', 'end', 'tsv_path', 'n_rows', 'sha256']


def write_manifest(
    manifest_path: Path,
    region_plans: List[RegionPlan],
    region_writers: List[RegionTsvWriter],
) -> None:
    """Write a TSV manifest listing all per-region TSV files.

    After all region workers complete, the parent calls this function with
    each worker's RegionTsvWriter (for sha256 + n_rows) and RegionPlan (for
    coordinates). Rows are sorted by region_id (which is already coord-ordered).

    Paths in the manifest are stored RELATIVE to the manifest file's directory
    so that moving the directory doesn't break the manifest.

    Args:
        manifest_path:  Destination manifest file. Parent directory must exist.
        region_plans:   RegionPlan objects (one per region, coord-sorted order).
        region_writers: Corresponding RegionTsvWriter objects (post-__exit__).
    """
    assert len(region_plans) == len(region_writers), (
        f"plan/writer count mismatch: {len(region_plans)} != {len(region_writers)}"
    )
    manifest_dir = manifest_path.parent
    rows = []
    for plan, writer in zip(region_plans, region_writers):
        rel_path = os.path.relpath(str(writer.path), str(manifest_dir))
        rows.append([
            plan.region_id,
            plan.chrom,
            str(plan.start),
            str(plan.end),
            rel_path,
            str(writer.n_rows),
            writer.sha256_hex,
        ])
    # Sort by region_id (already coord-sorted by construction, but be defensive)
    rows.sort(key=lambda r: r[0])
    with open(manifest_path, 'w') as fh:
        fh.write('\t'.join(_MANIFEST_HEADER) + '\n')
        for row in rows:
            fh.write('\t'.join(row) + '\n')


def load_manifest(manifest_path: Path) -> List[Dict]:
    """Read a manifest file and return a list of dicts (one per region).

    Verifies the header schema. Paths in the returned dicts are resolved
    to absolute paths relative to the manifest file's directory.

    Args:
        manifest_path:  Path to the manifest TSV.

    Returns:
        List of dicts with keys: region_id, chrom, start, end, tsv_path,
        n_rows, sha256. tsv_path is an absolute Path object.

    Raises:
        ValueError: If the header doesn't match the expected schema.
        FileNotFoundError: If manifest_path doesn't exist.
    """
    manifest_dir = manifest_path.parent
    records = []
    with open(manifest_path, 'r') as fh:
        header_line = fh.readline().rstrip('\n')
        if not header_line:
            raise ValueError(f"Empty manifest file: {manifest_path}")
        actual_cols = header_line.split('\t')
        for required in _MANIFEST_HEADER:
            if required not in actual_cols:
                raise ValueError(
                    f"Manifest {manifest_path} missing required column '{required}'. "
                    f"Found: {actual_cols}"
                )
        col_idx = {col: i for i, col in enumerate(actual_cols)}
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if not parts or all(p == '' for p in parts):
                continue
            rel_path = parts[col_idx['tsv_path']]
            abs_path = (manifest_dir / rel_path).resolve()
            records.append({
                'region_id': parts[col_idx['region_id']],
                'chrom': parts[col_idx['chrom']],
                'start': int(parts[col_idx['start']]),
                'end': int(parts[col_idx['end']]),
                'tsv_path': abs_path,
                'n_rows': int(parts[col_idx['n_rows']]),
                'sha256': parts[col_idx['sha256']],
            })
    return records
