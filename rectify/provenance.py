#!/usr/bin/env python3
"""
RECTIFY Provenance Tracking System

Provides lightweight, Snakemake-style provenance for the RECTIFY pipeline:
- SHA-256 hashing of inputs and outputs for data lineage
- Sidecar .provenance.json files alongside every output
- Step-skip logic: if inputs haven't changed, skip re-computation
- Human-readable headers in TSV/bedGraph output files

Usage:
    from rectify.provenance import ProvenanceTracker

    tracker = ProvenanceTracker(output_dir, args)
    if tracker.step_needed('correct', input_files=[bam_path], output_files=[corrected_tsv]):
        # ... run correction ...
        tracker.record_step('correct', input_files=[bam_path], output_files=[corrected_tsv])

Author: Kevin R. Roy
Date: 2026-03-28
"""

import hashlib
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Union

from . import __version__

# ---------------------------------------------------------------------------
# File hashing
# ---------------------------------------------------------------------------

def sha256_file(filepath: Union[str, Path], chunk_size: int = 1 << 20) -> str:
    """
    Compute SHA-256 hash of a file.

    Args:
        filepath: Path to file
        chunk_size: Read buffer size (default 1 MB)

    Returns:
        Hex digest string
    """
    h = hashlib.sha256()
    filepath = Path(filepath)
    with open(filepath, 'rb') as f:
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def sha256_string(s: str) -> str:
    """SHA-256 of a string (for hashing command-line flags, etc.)."""
    return hashlib.sha256(s.encode('utf-8')).hexdigest()


# ---------------------------------------------------------------------------
# Provenance header for output files
# ---------------------------------------------------------------------------

def generate_header_lines(
    step_name: str,
    input_files: Sequence[Union[str, Path]],
    command_string: str,
    version: str = __version__,
    extra: Optional[Dict[str, str]] = None,
    comment_char: str = '#',
) -> List[str]:
    """
    Generate provenance header lines for an output file.

    Args:
        step_name: Pipeline step name (e.g. 'correct', 'aggregate_3prime')
        input_files: Input file paths
        command_string: Full command invocation string
        version: RECTIFY version
        extra: Additional key=value metadata
        comment_char: Line prefix (default '#')

    Returns:
        List of header strings (each already prefixed with comment_char)
    """
    ts = datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
    lines = [
        f"{comment_char} RECTIFY v{version}",
        f"{comment_char} Step: {step_name}",
        f"{comment_char} Timestamp: {ts}",
        f"{comment_char} Command: {command_string}",
    ]

    for i, fp in enumerate(input_files):
        fp = Path(fp)
        if fp.exists():
            h = sha256_file(fp)
            lines.append(f"{comment_char} Input[{i}]: {fp.name}  sha256:{h}")
        else:
            lines.append(f"{comment_char} Input[{i}]: {fp.name}  (not found)")

    if extra:
        for k, v in extra.items():
            lines.append(f"{comment_char} {k}: {v}")

    return lines


# ---------------------------------------------------------------------------
# Sidecar provenance JSON
# ---------------------------------------------------------------------------

def _sidecar_path(output_file: Union[str, Path]) -> Path:
    """Return the .provenance.json sidecar path for a given output file."""
    p = Path(output_file)
    return p.parent / f".{p.name}.provenance.json"


def write_sidecar(
    output_file: Union[str, Path],
    step_name: str,
    input_files: Sequence[Union[str, Path]],
    command_string: str,
    version: str = __version__,
    extra: Optional[Dict[str, Any]] = None,
) -> Path:
    """
    Write a .provenance.json sidecar file next to an output file.

    The sidecar records everything needed to decide whether the output
    is still valid: input hashes, the command/version that produced it,
    and the output's own hash.

    Args:
        output_file: The output file this sidecar describes
        step_name: Pipeline step name
        input_files: Input file paths
        command_string: Full rectify invocation
        version: RECTIFY version string
        extra: Arbitrary extra metadata

    Returns:
        Path to sidecar file
    """
    output_file = Path(output_file)
    sidecar = _sidecar_path(output_file)

    record = {
        'rectify_version': version,
        'step': step_name,
        'timestamp': datetime.now(timezone.utc).isoformat(),
        'command': command_string,
        'python': sys.executable,
        'inputs': {},
        'output': {},
    }

    # Hash every input
    for fp in input_files:
        fp = Path(fp)
        if fp.exists():
            record['inputs'][str(fp)] = {
                'sha256': sha256_file(fp),
                'size': fp.stat().st_size,
                'mtime': fp.stat().st_mtime,
            }
        else:
            record['inputs'][str(fp)] = {'sha256': None, 'missing': True}

    # Hash the output itself
    if output_file.exists():
        record['output'] = {
            'path': str(output_file),
            'sha256': sha256_file(output_file),
            'size': output_file.stat().st_size,
        }

    # Config/flags hash — allows detecting if flags changed even if inputs didn't
    record['command_hash'] = sha256_string(f"{version}|{command_string}")

    if extra:
        record['extra'] = extra

    sidecar.parent.mkdir(parents=True, exist_ok=True)
    with open(sidecar, 'w') as f:
        json.dump(record, f, indent=2)

    return sidecar


def read_sidecar(output_file: Union[str, Path]) -> Optional[Dict]:
    """
    Read a sidecar provenance JSON if it exists.

    Returns:
        Parsed dict, or None if sidecar doesn't exist
    """
    sidecar = _sidecar_path(output_file)
    if not sidecar.exists():
        return None
    try:
        with open(sidecar, 'r') as f:
            return json.load(f)
    except (json.JSONDecodeError, OSError):
        return None


# ---------------------------------------------------------------------------
# Step-skip logic
# ---------------------------------------------------------------------------

def step_is_current(
    output_file: Union[str, Path],
    input_files: Sequence[Union[str, Path]],
    command_string: str,
    version: str = __version__,
    path_translator: Optional[Any] = None,
) -> bool:
    """
    Check whether a step's output is still valid (inputs + command unchanged).

    A step is current if ALL of the following are true:
    1. The output file exists
    2. A sidecar exists for it
    3. Every input file's SHA-256 matches the sidecar record
    4. The command+version hash matches (catches flag changes)

    Args:
        output_file: Primary output file for this step
        input_files: All input files for this step
        command_string: The exact command string (flags included)
        version: RECTIFY version

    Returns:
        True if the step can be skipped
    """
    output_file = Path(output_file)

    if not output_file.exists():
        return False

    record = read_sidecar(output_file)
    if record is None:
        return False

    # Check command + version hash
    expected_cmd_hash = sha256_string(f"{version}|{command_string}")
    if record.get('command_hash') != expected_cmd_hash:
        return False

    # Check each input hash
    recorded_inputs = record.get('inputs', {})
    _xlate = path_translator if callable(path_translator) else (lambda p: p)
    for fp in input_files:
        fp = Path(fp)
        fp_str = str(_xlate(fp))

        if fp_str not in recorded_inputs:
            return False

        rec = recorded_inputs[fp_str]
        if rec.get('missing'):
            return False

        if not fp.exists():
            return False

        current_hash = sha256_file(fp)
        if current_hash != rec.get('sha256'):
            return False

    return True


# ---------------------------------------------------------------------------
# ProvenanceTracker — high-level interface for run-all
# ---------------------------------------------------------------------------

class ProvenanceTracker:
    """
    High-level provenance tracker for a pipeline run.

    Accumulates step records and provides a single pipeline-level
    provenance manifest at the end.

    Usage:
        tracker = ProvenanceTracker(output_dir, command_string)

        if tracker.step_needed('correct', inputs, outputs):
            # ... do work ...
            tracker.record_step('correct', inputs, outputs)
        else:
            tracker.log_skip('correct')

        tracker.write_manifest()  # writes pipeline_provenance.json
    """

    def __init__(
        self,
        output_dir: Union[str, Path],
        command_string: str,
        version: str = __version__,
    ):
        self.output_dir = Path(output_dir)
        self.command_string = command_string
        self.version = version
        self.steps: List[Dict] = []
        self.skipped: List[str] = []
        self.start_time = time.monotonic()

    def step_needed(
        self,
        step_name: str,
        input_files: Sequence[Union[str, Path]],
        output_files: Sequence[Union[str, Path]],
    ) -> bool:
        """
        Check if a step needs to run.

        A step is skipped only if ALL output files are current
        (matching input hashes + command).
        """
        step_cmd = f"{self.command_string} [step={step_name}]"
        for out_f in output_files:
            if not step_is_current(out_f, input_files, step_cmd, self.version):
                return True
        return False

    def record_step(
        self,
        step_name: str,
        input_files: Sequence[Union[str, Path]],
        output_files: Sequence[Union[str, Path]],
        extra: Optional[Dict[str, Any]] = None,
    ) -> None:
        """Record provenance for a completed step (writes sidecars)."""
        step_cmd = f"{self.command_string} [step={step_name}]"

        for out_f in output_files:
            out_f = Path(out_f)
            if out_f.exists():
                write_sidecar(
                    out_f, step_name, input_files,
                    step_cmd, self.version, extra,
                )

        self.steps.append({
            'step': step_name,
            'status': 'completed',
            'inputs': [str(f) for f in input_files],
            'outputs': [str(f) for f in output_files],
            'timestamp': datetime.now(timezone.utc).isoformat(),
        })

    def log_skip(self, step_name: str) -> None:
        """Record that a step was skipped (outputs already current)."""
        self.skipped.append(step_name)
        self.steps.append({
            'step': step_name,
            'status': 'skipped',
            'timestamp': datetime.now(timezone.utc).isoformat(),
        })

    def write_manifest(self) -> Path:
        """
        Write a pipeline-level provenance manifest.

        Returns:
            Path to pipeline_provenance.json
        """
        elapsed = time.monotonic() - self.start_time
        manifest = {
            'rectify_version': self.version,
            'command': self.command_string,
            'started': self.steps[0]['timestamp'] if self.steps else None,
            'finished': datetime.now(timezone.utc).isoformat(),
            'elapsed_seconds': round(elapsed, 1),
            'steps': self.steps,
            'skipped_steps': self.skipped,
        }

        manifest_path = self.output_dir / 'pipeline_provenance.json'
        self.output_dir.mkdir(parents=True, exist_ok=True)
        with open(manifest_path, 'w') as f:
            json.dump(manifest, f, indent=2)

        return manifest_path

    def header_lines(
        self,
        step_name: str,
        input_files: Sequence[Union[str, Path]],
        comment_char: str = '#',
        extra: Optional[Dict[str, str]] = None,
    ) -> List[str]:
        """Generate provenance header lines for embedding in output files."""
        return generate_header_lines(
            step_name=step_name,
            input_files=input_files,
            command_string=self.command_string,
            version=self.version,
            extra=extra,
            comment_char=comment_char,
        )
