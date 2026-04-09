#!/usr/bin/env python3
"""
Provenance tracking for RECTIFY outputs.

Provides automatic generation and updating of provenance information:
- PROVENANCE.json: Machine-readable provenance with file hashes, timestamps, commands
- README.md: Human-readable documentation

Supports:
- Initial provenance generation on first run
- Incremental updates when reprocessing samples
- Command/parameter tracking for reproducibility

Author: Kevin R. Roy
Date: 2026-03-27
"""

import json
import hashlib
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Any
import sys
import os


def _is_relative_to(path: Path, base: Path) -> bool:
    """Python 3.8-compatible replacement for Path.is_relative_to() (added in 3.9)."""
    try:
        path.relative_to(base)
        return True
    except ValueError:
        return False


class ProvenanceTracker:
    """
    Track provenance for RECTIFY output directories.

    Usage:
        tracker = ProvenanceTracker(output_dir)
        tracker.set_command(sys.argv)
        tracker.set_config(config_dict)

        # ... process files ...
        tracker.add_output_file(output_path, source_files=[input1, input2])

        # At end:
        tracker.save()
    """

    def __init__(self, output_dir: Path, description: str = "RECTIFY output"):
        self.output_dir = Path(output_dir)
        self.description = description
        self.prov_file = self.output_dir / "PROVENANCE.json"
        self.readme_file = self.output_dir / "README.md"

        # Load existing provenance if present
        self.provenance = self._load_existing()

        # Track current run
        self.current_run = {
            'timestamp': datetime.now().isoformat(),
            'command': None,
            'config': None,
            'python_version': sys.version,
            'working_dir': os.getcwd(),
            'outputs': [],
        }

    def _load_existing(self) -> Dict:
        """Load existing provenance or create new."""
        if self.prov_file.exists():
            with open(self.prov_file) as f:
                return json.load(f)
        else:
            return {
                'directory': str(self.output_dir),
                'description': self.description,
                'created': datetime.now().isoformat(),
                'runs': [],
                'files': {},
            }

    def set_command(self, argv: List[str]):
        """Record the command that was run."""
        self.current_run['command'] = ' '.join(argv)

    def set_config(self, config: Dict):
        """Record configuration/parameters."""
        # Convert Path objects to strings for JSON serialization
        self.current_run['config'] = self._serialize_config(config)

    def _serialize_config(self, config: Any) -> Any:
        """Recursively convert config to JSON-serializable format."""
        if isinstance(config, dict):
            return {k: self._serialize_config(v) for k, v in config.items()}
        elif isinstance(config, (list, tuple)):
            return [self._serialize_config(v) for v in config]
        elif isinstance(config, Path):
            return str(config)
        elif hasattr(config, '__dict__'):
            return self._serialize_config(vars(config))
        else:
            try:
                json.dumps(config)
                return config
            except:
                return str(config)

    def add_output_file(self, filepath: Path, source_files: Optional[List[Path]] = None,
                        metadata: Optional[Dict] = None):
        """
        Record an output file with its provenance.

        Args:
            filepath: Path to output file
            source_files: List of input files used to generate this output
            metadata: Additional metadata (e.g., sample name, parameters)
        """
        filepath = Path(filepath)
        if not filepath.exists():
            return

        stat = filepath.stat()
        file_info = {
            'path': str(filepath),
            'name': filepath.name,
            'size_bytes': stat.st_size,
            'modified': datetime.fromtimestamp(stat.st_mtime).isoformat(),
            'md5': self._get_hash(filepath) if stat.st_size < 100_000_000 else 'large_file',
            'generated_by_run': len(self.provenance['runs']),  # Index of current run
        }

        if source_files:
            file_info['source_files'] = [
                {'path': str(p), 'md5': self._get_hash(p) if Path(p).exists() else 'not_found'}
                for p in source_files
            ]

        if metadata:
            file_info['metadata'] = metadata

        self.current_run['outputs'].append(file_info)

        # Update master file list
        rel_path = str(filepath.relative_to(self.output_dir)) if _is_relative_to(filepath, self.output_dir) else str(filepath)
        self.provenance['files'][rel_path] = file_info

    def _get_hash(self, filepath: Path, chunk_size: int = 8192) -> str:
        """Calculate MD5 hash of file."""
        filepath = Path(filepath)
        if not filepath.exists() or filepath.stat().st_size > 100_000_000:
            return 'skipped'

        h = hashlib.md5()
        try:
            with open(filepath, 'rb') as f:
                while chunk := f.read(chunk_size):
                    h.update(chunk)
            return h.hexdigest()
        except:
            return 'error'

    def save(self):
        """Save provenance to output directory."""
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Add current run to history
        self.provenance['runs'].append(self.current_run)
        self.provenance['last_updated'] = datetime.now().isoformat()

        # Write PROVENANCE.json atomically: write to .tmp then rename so a
        # crash mid-write never leaves a corrupt sidecar on disk.
        tmp_file = self.prov_file.with_suffix('.json.tmp')
        with open(tmp_file, 'w') as f:
            json.dump(self.provenance, f, indent=2)
        os.replace(tmp_file, self.prov_file)

        # Write README.md
        self._write_readme()

    def _write_readme(self):
        """Generate human-readable README."""
        with open(self.readme_file, 'w') as f:
            f.write(f"# {self.output_dir.name}\n\n")
            f.write(f"**Description:** {self.description}\n\n")
            f.write(f"**Created:** {self.provenance.get('created', 'Unknown')}\n")
            f.write(f"**Last updated:** {self.provenance.get('last_updated', 'Unknown')}\n")
            f.write(f"**Total runs:** {len(self.provenance['runs'])}\n\n")

            # File summary
            files = self.provenance.get('files', {})
            if files:
                f.write(f"## Files ({len(files)})\n\n")
                f.write("| File | Size | Modified |\n")
                f.write("|------|------|----------|\n")
                for name, info in sorted(files.items())[:50]:
                    size = self._format_size(info.get('size_bytes', 0))
                    mod = info.get('modified', '')[:10]
                    f.write(f"| {name} | {size} | {mod} |\n")
                if len(files) > 50:
                    f.write(f"\n*... and {len(files) - 50} more files*\n")
                f.write("\n")

            # Run history (last 5)
            runs = self.provenance.get('runs', [])
            if runs:
                f.write("## Recent Runs\n\n")
                for i, run in enumerate(runs[-5:]):
                    run_num = len(runs) - 5 + i if len(runs) > 5 else i
                    f.write(f"### Run {run_num + 1} ({run.get('timestamp', 'Unknown')[:10]})\n\n")
                    if run.get('command'):
                        f.write(f"```\n{run['command']}\n```\n\n")
                    if run.get('outputs'):
                        f.write(f"Outputs: {len(run['outputs'])} files\n\n")

            f.write("---\n")
            f.write("*Full provenance in `PROVENANCE.json`*\n")

    @staticmethod
    def _format_size(size_bytes: int) -> str:
        """Format bytes as human-readable."""
        for unit in ['B', 'KB', 'MB', 'GB']:
            if size_bytes < 1024:
                return f"{size_bytes:.1f} {unit}"
            size_bytes /= 1024
        return f"{size_bytes:.1f} TB"


def init_provenance(output_dir: Path, description: str = "RECTIFY output",
                    command: Optional[List[str]] = None,
                    config: Optional[Dict] = None) -> ProvenanceTracker:
    """
    Initialize provenance tracking for an output directory.

    Convenience function for common usage pattern.

    Args:
        output_dir: Output directory path
        description: Description of this output
        command: Command line arguments (default: sys.argv)
        config: Configuration dictionary

    Returns:
        ProvenanceTracker instance
    """
    tracker = ProvenanceTracker(output_dir, description)
    tracker.set_command(command or sys.argv)
    if config:
        tracker.set_config(config)
    return tracker
