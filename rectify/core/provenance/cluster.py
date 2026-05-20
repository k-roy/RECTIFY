"""Cluster detection. Used by sidecar.write to record ``cluster``; used by
skip_check to warn on cross-cluster mismatch.
"""

import os
import platform
import re
import socket
from typing import Literal

Cluster = Literal["sherlock", "hoffman2", "m1", "other"]


def detect_current_cluster() -> Cluster:
    """Heuristic detection of which cluster we're running on.

    Detection order matters: check more-specific markers first.

    Returns:
        One of ``"sherlock"``, ``"hoffman2"``, ``"m1"``, or ``"other"``.
    """
    # Sherlock: $SHERLOCK env var, or hostname matches sh<N>-<N>[<NN>]
    if os.environ.get("SHERLOCK"):
        return "sherlock"
    hostname = socket.gethostname().lower()
    if hostname.startswith("sh") and "stanford" in hostname:
        return "sherlock"
    # Some Sherlock compute nodes have hostnames like 'sh03-07n10' without
    # the .stanford suffix in $HOSTNAME — pattern-match.
    if re.match(r"^sh\d+-\d+", hostname):
        return "sherlock"

    # Hoffman2: hostname has .hoffman2.idre.ucla.edu suffix, or matches
    # n<NNNN> / login<N> compute/login nodes.
    if "hoffman2" in hostname or hostname.endswith(".idre.ucla.edu"):
        return "hoffman2"
    if re.match(r"^n\d{3,5}$", hostname) or re.match(r"^login\d+$", hostname):
        # These could be H2 OR another cluster; require an additional marker.
        # H2 has /u/project — check for it.
        if os.path.isdir("/u/project"):
            return "hoffman2"

    # M1 laptop: Darwin + arm64.
    # Python running under Rosetta 2 reports machine="x86_64" even on Apple
    # Silicon. Check sysctl hw.optional.arm64 as the authoritative indicator
    # (value "1" means the CPU is Apple Silicon, regardless of emulation layer).
    uname = platform.uname()
    if uname.system == "Darwin":
        if uname.machine == "arm64":
            return "m1"
        # Rosetta path: subprocess sysctl is too heavy; use ctypes-free approach
        # by asking sysctl via Python's ctypes binding or subprocess.
        # We keep it import-free with a stdlib-only call:
        import subprocess
        try:
            result = subprocess.run(
                ["sysctl", "-n", "hw.optional.arm64"],
                capture_output=True,
                text=True,
                timeout=2,
            )
            if result.returncode == 0 and result.stdout.strip() == "1":
                return "m1"
        except Exception:
            pass

    return "other"
