#!/usr/bin/env bash
# Launch IGV cleanly for the validation-read IGV review.
#
# Usage:
#   bash scripts/validation_data/launch_igv_review.sh <batch_file>
#
# What this script does:
#   1. Kills any running IGV instance (so the command port is free)
#   2. Re-stashes any autosaved sessions that may have accumulated since
#      AUTOLOAD_LAST_AUTOSAVE was disabled (belt and suspenders)
#   3. Launches IGV with the requested batch script
#
# Prerequisites (one-time, already done):
#   - AUTOLOAD_LAST_AUTOSAVE=false in ~/igv/prefs.properties (confirmed)
#   - Existing autosaves moved to ~/igv/autosave.disabled.<timestamp>/

set -euo pipefail

BATCH=${1:-}
if [[ -z "$BATCH" ]]; then
  echo "usage: $0 <batch_file>" >&2
  exit 1
fi
if [[ ! -f "$BATCH" ]]; then
  echo "batch file not found: $BATCH" >&2
  exit 1
fi

IGV_BIN=${IGV_BIN:-/Applications/IGV_2.19.2.app/Contents/MacOS/igv}
if [[ ! -x "$IGV_BIN" ]]; then
  echo "IGV not found at $IGV_BIN (set IGV_BIN env var to override)" >&2
  exit 1
fi

# 1. Kill any running IGV
pkill -f "IGV_" 2>/dev/null || true
pkill -f "igv.jar" 2>/dev/null || true
# Wait for processes to exit and port 60151 to free up
for _ in 1 2 3 4 5; do
  if ! lsof -i :60151 >/dev/null 2>&1; then
    break
  fi
  sleep 0.5
done

# 2. Re-stash any new autosaves
mkdir -p ~/igv/autosave.stash
mv ~/igv/autosave/*.xml ~/igv/autosave.stash/ 2>/dev/null || true

# 3. Launch with the batch file (suppress warnings)
"$IGV_BIN" -b "$BATCH" 2>&1 | grep -vE "WARNING|^openjdk|^OpenJDK|^SLF4J" &

echo "Launched IGV with batch: $BATCH"
echo "PID: $!"
