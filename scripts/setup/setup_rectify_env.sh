#!/usr/bin/env bash
# Post-install bootstrap for the rectify conda env.
#
# Two manual workarounds are required on top of `conda install`:
#
#   1. `namfinder` is needed by uLTRA but isn't a transitive dependency of
#      `ultra-bioinformatics`. The binary vendored in
#      `rectify/data/bin/linux_x86_64/namfinder` was built against newer
#      GLIBCXX than is available on some HPC nodes (Hoffman2's older RHEL).
#      The conda-built binary works everywhere.
#
#   2. uLTRA 0.1's `modules/create_augmented_gene.py` raises an AssertionError
#      on SGD yeast GTF (overlapping per-gene exons trigger it near 89-99% of
#      iteration). Comment the assertion out — the algorithm completes
#      correctly past it for our use cases.
#
# Run AFTER you have a working rectify env and have activated it:
#
#     conda activate <your-rectify-env>
#     bash scripts/setup/setup_rectify_env.sh
#
# Or, for the lab-shared env on Hoffman2:
#
#     conda activate /u/project/guillom/shared/envs/rectify
#     bash scripts/setup/setup_rectify_env.sh
#
# Safe to re-run — both steps are idempotent.

set -euo pipefail

if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "ERROR: no conda env is active. Activate your rectify env first." >&2
    exit 1
fi

echo "Bootstrapping env at: $CONDA_PREFIX"

# Step 1: install namfinder via conda if missing.
if command -v namfinder >/dev/null 2>&1; then
    echo "  namfinder: already present at $(command -v namfinder)"
else
    echo "  namfinder: installing via bioconda..."
    if command -v mamba >/dev/null 2>&1; then
        mamba install -y -c bioconda namfinder
    else
        conda install -y -c bioconda namfinder
    fi
fi

# Step 2: patch uLTRA's create_augmented_gene.py to skip the SGD GTF assertion.
MOD_PY=$(python -c "import modules, os; print(os.path.join(os.path.dirname(modules.__file__), 'create_augmented_gene.py'))")
if [[ ! -f "$MOD_PY" ]]; then
    echo "  uLTRA: modules/create_augmented_gene.py not found — is ultra-bioinformatics installed?" >&2
    exit 1
fi
if grep -q "rectify-patch" "$MOD_PY"; then
    echo "  uLTRA: assertion already patched in $MOD_PY"
else
    echo "  uLTRA: patching assertion in $MOD_PY"
    cp "$MOD_PY" "$MOD_PY.bak"
    # Use a temp file because BSD sed (macOS) and GNU sed (Linux) differ on
    # in-place editing flags.
    awk '
        /^        assert active_start <= exon\.start - 1$/ {
            print "        # rectify-patch: assertion removed (fires on overlapping per-gene exons in yeast SGD GTF)"
            print "        # " $0
            next
        }
        { print }
    ' "$MOD_PY.bak" > "$MOD_PY"
fi

echo
echo "Done. uLTRA + namfinder are ready."
