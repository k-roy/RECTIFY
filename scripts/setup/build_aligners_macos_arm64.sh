#!/usr/bin/env bash
# Build missing aligners from source on Apple Silicon (arm64).
#
# Why this script exists:
#   Apple Silicon Macs lack pre-built binaries for two of the five rectify
#   aligners. This script builds them with the macOS-specific patches needed:
#
#     1. namfinder (uLTRA dependency) — no osx-arm64 bioconda package.
#        Builds native arm64 from source via cmake. Needs the
#        CPLUS_INCLUDE_PATH workaround because some CLT installs leave
#        `/Library/Developer/CommandLineTools/usr/include/c++/v1/` empty
#        while the actual C++ headers live under the SDK.
#
#     2. deBGA + deSALT — bioconda's osx-arm64 deBGA ships a Linux ELF
#        (broken on macOS). Cross-compile both to x86_64 Mach-O via
#        `-arch x86_64` so they run under Rosetta. Three patches:
#          (a) deBGA uses `<omp.h>` only for `omp_get_wtime()` — replace
#              with a `clock_gettime` shim so no libomp linkage needed.
#          (b) deSALT's `desalt_index.c` calls `readlink("/proc/self/exe")`
#              to locate the deBGA binary — Linux-only. Patch to use
#              `_NSGetExecutablePath()` on macOS.
#          (c) deBGA's Makefile pattern rule swallows CFLAGS; the patch
#              adds `$(CFLAGS)` to the .c.o rule.
#
# Usage (from any cwd):
#   bash scripts/setup/build_aligners_macos_arm64.sh
#
# Output: ~/miniconda3/bin/namfinder, ~/miniconda3/envs/desalt/bin/{deBGA,deSALT}
# (and ~/.rectify/bin/deBGA symlink so deSALT can find it via its own
# argv[0]-based path discovery).
#
# Safe to re-run. Skips steps where the target binary already exists and
# is the expected arch.

set -euo pipefail

SOFTWARE_DIR="$HOME/software"
mkdir -p "$SOFTWARE_DIR"

# Locate the SDK with C++ headers. CLT install sometimes leaves the
# default v1/ path empty; the SDK's copy is reliable.
SDK="$(xcrun --show-sdk-path)"
CPP_HEADERS="$SDK/usr/include/c++/v1"
if [[ ! -f "$CPP_HEADERS/string" ]]; then
    echo "ERROR: C++ headers not found at $CPP_HEADERS" >&2
    echo "       Run: softwareupdate --install --all" >&2
    exit 1
fi

echo "Using SDK: $SDK"
echo "Using C++ headers: $CPP_HEADERS"

# ──────────────────────────────────────────────────────────────────────
# 1) namfinder (native arm64)
# ──────────────────────────────────────────────────────────────────────
NAMFINDER_BIN="$HOME/miniconda3/bin/namfinder"
if [[ -x "$NAMFINDER_BIN" ]] && "$NAMFINDER_BIN" --help >/dev/null 2>&1; then
    echo "namfinder: already installed at $NAMFINDER_BIN"
else
    echo "namfinder: building from source..."
    cd "$SOFTWARE_DIR"
    if [[ ! -d namfinder ]]; then
        git clone https://github.com/ksahlin/namfinder
    fi
    cd namfinder
    rm -rf build
    CPLUS_INCLUDE_PATH="$CPP_HEADERS" SDKROOT="$SDK" \
        cmake -B build \
            -DCMAKE_OSX_SYSROOT="$SDK" \
            -DCMAKE_C_FLAGS="-O3" \
            -DCMAKE_CXX_FLAGS="-O3"
    make -C build -j 2
    mkdir -p "$(dirname "$NAMFINDER_BIN")"
    cp build/namfinder "$NAMFINDER_BIN"
    echo "namfinder: installed -> $NAMFINDER_BIN"
fi

# ──────────────────────────────────────────────────────────────────────
# 2) deBGA + deSALT (cross-compile to x86_64 Mach-O, runs via Rosetta)
# ──────────────────────────────────────────────────────────────────────
DESALT_ENV_BIN="$HOME/miniconda3/envs/desalt/bin"
mkdir -p "$DESALT_ENV_BIN"

# 2a) deBGA
cd "$SOFTWARE_DIR"
if [[ ! -d deBGA ]]; then
    git clone --depth 1 https://github.com/HongzheGuo/deBGA.git
fi
cd deBGA

# Patch (a): replace <omp.h> + omp_get_wtime with clock_gettime shim.
if ! grep -q "omp.h removed for portability" seed_ali_core.c; then
    python3 - <<'PY'
import re
path = 'seed_ali_core.c'
s = open(path).read()
shim = (
    "/* omp.h removed for portability; omp_get_wtime replaced */\n"
    "#include <time.h>\n"
    "static inline double omp_get_wtime(void) {\n"
    "    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);\n"
    "    return ts.tv_sec + ts.tv_nsec * 1e-9;\n"
    "}"
)
s2 = re.sub(r'^#include\s*<omp\.h>\s*$', shim, s, count=1, flags=re.MULTILINE)
open(path, 'w').write(s2)
PY
fi

# Patch (c): Makefile pattern rule needs $(CFLAGS).
if ! grep -q 'CFLAGS) \$(INCLUDES)' Makefile; then
    sed -i.bak 's|\$(CC) -c \$(INCLUDES) \$< -o \$@|\$(CC) -c \$(CFLAGS) \$(INCLUDES) \$< -o \$@|' Makefile
fi

make clean >/dev/null 2>&1 || true
make CC="clang -arch x86_64" \
     CFLAGS="-O2 -isysroot $SDK -Wno-implicit-function-declaration -Wno-deprecated-non-prototype -ferror-limit=20" \
     LIBS="-arch x86_64 -lz -lpthread -lm" \
     -j 2
file ./deBGA | grep -q "Mach-O 64-bit executable x86_64" || {
    echo "ERROR: deBGA build did not produce expected Mach-O x86_64 binary" >&2
    exit 1
}
cp ./deBGA "$DESALT_ENV_BIN/deBGA"
echo "deBGA: installed -> $DESALT_ENV_BIN/deBGA"

# 2b) deSALT
cd "$SOFTWARE_DIR"
if [[ ! -d deSALT ]]; then
    git clone --depth 1 https://github.com/ydLiu-HIT/deSALT.git
fi
cd deSALT/src

# Patch (b): readlink("/proc/self/exe") -> _NSGetExecutablePath() on macOS.
if ! grep -q "_NSGetExecutablePath" desalt_index.c; then
    python3 - <<'PY'
path = 'desalt_index.c'
s = open(path).read()
header_addition = (
    "\n"
    "#if defined(__APPLE__)\n"
    "#include <mach-o/dyld.h>\n"
    "static int desalt_self_exe_path(char *buf, size_t bufsize) {\n"
    "    uint32_t size = bufsize;\n"
    "    if (_NSGetExecutablePath(buf, &size) != 0) return -1;\n"
    "    return (int)strlen(buf);\n"
    "}\n"
    "#else\n"
    "#include <unistd.h>\n"
    "static int desalt_self_exe_path(char *buf, size_t bufsize) {\n"
    "    int r = readlink(\"/proc/self/exe\", buf, bufsize);\n"
    "    if (r < 0 || (size_t)r >= bufsize) return -1;\n"
    "    buf[r] = '\\0';\n"
    "    return r;\n"
    "}\n"
    "#endif\n"
)
s = s.replace("int deBGA_index(", header_addition + "int deBGA_index(", 1)
s = s.replace(
    'r = readlink("/proc/self/exe", dir, 2048);\n'
    '    if (r < 0 || r >= 2048)\n'
    '    {\n'
    '        fprintf(stderr, "Failed, could not find the program of deSALT!\\n");\n'
    '    }\n'
    '    dir[r] = \'\\0\';',
    'r = desalt_self_exe_path(dir, sizeof(dir));\n'
    '    if (r < 0)\n'
    '    {\n'
    '        fprintf(stderr, "Failed, could not find the program of deSALT!\\n");\n'
    '        exit(1);\n'
    '    }'
)
open(path, 'w').write(s)
PY
fi

make clean >/dev/null 2>&1 || true
make CC="clang -arch x86_64" \
     CFLAGS="-O3 -isysroot $SDK -Wno-implicit-function-declaration -Wno-deprecated-non-prototype -Wno-unused-but-set-variable -Wno-unused-function -ferror-limit=20 -DHAVE_KALLOC" \
     LIBS="-arch x86_64 -lm -lz -lpthread" \
     -j 2
file ./deSALT | grep -q "Mach-O 64-bit executable x86_64" || {
    echo "ERROR: deSALT build did not produce expected Mach-O x86_64 binary" >&2
    exit 1
}
cp ./deSALT "$DESALT_ENV_BIN/deSALT"
echo "deSALT: installed -> $DESALT_ENV_BIN/deSALT"

# ──────────────────────────────────────────────────────────────────────
# 3) Symlinks for rectify's expected layout
# ──────────────────────────────────────────────────────────────────────
RECTIFY_BIN="$HOME/.rectify/bin"
mkdir -p "$RECTIFY_BIN"
ln -sf "$DESALT_ENV_BIN/deSALT" "$RECTIFY_BIN/deSALT"
ln -sf "$DESALT_ENV_BIN/deBGA"  "$RECTIFY_BIN/deBGA"
echo "Symlinks: $RECTIFY_BIN/{deSALT,deBGA}"

# ──────────────────────────────────────────────────────────────────────
# 4) Final smoke
# ──────────────────────────────────────────────────────────────────────
echo
echo "=== verification ==="
"$NAMFINDER_BIN" --help 2>&1 | head -1
"$DESALT_ENV_BIN/deBGA" 2>&1 | grep -E "^Program|^Version" | head -2
"$DESALT_ENV_BIN/deSALT" 2>&1 | grep -E "^Program|^Version" | head -2
echo
echo "Done. M1 should now be ready for 5/5 rectify align runs."
echo "Don't forget to also run scripts/setup/setup_rectify_env.sh once"
echo "(applies uLTRA assertion patch + any namfinder fallback)."
