#!/bin/bash
# =============================================================================
# SIRV absolute-truth error-structure measurement  (HANDOFF RESUME #1, the SIRV gate)
# -----------------------------------------------------------------------------
# WHAT: align real ONT DRS reads (WTC11, ENCODE) to the SIRV-Set 4 spike-in
#       reference, keep the reads that map (= the SIRV spike-ins, whose exact
#       sequence is KNOWN -> absolute truth, no RNA-mod / no biological-variant
#       confound), and measure their error structure with the SAME
#       `measure_error_structure` code used on real-yeast / pbsim / injector.
# WHY:  the magnitude of the real error clustering (sub5_gap_excess ~3 vs ~5) and
#       per-read over-dispersion (overdisp_v) is currently an UPPER BOUND from
#       read-vs-genome (conflates RNA-mod + minimap2 junction pile-up). SIRV is
#       junction-spanning AND mod-free -> it DECIDES the clean magnitude target
#       the injector re-fit needs, and re-checks guard #3 (alignment inflation)
#       in the junction-spanning case the gap-free test could not.
#
# OUTPUTS (under $WORK): sirv.sorted.bam(.bai) + measure-bam / autocorr-bam stdout
#       captured to sirv_errstruct.report.txt, and a .sirv_rc sentinel (rc inside).
#
# HEADER CAVEAT: confirm the SLURM header (account / partition) against the
#       `sherlock-sbatch` skill before submitting — larsms partition naming and
#       AVX-512 exclude are the documented Sherlock traps. The science body below
#       is the reviewed part.
# =============================================================================
#SBATCH --job-name=sirv_errstruct
#SBATCH --partition=normal              # owners (larsms) is faster if available
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/scratch/users/kevinroy/sirv_work/sirv_errstruct.%j.out
#SBATCH --error=/scratch/users/kevinroy/sirv_work/sirv_errstruct.%j.err
# AVX-512 SIGILL guard (documented trap; harmless if not honored):
#SBATCH --constraint="NO_GPU"

set -eo pipefail

# ---- paths -----------------------------------------------------------------
D=/home/groups/larsms/users/kevinroy/aligner_bench_live      # synced bench code
CB=/home/groups/larsms/users/kevinroy/anaconda3
P3=$CB/envs/pbsim3                                           # minimap2 + samtools
PY=$CB/envs/rectify/bin/python                               # full deps (import rectify)
WORK=/scratch/users/kevinroy/sirv_work
mkdir -p "$WORK"
cd "$WORK"

MM2=$P3/bin/minimap2
SAMTOOLS=$P3/bin/samtools

# ---- conda (activate BEFORE `set -u`: the profile references unbound vars) --
source "$CB/etc/profile.d/conda.sh"
conda activate "$CB/envs/rectify"
set -u

THREADS=${SLURM_CPUS_PER_TASK:-8}

# ---- 1. SIRV-Set 4 reference (tiny; CGL open mirror) ------------------------
REF_BASE=https://cgl.gi.ucsc.edu/data/LRGASP/data/references
if [[ ! -s sirv4.fasta ]]; then
  curl -fsSL -o sirv4.fasta.gz "$REF_BASE/lrgasp_sirv4.fasta.gz"
  gunzip -kf sirv4.fasta.gz
fi
if [[ ! -s sirv4.gtf ]]; then
  curl -fsSL -o sirv4.gtf.gz "$REF_BASE/lrgasp_sirvs4.gtf.gz"
  gunzip -kf sirv4.gtf.gz
fi

# ---- 2. real ONT DRS reads (WTC11, ENCODE; no auth) ------------------------
# One rep (~1.3 GiB) is plenty for error structure (yeast used ~30k reads;
# a DRS rep carries tens of thousands of SIRV spike-in reads). Add reps later.
REP=ENCFF155CFF                              # ENCSR392BGY rep 1 (DRS)
if [[ ! -s ${REP}.fastq.gz ]]; then
  curl -fsSL -o ${REP}.fastq.gz \
    "https://www.encodeproject.org/files/${REP}/@@download/${REP}.fastq.gz"
fi

# ---- 3. align DRS reads to the SIRV reference ------------------------------
# DRS native-RNA minimap2 flags (docs/aligners/minimap2.md): -ax splice -uf -k14.
# --MD is MANDATORY: events_from_bam_read derives substitutions from the MD tag
#   and RAISES if it is absent. U->T fixup (ENCODE DRS reads carry uracils;
#   minimap2 would treat U as N and fail to map) on sequence lines only (NR%4==2).
if [[ ! -s sirv.sorted.bam ]]; then
  zcat ${REP}.fastq.gz \
    | awk 'NR%4==2{gsub(/U/,"T");gsub(/u/,"t")}1' \
    | "$MM2" -ax splice -uf -k14 --MD --secondary=no -t "$THREADS" sirv4.fasta - \
    | "$SAMTOOLS" sort -@ "$THREADS" -o sirv.sorted.bam -
  "$SAMTOOLS" index sirv.sorted.bam
fi

MAPPED=$("$SAMTOOLS" view -c -F 0x904 sirv.sorted.bam)   # primary, mapped
echo "[sirv] primary-mapped (spike-in) reads = $MAPPED"

# ---- 4. measure the CLEAN error structure ----------------------------------
# Length-window [600,1000] matches the real-yeast measurement window so
# overdisp_v / sub5_gap_excess are length-comparable across sources (guard #2).
REPORT="$WORK/sirv_errstruct.report.txt"
VAL="$D/scripts/benchmark/error_realism_validate.py"
set +e                                  # capture rc even on failure (sentinel branch)
{
  echo "==== SIRV absolute-truth error structure ($(date -u +%FT%TZ)) ===="
  echo "reads file: $REP  mapped(primary)=$MAPPED  ref=sirv4.fasta"
  echo
  echo "---- measure-bam (full) ----"
  PYTHONPATH=$D "$PY" "$VAL" measure-bam sirv.sorted.bam
  echo
  echo "---- measure-bam (length-windowed [600,1000], matches yeast window) ----"
  PYTHONPATH=$D "$PY" "$VAL" measure-bam sirv.sorted.bam --len-min 600 --len-max 1000
  echo
  echo "---- autocorr-bam (PI-#2 global/local split; frac 0.3) ----"
  PYTHONPATH=$D "$PY" "$VAL" autocorr-bam sirv.sorted.bam --len-min 600 --len-max 1000 --frac 0.3
} 2>&1 | tee "$REPORT"
RC=${PIPESTATUS[0]}                      # rc of the brace group, not tee

# ---- 5. sentinel (rc inside, atomic) ---------------------------------------
printf '%s\n' "$RC" > "$WORK/.sirv_rc.tmp" && mv "$WORK/.sirv_rc.tmp" "$WORK/.sirv_rc"
echo "[sirv] DONE rc=$RC report=$REPORT"
exit $RC

# =============================================================================
# READING THE RESULT (what the re-fit needs):
#   overdisp_v          -> the real per-read over-dispersion target (yeast 0.70)
#   sub5_gap_excess     -> the CLEAN clustering magnitude (decides ~3 vs ~5;
#                          do NOT pre-bias toward either)
#   indel_run_frac_ge2  -> fat-tail indel-run target (yeast 0.39)
#   autocorr pearson_r  -> the global-vs-local split (yeast ~0.30)
# FOLLOW-ONS (separate steps, not this job):
#   * error TYPE split (frac_sub/ins/del, DRS deletion-dominant): run
#     empirical_cigar_error_profiler.py --isolation-flank 0 on sirv.sorted.bam.
#   * build_panel_from_gtf(sirv4.gtf, genome) -> per-junction recall/FDR on the
#     SAME SIRV BAM (junction-spanning guard-#3 cross-check).
#   * re-fit error_injector.calibrate_params to these SIRV targets.
# =============================================================================
