#!/bin/bash
#$ -N cdna_cal
#$ -cwd
#$ -j y
#$ -pe shared 8
#$ -l h_data=8G,h_rt=24:00:00
#$ -o /u/scratch/k/kevinroy/calibration_2026_05/cdna/slurm_logs/cdna_cal_$JOB_ID.log
#
# ONT PCR-cDNA (UMI-consensus) penalty-table calibration.
# Driver: scripts/calibration/empirical_cigar_error_profiler.py --protocol cdna
#
# Input: stage1_consensus.bam (UMI-collapsed cDNA reads, median ~1.9 kb) from
# /u/scratch/k/kevinroy/ont_cdna_v4/wt_rep2/out/. NOTE: this is a CONSENSUS BAM,
# not raw cDNA — the resulting tables are calibrated for consensus-grade inputs.
# Raw (pre-UMI-collapse) cDNA will have higher per-base error than this table
# expects; rectify should fall back to the heuristic when no UMI step ran.
#
# Wall budget: BAM->FASTQ ~10 min, split <5 min, align ~4-6 h serial-per-chunk
# even with --parallel-aligners (deSALT/uLTRA dominate), profile ~2-3 h.
# h_rt=24h is the safety margin; if alignment is the bottleneck convert the
# per-chunk loop to a `qsub -t 1-16` job array (one chunk per array task).
#
# Submit:   qsub run_profiler_cdna_h2.sh
# Monitor:  qstat -j <JOB_ID>  ;  tail -f <slurm_logs>/cdna_cal_<JOB_ID>.log
set -euo pipefail

# --- Env --------------------------------------------------------------------
# H2 batch shells need the modules init explicitly sourced; the default
# /etc/profile.d/modules.sh path does NOT exist on H2 compute nodes.
. /u/local/Modules/default/init/bash
module load conda/23.11.0
source activate rectify
module load samtools 2>/dev/null || true   # samtools may be on the env's bin

# --- Paths ------------------------------------------------------------------
WORK=/u/scratch/k/kevinroy/calibration_2026_05/cdna
CONS_BAM=/u/scratch/k/kevinroy/ont_cdna_v4/wt_rep2/out/stage1_consensus.bam
PROFILER=/u/scratch/k/kevinroy/calibration_2026_05/scripts/empirical_cigar_error_profiler.py
# Bundled reference: H2 ships only .fsa.gz (BGZF) with .fai + .gzi indices.
REF=/u/home/k/kevinroy/software/rectify/rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz
test -f "$REF" || { echo "ERROR: reference not found: $REF" >&2; exit 1; }
test -f "$REF.fai" -a -f "$REF.gzi" || { echo "ERROR: missing .fai/.gzi for $REF" >&2; exit 1; }
test -f "$CONS_BAM" || { echo "ERROR: consensus BAM not found: $CONS_BAM" >&2; exit 1; }

echo "=== cDNA calibration run ==="
echo "host:        $(hostname)"
echo "start:       $(date)"
echo "work:        $WORK"
echo "consensus:   $CONS_BAM ($(numfmt --to=iec $(stat -c %s "$CONS_BAM")))"
echo "profiler:    $PROFILER"
echo "ref:         $REF"

# --- Step 1: BAM -> FASTQ --------------------------------------------------
FQ="$WORK/wt_rep2_consensus.fastq.gz"
if [ -f "$FQ" ]; then
    echo "  [skip fastq] $FQ exists"
else
    echo "  [fastq] $CONS_BAM -> $FQ"
    samtools fastq -@ 4 "$CONS_BAM" 2>/dev/null | gzip -1 > "$FQ"
fi

# --- Step 2: split into 16 chunks ------------------------------------------
CHUNKS="$WORK/chunks"
if compgen -G "$CHUNKS/*chunk*.fastq.gz" >/dev/null; then
    echo "  [skip split] $CHUNKS already populated"
else
    echo "  [split] -n 16"
    mkdir -p "$CHUNKS"
    rectify split "$FQ" -n 16 -o "$CHUNKS"
fi

# --- Step 3: align per chunk with 5-aligner long-read panel ----------------
# rectify align writes BAMs to a flat dir; we move them into per-aligner subdirs
# to match the profiler's find_chunk_bams layout expectation.
ALIGN_DIR="$CHUNKS/align"
for a in minimap2 mapPacBio gapmm2 uLTRA deSALT; do
    mkdir -p "$ALIGN_DIR/$a"
done
for chunk_fq in "$CHUNKS"/*chunk*.fastq.gz; do
    chunk_id=$(basename "$chunk_fq" .fastq.gz)
    missing=0
    for a in minimap2 mapPacBio gapmm2 uLTRA deSALT; do
        [ -f "$ALIGN_DIR/$a/${chunk_id}.$a.bam" ] || missing=$((missing+1))
    done
    if [ "$missing" -eq 0 ]; then
        echo "  [skip align] $chunk_id (all 5 aligners present)"
        continue
    fi
    echo "  [align] $chunk_id ($missing aligners missing)"
    # uLTRA and deSALT live under --junction-aligners (they need an annotation);
    # --parallel-aligners runs the primaries concurrently inside one rectify
    # align call. With 8 H2 cores and 5 aligners that's ~2 cores each but well
    # below saturation.
    rectify align "$chunk_fq" \
        --aligners minimap2 mapPacBio gapmm2 \
        --junction-aligners uLTRA deSALT \
        --parallel-aligners \
        --no-consensus \
        --Scer \
        -o "$ALIGN_DIR" \
        --threads 8
    # Reorganize: <chunk_id>.<aligner>.bam[.bai] → align/<aligner>/
    for a in minimap2 mapPacBio gapmm2 uLTRA deSALT; do
        for f in "$ALIGN_DIR/${chunk_id}.${a}.bam" "$ALIGN_DIR/${chunk_id}.${a}.bam.bai"; do
            [ -e "$f" ] && mv "$f" "$ALIGN_DIR/$a/"
        done
    done
    rm -f "$ALIGN_DIR/${chunk_id}.rectified.bam" \
          "$ALIGN_DIR/${chunk_id}.rectified.bam.bai" \
          "$ALIGN_DIR/${chunk_id}_junctions.bed"
done

# --- Step 4: profiler ------------------------------------------------------
OUT="$WORK/error_profile_cdna_$(date +%Y%m%d)"
mkdir -p "$OUT"
echo "  [profile] -> $OUT"
python "$PROFILER" \
    --protocol cdna \
    --run-dir "$CHUNKS" \
    --reference "$REF" \
    --output-dir "$OUT" \
    --isolation-flank 10 \
    --union \
    --str-repeat \
    --workers 8 \
    2>&1 | tee "$OUT/profiler.log"

# --- Step 5: overhang calibrator -------------------------------------------
OVH_OUT="$WORK/overhang_table_cdna_$(date +%Y%m%d)"
mkdir -p "$OVH_OUT"
echo "  [overhang] -> $OVH_OUT"
python -m rectify.core.splice.calibrate_junction_overhang \
    --run-dir "$CHUNKS" \
    --reference "$REF" \
    --output-dir "$OVH_OUT" \
    2>&1 | tee "$OVH_OUT/calibrate.log"

echo "=== Done ==="
echo "end:      $(date)"
echo "profile:  $OUT/penalty_scores.tsv  +  str_penalty_scores.tsv"
echo "overhang: $OVH_OUT/junction_overhang_table.tsv"
echo
echo "Next steps (local M1):"
echo "  1. rsync $OUT $OVH_OUT M1:/path/to/rectify/scripts/calibration/runs/"
echo "  2. python scripts/calibration/bundle_protocol_tables.py \\"
echo "         --protocol cdna \\"
echo "         --inputs runs/error_profile_cdna_*/ \\"
echo "         --overhang runs/overhang_table_cdna_*/junction_overhang_table.tsv"
echo "  3. pytest tests/test_data_bundling.py tests/test_validation_reads_cdna.py"
