#!/usr/bin/env bash
# Load a DRS validation read into a running IGV session via the HTTP command port.
# Loads the THREE BAM versions of the read alongside the bedgraph + annotation:
#
#   1. mapped       (raw aligner output, pre-correction, pA-trimmed)
#   2. corrected    (post-walkback, still pA-trimmed)
#   3. pA-restored  (post-walkback, original poly-A added back as 3' soft-clip)
#
# IGV HTTP API supports only /load and /goto — no /new, /genome, or /snapshot.
# Old tracks accumulate; remove via right-click panel → Remove All Tracks when
# the view gets cluttered.
#
# Prerequisites:
#   - IGV running and listening on port 60151
#     (View → Preferences → Advanced → Enable Port)
#
# Usage:
#   bash scripts/validation_data/load_drs_read.sh <chrom> <start> <end> [aligner]
#
#   <aligner> defaults to minimap2.

set -euo pipefail

CHROM=${1:?chrom required, e.g. chrXIV}
START=${2:?start required, e.g. 10400}
END=${3:?end required, e.g. 10650}
ALIGNER=${4:-minimap2}

GENOME="/Users/kevinroy/work/rectify/rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz"
GFF="/Users/kevinroy/work/rectify/rectify/data/genomes/saccharomyces_cerevisiae/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz"
MAPPED_BAM="/Users/kevinroy/work/rectify/rectify/data/validation/aligners/validation_reads.${ALIGNER}.bam"
CORRECTED_BAM="/Users/kevinroy/work/rectify/rectify/data/validation/rectified/rectified_pA_tail_trimmed.bam"
PA_RESTORED_BAM="/Users/kevinroy/work/rectify/rectify/data/validation/rectified/rectified_pA_tail_soft_clipped.bam"
BG_PLUS="/Users/kevinroy/work/rectify/rectify/data/validation/rectified/corrected_3ends.plus.bedgraph"
BG_MINUS="/Users/kevinroy/work/rectify/rectify/data/validation/rectified/corrected_3ends.minus.bedgraph"

for f in "$GENOME" "$GFF" "$MAPPED_BAM" "$CORRECTED_BAM" "$PA_RESTORED_BAM" "$BG_PLUS" "$BG_MINUS"; do
    [[ -f "$f" ]] || { echo "missing: $f" >&2; exit 1; }
done

PORT=60151
HOST="http://localhost:$PORT"
if ! python3 -c "import socket; s=socket.socket(); s.settimeout(0.5); s.connect_ex(('127.0.0.1', $PORT)) == 0 or exit(1)" 2>/dev/null; then
    echo "IGV is not listening on port $PORT." >&2
    echo "Open IGV → Preferences → Advanced → Enable Port (default 60151), then re-run." >&2
    exit 1
fi

url_encode() {
    python3 -c "import urllib.parse, sys; print(urllib.parse.quote(sys.argv[1]))" "$1"
}

# Order matters for IGV's track display:
#   GFF → mapped → corrected → pA-restored → bedgraph+ → bedgraph−
FILES_CSV=""
for f in "$GFF" "$MAPPED_BAM" "$CORRECTED_BAM" "$PA_RESTORED_BAM" "$BG_PLUS" "$BG_MINUS"; do
    enc=$(url_encode "$f")
    if [[ -z "$FILES_CSV" ]]; then
        FILES_CSV="$enc"
    else
        FILES_CSV="${FILES_CSV},${enc}"
    fi
done
GENOME_ENC=$(url_encode "$GENOME")
LOCUS_ENC=$(url_encode "${CHROM}:${START}-${END}")

http_code() {
    curl -s -o /dev/null -w '%{http_code}' "$1"
}

LOAD_CODE=$(http_code "$HOST/load?file=${FILES_CSV}&genome=${GENOME_ENC}")
echo "load: HTTP $LOAD_CODE"
GOTO_CODE=$(http_code "$HOST/goto?locus=${LOCUS_ENC}")
echo "goto: HTTP $GOTO_CODE"

echo
echo "Loaded 3 versions of the read (aligner=${ALIGNER}):"
echo "  1. mapped:      validation_reads.${ALIGNER}.bam"
echo "  2. corrected:   rectified_pA_tail_trimmed.bam"
echo "  3. pA-restored: rectified_pA_tail_soft_clipped.bam"
echo "Plus annotation + plus/minus bedgraphs. Locus: ${CHROM}:${START}-${END}"
echo
echo "If the view is cluttered with prior tracks, right-click any panel → Remove All Tracks, then re-run."
