#!/usr/bin/env bash
# Drive IGV to a specific validation read, bring IGV to the foreground, save a
# screenshot of the IGV window. Writes the PNG to /tmp/igv_snapshots/<xv>.png.
#
# Usage:
#   bash scripts/validation_data/review_read.sh <xv_label> [aligner] [pad_bp]
#
# Examples:
#   bash scripts/validation_data/review_read.sh cat1_plus_1
#   bash scripts/validation_data/review_read.sh cat1_minus_1 minimap2
#   bash scripts/validation_data/review_read.sh cat5_plus_1 mapPacBio 100
#
# Reads coordinates from rectify/data/validation/corrected_reads.tsv (keyed by
# the read_id that the XV tag points at). The default pad is 50 bp on each side
# of the corrected_3prime position.

set -euo pipefail

XV=${1:?xv label required, e.g. cat1_plus_1}
ALIGNER=${2:-minimap2}
PAD=${3:-50}

REPO="/Users/kevinroy/work/rectify"
VAL="$REPO/rectify/data/validation"
TSV="$VAL/corrected_reads.tsv"
BAM_XV_MAP="$VAL/validation_reads.bam"
GENOME="$REPO/rectify/data/genomes/saccharomyces_cerevisiae/S288C_reference_sequence_R64-5-1_20240529.fsa.gz"
GFF="$REPO/rectify/data/genomes/saccharomyces_cerevisiae/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz"
MAPPED="$VAL/aligners/validation_reads.${ALIGNER}.bam"
CORRECTED="$VAL/rectified/rectified_pA_tail_trimmed.bam"
PA_REST="$VAL/rectified/rectified_pA_tail_soft_clipped.bam"
BG_PLUS="$VAL/rectified/corrected_3ends.plus.bedgraph"
BG_MINUS="$VAL/rectified/corrected_3ends.minus.bedgraph"
OUT="/tmp/igv_snapshots/${XV}.png"

mkdir -p /tmp/igv_snapshots

for f in "$GENOME" "$GFF" "$MAPPED" "$CORRECTED" "$PA_REST" "$BG_PLUS" "$BG_MINUS" "$TSV" "$BAM_XV_MAP"; do
    [[ -f "$f" ]] || { echo "missing: $f" >&2; exit 1; }
done

# Resolve XV label → full UUID via the validation BAM
UUID=$(python3 - <<PY
import pysam, sys
xv = "$XV"
with pysam.AlignmentFile("$BAM_XV_MAP", "rb") as bam:
    for r in bam:
        if r.is_secondary or r.is_supplementary or r.is_unmapped:
            continue
        try:
            if r.get_tag("XV") == xv:
                print(r.query_name); sys.exit(0)
        except KeyError:
            pass
sys.exit(1)
PY
) || { echo "XV label $XV not found in $BAM_XV_MAP" >&2; exit 1; }

# Pull corrected info from corrected_reads.tsv
INFO=$(python3 - <<PY
import csv
uuid = "$UUID"
with open("$TSV") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        if row["read_id"] == uuid:
            print(row["chrom"], row["strand"], row["alignment_start"],
                  row["alignment_end"], row["original_3prime"],
                  row["corrected_3prime"], row["correction_applied"],
                  row.get("five_prime_position", ""), sep="\t")
            break
PY
)
[[ -n "$INFO" ]] || { echo "UUID $UUID not in $TSV" >&2; exit 1; }
read -r CHROM STRAND ALN_START ALN_END ORIG_3P CORR_3P APPLIED FIVE_PRIME <<<"$INFO"

# Window centered on the corrected 3' end
LEFT=$(( CORR_3P - PAD ))
RIGHT=$(( CORR_3P + PAD ))
[[ $LEFT -lt 1 ]] && LEFT=1

PORT=60151
HOST="http://localhost:$PORT"
if ! python3 -c "import socket; s=socket.socket(); s.settimeout(0.5); s.connect_ex(('127.0.0.1',$PORT))==0 or exit(1)" 2>/dev/null; then
    echo "IGV is not listening on port $PORT. Launch IGV and retry." >&2
    exit 1
fi

url_encode() { python3 -c "import urllib.parse, sys; print(urllib.parse.quote(sys.argv[1]))" "$1"; }

FILES_CSV=""
for f in "$GFF" "$MAPPED" "$CORRECTED" "$PA_REST" "$BG_PLUS" "$BG_MINUS"; do
    enc=$(url_encode "$f")
    [[ -z "$FILES_CSV" ]] && FILES_CSV="$enc" || FILES_CSV="${FILES_CSV},${enc}"
done

curl -s -o /dev/null -w 'load: HTTP %{http_code}\n' \
  "$HOST/load?file=${FILES_CSV}&genome=$(url_encode "$GENOME")"
curl -s -o /dev/null -w 'goto: HTTP %{http_code}\n' \
  "$HOST/goto?locus=$(url_encode "${CHROM}:${LEFT}-${RIGHT}")"

# Front IGV (process name is "java" since IGV is a JVM app)
osascript -e 'tell application "System Events" to set frontmost of (first process whose name is "java") to true' >/dev/null 2>&1 || true
sleep 1.5
screencapture -x "$OUT"

echo
echo "  XV:                ${XV}"
echo "  read_id:           ${UUID}"
echo "  category:          $(python3 -c "
import pysam
with pysam.AlignmentFile('$BAM_XV_MAP','rb') as b:
    for r in b:
        if r.is_secondary or r.is_supplementary or r.is_unmapped: continue
        try:
            if r.get_tag('XV')=='$XV': print(r.get_tag('XG')); break
        except KeyError: pass
")"
echo "  locus:             ${CHROM}:${LEFT}-${RIGHT}  (corr_3p±${PAD})"
echo "  strand:            ${STRAND}"
echo "  alignment span:    [${ALN_START}, ${ALN_END})"
echo "  five_prime_pos:    ${FIVE_PRIME}"
echo "  original_3prime:   ${ORIG_3P}"
echo "  corrected_3prime:  ${CORR_3P}"
echo "  correction_applied: ${APPLIED}"
echo "  screenshot:        ${OUT}"
