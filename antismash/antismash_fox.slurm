#!/bin/bash
# antiSMASH SLURM array job script for Educloud Fox (UiO)
# Runs one .gbff file per array task.
#
# Usage:
#   1) Set INPUT_DIR/OUTPUT_ROOT and SBATCH values below.
#   2) Count files and submit with matching array range, for example:
#      N=$(find ./genbank -maxdepth 1 -type f -name '*.gbff' | wc -l)
#      sbatch --array=0-$((N-1)) antismash/antismash_fox_array.slurm

#SBATCH --account=ec000            # Replace with your Educloud project (ecXXX)
#SBATCH --job-name=antismash_arr
#SBATCH --output=logs/antismash_%A_%a.out
#SBATCH --error=logs/antismash_%A_%a.err
#SBATCH --time=24:00:00            # Adjust walltime (hh:mm:ss or d-hh:mm:ss)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -euo pipefail

# ── User-defined variables ────────────────────────────────────────────────────
INPUT_DIR="./genbank"              # Directory containing .gbff input files
OUTPUT_ROOT="./antismash"          # Root output directory for antiSMASH results
EXTRA_ARGS=""                      # Optional extra antiSMASH flags
# ─────────────────────────────────────────────────────────────────────────────

mkdir -p logs "$OUTPUT_ROOT"

# Load antiSMASH; check available versions with: module avail antismash
module purge
module load antiSMASH

echo "Host:      $(hostname)"
echo "Started:   $(date)"
echo "Array job: ${SLURM_ARRAY_JOB_ID:-N/A}"
echo "Task ID:   ${SLURM_ARRAY_TASK_ID:-N/A}"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: This script must be submitted as an array job (use sbatch --array=...)." >&2
    exit 1
fi

shopt -s nullglob
files=("$INPUT_DIR"/*.gbff)
if (( ${#files[@]} == 0 )); then
    echo "ERROR: No .gbff files found in $INPUT_DIR" >&2
    exit 1
fi

idx=$SLURM_ARRAY_TASK_ID
if (( idx < 0 || idx >= ${#files[@]} )); then
    echo "ERROR: SLURM_ARRAY_TASK_ID=$idx is out of range for ${#files[@]} file(s)." >&2
    echo "       Submit with --array=0-$(( ${#files[@]} - 1 ))" >&2
    exit 1
fi

filename="${files[$idx]}"
basepath="${filename%.*}"
basename="${basepath##*/}"

echo "Input:     $filename"
echo "Basename:  $basename"
echo "Output:    $OUTPUT_ROOT/$basename"

antismash \
    -t bacteria \
    --output-dir "$OUTPUT_ROOT/$basename" \
    --output-basename "$basename" \
    --fullhmmer \
    --clusterhmmer \
    --tigrfam \
    --asf \
    --cc-mibig \
    --cb-general \
    --cb-subclusters \
    --cb-knownclusters \
    --pfam2go \
    --rre \
    --smcog-trees \
    --tfbs \
    --genefinding-tool prodigal \
    $EXTRA_ARGS \
    "$basepath.gbff"

echo "Finished:  $(date)"
