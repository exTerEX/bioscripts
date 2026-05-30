#!/bin/bash
# Searches protein query sequences against a translated nucleotide database.
# Usage: sbatch tblastn_fox.slurm
# Edit the variables below before submitting.

#SBATCH --account= ### !!! Add project id
#SBATCH --job-name=tblastn
#SBATCH --output=logs/tblastn_%j.out
#SBATCH --error=logs/tblastn_%j.err
#SBATCH --time=08:00:00 ### !!! Add run time
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G

# ── User-defined variables ────────────────────────────────────────────────────
QUERY=""                           # Path to query protein FASTA file
DB=""                              # Path to BLAST nucleotide database (without extension)
OUT=""                             # Output file path
EVALUE="1e-5"
OUTFMT="16"                        # XML2 output
# OUTFMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
MAX_TARGET_SEQS=500
EXTRA_ARGS=""                      # Any additional tblastn arguments
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

if [[ -z "$QUERY" || -z "$DB" || -z "$OUT" ]]; then
    echo "ERROR: QUERY, DB, and OUT must be set before submitting." >&2
    exit 1
fi

mkdir -p logs "$(dirname "$OUT")"

# Load BLAST+; check available versions with: module avail BLAST+
module purge
module load BLAST+

echo "Host:     $(hostname)"
echo "Started:  $(date)"
echo "Query:    $QUERY"
echo "Database: $DB"
echo "Output:   $OUT"

tblastn \
    -query        "$QUERY" \
    -db           "$DB" \
    -out          "$OUT" \
    -evalue       "$EVALUE" \
    -outfmt       "$OUTFMT" \
    -max_target_seqs "$MAX_TARGET_SEQS" \
    -num_threads  "$SLURM_CPUS_PER_TASK" \
    $EXTRA_ARGS

echo "Finished: $(date)"
