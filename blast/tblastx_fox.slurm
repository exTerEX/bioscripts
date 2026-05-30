#!/bin/bash
# Searches translated nucleotide queries against a translated nucleotide database.
# NOTE: tblastx does NOT support -num_threads; each job runs single-threaded.
#       For many queries, consider splitting with --array (see comments below).
# Usage: sbatch tblastx_fox.slurm
# Edit the variables below before submitting.

#SBATCH --account= ### !!! Add project id
#SBATCH --job-name=tblastx
#SBATCH --output=logs/tblastx_%j.out
#SBATCH --error=logs/tblastx_%j.err
#SBATCH --time=08:00:00 ### !!! Add run time
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G

# ── User-defined variables ────────────────────────────────────────────────────
QUERY=""                           # Path to query nucleotide FASTA file
DB=""                              # Path to BLAST nucleotide database (without extension)
OUT=""                             # Output file path
EVALUE="1e-5"
OUTFMT="16"                        # XML2 output
# OUTFMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
MAX_TARGET_SEQS=500
EXTRA_ARGS=""                      # Any additional tblastx arguments
# ─────────────────────────────────────────────────────────────────────────────

# ── Array job usage (multiple queries) ───────────────────────────────────────
# If you have many query files in a directory, submit with:
#   sbatch --array=0-N tblastx_fox.slurm
# Then set QUERY_DIR below and uncomment the array block.
#
# QUERY_DIR=""
# QUERIES=("$QUERY_DIR"/*.fasta)
# QUERY="${QUERIES[$SLURM_ARRAY_TASK_ID]}"
# OUT="${QUERY%.fasta}.tblastx.tsv"
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

tblastx \
    -query        "$QUERY" \
    -db           "$DB" \
    -out          "$OUT" \
    -evalue       "$EVALUE" \
    -outfmt       "$OUTFMT" \
    -max_target_seqs "$MAX_TARGET_SEQS" \
    $EXTRA_ARGS

echo "Finished: $(date)"
