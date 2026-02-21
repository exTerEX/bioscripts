#!/usr/bin/env bash

set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

readonly REQUIRED_COMMANDS=("datasets" "dataformat" "unzip" "gzip")

# Define a function to display the help message
show_help() {
  echo "Usage: $0 [OPTIONS]"
  echo "Download genomes by taxon ID using NCBI datasets."
  echo ""
  echo "Options:"
  echo "  -h, --help          Display this help message"
  echo ""
  echo "  -t <taxon>          Taxon ID or comma-separated list of taxon IDs"
  echo ""
  echo "  -o <output>         Output directory for downloaded genomic data"
  echo ""
  echo "  -a <NCBI_API_KEY>   NCBI API key for higher rate limits (optional)."
  echo "                      Falls back to \$NCBI_API_KEY environment variable."
  exit 0
}

# Check that all required commands are available
check_dependencies() {
    local missing=()
    for cmd in "${REQUIRED_COMMANDS[@]}"; do
        if ! command -v "$cmd" &> /dev/null; then
            missing+=("$cmd")
        fi
    done

    if [[ ${#missing[@]} -gt 0 ]]; then
        echo "Error: the following required commands are not available:" >&2
        for cmd in "${missing[@]}"; do
            echo "  - $cmd" >&2
        done
        echo "" >&2
        echo "Install NCBI datasets: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/" >&2
        exit 1
    fi
}

# Check if no arguments are provided
if [[ $# -eq 0 ]]; then
    show_help
fi

# Process arguments
while getopts "ht:o:a:" opt; do
    case $opt in
        h)
            show_help
            ;;
        t)
            input_taxon="$OPTARG"
            ;;
        o)
            output_directory="$OPTARG"
            ;;
        a)
            NCBI_API_KEY="$OPTARG"
            ;;
        \?)
            echo "Error: Invalid option -$OPTARG" >&2
            show_help
            exit 1
            ;;
    esac
done

# Validate required arguments
if [[ ! -v input_taxon || ! -v output_directory ]]; then
    echo "Error: Both -t <taxon> and -o <output> are required." >&2
    exit 1
fi

check_dependencies

# Paths
filepath_ext="${output_directory}.zip"

# Create output directory if needed
mkdir -p "$output_directory"

# Download genomes (skip if zip already exists)
if [[ ! -e "$filepath_ext" ]]; then
    api_key_args=()
    if [[ -v NCBI_API_KEY ]]; then
        api_key_args=(--api-key "$NCBI_API_KEY")
    fi

    datasets download genome taxon "$input_taxon" --assembly-source RefSeq \
        "${api_key_args[@]}" \
        --dehydrated --include genome,protein,cds,gff3,gtf,gbff,seq-report \
        --filename "$filepath_ext"
fi

# Unzip for rehydration
unzip -o "$filepath_ext" -d "$output_directory"

# Rehydrate the genome data
datasets rehydrate --directory "$output_directory"

# Move files to final location
mv "$output_directory"/ncbi_dataset/data/* "$output_directory"

# Generate metadata files from assembly data report
report_file="$output_directory/assembly_data_report.jsonl"
if [[ -f "$report_file" ]]; then
    dataformat tsv genome --inputfile "$report_file" > "$output_directory/metadata.tsv"
    dataformat excel genome --inputfile "$report_file" --outputfile "$output_directory/metadata.xlsx"
fi

# Remove temporary files
rm -rf "$output_directory/ncbi_dataset" "$output_directory/README.md" "$output_directory/md5sum.txt" "$output_directory/assembly_data_report.jsonl"

# Rename the files within each genome directory to use the accession as prefix
for directory in "$output_directory"/*/; do
    [[ -d "$directory" ]] || continue

    accession_id=$(basename "$directory")

    [[ -f "$directory/cds_from_genomic.fna" ]] && mv "$directory/cds_from_genomic.fna" "$directory/${accession_id}_cds.fna"
    for fna in "$directory/${accession_id}"_*_genomic.fna; do
        [[ -f "$fna" ]] && mv "$fna" "$directory/${accession_id}_genomic.fna"
        break
    done
    [[ -f "$directory/genomic.gbff" ]]           && mv "$directory/genomic.gbff"           "$directory/${accession_id}.gbff"
    [[ -f "$directory/genomic.gff" ]]            && mv "$directory/genomic.gff"            "$directory/${accession_id}.gff"
    [[ -f "$directory/genomic.gtf" ]]            && mv "$directory/genomic.gtf"            "$directory/${accession_id}.gtf"
    [[ -f "$directory/protein.faa" ]]            && mv "$directory/protein.faa"            "$directory/${accession_id}.faa"
    [[ -f "$directory/sequence_report.jsonl" ]]  && mv "$directory/sequence_report.jsonl"  "$directory/${accession_id}.jsonl"

    # Compress genomic files
    for ext in gbff gff gtf faa fna; do
        for f in "$directory"/*."$ext"; do
            [[ -f "$f" ]] && gzip "$f"
        done
    done
done
