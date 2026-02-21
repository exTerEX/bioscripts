## Extract organism names from NCBI nucleotide accessions
#
## Usage:
#   $ python misc/extract_organism_from_nucleotides_accession.py -h
#   usage: extract_organism_from_nucleotides_accession.py [-h] --email EMAIL
#                                                         [--api-key API_KEY]
#                                                         [--column COLUMN]
#                                                         [--cache] [--threads THREADS]
#                                                         [--output OUTPUT]
#                                                         input
#
#   Extract organism names from NCBI nucleotide accessions
#
#   positional arguments:
#     input              Input: table file (csv/tsv/xlsx), single accession, or
#                        comma-separated list of accessions
#
#   options:
#     -h, --help         show this help message and exit
#     --email EMAIL      Email address for NCBI Entrez (required)
#     --api-key API_KEY  NCBI API key for higher rate limits
#     --column COLUMN    Column name containing accession numbers (default: Contig)
#     --cache            Cache downloaded FASTA files locally
#     --threads THREADS  Number of threads for parallel processing
#     --output OUTPUT    Output file path. If not given, outputs to stdout. No suffix defaults to CSV.

import argparse
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord


def fetch_fasta(
    accession: str,
    cache: bool = False,
    cache_dir: Path = Path(".cache/fna"),
) -> SeqRecord | None:
    """Fetch a FASTA record from NCBI nucleotide database.

    Args:
        accession: Nucleotide accession number (with or without version).
        cache: Whether to cache the downloaded record locally.
        cache_dir: Directory for cached FASTA files.

    Returns:
        SeqRecord if successful, None if fetch fails.
    """
    # Strip version number if present
    accession_base = accession.split(".")[0]
    cache_path = cache_dir / f"{accession_base}.fna"

    try:
        if cache_path.exists():
            return SeqIO.read(cache_path, "fasta")

        with Entrez.efetch(db="nuccore", id=accession_base, rettype="fasta", retmode="text") as handle:
            record = SeqIO.read(handle, "fasta")

        if cache:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            SeqIO.write(record, cache_path, "fasta")

        return record

    except Exception as e:
        print(f"Error fetching {accession}: {e}")
        return None


def extract_organism(accession: str, cache: bool = False) -> str:
    """Extract organism name from an NCBI nucleotide accession.

    Args:
        accession: NCBI nucleotide accession number.
        cache: Whether to cache downloaded records.

    Returns:
        Organism name or "-" if extraction fails.
    """
    record = fetch_fasta(accession, cache=cache)

    if record is None:
        return "-"

    description = record.description
    parts = description.split()

    # RefSeq accessions (NZ_, NC_) have accession as first word
    if len(parts) >= 3 and parts[0][:2] in ("NZ", "NC"):
        return f"{parts[1]} {parts[2]}"

    # GenBank accessions have organism as first two words
    if len(parts) >= 2:
        return f"{parts[0]} {parts[1]}"

    return "-"


def process_accessions(
    accessions: list[str],
    cache: bool = False,
    threads: int | None = None,
) -> list[str]:
    """Process multiple accessions in parallel using threads.

    Args:
        accessions: List of accession numbers to process.
        cache: Whether to cache downloaded records.
        threads: Number of threads (None for default).

    Returns:
        List of organism names in the same order as input accessions.
    """
    results: dict[int, str] = {}

    with ThreadPoolExecutor(max_workers=threads) as executor:
        # Submit all tasks with their index for ordering
        futures = {executor.submit(extract_organism, acc, cache): i for i, acc in enumerate(accessions)}

        for future in as_completed(futures):
            idx = futures[future]
            try:
                results[idx] = future.result()
            except Exception as e:
                print(f"Error processing accession at index {idx}: {e}")
                results[idx] = "-"

    # Return results in original order
    return [results[i] for i in range(len(accessions))]


def read_table_file(path: Path) -> pd.DataFrame:
    """Read a table file based on its extension.

    Args:
        path: Path to the table file.

    Returns:
        DataFrame with the table contents.

    Raises:
        ValueError: If file extension is not supported.
    """
    suffix = path.suffix.lower()

    match suffix:
        case ".csv":
            return pd.read_csv(path)
        case ".tsv" | ".txt":
            return pd.read_csv(path, sep="\t")
        case ".xlsx" | ".xls":
            return pd.read_excel(path)
        case _:
            raise ValueError(f"Unsupported file format: {suffix}")


def write_table_file(df: pd.DataFrame, path: Path):
    """Write a DataFrame to a table file based on its extension.

    Args:
        df: DataFrame to write.
        path: Output file path. No suffix defaults to CSV.

    Raises:
        ValueError: If file extension is not supported.
    """
    suffix = path.suffix.lower()

    match suffix:
        case ".csv" | "":  # No suffix defaults to CSV
            df.to_csv(path, index=False)
        case ".tsv" | ".txt":
            df.to_csv(path, sep="\t", index=False)
        case ".xlsx" | ".xls":
            df.to_excel(path, index=False)
        case _:
            raise ValueError(f"Unsupported file format: {suffix}")


def write_to_stdout(df: pd.DataFrame):
    """Write a DataFrame to stdout as tab-separated values.

    Args:
        df: DataFrame to write.
    """
    df.to_csv(sys.stdout, sep="\t", index=False)


def is_table_file(input_str: str) -> bool:
    """Check if input string looks like a table file path.

    Args:
        input_str: Input string to check.

    Returns:
        True if it looks like a table file, False otherwise.
    """
    path = Path(input_str)
    return path.exists() and path.suffix.lower() in (
        ".csv",
        ".tsv",
        ".txt",
        ".xlsx",
        ".xls",
    )


def process_table(
    input_path: Path,
    output_path: Path | None,
    column: str,
    cache: bool,
    threads: int | None,
):
    """Process a table file and add organism names.

    Args:
        input_path: Path to input table file.
        output_path: Path for output table file, or None for stdout.
        column: Column name containing accession numbers.
        cache: Whether to cache downloaded records.
        threads: Number of threads for parallel processing.
    """
    print(f"Reading input file: {input_path}", file=sys.stderr)
    df = read_table_file(input_path)

    if column not in df.columns:
        raise ValueError(f"Column '{column}' not found in input file")

    # Filter rows with valid accessions
    valid_mask = df[column].notna()
    accessions = df.loc[valid_mask, column].astype(str).tolist()

    print(f"Processing {len(accessions)} accessions...", file=sys.stderr)

    # Process accessions
    organisms = process_accessions(accessions, cache=cache, threads=threads)

    # Add results to dataframe
    df.loc[valid_mask, "Organism"] = organisms
    df["Organism"] = df["Organism"].fillna("-")

    # Write output
    if output_path:
        write_table_file(df, output_path)
        print(f"Output written to: {output_path}", file=sys.stderr)
    else:
        write_to_stdout(df)

    print(f"Total rows: {len(df)}", file=sys.stderr)


def process_accession_list(
    accessions: list[str],
    output_path: Path | None,
    cache: bool,
    threads: int | None,
):
    """Process a list of accessions and output results.

    Args:
        accessions: List of accession numbers.
        output_path: Optional output file path, or None for stdout.
        cache: Whether to cache downloaded records.
        threads: Number of threads for parallel processing.
    """
    print(f"Processing {len(accessions)} accessions...", file=sys.stderr)

    organisms = process_accessions(accessions, cache=cache, threads=threads)

    # Create result dataframe
    df = pd.DataFrame({"Accession": accessions, "Organism": organisms})

    if output_path:
        write_table_file(df, output_path)
        print(f"Output written to: {output_path}", file=sys.stderr)
    else:
        write_to_stdout(df)


def main(
    input_str: str,
    output_path: Path | None = None,
    email: str = "",
    api_key: str | None = None,
    column: str = "Contig",
    cache: bool = False,
    threads: int | None = None,
):
    """Main function to process accessions and extract organism names.

    Args:
        input_str: Input table file path, single accession, or comma-separated list.
        output_path: Path for output file.
        email: Email for NCBI Entrez.
        api_key: Optional NCBI API key for higher rate limits.
        column: Column name containing accession numbers (for table input).
        cache: Whether to cache downloaded records.
        threads: Number of threads for parallel processing.
    """
    # Configure Entrez
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    # Determine input type
    if is_table_file(input_str):
        process_table(Path(input_str), output_path, column, cache, threads)
    else:
        # Parse as accession(s)
        accessions = [acc.strip() for acc in input_str.split(",") if acc.strip()]
        if not accessions:
            raise ValueError("No valid accessions provided")
        process_accession_list(accessions, output_path, cache, threads)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract organism names from NCBI nucleotide accessions")

    parser.add_argument(
        "input",
        type=str,
        help="Input: table file (csv/tsv/xlsx), single accession, or comma-separated list",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output file path. If not given, outputs to stdout. No suffix defaults to CSV.",
    )
    parser.add_argument(
        "--email",
        type=str,
        required=True,
        help="Email address for NCBI Entrez (required)",
    )
    parser.add_argument(
        "--api-key",
        type=str,
        default=None,
        help="NCBI API key for higher rate limits",
    )
    parser.add_argument(
        "--column",
        type=str,
        default="Contig",
        help="Column name containing accession numbers (default: Contig)",
    )
    parser.add_argument(
        "--cache",
        action="store_true",
        help="Cache downloaded FASTA files locally",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Number of threads for parallel processing",
    )

    args = parser.parse_args()

    main(
        args.input,
        args.output,
        args.email,
        args.api_key,
        args.column,
        args.cache,
        args.threads,
    )
