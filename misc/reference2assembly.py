## Convert NCBI RefSeq/nucleotide accessions to assembly accessions
#
## Usage:
#   $ python misc/reference2assembly.py -h
#   usage: reference2assembly.py [-h] --email EMAIL [--api-key API_KEY]
#                                [--column COLUMN] [--output OUTPUT]
#                                [--threads THREADS] [--cache]
#                                input
#
#   Convert NCBI nucleotide accessions to assembly accessions
#
#   positional arguments:
#     input              Input: table file (csv/tsv/xlsx), single accession,
#                        or comma-separated list of accessions
#
#   options:
#     -h, --help         show this help message and exit
#     --email EMAIL      Email address for NCBI Entrez (required)
#     --api-key API_KEY  NCBI API key for higher rate limits
#     --column COLUMN    Column name containing accession numbers (default: Nucleotide_accession)
#     --output OUTPUT    Output file path. If not given, outputs to stdout
#     --threads THREADS  Number of threads for parallel processing
#     --cache            Cache intermediate results locally

import argparse
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
from Bio import Entrez


def refseq_to_assembly(accession: str, cache_dir: Path | None = None) -> str | None:
    """Convert an NCBI nucleotide accession to its assembly accession.

    Args:
        accession: NCBI nucleotide accession (e.g., 'NZ_AICP01000040.1').
        cache_dir: Optional directory for caching results.

    Returns:
        Assembly accession (e.g., 'GCF_023167545.1') or None if not found.
    """
    # Check cache first
    if cache_dir:
        cache_file = cache_dir / f"{accession.replace('.', '_')}.txt"
        if cache_file.exists():
            result = cache_file.read_text().strip()
            return result if result else None

    try:
        # Step 1: Search nuccore for the accession
        with Entrez.esearch(db="nuccore", term=accession, retmax=1) as handle:
            search_result = Entrez.read(handle)

        if not search_result["IdList"]:
            print(f"Warning: No nuccore ID found for: {accession}", file=sys.stderr)
            return None

        nuccore_id = search_result["IdList"][0]

        # Step 2: Link nuccore to assembly database
        with Entrez.elink(dbfrom="nuccore", db="assembly", id=nuccore_id) as handle:
            link_result = Entrez.read(handle)

        assembly_id = None
        if link_result and link_result[0].get("LinkSetDb"):
            for link_set in link_result[0]["LinkSetDb"]:
                if link_set["DbTo"] == "assembly" and link_set.get("Link"):
                    assembly_id = link_set["Link"][0]["Id"]
                    break

        if not assembly_id:
            print(f"Warning: No assembly ID found for: {accession}", file=sys.stderr)
            if cache_dir:
                cache_dir.mkdir(parents=True, exist_ok=True)
                cache_file.write_text("")
            return None

        # Step 3: Fetch assembly summary to get accession
        with Entrez.esummary(db="assembly", id=assembly_id) as handle:
            summary = Entrez.read(handle)

        doc_summary = summary.get("DocumentSummarySet", {}).get("DocumentSummary", [])
        if doc_summary:
            assembly_acc = doc_summary[0].get("AssemblyAccession", doc_summary[0].get("Accession"))
            if assembly_acc and cache_dir:
                cache_dir.mkdir(parents=True, exist_ok=True)
                cache_file.write_text(assembly_acc)
            return assembly_acc

    except Exception as e:
        print(f"Error processing {accession}: {e}", file=sys.stderr)

    return None


def process_accessions(
    accessions: list[str],
    cache: bool = False,
    threads: int | None = None,
) -> list[str | None]:
    """Process multiple accessions in parallel.

    Args:
        accessions: List of nucleotide accessions.
        cache: Whether to cache results.
        threads: Number of threads (None for default).

    Returns:
        List of assembly accessions in the same order as input.
    """
    cache_dir = Path(".cache/assembly") if cache else None
    results: dict[int, str | None] = {}

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(refseq_to_assembly, acc, cache_dir): i for i, acc in enumerate(accessions)}

        for future in as_completed(futures):
            idx = futures[future]
            try:
                results[idx] = future.result()
            except Exception as e:
                print(f"Error at index {idx}: {e}", file=sys.stderr)
                results[idx] = None

    return [results[i] for i in range(len(accessions))]


def read_table_file(path: Path) -> pd.DataFrame:
    """Read a table file based on its extension."""
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
    """Write a DataFrame to a table file based on its extension."""
    suffix = path.suffix.lower()
    match suffix:
        case ".csv" | "":
            df.to_csv(path, index=False)
        case ".tsv" | ".txt":
            df.to_csv(path, sep="\t", index=False)
        case ".xlsx" | ".xls":
            df.to_excel(path, index=False)
        case _:
            raise ValueError(f"Unsupported file format: {suffix}")


def write_to_stdout(df: pd.DataFrame):
    """Write a DataFrame to stdout as tab-separated values."""
    df.to_csv(sys.stdout, sep="\t", index=False)


def is_table_file(input_str: str) -> bool:
    """Check if input string looks like a table file path."""
    path = Path(input_str)
    return path.exists() and path.suffix.lower() in (".csv", ".tsv", ".txt", ".xlsx", ".xls")


def main(
    input_str: str,
    output_path: Path | None = None,
    email: str = "",
    api_key: str | None = None,
    column: str = "Nucleotide_accession",
    cache: bool = False,
    threads: int | None = None,
):
    """Main function to convert nucleotide accessions to assembly accessions.

    Args:
        input_str: Input table file path, single accession, or comma-separated list.
        output_path: Path for output file, or None for stdout.
        email: Email for NCBI Entrez.
        api_key: Optional NCBI API key for higher rate limits.
        column: Column name containing accession numbers (for table input).
        cache: Whether to cache results.
        threads: Number of threads for parallel processing.
    """
    # Configure Entrez
    Entrez.email = email
    Entrez.tool = "reference2assembly"
    if api_key:
        Entrez.api_key = api_key

    if is_table_file(input_str):
        # Process table file
        print(f"Reading input file: {input_str}", file=sys.stderr)
        df = read_table_file(Path(input_str))

        if column not in df.columns:
            raise ValueError(f"Column '{column}' not found in input file")

        valid_mask = df[column].notna()
        accessions = df.loc[valid_mask, column].astype(str).tolist()

        print(f"Processing {len(accessions)} accessions...", file=sys.stderr)
        assemblies = process_accessions(accessions, cache=cache, threads=threads)

        df.loc[valid_mask, "Assembly_accession"] = assemblies
        df["Assembly_accession"] = df["Assembly_accession"].fillna("-")

        if output_path:
            write_table_file(df, output_path)
            print(f"Output written to: {output_path}", file=sys.stderr)
        else:
            write_to_stdout(df)

        print(f"Total rows: {len(df)}", file=sys.stderr)

    else:
        # Process accession(s)
        accessions = [acc.strip() for acc in input_str.split(",") if acc.strip()]
        if not accessions:
            raise ValueError("No valid accessions provided")

        print(f"Processing {len(accessions)} accessions...", file=sys.stderr)
        assemblies = process_accessions(accessions, cache=cache, threads=threads)

        df = pd.DataFrame(
            {
                "Nucleotide_accession": accessions,
                "Assembly_accession": [a or "-" for a in assemblies],
            }
        )

        if output_path:
            write_table_file(df, output_path)
            print(f"Output written to: {output_path}", file=sys.stderr)
        else:
            write_to_stdout(df)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert NCBI nucleotide accessions to assembly accessions")

    parser.add_argument(
        "input",
        type=str,
        help="Input: table file (csv/tsv/xlsx), single accession, or comma-separated list",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output file path. If not given, outputs to stdout",
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
        default="Nucleotide_accession",
        help="Column name containing accession numbers (default: Nucleotide_accession)",
    )
    parser.add_argument(
        "--cache",
        action="store_true",
        help="Cache intermediate results locally",
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
