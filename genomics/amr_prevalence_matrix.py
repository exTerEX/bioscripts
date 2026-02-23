## Build AMR gene presence/absence matrix from GBFF files and AMR tool output
#
## Usage:
#   $ python genomics/amr_prevalence_matrix.py -h
#   usage: amr_prevalence_matrix.py [-h] [--format {amrfinder,rgi,abricate,tsv}]
#                                   [--gene-col GENE_COL] [--sample-col SAMPLE_COL]
#                                   [--threads THREADS] [--output OUTPUT]
#                                   amr_results gbff_dir
#
#   Build a presence/absence or count matrix of AMR genes per strain. Accepts
#   AMRFinderPlus, RGI, ABRicate, or generic TSV output and optionally enriches
#   with metadata from GBFF files (species, assembly level, isolation source).
#
#   positional arguments:
#     amr_results           AMR tool output file (TSV)
#     gbff_dir              Directory containing GBFF/GBK files for metadata
#
#   options:
#     -h, --help            show this help message and exit
#     --format              Input format: amrfinder, rgi, abricate, or tsv
#     --gene-col GENE_COL   Column name for gene identifier (for --format tsv)
#     --sample-col SAMPLE_COL Column name for sample identifier (for --format tsv)
#     --threads THREADS     Number of threads for parallel processing
#     --output OUTPUT       Output TSV file path (default: stdout)

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# ---------------------------------------------------------------------------
# Metadata extraction from GBFF
# ---------------------------------------------------------------------------


def extract_metadata(path: Path) -> dict[str, str]:
    """Extract strain-level metadata from a GBFF file.

    Args:
        path: Path to a GenBank file.

    Returns:
        Dictionary with organism, strain, isolation_source, host, and
        assembly_level fields.
    """
    meta: dict[str, str] = {
        "file": path.stem,
        "organism": "",
        "strain": "",
        "isolation_source": "",
        "host": "",
        "assembly_level": "",
    }

    try:
        record: SeqRecord = next(SeqIO.parse(path, "genbank"))
    except (StopIteration, Exception):
        return meta

    # Organism from annotations
    organism = record.annotations.get("organism", "")
    meta["organism"] = organism if isinstance(organism, str) else ""

    # Source feature qualifiers
    for feat in record.features:
        if feat.type != "source":
            continue
        quals = feat.qualifiers
        meta["strain"] = quals.get("strain", [""])[0]
        meta["isolation_source"] = quals.get("isolation_source", [""])[0]
        meta["host"] = quals.get("host", [""])[0]

        # Assembly level from structured comments
        for comment in record.annotations.get("structured_comment", {}).values():  # type: ignore
            if isinstance(comment, dict):
                meta["assembly_level"] = comment.get("Assembly Level", "")
        break

    return meta


# ---------------------------------------------------------------------------
# AMR result parsers
# ---------------------------------------------------------------------------


def _normalise_sample_id(raw: str, gbff_stems: set[str]) -> str:
    """Try to match a raw sample identifier to a GBFF file stem.

    Falls back to the raw string if no match is found.
    """
    # Exact match
    if raw in gbff_stems:
        return raw
    # Strip common extensions
    for ext in (".gbff", ".gbk", ".gb", ".fna", ".fasta", ".fa"):
        candidate = raw.removesuffix(ext)
        if candidate in gbff_stems:
            return candidate
    # Basename only
    candidate = Path(raw).stem
    if candidate in gbff_stems:
        return candidate
    return raw


def parse_amrfinder(path: Path) -> list[tuple[str, str]]:
    """Parse AMRFinderPlus output into (sample, gene) tuples."""
    pairs: list[tuple[str, str]] = []
    with path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample = row.get("Name", row.get("Contig id", ""))
            gene = row.get("Gene symbol", row.get("Sequence name", ""))
            if sample and gene:
                pairs.append((sample, gene))
    return pairs


def parse_rgi(path: Path) -> list[tuple[str, str]]:
    """Parse RGI main output into (sample, gene) tuples."""
    pairs: list[tuple[str, str]] = []
    with path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            # RGI ORF_ID contains the contig/sample info
            orf_id = row.get("ORF_ID", "")
            sample = orf_id.rsplit("_", 1)[0] if "_" in orf_id else orf_id
            gene = row.get("Best_Hit_ARO", row.get("AMR Gene Family", ""))
            if sample and gene:
                pairs.append((sample, gene))
    return pairs


def parse_abricate(path: Path) -> list[tuple[str, str]]:
    """Parse ABRicate output into (sample, gene) tuples."""
    pairs: list[tuple[str, str]] = []
    with path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample = row.get("#FILE", row.get("FILE", ""))
            gene = row.get("GENE", row.get("RESISTANCE", ""))
            if sample and gene:
                pairs.append((sample, gene))
    return pairs


def parse_generic_tsv(path: Path, sample_col: str, gene_col: str) -> list[tuple[str, str]]:
    """Parse a generic TSV with user-specified column names."""
    pairs: list[tuple[str, str]] = []
    with path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample = row.get(sample_col, "")
            gene = row.get(gene_col, "")
            if sample and gene:
                pairs.append((sample, gene))
    return pairs


# ---------------------------------------------------------------------------
# Matrix construction
# ---------------------------------------------------------------------------


def build_matrix(pairs: list[tuple[str, str]], gbff_stems: set[str]) -> pd.DataFrame:
    """Build a presence/absence matrix from (sample, gene) pairs.

    Args:
        pairs: List of (sample_identifier, gene_name) tuples.
        gbff_stems: Set of GBFF filename stems for sample normalisation.

    Returns:
        DataFrame with samples as rows and genes as columns (1/0 values).
    """
    counts: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))

    for raw_sample, gene in pairs:
        sample = _normalise_sample_id(raw_sample, gbff_stems)
        counts[sample][gene] += 1

    df = pd.DataFrame.from_dict(counts, orient="index").fillna(0).astype(int)
    df.index.name = "sample"
    df = df.reindex(sorted(df.columns), axis=1)
    df.sort_index(inplace=True)

    # Add total column
    df.insert(0, "total_genes", df.sum(axis=1))

    return df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main(
    amr_results: Path,
    gbff_dir: Path,
    fmt: str = "amrfinder",
    gene_col: str = "gene",
    sample_col: str = "sample",
    threads: int | None = None,
    output: Path | None = None,
) -> None:
    """Build AMR prevalence matrix with optional strain metadata.

    Args:
        amr_results: Path to AMR tool output TSV.
        gbff_dir: Directory containing GBFF files for metadata.
        fmt: Input format name.
        gene_col: Column name for gene (generic TSV mode).
        sample_col: Column name for sample (generic TSV mode).
        threads: Number of threads.
        output: Output TSV path, or None for stdout.
    """
    # 1. Parse AMR results
    match fmt:
        case "amrfinder":
            pairs = parse_amrfinder(amr_results)
        case "rgi":
            pairs = parse_rgi(amr_results)
        case "abricate":
            pairs = parse_abricate(amr_results)
        case "tsv":
            pairs = parse_generic_tsv(amr_results, sample_col, gene_col)
        case _:
            print(f"Error: Unknown format '{fmt}'", file=sys.stderr)
            sys.exit(1)

    if not pairs:
        print("Error: No AMR gene entries parsed from input.", file=sys.stderr)
        sys.exit(1)

    print(f"Parsed {len(pairs)} gene detections from {amr_results.name}", file=sys.stderr)

    # 2. Gather GBFF files for metadata
    gbff_files = list(gbff_dir.glob("*.gbff")) + list(gbff_dir.glob("*.gbk")) + list(gbff_dir.glob("*.gb"))
    gbff_stems = {p.stem for p in gbff_files}

    print(f"Found {len(gbff_files)} GenBank file(s) for metadata", file=sys.stderr)

    # 3. Build matrix
    matrix = build_matrix(pairs, gbff_stems)

    # 4. Extract metadata in parallel
    metadata_rows: dict[str, dict[str, str]] = {}

    if gbff_files:
        with ThreadPoolExecutor(max_workers=threads) as executor:
            futures = {executor.submit(extract_metadata, p): p for p in gbff_files}
            for future in as_completed(futures):
                try:
                    meta = future.result()
                    metadata_rows[meta["file"]] = meta
                except Exception as e:
                    print(f"  Warning: {e}", file=sys.stderr)

    # 5. Merge metadata
    if metadata_rows:
        meta_df = pd.DataFrame.from_dict(metadata_rows, orient="index")
        meta_df.index.name = "sample"
        meta_cols = ["organism", "strain", "isolation_source", "host", "assembly_level"]
        meta_df = meta_df[[c for c in meta_cols if c in meta_df.columns]]
        matrix = meta_df.join(matrix, how="right")

    # 6. Output
    if output:
        matrix.to_csv(output, sep="\t")
        print(f"\nMatrix written to: {output}", file=sys.stderr)
    else:
        matrix.to_csv(sys.stdout, sep="\t")

    n_samples = len(matrix)
    metadata_cols = ("total_genes", "organism", "strain", "isolation_source", "host", "assembly_level")
    n_genes = len([c for c in matrix.columns if c not in metadata_cols])
    print(f"Samples: {n_samples}, Unique genes: {n_genes}", file=sys.stderr)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Build a presence/absence or count matrix of AMR genes per strain. "
            "Accepts AMRFinderPlus, RGI, ABRicate, or generic TSV output and "
            "optionally enriches with metadata from GBFF files."
        ),
    )

    parser.add_argument("amr_results", type=Path, help="AMR tool output file (TSV)")
    parser.add_argument("gbff_dir", type=Path, help="Directory containing GBFF/GBK files for metadata")
    parser.add_argument(
        "--format",
        dest="fmt",
        type=str,
        default="amrfinder",
        choices=["amrfinder", "rgi", "abricate", "tsv"],
        help="Input format (default: amrfinder)",
    )
    parser.add_argument(
        "--gene-col",
        type=str,
        default="gene",
        help="Column name for gene identifier (for --format tsv)",
    )
    parser.add_argument(
        "--sample-col",
        type=str,
        default="sample",
        help="Column name for sample identifier (for --format tsv)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Number of threads for parallel processing. Defaults to number of CPUs.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output TSV file path (default: stdout)",
    )

    args = parser.parse_args()

    main(
        amr_results=args.amr_results,
        gbff_dir=args.gbff_dir,
        fmt=args.fmt,
        gene_col=args.gene_col,
        sample_col=args.sample_col,
        threads=args.threads,
        output=args.output,
    )
