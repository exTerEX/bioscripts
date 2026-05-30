## Fetch and merge eHOMD taxon/genome metadata, then export to Excel
#
## Usage:
#   $ python genomics/fetch_homd_metadata.py
#
#   Download taxon and genome metadata tables from eHOMD using the built-in
#   HOMD URLs, merge them on oral taxon IDs, derive a Gram-positive subset,
#   and export all tables to one Excel workbook.

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

DEFAULT_TAXON_URL = "https://www.homd.org/taxa/dld_table_all/browser"
DEFAULT_GENOME_URL = "https://www.homd.org/genome/dld_table_all/browser"
DEFAULT_OUTPUT = Path("data/eHOMD-metadata.xlsx")
DEFAULT_PHYLA = {"Actinobacteria", "Firmicutes"}

TAXON_DROP_COLUMNS = {
    "Warning",
    "Clone_count",
    "Clone_%",
    "Clone_rank",
    "Synonyms",
    "Genome_ID",
    "General_info",
    "Cultivability",
    "Phenotypic_characteristics",
    "Prevalence",
    "Disease",
    "References",
}

GENOME_DROP_COLUMNS = {
    "Oral Pathogen",
    "NCBI Genome-ID",
    "atcc_mn",
    "non_atcc_mn",
    "Genbank Acc no.",
    "16S rRNA",
    "16S rRNA Comment",
    "flag_id",
}

FINAL_COLUMNS = [
    "Oral Taxon-ID",
    "NCBI Taxon-ID",
    "Domain",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species",
    "Status",
    "Isolate Origin",
    "Body_site",
    "Genbank Assembly",
    "Genome-ID",
    "NCBI BioProject-ID",
    "NCBI BioSample-ID",
    "No. Contigs",
    "Sequencing Center",
    "Culture Collection",
    "Type_strain",
    "GC %",
    "16S_rRNA",
    "NCBI_pubmed_count",
    "NCBI_nucleotide_count",
    "NCBI_protein_count",
]


def load_homd_table(url: str) -> pd.DataFrame:
    """Load a HOMD browser TSV, skipping the banner row."""
    return pd.read_csv(url, sep="\t", skiprows=1, dtype=str)


def drop_known_columns(frame: pd.DataFrame, columns_to_drop: set[str]) -> pd.DataFrame:
    """Drop optional columns when present."""
    present = [column for column in columns_to_drop if column in frame.columns]
    return frame.drop(columns=present)


def build_merged_metadata(taxon: pd.DataFrame, genome: pd.DataFrame) -> pd.DataFrame:
    """Merge cleaned HOMD taxon and genome tables into one metadata table."""
    taxon_clean = drop_known_columns(taxon, TAXON_DROP_COLUMNS)
    genome_clean = drop_known_columns(genome, GENOME_DROP_COLUMNS)

    merged = genome_clean.merge(
        taxon_clean,
        how="outer",
        left_on="Oral_Taxon-ID",
        right_on="HMT_ID",
        suffixes=(".x", ".y"),
    )

    merged = merged.dropna(subset=["Oral_Taxon-ID", "HMT_ID"])

    merged = merged.drop(
        columns=[
            column
            for column in ["HMT_ID", "Genus.x", "Species.x", "Status.x", "Total Length"]
            if column in merged.columns
        ]
    )

    merged = merged.rename(
        columns={
            "Genus.y": "Genus",
            "Species.y": "Species",
            "Status.y": "Status",
            "Oral_Taxon-ID": "Oral Taxon-ID",
            "NCBI_taxon_id": "NCBI Taxon-ID",
        }
    )

    missing_columns = [column for column in FINAL_COLUMNS if column not in merged.columns]
    if missing_columns:
        raise ValueError(f"Merged HOMD table is missing expected column(s): {', '.join(missing_columns)}")

    return merged.loc[:, FINAL_COLUMNS]


def filter_gram_positive(frame: pd.DataFrame, phyla: set[str]) -> pd.DataFrame:
    """Return rows whose phylum belongs to the provided set."""
    return frame[frame["Phylum"].isin(phyla)].copy()


def write_workbook(
    taxon: pd.DataFrame,
    genome: pd.DataFrame,
    merged: pd.DataFrame,
    gram_positive: pd.DataFrame,
    output: Path,
) -> None:
    """Write all result tables to a single Excel workbook."""
    output.parent.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(output, engine="openpyxl") as writer:
        taxon.to_excel(writer, sheet_name="HOMD Taxon", index=False)
        genome.to_excel(writer, sheet_name="HOMD Genome", index=False)
        merged.to_excel(writer, sheet_name="HOMD Merged", index=False)
        gram_positive.to_excel(writer, sheet_name="HOMD Gram-positives", index=False)


def main() -> None:
    """Entry point."""
    if not DEFAULT_PHYLA:
        print("Error: DEFAULT_PHYLA must contain at least one phylum.", file=sys.stderr)
        sys.exit(1)

    try:
        taxon = load_homd_table(DEFAULT_TAXON_URL)
        genome = load_homd_table(DEFAULT_GENOME_URL)
        merged = build_merged_metadata(taxon, genome)
        gram_positive = filter_gram_positive(merged, DEFAULT_PHYLA)
        write_workbook(taxon, genome, merged, gram_positive, DEFAULT_OUTPUT)
    except Exception as exc:  # noqa: BLE001
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)

    print(f"Wrote workbook to {DEFAULT_OUTPUT}", file=sys.stderr)


if __name__ == "__main__":
    main()
