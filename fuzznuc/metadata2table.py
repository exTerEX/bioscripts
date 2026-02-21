## Given a fuzznuc TSV file, find nearby genes and add metadata
#
## Usage:
#   $ python fuzznuc/metadata2table.py -h
#   usage: metadata2table.py [-h] [--accession ACCESSION] [--distance DISTANCE]
#                            [--email EMAIL] [--api-key API_KEY] [--threads THREADS]
#                            input output
#
#   Given a fuzznuc TSV file, find nearby genes from NCBI GenBank records
#
#   positional arguments:
#     input                 Path to fuzznuc TSV file
#     output                Desired path/to/filename for the output TSV
#
#   options:
#     -h, --help            show this help message and exit
#     --accession ACCESSION NCBI accession to use (overrides SeqName from file)
#     --distance DISTANCE   Maximum distance in bps to search for genes (default: 200)
#     --email EMAIL         Email address for NCBI Entrez (required by NCBI)
#     --threads THREADS     Number of threads for parallel processing

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from Bio import Entrez, SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord


def fetch_genbank(accession: str, email: str) -> SeqRecord:
    """Fetch a GenBank record from NCBI Entrez.

    Args:
        accession: NCBI accession number.
        email: Email address for NCBI Entrez (required by NCBI).

    Returns:
        SeqRecord object containing the GenBank data.
    """
    with Entrez.efetch(db="nucleotide", id=accession, rettype="gbwithparts", retmode="text") as handle:
        return SeqIO.read(handle, "genbank")


def get_qualifier(feature: SeqFeature, key: str, default: str = "") -> str:
    """Safely extract a qualifier from a SeqFeature.

    Args:
        feature: BioPython SeqFeature object.
        key: Qualifier key to extract.
        default: Default value if key not found.

    Returns:
        The qualifier value or default.
    """
    values = feature.qualifiers.get(key, [default])
    return values[0] if values else default


def extract_gene_metadata(feature: SeqFeature, record: SeqRecord) -> dict[str, str | int]:
    """Extract metadata from a CDS feature.

    Args:
        feature: BioPython SeqFeature (CDS).
        record: Parent SeqRecord for sequence extraction.

    Returns:
        Dictionary containing gene metadata.
    """
    location = feature.location
    gene_seq = str(location.extract(record.seq))  # type: ignore

    return {
        "gene": get_qualifier(feature, "gene"),
        "gene_start": int(location.start) + 1,  # type: ignore ; Convert to 1-based
        "gene_end": int(location.end),  # type: ignore
        "gene_strand": "+" if location.strand == 1 else "-",  # type: ignore
        "locus_tag": get_qualifier(feature, "locus_tag"),
        "old_locus_tag": get_qualifier(feature, "old_locus_tag"),
        "EC_number": get_qualifier(feature, "EC_number"),
        "gene_sequence": gene_seq,
        "product": get_qualifier(feature, "product"),
        "translation": get_qualifier(feature, "translation"),
        "pseudo": "yes" if "pseudo" in feature.qualifiers else "no",
    }


def find_nearby_genes(
    record: SeqRecord,
    motif_start: int,
    motif_end: int,
    motif_strand: str,
    max_distance: int,
) -> list[dict[str, str | int]]:
    """Find genes on the same strand where motif is within distance of gene start.

    Args:
        record: SeqRecord containing the genome.
        motif_start: Start position of the motif (1-based).
        motif_end: End position of the motif (1-based).
        motif_strand: Strand of the motif ('+' or '-').
        max_distance: Maximum distance in bp from motif to gene start.

    Returns:
        List of dictionaries containing metadata for nearby genes.
    """
    nearby_genes: list[dict[str, str | int]] = []
    strand_int = 1 if motif_strand == "+" else -1

    # Convert to 0-based for BioPython comparison
    motif_start_0 = motif_start - 1
    motif_end_0 = motif_end

    for feature in record.features:
        if feature.type not in ("CDS"):
            continue

        # Check strand matches
        if feature.location.strand != strand_int:  # type: ignore
            continue

        # Gene start depends on strand
        # + strand: start is at location.start
        # - strand: start is at location.end (transcription goes backwards)
        if strand_int == 1:
            gene_start_pos = int(feature.location.start)  # type: ignore
            # Distance from motif end to gene start (motif should be upstream)
            distance = gene_start_pos - motif_end_0
        else:
            gene_start_pos = int(feature.location.end)  # type: ignore
            # Distance from gene start to motif start (motif should be downstream of gene end)
            distance = motif_start_0 - gene_start_pos

        # Only include if motif is upstream of gene start and within distance
        if 0 <= distance <= max_distance:
            nearby_genes.append(extract_gene_metadata(feature, record))

    return nearby_genes


def process_motif_row(
    row: dict[str, str],
    records_cache: dict[str, SeqRecord],
    accession_override: str | None,
    max_distance: int,
) -> list[dict[str, str | int]]:
    """Process a single motif row and find nearby genes.

    Args:
        row: Dictionary from the fuzznuc TSV file.
        records_cache: Cache of downloaded SeqRecords.
        accession_override: Optional accession to use instead of SeqName.
        max_distance: Maximum distance to search for genes.

    Returns:
        List of result dictionaries with motif and gene info combined.
        Empty list if no nearby genes are found.
    """
    accession = accession_override or row["SeqName"]
    record = records_cache[accession]

    motif_start = int(row["Start"])
    motif_end = int(row["End"])
    motif_strand = row["Strand"]

    # Handle reverse strand where Start > End
    if motif_start > motif_end:
        motif_start, motif_end = motif_end, motif_start

    # Extract motif sequence (0-based slicing)
    motif_seq = str(record.seq[motif_start - 1 : motif_end])  # type: ignore
    if motif_strand == "-":
        # Reverse complement for minus strand
        motif_seq = str(record.seq[motif_start - 1 : motif_end].reverse_complement())  # type: ignore

    nearby_genes = find_nearby_genes(record, motif_start, motif_end, motif_strand, max_distance)

    # Only return results if nearby genes are found
    if not nearby_genes:
        return []

    base_info = {
        "SeqName": row["SeqName"],
        "motif": motif_seq,
        "motif_start": motif_start,
        "motif_end": motif_end,
        "motif_strand": motif_strand,
        "score": row.get("Score", ""),
    }

    return [base_info | gene for gene in nearby_genes]


def main(
    input_path: Path,
    output_path: Path,
    accession: str | None = None,
    distance: int = 200,
    email: str = "user@example.com",
    api_key: str | None = None,
    threads: int | None = None,
):
    """Main function to process fuzznuc results and find nearby genes.

    Args:
        input_path: Path to fuzznuc TSV file.
        output_path: Path for output TSV file.
        accession: Optional accession override.
        distance: Maximum distance to search for genes.
        email: Email for NCBI Entrez.
        api_key: Optional NCBI API key for higher rate limits.
        threads: Number of threads for parallel processing.
    """
    # Configure Entrez
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    # Read input TSV
    with input_path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)

    if not rows:
        print("No motifs found in input file.")
        return

    # Determine which accessions to fetch
    accessions = {accession} if accession else {row["SeqName"] for row in rows}

    # Fetch GenBank records (parallelized)
    print(f"Fetching {len(accessions)} GenBank record(s) from NCBI...")
    records_cache: dict[str, SeqRecord] = {}

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(fetch_genbank, acc, email): acc for acc in accessions}
        for future in as_completed(futures):
            acc = futures[future]
            try:
                records_cache[acc] = future.result()
                print(f"  Downloaded: {acc}")
            except Exception as e:
                print(f"  Error fetching {acc}: {e}")

    # Process motif rows
    print(f"Processing {len(rows)} motif(s)...")
    all_results: list[dict[str, str | int]] = []

    for row in rows:
        results = process_motif_row(row, records_cache, accession, distance)
        all_results.extend(results)

    # Write output
    fieldnames = [
        "SeqName",
        "motif",
        "motif_start",
        "motif_end",
        "motif_strand",
        "score",
        "gene",
        "gene_start",
        "gene_end",
        "gene_strand",
        "locus_tag",
        "old_locus_tag",
        "EC_number",
        "product",
        "pseudo",
        "gene_sequence",
        "translation",
    ]

    with output_path.open("w") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(all_results)

    print(f"Output written to: {output_path}")
    print(f"Total results: {len(all_results)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Given a fuzznuc TSV file, find nearby genes from NCBI GenBank records"
    )

    parser.add_argument("input", type=Path, help="Path to fuzznuc TSV file")
    parser.add_argument("output", type=Path, help="Desired path/to/filename for the output TSV")
    parser.add_argument(
        "--accession",
        type=str,
        default=None,
        help="NCBI accession to use (overrides SeqName from file)",
    )
    parser.add_argument(
        "--distance",
        type=int,
        default=200,
        help="Maximum distance in bps to search for genes (default: 200)",
    )
    parser.add_argument(
        "--email",
        type=str,
        default="user@example.com",
        help="Email address for NCBI Entrez (required by NCBI)",
    )
    parser.add_argument(
        "--api-key",
        type=str,
        default=None,
        help="NCBI API key for higher rate limits",
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
        args.accession,
        args.distance,
        args.email,
        args.api_key,
        args.threads,
    )
