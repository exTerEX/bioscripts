## Extract upstream regions of genes from GenBank records
#
## Usage:
#   $ python misc/extract_upstream_gene_region.py -h
#   usage: extract_upstream_gene_region.py [-h] [--accession ACCESSION]
#                                          [--genbank GENBANK] [--locus-tags LOCUS_TAGS]
#                                          [--upstream UPSTREAM] [--include-start-codon]
#                                          [--output OUTPUT] [--email EMAIL]
#                                          [--api-key API_KEY] [--threads THREADS]
#
#   Extract upstream regions of genes from GenBank records
#
#   options:
#     -h, --help            show this help message and exit
#     --accession ACCESSION NCBI assembly accession (GCF_*/GCA_*) to download
#     --genbank GENBANK     Path to local GenBank file
#     --locus-tags LOCUS_TAGS
#                           Comma-separated list of locus_tags to extract
#     --upstream UPSTREAM   Number of base pairs upstream to extract (default: 200)
#     --include-start-codon Include the first 3 bp (start codon) of the gene
#     --output OUTPUT       Output FASTA file path. If not given, outputs to stdout
#     --email EMAIL         Email address for NCBI Entrez
#     --api-key API_KEY     NCBI API key for higher rate limits
#     --threads THREADS     Number of threads for parallel processing

import argparse
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord


def fetch_genbank(accession: str) -> SeqRecord:
    """Fetch a GenBank record from NCBI Entrez.

    Args:
        accession: NCBI nucleotide accession number.

    Returns:
        SeqRecord object containing the GenBank data.
    """
    with Entrez.efetch(
        db="nucleotide", id=accession, rettype="gbwithparts", retmode="text"
    ) as handle:
        return SeqIO.read(handle, "genbank")


def get_nucleotide_accessions_from_assembly(assembly_accession: str) -> list[str]:
    """Get all nucleotide accessions for an assembly accession.

    Args:
        assembly_accession: NCBI assembly accession (GCF_*/GCA_*).

    Returns:
        List of nucleotide accession numbers.
    """
    # Search for the assembly to get its UID
    with Entrez.esearch(db="assembly", term=assembly_accession, retmax=1) as handle:
        search_result = Entrez.read(handle)

    if not search_result["IdList"]:
        raise ValueError(f"Assembly not found: {assembly_accession}")

    assembly_uid = search_result["IdList"][0]

    # Link assembly to nucleotide database
    with Entrez.elink(dbfrom="assembly", db="nucleotide", id=assembly_uid, linkname="assembly_nuccore_refseq") as handle:
        link_result = Entrez.read(handle)

    nucleotide_ids: list[str] = []

    # Try RefSeq link first
    if link_result and link_result[0].get("LinkSetDb"):
        for link_set in link_result[0]["LinkSetDb"]:
            if link_set.get("Link"):
                nucleotide_ids = [link["Id"] for link in link_set["Link"]]
                break

    # If no RefSeq, try GenBank link
    if not nucleotide_ids:
        with Entrez.elink(dbfrom="assembly", db="nucleotide", id=assembly_uid, linkname="assembly_nuccore_insdc") as handle:
            link_result = Entrez.read(handle)

        if link_result and link_result[0].get("LinkSetDb"):
            for link_set in link_result[0]["LinkSetDb"]:
                if link_set.get("Link"):
                    nucleotide_ids = [link["Id"] for link in link_set["Link"]]
                    break

    if not nucleotide_ids:
        raise ValueError(f"No nucleotide records found for assembly: {assembly_accession}")

    # Fetch accession numbers from the UIDs
    with Entrez.efetch(db="nucleotide", id=",".join(nucleotide_ids), rettype="acc", retmode="text") as handle:
        accessions = [line.strip() for line in handle.read().strip().split("\n") if line.strip()]

    return accessions


def get_qualifier(feature: SeqFeature, key: str, default: str = "") -> str:
    """Safely extract a qualifier from a SeqFeature."""
    values = feature.qualifiers.get(key, [default])
    return values[0] if values else default


def extract_upstream_region(
    feature: SeqFeature,
    record: SeqRecord,
    upstream_bp: int,
    include_start_codon: bool,
) -> tuple[str, str, str, int, int] | None:
    """Extract the upstream region of a CDS feature.

    Args:
        feature: CDS SeqFeature.
        record: Parent SeqRecord.
        upstream_bp: Number of base pairs upstream to extract.
        include_start_codon: Whether to include the first 3 bp of the gene.

    Returns:
        Tuple of (sequence, gene_name, locus_tag, start, end) or None if extraction fails.
    """
    location = feature.location
    strand = location.strand
    seq_len = len(record.seq)

    # Get gene name and locus_tag
    gene_name = get_qualifier(feature, "gene")
    locus_tag = get_qualifier(feature, "locus_tag")
    
    # Need at least locus_tag to proceed
    if not locus_tag and not gene_name:
        return None

    if strand == 1:  # Forward strand
        gene_start = int(location.start)  # 0-based
        # Upstream region ends at gene start (or +3 if including start codon)
        end_pos = gene_start + 3 if include_start_codon else gene_start
        start_pos = max(0, gene_start - upstream_bp)
        
        upstream_seq = str(record.seq[start_pos:end_pos])
        # Report 1-based coordinates
        region_start = start_pos + 1
        region_end = end_pos

    else:  # Reverse strand (-1)
        gene_end = int(location.end)  # 0-based exclusive (so this is the "start" for - strand)
        # Upstream region starts at gene end (or -3 if including start codon)
        start_pos = gene_end - 3 if include_start_codon else gene_end
        end_pos = min(seq_len, gene_end + upstream_bp)
        
        upstream_seq = str(record.seq[start_pos:end_pos].reverse_complement())
        # Report 1-based coordinates (in genomic sense)
        region_start = start_pos + 1
        region_end = end_pos

    return upstream_seq, gene_name, locus_tag, region_start, region_end


def process_record(
    record: SeqRecord,
    locus_tags: set[str] | None,
    upstream_bp: int,
    include_start_codon: bool,
) -> list[SeqRecord]:
    """Process a GenBank record and extract upstream regions.

    Args:
        record: SeqRecord to process.
        locus_tags: Set of locus_tags to extract, or None for all CDS.
        upstream_bp: Number of base pairs upstream.
        include_start_codon: Whether to include start codon.

    Returns:
        List of SeqRecords with upstream sequences.
    """
    results: list[SeqRecord] = []
    accession = record.id

    for feature in record.features:
        if feature.type != "CDS":
            continue

        # Filter by locus_tag if specified
        if locus_tags:
            feature_locus = get_qualifier(feature, "locus_tag")
            if feature_locus not in locus_tags:
                continue

        extraction = extract_upstream_region(
            feature, record, upstream_bp, include_start_codon
        )

        if extraction is None:
            continue

        seq, gene_name, locus_tag, start, end = extraction

        # Create FASTA header
        # Format: >gene_name [locus_tag=*] [accession=*] [location=*]
        # Or:     >locus_tag [accession=*] [location=*]
        if gene_name:
            header = f"{gene_name} [locus_tag={locus_tag}] [accession={accession}] [location={start}-{end}]"
        else:
            header = f"{locus_tag} [accession={accession}] [location={start}-{end}]"

        fasta_record = SeqRecord(
            Seq(seq),
            id=header,
            description="",
        )
        results.append(fasta_record)

    return results


def process_nucleotide_accession(
    accession: str,
    locus_tags: set[str] | None,
    upstream_bp: int,
    include_start_codon: bool,
) -> list[SeqRecord]:
    """Fetch and process a single nucleotide accession.

    Args:
        accession: NCBI nucleotide accession to fetch.
        locus_tags: Set of locus_tags to extract, or None for all CDS.
        upstream_bp: Number of base pairs upstream.
        include_start_codon: Whether to include start codon.

    Returns:
        List of SeqRecords with upstream sequences.
    """
    try:
        record = fetch_genbank(accession)
        return process_record(record, locus_tags, upstream_bp, include_start_codon)
    except Exception as e:
        print(f"Error processing {accession}: {e}", file=sys.stderr)
        return []


def main(
    accession: str | None = None,
    genbank_path: Path | None = None,
    locus_tags_str: str | None = None,
    upstream_bp: int = 200,
    include_start_codon: bool = False,
    output_path: Path | None = None,
    email: str = "user@example.com",
    api_key: str | None = None,
    threads: int | None = None,
):
    """Main function to extract upstream regions.

    Args:
        accession: NCBI assembly accession (GCF_*/GCA_*) to download.
        genbank_path: Path to local GenBank file.
        locus_tags_str: Comma-separated locus_tags to extract.
        upstream_bp: Number of base pairs upstream.
        include_start_codon: Whether to include the first 3 bp.
        output_path: Output FASTA file path, or None for stdout.
        email: Email for NCBI Entrez.
        api_key: NCBI API key.
        threads: Number of threads for parallel processing.
    """
    # Configure Entrez
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    # Parse locus_tags
    locus_tags: set[str] | None = None
    if locus_tags_str:
        locus_tags = {tag.strip() for tag in locus_tags_str.split(",") if tag.strip()}

    all_results: list[SeqRecord] = []

    if genbank_path:
        # Process local GenBank file
        print(f"Reading GenBank file: {genbank_path}", file=sys.stderr)
        with genbank_path.open() as f:
            for record in SeqIO.parse(f, "genbank"):
                results = process_record(
                    record, locus_tags, upstream_bp, include_start_codon
                )
                all_results.extend(results)

    elif accession:
        # Get nucleotide accessions from assembly accession
        print(f"Fetching nucleotide accessions for assembly: {accession}", file=sys.stderr)
        try:
            nucleotide_accessions = get_nucleotide_accessions_from_assembly(accession)
        except Exception as e:
            print(f"Error: {e}", file=sys.stderr)
            return

        print(f"Found {len(nucleotide_accessions)} nucleotide record(s)", file=sys.stderr)

        if len(nucleotide_accessions) == 1:
            # Single contig, no threading needed
            print(f"Processing: {nucleotide_accessions[0]}\n", file=sys.stderr)
            all_results = process_nucleotide_accession(
                nucleotide_accessions[0], locus_tags, upstream_bp, include_start_codon
            )
        else:
            # Multiple contigs, use threading
            print(f"Processing {len(nucleotide_accessions)} contigs with threading...\n", file=sys.stderr)

            with ThreadPoolExecutor(max_workers=threads) as executor:
                futures = {
                    executor.submit(
                        process_nucleotide_accession, acc, locus_tags, upstream_bp, include_start_codon
                    ): acc
                    for acc in nucleotide_accessions
                }

                for future in as_completed(futures):
                    acc = futures[future]
                    try:
                        results = future.result()
                        all_results.extend(results)
                        print(f"  Processed: {acc} ({len(results)} sequences)", file=sys.stderr)
                    except Exception as e:
                        print(f"  Error processing {acc}: {e}", file=sys.stderr)

    else:
        raise ValueError("Either --accession or --genbank must be provided")

    # Output results
    if not all_results:
        print("No upstream regions extracted.", file=sys.stderr)
        return

    if output_path:
        with output_path.open("w") as f:
            SeqIO.write(all_results, f, "fasta")
        print(f"Output written to: {output_path}", file=sys.stderr)
    else:
        SeqIO.write(all_results, sys.stdout, "fasta")

    print(f"\nTotal sequences: {len(all_results)}", file=sys.stderr)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract upstream regions of genes from GenBank records"
    )

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--accession",
        type=str,
        help="NCBI assembly accession (GCF_*/GCA_*) to download",
    )
    input_group.add_argument(
        "--genbank",
        type=Path,
        help="Path to local GenBank file",
    )

    parser.add_argument(
        "--locus-tags",
        type=str,
        default=None,
        help="Comma-separated list of locus_tags to extract. If not given, extracts all CDS",
    )
    parser.add_argument(
        "--upstream",
        type=int,
        default=200,
        help="Number of base pairs upstream to extract (default: 200)",
    )
    parser.add_argument(
        "--include-start-codon",
        action="store_true",
        help="Include the first 3 bp (start codon) of the gene",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output FASTA file path. If not given, outputs to stdout",
    )
    parser.add_argument(
        "--email",
        type=str,
        default="user@example.com",
        help="Email address for NCBI Entrez",
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
        accession=args.accession,
        genbank_path=args.genbank,
        locus_tags_str=args.locus_tags,
        upstream_bp=args.upstream,
        include_start_codon=args.include_start_codon,
        output_path=args.output,
        email=args.email,
        api_key=args.api_key,
        threads=args.threads,
    )
