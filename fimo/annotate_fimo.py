## Annotate FIMO motif hits with nearby CDS metadata from GenBank records
#
## Usage:
#   $ python fimo/annotate_fimo.py -h
#   usage: annotate_fimo.py [-h] [--upstream-bp UPSTREAM_BP] [--email EMAIL]
#                           [--api-key API_KEY] [--tool TOOL] [--cache-dir CACHE_DIR]
#                           [--no-cache] [--threads THREADS] [--output OUTPUT]
#                           [--accession-col ACCESSION_COL] [--start-col START_COL]
#                           [--stop-col STOP_COL] [--strand-col STRAND_COL]
#                           [--keep-no-hit] [--no-deduplicate]
#                           input
#
#   Annotate FIMO TSV hits by finding CDS features on the same strand whose TSS
#   is within a configurable upstream window from each motif.

from __future__ import annotations

import argparse
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from http.client import HTTPConnection
from pathlib import Path

import pandas as pd
from Bio import Entrez, SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

# Work around occasional HTTP/1.1 issues with NCBI in some environments.
HTTPConnection._http_vsn = 10  # type: ignore
HTTPConnection._http_vsn_str = "HTTP/1.0"  # type: ignore

SIGN_TO_NUMBER = {"+": 1, "-": -1}


def normalise_description(description: str, description_regex: re.Pattern[str]) -> str:
    """Trim organism/contig description to a cleaner organism-like label."""
    first = description.split(",")[0]
    if description_regex.search(first):
        return " ".join(first.split(" ")[:-1]).strip()
    return first.strip()


def get_qualifier(feature: SeqFeature, key: str, default: str = "") -> str:
    """Safely read a feature qualifier value."""
    values = feature.qualifiers.get(key, [default])
    return values[0] if values else default


def configure_entrez(email: str, api_key: str | None, tool_name: str) -> None:
    """Configure Biopython Entrez client settings."""
    Entrez.email = email
    Entrez.tool = tool_name
    if api_key:
        Entrez.api_key = api_key
    Entrez.max_tries = 10


def fetch_genbank_record(
    accession: str,
    cache_dir: Path,
    use_cache: bool,
    retries: int = 3,
    retry_delay: float = 2.0,
) -> SeqRecord:
    """Fetch one GenBank record by accession, using optional local cache."""
    cache_path = cache_dir / f"{accession}.gb"

    if use_cache and cache_path.exists():
        return SeqIO.read(cache_path, "genbank")

    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            with Entrez.efetch(db="nuccore", id=accession, rettype="gbwithparts", retmode="text") as handle:
                record = SeqIO.read(handle, "genbank")

            if use_cache:
                cache_path.parent.mkdir(parents=True, exist_ok=True)
                SeqIO.write(record, cache_path, "genbank")

            return record
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            if attempt < retries:
                time.sleep(retry_delay)

    raise RuntimeError(f"Failed to fetch {accession}: {last_error}")


def check_motif_within_gene(record: SeqRecord, start_location: int, end_location: int) -> bool:
    """Return True if motif span is fully contained in a gene-like feature."""
    motif_start0 = min(start_location, end_location) - 1
    motif_end0 = max(start_location, end_location)

    for feature in record.features:
        if feature.type not in {"gene", "CDS"}:
            continue
        start = int(feature.location.start)
        end = int(feature.location.end)
        if motif_start0 >= start and motif_end0 <= end:
            return True

    return False


def find_candidate_cds(
    record: SeqRecord,
    motif_start: int,
    motif_stop: int,
    motif_strand: int,
    upstream_bp: int,
) -> SeqFeature | None:
    """Find first CDS on same strand with TSS in motif upstream search window."""
    start_1 = min(motif_start, motif_stop)
    stop_1 = max(motif_start, motif_stop)

    if motif_strand == -1:
        motif_space_start, motif_space_end = start_1 - upstream_bp, start_1
    else:
        motif_space_start, motif_space_end = stop_1, stop_1 + upstream_bp

    for feature in record.features:
        if feature.type != "CDS":
            continue
        if feature.location.strand != motif_strand:
            continue

        if motif_strand == 1:
            tss_location = int(feature.location.start) + 1
        else:
            tss_location = int(feature.location.end)

        if motif_space_start <= tss_location <= motif_space_end:
            return feature

    return None


def prefetch_records(
    accessions: set[str],
    cache_dir: Path,
    use_cache: bool,
    threads: int | None,
) -> tuple[dict[str, SeqRecord], set[str]]:
    """Download/cache records in parallel and return records plus failed IDs."""
    records: dict[str, SeqRecord] = {}
    failed: set[str] = set()

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {
            executor.submit(fetch_genbank_record, accession, cache_dir, use_cache): accession
            for accession in accessions
        }
        for future in as_completed(futures):
            accession = futures[future]
            try:
                records[accession] = future.result()
            except Exception as exc:  # noqa: BLE001
                print(f"Warning: could not load {accession}: {exc}", file=sys.stderr)
                failed.add(accession)

    return records, failed


def annotate_fimo_table(
    fimo: pd.DataFrame,
    records: dict[str, SeqRecord],
    accession_col: str,
    start_col: str,
    stop_col: str,
    strand_col: str,
    upstream_bp: int,
    keep_no_hit: bool,
    description_regex: re.Pattern[str],
) -> pd.DataFrame:
    """Annotate FIMO rows with nearby CDS metadata."""
    annotated_rows: list[dict[str, object]] = []

    for row in fimo.to_dict(orient="records"):
        accession = str(row[accession_col])
        if accession not in records:
            continue

        motif_start = int(row[start_col])
        motif_stop = int(row[stop_col])
        strand_str = str(row[strand_col])

        if strand_str not in SIGN_TO_NUMBER:
            print(f"Warning: unsupported strand '{strand_str}' for accession {accession}", file=sys.stderr)
            continue

        motif_strand = SIGN_TO_NUMBER[strand_str]
        record = records[accession]
        target_feature = find_candidate_cds(record, motif_start, motif_stop, motif_strand, upstream_bp)

        if target_feature is None and not keep_no_hit:
            continue

        description = normalise_description(record.description, description_regex)

        enriched = dict(row)
        enriched["organism"] = description
        enriched["in_gene"] = check_motif_within_gene(record, motif_start, motif_stop)

        if target_feature is None:
            enriched["gene_start"] = ""
            enriched["gene_stop"] = ""
            enriched["locus_tag"] = ""
            enriched["product"] = ""
            enriched["sequence"] = ""
            enriched["protein_sequence"] = ""
        else:
            enriched["gene_start"] = int(target_feature.location.start) + 1
            enriched["gene_stop"] = int(target_feature.location.end)
            enriched["locus_tag"] = get_qualifier(target_feature, "locus_tag")
            enriched["product"] = get_qualifier(target_feature, "product")
            enriched["sequence"] = str(target_feature.extract(record.seq))
            enriched["protein_sequence"] = get_qualifier(target_feature, "translation")

        annotated_rows.append(enriched)

    return pd.DataFrame(annotated_rows)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=("Annotate FIMO TSV hits by finding nearby CDS features from NCBI GenBank records.")
    )

    parser.add_argument("input", type=Path, help="Path to FIMO TSV file")
    parser.add_argument("--upstream-bp", type=int, default=100, help="Upstream window size in bp (default: 100)")
    parser.add_argument("--email", type=str, default="user@example.com", help="Email for NCBI Entrez")
    parser.add_argument("--api-key", type=str, default=None, help="NCBI API key")
    parser.add_argument("--tool", type=str, default="fimo-gene-annotator", help="Entrez tool name")
    parser.add_argument(
        "--cache-dir", type=Path, default=Path(".cache/entrez"), help="Directory for cached GenBank records"
    )
    parser.add_argument("--no-cache", action="store_true", help="Disable cache reads/writes")
    parser.add_argument("--threads", type=int, default=None, help="Worker count for accession prefetch")
    parser.add_argument(
        "--output", type=Path, default=None, help="Output TSV path (default: <input_stem>.annotated.tsv)"
    )

    parser.add_argument(
        "--accession-col",
        type=str,
        default="sequence_name",
        help="FIMO accession column name (default: sequence_name)",
    )
    parser.add_argument("--start-col", type=str, default="start", help="FIMO start column name (default: start)")
    parser.add_argument("--stop-col", type=str, default="stop", help="FIMO stop column name (default: stop)")
    parser.add_argument("--strand-col", type=str, default="strand", help="FIMO strand column name (default: strand)")
    parser.add_argument(
        "--description-regex",
        type=str,
        default=r"([Cc]ont|ctg|scaffold|chromosome|\.?\d{2,}_\d{1,}_\d{2,}\.?|_cov_|NODE|Scaf)",
        help="Regex used to detect contig-like description tails (default: current pattern)",
    )
    parser.add_argument("--keep-no-hit", action="store_true", help="Keep rows where no nearby CDS is found")
    parser.add_argument("--no-deduplicate", action="store_true", help="Do not deduplicate rows by locus_tag")

    return parser.parse_args()


def main() -> None:
    """Entry point."""
    args = parse_args()

    if not args.input.exists():
        print(f"Error: input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    configure_entrez(args.email, args.api_key, args.tool)

    fimo = pd.read_table(args.input, comment="#")

    required = {args.accession_col, args.start_col, args.stop_col, args.strand_col}
    missing = sorted(required - set(fimo.columns))
    if missing:
        print(f"Error: missing required column(s): {', '.join(missing)}", file=sys.stderr)
        sys.exit(1)

    try:
        description_re = re.compile(args.description_regex)
    except re.error as exc:
        print(f"Error: invalid description regex: {exc}", file=sys.stderr)
        sys.exit(1)

    accessions = {str(value) for value in fimo[args.accession_col].dropna().unique()}
    records, failed = prefetch_records(
        accessions=accessions,
        cache_dir=args.cache_dir,
        use_cache=not args.no_cache,
        threads=args.threads,
    )

    if not records:
        print("Error: no GenBank records were loaded.", file=sys.stderr)
        sys.exit(1)

    if failed:
        print(f"Warning: {len(failed)} accession(s) failed to load and will be skipped.", file=sys.stderr)

    annotated = annotate_fimo_table(
        fimo=fimo,
        records=records,
        accession_col=args.accession_col,
        start_col=args.start_col,
        stop_col=args.stop_col,
        strand_col=args.strand_col,
        upstream_bp=args.upstream_bp,
        keep_no_hit=args.keep_no_hit,
        description_regex=description_re,
    )

    if annotated.empty:
        print("Warning: no rows remained after annotation.", file=sys.stderr)

    if not args.no_deduplicate and "locus_tag" in annotated.columns and not annotated.empty:
        if "p-value" in annotated.columns:
            annotated = annotated.sort_values(by="p-value", ascending=True)
        elif "q-value" in annotated.columns:
            annotated = annotated.sort_values(by="q-value", ascending=True)
        annotated = annotated.drop_duplicates(subset=["locus_tag"], keep="first")

    output = args.output or args.input.with_suffix("").with_name(f"{args.input.stem}.annotated.tsv")
    annotated.to_csv(output, sep="\t", index=False)

    print(f"Wrote {len(annotated)} row(s) to {output}", file=sys.stderr)


if __name__ == "__main__":
    main()
