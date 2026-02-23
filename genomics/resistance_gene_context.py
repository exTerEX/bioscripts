## Extract genomic context around resistance genes from GenBank files
#
## Usage:
#   $ python genomics/resistance_gene_context.py -h
#   usage: resistance_gene_context.py [-h] [--window WINDOW] [--threads THREADS]
#                                     [--output OUTPUT]
#                                     gbff_dir locus_tags
#
#   Extract the genomic neighbourhood around resistance genes. Reports flanking
#   genes, mobile genetic element markers (transposases, integrases, IS elements),
#   and gene organisation within a configurable window.
#
#   positional arguments:
#     gbff_dir            Directory containing GBFF/GBK files
#     locus_tags          File with locus tags of interest (one per line), or
#                         comma-separated list
#
#   options:
#     -h, --help          show this help message and exit
#     --window WINDOW     Context window size in base pairs on each side (default: 10000)
#     --threads THREADS   Number of threads for parallel processing
#     --output OUTPUT     Output TSV file path (default: stdout)

from __future__ import annotations

import argparse
import csv
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

MGE_KEYWORDS: set[str] = {
    "transposase",
    "integrase",
    "recombinase",
    "resolvase",
    "insertion sequence",
    "is element",
    "conjugal transfer",
    "mobilization",
    "relaxase",
    "type iv secretion",
}

IS_ELEMENT_RE = re.compile(r"\bIS\d+\b|IS[A-Z][a-z]+\d*", re.IGNORECASE)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def get_qualifier(feature: SeqFeature, key: str, default: str = "") -> str:
    """Safely extract the first qualifier value from a SeqFeature."""
    values = feature.qualifiers.get(key, [default])
    return values[0] if values else default


def is_mge_feature(feature: SeqFeature) -> bool:
    """Check whether a feature looks like a mobile genetic element component.

    Args:
        feature: A BioPython SeqFeature.

    Returns:
        True if any qualifier suggests an MGE association.
    """
    for key in ("product", "note", "function"):
        text = get_qualifier(feature, key).lower()
        if any(kw in text for kw in MGE_KEYWORDS):
            return True
        if IS_ELEMENT_RE.search(text):
            return True

    return feature.type == "mobile_element"


def feature_dict(feature: SeqFeature, record_id: str, is_target: bool = False) -> dict[str, str]:
    """Convert a SeqFeature into a flat dictionary row.

    Args:
        feature: BioPython SeqFeature.
        record_id: Parent record accession.
        is_target: Whether this feature is the target resistance gene.

    Returns:
        Dictionary suitable for TSV output.
    """
    location = feature.location
    strand = location.strand  # type: ignore[union-attr]
    start = int(location.start) + 1  # type: ignore[union-attr]
    end = int(location.end)  # type: ignore[union-attr]

    return {
        "record_id": record_id,
        "locus_tag": get_qualifier(feature, "locus_tag"),
        "gene": get_qualifier(feature, "gene"),
        "product": get_qualifier(feature, "product"),
        "feature_type": feature.type,
        "start": str(start),
        "end": str(end),
        "strand": "+" if strand == 1 else "-" if strand == -1 else "?",
        "is_pseudo": str("pseudo" in feature.qualifiers or "pseudogene" in feature.qualifiers),
        "is_mge": str(is_mge_feature(feature)),
        "is_target": str(is_target),
    }


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------


def extract_context(
    record: SeqRecord,
    target_locus_tags: set[str],
    window: int,
) -> list[dict[str, str]]:
    """Extract genomic context around target locus tags in a single record.

    Args:
        record: GenBank SeqRecord.
        target_locus_tags: Set of locus tags to search for.
        window: Number of base pairs on each side of the target gene.

    Returns:
        List of row dictionaries for every feature in the context windows.
    """
    record_id = record.id or record.name or ""

    # Collect CDS / gene-level features with their locus_tags
    annotated_features: list[SeqFeature] = []
    target_features: list[SeqFeature] = []

    for feat in record.features:
        if feat.type not in ("CDS", "gene", "rRNA", "tRNA", "ncRNA", "mobile_element", "misc_feature"):
            continue
        annotated_features.append(feat)
        lt = get_qualifier(feat, "locus_tag")
        if lt in target_locus_tags:
            target_features.append(feat)

    if not target_features:
        return []

    # For each target, collect features within the window
    rows: list[dict[str, str]] = []
    seen_ranges: set[tuple[int, int, str]] = set()

    for target in target_features:
        t_start = int(target.location.start)  # type: ignore[union-attr]
        t_end = int(target.location.end)  # type: ignore[union-attr]

        win_start = max(0, t_start - window)
        win_end = min(len(record.seq), t_end + window)  # type: ignore

        range_key = (win_start, win_end, get_qualifier(target, "locus_tag"))
        if range_key in seen_ranges:
            continue
        seen_ranges.add(range_key)

        target_lt = get_qualifier(target, "locus_tag")

        for feat in annotated_features:
            f_start = int(feat.location.start)  # type: ignore[union-attr]
            f_end = int(feat.location.end)  # type: ignore[union-attr]

            # Check overlap with window
            if f_end < win_start or f_start > win_end:
                continue

            is_target = get_qualifier(feat, "locus_tag") == target_lt and feat.type == target.type
            row = feature_dict(feat, record_id, is_target=is_target)
            row["target_locus_tag"] = target_lt
            row["distance_to_target"] = str(
                min(abs(f_start - t_start), abs(f_end - t_end), abs(f_start - t_end), abs(f_end - t_start))
                if not is_target
                else 0
            )
            rows.append(row)

    return rows


def process_gbff(
    path: Path,
    target_locus_tags: set[str],
    window: int,
) -> list[dict[str, str]]:
    """Process a single GBFF file and return context rows.

    Args:
        path: Path to a GenBank file.
        target_locus_tags: Locus tags to search for.
        window: Context window in bp.

    Returns:
        List of context row dictionaries.
    """
    rows: list[dict[str, str]] = []
    try:
        for record in SeqIO.parse(path, "genbank"):
            rows.extend(extract_context(record, target_locus_tags, window))
    except Exception as e:
        print(f"Warning: Could not parse {path}: {e}", file=sys.stderr)
    return rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

_FIELDNAMES = [
    "target_locus_tag",
    "record_id",
    "locus_tag",
    "gene",
    "product",
    "feature_type",
    "start",
    "end",
    "strand",
    "distance_to_target",
    "is_target",
    "is_pseudo",
    "is_mge",
]


def main(
    gbff_dir: Path,
    locus_tags_input: str,
    window: int = 10_000,
    threads: int | None = None,
    output: Path | None = None,
) -> None:
    """Extract genomic context around resistance genes.

    Args:
        gbff_dir: Directory containing GBFF/GBK files.
        locus_tags_input: Path to file with locus tags or comma-separated string.
        window: Context window in bp on each side.
        threads: Number of threads.
        output: Output TSV path, or None for stdout.
    """
    # Parse locus tags
    locus_tags_path = Path(locus_tags_input)
    if locus_tags_path.is_file():
        target_locus_tags = {line.strip() for line in locus_tags_path.read_text().splitlines() if line.strip()}
    else:
        target_locus_tags = {tag.strip() for tag in locus_tags_input.split(",") if tag.strip()}

    if not target_locus_tags:
        print("Error: No locus tags provided.", file=sys.stderr)
        sys.exit(1)

    print(f"Searching for {len(target_locus_tags)} locus tag(s) with {window} bp window", file=sys.stderr)

    # Gather GBFF files
    gbff_files = list(gbff_dir.glob("*.gbff")) + list(gbff_dir.glob("*.gbk")) + list(gbff_dir.glob("*.gb"))
    if not gbff_files:
        print(f"Error: No GenBank files found in {gbff_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Processing {len(gbff_files)} file(s)...", file=sys.stderr)

    all_rows: list[dict[str, str]] = []

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(process_gbff, p, target_locus_tags, window): p for p in gbff_files}
        for future in as_completed(futures):
            path = futures[future]
            try:
                rows = future.result()
                all_rows.extend(rows)
                if rows:
                    print(f"  {path.name}: {len(rows)} features in context", file=sys.stderr)
            except Exception as e:
                print(f"  Error processing {path.name}: {e}", file=sys.stderr)

    # Sort by target then by position
    all_rows.sort(key=lambda r: (r["target_locus_tag"], r["record_id"], int(r["start"])))

    # Output
    out_handle = output.open("w", newline="") if output else sys.stdout
    writer = csv.DictWriter(out_handle, fieldnames=_FIELDNAMES, delimiter="\t", extrasaction="ignore")
    writer.writeheader()
    writer.writerows(all_rows)

    if output:
        out_handle.close()
        print(f"\nResults written to: {output}", file=sys.stderr)

    # Summary
    targets_found = {r["target_locus_tag"] for r in all_rows if r["is_target"] == "True"}
    mge_count = sum(1 for r in all_rows if r["is_mge"] == "True")
    print(f"Targets found: {len(targets_found)}/{len(target_locus_tags)}", file=sys.stderr)
    print(f"Total context features: {len(all_rows)}", file=sys.stderr)
    print(f"MGE-associated features: {mge_count}", file=sys.stderr)

    missing = target_locus_tags - targets_found
    if missing:
        print(f"Missing locus tags: {', '.join(sorted(missing))}", file=sys.stderr)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Extract the genomic neighbourhood around resistance genes. "
            "Reports flanking genes, mobile genetic element markers "
            "(transposases, integrases, IS elements), and gene organisation "
            "within a configurable window."
        ),
    )

    parser.add_argument("gbff_dir", type=Path, help="Directory containing GBFF/GBK files")
    parser.add_argument(
        "locus_tags",
        type=str,
        help="File with locus tags of interest (one per line), or comma-separated list",
    )
    parser.add_argument(
        "--window",
        type=int,
        default=10_000,
        help="Context window size in base pairs on each side (default: 10000)",
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
        gbff_dir=args.gbff_dir,
        locus_tags_input=args.locus_tags,
        window=args.window,
        threads=args.threads,
        output=args.output,
    )
