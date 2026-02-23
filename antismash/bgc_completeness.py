## Assess BGC completeness by comparing antiSMASH results against a reference gene list
#
## Usage:
#   $ python antismash/bgc_completeness.py -h
#   usage: bgc_completeness.py [-h] [--similarity SIMILARITY] [--threads THREADS]
#                              [--output OUTPUT]
#                              antismash_dir reference_genes
#
#   Check BGC completeness by comparing each detected cluster's gene content
#   against a user-supplied reference gene list (e.g. from MIBiG). Reports
#   which core biosynthetic genes are present, partial, or missing per cluster.
#
#   positional arguments:
#     antismash_dir         Directory containing antiSMASH result directories
#     reference_genes       TSV file with reference gene list (columns: gene, product, role)
#
#   options:
#     -h, --help            show this help message and exit
#     --similarity SIMILARITY Minimum name similarity (0-1) for fuzzy gene matching (default: 0.8)
#     --threads THREADS     Number of threads for parallel processing
#     --output OUTPUT       Output TSV file path (default: stdout)

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from difflib import SequenceMatcher
from pathlib import Path

# ---------------------------------------------------------------------------
# Reference gene loading
# ---------------------------------------------------------------------------


def load_reference_genes(path: Path) -> list[dict[str, str]]:
    """Load a reference gene list from a TSV file.

    Expected columns: gene (required), product (optional), role (optional).
    The role column can contain values like ``core``, ``accessory``,
    ``regulatory``, ``transport``, ``resistance``, ``tailoring``.

    Args:
        path: Path to the reference TSV file.

    Returns:
        List of dictionaries with gene, product, and role fields.
    """
    genes: list[dict[str, str]] = []
    with path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene = row.get("gene", "").strip()
            if not gene:
                continue
            genes.append(
                {
                    "gene": gene,
                    "product": row.get("product", "").strip(),
                    "role": row.get("role", "").strip(),
                }
            )

    if not genes:
        print(f"Error: No genes found in {path}", file=sys.stderr)
        sys.exit(1)

    return genes


# ---------------------------------------------------------------------------
# Gene matching
# ---------------------------------------------------------------------------


def _normalise(name: str) -> str:
    """Normalise a gene or product name for comparison."""
    return re.sub(r"[^a-z0-9]", "", name.lower())


def find_best_match(
    query: str,
    candidates: list[str],
    threshold: float,
) -> tuple[str, float] | None:
    """Find the best fuzzy match for *query* among *candidates*.

    Args:
        query: Gene name or product to search for.
        candidates: List of gene names / products from the BGC.
        threshold: Minimum similarity ratio (0-1).

    Returns:
        Tuple of (best_match, similarity) or None if nothing exceeds the
        threshold.
    """
    if not candidates:
        return None

    norm_query = _normalise(query)
    best_match = ""
    best_score = 0.0

    for cand in candidates:
        norm_cand = _normalise(cand)

        # Exact normalised match
        if norm_query == norm_cand:
            return cand, 1.0

        # Substring containment
        if norm_query in norm_cand or norm_cand in norm_query:
            score = max(len(norm_query), len(norm_cand)) / (len(norm_query) + len(norm_cand)) * 2
            score = max(score, 0.9)  # Substring match is high confidence
            if score > best_score:
                best_score = score
                best_match = cand
            continue

        # Fuzzy match
        score = SequenceMatcher(None, norm_query, norm_cand).ratio()
        if score > best_score:
            best_score = score
            best_match = cand

    if best_score >= threshold:
        return best_match, best_score

    return None


# ---------------------------------------------------------------------------
# antiSMASH JSON parsing
# ---------------------------------------------------------------------------


def extract_bgc_genes(record: dict) -> list[dict[str, list[dict[str, str]]]]:
    """Extract gene lists for each BGC region in an antiSMASH record.

    Args:
        record: A single antiSMASH record dictionary.

    Returns:
        List of dicts, one per region, each containing a ``genes`` list with
        gene/product/locus_tag/type info, plus region metadata.
    """
    regions: list[dict] = []

    if not record.get("areas"):
        return regions

    region_features = [f for f in record["features"] if f["type"] == "region"]

    for region_feat in region_features:
        quals = region_feat["qualifiers"]
        start_str, end_str = re.findall(r"\d+", region_feat["location"])
        region_start = int(start_str)
        region_end = int(end_str)

        # Collect CDS features within this region
        gene_list: list[dict[str, str]] = []
        for feat in record["features"]:
            if feat["type"] != "CDS":
                continue

            # Parse feature location
            loc_nums = re.findall(r"\d+", feat["location"])
            if len(loc_nums) < 2:
                continue
            feat_start = int(loc_nums[0])
            feat_end = int(loc_nums[-1])

            # Check if feature is within the region
            if feat_end < region_start or feat_start > region_end:
                continue

            fq = feat.get("qualifiers", {})
            gene_list.append(
                {
                    "gene": fq.get("gene", [""])[0],
                    "product": fq.get("product", [""])[0],
                    "locus_tag": fq.get("locus_tag", [""])[0],
                    "gene_kind": fq.get("gene_kind", [""])[0],
                }
            )

        regions.append(
            {
                "record_id": record["name"],
                "region": quals["region_number"][0],
                "start": start_str,
                "end": end_str,
                "product": " / ".join(quals.get("product", [])),
                "contig_edge": quals.get("contig_edge", [""])[0],
                "genes": gene_list,
            }
        )

    return regions


def parse_antismash_json(path: Path) -> list[dict]:
    """Parse an antiSMASH JSON and return all BGC regions with gene lists.

    Args:
        path: Path to the antiSMASH JSON file.

    Returns:
        List of region dictionaries.
    """
    with path.open() as f:
        data = json.load(f)

    all_regions: list[dict] = []
    for record in data["records"]:
        regions = extract_bgc_genes(record)
        for r in regions:
            r["file"] = path.stem  # type: ignore
        all_regions.extend(regions)

    return all_regions


# ---------------------------------------------------------------------------
# Completeness assessment
# ---------------------------------------------------------------------------


def assess_completeness(
    region: dict,
    reference_genes: list[dict[str, str]],
    similarity_threshold: float,
) -> list[dict[str, str]]:
    """Assess completeness of a BGC region against reference genes.

    Args:
        region: Region dictionary from antiSMASH parsing.
        reference_genes: List of reference gene dicts.
        similarity_threshold: Minimum similarity for matching.

    Returns:
        List of result rows, one per reference gene.
    """
    # Build candidate lists from BGC genes
    bgc_gene_names = [g["gene"] for g in region["genes"] if g["gene"]]
    bgc_products = [g["product"] for g in region["genes"] if g["product"]]

    # Map gene names/products to their full info
    gene_info_map: dict[str, dict[str, str]] = {}
    for g in region["genes"]:
        if g["gene"]:
            gene_info_map[g["gene"]] = g
        if g["product"]:
            gene_info_map[g["product"]] = g

    rows: list[dict[str, str]] = []
    matched_bgc_genes: set[str] = set()

    for ref in reference_genes:
        row: dict[str, str] = {
            "file": region["file"],
            "record_id": region["record_id"],
            "region": region["region"],
            "bgc_product": region["product"],
            "contig_edge": region["contig_edge"],
            "ref_gene": ref["gene"],
            "ref_product": ref["product"],
            "ref_role": ref["role"],
            "status": "missing",
            "matched_to": "",
            "match_similarity": "",
            "matched_locus_tag": "",
            "matched_gene_kind": "",
        }

        # Try matching by gene name first, then product
        match = find_best_match(ref["gene"], bgc_gene_names, similarity_threshold)
        if not match and ref["product"]:
            match = find_best_match(ref["product"], bgc_products, similarity_threshold)

        if match:
            matched_name, score = match
            row["status"] = "present" if score >= 0.95 else "partial"
            row["matched_to"] = matched_name
            row["match_similarity"] = f"{score:.2f}"
            matched_bgc_genes.add(matched_name)

            info = gene_info_map.get(matched_name, {})
            row["matched_locus_tag"] = info.get("locus_tag", "")
            row["matched_gene_kind"] = info.get("gene_kind", "")

        rows.append(row)

    # Summary row
    present = sum(1 for r in rows if r["status"] == "present")
    partial = sum(1 for r in rows if r["status"] == "partial")
    missing = sum(1 for r in rows if r["status"] == "missing")
    total = len(reference_genes)

    for row in rows:
        row["completeness"] = f"{present}/{total} ({present / total * 100:.0f}%)" if total else "0/0"
        row["present_count"] = str(present)
        row["partial_count"] = str(partial)
        row["missing_count"] = str(missing)

    return rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

_FIELDNAMES = [
    "file",
    "record_id",
    "region",
    "bgc_product",
    "contig_edge",
    "completeness",
    "present_count",
    "partial_count",
    "missing_count",
    "ref_gene",
    "ref_product",
    "ref_role",
    "status",
    "matched_to",
    "match_similarity",
    "matched_locus_tag",
    "matched_gene_kind",
]


def main(
    antismash_dir: Path,
    reference_genes_path: Path,
    similarity: float = 0.8,
    threads: int | None = None,
    output: Path | None = None,
) -> None:
    """Check BGC completeness across antiSMASH results.

    Args:
        antismash_dir: Directory containing antiSMASH result directories.
        reference_genes_path: Path to reference gene list TSV.
        similarity: Minimum similarity threshold for fuzzy matching.
        threads: Number of threads.
        output: Output TSV path, or None for stdout.
    """
    reference_genes = load_reference_genes(reference_genes_path)
    print(f"Loaded {len(reference_genes)} reference gene(s)", file=sys.stderr)

    json_files = list(antismash_dir.glob("*/*.json"))
    if not json_files:
        json_files = list(antismash_dir.glob("*.json"))

    if not json_files:
        print(f"Error: No antiSMASH JSON files found in {antismash_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Processing {len(json_files)} antiSMASH result(s)...", file=sys.stderr)

    all_rows: list[dict[str, str]] = []

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(parse_antismash_json, p): p for p in json_files}

        for future in as_completed(futures):
            path = futures[future]
            try:
                regions = future.result()
                for region in regions:
                    rows = assess_completeness(region, reference_genes, similarity)
                    all_rows.extend(rows)
                if regions:
                    print(f"  {path.name}: {len(regions)} region(s)", file=sys.stderr)
            except Exception as e:
                print(f"  Error processing {path.name}: {e}", file=sys.stderr)

    # Sort by file -> region -> status
    status_order = {"missing": 2, "partial": 1, "present": 0}
    all_rows.sort(key=lambda r: (r["file"], r["record_id"], r["region"], status_order.get(r["status"], 3)))

    # Output
    out_handle = output.open("w", newline="") if output else sys.stdout
    writer = csv.DictWriter(out_handle, fieldnames=_FIELDNAMES, delimiter="\t", extrasaction="ignore")
    writer.writeheader()
    writer.writerows(all_rows)

    if output:
        out_handle.close()
        print(f"\nResults written to: {output}", file=sys.stderr)

    # Summary
    regions_checked = len({(r["file"], r["record_id"], r["region"]) for r in all_rows})
    print(f"Regions assessed: {regions_checked}", file=sys.stderr)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Check BGC completeness by comparing each detected cluster's gene "
            "content against a user-supplied reference gene list (e.g. from "
            "MIBiG). Reports which core biosynthetic genes are present, "
            "partial, or missing per cluster."
        ),
    )

    parser.add_argument("antismash_dir", type=Path, help="Directory containing antiSMASH result directories")
    parser.add_argument(
        "reference_genes",
        type=Path,
        help="TSV file with reference gene list (columns: gene, product, role)",
    )
    parser.add_argument(
        "--similarity",
        type=float,
        default=0.8,
        help="Minimum name similarity (0-1) for fuzzy gene matching (default: 0.8)",
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
        antismash_dir=args.antismash_dir,
        reference_genes_path=args.reference_genes,
        similarity=args.similarity,
        threads=args.threads,
        output=args.output,
    )
