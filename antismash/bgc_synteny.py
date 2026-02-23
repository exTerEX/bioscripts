## Compare BGC gene arrangement across antiSMASH results for synteny analysis
#
## Usage:
#   $ python antismash/bgc_synteny.py -h
#   usage: bgc_synteny.py [-h] [--min-shared MIN_SHARED] [--similarity SIMILARITY]
#                         [--threads THREADS] [--output OUTPUT]
#                         antismash_dir
#
#   Compare the gene order and orientation of biosynthetic gene clusters
#   across antiSMASH results. Produces a pairwise synteny comparison matrix
#   showing gene conservation, inversions, and rearrangements.
#
#   positional arguments:
#     antismash_dir             Directory containing antiSMASH result directories
#
#   options:
#     -h, --help                show this help message and exit
#     --min-shared MIN_SHARED   Minimum shared genes to report a pair (default: 3)
#     --similarity SIMILARITY   Minimum similarity (0-1) for gene name matching (default: 0.8)
#     --threads THREADS         Number of threads for parallel processing
#     --output OUTPUT           Output TSV file path (default: stdout)

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from difflib import SequenceMatcher
from itertools import combinations
from pathlib import Path

# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------


@dataclass
class BGCGene:
    """A single gene within a BGC."""

    gene: str
    product: str
    locus_tag: str
    strand: int
    start: int
    end: int
    gene_kind: str


@dataclass
class BGCRegion:
    """A BGC region with ordered genes."""

    source_file: str
    record_id: str
    region_number: str
    product: str
    start: int
    end: int
    contig_edge: str
    genes: list[BGCGene] = field(default_factory=list)

    @property
    def label(self) -> str:
        """Short human-readable label for this region."""
        return f"{self.source_file}|{self.record_id}|r{self.region_number}"


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------


def _parse_strand(location: str) -> int:
    """Return +1 or -1 based on complement annotation."""
    return -1 if "complement" in location else 1


def _parse_location_bounds(location: str) -> tuple[int, int]:
    """Extract start and end positions from a feature location string."""
    nums = [int(n) for n in re.findall(r"\d+", location)]
    return (min(nums), max(nums)) if nums else (0, 0)


def parse_antismash_json(path: Path) -> list[BGCRegion]:
    """Parse an antiSMASH JSON file and extract BGC regions with ordered genes.

    Args:
        path: Path to the antiSMASH JSON file.

    Returns:
        List of BGCRegion objects with genes sorted by start position.
    """
    with path.open() as f:
        data = json.load(f)

    regions: list[BGCRegion] = []

    for record in data["records"]:
        region_features = [feat for feat in record["features"] if feat["type"] == "region"]

        for rfeat in region_features:
            quals = rfeat["qualifiers"]
            rstart, rend = _parse_location_bounds(rfeat["location"])

            region = BGCRegion(
                source_file=path.stem,
                record_id=record["name"],
                region_number=quals.get("region_number", ["?"])[0],
                product=" / ".join(quals.get("product", [])),
                start=rstart,
                end=rend,
                contig_edge=quals.get("contig_edge", [""])[0],
            )

            # Collect CDS features within this region
            for feat in record["features"]:
                if feat["type"] != "CDS":
                    continue
                fstart, fend = _parse_location_bounds(feat["location"])
                if fend < rstart or fstart > rend:
                    continue

                fq = feat.get("qualifiers", {})
                gene_name = fq.get("gene", [""])[0]
                product = fq.get("product", [""])[0]

                # Skip hypotheticals with no gene name
                identifier = gene_name or product
                if not identifier:
                    continue

                region.genes.append(
                    BGCGene(
                        gene=gene_name,
                        product=product,
                        locus_tag=fq.get("locus_tag", [""])[0],
                        strand=_parse_strand(feat["location"]),
                        start=fstart,
                        end=fend,
                        gene_kind=fq.get("gene_kind", [""])[0],
                    )
                )

            # Sort genes by genomic position
            region.genes.sort(key=lambda g: g.start)
            regions.append(region)

    return regions


# ---------------------------------------------------------------------------
# Gene matching
# ---------------------------------------------------------------------------


def _normalise(name: str) -> str:
    """Normalise a gene/product name for comparison."""
    return re.sub(r"[^a-z0-9]", "", name.lower())


def gene_similarity(a: BGCGene, b: BGCGene) -> float:
    """Compute a similarity score between two BGC genes.

    Compares gene name and product using fuzzy string matching.

    Args:
        a: First gene.
        b: Second gene.

    Returns:
        Similarity score between 0 and 1.
    """
    scores: list[float] = []

    # Compare gene names
    if a.gene and b.gene:
        na, nb = _normalise(a.gene), _normalise(b.gene)
        if na == nb:
            return 1.0
        scores.append(SequenceMatcher(None, na, nb).ratio())

    # Compare products
    if a.product and b.product:
        na, nb = _normalise(a.product), _normalise(b.product)
        if na == nb:
            return 1.0
        scores.append(SequenceMatcher(None, na, nb).ratio())

    # Cross-compare gene name with product
    if a.gene and b.product:
        scores.append(SequenceMatcher(None, _normalise(a.gene), _normalise(b.product)).ratio())
    if a.product and b.gene:
        scores.append(SequenceMatcher(None, _normalise(a.product), _normalise(b.gene)).ratio())

    return max(scores) if scores else 0.0


def build_gene_mapping(
    genes_a: list[BGCGene],
    genes_b: list[BGCGene],
    threshold: float,
) -> list[tuple[int, int, float]]:
    """Build a mapping of matching genes between two BGC regions.

    Uses greedy best-first matching.

    Args:
        genes_a: Genes from region A (positionally ordered).
        genes_b: Genes from region B (positionally ordered).
        threshold: Minimum similarity to accept a match.

    Returns:
        List of (index_a, index_b, similarity) tuples.
    """
    # Compute all pairwise similarities
    pairs: list[tuple[int, int, float]] = []
    for i, ga in enumerate(genes_a):
        for j, gb in enumerate(genes_b):
            sim = gene_similarity(ga, gb)
            if sim >= threshold:
                pairs.append((i, j, sim))

    # Greedy best-first assignment
    pairs.sort(key=lambda x: -x[2])
    used_a: set[int] = set()
    used_b: set[int] = set()
    mapping: list[tuple[int, int, float]] = []

    for i, j, sim in pairs:
        if i not in used_a and j not in used_b:
            mapping.append((i, j, sim))
            used_a.add(i)
            used_b.add(j)

    mapping.sort(key=lambda x: x[0])
    return mapping


# ---------------------------------------------------------------------------
# Synteny metrics
# ---------------------------------------------------------------------------


def compute_synteny_metrics(
    region_a: BGCRegion,
    region_b: BGCRegion,
    mapping: list[tuple[int, int, float]],
) -> dict[str, str]:
    """Compute synteny statistics for a pair of BGC regions.

    Args:
        region_a: First BGC region.
        region_b: Second BGC region.
        mapping: Gene mapping between the two regions.

    Returns:
        Dictionary of synteny metrics.
    """
    shared = len(mapping)
    total_a = len(region_a.genes)
    total_b = len(region_b.genes)

    # Check order conservation
    collinear = 0
    inversions = 0
    rearrangements = 0

    for k in range(len(mapping) - 1):
        ia, ib, _ = mapping[k]
        ia_next, ib_next, _ = mapping[k + 1]

        # Gene order: are they in the same order?
        if ib_next > ib:
            # Same direction — check strand conservation
            ga = region_a.genes[ia]
            gb = region_b.genes[ib]
            if ga.strand == gb.strand:
                collinear += 1
            else:
                inversions += 1
        else:
            rearrangements += 1

    # Jaccard index for gene repertoire
    jaccard = shared / (total_a + total_b - shared) if (total_a + total_b - shared) > 0 else 0.0

    # Synteny score = fraction of shared gene pairs that are collinear
    synteny_score = collinear / (shared - 1) if shared > 1 else (1.0 if shared == 1 else 0.0)

    # Average match similarity
    avg_sim = sum(m[2] for m in mapping) / shared if shared else 0.0

    # Gene order string for visual comparison
    order_a = " ".join(
        f"{'>' if region_a.genes[m[0]].strand == 1 else '<'}{region_a.genes[m[0]].gene or region_a.genes[m[0]].product}"
        for m in mapping
    )
    order_b = " ".join(
        f"{'>' if region_b.genes[m[1]].strand == 1 else '<'}{region_b.genes[m[1]].gene or region_b.genes[m[1]].product}"
        for m in mapping
    )

    return {
        "region_a": region_a.label,
        "bgc_product_a": region_a.product,
        "genes_a": str(total_a),
        "contig_edge_a": region_a.contig_edge,
        "region_b": region_b.label,
        "bgc_product_b": region_b.product,
        "genes_b": str(total_b),
        "contig_edge_b": region_b.contig_edge,
        "shared_genes": str(shared),
        "jaccard": f"{jaccard:.3f}",
        "collinear_pairs": str(collinear),
        "inversions": str(inversions),
        "rearrangements": str(rearrangements),
        "synteny_score": f"{synteny_score:.3f}",
        "avg_match_similarity": f"{avg_sim:.3f}",
        "gene_order_a": order_a,
        "gene_order_b": order_b,
    }


# ---------------------------------------------------------------------------
# Pairwise comparison
# ---------------------------------------------------------------------------


def compare_region_pair(
    pair: tuple[BGCRegion, BGCRegion],
    similarity_threshold: float,
    min_shared: int,
) -> dict[str, str] | None:
    """Compare a pair of BGC regions for synteny.

    Args:
        pair: Tuple of two BGCRegion objects.
        similarity_threshold: Gene matching threshold.
        min_shared: Minimum shared genes to report.

    Returns:
        Synteny metrics dictionary or None if insufficient shared genes.
    """
    region_a, region_b = pair

    mapping = build_gene_mapping(region_a.genes, region_b.genes, similarity_threshold)

    if len(mapping) < min_shared:
        return None

    return compute_synteny_metrics(region_a, region_b, mapping)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

_FIELDNAMES = [
    "region_a",
    "bgc_product_a",
    "genes_a",
    "contig_edge_a",
    "region_b",
    "bgc_product_b",
    "genes_b",
    "contig_edge_b",
    "shared_genes",
    "jaccard",
    "collinear_pairs",
    "inversions",
    "rearrangements",
    "synteny_score",
    "avg_match_similarity",
    "gene_order_a",
    "gene_order_b",
]


def main(
    antismash_dir: Path,
    min_shared: int = 3,
    similarity: float = 0.8,
    threads: int | None = None,
    output: Path | None = None,
) -> None:
    """Compare BGC synteny across antiSMASH results.

    Args:
        antismash_dir: Directory containing antiSMASH result directories.
        min_shared: Minimum shared genes to report a comparison.
        similarity: Minimum similarity for gene name matching.
        threads: Number of threads.
        output: Output TSV path, or None for stdout.
    """
    json_files = list(antismash_dir.glob("*/*.json"))
    if not json_files:
        json_files = list(antismash_dir.glob("*.json"))

    if not json_files:
        print(f"Error: No antiSMASH JSON files found in {antismash_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Processing {len(json_files)} antiSMASH result(s)...", file=sys.stderr)

    # Parse all regions
    all_regions: list[BGCRegion] = []
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(parse_antismash_json, p): p for p in json_files}
        for future in as_completed(futures):
            path = futures[future]
            try:
                regions = future.result()
                all_regions.extend(regions)
                if regions:
                    print(f"  {path.name}: {len(regions)} region(s)", file=sys.stderr)
            except Exception as e:
                print(f"  Error processing {path.name}: {e}", file=sys.stderr)

    if len(all_regions) < 2:
        print("Error: Need at least 2 BGC regions for comparison", file=sys.stderr)
        sys.exit(1)

    print(
        f"\nComparing {len(all_regions)} regions ({len(all_regions) * (len(all_regions) - 1) // 2} pairs)...",
        file=sys.stderr,
    )

    # Pairwise comparison
    pairs = list(combinations(all_regions, 2))
    all_rows: list[dict[str, str]] = []

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(compare_region_pair, pair, similarity, min_shared): pair for pair in pairs}
        for future in as_completed(futures):
            try:
                result = future.result()
                if result is not None:
                    all_rows.append(result)
            except Exception as e:
                print(f"  Comparison error: {e}", file=sys.stderr)

    # Sort by synteny score descending
    all_rows.sort(key=lambda r: float(r["synteny_score"]), reverse=True)

    # Output
    out_handle = output.open("w", newline="") if output else sys.stdout
    writer = csv.DictWriter(out_handle, fieldnames=_FIELDNAMES, delimiter="\t", extrasaction="ignore")
    writer.writeheader()
    writer.writerows(all_rows)

    if output:
        out_handle.close()
        print(f"\nResults written to: {output}", file=sys.stderr)

    print(f"Syntenic pairs found: {len(all_rows)}", file=sys.stderr)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Compare the gene order and orientation of biosynthetic gene "
            "clusters across antiSMASH results. Produces a pairwise synteny "
            "comparison matrix showing gene conservation, inversions, and "
            "rearrangements."
        ),
    )

    parser.add_argument("antismash_dir", type=Path, help="Directory containing antiSMASH result directories")
    parser.add_argument(
        "--min-shared",
        type=int,
        default=3,
        help="Minimum shared genes to report a pair (default: 3)",
    )
    parser.add_argument(
        "--similarity",
        type=float,
        default=0.8,
        help="Minimum similarity (0-1) for gene name matching (default: 0.8)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Number of threads for parallel processing",
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
        min_shared=args.min_shared,
        similarity=args.similarity,
        threads=args.threads,
        output=args.output,
    )
