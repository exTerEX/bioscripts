## Predict operons and compute per-operon co-expression statistics
#
## Usage:
#   $ python misc/operon_coexpression.py -h
#   usage: operon_coexpression.py [-h] [--max-gap MAX_GAP] [--min-genes MIN_GENES]
#                                 [--count-col COUNT_COL] [--threads THREADS]
#                                 [--output OUTPUT]
#                                 annotation counts
#
#   Predict operons from GBFF/GFF annotation (same-strand genes within a
#   maximum intergenic gap) and compute co-expression statistics per operon
#   from a read count or TPM matrix. Reports mean expression, coefficient
#   of variation, and Pearson correlation across samples for each operon.
#
#   positional arguments:
#     annotation            GBFF/GBK annotation file or GFF3 file
#     counts                Gene expression count matrix TSV (genes as rows,
#                           samples as columns)
#
#   options:
#     -h, --help            show this help message and exit
#     --max-gap MAX_GAP     Maximum intergenic distance (bp) for operon
#                           prediction (default: 50)
#     --min-genes MIN_GENES Minimum genes per operon to report (default: 2)
#     --count-col COUNT_COL Column name for gene identifiers in count matrix
#                           (default: locus_tag)
#     --threads THREADS     Number of threads for correlation computation
#     --output OUTPUT       Output TSV file path (default: stdout)

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path

from Bio import SeqIO

# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------


@dataclass
class GeneEntry:
    """A CDS feature with positional information."""

    locus_tag: str
    gene: str
    product: str
    start: int
    end: int
    strand: int  # +1 or -1
    record_id: str


@dataclass
class Operon:
    """A predicted operon (group of co-directional, closely-spaced genes)."""

    operon_id: str
    record_id: str
    strand: str
    genes: list[GeneEntry] = field(default_factory=list)

    @property
    def start(self) -> int:
        return min(g.start for g in self.genes)

    @property
    def end(self) -> int:
        return max(g.end for g in self.genes)

    @property
    def locus_tags(self) -> list[str]:
        return [g.locus_tag for g in self.genes]


# ---------------------------------------------------------------------------
# Annotation parsing
# ---------------------------------------------------------------------------


def parse_gbff_genes(path: Path) -> list[GeneEntry]:
    """Extract CDS features from GBFF/GBK file.

    Args:
        path: Path to the annotation file.

    Returns:
        List of GeneEntry objects sorted by (record_id, start).
    """
    genes: list[GeneEntry] = []

    for record in SeqIO.parse(path, "genbank"):
        for feature in record.features:
            if feature.type != "CDS":
                continue

            quals = feature.qualifiers
            locus_tag = quals.get("locus_tag", [""])[0]
            if not locus_tag:
                continue

            genes.append(
                GeneEntry(
                    locus_tag=locus_tag,
                    gene=quals.get("gene", [""])[0],
                    product=quals.get("product", [""])[0],
                    start=int(feature.location.start),  # type: ignore[union-attr]
                    end=int(feature.location.end),  # type: ignore[union-attr]
                    strand=feature.location.strand or 1,  # type: ignore[union-attr]
                    record_id=record.id,
                )
            )

    genes.sort(key=lambda g: (g.record_id, g.start))
    return genes


def parse_gff_genes(path: Path) -> list[GeneEntry]:
    """Extract CDS features from a GFF3 file.

    Args:
        path: Path to the GFF3 file.

    Returns:
        List of GeneEntry objects sorted by (record_id, start).
    """
    genes: list[GeneEntry] = []

    with path.open() as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9 or cols[2] != "CDS":
                continue

            attrs: dict[str, str] = {}
            for attr in cols[8].split(";"):
                if "=" in attr:
                    key, val = attr.split("=", 1)
                    attrs[key.strip()] = val.strip()

            locus_tag = attrs.get("locus_tag", attrs.get("ID", ""))
            if not locus_tag:
                continue

            strand = -1 if cols[6] == "-" else 1
            genes.append(
                GeneEntry(
                    locus_tag=locus_tag,
                    gene=attrs.get("gene", attrs.get("Name", "")),
                    product=attrs.get("product", ""),
                    start=int(cols[3]) - 1,  # GFF is 1-based
                    end=int(cols[4]),
                    strand=strand,
                    record_id=cols[0],
                )
            )

    genes.sort(key=lambda g: (g.record_id, g.start))
    return genes


# ---------------------------------------------------------------------------
# Operon prediction
# ---------------------------------------------------------------------------


def predict_operons(genes: list[GeneEntry], max_gap: int, min_genes: int) -> list[Operon]:
    """Predict operons based on intergenic distance and strand.

    Genes on the same strand with intergenic gaps <= *max_gap* bp are
    grouped into the same operon.

    Args:
        genes: Sorted list of gene entries.
        max_gap: Maximum intergenic distance for operon grouping.
        min_genes: Minimum number of genes per operon.

    Returns:
        List of predicted operons.
    """
    if not genes:
        return []

    operons: list[Operon] = []
    operon_counter = 0

    current_group: list[GeneEntry] = [genes[0]]

    for i in range(1, len(genes)):
        prev = genes[i - 1]
        curr = genes[i]

        # Same record, same strand, within gap threshold
        if curr.record_id == prev.record_id and curr.strand == prev.strand and curr.start - prev.end <= max_gap:
            current_group.append(curr)
        else:
            if len(current_group) >= min_genes:
                operon_counter += 1
                operons.append(
                    Operon(
                        operon_id=f"operon_{operon_counter:04d}",
                        record_id=current_group[0].record_id,
                        strand="+" if current_group[0].strand == 1 else "-",
                        genes=list(current_group),
                    )
                )
            current_group = [curr]

    # Don't forget the last group
    if len(current_group) >= min_genes:
        operon_counter += 1
        operons.append(
            Operon(
                operon_id=f"operon_{operon_counter:04d}",
                record_id=current_group[0].record_id,
                strand="+" if current_group[0].strand == 1 else "-",
                genes=list(current_group),
            )
        )

    return operons


# ---------------------------------------------------------------------------
# Expression statistics
# ---------------------------------------------------------------------------


def load_count_matrix(path: Path, id_col: str) -> tuple[list[str], dict[str, list[float]]]:
    """Load a gene expression count / TPM matrix.

    Args:
        path: Path to the count matrix TSV.
        id_col: Column name for gene identifiers.

    Returns:
        Tuple of (sample_names, gene_expression_dict) where
        gene_expression_dict maps locus_tag -> list of expression values.
    """
    with path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            print("Error: Empty count matrix", file=sys.stderr)
            sys.exit(1)

        if id_col not in reader.fieldnames:
            print(f"Error: Column '{id_col}' not found. Available: {', '.join(reader.fieldnames)}", file=sys.stderr)
            sys.exit(1)

        sample_names = [col for col in reader.fieldnames if col != id_col]
        expression: dict[str, list[float]] = {}

        for row in reader:
            gene_id = row[id_col].strip()
            if not gene_id:
                continue

            values: list[float] = []
            for sample in sample_names:
                try:
                    values.append(float(row[sample]))
                except (ValueError, KeyError):
                    values.append(0.0)

            expression[gene_id] = values

    return sample_names, expression


def _pearson_correlation(x: list[float], y: list[float]) -> float:
    """Compute Pearson correlation coefficient between two vectors.

    Args:
        x: First vector.
        y: Second vector.

    Returns:
        Pearson correlation coefficient, or 0.0 if undefined.
    """
    n = len(x)
    if n < 3:
        return 0.0

    mean_x = sum(x) / n
    mean_y = sum(y) / n

    num = sum((xi - mean_x) * (yi - mean_y) for xi, yi in zip(x, y, strict=False))
    den_x = math.sqrt(sum((xi - mean_x) ** 2 for xi in x))
    den_y = math.sqrt(sum((yi - mean_y) ** 2 for yi in y))

    if den_x == 0 or den_y == 0:
        return 0.0

    return num / (den_x * den_y)


def compute_operon_stats(
    operon: Operon,
    expression: dict[str, list[float]],
    sample_names: list[str],
) -> dict[str, str]:
    """Compute expression statistics for an operon.

    Args:
        operon: The operon to analyse.
        expression: Gene expression dictionary.
        sample_names: List of sample names.

    Returns:
        Dictionary with operon-level statistics.
    """
    locus_tags = operon.locus_tags
    gene_names = [g.gene or g.locus_tag for g in operon.genes]
    products = [g.product for g in operon.genes]

    # Collect expression profiles for operon genes
    profiles: list[list[float]] = []
    matched_tags: list[str] = []

    for lt in locus_tags:
        if lt in expression:
            profiles.append(expression[lt])
            matched_tags.append(lt)

    n_matched = len(profiles)

    # Per-gene mean expression
    gene_means: list[float] = []
    for profile in profiles:
        gene_means.append(sum(profile) / len(profile) if profile else 0.0)

    # Operon-level mean expression (average of gene means)
    operon_mean = sum(gene_means) / len(gene_means) if gene_means else 0.0

    # Coefficient of variation across gene means within the operon
    if len(gene_means) > 1 and operon_mean > 0:
        variance = sum((gm - operon_mean) ** 2 for gm in gene_means) / len(gene_means)
        cv = math.sqrt(variance) / operon_mean
    else:
        cv = 0.0

    # Mean pairwise Pearson correlation
    correlations: list[float] = []
    for i in range(len(profiles)):
        for j in range(i + 1, len(profiles)):
            r = _pearson_correlation(profiles[i], profiles[j])
            correlations.append(r)

    mean_corr = sum(correlations) / len(correlations) if correlations else 0.0
    min_corr = min(correlations) if correlations else 0.0

    return {
        "operon_id": operon.operon_id,
        "record_id": operon.record_id,
        "strand": operon.strand,
        "start": str(operon.start),
        "end": str(operon.end),
        "n_genes": str(len(operon.genes)),
        "n_matched": str(n_matched),
        "locus_tags": ", ".join(locus_tags),
        "gene_names": ", ".join(gene_names),
        "products": ", ".join(products),
        "mean_expression": f"{operon_mean:.2f}",
        "cv_within_operon": f"{cv:.3f}",
        "mean_pairwise_correlation": f"{mean_corr:.3f}",
        "min_pairwise_correlation": f"{min_corr:.3f}",
        "n_pairwise_comparisons": str(len(correlations)),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

_FIELDNAMES = [
    "operon_id",
    "record_id",
    "strand",
    "start",
    "end",
    "n_genes",
    "n_matched",
    "locus_tags",
    "gene_names",
    "products",
    "mean_expression",
    "cv_within_operon",
    "mean_pairwise_correlation",
    "min_pairwise_correlation",
    "n_pairwise_comparisons",
]


def main(
    annotation: Path,
    counts: Path,
    max_gap: int = 50,
    min_genes: int = 2,
    count_col: str = "locus_tag",
    threads: int | None = None,
    output: Path | None = None,
) -> None:
    """Predict operons and compute co-expression statistics.

    Args:
        annotation: Path to GBFF/GBK or GFF3 annotation file.
        counts: Path to gene expression count matrix TSV.
        max_gap: Maximum intergenic distance for operon prediction.
        min_genes: Minimum genes per operon.
        count_col: Column name for gene identifiers in count matrix.
        threads: Number of threads (reserved for future use).
        output: Output TSV path, or None for stdout.
    """
    # Detect annotation format
    suffix = annotation.suffix.lower()
    if suffix in {".gbff", ".gbk", ".gb"}:
        print(f"Parsing GBFF annotation: {annotation.name}", file=sys.stderr)
        genes = parse_gbff_genes(annotation)
    elif suffix in {".gff", ".gff3"}:
        print(f"Parsing GFF annotation: {annotation.name}", file=sys.stderr)
        genes = parse_gff_genes(annotation)
    else:
        print(f"Error: Unrecognised annotation format '{suffix}'. Use .gbff/.gbk/.gb or .gff/.gff3", file=sys.stderr)
        sys.exit(1)

    print(f"  {len(genes)} CDS features extracted", file=sys.stderr)

    # Predict operons
    operons = predict_operons(genes, max_gap, min_genes)
    print(f"  {len(operons)} operons predicted (max_gap={max_gap}, min_genes={min_genes})", file=sys.stderr)

    # Load expression data
    print(f"Loading count matrix: {counts.name}", file=sys.stderr)
    sample_names, expression = load_count_matrix(counts, count_col)
    print(f"  {len(expression)} genes × {len(sample_names)} samples", file=sys.stderr)

    # Compute statistics per operon
    rows: list[dict[str, str]] = []
    for operon in operons:
        row = compute_operon_stats(operon, expression, sample_names)
        rows.append(row)

    # Sort by mean correlation descending
    rows.sort(key=lambda r: float(r["mean_pairwise_correlation"]), reverse=True)

    # Output
    out_handle = output.open("w", newline="") if output else sys.stdout
    writer = csv.DictWriter(out_handle, fieldnames=_FIELDNAMES, delimiter="\t", extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)

    if output:
        out_handle.close()
        print(f"\nResults written to: {output}", file=sys.stderr)

    # Summary
    high_corr = sum(1 for r in rows if float(r["mean_pairwise_correlation"]) >= 0.8)
    print(f"Operons with mean correlation >= 0.8: {high_corr}/{len(rows)}", file=sys.stderr)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Predict operons from GBFF/GFF annotation (same-strand genes within "
            "a maximum intergenic gap) and compute co-expression statistics per "
            "operon from a read count or TPM matrix."
        ),
    )

    parser.add_argument("annotation", type=Path, help="GBFF/GBK annotation file or GFF3 file")
    parser.add_argument(
        "counts",
        type=Path,
        help="Gene expression count matrix TSV (genes as rows, samples as columns)",
    )
    parser.add_argument(
        "--max-gap",
        type=int,
        default=50,
        help="Maximum intergenic distance (bp) for operon prediction (default: 50)",
    )
    parser.add_argument("--min-genes", type=int, default=2, help="Minimum genes per operon to report (default: 2)")
    parser.add_argument(
        "--count-col",
        default="locus_tag",
        help="Column name for gene identifiers in count matrix (default: locus_tag)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Number of threads for correlation computation",
    )
    parser.add_argument("--output", type=Path, default=None, help="Output TSV file path (default: stdout)")

    args = parser.parse_args()

    main(
        annotation=args.annotation,
        counts=args.counts,
        max_gap=args.max_gap,
        min_genes=args.min_genes,
        count_col=args.count_col,
        threads=args.threads,
        output=args.output,
    )
