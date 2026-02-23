## Detect candidate genomic islands by GC deviation and structural features
#
## Usage:
#   $ python genomics/island_finder.py -h
#   usage: island_finder.py [-h] [--window WINDOW] [--step STEP]
#                           [--gc-threshold GC_THRESHOLD] [--min-length MIN_LENGTH]
#                           [--threads THREADS] [--output OUTPUT]
#                           gbff_dir
#
#   Scan GBFF files for candidate genomic islands using sliding-window GC
#   content deviation from the genome average. Regions with aberrant GC% are
#   inspected for flanking tRNA genes, integrases, direct repeats, and
#   mobile genetic element features. Outputs a TSV of candidate island
#   regions with gene content and structural evidence.
#
#   positional arguments:
#     gbff_dir              Directory containing GBFF/GBK genome files
#
#   options:
#     -h, --help            show this help message and exit
#     --window WINDOW       Sliding window size in bp (default: 10000)
#     --step STEP           Sliding window step in bp (default: 5000)
#     --gc-threshold GC_THRESHOLD
#                           GC deviation threshold in standard deviations
#                           (default: 2.0)
#     --min-length MIN_LENGTH
#                           Minimum island length in bp (default: 5000)
#     --threads THREADS     Number of threads for parallel processing
#     --output OUTPUT       Output TSV file path (default: stdout)

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------


@dataclass
class IslandCandidate:
    """A candidate genomic island region."""

    record_id: str
    organism: str
    start: int
    end: int
    length: int
    gc_content: float
    genome_gc: float
    gc_deviation: float
    gc_zscore: float
    flanking_trna_left: str
    flanking_trna_right: str
    has_integrase: bool
    has_transposase: bool
    n_cds: int
    gene_products: list[str] = field(default_factory=list)
    mge_evidence: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# GC content analysis
# ---------------------------------------------------------------------------


def gc_fraction(seq: str) -> float:
    """Compute GC fraction for a DNA sequence.

    Args:
        seq: DNA sequence string.

    Returns:
        GC fraction (0-1).
    """
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    total = len(seq) - seq.count("N")
    return gc / total if total > 0 else 0.0


def sliding_window_gc(
    sequence: str,
    window: int,
    step: int,
) -> list[tuple[int, int, float]]:
    """Calculate GC content over sliding windows.

    Args:
        sequence: DNA sequence.
        window: Window size in bp.
        step: Step size in bp.

    Returns:
        List of (start, end, gc_fraction) tuples.
    """
    results: list[tuple[int, int, float]] = []
    seq_len = len(sequence)

    for start in range(0, seq_len - window + 1, step):
        end = start + window
        gc = gc_fraction(sequence[start:end])
        results.append((start, end, gc))

    return results


def find_aberrant_regions(
    gc_windows: list[tuple[int, int, float]],
    genome_gc: float,
    gc_std: float,
    threshold: float,
    min_length: int,
) -> list[tuple[int, int, float]]:
    """Merge consecutive windows with aberrant GC into candidate regions.

    Args:
        gc_windows: Sliding window GC results.
        genome_gc: Genome-wide average GC.
        gc_std: Standard deviation of windowed GC values.
        threshold: Number of standard deviations for aberrant call.
        min_length: Minimum region length in bp.

    Returns:
        List of (start, end, mean_gc) for aberrant regions.
    """
    if gc_std == 0:
        return []

    aberrant: list[tuple[int, int, float]] = []

    # Identify aberrant windows
    for start, end, gc in gc_windows:
        zscore = abs(gc - genome_gc) / gc_std
        if zscore >= threshold:
            aberrant.append((start, end, gc))

    if not aberrant:
        return []

    # Merge consecutive/overlapping aberrant windows
    merged: list[tuple[int, int, list[float]]] = []
    curr_start, curr_end = aberrant[0][0], aberrant[0][1]
    curr_gcs: list[float] = [aberrant[0][2]]

    for start, end, gc in aberrant[1:]:
        if start <= curr_end:  # Overlapping or adjacent
            curr_end = max(curr_end, end)
            curr_gcs.append(gc)
        else:
            merged.append((curr_start, curr_end, curr_gcs))
            curr_start, curr_end = start, end
            curr_gcs = [gc]

    merged.append((curr_start, curr_end, curr_gcs))

    # Filter by minimum length
    regions: list[tuple[int, int, float]] = []
    for start, end, gcs in merged:
        if end - start >= min_length:
            mean_gc = sum(gcs) / len(gcs)
            regions.append((start, end, mean_gc))

    return regions


# ---------------------------------------------------------------------------
# Feature annotation within regions
# ---------------------------------------------------------------------------

_MGE_PATTERNS = re.compile(
    r"transpos|integrase|recombinas|insertion.?sequence|IS\d+|"
    r"phage|conjugat|mobiliz|resolvase|invertase",
    re.IGNORECASE,
)

_TRNA_PATTERN = re.compile(r"tRNA", re.IGNORECASE)


def annotate_island(
    record: SeqRecord,
    start: int,
    end: int,
    genome_gc: float,
    gc_std: float,
    region_gc: float,
) -> IslandCandidate:
    """Annotate a candidate island region with gene content and structural features.

    Args:
        record: BioPython SeqRecord.
        start: Island start coordinate.
        end: Island end coordinate.
        genome_gc: Genome-wide average GC.
        gc_std: Standard deviation of windowed GC.
        region_gc: Mean GC of the region.

    Returns:
        IslandCandidate object.
    """
    organism = record.annotations.get("organism", "")

    # Collect features within the region
    cds_products: list[str] = []
    has_integrase = False
    has_transposase = False
    mge_evidence: list[str] = []

    # Look for flanking tRNAs (within 2kb of boundaries)
    flanking_trna_left = ""
    flanking_trna_right = ""
    flank_window = 2000

    for feature in record.features:
        if feature.location is None:  # type: ignore[union-attr]
            continue

        feat_start = int(feature.location.start)  # type: ignore[union-attr]
        feat_end = int(feature.location.end)  # type: ignore[union-attr]

        # Check tRNAs flanking the island
        if feature.type == "tRNA":
            product = feature.qualifiers.get("product", [""])[0]
            if max(0, start - flank_window) <= feat_start <= start + flank_window:
                flanking_trna_left = product or "tRNA"
            if max(0, end - flank_window) <= feat_start <= end + flank_window:
                flanking_trna_right = product or "tRNA"

        # Features within the island
        if feat_end < start or feat_start > end:
            continue

        if feature.type == "CDS":
            quals = feature.qualifiers
            product = quals.get("product", ["hypothetical protein"])[0]
            gene = quals.get("gene", [""])[0]
            locus_tag = quals.get("locus_tag", [""])[0]

            label = gene or locus_tag
            cds_products.append(f"{label}: {product}" if label else product)

            # Check for MGE signatures
            searchable = f"{product} {gene}"
            if _MGE_PATTERNS.search(searchable):
                if "integrase" in searchable.lower() or "recombinas" in searchable.lower():
                    has_integrase = True
                    mge_evidence.append(f"integrase ({label})")
                if "transpos" in searchable.lower():
                    has_transposase = True
                    mge_evidence.append(f"transposase ({label})")
                if "phage" in searchable.lower():
                    mge_evidence.append(f"phage-related ({label})")
                if re.search(r"IS\d+", searchable):
                    mge_evidence.append(f"IS element ({label})")

    gc_deviation = region_gc - genome_gc
    gc_zscore = abs(gc_deviation) / gc_std if gc_std > 0 else 0.0

    return IslandCandidate(
        record_id=record.id,  # type: ignore
        organism=organism,  # type: ignore
        start=start,
        end=end,
        length=end - start,
        gc_content=region_gc,
        genome_gc=genome_gc,
        gc_deviation=gc_deviation,
        gc_zscore=gc_zscore,
        flanking_trna_left=flanking_trna_left,
        flanking_trna_right=flanking_trna_right,
        has_integrase=has_integrase,
        has_transposase=has_transposase,
        n_cds=len(cds_products),
        gene_products=cds_products,
        mge_evidence=mge_evidence,
    )


# ---------------------------------------------------------------------------
# Per-genome processing
# ---------------------------------------------------------------------------


def process_genome(
    path: Path,
    window: int,
    step: int,
    gc_threshold: float,
    min_length: int,
) -> list[dict[str, str]]:
    """Process a single GBFF genome file for candidate genomic islands.

    Args:
        path: Path to the GBFF/GBK file.
        window: Sliding window size.
        step: Sliding window step.
        gc_threshold: GC deviation threshold (SD).
        min_length: Minimum island length.

    Returns:
        List of result row dictionaries.
    """
    rows: list[dict[str, str]] = []

    for record in SeqIO.parse(path, "genbank"):
        seq_str = str(record.seq)
        seq_len = len(seq_str)

        if seq_len < window:
            continue

        # Genome-wide GC
        genome_gc = gc_fraction(seq_str)

        # Sliding window GC
        gc_windows = sliding_window_gc(seq_str, window, step)
        if not gc_windows:
            continue

        # Standard deviation of window GC values
        gc_values = [w[2] for w in gc_windows]
        mean_gc = sum(gc_values) / len(gc_values)
        variance = sum((g - mean_gc) ** 2 for g in gc_values) / len(gc_values)
        gc_std = math.sqrt(variance) if variance > 0 else 0.0

        # Find aberrant regions
        aberrant = find_aberrant_regions(gc_windows, genome_gc, gc_std, gc_threshold, min_length)

        for start, end, region_gc in aberrant:
            island = annotate_island(record, start, end, genome_gc, gc_std, region_gc)

            # Compute confidence score
            evidence_count = sum(
                [
                    bool(island.flanking_trna_left),
                    bool(island.flanking_trna_right),
                    island.has_integrase,
                    island.has_transposase,
                    island.gc_zscore >= 3.0,
                ]
            )

            rows.append(
                {
                    "file": path.name,
                    "record_id": island.record_id,
                    "organism": island.organism,
                    "start": str(island.start),
                    "end": str(island.end),
                    "length": str(island.length),
                    "gc_content": f"{island.gc_content:.4f}",
                    "genome_gc": f"{island.genome_gc:.4f}",
                    "gc_deviation": f"{island.gc_deviation:+.4f}",
                    "gc_zscore": f"{island.gc_zscore:.2f}",
                    "flanking_trna_left": island.flanking_trna_left,
                    "flanking_trna_right": island.flanking_trna_right,
                    "has_integrase": str(island.has_integrase),
                    "has_transposase": str(island.has_transposase),
                    "n_cds": str(island.n_cds),
                    "evidence_count": str(evidence_count),
                    "mge_evidence": "; ".join(island.mge_evidence),
                    "gene_content": " | ".join(island.gene_products[:30]),  # Truncate for readability
                }
            )

    return rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

_FIELDNAMES = [
    "file",
    "record_id",
    "organism",
    "start",
    "end",
    "length",
    "gc_content",
    "genome_gc",
    "gc_deviation",
    "gc_zscore",
    "flanking_trna_left",
    "flanking_trna_right",
    "has_integrase",
    "has_transposase",
    "n_cds",
    "evidence_count",
    "mge_evidence",
    "gene_content",
]


def main(
    gbff_dir: Path,
    window: int = 10000,
    step: int = 5000,
    gc_threshold: float = 2.0,
    min_length: int = 5000,
    threads: int | None = None,
    output: Path | None = None,
) -> None:
    """Scan genomes for candidate genomic islands.

    Args:
        gbff_dir: Directory containing GBFF/GBK genome files.
        window: Sliding window size in bp.
        step: Sliding window step in bp.
        gc_threshold: GC deviation threshold in standard deviations.
        min_length: Minimum island length in bp.
        threads: Number of threads.
        output: Output TSV path, or None for stdout.
    """
    gbff_files = list(gbff_dir.glob("**/*.gbff")) + list(gbff_dir.glob("**/*.gbk")) + list(gbff_dir.glob("**/*.gb"))

    if not gbff_files:
        print(f"Error: No GBFF/GBK files found in {gbff_dir}", file=sys.stderr)
        sys.exit(1)

    print(
        f"Scanning {len(gbff_files)} genome(s) (window={window}, step={step}, "
        f"threshold={gc_threshold} SD, min_length={min_length} bp)...",
        file=sys.stderr,
    )

    all_rows: list[dict[str, str]] = []

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(process_genome, p, window, step, gc_threshold, min_length): p for p in gbff_files}

        for future in as_completed(futures):
            path = futures[future]
            try:
                rows = future.result()
                all_rows.extend(rows)
                if rows:
                    print(f"  {path.name}: {len(rows)} candidate island(s)", file=sys.stderr)
                else:
                    print(f"  {path.name}: no candidates", file=sys.stderr)
            except Exception as e:
                print(f"  Error processing {path.name}: {e}", file=sys.stderr)

    # Sort by evidence count (descending), then GC z-score
    all_rows.sort(key=lambda r: (int(r["evidence_count"]), float(r["gc_zscore"])), reverse=True)

    # Output
    out_handle = output.open("w", newline="") if output else sys.stdout
    writer = csv.DictWriter(out_handle, fieldnames=_FIELDNAMES, delimiter="\t", extrasaction="ignore")
    writer.writeheader()
    writer.writerows(all_rows)

    if output:
        out_handle.close()
        print(f"\nResults written to: {output}", file=sys.stderr)

    # Summary
    high_conf = sum(1 for r in all_rows if int(r["evidence_count"]) >= 3)
    print(f"Total candidates: {len(all_rows)} ({high_conf} with >=3 evidence features)", file=sys.stderr)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Scan GBFF files for candidate genomic islands using sliding-window "
            "GC content deviation from the genome average. Regions with aberrant "
            "GC%% are inspected for flanking tRNA genes, integrases, direct repeats, "
            "and mobile genetic element features."
        ),
    )

    parser.add_argument("gbff_dir", type=Path, help="Directory containing GBFF/GBK genome files")
    parser.add_argument("--window", type=int, default=10000, help="Sliding window size in bp (default: 10000)")
    parser.add_argument("--step", type=int, default=5000, help="Sliding window step in bp (default: 5000)")
    parser.add_argument(
        "--gc-threshold",
        type=float,
        default=2.0,
        help="GC deviation threshold in standard deviations (default: 2.0)",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=5000,
        help="Minimum island length in bp (default: 5000)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Number of threads for parallel processing",
    )
    parser.add_argument("--output", type=Path, default=None, help="Output TSV file path (default: stdout)")

    args = parser.parse_args()

    main(
        gbff_dir=args.gbff_dir,
        window=args.window,
        step=args.step,
        gc_threshold=args.gc_threshold,
        min_length=args.min_length,
        threads=args.threads,
        output=args.output,
    )
