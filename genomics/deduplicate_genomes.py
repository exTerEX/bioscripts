#!/usr/bin/env python3
"""Remove duplicate genomes from a dehydrated NCBI datasets zip file.

Takes the duplicate report produced by detect_duplicate_strains.py and an NCBI
datasets zip file as input. For each cluster of duplicates the genome with the
fewest contigs is kept and the rest are removed from the archive metadata so
they will not be downloaded when the zip is rehydrated.

Modified files inside the zip:
    - ncbi_dataset/data/dataset_catalog.json  (assembly entries removed)
    - ncbi_dataset/data/assembly_data_report.jsonl  (lines removed)
    - ncbi_dataset/fetch.txt  (download lines removed)
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import json
import sys
import zipfile
from pathlib import Path

# ---------------------------------------------------------------------------
# Union-Find for clustering pairwise duplicates
# ---------------------------------------------------------------------------


class UnionFind:
    """Disjoint-set / union-find data structure."""

    def __init__(self) -> None:
        self.parent: dict[str, str] = {}
        self.rank: dict[str, int] = {}

    def find(self, x: str) -> str:
        if x not in self.parent:
            self.parent[x] = x
            self.rank[x] = 0
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: str, b: str) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1

    def get_clusters(self) -> dict[str, list[str]]:
        clusters: dict[str, list[str]] = {}
        for item in self.parent:
            root = self.find(item)
            clusters.setdefault(root, []).append(item)
        return clusters


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def source_id_to_accession(source_id: str) -> str:
    """Convert a source_id from detect_duplicate_strains to an accession.

    Source IDs from GBFF files look like ``GCF_008633005.1.gbff`` while those
    originating from accession lists are bare accessions such as
    ``GCF_008633005.1``.
    """
    # Strip common genomic file extensions
    for ext in (".gbff", ".gbk", ".gb", ".fna", ".fasta", ".fa"):
        if source_id.endswith(ext):
            return source_id[: -len(ext)]
    return source_id


def parse_duplicate_report(
    report_path: Path,
    confidence_filter: str | None = None,
) -> tuple[dict[str, list[str]], dict[str, int]]:
    """Parse the pairwise duplicate report and return clusters.

    Args:
        report_path: Path to duplicate_report.tsv.
        confidence_filter: If set, only consider matches at or above this
            confidence level (``high``, ``medium``, ``low``).

    Returns:
        Tuple of (clusters dict, contig_counts dict).
        Clusters maps a root accession to a list of member accessions.
        contig_counts maps accession -> number of contigs.
    """
    confidence_levels = {"high": 0, "medium": 1, "low": 2}
    max_level = confidence_levels.get(confidence_filter, 2) if confidence_filter else 2

    uf = UnionFind()
    contig_counts: dict[str, int] = {}

    with report_path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            level = confidence_levels.get(row["confidence"], 3)
            if level > max_level:
                continue

            acc1 = source_id_to_accession(row["strain_1"])
            acc2 = source_id_to_accession(row["strain_2"])

            uf.union(acc1, acc2)

            # Store contig counts when available
            with contextlib.suppress(ValueError, KeyError):
                contig_counts[acc1] = int(row["num_contigs_1"])
            with contextlib.suppress(ValueError, KeyError):
                contig_counts[acc2] = int(row["num_contigs_2"])

    # Only keep multi-member clusters
    clusters = {k: v for k, v in uf.get_clusters().items() if len(v) > 1}
    return clusters, contig_counts


def choose_representative(
    members: list[str],
    contig_counts_report: dict[str, int],
    contig_counts_jsonl: dict[str, int],
) -> str:
    """Pick the best representative genome (fewest contigs).

    Uses contig counts from the JSONL assembly report first, falling back to
    values from the duplicate report.

    Args:
        members: List of accessions in the cluster.
        contig_counts_report: Contig counts from duplicate_report.tsv.
        contig_counts_jsonl: Contig counts from assembly_data_report.jsonl.

    Returns:
        Accession of the representative genome.
    """

    def _contigs(acc: str) -> int:
        if acc in contig_counts_jsonl:
            return contig_counts_jsonl[acc]
        if acc in contig_counts_report:
            return contig_counts_report[acc]
        return sys.maxsize  # Unknown – never prefer

    return min(members, key=_contigs)


def parse_jsonl_metadata(content: str) -> dict[str, dict]:
    """Parse assembly_data_report.jsonl into per-accession metadata dicts.

    Returns:
        Mapping of accession to a dict with keys: num_contigs, organism,
        strain, total_length, gc_percent, biosample, assembly_level.
    """
    metadata: dict[str, dict] = {}
    for line in content.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            record = json.loads(line)
        except json.JSONDecodeError:
            continue
        accession = record.get("accession", "")
        if not accession:
            continue

        stats = record.get("assemblyStats", {})
        organism_info = record.get("organism", {})
        assembly_info = record.get("assemblyInfo", {})
        biosample_info = assembly_info.get("biosample", {})
        infraspecific = organism_info.get("infraspecificNames", {})

        metadata[accession] = {
            "num_contigs": _safe_int(stats.get("numberOfContigs")),
            "total_length": _safe_int(stats.get("totalSequenceLength")),
            "gc_percent": stats.get("gcPercent", ""),
            "contig_n50": _safe_int(stats.get("contigN50")),
            "organism": organism_info.get("organismName", ""),
            "strain": infraspecific.get("strain", ""),
            "biosample": biosample_info.get("accession", ""),
            "assembly_level": assembly_info.get("assemblyLevel", ""),
        }
    return metadata


def _safe_int(value) -> int | str:
    """Convert a value to int, returning empty string on failure."""
    try:
        return int(value)
    except (TypeError, ValueError):
        return ""


def parse_contig_counts_from_jsonl(content: str) -> dict[str, int]:
    """Extract accession -> numberOfContigs from assembly_data_report.jsonl."""
    counts: dict[str, int] = {}
    for line in content.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            record = json.loads(line)
        except json.JSONDecodeError:
            continue
        accession = record.get("accession", "")
        try:
            num_contigs = int(record.get("assemblyStats", {}).get("numberOfContigs", 0))
            counts[accession] = num_contigs
        except (TypeError, ValueError):
            pass
    return counts


def write_deduplication_report(
    removals_log: list[tuple[str, list[str], list[str]]],
    jsonl_metadata: dict[str, dict],
    contig_counts_report: dict[str, int],
    report_output: Path,
) -> None:
    """Write a TSV report detailing the deduplication results.

    Args:
        removals_log: List of (representative, removed_accessions, all_members).
        jsonl_metadata: Per-accession metadata parsed from JSONL.
        contig_counts_report: Contig counts from the duplicate report.
        report_output: Path for the output TSV.
    """
    with report_output.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(
            [
                "cluster_id",
                "cluster_size",
                "accession",
                "action",
                "organism",
                "strain",
                "num_contigs",
                "total_length",
                "gc_percent",
                "contig_n50",
                "biosample",
                "assembly_level",
            ]
        )

        for cluster_idx, (kept, _, members) in enumerate(sorted(removals_log, key=lambda r: -len(r[2])), start=1):
            for acc in sorted(members):
                meta = jsonl_metadata.get(acc, {})
                action = "keep" if acc == kept else "remove"

                num_contigs = meta.get("num_contigs", "")
                if num_contigs == "" and acc in contig_counts_report:
                    num_contigs = contig_counts_report[acc]

                writer.writerow(
                    [
                        cluster_idx,
                        len(members),
                        acc,
                        action,
                        meta.get("organism", ""),
                        meta.get("strain", ""),
                        num_contigs,
                        meta.get("total_length", ""),
                        meta.get("gc_percent", ""),
                        meta.get("contig_n50", ""),
                        meta.get("biosample", ""),
                        meta.get("assembly_level", ""),
                    ]
                )

    print(f"Deduplication report written to: {report_output}")


# ---------------------------------------------------------------------------
# Zip modification
# ---------------------------------------------------------------------------


def deduplicate_zip(
    zip_path: Path,
    report_path: Path,
    output_path: Path,
    report_output: Path | None = None,
    confidence_filter: str | None = None,
    dry_run: bool = False,
) -> None:
    """Create a new zip with duplicate genome entries removed.

    Args:
        zip_path: Path to the NCBI datasets zip file.
        report_path: Path to duplicate_report.tsv.
        output_path: Path for the deduplicated zip file.
        report_output: Path for the deduplication report TSV.
        confidence_filter: Minimum confidence level to consider.
        dry_run: If True, only report what would be removed.
    """
    # ------------------------------------------------------------------
    # 1. Parse duplicate report → clusters
    # ------------------------------------------------------------------
    clusters, contig_counts_report = parse_duplicate_report(report_path, confidence_filter)

    if not clusters:
        print("No duplicate clusters found – nothing to do.")
        return

    # ------------------------------------------------------------------
    # 2. Read assembly_data_report.jsonl from the zip for contig counts
    # ------------------------------------------------------------------
    jsonl_path = "ncbi_dataset/data/assembly_data_report.jsonl"

    with zipfile.ZipFile(zip_path, "r") as zf:
        try:
            jsonl_content = zf.read(jsonl_path).decode("utf-8")
        except KeyError:
            jsonl_content = ""
            print(
                f"Warning: {jsonl_path} not found in zip, contig counts from JSONL unavailable.",
                file=sys.stderr,
            )

    contig_counts_jsonl = parse_contig_counts_from_jsonl(jsonl_content)
    jsonl_metadata = parse_jsonl_metadata(jsonl_content)

    # ------------------------------------------------------------------
    # 3. Decide which accessions to remove
    # ------------------------------------------------------------------
    accessions_to_remove: set[str] = set()
    removals_log: list[tuple[str, list[str], list[str]]] = []  # (kept, [removed], [all])

    for _root, members in sorted(clusters.items()):
        representative = choose_representative(members, contig_counts_report, contig_counts_jsonl)
        to_remove = [m for m in members if m != representative]
        accessions_to_remove.update(to_remove)
        removals_log.append((representative, to_remove, members))

    # Summary
    print(f"Duplicate clusters:   {len(clusters)}")
    print(f"Genomes to remove:    {len(accessions_to_remove)}")
    print("Genomes remaining:    kept as representatives")
    print()

    for kept, removed_list, members in removals_log:
        contigs_kept = contig_counts_jsonl.get(kept, contig_counts_report.get(kept, "?"))
        print(f"  Cluster ({len(members)} members):")
        print(f"    Keep:   {kept}  (contigs: {contigs_kept})")
        for acc in removed_list:
            contigs_rm = contig_counts_jsonl.get(acc, contig_counts_report.get(acc, "?"))
            print(f"    Remove: {acc}  (contigs: {contigs_rm})")

    # Write deduplication report
    if report_output is not None:
        write_deduplication_report(removals_log, jsonl_metadata, contig_counts_report, report_output)

    if dry_run:
        print("\nDry-run mode – no files modified.")
        return

    # ------------------------------------------------------------------
    # 4. Write new zip excluding duplicates
    # ------------------------------------------------------------------
    catalog_path = "ncbi_dataset/data/dataset_catalog.json"
    fetch_path = "ncbi_dataset/fetch.txt"

    with (
        zipfile.ZipFile(zip_path, "r") as zf_in,
        zipfile.ZipFile(output_path, "w", compression=zipfile.ZIP_DEFLATED) as zf_out,
    ):
        for item in zf_in.infolist():
            raw = zf_in.read(item.filename)

            # --- dataset_catalog.json ---
            if item.filename == catalog_path:
                catalog = json.loads(raw.decode("utf-8"))
                original_count = len(catalog.get("assemblies", []))
                catalog["assemblies"] = [
                    entry
                    for entry in catalog.get("assemblies", [])
                    if entry.get("accession") not in accessions_to_remove
                ]
                filtered_count = len(catalog["assemblies"])
                raw = json.dumps(catalog, indent=2).encode("utf-8")
                print(f"\n{catalog_path}: {original_count} → {filtered_count} assemblies")

            # --- assembly_data_report.jsonl ---
            elif item.filename == jsonl_path:
                lines = raw.decode("utf-8").splitlines()
                kept_lines: list[str] = []
                removed_count = 0
                for line in lines:
                    stripped = line.strip()
                    if not stripped:
                        continue
                    try:
                        record = json.loads(stripped)
                    except json.JSONDecodeError:
                        kept_lines.append(line)
                        continue
                    if record.get("accession") in accessions_to_remove:
                        removed_count += 1
                    else:
                        kept_lines.append(line)
                raw = ("\n".join(kept_lines) + "\n").encode("utf-8")
                print(f"{jsonl_path}: removed {removed_count} entries")

            # --- fetch.txt ---
            elif item.filename == fetch_path:
                lines = raw.decode("utf-8").splitlines()
                kept_lines = []
                removed_count = 0
                for line in lines:
                    skip = False
                    for acc in accessions_to_remove:
                        # Lines contain paths like data/GCF_xxx/file
                        if f"data/{acc}/" in line or f"/{acc}/" in line:
                            skip = True
                            break
                    if skip:
                        removed_count += 1
                    else:
                        kept_lines.append(line)
                raw = ("\n".join(kept_lines) + "\n").encode("utf-8")
                print(f"{fetch_path}: removed {removed_count} lines")

            # --- Per-accession data directories ---
            else:
                # Skip files that belong to removed accessions
                skip = False
                for acc in accessions_to_remove:
                    if f"data/{acc}/" in item.filename or f"/{acc}/" in item.filename:
                        skip = True
                        break
                if skip:
                    continue

            zf_out.writestr(item, raw)

    print(f"\nDeduplicated zip written to: {output_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Remove duplicate genomes from an NCBI datasets zip file. "
            "Uses the duplicate report produced by detect_duplicate_strains.py "
            "to identify clusters and keeps the genome with the fewest contigs."
        )
    )

    parser.add_argument(
        "zip",
        type=Path,
        help="Path to the NCBI datasets zip file (dehydrated or hydrated)",
    )
    parser.add_argument(
        "report",
        type=Path,
        help="Path to duplicate_report.tsv produced by detect_duplicate_strains.py",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help=("Output path for the deduplicated zip file (default: <input>_deduplicated.zip)"),
    )
    parser.add_argument(
        "--confidence",
        choices=["high", "medium", "low"],
        default=None,
        help=(
            "Minimum confidence level to consider from the duplicate report. "
            "'high' only uses high-confidence matches, 'medium' uses high and "
            "medium, 'low' uses all. Default: use all matches."
        ),
    )
    parser.add_argument(
        "--report-output",
        type=Path,
        default=None,
        help=("Output path for the deduplication report TSV (default: <output_stem>_deduplication_report.tsv)"),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be removed without modifying anything",
    )

    args = parser.parse_args()

    # Validate inputs
    if not args.zip.exists():
        print(f"Error: zip file not found: {args.zip}", file=sys.stderr)
        sys.exit(1)
    if not args.report.exists():
        print(f"Error: duplicate report not found: {args.report}", file=sys.stderr)
        sys.exit(1)

    # Default output path
    output = args.output
    if output is None:
        output = args.zip.with_stem(args.zip.stem + "_deduplicated")

    # Avoid overwriting the input
    if output.resolve() == args.zip.resolve():
        print("Error: output path must differ from input zip.", file=sys.stderr)
        sys.exit(1)

    # Default report output path
    if args.report_output is not None:
        report_out = args.report_output
    else:
        report_out = output.with_name(output.stem + "_deduplication_report.tsv")

    deduplicate_zip(
        zip_path=args.zip,
        report_path=args.report,
        output_path=output,
        report_output=report_out,
        confidence_filter=args.confidence,
        dry_run=args.dry_run,
    )


if __name__ == "__main__":
    main()
