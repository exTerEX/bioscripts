## Annotate differential expression results with genomic feature information
#
## Usage:
#   $ python misc/expression_annotate.py -h
#   usage: expression_annotate.py [-h] [--id-col ID_COL] [--gbff-dir GBFF_DIR]
#                                 [--gff GFF] [--keywords KEYWORDS]
#                                 [--threads THREADS] [--output OUTPUT]
#                                 de_results
#
#   Annotate DESeq2/edgeR differential expression results with gene product
#   descriptions, COG categories, and biosynthetic/AMR/QS flags from GBFF
#   or GFF annotation files. Outputs an enriched TSV.
#
#   positional arguments:
#     de_results            DESeq2/edgeR results TSV (must contain a locus_tag
#                           or gene ID column)
#
#   options:
#     -h, --help            show this help message and exit
#     --id-col ID_COL       Column name for gene identifiers (default: locus_tag)
#     --gbff-dir GBFF_DIR   Directory containing GBFF/GBK annotation files
#     --gff GFF             GFF3 annotation file (alternative to --gbff-dir)
#     --keywords KEYWORDS   Comma-separated keywords for flagging genes of
#                           interest (default: built-in AMR/QS/BGC list)
#     --threads THREADS     Number of threads for parallel GBFF parsing
#     --output OUTPUT       Output TSV file path (default: stdout)

from __future__ import annotations

import argparse
import csv
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from Bio import SeqIO

# ---------------------------------------------------------------------------
# Default keyword sets for flagging
# ---------------------------------------------------------------------------

_DEFAULT_AMR_KEYWORDS = [
    "resistance",
    "beta-lactamase",
    "carbapenemase",
    "efflux",
    "aminoglycoside",
    "tetracycline",
    "chloramphenicol",
    "vancomycin",
    "methicillin",
    "oxacillin",
    "penicillin",
    "cephalosporin",
    "fluoroquinolone",
    "macrolide",
    "sulfonamide",
    "trimethoprim",
    "multidrug",
    "drug resistance",
    "antibiotic",
]

_DEFAULT_QS_KEYWORDS = [
    "quorum sensing",
    "autoinducer",
    "acyl-homoserine lactone",
    "homoserine lactone",
    "AHL",
    "luxI",
    "luxR",
    "lasI",
    "lasR",
    "rhlI",
    "rhlR",
    "AI-2",
    "lux",
    "signaling",
    "signal",
    "N-acyl",
    "synthase",
]

_DEFAULT_BGC_KEYWORDS = [
    "biosynthetic",
    "NRPS",
    "PKS",
    "polyketide",
    "nonribosomal",
    "terpene",
    "siderophore",
    "bacteriocin",
    "lanthipeptide",
    "thiopeptide",
    "RiPP",
    "arylpolyene",
    "homoserine",
]

_DEFAULT_MGE_KEYWORDS = [
    "transposase",
    "integrase",
    "insertion sequence",
    "IS element",
    "recombinase",
    "conjugal",
    "mobilization",
]


# ---------------------------------------------------------------------------
# Annotation extraction from GBFF
# ---------------------------------------------------------------------------


def parse_gbff_annotations(path: Path) -> dict[str, dict[str, str]]:
    """Parse a GBFF file and extract annotation per locus tag.

    Args:
        path: Path to the GBFF/GBK file.

    Returns:
        Dictionary mapping locus_tag -> annotation dict with gene, product,
        protein_id, ec_number, db_xref, note fields.
    """
    annotations: dict[str, dict[str, str]] = {}

    for record in SeqIO.parse(path, "genbank"):
        for feature in record.features:
            if feature.type != "CDS":
                continue

            quals = feature.qualifiers
            locus_tag = quals.get("locus_tag", [""])[0]
            if not locus_tag:
                continue

            # Extract location info
            strand = "+" if feature.location.strand == 1 else "-"  # type: ignore[union-attr]
            start = int(feature.location.start)  # type: ignore[union-attr]
            end = int(feature.location.end)  # type: ignore[union-attr]

            annotations[locus_tag] = {
                "gene": quals.get("gene", [""])[0],
                "product": quals.get("product", [""])[0],
                "protein_id": quals.get("protein_id", [""])[0],
                "ec_number": "; ".join(quals.get("EC_number", [])),
                "db_xref": "; ".join(quals.get("db_xref", [])),
                "note": "; ".join(quals.get("note", [])),
                "record_id": record.id,
                "organism": record.annotations.get("organism", ""),
                "strand": strand,
                "start": str(start),
                "end": str(end),
            }

    return annotations


def parse_gff_annotations(path: Path) -> dict[str, dict[str, str]]:
    """Parse a GFF3 file and extract annotation per locus tag / gene ID.

    Args:
        path: Path to the GFF3 file.

    Returns:
        Dictionary mapping ID -> annotation dict.
    """
    annotations: dict[str, dict[str, str]] = {}

    with path.open() as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            feat_type = parts[2]
            if feat_type not in {"CDS", "gene", "mRNA"}:
                continue

            attrs_str = parts[8]
            attrs: dict[str, str] = {}
            for attr in attrs_str.split(";"):
                if "=" in attr:
                    key, val = attr.split("=", 1)
                    attrs[key.strip()] = val.strip()

            # Determine ID
            locus_tag = attrs.get("locus_tag", attrs.get("ID", attrs.get("Name", "")))
            if not locus_tag:
                continue

            strand = parts[6]
            annotations[locus_tag] = {
                "gene": attrs.get("gene", attrs.get("Name", "")),
                "product": attrs.get("product", ""),
                "protein_id": attrs.get("protein_id", ""),
                "ec_number": attrs.get("ec_number", ""),
                "db_xref": attrs.get("Dbxref", ""),
                "note": attrs.get("Note", ""),
                "record_id": parts[0],
                "organism": "",
                "strand": strand,
                "start": parts[3],
                "end": parts[4],
            }

    return annotations


# ---------------------------------------------------------------------------
# Keyword flagging
# ---------------------------------------------------------------------------


def compile_keyword_patterns(keywords: list[str]) -> list[re.Pattern[str]]:
    """Compile keyword list into regex patterns.

    Args:
        keywords: List of keywords or regex patterns.

    Returns:
        List of compiled regex patterns for case-insensitive matching.
    """
    return [re.compile(re.escape(kw), re.IGNORECASE) for kw in keywords]


def flag_gene(
    annotation: dict[str, str],
    amr_patterns: list[re.Pattern[str]],
    qs_patterns: list[re.Pattern[str]],
    bgc_patterns: list[re.Pattern[str]],
    mge_patterns: list[re.Pattern[str]],
    custom_patterns: list[re.Pattern[str]],
) -> dict[str, str]:
    """Flag a gene based on its annotation.

    Args:
        annotation: Annotation dictionary for the gene.
        amr_patterns: AMR keyword regex patterns.
        qs_patterns: Quorum sensing keyword regex patterns.
        bgc_patterns: Biosynthetic keyword regex patterns.
        mge_patterns: Mobile genetic element keyword regex patterns.
        custom_patterns: User-defined keyword regex patterns.

    Returns:
        Dictionary with flag columns.
    """
    searchable = " ".join(
        [
            annotation.get("product", ""),
            annotation.get("gene", ""),
            annotation.get("note", ""),
            annotation.get("db_xref", ""),
        ]
    ).lower()

    flags: dict[str, str] = {
        "flag_amr": "",
        "flag_qs": "",
        "flag_bgc": "",
        "flag_mge": "",
        "flag_custom": "",
    }

    for pattern in amr_patterns:
        if pattern.search(searchable):
            flags["flag_amr"] = "AMR"
            break

    for pattern in qs_patterns:
        if pattern.search(searchable):
            flags["flag_qs"] = "QS"
            break

    for pattern in bgc_patterns:
        if pattern.search(searchable):
            flags["flag_bgc"] = "BGC"
            break

    for pattern in mge_patterns:
        if pattern.search(searchable):
            flags["flag_mge"] = "MGE"
            break

    for pattern in custom_patterns:
        if pattern.search(searchable):
            flags["flag_custom"] = pattern.pattern.replace("\\", "")
            break

    return flags


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

_ANNOTATION_FIELDS = [
    "gene",
    "product",
    "protein_id",
    "ec_number",
    "db_xref",
    "note",
    "record_id",
    "organism",
    "strand",
    "start",
    "end",
]

_FLAG_FIELDS = ["flag_amr", "flag_qs", "flag_bgc", "flag_mge", "flag_custom"]


def main(
    de_results: Path,
    id_col: str = "locus_tag",
    gbff_dir: Path | None = None,
    gff: Path | None = None,
    custom_keywords: list[str] | None = None,
    threads: int | None = None,
    output: Path | None = None,
) -> None:
    """Annotate DE results with genomic feature information.

    Args:
        de_results: Path to DESeq2/edgeR results TSV.
        id_col: Column name for gene identifiers.
        gbff_dir: Directory containing GBFF/GBK files.
        gff: Path to a GFF3 annotation file.
        custom_keywords: User-supplied keywords for flagging.
        threads: Number of threads.
        output: Output TSV path, or None for stdout.
    """
    if not gbff_dir and not gff:
        print("Error: Provide either --gbff-dir or --gff", file=sys.stderr)
        sys.exit(1)

    # Load DE results
    with de_results.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            print("Error: Empty results file", file=sys.stderr)
            sys.exit(1)
        if id_col not in reader.fieldnames:
            print(f"Error: Column '{id_col}' not found. Available: {', '.join(reader.fieldnames)}", file=sys.stderr)
            sys.exit(1)
        de_rows = list(reader)
        original_fields = list(reader.fieldnames)

    print(f"Loaded {len(de_rows)} DE results", file=sys.stderr)

    # Build annotation index
    annotations: dict[str, dict[str, str]] = {}

    if gbff_dir:
        gbff_files = list(gbff_dir.glob("**/*.gbff")) + list(gbff_dir.glob("**/*.gbk")) + list(gbff_dir.glob("**/*.gb"))
        if not gbff_files:
            print(f"Error: No GBFF/GBK files found in {gbff_dir}", file=sys.stderr)
            sys.exit(1)

        print(f"Parsing {len(gbff_files)} annotation file(s)...", file=sys.stderr)
        with ThreadPoolExecutor(max_workers=threads) as executor:
            futures = {executor.submit(parse_gbff_annotations, p): p for p in gbff_files}
            for future in as_completed(futures):
                path = futures[future]
                try:
                    result = future.result()
                    annotations.update(result)
                    print(f"  {path.name}: {len(result)} CDS features", file=sys.stderr)
                except Exception as e:
                    print(f"  Error parsing {path.name}: {e}", file=sys.stderr)
    elif gff:
        print(f"Parsing {gff.name}...", file=sys.stderr)
        annotations = parse_gff_annotations(gff)

    print(f"Annotation index: {len(annotations)} locus tags", file=sys.stderr)

    # Compile keyword patterns
    amr_patterns = compile_keyword_patterns(_DEFAULT_AMR_KEYWORDS)
    qs_patterns = compile_keyword_patterns(_DEFAULT_QS_KEYWORDS)
    bgc_patterns = compile_keyword_patterns(_DEFAULT_BGC_KEYWORDS)
    mge_patterns = compile_keyword_patterns(_DEFAULT_MGE_KEYWORDS)
    custom_patterns = compile_keyword_patterns(custom_keywords) if custom_keywords else []

    # Annotate each DE result
    out_fields = original_fields + _ANNOTATION_FIELDS + _FLAG_FIELDS
    annotated_rows: list[dict[str, str]] = []
    matched = 0

    for row in de_rows:
        gene_id = row.get(id_col, "").strip()
        anno = annotations.get(gene_id, {})

        if anno:
            matched += 1

        # Add annotation fields
        for field in _ANNOTATION_FIELDS:
            row[field] = anno.get(field, "")

        # Add flags
        flags = flag_gene(anno, amr_patterns, qs_patterns, bgc_patterns, mge_patterns, custom_patterns)
        row.update(flags)

        annotated_rows.append(row)

    # Output
    out_handle = output.open("w", newline="") if output else sys.stdout
    writer = csv.DictWriter(out_handle, fieldnames=out_fields, delimiter="\t", extrasaction="ignore")
    writer.writeheader()
    writer.writerows(annotated_rows)

    if output:
        out_handle.close()
        print(f"\nResults written to: {output}", file=sys.stderr)

    print(f"Matched: {matched}/{len(de_rows)} ({matched / len(de_rows) * 100:.1f}%)", file=sys.stderr)

    # Summary of flagged genes
    amr_count = sum(1 for r in annotated_rows if r.get("flag_amr"))
    qs_count = sum(1 for r in annotated_rows if r.get("flag_qs"))
    bgc_count = sum(1 for r in annotated_rows if r.get("flag_bgc"))
    mge_count = sum(1 for r in annotated_rows if r.get("flag_mge"))
    custom_count = sum(1 for r in annotated_rows if r.get("flag_custom"))

    print(
        f"Flagged: AMR={amr_count}, QS={qs_count}, BGC={bgc_count}, MGE={mge_count}, Custom={custom_count}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Annotate DESeq2/edgeR differential expression results with gene "
            "product descriptions, COG categories, and biosynthetic/AMR/QS "
            "flags from GBFF or GFF annotation files."
        ),
    )

    parser.add_argument("de_results", type=Path, help="DESeq2/edgeR results TSV (must contain a gene ID column)")
    parser.add_argument(
        "--id-col",
        default="locus_tag",
        help="Column name for gene identifiers (default: locus_tag)",
    )
    parser.add_argument("--gbff-dir", type=Path, default=None, help="Directory containing GBFF/GBK annotation files")
    parser.add_argument("--gff", type=Path, default=None, help="GFF3 annotation file (alternative to --gbff-dir)")
    parser.add_argument(
        "--keywords",
        default=None,
        help="Comma-separated keywords for flagging genes of interest (default: built-in AMR/QS/BGC list)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Number of threads for parallel GBFF parsing",
    )
    parser.add_argument("--output", type=Path, default=None, help="Output TSV file path (default: stdout)")

    args = parser.parse_args()

    custom_kw = [k.strip() for k in args.keywords.split(",")] if args.keywords else None

    main(
        de_results=args.de_results,
        id_col=args.id_col,
        gbff_dir=args.gbff_dir,
        gff=args.gff,
        custom_keywords=custom_kw,
        threads=args.threads,
        output=args.output,
    )
