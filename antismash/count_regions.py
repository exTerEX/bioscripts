## Given a bunch of antismash results, count the BGC regions
#
## Usage:
#   $ python antismash/count_regions.py -h
#   usage: count_regions.py [-h] [--by_contig] [--split_hybrids] directory output
#
#   Given a bunch of antismash results, count the BGC regions
#
#   positional arguments:
#     directory         Directory containing antiSMASH directories
#     output            Desired path/to/filename for the output TSV
#
#   options:
#     -h, --help        show this help message and exit
#     --by_contig       Count regions per each individual contig rather than per assembly
#     --split_hybrids   Count each hybrid region multiple times, once for each
#                       constituent BGC class. The total_count column is unaffected.
#     --threads THREADS Number of threads to use for parallel processing. Defaults to number of CPUs.

import argparse
import csv
import json
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def parse_json(path: Path) -> tuple[str, dict[str, list[list[str]]], dict[str, str]]:
    """Parse an antiSMASH JSON file and extract BGC region data.

    Args:
        path: Path to the antiSMASH JSON file to parse.

    Returns:
        A tuple containing:
            - str: The input file name from the JSON data.
            - dict[str, list[list[str]]]: A dictionary mapping contig names to
              lists of product lists for each area in the contig.
            - dict[str, str]: A dictionary mapping contig names to their
              descriptions (empty string if no description exists).
    """
    with path.open() as f:
        data = json.load(f)

    by_contig = {record["name"]: [area["products"] for area in record["areas"]] for record in data["records"]}

    descriptions = {record["name"]: record.get("description", "") for record in data["records"]}

    return data["input_file"], by_contig, descriptions


def tabulate(
    type_dict: dict[str, dict[str, list[list[str]]]],
    descriptions: dict[str, dict[str, str]],
    contig: bool = False,
    split_hybrids: bool = False,
) -> list[dict[str, str | int]]:
    """Tabulate BGC counts from parsed antiSMASH data.
    Converts a nested dictionary of genome/contig/region data into a list of
    dictionaries suitable for tabular output.

    Args:
        type_dict: Nested dictionary mapping genome names to contigs to lists of
            regions, where each region is a list of product type strings.
        descriptions: Nested dictionary mapping genome names to contigs to their
            description strings.
        contig: If True, produce one row per contig. If False, aggregate counts
            across all contigs for each genome. Defaults to False.
        split_hybrids: If True, count each product type in hybrid regions
            separately. If False, count multi-product regions as "hybrid".
            Defaults to False.

    Returns:
        A list of dictionaries, each containing:
            - Product type counts as key-value pairs
            - "record": The genome name (or "genome|contig" if contig=True)
            - "total_count": Total number of regions
            - "description": Description string for the record
    """
    table_list = []

    for genome, g_prods in type_dict.items():
        for cont, regions in g_prods.items():
            counts: Counter[str] = Counter()
            for region in regions:
                if len(region) > 1 and not split_hybrids:
                    counts["hybrid"] += 1
                else:
                    counts.update(region)

            if contig:
                table_list.append(
                    {
                        **counts,
                        "record": f"{genome}|{cont}",
                        "total_count": len(regions),
                        "description": descriptions[genome][cont],
                    }
                )

        if not contig:
            # Aggregate counts across all contigs for this genome
            genome_counts: Counter[str] = Counter()
            for _, regions in g_prods.items():
                for region in regions:
                    if len(region) > 1 and not split_hybrids:
                        genome_counts["hybrid"] += 1
                    else:
                        genome_counts.update(region)

            first_contig = next(iter(g_prods))
            num_contigs = len(g_prods)
            record_suffix = "s" if num_contigs > 1 else ""
            description_text = f"{descriptions[genome][first_contig]} [{num_contigs} total record{record_suffix}]"
            table_list.append(
                {
                    **genome_counts,
                    "record": genome,
                    "total_count": sum(len(r) for r in g_prods.values()),
                    "description": description_text,
                }
            )

    return table_list


def main(
    directory: Path,
    outpath: Path,
    contig: bool = False,
    split_hybrid: bool = False,
    threads: int | None = None,
):
    """Main function to process antiSMASH results and output BGC counts."""
    by_genome: dict[str, dict[str, list[list[str]]]] = {}
    descriptions: dict[str, dict[str, str]] = {}

    json_files = list(directory.glob("*/*.json"))

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(parse_json, path): path for path in json_files}

        for future in as_completed(futures):
            genome, types, description = future.result()
            by_genome[genome] = types
            descriptions[genome] = description

    table_list = tabulate(by_genome, descriptions, contig, split_hybrid)

    # Collect all product types from the data
    all_products = {k for d in table_list for k in d} - {
        "record",
        "total_count",
        "hybrid",
        "description",
    }

    fieldnames = ["record", "total_count", *sorted(all_products)]
    if not split_hybrid:
        fieldnames.append("hybrid")
    fieldnames.append("description")

    with outpath.open("w") as outf:
        writer = csv.DictWriter(outf, fieldnames=fieldnames, delimiter="\t", restval=0)
        writer.writeheader()
        writer.writerows(table_list)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Given a bunch of antismash results, count the BGC regions")

    parser.add_argument("directory", type=Path, help="Directory containing antiSMASH directories")
    parser.add_argument("output", type=Path, help="Desired path+name for the output TSV")
    parser.add_argument(
        "--by_contig",
        action="store_true",
        help="Count regions per each individual contig rather than per assembly",
    )
    parser.add_argument(
        "--split_hybrids",
        action="store_true",
        help=(
            "Count each hybrid region multiple times, once for each "
            "constituent BGC class. The total_count column is unaffected."
        ),
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Number of threads to use for parallel processing. Defaults to number of CPUs.",
    )

    args = parser.parse_args()

    main(args.directory, args.output, args.by_contig, args.split_hybrids, args.threads)
