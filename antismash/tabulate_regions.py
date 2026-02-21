## Given a bunch of antismash results, tabulate BGC regions
#
## Usage:
#   $ python antismash/tabulate_regions.py -h
#   usage: tabulate_regions.py [-h] [--threads THREADS] directory output
#
#   Given a bunch of antismash results, tabulate BGC regions
#
#   positional arguments:
#     directory          Directory containing antiSMASH directories
#     output             Desired path/to/filename for the output TSV
#
#   options:
#     -h, --help         show this help message and exit
#     --threads THREADS  Number of threads to use for parallel processing. Defaults to number of CPUs.

import argparse
import csv
import json
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def get_knownclusterblast_info(knownclusterblast: list | None, index: int) -> dict[str, str]:
    """Extract knownclusterblast hit info for a region.

    Args:
        knownclusterblast: The knownclusterblast results list, or None if not present.
        index: The index of the region to get hits for.

    Returns:
        Dictionary with knownclusterblast_hit, knownclusterblast_accession, and knownclusterblast_similarity keys.
    """
    empty_result = {"knownclusterblast_hit": "", "knownclusterblast_accession": "", "knownclusterblast_similarity": ""}

    if not knownclusterblast:
        return empty_result

    hits = knownclusterblast[index]["ranking"]
    if not hits:
        return empty_result

    similarity = hits[0][1]["similarity"]
    if similarity <= 15:
        return empty_result

    match similarity:
        case s if s > 75:
            sim_level = "high"
        case s if s > 50:
            sim_level = "medium"
        case _:
            sim_level = "low"

    return {
        "knownclusterblast_hit": hits[0][0]["description"],
        "knownclusterblast_accession": hits[0][0]["accession"],
        "knownclusterblast_similarity": sim_level,
    }


def extract_knownclusterblast(record: dict) -> list | None:
    """Safely extract knownclusterblast results from a record.
    Attempts to navigate the nested dictionary structure of an antiSMASH record
    to retrieve knownclusterblast results.

    Args:
        record: A dictionary containing antiSMASH analysis results with a nested
            structure containing modules and clusterblast information.

    Returns:
        A list of knownclusterblast results if found, or None if the required
        keys are missing or if any errors occur during extraction.
    """
    try:
        return record["modules"]["antismash.modules.clusterblast"]["knowncluster"]["results"]
    except (KeyError, TypeError, AttributeError):
        return None


def parse_json(path: Path) -> list[dict[str, str]]:
    """Parse an antiSMASH JSON file and extract region data.

    Args:
        path: Path to the antiSMASH JSON file.

    Returns:
        List of dictionaries containing region information.
    """
    with path.open() as f:
        data = json.load(f)

    result_list: list[dict[str, str]] = []

    for record in data["records"]:
        if not record["areas"]:
            continue

        regions = [feat for feat in record["features"] if feat["type"] == "region"]
        knownclusterblast = extract_knownclusterblast(record)

        for index, region in enumerate(regions):
            start, end = re.findall(r"\d+", region["location"])
            qualifiers = region["qualifiers"]

            region_dict = {
                "file": path.stem,
                "record_id": record["name"],
                "region": qualifiers["region_number"][0],
                "start": start,
                "end": end,
                "contig_edge": qualifiers["contig_edge"][0],
                "product": " / ".join(qualifiers["product"]),
                "record_desc": record["description"],
            } | get_knownclusterblast_info(knownclusterblast, index)

            result_list.append(region_dict)

    return result_list


def main(directory: Path, output: Path, threads: int | None = None):
    record_infos: list[dict[str, str]] = []

    json_files = list(directory.glob("*/*.json"))

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(parse_json, path): path for path in json_files}

        for future in as_completed(futures):
            record_infos.extend(future.result())

    fieldnames = [
        "file",
        "record_id",
        "region",
        "start",
        "end",
        "contig_edge",
        "product",
        "knownclusterblast_hit",
        "knownclusterblast_accession",
        "knownclusterblast_similarity",
        "record_desc",
    ]
    with output.open("w") as outf:
        writer = csv.DictWriter(outf, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(record_infos)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Given a bunch of antismash results, tabulate BGC regions")

    parser.add_argument("directory", type=Path, help="Directory containing antiSMASH directories")
    parser.add_argument("output", type=Path, help="Desired path/to/filename for the output TSV")
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Number of threads to use for parallel processing. Defaults to number of CPUs.",
    )

    args = parser.parse_args()

    main(args.directory, args.output, args.threads)
