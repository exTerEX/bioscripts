## Filter cblaster session JSON to keep S. anginosus group organisms
#
## Usage:
#   $ python cblaster/filter_session_com.py -h
#   usage: filter_session_com.py [-h] [--output OUTPUT]
#                                [--species-pattern SPECIES_PATTERN]
#                                [--deduplicate-scaffolds]
#                                input
#
#   Keep only organisms whose name matches an anginosus-group regex and
#   optionally deduplicate scaffold accessions with preference for RefSeq
#   (NZ_ prefixed) accessions.

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any


def load_session(path: Path) -> dict[str, Any]:
    """Load cblaster session JSON from disk."""
    with path.open() as fh:
        data = json.load(fh)

    if not isinstance(data, dict):
        raise ValueError("Top-level JSON value must be an object")
    if "organisms" not in data or not isinstance(data["organisms"], list):
        raise ValueError("JSON must contain an 'organisms' list")

    return data


def filter_organisms_by_species(
    organisms: list[dict[str, Any]],
    species_pattern: re.Pattern[str],
) -> list[dict[str, Any]]:
    """Keep only organisms whose 'name' matches species_pattern."""
    kept: list[dict[str, Any]] = []

    for index, organism in enumerate(organisms):
        name = str(organism.get("name", ""))
        if species_pattern.search(name):
            print(f"PRESERVE index {index}, containing species: {name}", file=sys.stderr)
            kept.append(organism)

    return kept


def deduplicate_scaffolds_prefer_refseq(organism: dict[str, Any]) -> dict[str, Any]:
    """Remove duplicate scaffold accessions while preferring NZ_ accessions."""
    scaffolds_raw = organism.get("scaffolds", [])
    if not isinstance(scaffolds_raw, list):
        return organism

    best_by_base: dict[str, dict[str, Any]] = {}

    for scaffold in scaffolds_raw:
        if not isinstance(scaffold, dict):
            continue
        accession = scaffold.get("accession")
        if not isinstance(accession, str) or not accession:
            continue

        base_accession = accession.removeprefix("NZ_")
        current = best_by_base.get(base_accession)

        # Prefer NZ_ version for duplicate accessions sharing the same base.
        if current is None:
            best_by_base[base_accession] = scaffold
            continue

        current_acc = str(current.get("accession", ""))
        if accession.startswith("NZ_") and not current_acc.startswith("NZ_"):
            best_by_base[base_accession] = scaffold

    keep_accessions = {str(s.get("accession", "")) for s in best_by_base.values()}
    kept_scaffolds = [
        scaffold
        for scaffold in scaffolds_raw
        if isinstance(scaffold, dict) and str(scaffold.get("accession", "")) in keep_accessions
    ]

    organism["scaffolds"] = kept_scaffolds
    return organism


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Filter cblaster session JSON to keep anginosus-group organisms and "
            "optionally deduplicate scaffolds with RefSeq preference."
        )
    )

    parser.add_argument("input", type=Path, help="Path to input session JSON")
    parser.add_argument("--output", type=Path, default=None, help="Output JSON path (default: overwrite input)")
    parser.add_argument(
        "--species-pattern",
        type=str,
        default=r"(anginosus|intermedius|constellatus|milleri)",
        help="Regex used to keep organisms by name (default: anginosus-group example)",
    )
    parser.add_argument(
        "--regex-pattern",
        dest="species_pattern",
        type=str,
        help="Alias for --species-pattern",
    )
    parser.add_argument(
        "--deduplicate-scaffolds",
        action="store_true",
        help="Deduplicate scaffold accessions with preference for NZ_ accessions",
    )

    return parser.parse_args()


def main() -> None:
    """Entry point."""
    args = parse_args()

    if not args.input.exists():
        print(f"Error: input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    try:
        species_re = re.compile(args.species_pattern, flags=re.IGNORECASE)
    except re.error as exc:
        print(f"Error: invalid species regex: {exc}", file=sys.stderr)
        sys.exit(1)

    try:
        session = load_session(args.input)
    except Exception as exc:  # noqa: BLE001
        print(f"Error: failed to load session JSON: {exc}", file=sys.stderr)
        sys.exit(1)

    organisms = session["organisms"]
    assert isinstance(organisms, list)

    filtered = filter_organisms_by_species(organisms, species_re)

    if args.deduplicate_scaffolds:
        filtered = [deduplicate_scaffolds_prefer_refseq(org) for org in filtered if isinstance(org, dict)]

    session["organisms"] = filtered

    output = args.output or args.input
    with output.open("w") as fh:
        json.dump(session, fh, indent=2)
        fh.write("\n")

    print(
        f"Wrote {len(filtered)} organism(s) to {output}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
