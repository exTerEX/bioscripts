## Merge a GenBank genome file with fuzznuc motif search results
#
## Usage:
#   $ python fuzznuc/merge_gbk.py -h
#   usage: merge_gbk.py [-h] [--descriptor DESCRIPTOR] genome fuzznuc output
#
#   Merge a GenBank genome file with fuzznuc motif search results
#
#   positional arguments:
#     genome                Path to the complete genome GenBank file
#     fuzznuc               Path to the fuzznuc results GenBank file
#     output                Desired path/to/filename for the merged GenBank file
#
#   options:
#     -h, --help            show this help message and exit
#     --descriptor DESCRIPTOR
#                           Optional descriptor to add to the note (before fuzznuc specification)

import argparse
import re
from pathlib import Path

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature


def parse_fuzznuc_features(fuzznuc_path: Path, descriptor: str | None = None) -> list[SeqFeature]:
    """Parse fuzznuc GenBank file and convert features to misc_binding.

    Args:
        fuzznuc_path: Path to the fuzznuc results GenBank file.
        descriptor: Optional descriptor to add to the note.

    Returns:
        List of SeqFeature objects with corrected locations and types.
    """
    features: list[SeqFeature] = []

    with fuzznuc_path.open() as f:
        content = f.read()

    # Parse misc_feature entries from the fuzznuc file
    # Pattern matches both regular and complement features
    feature_pattern = re.compile(
        r"misc_feature\s+(?:complement\()?(\d+)\.\.(\d+)\)?[^\n]*\n"
        r'\s+/note="([^"]*)"',
        re.MULTILINE,
    )

    for match in feature_pattern.finditer(content):
        pos1, pos2, _ = int(match.group(1)), int(match.group(2)), match.group(3)

        # Determine strand and ensure correct ordering (low..high)
        if pos1 > pos2:
            # Complement feature with reversed coordinates
            start, end = pos2 - 1, pos1  # Convert to 0-based
            strand = -1
        else:
            start, end = pos1 - 1, pos2  # Convert to 0-based
            # Check if it was a complement in the original
            strand = -1 if "complement" in match.group(0) else 1

        location = FeatureLocation(start, end, strand=strand)

        # Build the note
        note_parts = []
        if descriptor:
            note_parts.append(descriptor)
        note_parts.append("predicted using fuzznuc")

        new_note = "; ".join(note_parts)

        feature = SeqFeature(
            location=location,
            type="misc_binding",
            qualifiers={"note": [new_note]},
        )
        features.append(feature)

    return features


def merge_genbank(
    genome_path: Path,
    fuzznuc_path: Path,
    output_path: Path,
    descriptor: str | None = None,
):
    """Merge genome GenBank with fuzznuc features and write output.

    Args:
        genome_path: Path to the complete genome GenBank file.
        fuzznuc_path: Path to the fuzznuc results GenBank file.
        output_path: Path for the merged output GenBank file.
        descriptor: Optional descriptor to add to fuzznuc feature notes.
    """
    # Read the genome
    with genome_path.open() as f:
        genome_record = SeqIO.read(f, "genbank")

    # Parse fuzznuc features
    fuzznuc_features = parse_fuzznuc_features(fuzznuc_path, descriptor)

    # Combine all features
    all_features = list(genome_record.features) + fuzznuc_features

    # Sort features by start position, then by end position
    all_features.sort(key=lambda f: (int(f.location.start), int(f.location.end)))

    # Create new record with merged features
    genome_record.features = all_features

    # Write output
    with output_path.open("w") as f:
        SeqIO.write(genome_record, f, "genbank")

    print(f"Merged {len(fuzznuc_features)} fuzznuc features into genome")
    print(f"Total features: {len(all_features)}")
    print(f"Output written to: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Merge a GenBank genome file with fuzznuc motif search results")

    parser.add_argument("genome", type=Path, help="Path to the complete genome GenBank flat file")
    parser.add_argument("fuzznuc", type=Path, help="Path to the fuzznuc results GenBank flat file")
    parser.add_argument(
        "output",
        type=Path,
        help="Desired path/to/filename for the merged GenBank flat file",
    )
    parser.add_argument(
        "--descriptor",
        type=str,
        default=None,
        help="Optional descriptor to add to the note (before fuzznuc specification)",
    )

    args = parser.parse_args()

    merge_genbank(args.genome, args.fuzznuc, args.output, args.descriptor)


if __name__ == "__main__":
    main()
