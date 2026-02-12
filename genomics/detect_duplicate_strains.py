## Detect potential duplicate strains from GenBank files or assembly accessions
#
## Usage:
#   $ python misc/detect_duplicate_strains.py -h
#   usage: detect_duplicate_strains.py [-h] [--email EMAIL] [--api-key API_KEY]
#                                      [--threads THREADS] [--similarity SIMILARITY]
#                                      [--output OUTPUT] [--sequence-analysis] [--kmer-size KMER_SIZE]
#                                      [--sketch-size SKETCH_SIZE] [--ani-threshold ANI_THRESHOLD]
#                                      input
#
#   Detect potential duplicate strains from GBFF files or GenBank assembly IDs
#
#   positional arguments:
#     input                 Directory containing GBFF files, or file with assembly IDs (one per line)
#
#   options:
#     -h, --help            show this help message and exit
#     --email EMAIL         Email address for NCBI Entrez (required for fetching)
#     --api-key API_KEY     NCBI API key for higher rate limits
#     --threads THREADS     Number of threads for parallel processing
#     --similarity SIMILARITY
#                           Minimum sequence similarity threshold (0-1) for duplicate detection (default: 0.99)
#     --output OUTPUT       Output report file path (default: duplicate_report.tsv)
#     --sequence-analysis   Enable deep sequence analysis using MinHash sketching and k-mer comparison
#     --kmer-size KMER_SIZE
#                           K-mer size for sequence analysis (default: 21)
#     --sketch-size SKETCH_SIZE
#                           MinHash sketch size for fast similarity estimation (default: 1000)
#     --ani-threshold ANI_THRESHOLD
#                           ANI threshold for considering genomes as duplicates (default: 0.999)

import argparse
import csv
import gzip
import hashlib
import math
import re
import struct
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path
from typing import Iterator

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

# Constants for MinHash
MAX_HASH = (1 << 64) - 1
MERSENNE_PRIME = (1 << 61) - 1


class UnionFind:
    """Union-Find (Disjoint Set) data structure for clustering connected strains.

    Uses path compression and union by rank for near-constant time operations.
    """

    def __init__(self):
        self.parent: dict[str, str] = {}
        self.rank: dict[str, int] = {}

    def find(self, x: str) -> str:
        """Find the root representative of the set containing x.

        Args:
            x: Element to find.

        Returns:
            Root representative of x's set.
        """
        if x not in self.parent:
            self.parent[x] = x
            self.rank[x] = 0
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])  # Path compression
        return self.parent[x]

    def union(self, x: str, y: str) -> None:
        """Merge the sets containing x and y.

        Args:
            x: First element.
            y: Second element.
        """
        rx, ry = self.find(x), self.find(y)
        if rx == ry:
            return
        # Union by rank
        if self.rank[rx] < self.rank[ry]:
            rx, ry = ry, rx
        self.parent[ry] = rx
        if self.rank[rx] == self.rank[ry]:
            self.rank[rx] += 1

    def connected(self, x: str, y: str) -> bool:
        """Check if x and y are in the same set.

        Args:
            x: First element.
            y: Second element.

        Returns:
            True if x and y share the same root.
        """
        return self.find(x) == self.find(y)

    def get_clusters(self) -> dict[str, list[str]]:
        """Return all clusters as a mapping from root to members.

        Returns:
            Dictionary mapping cluster root to list of members.
        """
        clusters: dict[str, list[str]] = {}
        for item in self.parent:
            root = self.find(item)
            clusters.setdefault(root, []).append(item)
        return clusters


@dataclass
class StrainMetadata:
    """Container for strain metadata extracted from GenBank records."""

    source_id: str  # filename or assembly accession
    organism: str = ""
    strain: str = ""
    isolate: str = ""
    culture_collection: str = ""
    type_material: str = ""
    biosample: str = ""
    bioproject: str = ""
    assembly_accession: str = ""
    infraspecific_name: str = ""
    isolation_source: str = ""
    host: str = ""
    country: str = ""
    collection_date: str = ""
    total_length: int = 0
    num_contigs: int = 0
    gc_content: float = 0.0
    sequence_checksums: list[str] = field(default_factory=list)
    contig_lengths: list[int] = field(default_factory=list)
    # Sequence analysis fields
    minhash_sketch: list[int] = field(default_factory=list)
    kmer_frequencies: dict[str, int] = field(default_factory=dict)
    full_sequence: str = ""  # Stored only when sequence analysis is enabled


@dataclass
class DuplicateMatch:
    """Container for a potential duplicate match between two strains."""

    strain1: str
    strain2: str
    match_reasons: list[str]
    confidence: str  # "high", "medium", "low"
    metadata1: StrainMetadata
    metadata2: StrainMetadata
    # Sequence similarity metrics
    minhash_similarity: float = 0.0
    estimated_ani: float = 0.0
    shared_kmers_ratio: float = 0.0


def calculate_gc_content(sequence: str) -> float:
    """Calculate GC content of a sequence.

    Args:
        sequence: DNA sequence string.

    Returns:
        GC content as a fraction (0-1).
    """
    if not sequence:
        return 0.0
    sequence = sequence.upper()
    gc_count = sequence.count("G") + sequence.count("C")
    total = len(sequence)
    return gc_count / total if total > 0 else 0.0


def compute_sequence_hash(sequence: str) -> str:
    """Compute MD5 hash of a sequence for comparison.

    Args:
        sequence: DNA sequence string.

    Returns:
        MD5 hash string.
    """
    return hashlib.md5(sequence.upper().encode()).hexdigest()


# =============================================================================
# Sequence Analysis Functions (MinHash and k-mer based)
# =============================================================================


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence.

    Args:
        sequence: DNA sequence string.

    Returns:
        Reverse complement sequence.
    """
    complement = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(complement.get(base, "N") for base in reversed(sequence.upper()))


def canonical_kmer(kmer: str) -> str:
    """Return the canonical (lexicographically smaller) k-mer.

    Args:
        kmer: K-mer sequence.

    Returns:
        Canonical k-mer (original or reverse complement, whichever is smaller).
    """
    rc = reverse_complement(kmer)
    return min(kmer.upper(), rc)


def extract_kmers(sequence: str, k: int = 21) -> Iterator[str]:
    """Extract all canonical k-mers from a sequence.

    Args:
        sequence: DNA sequence string.
        k: K-mer size (default: 21).

    Yields:
        Canonical k-mers from the sequence.
    """
    sequence = sequence.upper()
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i : i + k]
        # Skip k-mers containing N
        if "N" not in kmer:
            yield canonical_kmer(kmer)


def murmurhash3_x64_128(data: bytes, seed: int = 42) -> int:
    """Simplified MurmurHash3-like hash for k-mer hashing.

    Uses Python's built-in hash with additional mixing for better distribution.

    Args:
        data: Bytes to hash.
        seed: Hash seed (default: 42).

    Returns:
        64-bit hash value.
    """
    # Use struct to create a more robust hash
    h = hashlib.md5(data + struct.pack("I", seed)).digest()
    # Convert first 8 bytes to integer
    return struct.unpack("<Q", h[:8])[0]


def kmer_hash(kmer: str, seed: int = 42) -> int:
    """Hash a k-mer to a 64-bit integer.

    Args:
        kmer: K-mer sequence.
        seed: Hash seed.

    Returns:
        64-bit hash value.
    """
    return murmurhash3_x64_128(kmer.encode(), seed)


def compute_minhash_sketch(
    sequence: str, k: int = 21, sketch_size: int = 1000, seed: int = 42
) -> list[int]:
    """Compute MinHash sketch of a genome sequence.

    MinHash sketching provides fast approximate Jaccard similarity estimation.
    The sketch consists of the smallest hash values from all k-mers.

    Args:
        sequence: Concatenated genome sequence.
        k: K-mer size (default: 21).
        sketch_size: Number of hash values in sketch (default: 1000).
        seed: Random seed for hashing.

    Returns:
        List of sketch_size smallest hash values (sorted).
    """
    if len(sequence) < k:
        return []

    # Collect all k-mer hashes
    hashes = set()
    for kmer in extract_kmers(sequence, k):
        h = kmer_hash(kmer, seed)
        hashes.add(h)

    # Return the smallest sketch_size hashes
    sorted_hashes = sorted(hashes)
    return sorted_hashes[:sketch_size]


def compute_kmer_frequencies(
    sequence: str, k: int = 21, sample_size: int = 10000
) -> dict[str, int]:
    """Compute k-mer frequency distribution (sampled for large genomes).

    Args:
        sequence: Genome sequence.
        k: K-mer size.
        sample_size: Maximum number of k-mers to sample.

    Returns:
        Dictionary of k-mer counts.
    """
    if len(sequence) < k:
        return {}

    # For large genomes, sample k-mers
    total_kmers = len(sequence) - k + 1
    step = max(1, total_kmers // sample_size)

    frequencies: Counter[str] = Counter()
    for i in range(0, total_kmers, step):
        kmer = sequence[i : i + k].upper()
        if "N" not in kmer:
            frequencies[canonical_kmer(kmer)] += 1

    return dict(frequencies)


def estimate_jaccard_from_minhash(sketch1: list[int], sketch2: list[int]) -> float:
    """Estimate Jaccard similarity from two MinHash sketches.

    Args:
        sketch1: First MinHash sketch.
        sketch2: Second MinHash sketch.

    Returns:
        Estimated Jaccard similarity (0-1).
    """
    if not sketch1 or not sketch2:
        return 0.0

    set1 = set(sketch1)
    set2 = set(sketch2)

    intersection = len(set1 & set2)
    union = len(set1 | set2)

    return intersection / union if union > 0 else 0.0


def estimate_ani_from_jaccard(jaccard: float, k: int = 21) -> float:
    """Estimate Average Nucleotide Identity (ANI) from Jaccard similarity.

    Uses the Mash distance formula: D = -1/k * ln(2J / (1+J))
    Then converts to ANI: ANI = 1 - D

    Args:
        jaccard: Jaccard similarity coefficient.
        k: K-mer size used.

    Returns:
        Estimated ANI (0-1), or 0 if calculation fails.
    """
    if jaccard <= 0:
        return 0.0
    if jaccard >= 1:
        return 1.0

    try:
        # Mash distance formula
        distance = -1.0 / k * math.log(2.0 * jaccard / (1.0 + jaccard))
        ani = max(0.0, 1.0 - distance)
        return min(1.0, ani)  # Clamp to [0, 1]
    except (ValueError, ZeroDivisionError):
        return 0.0


def calculate_shared_kmer_ratio(freq1: dict[str, int], freq2: dict[str, int]) -> float:
    """Calculate the ratio of shared k-mers between two genomes.

    Args:
        freq1: K-mer frequencies for first genome.
        freq2: K-mer frequencies for second genome.

    Returns:
        Ratio of shared k-mers to total unique k-mers.
    """
    if not freq1 or not freq2:
        return 0.0

    set1 = set(freq1.keys())
    set2 = set(freq2.keys())

    shared = len(set1 & set2)
    total = len(set1 | set2)

    return shared / total if total > 0 else 0.0


def calculate_kmer_frequency_correlation(
    freq1: dict[str, int], freq2: dict[str, int]
) -> float:
    """Calculate Pearson correlation of k-mer frequencies for shared k-mers.

    High correlation suggests same genome/strain even if not identical.

    Args:
        freq1: K-mer frequencies for first genome.
        freq2: K-mer frequencies for second genome.

    Returns:
        Pearson correlation coefficient (-1 to 1), or 0 if insufficient data.
    """
    shared_kmers = set(freq1.keys()) & set(freq2.keys())
    if len(shared_kmers) < 10:
        return 0.0

    x = [freq1[k] for k in shared_kmers]
    y = [freq2[k] for k in shared_kmers]

    n = len(x)
    mean_x = sum(x) / n
    mean_y = sum(y) / n

    # Calculate covariance and standard deviations
    cov = sum((x[i] - mean_x) * (y[i] - mean_y) for i in range(n))
    std_x = math.sqrt(sum((xi - mean_x) ** 2 for xi in x))
    std_y = math.sqrt(sum((yi - mean_y) ** 2 for yi in y))

    if std_x == 0 or std_y == 0:
        return 0.0

    return cov / (std_x * std_y)


def compute_sequence_analysis(
    metadata: StrainMetadata, k: int = 21, sketch_size: int = 1000
) -> None:
    """Compute sequence analysis metrics for a strain (modifies in place).

    Args:
        metadata: StrainMetadata object to update.
        k: K-mer size.
        sketch_size: MinHash sketch size.
    """
    if not metadata.full_sequence:
        return

    sequence = metadata.full_sequence.upper()

    # Compute MinHash sketch
    metadata.minhash_sketch = compute_minhash_sketch(sequence, k, sketch_size)

    # Compute k-mer frequencies (sampled)
    metadata.kmer_frequencies = compute_kmer_frequencies(sequence, k)


def compare_sequences(
    meta1: StrainMetadata, meta2: StrainMetadata, k: int = 21
) -> tuple[float, float, float, float]:
    """Compare two genomes using sequence analysis.

    Args:
        meta1: First strain metadata.
        meta2: Second strain metadata.
        k: K-mer size used.

    Returns:
        Tuple of (minhash_similarity, estimated_ani, shared_kmer_ratio, kmer_correlation).
    """
    minhash_sim = estimate_jaccard_from_minhash(
        meta1.minhash_sketch, meta2.minhash_sketch
    )
    ani = estimate_ani_from_jaccard(minhash_sim, k)
    shared_ratio = calculate_shared_kmer_ratio(
        meta1.kmer_frequencies, meta2.kmer_frequencies
    )
    correlation = calculate_kmer_frequency_correlation(
        meta1.kmer_frequencies, meta2.kmer_frequencies
    )

    return minhash_sim, ani, shared_ratio, correlation


def get_source_qualifier(record: SeqRecord, key: str, default: str = "") -> str:
    """Extract a qualifier from the source feature of a GenBank record.

    Args:
        record: BioPython SeqRecord object.
        key: Qualifier key to extract.
        default: Default value if key not found.

    Returns:
        The qualifier value or default.
    """
    for feature in record.features:
        if feature.type == "source":
            values = feature.qualifiers.get(key, [default])
            return values[0] if values else default
    return default


def get_dbxref_value(record: SeqRecord, db_name: str) -> str:
    """Extract a specific database cross-reference from source feature.

    Args:
        record: BioPython SeqRecord object.
        db_name: Database name to look for (e.g., "BioSample", "Assembly").

    Returns:
        The database reference value or empty string.
    """
    for feature in record.features:
        if feature.type == "source":
            db_xrefs = feature.qualifiers.get("db_xref", [])
            for xref in db_xrefs:
                if xref.startswith(f"{db_name}:"):
                    return xref.split(":", 1)[1]
    return ""


def extract_metadata_from_records(
    records: list[SeqRecord], source_id: str, store_sequence: bool = False
) -> StrainMetadata:
    """Extract strain metadata from a list of GenBank records.

    Args:
        records: List of BioPython SeqRecord objects.
        source_id: Identifier for the source (filename or accession).
        store_sequence: If True, store full sequence for analysis.

    Returns:
        StrainMetadata object with extracted information.
    """
    metadata = StrainMetadata(source_id=source_id)

    if not records:
        return metadata

    # Use first record for most metadata
    first_record = records[0]

    # Extract from annotations
    metadata.organism = first_record.annotations.get("organism", "")

    # Extract from source feature qualifiers
    metadata.strain = get_source_qualifier(first_record, "strain")
    metadata.isolate = get_source_qualifier(first_record, "isolate")
    metadata.culture_collection = get_source_qualifier(
        first_record, "culture_collection"
    )
    metadata.type_material = get_source_qualifier(first_record, "type_material")
    metadata.infraspecific_name = get_source_qualifier(
        first_record, "sub_species"
    ) or get_source_qualifier(first_record, "variety")
    metadata.isolation_source = get_source_qualifier(first_record, "isolation_source")
    metadata.host = get_source_qualifier(first_record, "host")
    metadata.country = get_source_qualifier(first_record, "country")
    metadata.collection_date = get_source_qualifier(first_record, "collection_date")

    # Extract database cross-references
    metadata.biosample = get_dbxref_value(first_record, "BioSample")
    metadata.bioproject = get_dbxref_value(first_record, "BioProject")
    metadata.assembly_accession = get_dbxref_value(first_record, "Assembly")

    # Calculate sequence statistics
    total_seq = ""
    for record in records:
        seq_str = str(record.seq)
        total_seq += seq_str
        metadata.sequence_checksums.append(compute_sequence_hash(seq_str))
        metadata.contig_lengths.append(len(seq_str))

    metadata.total_length = len(total_seq)
    metadata.num_contigs = len(records)
    metadata.gc_content = calculate_gc_content(total_seq)

    # Store full sequence if requested for sequence analysis
    if store_sequence:
        metadata.full_sequence = total_seq

    return metadata


def parse_gbff_file(path: Path, store_sequence: bool = False) -> StrainMetadata:
    """Parse a GBFF file and extract strain metadata.

    Args:
        path: Path to the GBFF file (can be gzipped).
        store_sequence: If True, store full sequence for analysis.

    Returns:
        StrainMetadata object with extracted information.
    """
    try:
        if path.suffix == ".gz":
            with gzip.open(path, "rt") as handle:
                records = list(SeqIO.parse(handle, "genbank"))
        else:
            with path.open() as handle:
                records = list(SeqIO.parse(handle, "genbank"))

        return extract_metadata_from_records(records, path.name, store_sequence)

    except Exception as e:
        print(f"Error parsing {path}: {e}", file=sys.stderr)
        return StrainMetadata(source_id=path.name)


def fetch_assembly_metadata(
    accession: str, email: str, api_key: str | None = None, store_sequence: bool = False
) -> StrainMetadata:
    """Fetch GenBank data for an assembly accession and extract metadata.

    Args:
        accession: Assembly accession (GCA_* or GCF_*).
        email: Email for NCBI Entrez.
        api_key: Optional NCBI API key.
        store_sequence: If True, store full sequence for analysis.

    Returns:
        StrainMetadata object with extracted information.
    """
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    try:
        # Search for the assembly in the assembly database
        with Entrez.esearch(db="assembly", term=accession) as handle:
            search_results = Entrez.read(handle)

        if not search_results["IdList"]:
            print(f"No assembly found for {accession}", file=sys.stderr)
            return StrainMetadata(source_id=accession)

        assembly_id = search_results["IdList"][0]

        # Get assembly summary for linked nucleotide records
        with Entrez.esummary(db="assembly", id=assembly_id) as handle:
            summary = Entrez.read(handle)

        doc_sum = summary["DocumentSummarySet"]["DocumentSummary"][0]

        # Get linked nucleotide sequences
        with Entrez.elink(
            dbfrom="assembly",
            db="nucleotide",
            id=assembly_id,
            linkname="assembly_nuccore_refseq",
        ) as handle:
            link_results = Entrez.read(handle)

        # Try RefSeq first, then GenBank
        nuc_ids = []
        if link_results[0]["LinkSetDb"]:
            nuc_ids = [link["Id"] for link in link_results[0]["LinkSetDb"][0]["Link"]]
        else:
            # Try GenBank nucleotide link
            with Entrez.elink(
                dbfrom="assembly",
                db="nucleotide",
                id=assembly_id,
                linkname="assembly_nuccore_insdc",
            ) as handle:
                link_results = Entrez.read(handle)
            if link_results[0]["LinkSetDb"]:
                nuc_ids = [
                    link["Id"] for link in link_results[0]["LinkSetDb"][0]["Link"]
                ]

        if not nuc_ids:
            print(f"No nucleotide records found for {accession}", file=sys.stderr)
            # Return what we can from assembly summary
            metadata = StrainMetadata(source_id=accession)
            metadata.organism = doc_sum.get("Organism", "")
            metadata.biosample = doc_sum.get("BioSampleAccn", "")
            metadata.assembly_accession = accession
            return metadata

        # Limit to first 10 records to avoid timeout
        nuc_ids = nuc_ids[:10]

        # Fetch GenBank records
        with Entrez.efetch(
            db="nucleotide",
            id=",".join(nuc_ids),
            rettype="gbwithparts",
            retmode="text",
        ) as handle:
            records = list(SeqIO.parse(handle, "genbank"))

        metadata = extract_metadata_from_records(records, accession, store_sequence)
        metadata.assembly_accession = accession

        return metadata

    except Exception as e:
        print(f"Error fetching {accession}: {e}", file=sys.stderr)
        return StrainMetadata(source_id=accession)


def normalize_strain_name(name: str) -> str:
    """Normalize a strain name for comparison.

    Args:
        name: Original strain name.

    Returns:
        Normalized strain name (lowercase, standardized separators).
    """
    if not name:
        return ""

    # Convert to lowercase
    normalized = name.lower()

    # Replace common separators with single space
    normalized = re.sub(r"[-_/\\]+", " ", normalized)

    # Remove multiple spaces
    normalized = re.sub(r"\s+", " ", normalized)

    # Remove common prefixes that vary
    prefixes_to_remove = [
        "strain ",
        "str ",
        "str. ",
        "isolate ",
        "iso ",
        "atcc ",
        "dsmz ",
        "dsm ",
        "nrrl ",
        "jcm ",
        "nbrc ",
        "ncimb ",
        "cip ",
        "lmg ",
        "ccug ",
        "mtcc ",
        "kctc ",
        "cgmcc ",
    ]
    for prefix in prefixes_to_remove:
        if normalized.startswith(prefix):
            normalized = normalized[len(prefix) :]

    return normalized.strip()


def extract_culture_collection_id(cc_string: str) -> tuple[str, str]:
    """Extract culture collection name and ID.

    Args:
        cc_string: Culture collection string (e.g., "ATCC:BAA-965").

    Returns:
        Tuple of (collection_name, id).
    """
    if ":" in cc_string:
        parts = cc_string.split(":", 1)
        return parts[0].upper(), parts[1]
    return "", cc_string


def calculate_contig_profile_similarity(
    lengths1: list[int], lengths2: list[int]
) -> float:
    """Calculate similarity between two contig length profiles.

    Args:
        lengths1: List of contig lengths for first strain.
        lengths2: List of contig lengths for second strain.

    Returns:
        Similarity score (0-1).
    """
    if not lengths1 or not lengths2:
        return 0.0

    # Sort lengths for comparison
    sorted1 = sorted(lengths1, reverse=True)
    sorted2 = sorted(lengths2, reverse=True)

    # Check if number of contigs is similar
    if len(sorted1) != len(sorted2):
        len_ratio = min(len(sorted1), len(sorted2)) / max(len(sorted1), len(sorted2))
        if len_ratio < 0.8:
            return len_ratio * 0.5

    # Compare lengths pairwise
    min_len = min(len(sorted1), len(sorted2))
    similarities = []
    for i in range(min_len):
        l1, l2 = sorted1[i], sorted2[i]
        if l1 == 0 and l2 == 0:
            sim = 1.0
        elif l1 == 0 or l2 == 0:
            sim = 0.0
        else:
            sim = min(l1, l2) / max(l1, l2)
        similarities.append(sim)

    return sum(similarities) / len(similarities) if similarities else 0.0


def compare_strains(
    meta1: StrainMetadata,
    meta2: StrainMetadata,
    similarity_threshold: float = 0.99,
    ani_threshold: float = 0.999,
    enable_sequence_analysis: bool = False,
    kmer_size: int = 21,
) -> DuplicateMatch | None:
    """Compare two strains for potential duplicate status.

    Args:
        meta1: Metadata for first strain.
        meta2: Metadata for second strain.
        similarity_threshold: Minimum sequence similarity for duplicate detection.
        ani_threshold: Minimum ANI threshold for sequence-based duplicate detection.
        enable_sequence_analysis: Whether to use MinHash/k-mer sequence analysis.
        kmer_size: K-mer size for ANI estimation.

    Returns:
        DuplicateMatch if potential duplicate found, None otherwise.
    """
    match_reasons: list[str] = []
    confidence_scores: list[int] = []

    # Sequence analysis metrics (computed if enabled)
    minhash_sim = 0.0
    estimated_ani = 0.0
    shared_kmer_ratio = 0.0

    # 1. Same BioSample (definitive match)
    if meta1.biosample and meta2.biosample and meta1.biosample == meta2.biosample:
        match_reasons.append(f"Same BioSample: {meta1.biosample}")
        confidence_scores.append(100)

    # 2. Same culture collection ID
    if meta1.culture_collection and meta2.culture_collection:
        coll1, id1 = extract_culture_collection_id(meta1.culture_collection)
        coll2, id2 = extract_culture_collection_id(meta2.culture_collection)
        if id1 and id2 and id1 == id2:
            match_reasons.append(f"Same culture collection ID: {id1}")
            confidence_scores.append(90)

    # 3. Identical sequences (check checksums)
    if meta1.sequence_checksums and meta2.sequence_checksums:
        common_checksums = set(meta1.sequence_checksums) & set(meta2.sequence_checksums)
        if common_checksums:
            ratio = len(common_checksums) / max(
                len(meta1.sequence_checksums), len(meta2.sequence_checksums)
            )
            if ratio >= similarity_threshold:
                match_reasons.append(
                    f"Identical sequences: {len(common_checksums)} contigs match "
                    f"({ratio:.1%} overlap)"
                )
                confidence_scores.append(95)

    # 4. Sequence-based analysis (MinHash + k-mer)
    if enable_sequence_analysis and meta1.minhash_sketch and meta2.minhash_sketch:
        minhash_sim, estimated_ani, shared_kmer_ratio, kmer_corr = compare_sequences(
            meta1, meta2, kmer_size
        )

        # Very high ANI suggests duplicate (>99.9% = essentially same strain)
        if estimated_ani >= ani_threshold:
            match_reasons.append(
                f"High sequence similarity: estimated ANI {estimated_ani:.4f} "
                f"(MinHash Jaccard: {minhash_sim:.4f})"
            )
            confidence_scores.append(95)
        elif estimated_ani >= 0.995:
            match_reasons.append(
                f"Very similar sequences: estimated ANI {estimated_ani:.4f} "
                f"(MinHash Jaccard: {minhash_sim:.4f})"
            )
            confidence_scores.append(85)
        elif estimated_ani >= 0.99:
            match_reasons.append(
                f"Similar sequences: estimated ANI {estimated_ani:.4f} "
                f"(MinHash Jaccard: {minhash_sim:.4f})"
            )
            confidence_scores.append(70)

        # High k-mer correlation with moderate ANI
        if kmer_corr >= 0.95 and estimated_ani >= 0.98:
            match_reasons.append(
                f"High k-mer frequency correlation: {kmer_corr:.4f} "
                f"(shared k-mers: {shared_kmer_ratio:.2%})"
            )
            confidence_scores.append(80)

    # 5. Same organism + similar strain names
    if meta1.organism and meta2.organism:
        if meta1.organism.lower() == meta2.organism.lower():
            norm_strain1 = normalize_strain_name(meta1.strain or meta1.isolate)
            norm_strain2 = normalize_strain_name(meta2.strain or meta2.isolate)

            if norm_strain1 and norm_strain2:
                # Check for exact match after normalization
                if norm_strain1 == norm_strain2:
                    match_reasons.append(
                        f"Same organism ({meta1.organism}) with matching strain names: "
                        f"'{meta1.strain or meta1.isolate}' ~ '{meta2.strain or meta2.isolate}'"
                    )
                    confidence_scores.append(85)
                # Check if one contains the other
                elif norm_strain1 in norm_strain2 or norm_strain2 in norm_strain1:
                    match_reasons.append(
                        f"Same organism ({meta1.organism}) with similar strain names: "
                        f"'{meta1.strain or meta1.isolate}' ~ '{meta2.strain or meta2.isolate}'"
                    )
                    confidence_scores.append(70)

    # 6. Very similar genome statistics
    if meta1.total_length > 0 and meta2.total_length > 0:
        length_ratio = min(meta1.total_length, meta2.total_length) / max(
            meta1.total_length, meta2.total_length
        )
        gc_diff = abs(meta1.gc_content - meta2.gc_content)

        if length_ratio >= 0.99 and gc_diff <= 0.001:
            # Check contig profile
            contig_sim = calculate_contig_profile_similarity(
                meta1.contig_lengths, meta2.contig_lengths
            )
            if contig_sim >= 0.95:
                match_reasons.append(
                    f"Near-identical genome profile: length ratio {length_ratio:.4f}, "
                    f"GC diff {gc_diff:.4f}, contig similarity {contig_sim:.2%}"
                )
                confidence_scores.append(80)

    # 7. Same type material designation
    if meta1.type_material and meta2.type_material:
        if meta1.type_material == meta2.type_material:
            match_reasons.append(f"Same type material: {meta1.type_material}")
            confidence_scores.append(85)

    # Only return match if we have reasons
    if not match_reasons:
        return None

    # Determine confidence level
    max_score = max(confidence_scores)
    if max_score >= 90:
        confidence = "high"
    elif max_score >= 70:
        confidence = "medium"
    else:
        confidence = "low"

    return DuplicateMatch(
        strain1=meta1.source_id,
        strain2=meta2.source_id,
        match_reasons=match_reasons,
        confidence=confidence,
        metadata1=meta1,
        metadata2=meta2,
        minhash_similarity=minhash_sim,
        estimated_ani=estimated_ani,
        shared_kmers_ratio=shared_kmer_ratio,
    )


def write_report(
    matches: list[DuplicateMatch],
    metadata_list: list[StrainMetadata],
    output: Path,
    include_sequence_analysis: bool = False,
):
    """Write the duplicate detection report.

    Args:
        matches: List of potential duplicate matches.
        metadata_list: List of all strain metadata.
        output: Output file path.
        include_sequence_analysis: Whether to include sequence analysis columns.
    """
    # Write matches report
    with output.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")

        # Header
        header = [
            "strain_1",
            "strain_2",
            "confidence",
            "match_reasons",
            "organism_1",
            "organism_2",
            "strain_name_1",
            "strain_name_2",
            "biosample_1",
            "biosample_2",
            "culture_collection_1",
            "culture_collection_2",
            "total_length_1",
            "total_length_2",
            "gc_content_1",
            "gc_content_2",
            "num_contigs_1",
            "num_contigs_2",
        ]
        if include_sequence_analysis:
            header.extend(
                [
                    "minhash_similarity",
                    "estimated_ani",
                    "shared_kmers_ratio",
                ]
            )
        writer.writerow(header)

        # Sort by confidence
        confidence_order = {"high": 0, "medium": 1, "low": 2}
        sorted_matches = sorted(matches, key=lambda m: confidence_order[m.confidence])

        for match in sorted_matches:
            m1, m2 = match.metadata1, match.metadata2
            row = [
                match.strain1,
                match.strain2,
                match.confidence,
                " | ".join(match.match_reasons),
                m1.organism,
                m2.organism,
                m1.strain or m1.isolate,
                m2.strain or m2.isolate,
                m1.biosample,
                m2.biosample,
                m1.culture_collection,
                m2.culture_collection,
                m1.total_length,
                m2.total_length,
                f"{m1.gc_content:.4f}",
                f"{m2.gc_content:.4f}",
                m1.num_contigs,
                m2.num_contigs,
            ]
            if include_sequence_analysis:
                row.extend(
                    [
                        f"{match.minhash_similarity:.4f}",
                        f"{match.estimated_ani:.4f}",
                        f"{match.shared_kmers_ratio:.4f}",
                    ]
                )
            writer.writerow(row)

    # Write metadata summary
    metadata_output = output.with_stem(output.stem + "_metadata")
    with metadata_output.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")

        writer.writerow(
            [
                "source_id",
                "organism",
                "strain",
                "isolate",
                "culture_collection",
                "type_material",
                "biosample",
                "bioproject",
                "assembly_accession",
                "total_length",
                "num_contigs",
                "gc_content",
                "isolation_source",
                "host",
                "country",
                "collection_date",
            ]
        )

        for meta in sorted(metadata_list, key=lambda m: m.organism):
            writer.writerow(
                [
                    meta.source_id,
                    meta.organism,
                    meta.strain,
                    meta.isolate,
                    meta.culture_collection,
                    meta.type_material,
                    meta.biosample,
                    meta.bioproject,
                    meta.assembly_accession,
                    meta.total_length,
                    meta.num_contigs,
                    f"{meta.gc_content:.4f}",
                    meta.isolation_source,
                    meta.host,
                    meta.country,
                    meta.collection_date,
                ]
            )

    print(f"Duplicate report written to: {output}")
    print(f"Metadata summary written to: {metadata_output}")


def write_cluster_report(
    clusters: dict[str, list[str]],
    matches: list[DuplicateMatch],
    metadata_lookup: dict[str, StrainMetadata],
    output: Path,
):
    """Write a cluster-based report grouping connected duplicate strains.

    Args:
        clusters: Mapping from cluster root to list of member source_ids.
        matches: List of pairwise duplicate matches.
        metadata_lookup: Mapping from source_id to StrainMetadata.
        output: Base output file path (cluster report uses a '_clusters' suffix).
    """
    cluster_output = output.with_stem(output.stem + "_clusters")

    # Build quick lookup: which matches involve a given source_id
    match_by_pair: dict[tuple[str, str], DuplicateMatch] = {}
    for match in matches:
        match_by_pair[(match.strain1, match.strain2)] = match
        match_by_pair[(match.strain2, match.strain1)] = match

    with cluster_output.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")

        writer.writerow([
            "cluster_id",
            "cluster_size",
            "source_id",
            "organism",
            "strain",
            "isolate",
            "culture_collection",
            "biosample",
            "assembly_accession",
            "total_length",
            "gc_content",
            "num_contigs",
            "linked_to",
            "link_confidence",
            "link_reasons",
        ])

        for cluster_idx, (root, members) in enumerate(
            sorted(clusters.items(), key=lambda kv: -len(kv[1])), start=1
        ):
            for member in sorted(members):
                meta = metadata_lookup.get(member)
                if meta is None:
                    continue

                # Find which other member in the cluster linked this one
                linked_to = ""
                link_confidence = ""
                link_reasons = ""
                for other in members:
                    if other == member:
                        continue
                    pair = (member, other)
                    if pair in match_by_pair:
                        m = match_by_pair[pair]
                        linked_to = other
                        link_confidence = m.confidence
                        link_reasons = " | ".join(m.match_reasons)
                        break  # Only need the first direct link

                writer.writerow([
                    cluster_idx,
                    len(members),
                    member,
                    meta.organism,
                    meta.strain,
                    meta.isolate,
                    meta.culture_collection,
                    meta.biosample,
                    meta.assembly_accession,
                    meta.total_length,
                    f"{meta.gc_content:.4f}",
                    meta.num_contigs,
                    linked_to,
                    link_confidence,
                    link_reasons,
                ])

    print(f"Cluster report written to: {cluster_output}")


def is_assembly_accession(text: str) -> bool:
    """Check if text looks like a GenBank assembly accession.

    Args:
        text: Text to check.

    Returns:
        True if it matches assembly accession pattern.
    """
    return bool(re.match(r"^GC[AF]_\d+(\.\d+)?$", text.strip()))


def main(
    input_path: Path,
    output: Path,
    email: str | None = None,
    api_key: str | None = None,
    threads: int | None = None,
    similarity_threshold: float = 0.99,
    enable_sequence_analysis: bool = False,
    kmer_size: int = 21,
    sketch_size: int = 1000,
    ani_threshold: float = 0.999,
):
    """Main function to detect duplicate strains.

    Args:
        input_path: Directory with GBFF files or file with assembly IDs.
        output: Output report file path.
        email: Email for NCBI Entrez.
        api_key: Optional NCBI API key.
        threads: Number of threads for parallel processing.
        similarity_threshold: Minimum similarity threshold.
        enable_sequence_analysis: Enable MinHash/k-mer sequence analysis.
        kmer_size: K-mer size for sequence analysis.
        sketch_size: MinHash sketch size.
        ani_threshold: ANI threshold for duplicate detection.
    """
    metadata_list: list[StrainMetadata] = []

    # Determine input type
    if input_path.is_dir():
        # Process GBFF files from directory
        gbff_files = (
            list(input_path.glob("*.gbff"))
            + list(input_path.glob("*.gbff.gz"))
            + list(input_path.glob("*.gb"))
            + list(input_path.glob("*.gb.gz"))
            + list(input_path.glob("*.gbk"))
            + list(input_path.glob("*.gbk.gz"))
        )

        if not gbff_files:
            print(f"No GenBank files found in {input_path}", file=sys.stderr)
            sys.exit(1)

        print(f"Processing {len(gbff_files)} GenBank files...")
        if enable_sequence_analysis:
            print("  (with sequence analysis enabled - this may take longer)")

        with ThreadPoolExecutor(max_workers=threads) as executor:
            futures = {
                executor.submit(parse_gbff_file, path, enable_sequence_analysis): path
                for path in gbff_files
            }

            for future in as_completed(futures):
                result = future.result()
                if result.organism:  # Only add if we got valid metadata
                    metadata_list.append(result)

    elif input_path.is_file():
        # Read assembly accessions from file
        with input_path.open() as f:
            lines = [line.strip() for line in f if line.strip()]

        # Check if lines are assembly accessions
        accessions = [line for line in lines if is_assembly_accession(line)]

        if accessions:
            if not email:
                print(
                    "Email is required for fetching from NCBI. Use --email option.",
                    file=sys.stderr,
                )
                sys.exit(1)

            print(f"Fetching metadata for {len(accessions)} assembly accessions...")
            if enable_sequence_analysis:
                print("  (with sequence analysis enabled - this may take longer)")
            Entrez.email = email
            if api_key:
                Entrez.api_key = api_key

            with ThreadPoolExecutor(max_workers=min(threads or 4, 4)) as executor:
                futures = {
                    executor.submit(
                        fetch_assembly_metadata,
                        acc,
                        email,
                        api_key,
                        enable_sequence_analysis,
                    ): acc
                    for acc in accessions
                }

                for future in as_completed(futures):
                    result = future.result()
                    if result.organism:
                        metadata_list.append(result)

        else:
            # Maybe the file contains paths to GBFF files
            gbff_files = [Path(line) for line in lines if Path(line).exists()]
            if gbff_files:
                print(f"Processing {len(gbff_files)} GenBank files from list...")
                if enable_sequence_analysis:
                    print("  (with sequence analysis enabled - this may take longer)")

                with ThreadPoolExecutor(max_workers=threads) as executor:
                    futures = {
                        executor.submit(
                            parse_gbff_file, path, enable_sequence_analysis
                        ): path
                        for path in gbff_files
                    }

                    for future in as_completed(futures):
                        result = future.result()
                        if result.organism:
                            metadata_list.append(result)
            else:
                print(
                    f"Could not parse input file: {input_path}. "
                    "Expected assembly accessions (GCA_*/GCF_*) or paths to GBFF files.",
                    file=sys.stderr,
                )
                sys.exit(1)

    else:
        print(f"Input path does not exist: {input_path}", file=sys.stderr)
        sys.exit(1)

    if len(metadata_list) < 2:
        print("Need at least 2 valid records to compare.", file=sys.stderr)
        sys.exit(1)

    print(f"Successfully loaded {len(metadata_list)} records.")

    # Compute sequence analysis if enabled
    if enable_sequence_analysis:
        print(
            f"Computing MinHash sketches (k={kmer_size}, sketch_size={sketch_size})..."
        )
        for i, meta in enumerate(metadata_list):
            if (i + 1) % 10 == 0:
                print(f"  Progress: {i + 1}/{len(metadata_list)} genomes analyzed")
            compute_sequence_analysis(meta, kmer_size, sketch_size)
        # Clear full sequences to free memory after analysis
        for meta in metadata_list:
            meta.full_sequence = ""

    print("Comparing strains for potential duplicates...")

    # Use Union-Find to cluster connected strains and skip redundant comparisons
    uf = UnionFind()
    matches: list[DuplicateMatch] = []
    total_comparisons = len(metadata_list) * (len(metadata_list) - 1) // 2
    comparison_count = 0
    skipped_count = 0

    for meta1, meta2 in combinations(metadata_list, 2):
        comparison_count += 1
        if comparison_count % 1000 == 0:
            print(
                f"  Progress: {comparison_count}/{total_comparisons} comparisons "
                f"({skipped_count} skipped as already connected)"
            )

        # Skip if both strains already belong to the same cluster
        if uf.connected(meta1.source_id, meta2.source_id):
            skipped_count += 1
            continue

        match = compare_strains(
            meta1,
            meta2,
            similarity_threshold,
            ani_threshold,
            enable_sequence_analysis,
            kmer_size,
        )
        if match:
            matches.append(match)
            # Merge their clusters so transitive duplicates are skipped
            uf.union(meta1.source_id, meta2.source_id)

    # Build cluster information
    clusters = uf.get_clusters()
    # Only keep clusters with more than one member
    duplicate_clusters = {k: v for k, v in clusters.items() if len(v) > 1}

    # Report results
    print(f"\nFound {len(matches)} unique duplicate pairs ({skipped_count} redundant comparisons skipped):")
    high_conf = sum(1 for m in matches if m.confidence == "high")
    med_conf = sum(1 for m in matches if m.confidence == "medium")
    low_conf = sum(1 for m in matches if m.confidence == "low")
    print(f"  High confidence: {high_conf}")
    print(f"  Medium confidence: {med_conf}")
    print(f"  Low confidence: {low_conf}")
    print(f"  Duplicate clusters: {len(duplicate_clusters)}")

    # Build metadata lookup for cluster report
    metadata_lookup = {meta.source_id: meta for meta in metadata_list}

    write_report(matches, metadata_list, output, enable_sequence_analysis)
    write_cluster_report(duplicate_clusters, matches, metadata_lookup, output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Detect potential duplicate strains from GBFF files or GenBank assembly IDs"
    )

    parser.add_argument(
        "input",
        type=Path,
        help="Directory containing GBFF files, or file with assembly IDs (one per line)",
    )
    parser.add_argument(
        "--email",
        type=str,
        help="Email address for NCBI Entrez (required for fetching)",
    )
    parser.add_argument(
        "--api-key",
        type=str,
        help="NCBI API key for higher rate limits",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Number of threads for parallel processing",
    )
    parser.add_argument(
        "--similarity",
        type=float,
        default=0.99,
        help="Minimum sequence similarity threshold (0-1) for duplicate detection (default: 0.99)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("duplicate_report.tsv"),
        help="Output report file path (default: duplicate_report.tsv)",
    )
    parser.add_argument(
        "--sequence-analysis",
        action="store_true",
        help="Enable deep sequence analysis using MinHash sketching and k-mer comparison",
    )
    parser.add_argument(
        "--kmer-size",
        type=int,
        default=21,
        help="K-mer size for sequence analysis (default: 21)",
    )
    parser.add_argument(
        "--sketch-size",
        type=int,
        default=1000,
        help="MinHash sketch size for fast similarity estimation (default: 1000)",
    )
    parser.add_argument(
        "--ani-threshold",
        type=float,
        default=0.999,
        help="ANI threshold for considering genomes as duplicates (default: 0.999)",
    )

    args = parser.parse_args()

    main(
        args.input,
        args.output,
        email=args.email,
        api_key=args.api_key,
        threads=args.threads,
        similarity_threshold=args.similarity,
        enable_sequence_analysis=args.sequence_analysis,
        kmer_size=args.kmer_size,
        sketch_size=args.sketch_size,
        ani_threshold=args.ani_threshold,
    )
