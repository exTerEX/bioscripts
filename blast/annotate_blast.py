## Annotate BLAST XML/XML2 results with gene-level information from GenBank records
#
## Usage:
#   $ python blast/annotate_blast.py -h
#   usage: annotate_blast.py [-h] [--gbff-dir GBFF_DIR] [--email EMAIL]
#                            [--api-key API_KEY] [--threads THREADS]
#                            [--strand STRAND] [--output OUTPUT]
#                            blast_xml
#
#   Parse BLAST XML/XML2 results and annotate hits with gene information
#   from GenBank records. Supports all BLAST program types (blastn, blastp,
#   blastx, tblastn, tblastx) and auto-detects the program from the XML.
#
#   For each HSP, annotation is extracted for both the hit (subject) and the
#   query sequence.  When --gbff-dir is provided, local files are tried first;
#   accessions not found locally are automatically fetched from NCBI as a
#   fallback.  A record_source column indicates where each record came from
#   ("local", "ncbi", or "not_found").
#
#   positional arguments:
#     blast_xml             Path to BLAST XML (outfmt 5) or XML2 (outfmt 16) file
#
#   options:
#     -h, --help            show this help message and exit
#     --gbff-dir GBFF_DIR   Directory containing GBFF/GBK files for local lookup
#     --email EMAIL         Email address for NCBI Entrez
#     --api-key API_KEY     NCBI API key for higher rate limits
#     --threads THREADS     Number of threads for parallel processing.
#                           Defaults to number of CPUs.
#     --strand STRAND       Filter by strand/frame: -3,-2,-1,1,2,3 for specific
#                           frame, or -/+ for all minus/plus strand frames
#     --output OUTPUT       Output TSV file path (default: stdout)

from __future__ import annotations

import argparse
import re
import sys
import time
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import pandas as pd
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIXML
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

BLAST_PROGRAMS: dict[str, dict[str, str]] = {
    "blastn": {"query_type": "nucleotide", "subject_type": "nucleotide"},
    "blastp": {"query_type": "protein", "subject_type": "protein"},
    "blastx": {"query_type": "nucleotide", "subject_type": "protein"},
    "tblastn": {"query_type": "protein", "subject_type": "nucleotide"},
    "tblastx": {"query_type": "nucleotide", "subject_type": "nucleotide"},
}

ANNOTATION_FEATURE_TYPES = {"CDS", "rRNA", "tRNA", "ncRNA", "tmRNA", "misc_RNA", "gene"}

FEATURE_TYPE_PRIORITY: dict[str, int] = {
    "CDS": 0,
    "rRNA": 1,
    "tRNA": 1,
    "ncRNA": 1,
    "tmRNA": 1,
    "misc_RNA": 2,
    "gene": 3,
}

NCBI_RETRY_ATTEMPTS = 3
NCBI_RETRY_DELAY = 2  # seconds

# Regex for common NCBI accession patterns (RefSeq + GenBank).
# Uses lookaround instead of \b to correctly capture versions in strings
# like "NC_000913.3_cds" where _ after the version defeats \b.
_ACCESSION_RE = re.compile(
    r"(?<![a-zA-Z0-9])("
    # RefSeq nucleotide / protein
    r"[ANXYW][CGMPRSTWZ]_\d{6,9}(?:\.\d+)?"
    r"|"
    # GenBank nucleotide  (1-2 letters + 5-8 digits)
    r"[A-Z]{1,2}\d{5,8}(?:\.\d+)?"
    r"|"
    # GenBank protein (3 letters + 5 digits)
    r"[A-Z]{3}\d{5}(?:\.\d+)?"
    r")(?![a-zA-Z0-9])"
)


def extract_accession(text: str) -> str | None:
    """Extract the first NCBI-like accession from *text*.

    Handles common BLAST identifiers such as ``lcl|…``, ``ref|…|``,
    ``gi|…|ref|…|``, and bare accession strings.

    Args:
        text: Query ID, description, or other text containing an accession.

    Returns:
        The accession string, or ``None`` if nothing recognisable is found.
    """
    if not text:
        return None
    m = _ACCESSION_RE.search(text)
    return m.group(1) if m else None


# ---------------------------------------------------------------------------
# BLAST Format Detection
# ---------------------------------------------------------------------------


def detect_blast_format(path: Path) -> str:
    """Detect whether file is BLAST XML (outfmt 5) or XML2 (outfmt 16).

    Args:
        path: Path to the BLAST output file.

    Returns:
        ``'xml'`` for BLAST XML or ``'xml2'`` for BLAST XML2.

    Raises:
        ValueError: If the format cannot be determined.
    """
    with path.open() as f:
        header = f.read(4096)

    if "BlastXML2" in header:
        return "xml2"
    if "BlastOutput" in header:
        return "xml"
    raise ValueError(f"Unrecognized BLAST output format in: {path}")


# ---------------------------------------------------------------------------
# BLAST XML (outfmt 5) Parser
# ---------------------------------------------------------------------------


def parse_blast_xml(path: Path) -> tuple[list[dict[str, Any]], str]:
    """Parse BLAST XML (outfmt 5) file into a list of HSP dictionaries.

    Args:
        path: Path to the BLAST XML file.

    Returns:
        Tuple of (list of HSP dictionaries, BLAST program name).
    """
    rows: list[dict[str, Any]] = []
    blast_program = ""

    with path.open() as f:
        for record in NCBIXML.parse(f):
            if not blast_program:
                blast_program = record.application.lower()

            query_id = getattr(record, "query_id", "") or ""
            query_desc = getattr(record, "query", "") or ""
            query_len = getattr(record, "query_length", 0) or 0

            for alignment in record.alignments:
                hit_accession = alignment.accession
                hit_desc = alignment.hit_def
                hit_len = alignment.length

                for hsp_idx, hsp in enumerate(alignment.hsps, 1):
                    query_frame, hit_frame = 0, 0
                    frame = getattr(hsp, "frame", None)
                    if frame:
                        if isinstance(frame, tuple):
                            query_frame = frame[0] if len(frame) > 0 else 0
                            hit_frame = frame[1] if len(frame) > 1 else 0
                        else:
                            query_frame = int(frame)

                    align_len = getattr(hsp, "align_length", 0) or len(hsp.query)
                    identities = getattr(hsp, "identities", 0) or 0
                    positives = getattr(hsp, "positives", None)
                    gaps = getattr(hsp, "gaps", 0) or 0

                    rows.append(
                        {
                            "query_id": query_id,
                            "query_description": query_desc,
                            "query_length": query_len,
                            "hit_accession": hit_accession,
                            "hit_description": hit_desc,
                            "hit_length": hit_len,
                            "hsp_num": hsp_idx,
                            "score": hsp.score,
                            "bit_score": hsp.bits,
                            "evalue": hsp.expect,
                            "identity": identities,
                            "identity_pct": round(identities / align_len * 100, 2) if align_len else 0.0,
                            "positives": positives,
                            "positives_pct": (
                                round(positives / align_len * 100, 2) if positives and align_len else None
                            ),
                            "gaps": gaps,
                            "alignment_length": align_len,
                            "query_start": hsp.query_start,
                            "query_end": hsp.query_end,
                            "query_frame": query_frame,
                            "hit_start": hsp.sbjct_start,
                            "hit_end": hsp.sbjct_end,
                            "hit_frame": hit_frame,
                            "query_seq": hsp.query,
                            "hit_seq": hsp.sbjct,
                            "midline": hsp.match,
                        }
                    )

    return rows, blast_program


# ---------------------------------------------------------------------------
# BLAST XML2 (outfmt 16) Parser
# ---------------------------------------------------------------------------


def _xml2_find(element: ET.Element, tag: str, ns: str) -> ET.Element | None:
    """Find a child element, with or without XML namespace."""
    return element.find(f"{{{ns}}}{tag}") if ns else element.find(tag)


def _xml2_findall(element: ET.Element, tag: str, ns: str) -> list[ET.Element]:
    """Find all child elements, with or without XML namespace."""
    return element.findall(f"{{{ns}}}{tag}") if ns else element.findall(tag)


def _xml2_text(element: ET.Element, tag: str, ns: str, default: str = "") -> str:
    """Get text content of a child element."""
    child = _xml2_find(element, tag, ns)
    return child.text if child is not None and child.text else default


def parse_blast_xml2(path: Path) -> tuple[list[dict[str, Any]], str]:
    """Parse BLAST XML2 (outfmt 16) file into a list of HSP dictionaries.

    Args:
        path: Path to the BLAST XML2 file.

    Returns:
        Tuple of (list of HSP dictionaries, BLAST program name).
    """
    tree = ET.parse(path)
    root = tree.getroot()

    # Detect namespace
    ns = ""
    if root.tag.startswith("{"):
        ns = root.tag.split("}")[0].lstrip("{")

    rows: list[dict[str, Any]] = []
    blast_program = ""

    # Handle both <BlastXML2> root and direct <BlastOutput2> root
    if root.tag.endswith("BlastXML2"):
        bo2_elements = _xml2_findall(root, "BlastOutput2", ns)
    elif root.tag.endswith("BlastOutput2"):
        bo2_elements = [root]
    else:
        raise ValueError(f"Unexpected root element in XML2: {root.tag}")

    for bo2 in bo2_elements:
        report_el = _xml2_find(bo2, "report", ns)
        if report_el is None:
            continue
        report = _xml2_find(report_el, "Report", ns)
        if report is None:
            continue

        if not blast_program:
            blast_program = _xml2_text(report, "program", ns).lower()

        results_el = _xml2_find(report, "results", ns)
        if results_el is None:
            continue
        results = _xml2_find(results_el, "Results", ns)
        if results is None:
            continue

        for search_el in _xml2_findall(results, "search", ns):
            search = _xml2_find(search_el, "Search", ns)
            if search is None:
                continue

            query_id = _xml2_text(search, "query-id", ns)
            query_title = _xml2_text(search, "query-title", ns)
            query_len = int(_xml2_text(search, "query-len", ns, "0"))

            hits_el = _xml2_find(search, "hits", ns)
            if hits_el is None:
                continue

            for hit in _xml2_findall(hits_el, "Hit", ns):
                hit_len = int(_xml2_text(hit, "len", ns, "0"))

                # Parse hit description(s) – take the first HitDescr
                desc_el = _xml2_find(hit, "description", ns)
                hit_accession = ""
                hit_desc = ""
                if desc_el is not None:
                    for hd in _xml2_findall(desc_el, "HitDescr", ns):
                        if not hit_accession:
                            hit_accession = _xml2_text(hd, "accession", ns)
                        if not hit_desc:
                            hit_desc = _xml2_text(hd, "title", ns)

                hsps_el = _xml2_find(hit, "hsps", ns)
                if hsps_el is None:
                    continue

                for hsp_idx, hsp in enumerate(_xml2_findall(hsps_el, "Hsp", ns), 1):
                    align_len = int(_xml2_text(hsp, "align-len", ns, "0"))
                    identities = int(_xml2_text(hsp, "identity", ns, "0"))
                    positives_val = _xml2_text(hsp, "positive", ns)
                    positives = int(positives_val) if positives_val else None
                    gaps = int(_xml2_text(hsp, "gaps", ns, "0"))

                    query_frame = int(_xml2_text(hsp, "query-frame", ns, "0"))
                    hit_frame = int(_xml2_text(hsp, "hit-frame", ns, "0"))

                    # If frame is 0, derive from strand element (blastn)
                    if query_frame == 0:
                        q_strand = _xml2_text(hsp, "query-strand", ns).lower()
                        if q_strand == "minus":
                            query_frame = -1
                        elif q_strand == "plus":
                            query_frame = 1

                    if hit_frame == 0:
                        h_strand = _xml2_text(hsp, "hit-strand", ns).lower()
                        if h_strand == "minus":
                            hit_frame = -1
                        elif h_strand == "plus":
                            hit_frame = 1

                    rows.append(
                        {
                            "query_id": query_id,
                            "query_description": query_title,
                            "query_length": query_len,
                            "hit_accession": hit_accession,
                            "hit_description": hit_desc,
                            "hit_length": hit_len,
                            "hsp_num": hsp_idx,
                            "score": int(_xml2_text(hsp, "score", ns, "0")),
                            "bit_score": float(_xml2_text(hsp, "bit-score", ns, "0")),
                            "evalue": float(_xml2_text(hsp, "evalue", ns, "0")),
                            "identity": identities,
                            "identity_pct": round(identities / align_len * 100, 2) if align_len else 0.0,
                            "positives": positives,
                            "positives_pct": (
                                round(positives / align_len * 100, 2) if positives and align_len else None
                            ),
                            "gaps": gaps,
                            "alignment_length": align_len,
                            "query_start": int(_xml2_text(hsp, "query-from", ns, "0")),
                            "query_end": int(_xml2_text(hsp, "query-to", ns, "0")),
                            "query_frame": query_frame,
                            "hit_start": int(_xml2_text(hsp, "hit-from", ns, "0")),
                            "hit_end": int(_xml2_text(hsp, "hit-to", ns, "0")),
                            "hit_frame": hit_frame,
                            "query_seq": _xml2_text(hsp, "qseq", ns),
                            "hit_seq": _xml2_text(hsp, "hseq", ns),
                            "midline": _xml2_text(hsp, "midline", ns),
                        }
                    )

    return rows, blast_program


# ---------------------------------------------------------------------------
# Unified BLAST Parsing
# ---------------------------------------------------------------------------


def parse_blast_results(path: Path) -> tuple[pd.DataFrame, str]:
    """Parse a BLAST output file (XML or XML2) into a DataFrame.

    Auto-detects the format, extracts all HSP data, and returns the BLAST
    program name discovered in the file header.

    Args:
        path: Path to the BLAST output file.

    Returns:
        Tuple of (DataFrame with all HSP data, BLAST program name).
    """
    fmt = detect_blast_format(path)

    if fmt == "xml":
        rows, program = parse_blast_xml(path)
    else:
        rows, program = parse_blast_xml2(path)

    if not rows:
        print("Warning: No BLAST results found in the input file.", file=sys.stderr)
        return pd.DataFrame(), program

    df = pd.DataFrame(rows)
    df.insert(0, "blast_program", program)

    print(f"Parsed {len(df)} HSPs from {fmt.upper()} file (program: {program})", file=sys.stderr)
    return df, program


# ---------------------------------------------------------------------------
# Strand / Frame Filtering
# ---------------------------------------------------------------------------


def filter_by_strand(df: pd.DataFrame, strand_filter: str, blast_program: str) -> pd.DataFrame:
    """Filter a BLAST results DataFrame by strand or reading frame.

    For translated searches (tblastn, blastx, tblastx), the filter applies to
    the translated reading frame column.  For blastn it applies to the hit
    strand.  Filtering is not supported for blastp and a warning is emitted.

    Args:
        df: DataFrame with BLAST results.
        strand_filter: One of ``'-3', '-2', '-1', '1', '2', '3', '-', '+'``.
        blast_program: The BLAST program used.

    Returns:
        Filtered DataFrame.
    """
    if df.empty:
        return df

    # Decide which column to filter on
    if blast_program in ("tblastn", "tblastx"):
        frame_col = "hit_frame"
    elif blast_program == "blastx":
        frame_col = "query_frame"
    elif blast_program == "blastn":
        frame_col = "hit_frame"
    else:
        print(
            f"Warning: Strand filtering is not applicable for {blast_program}.",
            file=sys.stderr,
        )
        return df

    if strand_filter == "+":
        mask = df[frame_col] > 0
    elif strand_filter == "-":
        mask = df[frame_col] < 0
    else:
        target_frame = int(strand_filter)
        mask = df[frame_col] == target_frame

    filtered = df[mask].copy()
    print(
        f"Strand filter '{strand_filter}' on {frame_col}: {len(df)} -> {len(filtered)} HSPs",
        file=sys.stderr,
    )
    return filtered


# ---------------------------------------------------------------------------
# Local GBFF Index
# ---------------------------------------------------------------------------


def build_local_gbff_index(directory: Path) -> dict[str, Path]:
    """Build an index mapping identifiers to GBFF file paths.

    Parses every GenBank-format file in *directory* and indexes each record
    by its accession (with and without version), the filename stem, and any
    BioSample / Assembly cross-references found in the record metadata.

    Args:
        directory: Directory containing GBFF / GBK files.

    Returns:
        Dictionary mapping identifier strings to file paths.
    """
    index: dict[str, Path] = {}
    gbff_files = list(directory.glob("*.gbff")) + list(directory.glob("*.gbk")) + list(directory.glob("*.gb"))

    if not gbff_files:
        print(f"Warning: No GBFF/GBK files found in {directory}", file=sys.stderr)
        return index

    print(f"Indexing {len(gbff_files)} GenBank file(s) from {directory}...", file=sys.stderr)

    for gbff_path in gbff_files:
        # Index by filename stem (e.g. GCF_000123456.1)
        index[gbff_path.stem] = gbff_path

        try:
            for record in SeqIO.parse(gbff_path, "genbank"):
                acc = record.id
                index[acc] = gbff_path

                # Without version
                if "." in acc:
                    index[acc.rsplit(".", 1)[0]] = gbff_path

                # By record name
                if record.name and record.name != acc:
                    index[record.name] = gbff_path

                # By dbxrefs (BioSample, Assembly, BioProject)
                for xref in getattr(record, "dbxrefs", []):
                    if ":" in xref:
                        db, value = xref.split(":", 1)
                        if db.strip().lower() in ("biosample", "assembly", "bioproject"):
                            index[value.strip()] = gbff_path
        except Exception as e:
            print(f"Warning: Could not parse {gbff_path}: {e}", file=sys.stderr)

    print(f"Indexed {len(index)} identifier(s) across {len(gbff_files)} file(s)", file=sys.stderr)
    return index


def find_local_gbff(accession: str, gbff_index: dict[str, Path]) -> Path | None:
    """Look up a local GBFF file for the given accession.

    Tries, in order: exact match, without version suffix, substring match.

    Args:
        accession: Sequence accession or identifier.
        gbff_index: Index produced by :func:`build_local_gbff_index`.

    Returns:
        Path to a local GBFF file, or ``None``.
    """
    if accession in gbff_index:
        return gbff_index[accession]

    acc_no_ver = accession.rsplit(".", 1)[0] if "." in accession else accession
    if acc_no_ver in gbff_index:
        return gbff_index[acc_no_ver]

    # Substring match as last resort
    for key, path in gbff_index.items():
        if accession in key or acc_no_ver in key:
            return path

    return None


# ---------------------------------------------------------------------------
# GenBank Record Fetching
# ---------------------------------------------------------------------------


def fetch_genbank_record(accession: str, db: str = "nucleotide") -> list[SeqRecord]:  # type: ignore
    """Fetch GenBank record(s) from NCBI Entrez with retries.

    Args:
        accession: NCBI accession to fetch.
        db: Entrez database (``'nucleotide'`` or ``'protein'``).

    Returns:
        List of SeqRecord objects, or empty list on failure.
    """
    rettype = "gbwithparts" if db == "nucleotide" else "gp"

    for attempt in range(NCBI_RETRY_ATTEMPTS):
        try:
            with Entrez.efetch(db=db, id=accession, rettype=rettype, retmode="text") as handle:
                records = list(SeqIO.parse(handle, "genbank"))
            return records
        except Exception as e:
            if attempt < NCBI_RETRY_ATTEMPTS - 1:
                time.sleep(NCBI_RETRY_DELAY * (attempt + 1))
            else:
                print(f"Error fetching {accession} from NCBI ({db}): {e}", file=sys.stderr)
                return []


def _load_local_records(
    accessions: list[str],
    gbff_index: dict[str, Path],
    record_cache: dict[str, list[SeqRecord]],
    record_source: dict[str, str],
) -> list[str]:
    """Load GenBank records from local GBFF files into the cache.

    Each file is parsed at most once even if multiple accessions reference it.

    Args:
        accessions: Accessions to resolve locally.
        gbff_index: Local GBFF index.
        record_cache: Cache to populate (maps accession -> records list).
        record_source: Tracks where each accession was loaded from.

    Returns:
        Accessions that could **not** be resolved locally.
    """
    loaded_files: dict[Path, list[SeqRecord]] = {}
    not_found: list[str] = []

    for acc in accessions:
        if acc in record_cache:
            continue

        path = find_local_gbff(acc, gbff_index)
        if path is None:
            not_found.append(acc)
            continue

        if path not in loaded_files:
            try:
                loaded_files[path] = list(SeqIO.parse(path, "genbank"))
            except Exception as e:
                print(f"Warning: Could not parse {path}: {e}", file=sys.stderr)
                not_found.append(acc)
                continue

        records = loaded_files[path]
        record_cache[acc] = records
        record_source[acc] = "local"

        # Also cache all record accessions from the same file
        for rec in records:
            record_cache[rec.id] = records  # type: ignore
            record_source[rec.id] = "local"  # type: ignore
            if "." in rec.id:  # type: ignore
                record_cache[rec.id.rsplit(".", 1)[0]] = records  # type: ignore
                record_source[rec.id.rsplit(".", 1)[0]] = "local"  # type: ignore

    return not_found


def _fetch_ncbi_records(
    accessions: list[str],
    db: str,
    record_cache: dict[str, list[SeqRecord]],
    record_source: dict[str, str],
    threads: int | None = None,
) -> None:
    """Fetch GenBank records from NCBI in parallel and populate the cache.

    Args:
        accessions: Accessions to fetch.
        db: Entrez database name.
        record_cache: Cache to populate.
        record_source: Tracks where each accession was loaded from.
        threads: Number of threads.
    """
    if not accessions:
        return

    # Cap concurrency to avoid NCBI rate limits
    max_threads = min(threads or 3, 3) if not Entrez.api_key else min(threads or 10, 10)

    print(f"Fetching {len(accessions)} record(s) from NCBI ({db})...", file=sys.stderr)

    def _worker(acc: str) -> tuple[str, list[SeqRecord]]:
        return acc, fetch_genbank_record(acc, db=db)

    with ThreadPoolExecutor(max_workers=max_threads) as executor:
        futures = {executor.submit(_worker, acc): acc for acc in accessions}
        completed = 0

        for future in as_completed(futures):
            acc = futures[future]
            try:
                accession, records = future.result()
                record_cache[accession] = records
                if records:
                    record_source[accession] = "ncbi"
                else:
                    record_source[accession] = "not_found"
                completed += 1
                if completed % 25 == 0 or completed == len(accessions):
                    print(f"  Fetched {completed}/{len(accessions)}", file=sys.stderr)
            except Exception as e:
                print(f"  Error fetching {acc}: {e}", file=sys.stderr)
                record_cache[acc] = []
                record_source[acc] = "not_found"


def prefetch_records(
    accessions: list[str],
    db: str,
    gbff_dir: Path | None,
    gbff_index: dict[str, Path],
    record_cache: dict[str, list[SeqRecord]],
    record_source: dict[str, str],
    threads: int | None = None,
) -> None:
    """Pre-fetch all required GenBank records (local first, then NCBI fallback).

    When *gbff_dir* is provided, local files are attempted first.  Any
    accession that cannot be resolved locally is automatically fetched from
    NCBI.  If neither source succeeds the accession is marked ``"not_found"``
    in *record_source*.

    Args:
        accessions: Unique accessions to resolve.
        db: Entrez database (``'nucleotide'`` or ``'protein'``).
        gbff_dir: Optional local GBFF directory.
        gbff_index: Local GBFF file index.
        record_cache: Cache to populate.
        record_source: Dict tracking provenance (``"local"``, ``"ncbi"``, ``"not_found"``).
        threads: Number of threads for NCBI fetches.
    """
    # Deduplicate, skipping already cached accessions
    remaining = [a for a in dict.fromkeys(accessions) if a and a not in record_cache]

    if not remaining:
        return

    # 1. Try local files first
    if gbff_dir and gbff_index:
        remaining = _load_local_records(remaining, gbff_index, record_cache, record_source)
        local_count = len(accessions) - len(remaining)
        if local_count:
            print(f"Loaded {local_count} accession(s) from local files", file=sys.stderr)

    # 2. Fetch anything still missing from NCBI
    if remaining:
        _fetch_ncbi_records(remaining, db, record_cache, record_source, threads)

    # 3. Mark any that still have no records
    for acc in accessions:
        if acc and acc not in record_source:
            if acc in record_cache and record_cache[acc]:
                record_source[acc] = "ncbi"
            else:
                record_source[acc] = "not_found"


# ---------------------------------------------------------------------------
# Gene Annotation Extraction
# ---------------------------------------------------------------------------


def get_qualifier(feature: SeqFeature, key: str, default: str = "") -> str:
    """Safely extract the first value of a qualifier from a SeqFeature.

    Args:
        feature: BioPython SeqFeature.
        key: Qualifier key.
        default: Default value if key is absent.

    Returns:
        First qualifier value, or *default*.
    """
    values = feature.qualifiers.get(key, [default])
    return values[0] if values else default


def find_overlapping_features(
    record: SeqRecord,
    start: int,
    end: int,
) -> list[SeqFeature]:
    """Find annotation features that overlap a coordinate range.

    Searches CDS, rRNA, tRNA, ncRNA, tmRNA, misc_RNA, and gene features.
    Results are sorted by overlap size (descending), with CDS preferred over
    gene features at equal overlap.

    Args:
        record: GenBank SeqRecord.
        start: Hit start coordinate (1-based, as reported by BLAST).
        end: Hit end coordinate (1-based, as reported by BLAST).

    Returns:
        List of overlapping SeqFeatures, best match first.
    """
    # Normalise to 0-based half-open interval for comparison
    q_start = min(start, end) - 1
    q_end = max(start, end)

    overlapping: list[tuple[int, int, SeqFeature]] = []

    for feature in record.features:
        if feature.type not in ANNOTATION_FEATURE_TYPES:
            continue

        feat_start = int(feature.location.start)  # type: ignore
        feat_end = int(feature.location.end)  # type: ignore

        overlap_start = max(q_start, feat_start)
        overlap_end = min(q_end, feat_end)
        overlap = overlap_end - overlap_start

        if overlap > 0:
            priority = FEATURE_TYPE_PRIORITY.get(feature.type, 99)
            overlapping.append((overlap, priority, feature))

    # Sort: largest overlap first, then by feature type priority (CDS first)
    overlapping.sort(key=lambda x: (-x[0], x[1]))
    return [feat for _, _, feat in overlapping]


def find_feature_by_protein_id(
    records: list[SeqRecord],
    protein_id: str,
) -> tuple[SeqFeature | None, SeqRecord | None]:
    """Find a CDS feature by protein_id across multiple GenBank records.

    Args:
        records: List of GenBank SeqRecords to search.
        protein_id: Protein accession to match.

    Returns:
        ``(feature, record)`` if found, otherwise ``(None, None)``.
    """
    pid_no_ver = protein_id.rsplit(".", 1)[0] if "." in protein_id else protein_id

    for record in records:
        for feature in record.features:
            if feature.type != "CDS":
                continue
            for pid in feature.qualifiers.get("protein_id", []):
                if pid == protein_id or pid.rsplit(".", 1)[0] == pid_no_ver:
                    return feature, record
    return None, None


def extract_gene_annotation(feature: SeqFeature, record: SeqRecord, prefix: str = "gene") -> dict[str, Any]:
    """Extract annotation information from a genomic feature.

    Args:
        feature: SeqFeature (CDS, rRNA, gene, etc.).
        record: Parent SeqRecord.
        prefix: Column name prefix (``"gene"`` for hit, ``"query_gene"`` for query).

    Returns:
        Dictionary with all annotation fields, keyed with *prefix*.
    """
    location = feature.location
    strand = location.strand  # type: ignore
    feat_start = int(location.start) + 1  # type: ignore ; 1-based inclusive
    feat_end = int(location.end)  # type: ignore ; 1-based inclusive

    gene_name = get_qualifier(feature, "gene")
    locus_tag = get_qualifier(feature, "locus_tag")
    old_locus_tag = get_qualifier(feature, "old_locus_tag")
    product = get_qualifier(feature, "product")
    protein_id = get_qualifier(feature, "protein_id")
    note = get_qualifier(feature, "note")

    # Pseudo detection
    is_pseudo = "pseudo" in feature.qualifiers or "pseudogene" in feature.qualifiers
    pseudo_reason = get_qualifier(feature, "pseudogene") if "pseudogene" in feature.qualifiers else ""

    # Nucleotide sequence
    try:
        nuc_seq = str(feature.location.extract(record.seq))  # type: ignore
    except Exception:
        nuc_seq = ""

    # Protein sequence
    protein_seq = get_qualifier(feature, "translation")
    if not protein_seq and feature.type == "CDS" and not is_pseudo:
        try:
            protein_seq = str(feature.location.extract(record.seq).translate(table=11, to_stop=True))  # type: ignore[union-attr]
        except Exception:
            protein_seq = ""

    # gene features normally lack protein info
    if feature.type == "gene":
        protein_id = ""
        protein_seq = ""

    return {
        f"{prefix}_name": gene_name,
        f"{prefix}_locus_tag": locus_tag,
        f"{prefix}_old_locus_tag": old_locus_tag,
        f"{prefix}_product": product,
        f"{prefix}_protein_id": protein_id,
        f"{prefix}_start": feat_start,
        f"{prefix}_end": feat_end,
        f"{prefix}_strand": "+" if strand == 1 else "-" if strand == -1 else "?",
        f"{prefix}_feature_type": feature.type,
        f"{prefix}_is_pseudo": is_pseudo,
        f"{prefix}_pseudo_reason": pseudo_reason,
        f"{prefix}_nucleotide_seq": nuc_seq,
        f"{prefix}_protein_seq": protein_seq,
        f"{prefix}_note": note,
    }


def _make_empty_annotation(prefix: str = "gene") -> dict[str, Any]:
    """Return an empty annotation dict with all expected keys."""
    return {
        f"{prefix}_name": "",
        f"{prefix}_locus_tag": "",
        f"{prefix}_old_locus_tag": "",
        f"{prefix}_product": "",
        f"{prefix}_protein_id": "",
        f"{prefix}_start": "",
        f"{prefix}_end": "",
        f"{prefix}_strand": "",
        f"{prefix}_feature_type": "",
        f"{prefix}_is_pseudo": "",
        f"{prefix}_pseudo_reason": "",
        f"{prefix}_nucleotide_seq": "",
        f"{prefix}_protein_seq": "",
        f"{prefix}_note": "",
    }


def _find_record_for_accession(
    accession: str,
    records: list[SeqRecord],
) -> SeqRecord | None:
    """Select the SeqRecord whose accession matches *accession*."""
    acc_no_ver = accession.rsplit(".", 1)[0] if "." in accession else accession

    for rec in records:
        if rec.id == accession or rec.name == accession:
            return rec
        rec_no_ver = rec.id.rsplit(".", 1)[0] if "." in rec.id else rec.id  # type: ignore
        if rec_no_ver == acc_no_ver:
            return rec

    # Fallback: single-record files
    return records[0] if len(records) == 1 else None


# ---------------------------------------------------------------------------
# Hit (subject) annotation
# ---------------------------------------------------------------------------


def _annotate_protein_subject(
    accession: str,
    records: list[SeqRecord],
    prefix: str = "gene",
) -> dict[str, Any]:
    """Annotate a protein subject hit from GenBank records.

    Tries to locate a CDS feature by ``protein_id``.  Falls back to basic
    annotation from protein record metadata.

    Args:
        accession: Hit accession (protein ID).
        records: GenBank records to search.
        prefix: Column-name prefix.

    Returns:
        Annotation dictionary.
    """
    feature, rec = find_feature_by_protein_id(records, accession)
    if feature and rec:
        return extract_gene_annotation(feature, rec, prefix=prefix)

    # Fallback: extract what we can from the protein record
    rec = records[0]
    ann = _make_empty_annotation(prefix)
    ann[f"{prefix}_protein_seq"] = str(rec.seq) if rec.seq else ""
    ann[f"{prefix}_product"] = rec.description or ""

    for feat in rec.features:
        if feat.type == "CDS":
            ann[f"{prefix}_name"] = get_qualifier(feat, "gene")
            ann[f"{prefix}_locus_tag"] = get_qualifier(feat, "locus_tag")
            coded_by = get_qualifier(feat, "coded_by")
            if coded_by:
                ann[f"{prefix}_note"] = f"coded_by={coded_by}"
                match = re.match(
                    r"(?:complement\()?([^:]+):<?(\d+)\.\.>?(\d+)\)?",
                    coded_by,
                )
                if match:
                    ann[f"{prefix}_start"] = int(match.group(2))
                    ann[f"{prefix}_end"] = int(match.group(3))
                    ann[f"{prefix}_strand"] = "-" if coded_by.startswith("complement") else "+"
            break
        elif feat.type == "Protein":
            ann[f"{prefix}_product"] = get_qualifier(feat, "name") or rec.description

    return ann


def annotate_hits(
    df: pd.DataFrame,
    blast_program: str,
    record_cache: dict[str, list[SeqRecord]],
    record_source: dict[str, str],
) -> pd.DataFrame:
    """Add gene-level annotation columns to every row of *df*.

    For nucleotide-subject programs (blastn, tblastn, tblastx) the gene
    overlapping the hit coordinates is located.  For protein-subject programs
    (blastp, blastx) the matching CDS is found by ``protein_id``.

    A ``hit_record_source`` column is added indicating where the GenBank
    record came from: ``"local"``, ``"ncbi"``, or ``"not_found"``.

    Args:
        df: DataFrame with BLAST results.
        blast_program: BLAST program name.
        record_cache: Pre-populated record cache.
        record_source: Provenance tracking dict.

    Returns:
        New DataFrame with annotation columns appended.
    """
    subject_type = BLAST_PROGRAMS.get(blast_program, {}).get("subject_type", "nucleotide")
    annotations: list[dict[str, Any]] = []

    for _, row in df.iterrows():
        accession = row["hit_accession"]
        records = record_cache.get(accession, [])
        source = record_source.get(accession, "not_found")

        base: dict[str, Any] = {"hit_record_source": source}

        if not records:
            note = "# record not found in local files or NCBI" if source == "not_found" else ""
            ann = _make_empty_annotation("gene")
            ann["gene_note"] = note
            annotations.append({**base, **ann})
            continue

        if subject_type == "nucleotide":
            target_record = _find_record_for_accession(accession, records)
            if target_record is None:
                annotations.append({**base, **_make_empty_annotation("gene")})
                continue

            features = find_overlapping_features(target_record, row["hit_start"], row["hit_end"])
            if features:
                annotations.append({**base, **extract_gene_annotation(features[0], target_record, prefix="gene")})
            else:
                ann = _make_empty_annotation("gene")
                ann["gene_product"] = "intergenic"
                annotations.append({**base, **ann})
        else:
            annotations.append({**base, **_annotate_protein_subject(accession, records, prefix="gene")})

    ann_df = pd.DataFrame(annotations)
    return pd.concat([df.reset_index(drop=True), ann_df.reset_index(drop=True)], axis=1)


# ---------------------------------------------------------------------------
# Query feature annotation
# ---------------------------------------------------------------------------


def _resolve_query_accessions(df: pd.DataFrame) -> dict[str, str]:
    """Build a mapping from query_id → resolved NCBI-like accession.

    Examines ``query_id`` and ``query_description`` columns and tries to
    extract a valid accession from each.

    Args:
        df: BLAST results DataFrame.

    Returns:
        Dict mapping ``query_id`` values to resolved accession strings.
        Entries where no accession could be extracted are omitted.
    """
    mapping: dict[str, str] = {}

    for _, grp in df.groupby("query_id", sort=False):
        row = grp.iloc[0]
        qid: str = row["query_id"]

        if qid in mapping:
            continue

        # Try to extract from query_id first, then from description
        acc = extract_accession(qid)
        if not acc:
            acc = extract_accession(str(row.get("query_description", "")))
        if acc:
            mapping[qid] = acc

    return mapping


def annotate_queries(
    df: pd.DataFrame,
    blast_program: str,
    record_cache: dict[str, list[SeqRecord]],
    record_source: dict[str, str],
) -> pd.DataFrame:
    """Add query-gene annotation columns to every row of *df*.

    Works analogously to :func:`annotate_hits`, but for the *query* side.
    For nucleotide-query programs (blastn, blastx, tblastx) the gene
    overlapping the query coordinates is located.  For protein-query programs
    (blastp, tblastn) the query protein is looked up by accession.

    Args:
        df: DataFrame (already annotated with hit info).
        blast_program: BLAST program name.
        record_cache: Pre-populated record cache (must include query records).
        record_source: Provenance tracking dict.

    Returns:
        New DataFrame with ``query_gene_*`` columns appended.
    """
    query_type = BLAST_PROGRAMS.get(blast_program, {}).get("query_type", "nucleotide")
    query_acc_map = _resolve_query_accessions(df)

    annotations: list[dict[str, Any]] = []

    for _, row in df.iterrows():
        qid: str = row["query_id"]
        accession = query_acc_map.get(qid)

        if not accession:
            ann = _make_empty_annotation("query_gene")
            ann["query_gene_note"] = "# query accession could not be resolved"
            annotations.append({"query_record_source": "not_found", **ann})
            continue

        records = record_cache.get(accession, [])
        source = record_source.get(accession, "not_found")

        base: dict[str, Any] = {"query_record_source": source}

        if not records:
            note = "# record not found in local files or NCBI" if source == "not_found" else ""
            ann = _make_empty_annotation("query_gene")
            ann["query_gene_note"] = note
            annotations.append({**base, **ann})
            continue

        if query_type == "nucleotide":
            target_record = _find_record_for_accession(accession, records)
            if target_record is None:
                annotations.append({**base, **_make_empty_annotation("query_gene")})
                continue

            features = find_overlapping_features(target_record, row["query_start"], row["query_end"])
            if features:
                annotations.append({**base, **extract_gene_annotation(features[0], target_record, prefix="query_gene")})
            else:
                ann = _make_empty_annotation("query_gene")
                ann["query_gene_product"] = "intergenic"
                annotations.append({**base, **ann})
        else:
            # Protein query (blastp, tblastn)
            annotations.append({**base, **_annotate_protein_subject(accession, records, prefix="query_gene")})

    ann_df = pd.DataFrame(annotations)
    return pd.concat([df.reset_index(drop=True), ann_df.reset_index(drop=True)], axis=1)


# ---------------------------------------------------------------------------
# Column Ordering
# ---------------------------------------------------------------------------

_COLUMN_ORDER = [
    "blast_program",
    # Query identifiers
    "query_id",
    "query_description",
    "query_length",
    # Query-gene annotation
    "query_record_source",
    "query_gene_name",
    "query_gene_locus_tag",
    "query_gene_old_locus_tag",
    "query_gene_product",
    "query_gene_protein_id",
    "query_gene_start",
    "query_gene_end",
    "query_gene_strand",
    "query_gene_feature_type",
    "query_gene_is_pseudo",
    "query_gene_pseudo_reason",
    # Hit identifiers
    "hit_accession",
    "hit_description",
    "hit_length",
    # HSP statistics
    "hsp_num",
    "evalue",
    "bit_score",
    "score",
    "identity",
    "identity_pct",
    "positives",
    "positives_pct",
    "gaps",
    "alignment_length",
    "query_start",
    "query_end",
    "query_frame",
    "hit_start",
    "hit_end",
    "hit_frame",
    # Hit-gene annotation
    "hit_record_source",
    "gene_name",
    "gene_locus_tag",
    "gene_old_locus_tag",
    "gene_product",
    "gene_protein_id",
    "gene_start",
    "gene_end",
    "gene_strand",
    "gene_feature_type",
    "gene_is_pseudo",
    "gene_pseudo_reason",
    # Sequences (at end for readability)
    "gene_nucleotide_seq",
    "gene_protein_seq",
    "gene_note",
    "query_gene_nucleotide_seq",
    "query_gene_protein_seq",
    "query_gene_note",
    "query_seq",
    "hit_seq",
    "midline",
]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main(
    blast_xml: Path,
    output: Path | None = None,
    gbff_dir: Path | None = None,
    email: str = "user@example.com",
    api_key: str | None = None,
    threads: int | None = None,
    strand_filter: str | None = None,
) -> None:
    """Parse BLAST results and annotate hits with gene information.

    Args:
        blast_xml: Path to BLAST XML / XML2 file.
        output: Output TSV path, or ``None`` for stdout.
        gbff_dir: Optional directory with local GBFF files.
        email: Email for NCBI Entrez.
        api_key: NCBI API key.
        threads: Number of threads.
        strand_filter: Strand / frame filter string.
    """
    # Configure Entrez
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    # 1. Parse BLAST results
    df, blast_program = parse_blast_results(blast_xml)
    if df.empty:
        print("No results to process.", file=sys.stderr)
        return

    # 2. Apply strand / frame filter
    if strand_filter:
        df = filter_by_strand(df, strand_filter, blast_program)
        if df.empty:
            print("No results remaining after strand filtering.", file=sys.stderr)
            return

    # 3. Build local GBFF index
    gbff_index: dict[str, Path] = {}
    if gbff_dir:
        if not gbff_dir.is_dir():
            print(f"Error: GBFF directory not found: {gbff_dir}", file=sys.stderr)
            sys.exit(1)
        gbff_index = build_local_gbff_index(gbff_dir)

    # 4. Pre-fetch GenBank records for *hits*
    record_cache: dict[str, list[SeqRecord]] = {}
    record_source: dict[str, str] = {}

    unique_hit_accessions = df["hit_accession"].unique().tolist()
    hit_subject_type = BLAST_PROGRAMS.get(blast_program, {}).get("subject_type", "nucleotide")
    hit_db = "nucleotide" if hit_subject_type == "nucleotide" else "protein"
    prefetch_records(unique_hit_accessions, hit_db, gbff_dir, gbff_index, record_cache, record_source, threads)

    # 5. Pre-fetch GenBank records for *queries*
    query_acc_map = _resolve_query_accessions(df)
    unique_query_accessions = list(dict.fromkeys(query_acc_map.values()))

    if unique_query_accessions:
        query_type = BLAST_PROGRAMS.get(blast_program, {}).get("query_type", "nucleotide")
        query_db = "nucleotide" if query_type == "nucleotide" else "protein"
        print(
            f"Resolving {len(unique_query_accessions)} query accession(s)...",
            file=sys.stderr,
        )
        prefetch_records(unique_query_accessions, query_db, gbff_dir, gbff_index, record_cache, record_source, threads)
    else:
        print("Warning: Could not resolve any query accessions from BLAST identifiers.", file=sys.stderr)

    # 6. Annotate hits
    print("Annotating hits with gene information...", file=sys.stderr)
    result = annotate_hits(df, blast_program, record_cache, record_source)

    # 7. Annotate queries
    print("Annotating queries with gene information...", file=sys.stderr)
    result = annotate_queries(result, blast_program, record_cache, record_source)

    # Re-order columns
    ordered_cols = [c for c in _COLUMN_ORDER if c in result.columns]
    extra_cols = [c for c in result.columns if c not in _COLUMN_ORDER]
    result = result[ordered_cols + extra_cols]

    # 8. Output
    if output:
        result.to_csv(output, sep="\t", index=False)
        print(f"\nResults written to: {output}", file=sys.stderr)
    else:
        result.to_csv(sys.stdout, sep="\t", index=False)

    # Summary
    total = len(result)
    hit_annotated = result["gene_locus_tag"].astype(str).ne("").sum()
    query_annotated = result["query_gene_locus_tag"].astype(str).ne("").sum()

    not_found_hits = (result["hit_record_source"] == "not_found").sum()
    not_found_queries = (result["query_record_source"] == "not_found").sum()

    print(f"Total annotated HSPs: {total}", file=sys.stderr)
    print(f"HSPs with hit gene annotation: {hit_annotated}", file=sys.stderr)
    print(f"HSPs with query gene annotation: {query_annotated}", file=sys.stderr)
    if not_found_hits:
        print(f"HSPs where hit record was not found: {not_found_hits}", file=sys.stderr)
    if not_found_queries:
        print(f"HSPs where query record was not found: {not_found_queries}", file=sys.stderr)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Parse BLAST XML/XML2 results and annotate hits with gene information "
            "from GenBank records. Supports blastn, blastp, blastx, tblastn, and "
            "tblastx. The BLAST program type is auto-detected from the XML file. "
            "Both hit (subject) and query sequences are annotated with the gene "
            "they fall within. When --gbff-dir is given but a record is not found "
            "locally, NCBI is used as fallback; records that cannot be resolved "
            "are flagged as 'not_found' in the output."
        ),
    )

    parser.add_argument(
        "blast_xml",
        type=Path,
        help="Path to BLAST XML (outfmt 5) or XML2 (outfmt 16) output file",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output TSV file path (default: stdout)",
    )
    parser.add_argument(
        "--gbff-dir",
        type=Path,
        default=None,
        help=(
            "Directory containing GBFF/GBK files for local lookup. "
            "Files are matched to BLAST hits by record accession, filename "
            "stem, BioSample, or assembly accession. Accessions not found "
            "locally are automatically fetched from NCBI as fallback."
        ),
    )
    parser.add_argument(
        "--email",
        type=str,
        default="user@example.com",
        help="Email address for NCBI Entrez (required when fetching from NCBI)",
    )
    parser.add_argument(
        "--api-key",
        type=str,
        default=None,
        help="NCBI API key for higher rate limits",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Number of threads for parallel processing. Defaults to number of CPUs.",
    )
    parser.add_argument(
        "--strand",
        type=str,
        default=None,
        choices=["-3", "-2", "-1", "1", "2", "3", "-", "+"],
        help=(
            "Filter results by strand/frame. Use -3,-2,-1,1,2,3 for a specific "
            "reading frame, or -/+ for all minus/plus strand frames. "
            "For tblastn/tblastx this filters the subject (hit) frame; "
            "for blastx it filters the query frame; for blastn it filters "
            "the hit strand."
        ),
    )

    args = parser.parse_args()

    main(
        blast_xml=args.blast_xml,
        output=args.output,
        gbff_dir=args.gbff_dir,
        email=args.email,
        api_key=args.api_key,
        threads=args.threads,
        strand_filter=args.strand,
    )
