# Bioscripts

A collection of Python scripts developed during my PhD for automating small, repetitive bioinformatics tasks. These tools primarily interface with NCBI/EMBL, process genomic data, and parse output from common bioinformatics tools.

## Requirements

- Python 3.10+
- BioPython
- Pandas

## Scripts

### antiSMASH (`antismash/`)

Tools for processing antiSMASH biosynthetic gene cluster (BGC) analysis results.

| Script | Description |
|--------|-------------|
| `count_regions.py` | Count BGC regions across multiple antiSMASH result directories |
| `tabulate_regions.py` | Tabulate BGC regions with detailed metadata including KnownClusterBlast hits |

### Fuzznuc (`fuzznuc/`)

Tools for working with EMBOSS fuzznuc motif search results.

| Script | Description |
|--------|-------------|
| `metadata2table.py` | Find genes near fuzznuc motif hits and extract metadata from NCBI GenBank |
| `merge_gbk.py` | Merge fuzznuc results into a GenBank genome file as misc_binding features |

### Miscellaneous (`misc/`)

General-purpose utilities for NCBI data retrieval and processing.

| Script | Description |
|--------|-------------|
| `get_organism_name_from_reference.py` | Extract organism names from NCBI nucleotide accessions |
| `reference2assembly.py` | Convert NCBI nucleotide accessions to assembly accessions |

## Usage

All scripts include a command-line interface. Use `-h` or `--help` for usage information:

```bash
python antismash/count_regions.py -h
python fuzznuc/metadata2table.py -h
```

Most scripts support multithreading via `--threads` and NCBI scripts accept `--email` and `--api-key` for Entrez queries.

## License

See [LICENSE](LICENSE) for details.
