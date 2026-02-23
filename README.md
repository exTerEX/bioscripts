# Bioscripts

A collection of scripts developed during my PhD for automating small, repetitive bioinformatics tasks. These tools primarily interface with NCBI/EMBL, process genomic data, and parse output from common bioinformatics tools.

## Requirements

- Python 3.10+
- BioPython
- Pandas
- [NCBI datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/) (for `genomics/` shell scripts)

## Scripts

### antiSMASH (`antismash/`)

Tools for processing antiSMASH biosynthetic gene cluster (BGC) analysis results.

| Script | Description |
|--------|-------------|
| `bgc_completeness.py` | Cross-reference antiSMASH BGCs against a reference gene list (e.g. MIBiG) and report which core biosynthetic genes are present, partial, or missing per cluster |
| `bgc_synteny.py` | Compare gene order and orientation of BGCs across antiSMASH results with pairwise synteny scoring, collinearity, and rearrangement detection |
| `count_regions.py` | Count BGC regions across multiple antiSMASH result directories |
| `tabulate_regions.py` | Tabulate BGC regions with detailed metadata including KnownClusterBlast hits |

### BLAST (`blast/`)

Tools for parsing, annotating, and filtering BLAST results.

| Script | Description |
|--------|-------------|
| `annotate_blast.py` | Parse BLAST XML/XML2 results and annotate both query and hit sequences with gene information (coordinates, product, locus tag, pseudo status, nucleotide & protein sequences) from local GBFF files or NCBI. Auto-detects BLAST program type and supports strand/frame filtering. |

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
| `expression_annotate.py` | Annotate DESeq2/edgeR differential expression results with gene product descriptions and AMR/QS/BGC/MGE flags from GBFF or GFF annotation files |
| `extract_upstream_gene_region.py` | Extract upstream regions of genes from GenBank records |
| `get_organism_name_from_reference.py` | Extract organism names from NCBI nucleotide accessions |
| `operon_coexpression.py` | Predict operons from GBFF/GFF annotation (same-strand, close-proximity genes) and compute per-operon co-expression statistics from count/TPM matrices |
| `reference2assembly.py` | Convert NCBI nucleotide accessions to assembly accessions |

### Genomics (`genomics/`)

Tools for downloading, rehydrating, and analysing NCBI genome assemblies.

| Script | Description |
|--------|-------------|
| `amr_prevalence_matrix.py` | Build a presence/absence matrix of AMR genes per strain from AMRFinderPlus, RGI, ABRicate, or generic TSV output, with optional GBFF metadata enrichment |
| `deduplicate_genomes.py` | Remove duplicate genomes from an NCBI datasets zip file using the report from `detect_duplicate_strains.py` |
| `detect_duplicate_strains.py` | Detect potential duplicate strains from GBFF files or assembly accessions using metadata, checksums, and optional MinHash/k-mer sequence analysis |
| `fetch_genomes_accession.sh` | Download genomes by GenBank/RefSeq accession using NCBI datasets |
| `fetch_genomes_taxon.sh` | Download genomes by taxon ID using NCBI datasets |
| `hydrate_genomes.sh` | Rehydrate a dehydrated NCBI datasets zip archive |
| `island_finder.py` | Scan GBFF genomes for candidate genomic islands by GC deviation with structural evidence (flanking tRNAs, integrases, transposases, MGE features) |
| `resistance_gene_context.py` | Extract the genomic neighbourhood around resistance genes from GBFF files, reporting flanking genes, distances, and mobile genetic element markers |

## Usage

All scripts include a command-line interface. Use `-h` or `--help` for usage information:

```bash
python antismash/count_regions.py -h
python fuzznuc/metadata2table.py -h
```

Most scripts support multithreading via `--threads` and NCBI scripts accept `--email` and `--api-key` for Entrez queries.

## License

See [LICENSE](LICENSE) for details.
