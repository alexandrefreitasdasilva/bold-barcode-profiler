# BOLD Barcode Profiler

[![Python Version](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![BLAST+](https://img.shields.io/badge/NCBI%20BLAST%2B-2.12%2B-brightgreen.svg)](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey.svg)]()

> **Automated end-to-end pipeline for host and vector DNA barcode identification (COI, COII, 16S, cytb) from metagenomic assemblies using BOLD Systems.**

---

## 📌 Overview

In metagenomic and virome studies, identifying the host organism, biological vector, or blood meal source is essential for contextualizing detected viruses and pathogens. **BOLD Barcode Profiler** provides an automated, reproducible workflow that accepts:
1. **Assembled metagenomic contigs** (`.fasta` from metaSPAdes, MEGAHIT, Trinity, Flye, etc.)
2. **Raw reference data exported from BOLD Systems** (`.tsv` for any taxonomic group, e.g., *Culicidae*, *Diptera*, *Chordata*, *Arthropoda*).

The pipeline automatically formats the database, performs local alignment using `blastn`, cross-references results with comprehensive taxonomic metadata, determines best hits per contig, and exports publication-ready tables and figures.

---

## 🔄 Workflow

```
  ┌───────────────────────────────┐     ┌───────────────────────────────┐
  │  Metagenomic Contigs (FASTA)  │     │   BOLD Systems Export (TSV)   │
  └──────────────┬────────────────┘     └──────────────┬────────────────┘
                 │                                     │
                 │                      ┌──────────────▼────────────────┐
                 │                      │   Database Sanitizer & Prep   │
                 │                      │  (Sanitize & build bold_db)   │
                 │                      └──────────────┬────────────────┘
                 │                                     │
                 └──────────────────────┬──────────────┘
                                        │
                         ┌──────────────▼──────────────┐
                         │   makeblastdb & blastn      │
                         │   (qcovs, pident, bitscore) │
                         └──────────────┬──────────────┘
                                        │
                         ┌──────────────▼──────────────┐
                         │  Metadata Cross-Referencing │
                         │  & Marker / Quality Filters │
                         └──────────────┬──────────────┘
                                        │
           ┌────────────────────────────┴───────────────────────────┐
           │                                                        │
┌──────────▼───────────┐                                 ┌──────────▼───────────┐
│   Tabular Summaries  │                                 │  Publication Figures │
│  - All filtered hits │                                 │  - 300 DPI PNG       │
│  - Best hits (TSV/CSV)                                 │  - Vector SVG        │
│  - 2-Column format   │                                 └──────────────────────┘
└──────────────────────┘
```

---

## ⚙️ Key Features

- **End-to-End Automation:** Runs the entire pipeline from raw inputs to final figures with a single command.
- **Robust BOLD Handling:** Automatically parses and sanitizes multi-gigabyte BOLD TSV dumps without hitting memory limits or `makeblastdb` local ID length restrictions (avoiding header length bugs).
- **Target Marker Filtering:** Easily restrict analysis to standard barcode markers like `COI-5P`, `COI`, `COII`, `16S`, `cytb`, or process all available markers.
- **Marker-Specific FASTA Extraction:** Automatically isolates identified contigs into separate FASTA files grouped by marker gene (e.g. `_contigs_COI-5P.fasta`, `_contigs_18S-3P.fasta`), annotated with taxonomic assignments for downstream phylogenetics and alignments.
- **NCBI Remote BLASTn Cross-Validation:** Optional one-command (`--remote-ncbi` or `--remote`) BLASTn verification against GenBank (`nt`) producing side-by-side BOLD vs. NCBI comparison tables.
- **Dual Output Formats:** Produces standard TSV/CSV tables plus a **2-column layout** optimized for copying directly into Word documents and Google Docs reports.
- **Publication Figures:** Generates multi-panel, publication-ready figures (PNG at 300 DPI and vector SVG) summarizing contig coverage, percent identity, bit score, and taxonomic assignments.

---

## 📦 Installation & Environment Setup

### Option 1: Using Conda / Mamba (Recommended)

Clone the repository and create the dedicated environment using the provided `environment.yml`:

```bash
git clone https://github.com/alexandrefreitasdasilva/bold-barcode-profiler.git
cd bold-barcode-profiler

# Create environment with Conda or Mamba
conda env create -f environment.yml

# Activate environment
conda activate bold-profiler
```

### Option 2: Using Micromamba

```bash
micromamba create -f environment.yml
micromamba activate bold-profiler
```

### Option 3: Manual Installation

Ensure **NCBI BLAST+** (`blastn` and `makeblastdb`) is installed in your system PATH, then install Python requirements:

```bash
pip install -r requirements.txt
```

---

## 🚀 Quick Start

### Basic Usage

Run the pipeline providing the contigs file and the raw BOLD TSV file:

```bash
python3 bold_barcode_profiler.py \
    --query my_contigs.fasta \
    --bold culicidae_bold_data.tsv \
    --outdir results/
```

### Advanced Usage with Filtering

Filter specifically for **COI** and **COII** barcodes with minimum 85% nucleotide identity and 60% query coverage:

```bash
python3 bold_barcode_profiler.py \
    --query /path/to/final.contigs.fasta \
    --bold /path/to/culicidae_bold.tsv \
    --outdir results_culicidae/ \
    --prefix sample126_culicidae \
    --marker "COI-5P,COII" \
    --min-identity 85.0 \
    --min-coverage 60.0 \
    --threads 8
```

### Cross-Validation with Remote NCBI BLASTn

Optionally cross-validate identified marker contigs against GenBank's full `nt` database to produce a side-by-side BOLD vs. NCBI comparison table:

```bash
python3 bold_barcode_profiler.py \
    --query /path/to/final.contigs.fasta \
    --bold /path/to/culicidae_bold.tsv \
    --outdir results_culicidae/ \
    --remote-ncbi
```

---

## 🧪 Testing with Example Data

A minimal example dataset is included in the `examples/` directory to verify your installation:

```bash
./examples/run_example.sh
```

Or run manually:

```bash
python3 bold_barcode_profiler.py \
    --query examples/example_contigs.fasta \
    --bold examples/example_bold.tsv \
    --outdir examples/example_results \
    --prefix demo_run \
    --marker "COI-5P,COII"
```

The output tables and figures will be generated in `examples/example_results/`.

---

## 📑 Command-Line Arguments

| Flag | Parameter | Type | Default | Description |
| :--- | :--- | :---: | :---: | :--- |
| `-q` | `--query` | `str` | *Required* | Path to input assembled contigs in FASTA format. |
| `-b` | `--bold` | `str` | *Required* | Path to raw BOLD TSV file exported from BOLD Systems. |
| `-o` | `--outdir` | `str` | `bold_results` | Directory to save all generated databases, alignments, and results. |
| `-p` | `--prefix` | `str` | *Basename* | Prefix for output filenames (defaults to query basename). |
| `-t` | `--threads` | `int` | `4` | Number of CPU threads to allocate for `blastn`. |
| `-m` | `--marker` | `str` | `all` | Specific barcode marker(s) to filter (e.g. `COI-5P`, `COII`, `16S`, or `all`). Multiple markers can be comma-separated. |
| | `--min-identity` | `float` | `80.0` | Minimum percentage nucleotide identity (`pident`). |
| | `--min-coverage` | `float` | `50.0` | Minimum query alignment coverage (`qcovs` in %). |
| | `--evalue` | `float` | `1e-5` | Maximum BLASTn E-value cutoff. |
| | `--max-target-seqs`| `int` | `5` | Maximum number of aligned target sequences to retain per query. |
| | `--skip-blast` | `flag` | `False` | Skip BLAST execution if raw BLAST results already exist in output dir. |
| | `--no-plot` | `flag` | `False` | Disable automatic figure generation. |
| | `--no-extract-fastas` | `flag` | `False` | Disable automatic extraction of marker-specific FASTA files. |
| | `--remote-ncbi`, `--remote` | `flag` | `False` | Execute remote NCBI BLASTn against GenBank (`nt`) for extracted marker contigs. |
| | `--ncbi-db` | `str` | `nt` | NCBI database for remote search (e.g., `nt`, `refseq_select_rna`). |
| | `--ncbi-evalue` | `float` | `1e-5` | Maximum E-value threshold for remote NCBI BLASTn. |
| | `--ncbi-max-target-seqs` | `int` | `5` | Maximum target sequences per query for remote NCBI BLASTn. |
| | `--ncbi-entrez-query` | `str` | `None` | Optional Entrez filter (e.g. `'Culicidae[Organism]'` or `'Diptera[Organism]'`). |

---

## 📂 Output Files

All outputs are organized in the designated `--outdir`:

| File | Description |
| :--- | :--- |
| `<prefix>_contigs_<marker>.fasta` | **Marker-partitioned contig FASTA** (e.g. `_contigs_COI-5P.fasta`, `_contigs_18S-3P.fasta`) with annotated headers. |
| `<prefix>_contigs_all_markers.fasta` | Consolidated FASTA with all identified marker contigs. |
| `<prefix>_bold_vs_ncbi_comparison.tsv` | Cross-validation table comparing BOLD identification vs. top NCBI GenBank hit (when `--remote-ncbi` is used). |
| `<prefix>_besthits.tsv` | Best hit per contig (highest bitscore and identity) with full BOLD taxonomic metadata. |
| `<prefix>_besthits.csv` | Comma-separated version of best hits for spreadsheet tools (Excel, LibreOffice). |
| `<prefix>_summary_2columns.tsv` | 2-column key-value format grouped per contig, ready to copy-paste into Google Docs / Word. |
| `<prefix>_all_filtered_hits.tsv` | All alignments passing user quality and marker thresholds. |
| `<prefix>_blastn_raw.tsv` | Raw BLASTn tabular output (`outfmt 6`). |
| `<prefix>_summary_figure.png` | **Multi-panel publication figure** at **300 DPI**. |
| `<prefix>_summary_figure.svg` | Vector format of the multi-panel summary figure for high-quality editing. |

### 🖼️ Example Summary Figure Panels:
* **Panel A:** Total contig length vs. aligned length against the BOLD reference.
* **Panel B:** Percentage nucleotide identity and BLASTn bit score.
* **Panel C:** Formatted metadata overview card with assigned species, genus, family, marker gene, alignment coverage, and BOLD Process IDs.

---

## 🌐 How to Download Data from BOLD Systems

1. Go to the [BOLD Systems Public Data Portal](http://www.boldsystems.org/).
2. In the Search bar, search for your target taxon (e.g., `Culicidae`, `Diptera`, `Chordata`, `Arthropoda`).
3. Under the **Download** tab, select **Combined (TSV)**.
4. Save the `.tsv` file and supply it directly to the `--bold` parameter.

---

## 🔬 Citation

If you use this workflow in your research, please cite:

1. **BOLD Systems:**  
   Ratnasingham, S., & Hebert, P. D. N. (2007). *BOLD: The Barcode of Life Data System (http://www.barcodinglife.org)*. Molecular Ecology Notes, 7(3), 355–364.
2. **NCBI BLAST+:**  
   Camacho, C., Coulouris, G., Avagyan, V., et al. (2009). *BLAST+: architecture and applications*. BMC Bioinformatics, 10, 421.

---

## 📄 License

This project is licensed under the [MIT License](LICENSE).
