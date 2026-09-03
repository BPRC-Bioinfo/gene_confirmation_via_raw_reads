# Gene confirmation via raw reads

This Snakemake workflow uses output excel files from VDJ-Insights or annotation refinement scripts to identify the raw reads coverage of genes.
It can be used to evaluate gene existence with raw reads coverage.


## Overview

The workflow:

- Finds gene coordinates from aligned BAM files (per-sample), using a region BED file and an annotation Excel file.
- Builds per-segment and per-region tables and outputs per-segment Excel summary files.

## Prerequisites

- Git (to clone the repository)
- Python 3 (>=3.8 recommended)
- Conda (Miniconda or Anaconda) — recommended
- Snakemake (installed via Conda/mamba or pip). This repository is tested with Snakemake that supports `--use-conda`.


## Installation (recommended)

1. Clone the repository:

```bash
git clone https://github.com/BPRC-Bioinfo/confirm_gene_via_raw_reads.git
cd confirm_gene_via_raw_reads
```

2. (Optional) Create a small Conda environment that contains Snakemake and Mamba (recommended):

```bash
# using conda
conda create -n snakemake -c conda-forge -c bioconda snakemake python=3.10 -y
conda activate snakemake

# or using mamba (faster)
conda create -n snakemake -c conda-forge mamba -y
mamba create -n snakemake -c conda-forge -c bioconda snakemake python=3.10 -y
conda activate snakemake
```

3. Verify snakemake is available:

```bash
snakemake --version
```

## Configure inputs

Edit `config.yaml` to point to your data. Example `config.yaml` values:

```yaml
bam_dir: /path/to/aligned_bam/   # directory containing sample BAMs (indexed .bam and .bai)
region_file: region_coordinates.bed
input_excel: your_annotation_file.xlsx
fasta_dir: /path/to/haps/        # directory with FASTA files referenced by the pipeline
```

Important notes:

- BAM files should be coordinate-sorted and have corresponding index files (`.bai`) in the same directory or alongside each BAM.
- `region_file` should be a BED file listing the genomic regions to analyze.

An example of a region_file is as follows:

```
contig_name  region_start_coord  region_end_coord  region_name
```

```
hap1_GCA_049350105.2_CM111665_2_f1	94425001	95632609	GCA_049350105.2_haps_hap1_TRA
hap1_T2TMmul8  89238939	90297198	T2T-MFA8v1.1_haps_hap1_TRB

```

**Region naming convention:** The `region_name` column should follow the pattern `name_(REGION_NAME)` where `REGION_NAME` is one of:
- `TRA` (T cell receptor alpha)
- `TRB` (T cell receptor beta)
- `TRG` (T cell receptor gamma)

- `input_excel` should be the excel output from either VDJ-Insight or annotation refinement script.


## Run the pipeline (local / single machine)

From the repository root, you can run the full workflow with Conda environment creation handled by Snakemake:

Basic run (uses Conda environments defined by rules):

```bash
snakemake -s coordinate_pipe4.smk --use-conda --cores 8
```

## Outputs

The workflow produces per-region output Excel files with the following columns:

| Column | Description |
|--------|-------------|
| Sample | Sample identifier |
| Short name | Short sample name |
| No reads PB | Number of PacBio reads |
| No reads ONT | Number of Oxford Nanopore reads |
| No reads PB 100% | Number of PacBio reads with 100% coverage |
| No reads ONT 100% | Number of Oxford Nanopore reads with 100% coverage |
| Region | Genomic region |
| Segment | Gene segment (TRA, TRB, or TRG) |
| Start coord | Start coordinate |
| End coord | End coordinate |
| Strand | DNA strand (+/-) |
| Target length | Length of target sequence |
| Target seq | Target sequence |
| SNP Positions PB | SNP positions in PacBio reads |
| SNP Positions ONT | SNP positions in Oxford Nanopore reads |

Output files are generated per segment (a, b, and g at present).

Check the Snakefile (`coordinate_pipe4.smk`) for the exact names and output directories used by each rule if you need to locate specific files.

## Example full workflow (local, 8 cores, mamba for conda)

```bash
conda activate snakemake            # activate your snakemake environment
cd confirm_gene_via_raw_reads
# edit config.yaml so paths point to your data
snakemake -s coordinate_pipe4.smk --use-conda --conda-frontend mamba --cores 8 --latency-wait 60
```
