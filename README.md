# Gene confirmation via raw reads

This Snakemake workflow use output excel file from VDJ-Insights or annotation refinement script to identify the raw reads coverage of the genes.
can be used to evaluate the gene existance with raw reads coverage.


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

An example of a region_file is as follow

```
contig_name  region_start_coord  region_end_coord  region_name
```

```
hap1_GCA_049350105.2_CM111665_2_f1	94425001	95632609	GCA_049350105.2_haps_hap1_TRA
hap1_T2TMmul8  89238939	90297198	T2T-MFA8v1.1_haps_hap1_TRA

```


- `input_excel` should be the excel output from either VDJ-Insight or annotation refinement script.


## Run the pipeline (local / single machine)

From the repository root, you can run the full workflow with Conda environment creation handled by Snakemake:

Basic run (uses Conda environments defined by rules):

```bash
snakemake -s coordinate_pipe.smk --use-conda --cores 8
```

## Outputs

The workflow produces per-region output files and combines them into per-segment Excel summary files. The main coordinate output filename expected by later steps is `rhesus_genes_coordinates_from_bam_bed.txt` (see Notes below) — do not rename it unless you adjust the Snakefile accordingly.

Check the Snakefile (`coordinate_pipe4.smk`) for the exact names and output directories used by each rule if you need to locate specific files.

## Troubleshooting

- Missing BAM/BAM index errors: ensure BAMs are coordinate-sorted and `.bai` indexes are present and named correctly.
- If the pipeline cannot find the coordinate output (expected filename), either regenerate it with the coordinate-finder step or edit the Snakefile to match your filename consistently.

## Example full workflow (local, 8 cores, mamba for conda)

```bash
conda activate snakemake            # activate your snakemake environment
cd confirm_gene_via_raw_reads
# edit config.yaml so paths point to your data
snakemake -s coordinate_pipe4.smk --use-conda --conda-frontend mamba --cores 8 --latency-wait 60
```
