# confirm_gene_via_raw_reads

This Snakemake workflow identifies gene coordinates from raw-read BAM files and produces segment-based summary Excel files.

## Overview

The workflow:

- Finds gene coordinates from aligned BAM files (per-sample), using a region BED file and an annotation Excel file.
- Builds per-segment and per-region tables and outputs per-segment Excel summary files.

Main files:

- `coordinate_pipe4.smk`: main Snakemake workflow
- `config.yaml`: configurable input paths for BAM files, the region BED file, the Excel annotation file, and the FASTA directory
- `region_coordinates.bed`: region coordinates used by the workflow
- `coordinate_study5.py`: script used by the coordinate_finder rule
- `bam_gene_matrix4.py`: script used for per-region analysis
- `gene_matrix.yaml`: Conda environment file used by parts of the workflow (used with `--use-conda`)

## Prerequisites

- Git (to clone the repository)
- Python 3 (>=3.8 recommended)
- Conda (Miniconda or Anaconda) — recommended
- Snakemake (installed via Conda/mamba or pip). This repository is tested with Snakemake that supports `--use-conda`.

Optional but recommended:

- Mamba (a faster Conda frontend) when creating/running Conda environments from Snakemake: `mamba`.
- A cluster scheduler (SLURM shown below) if you want to run region jobs on a cluster.

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

Note: You do NOT need to pre-create the environments listed in `gene_matrix.yaml`. Snakemake will create those automatically when using `--use-conda` (and `--conda-frontend mamba` if you have mamba installed).

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
- `input_excel` should match the expected annotation format used by `coordinate_study5.py` — check the script docstring or headers for the expected columns.

## Run the pipeline (local / single machine)

From the repository root, you can run the full workflow with Conda environment creation handled by Snakemake:

Basic run (uses Conda environments defined by rules):

```bash
snakemake -s coordinate_pipe4.smk --use-conda --cores 8
```

- `--cores 8` (or `-j 8`) tells Snakemake how many CPU cores to use in parallel.
- `--use-conda` instructs Snakemake to create and use the Conda environments defined in environment YAML files (e.g., `gene_matrix.yaml`).

If you have `mamba` installed and want Snakemake to use it for faster environment creation:

```bash
snakemake -s coordinate_pipe4.smk --use-conda --conda-frontend mamba --cores 8
```

Dry-run (no commands executed, shows what would run):

```bash
snakemake -s coordinate_pipe4.smk -n --use-conda
```

Run a single target (useful for debugging or re-running one step). Example: run the coordinate finder output file (replace with the actual expected output filename if different):

```bash
snakemake -s coordinate_pipe4.smk path/to/rhesus_genes_coordinates_from_bam_bed.txt --use-conda --cores 1
```

Run with verbose logging and a larger latency wait (helps on NFS filesystems):

```bash
snakemake -s coordinate_pipe4.smk --use-conda --cores 8 --verbose --latency-wait 60
```

## Run on an HPC cluster (SLURM example)

The per-region jobs are typically independent and suitable for cluster execution. Example using SLURM via `sbatch`:

```bash
snakemake -s coordinate_pipe4.smk --use-conda --jobs 200 --latency-wait 60 \
  --cluster "sbatch --time={resources.time} --cpus-per-task={threads} --mem={resources.mem} -A <your-account>"
```

- Adjust `--jobs` to the number of jobs you allow on the cluster.
- The `--cluster` string can be adapted to your cluster's SBATCH flags and the resources the rules declare.
- If your rules do not declare `resources.time` or `resources.mem`, you can hard-code cluster flags, e.g. `--cluster "sbatch -t 02:00:00 -c {threads} --mem=8G -A <acct>"`.

If you have a `cluster.yaml` configured for the project, use `--cluster-config cluster.yaml` or `--cluster "sbatch ... -p {params.partition}"` as needed.

## Outputs

The workflow produces per-region output files and combines them into per-segment Excel summary files. The main coordinate output filename expected by later steps is `rhesus_genes_coordinates_from_bam_bed.txt` (see Notes below) — do not rename it unless you adjust the Snakefile accordingly.

Check the Snakefile (`coordinate_pipe4.smk`) for the exact names and output directories used by each rule if you need to locate specific files.

## Troubleshooting

- Missing BAM/BAM index errors: ensure BAMs are coordinate-sorted and `.bai` indexes are present and named correctly.
- Conda environment creation fails: try using `--conda-frontend mamba` or pre-install required packages into a shared environment. Inspect the environment YAML (e.g., `gene_matrix.yaml`) for missing channels or packages.
- Permission / NFS latency issues: add `--latency-wait 60` (or higher) to your Snakemake command.
- If the pipeline cannot find the coordinate output (expected filename), either regenerate it with the coordinate-finder step or edit the Snakefile to match your filename consistently.

## Example full workflow (local, 8 cores, mamba for conda)

```bash
conda activate snakemake            # activate your snakemake environment
cd confirm_gene_via_raw_reads
# edit config.yaml so paths point to your data
snakemake -s coordinate_pipe4.smk --use-conda --conda-frontend mamba --cores 8 --latency-wait 60
```

## Notes

- The workflow expects the coordinate output file name to remain `rhesus_genes_coordinates_from_bam_bed.txt` unless you change it in the Snakefile.
- If you need help adapting the Snakefile for a cluster scheduler or adjusting resource requests for rules, inspect the rule definitions in `coordinate_pipe4.smk` and add/modify `threads`, `resources`, or `params` accordingly.

If you'd like, I can also:

- Add an example `config.example.yaml` filled with placeholder paths and sample values.
- Add a `Makefile` or wrapper script to simplify common runs.
- Provide a SLURM `cluster.yaml` suggestion based on the rule resources in `coordinate_pipe4.smk`.
