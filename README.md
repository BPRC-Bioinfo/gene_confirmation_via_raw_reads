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
git clone https://github.com/BPRC-Bioinfo/gene_confirmation_via_raw_reads.git
cd gene_confirmation_via_raw_reads
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

Edit `config/config.yaml` to point to your data. Example `config.yaml` values:

```yaml
bam_dir: /path/to/aligned_bam/   # top-level directory containing sample BAM folders (indexed .bam and .bai)
region_file: region_coordinates.bed
input_excel: your_annotation_file.xlsx
fasta_dir: /path/to/haps/        # directory with FASTA files referenced by the pipeline
```

Important notes and conventions

- BAM files should be coordinate-sorted and have corresponding index files (`.bai`) in the same directory or alongside each BAM.
- The pipeline searches for BAMs under `bam_dir` inside per-sample folders. It looks for BAM filenames using the pattern:

  {sample_id}_{tech}_merged_sorted_primary.bam

  where `{tech}` is `pacbio` or `nanopore` (case-insensitive). Example:

  sample42_pacbio_merged_sorted_primary.bam

  If files do not follow this exact name, the pipeline will still list all BAMs found under the matching sample folders (see the BAM presence output) but the automatic selection of the "preferred" BAM uses the above prefix.

- `region_file` should be a BED-like TSV listing the genomic regions to analyze. Example:

```
contig_name  region_start_coord  region_end_coord  region_name
```

```
hap1_GCA_049350105.2_CM111665_2_f1	94425001	95632609	GCA_049350105.2_haps_hap1_TRA
hap1_T2TMmul8	89238939	90297198	T2T-MFA8v1.1_haps_hap1_TRB
```

- Region naming convention: the `region_name` column should end with one of `TRA`, `TRB`, `TRG` (the workflow currently supports these three locus names).

- `input_excel` should be the excel output from either VDJ-Insight or an annotation refinement script. The expected columns used are: Sample, Short name, Region, Segment, Start coord, End coord, Strand.

## Run the pipeline (local / single machine)

From the repository root, you can run the full workflow with Conda environment creation handled by Snakemake.

Basic run (uses Conda environments defined by rules):

```bash
# NOTE: the Snakefile in this repository is `coordinate_pipe.smk` (not coordinate_pipe4.smk)
snakemake -s coordinate_pipe.smk --use-conda --cores 8
```

Example with mamba for conda:

```bash
conda activate snakemake            # activate your snakemake environment
cd gene_confirmation_via_raw_reads
# edit config/config.yaml so paths point to your data
snakemake -s coordinate_pipe.smk --use-conda --conda-frontend mamba --cores 8 --latency-wait 60
```

## What the pipeline counts (definitions)

- "No reads PB" / "No reads ONT": total number of reads (PacBio / Oxford Nanopore) that overlap the investigated region. These include reads that only partially cover the region (i.e. they do not have to cover the full target length).

- "No reads PB 100%" / "No reads ONT 100%": number of reads that satisfy both of the following conditions (the pipeline uses both checks):
  - percent identity vs reference == 100
  - percent gene covered == 100

  In other words, the read must both match the reference sequence exactly (100% identity on compared positions) and cover 100% of the target region.

## Outputs

The workflow produces per-region output Excel files (one per segment: V, D, and J) with the following columns (these are the final summary columns):

| Column | Description |
|--------|-------------|
| Sample | Sample identifier |
| Short name | Short sample name |
| No reads PB | Number of PacBio reads (any coverage) |
| No reads ONT | Number of Oxford Nanopore reads (any coverage) |
| No reads PB 100% | Number of PacBio reads with 100% identity AND 100% gene coverage |
| No reads ONT 100% | Number of Oxford Nanopore reads with 100% identity AND 100% gene coverage |
| Region | Genomic region |
| Segment | Gene segment (V, D, or J) |
| Start coord | Start coordinate |
| End coord | End coordinate |
| Strand | DNA strand (+/-) |
| Target length | Length of target sequence |
| Target seq | Target sequence |
| SNP Positions PB | SNP positions in PacBio reads (see interpretation below) |
| SNP Positions ONT | SNP positions in Oxford Nanopore reads (see interpretation below) |

The pipeline generates other intermediate and per-job files. Brief overview of commonly generated files and directories (useful when you inspect results):

- results/{segment_lower}_region_run_table.tsv
  - A run table listing each analysis job (one row per sample/haplotype/tech). This is the CSV created by the `make_region_run_table` rule.

- {segment_lower}_region_jobs/
  - A directory containing per-job TSV files (one file per job) created by the `split_region_run_table` checkpoint. Each file contains the job parameters (bam path, fasta path, coordinates, output prefix, etc.).

- {segment_lower}_region_results/
  - Per-job output directory where each job produces multiple files with the given output prefix. For a job with prefix `results/<prefix>` you will typically find:
    - <prefix>.gene_read_matrix.tsv  — matrix with one row per read and one column per gene position (REFERENCE row first).
    - <prefix>.read_vs_reference.tsv — per-read summary: covered positions, matches, mismatches, deletions, percent_identity_vs_reference, percent_gene_covered.
    - <prefix>.position_base_percent.tsv — per-position base counts and percentages (used to find variable positions).
    - <prefix>.variable_sites.txt — condensed list of variable positions (created by summarize_variable_positions.py).
    - <prefix>.done — small text file summarizing key results (used to build the final Excel summary).
    - <prefix>.log — execution log for the job.

- {segment_lower}_region_results/.complete
  - A marker file created by `mark_region_complete` when all jobs for the segment finished.

- {segment_lower}_region_summary.xlsx
  - The final Excel summary for the segment, created by `combine_v_region_done.py` by aggregating all `.done` files. This is the file most users will open for a quick overview.

You do not need to dig into all intermediate files for typical analyses, but they are available if you want to inspect read-level matrices, per-position base percentages, or job logs for debugging.

## Segments / regions

The pipeline currently operates on TRA/TRB/TRG regions (these are the locus names) and the per-segment summaries are produced for gene segments V, D and J (see `SEGMENTS = ["V","D","J"]` in the Snakefile). The README previously mentioned "a, b, and g" — that was a leftover phrasing and has been corrected.

Including C (constant) regions: the code filters target regions based on TARGET_REGIONS = ["TRA","TRB","TRG"] in `scripts/coordinate_study.py`. To include additional regions (for example C regions), you would need to add the region suffix to TARGET_REGIONS and ensure the input region_name values end with the expected suffix. You may also need to ensure the input Excel has appropriate Segment values.

## How SNP positions are reported (example)

The per-position summary (`<prefix>.position_base_percent.tsv`) contains one row per gene position with counts and percentages for each base and the major base percent. The `summarize_variable_positions.py` script condenses variable positions into a compact text string that is placed into the "SNP Positions" fields of the final summary.

Format produced by summarize_variable_positions.py

- Example content of `<prefix>.variable_sites.txt` (single-line):

  5:A90.00;T10.00 | 12:G100.00

Interpretation:
- Each entry before a `|` is a variable position encoded as `<pos>:<Base><percent>;<Base2><percent2>;...`.
- `pos` is the 1-based position inside the target/reference sequence (i.e. position 1 is the first base of the target sequence used for this job).
- For `5:A90.00;T10.00` — at gene position 5, 90.00% of reads had base A, 10.00% had T.
- For `12:G100.00` — at gene position 12, 100.00% of reads had base G (i.e. no variation at that position in the reads used).

Notes:
- The script encodes deletions as `DEL` (shown as `N` in the intermediate percentages then presented as `DEL_percent` in the position file and summarized accordingly).
- SNP positions shown in `SNP Positions PB` or `SNP Positions ONT` in the final Excel are only filled if there are reads from the respective technology for that sample/region.

## Common troubleshooting

- If the Snakefile or script cannot find the expected BAM or FASTA contig names, inspect the per-job log file in `{segment_lower}_region_results` and the `.gene_read_matrix.tsv` files to see which contig names were used; `bam_gene_matrix.py` prints helpful diagnostics when contig names do not match between BAM and FASTA.

- If your BAMs use a different naming convention, either rename them to the `{sample_id}_{tech}_merged_sorted_primary.bam` pattern or place only the intended BAMs in the sample folder and use the `--bam-presence-output` in `scripts/coordinate_study.py` to inspect what the pipeline found.


## Example full workflow (local, 8 cores, mamba for conda)

```bash
conda activate snakemake            # activate your snakemake environment
cd gene_confirmation_via_raw_reads
# edit config/config.yaml so paths point to your data
snakemake -s coordinate_pipe.smk --use-conda --conda-frontend mamba --cores 8 --latency-wait 60
```
