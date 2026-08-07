# confirm_gene_via_raw_reads

This Snakemake workflow identifies gene coordinates from raw-read BAM files and produces segment-based summary Excel files.

## Files

- coordinate_pipe4.smk: main Snakemake workflow
- config.yaml: configurable input paths for BAM files, the region BED file, the Excel annotation file, and the FASTA directory
- region_coordinates.bed: region coordinates used by the workflow
- coordinate_study5.py: script used by the coordinate_finder rule
- bam_gene_matrix4.py: script used for per-region analysis

## Prerequisites

- Python 3
- Snakemake
- Conda (for the environment defined in gene_matrix.yaml)
- The input files referenced in config.yaml must exist at the provided paths

## Configure inputs

Edit config.yaml and update the paths as needed:

```yaml
bam_dir: /path/to/aligned_bam/
region_file: region_coordinates.bed
input_excel: your_annotation_file.xlsx
fasta_dir: /path/to/haps/
```

## Run the pipeline

From this directory, run:

```bash
snakemake -s coordinate_pipe4.smk --use-conda
```

This will:

1. Run the coordinate finder step using the configured BAM directory, BED file, and Excel annotation file.
2. Build per-segment run tables.
3. Process each region job and create per-job output files.
4. Combine the results into per-segment Excel summaries.

## Notes

- The workflow expects the coordinate output file name to remain as rhesus_genes_coordinates_from_bam_bed.txt unless you change it in the Snakefile.
- If you want to dry-run the workflow first, use:

```bash
snakemake -s coordinate_pipe4.smk -n --use-conda
```
