# Gene confirmation via raw reads

This Snakemake workflow uses annotation Excel files produced by VDJ-Insights or annotation refinement scripts to investigate raw-read support for annotated genes.

The workflow evaluates whether annotated genes are supported by PacBio HiFi and/or Oxford Nanopore (ONT) reads by examining:

* The number of reads overlapping each annotated gene;
* The fraction of the target gene covered by each read;
* Sequence identity between each read and the target/reference sequence;
* Substitutions and deletions observed at individual positions;
* Agreement between PacBio and Oxford Nanopore observations.


The final output consists of summary Excel files for V, D, J, and constant-region gene segments, together with intermediate read-level and position-level files that can be used for manual inspection.


## Overview

The workflow:

* Identifies gene coordinates from aligned BAM files using a region BED file and an annotation Excel file.
* Determines PacBio and Oxford Nanopore read support for each annotated gene.
* Calculates read coverage and sequence identity relative to the target sequence.
* Reports positions at which the reads differ from the target sequence.
* Builds per-gene-segment summary tables for V, D, J genes and C constant segment.
* Produces Excel summary files as well as intermediate read-level and position-level result files.

## Prerequisites

* Git
* Python 3 (>=3.8 recommended)
* Conda (Miniconda or Anaconda)
* Snakemake


## Installation

1. Clone the repository:

```bash
git clone https://github.com/BPRC-Bioinfo/gene_confirmation_via_raw_reads.git
cd gene_confirmation_via_raw_reads
```

2. Create a Snakemake environment:

```bash
# using conda
conda create -n snakemake -c conda-forge -c bioconda snakemake python=3.10 -y
conda activate snakemake
```

3. Verify snakemake is available:

```bash
snakemake --version
```

## Generating primary BAM files

If coordinate-sorted primary BAM files have not already been generated, raw reads can be aligned to the haplotype assemblies using minimap2.

Both haplotypes should be included in the same alignment reference so that each read can map competitively to the most appropriate haplotype.

For example:


```
# Combine haplotypes into a single FASTA
cat sample_hap1.fasta sample_hap2.fasta > sample_haps.fasta

# Map raw reads and retain primary alignments
minimap2 -ax {preset} --secondary=no sample_haps.fasta raw_reads.fastq \
    | samtools view -u -F 0x900 \
    | samtools sort -o {sample}_{tech}_merged_sorted_primary.bam

# Create BAM index
samtools index {sample}_{tech}_merged_sorted_primary.bam
```

Use:


| Sequencing technology | minimap2 preset | {tech} value |
|-----------------------|-----------------|--------------|
| PacBio HiFi | map-hifi | pacbio |
| Oxford Nanopore | map-ont | nanopore |

--secondary=no prevents secondary alignments from being reported by minimap2.

The samtools view -F 0x900 step removes alignments carrying the secondary (0x100) or supplementary (0x800) SAM flags, leaving primary alignments for downstream analysis.


## Configure inputs

Edit `config.yaml` and provide the paths to the required input data.


Example:   

```yaml
bam_dir: /path/to/aligned_bam/   # directory containing sample BAMs (indexed .bam and .bai)
region_file: region_coordinates.bed
input_excel: your_annotation_file.xlsx # Excel output from VDJ-Insights
fasta_dir: /path/to/haps/        # directory with FASTA files referenced by the pipeline
```

### BAM files

`bam_dir` should point to the top-level directory containing the aligned BAM files for the samples.

Each BAM file must:

* be coordinate sorted;
* contain primary alignments;
* have a corresponding BAM index (.bai);
* use contig names that match the corresponding FASTA files.

The preferred BAM filename pattern is:

```
{sample_id}_{tech}_merged_sorted_primary.bam
{sample_id}_{tech}_merged_sorted_primary.bam.bai
```

where {tech} is either `pacbio` or `nanopore`.

For example:

```
sample42_pacbio_merged_sorted_primary.bam
sample42_pacbio_merged_sorted_primary.bam.bai
sample42_nanopore_merged_sorted_primary.bam
sample42_nanopore_merged_sorted_primary.bam.bai
```


### Region file

`region_file` should be a BED-like tab-separated file containing the interested genomic regions.


Expected structure:


```
contig_name    region_start_coord    region_end_coord    region_name
```

Example:


```
hap1_GCA_049350105.2_CM111665_2_f1    94425001    95632609    GCA_049350105.2_haps_hap1_TRA
hap1_T2TMmul8                           89238939    90297198    T2T-MFA8v1.1_haps_hap1_TRB
```

The region_name should end in one of the currently supported receptor loci:

* TRA — T cell receptor alpha
* TRB — T cell receptor beta
* TRG — T cell receptor gamma
* TRC - Gene segment constant

The contig names in this file must match the contig names present in both the BAM and FASTA files.


### Annotation Excel file

`input_excel` should be an Excel output file from either:

* VDJ-Insights
* Annotation refinement script

The annotation file provides the gene identities, coordinates, strand information, and target sequences used for raw-read confirmation.

### FASTA files

`fasta_dir` should contain the haplotypes FASTA sequences corresponding to the samples being analysed.

FASTA contig names must match those used in:

* the BAM header;
* the region file;
* the annotation coordinates.


## Run the pipeline

From the repository root:

```bash
snakemake -s coordinate_pipe.smk --use-conda --cores 8
```

## Outputs

The main user-facing outputs are the Excel summary files:


```
results/v_region_summary.xlsx
results/d_region_summary.xlsx
results/j_region_summary.xlsx
results/c_region_summary.xlsx
```

The workflow produces per-region output Excel files with the following columns:

| Column | Description |
|--------|-------------|
| Sample | Sample identifier |
| Short name | Short sample name |
| No. reads PB | Number of PacBio reads |
| No. reads ONT | Number of Oxford Nanopore reads |
| No. reads PB 100% | Number of PacBio reads with 100% identity and 100% target coverage |
| No. reads ONT 100% | Number of Oxford Nanopore reads with 100% identity and 100% target coverage |
| Region | Genomic region |
| Segment | Gene segment (V, D, J or C ) |
| Start coord | Start coordinate |
| End coord | End coordinate |
| Strand | DNA strand (+/-) |
| Target length | Length of target sequence |
| Target seq | Target sequence |
| SNP Positions PB | SNP positions in PacBio reads |
| SNP Positions ONT | SNP positions in Oxford Nanopore reads |

Output files are generated per segment (v, d, j or c).


## Intermediate output files

In addition to the Excel summary files, the workflow generates intermediate files including:

* run tables;
* per-job TSV files;
* gene/read matrices;
* read-vs-reference summaries;
* position-level base summaries;
* completion marker files;
* log files.

These files are useful for inspecting individual genes, evaluating read-level support, and troubleshooting unexpected results.


In addition to the Excel summary files, the pipeline generates intermediate files such as run tables, per-job TSV files, read-level result files, completion files, and log files. These files can be useful for inspecting individual genes, investigating read-level evidence, and troubleshooting unexpected results.


### Checking raw reads

For detailed inspection of an individual gene, the most useful files are:

```
results/*region_results/*gene_read_matrix.tsv
results/*region_results/*read_vs_reference.tsv
results/*region_results/*position_base_percent.tsv
```

#### Gene/read matrix

The file: `results/*region_results/*gene_read_matrix.tsv`

contains the target/reference sequence and the nucleotide call from every overlapping read at every target position.


```
read_id     gene_pos_1  gene_pos_2  gene_pos_3  gene_pos_4  gene_pos_5
REFERENCE   A           T           A           C           T
read1       A           T           A           C           A
read2       A           T           A           C           T
read3       A           T           A           -           T
read4       .           .           A           C           T
```

The symbols have the following meanings:

| Symbol | Meaning |
|--------|---------|
| A, C, G, T | Nucleotide observed in the read |
| N | Ambiguous nucleotide call in the read |
| - | Deletion relative to the reference at this position |
| . | The read does not cover this reference position |


In the above example:

* read1 contains a base (A) that does not match the reference

* read3 spans position 4, but the aligned read contains a deletion relative to the reference.

* positions 1 and 2 are outside the portion of the target covered by read4.


#### Read-vs-reference summary

The file: `results/*region_results/*read_vs_reference.tsv`

contains one row per retained read.

```
read_id    covered_positions    matches_to_reference    mismatches_to_reference    deletions    percent_identity_vs_reference    percent_gene_covered
read1      5                    4                       1                          0            80.00                             100.00
read2      5                    5                       0                          0            100.00                            100.00
read3      5                    4                       0                          1            80.00                             100.00
read4      3                    3                       0                          0            100.00                             60.00
```

The columns are interpreted as follows:

| Column | Description |
|--------|-------------|
| read_id | Read identifier |
| covered_positions | Number of target positions represented by a nucleotide or deletion in the alignment |
| matches_to_reference | Number of covered positions matching the reference |
| mismatches_to_reference | Number of covered positions containing a nucleotide different from the reference |
| deletions | Number of covered target positions represented by - |
| percent_identity_vs_reference | Percentage of covered positions that exactly match the reference |
| percent_gene_covered | Percentage of the complete target region covered by the read alignment |

Positions represented by . are not counted as covered positions.

Deletions represented by - are counted as covered positions, but they do not count as matches.

* read4 perfectly matches the reference sequence wherever it is aligned, so its identity is 100%. However, it covers only 3 of the 5 target positions, so its target coverage is 60%.

```
gene_pos  ref_base  depth_covered  A_count  C_count  G_count  T_count  N_count  DEL_count  A_percent  C_percent  G_percent  T_percent  N_percent  DEL_percent  major_base  major_base_percent
1         A         3              3        0        0        0        0        0          100.00     0.00       0.00       0.00       0.00       0.00         A           100.00
2         T         3              0        0        0        3        0        0            0.00     0.00       0.00     100.00       0.00       0.00         T           100.00
3         A         4              4        0        0        0        0        0          100.00     0.00       0.00       0.00       0.00       0.00         A           100.00
4         C         4              0        3        0        0        0        1            0.00    75.00       0.00       0.00       0.00      25.00         C            75.00
5         T         4              1        0        0        3        0        0           25.00     0.00       0.00      75.00       0.00       0.00         T            75.00
```


| Column | Description |
|--------|-------------|
| gene_pos | Position within the target sequence, starting at 1 |
| chrom | Reference contig |
| ref_pos_1based | Genomic reference coordinate |
| ref_base | Reference nucleotide |
| depth_covered | Number of reads covering the position |
| A_count | Number of reads containing A |
| C_count | Number of reads containing C |
| G_count | Number of reads containing G |
| T_count | Number of reads containing T |
| N_count | Number of reads containing N |
| DEL_count | Number of reads containing a deletion at the position |
| A_percent | Percentage of covered reads containing A |
| C_percent | Percentage of covered reads containing C |
| G_percent | Percentage of covered reads containing G |
| T_percent | Percentage of covered reads containing T |
| N_percent | Percentage of covered reads containing N |
| DEL_percent | Percentage of covered reads containing a deletion |
| major_base | Most frequently observed call at that position |
| major_base_percent | Percentage of covered reads containing the major call |



### Interpreting SNP positions

The SNP Positions PB and SNP Positions ONT columns provide a compact summary of base observations at positions of interest within the target sequence.

```
position:BASEpercent;BASEpercent
```


An example may look like:

```
5:A90.00;T10.00
```

Means that at position 5:

* 90% of the reads contain A
* 10% of the reads contain T

If the target/reference sequence contains A at this position, the T observations indicate a possible sequence difference in a subset of the reads. Multiple positions are separated by |

For example:

```
3:A88.89;G11.11 | 25:C77.78;T11.11;N11.11 | 26:A11.11;T88.89 | 27:C88.89;N11.11
```

#### Interpreting N

N represents an ambiguous nucleotide call. It means that due to no coverage the base sequence at that aligned position is reported as N rather than a specific nucleotide (A, C, G, or T).


### Interpreting PacBio and Oxford Nanopore evidence

PacBio HiFi and Oxford Nanopore sequencing technologies have different sequencing characteristics and error profiles.

When the same specific alternative nucleotide is independently observed in both PacBio and ONT reads, this provides stronger evidence that the observation represents a genuine biological sequence difference rather than a technology-specific sequencing error.

Interpretation should also consider:

* the number of reads supporting the alternative base;
* total read depth at the position;
* whether the same specific alternative occurs repeatedly in independent reads;
* read sequence identity;
* percentage of the target covered by the reads;
* presence of deletions or alignment ambiguity;
* the number of N calls at the position;
* whether the position lies in a difficult or repetitive genomic region;
* availability and depth of the other sequencing technology.

Low read depth should always be interpreted cautiously.


