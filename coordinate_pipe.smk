
# Ver 0.0.4

# Giang le

import os
import glob
from pathlib import Path

import pandas as pd


configfile: "config/config.yaml"

COORDINATE_TABLE = "rhesus_genes_coordinates_from_bam_bed.txt"

SEGMENTS = ["V", "D", "J"]

FASTA_DIR = config.get("fasta_dir", "/mnt/CGR_Analyse/HPRC_v2/raw-reads/raw_reads/haps/")

print (FASTA_DIR)

def segment_lower(segment):
    return segment.lower()


def run_table_for_segment(segment):
    return f"{segment_lower(segment)}_region_run_table.tsv"


def split_dir_for_segment(segment):
    return f"{segment_lower(segment)}_region_jobs"


def result_dir_for_segment(segment):
    return f"{segment_lower(segment)}_region_results"


def summary_xlsx_for_segment(segment):
    return f"{segment_lower(segment)}_region_summary.xlsx"


def collect_region_results(wildcards):
    segment = wildcards.segment
    job_ids = get_region_job_ids(wildcards)

    return expand(
        os.path.join(result_dir_for_segment(segment), "{job_id}.done"),
        job_id=job_ids,
    )


rule all:
    input:
#        "results/coordinates_from_bam_bed.txt"
        summaries = expand(
            "results/{segment_lower}_region_summary.xlsx",
            segment_lower=[segment.lower() for segment in SEGMENTS],
        ),
        done_files = expand(
            "results/{segment_lower}_region_results/.complete",
            segment_lower=[segment.lower() for segment in SEGMENTS],
        )


rule coordinate_finder:
    input:
        bed = config["region_file"],
        bam_dir = config["bam_dir"],
        excel = config["input_excel"],
    output:
        "results/coordinates_from_bam_bed.txt"
    shell:
        """
        python scripts/coordinate_study.py \
            --regions {input.bed} \
            --bam-dir {input.bam_dir} \
            --input-excel {input.excel} \
            --output {output}
        """


rule make_region_run_table:
    input:
        rules.coordinate_finder.output
    output:
        "results/{segment_lower}_region_run_table.tsv"
    params:
        fasta_dir = FASTA_DIR
    run:
        segment = wildcards.segment_lower.upper()

        if segment not in SEGMENTS:
            raise ValueError(
                f"Unsupported segment '{segment}'. Expected one of: {SEGMENTS}"
            )

        df = pd.read_csv(input[0], sep="\t")

        # In the coordinate output:
        # lowercase region = TRA/TRB/TRG
        # uppercase Region = Excel Segment column
        region_df = df[df["Segment"].astype(str) == segment].copy()

        rows = []

        for idx, row in region_df.reset_index(drop=True).iterrows():
            sample_id = str(row["sample_id"])
            haplotype = str(row["haplotype"])
            sample_excel = str(row["excel_sample"])
            shortname = str(row["Short name"])

            fasta = Path(params.fasta_dir) / f"{sample_id}_{haplotype}.fa"

            row_id = f"{segment.lower()}_{idx:06d}"

            base = {
                "region_row_id": row_id,
                "sample_id": sample_id,
                "haplotype": haplotype,
                "Short name": shortname,
                "sample_excel": sample_excel,
                "locus_region": row["region"],
                "segment": row["Segment"],
                "region_id": row["region_id"],
                "strand": row["Strand"],
                "start_coord": int(row["bam_invest_start_coord"]),
                "end_coord": int(row["bam_invest_end_coord"]),
                "fasta": str(fasta),
            }

            rows.append({
                **base,
                "job_id": f"{row_id}_pacbio",
                "tech": "pacbio",
                "bam": row["pacbio_bam"],
                "output": os.path.join(
                    result_dir_for_segment(segment),
                    f"{row_id}.pacbio.txt",
                ),
            })

            rows.append({
                **base,
                "job_id": f"{row_id}_nanopore",
                "tech": "nanopore",
                "bam": row["nanopore_bam"],
                "output": os.path.join(
                    result_dir_for_segment(segment),
                    f"{row_id}.nanopore.txt",
                ),
            })

        run_table = pd.DataFrame(rows)

        run_table.to_csv(
            output[0],
            sep="\t",
            index=False,
        )


checkpoint split_region_run_table:
    input:
        "{segment_lower}_region_run_table.tsv"
    output:
        directory("{segment_lower}_region_jobs")
    run:
        os.makedirs(output[0], exist_ok=True)

        df = pd.read_csv(input[0], sep="\t")

        for _, row in df.iterrows():
            job_id = row["job_id"]
            out_file = os.path.join(output[0], f"{job_id}.tsv")

            row.to_frame().T.to_csv(
                out_file,
                sep="\t",
                index=False,
            )


def get_region_job_ids(wildcards):
    checkpoint_output = checkpoints.split_region_run_table.get(
        segment_lower=wildcards.segment
    ).output[0]

    job_files = glob.glob(
        os.path.join(checkpoint_output, "*.tsv")
    )

    job_ids = [
        os.path.basename(job_file).replace(".tsv", "")
        for job_file in job_files
    ]

    return sorted(job_ids)


rule run_region:
    input:
        job = "{segment}_region_jobs/{job_id}.tsv"
    output:
        done = "{segment}_region_results/{job_id}.done"
    log:
        "{segment}_region_results/{job_id}.log"
    conda:
        "envs/gene_matrix.yaml"
    params:
        outdir = "{segment}_region_results"
    shell:
        r"""
        mkdir -p "{params.outdir}"

        bam=$(awk -F '\t' 'NR==2 {{print $15}}' "{input.job}")
        fasta=$(awk -F '\t' 'NR==2 {{print $12}}' "{input.job}")
        contig=$(awk -F '\t' 'NR==2 {{print $8}}' "{input.job}")
        start_coord=$(awk -F '\t' 'NR==2 {{print $10}}' "{input.job}")
        end_coord=$(awk -F '\t' 'NR==2 {{print $11}}' "{input.job}")
        output_file=$(awk -F '\t' 'NR==2 {{print "results/"$16}}' "{input.job}")
        sample_excel=$(awk -F '\t' 'NR==2 {{print $5}}' "{input.job}")
        shortname=$(awk -F '\t' 'NR==2 {{print $4}}' "{input.job}")
        strand=$(awk -F '\t' 'NR==2 {{print $9}}' "{input.job}")
        region=$(awk -F '\t' 'NR==2 {{print $6}}' "{input.job}")
        segment=$(awk -F '\t' 'NR==2 {{print $7}}' "{input.job}")


        echo "Job file: {input.job}" > "{log}"
        echo "BAM: $bam" >> "{log}"
        echo "FASTA: $fasta" >> "{log}"
        echo "Contig: $contig" >> "{log}"
        echo "Start coord: $start_coord" >> "{log}"
        echo "End coord: $end_coord" >> "{log}"
        echo "Output prefix: $output_file" >> "{log}"
        echo "python bam_gene_matrix4.py $bam $fasta $contig $start_coord $end_coord $output_file" >> "{log}"

        script_status="OK"

        python scripts/bam_gene_matrix.py \
            "$bam" \
            "$fasta" \
            "$contig" \
            "$start_coord" \
            "$end_coord" \
            "$output_file" >> "{log}" 2>> "{log}" || script_status="Error"

        read_summary_file="${{output_file}}.read_vs_reference.tsv"

        if [ "$script_status" = "OK" ] && [ -s "$read_summary_file" ]; then
            sequencer=$(echo "{wildcards.job_id}" | sed 's/.*_//g')
            total_reads=$(sed '1d' "$read_summary_file" | wc -l)
            qualify_reads=$(awk 'NR > 1 && $6 == 100 && $7 == 100' "$read_summary_file" | wc -l)

            ref=$(grep REF "${{output_file}}.gene_read_matrix.tsv" | cut -f2- | sed 's/\t//g')
            ref_len=$(echo "$ref" | wc -c | awk '{{print $0 - 1}}')

            position_summary_file="${{output_file}}.position_base_percent.tsv"
            variable_sites_file="${{output_file}}.variable_sites.txt"

            if [ -s "$position_summary_file" ]; then
                python scripts/summarize_variable_positions.py \
                    "$position_summary_file" \
                    "$variable_sites_file" >> "{log}" 2>&1

                if [ -s "$variable_sites_file" ]; then
                    positions=$(cat "$variable_sites_file")
                else
                    positions=""
                fi
            else
                positions="Error"
                echo "ERROR: position summary file was not created: $position_summary_file" >> "{log}"
            fi
        else
            sequencer=$(echo "{wildcards.job_id}" | sed 's/.*_//g')
            total_reads="Error"
            qualify_reads="Error"
            ref="Error"
            ref_len="Error"
            positions="Error"
            echo "ERROR: bam_gene_matrix4.py failed or output file was not created: $read_summary_file" >> "{log}"
        fi

        echo "Sample: $sample_excel" > "{output.done}"
        echo "Short name: $shortname" >> "{output.done}"
        echo "No. reads $sequencer : $total_reads" >> "{output.done}"
        echo "No. reads $sequencer 100% : $qualify_reads" >> "{output.done}"
        echo "Target ref: $ref" >> "{output.done}"
        echo "Start coord: $start_coord" >> "{output.done}"
        echo "End coord: $end_coord" >> "{output.done}"
        echo "Strand: $strand" >> "{output.done}"
        echo "Region: $region" >> "{output.done}"
        echo "Segment: $segment" >> "{output.done}"
        echo "Target length: $ref_len" >> "{output.done}"
        echo "SNP Positions: $positions" >> "{output.done}"
        """


rule mark_region_complete:
    input:
        done_files = collect_region_results
    output:
        complete = "{segment}_region_results/.complete"
    shell:
        """
        mkdir -p $(dirname "{output.complete}")
        touch "{output.complete}"
        """


rule combine_region_done:
    input:
        done_files = collect_region_results
    output:
        xlsx = "{segment}_region_summary.xlsx"
    conda:
        "envs/gene_matrix.yaml"
    script:
        "scripts/combine_v_region_done.py"