
# Ver 0.0.4
# Giang le


#!/usr/bin/env python3

import sys
import os
from collections import Counter
import difflib

import pysam


def usage():
    print(
        "\nUsage:\n"
        "  python bam_gene_matrix4.py <input.bam> <reference.fa> <chrom> "
        "<start_1based> <end_1based> <output_prefix>\n\n"
        "Example:\n"
        "  python bam_gene_matrix4.py rin14_hap1_TRA.bam RiN14_hap1_TRA_contig.fa "
        "hap1_RiN14_CM111665_2 94610322 94610597 rin14_TRA\n"
    )
    sys.exit(1)


def print_contig_diagnostics(chrom, bam, ref, bam_path, ref_fasta_path, start0, end0):
    bam_abs_path = os.path.realpath(bam_path)
    fasta_abs_path = os.path.realpath(ref_fasta_path)

    bam_contigs = list(bam.references)
    fasta_contigs = list(ref.references)

    print()
    print("ERROR: Requested contig was not found or could not be fetched.")
    print()

    print("Input files:")
    print(f"  BAM:   {bam_abs_path}")
    print(f"  FASTA: {fasta_abs_path}")
    print()

    print("Requested region:")
    print(f"  Contig:        {chrom}")
    print(f"  Start 1-based: {start0 + 1}")
    print(f"  End 1-based:   {end0}")
    print()

    print("Contig present in BAM?")
    print(f"  {'YES' if chrom in bam_contigs else 'NO'}")
    print()

    print("Contig present in FASTA?")
    print(f"  {'YES' if chrom in fasta_contigs else 'NO'}")
    print()

    bam_matches = difflib.get_close_matches(chrom, bam_contigs, n=20, cutoff=0.2)
    fasta_matches = difflib.get_close_matches(chrom, fasta_contigs, n=20, cutoff=0.2)

    print("Closest matching BAM contigs:")
    if bam_matches:
        for contig in bam_matches:
            print(f"  {contig}")
    else:
        print("  No close BAM contig matches found.")
    print()

    print("Closest matching FASTA contigs:")
    if fasta_matches:
        for contig in fasta_matches:
            print(f"  {contig}")
    else:
        print("  No close FASTA contig matches found.")
    print()

    print("All BAM contigs:")
    for contig in bam_contigs:
        print(f"  {contig}")
    print()

    print("All FASTA contigs:")
    for contig in fasta_contigs:
        print(f"  {contig}")
    print()

    sys.exit(2)


def validate_contig(chrom, bam, ref, bam_path, ref_fasta_path, start0, end0):
    bam_contigs = set(bam.references)
    fasta_contigs = set(ref.references)

    if chrom not in bam_contigs or chrom not in fasta_contigs:
        print_contig_diagnostics(
            chrom,
            bam,
            ref,
            bam_path,
            ref_fasta_path,
            start0,
            end0,
        )


def get_read_name(read):
    read_name = read.query_name

    if read.is_read1:
        read_name += "/1"
    elif read.is_read2:
        read_name += "/2"

    return read_name


def build_read_calls(read, seq, start0, end0, gene_len):
    calls = ["."] * gene_len

    for query_pos, ref_pos in read.get_aligned_pairs(matches_only=False):
        if ref_pos is None:
            continue

        if start0 <= ref_pos < end0:
            gene_index = ref_pos - start0

            if query_pos is None:
                calls[gene_index] = "-"
            else:
                calls[gene_index] = seq[query_pos].upper()

    return calls


def summarize_read(calls, ref_seq):
    covered_positions = 0
    matches = 0
    mismatches = 0
    deletions = 0

    for i, base in enumerate(calls):
        if base == ".":
            continue

        covered_positions += 1
        ref_base = ref_seq[i]

        if base == "-":
            deletions += 1
        elif base == ref_base:
            matches += 1
        else:
            mismatches += 1

    return covered_positions, matches, mismatches, deletions


def write_position_summary(
    output_path,
    all_read_calls,
    ref_seq,
    chrom,
    gene_start_1based,
    gene_len,
):
    with open(output_path, "w") as pos_fh:
        pos_fh.write(
            "gene_pos\t"
            "chrom\t"
            "ref_pos_1based\t"
            "ref_base\t"
            "depth_covered\t"
            "A_count\tC_count\tG_count\tT_count\tN_count\tDEL_count\t"
            "A_percent\tC_percent\tG_percent\tT_percent\tN_percent\tDEL_percent\t"
            "major_base\tmajor_base_percent\n"
        )

        for i in range(gene_len):
            ref_pos_1based = gene_start_1based + i
            ref_base = ref_seq[i]

            bases_at_position = [
                calls[i]
                for calls in all_read_calls
                if calls[i] != "."
            ]

            counts = Counter(bases_at_position)
            depth = len(bases_at_position)

            if depth > 0:
                a_pct = counts["A"] / depth * 100
                c_pct = counts["C"] / depth * 100
                g_pct = counts["G"] / depth * 100
                t_pct = counts["T"] / depth * 100
                n_pct = counts["N"] / depth * 100
                del_pct = counts["-"] / depth * 100

                major_base, major_count = counts.most_common(1)[0]
                major_pct = major_count / depth * 100
            else:
                a_pct = c_pct = g_pct = t_pct = n_pct = del_pct = 0.0
                major_base = "."
                major_pct = 0.0

            pos_fh.write(
                f"{i + 1}\t"
                f"{chrom}\t"
                f"{ref_pos_1based}\t"
                f"{ref_base}\t"
                f"{depth}\t"
                f"{counts['A']}\t"
                f"{counts['C']}\t"
                f"{counts['G']}\t"
                f"{counts['T']}\t"
                f"{counts['N']}\t"
                f"{counts['-']}\t"
                f"{a_pct:.2f}\t"
                f"{c_pct:.2f}\t"
                f"{g_pct:.2f}\t"
                f"{t_pct:.2f}\t"
                f"{n_pct:.2f}\t"
                f"{del_pct:.2f}\t"
                f"{major_base}\t"
                f"{major_pct:.2f}\n"
            )


def main():
    if len(sys.argv) != 7:
        usage()

    bam_path = sys.argv[1]
    ref_fasta_path = sys.argv[2]
    chrom = sys.argv[3]
    gene_start_1based = int(sys.argv[4])
    gene_end_1based = int(sys.argv[5])
    output_prefix = sys.argv[6]

    start0 = gene_start_1based
    end0 = gene_end_1based
    gene_len = end0 - start0

    if gene_len <= 0:
        print()
        print("ERROR: End coordinate must be greater than or equal to start coordinate.")
        print()
        print(f"Start coordinate: {gene_start_1based}")
        print(f"End coordinate:   {gene_end_1based}")
        print()
        sys.exit(1)

    matrix_out = f"{output_prefix}.gene_read_matrix.tsv"
    read_summary_out = f"{output_prefix}.read_vs_reference.tsv"
    position_summary_out = f"{output_prefix}.position_base_percent.tsv"

    all_read_calls = []

    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam, pysam.FastaFile(ref_fasta_path) as ref:
            validate_contig(
                chrom,
                bam,
                ref,
                bam_path,
                ref_fasta_path,
                start0,
                end0,
            )

            try:
                ref_seq = ref.fetch(chrom, start0, end0).upper()
            except ValueError:
                print_contig_diagnostics(
                    chrom,
                    bam,
                    ref,
                    bam_path,
                    ref_fasta_path,
                    start0,
                    end0,
                )

            with open(matrix_out, "w") as matrix_fh, open(read_summary_out, "w") as read_fh:
                matrix_fh.write(
                    "read_id\t"
                    + "\t".join([f"gene_pos_{i + 1}" for i in range(gene_len)])
                    + "\n"
                )

                matrix_fh.write("REFERENCE\t" + "\t".join(ref_seq) + "\n")

                read_fh.write(
                    "read_id\t"
                    "covered_positions\t"
                    "matches_to_reference\t"
                    "mismatches_to_reference\t"
                    "deletions\t"
                    "percent_identity_vs_reference\t"
                    "percent_gene_covered\n"
                )

                try:
                    reads = bam.fetch(chrom, start0, end0)
                except ValueError:
                    print_contig_diagnostics(
                        chrom,
                        bam,
                        ref,
                        bam_path,
                        ref_fasta_path,
                        start0,
                        end0,
                    )

                for read in reads:
                    if read.is_unmapped or read.query_sequence is None:
                        continue

                    read_name = get_read_name(read)

                    calls = build_read_calls(
                        read,
                        read.query_sequence,
                        start0,
                        end0,
                        gene_len,
                    )

                    matrix_fh.write(read_name + "\t" + "\t".join(calls) + "\n")
                    all_read_calls.append(calls)

                    covered, matches, mismatches, deletions = summarize_read(calls, ref_seq)

                    percent_identity = (
                        matches / covered * 100
                        if covered > 0
                        else 0.0
                    )

                    percent_gene_covered = covered / gene_len * 100

                    read_fh.write(
                        f"{read_name}\t"
                        f"{covered}\t"
                        f"{matches}\t"
                        f"{mismatches}\t"
                        f"{deletions}\t"
                        f"{percent_identity:.2f}\t"
                        f"{percent_gene_covered:.2f}\n"
                    )

        write_position_summary(
            position_summary_out,
            all_read_calls,
            ref_seq,
            chrom,
            gene_start_1based,
            gene_len,
        )

        print("SUCCESS")
        print(f"Wrote: {matrix_out}")
        print(f"Wrote: {read_summary_out}")
        print(f"Wrote: {position_summary_out}")

    except FileNotFoundError as e:
        print()
        print("ERROR: Input or output path was not found.")
        print()
        print(f"Details: {e}")
        print()
        print("Input files:")
        print(f"  BAM:   {os.path.realpath(bam_path)}")
        print(f"  FASTA: {os.path.realpath(ref_fasta_path)}")
        print()
        sys.exit(1)

    except OSError as e:
        print()
        print("ERROR: File or index problem.")
        print()
        print(f"Details: {e}")
        print()
        print("Possible causes:")
        print("  1. BAM file is missing.")
        print("  2. BAM index is missing.")
        print("  3. FASTA file is missing.")
        print("  4. FASTA index .fai is missing.")
        print("  5. Output directory does not exist.")
        print()
        print("Input files:")
        print(f"  BAM:   {os.path.realpath(bam_path)}")
        print(f"  FASTA: {os.path.realpath(ref_fasta_path)}")
        print()
        sys.exit(1)


if __name__ == "__main__":
    main()