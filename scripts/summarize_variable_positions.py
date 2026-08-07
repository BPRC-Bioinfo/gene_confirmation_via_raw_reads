#!/usr/bin/env python3

import sys
import pandas as pd


if len(sys.argv) != 3:
    sys.stderr.write(
        "\nUsage:\n"
        "  python summarize_variable_positions.py <position_base_percent.tsv> <output.txt>\n\n"
        "Example:\n"
        "  python summarize_variable_positions.py sample.position_base_percent.tsv sample.variable_sites.txt\n\n"
    )
    sys.exit(1)


input_tsv = sys.argv[1]
output_txt = sys.argv[2]

df = pd.read_csv(input_tsv, sep="\t")

variable_df = df[df["major_base_percent"] < 100].copy()

coded_positions = []

base_percent_columns = [
    ("A", "A_percent"),
    ("C", "C_percent"),
    ("G", "G_percent"),
    ("T", "T_percent"),
    ("N", "N_percent"),
    ("N", "DEL_percent"),  # encode deletions as N
]

for _, row in variable_df.iterrows():
    gene_pos = int(row["gene_pos"])
    base_parts = []

    for base, percent_col in base_percent_columns:
        percent = float(row[percent_col])

        if percent > 0:
            base_parts.append(f"{base}{percent:.2f}")

    coded_position = f"{gene_pos}:" + ";".join(base_parts)
    coded_positions.append(coded_position)

result = " | ".join(coded_positions)

with open(output_txt, "w") as out_fh:
    out_fh.write(result + "\n")