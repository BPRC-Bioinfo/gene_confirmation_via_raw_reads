import re
import pandas as pd


def clean_value(value):
    """
    Convert empty/Error values to n/a.
    """
    if value is None:
        return "n/a"

    value = value.strip()

    if value == "" or value.lower() == "error":
        return "n/a"

    return value


def is_real_value(value):
    """
    True only for values that should be allowed to overwrite existing data.
    """
    return value not in [None, "", "n/a"]


rows = {}

for done_file in snakemake.input.done_files:
    sample = None
    shortname = None
    target_seq = "n/a"
    target_length = "n/a"
    snp_positions = "n/a"
    region = "n/a"
    segment = "n/a"
    start_coord = "n/a"
    end_coord = "n/a"
    strand = "n/a"

    pb_reads = "n/a"
    ont_reads = "n/a"
    pb_reads_100 = "n/a"
    ont_reads_100 = "n/a"

    with open(done_file) as fh:
        for line in fh:
            line = line.strip()

            if line.startswith("Sample:"):
                sample = clean_value(line.split(":", 1)[1])

            elif line.startswith("Short name:"):
                shortname = clean_value(line.split(":", 1)[1])

            elif line.startswith("Target ref:"):
                target_seq = clean_value(line.split(":", 1)[1])

            elif line.startswith("Target length:"):
                target_length = clean_value(line.split(":", 1)[1])

            elif line.startswith("SNP Positions:"):
                snp_positions = clean_value(line.split(":", 1)[1])

            elif line.startswith("Region:"):
                region = clean_value(line.split(":", 1)[1])

            elif line.startswith("Segment:"):
                segment = clean_value(line.split(":", 1)[1])

            elif line.startswith("Start coord:"):
                start_coord = clean_value(line.split(":", 1)[1])

            elif line.startswith("End coord:"):
                end_coord = clean_value(line.split(":", 1)[1])

            elif line.startswith("Strand:"):
                strand = clean_value(line.split(":", 1)[1])

            elif line.startswith("No. reads"):
                match = re.match(
                    r"No\. reads\s+(\S+)\s*(100%)?\s*:\s*(.+)",
                    line
                )

                if match:
                    sequencer = match.group(1).lower()
                    is_100 = match.group(2) is not None
                    value = clean_value(match.group(3))

                    if sequencer in ["pacbio", "pb"]:
                        if is_100:
                            pb_reads_100 = value
                        else:
                            pb_reads = value

                    elif sequencer in ["nanopore", "ont"]:
                        if is_100:
                            ont_reads_100 = value
                        else:
                            ont_reads = value

    key = (sample, shortname)

    if key not in rows:
        rows[key] = {
            "Sample": sample,
            "Short name": shortname,
            "No reads PB": "n/a",
            "No reads ONT": "n/a",
            "No reads PB 100%": "n/a",
            "No reads ONT 100%": "n/a",
            "Region": "n/a",
            "Segment": "n/a",
            "Start coord": "n/a",
            "End coord": "n/a",
            "Strand": "n/a",
            "Target length": "n/a",
            "Target seq": "n/a",
            "SNP Positions PB": "n/a",
            "SNP Positions ONT": "n/a",
        }

    if is_real_value(pb_reads):
        rows[key]["No reads PB"] = pb_reads

    if is_real_value(ont_reads):
        rows[key]["No reads ONT"] = ont_reads

    if is_real_value(pb_reads_100):
        rows[key]["No reads PB 100%"] = pb_reads_100

    if is_real_value(ont_reads_100):
        rows[key]["No reads ONT 100%"] = ont_reads_100

    if is_real_value(region):
        rows[key]["Region"] = region

    if is_real_value(segment):
        rows[key]["Segment"] = segment

    if is_real_value(start_coord):
        rows[key]["Start coord"] = start_coord

    if is_real_value(end_coord):
        rows[key]["End coord"] = end_coord

    if is_real_value(strand):
        rows[key]["Strand"] = strand

    if is_real_value(target_length):
        rows[key]["Target length"] = target_length

    if is_real_value(target_seq):
        rows[key]["Target seq"] = target_seq

    if is_real_value(snp_positions):
        if is_real_value(pb_reads) or is_real_value(pb_reads_100):
            rows[key]["SNP Positions PB"] = snp_positions

        if is_real_value(ont_reads) or is_real_value(ont_reads_100):
            rows[key]["SNP Positions ONT"] = snp_positions


df = pd.DataFrame(rows.values())

df = df[
    [
        "Sample",
        "Short name",
        "No reads PB",
        "No reads ONT",
        "No reads PB 100%",
        "No reads ONT 100%",
        "Region",
        "Segment",
        "Start coord",
        "End coord",
        "Strand",
        "Target length",
        "Target seq",
        "SNP Positions PB",
        "SNP Positions ONT",
    ]
]

df.to_excel(snakemake.output.xlsx, index=False)