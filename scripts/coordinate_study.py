# ver 0.0.5
# Giang Le

#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


TARGET_REGIONS = ["TRA", "TRB", "TRG"]

GENE_COORD_COLUMNS = [
    "Sample",
    "Short name",
    "Region",
    "Segment",
    "Start coord",
    "End coord",
    "Strand",
]


def read_region_coordinates(path: str) -> pd.DataFrame:
    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=[
            "region_id",
            "region_start",
            "region_end",
            "region_name",
        ],
    )


def read_input_excel(path: str) -> pd.DataFrame:
    return pd.read_excel(path, usecols=GENE_COORD_COLUMNS)


def parse_region_id(region_id: str) -> tuple[str, str]:
    parts = str(region_id).split("_")

    if len(parts) < 2:
        raise ValueError(
            f"Could not parse region_id: {region_id}. "
            "Expected format like hap1_R04081_..."
        )

    haplotype = parts[0]
    sample_id = parts[1].lower()

    return haplotype, sample_id


def add_parsed_region_id_columns(regions: pd.DataFrame) -> pd.DataFrame:
    regions = regions.copy()

    parsed = regions["region_id"].apply(parse_region_id)
    regions["haplotype"] = parsed.apply(lambda x: x[0])
    regions["sample_id"] = parsed.apply(lambda x: x[1])

    return regions


def get_unique_region_ids(regions: pd.DataFrame) -> pd.DataFrame:
    unique_regions = regions[["region_id"]].drop_duplicates().reset_index(drop=True)
    return add_parsed_region_id_columns(unique_regions)


def extract_region_from_region_name(region_name: str) -> str | None:
    region_name = str(region_name)

    for region in TARGET_REGIONS:
        if region_name.endswith(f"_{region}"):
            return region

    return None


def find_matching_sample_folders(bam_dir: str, sample_id: str) -> list[Path]:
    bam_dir = Path(bam_dir)
    sample_id_lower = sample_id.lower()

    return sorted(
        folder
        for folder in bam_dir.iterdir()
        if folder.is_dir()
        and folder.name.lower().startswith(sample_id_lower)
    )


def find_bam_for_tech(
    matching_folders: list[Path],
    sample_id: str,
    tech: str,
) -> str | None:
    sample_id_lower = sample_id.lower()
    tech_lower = tech.lower()

    expected_prefix = f"{sample_id_lower}_{tech_lower}_merged_sorted_primary"

    matches = []

    for folder in matching_folders:
        for bam_file in folder.rglob("*.bam"):
            bam_name = bam_file.name.lower()

            if bam_name.startswith(expected_prefix):
                matches.append(bam_file)

    if not matches:
        return None

    return str(sorted(matches)[0])


def list_all_bams_in_folders(matching_folders: list[Path]) -> str:
    bam_files = []

    for folder in matching_folders:
        for bam_file in folder.rglob("*.bam"):
            bam_files.append(str(bam_file))

    return ";".join(sorted(bam_files))


def find_available_bams_for_sample(bam_dir: str, sample_id: str) -> dict:
    matching_folders = find_matching_sample_folders(bam_dir, sample_id)

    pacbio_bam = find_bam_for_tech(
        matching_folders=matching_folders,
        sample_id=sample_id,
        tech="pacbio",
    )

    nanopore_bam = find_bam_for_tech(
        matching_folders=matching_folders,
        sample_id=sample_id,
        tech="nanopore",
    )

    has_pacbio = pacbio_bam is not None
    has_nanopore = nanopore_bam is not None
    has_any_bam = has_pacbio or has_nanopore
    has_all_bams = has_pacbio and has_nanopore

    missing = []
    if not has_pacbio:
        missing.append("pacbio")
    if not has_nanopore:
        missing.append("nanopore")

    return {
        "sample_id": sample_id,
        "matching_folders": ";".join(str(folder) for folder in matching_folders),
        "all_bams_in_matching_folders": list_all_bams_in_folders(matching_folders),
        "pacbio_bam": pacbio_bam,
        "nanopore_bam": nanopore_bam,
        "has_pacbio": has_pacbio,
        "has_nanopore": has_nanopore,
        "has_any_bam": has_any_bam,
        "has_all_bams": has_all_bams,
        "missing": ";".join(missing),
    }


def build_bam_presence_table(
    unique_regions: pd.DataFrame,
    bam_dir: str,
) -> pd.DataFrame:
    rows = []

    for sample_id in sorted(unique_regions["sample_id"].unique()):
        rows.append(
            find_available_bams_for_sample(
                bam_dir=bam_dir,
                sample_id=sample_id,
            )
        )

    return pd.DataFrame(rows)


def build_sample_regions_table(
    regions: pd.DataFrame,
    samples_with_bams: pd.DataFrame,
) -> pd.DataFrame:
    regions = add_parsed_region_id_columns(regions)
    regions["region"] = regions["region_name"].apply(extract_region_from_region_name)

    sample_ids_with_bams = set(samples_with_bams["sample_id"])

    sample_regions = regions[
        regions["sample_id"].isin(sample_ids_with_bams)
        & regions["region"].isin(TARGET_REGIONS)
    ].copy()

    sample_regions = sample_regions[
        [
            "sample_id",
            "haplotype",
            "region",
            "region_id",
            "region_start",
            "region_end",
            "region_name",
        ]
    ]

    return sample_regions.sort_values(
        by=[
            "sample_id",
            "haplotype",
            "region",
        ]
    ).reset_index(drop=True)


def filter_input_excel_to_target_regions(input_excel: pd.DataFrame) -> pd.DataFrame:
    input_excel = input_excel.copy()

    input_excel = input_excel[
        input_excel["Region"].isin(TARGET_REGIONS)
    ].copy()

    return input_excel.sort_values(
        by=[
            "Sample",
            "Region",
            "Segment",
            "Start coord",
            "End coord",
        ]
    ).reset_index(drop=True)


def add_parsed_excel_sample_columns(input_excel: pd.DataFrame) -> pd.DataFrame:
    input_excel = input_excel.copy()
    sample_text = input_excel["Sample"].astype(str)

    haplotype = sample_text.str.extract(
        r"(hap[12])$",
        expand=False,
    )

    sample_id_trio_compact = sample_text.str.extract(
        r"^(.+?)triohap[12]_",
        expand=False,
    )

    sample_id_compact = sample_text.str.extract(
        r"^(.+?)hap[12]_",
        expand=False,
    )

    sample_id_underscore = sample_text.str.extract(
        r"^([^_]+)_",
        expand=False,
    )

    sample_id = (
        sample_id_trio_compact
        .fillna(sample_id_compact)
        .fillna(sample_id_underscore)
        .astype(str)
        .str.lower()
    )

    input_excel["haplotype"] = haplotype
    input_excel["sample_id"] = sample_id

    return input_excel


def match_regions_to_excel(
    sample_regions: pd.DataFrame,
    input_excel: pd.DataFrame,
) -> pd.DataFrame:
    excel = add_parsed_excel_sample_columns(input_excel)

    matched = sample_regions.merge(
        excel,
        left_on=[
            "sample_id",
            "haplotype",
            "region",
        ],
        right_on=[
            "sample_id",
            "haplotype",
            "Region",
        ],
        how="left",
    )

    return matched.sort_values(
        by=[
            "sample_id",
            "haplotype",
            "region",
            "Start coord",
            "End coord",
        ],
        na_position="last",
    ).reset_index(drop=True)


def add_bam_investigation_coordinates(matched: pd.DataFrame) -> pd.DataFrame:
    matched = matched.copy()

    matched = matched.rename(
        columns={
            "Sample": "excel_sample",
        }
    )

    matched["bam_invest_start_coord"] = (
        matched["region_start"] + matched["Start coord"] - 1
    )

    matched["bam_invest_end_coord"] = (
        matched["region_start"] + matched["End coord"] - 1
    )

    matched["bam_region_string"] = (
        matched["region_id"].astype(str)
        + ":"
        + matched["bam_invest_start_coord"].astype("Int64").astype(str)
        + "-"
        + matched["bam_invest_end_coord"].astype("Int64").astype(str)
    )

    return matched


def add_bam_paths_to_matched_table(
    matched: pd.DataFrame,
    samples_with_bams: pd.DataFrame,
) -> pd.DataFrame:
    bam_paths = samples_with_bams[
        [
            "sample_id",
            "pacbio_bam",
            "nanopore_bam",
            "has_pacbio",
            "has_nanopore",
            "missing",
        ]
    ].copy()

    return matched.merge(
        bam_paths,
        on="sample_id",
        how="left",
    )


def prepare_final_output_table(matched: pd.DataFrame) -> pd.DataFrame:
    final_table = matched.copy()

    final_table = final_table.rename(
        columns={
            "Region": "excel_region",
        }
    )

    final_columns = [
        "sample_id",
        "haplotype",
        "Short name",
        "region",
        "region_id",
        "region_start",
        "region_end",
        "region_name",
        "excel_sample",
        "excel_region",
        "Segment",
        "Start coord",
        "End coord",
        "Strand",
        "bam_invest_start_coord",
        "bam_invest_end_coord",
        "bam_region_string",
        "pacbio_bam",
        "nanopore_bam",
        "has_pacbio",
        "has_nanopore",
        "missing",
    ]

    existing_columns = [
        column for column in final_columns
        if column in final_table.columns
    ]

    return final_table[existing_columns].copy()


def print_missing_bam_summary(bam_presence: pd.DataFrame) -> None:
    samples_without_any_bam = bam_presence[
        ~bam_presence["has_any_bam"]
    ].copy()

    samples_with_partial_bams = bam_presence[
        bam_presence["has_any_bam"]
        & ~bam_presence["has_all_bams"]
    ].copy()

    if not samples_with_partial_bams.empty:
        print()
        print("Samples with only one BAM type available:")
        print(
            samples_with_partial_bams[
                [
                    "sample_id",
                    "has_pacbio",
                    "has_nanopore",
                    "missing",
                    "pacbio_bam",
                    "nanopore_bam",
                ]
            ].to_string(index=False)
        )

    if not samples_without_any_bam.empty:
        print()
        print("Samples removed because no BAM file was found:")
        print(
            samples_without_any_bam[
                [
                    "sample_id",
                    "matching_folders",
                    "all_bams_in_matching_folders",
                ]
            ].to_string(index=False)
        )


def print_unmatched_excel_rows(matched: pd.DataFrame) -> None:
    unmatched = matched[
        matched["excel_sample"].isna()
    ].copy()

    if unmatched.empty:
        return

    print()
    print("Region rows with no Excel match:")
    print(
        unmatched[
            [
                "sample_id",
                "haplotype",
                "region",
                "region_id",
                "region_name",
            ]
        ].drop_duplicates().to_string(index=False)
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare final BAM coordinate table from regions, Excel, and BAM directory. "
            "Samples are kept if at least one BAM is found: PacBio or Nanopore."
        )
    )

    parser.add_argument(
        "--regions",
        required=True,
        help="Regional coordinates TSV/BED-like file without header.",
    )

    parser.add_argument(
        "--input-excel",
        required=True,
        help="Excel file with gene coordinates.",
    )

    parser.add_argument(
        "--bam-dir",
        required=True,
        help="Top-level directory containing sample folders.",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output tab-delimited text file for the final matched table.",
    )

    parser.add_argument(
        "--bam-presence-output",
        default=None,
        help="Optional TSV output showing which BAMs were found for each sample.",
    )

    args = parser.parse_args()

    regions = read_region_coordinates(args.regions)
    input_excel = read_input_excel(args.input_excel)

    unique_regions = get_unique_region_ids(regions)

    bam_presence = build_bam_presence_table(
        unique_regions=unique_regions,
        bam_dir=args.bam_dir,
    )

    if args.bam_presence_output:
        bam_presence.to_csv(
            args.bam_presence_output,
            sep="\t",
            index=False,
        )

    samples_with_bams = bam_presence[
        bam_presence["has_any_bam"]
    ].copy()

    sample_regions = build_sample_regions_table(
        regions=regions,
        samples_with_bams=samples_with_bams,
    )

    input_excel_target_regions = filter_input_excel_to_target_regions(input_excel)

    matched_regions_to_excel = match_regions_to_excel(
        sample_regions=sample_regions,
        input_excel=input_excel_target_regions,
    )

    matched_regions_to_excel = add_bam_investigation_coordinates(
        matched_regions_to_excel
    )

    matched_regions_to_excel = add_bam_paths_to_matched_table(
        matched=matched_regions_to_excel,
        samples_with_bams=samples_with_bams,
    )

    final_table = prepare_final_output_table(matched_regions_to_excel)

    final_table.to_csv(
        args.output,
        sep="\t",
        index=False,
    )

    print(f"Wrote final matched table to: {args.output}")

    if args.bam_presence_output:
        print(f"Wrote BAM presence table to: {args.bam_presence_output}")

    print_missing_bam_summary(bam_presence)
    print_unmatched_excel_rows(matched_regions_to_excel)

    print()
    print(final_table.to_string(index=False))


if __name__ == "__main__":
    main()