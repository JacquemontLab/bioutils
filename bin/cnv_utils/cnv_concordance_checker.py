#!/usr/bin/env python3
"""
Compare two CNV datasets.

This script takes two CNV call tables (TSV format) and compares them on two levels:

1. Reciprocal CNV overlap (≥ 50% in both directions):
   - A CNV in dataset A is considered overlapping with one in dataset B if:
       overlap_length >= 50% of length_A AND
       overlap_length >= 50% of length_B
   - Requires matching SampleID, chromosome, and CNV type.

2. Gene-level match:
   - If both datasets contain a Gene_ID column, flags CNVs where
     SampleID, Type, and Gene_ID are all identical in the other dataset.

Outputs:
--------
Two TSV files (one for each input), containing:
    - CNV_based_overlap : 1 if reciprocal overlap found, else 0
    - gene_based_overlap    : 1 if gene match found, NA if Gene_ID not available

Typical use case:
-----------------
Useful for comparing CNV calls from different platforms (e.g. array vs WGS)
or different calling methods to assess concordance.

Each input file must:
    - Be a tab-separated file (.tsv)
    - Contain these mandatory columns:
        SampleID : string, sample identifier (must match between datasets)
        Chr      : string/int, chromosome name (same naming convention in both files)
        Start    : integer, start coordinate (consistent 0-based or 1-based)
        End      : integer, end coordinate (must be > Start)
        Type     : string, CNV type (e.g., DEL, DUP)
    - (Optional) Gene_ID : string, gene identifier (required for gene match comparison)

Example:
--------
python compare_cnv_calls.py \
    --file_a array_calls.tsv \
    --file_b wgs_calls.tsv \
    --out_a array_with_flags.tsv \
    --out_b wgs_with_flags.tsv
"""


import argparse
import pandas as pd
import duckdb
from pathlib import Path


def load_data(file_path: Path) -> pd.DataFrame:
    """Load TSV file into a DataFrame."""
    print(f"[INFO] Loading data from {file_path} ...")
    return pd.read_csv(file_path, sep="\t", low_memory=False)


def add_row_ids(df: pd.DataFrame) -> pd.DataFrame:
    """Add a unique row identifier column."""
    df["rowid"] = range(len(df))
    df
    return df


def flag_overlaps(con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    """Find reciprocal overlaps >= 50% between df_a and df_b."""
    query = """
    WITH overlaps_df AS (
      SELECT 
        df_a.rowid AS id_a,
        df_b.rowid AS id_b,
        df_a.SampleID AS SampleID_a,
        df_b.SampleID AS SampleID_b,
        df_a.Chr,
        GREATEST(df_a.Start, df_b.Start) AS overlap_start,
        LEAST(df_a.End, df_b.End) AS overlap_end,
        (LEAST(df_a.End, df_b.End) - GREATEST(df_a.Start, df_b.Start)) AS overlap_len,
        (df_a.End - df_a.Start) AS len_a,
        (df_b.End - df_b.Start) AS len_b
      FROM df_a
      JOIN df_b
        ON df_a.Chr = df_b.Chr
        AND df_a.SampleID = df_b.SampleID
        AND df_a.Type = df_b.Type
        AND df_a.End > df_b.Start
        AND df_a.Start < df_b.End
    )
    SELECT
      id_a,
      id_b,
      CAST(
        (overlap_len >= 0.5 * len_a AND overlap_len >= 0.5 * len_b) AS INTEGER
      ) AS reciprocal_overlap_flag
    FROM overlaps_df
    WHERE overlap_len > 0
    """
    return con.execute(query).df()


def flag_gene_matches(con: duckdb.DuckDBPyConnection, table_a: str, table_b: str, id_col: str) -> pd.DataFrame:
    """Check if (SampleID, Gene_ID) exists in the other table."""
    query = f"""
    SELECT {table_a}.rowid AS {id_col},
           CASE WHEN {table_b}.rowid IS NOT NULL THEN 1 ELSE 0 END AS gene_based_overlap
    FROM {table_a}
    LEFT JOIN {table_b}
      ON {table_a}.SampleID = {table_b}.SampleID
      AND {table_a}.Type = {table_b}.Type
      AND {table_a}.Gene_ID = {table_b}.Gene_ID
    """
    return con.execute(query).df()

def print_summary(df_a: pd.DataFrame, df_b: pd.DataFrame, name_a: str, name_b: str):
    """Print summary statistics."""
    print("\n=== SUMMARY ===")
    for name, df in [(name_a, df_a), (name_b, df_b)]:
        total = len(df)
        ovlp = df["CNV_based_overlap"].sum()
        gene = df["gene_based_overlap"].sum() if df["gene_based_overlap"].notna().any() else 0
        both = ((df["CNV_based_overlap"] == 1) & (df["gene_based_overlap"] == 1)).sum() if df["gene_based_overlap"].notna().any() else 0

        print(f"{name}: {total} CNVs")
        print(f"  CNV-overlap: {ovlp} ({ovlp / total * 100:.2f}%)")
        if df["gene_based_overlap"].notna().any():
            print(f"  Gene-overlap: {gene} ({gene / total * 100:.2f}%)")
            print(f"  Both overlap types: {both} ({both / total * 100:.2f}%)")
        print()


def main(file_a: Path, file_b: Path, out_a: Path, out_b: Path):
    # Load data
    df_a = add_row_ids(load_data(file_a))
    df_b = add_row_ids(load_data(file_b))

    # Connect to DuckDB
    con = duckdb.connect()
    con.register("df_a", df_a)
    con.register("df_b", df_b)

    # Flag overlaps
    print("[INFO] Computing reciprocal overlaps >= 50% ...")
    overlaps_df = flag_overlaps(con)
    df_a["CNV_based_overlap"] = 0
    df_b["CNV_based_overlap"] = 0
    df_a.loc[overlaps_df.loc[overlaps_df["reciprocal_overlap_flag"] == 1, "id_a"].unique(), "CNV_based_overlap"] = 1
    df_b.loc[overlaps_df.loc[overlaps_df["reciprocal_overlap_flag"] == 1, "id_b"].unique(), "CNV_based_overlap"] = 1

    # Flag gene matches if Gene_ID column exists in both
    if "Gene_ID" in df_a.columns and "Gene_ID" in df_b.columns:
        print("[INFO] Checking gene overlaps ...")
        gene_flag_a = flag_gene_matches(con, "df_a", "df_b", "id_a")
        df_a["gene_based_overlap"] = 0
        df_a.loc[gene_flag_a.loc[gene_flag_a["gene_based_overlap"] == 1, "id_a"], "gene_based_overlap"] = 1

        gene_flag_b = flag_gene_matches(con, "df_b", "df_a", "id_b")
        df_b["gene_based_overlap"] = 0
        df_b.loc[gene_flag_b.loc[gene_flag_b["gene_based_overlap"] == 1, "id_b"], "gene_based_overlap"] = 1
    else:
        print("[INFO] Skipping gene overlap check (Gene_ID column missing in one or both files)")
        df_a["gene_based_overlap"] = pd.NA
        df_b["gene_based_overlap"] = pd.NA

    # Summary
    print_summary(df_a, df_b, "File A", "File B")

    # Drop internal rowid before saving
    df_a = df_a.drop(columns=["rowid"])
    df_b = df_b.drop(columns=["rowid"])

    # Save results
    print(f"[INFO] Saving results to {out_a} and {out_b} ...")
    df_a.to_csv(out_a, sep="\t", index=False)
    df_b.to_csv(out_b, sep="\t", index=False)
    print("[INFO] Done.")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare two CNV call files.")
    parser.add_argument("--file_a", type=Path, required=True, help="Path to first TSV file")
    parser.add_argument("--file_b", type=Path, required=True, help="Path to second TSV file")
    parser.add_argument("--out_a", type=Path, default="comparison_a.tsv", help="Output TSV for file A results")
    parser.add_argument("--out_b", type=Path, default="comparison_b.tsv", help="Output TSV for file B results")

    args = parser.parse_args()

    main(args.file_a, args.file_b, args.out_a, args.out_b)