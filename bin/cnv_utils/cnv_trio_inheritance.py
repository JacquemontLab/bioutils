#!/usr/bin/env python3
 

# =============================================================================
# CNV Inheritance Analysis Pipeline
# Author: Florian Bénitière
# Date: 22/07/2025
# Implementation: Python (Polars)
#
# Description:
# This pipeline analyzes the inheritance of CNVs (Copy Number Variants) 
# in trio datasets. It filters CNVs based on pedigree information, 
# constructs parent–child mappings, performs overlap checks with bedtools, 
# and outputs an annotated child CNV table.
#
# Requirements:
# - Python packages: polars, pyarrow
# - External tools: zcat, awk, bedtools
#
# Input:
# - CNV file (.tsv.gz)
# - Pedigree file (tab-separated: SampleID, FatherID, MotherID)
#
# Output:
# - Annotated child CNV BED and TSV files
# =============================================================================


import os
import subprocess
import multiprocessing
import polars as pl
import argparse

# ------------------- Parse Arguments -------------------

parser = argparse.ArgumentParser(
    description="CNV Inheritance Analysis Pipeline (Polars implementation)"
)
parser.add_argument(
    "-c", "--cnv", required=True,
    help="Path to merged CNV file (.tsv.gz) (SampleID, Chr, Start, End)"
)
parser.add_argument(
    "-p", "--pedigree", required=True,
    help="Path to pedigree file (.tsv) (SampleID, FatherID, MotherID)"
)
parser.add_argument(
    "-o", "--output", default="annotated_child_cnv.tsv",
    help="Path to the output child CNVs"
)
parser.add_argument(
    "-t", "--type_col", default=None,
    help="Column name of CNV type (DEL, DUP, MIX etc...) for per-type annotation; leave empty to skip"
)
parser.add_argument(
    "-f", "--overlap", default="0.5",
    help="Comma-separated list of reciprocal overlap fractions (e.g., 0.5,0.1)"
)

args = parser.parse_args()

FILE_CNV = args.cnv
PEDIGREE_FILE = args.pedigree
OUTPUT = args.output
TYPE_COL = args.type_col
OVERLAPS = [float(x) for x in args.overlap.split(",")]

# ------------------- Filter CNV Calls -------------------

bash_cmd = f"""
if [[ "{FILE_CNV}" == *.gz ]]; then
  zcat {FILE_CNV} | awk 'FNR==NR {{ids[$1]; ids[$2]; ids[$3]; next}} FNR==1 || $1 in ids' {PEDIGREE_FILE} - > filtered.tsv
else
  awk 'FNR==NR {{ids[$1]; ids[$2]; ids[$3]; next}} FNR==1 || $1 in ids' {PEDIGREE_FILE} {FILE_CNV} > filtered.tsv
fi
"""

print("Running command:\n", bash_cmd)
subprocess.run(bash_cmd, shell=True, check=True)

# ------------------- Load Filtered Data -------------------

FILE = f"filtered.tsv"
data_used = pl.read_csv(FILE, separator="\t", infer_schema_length=10000)

# Ensure SampleID is string
data_used = data_used.with_columns(pl.col("SampleID").cast(pl.Utf8))

# Load pedigree file
family_info = pl.read_csv(PEDIGREE_FILE, separator="\t").with_columns([
    pl.col("SampleID").cast(pl.Utf8),
    pl.col("FatherID").cast(pl.Utf8),
    pl.col("MotherID").cast(pl.Utf8),
])

# Filter trios: child, father, mother all in CNV data
sample_ids = data_used["SampleID"].unique()
family_info = family_info.filter(
    (family_info["SampleID"].is_in(sample_ids.implode())) &
    (family_info["FatherID"].is_in(sample_ids.implode())) &
    (family_info["MotherID"].is_in(sample_ids.implode()))
)

# ------------------- Construct Parent → Children Lookup -------------------

parent_ids = pl.concat([
    family_info["FatherID"],
    family_info["MotherID"]
]).unique()

parent_of = {
    pid: ",".join(
        family_info.filter(
            (family_info["FatherID"] == pid) | (family_info["MotherID"] == pid)
        )["SampleID"].to_list()
    )
    for pid in parent_ids
}

# ------------------- Annotate CNVs with Children -------------------

n_cores = multiprocessing.cpu_count() - 1
chunk_size = (len(data_used) + n_cores - 1) // n_cores
chunks = [data_used[i:i+chunk_size] for i in range(0, len(data_used), chunk_size)]

def annotate_chunk(chunk: pl.DataFrame):
    return chunk.with_columns(
        pl.col("SampleID").map_elements(lambda sid: parent_of.get(sid, ""), return_dtype=pl.Utf8).alias("child")
    )

with multiprocessing.Pool(n_cores) as pool:
    result_chunks = pool.map(annotate_chunk, chunks)

data_used = pl.concat(result_chunks)

# ------------------- Reshape Data -------------------

cnv_ungrouped = data_used.with_columns(
    pl.col("child").str.split(",")
).explode("child")

child_cnv = cnv_ungrouped.filter(
    cnv_ungrouped["SampleID"].is_in(family_info["SampleID"].implode())
)

parents_cnv = cnv_ungrouped.filter(
    cnv_ungrouped["SampleID"].is_in(pl.concat([family_info["FatherID"], family_info["MotherID"]]).implode())
)

if TYPE_COL is None:
    # no type column, just SampleID and Chr
    child_cnv = child_cnv.with_columns(
        (pl.col("SampleID") + "_" + pl.col("Chr")).alias("Chr")
    )
    parents_cnv = parents_cnv.with_columns(
        (pl.col("child") + "_" + pl.col("Chr")).alias("Chr")
    )
else:
    # include the type column
    child_cnv = child_cnv.with_columns(
        (pl.col("SampleID") + "_" + pl.col("Chr") + "_" + pl.col(TYPE_COL)).alias("Chr")
    )
    parents_cnv = parents_cnv.with_columns(
        (pl.col("child") + "_" + pl.col("Chr") + "_" + pl.col(TYPE_COL)).alias("Chr")
    )

parents_cnv = parents_cnv.select(["Chr", "Start", "End"])
child_cnv = child_cnv.drop(["SampleID", "child"])

# ------------------- Export BED Files -------------------

child_bed = f"child_cnv.bed"
parent_bed = f"parents_cnv.bed"

child_cnv.write_csv(child_bed, separator="\t", include_header=False)
parents_cnv.write_csv(parent_bed, separator="\t", include_header=False)

# ------------------- Run bedtools intersect -------------------

for ovlap in OVERLAPS:
    bash_cmd = f"bedtools intersect -a {child_bed} -b {parent_bed} -f {ovlap} -r -wa > intersect_ovlap{ovlap}.bed"
    print("Running command:\n", bash_cmd)
    subprocess.run(bash_cmd, shell=True, check=True)

# ------------------- Save Final Annotated Output -------------------

# Re-import child BED with proper column names
vector_names = pl.read_csv(FILE, separator="\t", infer_schema_length=10000).columns[1:]  # skip first col
child_cnv = pl.read_csv(child_bed, separator="\t", has_header=False, infer_schema_length=10000)
child_cnv = child_cnv.rename({f"column_{i+1}": name for i, name in enumerate(vector_names)})

child_cnv = child_cnv.with_columns(
    (pl.col("Chr") + "_" + pl.col("Start").cast(pl.Utf8) + "_" + pl.col("End").cast(pl.Utf8)).alias("ID")
)

for ovlap in OVERLAPS:
    result = pl.read_csv(f"intersect_ovlap{ovlap}.bed", separator="\t", infer_schema_length=10000, has_header=False)
    result = result.rename({f"column_{i+1}": name for i, name in enumerate(vector_names)})
    
    result = result.with_columns(
        (pl.col("Chr") + "_" + pl.col("Start").cast(pl.Utf8) + "_" + pl.col("End").cast(pl.Utf8)).alias("ID")
    )
        
    child_cnv = child_cnv.with_columns(
        pl.col("ID").is_in(result["ID"].implode()).alias(f"observed_in_parent_ovlap{ovlap}")
    )

child_cnv = child_cnv.with_columns(
    pl.col("Chr").str.split("_").list.get(0).alias("SampleID"),
    pl.col("Chr").str.split("_").list.get(1).alias("Chr")
).drop("ID")
    

# Reorder so that SampleID is first
child_cnv = child_cnv.select(
    ["SampleID"] + [col for col in child_cnv.columns if col != "SampleID"]
)

child_cnv.write_csv(OUTPUT, separator="\t")

print("Pipeline completed successfully.")


print(f"💾 Results written to {OUTPUT}")


# --------------------- Cleanup Temporary Files ---------------------
files_to_remove = [
    parent_bed,
    child_bed,
    FILE
] + [f"intersect_ovlap{ov}.bed" for ov in OVERLAPS]

for f in files_to_remove:
    try:
        os.remove(f)
        print(f"🗑️ Removed {f}")
    except FileNotFoundError:
        print(f"⚠️ Skipped {f} (not found)")
