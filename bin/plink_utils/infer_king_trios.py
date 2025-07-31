#!/usr/bin/env python3

"""
Script Name: infer_king_trios.py
Author: Florian Bénitière
Date: 2025-07-29

Description:
------------
This script runs KING and PLINK to detect trios in genetic data. 
 Process KING and PLINK outputs with Polars to:
   - Detect children with exactly two parents.
   - Confirm parents are of opposite sexes.
   - Exclude cases where parents are related.
   - Assign Family IDs.

Usage:
------
    python infer_king_trios.py <plink_prefix> [output_tsv_file]
    
Output:
-------
- TSV with FamilyID, FatherID, and MotherID for valid trios.

"""

import subprocess
import sys
import os
import polars as pl


# --------------------- Argument Parsing ---------------------
if len(sys.argv) < 2:
    print("Usage: python script.py <plink_prefix> [output_tsv_file]")
    sys.exit(1)

plink_prefix = sys.argv[1]
OUTPUT_TSV = sys.argv[2] if len(sys.argv) > 2 else "king_trios.tsv"



# --------------------- Check KING ---------------------
if subprocess.call("command -v king >/dev/null 2>&1", shell=True) != 0:
    print("❌ Error: 'king' is not available in your PATH. Please install or load the module.")
    sys.exit(1)

# --------------------- Check PLINK ---------------------

def check_and_load_plink():
    """Check if plink is available, otherwise load via module."""
    result = subprocess.call("command -v plink >/dev/null 2>&1", shell=True)
    if result != 0:
        print("plink not found — loading module...")
        return "module load StdEnv/2020 && module load plink/1.9b_6.21-x86_64 &&"
    else:
        print("plink available.")
        return ""

load_module = check_and_load_plink()

# --------------------- Detect Number of CPUs ---------------------
cpus = os.environ.get("SLURM_CPUS_ON_NODE")
if cpus is None:
    cpus = str(os.cpu_count())
print(f"💻 Running KING with {cpus} cores...")

# --------------------- Run PLINK check-sex ---------------------
print("🔎 Running PLINK check-sex...")
# Make sure you have 'plink_prefix' defined somewhere before this
plink_command = (
    f"{load_module} plink --bfile {plink_prefix} --check-sex --threads {cpus} --out temp_plink_sexcheck"
)

subprocess.run(plink_command, shell=True, check=True)

# Convert space-separated to tab-separated
cmd = r"""awk 'BEGIN {OFS="\t"} {$1=$1; print}' temp_plink_sexcheck.sexcheck > temp_plink_sexcheck.sexcheck.tsv"""
subprocess.run(cmd, shell=True, check=True)

print("✅ PLINK check-sex finished successfully")

# --------------------- Run KING ---------------------
# First awk: copy col2 into col1 of .fam
subprocess.run(
    f"awk '{{ $1 = $2; print }}' OFS='\\t' {plink_prefix}.fam > tmp.fam",
    shell=True, check=True, executable="/bin/bash"
)

# Second awk: remove character chr of .bim
subprocess.run(
    f"awk '{{ sub(/^chr/, \"\", $1); print }}' OFS='\\t' {plink_prefix}.bim > tmp.bim",
    shell=True, check=True, executable="/bin/bash"
)

subprocess.run([
    "king",
    "-b", f"{plink_prefix}.bed",
    "--fam", "tmp.fam",
    "--bim", "tmp.bim",
    "--related",
    "--degree", "2",
    "--cpus", cpus,
    "--prefix", "temp_king_out"
], check=True)

print("✅ KING finished successfully")

# --------------------- Load KING and PLINK outputs ---------------------
dfs = []
for f in ["temp_king_out.kin", "temp_king_out.kin0"]:
    if os.path.exists(f):
        dfs.append(pl.read_csv(f, separator="\t", infer_schema_length=10000))

if len(dfs) == 0:
    raise FileNotFoundError("Neither temp_king_out.kin0 nor temp_king_out.kin found.")
elif len(dfs) == 1:
    family_king = dfs[0]
else:
    family_king = pl.concat(dfs, how="diagonal")

family_king = family_king.filter(pl.col("InfType") != "UN")

# PLINK sex check results
plink_sexcheck = (
    pl.read_csv(
        "temp_plink_sexcheck.sexcheck.tsv",
        separator="\t"
    )
    .select(["IID","SNPSEX"])   # IID and SNPSEX columns
)

# --------------------- Filter Parent-Offspring Pairs ---------------------
# Keep only parent-offspring relations
parents_relation = family_king.filter(pl.col("InfType") == "PO").select(["ID1", "ID2"])


# Add inverse relationships (child->parent)
inverse_relation = parents_relation.select(
    [pl.col("ID2").alias("ID1"), pl.col("ID1").alias("ID2")]
)

# Combine both directions
dt = pl.concat([parents_relation, inverse_relation])

# Add SNPSEX information for parent
dt = dt.join(plink_sexcheck.select(["IID", "SNPSEX"]),
             left_on="ID2", right_on="IID").rename({"SNPSEX": "sex_ID2"})

# --------------------- Identify Valid Trios ---------------------
# Group children by parents
grouped = dt.group_by("ID1").agg(pl.col("ID2").unique())

# Keep only children with exactly 2 unique parents
grouped = grouped.filter(pl.col("ID2").list.len() == 2)

# Extract parents
grouped = grouped.with_columns([
    pl.col("ID2").list.get(0).alias("parent_1"),
    pl.col("ID2").list.get(1).alias("parent_2")
])


# Keep only different-sex parents (1 + 2 = 3)
grouped = grouped.join(plink_sexcheck.rename({"IID": "parent_1", "SNPSEX": "sex_1"}), on="parent_1")
grouped = grouped.join(plink_sexcheck.rename({"IID": "parent_2", "SNPSEX": "sex_2"}), on="parent_2")


grouped = grouped.filter((pl.col("sex_1") + pl.col("sex_2")) == 3)


# Exclude related parents
existing_pairs = pl.concat([
    family_king.select(["ID1", "ID2"]),
    family_king.select([pl.col("ID2").alias("ID1"), pl.col("ID1").alias("ID2")])
]).unique()

check_relation = grouped.join(
    existing_pairs.rename({"ID1": "parent_1", "ID2": "parent_2"}),
    on=["parent_1", "parent_2"],
    how="anti"
)

# Get children with valid parents
iid_with_parents = check_relation.select("ID1").rename({"ID1": "IID"}).unique()

# --------------------- Assign Father and Mother ---------------------
# Fathers: sex = 1
paternal = dt.filter((pl.col("sex_ID2") == 1) & (pl.col("ID1").is_in(iid_with_parents["IID"].implode())))
# Mothers: sex = 2
maternal = dt.filter((pl.col("sex_ID2") == 2) & (pl.col("ID1").is_in(iid_with_parents["IID"].implode())))

df = iid_with_parents.join(paternal.select(["ID1","ID2"]).rename({"ID1":"IID","ID2":"PID"}), on="IID", how="left")
df = df.join(maternal.select(["ID1","ID2"]).rename({"ID1":"IID","ID2":"MID"}), on="IID", how="left")

# Rename final columns
df = df.rename({
    "IID": "SampleID",
    "PID": "FatherID",
    "MID": "MotherID"
})

print("✅ Trio assignment complete:")
print(df)

# --------------------- Assign Family IDs ---------------------
# Ensure IDs are strings
df = df.with_columns([
    pl.col("FatherID").cast(pl.Utf8),
    pl.col("MotherID").cast(pl.Utf8),
])

# Create a parent pair key
df = df.with_columns(
    (pl.col("FatherID").fill_null("NA") + "_" + 
     pl.col("MotherID").fill_null("NA")).alias("parent_pair")
)

# Assign FamilyID using dense ranking
df = df.with_columns(
    ("Family" + pl.col("parent_pair").rank(method="dense").cast(pl.Utf8)).alias("FamilyID")
)

# Drop helper column
df = df.drop("parent_pair")

# --------------------- Save Results ---------------------
df.write_csv(OUTPUT_TSV, separator="\t")
print(f"💾 Results written to {OUTPUT_TSV}")


# --------------------- Cleanup Temporary Files ---------------------
files_to_remove = [
    "temp_king_out.kin0",
    "temp_king_outX.kin0",
    "temp_king_out.log",
    "temp_king_outallsegs.txt",
    "temp_king_out.kin",
    "temp_king_outX.kin",
    "temp_plink_sexcheck.sexcheck",
    "temp_plink_sexcheck.nosex",
    "temp_plink_sexcheck.log",
    "temp_plink_sexcheck.hh",
    "temp_plink_sexcheck.sexcheck.tsv",
    "tmp.bim",
    "tmp.fam"
]

for f in files_to_remove:
    try:
        os.remove(f)
        print(f"🗑️ Removed {f}")
    except FileNotFoundError:
        print(f"⚠️ Skipped {f} (not found)")
