#!/bin/bash


# =============================================================================
# Script Name   : remove_dup_position_plink.sh
# Description   : Remove SNPs sharing the same chromosome+position,
#                 keeping only the SNP with the best call rate.
# Usage         : ./remove_dup_position_plink.sh <plink_prefix> [output_prefix]
# Requirements  : plink (v1.9 or v2), awk, sort, join, bash
# Author        : Florian Bénitière
# Date          : 2025-07-28
# =============================================================================


if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 <plink_prefix> [output_prefix]"
  exit 1
fi

plink_prefix="$1"
output_prefix="${2:-plink_no_dup_pos}"

# --------------------------------------
# Step 1: Detect CPU and memory resources
# --------------------------------------

cpus="${SLURM_CPUS_ON_NODE:-$(nproc)}"
echo "💻 Running with $cpus cores"

if [[ -n "${SLURM_MEM_PER_NODE:-}" ]]; then
  mem_MB=$(( SLURM_MEM_PER_NODE * 90 / 100 ))
  echo "Detected SLURM memory: $SLURM_MEM_PER_NODE MB (using 90%)"
else
  read -r total_mem used_mem free_mem shared_mem buff_cache available_mem < <(free -m | awk '/Mem:/ {print $2, $3, $4, $5, $6, $7}')
  mem_MB=$(( available_mem * 90 / 100 ))
  echo "Available memory (MB): $available_mem (using 90%)"
fi

# --------------------------------------
# Step 2: Ensure PLINK is available
# --------------------------------------

if ! command -v plink &> /dev/null; then
    echo "plink not found — loading module..."
    module load StdEnv/2020
    module load plink/1.9b_6.21-x86_64
else
    echo "plink available."
fi

# --------------------------------------
# Step 3: Prepare working directory
# --------------------------------------

tmpdir="./remove_dup_pos"
mkdir -p "$tmpdir"
echo "🔧 Using temporary directory: $tmpdir"

# --------------------------------------
# Step 4: Generate SNP position names
# --------------------------------------

awk -v OFS='\t' '{
        chr = $1
        sub(/^chr/, "", chr)
        chr = "chr" chr
        new_id = chr"_"$4
        print $2, new_id
    }' "${plink_prefix}.bim" > "${tmpdir}/name_position.tsv"

if ! awk '{print $2}' "${tmpdir}/name_position.tsv" | sort | uniq -d | grep -q .; then
  echo "No duplicate positions found. Skipping duplicate removal."
  plink --bfile "$plink_prefix" --make-bed --out "$output_prefix" --memory ${mem_MB} --threads ${cpus}
  rm -rf "$tmpdir" "$output_prefix.log" "$output_prefix.nosex"
  exit 0
fi

# --------------------------------------
# Step 5: Calculate per-SNP call rates
# --------------------------------------

echo "Calculating per-SNP call rates with PLINK..."
plink --bfile "$plink_prefix" --missing --out "${tmpdir}/stats_missing" --memory ${mem_MB} --threads ${cpus}

# --------------------------------------
# Step 6: Find duplicates and select best call rate SNPs
# --------------------------------------

echo "Processing call rates and identifying duplicate positions..."

awk -v OFS='\t' 'NR>1 {print $2, 1 - $5}' "${tmpdir}/stats_missing.lmiss" | sort -k1,1 > "${tmpdir}/snp_callrate.tsv"
sort -k1,1 "${tmpdir}/name_position.tsv" > "${tmpdir}/name_position_sorted.tsv"

# Join on SNP ID (column 1)
join -t $'\t' "${tmpdir}/snp_callrate.tsv" "${tmpdir}/name_position_sorted.tsv" > "${tmpdir}/merged_callrate_position.tsv"

# Sort by call rate descending (column 2)
sort -k2,2nr "${tmpdir}/merged_callrate_position.tsv" > "${tmpdir}/merged_sorted_desc.tsv"

# Print duplicates (all but first occurrence per position)
awk '
{
  pos = $3
  count[pos]++
  if(count[pos] > 1) print $1
}
' "${tmpdir}/merged_sorted_desc.tsv" > "${tmpdir}/duplicated_snps.txt"

dup_count=$(wc -l < "${tmpdir}/duplicated_snps.txt")
echo "Found $dup_count duplicate SNPs by position to exclude."

# --------------------------------------
# Step 7: Create filtered PLINK dataset
# --------------------------------------

echo "Removing duplicate SNPs and creating filtered dataset: $output_prefix"
plink --bfile "$plink_prefix" --exclude "${tmpdir}/duplicated_snps.txt" --make-bed --out "$output_prefix" --memory ${mem_MB} --threads ${cpus}

# --------------------------------------
# Step 8: Cleanup
# --------------------------------------

rm -rf "$tmpdir" "$output_prefix.log" "$output_prefix.nosex"

echo "✅ Done. Filtered dataset prefix: $output_prefix"
