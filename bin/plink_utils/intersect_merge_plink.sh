#!/bin/bash


# =============================================================================
# Script Name   : intersect_merge_plink.sh
# Description   : This script updates SNP IDs in multiple PLINK datasets 
#                 (array- and sequencing-based), identifies common SNPs across
#                 datasets, filters based on those common SNPs, and merges 
#                 the resulting PLINK files into a unified dataset.
#
# Input         : 
#   $1 - Comma-separated list of PLINK base paths (without .bed/.bim/.fam)
#   $2 - Output prefix for the merged dataset
# Requirements  : PLINK (v1.9 or v2), awk, bash, optional SLURM environment
# Author        : Florian Bénitière
# Date          : 2025-07-27
# =============================================================================

# --------------------------------------
# Step 0: Define input PLINK files
# --------------------------------------

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "❌ Error: Please provide:"
  echo "  1) comma-separated list of PLINK file prefixes"
  echo "  2) output prefix for merged dataset"
  echo "Usage: $0 \"path1,path2,path3,...\" merged_output_prefix"
  exit 1
fi

plink_files="$1"
merged_prefix="${2:-plink_merged}"

# --------------------------------------
# Step 1: Detect CPU and memory resources
# --------------------------------------
# Create a secure temporary working directory
tmpdir="./plinkmerge_tmp"
mkdir -p "$tmpdir"
echo "🔧 Created temporary working directory: $tmpdir"

# Use SLURM resources if available, otherwise fall back to local system info
cpus="${SLURM_CPUS_ON_NODE:-$(nproc)}"
echo "💻 Running with $cpus cores"

# Get memory safely: prefer SLURM allocated memory
if [[ -n "$SLURM_MEM_PER_NODE" ]]; then
  # Use 90% of SLURM allocated memory (as safety margin)
  mem_MB=$(( SLURM_MEM_PER_NODE * 90 / 100 ))
  echo "Detected SLURM memory: $SLURM_MEM_PER_NODE MB"
else
  # Fallback to checking system memory
  read total_mem used_mem free_mem shared_mem buff_cache available_mem <<< $(free -m | awk '/Mem:/ {print $2, $3, $4, $5, $6, $7}')
  echo "Available memory (MB): $available_mem"

  # Use 90% of available memory
  mem_MB=$(( available_mem * 90 / 100 ))
fi


# --------------------------------------
# Step 2: Ensure PLINK is available
# --------------------------------------

# Check if plink is available
if ! command -v plink &> /dev/null; then
    echo "plink not found — loading module..."
    module load StdEnv/2020
    module load plink/1.9b_6.21-x86_64
else
    echo "plink available."
fi

# --------------------------------------
# Step 3: Prepare SNP name updates and find intersection
# --------------------------------------

# Convert comma-separated string into array by replacing ',' with space
IFS=',' read -r -a plink_array <<< "$plink_files"

# For each PLINK dataset:
# - Generate new SNP IDs using chromosome_position_alleles format
# - Create an update file: <basename>_updated_SNPname.tsv
for base_path in "${plink_array[@]}"; do
  # Expand ~ to full home path
  base_path_expanded=$(eval echo "$base_path")

  # Define input .bim file and output .tsv with new SNP names
  plink_bim="${base_path_expanded}.bim"
  out_tsv="${tmpdir}/${base_path_expanded##*/}_updated_SNPname.tsv"  # filename only + _updated.tsv

  echo "Processing: $plink_bim"

  # Use chr, pos, allele1, allele2 to create a new unique SNP name
  # but only once per unique new SNP ID (to avoid duplicates).
    awk -v OFS='\t' '{
      chr = $1
      sub(/^chr/, "", chr)        # Remove leading "chr" if present
      chr = "chr"chr            # Add "chr" prefix explicitly
      new_id = chr"_"$4"_"$5"_"$6
      if (!seen[new_id]++) print $2, new_id
      }' "$plink_bim" > "$out_tsv"

    updated_files+=("$out_tsv")
done

# --------------------------------------
# Step 4: Compute SNP intersection across all datasets
# --------------------------------------

# Start with SNPs from first update file
cp "${updated_files[0]}" ${tmpdir}/intersect_tmp.tsv

# Iteratively intersect with SNPs from other files using new SNP names (column 2)
for ((i=1; i < ${#updated_files[@]}; i++)); do
  awk 'NR==FNR{a[$2]; next} $2 in a' ${tmpdir}/intersect_tmp.tsv "${updated_files[i]}" > intersect_next.tsv
  mv intersect_next.tsv ${tmpdir}/intersect_tmp.tsv
done

# Extract only SNP ids (2nd column) of final intersection
cut -f2 ${tmpdir}/intersect_tmp.tsv > ${tmpdir}/updated_intersection_SNP_all.txt

# Count total number of SNPs in intersection
total_snps=$(wc -l < ${tmpdir}/updated_intersection_SNP_all.txt)

echo "=== SNP Intersection Summary ==="
echo "Total SNP probes in intersection file: $total_snps"


# --------------------------------------
# Step 5: Apply SNP ID updates and filter by intersection
# --------------------------------------

# For each original dataset:
# - Apply updated SNP names using --update-name
# - Filter on intersection SNPs using --extract
# - Output reduced dataset with common SNPs
for base_path in "${plink_array[@]}"; do
  base_path_expanded=$(eval echo "$base_path")

  output_prefix="${tmpdir}/${base_path_expanded##*/}_plink_reduced"

  echo "Updating SNP names for: $base_path_expanded"
  plink --bfile "$base_path_expanded" \
         --update-name "${tmpdir}/${base_path_expanded##*/}_updated_SNPname.tsv" \
         --make-bed \
         --out "${tmpdir}/${base_path_expanded##*/}_update_name" --memory ${mem_MB} --threads ${cpus}

  echo "Filtering on shared SNPs..."
  plink --bfile "${tmpdir}/${base_path_expanded##*/}_update_name" \
         --extract "${tmpdir}/updated_intersection_SNP_all.txt" \
         --make-bed \
         --out "${output_prefix}" --memory ${mem_MB} --threads ${cpus}

  filtered_outputs+=("$output_prefix")
done


# --------------------------------------
# Step 6: Merge reduced PLINK datasets
# --------------------------------------
echo "Creating merge list file..."

# Add all but the first (reference) file to the merge list
for i in "${filtered_outputs[@]:1}"; do
  echo "${i}.bed ${i}.bim ${i}.fam" >> "${tmpdir}/merge_list.txt"
done


# Merge all datasets using the first one as the reference
echo "Merging filtered PLINK datasets..."

plink --bfile "${filtered_outputs[0]}" \
      --merge-list "${tmpdir}/merge_list.txt" \
      --make-bed \
      --out "$merged_prefix" \
      --memory ${mem_MB} --threads ${cpus}

# Clean up
rm -rf "$tmpdir" "${merged_prefix}.log" "${merged_prefix}.nosex"

echo "✅ Done. Merged dataset prefix: $merged_prefix"
