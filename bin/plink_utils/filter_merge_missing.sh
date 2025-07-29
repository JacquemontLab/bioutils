#!/bin/bash



##########################################################################################
#
# Script: filter_merge_missing.sh
# Author: Mame Seynabou Diop
# Date: 29/07/2025
#
# Description:
#   Applies per-chromosome SNP quality control filters and merges filtered datasets:
#     - Filters SNP and individual missingness, Hardy-Weinberg equilibrium, and minor allele frequency
#       using plink2 options provided by the user.
#     - Merges filtered per-chromosome binary PLINK files into one dataset.
#
# Notes:
#   Filtering parameters (e.g. --geno, --mind, --hwe, --maf) are optional and can be specified
#   as additional plink2 options when running the script. Defaults can be set within the script.
#
# Usage:
#   ./filter_merge_missing.sh <input_bfile_dir> <outputPath> [additional plink2 options]
#
# Example:
#   ./filter_merge_missing.sh data/genos output/ --geno 0.02 --mind 0.02 --hwe 1e-5 --maf 0.005
#
##########################################################################################

# ----------------------------- Argument checks ------------------------------------

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <input_bfile_prefix> <output_file> [additional plink2 options]"
    echo ""
    echo "Example:"
    echo "  $0 data/genos/geno output/ --geno 0.05 --mind 0.05 --hwe 1e-6 --maf 0.001"
    exit 1
fi

inputDir="$1"
outputPath="$2"
shift 2  # Remove inputPrefix and outputPath from $@

# Check input dir
if [[ ! -d "$inputDir" ]]; then
    echo "Error: input directory '$inputDir' does not exist."
    exit 1
fi


# ----------------------------- Find bfiles ----------------------------------------
# Get sorted list of all .bed files
bed_files=($(find "$inputDir" -maxdepth 1 -name "*.bed" | sort))


# Use the first .bed file as the base
base_bed="${bed_files[0]}"
base_prefix="${base_bed%.bed}"

# Load plink2 if needed
if ! command -v plink2 >/dev/null 2>&1; then
    echo "'plink2' not found — loading modules..."
    module load StdEnv/2020
    module load plink/2.00a3.6

    # Check again
    if ! command -v plink2 >/dev/null 2>&1; then
        echo "Error: plink2 is still not available after loading modules." >&2
        exit 1
    fi
else
    echo "plink2 is already available."
fi


# Get number of CPUs
cpus="${SLURM_CPUS_ON_NODE:-$(nproc)}"
echo "💻 Running with $cpus cores"

# Get memory safely: prefer SLURM allocated memory
if [[ -n "$SLURM_MEM_PER_NODE" ]]; then
  # Use 90% of SLURM allocated memory (as safety margin)
  plink_mem=$(( SLURM_MEM_PER_NODE * 90 / 100 ))
  echo "Detected SLURM memory: $SLURM_MEM_PER_NODE MB"
else
  # Fallback to checking system memory
  read total_mem used_mem free_mem shared_mem buff_cache available_mem <<< $(free -m | awk '/Mem:/ {print $2, $3, $4, $5, $6, $7}')
  echo "Available memory (MB): $available_mem"

  # Use 90% of available memory
  plink_mem=$(( available_mem * 90 / 100 ))
fi

echo "Setting PLINK memory to: $plink_mem MB"
echo "Setting PLINK threads to: $cpus"


# Define directories
filteredDir=filtered_chr

mkdir -p $filteredDir

# ----------------------------------------------
# Loop over all .bed files in the input directory
# ----------------------------------------------
for bed_file in "$inputDir"/*.bed; do
  base_name=$(basename "$bed_file" .bed)
  prefix="$inputDir/$base_name"

  echo "Filtering $base_name..."

  plink2 \
    --bfile "$prefix" \
    --threads "$cpus" \
    --memory "$plink_mem" \
    "$@" \
    --make-bed \
    --out "$filteredDir/${base_name}.filtered"
done


# Step 2: Create list.txt of all chromosomes except first vcf (for merge)
list_file="list.txt"
> "$list_file"


skip_first=true
for bed_file in "$filteredDir"/*.bed; do
  base_name=$(basename "$bed_file" .bed)
  prefix="$filteredDir/$base_name"

  if $skip_first; then
    first_vcf="$prefix"
    skip_first=false
    continue
  fi

  echo -e "${prefix}.bed\t${prefix}.bim\t${prefix}.fam" >> "$list_file"
done

# Step 3: Merge chromosomes (starting from first vcf)
plink2 \
  --bfile "$first_vcf" \
  --threads "$cpus" \
  --memory "$plink_mem" \
  --pmerge-list list.txt \
  --make-bed \
  --out "$outputPath"


rm -rf list.txt "$filteredDir"
