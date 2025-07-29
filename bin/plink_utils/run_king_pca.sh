#!/bin/bash


# ======================================================
# KING PCA Pipeline
# Author: Florian Bénitière
# Date: 2025-07-28
#
# Description:
#   This script performs ancestry projection PCA using KING 
#   by merging a target PLINK dataset with the 1000 Genomes
#   reference panel (GRCh37 or optionally lifted to GRCh38).
#   It ensures unique SNP IDs, harmonizes genome versions,
#   and produces PCA plots and output files for ancestry 
#   inference.
#
# Usage:
#   ./run_king_pca.sh <plink_prefix> [output_file] [genome_version]
#
# Arguments:
#   <plink_prefix>   Prefix of the input PLINK dataset (.bed/.bim/.fam)
#   [output_file]    Prefix for PCA output (default: king_ancestry_PC)
#   [genome_version] Genome build of the dataset, GRCh37 or GRCh38 
#                    (default: GRCh38)
#
# Example:
#   ./run_king_pca.sh my_dataset my_pca GRCh38
#
# Requirements:
#   - PLINK v1.9+
#   - KING v2.3+
#   - wget, awk, coreutils
#   - liftOver (if converting from GRCh37 to GRCh38)
#   - SLURM environment (optional, for automatic resource detection)
#
# Notes:
#   - The script automatically downloads the 1000 Genomes reference 
#     (GRCh37) from the KING website.
#   - If GRCh38 is requested, it performs liftOver conversion using
#     UCSC chain files.
#   - Intermediate files are handled in a temporary directory and
#     cleaned up at the end.
#   - The pipeline auto-detects CPU cores and available memory, and
#     adjusts PLINK resource usage accordingly.
#
# Output:
#   - KING PCA projection files (prefix: king_prefix)
#   - PCA plots (in .rplot and PDF formats if supported)
#
# ======================================================



set -e  # Exit on any error


# --------------------- Argument Parsing ---------------------
if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <plink_prefix> [output_file] [genome_version]"
  exit 1
fi

# ---------------------------
# Parse arguments
# ---------------------------
plink_prefix=${1:? "Error: You must provide the PLINK dataset prefix"}
output_file=${2:-"king_ancestry_PC.tsv"}
genome_version=${3:-"GRCh38"}

echo "➡ Starting pipeline with input: ${plink_prefix}, genome version: ${genome_version}"

# ---------------------------
# Check PLINK availability
# ---------------------------
if ! command -v plink &> /dev/null; then
    echo "plink not found — loading module..."
    module load StdEnv/2020
    module load plink/1.9b_6.21-x86_64 || { echo "❌ Failed to load PLINK module"; exit 1; }
else
    echo "✅ plink available."
fi

# --------------------- Check if KING is available ---------------------
if ! command -v king &> /dev/null; then
  echo "❌ Error: 'king' is not available in your PATH. Please install KING or load the module."
  exit 1
fi

# ---------------------------
# Resource detection
# ---------------------------
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

echo "Setting memory to: $mem_MB MB"
echo "Setting threads to: $cpus"


# --------------------------------------
# Prepare working directory
# --------------------------------------
tmpdir="./king_pca/"
mkdir -p "$tmpdir"
echo "🔧 Using temporary directory: $tmpdir"


if [ -f /lustre06/project/6008022/All_user_common_folder/RAW_DATA/Genetic/Reference_Data/king_ref/KGref_GRCh37_final.bed ]; then
    echo "✅ KGref files exist."
    cp /lustre06/project/6008022/All_user_common_folder/RAW_DATA/Genetic/Reference_Data/king_ref/KGref_GRCh37_final* $tmpdir
    cp /lustre06/project/6008022/All_user_common_folder/RAW_DATA/Genetic/Reference_Data/king_ref/KGref_GRCh38_final* $tmpdir
else
    echo "❌ KGref files is missing, needs to be prepared."
    # ---------------------------
    # Download 1000 Genomes reference (GRCh37)
    # ---------------------------
    wget https://www.kingrelatedness.com/ancestry/KGref.bed.xz
    wget https://www.kingrelatedness.com/ancestry/KGref.fam.xz
    wget https://www.kingrelatedness.com/ancestry/KGref.bim.xz

    ## Then command to uncompress :
    unxz *.xz

    for f in KGref.*; do
      base="${f%.*}"
      ext="${f##*.}"
      mv "$f" "${tmpdir}${base}_GRCh37_intermediate.${ext}"
    done

    # ---------------------------
    # Create unique SNP IDs for GRCh37 reference
    # ---------------------------
    awk -v OFS='\t' '{
      chr = $1
      sub(/^chr/, "", chr)      # Remove leading "chr" if present
      chr = "chr"chr            # Add "chr" prefix explicitly
      new_id = chr"_"$4"_"$5"_"$6
      print $2, new_id
      }' "${tmpdir}KGref_GRCh37_intermediate.bim" > "${tmpdir}name_to_update.tsv"

    plink --bfile ${tmpdir}KGref_GRCh37_intermediate --update-name ${tmpdir}name_to_update.tsv --make-bed --out ${tmpdir}KGref_GRCh37_final --memory ${mem_MB} --threads ${cpus}


    # ---------------------------
    # Optional: LiftOver to GRCh38
    # ---------------------------
    if [[ "$genome_version" == "GRCh38" ]]; then
        echo "➡ Converting 1000 Genomes reference from GRCh37 to GRCh38..."

        awk -v OFS='\t' '{print "chr"$1, $4, $4+1, $2, $3, $5, $6}' ${tmpdir}KGref_GRCh37_final.bim > ${tmpdir}KGref_GRCh37_map.bed

        wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

        # Check if liftOver is available
        if ! command -v liftOver &> /dev/null; then
            echo "liftOver not found — loading kentutils/453 module..."
            module load StdEnv/2020 kentutils/401 
        else
            echo "liftOver available."
        fi

        liftOver ${tmpdir}KGref_GRCh37_map.bed hg19ToHg38.over.chain.gz ${tmpdir}KGref_GRCh38.bed ${tmpdir}KGref_unmapped.bed
	rmhg19ToHg38.over.chain.gz
	
        awk -v OFS='\t' '{print $1}' ${tmpdir}KGref_GRCh37_map.bed | uniq -c
        awk -v OFS='\t' '{print $1, $4, $5, $2, $6, $7}' ${tmpdir}KGref_GRCh38.bed > ${tmpdir}KGref_GRCh38.bim
        awk -F'\t' '$1 ~ /^chr([0-9]+|X|Y|M)$/' ${tmpdir}KGref_GRCh38.bim  > ${tmpdir}KGref_GRCh38_filtered.bim


        #### Process reference genome (hg38)
        # Prepare chromosome and position mappings and SNP selections from hg38.
        cat ${tmpdir}KGref_GRCh38_filtered.bim | awk -v OFS='\t' '{print $2, $1}' > ${tmpdir}KGref_GRCh38_SNP_Chr_map.tsv
        cat ${tmpdir}KGref_GRCh38_filtered.bim | awk -v OFS='\t' '{print $2, $4}' > ${tmpdir}KGref_GRCh38_SNP_pos_map.tsv
        awk -v OFS='\t' '{print $2}' ${tmpdir}KGref_GRCh38_filtered.bim > ${tmpdir}KGref_GRCh38_SNP_selection.txt

        # Extract SNPs matching hg38 reference.
        plink --bfile ${tmpdir}KGref_GRCh37_final --extract ${tmpdir}KGref_GRCh38_SNP_selection.txt --make-bed --out ${tmpdir}KGref_GRCh38_filter --memory ${mem_MB} --threads ${cpus}

        # Update positions for hg38 SNPs.
        plink --bfile ${tmpdir}KGref_GRCh38_filter --update-map ${tmpdir}KGref_GRCh38_SNP_pos_map.tsv --make-bed --out ${tmpdir}KGref_GRCh38_filter_update_position --memory ${mem_MB} --threads ${cpus}

        # Update chromosomes for hg38 SNPs.
        plink --bfile ${tmpdir}KGref_GRCh38_filter_update_position --update-chr ${tmpdir}KGref_GRCh38_SNP_Chr_map.tsv --make-bed --out ${tmpdir}KGref_GRCh38_update_chr --memory ${mem_MB} --threads ${cpus}


        # Update name for hg38 SNPs.
        awk -v OFS='\t' '{
          chr = $1
          sub(/^chr/, "", chr)        # Remove leading "chr" if present
          chr = "chr"chr            # Add "chr" prefix explicitly
          new_id = chr"_"$4"_"$5"_"$6
          print $2, new_id
          }' "${tmpdir}KGref_GRCh38_update_chr.bim" > "${tmpdir}name_to_update.tsv"

        plink --bfile ${tmpdir}KGref_GRCh38_filter_update_position --update-name ${tmpdir}name_to_update.tsv --make-bed --out ${tmpdir}KGref_GRCh38_final --memory ${mem_MB} --threads ${cpus}
    fi
fi

# ---------------------------
# Update SNP names for input dataset
# ---------------------------
awk -v OFS='\t' '{
    chr = $1
    sub(/^chr/, "", chr)        # Remove leading "chr" if present
    chr = "chr"chr            # Add "chr" prefix explicitly
    new_id = chr"_"$4"_"$5"_"$6
    if (!seen[new_id]++) print $2, new_id
    }' "${plink_prefix}.bim" > "${tmpdir}update_SNPname.tsv"


plink --bfile "${plink_prefix}" \
        --update-name "${tmpdir}update_SNPname.tsv" \
        --make-bed \
        --out "${tmpdir}plink_updated_name" --memory ${mem_MB} --threads ${cpus}



# ---------------------------
# Prepare files for KING PCA
# ---------------------------
awk -v OFS='\t' '{print $2}' ${tmpdir}KGref_${genome_version}_final.bim > ${tmpdir}KGref_${genome_version}_SNP_selection.txt
plink --bfile "${tmpdir}plink_updated_name" --extract ${tmpdir}KGref_${genome_version}_SNP_selection.txt --make-bed --out ${tmpdir}plink_for_pcancestry --memory ${mem_MB} --threads ${cpus}

awk -v OFS='\t' '{print $2}' ${tmpdir}plink_for_pcancestry.bim > ${tmpdir}snp_to_select.txt
plink --bfile "${tmpdir}KGref_${genome_version}_final" --extract ${tmpdir}snp_to_select.txt --make-bed --out ${tmpdir}kingref_for_pcancestry --memory ${mem_MB} --threads ${cpus}



# ---------------------------
# Run KING PCA
# ---------------------------

# ---------------------------
# Check R availability (must be 3.6.x)
# ---------------------------
if ! command -v R &> /dev/null; then
    echo "❌ R not found — trying to load R 3.6.0..."
    module load nixpkgs/16.09 gcc/7.3.0
    module load openmpi/3.1.2
    module load r/3.6.0
    echo "✅ R 3.6.0 module loaded."
else
    R_VERSION=$(R --version | head -n 1 | grep -oE '[0-9]+\.[0-9]+')
    if [[ "$R_VERSION" == "3.6" ]]; then
        echo "✅ R $R_VERSION detected."
    else
        echo "⚠️ R version $R_VERSION detected — switching to R 3.6.0..."
        module load nixpkgs/16.09 gcc/7.3.0
        module load openmpi/3.1.2
        module load r/3.6.0
        echo "✅ R 3.6.0 module loaded."
    fi
fi

king -b ${tmpdir}kingref_for_pcancestry.bed,${tmpdir}plink_for_pcancestry.bed --pca --projection --rplot --prefix ${tmpdir}king_output --cpus ${cpus}



awk -v OFS='\t' 'NR==FNR {
        if (FNR==1) {next}              # skip ancestry header
        a[$2]=$9                        # store ancestry column by IID
        next
     }
     FNR==1 {
        # Print header: change IID → SampleID, add Ancestry
        $2="SampleID"
        print $2, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16 , "Ancestry"
        next
     }
     ($2 in a) {
        $2=$2                           # just to ensure SampleID prints correctly
        print $2, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16 , a[$2]
     }' OFS='\t' ${tmpdir}king_output_InferredAncestry.txt ${tmpdir}king_outputpc.txt > ${output_file}


rm -rf "$tmpdir"
echo "✅ KING PCA analysis completed."
