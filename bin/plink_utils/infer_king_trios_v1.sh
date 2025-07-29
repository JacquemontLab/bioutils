#!/bin/bash


# ------------------------------------------------------------------------------
# Script Name: extract_king_trios.sh
# Author: Florian Bénitière
# Date: 2025-07-28
#
# Description:
#   Runs KING on a PLINK binary dataset and parses the KING relationship log file
#   to extract trio information (FamilyID, SampleID, FatherID, MotherID)
#
# Usage:
#   ./extract_king_trios.sh plink_prefix [output_tsv_file]
#
# Arguments:
#   plink_prefix       Prefix of the PLINK BED dataset (e.g., data for data.bed)
#   output_tsv_file    (Optional) Output TSV file path. Default: king_trios.tsv
#
# Output:
#   A tab-separated file with the columns:
#     FamilyID   SampleID   FatherID   MotherID
# ------------------------------------------------------------------------------


set -e

# --------------------- Argument Parsing ---------------------
if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <plink_prefix> [output_tsv_file]"
  exit 1
fi

plink_bed="$1"
OUTPUT_TSV="${2:-king_trios.tsv}"

# --------------------- Check if KING is available ---------------------
if ! command -v king &> /dev/null; then
  echo "❌ Error: 'king' is not available in your PATH. Please install KING or load the module."
  exit 1
fi

# --------------------- Detect Number of CPUs ---------------------
cpus="${SLURM_CPUS_ON_NODE:-$(nproc)}"
echo "💻 Running KING with $cpus cores..."

# --------------------- Run KING ---------------------
king -b "${plink_bed}.bed" --build --degree 1 --cpus "$cpus" --prefix temp_king_out

INPUT_LOG="temp_king_outbuild.log"

# --------------------- Parse KING Output ---------------------
awk '
  {
    lines[NR] = $0
    if ($0 ~ /^Family KING[0-9]+:/) {
      family_at_line[NR] = $2
      sub(/:$/, "", family_at_line[NR])
    }
  }
  END {
    for (i = 2; i <= NR; i++) {
      if (lines[i] ~ /father/) {
        info_line = lines[i - 1]
        fam = ""
        for (j = i - 1; j >= 1; j--) {
          if (family_at_line[j] != "") {
            fam = family_at_line[j]
            break
          }
        }
        if (fam == "") fam = "NA"
        if (match(info_line, /\(P=([^,]+), O=([^,]+), R=([^)\)]+)\)/, arr)) {
          print fam "\t" arr[2] "\t" arr[3] "\t" arr[1]
        }
      }
    }
  }
' "$INPUT_LOG" | (
  echo -e "FamilyID\tSampleID\tFatherID\tMotherID"
  cat -
) > "$OUTPUT_TSV"

rm -f temp_king_outupdateparents.txt temp_king_outupdateids.txt temp_king_outallsegs.txt temp_king_outX.kin0 temp_king_outbuild.log temp_king_output_pcplot.Rout temp_king_output_ancestryplot.Rout

echo "✅ Trio table successfully written to: $OUTPUT_TSV"

