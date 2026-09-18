#!/usr/bin/env bash
set -uo pipefail

BASE=${PWD}

if [[ "${CC_CLUSTER:-}" == "rorqual" ]]; then
    COUNT_SCRIPT=/lustre09/project/6008022/LAB_WORKSPACE/SOFTWARE/bioutils/bin/manage_utils/count_sample_h5ad.py
elif [[ "$(hostname -s)" == "chusj-transcriptomic-server-1" ]]; then
    COUNT_SCRIPT=/mnt/chusj-transcriptomic-cephfs-1/LAB_WORKSPACE/SOFTWARE/bioutils/bin/manage_utils/count_sample_h5ad.py
else
    COUNT_SCRIPT=''
fi

shopt -s nullglob   # empty globs expand to nothing, not to the literal pattern

# ANSI styling (disabled when output isn't a terminal so pipes/files stay clean)
if [ -t 1 ]; then B=$'\033[1m'; R=$'\033[0m'; DIM=$'\033[2m'; else B=''; R=''; DIM=''; fi

# ---------------------------------------------------------------------
# One Python process for the whole tree.
# TSV columns (always 9):
#   DIR NFILES NSAMPLES N_TISSUE TISSUE_VALS N_REGION REGION_VALS N_SUBREGION SUBREGION_VALS
# A trailing "#TOTAL<TAB>N_TISSUE<TAB>N_REGION<TAB>N_SUBREGION" line carries
# exact tree-wide distinct counts.
# ---------------------------------------------------------------------
declare -A NFILES NSAMPLES
declare -A NTISSUE TISSUEVALS NREGION REGIONVALS NSUBREG SUBREGVALS
tot_tissue='?' ; tot_region='?' ; tot_subreg='?'

modality_dirs=( "$BASE"/*/expression/*/*/ )

if (( ${#modality_dirs[@]} > 0 )); then
    while IFS=$'\t' read -r c1 c2 c3 c4 c5 c6 c7 c8 c9; do
        if [[ "$c1" == "#TOTAL" ]]; then
            tot_tissue=$c2; tot_region=$c3; tot_subreg=$c4
            continue
        fi
        # c1=dir c2=nfiles c3=nsamples c4=ntissue c5=tissuevals
        # c6=nregion c7=regionvals c8=nsubreg c9=subregvals
        NFILES["$c1"]=$c2;   NSAMPLES["$c1"]=$c3
        NTISSUE["$c1"]=$c4;  TISSUEVALS["$c1"]=$c5
        NREGION["$c1"]=$c6;  REGIONVALS["$c1"]=$c7
        NSUBREG["$c1"]=$c8;  SUBREGVALS["$c1"]=$c9
    done < <(python "$COUNT_SCRIPT" --batch "${modality_dirs[@]}")
fi

echo
echo "🧬 ================================================================"
echo "🧬                 EXPRESSION DATABASE SUMMARY"
echo "🧬 ================================================================"
echo

total_h5ad=0
total_datasets=0

# emoji  label  count  values  -> "   │      emoji Label  N: v1, v2, …"
# Skips the line entirely when the field is NA/empty.
print_annot() {
    local emoji=$1 label=$2 count=$3 vals=$4
    [[ "$count" == "NA" || -z "$count" ]] && return
    printf "   │      %s %s%-10s%s %s: %s\n" "$emoji" "$DIM" "$label" "$R" "$count" "$vals"
}

for cohort_dir in "$BASE"/*/expression/; do
    [[ -d "$cohort_dir" ]] || continue

    cohort=$(basename "$(dirname "$cohort_dir")")
    echo "📦 ${B}${cohort}${R}"

    for project_dir in "$cohort_dir"*/; do
        [[ -d "$project_dir" ]] || continue
        project=$(basename "$project_dir")

        for modality_dir in "$project_dir"*/; do
            [[ -d "$modality_dir" ]] || continue
            modality=$(basename "$modality_dir")

            # Looked up from the batch pass; dirs with no h5ad are absent.
            n=${NFILES["$modality_dir"]:-0}
            (( n > 0 )) || continue

            case "$modality" in
                bulk_RNAseq) emoji="🧬" ;;
                snRNAseq)    emoji="🔬" ;;
                snATACseq)   emoji="🧪" ;;
                *)           emoji="📁" ;;
            esac

            n_samples=${NSAMPLES["$modality_dir"]:-?}

            printf "   ├── 📂 %-28s %s %-13s ✅ %3d h5ad / %s samples\n" \
                "$project" "$emoji" "$modality" "$n" "$n_samples"

            print_annot "" "Tissue"    "${NTISSUE[$modality_dir]:-NA}"  "${TISSUEVALS[$modality_dir]:-}"
            print_annot "" "Region"    "${NREGION[$modality_dir]:-NA}"  "${REGIONVALS[$modality_dir]:-}"
            print_annot "" "Subregion" "${NSUBREG[$modality_dir]:-NA}"  "${SUBREGVALS[$modality_dir]:-}"

            total_h5ad=$((total_h5ad + n))
            total_datasets=$((total_datasets + 1))
        done
    done
    echo
done

echo "🧬 ================================================================"
printf "📊 Total h5ad files:      %3d\n" "$total_h5ad"
printf "📂 Total data sets:       %3d\n" "$total_datasets"
printf " Distinct tissues:      %3s\n" "$tot_tissue"
printf " Distinct regions:      %3s\n" "$tot_region"
printf " Distinct subregions:   %3s\n" "$tot_subreg"
echo "🧬 ================================================================"
echo