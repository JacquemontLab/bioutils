#!/usr/bin/env bash
# summarize_database.sh — presentation-ready overview of the DATABASE_v2 tree.
#
# Aggregates UNIQUE SampleID counts across three data types per cohort/dataset:
#   • short variants  (ShortVariantsDB_curated.parquet)   → DuckDB
#   • CNV             (cnvDB.parquet)                      → DuckDB
#   • expression      (*.h5ad)                             → count_sample_h5ad.py
#
# Produces:
#   1. A per-dataset breakdown (grouped by cohort)
#   2. A cohort × datatype matrix (unique samples)
#   3. Grand totals
#
# Requirements: duckdb, python (with anndata) + count_sample_h5ad.py
#
# Usage:
#   ./summarize_database.sh [ROOT]
#   ./summarize_database.sh /mnt/.../DATABASE_v2
#
# Notes / pitfalls:
#   * "Unique samples" per datatype are counted independently; the same
#     individual may appear in several datatypes (this is expected and is a
#     union across datasets *within* a datatype, per cohort).
#   * konopka_abn encodes SampleIDs in the h5ad filename (ABN-XXXX_...),
#     matching the special-case in the original check_expression.sh.
#   * Counts are DISTINCT SampleID; NULLs are ignored.
set -uo pipefail

ROOT="${1:-.}"
COUNT_SCRIPT="${COUNT_SCRIPT:-/home/flben/bin/count_sample_h5ad.py}"

# ── ANSI styling (disabled when not a terminal so pipes/files stay clean) ────
if [ -t 1 ]; then
    B=$'\033[1m'; R=$'\033[0m'; DIM=$'\033[2m'; CY=$'\033[36m'; GR=$'\033[32m'
else
    B=''; R=''; DIM=''; CY=''; GR=''
fi

# ── helpers ─────────────────────────────────────────────────────────────────

# Count DISTINCT SampleID in a parquet dataset via DuckDB. Echoes an integer.
# Handles both a single .parquet file and a Hive-partitioned directory
# (e.g. .../ShortVariantsDB_curated.parquet/CHROM=chr1/*.parquet).
count_parquet_samples() {
    local pq="$1" src
    if [[ -d "$pq" ]]; then
        # Partitioned dataset: glob recursively, enable Hive partitioning.
        src="read_parquet('${pq}/**/*.parquet', hive_partitioning=true, union_by_name=true)"
    else
        src="read_parquet('$pq')"
    fi
    duckdb -noheader -list <<SQL 2>/dev/null
SELECT COUNT(DISTINCT SampleID)
FROM $src
WHERE SampleID IS NOT NULL;
SQL
}

# Count unique expression samples in a modality dir. Echoes an integer or NA.
count_expression_samples() {
    local modality_dir="$1" project="$2"
    if [[ "$project" == "konopka_abn" ]]; then
        # SampleID encoded in filename: ABN-XXXX_<region>_...
        find "$modality_dir" -maxdepth 1 -type f -name "*.h5ad" -printf "%f\n" \
            | sed -E 's/^(ABN-[^-_]+)_.*/\1/' | sort -u | wc -l
    else
        python "$COUNT_SCRIPT" "$modality_dir"/*.h5ad 2>/dev/null
    fi
}

# Emoji for an expression modality.
modality_emoji() {
    case "$1" in
        bulk_RNAseq) echo "🧬" ;;
        snRNAseq)    echo "🔬" ;;
        snATACseq)   echo "🧪" ;;
        *)           echo "📁" ;;
    esac
}

# ── accumulators (assoc arrays keyed by cohort) ─────────────────────────────
declare -A SV_BY_COHORT   # short variants unique samples
declare -A CNV_BY_COHORT  # cnv unique samples
declare -A EXP_BY_COHORT  # expression unique samples
declare -A DS_BY_COHORT   # dataset count
declare -A COHORT_SEEN

sv_total=0; cnv_total=0; exp_total=0; ds_total=0

add() { # $1=assoc-name $2=key $3=int  (safe add, treats non-int as 0)
    local -n arr="$1"; local k="$2" v="$3"
    [[ "$v" =~ ^[0-9]+$ ]] || v=0
    arr["$k"]=$(( ${arr["$k"]:-0} + v ))
}

echo
echo "🧬 ================================================================"
echo "🧬                  ${B}DATABASE — SAMPLE SUMMARY${R}"
echo "🧬 ================================================================"

# ── discover cohorts (top-level dirs under ROOT) ────────────────────────────
mapfile -t COHORTS < <(find "$ROOT" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort)

for cohort in "${COHORTS[@]}"; do
    printed_header=0
    cohort_dir="$ROOT/$cohort"

    _hdr() { (( printed_header )) || { echo; echo "📦 ${B}${cohort}${R}"; printed_header=1; }; }

    # ── short variants ──────────────────────────────────────────────────────
    while IFS= read -r cur; do
        # find already confirmed the path; use -e (not -f) so symlinks and
        # CephFS regular-file quirks don't silently drop a dataset.
        [[ -e "$cur" ]] || continue
        dir=$(dirname "$cur"); name=$(basename "$dir")
        n=$(count_parquet_samples "$cur"); [[ "$n" =~ ^[0-9]+$ ]] || n=0
        _hdr
        printf "   ├── 🧾 %-22s 🔬 %-13s ✅ %5s samples\n" "$name" "shortvariants" "$n"
        add SV_BY_COHORT "$cohort" "$n"
        add DS_BY_COHORT "$cohort" 1
        COHORT_SEEN["$cohort"]=1
        sv_total=$((sv_total + n)); ds_total=$((ds_total + 1))
    done < <(find "$cohort_dir" -path '*/sh*/*/ShortVariantsDB_curated.parquet' 2>/dev/null | sort)

    # ── CNV ─────────────────────────────────────────────────────────────────
    while IFS= read -r cnv; do
        [[ -e "$cnv" ]] || continue
        dir=$(dirname "$cnv"); name=$(basename "$dir")
        n=$(count_parquet_samples "$cnv"); [[ "$n" =~ ^[0-9]+$ ]] || n=0
        _hdr
        printf "   ├── 🧬 %-22s 🧬 %-13s ✅ %5s samples\n" "$name" "cnv" "$n"
        add CNV_BY_COHORT "$cohort" "$n"
        add DS_BY_COHORT "$cohort" 1
        COHORT_SEEN["$cohort"]=1
        cnv_total=$((cnv_total + n)); ds_total=$((ds_total + 1))
    done < <(find "$cohort_dir" -path '*/cn*/*/cnvDB.parquet' 2>/dev/null | sort)

    # ── expression ──────────────────────────────────────────────────────────
    if [[ -d "$cohort_dir/expression" ]]; then
        while IFS= read -r project_dir; do
            project=$(basename "$project_dir")
            while IFS= read -r modality_dir; do
                modality=$(basename "$modality_dir")
                n_h5ad=$(find "$modality_dir" -maxdepth 1 -type f -name "*.h5ad" | wc -l)
                (( n_h5ad > 0 )) || continue
                emoji=$(modality_emoji "$modality")
                n=$(count_expression_samples "$modality_dir" "$project")
                _hdr
                printf "   ├── 📂 %-22s %s %-13s ✅ %3d h5ad / %5s samples\n" \
                    "$project" "$emoji" "$modality" "$n_h5ad" "$n"
                add EXP_BY_COHORT "$cohort" "$n"
                add DS_BY_COHORT "$cohort" 1
                COHORT_SEEN["$cohort"]=1
                [[ "$n" =~ ^[0-9]+$ ]] && exp_total=$((exp_total + n))
                ds_total=$((ds_total + 1))
            done < <(find "$project_dir" -mindepth 1 -maxdepth 1 -type d | sort)
        done < <(find "$cohort_dir/expression" -mindepth 1 -maxdepth 1 -type d | sort)
    fi
done

# ── cohort × datatype matrix ────────────────────────────────────────────────
echo
echo "🧬 ================================================================"
echo "📊                 ${B}COHORT × DATATYPE  (unique samples)${R}"
echo "🧬 ================================================================"
{
    echo "Cohort|🔬 ShortVar|🧬 CNV|📂 Expression|📦 Datasets"
    echo "──────|──────────|─────|────────────|──────────"
    for cohort in "${COHORTS[@]}"; do
        [[ -n "${COHORT_SEEN[$cohort]:-}" ]] || continue
        printf "%s|%s|%s|%s|%s\n" \
            "${B}${cohort}${R}" \
            "${SV_BY_COHORT[$cohort]:-·}" \
            "${CNV_BY_COHORT[$cohort]:-·}" \
            "${EXP_BY_COHORT[$cohort]:-·}" \
            "${DS_BY_COHORT[$cohort]:-0}"
    done
    echo "──────|──────────|─────|────────────|──────────"
    printf "${B}TOTAL${R}|${B}%s${R}|${B}%s${R}|${B}%s${R}|${B}%s${R}\n" \
        "$sv_total" "$cnv_total" "$exp_total" "$ds_total"
} | column -t -s'|'

# ── grand totals ────────────────────────────────────────────────────────────
echo
echo "🧬 ================================================================"
printf "🔬 Short-variant samples : ${B}%5d${R}\n" "$sv_total"
printf "🧬 CNV samples           : ${B}%5d${R}\n" "$cnv_total"
printf "📂 Expression samples    : ${B}%5d${R}\n" "$exp_total"
printf "📦 Total datasets        : ${B}%5d${R}\n" "$ds_total"
echo "🧬 ================================================================"
echo