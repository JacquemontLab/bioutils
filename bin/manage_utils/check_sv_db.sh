#!/usr/bin/env bash
# check_sv_db.sh — count SampleID in curated/unfiltered/sampleDB per dataset
#                   create sampleDB.tsv (SampleID only) if absent
#                   + Cohort/Biobank breakdown from sampleDB when present
#   --fast : skip the unfiltered parquet (curated only)
set -euo pipefail

FAST=0
ARGS=()
for a in "$@"; do
    case "$a" in
        --fast) FAST=1 ;;
        *)      ARGS+=("$a") ;;
    esac
done
ROOT="${ARGS[0]:-.}"

# ANSI bold (disabled when output isn't a terminal so pipes/files stay clean)
if [ -t 1 ]; then B=$'\033[1m'; R=$'\033[0m'; else B=''; R=''; fi

# read_parquet() source list: curated (+ unfiltered unless --fast)
srclist() { # $1=cur $2=unf
    if (( FAST )); then echo "'$1'"; else echo "['$1', '$2']"; fi
}

for cur in "$ROOT"/*/sh*/*/ShortVariantsDB_curated.parquet; do
    dir=$(dirname "$cur")
    smp="$dir/sampleDB.tsv"
    unf="$dir/ShortVariantsDB_unfiltered.parquet"

    if [[ ! -f "$smp" ]]; then
        duckdb -noheader -list <<CREATE > "$smp"
SELECT 'SampleID';
SELECT DISTINCT SampleID
FROM read_parquet($(srclist "$cur" "$unf"))
WHERE SampleID IS NOT NULL
ORDER BY SampleID;
CREATE
        echo "🆕 created sampleDB.tsv → ${B}$(basename "$(dirname "$(dirname "$dir")")")${R}/$(basename "$dir")" >&2
    fi
done

# ── summary table ──────────────────────────────────────────────────────────
duckdb -noheader -list <<SQL | column -t -s'|'
$(
for cur in "$ROOT"/*/sh*/*/ShortVariantsDB_curated.parquet; do
    dir=$(dirname "$cur")
    smp="$dir/sampleDB.tsv"
    unf="$dir/ShortVariantsDB_unfiltered.parquet"
    name=$(basename "$dir")
    cohort=$(basename "$(dirname "$(dirname "$dir")")")
    src=$(srclist "$cur" "$unf")
    if (( FAST )); then unf_col=""; else
        unf_col="'|🧪 unfiltered=' || (SELECT COUNT(DISTINCT SampleID) FROM read_parquet('$unf')) ||"
    fi
    cat <<INNER
SELECT
  '📁 ${B}$cohort${R}/$name' ||
  '|🔬 curated=' || (SELECT COUNT(DISTINCT SampleID) FROM read_parquet('$cur')) ||
  $unf_col
  '|📋 sampleDB=' || s.n ||
  '|' || CASE WHEN m.missing = 0 THEN '✅ all in sampleDB'
              ELSE '⚠️ ' || m.missing || ' missing' END
FROM
  (SELECT COUNT(DISTINCT SampleID) AS n
   FROM read_csv('$smp', delim='\t', header=true)) s,
  (SELECT COUNT(DISTINCT SampleID) FILTER (
            WHERE CAST(SampleID AS VARCHAR) NOT IN (
              SELECT CAST(SampleID AS VARCHAR)
              FROM read_csv('$smp', delim='\t', header=true)
              WHERE SampleID IS NOT NULL)
          ) AS missing
   FROM read_parquet($src)) m;
INNER
done
)
SQL

# ── Cohort × Biobank breakdown (only where both columns exist) ──────────────
for cur in "$ROOT"/*/sh*/*/ShortVariantsDB_curated.parquet; do
    dir=$(dirname "$cur")
    smp="$dir/sampleDB.tsv"
    name=$(basename "$dir")
    cohort=$(basename "$(dirname "$(dirname "$dir")")")
    [[ -f "$smp" ]] || continue

    cols=$(duckdb -noheader -list -c \
        "SELECT column_name FROM (DESCRIBE SELECT * FROM read_csv('$smp', delim='\t', header=true));")

    miss=()
    grep -qxF Cohort  <<<"$cols" || miss+=(Cohort)
    grep -qxF Biobank <<<"$cols" || miss+=(Biobank)
    if (( ${#miss[@]} )); then
        echo
        echo "🗂️  ${B}$cohort${R}/$name — ❌ no breakdown ($(IFS=+; echo "${miss[*]}") not present)"
        continue
    fi

    echo
    echo "🗂️  ${B}$cohort${R}/$name — Cohort × Biobank"
    duckdb -noheader -list <<SQL | column -t -s'|'
SELECT '  ' || Cohort || '|' || Biobank || '|' || COUNT(*) AS row
FROM read_csv('$smp', delim='\t', header=true)
GROUP BY Cohort, Biobank
ORDER BY Cohort, Biobank;
SQL
done