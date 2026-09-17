#!/usr/bin/env bash
# check_cnv_db.sh — verify SampleID coverage + PC1..PC10 presence per dataset
#                   + Cohort/Biobank breakdown from sampleDB when present
set -euo pipefail

ROOT="${1:-.}"

# ANSI bold (disabled when not a terminal so pipes stay clean)
if [ -t 1 ]; then B=$'\033[1m'; R=$'\033[0m'; else B=''; R=''; fi

# ── summary table ──────────────────────────────────────────────────────────
duckdb -noheader -list <<SQL | column -t -s'|'
$(
for cnv in "$ROOT"/*/cn*/*/cnvDB.parquet; do
    dir=$(dirname "$cnv")
    smp="$dir/sampleDB.tsv"
    name=$(basename "$dir")
    cohort=$(basename "$(dirname "$(dirname "$dir")")")
    [[ -f "$smp" ]] || { echo "SELECT '📁 ${B}$cohort${R}/$name |❌ sampleDB.tsv missing| | | ';"; continue; }
    cat <<INNER
SELECT
  '📁 ${B}$cohort${R}/$name' ||
  '|🧬 cnvDB=' || c.n ||
  '|📋 sampleDB=' || s.n ||
  '|' || CASE WHEN c.missing = 0 THEN '✅ all in sampleDB'
              ELSE '⚠️ ' || c.missing || ' missing' END ||
  '|' || CASE WHEN pc.ok THEN '✅ PC ancestry' ELSE '❌ PC ancestry' END
FROM
  (SELECT COUNT(DISTINCT SampleID) AS n,
          COUNT(DISTINCT SampleID) FILTER (
            WHERE SampleID NOT IN (
              SELECT SampleID FROM read_csv('$smp', delim='\t', header=true)
              WHERE SampleID IS NOT NULL)
          ) AS missing
   FROM read_parquet('$cnv')) c,
  (SELECT COUNT(DISTINCT SampleID) AS n
   FROM read_csv('$smp', delim='\t', header=true)) s,
  (SELECT COUNT(*) = 10 AS ok
   FROM (DESCRIBE SELECT * FROM read_csv('$smp', delim='\t', header=true))
   WHERE column_name IN ('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')) pc;
INNER
done
)
SQL

# ── Cohort × Biobank breakdown (only where both columns exist) ──────────────
for cnv in "$ROOT"/*/cn*/*/cnvDB.parquet; do
    dir=$(dirname "$cnv")
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