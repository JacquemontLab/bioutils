#!/usr/bin/env bash
# Generate the three CNV documentation artefacts for every dataset under
#   $DB/<cohort>/cnv/<variant>/
# and place them in that dataset's docs/ directory:
#
#   1. cnvDB_columns_report.pdf   (pdf_columns_report.py on cnvDB.parquet)
#   2. geneDB_columns_report.pdf  (pdf_columns_report.py on geneDB.parquet)
#   3. cnv_dataset_qc.pdf         (cnv_dataset_report.Rmd via apptainer)
#
# Steps run sequentially: each already uses $CPUS cores / $TOTAL_MEMORY GB.
# Failures are recorded and the loop continues; exit status is non-zero if
# any step failed.
#
# Usage:
#   ./run_cnv_reports.sh [options]
#     --force            regenerate even if the docs/ output already exists
#     --dry-run          print what would run, execute nothing
#     --only PATTERN     restrict to datasets whose relative path matches
#                        PATTERN (bash glob, e.g. 'gtex/*' or '*wgs_dbgap*')
#     --verbose          stream tool output to the terminal as well as the log
#     --cpus N           override CPU count       (default 16)
#     --mem N            override memory in GB    (default 150)
#     --logdir DIR       override log directory
set -uo pipefail

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
if [[ "${CC_CLUSTER:-}" == "rorqual" ]]; then
    BIND=/lustre09/project/6008022
elif [[ "$(hostname -s)" == "chusj-transcriptomic-server-1" ]]; then
    BIND=/mnt/chusj-transcriptomic-cephfs-1
else
    echo "Unknown site: $(hostname -s)" >&2
    return 1 2>/dev/null || exit 1
fi
export BIND

BASE="$BIND/LAB_WORKSPACE"
DB="$BASE/DATABASE"
SIF="$BASE/SOFTWARE/Dockers/cnv_dataset_report_latest.sif"
RMD="$BASE/SOFTWARE/bioutils/bin/cnv_utils/cnv_dataset_report.Rmd"
COLUMNS_REPORT="$BASE/SOFTWARE/Git_pipeline/CNV-Annotation/bin/pdf_columns_report.py"
CNVANN_REPO="$BASE/SOFTWARE/Git_pipeline/CNV-Annotation/"
CNVCALLER_REPO="$BASE/SOFTWARE/Git_pipeline/CNV-Caller/"

CPUS=16
TOTAL_MEMORY=150
FORCE=0
DRYRUN=0
VERBOSE=0
ONLY='*'
LOGDIR=""

# --------------------------------------------------------------------------
# Argument parsing
# --------------------------------------------------------------------------
while [ $# -gt 0 ]; do
  case "$1" in
    --force)   FORCE=1; shift ;;
    --dry-run) DRYRUN=1; shift ;;
    --verbose) VERBOSE=1; shift ;;
    --only)    ONLY="${2:?--only needs a pattern}"; shift 2 ;;
    --cpus)    CPUS="${2:?--cpus needs a number}"; shift 2 ;;
    --mem)     TOTAL_MEMORY="${2:?--mem needs a number}"; shift 2 ;;
    --logdir)  LOGDIR="${2:?--logdir needs a path}"; shift 2 ;;
    -h|--help) sed -n '2,25p' "$0"; exit 0 ;;
    *) printf 'unknown option: %s\n' "$1" >&2; exit 2 ;;
  esac
done

RUN_TS="$(date +%Y%m%d_%H%M%S)"
LOGDIR="${LOGDIR:-$PWD/logs_cnv_reports_$RUN_TS}"

# --------------------------------------------------------------------------
# Pretty printing
# --------------------------------------------------------------------------
if [ -t 1 ]; then
  B=$'\033[1m'; DIM=$'\033[2m'; RED=$'\033[31m'; GRN=$'\033[32m'
  YLW=$'\033[33m'; CYN=$'\033[36m'; R=$'\033[0m'
else
  B=''; DIM=''; RED=''; GRN=''; YLW=''; CYN=''; R=''
fi

hr()      { printf '%s\n' "${DIM}$(printf '─%.0s' {1..78})${R}"; }
banner()  { hr; printf '%s%s%s\n' "$B" "$1" "$R"; hr; }
info()    { printf '   %sℹ%s  %s\n' "$CYN" "$R" "$1"; }
ok()      { printf '   %s✅%s %-34s %s%s%s\n' "$GRN" "$R" "$1" "$DIM" "$2" "$R"; }
skipped() { printf '   %s⏭%s  %-34s %s%s%s\n' "$YLW" "$R" "$1" "$DIM" "$2" "$R"; }
failed()  { printf '   %s❌%s %-34s %s%s%s\n' "$RED" "$R" "$1" "$DIM" "$2" "$R"; }

fmt_dur() {  # seconds -> MMmSSs
  local s="$1"
  printf '%dm%02ds' $((s / 60)) $((s % 60))
}

# --------------------------------------------------------------------------
# Pre-flight
# --------------------------------------------------------------------------
banner "CNV dataset reports — $RUN_TS"

preflight_fail=0
for f in "$SIF" "$RMD" "$COLUMNS_REPORT"; do
  if [ -e "$f" ]; then
    info "found $(basename "$f")"
  else
    failed "missing prerequisite" "$f"
    preflight_fail=1
  fi
done
[ -d "$DB" ] || { failed "missing DATABASE root" "$DB"; preflight_fail=1; }
[ "$preflight_fail" -eq 0 ] || exit 1

# Collect datasets: any <cohort>/cnv/<variant>/ holding a cnvDB.parquet.
datasets=()
shopt -s nullglob
for d in "$DB"/*/cnv/*/; do
  d="${d%/}"
  rel="${d#"$DB"/}"
  [ -f "$d/cnvDB.parquet" ] || continue
  # shellcheck disable=SC2254  # ONLY is intentionally a glob
  case "$rel" in $ONLY) datasets+=("$d") ;; esac
done
shopt -u nullglob

n_ds=${#datasets[@]}
if [ "$n_ds" -eq 0 ]; then
  failed "no datasets matched" "$DB/*/cnv/*/ (--only '$ONLY')"
  exit 1
fi

info "$n_ds dataset(s) to process · cpus=$CPUS · mem=${TOTAL_MEMORY}G"
info "logs: $LOGDIR"
[ "$DRYRUN" -eq 1 ] && info "${YLW}dry-run: nothing will be executed${R}"
[ "$DRYRUN" -eq 1 ] || mkdir -p "$LOGDIR"
echo

# --------------------------------------------------------------------------
# Step runner: run_step <label> <logfile> <cmd...>
# Returns the command's exit status; output goes to the log (and stdout if
# --verbose).
# --------------------------------------------------------------------------
run_step() {
  local label="$1" log="$2"; shift 2
  if [ "$DRYRUN" -eq 1 ]; then
    printf '   %s$%s %s\n' "$DIM" "$R" "$*"
    return 0
  fi
  if [ "$VERBOSE" -eq 1 ]; then
    "$@" 2>&1 | tee "$log"
    return "${PIPESTATUS[0]}"
  fi
  "$@" > "$log" 2>&1
}

show_log_tail() {
  local log="$1"
  [ -f "$log" ] || return 0
  printf '      %s--- last 15 lines of %s ---%s\n' "$DIM" "$log" "$R"
  tail -n 15 "$log" | sed "s/^/      ${DIM}| ${R}/"
}

# --------------------------------------------------------------------------
# Main loop
# --------------------------------------------------------------------------
results=()   # "rel|step|status|duration"
n_fail=0
n_done=0
n_skip=0
i=0

for ds in "${datasets[@]}"; do
  i=$((i + 1))
  rel="${ds#"$DB"/}"
  docs="$ds/docs"
  safe="${rel//\//_}"

  # dataset_name = <cohort>_<variant>, e.g. abn/cnv/wgs_sanders -> abn_wgs_sanders
  IFS='/' read -r p_cohort _p_datatype p_variant <<< "$rel"
  ds_name="${p_cohort}_${p_variant}"

  printf '%s[%d/%d]%s %s%s%s %s(%s)%s\n' \
    "$B" "$i" "$n_ds" "$R" "$B" "$rel" "$R" "$DIM" "$ds_name" "$R"
  [ "$DRYRUN" -eq 1 ] || mkdir -p "$docs"

  # ---- steps 1 & 2: column reports for cnvDB / geneDB -------------------
  for base in cnvDB geneDB; do
    parquet="$ds/$base.parquet"
    out="$base"_columns_report.pdf
    target="$docs/$out"
    label="$out"
    log="$LOGDIR/${safe}__${base}_columns_report.log"

    if [ ! -f "$parquet" ]; then
      skipped "$label" "no $base.parquet"
      results+=("$rel|$label|skip(no input)|-")
      n_skip=$((n_skip + 1))
      continue
    fi
    if [ -f "$target" ] && [ "$FORCE" -eq 0 ]; then
      skipped "$label" "exists (use --force)"
      results+=("$rel|$label|skip(exists)|-")
      n_skip=$((n_skip + 1))
      continue
    fi

    t0=$SECONDS
    if run_step "$label" "$log" \
         "$COLUMNS_REPORT" "$parquet" "$CPUS" "$TOTAL_MEMORY"; then
      # the script writes <base>_columns_report.pdf next to the parquet
      if [ "$DRYRUN" -eq 1 ]; then
        printf '   %s$%s mv %s %s\n' "$DIM" "$R" "$ds/$out" "$target"
        results+=("$rel|$label|dry-run|-")
      elif [ -f "$ds/$out" ] && mv -f "$ds/$out" "$target"; then
        d=$(fmt_dur $((SECONDS - t0)))
        ok "$label" "$d"
        results+=("$rel|$label|ok|$d")
        n_done=$((n_done + 1))
      else
        failed "$label" "script succeeded but $out not produced"
        show_log_tail "$log"
        results+=("$rel|$label|FAIL(no output)|-")
        n_fail=$((n_fail + 1))
      fi
    else
      failed "$label" "exit $? — see log"
      show_log_tail "$log"
      results+=("$rel|$label|FAIL|-")
      n_fail=$((n_fail + 1))
    fi
  done

  # ---- step 3: dataset QC report via apptainer ---------------------------
  label="cnv_dataset_qc.pdf"
  target="$docs/$label"
  log="$LOGDIR/${safe}__cnv_dataset_qc.log"

  if [ -f "$target" ] && [ "$FORCE" -eq 0 ]; then
    skipped "$label" "exists (use --force)"
    results+=("$rel|$label|skip(exists)|-")
    n_skip=$((n_skip + 1))
  else
    # Built as a variable so the dataset path is interpolated by bash, not R.
    r_expr=$(
      cat <<EOF
rmarkdown::render(
  "$RMD",
  params = list(
    path_dataset            = "$ds",
    dataset_name            = "$ds_name",
    path_CNVANNOTATION_repo = "$CNVANN_REPO",
    path_CNVCALLER_repo     = "$CNVCALLER_REPO"
  ),
  output_file       = "$target",
  intermediates_dir = tempdir(),
  quiet             = TRUE
)
EOF
    )

    t0=$SECONDS
    if run_step "$label" "$log" \
         apptainer exec -B "$BIND" "$SIF" Rscript -e "$r_expr" \
       && { [ "$DRYRUN" -eq 1 ] || [ -f "$target" ]; }; then
      d=$(fmt_dur $((SECONDS - t0)))
      if [ "$DRYRUN" -eq 1 ]; then
        results+=("$rel|$label|dry-run|-")
      else
        ok "$label" "$d"
        results+=("$rel|$label|ok|$d")
        n_done=$((n_done + 1))
      fi
    else
      failed "$label" "render failed or no output — see log"
      show_log_tail "$log"
      results+=("$rel|$label|FAIL|-")
      n_fail=$((n_fail + 1))
    fi
  fi

  echo
done

# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------
banner "Summary"
printf '%s%-46s %-30s %-16s %8s%s\n' "$B" "dataset" "artefact" "status" "time" "$R"
hr
for row in "${results[@]}"; do
  IFS='|' read -r r_rel r_step r_status r_dur <<< "$row"
  case "$r_status" in
    ok)     colour="$GRN" ;;
    FAIL*)  colour="$RED" ;;
    skip*)  colour="$YLW" ;;
    *)      colour="$DIM" ;;
  esac
  printf '%-46s %-30s %s%-16s%s %8s\n' \
    "$r_rel" "$r_step" "$colour" "$r_status" "$R" "$r_dur"
done
hr
printf '%s%d generated · %d skipped · %d failed%s   (total %s)\n' \
  "$B" "$n_done" "$n_skip" "$n_fail" "$R" "$(fmt_dur "$SECONDS")"
[ "$DRYRUN" -eq 1 ] || printf 'logs: %s\n' "$LOGDIR"

[ "$n_fail" -eq 0 ]