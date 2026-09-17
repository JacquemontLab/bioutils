#!/usr/bin/env bash
# Generate the ShortVariants column reports for every dataset under
#   $DB/<cohort>/shortvariants/<variant>/
# and place them in that dataset's docs/ directory:
#
#   1. ShortVariantsDB_curated_columns_report.pdf
#   2. ShortVariantsDB_unfiltered_columns_report.pdf
#
# Both come from ShortVariants-Annotation/bin/pdf_columns_report.py, which
# takes <parquet> <cpus> <mem_per_cpu_GB> and writes
# <basename>_columns_report.pdf next to its input.
#
# Steps run sequentially: each already uses $CPUS cores.
# Failures are recorded and the loop continues; exit status is non-zero if
# any step failed.
#
# Usage:
#   ./run_shortvariants_reports.sh [options]
#     --force            regenerate even if the docs/ output already exists
#     --dry-run          print what would run, execute nothing
#     --only PATTERN     restrict to datasets whose relative path matches
#                        PATTERN (bash glob, e.g. 'gtex/*' or '*wgs_dbgap*')
#     --verbose          stream tool output to the terminal as well as the log
#     --cpus N           override CPU count            (default 16)
#     --mem-per-cpu N    override memory per CPU in GB (default 4)
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
COLUMNS_REPORT="$BASE/SOFTWARE/Git_pipeline/ShortVariants-Annotation/bin/pdf_columns_report.py"

# Parquet basenames to report on, in order.
DBS=( ShortVariantsDB_curated ShortVariantsDB_unfiltered )

CPUS=16
MEM_PER_CPU=4
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
    --force)       FORCE=1; shift ;;
    --dry-run)     DRYRUN=1; shift ;;
    --verbose)     VERBOSE=1; shift ;;
    --only)        ONLY="${2:?--only needs a pattern}"; shift 2 ;;
    --cpus)        CPUS="${2:?--cpus needs a number}"; shift 2 ;;
    --mem-per-cpu) MEM_PER_CPU="${2:?--mem-per-cpu needs a number}"; shift 2 ;;
    --logdir)      LOGDIR="${2:?--logdir needs a path}"; shift 2 ;;
    -h|--help)     sed -n '2,28p' "$0"; exit 0 ;;
    *) printf 'unknown option: %s\n' "$1" >&2; exit 2 ;;
  esac
done

RUN_TS="$(date +%Y%m%d_%H%M%S)"
LOGDIR="${LOGDIR:-$PWD/logs_shortvariants_reports_$RUN_TS}"

# --------------------------------------------------------------------------
# Pretty printing
# --------------------------------------------------------------------------
if [ -t 1 ]; then
  B=$'\033[1m'; DIM=$'\033[2m'; RED=$'\033[31m'; GRN=$'\033[32m'
  YLW=$'\033[33m'; CYN=$'\033[36m'; R=$'\033[0m'
else
  B=''; DIM=''; RED=''; GRN=''; YLW=''; CYN=''; R=''
fi

hr()      { printf '%s\n' "${DIM}$(printf '─%.0s' {1..86})${R}"; }
banner()  { hr; printf '%s%s%s\n' "$B" "$1" "$R"; hr; }
info()    { printf '   %sℹ%s  %s\n' "$CYN" "$R" "$1"; }
ok()      { printf '   %s✅%s %-46s %s%s%s\n' "$GRN" "$R" "$1" "$DIM" "$2" "$R"; }
skipped() { printf '   %s⏭%s  %-46s %s%s%s\n' "$YLW" "$R" "$1" "$DIM" "$2" "$R"; }
failed()  { printf '   %s❌%s %-46s %s%s%s\n' "$RED" "$R" "$1" "$DIM" "$2" "$R"; }

fmt_dur() {  # seconds -> MMmSSs
  local s="$1"
  printf '%dm%02ds' $((s / 60)) $((s % 60))
}

# --------------------------------------------------------------------------
# Pre-flight
# --------------------------------------------------------------------------
banner "ShortVariants column reports — $RUN_TS"

preflight_fail=0
if [ -e "$COLUMNS_REPORT" ]; then
  info "found $(basename "$COLUMNS_REPORT") (ShortVariants-Annotation)"
else
  failed "missing prerequisite" "$COLUMNS_REPORT"
  preflight_fail=1
fi
[ -d "$DB" ] || { failed "missing DATABASE root" "$DB"; preflight_fail=1; }
[ "$preflight_fail" -eq 0 ] || exit 1

# Collect datasets: any <cohort>/shortvariants/<variant>/ holding at least
# one of the expected parquet files.
datasets=()
shopt -s nullglob
for d in "$DB"/*/shortvariants/*/; do
  d="${d%/}"
  rel="${d#"$DB"/}"
  has_any=0
  for base in "${DBS[@]}"; do
    # -e not -f: parquet may be a single file OR a partitioned directory
    [ -e "$d/$base.parquet" ] && { has_any=1; break; }
  done
  [ "$has_any" -eq 1 ] || continue
  # shellcheck disable=SC2254  # ONLY is intentionally a glob
  case "$rel" in $ONLY) datasets+=("$d") ;; esac
done
shopt -u nullglob

n_ds=${#datasets[@]}
if [ "$n_ds" -eq 0 ]; then
  failed "no datasets matched" "$DB/*/shortvariants/*/ (--only '$ONLY')"
  exit 1
fi

info "$n_ds dataset(s) to process · cpus=$CPUS · mem_per_cpu=${MEM_PER_CPU}G"
info "logs: $LOGDIR"
[ "$DRYRUN" -eq 1 ] && info "${YLW}dry-run: nothing will be executed${R}"
[ "$DRYRUN" -eq 1 ] || mkdir -p "$LOGDIR"
echo

# --------------------------------------------------------------------------
# Step runner: run_step <logfile> <cmd...>
# --------------------------------------------------------------------------
run_step() {
  local log="$1"; shift
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
results=()   # "rel|artefact|status|duration"
n_fail=0
n_done=0
n_skip=0
i=0

for ds in "${datasets[@]}"; do
  i=$((i + 1))
  rel="${ds#"$DB"/}"
  docs="$ds/docs"
  safe="${rel//\//_}"

  printf '%s[%d/%d]%s %s%s%s\n' "$B" "$i" "$n_ds" "$R" "$B" "$rel" "$R"
  [ "$DRYRUN" -eq 1 ] || mkdir -p "$docs"

  for base in "${DBS[@]}"; do
    parquet="$ds/$base.parquet"
    out="${base}_columns_report.pdf"
    target="$docs/$out"
    log="$LOGDIR/${safe}__${base}.log"

    if [ ! -e "$parquet" ]; then
      skipped "$out" "no $base.parquet"
      results+=("$rel|$out|skip(no input)|-")
      n_skip=$((n_skip + 1))
      continue
    fi
    if [ -f "$target" ] && [ "$FORCE" -eq 0 ]; then
      skipped "$out" "exists (use --force)"
      results+=("$rel|$out|skip(exists)|-")
      n_skip=$((n_skip + 1))
      continue
    fi

    t0=$SECONDS
    if run_step "$log" "$COLUMNS_REPORT" "$parquet" "$CPUS" "$MEM_PER_CPU"; then
      if [ "$DRYRUN" -eq 1 ]; then
        printf '   %s$%s mv %s %s\n' "$DIM" "$R" "$ds/$out" "$target"
        results+=("$rel|$out|dry-run|-")
      elif [ -f "$ds/$out" ] && mv -f "$ds/$out" "$target"; then
        d=$(fmt_dur $((SECONDS - t0)))
        ok "$out" "$d"
        results+=("$rel|$out|ok|$d")
        n_done=$((n_done + 1))
      else
        failed "$out" "script succeeded but $out not produced"
        show_log_tail "$log"
        results+=("$rel|$out|FAIL(no output)|-")
        n_fail=$((n_fail + 1))
      fi
    else
      failed "$out" "exit $? — see log"
      show_log_tail "$log"
      results+=("$rel|$out|FAIL|-")
      n_fail=$((n_fail + 1))
    fi
  done

  echo
done

# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------
banner "Summary"
printf '%s%-42s %-46s %-16s %8s%s\n' \
  "$B" "dataset" "artefact" "status" "time" "$R"
hr
for row in "${results[@]}"; do
  IFS='|' read -r r_rel r_step r_status r_dur <<< "$row"
  case "$r_status" in
    ok)     colour="$GRN" ;;
    FAIL*)  colour="$RED" ;;
    skip*)  colour="$YLW" ;;
    *)      colour="$DIM" ;;
  esac
  printf '%-42s %-46s %s%-16s%s %8s\n' \
    "$r_rel" "$r_step" "$colour" "$r_status" "$R" "$r_dur"
done
hr
printf '%s%d generated · %d skipped · %d failed%s   (total %s)\n' \
  "$B" "$n_done" "$n_skip" "$n_fail" "$R" "$(fmt_dur "$SECONDS")"
[ "$DRYRUN" -eq 1 ] || printf 'logs: %s\n' "$LOGDIR"

[ "$n_fail" -eq 0 ]