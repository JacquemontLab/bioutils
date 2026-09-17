#!/usr/bin/env bash
# Check that each cohort has README.md + metadata_sampleDB.tsv,
# and that each dataset's docs/ contains the required files.
#
# Expected layout: <cohort>/{README.md,metadata_sampleDB.tsv}
#                  <cohort>/<cnv|shortvariants>/<variant>/docs/
#
# README.md inside docs/ is NOT required; if present it's flagged
# (stray copy — you may want to remove it).
#
# Usage: ./check_docs.sh [ROOT]   (defaults to current dir)
# Exit status: 0 if nothing is missing, 1 otherwise.
set -uo pipefail

ROOT="${1:-.}"; ROOT="${ROOT%/}"   # trailing slash would break the rel/depth logic
missing_total=0
dirs_checked=0
dirs_with_gaps=0
cohorts_checked=0
cohorts_with_gaps=0

# ANSI bold (disabled if output is not a terminal)
if [ -t 1 ]; then B=$'\033[1m'; R=$'\033[0m'; else B=''; R=''; fi

# Files required at the cohort level.
cohort_required=( "README.md" "metadata_sampleDB.tsv" )

# Subdirectories whose presence identifies a top-level dir as a cohort.
cohort_markers=( cnv shortvariants variants expression phenotypic )

required_items() {
  case "$1" in
    cnv)
      printf '%s\n' \
        "exact:cnvDB_columns_report.pdf" \
        "exact:geneDB_columns_report.pdf" \
        "exact:cnv_dataset_qc.pdf" \
        "glob:*_release_note.pdf" \
        "exact:launch_report.txt"
      ;;
    shortvariants)
      printf '%s\n' \
        "exact:ShortVariantsDB_curated_columns_report.pdf" \
        "exact:ShortVariantsDB_unfiltered_columns_report.pdf" \
        "glob:*_release_note.pdf" \
        "exact:launch_report.txt"
      ;;
    variants)
      printf '%s\n' \
        "glob:*_release_note.pdf"
      ;;
  esac
}

# Build a display path with the first component (cohort) bold.
# e.g. cmc/cnv/array_called/docs -> <B>cmc<R>/cnv/array_called/docs
bold_path() {
  local rel="$1" head tail
  head="${rel%%/*}"
  tail="${rel#*/}"
  if [ "$head" = "$rel" ]; then
    printf '%s%s%s' "$B" "$rel" "$R"
  else
    printf '%s%s%s/%s' "$B" "$head" "$R" "$tail"
  fi
}

# ---------------------------------------------------------------------------
# Pass 1 — cohort-level files: <cohort>/{README.md,metadata_sampleDB.tsv}
# ---------------------------------------------------------------------------
while IFS= read -r cohortdir; do
  rel="${cohortdir#"$ROOT"/}"

  # A cohort is a top-level dir holding at least one known datatype dir;
  # this skips stray dirs (scripts/, logs/, ...) without hardcoding names.
  is_cohort=0
  for d in "${cohort_markers[@]}"; do
    [ -d "$cohortdir/$d" ] && { is_cohort=1; break; }
  done
  [ "$is_cohort" -eq 1 ] || continue

  cohorts_checked=$((cohorts_checked+1))
  missing_here=()
  for f in "${cohort_required[@]}"; do
    [ -f "$cohortdir/$f" ] || missing_here+=("$f")
  done

  disp="$(bold_path "$rel")"

  if [ ${#missing_here[@]} -eq 0 ]; then
    printf '\342\234\205 \360\237\223\246 %s\n' "$disp"
  else
    cohorts_with_gaps=$((cohorts_with_gaps+1))
    missing_total=$((missing_total+${#missing_here[@]}))
    printf '\342\235\214 \360\237\223\246 %s  (%d missing)\n' "$disp" "${#missing_here[@]}"
    for m in "${missing_here[@]}"; do
      printf '     - %s\n' "$m"
    done
  fi
done < <(find "$ROOT" -mindepth 1 -maxdepth 1 -type d ! -name '.*' | sort)

echo

# ---------------------------------------------------------------------------
# Pass 2 — dataset docs dirs at exactly <cohort>/<datatype>/<variant>/docs
# -maxdepth 4 keeps find from descending into nested docs
# (e.g. .../docs/cnv_caller/docs).
# ---------------------------------------------------------------------------
while IFS= read -r docdir; do
  rel="${docdir#"$ROOT"/}"

  IFS='/' read -r -a parts <<< "$rel"
  [ "${#parts[@]}" -eq 4 ] || continue
  [ "${parts[3]}" = "docs" ] || continue
  cohort="${parts[0]}"
  datatype="${parts[1]}"
  case "$datatype" in cnv|shortvariants|variants) ;; *) continue ;; esac

  dirs_checked=$((dirs_checked+1))
  missing_here=()
  release_notes=()          # actual release_note file(s) found (names vary)

  while IFS= read -r item; do
    [ -z "$item" ] && continue
    kind="${item%%:*}"; pat="${item#*:}"
    case "$kind" in
      exact) [ -f "$docdir/$pat" ] || missing_here+=("$pat") ;;
      glob)
        shopt -s nullglob
        matches=( "$docdir"/$pat )
        shopt -u nullglob
        [ ${#matches[@]} -ge 1 ] || missing_here+=("$pat")
        # keep the names so they can be reported (see below)
        case "$pat" in *_release_note.pdf) release_notes=( "${matches[@]}" ) ;; esac
        ;;
    esac
  done < <(required_items "$datatype")

  # README inside docs/ is optional: note it if present (candidate for removal)
  has_readme=0
  [ -f "$docdir/README.md" ] && has_readme=1

  disp="$(bold_path "$rel")"

  if [ ${#missing_here[@]} -eq 0 ]; then
    printf '\342\234\205 \360\237\223\201 %s\n' "$disp"
  else
    dirs_with_gaps=$((dirs_with_gaps+1))
    missing_total=$((missing_total+${#missing_here[@]}))
    printf '\342\235\214 \360\237\223\201 %s  (%d missing)\n' "$disp" "${#missing_here[@]}"
    for m in "${missing_here[@]}"; do
      printf '     - %s\n' "$m"
    done
  fi

  # Report the release note(s) actually present: the names are not
  # standardised, and they may be symlinks into biodoc.
  for n in "${release_notes[@]}"; do
    if [ -L "$n" ]; then
      printf '     \360\237\223\204 %s \342\206\222 %s\n' "${n##*/}" 
    else
      printf '     \360\237\223\204 %s\n' "${n##*/}"
    fi
  done
  [ ${#release_notes[@]} -gt 1 ] && \
    printf '     \342\232\240\357\270\217  %d release notes in one docs/\n' "${#release_notes[@]}"

  [ "$has_readme" -eq 1 ] && printf '     \342\204\271\357\270\217  README.md present (optional \342\200\224 remove?)\n'
done < <(
  find "$ROOT" -maxdepth 4 -type d -name docs \
       \( -path '*/cnv/*' -o -path '*/shortvariants/*' -o -path '*/variants/*' \) | sort
)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo
printf 'checked %d cohorts, %d with gaps\n' \
  "$cohorts_checked" "$cohorts_with_gaps"
printf 'checked %d docs dirs, %d with gaps\n' \
  "$dirs_checked" "$dirs_with_gaps"
printf '%d files missing total\n' "$missing_total"

[ "$missing_total" -eq 0 ]