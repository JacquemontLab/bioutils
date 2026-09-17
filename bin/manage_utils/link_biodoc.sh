#!/usr/bin/env bash
# sync_reports.sh -- copy docs/*_report.pdf from the DATABASE tree into biodoc
set -euo pipefail

if [[ "${CC_CLUSTER:-}" == "rorqual" ]]; then
    DST_ROOT=/lustre09/project/6008022/LAB_WORKSPACE/SOFTWARE/biodoc/DATABASE_CCDB
elif [[ "$(hostname -s)" == "chusj-transcriptomic-server-1" ]]; then
    DST_ROOT=/mnt/chusj-transcriptomic-cephfs-1/LAB_WORKSPACE/SOFTWARE/biodoc/DATABASE_SD4H
else
    DST_ROOT=''
fi

PATTERN='*_report.pdf|cnv_dataset_qc.pdf'
DRYRUN=0

usage() {
    cat >&2 <<EOF
Usage: ${0##*/} [-n] [-p PATTERN] [-d BIODOC_ROOT] <database_root>

  -n, --dry-run        show what would be copied, copy nothing
  -p, --pattern GLOB   file glob(s) inside docs/, '|'-separated
                       (default: $PATTERN)
  -d, --dest DIR       biodoc root (default: ${DST_ROOT:-<none: required on this host>})
  -h, --help           this message
EOF
    exit "${1:-2}"
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -n|--dry-run) DRYRUN=1; shift ;;
        -p|--pattern) PATTERN="${2:?-p needs a glob}"; shift 2 ;;
        -d|--dest)    DST_ROOT="${2:?-d needs a directory}"; shift 2 ;;
        -h|--help)    usage 0 ;;
        --)           shift; break ;;
        -*)           printf 'unknown option: %s\n' "$1" >&2; usage ;;
        *)            break ;;
    esac
done

[[ $# -eq 1 ]] || usage
SRC_ROOT="${1%/}"
DST_ROOT="${DST_ROOT%/}"
[[ -d $SRC_ROOT ]] || { printf 'not a directory: %s\n' "$SRC_ROOT" >&2; exit 1; }
[[ -n $DST_ROOT ]] || { printf 'unknown site %s: pass -d BIODOC_ROOT\n' "$(hostname -s)" >&2; exit 1; }
[[ -d $DST_ROOT ]] || { printf 'not a directory: %s\n' "$DST_ROOT" >&2; exit 1; }

# ---- expand PATTERN ('a|b|c') into a find -name test ------------------------
find_expr=()
IFS='|' read -r -a globs <<< "$PATTERN"
for g in "${globs[@]}"; do
    [[ -n $g ]] || continue
    (( ${#find_expr[@]} )) && find_expr+=(-o)
    find_expr+=(-name "$g")
done
(( ${#find_expr[@]} )) || { printf 'empty pattern: %s\n' "$PATTERN" >&2; exit 1; }

(( DRYRUN )) && printf '### DRY RUN -- nothing will be written ###\n\n' >&2

n_copied=0 n_skipped=0

# docs/ at exactly cohort/variant_type/source/docs -> avoids phenotypic/docs etc.
while IFS= read -r -d '' docs; do
    rel="${docs#"$SRC_ROOT"/}"; rel="${rel%/docs}"
    IFS=/ read -r cohort vtype src <<< "$rel"
    [[ $vtype == cnv || $vtype == shortvariants ]] || continue

    dst="$DST_ROOT/$rel"
    if [[ ! -d $dst ]]; then
        printf 'SKIP  no biodoc dir for %s\n' "$rel" >&2
        n_skipped=$((n_skipped+1)); continue
    fi

    found=0
    while IFS= read -r -d '' pdf; do
        found=1
        printf '%s  %s -> %s/\n' "$( (( DRYRUN )) && echo 'WOULD' || echo 'COPY ' )" \
               "${pdf#"$SRC_ROOT"/}" "$rel"
        (( DRYRUN )) || cp -pL --no-preserve=ownership -- "$pdf" "$dst/"
        n_copied=$((n_copied+1))
    done < <(find "$docs" -maxdepth 1 -type f \( "${find_expr[@]}" \) -print0)

    (( found )) || printf 'WARN  no %s in %s/docs\n' "$PATTERN" "$rel" >&2
done < <(find "$SRC_ROOT" -mindepth 4 -maxdepth 4 -type d -name docs -print0)

printf '\n%d file(s) %s, %d dir(s) skipped\n' \
       "$n_copied" "$( (( DRYRUN )) && echo 'would be copied' || echo copied )" "$n_skipped" >&2