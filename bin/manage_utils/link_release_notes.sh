#!/usr/bin/env bash
# link_release_notes.sh -- symlink biodoc release notes into the DATABASE docs/ dirs
set -euo pipefail


if [[ "${CC_CLUSTER:-}" == "rorqual" ]]; then
    BIND=/lustre09/project/6008022
elif [[ "$(hostname -s)" == "chusj-transcriptomic-server-1" ]]; then
    BIND=/mnt/chusj-transcriptomic-cephfs-1
else
    echo "Unknown site: $(hostname -s)" >&2
    return 1 2>/dev/null || exit 1
fi
export BIND


SRC_ROOT="$BIND/LAB_WORKSPACE/SOFTWARE/biodoc/DATABASE_SD4H"
PATTERN='*_release_note.pdf'
DRYRUN=0
RELATIVE=0
PRUNE=0

usage() {
    cat >&2 <<EOF
Usage: ${0##*/} [-n] [-r] [-R] [-p PATTERN] [-s BIODOC_ROOT] <database_root>

  -n, --dry-run        show what would be linked, change nothing
  -r, --relative       make relative symlinks (ln -rs) instead of absolute
  -R, --prune          delete pre-existing regular-file release notes in docs/
                       whose name differs from the link being created
  -p, --pattern GLOB   file glob inside each biodoc leaf (default: $PATTERN)
  -s, --source DIR     biodoc root (default: $SRC_ROOT)
  -h, --help           this message
EOF
    exit "${1:-2}"
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -n|--dry-run)  DRYRUN=1; shift ;;
        -r|--relative) RELATIVE=1; shift ;;
        -R|--prune)    PRUNE=1; shift ;;
        -p|--pattern)  PATTERN="${2:?-p needs a glob}"; shift 2 ;;
        -s|--source)   SRC_ROOT="${2:?-s needs a directory}"; shift 2 ;;
        -h|--help)     usage 0 ;;
        --)            shift; break ;;
        -*)            printf 'unknown option: %s\n' "$1" >&2; usage ;;
        *)             break ;;
    esac
done

[[ $# -eq 1 ]] || usage
DST_ROOT="${1%/}"
SRC_ROOT="${SRC_ROOT%/}"
[[ -d $SRC_ROOT ]] || { printf 'not a directory: %s\n' "$SRC_ROOT" >&2; exit 1; }
[[ -d $DST_ROOT ]] || { printf 'not a directory: %s\n' "$DST_ROOT" >&2; exit 1; }

(( DRYRUN )) && printf '### DRY RUN -- nothing will be written ###\n\n' >&2

n_linked=0 n_skipped=0 n_pruned=0
verb=$( (( DRYRUN )) && echo 'WOULD' || echo 'LINK ' )

# biodoc leaves are exactly cohort/variant_type/source
while IFS= read -r -d '' leaf; do
    rel="${leaf#"$SRC_ROOT"/}"
    IFS=/ read -r cohort vtype src <<< "$rel"
    [[ $vtype == cnv || $vtype == shortvariants ]] || continue

    docs="$DST_ROOT/$rel/docs"
    if [[ ! -d $docs ]]; then
        printf 'SKIP  no DATABASE docs/ for %s\n' "$rel" >&2
        n_skipped=$((n_skipped+1)); continue
    fi

    found=0
    while IFS= read -r -d '' note; do
        found=1
        base="${note##*/}"

        # drop stale, differently-named release notes that are real files
        if (( PRUNE )); then
            while IFS= read -r -d '' old; do
                [[ ${old##*/} == "$base" ]] && continue
                printf 'PRUNE %s/docs/%s\n' "$rel" "${old##*/}"
                (( DRYRUN )) || rm -f -- "$old"
                n_pruned=$((n_pruned+1))
            done < <(find "$docs" -maxdepth 1 -type f -name '*release_note*.pdf' -print0)
        fi

        printf '%s %s/docs/%s -> %s\n' "$verb" "$rel" "$base" "${note#"$SRC_ROOT"/}"
        if (( ! DRYRUN )); then
            if (( RELATIVE )); then
                ln -sfn --relative -- "$note" "$docs/$base"
            else
                ln -sfn -- "$note" "$docs/$base"
            fi
        fi
        n_linked=$((n_linked+1))
    done < <(find "$leaf" -maxdepth 1 -type f -name "$PATTERN" -print0)

    (( found )) || printf 'WARN  no %s in %s\n' "$PATTERN" "$rel" >&2
done < <(find "$SRC_ROOT" -mindepth 3 -maxdepth 3 -type d -print0)

printf '\n%d link(s) %s, %d pruned, %d dir(s) skipped\n' \
       "$n_linked" "$( (( DRYRUN )) && echo 'would be created' || echo created )" \
       "$n_pruned" "$n_skipped" >&2