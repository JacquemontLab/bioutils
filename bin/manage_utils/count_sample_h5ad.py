#!/usr/bin/env python3
"""Count unique sample identifiers and summarise annotation columns in .h5ad files.

Modes
-----
  count_sample_h5ad.py FILE...         -> |union of sample IDs|, one number
  count_sample_h5ad.py --batch DIR...  -> per dir, one TSV line (always 9 fields):
        DIR NFILES NSAMPLES N_TISSUE TISSUE_VALS N_REGION REGION_VALS N_SUBREGION SUBREGION_VALS
     plus a final tagged line with exact tree-wide distinct counts:
        #TOTAL<TAB>N_TISSUE<TAB>N_REGION<TAB>N_SUBREGION

Only /obs is touched: no anndata import, no access to X, var or unrelated obs
columns. For a categorical column, the level names are read from `categories`
(a few kB) and, by default, filtered against `codes` so that levels retained by
pandas after subsetting are not counted. --fast skips the `codes` read.

A present-but-blank column (single '' level) is reported as NA, not as 0.
Prints NA when no candidate column exists for a given field.
"""

from __future__ import annotations

import argparse
import os
import sys

import h5py
import numpy as np

# Checked in this order; first hit wins, for both the obs index and obs columns.
SAMPLE_COLUMNS: tuple[str, ...] = (
    "SampleID", "sampleID",
)

# Extra annotation columns to summarise in --batch mode. Each entry is a tuple
# of accepted aliases, tried in order; first hit wins. Order here drives the
# TSV field order, the #TOTAL order, and (mirrored in Bash) the print order.
ANNOTATION_COLUMNS: dict[str, tuple[str, ...]] = {
    "Tissue":    ("Tissue", "tissue"),
    "Region":    ("Region", "region"),
    "Subregion": ("Subregion", "subregion", "SubRegion"),
}

# Cap on distinct values inlined per column in the TSV, to keep lines readable.
MAX_VALS_INLINE = 5


def _decode(values) -> list[str]:
    """Normalise fixed-length bytes / object strings to str, dropping blanks."""
    out = []
    for v in values:
        s = v.decode() if isinstance(v, (bytes, np.bytes_)) else str(v)
        if s:                      # HDF5 stores missing strings as ""
            out.append(s)
    return out


def _values(node, parent: h5py.Group, name: str, exact: bool) -> set[str]:
    """Extract the observed values of one obs column, whatever its encoding."""
    # Modern categorical: group holding 'categories' + 'codes'.
    if isinstance(node, h5py.Group):
        cats = np.asarray(node["categories"][:])
        if not exact:
            return set(_decode(cats))
        codes = node["codes"][:]
        return set(_decode(cats[np.unique(codes[codes >= 0])]))

    data = node[:]

    # Legacy categorical: integer codes + /uns/<name>_categories.
    if np.issubdtype(data.dtype, np.integer):
        cat_path = f"uns/{name}_categories"
        root = parent.file
        if cat_path in root:
            cats = np.asarray(root[cat_path][:])
            return set(_decode(cats[np.unique(data[data >= 0])]))
        return {str(v) for v in np.unique(data)}   # genuine integer IDs

    return set(_decode(data))


def _column_values(
    obs, index_name: str | None, aliases: tuple[str, ...], exact: bool
) -> set[str] | None:
    """First matching column among `aliases`, resolved against a group `obs`."""
    for cand in aliases:
        if index_name == cand:
            return set(_decode(obs[index_name][:]))
        if cand in obs:
            return _values(obs[cand], obs, cand, exact)
    return None


def extract(path: str, exact: bool = True) -> dict[str, set[str] | None]:
    """Return {'SampleID': set|None, 'Tissue': set|None, ...} for one file.

    Keys are the logical names in ANNOTATION_COLUMNS plus 'SampleID'.
    A value of None means the corresponding column was absent.
    """
    result: dict[str, set[str] | None] = {"SampleID": None}
    for logical in ANNOTATION_COLUMNS:
        result[logical] = None

    with h5py.File(path, "r") as f:
        if "obs" not in f:
            raise KeyError("no /obs group")
        obs = f["obs"]

        # Very old layout: obs is a single compound dataset.
        if isinstance(obs, h5py.Dataset):
            fields = obs.dtype.names or ()
            for cand in SAMPLE_COLUMNS:
                if cand in fields:
                    result["SampleID"] = set(_decode(obs[cand]))
                    break
            for logical, aliases in ANNOTATION_COLUMNS.items():
                for cand in aliases:
                    if cand in fields:
                        result[logical] = set(_decode(obs[cand]))
                        break
            return result

        index_name = obs.attrs.get("_index")
        if isinstance(index_name, bytes):
            index_name = index_name.decode()

        result["SampleID"] = _column_values(obs, index_name, SAMPLE_COLUMNS, exact)
        for logical, aliases in ANNOTATION_COLUMNS.items():
            result[logical] = _column_values(obs, index_name, aliases, exact)

    return result


def sample_ids(path: str, exact: bool = True) -> set[str] | None:
    """Backward-compatible helper: sample IDs only."""
    return extract(path, exact)["SampleID"]


def union(paths, exact: bool) -> tuple[dict[str, set[str]], set[str]]:
    """Union of each field's values over `paths`.

    Returns (unions, found_fields) where `found_fields` is the set of logical
    names for which at least one file supplied the column (even if blank).
    """
    fields = ("SampleID", *ANNOTATION_COLUMNS)
    unions: dict[str, set[str]] = {k: set() for k in fields}
    found: set[str] = set()
    for path in paths:
        try:
            got = extract(path, exact)
        except (OSError, KeyError, IndexError) as exc:
            print(f"ERROR: {path}: {exc}", file=sys.stderr)
            continue
        for k in fields:
            if got[k] is not None:
                unions[k] |= got[k]
                found.add(k)
    return unions, found


def h5ads(directory: str) -> list[str]:
    """Non-recursive listing of .h5ad files, one syscall pass."""
    try:
        with os.scandir(directory) as it:
            return [e.path for e in it if e.name.endswith(".h5ad") and e.is_file()]
    except OSError as exc:
        print(f"ERROR: {directory}: {exc}", file=sys.stderr)
        return []


def _fmt_vals(vals: set[str]) -> str:
    """Sorted, comma-joined, tab-free value list with a +N overflow marker."""
    if not vals:
        return ""
    ordered = sorted(v.replace("\t", " ") for v in vals)
    shown = ordered[:MAX_VALS_INLINE]
    text = ", ".join(shown)
    if len(ordered) > MAX_VALS_INLINE:
        text += f", …(+{len(ordered) - MAX_VALS_INLINE})"
    return text


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("paths", nargs="+", metavar="FILE|DIR")
    p.add_argument("--batch", action="store_true",
                   help="treat arguments as directories; emit one TSV line each")
    p.add_argument("--fast", action="store_true",
                   help="trust categorical levels without reading codes")
    p.add_argument("--list", action="store_true",
                   help="print the sample IDs instead of the count (non-batch only)")
    args = p.parse_args()
    exact = not args.fast

    if args.batch:
        # Tree-wide distinct sets, accumulated exactly (not from truncated strings).
        totals: dict[str, set[str]] = {k: set() for k in ANNOTATION_COLUMNS}

        for directory in args.paths:
            files = h5ads(directory)
            if not files:
                continue
            unions, found = union(files, exact)
            n_samp = str(len(unions["SampleID"])) if "SampleID" in found else "NA"
            fields = [directory, str(len(files)), n_samp]
            for logical in ANNOTATION_COLUMNS:              # Tissue, Region, Subregion
                if logical in found and unions[logical]:    # present AND non-blank
                    fields.append(str(len(unions[logical])))
                    fields.append(_fmt_vals(unions[logical]))
                    totals[logical] |= unions[logical]
                else:                                        # absent or all-blank -> NA
                    fields.append("NA")
                    fields.append("")
            print("\t".join(fields))                         # always exactly 9 fields

        # Exact global distinct counts on stdout, tagged for the caller to filter.
        print("#TOTAL\t" + "\t".join(str(len(totals[k])) for k in ANNOTATION_COLUMNS))
        return 0

    unions, found = union(args.paths, exact)
    ids = unions["SampleID"]
    if args.list:
        print("\n".join(sorted(ids)))
    elif "SampleID" in found:
        print(len(ids))
    else:
        print("NA")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())