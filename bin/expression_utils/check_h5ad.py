#!/usr/bin/env python3
"""
check_h5ad.py — Quickly tell what state the data in an .h5ad file is in.

For a single-cell .h5ad, it answers two questions in two lines:
  1) Is .X raw counts, normalized, log-transformed, or scaled?
  2) What did .raw store (usually the untouched counts)?
And when .X looks log-transformed, it figures out WHICH log (log2 / ln / log10).

Why it's fast on huge files:
  It never loads the whole matrix. It reads only the first 500 cells (rows)
  straight from the HDF5 file with h5py. That's enough to characterize the data,
  and it stays instant even on files with millions of cells.

Usage:
  python check_h5ad.py myfile.h5ad
"""

import sys
import numpy as np
import h5py

# How many cells (rows) we sample to inspect the data. 500 is plenty to judge.
BLOCK = 500

# ---------------------------------------------------------------------------
# Terminal color helper (purely cosmetic).
# These are ANSI escape codes: g=green, y=yellow, c=cyan, d=dim/grey,
# b=bold, x=reset-back-to-normal.
# ---------------------------------------------------------------------------
C = dict(g="\033[32m", y="\033[33m", c="\033[36m", d="\033[90m", b="\033[1m", x="\033[0m")

def col(s, k):
    """Wrap string `s` in color `k` — but ONLY if we're printing to a real
    terminal. If output is piped to a file, return plain text so the file
    doesn't fill up with unreadable escape codes."""
    return f"{C[k]}{s}{C['x']}" if sys.stdout.isatty() else str(s)


# ---------------------------------------------------------------------------
# Read the first BLOCK rows (cells) of a matrix stored inside the .h5ad file,
# and return it as a dense NumPy array.
#
# A matrix can be stored dense or sparse (CSR/CSC). For the sparse cases we let
# anndata's `sparse_dataset` do the work: it wraps the on-disk sparse matrix so
# we can slice rows like a normal array, and it only reads the rows we ask for
# (it does NOT load the whole matrix). That saves us from hand-writing the
# fiddly index math for each sparse format.
# ---------------------------------------------------------------------------
try:                                    # import path moved between anndata versions
    from anndata.io import sparse_dataset
except ImportError:                     # older anndata
    from anndata.experimental import sparse_dataset


def _read_block(grp, n=BLOCK):
    """Read up to `n` rows of an anndata matrix as a dense NumPy array (default BLOCK)."""
    # Dense matrix: it's just an array on disk, slice the first n rows.
    if isinstance(grp, h5py.Dataset):
        return np.asarray(grp[:n])

    # Sparse (CSR or CSC): wrap it, slice the first n rows, densify.
    block = sparse_dataset(grp)[:n]
    return block.toarray()


# ---------------------------------------------------------------------------
# Look at a block of values and decide what kind of data it is.
# The logic is a series of checks, most-specific first.
# ---------------------------------------------------------------------------
def classify(sub):
    """Given a dense block `sub`, return an icon + label + color describing it."""
    # We only care about non-zero values. scRNA data is ~90% zeros, and zeros
    # tell us nothing about whether values are integer / decimal / negative.
    nz = sub[sub != 0]
    if nz.size == 0:
        return dict(icon="∅", label="empty", color="d", logged=False, vmin=None)

    # Are all non-zero values whole numbers? (allclose = tolerant of float noise.)
    # True  => looks like raw counts.
    # False => something transformed them into decimals.
    integer = np.allclose(nz, np.round(nz))

    vmin, vmax = float(nz.min()), float(nz.max())  # smallest / largest non-zero value

    # Per-cell total (sum across genes), then its "coefficient of variation"
    # = std / mean. If every cell was normalized to the same total (e.g. 10,000),
    # these sums are nearly identical, so CV ≈ 0. Raw counts vary a lot => CV big.
    sums = sub.sum(1)
    cv = sums.std() / sums.mean() if sums.mean() else np.nan

    # Heuristic: log1p squashes values hard. Even 20,000 counts becomes ~10 after
    # a log. So decimals with a max under ~30 are almost certainly log-transformed.
    logged = (not integer) and vmax < 30

    # Pick the verdict. Order matters: first match wins, most-specific at the top.
    if integer and vmin >= 0:      out = ("🔢", "raw counts", "g")            # whole numbers, no negatives
    elif vmin < 0:                 out = ("⚖️", "scaled (z-scored)", "y")     # negatives => scaled/centered
    elif logged and cv < 0.05:     out = ("📊", "normalized + log", "c")      # fixed total AND log scale
    elif logged:                   out = ("📉", "log-transformed", "c")       # log scale, totals still vary
    elif cv < 0.05:                out = ("📊", "normalized (linear)", "c")   # fixed total, but not logged
    else:                          out = ("🔸", "processed decimals", "y")    # decimals that fit no clean pattern

    return dict(icon=out[0], label=out[1], color=out[2], logged=logged, vmin=vmin)


# ---------------------------------------------------------------------------
# Read just the gene names from .var and .raw.var.
# We only need these when .X and .raw have different numbers of genes and we
# have to line them up. Gene names are stored in a fiddly "categorical" format
# on disk, so instead of parsing that by hand we let anndata read the (tiny)
# metadata. backed="r" means it does NOT load the big matrices — just the labels.
# ---------------------------------------------------------------------------
def _var_names(path):
    """Return (X gene names, raw gene names) — or (None, None) if unavailable."""
    try:
        import anndata as ad
        a = ad.read_h5ad(path, backed="r")            # loads metadata only, not .X
        xv = np.asarray(a.var_names)
        rv = np.asarray(a.raw.var_names) if a.raw is not None else None
        return xv, rv
    except Exception:
        return None, None


# ---------------------------------------------------------------------------
# Figure out the log base of .X (log2 vs ln vs log10) — BY MEASURING IT.
#
# We deliberately IGNORE any base recorded in .uns['log1p']: that record can be
# stale or wrong (e.g. someone re-logged the data and never updated it). The
# only thing we trust is the actual numbers, verified against .raw.
#
# How the measurement works:
#   A log-normalized .X is  X = log_base(shift + scale * raw_count), where:
#     - base  is 2, e, or 10          (which logarithm)
#     - shift is the pseudocount: 1 (log1p), 0.5, or 0 (plain log)
#     - scale is the unknown per-cell normalization factor (target_sum / lib_size)
#   We test the WHOLE GRID of (base × shift) candidates, not just the base.
#
#   For a given (base, shift): invert the log to recover  base**X - shift, which
#   should equal  scale * raw  — i.e. proportional to raw, a straight line through
#   the origin. The unknown per-cell `scale` is just that line's slope, so we don't
#   need target_sum. The candidate whose (base**X - shift) best matches raw wins.
#   Getting the shift WRONG bends that line (especially at low counts), dropping R²,
#   which is exactly how we tell log2(x) from log2(1+x) from log2(x+0.5).
# ---------------------------------------------------------------------------
BASES = {"log2": 2.0, "ln": np.e, "log10": 10.0}
SHIFTS = {"1p": 1.0, "0.5": 0.5, "0": 0.0}       # +1 (log1p), +0.5, or no shift


def _fit_score(X, R):
    """For each (base, shift) candidate, how well does base**X - shift ∝ raw?
    Returns (best_label, best_R2) using the median R² across sampled cells."""
    results = {}
    for bname, b in BASES.items():
        for sname, shift in SHIFTS.items():
            r2s = []
            for i in range(X.shape[0]):
                x, r = X[i], R[i]
                m = (x > 0) & (r > 0)              # genes expressed in both
                if m.sum() < 20:
                    continue
                # invert the candidate log: recovered should be proportional to raw
                recovered = np.power(b, x[m]) - shift     # = base**X - shift
                rr = r[m].astype(float)
                # R² of a straight line through the origin: recovered ≈ slope * raw
                slope = np.dot(recovered, rr) / np.dot(rr, rr)
                resid = recovered - slope * rr
                sst = np.sum((recovered - recovered.mean()) ** 2)
                if sst <= 0:
                    continue
                r2s.append(1 - np.sum(resid ** 2) / sst)
            if r2s:
                results[f"{bname}({sname})"] = float(np.median(r2s))
    if not results:
        return None, None
    best = max(results, key=results.get)
    return best, results[best]


def log_base(f, path):
    """Return (transform_name, status_tag), e.g. ('log2(1p)', '✔ verified'), or (None, reason)."""

    # No .raw counts => nothing to verify against. We refuse to guess.
    if "raw" not in f:
        return None, "no .raw to verify against"

    # A handful of cells is enough to identify the transform. Small read.
    READ_CELLS = 40
    X = _read_block(f["X"], READ_CELLS)           # first few cells of .X  (logged values)
    R = _read_block(f["raw"]["X"], READ_CELLS)    # first few cells of .raw (the counts)

    # If .X and .raw have different gene counts, align them by gene name so we're
    # comparing the same gene in both. (.raw usually has MORE genes than .X.)
    if X.shape[1] != R.shape[1]:
        xv, rv = _var_names(path)
        if xv is None or rv is None:
            return None, "can't align genes (.X vs .raw differ)"
        pos = {g: i for i, g in enumerate(rv)}                 # gene name -> column in raw
        keep = [(i, pos[g]) for i, g in enumerate(xv) if g in pos]  # genes present in both
        if len(keep) < 50:                                     # too few to trust the fit
            return None, "too few shared genes"
        xi = [k[0] for k in keep]; ri = [k[1] for k in keep]
        X, R = X[:, xi], R[:, ri]

    # .raw must actually be counts for this to mean anything. If .raw is itself
    # decimals/logged, the "log of raw" test is meaningless — bail honestly.
    rnz = R[R != 0]
    if rnz.size and not np.allclose(rnz, np.round(rnz)):
        return None, ".raw isn't counts, can't verify"

    # Test the full base × shift grid; take the best-fitting transform.
    name, r2 = _fit_score(X, R)
    if name is None:
        return None, "no usable cells to fit"

    # A good fit (high R²) means .X really is this exact transform of .raw.
    if r2 > 0.98:
        return name, "✔ verified using .raw"
    return None, f"⚠ .X is not a clean log of .raw (R²={r2:.2f})"


# ---------------------------------------------------------------------------
# Main: open the file, inspect .X, inspect .raw, print two tidy lines.
# ---------------------------------------------------------------------------
def main(path):
    with h5py.File(path, "r") as f:
        # Get the matrix dimensions for the header. Sparse matrices keep the
        # shape in an attribute; dense ones expose it directly.
        nrow = tuple(f["X"].attrs["shape"])[0] if "shape" in f["X"].attrs else f["X"].shape[0]
        ncol = tuple(f["X"].attrs["shape"])[1] if "shape" in f["X"].attrs else f["X"].shape[1]
        print(f"\n{col('▸ ' + path, 'b')}   {col(f'{nrow} cells × {ncol} genes', 'd')}")

        # --- .X line ---
        xs = classify(_read_block(f["X"]))
        line = f"  {xs['icon']}  {col('.X', 'b')}   {col(xs['label'], xs['color'])}"
        # Only bother computing the log base if .X actually looks logged.
        # We always verify against .raw (never trust a recorded base), so the
        # result is either "(log2, ✔ verified)" or the reason we couldn't verify.
        if xs["logged"]:
            base, tag = log_base(f, path)
            shown = f"{base}, {tag}" if base else tag   # tag holds the reason on failure
            line += col(f"   ({shown})", "d")
        print(line)

        # --- .raw line ---
        if "raw" not in f:
            print(f"  {col('·', 'd')}  {col('.raw', 'b')}  {col('absent', 'd')}")
        else:
            rs = classify(_read_block(f["raw"]["X"]))
            print(f"  {rs['icon']}  {col('.raw', 'b')} {col(rs['label'], rs['color'])}")
        print()


if __name__ == "__main__":
    # Expect exactly one argument: the path to the .h5ad file.
    if len(sys.argv) != 2:
        sys.exit("usage: python check_h5ad.py file.h5ad")
    main(sys.argv[1])