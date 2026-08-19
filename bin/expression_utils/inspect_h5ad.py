#!/usr/bin/env python3
"""
Inspect an AnnData .h5ad file with clear, color-coded output.

Usage:
    ./inspect_h5ad.py file.h5ad
    ./inspect_h5ad.py file.h5ad --max-cats 20 --no-x-stats
    ./inspect_h5ad.py file.h5ad --no-color        # or set NO_COLOR=1
"""
import argparse
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd


# --------------------------------------------------------------------------- #
# Color handling
# --------------------------------------------------------------------------- #
class C:
    """ANSI codes; blanked out when color is disabled."""
    RESET = "\033[0m"
    BOLD = "\033[1m"
    DIM = "\033[2m"
    RED = "\033[31m"
    GREEN = "\033[32m"
    YELLOW = "\033[33m"
    BLUE = "\033[34m"
    MAGENTA = "\033[35m"
    CYAN = "\033[36m"
    GREY = "\033[90m"

    @classmethod
    def disable(cls):
        for k, v in vars(cls).items():
            if isinstance(v, str) and v.startswith("\033"):
                setattr(cls, k, "")


def color_enabled(force_no_color: bool) -> bool:
    """Enable color only for a real TTY, unless disabled by flag/env."""
    if force_no_color or os.environ.get("NO_COLOR"):
        return False
    return sys.stdout.isatty()


# --------------------------------------------------------------------------- #
# Pretty-printing helpers
# --------------------------------------------------------------------------- #
def header(title: str) -> None:
    bar = "─" * (len(title) + 2)
    print(f"\n{C.BOLD}{C.CYAN}┌{bar}┐{C.RESET}")
    print(f"{C.BOLD}{C.CYAN}│ {title} │{C.RESET}")
    print(f"{C.BOLD}{C.CYAN}└{bar}┘{C.RESET}")


def kv(key: str, value, key_color: str = C.YELLOW) -> None:
    print(f"  {key_color}{key:<14}{C.RESET} {value}")


def human_bytes(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(n) < 1024.0:
            return f"{n:3.1f} {unit}"
        n /= 1024.0
    return f"{n:.1f} PB"


def fmt_list(values, color: str) -> str:
    return ", ".join(f"{color}{v!r}{C.RESET}" for v in values)


# --------------------------------------------------------------------------- #
# Section printers
# --------------------------------------------------------------------------- #
def print_overview(adata, input_file: str) -> None:
    header("OVERVIEW")
    try:
        size = os.path.getsize(input_file)
        kv("File", f"{input_file}  {C.DIM}({human_bytes(size)}){C.RESET}")
    except OSError:
        kv("File", input_file)
    n_obs, n_vars = adata.shape
    kv("Shape", f"{C.GREEN}{n_obs:,}{C.RESET} obs × "
                f"{C.GREEN}{n_vars:,}{C.RESET} vars")
    kv("X type", type(adata.X).__name__)
    try:
        kv("X dtype", adata.X.dtype)
    except AttributeError:
        pass


def print_obs(adata, max_cats: int) -> None:
    header(f"obs  —  cell metadata ({adata.obs.shape[1]} columns)")
    print(f"{C.DIM}First rows:{C.RESET}")
    print(adata.obs.head().to_string())
    print(f"\n{C.DIM}Per-column summary (showing up to {max_cats} values):{C.RESET}\n")

    for col in adata.obs.columns:
        vals = adata.obs[col]
        if isinstance(vals.dtype, pd.CategoricalDtype):
            cats = list(vals.cat.categories)
            more = f" {C.DIM}(+{len(cats) - max_cats} more){C.RESET}" if len(cats) > max_cats else ""
            tag = f"{C.MAGENTA}categorical{C.RESET}"
            print(f"  {C.BOLD}{col}{C.RESET} [{tag}, {C.GREEN}{len(cats)}{C.RESET} levels]: "
                  f"{fmt_list(cats[:max_cats], C.CYAN)}{more}")
        else:
            u = vals.unique()
            tag = f"{C.BLUE}{vals.dtype}{C.RESET}"
            more = f" {C.DIM}...{C.RESET}" if len(u) > max_cats else ""
            print(f"  {C.BOLD}{col}{C.RESET} [{tag}, n_unique={C.GREEN}{len(u):,}{C.RESET}]: "
                  f"{fmt_list(list(u[:max_cats]), C.CYAN)}{more}")


def print_var(adata) -> None:
    header(f"var  —  gene metadata ({adata.var.shape[1]} columns)")
    print(adata.var.head().to_string())


def print_slots(adata) -> None:
    header("SLOTS")
    for name, keys in (
        ("layers", list(adata.layers.keys())),
        ("obsm", list(adata.obsm.keys())),
        ("varm", list(adata.varm.keys())),
        ("obsp", list(adata.obsp.keys())),
        ("varp", list(adata.varp.keys())),
        ("uns", list(adata.uns.keys())),
    ):
        shown = fmt_list(keys, C.CYAN) if keys else f"{C.DIM}(empty){C.RESET}"
        kv(name, shown)


def print_x_stats(adata) -> None:
    header("X statistics")
    print(f"{C.DIM}Loading X into memory...{C.RESET}")
    X = adata.X[:]
    kv("type", type(X).__name__)
    kv("dtype", X.dtype)

    xmin, xmax = X.min(), X.max()
    kv("min / max", f"{C.GREEN}{xmin}{C.RESET} / {C.GREEN}{xmax}{C.RESET}")

    if hasattr(X, "nnz"):
        nnz = X.nnz
    else:
        nnz = int(np.count_nonzero(X))
    total = adata.shape[0] * adata.shape[1]
    sparsity = 100.0 * (1 - nnz / total) if total else 0.0
    kv("nnz", f"{nnz:,}  {C.DIM}({sparsity:.2f}% sparse){C.RESET}")

    looks_int = float(xmin).is_integer() and float(xmax).is_integer()
    kv("looks like", f"{C.GREEN}raw counts{C.RESET}" if looks_int
                     else f"{C.YELLOW}transformed/normalized{C.RESET}")


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def inspect_h5ad(input_file: str, max_cats: int, x_stats: bool) -> None:
    print(f"{C.BOLD}Reading:{C.RESET} {input_file} {C.DIM}(backed mode){C.RESET}")
    adata = ad.read_h5ad(input_file, backed="r")
    try:
        print_overview(adata, input_file)
        print_obs(adata, max_cats)
        print_var(adata)
        print_slots(adata)
        if x_stats:
            print_x_stats(adata)
        else:
            header("X statistics")
            print(f"  {C.DIM}skipped (--no-x-stats){C.RESET}")
    finally:
        if adata.isbacked:
            adata.file.close()
    print(f"\n{C.GREEN}{C.BOLD}✓ Done.{C.RESET}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Inspect an AnnData .h5ad file.")
    parser.add_argument("input", help="Input .h5ad file")
    parser.add_argument("--max-cats", type=int, default=10,
                        help="Max distinct values shown per obs column (default: 10)")
    parser.add_argument("--no-x-stats", action="store_true",
                        help="Skip loading X for min/max/nnz (faster on huge matrices)")
    parser.add_argument("--no-color", action="store_true",
                        help="Disable ANSI colors (also honors NO_COLOR env var)")
    args = parser.parse_args()

    if not color_enabled(args.no_color):
        C.disable()

    if not os.path.exists(args.input):
        sys.exit(f"{C.RED}Error:{C.RESET} file not found: {args.input}")

    inspect_h5ad(args.input, max_cats=args.max_cats, x_stats=not args.no_x_stats)


if __name__ == "__main__":
    main()
