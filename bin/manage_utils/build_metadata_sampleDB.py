#!/usr/bin/env python3
"""
Build metadata_sampleDB.tsv for a cohort.

Row set = UNION of every SampleID across phenotypic, cnv/*, shortvariants/*, and
expression/*.h5ad. Only samples with >=1 molecular modality are kept.

Cohort and Biobank are collected from EVERY source (phenotypic, each sampleDB,
each h5ad obs) and reconciled: a single value is emitted per sample, with a
warning printed on any cross-source disagreement. This ensures molecular-only
samples still get Cohort/Biobank instead of blanks.

Other phenotype columns (Age, Sex, Diagnosis, Ethnicity) come from phenotypic only.

Columns:
  SampleID, Cohort, Biobank, Age, Sex, Diagnosis, Ethnicity,
  in_phenotypic,
  has_cnv, cnv_dataset, has_shortvariants, shortvariants_dataset,
  has_expression, expression_type, Tissue

SampleIDs used as-is. Run from cohort root (dir with cnv/ expression/ phenotypic/ shortvariants/):
    python phenotypic/build_metadata_sampleDB.py
"""
from __future__ import annotations
import sys
from collections import defaultdict
from pathlib import Path
import pandas as pd

ROOT  = Path(".").resolve()
PHENO = ROOT / "phenotypic" / "phenotypic_sampleDB.tsv"
OUT   = ROOT / "metadata_sampleDB.tsv"

ID_CANDIDATES      = ("SampleID",)
COHORT_CANDIDATES  = ("Cohort",)
BIOBANK_CANDIDATES = ("Biobank", "Brain_bank", "BrainBank", "Brain_Bank")
REGION_CANDIDATES  = ("Region",)
TISSUE_CANDIDATES  = ("Tissue",)
ASSAY_TOKENS       = ("snATACseq", "snRNAseq", "snATAC", "snRNA", "bulk_RNAseq", "bulk_rna")
PHENO_ONLY_COLS    = ["Age", "Sex", "Diagnosis", "Ethnicity"]   # Cohort/Biobank handled separately

# global accumulators for reconciliation: sample -> {value: set(sources)}
cohort_obs:  dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
biobank_obs: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))

# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------
def find_col(df_or_obs, candidates) -> str | None:
    lower = {c.lower(): c for c in df_or_obs.columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    idx_name = getattr(df_or_obs.index, "name", None)
    if idx_name and idx_name.lower() in {c.lower() for c in candidates}:
        return "__index__"
    return None

def _series(df_or_obs, col: str) -> pd.Series:
    src = df_or_obs.index.to_series() if col == "__index__" else df_or_obs[col]
    return src.astype(str).str.strip()

def _clean(v: str) -> str | None:
    v = "" if v is None else str(v).strip()
    return None if v.lower() in {"", "nan", "na", "none", "<na>"} else v

def record_cohort_biobank(df_or_obs, ids: pd.Series, source: str) -> None:
    """Accumulate Cohort/Biobank values per sample from any table/obs, lowercased."""
    ccol = find_col(df_or_obs, COHORT_CANDIDATES)
    bcol = find_col(df_or_obs, BIOBANK_CANDIDATES)
    cser = _series(df_or_obs, ccol) if ccol else None
    bser = _series(df_or_obs, bcol) if bcol else None
    for i, sid in ids.items():
        if cser is not None:
            v = _clean(cser.get(i))
            if v: cohort_obs[sid][v.lower()].add(source)
        if bser is not None:
            v = _clean(bser.get(i))
            if v: biobank_obs[sid][v.lower()].add(source)

def read_sampledb_ids(path: Path, source: str) -> set[str]:
    """Read a sampleDB.tsv -> set of SampleIDs, and record Cohort/Biobank."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    col = find_col(df, ID_CANDIDATES)
    if col is None:
        print(f"  ! no SampleID column in {path} (cols={list(df.columns)})", file=sys.stderr)
        return set()
    ids = _series(df, col)
    record_cohort_biobank(df, ids, source)
    return set(ids.dropna().unique())

def tissue_from_filename(path: Path) -> set[str]:
    stem = path.stem
    for tok in ASSAY_TOKENS:
        idx = stem.find(tok + "_")
        if idx != -1:
            return {stem[idx + len(tok) + 1:]}
    return set()

def h5ad_scan(path: Path, source: str) -> tuple[set[str], set[str]]:
    """Unique SampleIDs + Tissue from h5ad obs; also record Cohort/Biobank."""
    import anndata as ad
    a = ad.read_h5ad(path, backed="r")
    obs = a.obs

    idcol = find_col(obs, ID_CANDIDATES)
    ids_src = _series(obs, idcol) if idcol is not None \
              else obs.index.to_series().astype(str).str.strip()
    record_cohort_biobank(obs, ids_src, source)
    ids = set(ids_src.dropna().unique())

    tissues: set[str] = set()
    rcol = find_col(obs, REGION_CANDIDATES)
    if rcol is not None and rcol != "__index__":
        tissues = set(obs[rcol].astype(str).str.strip().dropna().unique())
    if not tissues:
        tcol = find_col(obs, TISSUE_CANDIDATES)
        if tcol is not None and tcol != "__index__":
            vals = set(obs[tcol].astype(str).str.strip().dropna().unique())
            if vals - {"brain", "Brain"}:
                tissues = vals
    if not tissues:
        tissues = tissue_from_filename(path)
    tissues = {t for t in tissues if _clean(t)}

    try:
        a.file.close()
    except Exception:
        pass
    return ids, tissues

def collect_modality(base: Path) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    if not base.exists():
        return out
    for ds_dir in sorted(p for p in base.iterdir() if p.is_dir()):
        sdb = ds_dir / "sampleDB.tsv"
        if sdb.exists():
            out[ds_dir.name] = read_sampledb_ids(sdb, source=f"{base.name}/{ds_dir.name}")
        else:
            print(f"  ! {sdb} missing", file=sys.stderr)
    return out

def list_series(mod_map: dict[str, set[str]]) -> pd.Series:
    per: dict[str, list[str]] = {}
    for ds, ids in mod_map.items():
        for sid in ids:
            per.setdefault(sid, []).append(ds)
    return pd.Series({sid: ",".join(sorted(v)) for sid, v in per.items()}, dtype="string")

def reconcile(store: dict[str, dict[str, set[str]]], field: str,
              prefer: dict[str, str]) -> pd.Series:
    """Collapse per-sample multi-source values to one; warn on disagreement.
    `prefer` gives an authoritative value per sample (from phenotypic) to win ties."""
    out: dict[str, str] = {}
    for sid, valmap in store.items():
        vals = set(valmap.keys())
        if len(vals) == 1:
            out[sid] = next(iter(vals))
        else:
            chosen = prefer.get(sid) or sorted(vals)[0]
            srcs = {v: ",".join(sorted(s)) for v, s in valmap.items()}
            print(f"  ! {field} disagreement for {sid}: {srcs} -> using '{chosen}'",
                  file=sys.stderr)
            out[sid] = chosen
    return pd.Series(out, dtype="string")

# ----------------------------------------------------------------------
# 1) Read every source and collect IDs (+ Cohort/Biobank into accumulators)
# ----------------------------------------------------------------------
pheno_df = pd.read_csv(PHENO, sep="\t", dtype=str)
pheno_df["SampleID"] = pheno_df["SampleID"].astype(str).str.strip()
pheno_df = pheno_df.drop_duplicates("SampleID").set_index("SampleID")
pheno_ids = set(pheno_df.index)
record_cohort_biobank(pheno_df, pheno_df.index.to_series(), source="phenotypic")
print(f"phenotypic: {len(pheno_ids)} samples")

cnv_map = collect_modality(ROOT / "cnv")
sv_map  = collect_modality(ROOT / "shortvariants")
print(f"cnv datasets: {list(cnv_map)}")
print(f"shortvariants datasets: {list(sv_map)}")

expr_base   = ROOT / "expression"
expr_types:  dict[str, set[str]] = {}
expr_tissue: dict[str, set[str]] = {}
if expr_base.exists():
    for h5 in sorted(expr_base.glob("*/*/*.h5ad")):
        assay = h5.parent.name
        try:
            ids, tissues = h5ad_scan(h5, source=f"expression/{assay}")
        except Exception as e:
            print(f"  ! failed reading {h5}: {e}", file=sys.stderr)
            continue
        print(f"  {h5.relative_to(ROOT)}: {len(ids)} donors, assay={assay}, tissue={tissues}")
        for sid in ids:
            expr_types.setdefault(sid, set()).add(assay)
            if tissues:
                expr_tissue.setdefault(sid, set()).update(tissues)

# ----------------------------------------------------------------------
# 2) UNION roster
# ----------------------------------------------------------------------
all_ids = set(pheno_ids)
for m in (cnv_map, sv_map):
    for ids in m.values():
        all_ids |= ids
all_ids |= set(expr_types.keys())
print(f"\nUNION roster: {len(all_ids)} samples")

meta = pd.DataFrame(index=pd.Index(sorted(all_ids), name="SampleID"))

# ----------------------------------------------------------------------
# 3) Reconciled Cohort / Biobank (from all sources) + phenotype-only cols
# ----------------------------------------------------------------------
# phenotypic values act as the authoritative tie-breaker
pref_cohort  = {s: _clean(pheno_df.at[s, find_col(pheno_df, COHORT_CANDIDATES)]).lower()
                for s in pheno_ids
                if find_col(pheno_df, COHORT_CANDIDATES)
                and _clean(pheno_df.at[s, find_col(pheno_df, COHORT_CANDIDATES)])}
pref_biobank = {s: _clean(pheno_df.at[s, find_col(pheno_df, BIOBANK_CANDIDATES)]).lower()
                for s in pheno_ids
                if find_col(pheno_df, BIOBANK_CANDIDATES)
                and _clean(pheno_df.at[s, find_col(pheno_df, BIOBANK_CANDIDATES)])}

cohort_final  = reconcile(cohort_obs,  "Cohort",  pref_cohort)
biobank_final = reconcile(biobank_obs, "Biobank", pref_biobank)
meta["Cohort"]  = meta.index.map(cohort_final)
meta["Biobank"] = meta.index.map(biobank_final)

for c in PHENO_ONLY_COLS:
    if c in pheno_df.columns:
        meta[c] = meta.index.map(pheno_df[c])
meta["in_phenotypic"] = meta.index.isin(pheno_ids)

# ----------------------------------------------------------------------
# 4) Modality flags + dataset lists
# ----------------------------------------------------------------------
meta["cnv_dataset"]           = meta.index.map(list_series(cnv_map)).fillna("")
meta["has_cnv"]               = meta["cnv_dataset"].ne("")
meta["shortvariants_dataset"] = meta.index.map(list_series(sv_map)).fillna("")
meta["has_shortvariants"]     = meta["shortvariants_dataset"].ne("")

s_types  = pd.Series({s: ",".join(sorted(v)) for s, v in expr_types.items()},  dtype="string")
s_tissue = pd.Series({s: ",".join(sorted(v)) for s, v in expr_tissue.items()}, dtype="string")
meta["expression_type"] = meta.index.map(s_types).fillna("")
meta["Region"]          = meta.index.map(s_tissue).fillna("")
meta["has_expression"]  = meta["expression_type"].ne("")

# ----------------------------------------------------------------------
# 4b) Keep only samples with >=1 molecular modality
# ----------------------------------------------------------------------
keep = meta["has_cnv"] | meta["has_shortvariants"] | meta["has_expression"]
print(f"\ndropping {int((~keep).sum())} phenotype-only samples (no cnv/sv/expression)")
meta = meta[keep]

# ----------------------------------------------------------------------
# 5) Column order + write
# ----------------------------------------------------------------------
order = (["Cohort", "Biobank"] + [c for c in PHENO_ONLY_COLS if c in meta.columns]
         + ["in_phenotypic",
            "has_cnv", "cnv_dataset",
            "has_shortvariants", "shortvariants_dataset",
            "has_expression", "expression_type", "Region"])
meta = meta[order].reset_index()
meta.to_csv(OUT, sep="\t", index=False)
print(f"\nwrote {OUT}: {meta.shape[0]} rows x {meta.shape[1]} cols")

# ----------------------------------------------------------------------
# 6) Sanity summary
# ----------------------------------------------------------------------
print("\n-- coverage --")
for c in ["in_phenotypic", "has_cnv", "has_shortvariants", "has_expression"]:
    print(f"{c}: {int(meta[c].sum())}/{len(meta)}")
print("\nCohort:\n", meta["Cohort"].value_counts(dropna=False))
print("\nBiobank:\n", meta["Biobank"].value_counts(dropna=False))
print("\nmissing Cohort:", int(meta["Cohort"].isna().sum()),
      "| missing Biobank:", int(meta["Biobank"].isna().sum()))
mo = meta.loc[~meta["in_phenotypic"], "SampleID"].tolist()
print(f"\nmolecular-only (not in phenotypic): {len(mo)}")
if mo:
    print(mo[:20])