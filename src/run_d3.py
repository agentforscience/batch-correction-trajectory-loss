"""
D3 (experiment E5) - real multi-donor batches and a real-data null.

Two arms, because on real data there is no counterfactual:

**Arm A - real donors.** The HSPC -> erythrocyte lineage of ``Immune_ALL_human``
(the 2,629 cells scIB annotated with ``dpt_pseudotime``), spanning three real
donor batches whose cell-type composition is strongly imbalanced. No oracle
exists, so biology is scored two ways: variance attributable to the published
pseudotime, and **within-donor pseudotime preservation** - within one donor
there is no batch effect, so any reordering of that donor's own cells is
attributable to the method rather than to correction.

**Arm B - real-data null.** One donor only, split into random pseudobatches.
There is no batch effect at all, so the uncorrected embedding of those same
cells *is* an exact oracle and the full TVR estimator applies. This is
Antonsson & Melsted's calibration, scored on the D1/D2 biology axis.
"""
from __future__ import annotations

import argparse
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import methods as M
import metrics as MT
from run_d1 import score_embedding, dpt_reference, scib_trajectory_conservation, append

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
IMMUNE = ROOT / "datasets" / "scib" / "Immune_ALL_human.h5ad"


def within_batch_pseudotime(emb, batch, truetime, seed=0):
    """Mean |Spearman(DPT within a donor's own cells, published pseudotime)|.

    Computed donor by donor on the *corrected* embedding. Because a single
    donor contains no batch effect, a drop here is damage the method did, not
    batch structure it removed.
    """
    emb = np.asarray(emb)
    batch = np.asarray(batch)
    tt = np.asarray(truetime, dtype=float)
    out = []
    for b in np.unique(batch):
        m = batch == b
        if m.sum() < 60:
            continue
        sub = emb[m]
        root = int(np.argmin(tt[m]))
        try:
            dpt = MT.dpt_from_embedding(sub, root, seed=seed)
            ok = np.isfinite(dpt)
            if ok.sum() < 30:
                continue
            r = spearmanr(dpt[ok], tt[m][ok]).statistic
            if np.isfinite(r):
                out.append(abs(float(r)))
        except Exception:
            continue
    return float(np.mean(out)) if out else np.nan


def _prep(a, n_hvg=2000, seed=0):
    a = a.copy()
    a.X = a.layers["counts"].copy()
    a = a[:, np.asarray(a.X.sum(axis=0)).ravel() > 0].copy()
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    hv = sc.pp.highly_variable_genes(a, n_top_genes=min(n_hvg, a.n_vars),
                                     batch_key="Batch", inplace=False)
    a = a[:, np.asarray(hv["highly_variable"].values)].copy()
    return a


def load_lineage(arm, seed=0, donor="Oetjen_P"):
    import anndata as ad
    a = ad.read_h5ad(IMMUNE)
    a = a[np.isfinite(a.obs["dpt_pseudotime"].values)].copy()
    a.obs["TrueTime"] = a.obs["dpt_pseudotime"].astype(float).values
    a.obs["Group"] = a.obs["final_annotation"].astype(str).values
    if arm == "real":
        a.obs["Batch"] = a.obs["batch"].astype(str).values
    else:                                    # null: one donor, random pseudobatches
        a = a[a.obs["batch"].astype(str) == donor].copy()
        rng = np.random.default_rng(seed)
        a.obs["Batch"] = [f"Pseudo{i+1}" for i in rng.integers(0, 3, a.n_obs)]
    a.obs["Batch"] = pd.Categorical(a.obs["Batch"])
    a.obs["Group"] = pd.Categorical(a.obs["Group"])
    return _prep(a, seed=seed)


def run_unit(arm, seed, methods):
    t0 = time.time()
    a = load_lineage(arm, seed=seed)
    root = int(np.argmin(a.obs["TrueTime"].values))
    tt = np.asarray(a.obs["TrueTime"].values)
    grp = np.asarray(a.obs["Group"].astype(str).values)
    bat = np.asarray(a.obs["Batch"].astype(str).values)

    # Arm B has a genuine oracle (same cells, no batch effect); Arm A does not,
    # so its reference is the uncorrected embedding and TVR must be read as
    # "relative to unintegrated", not "relative to truth".
    oracle_emb, _ = M.run_uncorrected(a, "Batch", "Group", seed=seed)
    oracle_tv = MT.trajectory_variance(oracle_emb, tt, grp, seed=seed)
    pre_ref = dpt_reference(a, root, seed)

    from sklearn.metrics import normalized_mutual_info_score
    common = dict(arm=arm, seed=seed, n_cells=a.n_obs, n_genes=a.n_vars,
                  conf_nmi=round(float(normalized_mutual_info_score(bat, grp)), 4),
                  oracle_traj_var=oracle_tv,
                  has_true_oracle=(arm == "null"), tag="E5")
    rows = []
    for name in ["uncorrected"] + [m for m in methods if m != "uncorrected"]:
        try:
            emb, meta = M.METHODS[name](a, "Batch", "Group", seed=seed)
            sc_ = score_embedding(emb, a, oracle_emb, oracle_tv, root, seed)
            sc_["scib_traj_cons"] = scib_trajectory_conservation(pre_ref, emb, root, seed)
            sc_["within_batch_rho"] = within_batch_pseudotime(emb, bat, tt, seed)
            rows.append({**common, "method": name, **sc_,
                         "seconds": meta.get("seconds", np.nan),
                         "supervised": name in M.SUPERVISED})
        except Exception as e:
            traceback.print_exc()
            rows.append({**common, "method": name, "error": str(e)[:200]})
    print(f"  D3 arm={arm} seed={seed} cells={a.n_obs} -> {time.time()-t0:.0f}s", flush=True)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3])
    ap.add_argument("--out", default=str(OUT / "d3_raw.csv"))
    ap.add_argument("--slice", type=int, nargs=2, default=None)
    args = ap.parse_args()
    methods = list(M.METHODS)
    # Arm A has no seed-dependent randomisation in the data itself, but the
    # methods do, so it is still repeated across seeds.
    plan = [dict(arm=arm, seed=s) for arm in ["real", "null"] for s in args.seeds]
    if args.slice:
        i, n = args.slice
        plan = plan[i::n]
    print(f"D3: {len(plan)} units", flush=True)
    t0 = time.time()
    for i, cfg in enumerate(plan):
        print(f"[{i+1}/{len(plan)}]", flush=True)
        append(run_unit(methods=methods, **cfg), args.out)
    print(f"D3 done in {(time.time()-t0)/60:.1f} min -> {args.out}")


if __name__ == "__main__":
    main()
