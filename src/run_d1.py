"""
D1 - paired counterfactual simulation (experiments E1, E2, E3).

For each simulation unit we hold two matrices for the *same* cells: one with a
gene-wise batch effect and one without. The batch-free arm, uncorrected, is the
**oracle**: the embedding an ideal integration method would recover. Every
method is scored against it.

Three arms are run:

* ``batch``  - method applied to the data that has a batch effect (E1/E3)
* ``null``   - method applied to the *batch-free* data, which still carries
               batch *labels* (E2). Any biology lost here is pure artefact,
               because there was nothing to correct.
* reference rows ``oracle`` and ``uncorrected``.

Results are appended to results/d1_raw.csv after every unit, so a crash costs at
most one unit.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from simdata import simulate_pair, preprocess_pair, confounded_index, root_cell
import methods as M
import metrics as MT

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)


def batch_biology_confounding(batch, group):
    """Normalised mutual information between batch and branch labels.

    0 = batch composition identical across branches; 1 = batch determines
    branch. Reported as the *realised* confounding, since the sampling target is
    not always attainable.
    """
    from sklearn.metrics import normalized_mutual_info_score
    return float(normalized_mutual_info_score(np.asarray(batch), np.asarray(group)))


def score_embedding(emb, adata, oracle_emb, oracle_tv, root, seed):
    """All batch- and biology-axis metrics for one embedding."""
    bat = np.asarray(adata.obs["Batch"].values)
    grp = np.asarray(adata.obs["Group"].values)
    tt = np.asarray(adata.obs["TrueTime"].values)
    r = {}
    r.update(MT.batch_metrics(emb, bat, grp, seed=seed))
    tv = MT.trajectory_variance(emb, tt, grp, seed=seed)
    r["traj_var"] = tv
    r["TVR"] = tv / oracle_tv if oracle_tv > 0 else np.nan
    r["bio_removed"] = 1.0 - r["TVR"]
    r["step_var"] = MT.step_variance(emb, tt, seed=seed)
    r["rho_true"] = MT.pseudotime_recovery(emb, tt, root, seed=seed)
    r["knn_oracle"] = MT.knn_overlap(emb, oracle_emb, seed=seed)
    r["branch_ari"] = MT.branch_ari(emb, grp, seed=seed)
    r.update(MT.local_trajectory_metrics(emb, tt, grp))
    return r


def dpt_reference(adata_pre, root, seed=0):
    """Pre-integration AnnData carrying scIB's reference DPT.

    Computed once per unit and reused for every method, which is exactly what
    scIB does and saves a diffusion map per method.
    """
    import scanpy as sc
    pre = adata_pre.copy()
    pre.obsm["X_pca"] = M._pca(pre.X, seed=seed)
    sc.pp.neighbors(pre, use_rep="X_pca", random_state=seed)
    sc.tl.diffmap(pre, random_state=seed)
    pre.uns["iroot"] = int(root)
    sc.tl.dpt(pre)
    return pre


def scib_trajectory_conservation(pre, emb, root, seed=0):
    """scIB's published trajectory metric, for the metric audit (E6).

    Deliberately uses scIB's own semantics: DPT after integration compared with
    DPT on the *unintegrated* data, not with the truth.
    """
    import scanpy as sc
    import scib
    try:
        post = pre.copy()
        post.obsm["X_emb"] = np.asarray(emb, dtype=np.float32)
        sc.pp.neighbors(post, use_rep="X_emb", random_state=seed)
        sc.tl.diffmap(post, random_state=seed)
        post.uns["iroot"] = int(root)
        sc.tl.dpt(post)
        return float(scib.metrics.trajectory_conservation(
            pre, post, label_key="Group", pseudotime_key="dpt_pseudotime"))
    except Exception:
        return np.nan


def run_unit(seed, bfl, confounding, methods, arms=("batch",), n_keep=400,
             n_batches=3, pool=900, do_scib_traj=True, tag=""):
    t0 = time.time()
    a_with, a_no = simulate_pair(seed=seed, batch_fac_loc=bfl,
                                 n_batches=n_batches, pool_per_batch=pool)
    idx = confounded_index(a_with.obs, confounding, n_keep, seed=seed)
    a_with, a_no = a_with[idx].copy(), a_no[idx].copy()
    a_with, a_no = preprocess_pair(a_with, a_no)
    root = root_cell(a_with.obs)
    conf_nmi = batch_biology_confounding(a_with.obs["Batch"], a_with.obs["Group"])

    # ---- oracle: batch-free data, no correction
    oracle_emb, _ = M.run_uncorrected(a_no, "Batch", "Group", seed=seed)
    tt = np.asarray(a_with.obs["TrueTime"].values)
    grp = np.asarray(a_with.obs["Group"].values)
    oracle_tv = MT.trajectory_variance(oracle_emb, tt, grp, seed=seed)

    pre_ref = dpt_reference(a_with, root, seed) if do_scib_traj else None

    common = dict(seed=seed, batch_fac_loc=bfl, confounding=confounding,
                  conf_nmi=round(conf_nmi, 4), n_cells=a_with.n_obs,
                  n_genes=a_with.n_vars, oracle_traj_var=oracle_tv, tag=tag)
    rows = []

    rows.append({**common, "arm": "oracle", "method": "oracle",
                 **score_embedding(oracle_emb, a_with, oracle_emb, oracle_tv, root, seed),
                 "seconds": 0.0, "supervised": False})

    for arm in arms:
        adata = a_with if arm == "batch" else a_no
        for name in methods:
            try:
                emb, meta = M.METHODS[name](adata, "Batch", "Group", seed=seed)
                sc_ = score_embedding(emb, a_with, oracle_emb, oracle_tv, root, seed)
                if do_scib_traj and arm == "batch":
                    sc_["scib_traj_cons"] = scib_trajectory_conservation(pre_ref, emb, root, seed)
                rows.append({**common, "arm": arm, "method": name, **sc_,
                             "seconds": meta.get("seconds", np.nan),
                             "supervised": name in M.SUPERVISED})
            except Exception as e:
                traceback.print_exc()
                rows.append({**common, "arm": arm, "method": name,
                             "error": str(e)[:200], "supervised": name in M.SUPERVISED})
    print(f"  unit seed={seed} bfl={bfl} conf={confounding} arms={arms} "
          f"cells={a_with.n_obs} nmi={conf_nmi:.3f} -> {time.time()-t0:.0f}s", flush=True)
    return rows


def append(rows, path):
    df = pd.DataFrame(rows)
    df.to_csv(path, mode="a", header=not Path(path).exists(), index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3, 4, 5])
    ap.add_argument("--out", default=str(OUT / "d1_raw.csv"))
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--slice", type=int, nargs=2, default=None,
                    metavar=("I", "N"), help="run only units i, i+N, i+2N, ... "
                    "so the plan can be split across GPUs")
    args = ap.parse_args()

    methods = list(M.METHODS)
    seeds = args.seeds
    if args.smoke:
        seeds = seeds[:1]

    plan = []
    # E1 - batch-strength sweep, unconfounded
    for s in seeds:
        for bfl in ([0.2] if args.smoke else [0.0, 0.1, 0.2, 0.35, 0.5]):  # noqa: E501
            plan.append(dict(seed=s, bfl=bfl, confounding=0.0, arms=("batch",), tag="E1"))
    # E2 - null arm: methods applied to batch-free data (independent of bfl,
    # verified in src/validate_setup.py), so run once per seed
    for s in seeds:
        plan.append(dict(seed=s, bfl=0.0, confounding=0.0, arms=("null",), tag="E2"))
    # E3 - confounding sweep at a moderate batch effect
    if not args.smoke:
        for s in seeds:
            for c in [0.0, 0.35, 0.7]:
                plan.append(dict(seed=s, bfl=0.25, confounding=c, arms=("batch",), tag="E3"))

    if args.slice is not None:
        i, n = args.slice
        plan = plan[i::n]
    print(f"D1: {len(plan)} units x {len(methods)} methods", flush=True)
    t0 = time.time()
    for i, cfg in enumerate(plan):
        print(f"[{i+1}/{len(plan)}] {cfg['tag']}", flush=True)
        rows = run_unit(cfg["seed"], cfg["bfl"], cfg["confounding"], methods,
                        arms=cfg["arms"], tag=cfg["tag"])
        append(rows, args.out)
    print(f"D1 done in {(time.time()-t0)/60:.1f} min -> {args.out}")


if __name__ == "__main__":
    main()
