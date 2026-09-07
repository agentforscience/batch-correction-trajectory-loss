"""
D2 (experiment E4) - real-expression spike-in benchmark.

Same measurement as D1, but the expression is real: a real single-batch
trajectory dataset is split into pseudobatches and given a batch effect whose
per-gene factors are sampled from an actual protocol contrast (smartseq2 vs
inDrop3 in the scIB pancreas data). The uninjected version of the same cells is
the oracle.

Regimes (see ``src/spikein.py``):
  ``null``       strength 0 - the two arms are identical, so anything a method
                 changes is pure artefact (Antonsson & Melsted's calibration,
                 here on real expression)
  ``lfc``        gene-wise multiplicative only - the same functional form Splat
                 uses, i.e. the regime most favourable to D1's conclusions
  ``composite``  + depth shift + protocol-specific dropout
  ``celltype``   gene-wise shift that differs per cell type, which no global
                 gene x batch correction model can represent
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
import methods as M
import metrics as MT
from run_d1 import score_embedding, dpt_reference, scib_trajectory_conservation, append
from spikein import (estimate_real_batch_profile, assign_pseudobatches, inject,
                     load_trajectory_dataset)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"

#: (regime, strength, confounding) - the D2 grid
CONFIGS = [
    ("null", 0.0, 0.0),
    ("lfc", 1.0, 0.0),
    ("composite", 0.5, 0.0),
    ("composite", 1.0, 0.0),
    ("composite", 1.0, 0.8),
    ("celltype", 1.0, 0.0),
]


def _to_adata(X, template):
    import anndata as ad
    a = ad.AnnData(X=np.asarray(X, dtype=np.float32), obs=template.obs.copy(),
                   var=template.var.copy())
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    return a


def run_unit(dataset, regime, strength, confounding, seed, methods, profile,
             n_batches=3):
    t0 = time.time()
    base, _ = load_trajectory_dataset(dataset, seed=seed)
    rng = np.random.default_rng(1000 + seed)
    batches = assign_pseudobatches(base.obs["TrueTime"].values, n_batches,
                                   confounding, rng)
    x_inj, x_ora = inject(base.layers["counts"], batches,
                          base.obs["Group"].astype(str).values, profile, rng,
                          regime=regime, strength=strength, n_batches=n_batches)

    tpl = base.copy()
    tpl.obs["Batch"] = pd.Categorical([f"Batch{b+1}" for b in batches])
    a_inj, a_ora = _to_adata(x_inj, tpl), _to_adata(x_ora, tpl)
    keep = (np.asarray((a_inj.X > 0).sum(axis=0)).ravel() >= 3) & \
           (np.asarray((a_ora.X > 0).sum(axis=0)).ravel() >= 3)
    a_inj, a_ora = a_inj[:, keep].copy(), a_ora[:, keep].copy()

    root = int(np.argmin(a_inj.obs["TrueTime"].values))
    oracle_emb, _ = M.run_uncorrected(a_ora, "Batch", "Group", seed=seed)
    tt = np.asarray(a_inj.obs["TrueTime"].values)
    grp = np.asarray(a_inj.obs["Group"].astype(str).values)
    oracle_tv = MT.trajectory_variance(oracle_emb, tt, grp, seed=seed)
    pre_ref = dpt_reference(a_inj, root, seed)

    from sklearn.metrics import normalized_mutual_info_score
    conf_nmi = float(normalized_mutual_info_score(
        np.asarray(a_inj.obs["Batch"].astype(str)), grp))

    common = dict(dataset=dataset, regime=regime, strength=strength,
                  confounding=confounding, seed=seed, conf_nmi=round(conf_nmi, 4),
                  n_cells=a_inj.n_obs, n_genes=a_inj.n_vars,
                  oracle_traj_var=oracle_tv, tag="E4")
    rows = [{**common, "method": "oracle",
             **score_embedding(oracle_emb, a_inj, oracle_emb, oracle_tv, root, seed),
             "seconds": 0.0, "supervised": False}]
    for name in methods:
        try:
            emb, meta = M.METHODS[name](a_inj, "Batch", "Group", seed=seed)
            sc_ = score_embedding(emb, a_inj, oracle_emb, oracle_tv, root, seed)
            sc_["scib_traj_cons"] = scib_trajectory_conservation(pre_ref, emb, root, seed)
            rows.append({**common, "method": name, **sc_,
                         "seconds": meta.get("seconds", np.nan),
                         "supervised": name in M.SUPERVISED})
        except Exception as e:
            traceback.print_exc()
            rows.append({**common, "method": name, "error": str(e)[:200],
                         "supervised": name in M.SUPERVISED})
    print(f"  {dataset} {regime} s={strength} conf={confounding} seed={seed} "
          f"cells={a_inj.n_obs} nmi={conf_nmi:.3f} -> {time.time()-t0:.0f}s", flush=True)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3])
    ap.add_argument("--datasets", nargs="+", default=["setty", "pancreas"])
    ap.add_argument("--out", default=str(OUT / "d2_raw.csv"))
    ap.add_argument("--slice", type=int, nargs=2, default=None)
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()

    profile = estimate_real_batch_profile()
    methods = list(M.METHODS)
    configs = CONFIGS[:2] if args.smoke else CONFIGS
    seeds = args.seeds[:1] if args.smoke else args.seeds

    plan = [dict(dataset=d, regime=r, strength=s, confounding=c, seed=sd)
            for d in args.datasets for (r, s, c) in configs for sd in seeds]
    if args.slice:
        i, n = args.slice
        plan = plan[i::n]
    print(f"D2: {len(plan)} units x {len(methods)} methods", flush=True)
    t0 = time.time()
    for i, cfg in enumerate(plan):
        print(f"[{i+1}/{len(plan)}]", flush=True)
        append(run_unit(methods=methods, profile=profile, **cfg), args.out)
    print(f"D2 done in {(time.time()-t0)/60:.1f} min -> {args.out}")


if __name__ == "__main__":
    main()
