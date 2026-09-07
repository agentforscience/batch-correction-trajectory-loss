"""Pre-flight validation of the three assumptions the D1 design rests on.

1. The batch-free arm is invariant to the batch-strength knob (so the null arm
   needs to be run once per seed, not once per strength).
2. TVR responds monotonically to batch strength (the feasibility run's
   max|rho(PC, Step)| did not, which is why it was replaced).
3. Harmony actually does something at strong batch effects under the
   configuration used here (the unresolved item carried over from Phase 1).
"""
import sys, time
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from simdata import simulate_pair, preprocess_pair, confounded_index, root_cell
from metrics import trajectory_variance, batch_metrics, make_neighbors
import methods as M

t0 = time.time()
print("=" * 70)
print("[1] batch-free arm invariance to batch_fac_loc")
refs = {}
for bfl in [0.0, 0.2, 0.5]:
    aw, an = simulate_pair(seed=1, batch_fac_loc=bfl, pool_per_batch=300)
    refs[bfl] = an.X.copy()
    print(f"    bfl={bfl}: with-batch sum={aw.X.sum():.0f}  no-batch sum={an.X.sum():.0f}")
same = all(np.array_equal(refs[0.0], refs[b]) for b in refs)
print(f"    -> batch-free arms identical across strengths: {same}")

print("=" * 70)
print("[2] TVR sensitivity to batch strength (uncorrected vs oracle)")
rows = []
for bfl in [0.0, 0.1, 0.2, 0.35, 0.5]:
    aw, an = simulate_pair(seed=1, batch_fac_loc=bfl, pool_per_batch=400)
    aw, an = preprocess_pair(aw, an)
    e_w, _ = M.run_uncorrected(aw, "Batch", "Group", seed=0)
    e_o, _ = M.run_uncorrected(an, "Batch", "Group", seed=0)
    step, grp = aw.obs["Step"].values, aw.obs["Group"].values
    tv_w = trajectory_variance(e_w, step, grp)
    tv_o = trajectory_variance(e_o, step, grp)
    bm = batch_metrics(e_w, aw.obs["Batch"].values, grp)
    rows.append(dict(bfl=bfl, tv_uncorr=tv_w, tv_oracle=tv_o, TVR=tv_w / tv_o,
                     kbet=bm["kbet"], asw_batch=bm["asw_batch"]))
    print(f"    {rows[-1]}")
df = pd.DataFrame(rows)
print(f"    -> TVR range {df.TVR.min():.3f}..{df.TVR.max():.3f}  (needs to move)")

print("=" * 70)
print("[3] Harmony diagnostic at strong batch effect")
aw, an = simulate_pair(seed=1, batch_fac_loc=0.4, pool_per_batch=400)
aw, an = preprocess_pair(aw, an)
grp = aw.obs["Group"].values; bat = aw.obs["Batch"].values
base, _ = M.run_uncorrected(aw, "Batch", "Group", seed=0)
b0 = batch_metrics(base, bat, grp)
print(f"    uncorrected: kbet={b0['kbet']:.3f} asw_batch={b0['asw_batch']:.3f}")
hrows = []
for theta in [1.0, 2.0, 4.0]:
    for nclust in [10, 30, None]:
        e, meta = M.run_harmony(aw, "Batch", "Group", seed=0, theta=theta,
                                nclust=nclust, max_iter_harmony=30)
        bm = batch_metrics(e, bat, grp)
        hrows.append(dict(theta=theta, nclust=meta["nclust"], iters=meta["n_iter"],
                          kbet=bm["kbet"], asw_batch=bm["asw_batch"], sec=meta["seconds"]))
        print(f"    theta={theta} nclust={meta['nclust']} iters={meta['n_iter']} "
              f"kbet={bm['kbet']:.3f} asw={bm['asw_batch']:.3f}")
pd.DataFrame(hrows).to_csv("results/harmony_diagnostic.csv", index=False)
df.to_csv("results/tvr_sensitivity_check.csv", index=False)
print(f"\ntotal {time.time()-t0:.1f}s")
