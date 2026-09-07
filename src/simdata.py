"""
Paired counterfactual simulation (direction D1) and shared preprocessing.

``paired_simulation`` gives two count matrices for the *same cells* differing
only by the gene-wise batch factors. Everything here preserves that pairing:
subsampling for the confounding sweep uses one index vector applied to both
arms, and gene filtering uses one gene mask applied to both arms.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import scanpy as sc

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code" / "splatpy"))
from splatpy import SplatParams, paired_simulation, to_anndata  # noqa: E402


#: Branching-trajectory topology used throughout D1: a root path that splits
#: twice, giving four terminal states - the minimal structure that has both a
#: continuous axis (Step) and a discrete one (Group).
BASE = dict(
    n_genes=2000,
    group_prob=(0.35, 0.25, 0.25, 0.15),
    path_from=(0, 1, 1, 2),
    de_prob=0.25,
    de_fac_loc=0.35,
    de_fac_scale=0.4,
    path_n_steps=100,
    path_nonlinear_prob=0.1,
    path_sigma_fac=0.5,
    dropout_type="experiment",
    dropout_mid=1.5,
    dropout_shape=-1.0,
    lib_loc=9.5,
    lib_scale=0.2,
    bcv_common=0.1,
)


def simulate_pair(seed, batch_fac_loc, n_batches=3, pool_per_batch=900, **overrides):
    """Return (adata_with_batch, adata_without_batch) - identical cells.

    ``batch_fac_loc`` is the log-normal location of the per-gene batch factors,
    i.e. the batch-effect strength knob.
    """
    kw = dict(BASE)
    kw.update(overrides)
    p = SplatParams(
        batch_cells=tuple([pool_per_batch] * n_batches),
        batch_fac_loc=batch_fac_loc,
        batch_fac_scale=max(batch_fac_loc, 0.05),
        seed=seed,
        **kw,
    )
    with_b, no_b = paired_simulation(p, method="paths")
    a_with, a_no = to_anndata(with_b), to_anndata(no_b)
    # the pairing is the whole point - assert it rather than trust it
    assert (a_with.obs["Group"].values == a_no.obs["Group"].values).all()
    assert (a_with.obs["Step"].values == a_no.obs["Step"].values).all()
    assert (a_with.obs["Batch"].values == a_no.obs["Batch"].values).all()
    for a in (a_with, a_no):
        a.obs["TrueTime"] = _global_pseudotime(a.obs, kw["path_from"], kw["path_n_steps"])
    return a_with, a_no


def _global_pseudotime(obs, path_from, n_steps):
    """Splat's ``Step`` is position *within* a path; a cell on a child path is
    further along the lineage than its raw Step suggests. Global pseudotime adds
    the depth of the path in the ``path_from`` tree, giving a single ordering
    that DPT can legitimately be compared against."""
    depth = []
    for i in range(len(path_from)):
        d, cur = 0, path_from[i]
        while cur != 0:
            d += 1
            cur = path_from[cur - 1]
        depth.append(d)
    depth = np.asarray(depth, dtype=float)
    g = np.asarray([int(str(x).replace("Path", "")) for x in obs["Group"].values])
    return depth[g - 1] * n_steps + np.asarray(obs["Step"].values, dtype=float)


def confounded_index(obs, confounding, n_keep_per_batch, seed=0):
    """Subsample indices so that batch composition is confounded with branch.

    ``confounding`` = 0 gives each batch the same branch mix; = 1 gives each
    batch only its assigned branch. Intermediate values interpolate the
    sampling weights. This is the composition-imbalance axis of Maan et al.
    """
    rng = np.random.default_rng(seed)
    batches = np.asarray(obs["Batch"].values)
    groups = np.asarray(obs["Group"].values)
    ub, ug = np.unique(batches), np.unique(groups)
    keep = []
    for bi, b in enumerate(ub):
        pref = ug[bi % len(ug)]
        w = np.where(groups == pref,
                     (1 - confounding) / len(ug) + confounding,
                     (1 - confounding) / len(ug))
        w = np.where(batches == b, w, 0.0)
        if w.sum() <= 0:
            raise ValueError("no cells available for batch %s" % b)
        w = w / w.sum()
        n = min(n_keep_per_batch, int((batches == b).sum()))
        keep.append(rng.choice(len(batches), size=n, replace=False, p=w))
    idx = np.sort(np.concatenate(keep))
    return idx


def preprocess_pair(a_with, a_no, min_cells=3, n_top_genes=None):
    """Normalise both arms with one shared gene mask, keeping the cells paired."""
    mask = (np.asarray((a_with.X > 0).sum(axis=0)).ravel() >= min_cells) & \
           (np.asarray((a_no.X > 0).sum(axis=0)).ravel() >= min_cells)
    out = []
    for a in (a_with, a_no):
        a = a[:, mask].copy()
        a.layers["counts"] = a.X.copy()
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
        out.append(a)
    if n_top_genes is not None and n_top_genes < out[0].n_vars:
        # HVGs chosen on the *batch-free* arm so the choice cannot itself be a
        # batch artefact, then applied to both
        hv = sc.pp.highly_variable_genes(out[1], n_top_genes=n_top_genes, inplace=False)
        sel = np.asarray(hv["highly_variable"].values)
        out = [a[:, sel].copy() for a in out]
    return out[0], out[1]


def root_cell(obs):
    """Index of the earliest cell on the root branch - the DPT root."""
    groups = np.asarray(obs["Group"].values)
    steps = np.asarray(obs["Step"].values)
    root_group = sorted(np.unique(groups))[0]
    cand = np.where(groups == root_group)[0]
    return int(cand[np.argmin(steps[cand])])
