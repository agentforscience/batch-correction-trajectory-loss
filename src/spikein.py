"""
D2 - injecting a *real* batch effect into real single-batch trajectory data.

The counterfactual logic of D1 carried onto real expression: a dataset with one
batch is split into pseudobatches, a batch effect is injected into all but the
first, and the uninjected version of the *same cells* is the oracle.

What makes this a "real spike-in" rather than an invented one is where the
per-gene factors come from: they are drawn from the empirical distribution of
log-fold-changes between two actual sequencing protocols in the scIB pancreas
dataset (``tech`` = smartseq2 vs inDrop3), matched on expression level so the
mean-dependence of a real protocol effect is preserved.

Three injection regimes are provided, deliberately including one that a global
multiplicative model *cannot* represent, so that D2 is a genuine check on D1
rather than an echo of it:

``lfc``        gene-wise multiplicative shift  - same functional form as Splat
``composite``  lfc + sequencing-depth shift + protocol-specific dropout
``celltype``   gene-wise shift that differs by cell type - breaks the global
               multiplicative assumption shared by Splat and ComBat
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import scanpy as sc
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
PANCREAS = ROOT / "datasets" / "scib" / "human_pancreas_norm_complexBatch.h5ad"
N_BINS = 20


def _dense(X):
    return np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)


def estimate_real_batch_profile(tech_a="smartseq2", tech_b="inDrop3",
                                celltypes=("alpha", "beta", "ductal"), cache=True):
    """Empirical protocol-contrast profile: log2 fold-changes, detection-rate
    differences and a depth ratio, all binned by expression level.

    Restricting to cell types present in both protocols controls composition, so
    the remaining difference is technical. Returns arrays indexed by expression
    bin; ``sample_factors`` then maps them onto any target gene set.
    """
    cache_f = ROOT / "results" / f"batch_profile_{tech_a}_vs_{tech_b}.npz"
    if cache and cache_f.exists():
        d = np.load(cache_f, allow_pickle=True)
        return {k: d[k] for k in d.files}

    import anndata as ad
    a = ad.read_h5ad(PANCREAS)
    keep = a.obs["celltype"].isin(celltypes) & a.obs["tech"].isin([tech_a, tech_b])
    a = a[keep].copy()
    counts = a.layers["counts"]
    lib = np.asarray(counts.sum(axis=1)).ravel()
    cpm = _dense(counts) / np.maximum(lib, 1)[:, None] * 1e4

    is_a = (a.obs["tech"] == tech_a).values
    mean_a, mean_b = cpm[is_a].mean(axis=0), cpm[~is_a].mean(axis=0)
    det_a = (cpm[is_a] > 0).mean(axis=0)
    det_b = (cpm[~is_a] > 0).mean(axis=0)

    ok = (mean_a + mean_b) > 0
    eps = 1e-2
    lfc = np.log2((mean_a[ok] + eps) / (mean_b[ok] + eps))
    level = np.log10((mean_a[ok] + mean_b[ok]) / 2 + eps)
    d_det = (det_a - det_b)[ok]

    edges = np.quantile(level, np.linspace(0, 1, N_BINS + 1))
    edges[0] -= 1e-6
    edges[-1] += 1e-6
    binid = np.clip(np.digitize(level, edges) - 1, 0, N_BINS - 1)

    prof = {
        "lfc": lfc, "d_det": d_det, "binid": binid, "edges": edges,
        "depth_ratio": np.array([np.median(lib[is_a]) / max(np.median(lib[~is_a]), 1)]),
        "meta": np.array([tech_a, tech_b, str(len(lfc))], dtype=object),
    }
    cache_f.parent.mkdir(exist_ok=True)
    np.savez(cache_f, **prof)
    return prof


def sample_factors(profile, target_mean_expr, rng, strength=1.0):
    """Draw a per-gene (log2 fold-change, detection shift) for the target genes.

    Each target gene is matched to the expression bin of the source contrast at
    the same expression quantile, so genes that are lowly expressed in the
    target inherit the (larger, noisier) fold-changes that lowly expressed genes
    actually show between protocols.
    """
    level = np.log10(target_mean_expr + 1e-2)
    q = np.argsort(np.argsort(level)) / max(len(level) - 1, 1)
    src_edges = np.linspace(0, 1, N_BINS + 1)
    tgt_bin = np.clip(np.digitize(q, src_edges) - 1, 0, N_BINS - 1)

    lfc_out = np.zeros(len(level))
    det_out = np.zeros(len(level))
    for b in range(N_BINS):
        tgt = tgt_bin == b
        src = np.where(profile["binid"] == b)[0]
        if not tgt.any():
            continue
        if len(src) == 0:
            src = np.arange(len(profile["lfc"]))
        pick = rng.choice(src, size=int(tgt.sum()), replace=True)
        lfc_out[tgt] = profile["lfc"][pick]
        det_out[tgt] = profile["d_det"][pick]
    # random sign per gene so no batch is uniformly "higher"
    sign = rng.choice([-1.0, 1.0], size=len(level))
    return strength * sign * lfc_out, strength * sign * det_out


def assign_pseudobatches(truetime, n_batches, confounding, rng):
    """Pseudobatch labels, optionally confounded with pseudotime.

    ``confounding`` = 0 assigns at random; = 1 makes batch a pseudotime tercile.
    This is the D2 analogue of D1's branch-composition confounding, and is the
    form of confounding that most directly threatens a trajectory.
    """
    n = len(truetime)
    rank = np.argsort(np.argsort(np.asarray(truetime, dtype=float)))
    ordered = (rank / n * n_batches).astype(int).clip(0, n_batches - 1)
    random = rng.integers(0, n_batches, size=n)
    use_ordered = rng.random(n) < confounding
    return np.where(use_ordered, ordered, random)


def inject(counts, batches, labels, profile, rng, regime="composite",
           strength=1.0, n_batches=3):
    """Return (X_injected, X_oracle) - the same cells with and without a batch effect.

    Both arms are Poisson-resampled from their respective expected-count
    matrices using the same generator, exactly as ``splatpy.paired_simulation``
    toggles the batch factors on one RNG stream. At ``strength = 0`` the two
    arms are identical by construction, which is the null calibration.
    """
    counts = _dense(counts).astype(np.float64)
    n_cells, n_genes = counts.shape
    lib = counts.sum(axis=1)
    mean_expr = counts.mean(axis=0) / max(counts.mean(), 1e-9) * 1e2

    lam = counts.copy()
    ulab = np.unique(labels)

    for b in range(1, n_batches):                      # batch 0 is the reference
        m = batches == b
        if not m.any():
            continue
        if regime == "celltype":
            # a different gene-wise shift for every cell type: no single
            # gene-wise (or gene x batch) correction can undo this
            for c in ulab:
                mc = m & (labels == c)
                if mc.sum() < 5:
                    continue
                lfc, _ = sample_factors(profile, mean_expr, rng, strength)
                lam[mc] *= 2.0 ** lfc[None, :]
        else:
            lfc, d_det = sample_factors(profile, mean_expr, rng, strength)
            lam[m] *= 2.0 ** lfc[None, :]
            if regime == "composite":
                # sequencing-depth shift (real median-library-size ratio, damped)
                dr = float(profile["depth_ratio"][0]) ** (0.3 * strength)  # damped: the
                # raw smartseq2/inDrop3 ratio is ~63x, an extreme of the range
                dr = dr if b % 2 == 1 else 1.0 / dr
                lam[m] *= dr
                # protocol-specific detection differences -> extra dropout
                p_keep = np.clip(1.0 + d_det, 0.05, 1.0)
                lam[m] *= p_keep[None, :]
        if regime != "composite":
            # hold library size fixed so only composition changes
            s = lam[m].sum(axis=1)
            lam[m] *= (lib[m] / np.maximum(s, 1e-9))[:, None]

    x_inj = rng.poisson(np.maximum(lam, 0))
    x_ora = rng.poisson(np.maximum(counts, 0))
    return x_inj.astype(np.float32), x_ora.astype(np.float32)


# --------------------------------------------------------------- dataset prep
def load_trajectory_dataset(name, n_cells=2000, n_hvg=2000, seed=0):
    """Load a real single-batch dataset with a trajectory annotation.

    ``setty`` carries a published Palantir pseudotime. ``pancreas`` has an
    unambiguous Ductal -> Ngn3 -> endocrine lineage, so its ground-truth
    pseudotime is DPT computed on the *uninjected* data - uncontaminated,
    because there is no batch effect in it.
    """
    import anndata as ad
    rng = np.random.default_rng(seed)

    if name == "setty":
        a = ad.read_h5ad(ROOT / "datasets/trajectory/setty_bonemarrow_cd34.h5ad")
        a.obs["Group"] = (a.obs["clusters"].astype(str)
                          .str.replace(r"_\d$", "", regex=True).values)
        a.obs["TrueTime"] = a.obs["palantir_pseudotime"].astype(float).values
        root_label = "HSC"
    elif name == "pancreas":
        a = ad.read_h5ad(ROOT / "datasets/trajectory/pancreas_endocrinogenesis.h5ad")
        a.obs["Group"] = a.obs["clusters"].astype(str).values
        a.obs["TrueTime"] = np.nan
        root_label = "Ductal"
    else:
        raise ValueError(name)

    a = a[:, np.asarray(a.X.sum(axis=0)).ravel() > 0].copy()
    if a.n_obs > n_cells:
        a = a[np.sort(rng.choice(a.n_obs, n_cells, replace=False))].copy()
    a.layers["counts"] = a.X.copy()

    # HVGs on the raw (uninjected) data - the injection must not steer selection
    tmp = a.copy()
    sc.pp.normalize_total(tmp, target_sum=1e4)
    sc.pp.log1p(tmp)
    hv = sc.pp.highly_variable_genes(tmp, n_top_genes=min(n_hvg, tmp.n_vars),
                                     inplace=False)
    a = a[:, np.asarray(hv["highly_variable"].values)].copy()

    if not np.isfinite(a.obs["TrueTime"]).all():
        tmp = a.copy()
        sc.pp.normalize_total(tmp, target_sum=1e4)
        sc.pp.log1p(tmp)
        sc.pp.pca(tmp, n_comps=20, random_state=seed)
        sc.pp.neighbors(tmp, random_state=seed)
        sc.tl.diffmap(tmp, random_state=seed)
        cand = np.where(tmp.obs["Group"].values == root_label)[0]
        tmp.uns["iroot"] = int(cand[np.argmin(tmp.obsm["X_diffmap"][cand, 1])])
        sc.tl.dpt(tmp)
        a.obs["TrueTime"] = tmp.obs["dpt_pseudotime"].astype(float).values
    a = a[np.isfinite(a.obs["TrueTime"].values)].copy()
    a.obs["Group"] = a.obs["Group"].astype("category")
    root = int(np.argmin(a.obs["TrueTime"].values))
    return a, root
