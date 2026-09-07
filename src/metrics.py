"""
Metrics for the batch-correction biology/batch trade-off benchmark.

Two families:

* **Batch axis** - how well batches are mixed. kBET acceptance rate and
  ASW-batch are the two named in the hypothesis; graph iLISI and PCR-batch are
  reported alongside.
* **Biology axis** - how much of the *ground-truth* trajectory survives. The
  primary quantity, TVR (Trajectory Variance Retained), is oracle-referenced:
  it needs a batch-free counterfactual embedding of the *same cells*, which is
  what makes this benchmark different from scIB (see planning.md).

Design note on dimensionality: every embedding is reduced to the same
``n_pcs`` (default 10) before any variance-based metric is computed. Without
this, a low-dimensional latent space (scVI's 10 dims) would score higher than a
20-PC PCA purely because it carries fewer noise directions.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score
from sklearn.cluster import KMeans

import scib_metrics
from scib_metrics.nearest_neighbors import pynndescent

N_PCS_EVAL = 10       # uniform embedding dimensionality for variance metrics
N_NEIGHBORS = 25      # kBET / iLISI neighbourhood size (Tran et al. 2020 use 25)
KNN_OVERLAP_K = 30


# ------------------------------------------------------------------ utilities
def _dense(X):
    return np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)


def reduce_dims(X, n_pcs=N_PCS_EVAL, seed=0):
    """Top-``n_pcs`` PCs of an embedding, so all methods are compared at equal rank.

    Returns the scores matrix and the *unnormalised* explained variance of each
    component (needed as the PCR weights).
    """
    X = _dense(X).astype(np.float64)
    n_pcs = int(min(n_pcs, X.shape[1], X.shape[0] - 1))
    pca = PCA(n_components=n_pcs, random_state=seed)
    Z = pca.fit_transform(X)
    return Z, pca.explained_variance_


def _design_trajectory(step, group, degree=3):
    """Branch-specific polynomial design matrix for a branching trajectory.

    Columns: for each branch, an indicator interacted with 1, t, t^2, t^3 where
    t is the within-branch normalised pseudotime. This lets a PC be "explained
    by the trajectory" even though the trajectory is nonlinear and branching -
    a plain linear regression on Step badly under-counts that (this is exactly
    the insensitivity the feasibility run hit with max|rho(PC, Step)|).
    """
    step = np.asarray(step, dtype=float)
    t = (step - step.min()) / max(float(np.ptp(step)), 1e-12)
    groups = np.asarray(group)
    cols = []
    for g in np.unique(groups):
        m = (groups == g).astype(float)
        for d in range(degree + 1):
            cols.append(m * t**d)
    D = np.column_stack(cols)
    # drop constant/collinear columns
    keep = D.std(axis=0) > 1e-12
    D = D[:, keep]
    return np.column_stack([np.ones(len(D)), D])


def _design_onehot(labels):
    labels = np.asarray(labels)
    u = np.unique(labels)
    D = np.column_stack([(labels == g).astype(float) for g in u[1:]])
    return np.column_stack([np.ones(len(labels)), D]) if D.size else np.ones((len(labels), 1))


def _adj_r2(y, D):
    """Adjusted R^2 of OLS y ~ D (D already contains an intercept column)."""
    y = np.asarray(y, dtype=float)
    n, p = D.shape
    beta, *_ = np.linalg.lstsq(D, y, rcond=None)
    resid = y - D @ beta
    ss_res = float(resid @ resid)
    ss_tot = float(((y - y.mean()) ** 2).sum())
    if ss_tot <= 0:
        return 0.0
    r2 = 1.0 - ss_res / ss_tot
    if n - p <= 1:
        return max(r2, 0.0)
    adj = 1.0 - (1.0 - r2) * (n - 1) / (n - p)
    return float(np.clip(adj, 0.0, 1.0))


def pcr_covariate(X, design, n_pcs=N_PCS_EVAL, seed=0):
    """Principal-component regression: fraction of embedding variance explained
    by a covariate design matrix.

    This is scIB's PCR machinery (Luecken et al. 2022, Methods) generalised from
    a single batch covariate to an arbitrary design, so the same estimator can
    be pointed at the *biological* covariate. Returns
    ``sum_k adjR2_k * Var_k / sum_k Var_k`` over the top ``n_pcs`` components.
    """
    Z, var = reduce_dims(X, n_pcs=n_pcs, seed=seed)
    r2 = np.array([_adj_r2(Z[:, k], design) for k in range(Z.shape[1])])
    return float((r2 * var).sum() / var.sum())


# ------------------------------------------------------------- biology metrics
def trajectory_variance(X, step, group, n_pcs=N_PCS_EVAL, seed=0):
    """Absolute fraction of embedding variance attributable to the true trajectory."""
    return pcr_covariate(X, _design_trajectory(step, group), n_pcs=n_pcs, seed=seed)


def step_variance(X, step, n_pcs=N_PCS_EVAL, seed=0):
    """Same, but pooling branches: variance explained by a quadratic in Step alone."""
    t = np.asarray(step, dtype=float)
    t = (t - t.min()) / max(float(np.ptp(t)), 1e-12)
    D = np.column_stack([np.ones_like(t), t, t**2])
    return pcr_covariate(X, D, n_pcs=n_pcs, seed=seed)


def dpt_from_embedding(X, root_idx, n_neighbors=15, seed=0):
    """Diffusion pseudotime computed *on a given embedding*.

    Wraps scanpy's DPT so that pseudotime can be recovered from any method's
    output and compared against the ground-truth Step.
    """
    import anndata as ad
    import scanpy as sc

    a = ad.AnnData(np.asarray(_dense(X), dtype=np.float32))
    a.obsm["X_emb"] = a.X.copy()
    sc.pp.neighbors(a, use_rep="X_emb", n_neighbors=n_neighbors, random_state=seed)
    sc.tl.diffmap(a, random_state=seed)
    a.uns["iroot"] = int(root_idx)
    sc.tl.dpt(a)
    return np.asarray(a.obs["dpt_pseudotime"].values, dtype=float)


def pseudotime_recovery(X, step, root_idx, seed=0):
    """|Spearman| between DPT on the embedding and the ground-truth Step.

    Absolute value because DPT's direction is arbitrary up to the root choice.
    """
    try:
        dpt = dpt_from_embedding(X, root_idx, seed=seed)
        ok = np.isfinite(dpt)
        if ok.sum() < 10:
            return np.nan
        r = spearmanr(dpt[ok], np.asarray(step)[ok]).statistic
        return float(abs(r)) if np.isfinite(r) else np.nan
    except Exception:
        return np.nan


def knn_overlap(X, X_ref, k=KNN_OVERLAP_K, seed=0):
    """Mean fraction of each cell's k nearest neighbours shared with a reference
    embedding. A geometry-based check on TVR that makes no variance assumption."""
    from sklearn.neighbors import NearestNeighbors

    def nbrs(M):
        M = _dense(M)
        nn = NearestNeighbors(n_neighbors=k + 1).fit(M)
        return nn.kneighbors(M, return_distance=False)[:, 1:]

    A, B = nbrs(X), nbrs(X_ref)
    return float(np.mean([len(set(a) & set(b)) / k for a, b in zip(A, B)]))


def local_trajectory_metrics(X, truetime, group, k=KNN_OVERLAP_K):
    """Ground-truth-referenced *local* structure, complementary to the
    variance-based TVR.

    ``local_time_disp`` - median |true-pseudotime difference| to a cell's k
    nearest neighbours in the embedding, scaled by the standard deviation of
    the true pseudotime. Small means neighbours in the embedding really are
    neighbours in developmental time. ``nbr_purity`` - fraction of those
    neighbours on the same branch.

    These make no variance assumption, so if they track TVR the headline number
    is not an artefact of the PCR estimator.
    """
    from sklearn.neighbors import NearestNeighbors

    Xd = _dense(X)
    tt = np.asarray(truetime, dtype=float)
    grp = np.asarray(group)
    nn = NearestNeighbors(n_neighbors=k + 1).fit(Xd)
    idx = nn.kneighbors(Xd, return_distance=False)[:, 1:]
    disp = np.median(np.abs(tt[idx] - tt[:, None]), axis=1).mean() / max(tt.std(), 1e-12)
    purity = float((grp[idx] == grp[:, None]).mean())
    return {"local_time_disp": float(disp), "nbr_purity": purity}


def branch_ari(X, group, seed=0):
    """ARI of k-means on the embedding against the true branch labels."""
    Z, _ = reduce_dims(X, seed=seed)
    n_clusters = len(np.unique(group))
    km = KMeans(n_clusters=n_clusters, n_init=10, random_state=seed).fit(Z)
    return float(adjusted_rand_score(np.asarray(group), km.labels_))


# --------------------------------------------------------------- batch metrics
def make_neighbors(X, n_neighbors=N_NEIGHBORS, seed=0):
    return pynndescent(_dense(X).astype(np.float32), n_neighbors=n_neighbors,
                       random_state=seed, n_jobs=4)


def batch_metrics(X, batch, label, neighbors=None, seed=0):
    """kBET acceptance rate, ASW-batch, graph iLISI, PCR-batch.

    All four are oriented so that **higher = better mixed** except ``pcr_batch``,
    which is the fraction of variance attributable to batch (lower = better).
    """
    batch = np.asarray(batch)
    label = np.asarray(label)
    nn = neighbors if neighbors is not None else make_neighbors(X, seed=seed)
    out = {}
    try:
        # scib_metrics.kbet returns (acceptance_rate, stat, pvalue) despite its
        # annotation saying float; take the acceptance rate.
        kb = scib_metrics.kbet(nn, batch)
        out["kbet"] = float(kb[0] if isinstance(kb, tuple) else kb)
    except Exception as e:
        out["kbet"] = np.nan
        out["kbet_err"] = str(e)[:120]
    try:
        out["asw_batch"] = float(scib_metrics.silhouette_batch(
            _dense(X).astype(np.float32), label, batch))
    except Exception:
        out["asw_batch"] = np.nan
    try:
        out["ilisi"] = float(scib_metrics.ilisi_knn(nn, batch))
    except Exception:
        out["ilisi"] = np.nan
    try:
        out["pcr_batch"] = pcr_covariate(X, _design_onehot(batch), seed=seed)
    except Exception:
        out["pcr_batch"] = np.nan
    return out
