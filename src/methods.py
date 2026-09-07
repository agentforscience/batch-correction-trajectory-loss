"""
Integration-method wrappers.

Every wrapper has the same contract:

    run_<method>(adata, batch_key, label_key, seed, **kw) -> (X_emb, meta)

``adata.X`` is log1p-normalised expression, ``adata.layers['counts']`` holds raw
counts (needed by scVI/scANVI), and ``adata.obs[batch_key]`` / ``[label_key]``
hold the batch and biological-label annotations.

Methods whose native output is a corrected *matrix* (ComBat, Scanorama-matrix,
scGen) are reduced with PCA so that everything downstream sees an embedding of
the same kind; ``N_PCS`` is fixed across methods. BBKNN outputs only a kNN graph,
so its embedding is a diffusion map of that graph - a substitution that is
flagged in the results table.

Invocation conventions follow ``code/scib-pipeline`` (the pipeline that produced
the published scIB numbers) where they apply.
"""
from __future__ import annotations

import time
import warnings

import numpy as np
import scipy.sparse as sp

warnings.filterwarnings("ignore")

import logging
for _n in ("harmonypy", "scvi", "lightning", "pytorch_lightning", "scanpy", "jax"):
    logging.getLogger(_n).setLevel(logging.ERROR)

N_PCS = 20          # PCA rank used to build embeddings from corrected matrices
N_LATENT = 10       # latent dimensionality for the deep methods


def _dense(X):
    return np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)


def _pca(X, n_comps=N_PCS, seed=0):
    from sklearn.decomposition import PCA
    X = _dense(X).astype(np.float64)
    X = X - X.mean(axis=0, keepdims=True)
    n = int(min(n_comps, X.shape[1], X.shape[0] - 1))
    return PCA(n_components=n, random_state=seed).fit_transform(X)


def _timed(fn):
    def wrapper(*a, **kw):
        t0 = time.time()
        emb, meta = fn(*a, **kw)
        meta = dict(meta or {})
        meta["seconds"] = round(time.time() - t0, 2)
        return np.asarray(emb, dtype=np.float64), meta
    wrapper.__name__ = fn.__name__
    return wrapper


# ------------------------------------------------------------------ references
@_timed
def run_uncorrected(adata, batch_key, label_key, seed=0, **kw):
    """No correction: PCA of the log-normalised expression matrix."""
    return _pca(adata.X, seed=seed), {}


# ------------------------------------------------------------------ correctors
@_timed
def run_combat(adata, batch_key, label_key, seed=0, **kw):
    """Empirical-Bayes location/scale adjustment per gene (Johnson et al. 2007)."""
    import scanpy as sc
    a = adata.copy()
    a.X = _dense(a.X).astype(np.float64)
    sc.pp.combat(a, key=batch_key)
    return _pca(a.X, seed=seed), {}


@_timed
def run_harmony(adata, batch_key, label_key, seed=0, theta=2.0, nclust=None,
                max_iter_harmony=20, **kw):
    """Harmony: iterative soft-k-means clustering + linear batch correction in PC space.

    **PC standardisation is required here.** harmonypy's soft-k-means uses a
    fixed kernel width (``sigma``) on cosine distances. On raw (unstandardised)
    PCs, where PC1-2 carry an order of magnitude more variance than PC20, the
    assignments degenerate to one cluster per batch, the ridge design inside
    each cluster becomes rank-deficient, and Harmony returns its input unchanged
    ("Converged after 1 iteration"). This reproduced at every
    ``theta`` x ``nclust`` combination tried (``results/harmony_diagnostic.csv``)
    and is *not* a property of the method: dividing each PC by its standard
    deviation before the call restores normal behaviour (kBET 0.00 -> 0.96 at
    ``batch_fac_loc`` = 0.4). Harmony's correction is an additive offset in the
    supplied basis, so multiplying the columns back by the same standard
    deviations returns the result to the original PC scale and keeps the
    embedding comparable with every other method's.
    """
    import harmonypy
    pcs = _pca(adata.X, seed=seed)
    s = pcs.std(axis=0, keepdims=True)
    s[s == 0] = 1.0
    meta_df = adata.obs[[batch_key]].copy()
    if nclust is None:
        nclust = int(np.clip(adata.n_obs // 30, 5, 50))
    ho = harmonypy.run_harmony(
        pcs / s, meta_df, [batch_key], theta=theta, nclust=nclust,
        max_iter_harmony=max_iter_harmony, random_state=seed,
    )
    # harmonypy >= 2.0 returns Z_corr as (n_cells, n_dims) - do NOT transpose.
    Z = np.asarray(ho.Z_corr)
    if Z.shape[0] != adata.n_obs:
        Z = Z.T
    return Z * s, {"theta": theta, "nclust": nclust}


_SCANORAMA_CACHE = {}


def _scanorama(adata, batch_key, seed):
    """Run Scanorama once and cache both outputs.

    ``run_scanorama_embed`` and ``run_scanorama_matrix`` are two read-outs of
    one fit; the cache stops the expensive panorama search from running twice
    per unit. Keyed on data content, so it cannot leak across units.
    """
    key = (adata.n_obs, adata.n_vars, float(np.asarray(_dense(adata.X)).sum()), seed, batch_key)
    if key in _SCANORAMA_CACHE:
        return _SCANORAMA_CACHE[key]
    import scanorama
    a = adata.copy()
    a.X = _dense(a.X).astype(np.float32)
    order = np.argsort(a.obs[batch_key].astype(str).values, kind="stable")
    a = a[order].copy()
    splits = [a[a.obs[batch_key].astype(str) == b].copy()
              for b in a.obs[batch_key].astype(str).unique()]
    corrected = scanorama.correct_scanpy(splits, return_dimred=True, dimred=N_PCS)
    emb = np.vstack([_dense(c.obsm["X_scanorama"]) for c in corrected])
    mat = np.vstack([_dense(c.X) for c in corrected])
    inv = np.argsort(order)
    _SCANORAMA_CACHE.clear()
    _SCANORAMA_CACHE[key] = (emb[inv], mat[inv])
    return _SCANORAMA_CACHE[key]


@_timed
def run_scanorama_embed(adata, batch_key, label_key, seed=0, **kw):
    """Scanorama's low-dimensional panorama embedding (mutual-nearest-neighbour family)."""
    emb, _ = _scanorama(adata, batch_key, seed)
    return emb, {}


@_timed
def run_scanorama_matrix(adata, batch_key, label_key, seed=0, **kw):
    """Scanorama's corrected expression matrix, then PCA. scIB scores both variants."""
    _, mat = _scanorama(adata, batch_key, seed)
    return _pca(mat, seed=seed), {}


@_timed
def run_bbknn(adata, batch_key, label_key, seed=0, **kw):
    """BBKNN: batch-balanced kNN graph. Embedding = diffusion map of that graph."""
    import scanpy as sc
    import bbknn
    a = adata.copy()
    a.obsm["X_pca"] = _pca(a.X, seed=seed)
    n_per = int(a.obs[batch_key].value_counts().min())
    bbknn.bbknn(a, batch_key=batch_key, neighbors_within_batch=min(3, max(1, n_per // 10)))
    sc.tl.diffmap(a, n_comps=min(15, N_PCS), random_state=seed)
    # component 0 of a diffusion map is the constant vector - drop it
    return np.asarray(a.obsm["X_diffmap"])[:, 1:], {"graph_output": True}


@_timed
def run_scvi(adata, batch_key, label_key, seed=0, max_epochs=400, n_latent=N_LATENT,
             lr=1e-2, **kw):
    """scVI: conditional VAE on raw counts with batch as a covariate.

    ``max_epochs=400`` at ``lr=1e-2`` was chosen from an explicit convergence
    diagnostic (``src/diagnostics``, ``results/scvi_epoch_diagnostic.csv``):
    scvi-tools' defaults (200 epochs, lr 1e-3) leave the batch effect almost
    untouched here (kBET 0.001), which would have been misread as a method
    failure. Batch mixing keeps improving up to ~15k gradient steps while the
    trajectory-variance score stays flat at 0.45-0.52, so the biology loss
    reported for scVI is not an under-training artefact.
    """
    import scvi
    import torch
    scvi.settings.seed = seed
    a = adata.copy()
    a.X = a.layers["counts"].copy()
    scvi.model.SCVI.setup_anndata(a, batch_key=batch_key)
    m = scvi.model.SCVI(a, n_latent=n_latent, gene_likelihood="nb")
    m.train(max_epochs=max_epochs, early_stopping=False,
            accelerator="gpu" if torch.cuda.is_available() else "cpu",
            plan_kwargs={"lr": lr}, enable_progress_bar=False)
    return m.get_latent_representation(), {"n_latent": n_latent, "max_epochs": max_epochs}


@_timed
def run_scanvi(adata, batch_key, label_key, seed=0, max_epochs=400, n_latent=N_LATENT,
               lr=1e-2, **kw):
    """scANVI: scVI plus cell-label supervision (an oracle-ish advantage, noted in the report)."""
    import scvi
    import torch
    scvi.settings.seed = seed
    a = adata.copy()
    a.X = a.layers["counts"].copy()
    a.obs["_label"] = a.obs[label_key].astype(str).values
    acc = "gpu" if torch.cuda.is_available() else "cpu"
    scvi.model.SCVI.setup_anndata(a, batch_key=batch_key, labels_key="_label")
    base = scvi.model.SCVI(a, n_latent=n_latent, gene_likelihood="nb")
    base.train(max_epochs=max_epochs, early_stopping=False, accelerator=acc,
               plan_kwargs={"lr": lr}, enable_progress_bar=False)
    m = scvi.model.SCANVI.from_scvi_model(base, unlabeled_category="Unknown")
    m.train(max_epochs=max(20, max_epochs // 4), accelerator=acc, enable_progress_bar=False)
    return m.get_latent_representation(), {"n_latent": n_latent}


@_timed
def run_scgen(adata, batch_key, label_key, seed=0, max_epochs=100, **kw):
    """scGen: VAE + per-cell-type latent-space vector arithmetic. Label-supervised."""
    import scgen
    import scvi
    import torch
    scvi.settings.seed = seed
    a = adata.copy()
    a.X = _dense(a.X).astype(np.float32)
    a.obs["cell_type"] = a.obs[label_key].astype(str).values
    a.obs["batch"] = a.obs[batch_key].astype(str).values
    scgen.SCGEN.setup_anndata(a, batch_key="batch", labels_key="cell_type")
    m = scgen.SCGEN(a)
    m.train(max_epochs=max_epochs, batch_size=64, early_stopping=True,
            early_stopping_patience=25, enable_progress_bar=False,
            accelerator="gpu" if torch.cuda.is_available() else "cpu")
    corrected = m.batch_removal()
    return _pca(corrected.X, seed=seed), {"max_epochs": max_epochs}


METHODS = {
    "uncorrected": run_uncorrected,
    "ComBat": run_combat,
    "Harmony": run_harmony,
    "Scanorama-embed": run_scanorama_embed,
    "Scanorama-matrix": run_scanorama_matrix,
    "BBKNN": run_bbknn,
    "scVI": run_scvi,
    "scANVI": run_scanvi,
    "scGen": run_scgen,
}

#: Methods that use the biological labels during correction. They get an
#: advantage no unsupervised method has, so they are marked in every table.
SUPERVISED = {"scANVI", "scGen"}
