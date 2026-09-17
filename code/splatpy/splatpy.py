"""
splatpy: a faithful NumPy port of the Splat generative model from splatter
(Zappia, Phipson & Oshlack, Genome Biology 2017; https://github.com/Oshlack/splatter).

Ported because no R runtime is available in this workspace. Every step below
mirrors the corresponding function in splatter's R/splat-simulate.R (referenced
in the docstrings) so the simulation stays citable and auditable.

Supports method="single" | "groups" | "paths" and multi-batch simulation with
splatter's log-normal gene-wise batch factors, including `batch_rm_effect`,
which reproduces splatter's `batch.rmEffect=TRUE`: with the same seed it draws
the identical random stream but sets all batch factors to 1, yielding an exactly
matched batch-free counterfactual of the same cells. That pairing is what makes
"how much biology was removed" measurable rather than inferred.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np
from scipy.interpolate import CubicSpline


# ---------------------------------------------------------------- parameters
@dataclass
class SplatParams:
    """Defaults copied from splatter's SplatParams prototype (R/AllClasses.R)."""

    n_genes: int = 10000
    batch_cells: Sequence[int] = (100,)
    batch_fac_loc: float | Sequence[float] = 0.1
    batch_fac_scale: float | Sequence[float] = 0.1
    batch_rm_effect: bool = False

    mean_rate: float = 0.3
    mean_shape: float = 0.6

    lib_loc: float = 11.0
    lib_scale: float = 0.2
    lib_norm: bool = False

    out_prob: float = 0.05
    out_fac_loc: float = 4.0
    out_fac_scale: float = 0.5

    group_prob: Sequence[float] = (1.0,)
    de_prob: float | Sequence[float] = 0.1
    de_down_prob: float | Sequence[float] = 0.5
    de_fac_loc: float | Sequence[float] = 0.1
    de_fac_scale: float | Sequence[float] = 0.4

    bcv_common: float = 0.1
    bcv_df: float = 60.0

    dropout_type: str = "none"  # none | experiment | batch | group
    dropout_mid: float | Sequence[float] = 0.0
    dropout_shape: float | Sequence[float] = -1.0

    path_from: Sequence[int] = (0,)
    path_n_steps: int | Sequence[int] = 100
    path_skew: float | Sequence[float] = 0.5
    path_nonlinear_prob: float = 0.1
    path_sigma_fac: float = 0.8

    seed: int | None = None

    @property
    def n_groups(self) -> int:
        return len(self.group_prob)

    @property
    def n_batches(self) -> int:
        return len(self.batch_cells)

    @property
    def n_cells(self) -> int:
        return int(np.sum(self.batch_cells))


def _vec(x, n, dtype=float):
    a = np.asarray(x, dtype=dtype)
    return np.repeat(a, n) if a.ndim == 0 else a


# ------------------------------------------------------------------ helpers
def get_lnorm_factors(rng, n_facs, sel_prob, neg_prob, fac_loc, fac_scale):
    """splatter::getLNormFactors (R/splat-simulate.R:924)."""
    is_selected = rng.binomial(1, sel_prob, n_facs).astype(bool)
    n_selected = int(is_selected.sum())
    dir_selected = (-1.0) ** rng.binomial(1, neg_prob, n_selected)
    facs_selected = rng.lognormal(fac_loc, fac_scale, n_selected)
    # splatter reverses direction for factors < 1 so up/down stay balanced
    dir_selected[facs_selected < 1] = -1 * dir_selected[facs_selected < 1]
    factors = np.ones(n_facs)
    factors[is_selected] = facs_selected**dir_selected
    return factors


def bridge(rng, x, y, N=5, n=100, sigma_fac=0.8):
    """splatter's `bridge` (R/splat-simulate.R:993): smoothed Brownian bridge."""
    dt = 1.0 / (N - 1)
    t = np.linspace(0, 1, N)
    sigma2 = rng.uniform(0, sigma_fac * np.mean([x, y]))
    X = np.concatenate([[0.0], np.cumsum(rng.normal(0, sigma2, N - 1) * np.sqrt(dt))])
    BB = x + X - t * (X[N - 1] - y + x)
    # R's spline(..., method="natural") on equally spaced points
    spl = CubicSpline(np.arange(N), BB, bc_type="natural")
    out = spl(np.linspace(0, N - 1, n))
    out[out < 0] = 1e-6
    return out


def build_bridges(rng, xs, ys, n, sigma_facs):
    """Vectorised `bridge` over genes -> (n_genes, n) matrix."""
    return np.stack(
        [bridge(rng, xs[i], ys[i], N=5, n=n, sigma_fac=sigma_facs[i]) for i in range(len(xs))]
    )


def _logistic(x, x0, k):
    return 1.0 / (1.0 + np.exp(-k * (x - x0)))


def _path_order(path_from):
    """splatter::getPathOrder - topological order so parents precede children."""
    path_from = list(path_from)
    order, done = [], {0}
    remaining = list(range(1, len(path_from) + 1))
    while remaining:
        progressed = False
        for p in list(remaining):
            if path_from[p - 1] in done:
                order.append(p)
                done.add(p)
                remaining.remove(p)
                progressed = True
        if not progressed:
            raise ValueError("path_from contains a cycle or an unreachable node")
    return order


# ---------------------------------------------------------------- simulator
def splat_simulate(params: SplatParams, method: str = "paths", verbose: bool = False):
    """Port of splatter::splatSimulate. Returns a dict of counts + ground truth.

    The returned dict mirrors splatter's SingleCellExperiment slots:
      counts (cells x genes), Batch, Group, Step, ExpLibSize,
      GeneMean, BatchFac<b>, DEFacPath<p>/DEFacGroup<g>, TrueCounts, CellMeans.
    """
    assert method in ("single", "groups", "paths")
    p = params
    rng = np.random.default_rng(p.seed)
    n_genes, n_cells = p.n_genes, p.n_cells

    # --- cell metadata: batch assignment (splatSimulate setup)
    batch = np.concatenate([np.full(nc, i + 1) for i, nc in enumerate(p.batch_cells)])

    # --- 1. library sizes (splatSimLibSizes)
    if p.lib_norm:
        exp_lib = rng.normal(p.lib_loc, p.lib_scale, n_cells)
        exp_lib[exp_lib < 0] = exp_lib[exp_lib > 0].min() / 2
    else:
        exp_lib = rng.lognormal(p.lib_loc, p.lib_scale, n_cells)

    # --- 2. gene means + expression outliers (splatSimGeneMeans)
    base_means_gene = rng.gamma(shape=p.mean_shape, scale=1.0 / p.mean_rate, size=n_genes)
    outlier_facs = get_lnorm_factors(rng, n_genes, p.out_prob, 0, p.out_fac_loc, p.out_fac_scale)
    median_mean = np.median(base_means_gene)
    means_gene = base_means_gene.copy()
    is_outlier = outlier_facs != 1
    means_gene[is_outlier] = (median_mean * outlier_facs)[is_outlier]

    # --- 3. batch effects (splatSimBatchEffects)
    bfl = _vec(p.batch_fac_loc, p.n_batches)
    bfs = _vec(p.batch_fac_scale, p.n_batches)
    batch_facs = np.stack(
        [get_lnorm_factors(rng, n_genes, 1, 0.5, bfl[i], bfs[i]) for i in range(p.n_batches)]
    )  # (n_batches, n_genes); RNG is consumed even when the effect is removed
    if p.batch_rm_effect:
        batch_facs = np.ones_like(batch_facs)

    # (splatSimBatchCellMeans)
    batch_means_cell = batch_facs[batch - 1] * means_gene  # (n_cells, n_genes)

    group = np.ones(n_cells, dtype=int)
    step = np.zeros(n_cells, dtype=int)
    de_facs = None

    if method == "single":
        cell_facs = np.ones((n_cells, n_genes))
    else:
        # assign cells to groups/paths (splatSimulate)
        group = rng.choice(np.arange(1, p.n_groups + 1), size=n_cells, p=np.asarray(p.group_prob))
        de_prob = _vec(p.de_prob, p.n_groups)
        de_down = _vec(p.de_down_prob, p.n_groups)
        de_loc = _vec(p.de_fac_loc, p.n_groups)
        de_scale = _vec(p.de_fac_scale, p.n_groups)

        if method == "groups":
            # splatSimGroupDE + splatSimGroupCellMeans
            de_facs = np.stack(
                [
                    get_lnorm_factors(rng, n_genes, de_prob[i], de_down[i], de_loc[i], de_scale[i])
                    for i in range(p.n_groups)
                ]
            )
            cell_facs = de_facs[group - 1]
        else:
            # splatSimPathDE: DE factors compose along the path tree
            de_facs = np.ones((p.n_groups, n_genes))
            for path in _path_order(p.path_from):
                frm = p.path_from[path - 1]
                facs_from = np.ones(n_genes) if frm == 0 else de_facs[frm - 1]
                local = get_lnorm_factors(
                    rng, n_genes, de_prob[path - 1], de_down[path - 1],
                    de_loc[path - 1], de_scale[path - 1],
                )
                de_facs[path - 1] = local * facs_from

            # splatSimPathCellMeans
            n_steps = _vec(p.path_n_steps, p.n_groups, dtype=int)
            skew = _vec(p.path_skew, p.n_groups)
            sigma_facs = np.stack(
                [
                    np.where(rng.binomial(1, p.path_nonlinear_prob, n_genes).astype(bool),
                             p.path_sigma_fac, 0.0)
                    for _ in range(p.n_groups)
                ]
            )
            path_steps = []
            for idx in range(p.n_groups):
                frm = p.path_from[idx]
                facs_start = np.ones(n_genes) if frm == 0 else de_facs[frm - 1]
                path_steps.append(
                    build_bridges(rng, facs_start, de_facs[idx], int(n_steps[idx]), sigma_facs[idx])
                )  # (n_genes, n_steps)

            path_probs = []
            for idx in range(p.n_groups):
                pr = np.linspace(skew[idx], 1 - skew[idx], int(n_steps[idx]))
                path_probs.append(pr / pr.sum())

            step = np.array(
                [rng.choice(int(n_steps[g - 1]), p=path_probs[g - 1]) + 1 for g in group]
            )
            cell_facs = np.stack(
                [path_steps[group[i] - 1][:, step[i] - 1] for i in range(n_cells)]
            )

    # (splatSim*CellMeans) - renormalise to expected library size
    cell_means_gene = batch_means_cell * cell_facs
    cell_props = cell_means_gene / cell_means_gene.sum(axis=1, keepdims=True)
    base_means_cell = cell_props * exp_lib[:, None]

    # --- 4. BCV-adjusted means (splatSimBCVMeans)
    if np.isfinite(p.bcv_df):
        bcv = (p.bcv_common + 1.0 / np.sqrt(base_means_cell)) * np.sqrt(
            p.bcv_df / rng.chisquare(p.bcv_df, n_genes)
        )
    else:
        bcv = p.bcv_common + 1.0 / np.sqrt(base_means_cell)
    means_cell = rng.gamma(shape=1.0 / bcv**2, scale=base_means_cell * bcv**2)

    # --- 5. counts (splatSimTrueCounts)
    true_counts = rng.poisson(means_cell)

    # --- 6. dropout (splatSimDropout)
    counts = true_counts
    dropout = None
    if p.dropout_type != "none":
        if p.dropout_type == "experiment":
            mid = np.full(n_cells, float(np.asarray(p.dropout_mid).ravel()[0]))
            shape = np.full(n_cells, float(np.asarray(p.dropout_shape).ravel()[0]))
        elif p.dropout_type == "batch":
            mid = _vec(p.dropout_mid, p.n_batches)[batch - 1]
            shape = _vec(p.dropout_shape, p.n_batches)[batch - 1]
        elif p.dropout_type == "group":
            mid = _vec(p.dropout_mid, p.n_groups)[group - 1]
            shape = _vec(p.dropout_shape, p.n_groups)[group - 1]
        else:
            raise ValueError(f"unsupported dropout_type {p.dropout_type}")
        eta = np.log(means_cell)
        drop_prob = _logistic(eta, mid[:, None], shape[:, None])
        keep = rng.binomial(1, 1 - drop_prob)
        counts = true_counts * keep
        dropout = keep == 0

    return {
        "counts": counts.astype(np.int32),
        "true_counts": true_counts.astype(np.int32),
        "cell_means": means_cell,
        "dropout": dropout,
        "batch": batch,
        "group": group,
        "step": step,
        "exp_lib_size": exp_lib,
        "gene_mean": means_gene,
        "batch_facs": batch_facs,
        "de_facs": de_facs,
        "params": p,
    }


def to_anndata(sim):
    """Wrap a splat_simulate() result as AnnData with ground truth in obs/var."""
    import anndata as ad
    import pandas as pd

    n_cells, n_genes = sim["counts"].shape
    obs = pd.DataFrame(
        {
            "Batch": pd.Categorical([f"Batch{b}" for b in sim["batch"]]),
            "Group": pd.Categorical([f"Path{g}" for g in sim["group"]]),
            "Step": sim["step"],
            "ExpLibSize": sim["exp_lib_size"],
        },
        index=[f"Cell{i + 1}" for i in range(n_cells)],
    )
    var = pd.DataFrame({"GeneMean": sim["gene_mean"]},
                       index=[f"Gene{i + 1}" for i in range(n_genes)])
    for b in range(sim["batch_facs"].shape[0]):
        var[f"BatchFac{b + 1}"] = sim["batch_facs"][b]
    if sim["de_facs"] is not None:
        for g in range(sim["de_facs"].shape[0]):
            var[f"DEFac{g + 1}"] = sim["de_facs"][g]
    a = ad.AnnData(X=sim["counts"].astype(np.float32), obs=obs, var=var)
    a.layers["true_counts"] = sim["true_counts"].astype(np.float32)
    return a


def paired_simulation(params: SplatParams, method: str = "paths"):
    """Return (with_batch, without_batch) - identical cells, batch effect toggled.

    Both calls consume the same RNG stream (splatter's batch.rmEffect trick), so
    the two matrices differ *only* by the gene-wise batch factors. The no-batch
    arm is the ground-truth target any integration method should recover.
    """
    import dataclasses

    with_batch = splat_simulate(params, method=method)
    p_off = dataclasses.replace(params, batch_rm_effect=True)
    without_batch = splat_simulate(p_off, method=method)
    return with_batch, without_batch
