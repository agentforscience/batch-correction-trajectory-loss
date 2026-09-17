# Code Repositories

Everything below is either cloned from upstream or written here. Each entry says
what it is, whether it **runs in this environment**, and what it is for.

**Environment constraint that shapes this whole directory:** the container has
no R runtime, no root (`apt-get` fails), and conda is disallowed. Every R-only
tool below is documented but unusable. Python 3.12 + CUDA 13.0 on an RTX A6000.

---

## Written here

### `splatpy/` — NumPy port of the Splat generative model

`splatpy/splatpy.py`. Written because `splatter` is R-only. Follows
`code/splatter/R/splat-simulate.R` function by function; each Python function
names the R function it mirrors.

```python
import sys; sys.path.insert(0, "code/splatpy")
from splatpy import SplatParams, splat_simulate, paired_simulation, to_anndata

p = SplatParams(
    n_genes=2000, batch_cells=(600, 600),
    batch_fac_loc=0.2, batch_fac_scale=0.2,     # batch effect strength
    group_prob=(0.4, 0.3, 0.3), path_from=(0, 1, 1),  # branching trajectory
    de_prob=0.3, de_fac_loc=0.4, path_n_steps=100,
    dropout_type="experiment", dropout_mid=1.5, dropout_shape=-1.0,
    lib_loc=9.5, seed=1,
)
with_batch, no_batch = paired_simulation(p, method="paths")
adata = to_anndata(with_batch)   # obs: Batch, Group, Step, ExpLibSize
```

**Why this exists and not just any simulator.** `paired_simulation` calls the
model twice with the same seed, once with `batch_rm_effect=False` and once with
`True`. Splatter draws the batch factors either way and only then overwrites
them with 1, so the RNG stream is identical — the two outputs are the *same
cells* (verified: identical `Group`, `Step`, `ExpLibSize`) differing only by the
gene-wise batch factors. The batch-free arm is therefore the exact
counterfactual: whatever an integration method fails to recover, relative to it,
is biology it destroyed. No published benchmark has this axis (see
`papers/README.md` on scIB's trajectory metric).

Ground truth exposed: `Step` (position along the path = true pseudotime),
`Group` (which branch), `batch_facs` (per-gene, per-batch multiplier),
`de_facs` (per-gene, per-path DE factor), plus `true_counts` pre-dropout.

Verified: 2000 genes × 1200 cells ≈ 4 s. At `batch_fac_loc=0`, ASW-batch = 0.000
and Spearman(PC1, Step) = 0.86 on the main path; at `batch_fac_loc=0.4`,
ASW-batch = 0.408 and the same correlation collapses to 0.09. The knob works and
spans the interesting range.

**Not yet validated against real R splatter output.** The port was checked for
internal consistency and plausible marginals, not against `splatSimulate()`
row-for-row, because there is no R here. Read
`papers/crowell2023_shaky_foundations_simulation.pdf` before making realism
claims from it.

---

## Metrics

### `scib/` — github.com/theislab/scib (`cd67913`)
The reference implementation of every metric in Luecken et al. 2022.
**Installed and importable** (`scib` 1.1.7, via PyPI).

- `scib.metrics.trajectory_conservation(adata_pre, adata_post, label_key, pseudotime_key='dpt_pseudotime', batch_key=None)` — **verified importable**.
  This is the published trajectory metric. Note its semantics: it compares
  post-integration DPT to *pre-integration* DPT, so on real data it measures
  agreement with a corrupted reference, not with truth.
- Also: `nmi`, `ari`, `silhouette`, `silhouette_batch`, `graph_connectivity`,
  `isolated_labels`, `pcr_comparison`, `hvg_overlap`, `cell_cycle`,
  `ilisi_graph`, `clisi_graph`, `morans_i`.
- **`scib.metrics.kBET` needs rpy2 + the R kBET package — unusable here.**
  Use `scib_metrics.kbet` instead (below).

### `scib-metrics/` — github.com/YosefLab/scib-metrics (`792dbae`)
JAX reimplementation. **Installed and verified** (0.6.0). Pure Python, so it is
the kBET source for this project.

Available: `kbet`, `kbet_per_label`, `silhouette_batch` (= ASW-batch),
`silhouette_label`, `ilisi_knn`, `clisi_knn`, `lisi_knn`, `graph_connectivity`,
`isolated_labels`, `nmi_ari_cluster_labels_kmeans`, `nmi_ari_cluster_labels_leiden`,
`pcr_comparison`, `bras`, `sbee`, `nearest_neighbors.pynndescent`.

```python
import scib_metrics
nn = scib_metrics.nearest_neighbors.pynndescent(X, n_neighbors=25)
acceptance_rate, stat, pvalue = scib_metrics.kbet(nn, batch_labels)
```

**Gotcha:** jax here is CPU-only ("An NVIDIA GPU may be present… CUDA-enabled
jaxlib is not installed"). Metrics still run, just slower. `uv add
"jax[cuda12]"` would fix it if metric time becomes the bottleneck.

### `kBET/` — github.com/theislab/kBET (`afc5f43`)
Original R implementation. **R-only, unusable.** Kept as the definitional
reference for the χ² test on neighbourhood batch composition.

### `CellMixS/` — github.com/almutlue/CellMixS (`3dc033c`)
R/Bioconductor. **Unusable.** Its per-cell mixing score (CMS) is a good idea —
it localises *where* mixing failed rather than reporting one rejection rate —
and would need reimplementing in Python if wanted.

---

## Integration methods — verified working

Smoke-tested end to end on `splatpy` output (800 cells × 1200 genes) unless noted.

| Method | Package | Status | Output | Notes |
|---|---|---|---|---|
| ComBat | `scanpy.pp.combat` | **OK, 0.1 s** | corrected matrix | |
| BBKNN | `bbknn` 1.6.0 | **OK, 11 s** | kNN graph only | graph output limits which metrics apply |
| Scanorama | `scanorama` 1.7.4 | **OK, 0.3 s** | matrix + embedding | evaluate both separately, as scIB does |
| Harmony | `harmonypy` 2.0.0 | **OK, <0.1 s** | embedding | see API gotcha below |
| scVI | `scvi-tools` 1.5.0 | **OK, 8 s (GPU)** | embedding | needs raw counts |
| scANVI | `scvi-tools` 1.5.0 | **OK, 1.3 s (GPU)** | embedding | needs raw counts + labels |
| scGen | `code/scgen` (patched) | **OK, 1.5 s (GPU)** | corrected matrix + latent | needs labels; two patches, below |

**harmonypy 2.0 API change.** `run_harmony(...).Z_corr` is already
`(n_cells, n_dims)`. Older code (including several tutorials and the scIB
pipeline) does `.Z_corr.T`, which silently produces the wrong shape.
Do **not** transpose.

**Harmony convergence caveat, observed here.** On `splatpy` data with
`batch_fac_loc >= 0.2`, harmonypy logged "Converged after 1 iteration" and left
ASW-batch unchanged (0.207 → 0.202; 0.408 → 0.408) — it did essentially nothing.
At `batch_fac_loc = 0.1` it worked normally (0.076 → 0.002, kBET 0.101 → 0.954).
The likely cause is that with strongly separated batches each soft-k-means
cluster is captured by a single batch, so the diversity penalty has no gradient.
Before reporting "Harmony fails on strong batch effects", check `theta`,
`nclust` (default 40 is large for ~1000 cells), and `max_iter_harmony`.
This could be a real finding or a misconfiguration; it must be resolved, not
reported as-is.

### `scgen/` — github.com/theislab/scgen (`d79e1f0`, v2.1.1)
Installed editable from this clone. Upstream is unmaintained against
scvi-tools 1.5 and needed **two patches, both applied and verified**:

1. `scgen/_scgenvae.py`, `inference()` — scvi-tools ≥1.2 reads `qz`/`qzm`/`qzv`
   from the inference dict (`VAEMixin.get_latent_representation`), but scGen only
   emitted `qz_m`/`qz_v`, so `batch_removal()` crashed with
   `'NoneType' object has no attribute 'sqrt'`. Now emits both plus a
   `Normal(qz_m, qz_v.sqrt())`.
2. `scgen/_scgen.py`, `batch_removal()` — three `AnnData.concatenate(...)` calls,
   removed in anndata ≥0.11, replaced with `anndata.concat(..., label=...)`.

Also required `adjusttext` (an undeclared runtime import).
The PyPI `scgen` 2.1.0 release is **not** usable — it imports `scvi._compat`,
which no longer exists. Install from this clone.

### `scanorama/` — github.com/brianhie/scanorama (`4b4ea07`)
Cloned for reference; the PyPI package is what actually runs.

### `harmony-R/` — github.com/immunogenomics/harmony (`df19af2`)
R implementation. **Unusable.** Reference only; `harmonypy` is the port used.

### `scarches/` — github.com/theislab/scarches (`d35cde6`)
Cloned but **not installed and not smoke-tested.** Relevant if reference-mapping
(scPoli / scANVI surgery) becomes part of the design; otherwise ignore.

## Integration methods — NOT available

| Method | Blocker |
|---|---|
| MNN / fastMNN | `mnnpy` 0.1.9.5 fails to build on Python 3.12 (needs `longintrepr.h`, removed in 3.11). `scanpy.external.pp.mnn_correct` therefore raises `ModuleNotFoundError`. Scanorama and BBKNN are the MNN-family stand-ins. |
| LIGER | `pyliger` 0.2.4 pulls `louvain` 0.7.1, which fails to build. |
| Seurat v3 CCA/RPCA | R-only. |
| trVAE, DESC, SAUCIE, Conos | Not attempted — none is a top performer in scIB and each adds install risk. |

Losing MNN and LIGER matters: scIB places LIGER at the batch-removal extreme of
the trade-off, so the frontier measured here will under-sample that end. State
this as a limitation rather than presenting the frontier as complete.

---

## Prior benchmarks and designs to build on

### `cellanova/` — github.com/Janezjz/cellanova (`01ff1eb`)
**Installed and importable** (from the bundled wheel — `pip install -e .` fails
because the repo is flat-layout with `data/` and `figures/` at top level, which
confuses setuptools; use `code/cellanova/dist/cellanova-0.1.0-py3-none-any.whl`).
Exposes `cellanova.model`, `cellanova.utils`.

Two things to take from here beyond the method itself:
- `tutorials/` — worked notebooks for the global-distortion and gene-level
  distortion measurements described in `papers/README.md`.
- `data/t1d_example.h5ad` — **bonus dataset, 20,000 cells × 26,163 genes**,
  type-1-diabetes pancreatic islets. `obs` has `donor_id`, `disease_state`,
  `disease`, `assay`, `celltype`, and a `PseudoState` column. Real donor batches
  with a designated control pool, i.e. ready for the hold-out-control spike-in.
  No count layer, no obsm.
- `data/Simu1_SeuratObj.rds` — their simulation, R format, unreadable here.

### `scib-pipeline/` — github.com/theislab/scib-pipeline (`e97631a`)
Snakemake pipeline that reproduces Luecken 2022. Its conda env files pin the
exact method versions the published numbers came from. Useful as a spec for how
each method was invoked; **do not try to run it** (needs conda + R).

### `Iniquitate/` — github.com/hsmaan/Iniquitate (`cb20fe1`)
Snakemake workflow behind Maan et al., "the differential impacts of dataset
imbalance in single-cell data integration" (Nat Biotechnol 2024; PDF
unobtainable, abstract in `papers/abstracts_only/`). Also from the same author:
`balanced-clustering`, reworked clustering metrics for imbalanced settings —
worth pulling if cell-type imbalance is swept as a covariate.

### `batchbench/` — github.com/cellgeni/batchbench (`b458380`)
Nextflow + R pipeline from Chazarra-Gil et al. 2021. **Unusable as a pipeline.**
Read it for the entropy-of-mixing metric definitions.

### `Palantir/` — github.com/dpeerlab/Palantir (`db77fbb`)
Python, actively maintained. **Not yet installed.** This is the source of the
`palantir_pseudotime` and `palantir_branch_probs` annotations in
`datasets/palantir_cd34/`. Install (`uv add palantir`) if branch probabilities
need recomputing on *integrated* embeddings — which is exactly what measuring
trajectory destruction on that dataset requires.

## Simulators — reference only (all R)

- `splatter/` — `Oshlack/splatter` (`eeb1dfc`). The source `splatpy` was ported
  from. Key files: `R/splat-simulate.R`, `R/AllClasses.R` (default parameters).
- `dyngen/` — `dynverse/dyngen` (`620b3e8`). GRN-driven, gold-standard
  trajectories. Would be the better simulator if R were available.
