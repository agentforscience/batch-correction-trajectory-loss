# Resources Catalog

Everything gathered for *How Much Biology Do Batch Correction Methods Remove?*
All paths are relative to the workspace root.

- **Papers**: 25 PDFs + 4 abstract-only records → `papers/`, detail in `papers/README.md`
- **Datasets**: 8 usable `.h5ad` files, 2.9 GB → `datasets/`, detail in `datasets/README.md`
- **Code**: 15 cloned repositories + 1 written module, 1.4 GB → `code/`, detail in `code/README.md`
- **Synthesis**: `literature_review.md`; direction budget in `planning.md`

---

## Environment

Fresh `uv` venv at `.venv`, isolated by a `[tool.uv] package = false`
`pyproject.toml` so nothing touches a parent environment.

```bash
source .venv/bin/activate
```

Hardware: NVIDIA RTX A6000 (49 GB, **~42 GB already in use by another process** —
budget around 6 GB), 32 CPUs, 503 GB RAM, 737 GB free disk.

Key versions (verified installed): python 3.12.8, numpy 2.2.6, scipy 1.18.1,
scanpy 1.12.4, anndata 0.13.3, scvelo 0.3.4, torch 2.14.0+cu130 (**CUDA
available: True**), scvi-tools 1.5.0, scib 1.1.7, scib-metrics 0.6.0,
harmonypy 2.0.0, scanorama 1.7.4, bbknn 1.6.0, scgen 2.1.1 (editable, patched),
cellanova 0.1.0, jax 0.11.1 (**CPU only**).

**No R, no root, conda disallowed.** This is the single constraint that shaped
the most decisions — see `planning.md` R3 and `code/README.md`.

---

## Papers

25 downloaded, every one verified openable with `pypdf`. Full annotations in
`papers/README.md`.

| Title | Authors | Year | File | Key info |
|---|---|---|---|---|
| Benchmarking atlas-level data integration in single-cell genomics | Luecken et al. | 2022 | `papers/luecken2022_scib_atlas_integration_benchmark.pdf` | **Deep-read.** The reference benchmark. All 14 metric definitions; the batch-vs-bio trade-off; trajectory metric is relative to *unintegrated* data, not truth |
| Signal recovery in single cell batch integration (CellANOVA) | Zhang et al. | 2023 | `papers/zhang2023_signal_recovery_batch_integration.pdf` | **Deep-read.** Global + gene-level distortion metrics; hold-out-control spike-in design |
| Batch correction methods … are often poorly calibrated | Antonsson & Melsted | 2025 | `papers/antonsson2025_batch_correction_poorly_calibrated.pdf` | **Deep-read.** Null-pseudobatch calibration design; Harmony/ComBat best, scVI/LIGER/MNN worst |
| Trade-off between conservation of biological variation and batch effect removal | Li et al. | 2022 | `papers/li2022_pareto_tradeoff_scvi.pdf` | **Deep-read.** Explicit Pareto front, but within scVI only |
| A benchmark of batch-effect correction methods for scRNA-seq | Tran et al. | 2020 | `papers/tran2020_benchmark_batch_correction.pdf` | **Deep-read.** kBET/LISI/ASW/ARI; recommends Harmony; kBET k-sweep convention |
| Splatter: simulation of scRNA-seq data | Zappia et al. | 2017 | `papers/zappia2017_splatter.pdf` | Generative model ported to `code/splatpy/` |
| dyngen: a multi-modal simulator of single cells | Cannoodt et al. | 2021 | `papers/cannoodt2021_dyngen.pdf` | Gold-standard trajectories; R-only |
| Simulating multiple faceted variability (SymSim) | Zhang et al. | 2019 | `papers/zhang2019_symsim.pdf` | Batch effects on a true cell tree; R-only |
| scDesign2 | Sun et al. | 2021 | `papers/sun2021_scdesign2.pdf` | Copula simulator; realism reference; R-only |
| The shaky foundations of simulating scRNA-seq data | Crowell et al. | 2023 | `papers/crowell2023_shaky_foundations_simulation.pdf` | Simulator-realism critique; read before claiming realism for `splatpy` |
| Batch effects … corrected by matching mutual nearest neighbors | Haghverdi et al. | 2018 | `papers/haghverdi2018_mnn_correct.pdf` | MNN; `mnnpy` unbuildable here |
| BBKNN: fast batch alignment | Polański et al. | 2020 | `papers/polanski2020_bbknn.pdf` | Graph-output method |
| Comprehensive Integration of Single-Cell Data (Seurat v3) | Stuart et al. | 2019 | `papers/stuart2019_seurat_v3_integration.pdf` | R-only |
| Probabilistic harmonization … (scANVI) | Xu et al. | 2021 | `papers/xu2021_scanvi_probabilistic_harmonization.pdf` | Label-supervised; top bio-conservation in scIB |
| scGen predicts single-cell perturbation responses | Lotfollahi et al. | 2019 | `papers/lotfollahi2019_scgen.pdf` | bioRxiv version; label-supervised |
| Flexible comparison … using BatchBench | Chazarra-Gil et al. | 2021 | `papers/chazarragil2021_batchbench.pdf` | Entropy-of-mixing metrics |
| CellMixS | Lütge et al. | 2021 | `papers/lutge2021_cellmixs.pdf` | Per-cell mixing score; localises failures |
| pipeComp | Germain et al. | 2020 | `papers/germain2020_pipecomp.pdf` | Whole-pipeline evaluation |
| Influence of integration on DGE performance | Nguyen et al. | 2022 | `papers/nguyen2022_integration_influence_dge.pdf` | Integration changes downstream DE |
| The specious art of single-cell genomics | Chari & Pachter | 2023 | `papers/chari2023_specious_art_single_cell_genomics.pdf` | 2-D embeddings distort distances |
| PAGA | Wolf et al. | 2019 | `papers/wolf2019_paga.pdf` | Trajectory/topology inference |
| Best practices for single-cell analysis | Heumos et al. | 2023 | `papers/heumos2023_best_practices_singlecell.pdf` | Pipeline defaults |
| Eleven grand challenges in single-cell data science | Lähnemann et al. | 2020 | `papers/lahnemann2020_grand_challenges.pdf` | Frames integration as an open problem |
| Assessing and mitigating batch effects in large-scale omics | Zhou et al. | 2024 | `papers/zhou2024_assessing_mitigating_batch_effects.pdf` | Cross-omics view |
| Integrating scRNA-seq datasets with substantial batch effects | Kang et al. | 2025 | `papers/kang2025_integrating_substantial_batch_effects.pdf` | Strong-batch-effect regime |

### Abstract only (4) — `papers/abstracts_only/README.md`

`www.biorxiv.org` returns HTTP 429 / Cloudflare error 1015 for every request from
this workspace's IP, and none of these four has a PMC or OA copy.

| Title | Why it matters | Substitute |
|---|---|---|
| Erasure of Biologically Meaningful Signal … (Tyler et al.) | Closest prior statement of this hypothesis | Abstract is detailed enough to position the work |
| The differential impacts of dataset imbalance … (Maan et al.) | Imbalance drives biology loss → sweep it | Code cloned: `code/Iniquitate/` |
| A comparison of single-cell trajectory inference methods (Saelens et al.) | Source of the trajectory-accuracy metric family | Metrics restated in Luecken 2022 |
| A test metric for assessing scRNA-seq batch correction (kBET, Büttner et al.) | kBET definition | Fully restated in Luecken 2022 Methods (read); implemented in `scib_metrics.kbet` |

Also unobtainable as PDFs: **Harmony** (Korsunsky 2019) and **Scanorama**
(Hie 2019) — PMC author manuscripts outside the OA subset, publisher PDFs return
HTML. Both methods are installed and verified; both are described in Tran 2020
and Luecken 2022.

---

## Datasets

Total 2.9 GB, excluded from git via `datasets/.gitignore`. Reproduce with
`bash datasets/download.sh`. Every file was loaded and inspected post-download;
full obs/label breakdowns in `datasets/README.md`.

| Name | Source | Size | Batches | Trajectory ground truth | Location |
|---|---|---|---|---|---|
| scIB pancreas | figshare 12420968 | 16,382 × 19,093 | 9 protocols (`tech`) | none (discrete types) | `datasets/scib/human_pancreas_norm_complexBatch.h5ad` |
| scIB immune human | figshare 12420968 | 33,506 × 12,303 | 10 donors, 5 studies (`batch`) | **`obs['dpt_pseudotime']`**, HSPC→erythrocyte | `datasets/scib/Immune_ALL_human.h5ad` |
| scIB lung atlas | figshare 12420968 | 32,472 × 15,148 | 16 donors (`batch`), 2 protocols, 2 locations | none; **batch/biology confounded** | `datasets/scib/Lung_atlas_public.h5ad` |
| Palantir CD34+ marrow rep1/2/3 | dp-lab-data-public S3 | 5,780 / 6,501 / 12,046 cells | 3 real donors | **`palantir_pseudotime` + `obsm['palantir_branch_probs']`** | `datasets/palantir_cd34/` |
| Paul15 hematopoiesis | `sc.datasets.paul15()` | 2,730 × 3,451 | single | `paul15_clusters` (MEP/GMP branches) | `datasets/trajectory/paul15_hematopoiesis.h5ad` |
| Pancreas endocrinogenesis | `scv.datasets.pancreas()` | 3,696 × 27,998 | single | Ductal→Ngn3→Pre-endocrine→α/β/δ/ε | `datasets/trajectory/pancreas_endocrinogenesis.h5ad` |
| Setty marrow (rep1, w/ splicing) | `scv.datasets.bonemarrow()` | 5,780 × 14,319 | single | `palantir_pseudotime` | `datasets/trajectory/setty_bonemarrow_cd34.h5ad` |
| Dentate gyrus | `scv.datasets.dentategyrus()` | 2,930 × 13,913 | single | neurogenesis, 14 clusters | `datasets/trajectory/dentategyrus.h5ad` |
| T1D islets (bonus) | bundled with cellanova | 20,000 × 26,163 | donors + disease state | `PseudoState` | `code/cellanova/data/t1d_example.h5ad` |
| Synthetic | `code/splatpy` | generated | arbitrary | **`Step` = exact pseudotime, `Group` = branch, plus a matched batch-free counterfactual** | generated on demand |

**Downloaded but unusable**: `moignard15_blood_dev.h5ad` — qPCR, 42 genes,
`.X` dtype `object`.

**Deliberately not downloaded**: the 1M-cell mouse brain task and 6 scATAC tasks
from the same figshare article (out of scope); `Immune_ALL_hum_mou.h5ad`
(4.3 GB cross-species — batch and biology are inseparable there, see
`planning.md` R4).

---

## Code repositories

15 cloned + 1 written. Full detail, including what runs and what does not, in
`code/README.md`.

| Name | URL | Purpose | Location | Status |
|---|---|---|---|---|
| **splatpy** | *written here* | NumPy port of the Splat model; `paired_simulation()` gives the batch-free counterfactual | `code/splatpy/` | **Working, smoke-tested** |
| scib | github.com/theislab/scib | Reference metrics incl. `trajectory_conservation` | `code/scib/` | Installed (1.1.7); `kBET` submodule needs R |
| scib-metrics | github.com/YosefLab/scib-metrics | Pure-Python kBET, ASW-batch, LISI, etc. | `code/scib-metrics/` | Installed (0.6.0), verified; jax is CPU-only |
| cellanova | github.com/Janezjz/cellanova | Signal recovery + distortion metrics + T1D dataset | `code/cellanova/` | Installed from bundled wheel |
| scgen | github.com/theislab/scgen | Label-supervised integration | `code/scgen/` | **Patched twice**, editable install, verified |
| scanorama | github.com/brianhie/scanorama | MNN-family integration | `code/scanorama/` | PyPI version verified |
| scarches | github.com/theislab/scarches | Reference mapping / scPoli | `code/scarches/` | Cloned, not installed, not tested |
| Palantir | github.com/dpeerlab/Palantir | Source of the CD34+ pseudotime/branch annotations | `code/Palantir/` | Cloned; install if branch probs must be recomputed |
| Iniquitate | github.com/hsmaan/Iniquitate | Dataset-imbalance analysis workflow | `code/Iniquitate/` | Reference (Snakemake) |
| scib-pipeline | github.com/theislab/scib-pipeline | Reproduces scIB; pins exact method invocations | `code/scib-pipeline/` | Reference only (needs conda + R) |
| splatter | github.com/Oshlack/splatter | Source of the ported model | `code/splatter/` | Reference (R) |
| dyngen | github.com/dynverse/dyngen | GRN simulator, gold-standard trajectories | `code/dyngen/` | Reference (R) |
| kBET | github.com/theislab/kBET | Original kBET | `code/kBET/` | Reference (R) |
| CellMixS | github.com/almutlue/CellMixS | Per-cell mixing scores | `code/CellMixS/` | Reference (R) |
| harmony (R) | github.com/immunogenomics/harmony | Harmony reference impl. | `code/harmony-R/` | Reference (R) |
| batchbench | github.com/cellgeni/batchbench | Entropy-of-mixing benchmark | `code/batchbench/` | Reference (Nextflow + R) |

### Integration methods verified end-to-end

Smoke-tested on `splatpy` output (800 cells × 1,200 genes):
ComBat 0.1 s · BBKNN 11 s · Scanorama 0.3 s · Harmony <0.1 s · scVI 8 s (GPU) ·
scANVI 1.3 s (GPU) · scGen 1.5 s (GPU).

**Unavailable**: MNN/fastMNN (`mnnpy` won't build on Python 3.12 — needs
`longintrepr.h`), LIGER (`pyliger` → `louvain` build failure), Seurat v3 (R).

### Feasibility check already run

A batch-strength sweep with Harmony confirmed the whole chain works and that the
effect is measurable (`splatpy` → PCA → Harmony → kBET/ASW-batch/trajectory):

| `batch_fac_loc` | method | ASW-batch | kBET | max\|ρ(PC, Step)\| |
|---|---|---|---|---|
| 0.0 | uncorrected | −0.000 | 0.965 | 0.849 |
| 0.1 | uncorrected | 0.076 | 0.101 | 0.788 |
| 0.1 | harmony | 0.002 | 0.954 | 0.837 |
| 0.2 | uncorrected | 0.207 | 0.000 | 0.851 |
| 0.2 | harmony | 0.202 | 0.001 | 0.850 |
| 0.4 | uncorrected | 0.408 | 0.000 | 0.842 |
| — | batch-free oracle | 0.000 | 0.941 | 0.844 |

Two things to carry forward. **Harmony did nothing at `batch_fac_loc ≥ 0.2`**
("Converged after 1 iteration") — likely a `theta`/`nclust` configuration issue
on 1,200 cells rather than a real failure; resolve before reporting.
And **`max|ρ(PC, Step)|` is too insensitive** to serve as the trajectory metric —
it barely moved while the batch effect went from 0 to 0.408. Use variance
attributable to `Step` across the whole embedding, or pseudotime recovery
against the oracle, instead.

---

## Notes on resource gathering

### Search strategy

The paper-finder service was down (`localhost:8000` refused connections), so
retrieval was manual and scripted:

- `tools/search_epmc.py` — Europe PMC REST search, `resultType=core`,
  citation-sorted, emitting JSONL with abstracts. 18 queries across two rounds,
  logged in `logs/search/`.
- Semantic Scholar Graph API and Unpaywall for resolving PDFs by DOI.
- `tools/fetch_paper.py` — multi-source downloader (Europe PMC render → PMC
  direct → Unpaywall OA locations → Semantic Scholar `openAccessPdf` → bioRxiv),
  validating the `%PDF-` magic bytes and a minimum size so HTML error pages
  aren't saved as `.pdf`.

Round 1 queries were broad (benchmarking, kBET, over-correction, simulation,
trajectory, metrics, spike-in); round 2 targeted specific known works and
recovered the papers that matter most here (CellANOVA, Antonsson, Splatter,
dyngen, SymSim). Citation-sorted Europe PMC results carry substantial off-topic
noise, so all candidates were screened by title and abstract before download.

### Selection criteria

Prioritised, in order: papers that measure the batch-vs-biology trade-off;
benchmark papers whose metrics or datasets are directly reusable; simulators
that can generate ground-truth trajectories; and the methods being benchmarked.
Deprioritised new-method papers with no benchmarking content, and anything on
scATAC, spatial, or multi-omic integration.

### Challenges encountered

1. **paper-finder unavailable** → replaced with a scripted Europe PMC workflow.
2. **bioRxiv blocked** (HTTP 429 / Cloudflare 1015, persistent across retries
   and user agents) → 4 papers reduced to abstracts, retrieved via the bioRxiv
   *API* (`api.biorxiv.org`, which is not blocked). None is blocking: kBET's
   definition is fully restated in Luecken 2022, and Maan's code is cloned.
3. **Harmony and Scanorama PDFs unobtainable** — PMC author manuscripts outside
   the OA subset; `nature.com` and PMC both return HTML. Not blocking.
4. **No R and no root** → `apt-get install r-base-core` fails; conda is
   disallowed. This eliminated splatter, dyngen, SymSim, scDesign2, Seurat, the
   R kBET, CellMixS and BatchBench. Mitigated by porting the Splat model to
   `code/splatpy/` directly from the R source, and by using `scib_metrics` for
   pure-Python kBET.
5. **`uv add` failed initially** — hatchling could not find a package to build.
   Fixed with `[tool.uv] package = false` in `pyproject.toml`.
6. **Package incompatibilities**: `scgen` 2.1.0 on PyPI imports the removed
   `scvi._compat`; `scgen` master then broke twice against scvi-tools 1.5 and
   anndata 0.13. Both patched (documented in `code/README.md`). `mnnpy` and
   `pyliger` could not be built at all.
7. **`harmonypy` 2.0 changed `Z_corr` orientation** — it is now
   `(n_cells, n_dims)`. Transposing it, as older tutorials do, fails silently
   in shape-tolerant code.

### Gaps and workarounds

| Gap | Workaround | Residual risk |
|---|---|---|
| No R simulator | `code/splatpy` ported from splatter's R source | Not validated against R output; realism unverified (Crowell 2023) |
| MNN, LIGER unavailable | Scanorama + BBKNN as MNN-family stand-ins | The batch-removal extreme of the frontier is under-sampled — must be stated |
| 4 papers abstract-only | Content recovered from abstracts + citing papers | Low; kBET and Saelens metrics are restated in Luecken 2022 |
| jax CPU-only | Metrics still run | Slower; `uv add "jax[cuda12]"` if it becomes the bottleneck |
| ~42 GB of the GPU already occupied | Batch sizes must stay modest | scVI/scANVI/scGen all ran fine in the smoke test |

---

## Recommendations for experiment design

1. **Primary datasets** — synthetic paired simulations from `code/splatpy`
   (the only source of a true counterfactual); `Immune_ALL_human.h5ad` and
   `datasets/palantir_cd34/` for real trajectories under real batches;
   `paul15` / `pancreas_endocrinogenesis` / `setty_bonemarrow_cd34` as
   injection substrates; `human_pancreas` (easy) and `Lung_atlas_public`
   (confounded) as the two calibration extremes.
2. **Baselines** — ComBat, BBKNN, Scanorama (matrix and embedding separately),
   Harmony, scVI, scANVI, scGen, plus two essential reference points: the
   **unintegrated** data and, in the synthetic arm, the **batch-free oracle**.
3. **Metrics** — batch axis: kBET (swept over k, per Tran 2020) and ASW-batch,
   both named in the hypothesis, plus graph iLISI and PCR-batch. Biology axis:
   ground-truth trajectory variance retained (D1/D2) as primary, with
   `scib.metrics.trajectory_conservation` reported alongside for comparability,
   plus CellANOVA's global and gene-level distortion. Report the axes
   separately; the scIB 0.6/0.4 composite is a secondary summary only.
4. **Code to reuse** — `scib_metrics` for all pure-Python metrics; `scib` for
   `trajectory_conservation`, `hvg_overlap` and `cell_cycle`; `code/splatpy`
   for simulation; `code/cellanova/tutorials/` for the distortion measurements;
   `code/scib-pipeline/` as the spec for how each method should be invoked.
5. **Design discipline** — multiple seeds with error bars (non-domination is
   meaningless without variance); sweep batch strength, batch–biology
   confounding, and cell-type composition imbalance (Maan 2024); include the
   no-batch-effect null (Antonsson 2025) where the correct amount of change is
   exactly zero.
