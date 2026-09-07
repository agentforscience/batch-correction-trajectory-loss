# Literature Review

**Topic.** How much biology do batch-correction methods remove? A controlled
benchmark using synthetic and real scRNA-seq spike-ins.

25 PDFs in `papers/` (5 read in full via chunked PDFs, the rest skimmed from
abstracts and targeted sections), 4 further papers as abstracts in
`papers/abstracts_only/`. Search was manual: the paper-finder service at
`localhost:8000` was not running, so retrieval used the Europe PMC REST API
(`tools/search_epmc.py`), plus Semantic Scholar and Unpaywall for PDF resolution.

---

## 1. Research area overview

Batch effects are the central obstacle to combining scRNA-seq datasets, and at
least 49 integration methods existed by late 2020 (Zappia's scRNA-tools census,
via Luecken 2022). The field's evaluation culture formed around a single
framing: **remove batch effect, preserve cell-type structure**. Early benchmarks
(Tran 2020) scored methods on kBET, LISI, ASW and ARI, all of which reduce
"biology" to discrete cell-type labels. Luecken 2022 broadened this with
label-free metrics — cell-cycle variance, HVG overlap, and trajectory
conservation — and, crucially, framed the result as a *trade-off*.

Since ~2021 a distinct, smaller literature has argued that the trade-off is not
merely a nuisance but a systematic failure: correction removes real biology.
Tyler et al. (2021/2023) titled their preprint *Erasure of Biologically
Meaningful Signal*. Zhang et al. (2023, CellANOVA) opened with "how much
biological signal is erased during integration?" — the same question as this
project — and answered it by building a recovery method. Antonsson & Melsted
(2025) showed that most methods are not even *calibrated*: they alter data that
contains no batch effect at all.

**The synthesis: everyone agrees the trade-off exists; nobody has measured it
against a ground truth.** That is the gap. See `planning.md`.

---

## 2. Key papers

### Luecken, Büttner, …, Colomé-Tatché & Theis (2022) — scIB
*Benchmarking atlas-level data integration in single-cell genomics.*
Nature Methods 19:41-50. `papers/luecken2022_scib_atlas_integration_benchmark.pdf`

- **Contribution.** 68 method+preprocessing combinations, 85 batches, >1.2M
  cells, 13 tasks, 14 metrics. The field's reference benchmark.
- **Methodology.** Metrics in two groups. *Batch removal*: PCR-batch, batch ASW,
  graph iLISI, graph connectivity, kBET. *Bio-conservation*, split into
  label-based (NMI, ARI, cell-type ASW, isolated-label F1 and silhouette,
  graph cLISI) and label-free (cell-cycle conservation, HVG conservation,
  trajectory conservation). Aggregate = 0.6·bio + 0.4·batch, each metric
  min–max scaled within a task. New extensions of kBET and LISI to graph outputs
  so that graph-, embedding-, and matrix-output methods are comparable.
- **Datasets.** Pancreas (16,382 cells / 9 batches), Lung (32,472 / 16),
  Immune human (33,506 / 10), Immune human+mouse (97,952 / 23), mouse brain
  (978,734 / 4), 6 ATAC tasks, and **2 Splatter simulations**
  (12,097 cells / 6 batches — composition variation; 19,318 / 16 — nested batch
  effects). The first four are in `datasets/scib/`.
- **Results.** scANVI, Scanorama, scVI and scGen perform well overall.
  "BBKNN and Seurat v3 tended to remove batch variation, scANVI and scGen
  prioritized bio-conservation." Harmony works on simple tasks and drops out of
  the top three on complex ones. Scaling shifts everything toward batch removal
  at the cost of bio-conservation; HVG selection helps overall.
- **The limitation that matters here.** Trajectory conservation =
  `(spearman(dpt_post, dpt_pre) + 1) / 2`, where `dpt_pre` is diffusion
  pseudotime on the **unintegrated** data. On real data with a real batch
  effect, `dpt_pre` is corrupted, so a method that correctly removes the batch
  effect and a method that destroys the trajectory both score low. No percentage
  of "biology removed" can be derived from it.
- Code: `code/scib/`, `code/scib-pipeline/`. Data: `datasets/scib/`.

### Zhang, Mathew, …, Ma & Zhang (2023) — CellANOVA
*Signal recovery in single cell batch integration.* bioRxiv 2023.05.05.539614.
`papers/zhang2023_signal_recovery_batch_integration.pdf`

- **Contribution.** "Current paradigms for single cell data integration are
  unnecessarily aggressive, removing biologically meaningful variation." Presents
  CellANOVA, which uses a **pool-of-controls** design to learn a latent batch
  basis `V` from samples that should not differ, then recovers any variation
  orthogonal to `V` in the non-control samples.
- **Two distortion measures, both directly reusable.**
  1. *Global distortion* — per-cell Pearson correlation between pre- and
     post-integration expression vectors. For Harmony, Seurat RPCA and LIGER
     this averages **below 0.5**; after CellANOVA recovery it exceeds 0.9.
  2. *Gene-specific distortion* — hold out one sample, run DE between cell types
     **within that single sample** before and after correction, and correlate the
     adjusted p-values. All cells come from one batch, so pre-correction p-values
     are unconfounded and any shift is artifact. Existing methods "artificially
     reduce the p-values", inflating type-1 error.
- **A real-data spike-in protocol.** Relabel one control sample as a fake
  treatment sample; any signal recovered for it is false by construction. Applies
  to the T1D islet data (11 healthy / 5 T1D / 8 AAB+), an NSCLC immunotherapy
  longitudinal cohort, a mouse irradiation study, and a scRNA-vs-snRNA kidney study.
- **Recovery readout.** For each cell take its 30 *out-of-batch* nearest
  neighbours and ask what fraction share its condition label. Since the condition
  label is never used during recovery, enrichment is evidence of preserved signal.
- Code: `code/cellanova/` (installed). Bundles `data/t1d_example.h5ad`
  (20,000 cells × 26,163 genes, real donors + disease state).

### Antonsson & Melsted (2025)
*Batch correction methods used in scRNA-seq analyses are often poorly calibrated.*
Genome Research 35:1832-1841. `papers/antonsson2025_batch_correction_poorly_calibrated.pdf`

- **Contribution.** The **null-calibration** design: split a real dataset into
  two random pseudobatches — no batch effect exists — and measure what each
  method changes anyway. "Ideally, the application of batch effect correction
  should not correct the data at all… Under this null hypothesis, any significant
  change can be classified as an artifact."
- **Methodology.** 25 random splits per dataset (PBMC3K, PBMC4K, mouse brain,
  human jejunum, mouse heart); measure change at the kNN graph, at Leiden
  clustering (ARI vs uncorrected), and at DE (MAST, Bonferroni α=0.05).
  Baseline: downsample one batch's reads 50%.
- **Results.** Three tiers. ComBat and Harmony best (ARI 0.93 / 0.92 on PBMC3K,
  both 0.78 on mouse brain), Seurat middle (0.82 / 0.71), and MNN, scVI, LIGER
  clearly worst — LIGER and scVI created clusters that were pure batch artifacts.
  Only ComBat and Harmony perturbed the data less than the 50%-downsampling
  baseline. On DE between B cells and CD8+ T cells, MNN and Seurat reported
  **>800 DE genes where the uncorrected data had 179**.
  Recommendation: Harmony only.
- **Note the tension with scIB**, which ranks Harmony outside the top three on
  complex tasks while scVI ranks high. The two studies optimise different things
  — scIB rewards batch removal at 40% weight, this paper penalises any change at
  all — and that tension is itself evidence for the Pareto framing.

### Li, McCarthy, Shim & Wei (2022)
*Trade-off between conservation of biological variation and batch effect removal
in deep generative modeling.* BMC Bioinformatics 23:460.
`papers/li2022_pareto_tradeoff_scvi.pdf`

- The only paper that constructs an explicit **Pareto front** between the two
  objectives, using Pareto multi-task learning on scVI's ZINB objective, and
  showing it beats naive scalarisation. Introduces MINE (mutual information
  neural estimation) as a batch-effect measure, arguing it improves on MMD.
- **Scope limit**: one method's achievable frontier, not the cross-method
  empirical frontier the present hypothesis concerns. Complementary, and the
  reason `planning.md` folds per-method hyperparameter sweeps into D1 — comparing
  methods at single default points would be a weaker claim than comparing curves.

### Tran, Ang, Chevrier, Zhang et al. (2020)
*A benchmark of batch-effect correction methods for scRNA-seq data.*
Genome Biology 21:12. `papers/tran2020_benchmark_batch_correction.pdf`

- 14 methods, 10 datasets, 5 scenarios (same cell types / different protocols;
  non-identical cell types; >2 batches; >500k cells; simulation), scored on
  kBET, LISI, ASW, ARI, plus a DE F-score.
- Recommends **Harmony first** (runtime), with LIGER and Seurat 3 as alternatives.
- Two conventions worth adopting: kBET rejection rates computed at local sample
  sizes of 5/10/15/20/25% (kBET is strongly k-dependent, so a single k is
  fragile), and the use of **Splatter** simulations specifically to test the
  effect of correction on downstream DE recovery.
- Explicitly warns that ASW/ARI alone "may not give a full picture."

### Tyler, Guccione & Schadt (2021/2023) — abstract only
*Erasure of Biologically Meaningful Signal by Unsupervised scRNAseq
Batch-correction Methods.* 10.1101/2021.11.15.468733.

The closest prior statement of this project's hypothesis. Introduces two
measures of batch effect — at the level of cluster composition, and along
overlapping topologies — and finds that standard correction "erase[s] cell-type
and cell-state variation in real-world biological datasets, single cell gene
expression atlases, and in silico experiments". Warns these erasures "may create
the artefactual appearance of external validation/replication of findings".
Concludes that either biological effects must be balanced across batches (as in
bulk experimental design) or technical effects must be explicitly modelled.
**PDF unobtainable** — bioRxiv returns HTTP 429 from this workspace.

### Maan, Zhang, …, Wang (2022/2024) — abstract only
*The differential impacts of dataset imbalance in single-cell data integration.*
Nat Biotechnol 2024. Cell-type imbalance across batches drives how much biology
is destroyed — which makes composition imbalance a required covariate to sweep,
not an incidental detail. Code: `code/Iniquitate/` (cloned).

---

## 3. Simulation literature

- **Zappia, Phipson & Oshlack (2017), Splatter** — `papers/zappia2017_splatter.pdf`.
  Gamma-Poisson hierarchical model: gene means ~ Gamma, library sizes ~ log-normal,
  BCV-inflated means via Gamma, counts via Poisson, optional logistic dropout.
  Batch effects enter as **per-gene log-normal multiplicative factors**
  (`batch.facLoc`, `batch.facScale`); `method="paths"` places cells along
  Brownian-bridge trajectories between DE endpoints. `batch.rmEffect=TRUE` yields
  a matched batch-free counterfactual. Ported to `code/splatpy/`.
- **Cannoodt et al. (2021), dyngen** — GRN-driven, gold-standard trajectories,
  multi-modal. The better tool; R-only, so unusable here (see `planning.md` R3).
- **Zhang, Xu & Yosef (2019), SymSim** — explicitly models "multiple faceted
  variability" including batch effects on a true cell tree. R-only.
- **Sun et al. (2021), scDesign2** — copula-based, fitted to real data; the
  standard for "does my simulation look real". R-only.
- **Crowell, Leonardo, Soneson & Robinson (2023), *The shaky foundations of
  simulating scRNA-seq data*** — `papers/crowell2023_shaky_foundations_simulation.pdf`.
  Systematic critique of simulator realism. **Read before making any realism
  claim about `splatpy`**; it also supplies the criteria for how to check.

---

## 4. Standard metrics (definitions as used in scIB unless noted)

**Batch removal**

| Metric | Definition | Implementation here |
|---|---|---|
| kBET | χ² test that the batch composition of each cell's k-NN neighbourhood matches the global composition; reported as `1 −` rejection rate, computed per cell-type label and averaged. k defaults to the median per-batch cell count within a label, clipped to [10,100]. | `scib_metrics.kbet` (pure Python) |
| ASW-batch | Per-cell `\|silhouette(batch)\|`, `1 − ` that, averaged within each cell-type label then across labels. 1 = ideal mixing. | `scib_metrics.silhouette_batch` |
| graph iLISI | Inverse Simpson's index over batch labels in graph-distance neighbourhoods, median over cells, rescaled to [0,1]. | `scib_metrics.ilisi_knn` |
| PCR-batch | Σ over PCs of `Var(C\|PC_i) × R²(PC_i \| batch)` — total variance attributable to batch. | `scib_metrics.pcr_comparison` |
| graph connectivity | Mean over cell types of `\|largest connected component\| / \|cells\|` in the label-subset kNN graph. | `scib_metrics.graph_connectivity` |

**Bio-conservation**

| Metric | Definition | Implementation here |
|---|---|---|
| NMI / ARI | vs resolution-optimised Louvain/Leiden clustering (scIB sweeps resolution 0.1–2.0 in steps of 0.1 and takes the best NMI). | `scib_metrics.nmi_ari_cluster_labels_leiden` |
| cell-type ASW | `(ASW(labels) + 1)/2`. | `scib_metrics.silhouette_label` |
| graph cLISI | iLISI's dual on cell-type labels. | `scib_metrics.clisi_knn` |
| isolated-label F1 / silhouette | Scores restricted to labels present in the fewest batches. | `scib_metrics.isolated_labels` |
| HVG conservation | Overlap coefficient of the top-500 HVGs per batch, pre vs post. | `scib.metrics.hvg_overlap` |
| cell-cycle conservation | `1 − \|Var_after − Var_before\| / Var_before` for S/G2M score variance, per batch. | `scib.metrics.cell_cycle` |
| **trajectory conservation** | `(spearman(dpt_post, dpt_pre) + 1)/2`; 0 if the graph is too disconnected to compute DPT. | `scib.metrics.trajectory_conservation` |

**Distortion (from CellANOVA, not in scIB)** — per-cell pre/post expression
correlation; within-batch DE p-value correlation pre vs post; out-of-batch 30-NN
condition enrichment.

**Aggregation.** scIB: `0.6·bio + 0.4·batch`, min–max scaled per task. Note that
this weighting is a *choice*, and a hypothesis about a Pareto frontier should
report the two axes separately rather than collapse them.

---

## 5. Gaps and opportunities

1. **No ground-truth denominator.** Every published trajectory-preservation
   number is relative to the unintegrated data, not to truth. A paired
   counterfactual makes "% of trajectory variance removed" computable.
2. **Trajectories are under-measured.** scIB has one trajectory metric on one
   lineage in one task. The trajectory-inference field (Saelens 2019) has a
   richer metric family — topology (HIM), branch assignment (F1_branches),
   geodesic-distance correlation — none of which has been applied to integration.
3. **Methods are compared at single default points.** A "no method dominates"
   claim needs per-method curves over their own hyperparameters (Li 2022 did this
   for scVI alone).
4. **Benchmarks disagree, and nobody has reconciled them.** Tran 2020 recommends
   Harmony/LIGER/Seurat; scIB ranks scANVI/Scanorama/scVI/scGen; Antonsson 2025
   recommends Harmony only and puts scVI and LIGER at the bottom. These are
   consistent if the methods sit at different points on a frontier and the three
   studies weight the axes differently — which is testable.
5. **Composition imbalance is a first-order confounder** (Maan 2024) and is
   rarely swept.

## 6. Recommendations for the experiment

- **Datasets.** Synthetic: `code/splatpy` paired simulations (D1). Real spike-in:
  `paul15_hematopoiesis`, `pancreas_endocrinogenesis`, `setty_bonemarrow_cd34`
  with injected effects estimated from real protocol contrasts (D2). Real
  batches: `datasets/palantir_cd34/` (3 donors + branch probabilities) and
  `Immune_ALL_human.h5ad` (10 batches + `dpt_pseudotime`); `human_pancreas` as
  the easy calibration task and `Lung_atlas_public` as the confounded task (D3).
- **Methods.** ComBat, BBKNN, Scanorama (matrix and embedding scored
  separately), Harmony, scVI, scANVI, scGen — all verified working. MNN and
  LIGER are unavailable (see `code/README.md`); since scIB places LIGER at the
  batch-removal extreme, the measured frontier under-samples that end and this
  must be stated.
- **Metrics.** Batch axis: kBET and ASW-batch (named in the hypothesis), plus
  graph iLISI and PCR-batch. Biology axis: the ground-truth trajectory measures
  from D1/D2 as the primary readout, with `scib.metrics.trajectory_conservation`
  reported alongside for comparability to published numbers. Add CellANOVA's
  global and gene-level distortion measures. Report the two axes separately;
  give the scIB 0.6/0.4 composite only as a secondary summary.
- **Design.** Several seeds per configuration with error bars — a Pareto-frontier
  claim is a claim about non-domination, which is meaningless without variance.
  Sweep batch strength and batch–biology confounding; include the
  no-batch-effect null (Antonsson) as the anchor point where the correct amount
  of change is zero, and the unintegrated data as the anchor at the other end.
- **Watch out for.** The Harmony non-convergence observed during setup
  (`code/README.md`) — resolve it before reporting it. kBET's k-sensitivity —
  follow Tran 2020 and sweep k. `harmonypy` 2.0's un-transposed `Z_corr`.
  And avoid reading any metric off a 2-D embedding (Chari & Pachter 2023).
