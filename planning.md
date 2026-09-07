# Planning — *How Much Biology Do Batch Correction Methods Remove?*

> Phase 0 (Motivation & Novelty) and Phase 1 (Experimental Design) were written by
> the `experiment_runner` agent on 2026-09-06. The **direction budget** section
> below (D1/D2/D3 kept, R1–R7 pruned) was produced by the `resource_finder`
> agent and is retained verbatim; the ranking was **not** revised.

---

## Motivation & Novelty Assessment

### Why This Research Matters

Batch correction is a mandatory, load-bearing step in essentially every
multi-sample single-cell atlas: the corrected embedding is what clustering,
cell-type annotation, trajectory inference and differential expression are all
computed on. Method choice is currently made by consulting benchmark
leaderboards (scIB being the standard), which score methods on a composite of
"batch removal" and "bio conservation". But the bio-conservation side of that
composite is measured *relative to the uncorrected data*, which under a real
batch effect is itself corrupted. If the leaderboard's notion of "biology
preserved" is systematically biased, then thousands of atlases are being built
on a method ranking that nobody has validated against a known truth. The people
who benefit from resolving this are anyone building or reading a
multi-batch atlas — and specifically anyone drawing a developmental trajectory
through one, because trajectories are continuous, low-variance structures and
therefore the biology most at risk from over-correction.

### Gap in Existing Work

From `literature_review.md`, three papers approach this claim and none closes it:

- **Luecken et al. 2022 (scIB)** reports the trade-off qualitatively and its
  trajectory metric is `(spearman(dpt_post, dpt_pre) + 1)/2` — the reference is
  the *unintegrated* pseudotime. It cannot separate "removed biology" from
  "removed batch", and cannot yield a percentage.
- **Zhang et al. 2023 (CellANOVA)** quantifies distortion as deviation from the
  original data, but distortion-from-original is not distortion-from-truth.
- **Antonsson & Melsted 2025** measures what correction does when there is *no*
  batch effect. That isolates the artifact but removes the trade-off, since with
  no batch effect there is nothing to trade against.
- **Li et al. 2022** draws an explicit Pareto front, but only *within* scVI by
  varying one loss weight — not *across* methods.

What is missing is a design in which the batch-free truth is known **for the
same cells**, so that "fraction of trajectory variance removed" is a computed
number rather than a proxy, and in which that number can be placed on the same
axes as kBET and ASW-batch to test whether the frontier is genuinely
non-dominated.

### Our Novel Contribution

Four things, in decreasing order of novelty:

1. **A counterfactual ground truth for trajectory destruction.** `splatpy`'s
   `paired_simulation()` runs the Splat model twice on one RNG stream with the
   gene-wise batch factors toggled, producing the *same cells* differing only by
   batch. The batch-free arm is an exact oracle. We define
   **TVR (Trajectory Variance Retained)** = the fraction of embedding variance
   attributable to the true pseudotime, normalised by the oracle's — so
   "biology removed" = 1 − TVR is directly the quantity the hypothesis names a
   range for (15–40%).
2. **A four-way decomposition that no prior benchmark reports together**:
   `oracle` (batch-free, uncorrected), `uncorrected` (batch present, no
   correction), `method-on-batch` (the usual measurement), and
   **`method-on-oracle`** — running each correction method on data with batch
   *labels* but no batch *effect*. The last isolates pure artifact, extending
   Antonsson's null design to a setting where the trade-off still exists.
3. **A direct audit of the published metric.** We compute scIB's
   `trajectory_conservation` alongside our oracle-referenced TVR on the same
   embeddings. If the two disagree, the scIB trajectory ranking is measurably
   biased, and we can say in which direction and by how much.
4. **Statistical Pareto analysis rather than eyeballing.** Non-domination is
   tested with bootstrap CIs over seeds, not asserted from a scatter plot.

### Experiment Justification

- **E1 — Batch-strength sweep (D1).** *Why:* the hypothesis states a percentage
  range; a percentage is only meaningful as a function of how large the batch
  effect is. Sweeping `batch_fac_loc` ∈ {0, 0.1, 0.2, 0.35, 0.5} tells us
  whether 15–40% is the right ballpark and where in the regime it holds.
  *What we learn:* the primary numbers, plus whether the ordering of methods is
  stable across regimes.
- **E2 — Null / artifact arm (D1, `method-on-oracle`).** *Why:* if a method
  removes 25% of trajectory variance when there is no batch effect at all, then
  its 25% under a real batch effect is not a trade-off — it is unconditional
  damage. This separates "the price of correction" from "the price of the
  method". *What we learn:* the intercept of each method's damage curve.
- **E3 — Batch–biology confounding sweep (D1).** *Why:* Maan et al. argue
  composition imbalance drives biology loss; the hypothesis predicts the
  frontier *widens* under confounding. This is the sharpest available test of
  the "no method dominates" clause. *What we learn:* whether non-domination is a
  property of the methods or of the regime.
- **E4 — Real-expression spike-in (D2).** *Why:* `splatpy` is unvalidated
  against R splatter and Crowell et al. 2023 argue simulators are unrealistic.
  E4 keeps the counterfactual logic but replaces simulated expression with real
  `paul15` / `pancreas_endocrinogenesis` expression, injecting a batch effect
  whose per-gene factors are **estimated from a real protocol contrast** in the
  scIB pancreas data. *What we learn:* whether E1's percentages survive on real
  expression — the credibility check without which E1 stands alone.
- **E5 — Real multi-donor trajectories (D3).** *Why:* no ground truth exists
  here, but this is the setting practitioners actually face. We report
  within-donor pseudotime preservation (within one donor there is no batch
  effect, so any reordering is method artifact) on Palantir CD34+ and
  Immune_ALL_human. *What we learn:* whether the frontier shape found under
  controlled conditions is visible in the wild.
- **E6 — Metric audit.** *Why:* contribution 3 above. Runs at zero extra cost
  on E1/E4 embeddings.

---

## Research Question

Do single-cell batch-correction methods that score best on batch-mixing metrics
(kBET, ASW-batch) systematically destroy more of a known differentiation
trajectory than methods that score worse — specifically 15–40% of trajectory
variance — and is the resulting set of (batch-mixing, biology-preserved) points
a Pareto frontier on which no method dominates?

## Hypothesis Decomposition

| ID | Sub-hypothesis | Refuted if |
|---|---|---|
| **H1** | Batch correction removes measurable ground-truth trajectory variance. | 1 − TVR ≈ 0 for all methods (CI includes 0). |
| **H2** | The amount removed is in the 15–40% range. | Point estimates cluster well outside [0.15, 0.40] at realistic batch strength. |
| **H3** | Biology removed correlates *positively* with batch-mixing performance across methods. | Spearman(batch score, 1−TVR) ≤ 0 across methods. |
| **H4** | No method dominates: the (batch, bio) point set has ≥2 Pareto-optimal methods with non-overlapping CIs. | One method is simultaneously best on both axes. |
| **H5** | Confounding widens the frontier / increases biology loss. | Loss is flat in the confounding sweep. |
| **H6** | scIB's `trajectory_conservation` disagrees with oracle-referenced TVR. | The two rank methods identically (Spearman ≈ 1). |

Independent variables: method, batch-effect strength, batch–biology confounding,
data source (synthetic / real-expression spike-in / real batches), seed.
Dependent variables: batch axis (kBET acceptance rate, ASW-batch, iLISI,
PCR-batch) and biology axis (TVR, pseudotime recovery ρ vs *true* Step, kNN
overlap with oracle, branch ARI, scIB `trajectory_conservation`).

## Methodology

### Metrics (definitions)

**Trajectory Variance Retained (TVR)** — primary, novel. For an embedding `E`,
let `PCR_step(E)` be principal-component regression of the *true* pseudotime on
`E`: PCA `E`, regress each PC on `Step`, and take the variance-weighted mean R²
(this is scIB's PCR machinery applied to the biological covariate rather than
the batch covariate; it is scale- and rotation-invariant). Then

```
TVR(method) = PCR_step(E_method) / PCR_step(E_oracle)
biology removed = 1 - TVR
```

`E_oracle` = PCA of the batch-free counterfactual. Values > 1 are possible and
meaningful (a method can concentrate trajectory variance).

**Supporting biology metrics.** (a) `rho_true` — Spearman correlation of
diffusion pseudotime computed on `E_method` against the true `Step`; (b)
`knn_oracle` — mean overlap of each cell's k=30 neighbours between `E_method`
and `E_oracle`; (c) branch ARI from k-means on `E_method` against true `Group`;
(d) scIB `trajectory_conservation` for comparability.

**Batch metrics.** kBET acceptance rate (`scib_metrics.kbet`, k=25), ASW-batch
(`scib_metrics.silhouette_batch`), graph iLISI, PCR-batch. Both metrics named in
the hypothesis (kBET, ASW-batch) are primary.

### Methods compared

`uncorrected`, ComBat, Harmony, Scanorama-embed, Scanorama-matrix, BBKNN
(diffusion-map embedding of its graph), scVI, scANVI, scGen, plus the `oracle`
reference. **MNN and LIGER are unavailable in this environment** (build
failures, documented in `code/README.md`); since scIB places LIGER at the
batch-removal extreme, the frontier reported here under-samples that end.

### Statistical analysis plan (pre-registered)

- 5 seeds per configuration; all point estimates reported as mean ± bootstrap
  95% CI (10,000 resamples, percentile method).
- H1: one-sample t-test / bootstrap CI on `1 − TVR` vs 0, per method, with
  Benjamini–Hochberg FDR across methods at q = 0.05.
- H3: Spearman correlation across method means between the batch axis and
  `1 − TVR`, with a permutation p-value (10,000 permutations).
- H4: Pareto set computed on method means; a method is declared non-dominated
  robustly if it remains on the frontier in ≥95% of 10,000 bootstrap resamples
  over seeds.
- H5: two-way analysis over (method × confounding) on `1 − TVR`, reporting
  effect sizes (partial η²) rather than relying on p-values alone.
- H6: Spearman between method rankings under TVR and under
  `trajectory_conservation`, plus a Bland–Altman-style agreement plot.
- No metric or method is dropped after seeing results; all configurations run
  are reported.

### Expected outcomes

Support for the hypothesis looks like: `1 − TVR` significantly > 0 for most
methods, several methods in [0.15, 0.40] at moderate batch strength, a positive
batch-vs-biology Spearman across methods, and ≥2 robustly non-dominated methods.
Refutation looks like: near-zero or negative `1 − TVR`, no correlation between
the axes, or a single dominating method.

### Timeline

| Milestone | Est. |
|---|---|
| Harness + metric implementation, validated on a tiny config | 45 min |
| E1/E2 batch-strength + null sweep | 60 min |
| E3 confounding sweep | 30 min |
| E4 real-expression spike-in | 45 min |
| E5 real multi-donor | 30 min |
| Analysis, figures, statistics | 45 min |
| REPORT.md / README.md | 30 min |

### Potential challenges

1. **Harmony "converged after 1 iteration" at strong batch effects** (observed
   in the feasibility run). Resolved before the main sweep by testing
   `theta` / `nclust` / `max_iter_harmony` explicitly — see
   `results/harmony_diagnostic.csv`.
2. **`max|ρ(PC, Step)|` was too insensitive** as a trajectory metric. Replaced
   by TVR (variance-weighted PCR over the whole embedding), which responds to
   the batch-strength knob.
3. **splatpy realism.** Mitigated by E4/E5; stated as a limitation regardless.
4. **GPU contention** (2 of 4 A6000s are ~87% occupied). Pinned to device 1.
5. **BBKNN emits only a graph.** Handled by taking a diffusion map of its graph
   as the embedding, and flagging the choice.

### Success criteria

The study succeeds if every sub-hypothesis H1–H6 receives a clear
supported/refuted/ambiguous verdict backed by a CI or a test, whether or not the
verdicts favour the original hypothesis, and if the D1 and D2 arms are compared
head-to-head so that a disagreement between simulation and real expression is
visible rather than hidden.

---

## Direction budget (from `resource_finder`)

**Hypothesis.** Methods that achieve the best batch-mixing scores (kBET,
ASW-batch) systematically sacrifice biological signal, removing 15–40% of the
variance in known differentiation trajectories, while methods that preserve
biology best leave measurable residual batch effects — creating a Pareto
frontier that no current method dominates.

## The gap the literature leaves open

Three papers state a version of this claim; none measures it against a ground
truth.

- **Luecken et al. 2022 (scIB)** reports the trade-off qualitatively and plots
  batch-correction vs bio-conservation (Extended Data Fig. 4). But its
  trajectory-conservation metric is `(spearman(dpt_post, dpt_pre)+1)/2` — the
  reference is the *unintegrated* pseudotime, which under a real batch effect is
  itself corrupted. The metric therefore cannot separate "removed biology" from
  "removed batch", and cannot produce a percentage.
- **Zhang et al. 2023 (CellANOVA)** shows integration is "unnecessarily
  aggressive" and quantifies distortion (per-cell pre/post expression
  correlation drops below 0.5). But distortion-from-original is not
  distortion-from-truth, and their focus is a recovery method.
- **Antonsson & Melsted 2025** measures what correction does when there is *no*
  batch effect (random pseudobatches). That isolates the numerator of our
  question but not the trade-off, since with no batch effect there is nothing to
  trade against.

What is missing is a design where the batch-free truth is *known for the same
cells*, so that "fraction of trajectory variance removed" is a computable number
rather than a proxy. `code/splatpy/paired_simulation()` supplies exactly that
(same RNG stream, batch factors toggled), and the two spike-in designs below
extend the same logic to real data.

## Directions considered

Scored 1–5 on: **Ev** evidence from the literature that the direction is live,
**Rel** relevance to the hypothesis as stated, **IG** expected information gain,
**Feas** implementation feasibility in *this* environment (no R, no root).

| # | Direction | Ev | Rel | IG | Feas | Verdict |
|---|---|---|---|---|---|---|
| D1 | Paired counterfactual simulation | 5 | 5 | 5 | 5 | **KEEP** |
| D2 | Injected-batch spike-in into real trajectory data | 4 | 5 | 5 | 4 | **KEEP** |
| D3 | Real multi-donor trajectory + null-pseudobatch calibration | 5 | 4 | 4 | 5 | **KEEP** |
| R1 | scATAC / multi-omic integration | 4 | 1 | 2 | 2 | reject |
| R2 | Build a new signal-recovery method | 5 | 2 | 3 | 2 | reject |
| R3 | dyngen / SymSim GRN-mechanistic simulation | 4 | 4 | 4 | 1 | reject |
| R4 | Cross-species integration | 4 | 1 | 2 | 3 | reject |
| R5 | Within-method hyperparameter Pareto sweep | 5 | 3 | 2 | 4 | reject as a direction; fold into D1 |
| R6 | Spatial transcriptomics batch correction | 3 | 1 | 2 | 2 | reject |
| R7 | Disease/perturbation signal preservation | 5 | 2 | 4 | 4 | reject as a direction; fold into D3 |

---

## Kept: D1 — Paired counterfactual simulation

Simulate a branching differentiation trajectory with `code/splatpy`, twice per
seed: once with a batch effect, once without, same cells. Integrate the
batch arm with each method; score against the batch-free arm.

- **Ground truth available**: `Step` (true pseudotime), `Group` (true branch),
  per-gene `batch_facs` and `de_facs`, `true_counts` pre-dropout, and the
  matched no-batch expression matrix.
- **The headline number becomes computable.** "Fraction of trajectory variance
  removed" = 1 − (variance in the corrected embedding attributable to `Step`) /
  (same quantity in the batch-free oracle embedding). This is the 15–40% claim,
  measured rather than proxied.
- **Sweep** (this is where R5 lives): batch strength `batch_fac_loc/scale`
  ∈ {0, 0.05, 0.1, 0.2, 0.4}; batch–biology confounding (branch composition
  imbalance across batches, per Maan et al.); number of batches; cell count;
  and per-method key hyperparameters (Harmony `theta`, scVI `n_latent`) so that
  each method contributes a *curve*, not a point — a method is only shown to be
  non-dominated if it is non-dominated across its own tuning range.
- **Cost**: ~4 s/simulation, seconds-to-minutes per integration on the A6000.
  Cheap enough to run many seeds, which the hypothesis needs (a Pareto frontier
  claim requires error bars).
- **Known risk**: `splatpy` is unvalidated against R splatter, and Crowell et al.
  2023 argues simulators are generally unrealistic. Hence D2 and D3 are not
  optional — the simulation result is only credible if the real-data arms agree.

## Kept: D2 — Injected-batch spike-in into real trajectory data

Take a real single-batch dataset with a well-characterised trajectory
(`paul15_hematopoiesis`, `pancreas_endocrinogenesis`, `setty_bonemarrow_cd34`),
split it into pseudobatches, and inject a batch effect. The uninjected data is
the ground truth for the same cells — the counterfactual logic of D1 on real
expression.

- **Make the injected effect realistic rather than invented**: estimate
  per-gene batch factors from an actual technology contrast in
  `datasets/scib/human_pancreas_norm_complexBatch.h5ad` (9 protocols) or from
  the v2/v3 10x contrast in `Immune_ALL_human.h5ad`, then apply those factors.
  This is the "real spike-in" in the project title, as distinct from D1's
  synthetic one.
- **Sweep the confounding deliberately**: assign pseudobatches at random
  (unconfounded), by pseudotime quantile (fully confounded with the trajectory),
  and at intermediate mixtures. The hypothesis predicts the Pareto frontier
  widens as confounding increases; this is the sharpest available test of it.
- **Risk**: a multiplicative gene-wise injection is the same functional form
  `splatpy` uses, so D2 could inherit D1's bias rather than check it. Mitigate by
  also injecting non-multiplicative effects (library-size shift, dropout-rate
  shift, cell-type-specific effects) and by reporting whether conclusions differ.

## Kept: D3 — Real multi-donor trajectory, plus null calibration

Two complementary arms on genuinely real batch structure, where no ground truth
exists and the measurement has to be indirect.

- **Real batches**: `datasets/palantir_cd34/` (3 donors, `palantir_pseudotime`
  and `palantir_branch_probs` per cell) and `Immune_ALL_human.h5ad` (10 batches,
  HSPC→erythrocyte lineage with the precomputed `dpt_pseudotime` scIB used).
  Score with `scib.metrics.trajectory_conservation` plus per-donor
  within-batch pseudotime preservation — *within* a single donor there is no
  batch effect, so any change to that donor's internal trajectory ordering is
  attributable to the method, not to correction of a real effect. This
  recovers a ground-truth-like contrast without a counterfactual.
- **Null calibration** (Antonsson & Melsted's design): random pseudobatches over
  the same real data, no batch effect present. Everything a method changes here
  is pure artifact. This anchors the "biology removed" axis at a point where the
  correct answer is exactly zero, and cross-checks whether the percentages from
  D1/D2 are in a believable range.
- **R7 folds in here**: the CellANOVA hold-out-control trick (relabel one control
  donor as a fake treatment; any "signal" found is false by construction) applies
  directly to `code/cellanova/data/t1d_example.h5ad` and to the Palantir
  replicates, and is a cheap second null.

---

## Rejected, with reasons

- **R1 — scATAC / multi-omic.** The hypothesis is about differentiation
  trajectories in scRNA-seq. scIB shows feature-space choice dominates method
  choice for ATAC, so results would answer a different question. Also skipped
  the ATAC files during download.
- **R2 — build a new recovery method.** The hypothesis is a measurement claim
  ("no method dominates"), not a method proposal. CellANOVA already occupies
  this space and is cloned; competing with it would replace the question rather
  than answer it.
- **R3 — dyngen / SymSim.** Mechanistically better simulators with gold-standard
  trajectories, and genuinely the right tool. **Rejected purely on
  feasibility**: both are R-only, and there is no R runtime and no root
  (`apt-get install r-base-core` → permission denied; conda disallowed). This is
  the most costly rejection in the list — `splatpy` is the substitute, and its
  unvalidated status must be stated as a limitation. *If an R runtime becomes
  available, reopen this immediately.*
- **R4 — cross-species.** Species differences are simultaneously the batch and
  the biology; there is no ground truth that separates them, so it cannot serve
  a quantitative claim. (`Immune_ALL_hum_mou.h5ad`, 4.3 GB, deliberately not
  downloaded.)
- **R5 — within-method Pareto sweep.** Li et al. 2022 already traced scVI's
  Pareto front via multi-task learning. As a standalone direction it duplicates
  published work; as a *component* of D1 it is essential, because a
  "no method dominates" claim is weak if each method is represented by one
  default-parameter point. Folded into D1.
- **R6 — spatial.** Different modality, off-hypothesis.
- **R7 — disease/perturbation signal preservation.** Well-motivated (it is
  CellANOVA's and Tyler et al.'s focus) but the hypothesis names differentiation
  trajectories specifically. Retained only as the hold-out-control null inside D3.

## Revision rule

This ranking is fixed unless new evidence invalidates it. The two triggers that
would justify revisiting: (a) an R runtime becomes available → reopen R3; (b) D1
and D2/D3 disagree materially on the magnitude of biology removed → the
simulation model itself becomes the object of study, and the ranking must be
rewritten in `STATE.md` with the reason.
