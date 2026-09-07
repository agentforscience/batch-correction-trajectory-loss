# Research State

- Current phase: `None`
- Pipeline completed: `False`

## Previous phases

resource_finder (succeeded), experiment_runner (failed)

## Current phase context

- Phase: `experiment_runner`
- Status: `failed`
- Started: `2026-09-06T08:32:37.839074Z`
- Next steps:
  - Validate the report and experimental artifacts before finalizing.

## Workspace check

- Expected: `/workspaces/how_much_biology_do_batch_corr_20260906_074942_ba309bd5`
- Actual: `/app`
- Directory usable: `True`
- Current process matches workspace: `False`

## Output validation

- Valid: `False`
- Expected: `REPORT.md`
- Missing: `REPORT.md`
- Outside workspace: None

## Agent notes

<!-- NEURICO_AGENT_NOTES_START -->
### resource_finder
<!-- NEURICO_AGENT_NOTES_START:resource_finder -->
### Phase 1 (resource_finder) — COMPLETE

**Environment.** Fresh `uv` venv at `.venv`, isolated with `[tool.uv] package = false`
in `pyproject.toml`. Activate with `source .venv/bin/activate`. GPU: RTX A6000,
torch 2.14+cu130 CUDA=True, but **~42 GB of 49 GB already used by another process**.
32 CPUs, 503 GB RAM, 737 GB disk.

**Hard constraint that shaped everything: no R runtime, no root, conda disallowed.**
This eliminated splatter, dyngen, SymSim, scDesign2, Seurat v3, the R kBET,
CellMixS and BatchBench. Do not spend time retrying these.

**Artifacts on disk**
- `papers/` — 25 PDFs, all verified openable; `papers/README.md` annotates each;
  `papers/pages/` has 3-page chunks for the 5 deeply-read papers;
  `papers/abstracts_only/README.md` covers 4 unreachable ones.
- `datasets/` — 2.9 GB, 8 usable `.h5ad`, git-ignored, `datasets/download.sh`
  reproduces all of it, `datasets/README.md` has verified obs/label breakdowns.
- `code/` — 15 cloned repos + `code/splatpy/` (written here); `code/README.md`.
- `literature_review.md`, `resources.md`, `planning.md`.
- `tools/` — `search_epmc.py`, `fetch_paper.py` (reusable retrieval scripts).

**Key finding that defines the experiment.** Every published measurement of
trajectory preservation is relative to the *unintegrated* data, not to truth.
scIB's `trajectory_conservation` = `(spearman(dpt_post, dpt_pre)+1)/2`, and under
a real batch effect `dpt_pre` is itself corrupted — so no percentage of "biology
removed" can be derived from it. `code/splatpy/paired_simulation()` closes this:
it runs the Splat model twice on the same RNG stream with the batch factors
toggled, producing the *same cells* (verified: identical Group, Step, ExpLibSize)
differing only by batch. The batch-free arm is an exact counterfactual, which
makes "% of trajectory variance removed" computable rather than proxied.

**Direction budget (3 kept, 7 pruned — full scoring in `planning.md`).**
- D1 paired counterfactual simulation (`splatpy`), sweeping batch strength,
  confounding, imbalance, and per-method hyperparameters.
- D2 injected-batch spike-in into real trajectory data (paul15 /
  pancreas_endocrinogenesis / setty_bonemarrow), with the injected effect
  estimated from a real protocol contrast — the "real spike-in" of the title.
- D3 real multi-donor trajectories (Palantir CD34+ 3 donors, Immune_ALL_human
  10 batches) plus the Antonsson null-pseudobatch calibration.
- Pruned: scATAC/multi-omic, building a new recovery method, dyngen/SymSim
  (feasibility only — reopen if R appears), cross-species, spatial; and two
  folded in rather than dropped (within-method Pareto sweep → D1;
  disease-signal preservation → D3).

**Methods verified end-to-end** on splatpy output: ComBat 0.1s, BBKNN 11s,
Scanorama 0.3s, Harmony <0.1s, scVI 8s (GPU), scANVI 1.3s (GPU), scGen 1.5s (GPU).
**Unavailable: MNN/fastMNN** (`mnnpy` won't build on py3.12) and **LIGER**
(`pyliger`→`louvain` build failure). scIB places LIGER at the batch-removal
extreme, so the measured frontier under-samples that end — state as a limitation.

**Traps to avoid (each cost time here)**
1. `harmonypy` 2.0 returns `Z_corr` as `(n_cells, n_dims)` — do **not** transpose.
2. `scib.metrics.kBET` needs rpy2+R. Use `scib_metrics.kbet` (pure Python).
3. PyPI `scgen` 2.1.0 imports the removed `scvi._compat`. Use the editable
   install from `code/scgen`, which carries two applied patches (inference dict
   keys for scvi-tools ≥1.2; `AnnData.concatenate` → `anndata.concat`).
4. `cellanova` will not `pip install -e .` (flat layout). Install its bundled
   wheel: `code/cellanova/dist/cellanova-0.1.0-py3-none-any.whl`.
5. Palantir replicates have different gene sets and **no raw counts layer** —
   intersect `var_names`, and scVI/scANVI cannot run on them as-is.
6. jax is CPU-only; metrics run but slowly. `uv add "jax[cuda12]"` if needed.

**Two unresolved items for Phase 2, both surfaced by the feasibility run**
(numbers in `resources.md`):
- **Harmony did nothing at `batch_fac_loc >= 0.2`** — logged "Converged after 1
  iteration", ASW-batch unchanged (0.408 → 0.408), while it worked normally at
  0.1 (0.076 → 0.002). Probably `theta`/`nclust`/`max_iter_harmony` on only
  1,200 cells, not a real method failure. **Resolve before reporting.**
- **`max|spearman(PC_k, Step)|` is too insensitive** as a trajectory metric — it
  stayed ~0.85 across the full batch-strength sweep. Replace with variance
  attributable to `Step` over the whole embedding, or pseudotime recovery
  against the oracle embedding.

**Known risk.** `splatpy` is faithful to splatter's R source function-by-function
but has **not** been validated against real `splatSimulate()` output (no R).
Read `papers/crowell2023_shaky_foundations_simulation.pdf` before making realism
claims. This is why D2 and D3 are not optional: the simulation result is only
credible if the real-data arms agree with it.

**Next phase: experiment_runner.** Concrete next steps:
1. Read `planning.md` (directions), then `code/README.md` (what runs, what does not).
2. Resolve the Harmony convergence question and settle the trajectory metric.
3. Build the D1 harness: paired simulation → 7 methods → batch + biology metrics,
   several seeds, sweeping batch strength and confounding.
4. Port the same scoring to D2 and D3; report the two Pareto axes separately.
<!-- NEURICO_AGENT_NOTES_END:resource_finder -->

### experiment_runner
<!-- NEURICO_AGENT_NOTES_START:experiment_runner -->
Update this section at the end of the `experiment_runner` phase.
<!-- NEURICO_AGENT_NOTES_END:experiment_runner -->

<!-- NEURICO_AGENT_NOTES_END -->
