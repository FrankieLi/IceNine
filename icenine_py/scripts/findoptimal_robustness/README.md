# FindOptimal robustness study (2026-10-06)

Why does the multi-level reconstruction (`AdaptiveVoxelReconstructor.reconstruct_voxel`) return a
CSL relative of the truth, and what fixes it? Report: `docs/findoptimal_robustness_report.md`;
results: `benchmarks/findoptimal_robustness/`; plan and summary: `MIGRATION_HISTORY.md`.

All commands run from `icenine_py/`, with `uv run python`. Raw per-run caches go to
`scripts/findoptimal_robustness/cache/` (gitignored, about 300 MB); every step is resumable.

| Script | What it does |
|---|---|
| `common.py` | voxel set, per-case detector images (as `findoptimal_sweep.py` experiment A), worker setup, Wilson intervals |
| `csl.py` | cubic CSL table (Sigma <= 29, 21 entries), distinct relatives of an orientation, Brandon classification, shared reflections |
| `e0_run.py` | `select` the 200 voxels (+ spares), `run` instrumented full reconstructions (recorder hook), 3 seeds, clean and realistic |
| `analyze_e0.py`, `analyze_dependence.py` | where the truth is lost (S1/S2/S3), wrong rate per seed/voxel, dependence on r_perp, boundary, shared reflections |
| `f1_run.py` | F1: CSL relatives of the final answer, quick MC, FindOptimal on the best (post hoc; `--source` for another run set) |
| `fixes_run.py` | F1b, F2a/b, F3a/b/c: full re-runs with a search knob changed (`--subset-wr`: wrong runs + equal right runs) |
| `summarize_fixes.py`, `f1_remaining.py` | wrong rate (Wilson CI), cost and what remains wrong |
| `features.py` | one-pass candidate features (hit fractions by family/detector, CSL-shared vs non-shared reflections) |
| `e2_dataset.py`, `e2_models.py` | candidate dataset (harvested + synthetic), LR / gradient-boosted trees, voxel-disjoint folds, ROC, learning curve |
| `e2_endtoend.py`, `e2_summary.py` | classifier as pruning reranker (`rerank`) and as final chooser (`final`) |
| `plots.py` | the four plots |
| `run_chain*.sh` | the job chains that were run |

Reproduction order: `e0_run.py select`, `e0_run.py run`, `analyze_e0.py`, `f1_run.py run`,
`e2_dataset.py run`, `analyze_dependence.py`, `fixes_run.py run --fix F1b` (and `--fix F2a ... --subset-wr`),
`e2_models.py run --save-models`, `e2_endtoend.py final|rerank`, `f1_run.py run --source e2_rerank_GBT`,
`summarize_fixes.py`, `e2_summary.py`, `plots.py`.

Test: `tests/test_findoptimal_robustness.py` (CSL generator, folds, feature extractor) and the
bit-identity test of the recorder hook in `tests/test_findoptimal_refactor.py`.
