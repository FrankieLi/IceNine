# RUN_LOG: Task C, candidate overlap and de-duplication (2026-10-14)

Task branch `feature/eval-reduction-C` (parent `feature/eval-reduction`). Plan: `~/.claude/plans/eval-reduction.md`,
Task C. All Python through `uv run` from `icenine_py/`. Clean images only, no realistic images, no extra seeds. The
library default path (every new key off) is the baseline; `BatchedQuickMC` and `SearchDevice` stay off in every arm.
The raw per-run arrays are in `benchmarks/eval_reduction/C/cache/` (gitignored); everything quoted comes from the
`summary_*.json` and `tables_*.md` files beside this log.

## 1. What changed in the repository, and why

| File | Change | Why |
|---|---|---|
| `icenine/candidate_dedup.py` (new) | `neighbour_lists`, `greedy_nms`, `component_labels`, `merged_clique_labels` (symmetry-reduced cKDTree clustering); `DedupSpec`; `keep_before_quick_mc` (C0, C1), `distinct_order` (C2) | the fast clustering the plan asks for, and the three stage functions |
| `tests/test_candidate_dedup.py` (new) | clustering against `csl.dedup_orientations` and a brute-force reduced misorientation (random sets, cubic FZ-boundary sets, a 4-fold group, no symmetry), the spec parser, the stage functions | correctness |
| `scripts/eval_reduction/candidate_probe.py`, `candidate_overlap.py`, `candidate_overlap_figure.py` (new) | step 1 measurement from the existing recorder hooks, no library change | step 1 |
| `scripts/eval_reduction/candidate_arms.py` (new) | step 3, E0 arms | step 3 |

## 2. Step 1: runs (commands from `icenine_py/`)

| What | Command | Wall | Label |
|---|---|---|---|
| pilot, voxel 15901, C++ rule | `candidate_overlap.py run-full --voxel 15901 --rule cpp` | 95 s, 237,739 evaluations | single process, load 1.2 |
| 6 full-sample voxels x lenient and C++ rule | `benchmarks/eval_reduction/C/run_step1.sh` (6 processes at once) | 83 to 157 s each | contended (counts only; no timing claim) |
| first 30 E0 voxels x lenient and C++ rule | same script, `run-e0 --rule R --workers 10` | 3 to 4 s per voxel | contended multi-worker |
| analysis | `candidate_overlap.py analyze --set full` and `--set e0` | seconds | |

The clique labels recovered from the "fz" recorder events reproduce the recorded candidate list in every run (the
`clique_ok` flag is true for all levels of all runs). The evaluation counts of the pilot match the plan's ledger numbers
(8,110 / 2,487 / 697 / 204 candidates per level; 137,976 quick-MC evaluations over four levels).

## 3. Step 1: decision point (non-blocking, default taken)

* C4 (a shorter quick MC on duplicates) is built only if over half of the quick-MC outputs are duplicates. They are not:
  within 0.5 x the level diameter 23.1% (C++ rule) and 25.6% (lenient) of the full-sample outputs at levels 0 to 2, 13.0% and
  16.3% of the E0 outputs; within 0.1 deg only 0.06% and 0.03%. **C4 is not built or run.**
* The merged clique (definition ii) as the plan words it (link radius 1 x diameter) is degenerate at level 0: 99.5%
  (C++ rule) and 100% (lenient) of the cliques overlap another clique and the largest component holds every clique of the
  run. The arms with `merged` are therefore run twice: with the literal link radius 1 x diameter (expected to fail, kept
  so that both definitions are reported) and with link radius 0.5 x diameter (largest component 43 to 61 cliques).
* C0 to C3 are built as the plan says.
