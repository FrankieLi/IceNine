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
  so that both definitions are reported) and with link radius 0.5 x diameter (largest component 18 to 61 cliques).
* C0 to C3 are built as the plan says.

## 4. Step 2: what changed (in addition to the table in section 1)

| File | Change | Why |
|---|---|---|
| `icenine/orientation_search.py` | `SearchCandidate.clique` (FZ point index, default -1, set in `run_discrete_search_spaced`); `SearchParameters.candidate_dedup` (default `"off"`, parsed in `__post_init__`, `from_config`) | the clique of a candidate; the key |
| `icenine/config_file.py` | `candidate_dedup` field, key `CandidateDedup` (validated with `DedupSpec.parse`) | config |
| `icenine/reconstructor.py` | `dedup = DedupSpec.parse(...)` once per voxel; C0/C1 on `candidates` after the "discrete" record and before the quick-MC branch (`_dedup_before_quick_mc`); C2 after the sort and before `n_keep`; `n_keep` from the discrete list (`quarter=orig`) or the shorter list; additive keys `dd_*` in the "quick_mc" record with the key on. With the key off every added branch is skipped | hooks |
| `scripts/eval_reduction/eval_ledger.py` | `--opt NAME=VALUE` splits on the first `=` only | a spec contains `=` |
| `.gitignore` | `icenine_py/benchmarks/eval_reduction/C/cache/` | raw arrays (about 200 MB) |
| `tests/test_candidate_dedup_levels.py` (new) | config key, clique index, key off bit-identical and no `dd_` keys, C0 / C1 / C2 / batched quick MC on the ThreeVoxels toy problem with eight extra FZ points | level loop |

Untouched, as the brief says: `_batched_quick_mc`, `batched_cost.py`, `batched_quick_mc.py`, `cost_functions.py`, `_rasterize.c`, all physics code, every golden test.

## 5. Step 3: runs (commands from `icenine_py/`)

| What | Command | Wall | Label |
|---|---|---|---|
| E0 clean, 14 arms (`arms.txt`), 200 voxels each | `run_e0_arms.sh` (`candidate_arms.py run --arm A --spec S --workers 10`) | about 95 s per arm | contended multi-worker; the saved `run_meta.json` of `base` and `c0a_k1` was overwritten by a later no-op call |
| six-voxel ledger, 14 arms x 6 voxels = 84 runs | `run_ledger_arms.sh` (`eval_ledger.py run --voxel V --tag A --opt candidate_dedup=S`) | 8 to 18 min per arm | one process at a time; the preflight inside the tool passed in 15 of the 84 runs (the 1-minute load was above 1.5 from the previous run), so 69 are labelled contended; no timing claim |
| 24 more voxels of the 30-voxel set, 6 arms (144 runs) | `run_full30.sh` | about 50 min | 8 processes at once, contended |
| E0 sanity check: the uninstrumented `e0_cpprule.py --arm off --msin 0 --tag c_ref --variant clean` against the instrumented `base` arm | one run, 4 workers | 18 min (ledger running; `step3/run_meta_e0_uninstrumented_off.json`) | identical error and evaluation count on 200 of 200 voxels |
| full suite | `uv run pytest tests/ --junitxml=...` | 306 s | 996 passed, 37 skipped, 1 deselected |
| summaries | `candidate_arms.py summary`, `candidate_summary.py`, `candidate_arms_figure.py`, `candidate_clique_truth.py`, `collect_tables.py --root benchmarks/eval_reduction/C --prefix evalC_` | seconds | |

Preflight files: `step3/preflight_e0.json`, `ledger/preflight_start.json`, `full30/preflight_start.json` (the E0 one was taken while the machine was still loaded by the previous tool call; load 5.2 to 9.1).

## 6. Decision points (non-blocking, taken by the plan's default)

* After Step 1: C0 to C3 built; C4 not (outputs within 0.5 D: 23.1% and 25.6% of the full-sample outputs, not over half). The one-diameter merged clique was run as the plan words it
  and also at 0.5 D, because it is degenerate at level 0 (section 3).
* Extra arms beyond the plan's list, so that the plan's design choices can be read: `quarter=filtered` for C1 at 0.5 D and for C3 (the plan does not say whether the freed slots are refilled for C0
  and C1; `quarter=orig` refills them and is the default of the key), and the 30-voxel extension of section 5.
* After Step 3: the arm named for Task G by the plan's rule (fewest evaluations with no broken E0 voxel and no more wrong voxels than the baseline on the six) is **the baseline**: no arm has zero broken E0 voxels
  (`step3/summary.json`, `_choice`). Evidence the owner may weigh against that rule is in the MIGRATION_HISTORY subsection (`clique=cpp:1` saves 23.9% of the evaluations on 30 full-sample voxels with 2 wrong against 1; `pre=0.5`
  saves 8.9% with 0 wrong; on E0 they break 24 and 30 voxels).

## 7. What went wrong, and what was redone

* `eval_ledger.py --opt` split on every `=`; a spec such as `pre=0.5` failed (one pilot run lost). Fixed to split once.
* The first ledger chain waited up to 90 s per run for a quiet load and was killed after four base runs and restarted with a 16 s wait (the finished runs were kept; the runs are deterministic given the code and seed).
* A pre-commit check rejected files with lines over 100 characters inside `# fmt: skip` blocks (twice); lines wrapped by hand.
* `audit_numbers.py` leaves 52 numbers of the MIGRATION_HISTORY subsection unmatched: all are thousands-separated cells ("46,031") that the generated tables contain verbatim (the tool does not read comma-formatted
  source numbers); they were synced from `doc_tables_all.md`, not typed. The other numbers match (`step1/summary_*.json` has the pooled shares under `shares`, `step3/summary.json` the derived savings under `_derived`).
* The plan's findings-document deliverable was not written as a separate file: the harness refused to write a findings `.md` for this subagent. The findings are in the MIGRATION_HISTORY subsection (with the generated tables) and the report.

## 8. Could be done better

* The paired "broken / fixed" counts mix the effect of an arm with the effect of any change of the random draw (a different number of quick-MC runs shifts the later draws). Separating them needs extra seeds, which are not allowed; a
  per-candidate random stream (as `BatchedQuickMC` draw option B has) would remove the shift for the filters before the quick MC.
* E0 hardly measures evaluations (92.2% of them are the screens and re-costs); a full-sample accuracy set larger than 30 voxels, or the BFS pilot, would decide between the baseline and `clique=cpp:1` / `pre=0.5`.
* The clique filter keeps the best-scored candidate of a clique; Step 1 shows the truth-nearest one is not always that (9 of 30 E0 runs at level 0). A score that ranks within a clique better than the discrete cost was not tried.
* Evaluations are counted, not timed: a quiet-machine timing of the arm that is finally chosen, with and without `BatchedQuickMC`, is open.
