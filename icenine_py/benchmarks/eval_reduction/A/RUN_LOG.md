# RUN_LOG: Task A, evaluation ledger and the residual C++ vs Python gap (2026-10-13)

Task branch `feature/eval-reduction-A` (parent `feature/eval-reduction`). Plan: `~/.claude/plans/eval-reduction.md`,
Task A. All Python through `uv run` from `icenine_py/`; all C++ in a scratch copy of the repository (`Src/` in the
repository is never touched; the scratch directory is not kept, the patches are). Clean full-sample images only
(Phase D clean stack, Q-max 8 configuration), no realistic images, no extra seeds. Timings are labelled
"single-worker" (one process, torch threads 1, `preflight.require_quiet` passed at its start) or "contended".

## 1. What changed in the repository, and why

| File | Change | Why |
|---|---|---|
| `icenine/cost_functions.py`, `icenine/_rasterize.c` | opt-in `cpp_edge_rules` (config key `CppEdgeRules`, default off): omega below the first frame has no frame, a negative detector coordinate is pixel -1, a pixel-less triangle with an in-bounds vertex lights the detector; helpers `_omega_bins`, `_pixel_coordinate`; C extension constant `HAS_VERTEX_LIT` and bit 2 of the `cpp_rule` argument | the three differences Step 0 found between the Python and C++ cost for one fixed orientation |
| `icenine/batched_screen.py` (new) | `ScreenEvaluator`: batched exact screen (`overlaps_exact`) and C++ omega interpolation (`overlaps_interp`, `interp_peaks`) | Step 3: `BatchedScreen`, `CppOmegaInterp` |
| `icenine/orientation_search.py`, `icenine/reconstructor.py`, `icenine/config_file.py` | `run_discrete_search_spaced(screen=, screen_interp=, screen_chunk=)`; `_screen_kwargs`; `_cost_flags` carries `cpp_edge_rules`; config keys `CppEdgeRules`, `CppOmegaInterp`, `BatchedScreen`; `SearchParameters` fields | wiring; every switch off = the old call, byte for byte |
| `scripts/eval_reduction/` (new) | `parity_inputs.py`, `parity_common.py`, `parity_compare.py`, `parity_summary.py` (Step 0); `eval_ledger.py` (Step 1); `cpp_ledger_parse.py`, `ledger_attribution.py` (Step 2); `screen_check.py`, `interp_vs_cpp.py`, `screen_timing.py`, `e0_interp_summary.py` (Step 3) | the measurements |
| `scripts/phase_d/e0_cpprule.py` | arms `edge`, `interp`, `lenint` | E0 clean for the new switches |
| `tests/test_cpp_edge_rules.py`, `tests/test_batched_screen.py`, `tests/test_eval_reduction_tools.py` (new) | 4 + 7 + 2 tests | switches off = old behaviour; switches on = what they say |
| `benchmarks/eval_reduction/A/` | patches of the scratch C++ builds, results, this log | evidence |

Not committed: the C++ scratch builds. Patches: `step0/cpp_parity_patch.py.txt` (fixed-orientation driver and per-peak
dump), `step2/cpp_ledger_patch.py.txt` (stage counters of the C++ no-start search). Both apply to a fresh copy of
`CMakeLists.txt`, `Src/`, `XDM++/`, `CostFunctionBenchmark/` with `python patch.py` from the scratch root and build
with `cmake -DCMAKE_BUILD_TYPE=Release . && make -j8` (about 11 s).

## 2. Runs (commands from `icenine_py/`, all clean images)

| Step | What | Command | Wall | Label |
|---|---|---|---|---|
| 0 | C++ scratch build with the parity driver | `python patch.py && cmake -DCMAKE_BUILD_TYPE=Release . && make -j8` in the scratch copy | 11 s | |
| 0 | inputs: 2,230 orientations of three voxels | `parity_inputs.py --out DIR` | 2 s | |
| 0 | C++ driver, 1,077 / 577 / 576 records | `PARITY_IN=... PARITY_OUT=... IceNine r recon_V.config` from `Examples/Example2.ManyGrains` | 14 to 22 s each | single process |
| 0 | Python comparison, 4 Python variants per orientation | `parity_compare.py --dir DIR --out benchmarks/eval_reduction/A/step0` (absolute `--out`: the script changes directory) | about 60 s | |
| 0 | summary | `parity_summary.py`, `pixelless_check.py` | seconds | |
| 1 | pilot, `cpp` arm, voxel 15901 | `eval_ledger.py run --voxel 15901 --tag pilot_cpp --opt cpp_overlap_rule=1` | 117.7 s | single-worker |
| 1 | ledger, 6 voxels x 4 arms (base, cpp, r2, cpp_edge) | `ledger/run_all.sh` | about 62 min | back to back, 16 of 24 labelled contended (load 1.5 to 1.6) |
| 2 | C++ counters, 6 voxels | `LEDGER=1 IceNine r recon_V.config` (scratch build with `cpp_ledger_patch.py.txt`) | about 20 s each | single process |
| 2 | attribution | `cpp_ledger_parse.py`, `ledger_attribution.py`, `hit_ratio_check.py` | seconds | |
| 3 | pass sets, voxel 15901, four ways, edge rules off and on | `screen_check.py --voxel 15901 [--edge]` | 40 s each | contended (ledger running) |
| 3 | Python interpolation against the C++ fixed pairs | `interp_vs_cpp.py` | 20 s | |
| 3 | ledger with the switches, 6 voxels x 4 arms | `ledger/run_switches.sh` | about 50 min | back to back, 8 of 24 labelled contended |
| 3 | E0 clean, arms interp, edge, lenint, 10 workers | `ledger/../e0/run_e0.sh` (`e0_cpprule.py run --arm A --workers 10 --msin 0 --tag a --variant clean`) | about 2 to 5 min per arm | contended multi-worker |
| 3 | level-0 screen microbenchmark | `screen_timing.py --out step3` | 16 s | single-worker, preflight passed (load 1.22) |
| 3 | quiet repeat of two arms on two voxels | `ledger/run_timing.sh` | 8 min | single-worker, all four passed |
| all | full suite, switches off | `uv run pytest tests/ --junitxml=...` | 303 s | 889 passed, 37 skipped, 1 deselected |
| all | full suite with the defaults flipped temporarily (twice), then restored | same | 336 s, 303 s | `default_flip_would_change.md` |

## 3. What went wrong, and what was redone

* The first `parity_compare` run matched too few peaks: reflections of one family (hkl and 2h 2k 2l) share a direction,
  and the C++ output carried no direction for peaks without a frame. Both sides now carry direction and |G|, matching
  uses all three, and every Python peak is matched (0 unmatched).
* The first Python probe returned empty results because the phase index was 0; the voxels carry phase 1.
* `parity_compare` and `eval_ledger` change the working directory to `Examples/Example2.ManyGrains`; a relative
  `--out` or `--dir` wrote into that directory (one ledger file, removed and moved; the scripts now resolve the path
  first). The first ledger chain was killed after one run and restarted for this reason.
* `sed -i` is not the GNU one on this machine (several edits were done with Python instead).
* The ledger runs stored every quick-MC call (about 2.3 MB per run, 195 MB in all); `eval_ledger.py shrink` stores the
  quick calls as a summary (min, median, max evaluations, blocks, stop codes) and the FindOptimal and variance calls in
  full (688 KB in all). The raw C++ driver output (21 MB) is not kept.
* `run_timing.sh` first waited for a load below 1.0, which this machine never reaches while the editor and the agent
  run (idle load 1.2 to 1.4); it was restarted with 1.4.
* The Write tool refused to create `docs/findings_eval_reduction_2026-10.md` (a rule against report files); the
  findings are in MIGRATION_HISTORY "Task A" and the two open questions of
  `docs/findings_finisher_mc_study_2026-10.md` instead. Task G creates the consolidated findings document.

## 4. Could be done better

* The back-to-back ledger runs read a 1-minute load of 1.5 to 1.6 from their own previous run and are labelled
  contended; a cooldown between runs would have given clean labels at the cost of about 40 minutes. Only the two
  voxels of the quiet repeat carry a clean label for the whole-search timings.
* FindOptimal and VarianceMinimizing are compared on six voxels, one seed each; their counts are random, so only the
  mean is comparable.
* Which of the three edge conventions makes FindOptimal converge on the first candidate was not separated
  (the hit ratio at the exact truth is 1.0 in every variant).
* The batched screen is checked against the unbatched screen on 43,974 pairs of one voxel (identical) and on the
  fixed test set; the whole searches of all six voxels have identical evaluations and final errors with and without
  BatchedScreen (the exact screen), and with the interpolation at one and at 32 FZ points per call.
* The E0 comparison used the Task 0 baselines (`on`, `off`) from the cache of the centroid run; the default code path
  is unchanged (checked on voxel 15901 and by the golden tests), so they were not re-run.

## 5. Task A fixes (2026-10-13, `feature/eval-reduction-A-fixes`)

| What | Command | Wall | Label |
|---|---|---|---|
| attribution re-run (adds `random_ratio_range_*`, `l0_candidates_cpp_vs_cpp_edge`; existing keys unchanged) | `ledger_attribution.py --ledger benchmarks/eval_reduction/A/ledger --cpp benchmarks/eval_reduction/A/step2 --arms cpp cpp_edge cpp_edge_interp cpp_edge_interp_batched --out benchmarks/eval_reduction/A/step2` | seconds | |
| Step 0 summary re-run (99.96% formatting; `summary.json` unchanged) | `parity_summary.py --dir <parity inputs> --out benchmarks/eval_reduction/A/step0` | seconds | |
| peak RSS, voxel 15901, `BatchedScreen` off then on | `screen_memory.py run --voxel 15901 --batched 0` / `--batched 1`, then `summary` | 173 s, 145 s | one process, contended (load 4.7, 3.0; the second overlapped the suite) |
| full suite with only `BatchedScreen` flipped | `uv run pytest tests/ --junitxml=...` | 282 s | 887 passed, 2 failed (default assertion, `test_profiling.py` instrument count), 37 skipped, 1 deselected; flip reverted |

## 6. Task A-exact (2026-10-13, `feature/eval-reduction-A-exact`)

`CppEdgeRules` lights a detector exactly as C++ (overlap > 0 or an in-bounds vertex); extension rebuilt
(`HAS_VERTEX_LIT` = 2). Outputs in `a_exact/`. Commands from `icenine_py/`; `DIR` is the scratch directory of the Step 0
inputs (`parity_inputs.py --out DIR`) and C++ outputs; logs in the scratch directory, not kept.

| What | Command | Wall | Result |
|---|---|---|---|
| C++ parity driver rebuilt in a scratch copy | `python patch.py && cmake -DCMAKE_BUILD_TYPE=Release . && make -j8` (patch `step0/cpp_parity_patch.py.txt`) | not recorded | first A-exact commit; C++ `Src/` unchanged by this task |
| C++ driver rerun, 1,077 / 577 / 576 records | `PARITY_IN=... PARITY_OUT=... IceNine r recon_V.config` from `Examples/Example2.ManyGrains` | not recorded | first A-exact commit; outputs reused by the rows below |
| Python comparison at the A-exact code | `parity_compare.py --dir DIR --out <scratch>/rerun` (absolute `--out`) | 41 s | `results.json.gz` byte-identical to the committed `step0/results.json.gz` (`cmp`; same SHA-256, `a_exact/parity_rerun.json`), so `step0/summary.json` and the Task A tables stand |
| slow path against the C path, dense images | `slow_path_check.py --out <abs>/benchmarks/eval_reduction/A/a_exact/slow_path_check.json` | 4 s | key `plain`: 48 evaluations, 0 differ in any field |
| same, `--stub` (dark images, slow-path triangle counter returns 0 overlap, 5 lit pixels) | `slow_path_check.py --stub --out <same file>` | 1 s | key `stub`: complete cost with the switch, peak counts equal 12 of 12; without the switch 0 of 12 (peak-on-detector sums 1263 slow against 766 C); the screen never calls the triangle counter, so it does not test the lit rule |

Review fixes (`.claude/reports/20261013-2304-review-A-exact.md`): these outputs were not saved in the first A-exact
commit; the screen stub sum was wrongly cited as evidence for the lit rule in MIGRATION_HISTORY (rewritten); a unit
test of the slow-path rule (`test_slow_path_lit_rule_with_a_stubbed_triangle_counter`) was added.
