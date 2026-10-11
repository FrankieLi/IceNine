# RUN_LOG: level-0 candidates, Step 1 (trace), 2026-10-10

Task branch `feature/finisher-mc-study-level0-candidates` (parent `feature/finisher-mc-study`).
Plan: "Phase D task: level-0 candidates" in the finisher/MC study plan. This log covers Step 1 only
(no fix to `icenine/` search behaviour; no E0; no multi-worker job). Everything ran single-process
on the shared laptop (load average 1.4 to 2.1 at the start): timings below are "contended" and are
not benchmark timings.

## 1. What changed in the repository, and why

| File | Change | Why |
|---|---|---|
| `icenine/orientation_search.py` `run_discrete_search_spaced` | new optional argument `recorder` (default `None`) | per-FZ-point trace of the screen, re-cost and spacing filter; with `None` the code path is the old one (no extra work beyond `is not None` tests) |
| `icenine/reconstructor.py` `AdaptiveVoxelReconstructor` | attribute `discrete_recorder` (default `None`), passed only when set | hands the hook to the search without touching the signature seen by callers that monkeypatch `run_discrete_search_spaced` (`seed_diag.py`) |
| `tests/test_discrete_recorder.py` | new, 2 tests | default and recorded results identical; records consistent |
| `scripts/phase_d/level0_trace.py` | new: `py`, `cpp`, `compare`, `s3` | Python trace, C++ log parser, S4/S1/S2/S3 comparison, S3 on C++'s own inputs; `py --cpp-rule --centroid` is a read-only probe that swaps in the C++ rules at run time |
| `scripts/phase_d/level0_summary.py` | new | `summary.json` and `tables.md` (doc_tables) |
| `benchmarks/phase_d_seed_diag/level0/` | traces, comparison JSONs, tables, C++ patch as `.txt` | evidence |

The library behaviour is unchanged. The golden-style tests and the full suite were run (section 5).

## 2. C++ side (manual steps)

1. The earlier session's scratch directory still existed (not kept by design, per MIGRATION_HISTORY); I
   copied only the one-voxel mic and the config for 15901 from it
   (`run/one_15901.mic`, `run/recon_15901.config`) and rewrote the two absolute paths in the config with `sed`.
   This is the Q-max 8 config used by `seed_diag` (`MaxQ 8`, `MaxLocalResolution 3`, 200 MC steps, 2 restarts,
   30 candidates, `EtaLimit 86`, `SampleSymmetry Cubic`), clean full-sample images
   (`ScatteringData_PhaseD/full/clean/500Grains.sim`).
2. Scratch copy of the repo sources outside the tree: `CMakeLists.txt`, `Src/`, `XDM++/`,
   `CostFunctionBenchmark/` into the session scratchpad (`cpp_l0/`). `Src/` in the repo was not touched.
3. Patch 1 (`cpp_debug_patch_level0.py.txt`, applied with `python patch.py` from the scratch root): prints for
   the first `GetSpacedCandidates` call (level 0) per FZ point the peak_overlap of every local-grid point
   (`L0O`), the screened indices, their re-cost 1 - confidence and the kept indices (`L0P`), and per call the
   header and totals (`L0H`, `L0S`). Gated by env `L0TRACE=1`.
4. Patch 2 (hand edits in the scratch copy; the full result is `cpp_debug_full_scratch.diff.txt`, a diff of the
   scratch `Src/` against the repo `Src/`): per-peak dump for the first three FZ points at j = 0 (`L0Y` per
   peak, `L0W` window counts inside `PixelBasedPeakOverlapCounter`), the centre and local matrices (`L0M`), and
   exact versus interpolated omegas for j = 0..2 (`L0X`). Gated by env `L0X` plus a flag that is true only inside
   the first `GetSpacedCandidates` call.
5. Build: `cmake -DCMAKE_BUILD_TYPE=Release .` (1.3 s) and `make -j8` (10.7 s wall, 17:19) in the scratch copy.
   Seven builds in all (patch, one compile error fixed, then rebuilds after added prints), each about 10 s.
6. Run: from `Examples/Example2.ManyGrains`, `L0TRACE=1 <scratch>/IceNine r <scratch>/run/recon_15901.config`,
   single process, 28.4 s wall (17:19:48 to 17:20:17, contended). The C++ level-0 candidate count 8,151 equals the
   stored `mc_parity_stages` value, so the scratch build reproduces the earlier one.
7. Parsed into `cpp_v15901.json` with `level0_trace.py cpp --log <scratch log>`.

## 3. Python side

| Run | Command (from `icenine_py/`) | Wall | Output |
|---|---|---|---|
| level 0 only, library default | `uv run python scripts/phase_d/level0_trace.py py --voxel 15901` | 39.6 s | (superseded by the full run) |
| full no-start search, library default (`--full`) | `... py --voxel 15901 --full` | 169.5 s (overlapped with a C++ build) | `py_v15901.json`, `logs/py_full.log` |
| level 0, C++ overlap rule | `... py --voxel 15901 --cpp-rule` | 23.4 s | `py_cpprule_v15901_L0.json` |
| full, C++ rule + centroid screen | `... py --voxel 15901 --full --cpp-rule --centroid` | 120.5 s | `py_cpprule_centroid_v15901.json` |
| level 0, C++ rule + centroid, 5 other voxels | same, `--voxel V` for 17034 18558 19889 22867 5242 | 22.4 to 26.2 s each | `py_cpprule_centroid_vV_L0.json` (counts only; per-FZ records stripped) |
| S3 on C++'s inputs | `... s3 --voxel 15901` | short (not timed) | `s3_on_cpp_inputs_v15901.json` |

Seed 1, `rng = default_rng(1)`; level 0 has no random draws. The three `compare` calls are
`... compare --voxel 15901 --py-tag "" | _cpprule | _cpprule_centroid`; `level0_summary.py` writes
`summary.json` and `tables.md`. The probe builds a copy of `icenine/_rasterize.c` whose
`count_qualified_peaks_c` is a literal port of C++ `CountQualifiedPeaks` (compiled with `cc` into a temporary
directory at run time; the library C file and `cost_functions.py` are unchanged).

## 4. What was found (see `tables.md`, `summary.json`)

- S4 settings: identical at all four levels (N_fz 4,886 at level 0, N_local 9, number of reciprocal vectors
  26/50/64/112, diameter 5.000/3.333/2.222/1.481 deg, 24 symmetry operators). The C++ nQMax is 5..8; the Python
  hook does not expose nQMax (the cost function does not keep `max_q`), equal reciprocal-vector counts are the
  evidence that the cut is the same.
- S1 screen differs. Pairs (FZ point, local-grid point) passing the screen: C++ 12,326, Python 35,267; the
  pass set is identical on only 7 of 4,886 FZ points; peak_overlap is equal on 27.3 % of the 43,974 pairs.
  The difference is the same for every local-grid index j (no dependence on j).
- Cause located by reading and by a per-peak C++ dump: C++ `CountQualifiedPeaks` (`Src/OverlapInfo.tmpl.cpp`
  lines 107-131, loop 113-131) invalidates the overlap of a peak when a lit detector shows no overlap; the Python
  `count_qualified_peaks` (`icenine/cost_functions.py` lines 191-213) and the C port used on the hot path
  (`icenine/_rasterize.c` lines 334-355) do not. In the 3 dumped FZ points (76 peaks), 6 peaks have overlap on
  detector 0 and a lit, non-overlapping detector 1: C++ counts 0 overlap for all 6 (Python's rule counts them).
- Probe, same trace with the C++ rule: screened pairs 12,485 (C++ 12,326), peak_overlap equal on 94.0 %, level-0
  candidates 8,206 (C++ 8,151). Adding the centroid for the radius-3 screen (C++ `VoxelToCenter`; the library
  projects vertex 0): screened pairs 12,338, pass set identical on 4,466 of 4,886 FZ points, peak_overlap equal
  on 98.6 %, candidates 8,110. Over the six interior voxels the probe gives 0.995 to 1.003 times the C++
  level-0 count (library default 1.24 to 1.53).
- S3: the library `_spacing_filter` applied to the C++ trace's own screened lists gives the C++ kept lists on
  3,537 of 3,537 FZ points with at least 2 candidates (same set and same order). S3 is not a source of the gap.
- S2: on pairs screened by both, the re-cost differs by a median of -0.030 (Python lower) with the library
  default and by 0 with the C++ rule (85 % within 1e-6, 99.5 % within 1e-3). This is the same rule.

## 5. Problems, redone work, slow steps

- First patch script failed three times on whitespace in multi-line anchors; rewritten with whitespace-tolerant
  regex. The C++ build itself was fast (10 s).
- The first heredoc attempt to create the patch file was blocked by the Bash hook (the command text contained a
  bare `python3`); the file had to be re-created with the Write tool.
- `timeout` does not exist on macOS (one wasted run); the C++ run takes about 28 s so no timeout is needed.
- The first per-peak dump was gated wrongly: the shared `LocalAngularInterpProcess::operator()` is also used by
  the basic reconstructor that runs before the adaptive one ("Normal Reconstruction", pixel radius 1), so the
  first three FZ points dumped were its, not level 0's. It showed `radius 1` in the counter print, which exposed
  the mistake; the dump was re-gated to the first `GetSpacedCandidates` call and redone (radius 3 confirmed). The
  conclusion from the L0O/L0P data, which were gated correctly, did not change.
- First `compare` showed `n_q_max` Python as null (the cost function does not store `max_q`); the hook no longer
  reports it.
- The Python full run overlapped a C++ rebuild (169 s against 120 s for the later full run): contended.

## 6. Tests

`uv run pytest tests/test_discrete_recorder.py`: 2 passed. Full suite `uv run pytest tests/ --junitxml=<scratch>/junit.xml`
(no extra `-q`, 17:41 to 17:46, 304 s, contended by nothing else of mine): 846 passed, 34 skipped, 1 deselected.
The existing golden tests (`test_findoptimal_refactor`, `test_cma_optimizer`, `test_reconstruction_validation`,
`test_mc_cpp_parity`, `test_variance_stage_parity`) are in that run and pass unchanged, i.e. the default path is unchanged.

## 7. Could be done better

- A permanent, gated trace switch in the C++ source (env var) would remove the patch-and-rebuild routine, but
  changes `Src/`; as a scratch patch it needs `rep` anchors that survive whitespace, so keep patches as unified
  diffs (`cpp_debug_full_scratch.diff.txt`) applied with `patch -p0`, not as string replacements.
- Gate every debug print on the call site, not on shared code: the basic and the adaptive reconstructor share the
  search-point processor. Print one self-describing header per call so a log shows which path produced it.
- Keep the one-voxel C++ configs and mics in the repo (a `configs_cpp/` next to `scripts/phase_d/configs/`) with
  a relative output path; today they live in a session scratch directory and had to be recovered by hand.
- A generic "per-item recorder" convention (event name + dict) for search functions would let parity traces
  reuse one comparison script; `level0_trace.py compare` is already generic over (j, cost, kept) records.
- The probe (`--cpp-rule --centroid`) is a runtime monkeypatch; once the fix is decided, its checks should become a
  test that pins the per-FZ-point counts (Step 2).
- Per-FZ JSON is large (0.6 to 1.2 MB per voxel); store arrays as `npz` or drop the all-j peak_overlap arrays
  for non-15901 voxels (done here by stripping).
- Scripted running of the five other voxels in one shell loop worked (about 2.3 min total); the same loop could
  also launch the C++ side if the configs were kept.
