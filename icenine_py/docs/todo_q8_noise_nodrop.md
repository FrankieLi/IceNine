---
title: "TODO: realistic images without the random spot drop (ROADMAP Q8)"
subtitle: "Plan: a no-drop noise variant for the Phase D full-sample images and the E0 per-voxel images"
date: "2026-10-13"
geometry: margin=1in
fontsize: 11pt
---

# Status

**Approved 2026-10-13; Step 0 done; Steps 1, 2, 4 and 5 done 2026-10-14, Step 3 not needed** (branch `feature/python-speed-parity-q8-nodrop`, not merged; results and the owner's open decisions are in MIGRATION_HISTORY "Q8"; the Phase D pilot row of Step 4 was not run) (owner: remove the random spot drop, keep jitter, hot pixels and
blobs; Ilium Agent does the render; generate the Phase D data on Ilium Agent's machine; Step 4 on the
200-voxel set only). Runs after S3 (same checkout, and its renders and re-runs would disturb S3's
timings). Branch: task branch
`feature/python-speed-parity-q8-nodrop` off `feature/python-speed-parity` (the ROADMAP lives only on
the agents' feature branches, so the Q8 row can be updated in the same branch). An `implementer`
agent builds it, `/review-task` when done, fixes to a fresh implementer. Not during S3's timing steps
(S3 Steps 0 and 4 need a quiet machine).

# Why

Both "realistic" noise models erase each spot at random with probability 0.1, regardless of
intensity. A strong spot does not vanish, and 10% is far too high. It penalises every rule that asks
detector 1 to confirm a peak (C++ and ray rules reject about 9% of real peaks at the truth on the
realistic images), so every realistic-image comparison so far is tilted against the strict rules
(findings_finisher_mc_study §7 item 8). Of the options listed there (no drop; intensity-based drop
near a threshold; about 1% fixed), the owner chose **no drop**.

# Two noise models in scope

| images | where the noise is | miss draw | used by |
|---|---|---|---|
| Phase D full-sample (`ScatteringData_PhaseD/<tag>/realistic/`) | `scripts/phase_d/noise.py` (`NoiseParams`, per connected component of a frame) | always drawn (`noise.py:74`), so `p_miss 0` leaves every other draw unchanged | Phase D pilot and reseed, `sanity_checks.py`, the `*_realistic*.config` arms |
| E0 per-voxel (`scripts/findoptimal_robustness`, variant "all") | library `icenine/orientation_eval.py` `RealismConfig` (per spot window) | drawn only when `p_miss > 0`, as the last draw of `make_realistic_windows` (`orientation_eval.py:682`) | E0 200-voxel set: level-0 study (Q3), ray-overlap study, FindOptimal robustness |

# Steps

**Step 0 (done 2026-10-13): Phase D data on Ilium Agent's machine.** `render_full.py --workers 6`
(55 s; 5.1 GB summed over the workers, Newt Agent's 10-worker run summed 12.4 GB). Every
deterministic field of the run record equals Newt Agent's `scripts/phase_d/results/render_full.json`
(57,389 spots, 5,743 missed, 7,414,841 lit clean pixels; only timing and memory differ), so the data
is identical. Output (gitignored): `Examples/Example2.ManyGrains/ScatteringData_PhaseD/full/{clean,
realistic}/`. The sample `.mic` files were already in the repo. Not rendered: `full_q16` (Q-max 16;
Newt Agent's 10-worker run summed 19.7 GB, so use 4 workers here) and `pilot`; render them only if a
Step 4 re-run needs them. The render overwrites the tracked `results/render_full.json` and
`results/preflight_render_full.json`; restore them with `git checkout` (done).

**Step 1: no-drop variants (opt-in; defaults and existing images unchanged).**

- Phase D: `NoiseParams` gets named variants `toy` (today's values) and `nodrop` (`p_miss 0`).
  `render_full.py --noise nodrop --variants realistic` writes `<tag>/realistic_nodrop/` next to
  `realistic/` (kept for paired comparisons), and records the variant in its run record under a new
  name (do not overwrite the tracked `render_full.json`). `make_configs.py` adds
  `*_realistic_nodrop.config` for the existing arms (`mc`, `cma`, `cma_noretry`).
- E0: `RealismConfig.named` gets `noise_nodrop` and `all_nodrop` (`p_miss 0`); the E0 harness
  (`scripts/findoptimal_robustness/common.py`) accepts the variant `realistic_nodrop`.
- Tests (`tests/test_phase_d.py`, `tests/test_orientation_eval.py` or the existing realism tests):
  - `nodrop` frame equals the `toy` frame except on the components `toy` missed and their hot pixels
    and blobs;
  - no clean spot is missing in `nodrop`;
  - E0: each window equals variant `all` wherever `all` did not miss. This holds only if each call to
    `make_realistic_windows` has its own generator (E0 seeds one per voxel, `common.py:103`), since
    the miss draw comes last. If the generator is shared, draw and ignore the miss instead, and say so.

**Step 2: render and sanity checks.** Render `full/realistic_nodrop` (about 1 min, 6 workers). Run
`sanity_checks.py --variant realistic_nodrop` (writes `results/sanity_full_realistic_nodrop.json`)
and put three variants side by side (clean, realistic, realistic_nodrop): cost at the truth, and the
fraction of real peaks the C++ and ray rules reject at the truth (about 9% on `realistic` today).
Expected: `realistic_nodrop` close to clean on rejections, and on cost apart from the jitter.

**Step 3 (only if Step 4 needs it): Q-max 16.** `render_full.py --max-q 16 --tag full_q16 --workers
4 --noise nodrop --variants realistic`.

**Step 4: re-runs on the 200-voxel set (owner 2026-10-13: the two E0 rows; the Phase D pilot row
waits with X1). Costs are estimates.**

| re-run | what it revisits | images | estimated cost |
|---|---|---|---|
| E0 `mc` arm, `CppOverlapRule` off and on, seed 0 | Q3 (realistic half: rule off 25 wrong of 200 on `all`) | E0 `all_nodrop` | about 1 h per arm on 10 workers (scaled from the original E0 budget) |
| ray-overlap study, realistic columns (lenient 18 vs ray 29-32 wrong of 200) | lenient vs ray rule | E0 `all_nodrop` | similar per arm |
| Phase D pilot and reseed, realistic arms | X1 / X2 | `full/realistic_nodrop` (+ `full_q16`) | 2-5 h per run; the owner deferred long runs (X1 later) |

Each re-run quotes wrong counts and evaluations side by side with the clean and the old realistic
runs, paired per voxel, no p-values.

**Step 5: docs.** MIGRATION_HISTORY section; findings_finisher_mc_study §7 item 8 status; caveat lines in
the findings that quote realistic results (pointing to the no-drop figures); ROADMAP Q8 row and the
"How to use this file" rule "clean images only until Q8" (a shared section: the owner decides when to
lift it).

# Rules

No default flips (the defaults stay `toy` and `all`; switching any default to no-drop is listed for
the owner); no golden regeneration; existing images and run records untouched; clean-image results
are not affected.
