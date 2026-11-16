# AR6 (`grainopt`) state

Rewritten at every step close; at most 120 lines. How: `RUNBOOK.md`. Ideas: `IDEAS.md`. Plan:
`PLAN.md`. Results: `findings.md`, `benchmarks/grain_opt/`. History: `STATE_HISTORY.md`.

## Now

- **GO1 closed** (2026-11-15 22:00 PST): only the strict cost (grain mean) puts its optimum within
  0.01 deg of the ground truth; Sinkhorn is flat within 0.03-0.04 deg. `findings.md` "GO1".
- **GO2 (optimizer Pareto on G20) running.** Strict-cost arms done (`benchmarks/grain_opt/go2/
  tables.md`): field CMA on the strict cost at 100 evaluations leaves 3 of 853 voxels above 0.1
  deg (the 3 wrong paint voxels) against 18 for AR3 C. Running: `go2_Ask150_As100_V20` (20
  grains, 2-4 min per grain), then `go2_smooth_arms2.sh` (intensity-Sinkhorn and LS basin arms
  plus a strict reference on the 10-grain subset).
- **GO3 done** (samples `rand`, `corr3`). **I31 done** (C++ strain survey).
- **Next action:** when the smooth arms end: `summarize.py` on all go2 runs (`--out go2/tables`),
  write the GO2 findings section (Pareto frontier, wall time per voxel, the Sinkhorn arms' gate
  outcome), add Pareto rows, close GO2. Then plan GO4: robustness of the GO2 frontier pipelines
  on `rand` and `corr3` (`--sample rand --start centre --group true`) and on `grad`.

## Running jobs

- `run_pipeline.py --name go2_Ask150_As100_V20` (started about 21:20 PST 2026-11-15).
- `go2_smooth_arms2.sh` (waits for the run above, then 3 runs on 10 grains); log
  `benchmarks/grain_opt/runs_go2_smooth.log`.

## Data (git-ignored, in this worktree; link, do not copy)

- `Examples/Example2.ManyGrains/ScatteringData_PhaseD/stacks/grad_clean.npy`: `grad` binary
  stack (rebuilt 2026-11-15 from a fresh render; equals the stored record).
- `.../ScatteringData_PhaseD/grad/intensity_clean.npz`: clean intensities (sparse CSR by frame).
- `.../ScatteringData_PhaseD/grad/clean/`: the ASCII render.
- `big`, `held`, `mid` stacks: not present on this machine's worktrees (STATE.md of Newt lists
  removed worktrees). Re-render with `render_full.py --mic <truth>` when needed (about 1 min).
- Committed ground truth and starts: `benchmarks/autoresearch/ar3/gradient/sample/` (ground truth .mic,
  sample.json), `.../gradient/runs/P/maps.npz` (painted start "P"),
  `.../ar3/fine/runs/C/refine.json.gz` (AR3 step 4 per-voxel CMA, reference).

## Code

- `icenine/grain_fit.py` (SpotModel, regions, Sinkhorn / LS objectives, FieldA, GrainProblem);
  tests `tests/test_grain_fit.py` (16 pass).
- `scripts/grain_opt/`: `common.py` (Sample loader), `run_grains.py` (stage A/B on G20),
  `summarize.py` (tables), `probe_*.py`, `check_render.py`, `build_store.py`, `build_page.py`.

## Tracking page and wake job

- Page: https://claude.ai/artifact/B2wANBjWDeTN86JNgmqBKb (publish `benchmarks/grain_opt/page/index.html`).
- Wake job: session cron `8820401c`, every 3 h at :17, created 2026-11-15 20:55 PST, expires
  2026-11-22; recreate on 2026-11-21.

## Owner items

- Style (owner, 2026-11-16): write "ground truth", not "truth" (RUNBOOK "Wording"). Existing
  docs and the page get the wording at the GO1/GO2 close.

- None open. Merges to develop: owner only.
