# AR6 (`grainopt`) state

Rewritten at every step close; at most 120 lines. How: `RUNBOOK.md`. Ideas: `IDEAS.md`. Plan:
`PLAN.md`. Results: `findings.md`, `benchmarks/grain_opt/`. History: `STATE_HISTORY.md`.

## Now

- **Step GO1 (objective bias, ideas I1 I2 I17 I21) is running.** First run (20:24) was cut by its
  40-minute timeout and its rows were lost in a pipe. Rerun 20:56 PST 2026-11-15, one process per
  grain: `probe_bias.py --grains 381 --evals 120 --tag bias_381` (and 302); rows are appended to
  `benchmarks/grain_opt/probe/bias_<grain>.jsonl` as they finish (7 arms each).
- **Next action:** read the GO1 rows; pick the objective with the smallest offset (rule in PLAN
  GO1: below 0.01 deg on both grains). If one qualifies: plan GO2 (optimizer Pareto on G20). If
  none: plan a cause study (I3 soft pixel-touched model, I4 noise floor) first.
- Open issue: the objective's minimum sat about 0.05 deg from the truth in the smoke test; the
  strict-cost gate rejected both updates.

## Running jobs

- `probe_bias.py` x2 (GO1, grains 381 and 302), started 20:56 PST 2026-11-15, no timeout; logs
  `benchmarks/grain_opt/probe/bias_<grain>.log`.

## Data (git-ignored, in this worktree; link, do not copy)

- `Examples/Example2.ManyGrains/ScatteringData_PhaseD/stacks/grad_clean.npy`: `grad` binary
  stack (rebuilt 2026-11-15 from a fresh render; equals the stored record).
- `.../ScatteringData_PhaseD/grad/intensity_clean.npz`: clean intensities (sparse CSR by frame).
- `.../ScatteringData_PhaseD/grad/clean/`: the ASCII render.
- `big`, `held`, `mid` stacks: not present on this machine's worktrees (STATE.md of Newt lists
  removed worktrees). Re-render with `render_full.py --mic <truth>` when needed (about 1 min).
- Committed truth and starts: `benchmarks/autoresearch/ar3/gradient/sample/` (truth .mic,
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

- None open. Merges to develop: owner only.
