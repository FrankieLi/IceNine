# AR6 (`grainopt`) state

Rewritten at every step close; at most 120 lines. How: `RUNBOOK.md`. Ideas: `IDEAS.md`. Plan:
`PLAN.md`. Results: `findings.md`, `benchmarks/grain_opt/`. History: `STATE_HISTORY.md`.

## Now

- **Step GO1 (objective bias, ideas I1 I2 I17 I21) is running.** `probe_bias.py --grains 381 302`
  (7 arms, CMA over the grain rotation from the truth). Started 2026-11-15 20:24 PST under a
  40-minute `timeout`; rows go to the console only (rerun the missing arms with
  `--grains ... --variants ...` if it was cut off), and `benchmarks/grain_opt/probe/bias.json` at
  the end.
- **Next action:** read the GO1 rows; pick the objective with the smallest offset (rule in PLAN
  GO1: below 0.01 deg on both grains). If one qualifies: plan GO2 (optimizer Pareto on G20). If
  none: plan a cause study (I3 soft pixel-touched model, I4 noise floor) first.
- Open issue: the objective's minimum sat about 0.05 deg from the truth in the smoke test; the
  strict-cost gate rejected both updates.

## Running jobs

- `scripts/grain_opt/probe_bias.py` (GO1), one process, started 20:24 PST 2026-11-15.

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
