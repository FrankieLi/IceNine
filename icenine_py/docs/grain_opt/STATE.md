# AR6 (`grainopt`) state

Rewritten at every step close; at most 120 lines. How: `RUNBOOK.md`. Ideas: `IDEAS.md`. Plan:
`PLAN.md`. Results: `findings.md`, `benchmarks/grain_opt/`. History: `STATE_HISTORY.md`.

## Now

- **GO1 (objective bias) nearly done.** Rows in `benchmarks/grain_opt/probe/bias_*.jsonl`.
  So far every Sinkhorn variant has a flat bottom: CMA from the truth ends 0.031 to 0.041 deg
  away (balanced 0.013 to 0.095) with the objective only 0.2 to 0.7% lower; the strict cost
  (grain mean) ends 0.006 to 0.008 deg away (passes the 0.01 rule). Left: intensity and LS arms,
  dilation arms on grain 302.
- **GO2 (optimizer Pareto on G20) started early on the strict-cost arms** (PLAN GO2 arms 1 to 4)
  via `scripts/grain_opt/go2_strict_arms.sh`; Sinkhorn / LS arms (5, 6) start after GO1 ends.
- **GO3 done** (samples `rand`, `corr3`, commit 52cab1b7).
- **Next action:** when GO1's last rows land, write the GO1 findings section, IDEAS statuses,
  progress.json step row; start GO2 arms 5 and 6; when all GO2 arms end, run
  `summarize.py go2_*`, add Pareto rows, close GO2 (runbook "Close").

## Running jobs

- `probe_bias.py` (GO1): grains 381 and 302 (tags bias_381, bias_302), started 20:56 PST.
- `go2_strict_arms.sh` (GO2 arms 1-4, 10 runs in sequence), started about 21:15 PST; log
  `benchmarks/grain_opt/runs_go2_strict.log`, per-run logs `benchmarks/grain_opt/runs/*.log`.

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
