# AR6 (`grainopt`) state

Rewritten at every step close; at most 120 lines. How: `RUNBOOK.md`. Ideas: `IDEAS.md`. Plan:
`PLAN.md`. Results: `findings.md`, `benchmarks/grain_opt/`. History: `STATE_HISTORY.md`.

## Now

- **GO4b closed** (2026-11-15 22:55 PST): per-voxel step ladder L/183 recovers corr3 up to 2 deg
  (9/11/47% above 0.1 deg at 0.5/1/2); rand is a search problem. `findings.md` "GO4b".
- **GO4c closed** (23:25 PST): L/366/3 recovers corr3 to 3/9/31% above 0.1 deg at 0.5/1/2 deg;
  rand about half; 5 deg needs a global search (I36 blocked on ray-index cost).
- **GO5 (strain) delegated** to an implementer agent, branch `feature/grain-opt-strain`. On
  hand-back: verify tests, review the diff, merge `--no-ff` into `feature/grain-opt`, place its
  findings text, close GO5.
- **Next action:** wait for the GO5 hand-back (merge after tests and review). Then plan GO6 (objective variants: I18 Poisson, I19 moments, I3 soft pixel-touched model)
  or GO4d (cached ray index for I36), whichever the results point to.

## Running jobs

- Implementer agent for GO5 (branch `feature/grain-opt-strain`). If the session died: check that
  branch's commits and resume from them; never start a second agent on the branch.

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
