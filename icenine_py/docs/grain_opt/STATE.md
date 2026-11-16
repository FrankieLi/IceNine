# AR6 (`grainopt`) state

Rewritten at every step close; at most 120 lines. How: `RUNBOOK.md`. Ideas: `IDEAS.md`. Plan:
`PLAN.md`. Results: `findings.md`, `benchmarks/grain_opt/`. History: `STATE_HISTORY.md`.

## Now

- **GO4 closed** (2026-11-15 22:45 PST): from a grain-centre start, local fits recover only voxels
  within about 0.5 deg of their start. `findings.md` "GO4".
- **GO4b (far starts: neighbour continuation N, step ladder L) running**: `go4b_arms.sh` (10 runs),
  started 22:36 PST. Smoke: corr3 grain 20 (1 deg) median 0.923 to 0.163 deg at 103 eval/voxel.
- **GO5 (strain) delegated** to an implementer agent in its own worktree, branch
  `feature/grain-opt-strain` (PLAN "GO5"). On hand-back: verify tests, review the diff (code under
  `icenine/`), merge `--no-ff` into `feature/grain-opt`, place its findings text.
- **Next action:** when GO4b ends: summarize (`--out go4b/<sample>`), per-level table against
  GO4's best arm, close GO4b. When GO5 hands back: check and merge as above.

## Running jobs

- `go4b_arms.sh` (GO4b, 10 runs); log `benchmarks/grain_opt/runs_go4b.log`.
- Implementer agent for GO5 (background; worktree under `.claude/worktrees/agent-*`, branch
  `feature/grain-opt-strain`). If the session died: check that branch's commits and resume the
  step from them; never start a second agent on the branch.

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
