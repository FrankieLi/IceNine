# AR6 (`grainopt`) state

Rewritten at every step close; at most 120 lines. How: `RUNBOOK.md`. Ideas: `IDEAS.md`. Plan:
`PLAN.md`. Results: `findings.md`, `benchmarks/grain_opt/`. History: `STATE_HISTORY.md`.

## Now

- **GO2 closed** (2026-11-15 22:35 PST): the frontier is the grain field on the strict cost (100
  objective evaluations per voxel: 3 of 853 voxels > 0.1 deg against 18 for AR3 C); smooth basin
  stages add nothing from the painted start. `findings.md` "GO2", `benchmarks/grain_opt/go2/`.
- **GO4 (robustness on `rand`, `corr3`) running**: `scripts/grain_opt/go4_arms.sh` (10 runs),
  started 22:40 PST; PLAN.md "GO4".
- **Next action:** when GO4 ends: `summarize.py go4_rand_* --out go4/rand` (and corr3), per-level
  tables into a GO4 findings section, Pareto rows per sample (new chart group or note), close GO4.
  Then plan GO5 (strain: SpotModel opt-in strained g, pure-dilation sign test, strained sample
  render; sources `findings_cpp_strain.md`, HEXRD reference in shmem "ref-hexrd").

## Running jobs

- `go4_arms.sh` (GO4, 10 runs in sequence); log `benchmarks/grain_opt/runs_go4.log`.

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
