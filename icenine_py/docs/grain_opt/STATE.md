# AR6 (`grainopt`) state

Rewritten at every step close; at most 120 lines. How: `RUNBOOK.md`. Ideas: `IDEAS.md`. Plan:
`PLAN.md`. Results: `findings.md`, `benchmarks/grain_opt/`. History: `STATE_HISTORY.md`.

## Now

- **GO7g closed** (2026-11-19 22:55 PST, G2b, I68 yes; owner card "Rocking-aware threshold"):
  a strict-or-Poisson acceptance gate (`run_pipeline.py --gate either --store grad_rock10t0
  --model_rocking 1.0`) on the alpha 1.0 rocking images accepts grain 197 and leaves 5 of 853 C
  voxels above 0.1 deg (clean 3, strict gate 20), median 0.0134; on clean images 3 of 853, 0.0144.
  Lower thresholds (x 0.5, x 0.25, any) restore the ground-truth strict cost and the unlit-peak share
  (0.3%) but grain 197 stays rejected: the strict gate on it is a knife edge (A-stage field 0.0927
  against the painted start 0.0877). The 3-frame sum widens clean spots (2.67 x) and fails.
- **Pipeline default changed** (2026-11-19 22:56 PST, scripts only, at the coordinator's request
  after GO7g): `run_pipeline.py --gate either` (A / B stages also accepted when the group's Poisson
  deviance does not rise) is now the default; `--gate strict` turns it off and reproduces earlier
  runs. Checked on clean G20: 3 of 853 C voxels above 0.1 deg, median 0.0144 (same as CP2). Not in
  the library (would need the owner's approval in their own words).
- **GO7f closed** (D-20261119-I1 library part merged df583c97): a, b, c no; d: rocking damage is
  threshold dilution. **GO9c / GO9d** (I20, I19) no. **GO1b** (I4): no clear bias, offsets within
  scan noise. **GO1c** (I3): soft "pixel touched" surrogate (tau 0.25 px) qualifies, 0.0045 /
  0.0017 deg; **GO1d no** (from displaced starts it ends 0.005-0.02 deg away, 3 of 6 within 0.01 at
  0.03 deg; the strict cost ends at 0.002 deg and is faster). I3 closed.
- **Open choices for the owner**: I66 (retrain GNLayerNet on grain-like windows), I35 (G1b from a
  real painted start on rand / corr3), the GO7f defaults (renderer attributes vs config keys, PSF lit
  threshold 1e-3, |sin eta| floor 0.01, high500 freeze test policy), and whether the strict-or-Poisson
  gate should move into the library (needs approval; it is the scripts default now).
- **Next action:** the proposed-idea queue is down to I8 (blur annealing for the Sinkhorn basin
  finder), low expected value since the GN fast fit and the G1b pipeline replaced Sinkhorn; the
  next step needs the owner's pick among the open choices above (or I8). Wake cron as fallback.

## Running jobs

- none.

## Data (git-ignored, in this worktree; link, do not copy)

- `strain` sample images are symlinks into the GO5 agent worktree (branch feature/grain-opt-strain)
  `.claude/worktrees/agent-a58fece2962bad42f` (stack, `strain/` folder, sample `truth.mic` and
  `grid.mic`): do not remove that worktree.
- `rand`, `corr3` stacks and stores: `.../stacks/<tag>_clean.npy`, `.../<tag>/intensity_clean.npz`.
- Rocking stacks (GO7): `.../stacks/grad_rock03_clean.npy`, `grad_rock10_clean.npy` (and stores in
  `grad_rock03/`, `grad_rock10/`); rebuild with `rock_stack.py --alpha A --out <name>` (about 3 min).
- PSF stacks (GO7b): `grad_psf05`, `grad_psf10`, `strain_psf05`, `strain_psf10` (stacks and stores);
  rebuild with `psf_stack.py --tag T --sigma S --out <name>` (about 2 min).
- Ray caches (GO4d): `.../corr3/observed_rays.npz`, `.../rand/observed_rays.npz` (paired observed
  rays; rebuild with `ray_cache.py --sample <s>`, 19 / 48 min).
- `strainrn` (GO5m): noisy `strain` images, `.../stacks/strainrn_clean.npy` and
  `.../strainrn/intensity_clean.npz`, built in this worktree with `build_store.py --tag strain
  --variant realistic_nodrop --out_tag strainrn` (reads the GO5 worktree's render; 3 s).

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

- `scripts/grain_opt/fit_strain.py`: strain fit (GO5b arms E@GT, FE@GT, FE@P; GO5c J@P, S@P, SR@P, S@GT; GO5d E@F, O@P, Ost@P; GO5e-g `--cost sd|int|cen|ch|hyb`, `--lam`, `--seed`; GO5g Sx@P; GO5i GN@P, GNH@P with `gn_fit`; GO5j GA1@P, GA2@P; GO5k GF2@P, GFF@P; GO5l GP4@P, GP6@P).
- `icenine/grain_fit.py` (SpotModel, regions, Sinkhorn / LS objectives, FieldA, GrainProblem);
  tests `tests/test_grain_fit.py` (16 pass).
- `scripts/grain_opt/`: `common.py` (Sample loader), `run_grains.py` (stage A/B on G20),
  `summarize.py` (tables), `probe_*.py`, `check_render.py`, `build_store.py`, `build_page.py`.

## Tracking page and wake job

- Page: https://claude.ai/artifact/B2wANBjWDeTN86JNgmqBKb (publish `benchmarks/grain_opt/page/index.html`).
- Wake job: session cron `dc9287aa`, every 3 h at :17, created 2026-11-16 01:50 PST (after a
  machine restart killed `8820401c`), expires 2026-11-23; recreate on 2026-11-22.

## Owner items

- Standing rule (owner via the coordinator, 2026-11-16): no idling between steps; after a
  step closes, plan and start the next one at once (RUNBOOK "On every wake"). The wake job is
  a fallback.
- Style (owner, 2026-11-16): write "ground truth", not "truth" (RUNBOOK "Wording"). Existing
  docs and the page get the wording at the GO1/GO2 close.

- None open. Merges to develop: owner only.
