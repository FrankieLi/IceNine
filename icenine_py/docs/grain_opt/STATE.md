# AR6 (`grainopt`) state

Rewritten at every step close; at most 120 lines. How: `RUNBOOK.md`. Ideas: `IDEAS.md`. Plan:
`PLAN.md`. Results: `findings.md`, `benchmarks/grain_opt/`. History: `STATE_HISTORY.md`.

## Now

- **GO7f library part merged** (df583c97, 2026-11-19 21:15 PST, D-20261119-I1): opt-in
  `SpotModel(psf_sigma_px, rocking_deg)`, renderer PSF / rocking, strict-cost `rocking_deg`
  window; defaults bit-identical; three code reviews. Script flags `--model_psf` (fit_strain) and
  `--rocking_deg` (run_pipeline) in b20afe48.
- **GO7f-a no** (model PSF 1 px, ratio 2.5, `strain_psf10`): strain medians 2.16e-4 / 4.15e-4 at
  5e-4 / 1e-3 (GO7e PSF arm 2.12e-4 / 3.49e-4), orientation medians 0.042 / 0.086 deg.
- **GO7f-b no** (strict-cost rocking window 1.0 on `grad_rock10`): grain 197 still rejected (and
  493 now too); C voxels above 0.1 deg 21 of 853 (no window 20, clean 3); median 0.0138. The window
  barely lowers the strict cost; torch path 3.2 s per voxel.
- **Open choices for the owner** (unchanged): I66 (retrain GNLayerNet on grain-like windows), I35
  (G1b from a real painted start on rand / corr3, paint likely over 2 h). Owner decisions on GO7f
  defaults pending: renderer options as attributes vs config keys, PSF lit threshold 1e-3 of the
  frame maximum, |sin eta| floor 0.01, high500 freeze test policy.
- **Next action:** close GO7f-c (model PSF without the ratio filter) and GO7f-d (rocking damage
  diagnosis) by the PLAN rules; then pick the next G-goal step per the program balance.

## Running jobs

- GO7f-c: `fit_strain.py --arm GP4@P --cost hyb --lam 8 --per_level 12 --stack
  strain_psf10_clean.npy --model_psf 1 --out go7f` (started 2026-11-19 22:15 PST, log
  `benchmarks/grain_opt/go7f/go7f_c.log`).
- GO7f-d: `scripts/grain_opt/go7f_diag.py` (started 22:19 PST, log `.../go7f/go7f_d.log`).

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
