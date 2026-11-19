# AR6 (`grainopt`) state

Rewritten at every step close; at most 120 lines. How: `RUNBOOK.md`. Ideas: `IDEAS.md`. Plan:
`PLAN.md`. Results: `findings.md`, `benchmarks/grain_opt/`. History: `STATE_HISTORY.md`.

## Now

- **GO5s closed** (2026-11-18 23:50 PST, I61 yes): sample `straingq` (mean + linear + quadratic +
  smooth random field; oracle linear share 0.33). Linear fits keep the mean but explain almost none
  of the intragranular field (rel. error 0.98-1.10; linear floor 0.82).
- GO5r: Gauss-Newton + spot spreads, I60 no (best mean at 15 s per grain). I65 logged.
- **Approval** D-20261118-I1: strain queue I60-I64 (literature review, sibling thread "Review neural
  methods for grain gradients"), then I25 at low priority.
- **Checkpoint 2026-11-19 00:16 PST** (owner shutting the machine down): GO5t (I62) was stopped at
  20 / 36 (E@GT) and 21 / 36 (GP4@P) grains, before the 1e-3 level; `fit_strain.py` writes its
  result file only at the end, so nothing is judged yet. The partial per-grain lines are kept in
  `benchmarks/grain_opt/go5t/*.log.partial.jsonl` (preview, not judged: relative intragranular
  error at 5e-4 about 1.1-1.2 on 8-9 grains). The plan and rule (PLAN "GO5t") stand.
- **Next action (resume first):** rerun GO5t from scratch with the two commands below (about 65
  minutes, 2 processes), judge by the GO5t rule, close the step; then I63 only if the rule calls
  for it, I64 (GP / Laplacian prior), I25, I65.
  `uv run python scripts/grain_opt/fit_strain.py --arm E@GT --cost hyb --lam 8 --per_level 12
  --sample straingq --ratio_hi 2.5 --grad --quad --evals_e 1500 --out go5t`;
  `uv run python scripts/grain_opt/fit_strain.py --arm GP4@P --cost hyb --lam 8 --per_level 12
  --sample straingq --ratio_hi 2.5 --grad --quad --evals_g 1000 --out go5t`.

## Running jobs

- None (GO5t stopped at the checkpoint; rerun it first).

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
