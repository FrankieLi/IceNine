# AR6 (`grainopt`) state

Rewritten at every step close; at most 120 lines. How: `RUNBOOK.md`. Ideas: `IDEAS.md`. Plan:
`PLAN.md`. Results: `findings.md`, `benchmarks/grain_opt/`. History: `STATE_HISTORY.md`.

## Now

- **GO5f closed** (2026-11-16 15:20 PST): per-spot centroid match (`--cost cen`) is the most
  consistent strain cost so far (repeat scatter 5e-6; 1.3 / 2.2 / 1.7e-4 at the ground-truth
  orientation) but sharper than coverage only at 1e-3 (I42 no by the rule). From the painted
  start (staged) it gives 1.4 / 2.2 / 3.7e-4 (coverage 7.1 / 3.4 / 3.0e-4) while the orientation
  drifts to 0.06 deg (coverage 0.012). Chamfer (`ch`) is biased. `findings.md` "GO5f".
- GO5e: one-sided sub-pixel costs biased. GO5d: strain-free field fits take up strain as rotation.
- **Next action:** plan and run GO5g (I43): hybrid cost, coverage for the orientation field and
  centroids for strain (weighted sum, or staged: field on coverage, then strain and mean
  orientation on centroids), merged spots down-weighted; painted start, judged against `cen`
  `E@GT` and coverage `S@P`, plus orientation error. Then GO6 or the cached ray index.

## Running jobs

- None.

## Data (git-ignored, in this worktree; link, do not copy)

- `strain` sample images are symlinks into the GO5 agent worktree (branch feature/grain-opt-strain)
  `.claude/worktrees/agent-a58fece2962bad42f` (stack, `strain/` folder, sample `truth.mic` and
  `grid.mic`): do not remove that worktree.
- `rand`, `corr3` stacks and stores: `.../stacks/<tag>_clean.npy`, `.../<tag>/intensity_clean.npz`.

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

- `scripts/grain_opt/fit_strain.py`: strain fit (GO5b arms E@GT, FE@GT, FE@P; GO5c J@P, S@P, SR@P, S@GT; GO5d E@F, O@P, Ost@P; GO5e/f `--cost sd|int|cen|ch`, `--seed`).
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
