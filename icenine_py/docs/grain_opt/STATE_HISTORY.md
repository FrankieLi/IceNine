# AR6 (`grainopt`) state history

Finished steps and old notes, newest first. Not read on a wake.

## GO4, GO4b, GO4c, GO5 (2026-11-15/16)

GO4/4b/4c tables in `benchmarks/grain_opt/go4*/`; GO5 built by an implementer agent (branch
feature/grain-opt-strain, merge d7b06158), reviewed by a code-reviewer agent: one latent bug (strain
not rotated with symmetry swaps) fixed, one test gap (frame of M) closed. The gated global search
(I36) was tried and blocked on cost.

## GO2 (2026-11-15): optimizer Pareto on G20

Tables `benchmarks/grain_opt/go2/`; findings "GO2". Sinkhorn arms cut to one 20-grain arm plus
two 10-grain arms to stay inside the 2 h step cap (2-4 min per grain for Sinkhorn stages).

## GO1 (2026-11-15): objective bias

Rows `benchmarks/grain_opt/probe/bias_*.jsonl`; table in `findings.md` "GO1". The first run was
cut by a 40-minute timeout and lost its rows in a pipe; the rerun writes rows as it goes.
GO3 (samples) and I31 (C++ strain survey) were done while GO1 ran.

## GO0 (2026-11-15): work item (1), commit 5a3162d2

Grouping, regions, coverage renderer, objectives. Numbers in `findings.md` "Work item (1)".
The `grad` stack was rebuilt from a fresh render. The coordinator's brief (work items 1 to 4)
became AR6 on 2026-11-16 at the owner's request.

## GO5 close (2026-11-16 00:00 PST): STATE "Now" at the time

- **GO5 closed** (2026-11-16 00:00 PST): opt-in strain in SpotModel and the stock renderer
  (g = O M g0), merged after tests and review (fixes: strain follows symmetry swaps in
  GrainProblem; crystal-frame test). Sample `strain` built. `findings.md` "GO5".
- GO4c closed: corr3 recovered to 2 deg by L/366/3; rand and 5 deg need a global search (I36
  blocked on ray-index cost).
- **Next action:** plan and run GO5b: per-grain strain fit on `strain` (I38): grain field +
  6 strain numbers, CMA on the coverage or strict-like cost, from the ground-truth orientation
  first, then from the painted start; recovered vs ground-truth strain per level. Then GO6
  (objective variants) or the cached ray index for I36.

## GO5b close (2026-11-16 02:20 PST): STATE "Now" at the time

- **GO5b closed** (2026-11-16 02:20 PST): per-grain strain fit on `strain` (S12, coverage cost,
  CMA). Strain only at the ground-truth orientation reaches about 2e-4 at every level (5e-4 and
  1e-3 recovered, 1e-4 below the floor). Field + strain: 3.4-3.8e-4 from the ground truth,
  5.3-6.7e-4 from the painted start; 19 of 24 joint fits stop above the ground truth's cost
  (optimizer-limited). Robust start field added (painted outliers). `findings.md` "GO5b".
- GO4c closed: corr3 recovered to 2 deg by L/366/3; rand and 5 deg need a global search (I36
  blocked on ray-index cost).
- **Next action:** plan and run GO5c (I39): a joint strain fit that reaches the cost optimum
  (staged legs field -> strain only -> joint, CMA restarts or a bigger budget) from the painted
  start, target the 2e-4 floor; then I40 (sharper sub-pixel cost for 1e-4). Runs take minutes
  (`fit_strain.py`, about 30 ms per evaluation). Then GO6 or the cached ray index for I36.

## Now as of GO5c close (2026-11-16 11:45 PST)

- **GO5c closed** (2026-11-16 11:45 PST): reaching the cost optimum does not recover strain.
  `J@P` (2400 joint evaluations) ends below the ground truth's cost in 10 of 12 grains but keeps
  4.5-6.5e-4 strain error: with the orientation field free the coverage cost is
  objective-limited (0.01 deg rotation trades with 2-4e-4 strain). Staged `S@P` gives 3.4e-4 /
  3.0e-4 at 5e-4 / 1e-3 (floor 1.9e-4 / 2.0e-4). `findings.md` "GO5c".
- GO5b closed: strain-only at the ground-truth orientation reaches ~2e-4 (1e-4 below the floor).
- **Next action:** plan and run GO5d (I41): fix the orientation field from a strain-insensitive
  fit (e.g. field leg without strain, or the strict cost), then fit strain alone; compare with
  `S@P` and the `E@GT` floor. Then I40 (sharper sub-pixel cost). Runs take minutes
  (`fit_strain.py`). Then GO6 or the cached ray index for I36.
