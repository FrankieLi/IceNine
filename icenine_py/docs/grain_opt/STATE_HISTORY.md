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

## Now as of GO5d close (2026-11-16 13:30 PST)

- **GO5d closed** (2026-11-16 13:30 PST): fixing the field first, then strain alone, does not
  help (I41 no). Field fits without a strain model take up the strain as rotation (0.009 / 0.022
  / 0.040 deg at 1e-4 / 5e-4 / 1e-3); strain after them 1.5 / 5.3 / 6.0e-4, worse than the staged
  joint `S@P` (7.1 / 3.4 / 3.0e-4) above 1e-4. At a fixed exact orientation the strain fit
  scatters 1-3e-4 run to run: the coverage cost is flat there. `findings.md` "GO5d".
- GO5c: reaching the optimum does not recover strain (orientation-strain trade-off).
- **Next action:** plan and run I40 (GO5e): a sharper sub-pixel cost (intensity-weighted
  coverage or distance to lit pixels) so 1e-4 strain moves the cost. Judge first at the
  ground-truth orientation (`E@GT`, floor band 1-3e-4), then in the staged joint fit from the
  painted start (`S@P`). Runs take minutes (`fit_strain.py`). Then GO6 or the cached ray index
  for I36.

## Now as of GO5e close (2026-11-16 15:15 PST)

- **GO5e closed** (2026-11-16 15:15 PST): one-sided sub-pixel costs (signed distance to the lit
  set, clean intensity) are sharp (repeat scatter 6e-5) but biased: at the ground-truth
  orientation the strain error is 0.9-2.5e-3 (coverage 1.6-2.7e-4); they reward pulling the
  predicted footprint into the spot. I40 no; painted-start fit not run. `findings.md` "GO5e".
- GO5d: a field fit without strain takes up strain as rotation (I41 no). GO5c: the joint
  coverage fit is objective-limited.
- **Next action:** plan and run GO5f (I42): a two-sided cost, per-spot centroid match (predicted
  mass centroid vs observed intensity centroid, HEXRD-style) and/or symmetric signed distance;
  judge first at the ground-truth orientation (`E@GT`, two seeds) against the coverage band,
  then the staged painted-start fit. `fit_strain.py --cost`. Then GO6 or the cached ray index.

## Now as of GO5f close (2026-11-16 15:20 PST)

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

## Now as of GO5g close (2026-11-16 15:35 PST)

- **GO5g closed** (2026-11-16 15:35 PST): the hybrid cost `cen` + 8 x coverage (`--cost hyb
  --lam 8`), staged from the painted start, recovers strain at 1.9 / 1.5 / 2.2e-4 (1e-4 / 5e-4 /
  1e-3) with orientation 0.008-0.012 deg; coverage alone 7.1 / 3.4 / 3.0e-4. I43 yes.
  `findings.md` "GO5g"; costs table `COSTS.md`.
- GO5f: centroids alone good strain, orientation drifts. GO5c-e: coverage objective-limited,
  one-sided sub-pixel costs biased.
- **Next action:** plan and run GO5h (I44): confirm the hybrid on 12 grains per level (S36),
  lambda 8 and 16, painted start, judged against coverage `S@P` on the same grains; then
  Gauss-Newton on centroid residuals (I7 / I19). Then GO6 or the cached ray index.

## Now as of GO5h close (2026-11-16 16:10 PST)

- **GO5h closed** (2026-11-16 16:10 PST): the hybrid strain fit holds on 36 grains (12 per
  level, painted start, staged): lambda 8 gives 2.3 / 1.5 / 2.1e-4 against coverage 4.6 / 2.8 /
  6.3e-4 on the same grains, orientation 0.009-0.017 deg; lambda 16 equal within scatter.
  Remaining: a tail of 1e-3 grains at 4-9e-4; about 50 s per grain with CMA. `findings.md` "GO5h".
- GO5g: hybrid introduced (I43 yes). GO5c-f: coverage objective-limited, one-sided sub-pixel
  costs biased, centroids alone let the orientation drift.
- **Next action:** plan and run GO5i (I44 second part): Gauss-Newton (or LM) on per-spot centroid
  residuals for strain and mean rotation, the field gradient fixed from the coverage field leg;
  judge on S36 against `hyb` lambda 8 (strain, orientation, wall time). Then GO6 or the cached
  ray index.

## Now as of GO5i close (2026-11-16 16:30 PST)

- **GO5i closed** (2026-11-16 16:30 PST): Gauss-Newton on centroid residuals (`--arm GN@P`,
  Tukey weights 3 / 1.5 / 1 px) gives strain 1.2 / 1.3 / 2.9e-4 on S36 at 10 s per grain (CMA
  hybrid 2.3 / 1.5 / 2.1e-4 at 53 s), but orientation 0.021 / 0.039 deg at 5e-4 / 1e-3: the
  field gradient comes from the zero-strain coverage leg. I44 second part no. `findings.md`.
- GO5h: CMA hybrid confirmed on 36 grains. GO5g: hybrid introduced.
- **Next action:** plan and run GO5j (I45): alternate a field-gradient refit on coverage at the
  Gauss-Newton strain with Gauss-Newton (two rounds); S36, judged like GO5i against the CMA
  hybrid (strain, orientation at most 0.02 deg, half the time). Then GO6 or the cached ray index.

## Now as of GO5j close (2026-11-16 16:45 PST)

- **GO5j closed** (2026-11-16 16:45 PST): Gauss-Newton plus gradient refits (`--arm GA2@P`)
  gives strain 1.2 / 1.2 / 3.1e-4 and orientation 0.009 / 0.014 / 0.034 deg on S36 at 24 s per
  grain (CMA hybrid 2.3 / 1.5 / 2.1e-4, 0.011 / 0.009 / 0.017 deg, 53 s). 5e-4 fixed; at 1e-3
  Gauss-Newton escapes past the Tukey cut-off on grains with few spots (452, 0, 485, 349).
  I45 no. `findings.md` "GO5j".
- GO5i: Gauss-Newton 5x faster than CMA, gradient biased. GO5h: CMA hybrid confirmed.
- **Next action:** plan and run GO5k (I46): freeze inliers and weights per Gauss-Newton stage,
  at least half the spots inliers, fall back to the CMA hybrid for grains with few spots (and/or
  rebuild windows at the fitted strain); S36, judged like GO5i/j. Then GO6 or the cached ray index.

## Now as of GO5k close (2026-11-16 17:00 PST)

- **GO5k closed** (2026-11-16 17:00 PST): frozen Gauss-Newton weights plus a CMA-hybrid fallback
  for grains with fewer than 40 spots (`--arm GFF@P`) give strain 1.3 / 1.3 / 2.4e-4 on S36 at
  19 s per grain, orientation 0.009 / 0.015 / 0.028 deg (limit 0.02 at 1e-3). Two grains stall
  (all stages skipped, strain left at 0). I46 no. `findings.md` "GO5k".
- GO5j: gradient refit fixes 5e-4. GO5h: CMA hybrid confirmed (2.3 / 1.5 / 2.1e-4, 53 s).
- **Next action:** plan and run GO5l (I47): `GFF@P` plus fallback when Gauss-Newton skips its
  stages, plus a short joint hybrid CMA polish only for grains with fitted strain RMS > 6e-4;
  S36, judged like GO5i-k. Then GO6 or the cached ray index.

## Now as of GO5l close (2026-11-16 17:15 PST)

- **GO5l closed** (2026-11-16 17:15 PST): the fast strain fit `--arm GP4@P --cost hyb --lam 8`
  (Gauss-Newton on centroids, frozen Tukey weights, two coverage gradient refits, CMA-hybrid
  fallback for few spots or stalls, hybrid polish for strain above 6e-4) matches the CMA hybrid
  on S36: strain 1.3 / 1.3 / 2.1e-4, orientation 0.009 / 0.015 / 0.020 deg, 25 s per grain
  (CMA 53 s). I47 yes. `findings.md` "GO5l"; `COSTS.md`.
- GO5g-k: hybrid cost, confirmation on 36 grains, Gauss-Newton and its failure modes.
- **Next action:** plan and run GO5m (I48): the strain fits on noisy data (Poisson counts plus
  background on the clean intensity store, thresholded into a new binary stack and store); S36,
  `GP4@P` and the CMA hybrid against their clean results. Then GO6 or the cached ray index.

## Now as of GO5m close (2026-11-16 17:50 PST)

- **GO5m closed** (2026-11-16 17:50 PST): on images with the detector noise model
  (`--stack strainrn_clean.npy`, built from the `strain` sample's `realistic_nodrop` render)
  the CMA hybrid is robust (2.5 / 1.9 / 2.8e-4, clean 2.3 / 1.5 / 2.1e-4) and the fast fit is
  not (1.8 / 2.7 / 3.7e-4, clean 1.3 / 1.3 / 2.1e-4). I48 no. `findings.md` "GO5m".
- GO5l: fast fit `GP4@P` matches the CMA hybrid on clean images in half the time.
- **Next action:** plan and run GO5n (I49): noise-robust observed centroids for Gauss-Newton
  (binary weights; isolated lit pixels removed), fast fit on the noisy images against the noisy
  CMA hybrid. Then GO6 or the cached ray index.

## Now as of GO5n close (2026-11-16 18:15 PST)

- **GO5n closed** (2026-11-16 18:15 PST): binary-weight observed centroids (`--cen_w bin` /
  `binclean`) make the fast fit worse on noisy images (strain 4.6-6.7e-4 against 1.8 / 2.7 /
  3.7e-4 with intensity weights): the intensities carry the sub-pixel spot position. I49 no.
- GO5m: CMA hybrid robust to detector noise (2.5 / 1.9 / 2.8e-4), fast fit not.
- **Next action:** plan and run GO5o (I50): intensity weights with isolated lit pixels removed
  (`--cen_w intclean`), fast fit on noisy and on clean images, against the CMA hybrid on each.
  Then GO6 or the cached ray index.

## Now as of GO5o close (2026-11-16 18:35 PST)

- **GO5o closed** (2026-11-16 18:35 PST): intensity-weighted centroids with isolated lit pixels
  removed (`--cen_w intclean`) bring the noisy fast fit to 2.2 / 2.3 / 3.3e-4 (clean unchanged),
  but at 26.4 s per grain against the 25.3 s limit: I50 no, on time alone. Strain line
  time-boxed. Strain chart on the page (`strain_pareto.py` -> `strain_pareto.json`).
- Strain summary: CMA hybrid `cen` + 8 x coverage (clean 2.2 / 1.5 / 2.1e-4, noisy 2.5 / 1.9 /
  2.8e-4, 51 s); fast fit `GP4@P` (clean 1.3 / 1.3 / 2.1e-4, 25 s; noisy with `intclean` 2.2 /
  2.3 / 3.3e-4, 26 s). `COSTS.md`.
- **Program order** (owner via the coordinator, 2026-11-16; PLAN "GO6 onwards"): (1) G1b, (2) G2b
  (I23, I24), (3) G2c (I25-I27), (4) G2a (I18, I22). At most one more strain step unless a
  result demands it.
- **Next action:** plan and run GO4d (G1b, I36): cache the ray index per sample so the
  ray-indexing seed path is cheap, then the gated global search for voxels whose strict cost
  stays high after the ladder, on `rand` / `corr3` at 2 and 5 deg (GO4c: L/366/3 leaves 31% and
  72% of corr3 voxels above 0.1 deg there).

## Now as of GO4d close (2026-11-16 19:40 PST)

- **GO4d closed** (2026-11-16 19:40 PST, G1b): the ray pairing is cached per sample
  (`ray_cache.py`; `observed_rays.npz` in `corr3/` and `rand/`); `L/366/3,G/15/0.15` takes the
  share of voxels above 0.1 deg at 0.5 / 1 / 2 / 5 deg to 3 / 5 / 9 / 36% on corr3 (was 3 / 9 /
  31 / 72%) and 26 / 15 / 52 / 44% on rand (was 48 / 53 / 77 / 83%), 0.3 / 0.6 s per voxel. I36
  yes. `findings.md` "GO4d".
- Strain line time-boxed after GO5o (CMA hybrid; fast fit `GP4@P`; `COSTS.md`, strain chart).
- **Program order** (owner via the coordinator, 2026-11-16; PLAN "GO6 onwards"): (1) G1b done
  for now (I51 open), (2) G2b (I23, I24), (3) G2c (I25-I27), (4) G2a (I18, I22).
- **Next action:** plan and run GO7 (G2b, I23): rocking curve alpha > 0, spots spread over
  frames; first find whether the renderer and `SpotModel` support it, then render a sample with
  alpha > 0 and measure what the grain fit and the strict cost lose and whether sub-frame omega
  from frame ratios helps.

## Now as of GO7 close (2026-11-16 19:50 PST)

- **GO7 closed** (2026-11-16 19:50 PST, G2b I23 part 1): emulated rocking curve on the `grad`
  images (`rock_stack.py`, alpha 0.3 / 1.0 deg; `run_pipeline.py --stack`): the strict-field fit
  keeps its median (0.0149 / 0.0134 deg against 0.0146) but one grain (197) is rejected by the
  acceptance gate: 20 of 853 voxels above 0.1 deg against 3. Harm by the rule. Part 2 (model with
  a rocking curve) needs owner approval (library change); low priority.
- GO4d (G1b): cached ray pairing + gated seed path, I36 yes.
- **Program order** (owner via the coordinator, 2026-11-16): (1) G1b done for now, (2) G2b (I23
  part 1 done; I24 next), (3) G2c (I25-I27), (4) G2a (I18, I22).
- **Next action:** plan and run GO7b (G2b, I24): detector point-spread function on the images
  (emulated: blur the clean intensities with a Gaussian PSF, sigma 0.5 and 1 px, threshold), the
  same `Astrict/100` G20 test, plus the fast strain fit on the `strain` sample with the PSF.

## Now as of GO7b close (2026-11-16 20:00 PST)

- **GO7b closed** (2026-11-16 20:00 PST, G2b I24 images part): Gaussian PSF on the images
  (`psf_stack.py`): 0.5 px harmless to the fast strain fit; 1.0 px breaks it (strain 2.3 / 2.7 /
  5.7e-4, orientation up to 0.06 deg) and raises the orientation-fit median 1.3 x. Grain 197
  again lost to the strict-cost gate. `findings.md` "GO7b".
- GO7: rocking-curve images cost one grain via the gate. GO4d (G1b): I36 yes.
- **Program order** (owner via the coordinator, 2026-11-16): (2) G2b in progress (I52 next, then
  model-side I23/I24 parts need owner approval), (3) G2c (I25-I27), (4) G2a (I18, I22).
- **Next action:** plan and run GO7c (I52): PSF-aware predicted coverage in the strain costs
  (`fit_strain.py`, predicted coverage points spread by the PSF Gaussian, Monte Carlo copies),
  fast fit on `strain_psf10`, judged against the clean fast fit. Then G2c.

## Now as of GO7c close (2026-11-16 20:10 PST)

- **GO7c closed** (2026-11-16 20:10 PST, G2b I52): PSF-aware predicted coverage does not rescue
  the strain fit on 1 px PSF images (2.3 / 2.6 / 6.8e-4); Gauss-Newton-only grains degrade as
  much, so merged spots (contaminated centroids) are the cause (I53 open). G2b time-boxed here.
- GO7b: 1 px PSF breaks the fast strain fit, 0.5 px harmless. GO7: rocking costs one grain via
  the gate. GO4d (G1b): I36 yes.
- **Program order** (owner via the coordinator, 2026-11-16): (3) G2c (I25-I27) next, then (4)
  G2a (I18, I22). Model-side rocking / PSF (library changes) wait for owner approval.
- **Next action:** scope and run the first G2c step (GO8). I25: `GNLayerNet`
  (`docs/orientation_nn_design.md`; learned per-voxel Gauss-Newton refiner, offsets up to 1 deg;
  no trained weights in this worktree) as a per-voxel refiner inside grains on `rand` / `corr3`
  at 0.5 and 1 deg after the ladder; first check where the model code and training data live,
  and the training cost (each run under 2 h).

## Now as of GO9b close (2026-11-16 20:35 PST)

- **GO9 / GO9b closed** (2026-11-16 20:35 PST, G2a): Poisson deviance with free per-reflection
  scales is the closest smooth objective at the ground truth (0.0101 / 0.0060 deg on grains 381 /
  302; strict 0.0060 / 0.0082), missing the rule by 0.0001 deg (I18 no); tying the scales to the
  structure factor makes grain 302 worse (0.030; I22 no). `probe_poisson.py`.
- GO7-GO7c (G2b): rocking costs a grain via the gate; 1 px PSF breaks the strain fit (merged
  spots). GO4d (G1b): cached ray pairing, I36 yes.
- **Owner item:** G2c (I25 GNLayerNet) needs a dataset, training and an implementer (approval in
  DECISIONS.md); asked 2026-11-16 20:15 PST, options sent. Waiting.
- **Next action:** while G2c waits, plan and run GO4e (G1b, I51): after the gated seed path, a
  neighbour continuation from the improved voxels (`L/366/3,G/15/0.15,N/30/0.3`) and/or a gate on
  neighbour disagreement, on corr3 / rand (GO4d left rand 2 deg at 52%).

## Now as of GO4e close (2026-11-16 20:35 PST)

- **GO4e closed** (2026-11-16 20:35 PST, G1b I51): a neighbour pass after the gated seed path
  gains 1-4 points everywhere (best pipeline now `L/366/3,G/15/0.15,N/30/0.3`: corr3 2 / 4 / 9 /
  33%, rand 25 / 13 / 48 / 41% above 0.1 deg at 0.5 / 1 / 2 / 5 deg); a second seed pass adds
  nothing. I51 no.
- GO9/GO9b (G2a): Poisson deviance best smooth objective, I18 / I22 no. GO7-7c (G2b).
- **Owner item:** G2c (I25 GNLayerNet) needs a dataset, training and an implementer (approval in
  DECISIONS.md); asked 2026-11-16 20:15 PST. Waiting.
- **Next action:** GO4f (G1b diagnosis): on rand at 2 deg after the best pipeline, for the voxels
  still wrong: were they gated, did the seed path return a candidate, and is the strict cost at
  the ground truth lower than at their end (search failure or objective failure)?

## Now as of GO4f close (2026-11-16 20:45 PST)

- **GO4f closed** (2026-11-16 20:45 PST, G1b diagnosis): every voxel still wrong after the best
  G1b pipeline was gated, and the strict cost at the ground truth is far lower than at its end
  (rand 2 deg: 0.079 against 0.393): search failures; the seed path's candidates miss the ground
  truth. `go4f_diag.py`.
- GO4e: neighbour pass after the seed path, +1-4 points. GO9/9b (G2a), GO7-7c (G2b) closed.
- **Owner item:** G2c (I25 GNLayerNet) waits for the owner (asked 2026-11-16 20:15 PST).
- **Next action:** plan and run GO4g (I54): for gated voxels the seed path leaves above the gate,
  a wider search: all 20 indexing candidates through the strict cost and the finisher on 5, then
  the reconstructor's grid-search fallback (`hooks.grid_seed`); rand / corr3, judged against
  GO4e's best pipeline.

## Now as of GO4g close (2026-11-16 21:55 PST)

- **GO4g closed** (2026-11-16 21:55 PST, G1b I54 yes): after `L/366/3,G/15/0.15,N/30/0.3`, the
  all-candidates seed path `W/20/0.15` takes rand to 25 / 5 / 21 / 14% and corr3 to 2 / 3 / 3 /
  21% above 0.1 deg (0.5 / 1 / 2 / 5 deg); adding the grid search `X/1/0.15` (rand only) gives
  17 / 2 / 7 / 6%, wrong voxels 229 -> 25 of 1,043, at 4 s per voxel. `findings.md` "GO4g".
- GO4f: remaining failures were search failures. GO9 / GO9b (G2a), GO7-7c (G2b) closed.
- **Owner item:** G2c (I25 GNLayerNet) waits for the owner (asked 2026-11-16 20:15 PST).
- **Next action:** plan and run GO4h (I55): `X` on corr3 and a per-voxel strict polish after
  `W,X` (e.g. `V/61/0.02`) on both samples, aimed at rand 0.5 deg (17%, mostly 0.1-1 deg);
  judged against GO4g.

## Now as of GO4h close (2026-11-18 16:45 PST)

- **GO4h closed** (2026-11-18 16:45 PST, G1b I55 no): full search + per-voxel polish ends corr3 at
  1 / 1 / 0.4 / 7% and rand at 16 / 1.4 / 5 / 4% above 0.1 deg (0.5 / 1 / 2 / 5 deg), medians
  0.005-0.009 deg; rand's 0.5 deg near misses (0.1-0.75 deg) sit below the strict-cost gate.
  Best G1b pipeline: `L/366/3,G/15/0.15,N/30/0.3,W/20/0.15,X/1/0.15,V/61/0.02` (2-4 s per voxel).
- Machine restart on 2026-11-17; GO4h had finished. Wake job `dc9287aa` survived.
- **Program order** (owner via the coordinator, 2026-11-18): I53 (split merged spots) next, then
  I25 GNLayerNet at low priority; the rest of G2c (I26, I27) skipped. The owner's I25 choice
  came relayed; it is not in DECISIONS.md (the logging script was refused for a relayed
  approval), so ask the owner to confirm in their own words before launching an implementer.
- **Next action:** plan and run GO7d (I53): split observed connected components that a PSF merged,
  using the predicted spot positions (assign each lit pixel to the nearest predicted spot of the
  grain or its neighbours) before taking centroids; fast strain fit on `strain_psf10`, judged
  against GO7b (2.3 / 2.7 / 5.7e-4) and clean (1.3 / 1.3 / 2.1e-4).

## Now as of GO7e close (2026-11-18 17:30 PST)

- **GO7d / GO7e closed** (2026-11-18 17:30 PST, G2b): splitting merged spots does not fix the 1 px
  PSF damage (I53 no); the cause is the region area-ratio filter (0.8-1.6) dropping blurred
  spots. Widening it to 2.5 (`--ratio_hi 2.5`) helps on PSF images (2.1 / 2.1 / 3.5e-4) but misses
  the limits (I56 no), and improves the clean fast fit to 1.3 / 1.2 / 1.4e-4, orientation <= 0.013
  deg (was 2.1e-4 / 0.020 at 1e-3), 35 s per grain. G2b stops here.
- GO4h (G1b): full search + polish, corr3 1 / 1 / 0.4 / 7%, rand 16 / 1.4 / 5 / 4%.
- **Program order** (owner via the coordinator, 2026-11-18): I53 done; I25 GNLayerNet next at low
  priority (the relayed choice is not in DECISIONS.md; confirmation in the owner's own words was
  asked for 2026-11-18 16:45 PST before any implementer). G2c otherwise skipped.
- **Next action:** GO5p (I57, one strain check the GO7e result calls for): the fast fit with
  `--ratio_hi 2.5` on the noisy images (`strainrn`), to adopt it as the default; then I25 scoping
  in this session (dataset and training cost; no implementer until the owner confirms).

## Now as of GO5p close (2026-11-18 17:50 PST)

- **GO5p closed** (2026-11-18 17:50 PST, I57 no, narrowly): on noisy images the fast fit with the
  ratio filter at 2.5 gives 1.6 / 2.6 / 3.1e-4 (CMA hybrid 2.5 / 1.9 / 2.8e-4), missing 1.25 x
  the CMA hybrid at 5e-4. Strain line closed; `--ratio_hi 2.5` recommended for clean data.
- GO7d / GO7e (G2b) closed. GO4h (G1b): corr3 1 / 1 / 0.4 / 7%, rand 16 / 1.4 / 5 / 4%.
- **I25 scoping** (2026-11-18): no GNLayerNet weights or datasets exist on this machine; I25 needs
  the 3.2 GB dataset regenerated and a 60-epoch training. The niche it might fill: the
  realism-trained net was 2-3 x better than robust Gauss-Newton on windows with overlapping
  neighbours (`docs/orientation_nn_design.md`), the situation of rand's near misses. Waiting for
  the owner's own confirmation (asked 2026-11-18 16:45 PST) before an implementer.
- **Next action:** GO4i (I58) running; close it, then I25 if confirmed, else report and pick the
  next G1b or G2a idea from IDEAS.
