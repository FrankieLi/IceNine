# AR6 (`grainopt`): full-grain fitting research program

Owner's request (Frankie, 2026-11-16), in the owner's words where possible. Branch
`feature/grain-opt`, worktree `.claude/worktrees/grain-opt`. Separate from Newt's autoresearch
(`docs/autoresearch/`): never edit those files or the main checkout.

## Goals

- **G1. Find an algorithm that reaches the optimum of the objectives we have.**
  - **G1a.** For a fixed objective, sketch the Pareto frontier of cost (evaluations) against
    quality, across algorithms and settings.
  - **G1b.** Robustness against the kind of gradient: from linear gradients (easy) to large,
    fairly random voxel-to-voxel misorientations (0 < d < 5 deg; plastically deformed samples
    with geometrically necessary dislocations).
  - **G1c.** Sensitivity to elastic strain (lattice deformation): can it be captured and
    recovered? Reference to look up (owner, 2026-11-16): HEXRD far-field HEDM strain methods,
    https://hexrd.github.io (local repo `~/Research/github/hexrd`, `fitgrains.py`). Far field
    puts the detector far away, so a grain covers only a few pixels; strain resolution depends
    only on the ratio of detector pixel size to sample-detector distance.
- **G2. Better objectives to guide the reconstruction.**
  - **G2a.** Different ways to use intensity.
  - **G2b.** More sophisticated diffraction simulation (rocking curve alpha > 0, partial-volume
    coverage, strain in the forward model, detector response).
  - **G2c.** More sophisticated models, including neural networks (reuse
    `docs/orientation_nn_design.md`, `docs/nn_inverse_problem_formulation.md`,
    `docs/research_ideas_nerf.md`, `docs/todo_intensity_and_distribution_losses.md`,
    `docs/todo_hybrid_nn_refinement.md`, and shmem).

## Standing judging rule (applies to every step unless its section overrides it)

Written before any run of the step.

- **Benchmark G20**: the 20 `grad` grains picked by `run_grains.pick_true_grains(T, 20)` (4 per
  gradient level 0, 0.1, 0.2, 0.5, 1.0 deg; sizes at quantiles 0.1 to 0.95 of the level). The
  fitted group of a true grain is the start-map flood-fill group (2 deg) holding most of its
  voxels. Start map: the painted map (AR3 step 3, `runs/P/maps.npz` "P").
- **Quality**: per-voxel symmetry-reduced misorientation to the voxel's own ground truth: median, 95th
  percentile, share above 0.1 deg, wrong (above 1 deg), with Wilson 95% intervals. Boundary
  voxels (a neighbour outside the group) reported separately. The strict binary cost
  (`VoxelCostFunction`, local, hard) is the referee and is reported per voxel.
- **Cost**: objective evaluations per grain (one evaluation renders every voxel of the group
  once) and wall seconds (labelled single-worker or contended by `preflight_label`). For a
  per-voxel method, strict evaluations per voxel. The Pareto x axis is "objective evaluations
  per voxel": one evaluation of a grain objective scores every voxel of the grain once (counts 1
  per voxel); a per-voxel fit counts its strict-cost calls per voxel; chained stages add up.
- **Reference points**: start (0), AR3 step 4 CMA per-voxel fit `C` (61 strict evaluations per
  voxel; on the painted voxels 334 of 19,013 above 0.1 deg, 1.6 to 2.0%), on the same voxels.
- **Comparisons**: a gain or loss only when Wilson intervals do not overlap. No p-values, no
  extra seeds unless the step says so. One run per arm.
- **Acceptance inside a fit**: a grain update is kept only if the objective and the mean strict
  cost of its voxels do not get worse.
- **Stop rule**: a run stops at 2x its time budget; a step at 2 h of compute unless a decision
  entry raises it. Check `uptime` before runs (Newt's autoresearch shares the machine; at most
  4 heavy processes for this track).

## Steps

Each step gets its own section below before it runs: question, arms, judging rule overrides,
budget, what is committed. Ideas live in `IDEAS.md`; a step takes one or more ideas.

### GO1: objective bias (ideas I1, I2)

Question: where does each objective's minimum sit relative to the ground truth? Level-0 grains of
`grad` (one true orientation each): grains 381 and 302. CMA over the grain rotation only
(3 numbers), started AT the ground truth, 150 evaluations, sigma 0.03 deg. Arms: unbalanced Sinkhorn
(reach 3 px) at blur 1 and 0.5 px, balanced at blur 1, with the pixel-box dilation of the
prediction at blur 1 and 0.5, intensity-weighted, least squares. Reported: median offset of the
optimum from the ground truth (deg), objective and mean strict cost at ground truth and optimum. Rule: an
objective is "unbiased enough" if its offset is below 0.01 deg on both grains (the start's median
error on level-0 grains is 0.016 deg, on all painted voxels 0.041 deg). If no arm qualifies, the
next step studies the cause (I3, I4) before any optimizer study.

### GO2: optimizer Pareto on G20 (ideas I5, I6, I9, I10, I28, I29)

Question: for each objective, how many objective evaluations per voxel buy how much accuracy, and which chain
of stages reaches the strict-cost optimum cheapest? GO1 found every Sinkhorn variant flat within
0.03 to 0.04 deg of the ground truth, and the strict cost (grain mean) within 0.006 to 0.008 deg. So
Sinkhorn is a basin finder at best, and the finish must use the strict cost (or a sharper smooth
cost, I3).

Benchmark G20 from the painted start, start-map flood-fill groups. Arms (pipelines of
`run_pipeline.py`; A = grain-field CMA, B = per-voxel Adam with smoothness, V = per-voxel CMA on
the strict cost):

1. Field CMA on the strict cost: `Astrict/50`, `/100`, `/200`, `/400` (I28).
2. The same with the quadratic field: `Astrict/200/o2`.
3. Field then per-voxel polish: `Astrict/100,V/10`, `Astrict/100,V/30`, `Astrict/100,V/61` (I10).
4. Per-voxel only, from the start: `V/61` (AR3 C re-run on the same voxels, I9) and `V/20`.
5. Sinkhorn first: `Asink/150,Astrict/100,V/20`; and `Asink/150,Bsink/40,Astrict/50` (I6).
6. Least squares first: `Als/150,Astrict/100,V/20` (I21 as basin finder).

Judging: the standing rule. Pareto axes: objective evaluations per voxel (x) against share above 0.1 deg and
median error (y), on the voxels of the 20 grains with a C record; also wall seconds per voxel.
The frontier is the set of arms no other arm beats on both axes (with interval overlap counted as
a tie). Budget: strict-cost arms take seconds per grain; Sinkhorn arms about 1 min per grain;
whole step under 1 h of compute. Committed: run.json and voxels.json.gz per arm (each under
1 MB), tables.md, Pareto rows.

### GO3: deformed samples (ideas I11, I12) -- built

`rand` (white-noise voxel rotations) and `corr3` (Gaussian-smoothed, length 3 voxel sides), RMS
levels 0.5, 1, 2, 5 deg drawn per grain (seed 0). Rendered clean with the stock forward model;
coverage mass on lit pixels at the ground truth 97 to 99% on 10 grains of each. In `rand`, 69% (2 deg)
and 97% (5 deg) of neighbour pairs exceed 2 deg, so flood-fill grouping fragments those grains:
GO4 uses true grain ids (`--group true`) and the grain-centre start (`--start centre`).

### GO4: robustness across kinds of orientation variation (ideas I13, I6)

Question: does the GO2 frontier (grain field on the strict cost, then a per-voxel polish) hold
when the voxel-to-voxel variation is random or spatially correlated, from 0.5 to 5 deg RMS?
Samples `rand` (white noise) and `corr3` (correlation length 3 voxel sides), GO3. Start: each
voxel at its grain's orientation before deformation (`--start centre`, a perfect paint that knows
nothing of the variation); membership: true grain ids (`--group true`; flood fill fails at 2 and
5 deg). Benchmark: 20 grains per sample by `pick_true_grains` (5 per level 0.5, 1, 2, 5 deg).

Arms (per sample): `Astrict/100`; `Astrict/100,V/61`; `Astrict/100,V/61/0.3` (per-voxel step
0.3 deg); `V/61/0.3`; `Als/150,Astrict/100,V/61/0.3` (least squares as basin finder).
Judging: the standing rule, reported per level (median, p95, above 0.1 deg, wrong above 1 deg,
Wilson intervals); no C reference exists on these samples, so all voxels of the groups count.
Expected: the field cannot follow white noise; the question is which polish recovers which level.
Budget: strict arms seconds per grain, the LS arm about 10 s per grain; under 1 h in total.

### GO4b: reaching far starts (ideas I33, I34)

Question: which per-voxel strategy recovers deformed grains whose voxels start 0.5 to 5 deg from
their ground truth? Same samples, start, membership and 20-grain benchmarks as GO4. New stages in
`run_pipeline.py`:

- `N/<evals>/<sigma>` neighbour continuation (I33): voxels ordered by strict cost (best first);
  each voxel, in breadth-first order from the best 10%, scores its own orientation and the current
  orientations of its fitted neighbours in the grain (one strict evaluation each), keeps the best,
  then runs per-voxel CMA (`evals`, step `sigma`). Two sweeps. Every evaluation counts.
- `L/<evals>` step ladder (I34): per-voxel CMA in three legs of `evals`/3 each at steps 1.0, 0.3
  and 0.05 deg, each leg from the best so far.

Arms per sample: `L/61`, `L/183`, `N/30/0.3`, `L/61,N/30/0.3`, `Als/150,Astrict/100,L/61,N/30/0.3`.
Judging: the standing rule per level, against GO4's best arm per level. Budget: under 1 h.

### GO4c: larger ladders (idea I37; I36 blocked)

Same setup as GO4b. Arms per sample: `L/244/3` (legs 3, 1, 0.3, 0.05 deg, 61 each), `L/366/3`
(91 each), `L/183,N/30/0.3`. Judging as GO4b, against L/183. The gated global search (I36, stage
`G`) was tried first: the reconstructor's seed path did not finish one 8-voxel `rand` grain in
15 minutes (it pairs 1.56M spot fragments across the detectors, about 25% of memory), so it is
blocked until the ray index is cached per sample.

### GO5: strain in the forward model (ideas I14, I15, I16, I30)

Implemented by an `implementer` agent in its own worktree (branch `feature/grain-opt-strain` off
`feature/grain-opt`), reviewed before merge. Conventions from `findings_cpp_strain.md`: crystal-
frame symmetric strain eps; M = (I + eps)^-1; g = O M g0; Bragg angle and omega from the strained
|g|; .mic columns 13-18 hold M (19-column layout), written with at least 9 significant digits.

1. `SpotModel` (grain_fit.py): optional per-voxel strain (n, 6) [e11 e22 e33 e23 e13 e12],
   default None, bit-identical output when None. Differentiable in eps.
2. Tests: pure dilation eps = d I gives |g| = |g0|/(1+d) and smaller 2-theta; a small shear
   moves spots as finite differences predict; None equals the current output exactly.
3. Stock render path (I15): opt-in per-voxel deformation from the .mic in
   `forward_simulation.py`, default off and bit-identical; strained-g Bragg filter uses the
   strained |g| (the C++ bug is not copied). Test against SpotModel spot positions.
4. Sample generator `make_strained_sample.py`: the `grad` truth plus a per-grain constant strain
   tensor of RMS 1e-4, 5e-4 or 1e-3 (levels drawn per grain, seed 0); render one sample
   (`strain`); store as GO3.
5. Sensitivity (I16, I30): on 5 grains of `strain`, curvature of the mean strict cost and of the
   intensity-weighted Sinkhorn in each strain component at the ground truth (finite differences),
   expressed as the strain step that raises the cost by 10%; compare with pixel/distance
   (4.4e-4, 2.7e-4 for the two detectors).

Budget: implementation about 2 to 3 h of agent time; runs under 30 min. Committed: code, tests,
the generator, sample records (not the .mic), a short findings section.

### GO5b: per-grain strain fit (idea I38)

Question: can a per-grain fit recover the elastic strain of sample `strain` (GO5) from clean
images, at tensor RMS 1e-4, 5e-4 and 1e-3, and does it still work from the painted start?

Setup. Benchmark S12: 12 grains of `strain`, 4 per strain level, sizes at quantiles 0.2, 0.4,
0.6, 0.8 of the level (true grain ids as membership). Cost: the coverage cost of GO5 (1 minus the
share of predicted coverage mass on lit pixels, mean over the grain's spots), the sharper
strict-like cost (`VoxelCostFunction` has no strain). Parameters: the grain field (`FieldA`
order 1, 9 numbers) and one crystal-frame strain per grain (6 numbers, [e11 e22 e33 e23 e13
e12]). The spot windows are built without strain (no ground-truth strain leaks in). Optimizer:
CMA, per-coordinate steps 0.05 deg (field) and 5e-4 (strain), strain starting at 0.

Arms (script `fit_strain.py`, output `benchmarks/grain_opt/go5b/`):

- `E@GT`: orientations held at the ground truth; CMA on the 6 strain numbers, 300 evaluations
  (can the cost see strain at all).
- `FE@GT`: field + strain from the field fitted to the ground-truth orientations, 600 evaluations.
- `FE@P`: from the painted start "P": a field-only CMA on the coverage cost (200 evaluations, step
  0.05 deg), then field + strain (600 evaluations).

Measures per grain: strain error = tensor RMS of (fitted minus ground-truth strain), with the
ground truth re-expressed in the fitted symmetry representative; per-component errors; cost at
the fit and at the ground truth (ground-truth orientations and strain); per-voxel misorientation
to the ground truth (median, p95) for the field arms.

Judging rule (written before any run). Per level and arm, the median strain error over its 4
grains: "recovered" when at most 0.5 x the level (twice better than answering zero), "partial"
between 0.5 and 1.0 x, "not recovered" at or above the level. A grain whose cost at the fit is
above its cost at the ground truth is an optimizer shortfall; one with a lower cost but a large
error (above the level) is an objective shortfall (the cost optimum is not at the ground truth).
`FE@P` against `FE@GT`: worse at a level only when its median error is above 1.5 x and at least
2 of the 4 grains are worse. Reference scale: pixel / distance 2.7e-4 to 4.4e-4 (GO5: 1e-3 strain
moves spots 0.2 to 0.5 px). One run per arm, seed 0; no p-values.

Budget: under 1 h of compute (measured per evaluation before the full run); at most 3
processes (one per arm). Committed: the script, `go5b/*.json` (per-grain rows), findings.

### GO5c: joint strain fit that reaches the optimum (idea I39)

Question: GO5b's joint fits (field + strain, 15 numbers) stopped above the ground truth's cost in
19 of 24 grains. Does a staged or longer search reach the optimum and bring the painted-start
strain error down to the strain-only floor of `E@GT` (median 1.7e-4 / 1.9e-4 / 2.0e-4)?

Same benchmark S12, cost, robust painted start and measures as GO5b (`fit_strain.py`, output
`benchmarks/grain_opt/go5c/`). Arms:

- `J@P`: painted start, field leg 200, then one joint leg of 2400 evaluations (budget control).
- `S@P`: painted start, field leg 200, strain-only leg 300 (field fixed), joint leg 600 with
  smaller steps (0.02 deg, 2e-4): 1100 evaluations.
- `SR@P`: `S@P`, then two more joint legs of 300 from the best point, steps halved each leg:
  1700 evaluations.
- `S@GT`: `S@P`'s strain and joint legs from the field fitted to the ground truth (no field leg).

Judging rule (written before any run). Per level and arm, median strain error over 4 grains.
"Reaches the floor" at a level when the median is at most 1.5 x `E@GT`'s median at that level.
Better than GO5b `FE@P` at a level when the median is below `FE@P`'s / 1.5 and at least 3 of the
4 grains are better. "Short" (final cost above the ground truth's) is counted per arm; the step
answers I39 when an arm has at most 4 of 12 short and reaches the floor at 5e-4 and 1e-3. The
1e-4 level is reported, not judged (below the floor, I40). One run per arm, seed 0.

Budget: about 30 ms per evaluation, under 30 min in total; 4 processes at most.

### GO5d: fix the orientation, then strain alone (idea I41)

Question: GO5c showed the coverage cost trades a 0.01 deg rotation of the grain field against
2-4e-4 strain, so a joint fit lands at a lower cost with the wrong strain. If the field is fixed
first from a fit that does not fit strain, and only strain is then fitted, does the painted-start
strain error come down to the strain-only floor of `E@GT` (median 1.7e-4 / 1.9e-4 / 2.0e-4)?

Same benchmark S12, coverage cost, robust painted start and measures as GO5b/GO5c
(`fit_strain.py`, output `benchmarks/grain_opt/go5d/`). In every arm the last leg is strain only
(300 evaluations, step 5e-4, from 0) at a fixed field; `mis_fit_med` is the fixed field's
orientation error. Arms:

- `E@F`: the field fitted to the ground-truth orientations (least squares, no CMA), then strain
  only. Control: the floor when the orientation is limited to an order-1 field.
- `O@P`: painted start, field leg on the coverage cost with strain fixed at 0 (600 evaluations,
  step 0.05 deg), then strain only.
- `Ost@P`: painted start, field leg on the strict cost (no strain model; 200 evaluations, step
  0.05 deg), then strain only on the coverage cost.

Judging rule (written before any run). Per level and arm, median strain error over 4 grains.
"Reaches the floor" at a level when the median is at most 1.5 x `E@GT`'s median at that level.
Better than GO5c `S@P` at a level when the median is below `S@P`'s and at least 3 of the 4 grains
are better. I41 is answered yes when a painted-start arm reaches the floor at 5e-4 and 1e-3. If
`E@F` itself misses the floor, the order-1 field is the limit and the painted arms are judged
against `E@F` instead (same 1.5 x). The 1e-4 level is reported, not judged. One run per arm,
seed 0.

Budget: about 30 ms per coverage evaluation, strict about 25 ms per voxel; under 30 min in
total; 3 processes.

### GO5e: a sharper sub-pixel strain cost (idea I40)

Question: the coverage cost is flat over about 1e-4 to 3e-4 strain (GO5d: at the exact
orientation, repeat strain-only fits scatter in that band). Does a cost that sees sub-pixel
shifts of the predicted coverage bring the strain-only error below that band, and then the
painted-start staged fit with it?

Costs, on the same predicted coverage points (continuous column, row and frame positions, mass)
as the coverage cost, per spot in its window, mean over spots:

- `sd` (signed distance): a 3D signed distance map of the lit pixels in the window (frames,
  rows, columns; distance in pixel or frame units, positive outside the lit set, negative inside,
  from Euclidean distance transforms), read by trilinear interpolation at each point, clipped to
  [-3, 5]; cost = mass-weighted mean. Lowest when the predicted mass sits deepest inside the
  observed spot.
- `int` (intensity): the clean intensity of the lit pixels (0 elsewhere), read by trilinear
  interpolation; cost = 1 - mass-weighted mean intensity / the spot's top intensity.

Benchmark S12 and measures as GO5b to GO5d (`fit_strain.py --cost`, output
`benchmarks/grain_opt/go5e/`). Arms:

- `E@GT` (ground-truth orientations, strain only, 300 evaluations) with `sd` and with `int`,
  each with seeds 0 and 1 (repeat scatter).
- `S@P` (painted start: field leg 200 on the same cost at zero strain, strain only 300, joint
  600 at steps 0.02 deg and 2e-4) with the cost that wins the first part, seed 0. Not run if
  neither cost is sharper.

Judging rule (written before any run). Reference: the coverage cost's floor, the mean of the
`E@GT` (GO5b) and `E@F` (GO5d) medians: 1.66e-4 / 1.57e-4 / 2.67e-4 at 1e-4 / 5e-4 / 1e-3; its
repeat scatter (median over grains of the absolute difference between the two runs) 3.5e-5.

- A cost is "sharper" at a level when the mean of its two seeds' medians is at most the
  reference / 1.5 at that level. I40 is answered yes for a cost that is sharper at 2 or more
  levels, one of them 1e-4. The winner is the cost sharper at more levels (ties: lower mean of
  medians over levels). The repeat scatter is reported, not judged.
- The staged fit with the winner "reaches the new floor" at a level when its median is at most
  1.5 x the winner's `E@GT` mean of medians there, and beats the coverage `S@P` (7.1e-4 /
  3.4e-4 / 3.0e-4) at a level when its median is lower and 3 of 4 grains are better. Judged at
  5e-4 and 1e-3; 1e-4 reported.

Budget: the maps are built once per grain; an evaluation adds one interpolation per spot (a few
ms over the coverage cost's 30 ms). Five runs of 2 to 6 minutes; 3 processes at most.

### GO5f: a two-sided sub-pixel strain cost (idea I42)

Question: GO5e's one-sided sub-pixel costs were sharp but biased (they reward shrinking the
predicted footprint into the spot). Does a cost that compares the predicted spot with the
observed one from both sides bring the strain-only error at the ground-truth orientation below
the coverage band, and then the painted-start staged fit?

Costs, on the same predicted coverage points as before (continuous column, row and frame; frame
in frame units), per spot, mean over spots:

- `cen` (centroid): squared distance between the predicted mass centroid and the observed
  clean-intensity centroid of the spot's lit pixels (pixel centres), truncated at 4 (2 px) so a
  spot merged with a neighbour's cannot dominate. Size changes do not pay; position does.
- `ch` (chamfer): the positive part of GO5e's signed distance, averaged over the predicted mass
  (predicted mass off the lit set), plus the intensity-weighted mean distance from each lit pixel
  centre to the nearest predicted point minus its value at the ground truth's point spacing
  (observed pixels the prediction leaves uncovered), both clipped at 5.

Benchmark S12 and measures as GO5b to GO5e (`fit_strain.py --cost`, output
`benchmarks/grain_opt/go5f/`). Arms: `E@GT` with `cen` and with `ch`, seeds 0 and 1; then `S@P`
(field leg 200 at zero strain, strain 300, joint 600 at 0.02 deg and 2e-4) with the winner,
seed 0, only if a cost is sharper.

Judging rule (written before any run; the same as GO5e). Reference: coverage `E@GT` / `E@F`
mean of medians 1.66e-4 / 1.57e-4 / 2.67e-4 at 1e-4 / 5e-4 / 1e-3.

- A cost is "sharper" at a level when the mean of its two seeds' medians is at most the
  reference / 1.5. I42 is answered yes for a cost sharper at 2 or more levels, one of them 1e-4.
  Winner: sharper at more levels (ties: lower mean of medians over levels). Repeat scatter and
  "short" are reported, not judged.
- The staged fit with the winner "reaches the new floor" at a level when its median is at most
  1.5 x the winner's `E@GT` mean of medians, and beats the coverage `S@P` (7.1e-4 / 3.4e-4 /
  3.0e-4) at a level when its median is lower and 3 of 4 grains are better. Judged at 5e-4 and
  1e-3; 1e-4 reported.

Budget: `cen` costs about the coverage cost's 30 ms per evaluation; `ch` adds a nearest-point
search per spot. Five runs of 2 to 8 minutes; 3 processes at most.

### GO5g: hybrid cost, coverage for orientation and centroids for strain (idea I43)

Question: GO5f's centroid cost recovers strain well from the painted start at 1e-4 and 5e-4 but
lets the orientation field drift (0.06 deg; coverage keeps 0.012 deg). Does a hybrid that keeps
coverage for the orientation field and centroids for strain get both: strain near the `cen`
floor and orientation near coverage's?

Cost `hyb`: `cen` + lambda x coverage (GO5f's `cen` already truncates merged spots at 4 px^2; no
further down-weighting in this step). At the ground truth a 2e-4 strain step changes `cen` by
about 0.01 and coverage by about 0.005, so lambda = 2 weighs them about equally. Benchmark S12
and measures as before (`fit_strain.py --cost hyb --lam`, output `benchmarks/grain_opt/go5g/`).
Arms, seed 0:

- `E@GT` `hyb` lambda 2: strain only at the ground-truth orientation (the hybrid's floor).
- `S@P` `hyb` lambda 2: GO5c staged protocol, every leg on the hybrid.
- `S@P` `hyb` lambda 8: the same, coverage weighted more.
- `Sx@P` (by role): field leg 200 on coverage (zero strain), strain-only leg 300 on `cen`,
  joint leg 600 on `hyb` lambda 2 (steps 0.02 deg, 2e-4).

Judging rule (written before any run). Per level and arm, median strain error over 4 grains and
median voxel orientation error. I43 is answered yes for a painted-start arm whose strain median
is at most 1.5 x `cen` `E@GT` (GO5f mean of medians 1.27e-4 / 2.22e-4 / 1.68e-4, so 1.9e-4 /
3.3e-4 / 2.5e-4) at 5e-4 and 1e-3, and whose orientation median is at most 0.02 deg at every
level. Better than GO5f `cen` `S@P` at 1e-3 when its median there is below 3.7e-4 and 3 of 4
grains are better. 1e-4 strain is reported, not judged. Best arm: the one meeting the rule with
the lowest mean strain median over 5e-4 and 1e-3.

Budget: about 60 ms per hybrid evaluation; four runs of 2 to 8 minutes; 3 processes at most.

### GO5h: confirm the hybrid strain fit on more grains (idea I44)

Question: GO5g's hybrid (`cen` + 8 x coverage, staged from the painted start) met its rule on 4
grains per level with a small margin at 1e-3. Does it hold on 12 grains per level, and is a
heavier coverage weight (lambda 16) better?

Benchmark S36: 12 grains per strain level, sizes at quantiles 1/13 to 12/13 of the level
(`fit_strain.py --per_level 12`; it overlaps S12). Painted start, robust start field, GO5c staged
protocol (field leg 200 at zero strain, strain only 300, joint 600 at 0.02 deg and 2e-4), seed 0.
Output `benchmarks/grain_opt/go5h/`. Arms:

- `S@P` coverage (baseline on the same grains).
- `S@P` `hyb` lambda 8.
- `S@P` `hyb` lambda 16.

Judging rule (written before any run). Per level and arm, median strain error over 12 grains and
median voxel orientation error. A hybrid arm "confirms" when at 5e-4 and at 1e-3 its strain
median is at most the coverage arm's median / 1.5 on the same grains with at least 8 of 12 grains
better, and its orientation median is at most 0.02 deg at every level. I44 (first part) is
answered yes when a hybrid arm confirms; the better lambda is the confirming arm with the lower
mean strain median over 5e-4 and 1e-3. 1e-4 strain is reported, not judged. Grains whose
orientation ends above 0.05 deg are counted and named.

Budget: about 60 ms per hybrid evaluation, 1,100 evaluations per grain: about 25-35 minutes per
hybrid arm, 15 for coverage; 3 processes in parallel, under 45 minutes.

### GO6 onwards

Planned order (each gets a section when it starts): GO6 objective variants (I17 to I22); GO7 richer
simulation (I23, I24); GO8 learned models (I25 to I27).
