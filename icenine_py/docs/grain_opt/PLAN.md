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

### GO5i: Gauss-Newton on centroid residuals (idea I44, second part)

Question: the hybrid strain fit (GO5h) costs about 50 s per grain with CMA (1,100 evaluations).
Per-spot centroid residuals are smooth, low-dimensional least squares: does Gauss-Newton on
them, with the orientation field's gradient taken from the coverage field leg, match the hybrid's
accuracy at a fraction of the time?

Method: residuals r = predicted mass centroid - observed intensity centroid per spot (column,
row, frame; 3 per spot), over the spots that have predicted mass at the start. Parameters: the
field's mean rotation (`theta[:3]`) and one strain (6); the field's gradient (`theta[3:]`) is
fixed from the coverage field leg. Levenberg-Marquardt (Marquardt diagonal scaling, damping x10
on a rejected step, /3 on an accepted one), finite-difference Jacobian (steps 1e-5 rad and
1e-5, 9 extra predictions per iteration), per-spot weights, at most 15 iterations per stage, stop
when the step is below 1e-7. Weights: Tukey biweight with the cut-off annealed over three stages,
3, 1.5 and 1 px (amendment before the judged runs: a one-grain smoke test with the first
choice, Huber at 2 px, ended at a centroid cost below the ground truth's with 6e-4 strain error;
spots merged with a neighbour's still pull under Huber; Tukey at 1 px gave 1.5-2.4e-4 on grains
47, 111 and 237).

Benchmark S36, painted start, robust start field; output `benchmarks/grain_opt/go5i/`. Arms
(`--cost hyb --lam 8` so that `f_fit` and `f_gt` are hybrid costs):

- `GN@P`: coverage field leg 200 at zero strain, then Gauss-Newton.
- `GNH@P`: `GN@P`, then one joint CMA leg of 300 on `hyb` lambda 8 at steps 0.01 deg and 1e-4.

Judging rule (written before any run). Reference: GO5h `hyb` lambda 8 on the same grains
(strain 2.3e-4 / 1.5e-4 / 2.1e-4, orientation 0.011 / 0.009 / 0.017 deg, mean wall time per
grain). An arm "matches" when at 5e-4 and 1e-3 its strain median is at most 1.25 x the
reference's, its orientation median is at most 0.02 deg at every level, and its mean wall time
per grain is at most half the reference's (same machine load class; the load is reported). I44
(second part) is answered yes when an arm matches; the faster matching arm is preferred. 1e-4
strain is reported, not judged.

Budget: about 6 s field leg plus 3-5 s Gauss-Newton per grain (`GN@P`), plus about 18 s for
`GNH@P`; both arms under 20 minutes; 2 processes.

### GO5j: alternate field-gradient refit and Gauss-Newton (idea I45)

Question: GO5i's Gauss-Newton left the orientation field's gradient as the zero-strain coverage
leg fitted it, which is biased in proportion to the strain (0.039 deg at 1e-3). If the gradient
is refitted on coverage at the Gauss-Newton strain and mean rotation, then Gauss-Newton is run
again, does the fit match the CMA hybrid on strain and orientation at half its time or less?

Method: GO5i's `GN@P` (coverage field leg 200 at zero strain, Gauss-Newton with Tukey cut-off
3 / 1.5 / 1 px), then per round: CMA of 100 evaluations on coverage over the gradient only
(`theta[3:]`, 6 numbers, step 0.01 deg) with the mean rotation and strain held at the
Gauss-Newton values, then Gauss-Newton again from there (same annealing). Benchmark S36,
painted start, `--cost hyb --lam 8` for `f_fit`; output `benchmarks/grain_opt/go5j/`. Arms:

- `GA1@P`: one round.
- `GA2@P`: two rounds.

Judging rule (written before any run; the same as GO5i). Reference: GO5h `hyb` lambda 8 on the
same grains (strain 2.3e-4 / 1.5e-4 / 2.1e-4, orientation 0.011 / 0.009 / 0.017 deg, 53 s per
grain). An arm "matches" when at 5e-4 and 1e-3 its strain median is at most 1.25 x the
reference's, its orientation median is at most 0.02 deg at every level, and its mean wall time
per grain is at most half the reference's (26.7 s; load reported). I45 is answered yes when an
arm matches; the faster matching arm is preferred. 1e-4 strain is reported, not judged; grains
above 0.05 deg are named.

Budget: about 10 s (`GN@P`) plus 3 s coverage and 3-5 s Gauss-Newton per round, per grain:
both arms under 15 minutes; 2 processes.

### GO5k: stop the Gauss-Newton escape (idea I46)

Question: GO5j's Gauss-Newton failed at 1e-3 on grains with few spots, where it pushed
residuals past the Tukey cut-off (weights recomputed every iteration) and the weighted cost fell
to about 0. With the inlier set and weights frozen per stage, and a CMA-hybrid fallback for
grains with few spots, does the fast fit match the CMA hybrid at every level?

Method changes to `GA2@P` (two gradient-refit rounds):

- Frozen weights: in each Gauss-Newton stage (Tukey cut-off 3, 1.5, 1 px) the per-spot weights
  are computed once from the residuals at the stage start and kept for all its iterations. A
  stage whose inliers (weight > 0) are fewer than half the active spots is skipped (its start
  is kept).
- Fallback: a grain with fewer than 40 active spots at the start (centroid residuals defined)
  skips Gauss-Newton and runs GO5h's CMA hybrid protocol (strain only 300, joint 600 at 0.02 deg
  and 2e-4, `hyb` lambda 8) from the coverage field leg.

Benchmark S36, painted start, `--cost hyb --lam 8`; output `benchmarks/grain_opt/go5k/`. Arms:

- `GF2@P`: frozen weights, no fallback.
- `GFF@P`: frozen weights and fallback.

Judging rule (written before any run; the same as GO5i/j). Reference GO5h `hyb` lambda 8 on the
same grains (strain 2.3e-4 / 1.5e-4 / 2.1e-4, orientation 0.011 / 0.009 / 0.017 deg, 53 s per
grain). An arm "matches" when at 5e-4 and 1e-3 its strain median is at most 1.25 x the
reference's, its orientation median is at most 0.02 deg at every level, and its mean wall time
per grain is at most half the reference's (26.7 s; load reported). I46 is answered yes when an
arm matches; the faster matching arm is preferred. Fallback grains are counted and named; 1e-4
strain is reported, not judged.

Budget: about 24 s per Gauss-Newton grain, 53 s per fallback grain; both arms under 20 minutes;
2 processes.

### GO5l: fallback on stalls, polish large strains (idea I47)

Question: GO5k's `GFF@P` met the strain and time limits but left the orientation at 1e-3 at
0.028 deg (limit 0.02), and two grains stalled (every Gauss-Newton stage skipped, strain left
at 0). With a fallback on stalls and a short joint hybrid polish only where the fitted strain is
large, does the fast fit match the CMA hybrid at every level within half its time?

Method: `GFF@P`, plus:

- Stall fallback: if every Gauss-Newton stage of the first round is skipped, the grain runs
  GO5h's CMA hybrid protocol from the coverage field leg (as the spot-count fallback does).
- Polish: if the Gauss-Newton strain has tensor RMS above 6e-4, one joint CMA leg on `hyb`
  lambda 8 over field and strain (15 numbers) from the Gauss-Newton result, steps 0.02 deg and
  1e-4.

Benchmark S36, painted start, `--cost hyb --lam 8`; output `benchmarks/grain_opt/go5l/`. Arms:

- `GP4@P`: polish leg of 400 evaluations.
- `GP6@P`: polish leg of 600 evaluations.

Judging rule (written before any run; the same as GO5i-k). Reference GO5h `hyb` lambda 8 on the
same grains (strain 2.3e-4 / 1.5e-4 / 2.1e-4, orientation 0.011 / 0.009 / 0.017 deg, 53 s per
grain). An arm "matches" when at 5e-4 and 1e-3 its strain median is at most 1.25 x the
reference's, its orientation median is at most 0.02 deg at every level, and its mean wall time
per grain is at most half the reference's (26.7 s; load reported). I47 is answered yes when an
arm matches; the faster matching arm is preferred. Fallback and polished grains are counted and
named; 1e-4 strain is reported, not judged.

Budget: about 19 s per grain plus about 25-36 s per polished grain; both arms under 25
minutes; 2 processes.

### GO5m: strain fits on noisy images (idea I48)

Question: every strain result so far is on clean renders. Do the CMA hybrid (GO5h) and the fast
fit (`GP4@P`, GO5l) keep their accuracy on images with the project's detector noise model?

Data: the `strain` sample's existing `realistic_nodrop` render (`scripts/phase_d/noise.py`,
variant "nodrop", per connected spot of a frame: each lit pixel dropped and each unlit
4-neighbour lit with probability 0.05, one hot pixel with probability 0.05, a 2-4 px blob with
probability 0.1, all within 32 px of the spot; no spot drop). Built into a binary stack and an
intensity store in this worktree (`build_store.py --out_tag strainrn`; the render folder is
read, nothing is written into the GO5 worktree). Ground truth, painted start, grain membership
and spot windows as before; only the images change (`fit_strain.py --stack strainrn_clean.npy`).

Benchmark S36, painted start, `--cost hyb --lam 8`; output `benchmarks/grain_opt/go5m/`. Arms:

- `S@P` (CMA hybrid, GO5h protocol) on the noisy images.
- `GP4@P` (fast fit) on the noisy images.

Judging rule (written before any run). Clean references on the same grains: CMA hybrid
2.3e-4 / 1.5e-4 / 2.1e-4, orientation 0.011 / 0.009 / 0.017 deg; `GP4@P` 1.3e-4 / 1.3e-4 /
2.1e-4, orientation 0.009 / 0.015 / 0.020 deg. A fit is "robust to noise" when at 5e-4 and 1e-3
its noisy strain median is at most 1.5 x its clean median and its noisy orientation median is at
most 0.03 deg at every level. I48 is answered yes when both fits are robust; the fast fit "holds"
when it is robust and its noisy strain median is at most 1.25 x the noisy CMA hybrid's at 5e-4
and 1e-3. Fallback and polished grains are counted; 1e-4 strain is reported, not judged.

Budget: building the stack about 3 minutes; CMA hybrid about 32 minutes, fast fit about 15
minutes, in parallel; under 45 minutes.

### GO5n: noise-robust observed centroids (idea I49)

Question: on noisy images the fast fit's strain doubled (GO5m), most likely because the
intensity-weighted observed centroid moves with added edge pixels (which copy a neighbour's
intensity), hot pixels and blobs. Does a more robust observed centroid restore it?

Centroid variants (`fit_strain.py --cen_w`, applied to every centroid use: Gauss-Newton and the
hybrid's `cen` term):

- `bin`: binary weights (every kept lit pixel counts 1).
- `binclean`: binary weights after removing isolated lit pixels (fewer than 2 lit 8-neighbours
  among the spot's pixels in the same frame).

Benchmark S36, painted start, noisy images (`--stack strainrn_clean.npy`), `GP4@P`, `hyb`
lambda 8; output `benchmarks/grain_opt/go5n/`. Arms: `GP4@P --cen_w bin`, `GP4@P --cen_w
binclean`.

Judging rule (written before any run). Reference: GO5m's noisy CMA hybrid on the same grains
(strain 2.5e-4 / 1.9e-4 / 2.8e-4, 51 s per grain). An arm "holds" when at 5e-4 and 1e-3 its
strain median is at most 1.25 x the reference's, its orientation median is at most 0.03 deg at
every level, and its mean wall time per grain is at most half the reference's (25.5 s; load
reported). I49 is answered yes when an arm holds; the arm with the lower mean strain median over
5e-4 and 1e-3 is preferred. 1e-4 strain is reported, not judged; GO5m's noisy `GP4@P` with
intensity weights (1.8e-4 / 2.7e-4 / 3.7e-4) is the before.

Budget: about 27 s per grain, two arms in parallel, under 25 minutes.

### GO5o: intensity-weighted centroids without isolated pixels (idea I50)

Question: GO5n showed the intensity weights are what make the observed centroid sub-pixel. If
they are kept and only isolated lit pixels (fewer than 2 lit 8-neighbours among the spot's
pixels in the frame: hot pixels, stray flips) are removed, does the fast fit hold on noisy
images without losing its clean accuracy?

Centroid variant `intclean` (`fit_strain.py --cen_w intclean`): clean-intensity weights, isolated
lit pixels removed, in every centroid use. Benchmark S36, painted start, `GP4@P`, `hyb` lambda
8; output `benchmarks/grain_opt/go5o/`. Arms: `GP4@P --cen_w intclean` on the noisy images
(`--stack strainrn_clean.npy`) and on the clean images.

Judging rule (written before any run). Noisy: the arm "holds" when at 5e-4 and 1e-3 its strain
median is at most 1.25 x GO5m's noisy CMA hybrid (1.9e-4, 2.8e-4), its orientation median is at
most 0.03 deg at every level, and its mean wall time per grain is at most half the noisy CMA
hybrid's (25.5 s). Clean: no loss when at 5e-4 and 1e-3 its strain median is at most 1.25 x
GO5l's clean `GP4@P` (1.3e-4, 2.1e-4) and its orientation median at most 0.02 deg at every
level. I50 is answered yes when both hold. 1e-4 strain is reported, not judged.

Budget: about 27 s per grain, two arms in parallel, under 25 minutes.

### GO4d: cached ray index and the gated global search (goal G1b, idea I36)

Question: GO4c left 31% / 72% of `corr3` voxels and 77% / 83% of `rand` voxels more than 0.1
deg from the ground truth at 2 / 5 deg (`L/366/3`); those need a global search per voxel. The
reconstructor's ray-indexing seed path (index, strict cost on 5 candidates, finisher on 2) was
blocked on cost: its one-time pairing of all spots of the two detectors (about 1.56M fragments
on `rand`) took more than 15 minutes per process. With the paired rays cached per sample, does
the gated seed path recover the voxels the ladder leaves wrong?

Part A (unblock, infrastructure): `scripts/grain_opt/ray_cache.py` pairs the observed rays of a
sample once (`SmoothCost.observed_rays`, the reconstructor's own cached `SmoothCost`) and saves
them (`<sample>/observed_rays.npz`, git-ignored); later processes load them into the same
object, so `index_seed` skips the pairing. Check: the loaded arrays equal the built ones; the
pairing time and the load time are reported.

Part B (judged): `run_pipeline.py --start centre --group true` on `corr3` and `rand` (GO4
benchmark, all levels 0.5 / 1 / 2 / 5 deg), output `benchmarks/grain_opt/go4d/`. Arms:

- `L/366/3,G/15/0.15`: the ladder, then the seed path on voxels whose strict cost is above 0.15.
- `L/366/3,G/10/0.10`: the same with threshold 0.10.

(GO4b: at the ground truth the median strict cost is 0.065-0.096, at a wrong end 0.10-0.71.)

Judging rule (written before any run). Reference: GO4c `L/366/3` on the same sample. Per level,
share of voxels more than 0.1 deg from the ground truth (Wilson 95% intervals) and share wrong
(above 1 deg). I36 is answered yes for an arm when on both samples the share above 0.1 deg at
5 deg falls by at least a third against the reference (corr3 72% -> at most 48%, rand 83% -> at
most 55%) and no level gets worse by more than 2 points. The cheaper passing arm (ms per voxel)
is preferred. Seed-path calls, improved voxels and wall time per voxel are reported.

Budget: Part A one pairing per sample (15-30 minutes each, in parallel); Part B the gated voxels
times the seed path (to be measured on a smoke grain before the runs); each run under 2 hours.

### GO7: rocking curve on the images (goal G2b, idea I23, part 1)

Question: the renderer puts every spot in exactly one omega frame; real spots spread over
frames, by Delta omega = alpha / |sin eta| (`docs/omega_peak_width_derivation.md`, eta the
azimuth from the rotation axis, alpha the Bragg window). How much does the orientation fit lose
when the images have a rocking curve but the forward model does not?

Data (emulated, no library change): `scripts/grain_opt/rock_stack.py` reads the clean `grad`
intensity store and spreads each lit pixel's intensity over the frames that the window
[omega_c - Delta omega / 2, omega_c + Delta omega / 2] overlaps (omega_c the frame centre: the
render does not keep the sub-frame omega), with Delta omega from the pixel's eta (lab direction
from the sample origin, `spot_pairing.pixel_to_lab`), capped at 8 frames, and thresholds the
spread intensity at the 5th percentile of the clean intensities (a fixed detector threshold:
weak spots spread thin fall below it). Stacks `stacks/grad_rockA_clean.npy` for alpha = 0.3 and
1.0 deg (frame width 1 deg); lit-pixel counts and spots lost are reported.

Arms (`run_pipeline.py`, benchmark G20, painted start, `--stack`): `Astrict/100` (the GO2
frontier arm) on clean, alpha 0.3 and alpha 1.0 images. Output `benchmarks/grain_opt/go7/`.

Judging rule (written before any run). Per arm, voxels more than 0.1 deg from the ground truth
(Wilson 95%) and median error, on the 853 voxels with a C record. The rocking curve "harms" the
fit at an alpha when the share above 0.1 deg has a Wilson interval entirely above the clean
arm's, or the median is more than 1.5 x the clean median. I23 part 1 is answered by the harm per
alpha. A forward model with a rocking curve (part 2) changes `icenine/` and needs the owner's
approval before it is planned.

Budget: building a stack about 2 minutes; each arm about 2 minutes (GO2: 5.2 ms per voxel at
200 evaluations); 3 processes.

### GO7b: detector point-spread function on the images (goal G2b, idea I24)

Question: real detectors blur each spot (point-spread function) before the threshold; the
forward model has no PSF. How much do the orientation fit and the fast strain fit lose?

Data (emulated, no library change): `scripts/grain_opt/psf_stack.py` blurs each frame's clean
intensities with a Gaussian PSF (sigma 0.5 or 1.0 px, `scipy.ndimage.gaussian_filter`) and
thresholds at the 5th percentile of the clean intensities; stack and intensity store written in
this worktree (`grad_psf05`, `grad_psf10`, `strain_psf05`, `strain_psf10`).

Arms: `run_pipeline.py --pipeline Astrict/100` (G20, painted start) on `grad_psf05` and
`grad_psf10` (clean reference: GO7's `go7_clean_As100`); `fit_strain.py --arm GP4@P --cost hyb
--lam 8 --per_level 12` on `strain_psf05` and `strain_psf10` (clean reference: GO5l `GP4@P`,
1.3e-4 / 1.3e-4 / 2.1e-4). Output `benchmarks/grain_opt/go7b/`.

Judging rule (written before any run). Orientation: as GO7 (harm at a sigma when the share of
the 853 C voxels above 0.1 deg has a Wilson interval entirely above the clean arm's, or the
median is above 1.5 x the clean median). Strain: harm at a sigma when the strain median at 5e-4
or 1e-3 is above 1.5 x the clean median, or the orientation median is above 0.03 deg at a level.
I24 (images part) is answered by the harm per sigma and fit; a PSF in the forward model is
library work for the owner to approve.

Budget: stacks about 2 minutes each; orientation arms 2 minutes; strain arms about 15 minutes;
4 processes at most.

### GO7c: PSF-aware predicted coverage in the strain costs (goal G2b, idea I52)

Question: with a 1 px PSF on the images the fast strain fit breaks (GO7b: 2.3 / 2.7 / 5.7e-4,
orientation up to 0.06 deg), most likely because the observed spots grow by about a pixel while
the predicted coverage stays sharp, flattening the coverage term. If the predicted coverage is
blurred by the same PSF before it is scored, does the fit recover?

Method (`fit_strain.py --psf_sigma S`, scripts only): in the coverage cost every predicted
coverage point is replaced by K = 4 copies displaced in column and row by Gaussian offsets of
sigma S px (fixed seed, so the cost is deterministic), each with a quarter of the mass. The
centroid cost is unchanged (a symmetric blur keeps the centroid). This affects the coverage leg,
the gradient refits and the hybrid polish and fallback of `GP4@P`.

Arm: `GP4@P --cost hyb --lam 8 --per_level 12 --stack strain_psf10_clean.npy --psf_sigma 1.0`;
output `benchmarks/grain_opt/go7c/`.

Judging rule (written before any run). Reference: GO5l clean `GP4@P` (1.3e-4 / 1.3e-4 /
2.1e-4) and GO7b `GP4@P` on the same 1 px images without the PSF model (2.3e-4 / 2.7e-4 /
5.7e-4). I52 is answered yes when the PSF-aware arm has no harm by GO7b's rule (strain median at
5e-4 and 1e-3 at most 1.5 x clean, orientation median at most 0.03 deg at every level), and
partly when it improves both strain medians at 5e-4 and 1e-3 by at least 1.5 x against GO7b
without meeting the no-harm limits. 1e-4 is reported, not judged; wall time per grain reported.

Budget: about 30-40 s per grain (four times the coverage points), under 25 minutes.

### GO9: Poisson likelihood on pixel counts as an orientation objective (goal G2a, idea I18)

Question: GO1 found every smooth objective biased (optimum 0.012-0.095 deg from the ground truth)
and only the strict cost unbiased (0.006 / 0.008 deg). Does a Poisson likelihood on the observed
pixel counts, with a per-reflection scale and a background, have its optimum at the ground
truth?

Objective (`scripts/grain_opt/probe_poisson.py`, no library change): per spot window, predicted
counts mu = s_j m + b over the union of observed lit pixels and pixels with predicted coverage
mass m (grain plus fixed neighbour shell, as GO1's objectives), observed counts y = clean
intensity / gain (gain = 5th-percentile intensity / 50, so the weakest lit pixels hold about 50
counts); s_j the closed-form scale sum(y) / sum(m); b = 0.5 counts per pixel; objective = Poisson
deviance summed over spots. Variants: `pois` as above, `pois_b5` with b = 5.

Protocol and rule: GO1's (`grad`, level-0 grains 381 and 302, CMA over the grain rotation only,
sigma 0.03 deg, 150 evaluations, started at the ground truth). An objective qualifies if its
optimum is below 0.01 deg from the ground truth on both grains (written before any run). GO1's
numbers are the reference (strict 0.0060 / 0.0082; least squares 0.0117 / 0.0423;
intensity-weighted Sinkhorn 0.0074 / 0.0263). Output `benchmarks/grain_opt/go9/`.

Budget: a few minutes.

### GO9b: Poisson likelihood with structure-factor scales (goal G2a, idea I22)

Question: GO9's Poisson deviance missed the 0.01 deg rule by 0.0001 deg on grain 381, with a
free scale per reflection. The renderer gives each triangle the intensity
I = F / (|sin eta| sin 2 theta) (F the reflection's form intensity, `peak_filters` /
`forward_simulation`), so the scales are known up to one constant. With s_j = c I_j (one c per
grain, closed form c = sum(y) / sum(I_j m)), does the optimum move to the ground truth?

Objective (`probe_poisson.py --variant pois_sf`): GO9's Poisson deviance (b = 0.5) with s_j =
c I_j, where F comes from the structure's reflection list (filtered by Q-max as `SmoothCost`
does, indexed by the spot's reflection), and eta and 2 theta from the lab direction of the
spot's predicted centroid (sample origin, `spot_pairing.pixel_to_lab`; eta from the rotation
axis around the beam). Arms: `pois_sf`, and `pois` rerun as the control.

Protocol and rule: GO1's (grains 381 and 302, CMA over the rotation, sigma 0.03 deg, 150
evaluations, from the ground truth); qualifies if the optimum is below 0.01 deg from the ground
truth on both grains (written before any run). Output `benchmarks/grain_opt/go9b/`.

Budget: a few minutes.

### GO4e: neighbour pass after the gated seed path (goal G1b, idea I51)

Question: GO4d's ladder + gated seed path (`L/366/3,G/15/0.15`) left rand at 2 deg with 52% of
voxels above 0.1 deg (and 44% at 5 deg), many near but outside the basin with a strict cost below
the gate, others seeded from a neighbouring grain. Do the voxels the seed path fixed help their
neighbours through a neighbour continuation (candidates = own and fitted in-grain neighbours'
orientations, then per-voxel CMA), and does a second seed pass afterwards pick up more?

Arms (`run_pipeline.py --start centre --group true`, corr3 and rand, output
`benchmarks/grain_opt/go4e/`):

- `L/366/3,G/15/0.15,N/30/0.3`.
- `L/366/3,G/15/0.15,N/30/0.3,G/15/0.15`.

Judging rule (written before any run). Reference: GO4d `L/366/3,G/15/0.15` on the same sample.
I51 is answered yes for an arm when on rand the share above 0.1 deg at 2 deg falls by at least a
third (52% -> at most 35%) and no level on either sample gets worse by more than 2 points; the
cheaper passing arm is preferred. Shares, wrong (above 1 deg) and ms per voxel reported.

Budget: GO4d took 5-12 minutes per run; N/30 adds about 0.01 s per voxel and a second seed pass
up to the first's cost; 4 runs in parallel, under 30 minutes.

### GO4f: why voxels stay wrong after the G1b pipeline (goal G1b, diagnosis)

Diagnosis only, no arms and no judging rule (this section was written at the step's close; the
question was set in STATE.md before the run): for the voxels still more than 1 deg from the
ground truth after `L/366/3,G/15/0.15,N/30/0.3` (GO4e), were they gated by the seed path, did it
lower their strict cost, and is the strict cost at the ground truth lower than at their end
(search failure) or not (objective failure)? `scripts/grain_opt/go4f_diag.py`, output
`benchmarks/grain_opt/go4f/`.

### GO4g: wider search for the voxels the seed path misses (goal G1b, idea I54)

Question: GO4f showed every remaining wrong voxel is a search failure (the strict cost at the
ground truth is far lower than at the end; the seed path's candidates miss it). Does a wider
search on the gated voxels still above the gate recover them?

Stages added to `run_pipeline.py` (each runs only on voxels whose strict cost is still above the
threshold and keeps a result only if it lowers the strict cost):

- `W/<n>/<threshold>`: the seed path with every indexing candidate costed (`n_cost` 20 = the
  indexer's top 20) and the finisher on the 5 lowest (`n_refine` 5; GO4d used 5 and 2).
- `X/<n>/<threshold>`: the reconstructor's per-voxel grid search (`reconstruct_voxel`, the
  painting fallback, `grid_seed`).

Arms (corr3 and rand, `--start centre --group true`, output `benchmarks/grain_opt/go4g/`):

- `L/366/3,G/15/0.15,N/30/0.3,W/20/0.15`.
- `L/366/3,G/15/0.15,N/30/0.3,W/20/0.15,X/1/0.15`.

Judging rule (written before any run). Reference: GO4e `L/366/3,G/15/0.15,N/30/0.3` on the same
sample. I54 is answered yes for an arm when the share above 0.1 deg falls by at least a third at
2 deg on rand (48% -> at most 32%) and at 5 deg on both samples (corr3 33% -> at most 22%, rand
41% -> at most 27%), with no level worse by more than 2 points; the cheaper passing arm is
preferred. Voxels searched and improved per stage, and wall time per voxel, are reported. A run
over 2 hours is stopped and reported as over budget.

Budget: `W` about 2-3 times the seed path per searched voxel (about 3 s); `X` unknown (measured
on a smoke grain first; if above 10 s per voxel, `X` runs on rand only).

### GO4h: grid search on corr3 and a final per-voxel polish (goal G1b, idea I55)

Question: GO4g's full search (`W,X`) left rand at 0.5 deg at 17% above 0.1 deg, nearly all between
0.1 and 1 deg (a precision problem), and `X` was not run on corr3. Does a final per-voxel strict
polish fix the near misses, and what does `X` do on corr3?

Arm, both samples (`run_pipeline.py --start centre --group true`, output
`benchmarks/grain_opt/go4h/`): `L/366/3,G/15/0.15,N/30/0.3,W/20/0.15,X/1/0.15,V/61/0.02` (the
per-stage rows give the result after `X` and after the polish).

Judging rule (written before any run). Reference: GO4g `+W,X` on rand (17 / 2 / 7 / 6% above 0.1
deg at 0.5 / 1 / 2 / 5 deg) and `+W` on corr3 (2 / 3 / 3 / 21%). I55 is answered yes when on
rand the share above 0.1 deg at 0.5 deg after the polish falls by at least a third (17% -> at most
11%) with no level worse by more than 2 points on either sample. corr3 after `X` is reported
against `+W`. Wall time per voxel reported; a run over 2 hours is stopped and reported.

Budget: rand about 75 minutes (GO4g `W,X` 70 minutes; the polish about 0.1 s per voxel); corr3
about 25 minutes (`X` on the voxels left above the gate after `W`, about 80 x 16 s).

### GO7d: split merged spots before taking centroids (goal G2b, idea I53)

Question: on 1 px PSF images the fast strain fit degrades even on its Gauss-Newton-only grains
(GO7c), most likely because the blur merges nearby spots into one connected component and the
observed centroid mixes a neighbour's light. If each observed lit pixel in a spot window is
assigned to the grain or to a neighbour by the nearest predicted point, does the fit recover?

Method (`fit_strain.py --cen_split`, scripts only): for each spot, the grain's predicted coverage
points at the start orientation (no strain) and the fixed shell (neighbour) prediction stored
with the region; an observed lit pixel enters the observed centroid only if its nearest grain
point is no farther than its nearest shell point (column and row in pixels, frame in frames).
Spots without shell points in the window are unchanged.

Arm: `GP4@P --cost hyb --lam 8 --per_level 12 --stack strain_psf10_clean.npy --cen_split`;
control `GP4@P --cen_split` on clean images. Output `benchmarks/grain_opt/go7d/`.

Judging rule (written before any run). References: GO5l clean `GP4@P` (1.3e-4 / 1.3e-4 /
2.1e-4, orientation up to 0.020 deg) and GO7b `GP4@P` on the 1 px images (2.3e-4 / 2.7e-4 /
5.7e-4). I53 is answered yes when the PSF arm has no harm by GO7b's rule (strain median at 5e-4
and 1e-3 at most 1.5 x clean, orientation median at most 0.03 deg at every level), and partly when
it improves both strain medians at 5e-4 and 1e-3 by at least 1.5 x against GO7b. The clean
control must not lose more than 1.25 x at 5e-4 and 1e-3. 1e-4 reported, not judged.

Budget: about 20-25 s per grain; two runs in parallel, under 25 minutes.

### GO7e: widen the region area-ratio filter for blurred images (goal G2b, idea I56)

Question: GO7d found that the region build's area-ratio filter (observed lit pixels / pixels
touched by the start prediction, kept in 0.8-1.6) drops the blur-enlarged spots of the 1 px PSF
images (a grain keeps 58% of its spots). With the upper bound raised, does the fast strain fit
recover on those images?

Method: `fit_strain.py --ratio_hi H` passes `region_kw=dict(ratio=(0.8, H))` to `GrainProblem`
(scripts only). Arms: `GP4@P --cost hyb --lam 8 --per_level 12 --ratio_hi 2.5` on
`strain_psf10_clean.npy`, and the same on clean images (control). Output
`benchmarks/grain_opt/go7e/`. Spot counts per grain are reported.

Judging rule (written before any run; GO7b's). References: GO5l clean `GP4@P` (1.3e-4 / 1.3e-4 /
2.1e-4) and GO7b on the 1 px images (2.3e-4 / 2.7e-4 / 5.7e-4). I56 is answered yes when the PSF
arm has no harm (strain median at 5e-4 and 1e-3 at most 1.5 x clean, orientation median at most
0.03 deg at every level), partly when both strain medians at 5e-4 and 1e-3 improve by at least
1.5 x against GO7b. The clean control must stay within 1.25 x of GO5l at 5e-4 and 1e-3 (a wider
filter also admits spots overlapping other grains). 1e-4 reported, not judged.

Budget: about 20-30 s per grain (more spots); two runs in parallel, under 30 minutes.

### GO5p: the wider ratio filter on noisy images (idea I57)

Question: GO7e found that widening the region area-ratio filter to 2.5 improves the clean fast
strain fit (1.3 / 1.2 / 1.4e-4). Does it hold on the noisy images (GO5m's `strainrn`), so that it
can become the fast fit's default?

Arms (S36, painted start, `GP4@P --cost hyb --lam 8 --stack strainrn_clean.npy --ratio_hi 2.5`,
output `benchmarks/grain_opt/go5p/`): plain; and with `--cen_w intclean` (GO5o).

Judging rule (written before any run). References on the noisy images: GO5m fast fit (1.8e-4 /
2.7e-4 / 3.7e-4) and noisy CMA hybrid (2.5e-4 / 1.9e-4 / 2.8e-4, 51 s per grain). The wider filter
becomes the default when an arm's noisy strain medians at 5e-4 and 1e-3 are at most GO5m's fast
fit and at most 1.25 x the noisy CMA hybrid's (2.4e-4, 3.5e-4), with orientation medians at most
0.03 deg at every level; the faster passing arm is preferred. 1e-4 reported, not judged.

Budget: about 35 s per grain; two runs in parallel, under 30 minutes.

### GO4i: a lower gate for rand's near misses (goal G1b, idea I58)

Question: after GO4h, rand at 0.5 deg keeps 16% of voxels above 0.1 deg, nearly all between 0.1
and 1 deg, with strict costs below the 0.15 gate (median 0.133, 10th percentile 0.110; the
voxels within 0.1 deg have a 90th percentile of 0.110). Does a second wide-search pass with the
gate at 0.11 recover them without disturbing the good voxels?

Arm (rand, `run_pipeline.py --start centre --group true`, output `benchmarks/grain_opt/go4i/`):
`L/366/3,G/15/0.15,N/30/0.3,W/20/0.15,X/1/0.15,W/20/0.11,X/1/0.11,V/61/0.02` (GO4h's pipeline with
a second `W,X` pass at gate 0.11 before the polish). corr3 is not run (0.5 deg already at 1%).

Judging rule (written before any run). Reference: GO4h rand (16 / 1.4 / 5 / 4% above 0.1 deg at
0.5 / 1 / 2 / 5 deg). I58 is answered yes when rand at 0.5 deg falls by at least a third (to at
most 11%) and no level gets worse by more than 2 points. Searched and improved voxels per stage
and wall time per voxel reported; stopped and reported if over 2 hours.

Budget: GO4h rand took 71 minutes; the second pass adds about 60 voxels x (1 + 16) s, about
17 minutes.

### GO5q: intragranular strain gradient (goal G1c, idea I59; owner request 2026-11-19)

Question: the `strain` sample has one strain tensor per grain. With a mean strain plus a linear
strain gradient across each grain, do the strain fits recover the mean and the gradient?

Sample `straing` (`scripts/grain_opt/make_strain_gradient_sample.py`, scripts only; the renderer's
opt-in per-voxel strain path, `.mic` columns 13-18): the `grad` orientations; per grain the same
mean strain as `strain` (same seed, same levels 1e-4 / 5e-4 / 1e-3) plus an in-plane linear
gradient eps(v) = eps_mean + D d_v, d_v = (xy_v - c) / r (c the grain centroid, r the RMS voxel
distance from it), D a 6 x 2 Gaussian matrix scaled so the gradient part has tensor RMS over the
grain's voxels equal to the grain's level (so the strain varies across a grain about as much as
its mean). The layer is one voxel thick, so only the in-plane gradient (6 x 2 = 12 numbers) is
observable, not 6 x 3. No voxel noise in this step. Rendered clean (`render_full.py --strain`),
stack and store built in this worktree; per-voxel ground-truth strain in `strain_vox.npy`.

Fits (`fit_strain.py --grad`: per-grain strain = 6 mean + 12 gradient numbers; mean error = tensor
RMS of the fitted minus ground-truth grain-mean strain; gradient error = RMS over the grain's
voxels of the tensor RMS of the fitted minus ground-truth gradient part, reported also relative to
the ground-truth gradient RMS). S36, output `benchmarks/grain_opt/go5q/`. Arms:

- `E@GT` mean only (6) and with gradient (18): strain only at the ground-truth orientation,
  `hyb` lambda 8, 300 / 900 evaluations.
- `S@P --grad`: the CMA hybrid staged from the painted start (field 200, strain 18 numbers 900,
  joint 900).
- `GP4@P --grad --ratio_hi 2.5`: the fast fit for the mean (Gauss-Newton on centroids does not see
  a linear gradient to first order), then a strain leg over the 18 numbers on `hyb` lambda 8 (600
  evaluations, field fixed).

Judging rule (written before any run). The gradient counts as recovered by an arm at a level when
the median relative gradient error is at most 0.5 (at least half the gradient explained). I59 is
answered yes for an arm that recovers the gradient at 5e-4 and 1e-3 and whose mean-strain median at
those levels is at most 1.5 x the same method's uniform-strain result on `strain` (CMA hybrid 1.5e-4
/ 2.1e-4; fast fit with filter 2.5 1.2e-4 / 1.4e-4). 1e-4 reported, not judged; orientation reported.

Budget: sample and render about 5 minutes; `E@GT` arms about 10-25 minutes; `S@P --grad` about 70
minutes; fast fit about 40 minutes; at most 3 processes, each under 2 hours.

### GO5r: Gauss-Newton gradient fit with spot-spread residuals (goal G1c, idea I60)

Question: GO5q found the linear strain gradient only partly recovered (relative error 0.58-0.97
from the painted start with the CMA hybrid), because a linear gradient changes each spot's shape,
not its centroid. Does Gauss-Newton on per-spot centroids plus second moments (the spread)
recover it better, and faster?

Method (`fit_strain.py --grad --gnm`, scripts only): residuals per spot = predicted minus observed
intensity-weighted centroid (column, row, frame) and second moments (variance in column, in row,
and the column-row covariance, in px^2; the predicted moments add the pixel box variance 1/12 px^2
per axis); amended before the judged runs: the column variance is dropped (at the ground
truth on a smoke grain its residual has median -4.2 px^2, larger than the gradient's effect of
2.1 px^2, likely the spot's smear along the column during a 1 deg frame, which the single-omega
prediction lacks; the row variance and covariance biases are 0.07 and 0.06 px^2 against
effects of 0.35 and 0.98 px^2); Levenberg-Marquardt with a finite-difference Jacobian, frozen Tukey weights per stage
on the centroid distance (3 / 1.5 / 1 px, as GO5k). Arms on `straing`, S36, ratio filter 2.5,
`hyb` lambda 8 for `f_fit`; output `benchmarks/grain_opt/go5r/`:

- `E@GT --grad --gnm`: strain only (18 numbers) at the ground-truth orientation.
- `GP4@P --grad --gnm`: the fast fit for the field and mean (GO5l), then Gauss-Newton over the
  mean rotation and the 18 strain numbers on centroid + moment residuals (replacing GO5q's CMA
  gradient leg). Amendment before the judged runs: when every Gauss-Newton stage is skipped
  (fewer than half the spots within the Tukey cut-off; on a smoke grain this came from painted
  outlier voxels putting spurious spots into the region), GO5q's CMA gradient leg runs instead.

Judging rule (written before any run). References from GO5q on the same grains (relative gradient
error at 5e-4 / 1e-3): strain only at the ground-truth orientation 0.77 / 0.45; CMA hybrid from the
painted start 0.97 / 0.58 (123 s per grain); fast fit + CMA gradient leg 0.79 / 0.66 (76 s). I60 is
answered yes when the painted-start arm reaches relative gradient error at most 0.5 at 5e-4 and
1e-3 with mean-strain medians within GO5q's limits (2.2e-4, 3.2e-4); partly when it beats the CMA
hybrid's gradient error by at least 1.5 x at both levels. Time per grain and orientation reported;
the 1e-4 level reported, not judged.

Budget: about 30-60 s per grain (21-parameter Jacobians); two runs in parallel, under 60 minutes.

### GO5s: non-linear strain-gradient sample and the linear fit on it (goal G1c, idea I61)

Question: a smooth per-grain strain field richer than linear is the case where a quadratic field
(I62) or a prior (I64) could pay off. Build that sample, and measure what the linear (6 + 12) fits
recover on it, as the baseline for I62 and I64.

Sample `straingq` (`make_strain_gradient_sample.py --mode nonlinear`, scripts only): the `strain`
grain means (same seed and draws), plus per grain three intragranular parts of equal RMS: an
in-plane linear term, a quadratic term (d_x^2, d_x d_y, d_y^2, each centred over the grain) and a
smooth Gaussian random field per strain component (squared-exponential kernel, correlation length
3 voxel spacings), each scaled to tensor RMS level / sqrt(3) over the grain's voxels (amended
before any judged run: the quadratic and random parts are made orthogonal to the linear basis per
grain, since without that the oracle linear share was 0.60, above the rule's 0.5), so the total
intragranular part has RMS about the grain's level (as `straing`). Rendered clean and stored in
this worktree; the per-voxel ground truth in `strain_vox.npy`. The linear share of the
intragranular field (the least-squares linear projection of the ground truth, per grain) is
reported as the best a linear fit could do.

Arms (S36, `--sample straingq --grad --ratio_hi 2.5`, `hyb` lambda 8; output
`benchmarks/grain_opt/go5s/`): strain only with the linear field at the ground-truth orientation
(CMA, 900), and the fast fit + CMA gradient leg from the painted start (GO5q's arms).

Judging rule (written before any run). I61 is answered (sample delivered) when the sample renders
and the oracle linear projection explains at most 50% of the intragranular variance (median over
grains; so a non-linear field has room to gain). The linear fits are judged against their `straing`
results (GO5q): they hold when their mean-strain medians at 5e-4 and 1e-3 are within 1.5 x of
GO5q's; their relative error against the full intragranular field is reported next to the oracle
linear floor.

Budget: sample and render about 5 minutes; the two arms about 40 and 50 minutes in parallel.

### GO5t: quadratic per-grain strain field (goal G1c, idea I62)

Question: on `straingq` (linear + quadratic + smooth random intragranular strain), does a quadratic
per-grain strain field (mean 6 + linear 12 + quadratic 18 = 36 numbers; centred quadratic terms
d_x^2, d_x d_y, d_y^2) recover more of the intragranular field than the linear field (GO5s)?

Method: `fit_strain.py --sample straingq --grad --quad` (scripts only), `hyb` lambda 8, ratio
filter 2.5, S36; output `benchmarks/grain_opt/go5t/`. Arms: strain only at the ground-truth
orientation (CMA over 36 numbers, 1,500 evaluations); fast fit + quadratic CMA leg from the painted
start (1,000 evaluations).

Judging rule (written before any run). References: GO5s's linear arms on the same grains (relative
error against the full intragranular field 1.05 / 0.98 at 5e-4 / 1e-3 at the ground-truth
orientation; 1.10 / 0.98 from the painted start); the oracle floor for linear + quadratic is about
0.58. I62 is answered yes for an arm whose relative error is at least 0.1 lower than its linear
counterpart at 5e-4 and at 1e-3, with mean-strain medians within 1.5 x of the linear arm's. I63 (a
coordinate network) is run only if the quadratic arm is yes and still at least 0.2 above the 0.58
floor at 1e-3. The 1e-4 level is reported, not judged.

Budget: about 100 s and 110 s per grain; two runs in parallel, about 65 minutes.

### GO5u: GP prior on the intragranular strain residual (goal G1c, idea I64)

Question: GO5t's quadratic field helped at 1e-3 but fitted noise at 5e-4: the right amount of
intragranular structure depends on the strain level. Does a GP prior on a smooth residual, on top
of the mean and linear field, let the data decide, and do better than the quadratic at 5e-4 without
losing at 1e-3?

Method (`fit_strain.py --gpk 5 --gp_beta B`, scripts only): per grain, the residual is
sum_k c_k phi_k(v) per strain component, phi_k the K = 5 leading eigenvectors of a
squared-exponential kernel over the grain's voxels (correlation length 3 voxel spacings), made
orthogonal to the mean and linear basis and scaled to RMS 1, with eigenvalues lambda_k normalised to
lambda_1 = 1; strain numbers per grain 6 + 12 + 30 = 48. The prior adds
B x mean_k,j (c_jk^2 / (lambda_k sigma^2)) to the cost, sigma = 5e-4 (a GP prior in the reduced
basis; B = 0 is the unregularised 48-number fit). On `straingq`, S36, `hyb` lambda 8, ratio filter
2.5; output `benchmarks/grain_opt/go5u/`. Arms: strain only at the ground-truth orientation (CMA,
1,500 evaluations) with B = 0 and B = 0.05; the fast fit + GP leg from the painted start (1,000
evaluations) with B = 0.05.

Judging rule (written before any run). References: GO5t's quadratic arms (relative intragranular
error 1.24 / 0.72 at 5e-4 / 1e-3 at the ground-truth orientation; 1.09 / 0.74 from the painted
start). I64 is answered yes for a prior arm (B > 0) whose relative error at 5e-4 is at least 0.1
lower than its quadratic counterpart's and at 1e-3 at most 0.05 higher, with mean-strain medians
within 1.5 x of GO5s's linear arms. The B = 0 arm shows what the prior itself contributes. The
1e-4 level is reported, not judged.

Budget: about 100-110 s per grain; three runs in parallel, about 70 minutes.

### GO5v: build the region from the outlier-cleaned start (idea I65)

Question: GO5r found that painted outlier voxels put spurious spots into a grain's region (the spot
windows are built from the painted map; grain 237: centroid residuals of 46 px even at the ground
truth). In GO7e's clean fast fit the 9 stalled grains, the 9 CMA fallbacks and the 9 grains with
painted outliers are the same count. If the region is built from the start map with each grain's
outlier voxels set to its medoid orientation (as `robust_theta0` does for the start field), do the
stalls disappear and does the fit improve?

Method (`fit_strain.py --clean_start`, scripts only): before `GrainProblem`, the members' painted
orientations farther than 1.5 deg from the medoid voxel are replaced by the medoid's in the start
map passed to it (only the grain's own members; the shell is unchanged). Arm: `GP4@P --cost hyb
--lam 8 --per_level 12 --ratio_hi 2.5 --clean_start` on clean `strain`; output
`benchmarks/grain_opt/go5v/`.

Judging rule (written before any run). Reference: GO7e clean `GP4@P --ratio_hi 2.5` (strain 1.3e-4 /
1.2e-4 / 1.4e-4, orientation at most 0.013 deg, 9 stalls and 9 fallbacks, 35 s per grain). I65 is
answered yes when the stalls fall to at most 3 and the strain medians at 5e-4 and 1e-3 are at most
1.1 x the reference's with orientation medians at most 0.02 deg. Wall time per grain is reported.

Budget: about 35 s per grain; one run, about 25 minutes.

### GO8a: retrain GNLayerNet (goal G2c, idea I25, part 1; owner approval D-20261118-I1)

Question: no GNLayerNet weights or datasets exist on this machine (scoped 2026-11-19). Before using
the network inside grains (GO8b), can the documented realism-trained network be retrained here to
its documented accuracy?

Method: the recipe in `docs/orientation_nn_design.md` section 9 (scripts only): datasets with
`generate_toy_orientation_dataset.py` (ManyGrains, 30 voxels by radius, 24 trained and 6 held out,
500 offsets per voxel within 1 deg; the distractor layer `--neighbors 2 --twin`, tag `arch_dis`; a
`stage3_multi` test file only for the aux sidecar), `make_dataset_aux.py`, then
`train_toy_orientation_nn.py --arch gn --gn-iters 3 --loss decoupled --device mps --lr 1e-4
--clip 1 --cosine --batch-size 64 --epochs 60 --ema 0.998 --checkpoint ema --val-voxels 4
--realistic-train all --seed 0 --save-model`. Datasets and weights under `scripts/` (git-ignored);
results JSON in `benchmarks/grain_opt/go8/`. Timing: a 720-sample smoke dataset took 10 s.

Judging rule (written before any run). Reference (design doc table, realism-trained net, held-out
median error): 0.0193 deg on clean windows and 0.0821 deg on `all` windows (2-seed means; seed 0
0.0176 and 0.0827). GO8a passes (the retrained network reproduces) when its held-out median error
is at most 1.25 x the seed-0 reference on both variants (0.022 and 0.103 deg). If it fails, GO8b is
not started and the gap is reported.

Budget: datasets about 5 minutes and 7 GB of disk; training about 15-20 minutes on the GPU.

### GO8b: GNLayerNet inside grains (goal G2c, idea I25, part 2)

Question: GO8a retrained the realism-trained GNLayerNet (a learned per-voxel Gauss-Newton refiner,
offsets up to 1 deg, trained on simulated windows with neighbour overlap and detector noise). On the
observed `rand` / `corr3` images, does it refine the voxels a local search leaves off, after the
step ladder?

Method (scripts only): a new `run_pipeline.py` stage `Q/<iterations>/<threshold>`: for each voxel
of a grain whose strict cost is above the threshold (0 = every voxel), up to `iterations` rounds of
(nominal = the voxel's current orientation; peak set, window placement, context and nominal offsets
from `perturbation_sweep.prepare_nominal` on the ManyGrains geometry, whose voxel grid equals the
samples' grid; observed windows cut from the sample's binary stack, 32 x 32 px at each peak's
detector, frames frame0 +- 4, each lit pixel coded by its frame offset, the lit frame nearest frame0
when several are lit, as the training windows; the network's offset applied as exp([delta]x) R_nom),
keeping a round only if it lowers the voxel's strict cost. Weights
`scripts/toy_orientation_sweep_model_realistic_s0.pt`. Benchmark GO4 (`--start centre --group
true`, 1,043 voxels per sample, levels 0.5 / 1 / 2 / 5 deg); output `benchmarks/grain_opt/go8/`.
Arms on rand and corr3: `L/366/3,Q/3/0` (network after the ladder, every voxel) and `Q/3/0` (network
from the grain-centre start, to see its reach).

Judging rule (written before any run). Reference: GO4c `L/366/3` on the same sample (rand 48 / 53 /
77 / 83%, corr3 3 / 9 / 31 / 72% above 0.1 deg). I25 is answered yes when `L/366/3,Q/3/0` lowers the
share above 0.1 deg on rand by at least a third at 0.5 and at 1 deg (to at most 32% and 35%) with no
level worse by more than 2 points on either sample; partly when it lowers rand at 0.5 and 1 deg by
at least 5 points each. The `Q/3/0` arm and the GO4i pipeline (rand 4.5 / 0.4 / 3.3 / 3.2%) are
reported for context. Wall time per voxel reported; a run over 2 hours is stopped and reported.

Budget: preparing a voxel nominal (peak set, observer) and one network call per round; about
0.3-1 s per voxel per round (measured on a smoke grain first); 4 runs, each under 2 hours.

### GO2b: boundary pass for paint errors (goal G1a, idea I32)

Question: after the G20 strict-field fit (`Astrict/100`, painted start, flood groups), 6 of the
1,082 voxels stay wrong (39-57 deg; 3 of them among the 853 voxels with a C record): boundary voxels
painted with the neighbouring grain's orientation, so flood-fill puts them in the pick's group while
their ground-truth grain is the neighbour. Their strict costs (0.13-0.26) sit above the fitted
voxels' (about 0.08). Does a boundary pass that lets each high-cost voxel try its neighbours'
orientations from other groups recover them?

Method: a new `run_pipeline.py` stage `P/<evals>/<threshold>`: for each group member whose strict cost
is above the threshold, candidates = the painted orientations of its neighbours (within the
neighbour radius) that belong to other flood groups plus the current orientations of its in-group
neighbours; the candidate with the lowest strict cost is taken if it beats the voxel's, then a
per-voxel strict CMA polish of `evals` evaluations (sigma 0.02 deg) as the `V` stage. Arm on G20:
`Astrict/100,P/30/0.12`; output `benchmarks/grain_opt/go2b/`.

Judging rule (written before any run). Reference: GO7's clean `Astrict/100` (853 C voxels: 3 above
0.1 deg, median 0.0146 deg; all 1,082 voxels: 6 wrong above 1 deg). I32 is answered yes when the
wrong voxels among all 1,082 fall from 6 to at most 2 and the C-voxel share above 0.1 deg falls from
3 to at most 1, with the C-voxel median within 5% of 0.0146. Voxels tried and changed, and time, are
reported.

Budget: G20 `Astrict/100` takes about 3 minutes; the pass adds a few strict evaluations per gated
voxel.

### GO2c: boundary pass with polished candidates (goal G1a, idea I32, second design)

Question: GO2b's pass compared unpolished candidates and lost the right ones (the strict cost at a
candidate 0.2-0.4 deg off is worse than at the wrong painted orientation). If each gated voxel's
distinct candidates (those more than 0.5 deg apart; at most the 3 with the lowest strict cost) are
first polished with a short step ladder (60 evaluations, steps 1 / 0.3 / 0.05 deg), and the best
polished one is taken if it beats the voxel, are the paint errors recovered?

Method: `run_pipeline.py` stage `P` with a polish-first mode (`P2/<evals>/<threshold>`); arm on G20
`Astrict/100,P2/60/0.12`; output `benchmarks/grain_opt/go2b/`.

Judging rule (written before any run; GO2b's): I32 is answered yes when the wrong voxels among all
1,082 fall from 6 to at most 2 and the C-voxel share above 0.1 deg from 3 to at most 1, with the
C-voxel median within 5% of 0.0146.

Budget: about 3 x 60 strict evaluations per gated voxel (about 70 voxels): a few minutes.

### GO6 onwards

Strain queue (owner, 2026-11-19, from the literature review in the sibling thread "Review neural
methods for grain gradients"): after GO5q, I60 (Gauss-Newton gradient fit), I61 (non-linear
gradient sample), I62 (quadratic field), I63 (coordinate network, only if I62 leaves a residual),
I64 (GP / Laplacian prior on a per-voxel residual). I25 stays queued at low priority.

Program balance (owner via the coordinator, 2026-11-16): GO5b to GO5o were all G1c strain. After
GO5o the strain line is time-boxed (at most one more strain step unless a result demands it),
and the next steps rotate to goals with no steps yet, in this order: (1) G1b unblock: cache the
ray index per sample so the gated global search (I36) can run on `rand` / `corr3` at 2-5 deg;
(2) G2b richer simulation (I23 rocking curve alpha > 0, I24 partial volume / PSF); (3) G2c learned
models (I25 to I27); (4) G2a remaining intensity uses (I18 Poisson, I22 structure factor).

Planned order (each gets a section when it starts): GO6 objective variants (I17 to I22); GO7 richer
simulation (I23, I24); GO8 learned models (I25 to I27).
