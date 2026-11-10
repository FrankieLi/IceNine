# AR3 step 1: painting grains on the pilot region

Step 1 of "AR3: index the layer and paint grains" in `docs/autoresearch/PLAN.md` (2026-11-10),
under D-20261110-N1. This is a study: `icenine/grain_paint.py` runs outside the reconstructor, and
no config key switches it on. Scripts: `scripts/autoresearch/ar3/`. Results:
`benchmarks/autoresearch/ar3/` (`summary.json`, `tables.md`, `ar3_paint.png`).

How this page was made. The implementer's write of a new `.md` file is refused by the tool, so it
returned the text and the main session wrote the file. The tables are copied from `tables.md` by
`scripts/dev/sync_doc_tables.py`, not typed.

## The result in one picture

Reading the picture. Left: distance of the painted voxels to the truth, painted without a fit
(orange) and fitted one by one by the key-on BFS (blue), over the 1,978 voxels that are right in
both runs; the two curves lie on top of each other. Middle: the side arm (a local fit on contested
voxels) against wrong voxels; every point is at 0. Right: voxels where a grain passes outside its
own true grain, by distance from that grain.

![Accuracy without a fit, the side arm, and how far a wrong grain floods](../../benchmarks/autoresearch/ar3/ar3_paint.png)

- Painting the 2,000-voxel pilot costs 156,766 evaluations. The key-on BFS costs 2,819,906, the
  key-off BFS 25,547,935. Painting is 5.6% of the first and 0.6% of the second.
- Painting has 0 wrong voxels of 2,000 (0.0 to 0.2%). Key on has 22, key off 37. Against key on:
  22 fixed, 0 broken.
- Of the 126 voxels where two or more grains pass, the lowest strict cost names the true grain in
  126 of 126 (97 to 100%). 44 of them lie between grains of different orientation; there it is 44
  of 44 (92 to 100%).
- Without a fit, painted voxels are as close to the truth as fitted ones: median 0.0169 degrees
  against 0.0159.
- The gate for step 2 is met: 0 breaks against key on (limit 20), and 5.6% of key on's evaluations
  (limit 50%).
- Main caveat: the truth of this sample has exactly one orientation per grain. So the premise of
  painting is true by construction here.

## Terms

- **Grain entry**: one accepted seed (voxel, orientation, hit ratio).
- **Hit ratio**: the share of predicted spots that fall on lit pixels.
- **Score** a grain on a voxel: one strict evaluation of the grain's orientation there, with no
  fit. It gives a cost and a hit ratio.
- A grain **passes** on a voxel by the BFS neighbour test, unchanged: hit ratio above 0.9 times the
  best seen so far in that grain's flood.
- **Flood**: score a grain on the neighbours of its seed, then on the neighbours of every voxel
  where it passes. A flood does not depend on other grains.
- **Contested**: two or more grains pass. **Unpainted**: none passes.
- **Paint**: give the voxel the passing grain with the lowest strict cost.
- **Wrong**: over 1 degree from the truth (crystal symmetry counted).
- **Edge-sharing**, **corner-only**, **interior**: another true grain shares an edge with the
  voxel, only a corner, or no vertex.
- **Key on**, **key off**: the BFS run with and without the indexing seed of AR2.

## What was built and run

`grain_paint.paint_grains` takes the voxels in the shuffled order of the AR2 BFS runs (no new
draw). For the next voxel that is unpainted and not settled:

1. If a grain was scored on it, one local fit starts from the lowest-cost such grain. If the fit
   passes the neighbour test, the voxel is settled by a fit.
2. Otherwise the seed path runs (index, strict cost on 5 candidates, finisher on 2). If the result
   passes the seed test, it becomes a grain entry and is flooded.
3. Otherwise the grid search runs (the fallback). A voxel it rejects is failed.

Then each voxel takes the passing grain of lowest cost. The side arm fits each contested voxel from
its two lowest-cost grains and keeps the lower final cost. Clean images, Q-max 8, FusedCost on, one
worker.

Baseline reproduced: the key-on profile rerun made 2,816,897 evaluations, equal to the stored
FusedCost key-on run. The stored key-off numbers were recomputed: 37 wrong, 25,547,935
evaluations, 69 full searches.

## Evaluations

Reading the table. One row per arm, with the evaluations by stage. The two stored key-on rows
differ only in FusedCost. The last row adds the side arm to paint.

<!-- table:ar3_evals -->
| arm | total | seed paths | scoring | local fits | grid searches | wrong |
|---|---|---|---|---|---|---|
| key off (stored rerun) | 25,547,935 | 0 | 0 | 2,665,070 | 22,882,865 | 37 of 2000 (1.3 to 2.5%) |
| key on (stored) | 2,819,906 | 122,198 | 0 | 2,697,708 | 0 | 22 of 2000 (0.7 to 1.7%) |
| key on, FusedCost (stored) | 2,816,897 | 122,154 | 0 | 2,694,743 | 0 | 22 of 2000 (0.7 to 1.7%) |
| paint (FusedCost) | 156,766 | 119,033 | 6,697 | 31,036 | 0 | 0 of 2000 (0.0 to 0.2%) |
| paint + fit the contested (side arm) | 324,092 | 119,033 | 6,697 | 198,362 | 0 | 0 of 2000 (0.0 to 0.2%) |
<!-- /table:ar3_evals -->

- Seed paths are the largest share of paint: 119,033 of 156,766 (76%), for 62 seeds (1,920 each).
- 48 local fits were rejected before their voxels became seeds: 31,036 evaluations (20%, 647
  each). None was accepted, so "settled by a fit" ended no voxel.
- Scoring: 6,697 evaluations (4.3%). That is about 108 scores per entry, and 33 passes on average.
- The grid search ran on no voxel: all 62 seeds passed the seed test.
- The one-time pairing (17.2 s) is outside the count.

## Wrong voxels, paired

Reading the table. Each row compares paint with a stored run, voxel by voxel. Fixed: wrong there,
right in paint. Broken: the opposite.

<!-- table:ar3_paired -->
| vs | wrong there | wrong in paint | fixed | broken | wrong in both |
|---|---|---|---|---|---|
| key on (stored) | 22 | 0 | 22 | 0 | 0 |
| key on, FusedCost (stored) | 22 | 0 | 22 | 0 | 0 |
| key off (stored) | 37 | 0 | 37 | 0 | 0 |
<!-- /table:ar3_paired -->

Reading the table. The wrong voxels of paint by kind. All rows are 0.

<!-- table:ar3_wrong -->
| what | n |
|---|---|
| wrong voxels of paint | 0 |
|   of which crossing a boundary, sharing an edge | 0 |
|   of which crossing a boundary, corner only | 0 |
|   of which interior voxels | 0 |
<!-- /table:ar3_wrong -->

- The 22 voxels wrong in key on all share an edge with another grain. In paint, 18 of them are
  contested between two grains, 1 between three, and 3 have one passing grain. All 22 are painted
  right.
- Paint has 0 wrong voxels, so there is no break to trace and no boundary crossing to split by
  edge or corner.
- Two key-off runs alone differ by 20 voxels one way and 7 the other (AR2). 22 against 0 is
  larger. But paint and the BFS differ in two ways at once (the route and the per-voxel fit), and
  neither is credited alone.
- Entries: 62 for 60 true grains. Every true grain has one. 2 grains have two (0.055 and 0.041
  degrees apart, kept as their own entries). No entry is wrong.

## The direct test of the boundary failure

Reading the table. Each row is a set of contested voxels. The columns say how often the true grain
passes, how often the lowest cost names it, how often the highest hit ratio names it, and how often
those two rules name different grains. The rows that say "different orientations" are the
contested voxels whose passing grains are over 1 degree apart.

<!-- table:ar3_direct -->
| set | n | true grain passes | lowest cost is the true grain | highest hit is the true grain | lowest cost and highest hit differ |
|---|---|---|---|---|---|
| all: all kinds | 126 | 126 of 126 (97 to 100%) | 126 of 126 (97 to 100%) | 126 of 126 (97 to 100%) | 0 of 126 (0 to 3%) |
| all: edge | 63 | 63 of 63 (94 to 100%) | 63 of 63 (94 to 100%) | 63 of 63 (94 to 100%) | 0 of 63 (0 to 6%) |
| all: corner | 32 | 32 of 32 (89 to 100%) | 32 of 32 (89 to 100%) | 32 of 32 (89 to 100%) | 0 of 32 (0 to 11%) |
| all: interior | 31 | 31 of 31 (89 to 100%) | 31 of 31 (89 to 100%) | 31 of 31 (89 to 100%) | 0 of 31 (0 to 11%) |
| lowest cost and highest hit different orientations: all kinds | 44 | 44 of 44 (92 to 100%) | 44 of 44 (92 to 100%) | 44 of 44 (92 to 100%) | 0 of 44 (0 to 8%) |
| lowest cost and highest hit different orientations: edge | 44 | 44 of 44 (92 to 100%) | 44 of 44 (92 to 100%) | 44 of 44 (92 to 100%) | 0 of 44 (0 to 8%) |
<!-- /table:ar3_direct -->

- Of the 126 contested voxels, 44 lie between grains of different orientation, and all 44 are
  edge-sharing. The other 82 are voxels of the two grains that have two entries.
- On the 44, the true grain passes and has the lowest cost in 44 of 44 (92 to 100%). So the cost
  tells touching grains apart once both are scored, whatever the BFS route.
- The lowest cost and the highest hit ratio never name different grains here: 0 of 126 (0 to 3%).

## How far a wrong grain floods

Reading the table. Voxels where a grain passes outside its own true grain, by distance from that
grain (in voxel sides).

<!-- table:ar3_flood -->
| distance from the grain's true region (voxel sides) | n |
|---|---|
| 0 to 1 | 44 |
| 1 to 2 | 0 |
| 2 to 3 | 0 |
| 3 to 4 | 0 |
| 4 or more | 0 |
<!-- /table:ar3_flood -->

- 44 passes lie outside the true grain, against 2,021 inside. All 44 are within 1 side and
  edge-sharing. 26 of the 62 entries have at least one; the most is 6.
- No wrong grain passes beyond its neighbours on this sample (same relative test as the BFS).

## Accuracy without a fit

Reading the table. Distance to the truth (degrees) over the 1,978 voxels painted by a grain that
are right in both runs. Painted voxels carry their grain's orientation unchanged.

<!-- table:ar3_accuracy -->
| arm | median | 90th | 99th | largest |
|---|---|---|---|---|
| paint (no fit) | 0.0169 | 0.0362 | 0.0594 | 0.0594 |
| key on (fitted) | 0.0159 | 0.0384 | 0.0577 | 0.1085 |
<!-- /table:ar3_accuracy -->

- Without a fit, the painted orientation is as accurate as the fitted one on this subset. This is
  consistent with no final fit per voxel or per grain being needed (goal 4). It is not shown for
  samples whose grains hold several orientations.

## The side arm

Reading the table. The rule: fit a contested voxel when its two lowest costs differ by less than
m. Columns: voxels fitted, fits, evaluations, fixed, broken, wrong after.

<!-- table:ar3_margin -->
| fit if the two lowest costs differ by less than | voxels | fits | evals | fixed | broken | wrong |
|---|---|---|---|---|---|---|
| 0.001 | 8 | 16 | 11,094 | 0 | 0 | 0 |
| 0.003 | 12 | 24 | 16,390 | 0 | 0 | 0 |
| 0.01 | 20 | 40 | 26,982 | 0 | 0 | 0 |
| 0.03 | 66 | 132 | 87,886 | 0 | 0 | 0 |
| 0.1 | 120 | 240 | 159,382 | 0 | 0 | 0 |
| 0.3 | 126 | 252 | 167,326 | 0 | 0 | 0 |
| 1 | 126 | 252 | 167,326 | 0 | 0 | 0 |
| all | 126 | 252 | 167,326 | 0 | 0 | 0 |
<!-- /table:ar3_margin -->

- 252 fits on 126 contested voxels cost 167,326 evaluations and 11.0 s (contended). 0 fixed, 0
  broken: paint has nothing to fix on this run. Paint plus the side arm is 324,092 evaluations,
  11.5% of key on.
- The curve picks no margin: every m gives 0 wrong. A margin rule is tested on voxels outside the
  pilot (step 2).

## Where the time goes

Reading the table. The key-on pilot (FusedCost, one worker, timers only, contended), split by
part. Under cProfile the same run took 235.6 s with the same evaluations. Its top functions by own
time are the fused evaluation (106 s) and the orientation sampling in the local fits (`slerp`,
39 s).

<!-- table:ar3_profile -->
| part | s | calls or evaluations | share |
|---|---|---|---|
| local fits, cost evaluations | 121.2 | 2,694,743 | 54.1% |
| local fits, optimizer's own work | 65.5 | 4,095 | 29.2% |
| seed path: indexing | 11.9 | 64 | 5.3% |
| seed path: finisher, cost evaluations | 4.3 | 121,834 | 1.9% |
| seed path: finisher, optimizer's own work | 3.1 | 128 | 1.4% |
| seed path: strict costs on candidates and the rest of the seed path | 0.0 | 64 | 0.0% |
| overlap evaluations (the BFS's own) | 0.0 | 64 | 0.0% |
| grid search | 0.0 | 0 | 0.0% |
| one-time pairing (inside the first seed path) | 17.6 | n/a | 7.8% |
| everything else (BFS bookkeeping, neighbours, records) | 0.5 | n/a | 0.2% |
<!-- /table:ar3_profile -->

- Cost evaluations inside local fits are 54% of the key-on time. The optimizer's own work in them
  is 29%. So batching the evaluations alone could save at most a little over half.
- Indexing plus finisher is 9%, and the pairing 8%.

Reading the table. Wall times of this session, in seconds. All are contended, one worker each.

<!-- table:ar3_wall -->
| what | s |
|---|---|
| paint run (FusedCost, unbatched loop) | 37.9 |
|   of which the one-time pairing | 17.2 |
| key on, FusedCost, timers only, this session | 224.3 |
| key on, FusedCost, stored (AR2 step 2) | 213.1 |
| key on under cProfile (not a timing) | 235.6 |
| side arm (fits on the contested voxels) | 11.0 |
<!-- /table:ar3_wall -->

- The paint loop scores one voxel at a time in Python. Its time is a first reading of an unbatched
  loop, not a tested wall time.
- The timing preflight refused every run (an editor process, cpptools, at about 98% CPU, as in
  AR2). So every time here is labelled contended.

## Is the premise true on this sample

Reading the table. For the pilot region and the whole layer: the number of true grains and voxels,
the largest spread of orientation inside one true grain, and how many grains spread over 0.001 and
over 1 degree.

<!-- table:ar3_premise -->
| scope | n_grains | n_voxels | max_spread_deg | n_grains_over_1e3_deg | n_grains_over_1_deg |
|---|---|---|---|---|---|
| region | 60 | 2000 | 0.00000 | 0 | 0 |
| layer | 497 | 24570 | 0.00000 | 0 | 0 |
<!-- /table:ar3_premise -->

- The truth holds exactly one orientation per grain (largest spread 0.0 degrees). This says
  nothing for samples with several orientations per grain.

## Criteria

- Evaluations by stage against key on and key off: met.
- Wrong voxels with intervals, paired, breaks traced: met (none to trace).
- Direct test with intervals: met (44 of 44).
- Gate for step 2: met.

## What was not tested

- The whole layer (a larger hull pairs rays differently).
- Samples with several orientations per grain.
- Realistic images, another Q-max, another voxel order (one fixed order was used).
- The fallback and the "settled by a fit" branch (they ended no voxel).
- A margin rule (none was needed here).
- A tested wall time, and a batched version.
- Why 2 grains got two entries (not traced).
