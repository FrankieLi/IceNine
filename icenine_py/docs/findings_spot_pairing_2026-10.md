# Findings: matching spots across the two detectors (Task D, 2026-10-14)

## What we asked

Each diffracted ray hits both detectors, the near one (detector 0) and the far one (detector 1). So one peak makes
two spots, one on each detector, in the same omega frame. A straight line through those two spots is the ray. If we
follow that line back to the sample layer, we find where in the sample the peak came from. We call that point the
ray's **origin**.

If we can do this from the images alone, with no orientation search, we get two things:

1. For each voxel, a list of the peaks that came from it, each with an exact direction. A cost function or a
   Gauss-Newton step could then match peaks smoothly, instead of counting overlapping pixels (Tasks E and F).
2. A way to throw out spots that cannot belong to a voxel, before trying any orientation.

We tested this on the clean simulated images of the full sample (24,570 voxels). The simulator tells us which voxel
drew which spot, so every answer can be checked against the truth.

![How a spot pair gives a ray and its origin on the sample layer](../benchmarks/eval_reduction/D/figures/pairing_geometry.png)

## What we found, in short

1. **A spot almost never comes from one voxel.** All the voxels of a grain diffract into the same spot, so a typical
   spot is drawn by about 44 voxels at once. Out of 57,389 spots, only 30 come from a single voxel. About 2,600 spots
   mix several grains.

2. **We can match spots across the two detectors, as long as we also check the diffraction angle.**
   - The line through two spots must cross the sample layer inside the sample. That test alone is not enough:
     most spots have about three possible partners, and fewer than half of the matches are right.
   - We added a second test: the angle between the line and the beam must be one of the copper diffraction angles.
     This test needs only the crystal structure, not the orientation.
   - With both tests, we find 99% of the spot pairs that come from one grain. 85% of all the matches we make are
     right.
   - Spots that mix several grains are matched only 35% of the time. Their centre is an average over several grains,
     so it does not lie on any one ray.

3. **The origin is accurate, except for rays that leave the sample almost flat.** The typical error is 2.7 µm, less
   than a third of a voxel (9.4 µm). For the flattest rays it grows to about 31 µm, because a small error in the line
   moves its crossing point a long way along the layer.

4. **But the origin points to the grain, not to the voxel.** A spot is drawn by the whole grain, so its origin is the
   middle of the grain. To find a voxel's peaks, we take every ray whose origin is within some distance of the voxel.
   - Within 30 µm, we find two thirds of the voxel's real peaks, but also pick up about 45 wrong ones.
   - Within 50 µm, we find 82% of them, with about 171 wrong ones.
   There is no distance that gives most of the real peaks and few wrong ones.

5. **A cost built on these peak lists does not tell right from wrong orientations any better than the pixel cost.**
   - It is gentler, not sharper. Its score halves 0.38° away from the true orientation; the pixel cost's score
     halves at 0.26°.
   - Both reject orientations 1° or more away, as long as the peak list is good.
   - With the short list the plan suggested (rays within 15.7 µm of the voxel), the true orientation scores only
     8%, because most of the voxel's peaks are missing from the list.

6. **We did not build the classifier (Step 4).** The plan said to build it only if fewer than 95% of single-voxel
   spots were matched. That number does not exist here (there are almost no single-voxel spots). For one-grain spots
   we match 99%, so we skipped it. Two problems remain open: 15% of our matches are wrong, and spots that mix several
   grains are mostly unmatched.

## What this means for the next tasks

- **Tasks E and F can use the matched rays.** For one-grain spots we have the ray's direction and, for most rays, a
  good origin. Their value is that a ray gives a smooth error to minimise. As a plain hit-or-miss score they are no
  better than the pixel cost.
- **We cannot assign peaks to single voxels from the images alone on this sample.** Almost every spot is shared by a
  whole grain. The origin finds the grain, not the voxel.
- **As a filter before the search, this would lose too many real peaks.** Any distance that keeps the wrong peaks
  few also drops a third or more of the real ones.

## The details

### Step 1: which voxels drew each spot

We re-ran the simulator on the same clean images and recorded, for every lit pixel, which voxel and which
reflection drew it. The re-run gives exactly the same 360 images and the same 57,389 spots as before, so the record
matches the images the other studies use. (It took 35 s on 10 workers.)

Two spots form a **true pair** when they share a voxel and reflection. There are 27,888 true pairs, together
carrying 1,375,053 peaks. 3,890 spots on detector 0 have no partner on detector 1; we did not find out why.

The table counts spots by how many sources they have. "Sources" means (voxel, reflection) pairs.

<!-- table:evalD_step1_components -->
| detector | class | components | with a true partner | median pixels | p90 pixels | median sources | max sources |
|---|---|---|---|---|---|---|---|
| 0 | single-source | 10 | 0 | 10 | 11 | 1 | 1 |
| 0 | one grain, several voxels | 28411 | 24671 | 97 | 223 | 44 | 185 |
| 0 | several grains | 1770 | 1630 | 230 | 462 | 106 | 296 |
| 1 | single-source | 20 | 20 | 3 | 10 | 1 | 1 |
| 1 | one grain, several voxels | 26355 | 26355 | 103 | 235 | 44 | 185 |
| 1 | several grains | 823 | 823 | 249 | 472 | 103 | 296 |
<!-- /table:evalD_step1_components -->

### Step 2: matching spots across the detectors

For each frame we try every spot on detector 0 against every spot on detector 1. A frame has at most a few hundred
spots per detector, so this is cheap. We tried several sets of tests, named in the table as follows:

- `hull`: the line must cross the sample layer inside the sample.
- `ring`: also, the angle to the beam must match a copper diffraction angle.
- `shape`: also, the two spots must have similar size and shape.
- `tight`: each spot's centre may be off by 1 pixel.
- `plan`: adds the looser tolerances the plan listed (10 µm detector registration and 5.8 µm from the voxel size).
  The simulated images have no registration error, so this only makes the test looser.
- `c0.5`, `c2`: centre tolerance of 0.5 or 2 pixels instead of 1.
- `binary`: spot centre computed without weighting by intensity.

Our main choice is `ring_tight`. In the table, **precision** is the share of our matches that are right, and
**recall** is the share of true pairs we found. The brackets are 95% intervals.

<!-- table:evalD_step2_presets -->
| preset | predicted pairs | precision | link recall | peak recall | recall, one grain | recall, several grains |
|---|---|---|---|---|---|---|
| hull_tight | 56175 | 0.463 [0.459, 0.467] (26023/56175) | 0.933 [0.930, 0.936] (26023/27888) | 0.947 | 0.942 [0.939, 0.945] (22141/23502) | 0.885 [0.875, 0.894] (3882/4386) |
| hull_plan | 69314 | 0.376 [0.372, 0.379] (26034/69314) | 0.934 [0.931, 0.936] (26034/27888) | 0.947 | 0.942 [0.939, 0.945] (22141/23502) | 0.888 [0.878, 0.897] (3893/4386) |
| ring_tight | 29029 | 0.854 [0.850, 0.858] (24788/29029) | 0.889 [0.885, 0.892] (24788/27888) | 0.887 | 0.990 [0.988, 0.991] (23256/23502) | 0.349 [0.335, 0.364] (1532/4386) |
| ring_plan | 56027 | 0.487 [0.483, 0.491] (27283/56027) | 0.978 [0.977, 0.980] (27283/27888) | 0.980 | 0.998 [0.998, 0.999] (23459/23502) | 0.872 [0.862, 0.881] (3824/4386) |
| shape_tight | 25798 | 0.909 [0.905, 0.913] (23452/25798) | 0.841 [0.837, 0.845] (23452/27888) | 0.844 | 0.966 [0.964, 0.968] (22705/23502) | 0.170 [0.159, 0.182] (747/4386) |
| shape_plan | 44653 | 0.545 [0.540, 0.550] (24338/44653) | 0.873 [0.869, 0.877] (24338/27888) | 0.884 | 0.974 [0.972, 0.976] (22897/23502) | 0.329 [0.315, 0.343] (1441/4386) |
| ring_tight_c0.5 | 25910 | 0.927 [0.923, 0.930] (24007/25910) | 0.861 [0.857, 0.865] (24007/27888) | 0.856 | 0.985 [0.983, 0.986] (23139/23502) | 0.198 [0.186, 0.210] (868/4386) |
| ring_tight_c2 | 35907 | 0.723 [0.718, 0.727] (25954/35907) | 0.931 [0.928, 0.934] (25954/27888) | 0.932 | 0.993 [0.992, 0.994] (23344/23502) | 0.595 [0.580, 0.610] (2610/4386) |
| shape_tight_c0.5 | 24090 | 0.957 [0.955, 0.960] (23058/24090) | 0.827 [0.822, 0.831] (23058/27888) | 0.826 | 0.961 [0.959, 0.964] (22588/23502) | 0.107 [0.098, 0.117] (470/4386) |
| shape_tight_c2 | 29689 | 0.807 [0.802, 0.811] (23948/29689) | 0.859 [0.855, 0.863] (23948/27888) | 0.868 | 0.970 [0.967, 0.972] (22788/23502) | 0.264 [0.252, 0.278] (1160/4386) |
| ring_tight_binary | 28200 | 0.849 [0.845, 0.853] (23942/28200) | 0.859 [0.854, 0.863] (23942/27888) | 0.852 | 0.957 [0.954, 0.959] (22488/23502) | 0.332 [0.318, 0.346] (1454/4386) |
| shape_tight_binary | 25027 | 0.906 [0.903, 0.910] (22684/25027) | 0.813 [0.809, 0.818] (22684/27888) | 0.812 | 0.935 [0.931, 0.938] (21968/23502) | 0.163 [0.153, 0.174] (716/4386) |
<!-- /table:evalD_step2_presets -->

Where `ring_tight` goes wrong:

- Of its 4,241 wrong matches, 4,149 join two spots that each have a different true partner.
- Of the true pairs it misses, 3,036 fail the angle test. Most of these (2,813) involve a spot that mixes several
  grains. Only 36 misses fail the inside-the-sample test.

**How far off the origin is.** The next table is for correct matches of one-grain spots. It groups them by how
steeply the ray leaves the sample layer (the z part of the ray's direction). The true origin is the average position
of the spot's voxels. The last column is the error the plan's tolerances predict. It is 20 to 30 times the real
error, so those tolerances work as a loose limit, not as an error estimate.

<!-- table:evalD_step2_origin_error -->
| ray z-slope (abs u_z) | pairs | lab-frame error median (um) | lab-frame error p90 (um) | sample-frame error median (um) | sample-frame error p90 (um) | predicted sigma median (um) |
|---|---|---|---|---|---|---|
| 0 to 0.02 | 2131 | 30.8 | 92.8 | 30.7 | 92.8 | 364.8 |
| 0.02 to 0.05 | 5569 | 4.2 | 19.7 | 4.7 | 19.8 | 180.7 |
| 0.05 to 0.1 | 5838 | 1.6 | 4.3 | 2.4 | 5.2 | 74.7 |
| 0.1 to 1 | 9714 | 1.0 | 2.6 | 1.9 | 4.0 | 37.0 |
<!-- /table:evalD_step2_origin_error -->

The map below shows the origin error across the sample. It has no spatial pattern.

![Origin error over the sample](../benchmarks/eval_reduction/D/figures/origin_error_map.png)

### Step 3: assigning peaks to a voxel, and does it help?

We tried two ways to decide which matched rays belong to a voxel:

- **By distance** (the plan's way): keep the rays whose origin is within r of the voxel's centre.
- **By footprint** (added): start a ray at the voxel's centre, in the matched pair's direction, and keep the pair if
  that ray hits both spots.

We checked these on 30 voxels from the earlier ray study, which include the six full-sample voxels used before. (The
plan asked for the six plus 30 others; we used the 30.) A run on 300 voxels spread over the sample gives the same
picture.

In the table, "true peaks found" is the share of the voxel's real peaks that we keep, and "false attached" is how
many wrong peaks we keep per voxel. "Interior" and "boundary" mean voxels inside a grain or next to a grain boundary.
Even a perfect rule could find at most 89% of the real peaks, because 11% of them sit in spot pairs we never matched.

<!-- table:evalD_step3_attach_voxels30 -->
| attachment | true peaks found | attached pairs per voxel | false attached per voxel | share of attached that are false | recall, interior | recall, boundary |
|---|---|---|---|---|---|---|
| A r=2 um | 0.004 [0.002, 0.008] (7/1714) | 0.3 | 0.1 | 0.222 | 0.008 | 0.002 |
| A r=5 um | 0.020 [0.014, 0.028] (34/1714) | 1.7 | 0.5 | 0.320 | 0.049 | 0.003 |
| A r=10 um | 0.127 [0.112, 0.144] (218/1714) | 9.6 | 2.3 | 0.240 | 0.320 | 0.016 |
| A r=15 um | 0.231 [0.212, 0.252] (396/1714) | 19.7 | 6.5 | 0.330 | 0.494 | 0.079 |
| A r=20 um | 0.372 [0.349, 0.395] (637/1714) | 35.8 | 14.6 | 0.407 | 0.546 | 0.271 |
| A r=30 um | 0.663 [0.640, 0.685] (1136/1714) | 82.8 | 44.9 | 0.542 | 0.670 | 0.658 |
| A r=50 um | 0.821 [0.802, 0.838] (1407/1714) | 218.1 | 171.2 | 0.785 | 0.817 | 0.823 |
| A r=100 um | 0.863 [0.846, 0.878] (1479/1714) | 817.5 | 768.2 | 0.940 | 0.858 | 0.866 |
| B footprint | 0.851 [0.833, 0.867] (1458/1714) | 386.8 | 338.2 | 0.874 | 0.852 | 0.850 |
<!-- /table:evalD_step3_attach_voxels30 -->

![Real peaks found against wrong peaks kept, for each rule](../benchmarks/eval_reduction/D/figures/attribution.png)

**Does a peak-list score tell right from wrong orientations?** For each voxel we scored the true orientation and
many wrong ones. The wrong ones were turned away from the truth by 0.05° to 20° about random axes, plus larger turns
about the <100>, <110> and <111> axes. We compared two scores:

- the **pixel score**: 1 minus the usual pixel-overlap cost (strict C++ rule);
- the **peak-list score**: the share of the peaks this orientation predicts that match a ray in the voxel's list.
  A match means the same frame and a direction within 0.05°, about one pixel of error on each detector.

"Oracle" means the voxel's list is exactly its true peaks: the best any assignment rule could do.

The first table gives each score at the true orientation. The "half-score radius" is how far you can turn away
from the truth before the score falls to half its value at the truth.

<!-- table:evalD_step3_truth_scores -->
| score | median at the truth | minimum at the truth | half-score radius (deg) |
|---|---|---|---|
| pixel overlap (C++ rule) | 0.936 | 0.894 | 0.258 |
| attributed, A r=15.7 um | 0.078 | 0.000 | 0.384 |
| attributed, A r=30 um | 0.719 | 0.018 | 0.376 |
| attributed, B footprint | 0.792 | 0.696 | 0.377 |
| attributed, oracle (true pairs) | 0.789 | 0.696 | 0.378 |
<!-- /table:evalD_step3_truth_scores -->

The next table shows how each score falls as the orientation turns away from the truth. Each score is divided by its
value at the truth.

<!-- table:evalD_step3_sharpness -->
| misorientation (deg) | pixel overlap (C++ rule) | attributed, A r=15.7 um | attributed, A r=30 um | attributed, B footprint | attributed, oracle (true pairs) |
|---|---|---|---|---|---|
| 0.05 | 0.927 | 0.976 | 0.967 | 0.958 | 0.958 |
| 0.1 | 0.819 | 0.944 | 0.922 | 0.921 | 0.920 |
| 0.25 | 0.512 | 0.764 | 0.732 | 0.731 | 0.733 |
| 0.5 | 0.235 | 0.338 | 0.338 | 0.342 | 0.342 |
| 1 | 0.061 | 0.072 | 0.085 | 0.087 | 0.088 |
| 2 | 0.017 | 0.012 | 0.017 | 0.020 | 0.018 |
| 5 | 0.004 | 0.000 | 0.003 | 0.003 | 0.003 |
| 10 | 0.001 | 0.000 | 0.000 | 0.000 | 0.000 |
| 20 | 0.003 | 0.000 | 0.000 | 0.000 | 0.000 |
<!-- /table:evalD_step3_sharpness -->

![Score against distance from the true orientation](../benchmarks/eval_reduction/D/figures/sharpness.png)

The last table asks whether any wrong orientation scores close to the truth. "Best wrong / truth" is the highest
wrong score divided by the true score. The last two columns count voxels where a wrong orientation scores at least
half of, or at least as much as, the truth.

<!-- table:evalD_step3_separation -->
| wrong orientations | score | median best wrong / truth | max over voxels | voxels >= 0.5 truth | voxels >= truth |
|---|---|---|---|---|---|
| random axis, 1 degree and more | pixel overlap (C++ rule) | 0.110 | 0.173 | 0 | 0 |
| random axis, 1 degree and more | attributed, A r=15.7 um | 0.178 | 1.035 | 5 | 1 |
| random axis, 1 degree and more | attributed, A r=30 um | 0.138 | 0.983 | 2 | 0 |
| random axis, 1 degree and more | attributed, B footprint | 0.139 | 0.260 | 0 | 0 |
| random axis, 1 degree and more | attributed, oracle (true pairs) | 0.139 | 0.260 | 0 | 0 |
| random axis, 0.25 to 0.5 degree | pixel overlap (C++ rule) | 0.640 | 0.818 | 28 | 0 |
| random axis, 0.25 to 0.5 degree | attributed, A r=15.7 um | 0.962 | 1.035 | 23 | 11 |
| random axis, 0.25 to 0.5 degree | attributed, A r=30 um | 0.849 | 1.526 | 30 | 3 |
| random axis, 0.25 to 0.5 degree | attributed, B footprint | 0.830 | 0.921 | 30 | 0 |
| random axis, 0.25 to 0.5 degree | attributed, oracle (true pairs) | 0.830 | 0.921 | 30 | 0 |
| rotations about <100>, <110>, <111> | pixel overlap (C++ rule) | 0.310 | 0.350 | 0 | 0 |
| rotations about <100>, <110>, <111> | attributed, A r=15.7 um | 0.316 | 1.000 | 4 | 1 |
| rotations about <100>, <110>, <111> | attributed, A r=30 um | 0.309 | 0.667 | 2 | 0 |
| rotations about <100>, <110>, <111> | attributed, B footprint | 0.310 | 0.435 | 0 | 0 |
| rotations about <100>, <110>, <111> | attributed, oracle (true pairs) | 0.307 | 0.435 | 0 | 0 |
<!-- /table:evalD_step3_separation -->

**The 0.05° match tolerance sets how sharp the peak-list score is.** With 0.02° the score halves at 0.18°; with 0.2°
it halves at 0.72°. We chose 0.05° from the pixel size before running this test, so it was not tuned. The full scan is
in `step3/tables_separate.md`. We did not test whether a smooth direction error helps a search.

## Not tested, or left open

- Splitting spots that mix several grains.
- Noisy images.
- A distance limit set per ray from its predicted origin error.
- Using the peak lists in a search (Tasks E and F).
- Why 3,890 detector-0 spots have no partner.
- A sample layer with real thickness (the simulated layer is a plane at z = 0).
- Whether a classifier would cut the 15% wrong matches or recover the several-grain spots.

## Files

- Library: `icenine/spot_pairing.py`, tests in `tests/test_spot_pairing.py`.
- Simulator change: the opt-in `ForwardSimulation.peak_hook`. It is off by default and then changes nothing.
- Recording script: `scripts/phase_d/render_provenance.py`. It sits beside `render_full.py`, which is unchanged.
- Study scripts: `scripts/eval_reduction/pairing_*.py`.
- Run record: `benchmarks/eval_reduction/D/RUN_LOG.md`.
- Results and figures: `benchmarks/eval_reduction/D/`. The tables are generated (`doc_tables_all.md`) and copied in
  with `scripts/dev/sync_doc_tables.py`.

All runs use clean images, the Q-max 8 configuration and sample seed 0. No golden test or default was changed. The
same results are in the "Task D" section of `MIGRATION_HISTORY.md`.
