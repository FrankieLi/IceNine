# Findings: spot pairing across detectors and peak-to-voxel attribution (Task D, 2026-10-14)

Branch `feature/eval-reduction-D`. Run record: `benchmarks/eval_reduction/D/RUN_LOG.md`. Library:
`icenine/spot_pairing.py` (pure geometry, numpy, `cKDTree` for the attribution; tests `tests/test_spot_pairing.py`),
plus the opt-in `ForwardSimulation.peak_hook` (default `None`, no effect when off). Provenance writer:
`scripts/phase_d/render_provenance.py` (a sibling of `render_full.py`, which is unchanged). Scripts:
`scripts/eval_reduction/pairing_*.py`. Tables below are generated (`benchmarks/eval_reduction/D/doc_tables_all.md`,
synced with `scripts/dev/sync_doc_tables.py`); figures are in `benchmarks/eval_reduction/D/figures/`. Everything is
the clean full-sample render (24,570 voxels, Q-max 8 configuration, sample seed 0); no realistic images, no extra
seeds, no golden test changed, no default flipped. The same results, with the same tables, are in the "Task D"
subsection of `MIGRATION_HISTORY.md`.

## The question

The two detector distances see the same diffracted ray twice. A ray is a line, so two spots in one omega frame that
lie on one line through the illuminated layer are one peak, and the line's crossing of the layer is where the peak
came from. If that works once per image set, with no search over orientations, it gives (i) per voxel an observed
peak list (direction, frame, origin) that a cost or a Gauss-Newton step can match continuously (Tasks E, F), and
(ii) a way to drop spots that cannot belong to a voxel before any orientation is tried.

![Pairing geometry: two detectors, the ray, its origin on the layer, the tolerance terms](../benchmarks/eval_reduction/D/figures/pairing_geometry.png)

## Short answer

1. **Almost no spot has one source.** Of 30,191 detector-0 and 27,198 detector-1 spots (4-connected components of lit
   pixels), only 10 and 20 are drawn by a single (voxel, reflection) pair. The others are unions of all the voxels of
   one grain that diffract into the same spot (median 44 sources, about 100 pixels) or of several grains (1,770 on
   detector 0). No true pair has two single-source spots (0 of 27,888), so the plan's "single-source recall" cannot
   be measured on this sample; the classes used below are single-source, one grain, several grains.
2. **Geometry pairs the one-grain spots, but only with the ring angle.** Rule: the line through the two spot
   centroids meets the layer z = 0 inside the sample outline (widened by the origin's own tolerance), and its
   scattering angle matches one of the eight Cu ring angles; 1 pixel centroid tolerance; greedy assignment.
   Recall 0.990 on one-grain pairs (23,256 of 23,502), 0.349 on pairs involving a union of several grains
   (1,532 of 4,386); precision 0.854 over all predicted pairs (24,788 of 29,029). With the outline test alone,
   precision is 0.463 and 80% of the detector-0 spots that have a partner have more than one candidate (median 3).
   The ring gate uses the crystal structure, not the orientation, so it is still a search-free test.
3. **The origin is good except for rays nearly in the layer.** Correct one-grain pairs: median error 2.7 µm,
   90th percentile 15.7 µm, 99th percentile 88 µm (voxel side 9.4 µm). The median is 31 µm when the ray's slope out of
   the layer is below 0.02 and 1.9 µm above 0.1. The map over the sample plane shows no spatial pattern.
4. **Attribution by origin alone is weak.** A union spot's origin is the centre of the grain's footprint, not the
   voxel's, so a radius around the voxel centroid trades recall for false pairs. On the 30-voxel set, r = 30 µm finds
   66.3% of a voxel's true peaks with about 45 false pairs attached per voxel; r = 50 µm finds 82.1% with about 171.
   The ceiling is 0.887 (the share of true peaks whose spot pair is found at all).
5. **The attributed-peak cost is not sharper than the pixel cost.** It is smoother: with a 0.05 degree direction
   tolerance its score falls to half at 0.38 degrees from the truth, against 0.26 degrees for the pixel-overlap cost.
   At 1 degree and more both reject wrong orientations when the attributed list is good (best wrong over truth 0.139
   against 0.110). The short list from the plan's radius (r = 15.7 µm) scores only 0.078 at the truth (median).
6. **Step 4 (classifier) did not run.** The plan's rule runs it when single-source recall is below 0.95; that recall
   is undefined here, and the nearest class (one grain) has 0.990. Precision (0.854) and the several-grain spots are
   left open.

## Step 1: ground truth from the simulator

`render_provenance.py` re-renders the clean images through a recording image (35 s on 10 workers, contended). All 360
frame images equal the stock `clean/` files and the 57,389 components equal the render's spot count, so the
provenance is of the same images the studies use. 27,888 true pairs (spot pairs sharing a source) carry 1,375,053
peaks. 3,890 detector-0 spots have no partner on detector 1; the cause was not traced.

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

## Step 2: geometric pairing

For each frame, every detector-0 spot is tested against every detector-1 spot (a broadcast; at most a few hundred
spots per detector per frame). Presets: `hull_*` (outline gate only), `ring_*` (outline and ring angle), `shape_*`
(adds a rigid pixel-count and bounding-box gate). `tight` is a 1 pixel centroid tolerance on both detectors; `plan`
adds the registration (10 µm) and voxel size x (L1/L0 - 1) (5.8 µm) on detector 1, the ray rule's arm A terms read as
half-widths. The simulated data has no registration error, so `plan` only loosens the gate. `_c0.5` and `_c2` are
centroid tolerances of 0.5 and 2 pixels; `_binary` is the unweighted centroid. The main pairing is `ring_tight`.
Rates carry Wilson 95% intervals; peak recall is the share of the 1,375,053 true peaks whose spot pair is found.

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

Of the 4,241 false pairs of `ring_tight`, 4,149 join two spots that both have another true partner, and 92 involve a
spot with no partner. Of the true pairs it misses, 3,036 fail the ring gate (223 one-grain, 2,813 several-grain) and
36 the outline gate; none fails the ray or origin-tolerance gate. That several-grain pairs mostly fail the ring gate
fits the picture of a centroid that mixes different grains and so is not on any one ray; the unions were not split
to test it.

**Origin error by ray slope** (correct one-grain pairs; the truth is the mean of the source voxel centroids). The
`plan` tolerance predicts a sigma 20 to 30 times the median error: it works as a gate, not as an error estimate.

<!-- table:evalD_step2_origin_error -->
| ray z-slope (abs u_z) | pairs | lab-frame error median (um) | lab-frame error p90 (um) | sample-frame error median (um) | sample-frame error p90 (um) | predicted sigma median (um) |
|---|---|---|---|---|---|---|
| 0 to 0.02 | 2131 | 30.8 | 92.8 | 30.7 | 92.8 | 364.8 |
| 0.02 to 0.05 | 5569 | 4.2 | 19.7 | 4.7 | 19.8 | 180.7 |
| 0.05 to 0.1 | 5838 | 1.6 | 4.3 | 2.4 | 5.2 | 74.7 |
| 0.1 to 1 | 9714 | 1.0 | 2.6 | 1.9 | 4.0 | 37.0 |
<!-- /table:evalD_step2_origin_error -->

![Origin error over the sample plane](../benchmarks/eval_reduction/D/figures/origin_error_map.png)

## Step 3: peak-to-voxel attribution and the discriminating test

Rule A (the plan's): a pair's origin within r of the voxel centroid. Rule B (added): the ray from the voxel centroid
along the pair's direction hits the bounding box of both spots (tolerance 1 pixel plus the arc motion). A voxel's true
peaks are those it drew on both detectors; an attached pair is true when the voxel is a source of both spots with a
shared reflection. The 30 voxels are the ray-study set, which contains the six interior full-sample voxels (the plan's
"six and 30 others" is 30 here, not 36). `step3/attach_summary.json` also has 300 evenly spread voxels with the same
picture (r = 30 µm 70.2%, r = 50 µm 83.0%).

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

![True peaks found against false pairs attached, by rule](../benchmarks/eval_reduction/D/figures/attribution.png)

**The discriminating test.** At the true orientation and at wrong ones (8 random axes at each of 0.05 to 20 degrees,
and rotations of 10 to 60 degrees about <100>, <110>, <111>), two scores per voxel. `pixel` is 1 - cost with the
strict C++ rule, radius 0. The attributed score is the share of the orientation's predicted observable peaks that
match an attributed observed peak (same frame, ray direction within 0.05 degrees, which is one pixel of centroid
error on both detectors over the 2 mm detector separation). The oracle list is the voxel's own true pairs. Scores
are normalised by the voxel's score at the truth.

<!-- table:evalD_step3_truth_scores -->
| score | median at the truth | minimum at the truth | half-score radius (deg) |
|---|---|---|---|
| pixel overlap (C++ rule) | 0.936 | 0.894 | 0.258 |
| attributed, A r=15.7 um | 0.078 | 0.000 | 0.384 |
| attributed, A r=30 um | 0.719 | 0.018 | 0.376 |
| attributed, B footprint | 0.792 | 0.696 | 0.377 |
| attributed, oracle (true pairs) | 0.789 | 0.696 | 0.378 |
<!-- /table:evalD_step3_truth_scores -->

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

![Score against misorientation, pixel cost and attributed costs](../benchmarks/eval_reduction/D/figures/sharpness.png)

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

The attributed cost's sharpness is set by its direction tolerance: 0.02 degrees gives a half-score radius of 0.18
degrees (completeness at the truth 0.736), 0.2 degrees gives 0.72 degrees (0.825); a one-frame window widens it again
(`step3/tables_separate.md`, `step3_tolerance`). The 0.05 degree default was set from the pixel tolerance before the
scan, so the table is evidence, not a tuned result. Whether a smooth direction residual helps a search was not tested.

## What this means for the next tasks

* For Tasks E and F, an observed peak list is available for one-grain spots (recall 0.990, origin error a few µm on
  steep rays). Its value is continuity (a direction residual), not discrimination: as a hit-or-miss score it is no
  sharper than the pixel cost.
* Attributing peaks to a single voxel from the images alone does not work well on this sample, because nearly every
  spot is a union of a grain's voxels. The origin locates the grain, not the voxel.
* As a pre-search filter, the attributed list from a radius rule would drop many true peaks at any radius that keeps
  the false ones few (66% found at about 45 false per voxel).

## Not tested / left open

Splitting unions of several grains (their centroid is not on any one ray); noise; a per-pair attribution radius from
the predicted origin tolerance; searches with the attributed cost (Tasks E, F); the cause of the 3,890 unpartnered
detector-0 spots; a real layer thickness (the simulated layer is exact, z = 0); whether a classifier (Step 4) would
raise the 0.854 precision or recover the several-grain pairs.
