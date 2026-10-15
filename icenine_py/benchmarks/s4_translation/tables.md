<!-- table:s4_equality -->
| set | table | pairs | peak_overlap_differs | pass_differs | peak_coordinate_flips | unexplained | largest_boundary_distance |
|---|---|---|---|---|---|---|---|
| Phase D (20 voxels) | float64 | 879480 | 56 | 15 | 10459 | 0 | 2.4e-04 |
| Phase D (20 voxels) | float32 | 879480 | 35 | 7 | 6662 | 0 | 2.4e-04 |
| ThreeVoxels (3 voxels) | float64 | 131922 | 0 | 0 | 1625 | 0 | n/a |
| ThreeVoxels (3 voxels) | float32 | 131922 | 0 | 0 | 1044 | 0 | n/a |
<!-- /table:s4_equality -->

<!-- table:s4_full_layer -->
| item | value |
|---|---|
| ordering | O2 counts, u8, 10 threads, voxel tiles of 4096 |
| voxels x orientations | 24570 x 43974 |
| wall time (s), kernel + pass reduction | 6.58 (6.00 + 0.58) |
| BatchedScreen for the same work, single worker (h) | 3.94 |
| passers (voxel, orientation) pairs | 848266117 (78.5%) |
| passers per voxel, median (min-max) | 34466 (32991-35912) |
| passer list as int32 pairs / as a bitset / counts as uint8 (GB) | 6.79 / 0.135 / 1.08 |
| peak RSS of the run (GB) | 1.99 |
| label | contended |
<!-- /table:s4_full_layer -->

<!-- table:s4_linearity -->
| arm | s_500 | s_1000 | s_2000 | us_per_voxel_500 | us_per_voxel_1000 | us_per_voxel_2000 |
|---|---|---|---|---|---|---|
| O1 pass packed | 1.99 | 3.98 | 7.97 | 3.98e+03 | 3.98e+03 | 3.99e+03 |
| O2 counts u8 b=32FZ | 0.661 | 1.26 | 2.42 | 1.32e+03 | 1.26e+03 | 1.21e+03 |
| O3 counts u8 tile=256 | 0.586 | 1.18 | 2.38 | 1.17e+03 | 1.18e+03 | 1.19e+03 |
| O5 counts u8 tile=256 | 0.721 | 1.51 | 3.13 | 1.44e+03 | 1.51e+03 | 1.57e+03 |
| O2 counts u8 b=32FZ th=10 | 0.113 | 0.215 | 0.435 | 227 | 215 | 218 |
<!-- /table:s4_linearity -->

<!-- table:s4_o4 -->
| quantity | model | measured |
|---|---|---|
| bright rate per O2 test (patch, sampled peaks) | 0.0654 | n/a |
| voxels per distinct pixel hit | 1.59 | n/a |
| bright pixels visited per peak | 86.4 | 225 |
| candidate voxels tested per bright pixel | 1.59 | 6.97 |
| cost per bright pixel (ns) | 5 / 10 / 20 (assumed) | 33.2 |
| O2 time / O4 time, 1 thread | 5.06 / 2.57 / 1.30 | 0.3 |
| O4 / O2 wall, 1 thread (s) | n/a | 7.97 / 2.39 |
| O4 / O2 wall, 10 threads (s) | n/a | 1.13 / 0.41 |
<!-- /table:s4_o4 -->

<!-- table:s4_orderings -->
| arm | threads | median_s | ns_per_test | pairs_per_s | full_layer_s | equal |
|---|---|---|---|---|---|---|
| O1 counts u8 | 1 | 19 | 8.91 | 4.62e+06 | 234 | True |
| O1 pass u8 | 1 | 10.3 | 4.84 | 8.51e+06 | 127 | True |
| O2 counts u8 b=32FZ | 1 | 2.42 | 1.13 | 3.64e+07 | 29.7 | True |
| O2 pass u8 b=32FZ | 1 | 3.18 | 1.49 | 2.77e+07 | 39.1 | True |
| O3 counts u8 tile=256 | 1 | 2.38 | 1.12 | 3.69e+07 | 29.3 | True |
| O3 counts u8 tile=2000 | 1 | 2.17 | 1.02 | 4.05e+07 | 26.7 | True |
| O5 counts u8 tile=64 | 1 | 3.4 | 1.59 | 2.59e+07 | 41.8 | True |
| O5 counts u8 tile=256 | 1 | 3.13 | 1.47 | 2.81e+07 | 38.5 | True |
| O5 counts u8 tile=1024 | 1 | 2.61 | 1.22 | 3.37e+07 | 32.1 | True |
| O1 counts packed | 1 | 15.2 | 7.09 | 5.8e+06 | 186 | True |
| O1 pass packed | 1 | 7.97 | 3.73 | 1.1e+07 | 97.9 | True |
| O2 counts packed b=32FZ | 1 | 2.4 | 1.12 | 3.67e+07 | 29.5 | True |
| O2 pass packed b=32FZ | 1 | 3.32 | 1.55 | 2.65e+07 | 40.8 | True |
| O3 counts packed tile=256 | 1 | 2.54 | 1.19 | 3.46e+07 | 31.2 | True |
| O3 counts packed tile=2000 | 1 | 2.37 | 1.11 | 3.72e+07 | 29.1 | True |
| O5 counts packed tile=64 | 1 | 2.62 | 1.23 | 3.36e+07 | 32.2 | True |
| O5 counts packed tile=256 | 1 | 2.51 | 1.17 | 3.5e+07 | 30.8 | True |
| O5 counts packed tile=1024 | 1 | 2.41 | 1.13 | 3.65e+07 | 29.6 | True |
| O2 counts u8 b=1FZ | 1 | 2.43 | 1.14 | 3.62e+07 | 29.8 | True |
| O2 counts u8 b=256FZ | 1 | 2.42 | 1.13 | 3.64e+07 | 29.7 | True |
<!-- /table:s4_orderings -->

<!-- table:s4_projection -->
| scenario | seeds | total_h | share_screen | share_recost | share_refine | share_neighbours |
|---|---|---|---|---|---|---|
| (a) today's BFS | 350 | 7.65 | 0.007 | 0.109 | 0.819 | 0.064 |
| (b) global screen and re-cost, per-voxel refinement | 350 | 447 | 0.000 | 0.015 | 0.985 | 0.000 |
| (c) hybrid: BFS, seeds take the global level-0 list | 350 | 6.86 | 0.001 | 0.014 | 0.913 | 0.072 |
| (a) today's BFS | 395 | 8.57 | 0.007 | 0.110 | 0.826 | 0.057 |
| (b) global screen and re-cost, per-voxel refinement | 395 | 447 | 0.000 | 0.015 | 0.985 | 0.000 |
| (c) hybrid: BFS, seeds take the global level-0 list | 395 | 7.68 | 0.001 | 0.014 | 0.921 | 0.064 |
| (a) today's BFS | 504 | 10.8 | 0.007 | 0.111 | 0.836 | 0.045 |
| (b) global screen and re-cost, per-voxel refinement | 504 | 447 | 0.000 | 0.015 | 0.985 | 0.000 |
| (c) hybrid: BFS, seeds take the global level-0 list | 504 | 9.66 | 0.001 | 0.014 | 0.934 | 0.051 |
<!-- /table:s4_projection -->

<!-- table:s4_q1 -->
| check | n | result | exceptions |
|---|---|---|---|
| 1a float64: u(x + dx) - (u(x) + A dx), largest | 1760 | 1.8e-12 px | none |
| 1a float32 (evaluate's stage-D input) against the float64 law, largest | 1760 | 8.5e-04 px | float32 alone at x: 8.0e-04 px |
| 1a peak sets (omegas, normals, frames, plane hits) identical | 1760 | 1760 of 1760 moves | none |
| 1b simulator peak sets identical (count and frames) | 1440 | 1440 of 1440 moves | none |
| 1b(i) integer shift: S(x') = S(x) + m (float32 simulator) | 480 | 470 equal, 0 failures | 10 at the detector edge (3 off it), 0 float flips |
| 1b(i) the same in float64 vertices | 480 | 470 equal, 0 failures | 10 at the detector edge, 0 float flips |
| 1b(ii) phase re-rasterised: S(x') = R(T + phi) + m (float32 simulator) | 201240 | 116531 equal, 0 failures | 84644 at the detector edge (84439 off it), 65 float flips (largest distance 2.3e-04 px) |
| 1b(ii) the same in float64 vertices | 201240 | 116596 equal, 0 failures | 84644 at the detector edge, 0 float flips |
| 1b(iii) pure pixel-set translation S(x) + round(A dx) = S(x') | 114695 | 30746 equal (26.8%, Wilson 0.266-0.271) | 86545 edge spots excluded |
| 1e dilated bitmap and box bounds = pixel_radius_overlap_uint8 | 33882912 | equal | none (4 frames x 2 detectors, border band included) |
<!-- /table:s4_q1 -->

<!-- table:s4_q1_cost -->
| measure | value |
|---|---|
| passers re-costed, exact and translated (20 voxels) | 691539 |
| passers with any OverlapInfo field changed | 690613 (99.9%) |
| passers with the re-cost (1 - confidence) changed | 203241 (29.4%) |
| largest re-cost change | 0.082 |
| largest quality change (final cost = 1 - quality) | 0.042 |
| cliques whose spacing-filtered set is unchanged | 69881 of 97720 (71.5%) |
| voxels whose best filtered candidate is unchanged | 13 of 20 (Wilson 0.433-0.819) |
| top-30 filtered candidates in common, median over voxels | 19.5 of 30 |
<!-- /table:s4_q1_cost -->

<!-- table:s4_q1_pure -->
| stat | value |
|---|---|
| spots compared (both on the detector, unclipped) | 114695 |
| exactly equal | 30746 |
| symmetric difference, median / 95th pct / max (px) | 3 / 12 / 17 |
| symmetric difference / spot size, median / 95th pct | 0.30 / 1.57 |
| spots with every differing pixel within 1 px (8-neighbour) of the other set | 114695 |
| spots with every differing pixel on its set's 4-boundary | 114607 |
| spots leaving / entering the detector over the moves | 1676 / 1717 |
<!-- /table:s4_q1_pure -->

<!-- table:s4_recost -->
| item | value | n |
|---|---|---|
| today: evaluate per (orientation, voxel) (us) | 266 | 4000 |
| orientation-major: per pair, float64 vertices + stage-D C call (us) | 28.9 | 321456 |
|   the same with one vertex product per orientation (us) | 29.8 | 321456 |
|   set-up per orientation: peaks, A and c (us) | 458 | 200 |
| speed-up per pair | 9.2x | 321456 |
| equality: pairs compared with the Step 0 per-voxel evaluate | 691539 | 20 |
|   pairs with any OverlapInfo field different | 33962 (4.9%) | 691539 |
|   pairs with the re-cost (1 - confidence) different | 85 | 691539 |
|   largest quality difference | 0.0092 | 691539 |
|   sampled differing pairs explained by an integer-vertex flip | 400 of 400 (largest distance 4.1e-04 px) | 400 |
<!-- /table:s4_recost -->

<!-- table:s4_step0 -->
| part | value | unit | n |
|---|---|---|---|
| BatchedScreen, whole level-0 screen of one voxel (chunk 32) | 0.577 | s | 10 |
|   of which torch geometry (Bragg solve, eta, stages B-C) | 0.401 | s | 10 |
|   of which the per-orientation stage-D loop | 0.212 | s | 10 |
| screen passers per voxel (of 43,974), median | 3.45e+04 | orientations | 10 |
| re-cost of the passers (local evaluate), per voxel | 8.56 | s | 3 |
| re-cost, per evaluate | 251 | us | 3 |
| no-start search (reconstruct_voxel), wall | 73.6 | s | 3 |
|   screens of all levels | 1.56 | s | 3 |
|   re-costs of all levels | 17.4 | s | 3 |
|   quick MC, all levels (level 0) | 50.9 (34.1) | s | 3 |
|   rest of the discrete stage (spacing filter, loop) | 3.07 | s | 3 |
|   FindOptimal, VarianceMinimizing, final | 0.265 | s | 3 |
| BFS neighbour step (CMA, 250 evaluations) | 0.0732 | s | 10 |
<!-- /table:s4_step0 -->

<!-- table:s4_threads -->
| arm | 1 | 2 | 4 | 6 | 8 | 10 | speedup_10 |
|---|---|---|---|---|---|---|---|
| O2 counts u8 b=32FZ | 2.42 | 1.23 | 0.655 | 0.555 | 0.483 | 0.435 | 5.56 |
| O2 pass u8 b=32FZ | 3.18 | 1.61 | 0.896 | 0.715 | 0.604 | 0.527 | 6.03 |
| O1 pass packed | 7.97 | 4.31 | 2.39 | 2.83 | 1.93 | 1.78 | 4.47 |
| O3 counts u8 tile=256 | 2.38 | 1.26 | 0.71 | 0.761 | 0.521 | 0.536 | 4.44 |
| O5 counts u8 tile=256 | 3.13 | 1.82 | 1.07 | 0.9 | 0.783 | 0.692 | 4.53 |
<!-- /table:s4_threads -->

<!-- table:s4_variants -->
| arm | median_s | peaks_s | kernel_s | frames_changed | counts_differ | pass_differ |
|---|---|---|---|---|---|---|
| O2 streamed (no table), b=1FZ | 3.89 | 1.39 | 2.48 | 1 | some (patch) | n/a |
| O2 streamed (no table), b=32FZ | 2.78 | 0.346 | 2.43 | 0 | 0 | n/a |
| O2 streamed (no table), b=256FZ | 2.71 | 0.285 | 2.42 | 0 | 0 | n/a |
| O2 on MPS (float32 offsets), block=288 | 8.68 | n/a | n/a | n/a | 10 | 1 |
| O2 on MPS (float32 offsets), block=576 | 8.27 | n/a | n/a | n/a | 10 | 1 |
| O2 on MPS (float32 offsets), block=1152 | 25 | n/a | n/a | n/a | 10 | 1 |
<!-- /table:s4_variants -->
