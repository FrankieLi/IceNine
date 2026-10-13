<!-- table:phase_d_cost_check -->
| images | voxels | q_true | q_0.25 | q_0.5 | q_1.0_max | min_gap_true_minus_1deg |
|---|---|---|---|---|---|---|
| clean | 6 | 0.90 - 0.95 | 0.46 - 0.62 | 0.18 - 0.34 | 0.06 - 0.13 | 0.80 |
| realistic | 6 | 0.69 - 0.81 | 0.42 - 0.51 | 0.16 - 0.29 | 0.05 - 0.13 | 0.60 |
| realistic_q16 | 6 | 0.69 - 0.77 | 0.38 - 0.51 | 0.16 - 0.28 | 0.05 - 0.12 | 0.57 |
<!-- /table:phase_d_cost_check -->

<!-- table:phase_d_frames -->
| quantity | value |
|---|---|
| frames (180 omega x 2 detectors) | 360 |
| lit pixels per frame, median | 20484.5 |
| lit pixels per frame, min | 13413 |
| lit pixels per frame, max | 26295 |
| lit pixels, clean total | 7414841 |
| lit pixels, realistic total | 6592953 |
| spots (connected components), all frames | 57389 |
| spots missed | 5743 |
| hot pixels added | 2293 |
| blobs added | 5148 |
<!-- /table:phase_d_frames -->

<!-- table:phase_d_memmap -->
| images | evaluations | identical | dense_rss_gb | memmap_rss_gb | dense_load_s | memmap_load_s |
|---|---|---|---|---|---|---|
| clean | 40 | True | 7.88 | 0.37 | 8.2 | 0.4 |
| realistic | 40 | True | 7.67 | 0.37 | 7.7 | 0.4 |
| realistic_q16 | 40 | True | 7.04 | 0.36 | 57.5 | 0.4 |
<!-- /table:phase_d_memmap -->

<!-- table:phase_d_probe -->
| quantity | value |
|---|---|
| no-start seed voxel (reconstruct_voxel), n | 4 |
|   wall time per seed (s) | 177, 376, 711, 204 |
|   final error (deg) | 0.005, 0.013, 0.012, 0.005 |
| neighbour step (local_optimization from 0.3 deg off), n | 20 |
|   median / mean wall time (s) | 0.4 / 0.4 |
|   median final error (deg) | 0.300 |
<!-- /table:phase_d_probe -->

<!-- table:phase_d_q16 -->
| quantity | value |
|---|---|
| render wall time, 10 workers (s) | 245.7 |
| max worker sim time (s) | 232.3 |
| max worker RSS (GB) | 2.25 |
| lit pixels, clean total | 58842023 |
| lit pixels, realistic total | 51498581 |
| lit fraction of all pixels (Q-max 8: 0.00491) | 0.039 |
| spots, all frames | 220381 |
| spot size px: median / p95 / max | 190 / 753 / 5194 |
| spots larger than 2000 px | 277 |
<!-- /table:phase_d_q16 -->

<!-- table:phase_d_render -->
| run | voxels | workers | wall_s | max_worker_sim_s | max_worker_write_noise_s | max_worker_rss_gb |
|---|---|---|---|---|---|---|
| pilot (500 voxels) | 500 | 10 | 13.8 | 0.8 | 10.3 | 1.12 |
| full | 24570 | 10 | 48.7 | 35.8 | 11.0 | 1.39 |
| full Q-max 16 | 24570 | 10 | 245.7 | 232.3 | 11.4 | 2.25 |
<!-- /table:phase_d_render -->

<!-- table:phase_d_sample -->
| quantity | value |
|---|---|
| voxels | 24570 |
| grains (distinct orientations) | 497 |
| voxels per grain min / median / max | 4 / 44 / 185 |
| adjacent grain pairs | 1469 |
| redrawn grains (of 497) | 5 |
| min new vs any old (deg) | 1.29 |
| min new vs adjacent new (deg) | 3.9 |
| min old vs old (deg) | 7.25 |
| mic round-trip max error (deg) | 2.96e-06 |
<!-- /table:phase_d_sample -->

<!-- table:phase_d_sizes -->
| quantity | value |
|---|---|
| spot size px: median / p95 / max (Q-max 8) | 103 / 308.6 / 1072 |
| lit fraction of all pixels (Q-max 8) | 0.00491 |
<!-- /table:phase_d_sizes -->
