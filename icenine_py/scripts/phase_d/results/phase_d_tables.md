<!-- table:phase_d_cost_check -->
| images | voxels | q_true | q_0.25 | q_0.5 | q_1.0_max | min_gap_true_minus_1deg |
|---|---|---|---|---|---|---|
| clean | 6 | 0.90 - 0.95 | 0.46 - 0.62 | 0.18 - 0.34 | 0.06 - 0.13 | 0.80 |
| realistic | 6 | 0.69 - 0.81 | 0.42 - 0.51 | 0.16 - 0.29 | 0.05 - 0.13 | 0.60 |
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

<!-- table:phase_d_render -->
| run | voxels | workers | wall_s | max_worker_sim_s | max_worker_write_noise_s | max_worker_rss_gb |
|---|---|---|---|---|---|---|
| pilot (500 voxels) | 500 | 10 | 13.8 | 0.8 | 10.3 | 1.12 |
| full | 24570 | 10 | 47.5 | 35.6 | 9.8 | 1.49 |
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
