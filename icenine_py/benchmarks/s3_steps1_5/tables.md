<!-- table:s3_ab_evals -->
| arm | candidates per level, median | local evaluations per run, mean | global evaluations per run, mean | quick-MC evaluations per run, mean | run time s, median (10 workers, contended) |
|---|---|---|---|---|---|
| key off | 178/47/12/3 | 4125 | 44500 | 2896 | 3.8 |
| key on (CPU) | 178/47/12/3 | 4120 | 44501 | 2896 | 2.3 |
| key on (MPS; seed 0, voxels 0-49) | 176/47/12/3 | 4121 | 44499 | 2891 | n/a |
<!-- /table:s3_ab_evals -->

<!-- table:s3_ab_paired -->
| pair (a vs b) | scope | paired runs | wrong a | wrong b | fixed (wrong a, right b) | broken (right a, wrong b) | both wrong | both right |
|---|---|---|---|---|---|---|---|---|
| off vs on_cpu | all seeds | 600 | 111 | 107 | 60 | 56 | 51 | 433 |
| off vs on_cpu | seed 0 | 200 | 35 | 35 | 18 | 18 | 17 | 147 |
| off vs on_cpu | seed 1 | 200 | 37 | 40 | 15 | 18 | 22 | 145 |
| off vs on_cpu | seed 2 | 200 | 39 | 32 | 27 | 20 | 12 | 141 |
| on_cpu vs on_mps | all seeds | 50 | 8 | 8 | 0 | 0 | 8 | 42 |
| on_cpu vs on_mps | seed 0 | 50 | 8 | 8 | 0 | 0 | 8 | 42 |
| off_2w vs on_mps | all seeds | 50 | 12 | 8 | 8 | 4 | 4 | 34 |
| off_2w vs on_mps | seed 0 | 50 | 12 | 8 | 8 | 4 | 4 | 34 |
| on_cpu_2w vs on_mps | all seeds | 50 | 8 | 8 | 0 | 0 | 8 | 42 |
| on_cpu_2w vs on_mps | seed 0 | 50 | 8 | 8 | 0 | 0 | 8 | 42 |
| key off, seed 0 vs seed 1 | same voxels | 200 | 35 | 37 | 20 | 22 | 15 | 143 |
| key off, seed 0 vs seed 2 | same voxels | 200 | 35 | 39 | 14 | 18 | 21 | 147 |
| key off, seed 1 vs seed 2 | same voxels | 200 | 37 | 39 | 21 | 23 | 16 | 140 |
<!-- /table:s3_ab_paired -->

<!-- table:s3_ab_timing_2w -->
| arm | runs | wrong (> 1 deg) | run time s, median | run time s, mean | run time s, sum |
|---|---|---|---|---|---|
| key off (2 workers) | 50 | 12 | 2.16 | 2.14 | 106.8 |
| key on (CPU, 2 workers) | 50 | 8 | 1.28 | 1.28 | 64.2 |
| key on (MPS; seed 0, voxels 0-49) | 50 | 8 | 1.45 | 1.48 | 74.0 |
<!-- /table:s3_ab_timing_2w -->

<!-- table:s3_ab_wrong -->
| arm | scope | wrong (> 1 deg) | runs | rate | Wilson 95% low | Wilson 95% high |
|---|---|---|---|---|---|---|
| key off | seed 0 | 35 | 200 | 0.175 | 0.129 | 0.234 |
| key off | seed 1 | 37 | 200 | 0.185 | 0.137 | 0.245 |
| key off | seed 2 | 39 | 200 | 0.195 | 0.146 | 0.255 |
| key off | pooled | 111 | 600 | 0.185 | 0.156 | 0.218 |
| key on (CPU) | seed 0 | 35 | 200 | 0.175 | 0.129 | 0.234 |
| key on (CPU) | seed 1 | 40 | 200 | 0.200 | 0.150 | 0.261 |
| key on (CPU) | seed 2 | 32 | 200 | 0.160 | 0.116 | 0.217 |
| key on (CPU) | pooled | 107 | 600 | 0.178 | 0.150 | 0.211 |
| key on (MPS; seed 0, voxels 0-49) | seed 0 | 8 | 50 | 0.160 | 0.083 | 0.285 |
| key on (MPS; seed 0, voxels 0-49) | pooled | 8 | 50 | 0.160 | 0.083 | 0.285 |
<!-- /table:s3_ab_wrong -->

<!-- table:s3_crossover -->
| batch B | CPU compact us per orientation | CPU fixed-shape us per orientation | MPS fixed-shape us per orientation |
|---|---|---|---|
| 1 | 258 | 290 | 3415 |
| 2 | 159 | 198 | 1560 |
| 4 | 110 | 145 | 831 |
| 8 | 83 | 119 | 420 |
| 16 | 70 | 110 | 220 |
| 32 | 62 | 103 | 173 |
| 48 | 60 | 100 | 111 |
| 64 | 59 | 98 | 108 |
| 96 | 57 | 97 | 41 |
| 128 | 56 | 96 | 40 |
| 192 | 56 | 97 | 35 |
| 256 | 56 | 96 | 32 |
| 512 | 55 | 96 | 30 |
<!-- /table:s3_crossover -->

<!-- table:s3_timing_e0 -->
| arm (off = sequential quick MC) | searches | wall s, median | quick MC s, median | quick MC us per evaluation, median | candidates per level, median | evaluations per search, median | wall speedup (sum of searches) | quick-MC us per evaluation speedup |
|---|---|---|---|---|---|---|---|---|
| off | 20 | 1.97 | 0.98 | 359 | 168/44/12/2 | 48362 | n/a | n/a |
| cpu | 20 | 1.19 | 0.21 | 76 | 168/44/10/2 | 48312 | 1.69 | 4.71 |
| mps | 20 | 1.16 | 0.18 | 67 | 168/44/10/2 | 48312 | 1.73 | 5.33 |
<!-- /table:s3_timing_e0 -->

<!-- table:s3_timing_three -->
| arm (off = sequential quick MC) | searches | wall s, median | quick MC s, median | quick MC us per evaluation, median | candidates per level, median | evaluations per search, median | wall speedup (sum of searches) | quick-MC us per evaluation speedup |
|---|---|---|---|---|---|---|---|---|
| off | 6 | 3.78 | 2.24 | 348 | 397/105/27/6 | 53853 | n/a | n/a |
| cpu | 6 | 2.01 | 0.49 | 77 | 397/103/26/6 | 53924 | 1.89 | 4.50 |
| mps | 6 | 1.94 | 0.41 | 65 | 397/103/26/6 | 53924 | 1.95 | 5.36 |
<!-- /table:s3_timing_three -->
