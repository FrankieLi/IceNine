<!-- table:evalE_timing_ms -->
| cost | full-sample: ms per evaluation, batch of 256 | full-sample: ms per evaluation, one per call | full-sample: x strict | E0: ms per evaluation, batch of 256 | E0: ms per evaluation, one per call | E0: x strict |
|---|---|---|---|---|---|---|
| strict cost (one call per orientation) | 0.381 | 0.381 | 1.000 | 0.392 | 0.392 | 1.000 |
| S | 0.470 | 0.809 | 1.233 | 0.228 | 0.601 | 0.581 |
| E1 | 1.332 | 1.649 | 3.492 | 0.240 | 0.634 | 0.613 |
| E2 | 2.080 | 2.675 | 5.453 | 0.886 | 1.590 | 2.261 |
| E3 | 0.473 | 0.844 | 1.240 | 0.227 | 0.640 | 0.580 |
| E4 | 1.186 | 2.153 | 3.109 | 0.313 | 1.248 | 0.798 |
| E1+E2 | 6.663 | 7.268 | 17.466 | 0.950 | 1.710 | 2.424 |
| E1+E2+E3 | 6.728 | 7.485 | 17.637 | 0.952 | 1.747 | 2.430 |
| E2+E3 | 2.136 | 2.792 | 5.599 | 0.888 | 1.627 | 2.265 |
| E3+E4 | 1.639 | 2.811 | 4.297 | 0.495 | 1.632 | 1.262 |
| E1+E2+E3+E4 | 7.903 | 9.499 | 20.718 | 1.223 | 2.745 | 3.120 |
<!-- /table:evalE_timing_ms -->

<!-- table:evalE_timing_per_voxel -->
| voxel set | arm | modelled strict cost time per voxel, baseline (s) | modelled strict cost time per voxel, arm (s) | smooth evaluations per voxel | modelled smooth cost time per voxel (s) | SmoothCost build (s) | extra time as shipped (a build per voxel) | projection: one build per image set (cache not built) |
|---|---|---|---|---|---|---|---|---|
| six full-sample | E3+E4,tau=0.25 | 126.7 | 127.6 | 15030 | 24.6 | 16.6 | +33.3% | +20.2% |
| six full-sample | E1+E2+E3+E4 | 126.7 | 132.8 | 15755 | 124.5 | 16.6 | +116.2% | +103.1% |
<!-- /table:evalE_timing_per_voxel -->

<!-- table:evalE_timing_voxel_wall -->
| voxel set | voxel | timing | key off (s) | E3+E4,tau=0.25 (s) | E3+E4,tau=0.25: extra time | E1+E2+E3+E4 (s) | E1+E2+E3+E4: extra time | same error as the end-to-end run |
|---|---|---|---|---|---|---|---|---|
| E0 | 7336 | single-worker, median of 2 | 2.1 | 6.0 | +179% | 6.4 | +201% | yes |
| E0 | 3108 | single-worker, median of 2 | 2.4 | 6.5 | +173% | 7.1 | +195% | yes |
| E0 | 2910 | single-worker, median of 2 | 3.0 | 7.0 | +130% | 7.7 | +153% | yes |
| E0 | 16182 | single-worker, median of 2 | 2.5 | 6.5 | +164% | 6.9 | +179% | yes |
| full-sample | 15901 | contended, median of 1 | 117.4 | 159.6 | +36% | 263.9 | +125% | yes |
<!-- /table:evalE_timing_voxel_wall -->
