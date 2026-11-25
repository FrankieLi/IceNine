<!-- table:tex1_cost_wall -->
| arm | workers | evals | wall_s | wall_h |
|---|---|---|---|---|
| CSP | 12 workers, 13 ranks on 12 cores | 4,025,337,766 | 32,052 | 8.90 |
| RCP | 1 worker | 284,235,185 | 22,573 | 6.27 |
<!-- /table:tex1_cost_wall -->

<!-- table:tex1_noise_grains -->
| arm | cls | grains | recovered | pooled_median | wrong_in_recovered | size_recovered | size_other |
|---|---|---|---|---|---|---|---|
| CSP | noise | 24 | 0 (0%, 14%) | n/a | none recovered | none | [65, 1705] |
| CSP | both | 26 | 0 (0%, 13%) | n/a | none recovered | none | [108, 1427] |
| RCP | noise | 24 | 11 (28%, 65%) | 0.053 | 354 of 1975 (17.9%) | [65, 270] | [249, 1705] |
| RCP | both | 26 | 25 (81%, 99%) | 0.029 | 1339 of 11142 (12.0%) | [108, 1427] | [522, 522] |
<!-- /table:tex1_noise_grains -->

<!-- table:tex1_pairs -->
| subset | both_wrong | only_a_wrong | only_b_wrong | neither |
|---|---|---|---|---|
| all_voxels | 2,482 | 9,078 | 1,510 | 85,210 |
| sharp | 8 | 1,043 | 7 | 36,355 |
| noise | 1,247 | 2,245 | 883 | 5,473 |
| gradient | 305 | 2,074 | 86 | 36,890 |
| both | 922 | 3,716 | 534 | 6,492 |
<!-- /table:tex1_pairs -->

<!-- table:tex1_per_class -->
| arm | cls | voxels | wrong | share | no_output | median | p90 | trivial_median | trivial_wrong | closer |
|---|---|---|---|---|---|---|---|---|---|---|
| CSP | all | 98,280 | 11,560 | 11.762% (11.562, 11.965) | 0 | 0.0744 | 1.1153 | 0.6859 | 35,520 | n/a (pools sharp voxels) |
| CSP | sharp | 37,413 | 1,051 | 2.809% (2.647, 2.982) | 0 | 0.0229 | 0.0546 | 0.0000 | 0 | 0.0% |
| CSP | noise | 9,848 | 3,492 | 35.459% (34.520, 36.409) | 0 | 0.8353 | 1.5108 | 0.7947 | 0 | 42.9% |
| CSP | gradient | 39,355 | 2,379 | 6.045% (5.814, 6.285) | 0 | 0.1896 | 0.4327 | 1.4813 | 26,558 | 92.5% |
| CSP | both | 11,664 | 4,638 | 39.763% (38.879, 40.655) | 0 | 0.8254 | 1.6459 | 1.6647 | 8,962 | 76.9% |
| RCP | all | 98,280 | 3,992 | 4.062% (3.940, 4.187) | 0 | 0.0303 | 0.3923 | 0.6859 | 35,520 | n/a (pools sharp voxels) |
| RCP | sharp | 37,413 | 15 | 0.040% (0.024, 0.066) | 0 | 0.0187 | 0.0410 | 0.0000 | 0 | 0.0% |
| RCP | noise | 9,848 | 2,130 | 21.629% (20.827, 22.453) | 0 | 0.6580 | 1.2002 | 0.7947 | 0 | 58.4% |
| RCP | gradient | 39,355 | 391 | 0.994% (0.900, 1.096) | 0 | 0.0570 | 0.2418 | 1.4813 | 26,558 | 98.4% |
| RCP | both | 11,664 | 1,456 | 12.483% (11.895, 13.095) | 0 | 0.0300 | 1.1147 | 1.6647 | 8,962 | 91.3% |
<!-- /table:tex1_per_class -->

<!-- table:tex1_spread_ratio -->
| arm | grains | median | q25 | q75 | within |
|---|---|---|---|---|---|
| CSP | 125 | 1.134 | 1.060 | 1.501 | 87 of 125 |
| RCP | 125 | 1.003 | 0.988 | 1.023 | 125 of 125 |
<!-- /table:tex1_spread_ratio -->
