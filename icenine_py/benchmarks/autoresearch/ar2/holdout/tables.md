<!-- table:ar2a_arms -->
| vset | arm | wrong | mean_evals | median_evals | median_err |
|---|---|---|---|---|---|
| boundary | best 1 | 5 of 150 (1.4 to 7.6%) | 1107 | 1121 | 0.0154 |
| boundary | best 2 | 1 of 150 (0.1 to 3.7%) | 1978 | 1931 | 0.0146 |
| boundary | best 3 | 1 of 150 (0.1 to 3.7%) | 2581 | 2551 | 0.0146 |
| interior | best 1 | 0 of 150 (0.0 to 2.5%) | 1497 | 1424 | 0.0255 |
| interior | best 2 | 0 of 150 (0.0 to 2.5%) | 2114 | 2028 | 0.0255 |
| interior | best 3 | 0 of 150 (0.0 to 2.5%) | 2687 | 2568 | 0.0255 |
| both | best 1 | 5 of 300 (0.7 to 3.8%) | 1302 | 1226 | 0.0191 |
| both | best 2 | 1 of 300 (0.1 to 1.9%) | 2046 | 1984 | 0.0189 |
| both | best 3 | 1 of 300 (0.1 to 1.9%) | 2634 | 2562 | 0.0189 |
<!-- /table:ar2a_arms -->

<!-- table:ar2a_explained -->
| bin | n | at1 | at5 | at20 |
|---|---|---|---|---|
| 35 or fewer | 21 | 8 of 21 | 21 of 21 | 21 of 21 |
| 36 to 40 | 71 | 33 of 71 | 71 of 71 | 71 of 71 |
| 41 to 45 | 117 | 87 of 117 | 117 of 117 | 117 of 117 |
| 46 to 50 | 78 | 68 of 78 | 78 of 78 | 78 of 78 |
| 51 or more | 13 | 11 of 13 | 13 of 13 | 13 of 13 |
<!-- /table:ar2a_explained -->

<!-- table:ar2a_recall -->
| vset | tol | at1 | at5 | at20 |
|---|---|---|---|---|
| boundary | 0.5 degrees | 78 of 150 (44.1 to 59.8%) | 150 of 150 (97.5 to 100.0%) | 150 of 150 (97.5 to 100.0%) |
| boundary | 2 degrees | 78 of 150 (44.1 to 59.8%) | 150 of 150 (97.5 to 100.0%) | 150 of 150 (97.5 to 100.0%) |
| boundary | 5 degrees | 78 of 150 (44.1 to 59.8%) | 150 of 150 (97.5 to 100.0%) | 150 of 150 (97.5 to 100.0%) |
| interior | 0.5 degrees | 129 of 150 (79.5 to 90.7%) | 150 of 150 (97.5 to 100.0%) | 150 of 150 (97.5 to 100.0%) |
| interior | 2 degrees | 129 of 150 (79.5 to 90.7%) | 150 of 150 (97.5 to 100.0%) | 150 of 150 (97.5 to 100.0%) |
| interior | 5 degrees | 129 of 150 (79.5 to 90.7%) | 150 of 150 (97.5 to 100.0%) | 150 of 150 (97.5 to 100.0%) |
| both | 0.5 degrees | 207 of 300 (63.6 to 74.0%) | 300 of 300 (98.7 to 100.0%) | 300 of 300 (98.7 to 100.0%) |
| both | 2 degrees | 207 of 300 (63.6 to 74.0%) | 300 of 300 (98.7 to 100.0%) | 300 of 300 (98.7 to 100.0%) |
| both | 5 degrees | 207 of 300 (63.6 to 74.0%) | 300 of 300 (98.7 to 100.0%) | 300 of 300 (98.7 to 100.0%) |
<!-- /table:ar2a_recall -->

<!-- table:ar2a_time -->
| vset | n | median_s | mean_s | p90_s | max_s |
|---|---|---|---|---|---|
| boundary | 150 | 0.19 | 0.19 | 0.25 | 0.32 |
| interior | 150 | 0.15 | 0.15 | 0.19 | 0.26 |
| both | 300 | 0.17 | 0.17 | 0.23 | 0.32 |
<!-- /table:ar2a_time -->

<!-- table:ar2a_trace -->
| voxel | vset | wrong_in_best | truth_rank | truth_in_first5 | truth_refined_in_best | taken | taken_err | taken_cost | truth_cost | truth_err |
|---|---|---|---|---|---|---|---|---|---|---|
| 1687 | boundary | 1 | 1 | yes | 2 | grain 110, touches | 53.166 | 0.137 to 0.137 | 0.165 to 0.061 | 0.009 |
| 8343 | boundary | 1 | 1 | yes | 2 | grain 379, touches | 29.542 | 0.184 to 0.184 | 0.321 to 0.087 | 0.013 |
| 11568 | boundary | 1 | 1 | yes | 2 | grain 257, touches | 22.196 | 0.218 to 0.141 | 0.223 to 0.122 | 0.013 |
| 17361 | boundary | 1,2,3 | 1 | yes | 2 | grain 254, touches | 57.305 | 0.133 to 0.133 | 0.162 to 0.162 | 0.044 |
| 21535 | boundary | 1 | 2 | yes | 2 | grain 165, touches | 47.622 | 0.299 to 0.279 | 0.346 to 0.081 | 0.000 |
<!-- /table:ar2a_trace -->
