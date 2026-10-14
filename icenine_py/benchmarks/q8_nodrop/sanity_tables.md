<!-- table:q8_sanity_cost -->
| variant | six_voxel_cost_median | old | cpp | f0_R2 | f0_A |
|---|---|---|---|---|---|
| clean | 0.0674 | 0.0643 | 0.0643 | 0.0643 | 0.0643 |
| realistic | 0.2586 | 0.2646 | 0.2831 | 0.2831 | 0.2825 |
| realistic_nodrop | 0.1168 | 0.1162 | 0.1162 | 0.1162 | 0.1162 |
<!-- /table:q8_sanity_cost -->

<!-- table:q8_sanity_rejected -->
| variant | setting | det0_overlaps | cpp | cpp_ci | f0_R2 | f0_R2_ci | f0_A | f0_A_ci |
|---|---|---|---|---|---|---|---|---|
| clean | screen | 524 | 3 of 524 (0.6%) | 0.2 to 1.7% | 3 of 524 (0.6%) | 0.2 to 1.7% | 3 of 524 (0.6%) | 0.2 to 1.7% |
| clean | local | 521 | 0 of 521 (0.0%) | 0.0 to 0.7% | 0 of 521 (0.0%) | 0.0 to 0.7% | 0 of 521 (0.0%) | 0.0 to 0.7% |
| realistic | screen | 475 | 43 of 475 (9.1%) | 6.8 to 12.0% | 43 of 475 (9.1%) | 6.8 to 12.0% | 39 of 475 (8.2%) | 6.1 to 11.0% |
| realistic | local | 471 | 41 of 471 (8.7%) | 6.5 to 11.6% | 40 of 471 (8.5%) | 6.3 to 11.4% | 36 of 471 (7.6%) | 5.6 to 10.4% |
| realistic_nodrop | screen | 524 | 3 of 524 (0.6%) | 0.2 to 1.7% | 3 of 524 (0.6%) | 0.2 to 1.7% | 3 of 524 (0.6%) | 0.2 to 1.7% |
| realistic_nodrop | local | 521 | 0 of 521 (0.0%) | 0.0 to 0.7% | 0 of 521 (0.0%) | 0.0 to 0.7% | 0 of 521 (0.0%) | 0.0 to 0.7% |
<!-- /table:q8_sanity_rejected -->
