<!-- table:q8_sanity_cost -->
| variant | 6 voxels, lenient (sanity_checks) | lenient | C++ rule | ray 2px | ray geometric |
|---|---|---|---|---|---|
| clean | 0.0674 | 0.0643 | 0.0643 | 0.0643 | 0.0643 |
| realistic | 0.2586 | 0.2646 | 0.2831 | 0.2831 | 0.2825 |
| realistic_nodrop | 0.1168 | 0.1162 | 0.1162 | 0.1162 | 0.1162 |
<!-- /table:q8_sanity_cost -->

<!-- table:q8_sanity_rejected -->
| variant | setting | det0_overlaps | C++ rule | C++ rule 95% interval | ray 2px | ray 2px 95% interval | ray geometric | ray geometric 95% interval |
|---|---|---|---|---|---|---|---|---|
| clean | screen | 524 | 3 of 524 (0.6%) | 0.19 to 1.67% | 3 of 524 (0.6%) | 0.19 to 1.67% | 3 of 524 (0.6%) | 0.19 to 1.67% |
| clean | local | 521 | 0 of 521 (0.0%) | 0.00 to 0.73% | 0 of 521 (0.0%) | 0.00 to 0.73% | 0 of 521 (0.0%) | 0.00 to 0.73% |
| realistic | screen | 475 | 43 of 475 (9.1%) | 6.79 to 11.97% | 43 of 475 (9.1%) | 6.79 to 11.97% | 39 of 475 (8.2%) | 6.06 to 11.03% |
| realistic | local | 471 | 41 of 471 (8.7%) | 6.48 to 11.60% | 40 of 471 (8.5%) | 6.30 to 11.36% | 36 of 471 (7.6%) | 5.57 to 10.40% |
| realistic_nodrop | screen | 524 | 3 of 524 (0.6%) | 0.19 to 1.67% | 3 of 524 (0.6%) | 0.19 to 1.67% | 3 of 524 (0.6%) | 0.19 to 1.67% |
| realistic_nodrop | local | 521 | 0 of 521 (0.0%) | 0.00 to 0.73% | 0 of 521 (0.0%) | 0.00 to 0.73% | 0 of 521 (0.0%) | 0.00 to 0.73% |
<!-- /table:q8_sanity_rejected -->
