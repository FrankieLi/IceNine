<!-- table:q8_e0_truth_cost -->
| msin | images | lenient | C++ rule | ray 2px | ray geometric |
|---|---|---|---|---|---|
| 0.3 | clean | 0.0682 | 0.0682 | 0.0682 | 0.0682 |
| 0.3 | realistic (10% drop) | 0.2627 | 0.2794 | 0.2793 | 0.2793 |
| 0.3 | realistic, no drop | 0.1148 | 0.1148 | 0.1148 | 0.1148 |
| 0.0 | clean | 0.0662 | 0.0662 | 0.0662 | 0.0662 |
| 0.0 | realistic (10% drop) | 0.2531 | 0.2681 | 0.2677 | 0.2677 |
| 0.0 | realistic, no drop | 0.1120 | 0.1120 | 0.1120 | 0.1120 |
<!-- /table:q8_e0_truth_cost -->

<!-- table:q8_e0_truth_rejected -->
| msin | images | det0_overlaps | C++ rule | C++ rule 95% interval | ray 2px | ray 2px 95% interval | ray geometric | ray geometric 95% interval |
|---|---|---|---|---|---|---|---|---|
| 0.3 | clean | 12439 | 0 of 12439 (0.0%) | 0.00 to 0.03% | 0 of 12439 (0.0%) | 0.00 to 0.03% | 0 of 12439 (0.0%) | 0.00 to 0.03% |
| 0.3 | realistic (10% drop) | 11174 | 974 of 11174 (8.7%) | 8.21 to 9.25% | 965 of 11174 (8.6%) | 8.13 to 9.17% | 964 of 11174 (8.6%) | 8.12 to 9.16% |
| 0.3 | realistic, no drop | 12439 | 0 of 12439 (0.0%) | 0.00 to 0.03% | 0 of 12439 (0.0%) | 0.00 to 0.03% | 0 of 12439 (0.0%) | 0.00 to 0.03% |
| 0 | clean | 12932 | 0 of 12932 (0.0%) | 0.00 to 0.03% | 0 of 12932 (0.0%) | 0.00 to 0.03% | 0 of 12932 (0.0%) | 0.00 to 0.03% |
| 0 | realistic (10% drop) | 11667 | 975 of 11667 (8.4%) | 7.87 to 8.87% | 966 of 11667 (8.3%) | 7.79 to 8.79% | 965 of 11667 (8.3%) | 7.79 to 8.78% |
| 0 | realistic, no drop | 12932 | 0 of 12932 (0.0%) | 0.00 to 0.03% | 0 of 12932 (0.0%) | 0.00 to 0.03% | 0 of 12932 (0.0%) | 0.00 to 0.03% |
<!-- /table:q8_e0_truth_rejected -->
