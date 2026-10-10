<!-- table:phase_d_reseed_compare -->
| run | pilot | wrong_pilot | wrong_reseed | unres_pilot | unres_reseed | lost_pilot | lost_reseed | acc_wrong_pilot | acc_wrong_reseed | n_clusters_differing | p |
|---|---|---|---|---|---|---|---|---|---|---|---|
| mc_clean | mc_clean | 740 | 633 | 710 | 601 | 28 | 19 | 32 | 32 | 10 | 0.00195 |
| cma_noretry_clean | cma_noretry_clean | 734 | 633 | 677 | 576 | 28 | 19 | 57 | 57 | 10 | 0.00195 |
| mc_realistic_q16 | mc_realistic_q16 | 709 | 608 | 691 | 590 | 27 | 18 | 19 | 19 | 9 | 0.00391 |
| cma_noretry_realistic_q16 | cma_realistic_q16 | 617 | 568 | 580 | 533 | 26 | 19 | 37 | 35 | 12 | 0.533 |
<!-- /table:phase_d_reseed_compare -->

<!-- table:phase_d_reseed_counts -->
| run | seeds | seed_rej | reseeds | accepted | rejected | median_cluster | final_voxels | reseed_min | full_search_min | share_pct | ev_share_pct | unres_clusters | largest_unres | largest_grains | largest_lost_inside | largest_has_rej_seed | wall_h | wall_pilot_h | ev_reseed |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mc_clean | 34 | 2 | 11 | 11 | 0 | 6.0 | 0 | 95.7 | 87.3 | 44.6 | 21.8 | 1 | 601 | 20 | 19 | True | 3.57 | 2.27 | 3,723,137 |
| cma_noretry_clean | 33 | 1 | 10 | 10 | 0 | 7.0 | 0 | 92.9 | 92.7 | 46.6 | 22.1 | 1 | 576 | 19 | 19 | True | 3.32 | 2.06 | 3,403,584 |
| mc_realistic_q16 | 36 | 4 | 10 | 9 | 1 | 7.0 | 1 | 98.5 | 97.9 | 39.5 | 19.4 | 2 | 589 | 21 | 18 | True | 4.15 | 2.41 | 3,442,938 |
| cma_noretry_realistic_q16 | 43 | 10 | 10 | 9 | 1 | 7.5 | 4 | 96.9 | 96.7 | 31.8 | 18.2 | 2 | 529 | 18 | 18 | True | 5.09 | 2.40 | 3,638,977 |
<!-- /table:phase_d_reseed_counts -->

<!-- table:phase_d_reseed_log -->
| run | k | cluster | voxel | accepted | left_refit | seconds | evals |
|---|---|---|---|---|---|---|---|
| mc_clean | 1 | 21 | 1 | yes | 4 | 181 | 345,600 |
| mc_clean | 2 | 16 | 67 | yes | 6 | 815 | 339,345 |
| mc_clean | 3 | 6 | 227 | yes | 0 | 1273 | 340,846 |
| mc_clean | 4 | 59 | 854 | yes | 17 | 344 | 358,345 |
| mc_clean | 5 | 4 | 210 | yes | 0 | 1537 | 333,165 |
| mc_clean | 6 | 8 | 226 | yes | 1 | 699 | 334,393 |
| mc_clean | 7 | 4 | 238 | yes | 0 | 170 | 329,932 |
| mc_clean | 8 | 17 | 528 | yes | 3 | 181 | 341,184 |
| mc_clean | 9 | 3 | 961 | yes | 0 | 179 | 331,783 |
| mc_clean | 10 | 1 | 1888 | yes | 0 | 180 | 334,067 |
| mc_clean | 11 | 1 | 1938 | yes | 0 | 180 | 334,477 |
| cma_noretry_clean | 1 | 21 | 1 | yes | 4 | 176 | 330,144 |
| cma_noretry_clean | 2 | 16 | 67 | yes | 6 | 177 | 332,749 |
| cma_noretry_clean | 3 | 6 | 227 | yes | 0 | 191 | 360,732 |
| cma_noretry_clean | 4 | 55 | 203 | yes | 15 | 180 | 341,070 |
| cma_noretry_clean | 5 | 8 | 226 | yes | 1 | 174 | 330,645 |
| cma_noretry_clean | 6 | 4 | 238 | yes | 0 | 787 | 328,035 |
| cma_noretry_clean | 7 | 14 | 528 | yes | 0 | 945 | 332,111 |
| cma_noretry_clean | 8 | 1 | 861 | yes | 0 | 563 | 331,954 |
| cma_noretry_clean | 9 | 1 | 1888 | yes | 0 | 1666 | 358,407 |
| cma_noretry_clean | 10 | 1 | 1938 | yes | 0 | 714 | 357,737 |
| mc_realistic_q16 | 1 | 21 | 1 | yes | 4 | 574 | 349,205 |
| mc_realistic_q16 | 2 | 16 | 67 | yes | 6 | 1655 | 344,685 |
| mc_realistic_q16 | 3 | 6 | 227 | yes | 0 | 701 | 340,417 |
| mc_realistic_q16 | 4 | 56 | 203 | yes | 15 | 190 | 365,918 |
| mc_realistic_q16 | 5 | 8 | 226 | yes | 1 | 182 | 342,250 |
| mc_realistic_q16 | 6 | 4 | 238 | yes | 0 | 184 | 341,646 |
| mc_realistic_q16 | 7 | 14 | 528 | yes | 0 | 188 | 348,396 |
| mc_realistic_q16 | 8 | 1 | 861 | no | 1 | 532 | 336,641 |
| mc_realistic_q16 | 9 | 1 | 1888 | yes | 0 | 385 | 336,459 |
| mc_realistic_q16 | 10 | 1 | 1938 | yes | 0 | 1320 | 337,321 |
| cma_noretry_realistic_q16 | 1 | 21 | 1 | yes | 4 | 197 | 365,199 |
| cma_noretry_realistic_q16 | 2 | 16 | 67 | yes | 6 | 547 | 363,104 |
| cma_noretry_realistic_q16 | 3 | 6 | 227 | yes | 0 | 399 | 363,442 |
| cma_noretry_realistic_q16 | 4 | 55 | 203 | yes | 15 | 1338 | 371,351 |
| cma_noretry_realistic_q16 | 5 | 9 | 226 | yes | 1 | 2414 | 363,728 |
| cma_noretry_realistic_q16 | 6 | 4 | 238 | no | 4 | 182 | 362,593 |
| cma_noretry_realistic_q16 | 7 | 14 | 528 | yes | 0 | 184 | 366,173 |
| cma_noretry_realistic_q16 | 8 | 1 | 861 | yes | 0 | 182 | 360,639 |
| cma_noretry_realistic_q16 | 9 | 1 | 1888 | yes | 0 | 184 | 359,393 |
| cma_noretry_realistic_q16 | 10 | 1 | 1938 | yes | 0 | 187 | 363,355 |
<!-- /table:phase_d_reseed_log -->
