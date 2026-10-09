<!-- table:phase_d_pilot_axis -->
| run | subset | n | wrong | rate | median_right |
|---|---|---|---|---|---|
| mc_clean | near_axis | 416 | 164 | 39.42 | 0.0119 |
| mc_clean | far_axis | 1584 | 576 | 36.36 | 0.0122 |
| cma_clean | near_axis | 416 | 164 | 39.42 | 0.0062 |
| cma_clean | far_axis | 1584 | 578 | 36.49 | 0.0059 |
| cma_noretry_clean | near_axis | 416 | 164 | 39.42 | 0.0055 |
| cma_noretry_clean | far_axis | 1584 | 570 | 35.98 | 0.0053 |
| mc_realistic_q16 | near_axis | 416 | 162 | 38.94 | 0.0026 |
| mc_realistic_q16 | far_axis | 1584 | 547 | 34.53 | 0.0167 |
| cma_realistic_q16 | near_axis | 416 | 163 | 39.18 | 0.0153 |
| cma_realistic_q16 | far_axis | 1584 | 454 | 28.66 | 0.0140 |
<!-- /table:phase_d_pilot_axis -->

<!-- table:phase_d_pilot_boundary -->
| run | subset | n | wrong | rate | median_right |
|---|---|---|---|---|---|
| mc_clean | boundary | 1411 | 568 | 40.26 | 0.0122 |
| mc_clean | interior | 589 | 172 | 29.20 | 0.0109 |
| cma_clean | boundary | 1411 | 570 | 40.40 | 0.0057 |
| cma_clean | interior | 589 | 172 | 29.20 | 0.0064 |
| cma_noretry_clean | boundary | 1411 | 562 | 39.83 | 0.0055 |
| cma_noretry_clean | interior | 589 | 172 | 29.20 | 0.0055 |
| mc_realistic_q16 | boundary | 1411 | 541 | 38.34 | 0.0161 |
| mc_realistic_q16 | interior | 589 | 168 | 28.52 | 0.0100 |
| cma_realistic_q16 | boundary | 1411 | 483 | 34.23 | 0.0095 |
| cma_realistic_q16 | interior | 589 | 134 | 22.75 | 0.0249 |
<!-- /table:phase_d_pilot_boundary -->

<!-- table:phase_d_pilot_cost -->
| run | seeds | seed_rej | nb_fits | nb_first | retry acc/att | revisit acc/att | ev_seed | ev_nb | ev_rev | seed (min) | neighbour (min) | revisit (min) | s / seed | s / neighbour fit |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mc_clean | 34 | 2 | 1966 | 770 | 0/0 | 488/1216 | 11,274,936 | 1,290,496 | 795,595 | 114.3 | 14.0 | 7.7 | 201.8 | 0.426 |
| cma_clean | 34 | 2 | 1966 | 783 | 2/1895 | 502/1214 | 11,572,106 | 790,399 | 483,426 | 121.1 | 9.3 | 5.5 | 213.7 | 0.283 |
| cma_noretry_clean | 33 | 1 | 1967 | 783 | 0/0 | 508/1215 | 11,177,832 | 493,717 | 304,965 | 113.8 | 6.4 | 3.1 | 206.9 | 0.195 |
| mc_realistic_q16 | 36 | 4 | 1964 | 755 | 0/0 | 522/1272 | 12,161,199 | 1,291,148 | 834,615 | 120.9 | 15.0 | 8.6 | 201.6 | 0.457 |
| cma_realistic_q16 | 35 | 1 | 1965 | 780 | 4/1978 | 603/1395 | 12,654,234 | 790,650 | 549,188 | 128.9 | 9.4 | 5.6 | 220.9 | 0.287 |
<!-- /table:phase_d_pilot_cost -->

<!-- table:phase_d_pilot_diag -->
| run | n_grains | grains_with_a_seed | lost_grains | lost_grains_with_a_seed | voxels_in_lost_grains | unresolved_in_lost_grains | unresolved_total | accepted_total | accepted_wrong | accepted_wrong_cross_boundary_like |
|---|---|---|---|---|---|---|---|---|---|---|
| mc_clean | 60 | 32 | 28 | 0 | 722 | 696 | 710 | 1290 | 32 | 32 |
| cma_clean | 60 | 32 | 28 | 0 | 722 | 673 | 681 | 1319 | 61 | 61 |
| cma_noretry_clean | 60 | 32 | 28 | 0 | 722 | 676 | 677 | 1323 | 57 | 57 |
| mc_realistic_q16 | 60 | 32 | 27 | 0 | 696 | 677 | 691 | 1309 | 19 | 19 |
| cma_realistic_q16 | 60 | 34 | 26 | 0 | 608 | 579 | 580 | 1420 | 37 | 37 |
<!-- /table:phase_d_pilot_diag -->

<!-- table:phase_d_pilot_main -->
| run | n | wrong | wrong % | Wilson 95% (%) | grain-bootstrap 95% (%) | median right err (deg) | unresolved | found | partial | lost | frag | wall (h, contended) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mc_clean | 2000 | 740 | 37.00 | 34.9-39.1 | 24.5-51.5 | 0.0119 | 710 | 32 | 0 | 28 | 29 | 2.27 |
| cma_clean | 2000 | 742 | 37.10 | 35.0-39.2 | 24.6-51.6 | 0.0061 | 681 | 32 | 0 | 28 | 32 | 2.27 |
| cma_noretry_clean | 2000 | 734 | 36.70 | 34.6-38.8 | 24.2-51.3 | 0.0055 | 677 | 32 | 0 | 28 | 31 | 2.06 |
| mc_realistic_q16 | 2000 | 709 | 35.45 | 33.4-37.6 | 22.7-50.2 | 0.0124 | 691 | 32 | 1 | 27 | 25 | 2.41 |
| cma_realistic_q16 | 2000 | 617 | 30.85 | 28.9-32.9 | 19.2-44.2 | 0.0141 | 580 | 34 | 0 | 26 | 26 | 2.40 |
<!-- /table:phase_d_pilot_main -->

<!-- table:phase_d_pilot_pairs -->
| pair | wrong_a | wrong_b | n_clusters | n_clusters_differing | p |
|---|---|---|---|---|---|
| mc_clean vs cma_clean | 740 | 742 | 60 | 6 | 0.8750 |
| cma_clean vs cma_noretry_clean | 742 | 734 | 60 | 2 | 0.5000 |
| mc_realistic_q16 vs cma_realistic_q16 | 709 | 617 | 60 | 11 | 0.3682 |
<!-- /table:phase_d_pilot_pairs -->

<!-- table:phase_d_pilot_pairs_accepted -->
| pair | wrong_a | wrong_b | n_clusters | n_clusters_differing | p |
|---|---|---|---|---|---|
| mc_clean vs cma_clean | 32 | 61 | 60 | 19 | 0.0001 |
| cma_clean vs cma_noretry_clean | 61 | 57 | 60 | 5 | 0.3125 |
| mc_realistic_q16 vs cma_realistic_q16 | 19 | 37 | 60 | 15 | 0.0007 |
<!-- /table:phase_d_pilot_pairs_accepted -->

<!-- table:phase_d_pilot_projection -->
| run | proj_full_h_voxel_linear | proj_full_h_seeds_per_grain | proj_full_h_every_piece_seeded |
|---|---|---|---|
| mc_clean | 27.9 | 20.2 | 38.7 |
| cma_clean | 27.8 | 19.8 | 39.4 |
| cma_noretry_clean | 25.3 | 17.7 | 37.1 |
| mc_realistic_q16 | 29.6 | 21.5 | 39.1 |
| cma_realistic_q16 | 29.5 | 20.9 | 40.6 |
<!-- /table:phase_d_pilot_projection -->

<!-- table:phase_d_pilot_source -->
| run | source | n | wrong | rate | median_right |
|---|---|---|---|---|---|
| mc_clean | seed | 32 | 0 | 0.00 | 0.0126 |
| mc_clean | neighbor | 770 | 28 | 3.64 | 0.0065 |
| mc_clean | revisit | 488 | 4 | 0.82 | 0.0136 |
| mc_clean | unresolved | 710 | 708 | 99.72 | 0.0657 |
| cma_clean | seed | 32 | 0 | 0.00 | 0.0100 |
| cma_clean | neighbor | 783 | 42 | 5.36 | 0.0052 |
| cma_clean | neighbor_retry | 2 | 2 | 100.00 | n/a |
| cma_clean | revisit | 502 | 17 | 3.39 | 0.0070 |
| cma_clean | unresolved | 681 | 681 | 100.00 | n/a |
| cma_noretry_clean | seed | 32 | 0 | 0.00 | 0.0076 |
| cma_noretry_clean | neighbor | 783 | 42 | 5.36 | 0.0042 |
| cma_noretry_clean | revisit | 508 | 15 | 2.95 | 0.0077 |
| cma_noretry_clean | unresolved | 677 | 677 | 100.00 | n/a |
| mc_realistic_q16 | seed | 32 | 0 | 0.00 | 0.0164 |
| mc_realistic_q16 | neighbor | 755 | 14 | 1.85 | 0.0083 |
| mc_realistic_q16 | revisit | 522 | 5 | 0.96 | 0.0196 |
| mc_realistic_q16 | unresolved | 691 | 690 | 99.86 | 0.2167 |
| cma_realistic_q16 | seed | 34 | 0 | 0.00 | 0.0109 |
| cma_realistic_q16 | neighbor | 780 | 19 | 2.44 | 0.0143 |
| cma_realistic_q16 | neighbor_retry | 3 | 3 | 100.00 | n/a |
| cma_realistic_q16 | revisit | 603 | 15 | 2.49 | 0.0140 |
| cma_realistic_q16 | unresolved | 580 | 580 | 100.00 | n/a |
<!-- /table:phase_d_pilot_source -->

<!-- table:phase_d_pilot_unresolved -->
| run | unresolved | wrong | median_err | max_err | wrong_all | wrong_boundary | cross |
|---|---|---|---|---|---|---|---|
| mc_clean | 710 | 708 | 41.475 | 58.629 | 740 | 568 | 497 |
| cma_clean | 681 | 681 | 41.662 | 58.666 | 742 | 570 | 505 |
| cma_noretry_clean | 677 | 677 | 41.703 | 58.669 | 734 | 562 | 494 |
| mc_realistic_q16 | 691 | 690 | 42.042 | 61.073 | 709 | 541 | 485 |
| cma_realistic_q16 | 580 | 580 | 43.500 | 58.651 | 617 | 483 | 452 |
<!-- /table:phase_d_pilot_unresolved -->
