<!-- table:phase_d_reseed2_compare -->
| run | pilot | wrong_pilot | wrong_reseed | unres_pilot | unres_reseed | lost_pilot | lost_v1 | lost_reseed | wrong_v1 | unres_v1 | acc_wrong_pilot | acc_wrong_reseed | n_clusters_differing | p |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mc_clean | mc_clean | 740 | 37 | 710 | 1 | 28 | 19 | 0 | 633 | 601 | 32 | 36 | 30 | <1e-4 (MC floor) |
| cma_noretry_clean | cma_noretry_clean | 734 | 60 | 677 | 1 | 28 | 19 | 0 | 633 | 576 | 57 | 59 | 29 | <1e-4 (MC floor) |
| mc_realistic_q16 | mc_realistic_q16 | 709 | 23 | 691 | 2 | 27 | 18 | 0 | 608 | 590 | 19 | 21 | 30 | <1e-4 (MC floor) |
| cma_noretry_realistic_q16 | cma_realistic_q16 | 617 | 38 | 580 | 2 | 26 | 19 | 0 | 568 | 533 | 37 | 36 | 27 | <1e-4 (MC floor) |
<!-- /table:phase_d_reseed2_compare -->

<!-- table:phase_d_reseed2_counts -->
| run | seeds | seed_rej | reseeds | accepted | rejected | median_cluster | final_voxels | lost_failed | lost_nosearch | reseed_min | full_search_min | share_pct | ev_share_pct | unres_clusters | largest_unres | largest_grains | largest_lost_inside | largest_has_rej_seed | wall_h | wall_pilot_h | ev_reseed |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mc_clean | 34 | 2 | 35 | 32 | 3 | 22.0 | 3 | 0 | 0 | 109.2 | 103.6 | 47.7 | 47.7 | 1 | 1 | 1 | 0 | True | 3.82 | 2.27 | 12,201,531 |
| cma_noretry_clean | 33 | 1 | 30 | 29 | 1 | 22.5 | 1 | 0 | 0 | 93.7 | 91.6 | 46.6 | 46.4 | 1 | 1 | 1 | 0 | True | 3.35 | 2.06 | 10,364,921 |
| mc_realistic_q16 | 36 | 4 | 44 | 34 | 10 | 15.5 | 10 | 0 | 0 | 141.1 | 135.6 | 51.2 | 51.9 | 2 | 1 | 1 | 0 | False | 4.59 | 2.41 | 15,396,455 |
| cma_noretry_realistic_q16 | 43 | 10 | 42 | 28 | 14 | 15.0 | 14 | 0 | 0 | 135.3 | 133.4 | 47.8 | 48.4 | 2 | 1 | 1 | 0 | True | 4.71 | 2.40 | 15,357,913 |
<!-- /table:phase_d_reseed2_counts -->

<!-- table:phase_d_reseed2_log -->
| run | k | cluster | voxel | accepted | left_refit | seconds | evals |
|---|---|---|---|---|---|---|---|
| mc_clean | 1 | 21 | 1 | yes | 4 | 186 | 345,600 |
| mc_clean | 2 | 600 | 715 | yes | 559 | 209 | 385,157 |
| mc_clean | 3 | 559 | 1080 | yes | 496 | 207 | 382,309 |
| mc_clean | 4 | 496 | 706 | yes | 469 | 195 | 360,266 |
| mc_clean | 5 | 469 | 38 | yes | 413 | 206 | 380,428 |
| mc_clean | 6 | 282 | 731 | yes | 247 | 199 | 368,940 |
| mc_clean | 7 | 121 | 34 | yes | 66 | 206 | 380,603 |
| mc_clean | 8 | 14 | 4 | yes | 0 | 187 | 345,285 |
| mc_clean | 9 | 131 | 634 | yes | 76 | 210 | 381,807 |
| mc_clean | 10 | 21 | 10 | yes | 0 | 191 | 348,383 |
| mc_clean | 11 | 52 | 669 | no | 52 | 182 | 332,773 |
| mc_clean | 12 | 51 | 1679 | yes | 26 | 194 | 354,323 |
| mc_clean | 13 | 26 | 665 | yes | 23 | 188 | 345,043 |
| mc_clean | 14 | 23 | 668 | no | 23 | 182 | 332,682 |
| mc_clean | 15 | 22 | 1675 | yes | 6 | 190 | 348,043 |
| mc_clean | 16 | 6 | 654 | yes | 0 | 182 | 333,539 |
| mc_clean | 17 | 55 | 626 | yes | 0 | 199 | 364,974 |
| mc_clean | 18 | 125 | 1345 | yes | 102 | 189 | 349,263 |
| mc_clean | 19 | 45 | 189 | yes | 29 | 191 | 351,250 |
| mc_clean | 20 | 22 | 504 | yes | 0 | 186 | 342,733 |
| mc_clean | 21 | 7 | 62 | no | 7 | 181 | 334,268 |
| mc_clean | 22 | 16 | 67 | yes | 6 | 183 | 337,161 |
| mc_clean | 23 | 6 | 227 | yes | 0 | 184 | 339,709 |
| mc_clean | 24 | 57 | 1056 | yes | 20 | 195 | 358,715 |
| mc_clean | 25 | 59 | 854 | yes | 17 | 196 | 361,194 |
| mc_clean | 26 | 4 | 210 | yes | 0 | 176 | 332,463 |
| mc_clean | 27 | 20 | 215 | yes | 0 | 181 | 343,014 |
| mc_clean | 28 | 8 | 226 | yes | 1 | 176 | 334,363 |
| mc_clean | 29 | 4 | 238 | yes | 0 | 171 | 328,491 |
| mc_clean | 30 | 17 | 528 | yes | 3 | 176 | 338,912 |
| mc_clean | 31 | 1 | 536 | yes | 0 | 170 | 328,817 |
| mc_clean | 32 | 6 | 808 | yes | 0 | 171 | 332,540 |
| mc_clean | 33 | 3 | 961 | yes | 0 | 171 | 330,969 |
| mc_clean | 34 | 1 | 1888 | yes | 0 | 173 | 334,520 |
| mc_clean | 35 | 1 | 1938 | yes | 0 | 172 | 332,994 |
| cma_noretry_clean | 1 | 21 | 1 | yes | 4 | 177 | 330,144 |
| cma_noretry_clean | 2 | 575 | 714 | yes | 534 | 206 | 379,652 |
| cma_noretry_clean | 3 | 352 | 684 | yes | 325 | 199 | 367,890 |
| cma_noretry_clean | 4 | 325 | 1343 | yes | 262 | 188 | 349,110 |
| cma_noretry_clean | 5 | 262 | 680 | yes | 237 | 184 | 342,181 |
| cma_noretry_clean | 6 | 237 | 731 | yes | 205 | 184 | 341,370 |
| cma_noretry_clean | 7 | 93 | 686 | yes | 39 | 185 | 344,270 |
| cma_noretry_clean | 8 | 14 | 4 | yes | 0 | 196 | 362,465 |
| cma_noretry_clean | 9 | 182 | 300 | yes | 128 | 189 | 350,216 |
| cma_noretry_clean | 10 | 74 | 144 | yes | 19 | 185 | 343,062 |
| cma_noretry_clean | 11 | 19 | 10 | yes | 0 | 197 | 363,897 |
| cma_noretry_clean | 12 | 25 | 1675 | no | 25 | 192 | 356,066 |
| cma_noretry_clean | 13 | 24 | 668 | yes | 7 | 198 | 361,896 |
| cma_noretry_clean | 14 | 7 | 654 | yes | 0 | 180 | 329,909 |
| cma_noretry_clean | 15 | 54 | 575 | yes | 0 | 188 | 342,723 |
| cma_noretry_clean | 16 | 111 | 1345 | yes | 92 | 183 | 334,315 |
| cma_noretry_clean | 17 | 36 | 357 | yes | 15 | 185 | 338,103 |
| cma_noretry_clean | 18 | 16 | 67 | yes | 6 | 182 | 332,771 |
| cma_noretry_clean | 19 | 6 | 227 | yes | 0 | 198 | 360,109 |
| cma_noretry_clean | 20 | 56 | 1056 | yes | 20 | 185 | 338,810 |
| cma_noretry_clean | 21 | 15 | 799 | yes | 0 | 198 | 361,753 |
| cma_noretry_clean | 22 | 55 | 203 | yes | 15 | 184 | 339,161 |
| cma_noretry_clean | 23 | 20 | 215 | yes | 0 | 181 | 334,522 |
| cma_noretry_clean | 24 | 8 | 226 | yes | 1 | 178 | 330,497 |
| cma_noretry_clean | 25 | 4 | 238 | yes | 0 | 176 | 326,963 |
| cma_noretry_clean | 26 | 14 | 528 | yes | 0 | 180 | 332,384 |
| cma_noretry_clean | 27 | 1 | 536 | yes | 0 | 177 | 327,499 |
| cma_noretry_clean | 28 | 1 | 861 | yes | 0 | 178 | 329,025 |
| cma_noretry_clean | 29 | 1 | 1888 | yes | 0 | 194 | 356,919 |
| cma_noretry_clean | 30 | 1 | 1938 | yes | 0 | 194 | 357,239 |
| mc_realistic_q16 | 1 | 21 | 1 | yes | 4 | 198 | 349,205 |
| mc_realistic_q16 | 2 | 586 | 715 | yes | 543 | 227 | 397,595 |
| mc_realistic_q16 | 3 | 543 | 1080 | yes | 478 | 222 | 388,620 |
| mc_realistic_q16 | 4 | 449 | 38 | yes | 390 | 222 | 383,817 |
| mc_realistic_q16 | 5 | 257 | 731 | yes | 222 | 218 | 376,537 |
| mc_realistic_q16 | 6 | 97 | 686 | yes | 41 | 222 | 382,777 |
| mc_realistic_q16 | 7 | 15 | 5 | yes | 0 | 199 | 345,349 |
| mc_realistic_q16 | 8 | 133 | 634 | yes | 78 | 225 | 388,383 |
| mc_realistic_q16 | 9 | 21 | 10 | yes | 0 | 204 | 353,844 |
| mc_realistic_q16 | 10 | 26 | 1675 | no | 26 | 193 | 335,489 |
| mc_realistic_q16 | 11 | 25 | 668 | no | 25 | 191 | 335,464 |
| mc_realistic_q16 | 12 | 24 | 1676 | no | 24 | 182 | 334,721 |
| mc_realistic_q16 | 13 | 23 | 667 | yes | 14 | 190 | 348,659 |
| mc_realistic_q16 | 14 | 57 | 625 | yes | 3 | 205 | 374,969 |
| mc_realistic_q16 | 15 | 125 | 177 | yes | 89 | 200 | 366,707 |
| mc_realistic_q16 | 16 | 59 | 1184 | no | 59 | 184 | 337,659 |
| mc_realistic_q16 | 17 | 58 | 1803 | yes | 43 | 195 | 357,102 |
| mc_realistic_q16 | 18 | 21 | 43 | yes | 0 | 191 | 349,957 |
| mc_realistic_q16 | 19 | 21 | 504 | yes | 0 | 192 | 352,558 |
| mc_realistic_q16 | 20 | 16 | 67 | yes | 6 | 187 | 342,944 |
| mc_realistic_q16 | 21 | 6 | 227 | yes | 0 | 187 | 342,730 |
| mc_realistic_q16 | 22 | 14 | 670 | no | 14 | 184 | 337,483 |
| mc_realistic_q16 | 23 | 13 | 1678 | no | 13 | 183 | 335,392 |
| mc_realistic_q16 | 24 | 12 | 671 | no | 12 | 181 | 334,103 |
| mc_realistic_q16 | 25 | 11 | 152 | no | 11 | 183 | 336,365 |
| mc_realistic_q16 | 26 | 10 | 1681 | no | 10 | 182 | 334,806 |
| mc_realistic_q16 | 27 | 9 | 1147 | yes | 0 | 189 | 347,655 |
| mc_realistic_q16 | 28 | 29 | 700 | yes | 2 | 194 | 356,009 |
| mc_realistic_q16 | 29 | 1 | 156 | yes | 0 | 182 | 335,296 |
| mc_realistic_q16 | 30 | 56 | 203 | yes | 15 | 197 | 365,925 |
| mc_realistic_q16 | 31 | 30 | 541 | yes | 10 | 191 | 354,801 |
| mc_realistic_q16 | 32 | 8 | 226 | yes | 1 | 185 | 344,221 |
| mc_realistic_q16 | 33 | 4 | 238 | yes | 0 | 183 | 341,156 |
| mc_realistic_q16 | 34 | 14 | 528 | yes | 0 | 187 | 347,571 |
| mc_realistic_q16 | 35 | 10 | 537 | yes | 1 | 185 | 344,073 |
| mc_realistic_q16 | 36 | 1 | 536 | yes | 0 | 182 | 338,677 |
| mc_realistic_q16 | 37 | 1 | 789 | yes | 0 | 182 | 339,774 |
| mc_realistic_q16 | 38 | 1 | 861 | no | 1 | 181 | 337,040 |
| mc_realistic_q16 | 39 | 1 | 1296 | yes | 0 | 181 | 336,855 |
| mc_realistic_q16 | 40 | 1 | 1443 | yes | 0 | 179 | 336,139 |
| mc_realistic_q16 | 41 | 1 | 1584 | yes | 0 | 179 | 334,805 |
| mc_realistic_q16 | 42 | 1 | 1703 | yes | 0 | 180 | 336,969 |
| mc_realistic_q16 | 43 | 1 | 1888 | yes | 0 | 181 | 336,991 |
| mc_realistic_q16 | 44 | 1 | 1938 | yes | 0 | 181 | 339,263 |
| cma_noretry_realistic_q16 | 1 | 21 | 1 | yes | 4 | 201 | 365,199 |
| cma_noretry_realistic_q16 | 2 | 526 | 714 | yes | 482 | 210 | 380,244 |
| cma_noretry_realistic_q16 | 3 | 294 | 683 | yes | 267 | 209 | 377,521 |
| cma_noretry_realistic_q16 | 4 | 240 | 731 | yes | 207 | 206 | 372,773 |
| cma_noretry_realistic_q16 | 5 | 97 | 686 | yes | 41 | 209 | 378,138 |
| cma_noretry_realistic_q16 | 6 | 15 | 5 | yes | 0 | 200 | 363,331 |
| cma_noretry_realistic_q16 | 7 | 188 | 300 | yes | 132 | 211 | 381,470 |
| cma_noretry_realistic_q16 | 8 | 75 | 144 | yes | 20 | 205 | 373,638 |
| cma_noretry_realistic_q16 | 9 | 20 | 10 | yes | 0 | 200 | 365,979 |
| cma_noretry_realistic_q16 | 10 | 26 | 1675 | no | 26 | 196 | 360,039 |
| cma_noretry_realistic_q16 | 11 | 25 | 668 | no | 25 | 196 | 359,116 |
| cma_noretry_realistic_q16 | 12 | 24 | 1676 | no | 24 | 195 | 358,982 |
| cma_noretry_realistic_q16 | 13 | 23 | 667 | yes | 15 | 198 | 364,377 |
| cma_noretry_realistic_q16 | 14 | 57 | 575 | yes | 0 | 203 | 371,928 |
| cma_noretry_realistic_q16 | 15 | 109 | 1345 | yes | 89 | 200 | 366,872 |
| cma_noretry_realistic_q16 | 16 | 34 | 1356 | yes | 14 | 202 | 372,501 |
| cma_noretry_realistic_q16 | 17 | 16 | 67 | yes | 6 | 192 | 362,290 |
| cma_noretry_realistic_q16 | 18 | 6 | 227 | yes | 0 | 192 | 362,175 |
| cma_noretry_realistic_q16 | 19 | 15 | 670 | no | 15 | 190 | 361,814 |
| cma_noretry_realistic_q16 | 20 | 14 | 1678 | no | 14 | 186 | 358,909 |
| cma_noretry_realistic_q16 | 21 | 13 | 671 | no | 13 | 186 | 358,291 |
| cma_noretry_realistic_q16 | 22 | 12 | 1681 | no | 12 | 186 | 359,013 |
| cma_noretry_realistic_q16 | 23 | 11 | 669 | no | 11 | 186 | 359,249 |
| cma_noretry_realistic_q16 | 24 | 10 | 152 | no | 10 | 188 | 362,189 |
| cma_noretry_realistic_q16 | 25 | 27 | 700 | yes | 0 | 191 | 367,667 |
| cma_noretry_realistic_q16 | 26 | 55 | 1056 | yes | 20 | 192 | 373,575 |
| cma_noretry_realistic_q16 | 27 | 14 | 799 | no | 14 | 184 | 360,656 |
| cma_noretry_realistic_q16 | 28 | 13 | 1809 | no | 13 | 185 | 361,736 |
| cma_noretry_realistic_q16 | 29 | 12 | 794 | no | 12 | 186 | 361,846 |
| cma_noretry_realistic_q16 | 30 | 11 | 1806 | no | 11 | 185 | 361,441 |
| cma_noretry_realistic_q16 | 31 | 10 | 188 | yes | 0 | 187 | 364,554 |
| cma_noretry_realistic_q16 | 32 | 55 | 203 | yes | 15 | 191 | 372,708 |
| cma_noretry_realistic_q16 | 33 | 20 | 215 | yes | 0 | 188 | 367,350 |
| cma_noretry_realistic_q16 | 34 | 9 | 226 | yes | 1 | 188 | 364,376 |
| cma_noretry_realistic_q16 | 35 | 4 | 238 | yes | 0 | 187 | 364,053 |
| cma_noretry_realistic_q16 | 36 | 14 | 528 | yes | 0 | 189 | 366,427 |
| cma_noretry_realistic_q16 | 37 | 1 | 536 | yes | 0 | 186 | 361,834 |
| cma_noretry_realistic_q16 | 38 | 9 | 1148 | no | 9 | 186 | 361,100 |
| cma_noretry_realistic_q16 | 39 | 3 | 674 | yes | 0 | 188 | 367,195 |
| cma_noretry_realistic_q16 | 40 | 1 | 861 | yes | 0 | 186 | 361,029 |
| cma_noretry_realistic_q16 | 41 | 1 | 1888 | yes | 0 | 185 | 361,040 |
| cma_noretry_realistic_q16 | 42 | 1 | 1938 | yes | 0 | 187 | 363,288 |
<!-- /table:phase_d_reseed2_log -->

<!-- table:phase_d_reseed2_three -->
| run | wrong_pilot | wrong_v1 | wrong_v2 | unres_pilot | unres_v1 | unres_v2 | lost_pilot | lost_v1 | lost_v2 |
|---|---|---|---|---|---|---|---|---|---|
| mc_clean | 740 | 633 | 37 | 710 | 601 | 1 | 28 | 19 | 0 |
| cma_noretry_clean | 734 | 633 | 60 | 677 | 576 | 1 | 28 | 19 | 0 |
| mc_realistic_q16 | 709 | 608 | 23 | 691 | 590 | 2 | 27 | 18 | 0 |
| cma_noretry_realistic_q16 | 617 | 568 | 38 | 580 | 533 | 2 | 26 | 19 | 0 |
<!-- /table:phase_d_reseed2_three -->
