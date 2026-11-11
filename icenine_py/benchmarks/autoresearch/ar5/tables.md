<!-- table:p1_pick -->
| setting | evals | cut | gain | lost | ok |
|---|---|---|---|---|---|
| full pass | 4,821,806 | 0.0 | 190 | 0.0 | yes |
| P1 d=0 | 84,963 | 98.2 | 76 | 60.0 | no |
| P1 d=0.02 | 87,037 | 98.2 | 88 | 53.7 | no |
| P1 d=0.05 | 91,429 | 98.1 | 102 | 46.3 | no |
| P1 d=0.1 | 115,036 | 97.6 | 125 | 34.2 | no |
| P1 d=0.2 | 205,804 | 95.7 | 156 | 17.9 | no |
| P2 m=1 | 2,164,026 | 55.1 | 181 | 4.7 | yes |
| P2 m=2 | 3,521,825 | 27.0 | 190 | 0.0 | yes |
| P3 d=0.2 m=1 | 204,401 | 95.8 | 156 | 17.9 | no |
| P3 d=0.2 m=2 | 205,438 | 95.7 | 156 | 17.9 | no |
<!-- /table:p1_pick -->

<!-- table:p1_rules -->
| setting | sample | kept | of | evals | wrong | fixed | broken | gain | lost |
|---|---|---|---|---|---|---|---|---|---|
| full | big | 18,665 | 18,665 | 1,138,565 | 24 | 0 | 0 | 43 | 0.0 |
| full | first gradient sample | 60,381 | 60,381 | 3,683,241 | 63 | 0 | 0 | 147 | 0.0 |
| full | mid | 50,099 | 50,099 | 3,056,039 | 31 | 0 | 0 | 53 | 0.0 |
| P1 d=0 | big | 20 | 18,665 | 19,885 | 53 | 2 | 31 | 14 | 67.4 |
| P1 d=0 | first gradient sample | 77 | 60,381 | 65,078 | 148 | 7 | 92 | 62 | 57.8 |
| P1 d=0 | mid | 35 | 50,099 | 52,234 | 54 | 1 | 24 | 30 | 43.4 |
| P1 d=0.02 | big | 32 | 18,665 | 20,617 | 49 | 1 | 26 | 18 | 58.1 |
| P1 d=0.02 | first gradient sample | 99 | 60,381 | 66,420 | 140 | 6 | 83 | 70 | 52.4 |
| P1 d=0.02 | mid | 43 | 50,099 | 52,722 | 53 | 1 | 23 | 31 | 41.5 |
| P1 d=0.05 | big | 56 | 18,665 | 22,081 | 45 | 0 | 21 | 22 | 48.8 |
| P1 d=0.05 | first gradient sample | 147 | 60,381 | 69,348 | 130 | 4 | 71 | 80 | 45.6 |
| P1 d=0.05 | mid | 80 | 50,099 | 54,979 | 51 | 1 | 21 | 33 | 37.7 |
| P1 d=0.1 | big | 171 | 18,665 | 29,096 | 42 | 0 | 18 | 25 | 41.9 |
| P1 d=0.1 | first gradient sample | 419 | 60,381 | 85,940 | 110 | 2 | 49 | 100 | 32.0 |
| P1 d=0.1 | mid | 203 | 50,099 | 62,482 | 49 | 0 | 18 | 35 | 34.0 |
| P1 d=0.2 | big | 855 | 18,665 | 70,820 | 36 | 0 | 12 | 31 | 27.9 |
| P1 d=0.2 | first gradient sample | 1,223 | 60,381 | 134,984 | 85 | 1 | 23 | 125 | 15.0 |
| P1 d=0.2 | mid | 1,106 | 50,099 | 117,565 | 40 | 0 | 9 | 44 | 17.0 |
| P2 m=1 | big | 11,054 | 18,665 | 692,959 | 24 | 0 | 0 | 43 | 0.0 |
| P2 m=1 | first gradient sample | 23,126 | 60,381 | 1,471,067 | 72 | 0 | 9 | 138 | 6.1 |
| P2 m=1 | mid | 17,859 | 50,099 | 1,139,498 | 36 | 0 | 5 | 48 | 9.4 |
| P2 m=2 | big | 15,257 | 18,665 | 949,342 | 24 | 0 | 0 | 43 | 0.0 |
| P2 m=2 | first gradient sample | 41,182 | 60,381 | 2,572,483 | 63 | 0 | 0 | 147 | 0.0 |
| P2 m=2 | mid | 29,774 | 50,099 | 1,866,313 | 31 | 0 | 0 | 53 | 0.0 |
| P3 d=0.2 m=1 | big | 841 | 18,665 | 69,966 | 36 | 0 | 12 | 31 | 27.9 |
| P3 d=0.2 m=1 | first gradient sample | 1,214 | 60,381 | 134,435 | 85 | 1 | 23 | 125 | 15.0 |
| P3 d=0.2 m=1 | mid | 1,068 | 50,099 | 115,247 | 40 | 0 | 9 | 44 | 17.0 |
| P3 d=0.2 m=2 | big | 854 | 18,665 | 70,759 | 36 | 0 | 12 | 31 | 27.9 |
| P3 d=0.2 m=2 | first gradient sample | 1,218 | 60,381 | 134,679 | 85 | 1 | 23 | 125 | 15.0 |
| P3 d=0.2 m=2 | mid | 1,103 | 50,099 | 117,382 | 40 | 0 | 9 | 44 | 17.0 |
<!-- /table:p1_rules -->

<!-- table:p2_failed -->
| sample | voxel | kind | cand | rank | c20 | r5 | s8best | gridhit | griderr |
|---|---|---|---|---|---|---|---|---|---|
| first gradient sample | 19058 | first seed | 45 | 1 | yes | yes | 0.792 | 0.971 | 0.02 |
| first gradient sample | 367 | reseed | 48 | 1 | yes | yes | 0.783 | 0.980 | 0.00 |
| first gradient sample | 18267 | reseed | 49 | 1 | yes | yes | 0.774 | 1.000 | 0.00 |
| first gradient sample | 20136 | reseed | 53 | 3 | yes | yes | 0.721 | 0.984 | 0.01 |
| first gradient sample | 1672 | reseed | 38 | 1 | yes | yes | 0.786 | 0.984 | 0.00 |
| first gradient sample | 3957 | reseed | 46 | 1 | yes | yes | 0.795 | 0.993 | 0.01 |
| first gradient sample | 8146 | reseed | 49 | 1 | yes | yes | 0.758 | 0.990 | 0.01 |
| first gradient sample | 2740 | reseed | 38 | 1 | yes | yes | 0.784 | 0.996 | 0.00 |
| first gradient sample | 2744 | reseed | 54 | 1 | yes | yes | 0.784 | 0.966 | 0.01 |
| first gradient sample | 3494 | reseed | 36 | 1 | yes | yes | 0.689 | 0.401 | 1.40 |
| first gradient sample | 19063 | reseed | 44 | 1 | yes | yes | 0.791 | 0.988 | 0.03 |
| first gradient sample | 6768 | reseed | 51 | 1 | yes | yes | 0.800 | 0.984 | 0.00 |
| first gradient sample | 11536 | reseed | 39 | 1 | yes | yes | 0.784 | 1.000 | 0.00 |
| first gradient sample | 11654 | reseed | 58 | 1 | yes | yes | 0.791 | 0.997 | 0.00 |
| first gradient sample | 24194 | reseed | 32 | 1 | yes | yes | 0.656 | 0.995 | 0.01 |
| first gradient sample | 23825 | reseed | 45 | 1 | yes | yes | 0.756 | 0.998 | 0.01 |
| mid | 17538 | first seed | 40 | 1 | yes | yes | 0.779 | 0.996 | 0.01 |
| mid | 7134 | reseed | 47 | 1 | yes | yes | 0.711 | 0.967 | 0.01 |
| mid | 7167 | reseed | 42 | 1 | yes | yes | 0.723 | 0.988 | 0.00 |
| mid | 20676 | reseed | 42 | 1 | yes | yes | 0.734 | 1.000 | 0.00 |
| mid | 1339 | reseed | 47 | 1 | yes | yes | 0.747 | 0.751 | 34.96 |
| mid | 1158 | reseed | 36 | 1 | yes | yes | 0.799 | 0.969 | 0.01 |
| mid | 17608 | reseed | 53 | 1 | yes | yes | 0.787 | 1.000 | 0.00 |
| mid | 15577 | reseed | 48 | 1 | yes | yes | 0.763 | 0.853 | 0.10 |
| mid | 13439 | reseed | 39 | 1 | yes | yes | 0.771 | 0.986 | 0.03 |
| mid | 8425 | reseed | 44 | 1 | yes | yes | 0.741 | 0.997 | 0.00 |
| mid | 10395 | reseed | 61 | 2 | yes | yes | 0.796 | 0.970 | 0.01 |
| mid | 10603 | reseed | 62 | 2 | yes | yes | 0.783 | 1.000 | 0.01 |
| mid | 12282 | reseed | 47 | 1 | yes | yes | 0.785 | 0.977 | 0.00 |
| mid | 16731 | reseed | 57 | 1 | yes | yes | 0.742 | 0.992 | 0.01 |
| mid | 20251 | reseed | 55 | 2 | yes | yes | 0.769 | 0.734 | 45.70 |
<!-- /table:p2_failed -->

<!-- table:p2_rungs -->
| sample | voxel | r1hit | r1ev | r1err | r3hit | r3ev | r3err | ends | ladev | gridev |
|---|---|---|---|---|---|---|---|---|---|---|
| first gradient sample | 19058 | 0.792 | 4,350 | 0.26 | 0.797 | 210,151 | 0.21 | grid | 557,093 | 342,592 |
| first gradient sample | 367 | 0.783 | 3,561 | 0.33 | 0.620 | 209,558 | 47.71 | grid | 557,965 | 344,846 |
| first gradient sample | 18267 | 0.774 | 5,218 | 52.54 | 0.864 | 212,299 | 0.11 | R3 | 217,517 | 343,763 |
| first gradient sample | 20136 | 0.721 | 4,773 | 0.21 | 0.673 | 209,577 | 24.34 | grid | 560,222 | 345,872 |
| first gradient sample | 1672 | 0.786 | 5,317 | 25.88 | 0.967 | 208,144 | 0.03 | R3 | 213,461 | 347,730 |
| first gradient sample | 3957 | 0.795 | 3,762 | 0.43 | 1.000 | 207,507 | 0.00 | R3 | 211,269 | 346,279 |
| first gradient sample | 8146 | 0.752 | 5,495 | 18.02 | 0.757 | 209,962 | 18.02 | grid | 562,048 | 346,591 |
| first gradient sample | 2740 | 0.784 | 5,483 | 0.07 | 0.877 | 211,540 | 0.05 | R3 | 217,023 | 345,796 |
| first gradient sample | 2744 | 0.784 | 4,288 | 0.06 | 0.726 | 208,260 | 0.19 | grid | 558,590 | 346,042 |
| first gradient sample | 3494 | 0.689 | 4,092 | 0.57 | 0.997 | 207,132 | 0.00 | R3 | 211,224 | 343,279 |
| first gradient sample | 19063 | 0.791 | 5,469 | 0.18 | 0.926 | 210,505 | 0.08 | R3 | 215,974 | 347,701 |
| first gradient sample | 6768 | 0.800 | 4,415 | 0.10 | 0.802 | 211,105 | 0.10 | R3 | 215,520 | 346,340 |
| first gradient sample | 11536 | 0.784 | 4,715 | 0.15 | 0.999 | 207,049 | 0.00 | R3 | 211,764 | 339,153 |
| first gradient sample | 11654 | 0.791 | 4,342 | 0.25 | 0.980 | 207,732 | 0.01 | R3 | 212,074 | 345,632 |
| first gradient sample | 24194 | 0.656 | 3,590 | 0.58 | 0.655 | 209,128 | 0.58 | grid | 557,060 | 344,342 |
| first gradient sample | 23825 | 0.756 | 4,387 | 0.16 | 0.845 | 209,810 | 0.09 | R3 | 214,197 | 342,629 |
| mid | 17538 | 0.752 | 7,630 | 0.33 | 0.758 | 202,240 | 0.33 | grid | 542,616 | 332,746 |
| mid | 7134 | 0.711 | 8,320 | 37.59 | 0.948 | 205,521 | 0.01 | R3 | 213,841 | 347,951 |
| mid | 7167 | 0.723 | 6,980 | 0.08 | 0.757 | 208,501 | 0.10 | grid | 556,571 | 341,090 |
| mid | 20676 | 0.734 | 5,813 | 0.58 | 0.809 | 204,078 | 0.07 | R3 | 209,891 | 332,369 |
| mid | 1339 | 0.747 | 6,692 | 34.96 | 0.729 | 207,345 | 35.02 | grid | 554,245 | 340,208 |
| mid | 1158 | 0.799 | 6,061 | 0.24 | 0.707 | 203,782 | 0.33 | grid | 551,380 | 341,537 |
| mid | 17608 | 0.787 | 7,442 | 0.20 | 0.863 | 211,212 | 0.15 | R3 | 218,654 | 340,577 |
| mid | 15577 | 0.763 | 5,089 | 36.67 | 0.747 | 203,739 | 36.69 | grid | 544,411 | 335,583 |
| mid | 13439 | 0.771 | 4,868 | 0.05 | 0.792 | 203,631 | 0.04 | grid | 548,421 | 339,922 |
| mid | 8425 | 0.696 | 7,772 | 56.26 | 0.702 | 204,150 | 56.30 | grid | 553,968 | 342,046 |
| mid | 10395 | 0.796 | 5,214 | 0.11 | 0.974 | 201,128 | 0.01 | R3 | 206,342 | 339,552 |
| mid | 10603 | 0.783 | 4,234 | 0.28 | 0.999 | 201,648 | 0.01 | R3 | 205,882 | 335,381 |
| mid | 12282 | 0.785 | 4,921 | 0.11 | 0.516 | 200,675 | 0.59 | grid | 540,868 | 335,272 |
| mid | 16731 | 0.731 | 6,145 | 53.41 | 0.859 | 210,412 | 0.15 | R3 | 216,557 | 342,091 |
| mid | 20251 | 0.769 | 5,853 | 0.09 | 0.989 | 201,643 | 0.02 | R3 | 207,496 | 341,439 |
<!-- /table:p2_rungs -->

<!-- table:p2_totals -->
| group | n | r1pass | r1right | r3pass | r3right | ends | ladev | gridev | share | rl | rg |
|---|---|---|---|---|---|---|---|---|---|---|---|
| first gradient sample | 16 | 0 | 13 | 10 | 13 | 0 / 10 / 6 | 5,493,001 | 5,518,587 | 99.5 | 16 | 15 |
| mid | 15 | 0 | 10 | 7 | 12 | 0 / 7 / 8 | 5,871,143 | 5,087,764 | 115.4 | 14 | 13 |
| both | 31 | 0 | 23 | 17 | 25 | 0 / 17 / 14 | 11,364,144 | 10,606,351 | 107.1 | 30 | 28 |
<!-- /table:p2_totals -->

<!-- table:p3_flag -->
| sample | thr | flagged | share | wf | wu | nw |
|---|---|---|---|---|---|---|
| big | 0.50 | 0 | 0.0 | 0 | 24 | 24 |
| big | 0.55 | 0 | 0.0 | 0 | 24 | 24 |
| big | 0.60 | 0 | 0.0 | 0 | 24 | 24 |
| big | 0.65 | 0 | 0.0 | 0 | 24 | 24 |
| big | 0.70 | 0 | 0.0 | 0 | 24 | 24 |
| big | 0.75 | 0 | 0.0 | 0 | 24 | 24 |
| big | 0.80 | 0 | 0.0 | 0 | 24 | 24 |
| big | 0.85 | 1 | 0.0 | 0 | 24 | 24 |
| big | 0.90 | 80 | 0.3 | 0 | 24 | 24 |
| big | 0.95 | 1,532 | 6.2 | 24 | 0 | 24 |
| first gradient sample | 0.50 | 0 | 0.0 | 0 | 63 | 63 |
| first gradient sample | 0.55 | 0 | 0.0 | 0 | 63 | 63 |
| first gradient sample | 0.60 | 0 | 0.0 | 0 | 63 | 63 |
| first gradient sample | 0.65 | 0 | 0.0 | 0 | 63 | 63 |
| first gradient sample | 0.70 | 1 | 0.0 | 0 | 63 | 63 |
| first gradient sample | 0.75 | 1 | 0.0 | 0 | 63 | 63 |
| first gradient sample | 0.80 | 3 | 0.0 | 0 | 63 | 63 |
| first gradient sample | 0.85 | 34 | 0.1 | 7 | 56 | 63 |
| first gradient sample | 0.90 | 294 | 1.2 | 32 | 31 | 63 |
| first gradient sample | 0.95 | 2,785 | 11.3 | 63 | 0 | 63 |
| mid | 0.50 | 0 | 0.0 | 0 | 31 | 31 |
| mid | 0.55 | 0 | 0.0 | 0 | 31 | 31 |
| mid | 0.60 | 0 | 0.0 | 0 | 31 | 31 |
| mid | 0.65 | 0 | 0.0 | 0 | 31 | 31 |
| mid | 0.70 | 0 | 0.0 | 0 | 31 | 31 |
| mid | 0.75 | 2 | 0.0 | 0 | 31 | 31 |
| mid | 0.80 | 7 | 0.0 | 2 | 29 | 31 |
| mid | 0.85 | 54 | 0.2 | 5 | 26 | 31 |
| mid | 0.90 | 484 | 2.0 | 13 | 18 | 31 |
| mid | 0.95 | 3,709 | 15.1 | 31 | 0 | 31 |
<!-- /table:p3_flag -->

<!-- table:p3_left -->
| sample | n | flagged | share | wrong | caught | left | unflagged | lo | hi |
|---|---|---|---|---|---|---|---|---|---|
| big | 24,570 | 1,532 | 6.2 | 24 | 24 | 0 | 23,038 | 0.000 | 0.017 |
| first gradient sample | 24,570 | 2,785 | 11.3 | 63 | 63 | 0 | 21,785 | 0.000 | 0.018 |
| mid | 24,570 | 3,709 | 15.1 | 31 | 31 | 0 | 20,861 | 0.000 | 0.018 |
<!-- /table:p3_left -->

<!-- table:p3_source -->
| sample | group | n | flagged | wrong | wf |
|---|---|---|---|---|---|
| big | seed | 37 | 0 | 0 | 0 |
| big | neighbor | 19,135 | 768 | 13 | 13 |
| big | revisit | 5,106 | 757 | 11 | 11 |
| big | reseed | 292 | 7 | 0 | 0 |
| big | boundary | 11,054 | 1,348 | 24 | 24 |
| big | all | 24,570 | 1,532 | 24 | 24 |
| first gradient sample | seed | 347 | 14 | 0 | 0 |
| first gradient sample | neighbor | 9,511 | 687 | 20 | 20 |
| first gradient sample | revisit | 13,780 | 1,928 | 23 | 23 |
| first gradient sample | reseed | 931 | 155 | 20 | 20 |
| first gradient sample | boundary | 23,126 | 2,766 | 63 | 63 |
| first gradient sample | all | 24,570 | 2,785 | 63 | 63 |
| mid | seed | 156 | 4 | 0 | 0 |
| mid | neighbor | 14,205 | 1,100 | 12 | 12 |
| mid | revisit | 8,754 | 2,378 | 10 | 10 |
| mid | reseed | 1,453 | 225 | 8 | 8 |
| mid | boundary | 17,859 | 3,621 | 31 | 31 |
| mid | all | 24,570 | 3,709 | 31 | 31 |
<!-- /table:p3_source -->

<!-- table:p3_v1 -->
| thr | flagged | done | ev | mean | share | kept | wrong | fixed | broken | after |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.95 | 3,709 | 3,709 | 27,882,533 | 7,518 | 105.4 | 3,485 | 31 | 13 | 0 | 18 |
<!-- /table:p3_v1 -->

<!-- table:p3_votes -->
| group | nw | nr | mw | mr | cut | w | r | sep |
|---|---|---|---|---|---|---|---|---|
| big | 24 | 1508 | 2516 | 1476 | 2446 | 13 of 24 | 173 of 1508 | yes |
| first gradient sample | 63 | 2722 | 2707 | 1714 | 3035 | 31 of 63 | 324 of 2722 | yes |
| mid | 31 | 3678 | 2553 | 1412 | 1975 | 25 of 31 | 955 of 3678 | yes |
| all three | 118 | 7908 | 2573 | 1535 | 2369 | 67 of 118 | 1410 of 7908 | yes |
<!-- /table:p3_votes -->

<!-- table:predictions -->
| n | what | result | numbers |
|---|---|---|---|
| 1 | a try-test setting cuts 70% of the pass and loses at most 5% of its gain | does not agree | best cut within the 5% loss: 55.1% |
| 2 | a candidate near the truth is in the first 100 and not the 20 costed, in half the cases | does not agree | 0 of 31 |
| 3 | the ladder passes its rule | does not agree | 107.1% of the grid evaluations |
| 4 | a hit-ratio threshold flags 80% of wrong voxels and at most 5% of voxels | does not agree | first part holds (100% of wrong voxels caught); second part fails (8.8% of voxels flagged on the pick samples, 15.1% on mid) |
| 5 | the vote count does not separate wrong from right among flagged voxels | does not agree | a weak separation exists; cut 2369: 67 of 118 wrong, 1410 of 7908 right at or above |
| 6 | T1 needs under 14 M (gradient) and 22 M (mid) evaluations | not evaluated | no rule passed, so T1 was not run |
<!-- /table:predictions -->

<!-- table:rules -->
| rule | what | verdict | numbers |
|---|---|---|---|
| try test (Part 1) | on mid: at least half of the pass cut, at most 5% of the gain lost | not met | 62.7% cut, 9.4% of the gain lost (48 of 53) |
| ladder (Part 2) | at most 2/3 of the grid evaluations, at most one case fewer right | not met | 107.1% of the evaluations; right 30 against 28 |
| V1 (Part 3) | fixes more than it breaks, fixes a quarter, costs under 5% of the run | not met | fixed 13, broken 0, 31 wrong processed; about 105% of the run |
<!-- /table:rules -->
