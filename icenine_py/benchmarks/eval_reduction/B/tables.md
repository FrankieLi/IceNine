<!-- table:evalB_denominator_e0 -->
| cost | rule | under 0.05 deg | 0.05 to 0.5 deg | 0.5 to 5 deg |
|---|---|---|---|---|
| complete cost | lenient | 11 of 50 | 225 of 460 | 167 of 173 |
| complete cost | C++ | 12 of 50 | 279 of 460 | 168 of 173 |
| complete cost | ray rule R2 | 11 of 50 | 225 of 460 | 167 of 173 |
| level-0 screen | lenient | 0 of 0 | 63 of 123 | 30 of 36 |
| level-0 screen | C++ | 0 of 0 | 64 of 123 | 30 of 36 |
| level-0 screen | ray rule R2 | 0 of 0 | 63 of 123 | 30 of 36 |
<!-- /table:evalB_denominator_e0 -->

<!-- table:evalB_denominator_full -->
| cost | rule | under 0.05 deg | 0.05 to 0.5 deg | 0.5 to 5 deg |
|---|---|---|---|---|
| complete cost | lenient | 0 of 0 | 42 of 84 | 810 of 909 |
| complete cost | C++ | 0 of 0 | 44 of 84 | 868 of 909 |
| complete cost | ray rule R2 | 0 of 0 | 42 of 84 | 844 of 909 |
| level-0 screen | lenient | 0 of 0 | 3 of 6 | 221 of 264 |
| level-0 screen | C++ | 0 of 0 | 3 of 6 | 240 of 264 |
| level-0 screen | ray rule R2 | 0 of 0 | 3 of 6 | 229 of 264 |
<!-- /table:evalB_denominator_full -->

<!-- table:evalB_edges_per_deg_e0 -->
| event type | kind | complete cost, under 0.05 deg | complete cost, 0.05 to 0.5 deg | complete cost, 0.5 to 5 deg | level-0 screen, under 0.05 deg | level-0 screen, 0.05 to 0.5 deg | level-0 screen, 0.5 to 5 deg |
|---|---|---|---|---|---|---|---|
| bragg_limit | geometry, hard cut | 0 of 0 | 0 of 0.15 | 0 of 0.05 | 0 of 0 | 0 of 0.07 | 0 of 0 |
| eta_cut | geometry, hard cut | 0 of 0.60 | 0 of 1.00 | 0 of 1.1 | 0 of 0 | 0.23 of 0.40 | 0.05 of 0.23 |
| no_frame | geometry, hard cut | 0.67 of 1.8 | 0.12 of 1.0 | 0.03 of 0.79 | 0 of 0 | 0.03 of 0.26 | 0.02 of 0.13 |
| frame | physics, hard cut | 43 of 82 | 26 of 86 | 3.0 of 83 | 11 of 21 | 11 of 22 | 2.6 of 19 |
| detector | geometry, hard cut | 0 of 2.1 | 0.01 of 2.1 | 0 of 2.8 | 0.27 of 0.27 | 0.04 of 0.12 | 0.02 of 0.24 |
| pixel | artefact | 2,066 of 4,876 | 750 of 4,878 | 31 of 4,814 | 0 of 296 | 1.7 of 282 | 1.8 of 266 |
| subpixel | artefact | 0 of 0 | 0 of 0 | 0 of 0 | 0 of 0 | 0 of 0 | 0 of 0 |
| rule | rule design | 0 of 0 | 0 of 0 | 0 of 0 | 0 of 0 | 0 of 0 | 0 of 0 |
<!-- /table:evalB_edges_per_deg_e0 -->

<!-- table:evalB_edges_per_deg_full -->
| event type | kind | complete cost, under 0.05 deg | complete cost, 0.05 to 0.5 deg | complete cost, 0.5 to 5 deg | level-0 screen, under 0.05 deg | level-0 screen, 0.05 to 0.5 deg | level-0 screen, 0.5 to 5 deg |
|---|---|---|---|---|---|---|---|
| bragg_limit | geometry, hard cut | 0 of 0 | 0 of 0 | 0 of 0.09 | 0 of 0 | 0 of 0 | 0 of 0.01 |
| eta_cut | geometry, hard cut | 0 of 0 | 0 of 0.96 | 0 of 0.97 | 0 of 0 | 0 of 0 | 0 of 0.25 |
| no_frame | geometry, hard cut | 0 of 0 | 0.02 of 1.2 | 0.02 of 1.0 | 0 of 0 | 0.02 of 0.09 | 0.01 of 0.31 |
| frame | physics, hard cut | 51 of 81 | 39 of 86 | 2.7 of 86 | 7.8 of 13 | 11 of 20 | 1.8 of 20 |
| detector | geometry, hard cut | 0 of 2.0 | 0 of 2.8 | 0 of 3.2 | 0 of 0 | 0.15 of 0.20 | 0.01 of 0.20 |
| pixel | artefact | 150 of 5,091 | 372 of 4,982 | 26 of 4,975 | 0 of 308 | 0.91 of 272 | 1.4 of 277 |
| subpixel | artefact | 0 of 0 | 0 of 0 | 0 of 0 | 0 of 0 | 0 of 0 | 0 of 0 |
| rule | rule design | 0 of 0 | 0 of 0 | 0 of 0 | 0 of 0 | 0 of 0 | 0 of 0 |
<!-- /table:evalB_edges_per_deg_full -->

<!-- table:evalB_flips_e0 -->
| cost | rule | cause: frame | cause: pixel | cause: detector | cause: eta_cut | cause: no_frame | cause: rule | cause: unattributed | per degree |
|---|---|---|---|---|---|---|---|---|---|
| complete cost | lenient | 10535 | 22052 | 1 | 14 | 62 | 0 | 0 | 78 |
| level-0 screen | lenient | 3692 | 67 | 0 | 67 | 10 | 0 | 0 | 9.1 |
| complete cost | C++ | 8720 | 26467 | 1 | 0 | 56 | 0 | 0 | 84 |
| level-0 screen | C++ | 3658 | 663 | 22 | 67 | 10 | 0 | 0 | 11 |
| complete cost | ray rule R2 | 10533 | 22058 | 1 | 14 | 61 | 0 | 0 | 78 |
| level-0 screen | ray rule R2 | 3692 | 67 | 0 | 67 | 10 | 0 | 0 | 9.1 |
<!-- /table:evalB_flips_e0 -->

<!-- table:evalB_flips_full -->
| cost | rule | cause: frame | cause: pixel | cause: detector | cause: eta_cut | cause: no_frame | cause: rule | cause: unattributed | per degree |
|---|---|---|---|---|---|---|---|---|---|
| complete cost | lenient | 5972 | 5729 | 0 | 4 | 19 | 0 | 0 | 20 |
| level-0 screen | lenient | 2538 | 857 | 4 | 10 | 14 | 0 | 0 | 5.7 |
| complete cost | C++ | 3878 | 4736 | 0 | 0 | 11 | 0 | 0 | 14 |
| level-0 screen | C++ | 1610 | 797 | 14 | 0 | 5 | 0 | 0 | 4.0 |
| complete cost | ray rule R2 | 4524 | 4147 | 0 | 0 | 15 | 0 | 0 | 14 |
| level-0 screen | ray rule R2 | 1779 | 633 | 4 | 6 | 7 | 0 | 0 | 4.0 |
<!-- /table:evalB_flips_full -->

<!-- table:evalB_gross_share_e0 -->
| cost | regime | net variation | gross change | bragg_limit | eta_cut | no_frame | frame | detector | pixel | subpixel | rule | denominator | unattributed |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| complete cost | under 0.05 deg | 134 | 134 | 0.0% | 0.0% | 0.2% | 12.8% | 0.0% | 86.5% | 0.0% | 0.0% | 0.5% | 0.0% |
| complete cost | 0.05 to 0.5 deg | 404 | 405 | 0.0% | 0.0% | 0.1% | 17.0% | 0.0% | 82.4% | 0.0% | 0.0% | 0.5% | 0.0% |
| complete cost | 0.5 to 5 deg | 8.1 | 8.1 | 0.0% | 0.0% | 0.2% | 34.9% | 0.0% | 64.2% | 0.0% | 0.0% | 0.8% | 0.0% |
| level-0 screen | under 0.05 deg | 18 | 18 | 0.0% | 0.0% | 0.0% | 98.7% | 1.3% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% |
| level-0 screen | 0.05 to 0.5 deg | 197 | 205 | 0.0% | 1.9% | 0.2% | 81.4% | 0.2% | 12.6% | 0.0% | 0.0% | 3.7% | 0.0% |
| level-0 screen | 0.5 to 5 deg | 29 | 30 | 0.0% | 1.3% | 0.4% | 56.5% | 0.2% | 39.7% | 0.0% | 0.0% | 1.9% | 0.0% |
<!-- /table:evalB_gross_share_e0 -->

<!-- table:evalB_gross_share_full -->
| cost | regime | net variation | gross change | bragg_limit | eta_cut | no_frame | frame | detector | pixel | subpixel | rule | denominator | unattributed |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| complete cost | under 0.05 deg | 5.8 | 5.8 | 0.0% | 0.0% | 0.0% | 73.7% | 0.0% | 26.3% | 0.0% | 0.0% | 0.0% | 0.0% |
| complete cost | 0.05 to 0.5 deg | 75 | 75 | 0.0% | 0.0% | 0.0% | 39.3% | 0.0% | 59.5% | 0.0% | 0.0% | 1.2% | 0.0% |
| complete cost | 0.5 to 5 deg | 43 | 43 | 0.0% | 0.0% | 0.2% | 42.9% | 0.0% | 55.6% | 0.0% | 0.0% | 1.3% | 0.0% |
| level-0 screen | under 0.05 deg | 2.3 | 2.3 | 0.0% | 0.0% | 0.0% | 100.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% |
| level-0 screen | 0.05 to 0.5 deg | 35 | 35 | 0.0% | 0.0% | 0.2% | 90.9% | 0.6% | 7.8% | 0.0% | 0.0% | 0.6% | 0.0% |
| level-0 screen | 0.5 to 5 deg | 88 | 89 | 0.0% | 0.0% | 0.3% | 53.6% | 0.2% | 43.5% | 0.0% | 0.0% | 2.4% | 0.0% |
<!-- /table:evalB_gross_share_full -->

<!-- table:evalB_jumps_e0 -->
| type | under 0.05 deg | 0.05 to 0.5 deg | 0.5 to 5 deg |
|---|---|---|---|
| frame | 0.0237 / 0.0237 (n=1277) | 0.0133 / 0.0237 (n=7078) | 0.0133 / 0.0237 (n=365) |
| pixel | 0.0024 / 0.0075 (n=61988) | 0.0013 / 0.0075 (n=202622) | 0.0013 / 0.0075 (n=3765) |
| detector | none | 0.0075 / 0.0075 (n=4) | none |
| eta_cut | none | none | none |
| no_frame | 0.0237 / 0.0237 (n=20) | 0.0237 / 0.0237 (n=32) | 0.0042 / 0.0133 (n=4) |
| bragg_limit | none | none | none |
<!-- /table:evalB_jumps_e0 -->

<!-- table:evalB_jumps_full -->
| type | under 0.05 deg | 0.05 to 0.5 deg | 0.5 to 5 deg |
|---|---|---|---|
| frame | 0.0237 / 0.0237 (n=308) | 0.0237 / 0.0237 (n=2109) | 0.0237 / 0.0237 (n=1472) |
| pixel | 0.0024 / 0.0075 (n=900) | 0.0024 / 0.0133 (n=20081) | 0.0013 / 0.0075 (n=14018) |
| detector | none | none | none |
| eta_cut | none | none | none |
| no_frame | none | 0.0237 / 0.0237 (n=1) | 0.0133 / 0.0237 (n=10) |
| bragg_limit | none | none | none |
<!-- /table:evalB_jumps_full -->

<!-- table:evalB_minima_e0 -->
| cost | rule | walk | prominence | rays with a minimum | minima | median depth | median width (deg) | median position (deg) | cause of the dip (top 3) |
|---|---|---|---|---|---|---|---|---|---|
| complete cost | lenient | dense | 0.001 | 524 of 600 | 1727 | 0.0016 | 0.000 | 0.19 | pixel 1520, denominator 199, frame 7 |
| complete cost | lenient | dense | 0.01 | 13 of 600 | 17 | 0.0105 | 0.002 | 0.14 | denominator 12, frame 5 |
| complete cost | lenient | dense | 0.05 | 0 of 600 | 0 | n/a | n/a | n/a | none |
| complete cost | lenient | profile | 0.001 | 189 of 600 | 272 | 0.0021 | 0.004 | 0.00 | none |
| complete cost | lenient | profile | 0.01 | 1 of 600 | 1 | 0.0106 | 0.004 | 0.04 | none |
| complete cost | lenient | profile | 0.05 | 0 of 600 | 0 | n/a | n/a | n/a | none |
| complete cost | C++ | dense | 0.001 | 525 of 600 | 1780 | 0.0016 | 0.000 | 0.19 | pixel 1569, denominator 203, frame 7 |
| complete cost | C++ | dense | 0.01 | 12 of 600 | 16 | 0.0119 | 0.002 | 0.13 | denominator 11, frame 5 |
| complete cost | C++ | dense | 0.05 | 0 of 600 | 0 | n/a | n/a | n/a | none |
| complete cost | C++ | profile | 0.001 | 197 of 600 | 285 | 0.0021 | 0.004 | 0.00 | none |
| complete cost | C++ | profile | 0.01 | 2 of 600 | 2 | 0.0102 | 0.004 | 0.03 | none |
| complete cost | C++ | profile | 0.05 | 0 of 600 | 0 | n/a | n/a | n/a | none |
| complete cost | ray rule R2 | dense | 0.001 | 563 of 600 | 2200 | 0.0015 | 0.001 | 0.22 | pixel 1993, denominator 199, frame 7 |
| complete cost | ray rule R2 | dense | 0.01 | 11 of 600 | 15 | 0.0105 | 0.002 | 0.14 | denominator 10, frame 5 |
| complete cost | ray rule R2 | dense | 0.05 | 0 of 600 | 0 | n/a | n/a | n/a | none |
| complete cost | ray rule R2 | profile | 0.001 | 196 of 600 | 281 | 0.0020 | 0.004 | 0.00 | none |
| complete cost | ray rule R2 | profile | 0.01 | 1 of 600 | 1 | 0.0107 | 0.004 | 0.04 | none |
| complete cost | ray rule R2 | profile | 0.05 | 0 of 600 | 0 | n/a | n/a | n/a | none |
| level-0 screen | lenient | dense | 0.001 | 38 of 600 | 46 | 0.0386 | 0.041 | 0.23 | denominator 28, detector 13, frame 5 |
| level-0 screen | lenient | dense | 0.01 | 37 of 600 | 45 | 0.0404 | 0.042 | 0.23 | denominator 27, detector 13, frame 5 |
| level-0 screen | lenient | dense | 0.05 | 6 of 600 | 11 | 0.0714 | 0.015 | 0.23 | denominator 6, frame 5 |
| level-0 screen | lenient | profile | 0.001 | 159 of 600 | 171 | 0.0087 | 0.006 | 1.69 | none |
| level-0 screen | lenient | profile | 0.01 | 78 of 600 | 80 | 0.0207 | 0.006 | 0.60 | none |
| level-0 screen | lenient | profile | 0.05 | 6 of 600 | 6 | 0.0646 | 0.020 | 0.18 | none |
| level-0 screen | C++ | dense | 0.001 | 37 of 600 | 45 | 0.0404 | 0.042 | 0.23 | denominator 27, detector 13, frame 5 |
| level-0 screen | C++ | dense | 0.01 | 37 of 600 | 45 | 0.0404 | 0.042 | 0.23 | denominator 27, detector 13, frame 5 |
| level-0 screen | C++ | dense | 0.05 | 6 of 600 | 11 | 0.0714 | 0.015 | 0.23 | denominator 6, frame 5 |
| level-0 screen | C++ | profile | 0.001 | 158 of 600 | 170 | 0.0107 | 0.007 | 1.63 | none |
| level-0 screen | C++ | profile | 0.01 | 84 of 600 | 86 | 0.0283 | 0.007 | 0.70 | none |
| level-0 screen | C++ | profile | 0.05 | 8 of 600 | 8 | 0.0554 | 0.020 | 0.45 | none |
| level-0 screen | ray rule R2 | dense | 0.001 | 23 of 600 | 29 | 0.0460 | 0.063 | 0.23 | denominator 24, frame 5 |
| level-0 screen | ray rule R2 | dense | 0.01 | 23 of 600 | 29 | 0.0460 | 0.063 | 0.23 | denominator 24, frame 5 |
| level-0 screen | ray rule R2 | dense | 0.05 | 6 of 600 | 10 | 0.0714 | 0.008 | 0.22 | denominator 5, frame 5 |
| level-0 screen | ray rule R2 | profile | 0.001 | 151 of 600 | 159 | 0.0104 | 0.007 | 1.83 | none |
| level-0 screen | ray rule R2 | profile | 0.01 | 79 of 600 | 81 | 0.0250 | 0.007 | 0.91 | none |
| level-0 screen | ray rule R2 | profile | 0.05 | 8 of 600 | 8 | 0.0611 | 0.014 | 0.48 | none |
<!-- /table:evalB_minima_e0 -->

<!-- table:evalB_minima_full -->
| cost | rule | walk | prominence | rays with a minimum | minima | median depth | median width (deg) | median position (deg) | cause of the dip (top 3) |
|---|---|---|---|---|---|---|---|---|---|
| complete cost | lenient | dense | 0.001 | 120 of 120 | 3641 | 0.0018 | 0.019 | 1.98 | pixel 3120, frame 414, denominator 104 |
| complete cost | lenient | dense | 0.01 | 97 of 120 | 171 | 0.0138 | 0.204 | 3.24 | pixel 140, frame 24, denominator 7 |
| complete cost | lenient | dense | 0.05 | 0 of 120 | 0 | n/a | n/a | n/a | none |
| complete cost | C++ | dense | 0.001 | 118 of 120 | 1040 | 0.0022 | 0.011 | 0.97 | pixel 880, denominator 99, frame 59 |
| complete cost | C++ | dense | 0.01 | 72 of 120 | 109 | 0.0132 | 0.136 | 2.96 | pixel 90, frame 14, denominator 5 |
| complete cost | C++ | dense | 0.05 | 0 of 120 | 0 | n/a | n/a | n/a | none |
| complete cost | ray rule R2 | dense | 0.001 | 120 of 120 | 2862 | 0.0016 | 0.006 | 0.83 | pixel 2663, denominator 111, frame 82 |
| complete cost | ray rule R2 | dense | 0.01 | 74 of 120 | 131 | 0.0132 | 0.187 | 3.14 | pixel 107, frame 13, denominator 10 |
| complete cost | ray rule R2 | dense | 0.05 | 0 of 120 | 0 | n/a | n/a | n/a | none |
| level-0 screen | lenient | dense | 0.001 | 120 of 120 | 577 | 0.0156 | 0.140 | 2.54 | frame 341, pixel 169, denominator 65 |
| level-0 screen | lenient | dense | 0.01 | 120 of 120 | 536 | 0.0190 | 0.142 | 2.55 | frame 331, pixel 169, denominator 34 |
| level-0 screen | lenient | dense | 0.05 | 71 of 120 | 92 | 0.0614 | 0.386 | 2.96 | pixel 49, frame 42, denominator 1 |
| level-0 screen | C++ | dense | 0.001 | 103 of 120 | 210 | 0.0394 | 0.291 | 2.50 | pixel 80, denominator 66, frame 64 |
| level-0 screen | C++ | dense | 0.01 | 92 of 120 | 180 | 0.0526 | 0.303 | 2.54 | pixel 79, frame 61, denominator 40 |
| level-0 screen | C++ | dense | 0.05 | 68 of 120 | 98 | 0.0556 | 0.333 | 2.68 | pixel 69, frame 28, denominator 1 |
| level-0 screen | ray rule R2 | dense | 0.001 | 108 of 120 | 215 | 0.0385 | 0.299 | 2.57 | frame 76, pixel 72, denominator 66 |
| level-0 screen | ray rule R2 | dense | 0.01 | 100 of 120 | 184 | 0.0526 | 0.300 | 2.59 | frame 73, pixel 72, denominator 38 |
| level-0 screen | ray rule R2 | dense | 0.05 | 74 of 120 | 96 | 0.0556 | 0.394 | 2.72 | pixel 59, frame 37 |
<!-- /table:evalB_minima_full -->

<!-- table:evalB_rules_e0 -->
| cost | rule | under 0.05 deg: cost-changing events/deg | under 0.05 deg: steps with a cost change | under 0.05 deg: steps on the cost = 1 plateau | 0.05 to 0.5 deg: cost-changing events/deg | 0.05 to 0.5 deg: steps with a cost change | 0.05 to 0.5 deg: steps on the cost = 1 plateau | 0.5 to 5 deg: cost-changing events/deg | 0.5 to 5 deg: steps with a cost change | 0.5 to 5 deg: steps on the cost = 1 plateau |
|---|---|---|---|---|---|---|---|---|---|---|
| complete cost | lenient | 2,140 | 38931 of 60000 | 0 | 1,107 | 214581 of 540000 | 0 | 71 | 8039 of 240000 | 113878 |
| complete cost | C++ | 2,109 | 38661 of 60000 | 0 | 777 | 157989 of 540000 | 0 | 34 | 4044 of 240000 | 135130 |
| complete cost | C++ + CppEdgeRules | 2,109 | 38650 of 60000 | 0 | 776 | 157736 of 540000 | 0 | 34 | 4025 of 240000 | 135530 |
| complete cost | ray rule R2 | 2,138 | 38920 of 60000 | 0 | 1,105 | 214358 of 540000 | 0 | 71 | 8024 of 240000 | 113878 |
| level-0 screen | lenient | 11 | 335 of 60000 | 0 | 13 | 3585 of 540000 | 0 | 5.1 | 630 of 240000 | 96690 |
| level-0 screen | C++ | 11 | 335 of 60000 | 0 | 13 | 3574 of 540000 | 0 | 4.5 | 549 of 240000 | 110941 |
| level-0 screen | C++ + CppEdgeRules | 11 | 335 of 60000 | 0 | 13 | 3569 of 540000 | 0 | 4.4 | 546 of 240000 | 110941 |
| level-0 screen | ray rule R2 | 11 | 327 of 60000 | 0 | 11 | 3126 of 540000 | 0 | 3.4 | 429 of 240000 | 96690 |
<!-- /table:evalB_rules_e0 -->

<!-- table:evalB_rules_full -->
| cost | rule | under 0.05 deg: cost-changing events/deg | under 0.05 deg: steps with a cost change | under 0.05 deg: steps on the cost = 1 plateau | 0.05 to 0.5 deg: cost-changing events/deg | 0.05 to 0.5 deg: steps with a cost change | 0.05 to 0.5 deg: steps on the cost = 1 plateau | 0.5 to 5 deg: cost-changing events/deg | 0.5 to 5 deg: steps with a cost change | 0.5 to 5 deg: steps on the cost = 1 plateau |
|---|---|---|---|---|---|---|---|---|---|---|
| complete cost | lenient | 210 | 1064 of 12000 | 0 | 645 | 29171 of 108000 | 0 | 141 | 72628 of 1080000 | 46723 |
| complete cost | C++ | 201 | 1026 of 12000 | 0 | 411 | 19898 of 108000 | 0 | 29 | 15575 of 1080000 | 405232 |
| complete cost | C++ + CppEdgeRules | 201 | 1026 of 12000 | 0 | 411 | 19891 of 108000 | 0 | 29 | 15487 of 1080000 | 410844 |
| complete cost | ray rule R2 | 203 | 1032 of 12000 | 0 | 619 | 28222 of 108000 | 0 | 64 | 33332 of 1080000 | 285050 |
| level-0 screen | lenient | 7.8 | 47 of 12000 | 0 | 13 | 680 of 108000 | 0 | 6.4 | 3638 of 1080000 | 81706 |
| level-0 screen | C++ | 7.8 | 47 of 12000 | 0 | 12 | 666 of 108000 | 0 | 3.2 | 1879 of 1080000 | 343775 |
| level-0 screen | C++ + CppEdgeRules | 7.8 | 47 of 12000 | 0 | 12 | 669 of 108000 | 0 | 3.2 | 1838 of 1080000 | 343775 |
| level-0 screen | ray rule R2 | 7.8 | 47 of 12000 | 0 | 11 | 612 of 108000 | 0 | 3.3 | 1939 of 1080000 | 278639 |
<!-- /table:evalB_rules_full -->

<!-- table:evalB_truth_lowest_e0 -->
| cost | rule | rays | rays whose lowest cost is below the truth's |
|---|---|---|---|
| complete cost | lenient | 600 | 38 |
| complete cost | C++ | 600 | 38 |
| complete cost | ray rule R2 | 600 | 38 |
| level-0 screen | lenient | 600 | 13 |
| level-0 screen | C++ | 600 | 13 |
| level-0 screen | ray rule R2 | 600 | 4 |
<!-- /table:evalB_truth_lowest_e0 -->

<!-- table:evalB_truth_lowest_full -->
| cost | rule | rays | rays whose lowest cost is below the truth's |
|---|---|---|---|
| complete cost | lenient | 120 | 0 |
| complete cost | C++ | 120 | 0 |
| complete cost | ray rule R2 | 120 | 0 |
| level-0 screen | lenient | 120 | 0 |
| level-0 screen | C++ | 120 | 0 |
| level-0 screen | ray rule R2 | 120 | 0 |
<!-- /table:evalB_truth_lowest_full -->
