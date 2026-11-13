<!-- table:patience_changed -->
| sample | changed | was at 0.95 | now at 0.95 | median hit G | median hit S | median cost G | median cost S | median error G (deg) | median error S (deg) | max error S (deg) |
|---|---|---|---|---|---|---|---|---|---|---|
| held | 41 | 25 | 20 | 0.959 | 0.938 | 0.0910 | 0.1086 | 0.0284 | 0.0487 | 0.2142 |
| mid | 20 | 18 | 16 | 0.963 | 0.973 | 0.0947 | 0.0882 | 0.0375 | 0.0249 | 0.1864 |
| grad | 1 | 0 | 0 | 0.931 | 0.917 | 0.1115 | 0.1161 | 0.0431 | 0.1091 | 0.1091 |
| big | 4 | 4 | 4 | 0.969 | 0.981 | 0.0905 | 0.0821 | 0.0183 | 0.0293 | 0.0471 |
<!-- /table:patience_changed -->

<!-- table:patience_cuts -->
| sample | MC refinements | cut | evaluations in cut ones | cut end under 0.95 | cut end at least 0.95 | taken seed ends changed | median angle to G (deg) | S wrong where G right, matched steps (n) | S-only steps wrong where G's map right (n) |
|---|---|---|---|---|---|---|---|---|---|
| held | 5,488 | 763 | 2,872,706 | 420 | 343 | 41 | 0.0368 | 0 (3,640) | 0 (59) |
| mid | 1,645 | 348 | 1,296,160 | 151 | 197 | 20 | 0.0477 | 0 (1,478) | 1 (134) |
| grad | 994 | 15 | 48,986 | 3 | 12 | 1 | 0.1223 | 0 (1,280) | 0 (4) |
| big | 162 | 55 | 206,422 | 23 | 32 | 4 | 0.0193 | 0 (281) | 0 (35) |
<!-- /table:patience_cuts -->

<!-- table:patience_gate -->
| gate item | mid | grad | big |
|---|---|---|---|
| seed steps replayed (all with an MC refinement) | 770 | 478 | 99 |
| (g1) equal to G where uncut or cut after the last gain | yes | yes | yes |
| (g2) taken ends changed (limit) | 25 (30) | 1 (3) | 4 (6) |
| (g3) changed end wrong where G right | 0 | 0 | 0 |
| (g3) failed seed tests G / replay | 12 / 12 | 11 / 11 | 0 / 0 |
| (g4) evaluations saved, share of G's run (mid at least 8%) | 1,431,626 (11.01%) | 20,390 (0.20%) | 230,137 (6.04%) |
| new rescues / new grid searches | 0 / 0 | 0 / 0 | 0 / 0 |
<!-- /table:patience_gate -->

<!-- table:patience_layer -->
| sample | arm | evaluations with the pass | seed step | wrong after pass | Wilson 95% (%) | fixed / broken vs G | fixed / broken vs S8 | broken minus fixed vs S8 |
|---|---|---|---|---|---|---|---|---|
| held | G | 30,768,867 | 14,150,515 | 69 | 0.22 to 0.36 | n/a | 25 / 24 | -1 |
| held | S | 29,010,916 | 12,039,131 | 69 | 0.22 to 0.36 | 0 / 0 | 25 / 24 | -1 |
| mid | G | 13,004,457 | 5,417,024 | 35 | 0.10 to 0.20 | n/a | 9 / 13 | 4 |
| mid | S | 11,626,275 | 4,044,515 | 37 | 0.11 to 0.21 | 0 / 2 | 8 / 14 | 6 |
| grad | G | 10,196,692 | 1,528,284 | 56 | 0.18 to 0.30 | n/a | 31 / 24 | -7 |
| grad | S | 10,181,935 | 1,512,978 | 55 | 0.17 to 0.29 | 1 / 0 | 32 / 24 | -8 |
| big | G | 3,807,747 | 687,902 | 21 | 0.06 to 0.13 | n/a | 8 / 5 | -3 |
| big | S | 3,616,809 | 493,914 | 21 | 0.06 to 0.13 | 0 / 0 | 8 / 5 | -3 |
<!-- /table:patience_layer -->

<!-- table:patience_paired -->
| sample | pass | fixed vs G | broken vs G | both wrong vs G | fixed vs S8 | broken vs S8 | both wrong vs S8 |
|---|---|---|---|---|---|---|---|
| held | before | 0 | 0 | 101 | 26 | 32 | 69 |
| held | after | 0 | 0 | 69 | 25 | 24 | 45 |
| mid | before | 0 | 3 | 75 | 21 | 15 | 63 |
| mid | after | 0 | 2 | 35 | 8 | 14 | 23 |
| grad | before | 0 | 0 | 206 | 63 | 59 | 147 |
| grad | after | 1 | 0 | 55 | 32 | 24 | 31 |
| big | before | 0 | 0 | 74 | 10 | 17 | 57 |
| big | after | 0 | 0 | 21 | 8 | 5 | 16 |
<!-- /table:patience_paired -->

<!-- table:patience_predictions -->
| prediction | held | mid | grad | big |
|---|---|---|---|---|
| 1 replay (gate) | n/a | True | False | True |
| 2 replay saving | n/a | True | True | True |
| 3 held seed 10-30% / total 4-12% fewer | True (14.9%; 5.7%) | n/a | n/a | n/a |
| 4 total fewer (mid 8-12, big 4-7, grad -1..1) | n/a | True (10.6%) | True (0.14%) | True (5.0%) |
| 5 broken - fixed vs G in -4..4, median shift <= 0.0005 | True (0) | True (2) | True (-1) | True (0) |
<!-- /table:patience_predictions -->

<!-- table:patience_rules -->
| sample | (1) total fewer vs G (%) | (1) | (2) broken - fixed vs G | (2) | (3) median shift (deg) | (3) | (4) extra wrong seeds; grid S/G | (4) | (5) vs S8: fewer %; broken - fixed; median shift | (5) | G passes (5) | S passes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| held | 5.71; seed 14.92 | pass | 0 | pass | 0.00002 | pass | 0; 20/19 | pass | 51.6; -1; -0.00258 | pass | yes | pass |
| mid | 10.60 | pass | 2 | pass | -0.00003 | pass | 0; 4/4 | pass | 56.0; 6; -0.00114 | fail | yes | FAIL |
| grad | 0.14 | pass | -1 | pass | -0.00000 | pass | 0; 5/5 | pass | 43.4; -8; -0.00130 | pass | yes | pass |
| big | 5.01 | pass | 0 | pass | -0.00000 | pass | 0; 0/0 | pass | 56.6; -3; -0.00116 | pass | yes | pass |
<!-- /table:patience_rules -->

<!-- table:patience_shares -->
| sample | taken seed ends changed, of matched steps (Wilson 95%) | cut refinements whose end is under hit ratio 0.95 (Wilson 95%) |
|---|---|---|
| held | 41 of 3,640 (0.8% to 1.5%) | 420 of 763 (51.5% to 58.5%) |
| mid | 20 of 1,478 (0.9% to 2.1%) | 151 of 348 (38.3% to 48.6%) |
| grad | 1 of 1,280 (0.0% to 0.4%) | 3 of 15 (7.0% to 45.2%) |
| big | 4 of 281 (0.6% to 3.6%) | 23 of 55 (29.7% to 55.0%) |
| all four | n/a | 597 of 1,181 (47.7% to 53.4%) |
<!-- /table:patience_shares -->

<!-- table:patience_stages -->
| sample | arm | seed step | rescue | grid | neighbour | revisit | pass | total |
|---|---|---|---|---|---|---|---|---|
| held | G | 14,150,515 | 167,167 | 6,563,817 | 1,480,653 | 3,210,857 | 5,195,858 | 30,768,867 |
| held | S | 12,039,131 | 175,175 | 6,910,035 | 1,480,653 | 3,209,271 | 5,196,651 | 29,010,916 |
| mid | G | 5,417,024 | 32,032 | 1,334,072 | 1,488,949 | 1,678,659 | 3,053,721 | 13,004,457 |
| mid | S | 4,044,515 | 32,032 | 1,334,072 | 1,488,949 | 1,676,951 | 3,049,756 | 11,626,275 |
| grad | G | 1,528,284 | 34,034 | 1,697,698 | 1,477,542 | 1,779,614 | 3,679,520 | 10,196,692 |
| grad | S | 1,512,978 | 34,034 | 1,697,698 | 1,477,542 | 1,779,675 | 3,680,008 | 10,181,935 |
| big | G | 687,902 | 0 | 0 | 1,496,513 | 488,000 | 1,135,332 | 3,807,747 |
| big | S | 493,914 | 0 | 0 | 1,496,513 | 490,135 | 1,136,247 | 3,616,809 |
<!-- /table:patience_stages -->

<!-- table:patience_unmatched -->
| sample | seed steps G / S | matched | only in S | only in G | S-only with a cut | S-only wrong where G's map right | G-only wrong where S's map right | wrong seeds, matched + S-only (limit 1) |
|---|---|---|---|---|---|---|---|---|
| held | 3,705 / 3,699 | 3,640 | 59 | 65 | 16 | 0 | 0 | 0 |
| mid | 1,607 / 1,612 | 1,478 | 134 | 129 | 52 | 1 | 0 | 1 |
| grad | 1,280 / 1,284 | 1,280 | 4 | 0 | 0 | 0 | 0 | 0 |
| big | 317 / 316 | 281 | 35 | 36 | 15 | 0 | 0 | 0 |
<!-- /table:patience_unmatched -->

<!-- table:patience_validity -->
| sample | finished | voxels processed | aligned draws | seed steps G / S | matched steps | fits compared (matched steps) | fits differ | of them cut MC fits | other |
|---|---|---|---|---|---|---|---|---|---|
| held | True | 24,570 | True | 3,705 / 3,699 | 3,640 | 12,817 | 733 | 733 | 0 |
| mid | True | 24,570 | True | 1,607 / 1,612 | 1,478 | 3,920 | 280 | 280 | 0 |
| grad | True | 24,570 | True | 1,280 / 1,284 | 1,280 | 3,027 | 15 | 15 | 0 |
| big | True | 24,570 | True | 317 / 316 | 281 | 512 | 35 | 35 | 0 |
<!-- /table:patience_validity -->

<!-- table:patience_wrong -->
| sample | arm | failed seed tests | rescued | grid searches | wrong after pass | Wilson 95% (%) | median error (deg) |
|---|---|---|---|---|---|---|---|
| held | G | 65 | 46 | 19 | 69 | 0.22 to 0.36 | 0.02599 |
| held | S | 67 | 47 | 20 | 69 | 0.22 to 0.36 | 0.02601 |
| mid | G | 12 | 8 | 4 | 35 | 0.10 to 0.20 | 0.02861 |
| mid | S | 12 | 8 | 4 | 37 | 0.11 to 0.21 | 0.02859 |
| grad | G | 11 | 6 | 5 | 56 | 0.18 to 0.30 | 0.02134 |
| grad | S | 11 | 6 | 5 | 55 | 0.17 to 0.29 | 0.02134 |
| big | G | 0 | 0 | 0 | 21 | 0.06 to 0.13 | 0.03044 |
| big | S | 0 | 0 | 0 | 21 | 0.06 to 0.13 | 0.03044 |
<!-- /table:patience_wrong -->
