<!-- table:stallhit_changed -->
| sample | seed steps G / H | matched | only in H | only in G | taken ends changed (matched) | median angle to G (deg) | changed, now at 0.95 or more | changed, max error H (deg) | H wrong where G right: matched + only-H |
|---|---|---|---|---|---|---|---|---|---|
| held2 | 3,631 / 3,632 | 3,600 | 32 | 31 | 12 | 0.0622 | 12 | 0.0986 | 0 + 0 |
| held | 3,705 / 3,700 | 3,687 | 13 | 18 | 3 | 0.0192 | 3 | 0.0634 | 0 + 0 |
| mid | 1,607 / 1,607 | 1,582 | 25 | 25 | 7 | 0.0331 | 7 | 0.1303 | 0 + 0 |
| grad | 1,280 / 1,280 | 1,280 | 0 | 0 | 0 | 0.0000 | 0 | 0.0000 | 0 + 0 |
| big | 317 / 317 | 317 | 0 | 0 | 1 | 0.0317 | 1 | 0.0449 | 0 + 0 |
<!-- /table:stallhit_changed -->

<!-- table:stallhit_cuts -->
| sample | MC refinements | stall checks | stopped (cuts) | went on | hit ratio of the checks that went on (min / median / max) | lowest hit ratio at a stop | evaluations in cut ones |
|---|---|---|---|---|---|---|---|
| held2 | 5,652 | 833 | 357 | 476 | 0.446 / 0.902 / 0.949 | 0.9500 | 1,324,060 |
| held | 5,443 | 789 | 343 | 446 | 0.491 / 0.903 / 0.950 | 0.9503 | 1,266,396 |
| mid | 1,606 | 346 | 197 | 149 | 0.672 / 0.911 / 0.949 | 0.9502 | 781,191 |
| grad | 991 | 15 | 12 | 3 | 0.884 / 0.917 / 0.926 | 0.9514 | 37,253 |
| big | 153 | 45 | 27 | 18 | 0.525 / 0.895 / 0.950 | 0.9617 | 117,288 |
<!-- /table:stallhit_cuts -->

<!-- table:stallhit_fallbacks -->
| sample | arm | failed seed tests | rescued | grid searches |
|---|---|---|---|---|
| held2 | G | 64 | 43 | 21 |
| held2 | H | 64 | 43 | 21 |
| held | G | 65 | 46 | 19 |
| held | H | 65 | 46 | 19 |
| mid | G | 12 | 8 | 4 |
| mid | H | 12 | 8 | 4 |
| grad | G | 11 | 6 | 5 |
| grad | H | 11 | 6 | 5 |
| big | G | 0 | 0 | 0 |
| big | H | 0 | 0 | 0 |
<!-- /table:stallhit_fallbacks -->

<!-- table:stallhit_layer -->
| sample | arm | evaluations with the pass | seed step | wrong after pass | Wilson 95% (%) | median error (deg) | fixed / broken vs G | fixed / broken vs S8 |
|---|---|---|---|---|---|---|---|---|
| held2 | G | 32,554,555 | 15,269,080 | 67 | 0.21 to 0.35 | 0.02648 | n/a | 19 / 15 |
| held2 | H | 31,143,763 | 13,862,375 | 67 | 0.21 to 0.35 | 0.02649 | 0 / 0 | 19 / 15 |
| held | G | 30,768,867 | 14,150,515 | 69 | 0.22 to 0.36 | 0.02599 | n/a | 25 / 24 |
| held | H | 29,588,170 | 12,971,953 | 69 | 0.22 to 0.36 | 0.02600 | 0 / 0 | 25 / 24 |
| mid | G | 13,004,457 | 5,417,024 | 35 | 0.10 to 0.20 | 0.02861 | n/a | 9 / 13 |
| mid | H | 11,950,397 | 4,363,940 | 35 | 0.10 to 0.20 | 0.02861 | 0 / 0 | 9 / 13 |
| grad | G | 10,196,692 | 1,528,284 | 56 | 0.18 to 0.30 | 0.02134 | n/a | 31 / 24 |
| grad | H | 10,180,191 | 1,511,783 | 56 | 0.18 to 0.30 | 0.02134 | 0 / 0 | 31 / 24 |
| big | G | 3,807,747 | 687,902 | 21 | 0.06 to 0.13 | 0.03044 | n/a | 8 / 5 |
| big | H | 3,635,256 | 515,411 | 21 | 0.06 to 0.13 | 0.03045 | 0 / 0 | 8 / 5 |
<!-- /table:stallhit_layer -->

<!-- table:stallhit_paired -->
| sample | pass | fixed vs G | broken vs G | both wrong vs G | fixed vs S8 | broken vs S8 | both wrong vs S8 |
|---|---|---|---|---|---|---|---|
| held2 | before | 0 | 0 | 97 | 28 | 17 | 80 |
| held2 | after | 0 | 0 | 67 | 19 | 15 | 52 |
| held | before | 0 | 0 | 101 | 26 | 32 | 69 |
| held | after | 0 | 0 | 69 | 25 | 24 | 45 |
| mid | before | 0 | 0 | 75 | 22 | 13 | 62 |
| mid | after | 0 | 0 | 35 | 9 | 13 | 22 |
| grad | before | 0 | 0 | 206 | 63 | 59 | 147 |
| grad | after | 0 | 0 | 56 | 31 | 24 | 32 |
| big | before | 0 | 0 | 74 | 10 | 17 | 57 |
| big | after | 0 | 0 | 21 | 8 | 5 | 16 |
<!-- /table:stallhit_paired -->

<!-- table:stallhit_predictions -->
| prediction | held2 | held | mid | grad | big |
|---|---|---|---|---|---|
| 1 replay (gate) | n/a | n/a | True (8 changed) | True (0 changed) | True (1 changed) |
| 2 replay saving | n/a | n/a | True (7.86%) | True (0.16%) | True (4.53%) |
| 3 total fewer (mid 6-9, big 3.5-5.5, grad -0.5..0.5, held 2-6, held2 2-8; held2 seed 6-18) | True (4.33%) | True (3.84%) | True (8.11%) | True (0.16%) | True (4.53%) |
| 4 broken - fixed vs G in -3..3, median shift <= 0.0005 | True (0) | True (0) | True (0) | True (0) | True (0) |
| 5 mid passes rule 5 | n/a | n/a | True | n/a | n/a |
<!-- /table:stallhit_predictions -->

<!-- table:stallhit_rules -->
| sample | (1) total fewer vs G (%) | (1) | (2) broken - fixed vs G | (2) | (3) median shift (deg) | (3) | (4) extra wrong seeds (matched + unmatched); grid H/G | (4) | (5) vs S8: fewer %; broken - fixed; median shift | (5) | G passes (5) | H passes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| held2 | 4.33; seed 9.21 | pass | 0 | pass | 0.00000 | pass | 0; 21/21 | pass | 49.1; -4; -0.00131 | pass | yes | pass |
| held | 3.84 | pass | 0 | pass | 0.00001 | pass | 0; 19/19 | pass | 50.6; -1; -0.00259 | pass | yes | pass |
| mid | 8.11 | pass | 0 | pass | -0.00000 | pass | 0; 4/4 | pass | 54.8; 4; -0.00111 | pass | yes | pass |
| grad | 0.16 | pass | 0 | pass | 0.00000 | pass | 0; 5/5 | pass | 43.4; -7; -0.00129 | pass | yes | pass |
| big | 4.53 | pass | 0 | pass | 0.00001 | pass | 0; 0/0 | pass | 56.4; -3; -0.00116 | pass | yes | pass |
<!-- /table:stallhit_rules -->

<!-- table:stallhit_shares -->
| sample | taken seed ends changed, of matched steps (Wilson 95%) |
|---|---|
| held2 | 12 of 3,600 (0.2% to 0.6%) |
| held | 3 of 3,687 (0.0% to 0.2%) |
| mid | 7 of 1,582 (0.2% to 0.9%) |
| grad | 0 of 1,280 (0.0% to 0.3%) |
| big | 1 of 317 (0.1% to 1.8%) |
<!-- /table:stallhit_shares -->

<!-- table:stallhit_stages -->
| sample | arm | seed step | rescue | grid | neighbour | revisit | pass | total |
|---|---|---|---|---|---|---|---|---|
| held2 | G | 15,269,080 | 176,176 | 7,252,187 | 1,480,531 | 3,147,539 | 5,229,042 | 32,554,555 |
| held2 | H | 13,862,375 | 176,176 | 7,252,187 | 1,480,592 | 3,144,916 | 5,227,517 | 31,143,763 |
| held | G | 14,150,515 | 167,167 | 6,563,817 | 1,480,653 | 3,210,857 | 5,195,858 | 30,768,867 |
| held | H | 12,971,953 | 167,167 | 6,563,817 | 1,480,653 | 3,209,027 | 5,195,553 | 29,588,170 |
| mid | G | 5,417,024 | 32,032 | 1,334,072 | 1,488,949 | 1,678,659 | 3,053,721 | 13,004,457 |
| mid | H | 4,363,940 | 32,032 | 1,334,072 | 1,488,949 | 1,679,147 | 3,052,257 | 11,950,397 |
| grad | G | 1,528,284 | 34,034 | 1,697,698 | 1,477,542 | 1,779,614 | 3,679,520 | 10,196,692 |
| grad | H | 1,511,783 | 34,034 | 1,697,698 | 1,477,542 | 1,779,614 | 3,679,520 | 10,180,191 |
| big | G | 687,902 | 0 | 0 | 1,496,513 | 488,000 | 1,135,332 | 3,807,747 |
| big | H | 515,411 | 0 | 0 | 1,496,513 | 488,000 | 1,135,332 | 3,635,256 |
<!-- /table:stallhit_stages -->

<!-- table:stallhit_validity -->
| sample | finished | voxels processed | aligned draws | fits compared (matched steps) | fits differ | cut MC fits | only one counted check | other | every stop at 0.95 or more |
|---|---|---|---|---|---|---|---|---|---|
| held2 | True | 24,570 | True | 13,242 | 780 | 353 | 427 | 0 | True |
| held | True | 24,570 | True | 12,949 | 746 | 342 | 404 | 0 | True |
| mid | True | 24,570 | True | 4,210 | 312 | 188 | 124 | 0 | True |
| grad | True | 24,570 | True | 3,027 | 15 | 12 | 3 | 0 | True |
| big | True | 24,570 | True | 574 | 43 | 27 | 16 | 0 | True |
<!-- /table:stallhit_validity -->
