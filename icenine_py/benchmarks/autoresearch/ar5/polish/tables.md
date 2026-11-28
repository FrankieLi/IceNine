<!-- table:polish_ends -->
| sample | G2 accepted (of failed rescues; Wilson 95%) | accepted ends wrong | largest error of an accepted end (degree) | median angle to H's map (degree) | largest angle to H's map (degree) | G2 evaluations |
|---|---|---|---|---|---|---|
| held2 | 12 of 21 (37% to 76%) | 0 | 0.01553 | 0.0042 | 42.585 | 60,517 |
| held | 15 of 19 (57% to 91%) | 0 | 0.00915 | 0.0076 | 0.722 | 42,960 |
| mid | 2 of 4 (15% to 85%) | 0 | 0.00303 | 0.0193 | 0.024 | 12,856 |
| grad | 3 of 5 (23% to 88%) | 0 | 0.00235 | 0.0050 | 0.552 | 12,856 |
| big | 0 of 0 | 0 | n/a | n/a | n/a | 0 |
| all five | 32 of 49 (51% to 77%) | 0 | 0.01553 | 0.0064 | 42.585 | 129,189 |
<!-- /table:polish_ends -->

<!-- table:polish_evals -->
| sample | J total with the pass | H total with the pass | % fewer (J) | J before the pass | H before the pass | % fewer before the pass | grid searches left (H) | G2 accepted |
|---|---|---|---|---|---|---|---|---|
| held2 | waits (pass) | 31,143,763 | waits | 21,834,061 | 25,916,246 | 15.75 | 9 (H 21) | 12 of 21 |
| held | 24,457,773 | 29,588,170 | 17.34 | 19,259,536 | 24,392,617 | 21.04 | 4 (H 19) | 15 of 19 |
| mid | 11,295,158 | 11,950,397 | 5.48 | 8,242,901 | 8,898,140 | 7.36 | 2 (H 4) | 2 of 4 |
| grad | 9,175,867 | 10,180,191 | 9.87 | 5,495,005 | 6,500,671 | 15.47 | 2 (H 5) | 3 of 5 |
| big | 3,635,256 | 3,635,256 | 0.00 | 2,499,924 | 2,499,924 | 0.00 | 0 (H 0) | 0 of 0 |
<!-- /table:polish_evals -->

<!-- table:polish_held2_now -->
| held2 number | value | status |
|---|---|---|
| run evaluations before the pass, J | 21,834,061 | final |
| run evaluations before the pass, H | 25,916,246 | final |
| % fewer before the pass | 15.75 | final |
| grid searches left, J (H) | 9 (H 21) | final |
| G2 accepted | 12 of 21 | final |
| wrong voxels before the pass, J | 95 (0.32% to 0.47%) | before the pass only |
| wrong voxels before the pass, H | 97 (0.32% to 0.48%) | before the pass only |
| fixed / broken by J vs H, before the pass | 2 / 0 | before the pass only |
| extra wrong seeds (matched + unmatched; the seed step's own end) | 0 | final |
| total with the pass, wrong after the pass, rules 1, 2, 3, 5, predictions 4 and 5 | not computed | waits for step 13's held2 result |
<!-- /table:polish_held2_now -->

<!-- table:polish_house -->
| sample | % fewer evaluations than S8 (with the pass) | fixed vs S8 | broken vs S8 | median error minus S8's (degree) |
|---|---|---|---|---|
| held | 59.2 | 25 | 23 | -0.00262 |
| mid | 57.3 | 9 | 13 | -0.00111 |
| grad | 49.0 | 31 | 24 | -0.00130 |
| big | 56.4 | 8 | 5 | -0.00116 |
<!-- /table:polish_house -->

<!-- table:polish_rules -->
| sample | (1) total with the pass | (2) broken minus fixed at most 5 | (3) median error | (4) wrong seeds (the seed step's own end, before the fallback; G2 cannot change it) | G2 ends accepted, wrong, largest error (degree) | (5) house rule against S8 | J passes |
|---|---|---|---|---|---|---|---|
| held2 | waits | waits | waits | True (seeds only) | 12 accepted, 0 wrong, largest 0.01553 | waits | waits for step 13 |
| held | True (17.34% fewer) | True | True | True | 15 accepted, 0 wrong, largest 0.00915 | True | True |
| mid | True (5.48% fewer) | True | True | True | 2 accepted, 0 wrong, largest 0.00303 | True | True |
| grad | True (9.87% fewer) | True | True | True | 3 accepted, 0 wrong, largest 0.00235 | True | True |
| big | True (0.00% fewer) | True | True | True | 0 accepted, 0 wrong, largest n/a | True | True |
<!-- /table:polish_rules -->

<!-- table:polish_validity -->
| sample | J finished | aligned draws | seed steps | fits that differ in matched steps | grid searches both arms ran with equal evaluations |
|---|---|---|---|---|---|
| held2 | True | True | 3633 (H 3632); matched 3610 | 0 | 9 of 9 |
| held | True | True | 3700 (H 3700); matched 3700 | 0 | 4 of 4 |
| mid | True | True | 1607 (H 1607); matched 1607 | 0 | 2 of 2 |
| grad | True | True | 1281 (H 1280); matched 1280 | 0 | 2 of 2 |
| big | True | True | 317 (H 317); matched 317 | 0 | 0 of 0 |
<!-- /table:polish_validity -->

<!-- table:polish_wrong -->
| sample | wrong after the pass, J (Wilson 95%) | wrong, H | wrong, S8 | fixed by J vs H | broken by J vs H | wrong in both | median error J (degree) | median error H | extra wrong seeds (matched + unmatched; the seed step's own end) |
|---|---|---|---|---|---|---|---|---|---|
| held | 68 (0.22% to 0.35%) | 69 (0.22% to 0.36%) | 70 (0.23% to 0.36%) | 1 | 0 | 68 | 0.02597 | 0.02600 | 0 |
| mid | 35 (0.10% to 0.20%) | 35 (0.10% to 0.20%) | 31 (0.09% to 0.18%) | 0 | 0 | 35 | 0.02861 | 0.02861 | 0 |
| grad | 56 (0.18% to 0.30%) | 56 (0.18% to 0.30%) | 63 (0.20% to 0.33%) | 0 | 0 | 56 | 0.02134 | 0.02134 | 0 |
| big | 21 (0.06% to 0.13%) | 21 (0.06% to 0.13%) | 24 (0.07% to 0.15%) | 0 | 0 | 21 | 0.03045 | 0.03045 | 0 |
| held2 | waits | waits | waits | waits | waits | waits | waits | waits | 0 |
<!-- /table:polish_wrong -->
