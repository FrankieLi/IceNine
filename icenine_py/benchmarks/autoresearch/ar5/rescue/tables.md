<!-- table:accuracy -->
Reading line: wrong voxels out of 24,570 (more than 1 degree from the truth), before and after the boundary pass, and the median error after the pass in degrees.

| sample | arm | wrong before pass (Wilson 95%) | wrong after pass (Wilson 95%) | median error after pass |
|---|---|---|---|---|
| held | S8 | 95 (0.32% to 0.47%) | 70 (0.23% to 0.36%) | 0.02859 |
| held | A | 104 (0.35% to 0.51%) | 73 (0.24% to 0.37%) | 0.02607 |
| held | G | 101 (0.34% to 0.50%) | 69 (0.22% to 0.36%) | 0.02599 |
| mid | S8 | 84 (0.28% to 0.42%) | 31 (0.09% to 0.18%) | 0.02972 |
| mid | A | 75 (0.24% to 0.38%) | 35 (0.10% to 0.20%) | 0.02862 |
| mid | G | 75 (0.24% to 0.38%) | 35 (0.10% to 0.20%) | 0.02861 |
| grad | S8 | 210 (0.75% to 0.98%) | 63 (0.20% to 0.33%) | 0.02264 |
| grad | A | 206 (0.73% to 0.96%) | 56 (0.18% to 0.30%) | 0.02135 |
| grad | G | 206 (0.73% to 0.96%) | 56 (0.18% to 0.30%) | 0.02134 |
| big | S8 | 67 (0.21% to 0.35%) | 24 (0.07% to 0.15%) | 0.03161 |
| big | A | 74 (0.24% to 0.38%) | 21 (0.06% to 0.13%) | 0.03044 |
| big | G | 74 (0.24% to 0.38%) | 21 (0.06% to 0.13%) | 0.03044 |
<!-- /table:accuracy -->

<!-- table:evals -->
Reading line: evaluations by stage, run plus boundary pass, arm A (stored step 8) and arm G; `fewer` is G's total against A's. Wall times are not given (contended).

| sample | A total | G total | fewer | A seed step | G seed step | G rescue | A grid | G grid | neighbour+revisit (G) | pass A / G |
|---|---|---|---|---|---|---|---|---|---|---|
| held | 46,516,396 | 30,768,867 | 33.9% | 14,173,739 | 14,150,515 | 167,167 | 22,458,034 | 6,563,817 | 4,691,510 | 5,192,137 / 5,195,858 |
| mid | 15,645,472 | 13,004,457 | 16.9% | 5,417,024 | 5,417,024 | 32,032 | 4,007,241 | 1,334,072 | 3,167,608 | 3,053,599 / 3,053,721 |
| grad | 12,189,777 | 10,196,692 | 16.4% | 1,527,842 | 1,528,284 | 34,034 | 3,726,723 | 1,697,698 | 3,257,156 | 3,679,276 / 3,679,520 |
| big | 3,807,747 | 3,807,747 | 0.0% | 687,902 | 687,902 | 0 | 0 | 0 | 1,984,513 | 1,135,332 / 1,135,332 |
<!-- /table:evals -->

<!-- table:paired -->
Reading line: the base is the arm named after `vs`; fixed = wrong in the base and right in the first-named arm, broken = the reverse, counted by voxel.

| sample | pair | before pass fixed / broken | after pass fixed / broken | both wrong after | taken end differs (matched steps) |
|---|---|---|---|---|---|
| held | G vs A | 3 / 0 | 4 / 0 | 69 | 0 of 3702 |
| held | G vs S8 | 26 / 32 | 25 / 24 | 45 | - |
| held | A vs S8 | 26 / 35 | 25 / 28 | 45 | - |
| mid | G vs A | 0 / 0 | 0 / 0 | 35 | 0 of 1607 |
| mid | G vs S8 | 22 / 13 | 9 / 13 | 22 | - |
| mid | A vs S8 | 22 / 13 | 9 / 13 | 22 | - |
| grad | G vs A | 0 / 0 | 0 / 0 | 56 | 0 of 1277 |
| grad | G vs S8 | 63 / 59 | 31 / 24 | 32 | - |
| grad | A vs S8 | 63 / 59 | 31 / 24 | 32 | - |
| big | G vs A | 0 / 0 | 0 / 0 | 21 | 0 of 317 |
| big | G vs S8 | 10 / 17 | 8 / 5 | 16 | - |
| big | A vs S8 | 10 / 17 | 8 / 5 | 16 | - |
<!-- /table:paired -->

<!-- table:predictions -->
Reading line: the five predictions written before any run.

| prediction | met | numbers |
|---|---|---|
| 1 | met | cost check: grid end cheaper than every kept end in 23 of 23 (needed 21) |
| 2 | NOT met | replay: rescued 14 of 23 (12 to 21 met), all right (met), median error 0.0030 (met), mean tries 2.87 (at most 2.5 NOT met) |
| 3 | met | held: grid fallbacks 65 to 19 (70.8% fewer, predicted 40 to 85), total 33.9% fewer (predicted 20 to 45) |
| 4 | met | mid 16.9% fewer (12 to 21), grad 16.4% (15 to 28), big identical: True |
| 5 | NOT met | broken minus fixed against A (after pass) held -4, mid 0, grad 0, big 0 (-3 to 5); median error moves held -0.00008, mid -0.00000, grad -0.00001, big +0.00000 (at most 0.001); held's -4 is outside the range, on the favourable side |
<!-- /table:predictions -->

<!-- table:rescues -->
Reading line: what rule G1 did at the failed seed steps of arm G (errors in degrees against the truth, scoring only; wrong means over 1 degree).

| sample | failed seed steps (A / G) | rescues tried | rescued | tries | grid fallbacks left | rescued hit ratio min | rescued error median / max | wrong rescued seeds |
|---|---|---|---|---|---|---|---|---|
| held | 65 / 65 | 65 | 46 | 167 | 19 | 0.963 | 0.0016 / 0.0069 | 0 |
| mid | 12 / 12 | 12 | 8 | 32 | 4 | 0.967 | 0.0035 / 0.0108 | 0 |
| grad | 11 / 11 | 11 | 6 | 34 | 5 | 0.984 | 0.0027 / 0.0062 | 0 |
| big | 0 / 0 | 0 | 0 | 0 | 0 | - | - | 0 |
<!-- /table:rescues -->

<!-- table:rules -->
Reading line: the plan's five frozen rules for arm G against arm A (and the house rule against S8), per sample; a rule that is not met says NO.

| sample | (1) fallbacks / total | (2) net vs A | (3) median vs A | (4) wrong rescued | (5) house vs S8 (A passes it?) | G passes |
|---|---|---|---|---|---|---|
| held | yes: grid 65 to 19, total 33.9% fewer | yes: -4 | yes: -0.00008 | yes: 0 | yes (A: yes) | PASS |
| mid | yes: grid 12 to 4, total 16.9% fewer | yes: 0 | yes: -0.00000 | yes: 0 | yes (A: yes) | PASS |
| grad | yes: grid 11 to 5, total 16.4% fewer | yes: 0 | yes: -0.00001 | yes: 0 | yes (A: yes) | PASS |
| big | yes: grid 0 to 0, total 0.0% fewer | yes: 0 | yes: +0.00000 | yes: 0 | yes (A: yes) | PASS |
<!-- /table:rules -->

<!-- table:verdict -->
Reading line: the plan's verdict: the step passes if G passes on held and on all disclosed samples; if it fails on held there is no key.

| sample | G passes |
|---|---|
| held | PASS |
| mid | PASS |
| grad | PASS |
| big | PASS |

Step verdict: PASS
<!-- /table:verdict -->
