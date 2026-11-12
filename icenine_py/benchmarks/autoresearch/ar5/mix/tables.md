<!-- table:accuracy -->
| sample | run | wrong before the pass | 95% interval | wrong after the pass | 95% interval | S8 before the pass | S8 after the pass | S10 before the pass | S10 after the pass |
|---|---|---|---|---|---|---|---|---|---|
| big | S11 | 70 | 0.23 to 0.36% | 18 | 0.05 to 0.12% | 67 | 24 | 77 | 22 |
| first gradient sample | S11 | 207 | 0.74 to 0.96% | 62 | 0.20 to 0.32% | 210 | 63 | 198 | 67 |
| mid | S11 | 77 | 0.25 to 0.39% | 36 | 0.11 to 0.20% | 84 | 31 | 73 | 35 |
<!-- /table:accuracy -->

<!-- table:check -->
| seed steps checked | with MC run on a candidate | identical to Part 1 | evaluations, script layer | evaluations, Part 1 |
|---|---|---|---|---|
| 20 | 10 | 20 | 121659 | 121659 |
<!-- /table:check -->

<!-- table:fallback_cost -->
| sample | S8 grid-fallback evaluations | S8 fallbacks | S8 per fallback | S10 grid-fallback evaluations | S10 fallbacks | S10 per fallback | S11 grid-fallback evaluations | S11 fallbacks | S11 per fallback |
|---|---|---|---|---|---|---|---|---|---|
| big | 0 | 0 | n/a | 0 | 0 | n/a | 0 | 0 | n/a |
| first gradient sample | 5449501 | 16 | 340594 | 5108107 | 15 | 340540 | 3741408 | 11 | 340128 |
| mid | 5002750 | 15 | 333517 | 5008366 | 15 | 333891 | 4004585 | 12 | 333715 |
<!-- /table:fallback_cost -->

<!-- table:nogrid -->
| sample | S8 without grid fallbacks | S10 without grid fallbacks | S11 without grid fallbacks | S11 fewer than S8 (%) | S11 fewer than S10 (%) | grid-fallback evaluations S11 saves against S8 | all evaluations S11 saves against S8 | grid-fallback evaluations S11 saves against S10 | all evaluations S11 saves against S10 |
|---|---|---|---|---|---|---|---|---|---|
| big | 8340428 | 5898258 | 4015342 | 51.9 | 31.9 | 0 | 4325086 | 0 | 1882916 |
| first gradient sample | 12550729 | 10059891 | 8566586 | 31.7 | 14.8 | 1708093 | 5692236 | 1366699 | 2860004 |
| mid | 21440324 | 17660338 | 12078111 | 43.7 | 31.6 | 998165 | 10360378 | 1003781 | 6586008 |
<!-- /table:nogrid -->

<!-- table:paired -->
| sample | run | fixed vs S8, before the pass | broken vs S8, before the pass | fixed vs S8, after the pass | broken vs S8, after the pass | wrong in both, after the pass | fixed vs S10, before the pass | broken vs S10, before the pass | fixed vs S10, after the pass | broken vs S10, after the pass |
|---|---|---|---|---|---|---|---|---|---|---|
| big | S11 | 10 | 13 | 13 | 7 | 11 | 16 | 9 | 9 | 5 |
| first gradient sample | S11 | 54 | 51 | 30 | 29 | 33 | 53 | 62 | 39 | 34 |
| mid | S11 | 25 | 18 | 10 | 15 | 21 | 24 | 28 | 12 | 13 |
<!-- /table:paired -->

<!-- table:part1 -->
| sample | arm | steps | seed-test passes | passes with M only | passes with the arm only | ends within 1 degree | M only | arm only | seed-step evaluations | fewer than M (%) | candidates where MC ran (%) | steps with a missing MC end |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| big | M | 329 | 329 | n/a | n/a | 329 | n/a | n/a | 2856712 | 0.0 | 100.0 | 0 |
| big | X60 | 329 | 329 | 0 | 0 | 329 | 0 | 0 | 2251212 | 21.2 | 83.9 | 0 |
| big | X200 | 329 | 329 | 0 | 0 | 329 | 0 | 0 | 938875 | 67.1 | 37.5 | 0 |
| first gradient sample | M | 1279 | 1263 | n/a | n/a | 1254 | n/a | n/a | 3008757 | 0.0 | 100.0 | 0 |
| first gradient sample | X60 | 1279 | 1265 | 0 | 2 | 1255 | 0 | 1 | 2138986 | 28.9 | 74.1 | 1 |
| first gradient sample | X200 | 1279 | 1267 | 0 | 4 | 1260 | 1 | 7 | 1563116 | 48.0 | 49.3 | 0 |
| mid | M | 1611 | 1596 | n/a | n/a | 1595 | n/a | n/a | 12007352 | 0.0 | 100.0 | 0 |
| mid | X60 | 1611 | 1596 | 0 | 0 | 1595 | 0 | 0 | 10600545 | 11.7 | 90.2 | 0 |
| mid | X200 | 1611 | 1600 | 0 | 4 | 1599 | 0 | 4 | 5643575 | 53.0 | 59.2 | 0 |
<!-- /table:part1 -->

<!-- table:pick -->
| arm | seed-test passes | ends within 1 degree | seed-step evaluations | fewer than M (%) | candidates where MC ran (%) | steps with a missing MC end | rule (a) | rule (b) |
|---|---|---|---|---|---|---|---|---|
| M | 1592 | 1583 | 5865469 | 0.0 | 100.0 | 0 | n/a | n/a |
| X60 | 1594 | 1584 | 4390198 | 25.2 | 75.7 | 1 | True | True |
| X200 | 1596 | 1589 | 2501991 | 57.3 | 47.3 | 0 | True | True |
<!-- /table:pick -->

<!-- table:predictions -->
| sample | S11 total with the pass | prediction 2: below | S11 grid fallbacks | S8 grid fallbacks | S11 seed-step evaluations | Part 1 count, same sample | difference (%) | candidates where MC ran (%) | candidates refined |
|---|---|---|---|---|---|---|---|---|---|
| big | 4015342 | 4400000 | 0 | 0 | 888360 | 938875 | -5.4 | 39.0 | 436 |
| first gradient sample | 12307994 | 14200000 | 11 | 16 | 1626738 | 1563116 | 4.1 | 48.8 | 2158 |
| mid | 16082696 | 17500000 | 12 | 15 | 5825184 | 5643575 | 3.2 | 60.7 | 2814 |
<!-- /table:predictions -->

<!-- table:rule -->
| sample | run | broken minus fixed vs S8 (after the pass) | median error, S11 (deg) | median error, S8 (deg) | difference (deg) | rule: at least 20% fewer evaluations | rule: broken minus fixed at most 5 | rule: median at most 0.005 degree above | passes | fewer evaluations than S8 (%) |
|---|---|---|---|---|---|---|---|---|---|---|
| big | S11 | -6 | 0.03077 | 0.03161 | -0.00084 | True | True | True | True | 51.9 |
| first gradient sample | S11 | -1 | 0.02104 | 0.02264 | -0.00160 | True | True | True | True | 31.6 |
| mid | S11 | 5 | 0.02913 | 0.02972 | -0.00060 | True | True | True | True | 39.2 |
<!-- /table:rule -->

<!-- table:shares -->
| sample | against | evaluations S11 saves | seed step (%) | grid fallback (%) | neighbour fits (%) | revisits (%) | boundary pass (%) |
|---|---|---|---|---|---|---|---|
| big | S8 | 4325086 | 45.5 | 0.0 | 46.1 | 8.3 | 0.1 |
| first gradient sample | S8 | 5692236 | 24.3 | 30.0 | 23.2 | 22.6 | -0.0 |
| mid | S8 | 10360378 | 59.7 | 9.6 | 22.9 | 7.7 | 0.1 |
| big | S10 | 1882916 | 101.0 | 0.0 | -0.0 | -0.8 | -0.2 |
| first gradient sample | S10 | 2860004 | 52.3 | 47.8 | -0.0 | 0.1 | -0.2 |
| mid | S10 | 6586008 | 85.1 | 15.2 | -0.0 | -0.4 | 0.1 |
<!-- /table:shares -->

<!-- table:stages -->
| sample | run | seed step | grid fallback | neighbour fits | revisits | before the pass | boundary pass | total with the pass | fewer than S8 (%) | fewer than S10 (%) | grid fallbacks (count) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| big | S8 | 2856712 | 0 | 3491437 | 853714 | 7201863 | 1138565 | 8340428 | 0.0 | n/a | 0 |
| big | S10 | 2790857 | 0 | 1496574 | 479704 | 4767135 | 1131123 | 5898258 | 29.3 | 0.0 | 0 |
| big | S11 | 888360 | 0 | 1496635 | 494710 | 2879705 | 1135637 | 4015342 | 51.9 | 31.9 | 0 |
| first gradient sample | S8 | 3008757 | 5449501 | 2798001 | 3060730 | 14316989 | 3683241 | 18000230 | 0.0 | n/a | 16 |
| first gradient sample | S10 | 3122117 | 5108107 | 1477481 | 1779065 | 11486770 | 3681228 | 15167998 | 15.7 | 0.0 | 15 |
| first gradient sample | S11 | 1626738 | 3741408 | 1477542 | 1776503 | 8622191 | 3685803 | 12307994 | 31.6 | 18.9 | 11 |
| mid | S8 | 12007352 | 5002750 | 3859661 | 2517272 | 23387035 | 3056039 | 26443074 | 0.0 | n/a | 15 |
| mid | S10 | 11431628 | 5008366 | 1489254 | 1687748 | 19616996 | 3051708 | 22668704 | 14.3 | 0.0 | 15 |
| mid | S11 | 5825184 | 4004585 | 1489620 | 1715137 | 13034526 | 3048170 | 16082696 | 39.2 | 29.1 | 12 |
<!-- /table:stages -->

<!-- table:uncounted -->
| sample | uncounted hit-ratio checks (one per mixed candidate) | counted evaluations of S11 with the pass | uncounted as a share of counted (%) |
|---|---|---|---|
| big | 436 | 4015342 | 0.0109 |
| first gradient sample | 2158 | 12307994 | 0.0175 |
| mid | 2814 | 16082696 | 0.0175 |
<!-- /table:uncounted -->

<!-- table:whatif -->
| sample | replay step | candidate | CMA end hit ratio | CMA end cost | CMA evaluations | MC end hit ratio | MC end cost | MC evaluations | next candidate: CMA evaluations | next candidate: CMA end hit ratio | X200 seed-step evaluations, this page | the plan's what-if |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mid | 442 | 1 | 0.94752 | 0.06939 | 201 | 0.95300 | 0.07189 | 5983 | 201 | 0.979 | 5643575 | 5643374 |
<!-- /table:whatif -->
