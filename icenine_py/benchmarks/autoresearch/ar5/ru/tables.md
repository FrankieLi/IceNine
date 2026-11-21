<!-- table:ru_gates -->
| gate | frozen bar | measured | verdict |
|---|---|---|---|
| r1 evaluations | at least 15% fewer than 23,402,740 | 28.9% fewer (16,648,743) | met |
| r2 mid | broken minus fixed at most 3; median at most +0.001; grid at most +2 | net -1; median -0.00044; grid +0 | met |
| r2 held | broken minus fixed at most 3; median at most +0.001; grid at most +2 | net -1; median -0.00067; grid +0 | met |
| r3 wall (contended) | at least 15% lower than 4,545 s | 31.7% lower (3,104 s) | met |
<!-- /table:ru_gates -->

<!-- table:ru_held8_arms -->
| arm | evals | wrong | median error (deg) | wall s | label |
|---|---|---|---|---|---|
| S8 | 58,102,056 | 49 | 0.02756 | 8390 | contended |
| PXF J5SD | 16,226,609 | 26 | 0.02267 | 2980 | contended |
| RU J5SD | 10,905,900 | 26 | 0.02206 | 2008 | contended |
<!-- /table:ru_held8_arms -->

<!-- table:ru_held8_gates -->
| gate | frozen bar | measured | verdict |
|---|---|---|---|
| h1 RU vs S8 | at least 20% fewer; net at most 5; median at most +0.005 | 81.2% fewer; fixed 30 broken 7; median -0.00549 | met |
| h2 RU vs PXF | at least 10% fewer; net at most 5; median at most +0.001 | 32.8% fewer; fixed 9 broken 9; median -0.00060 | met |
<!-- /table:ru_held8_gates -->

<!-- table:ru_held8_tries -->
| seed steps | ended by a try | Wilson 95% | failed tries | try evals |
|---|---|---|---|---|
| 3,532 | 2,564 (72.6%) | 71.1% to 74.0% | 967 | 1,024,404 |
<!-- /table:ru_held8_tries -->

<!-- table:ru_live_evals -->
| sample | PXF evals | RU evals | fewer | wrong PXF / RU | fixed | broken | median higher (deg) | grid PXF / RU |
|---|---|---|---|---|---|---|---|---|
| mid | 6,822,192 | 5,177,331 | 24.1% | 16 / 15 | 5 | 4 | -0.00044 | 0 / 0 |
| held | 16,580,548 | 11,471,412 | 30.8% | 40 / 39 | 11 | 10 | -0.00067 | 2 / 2 |
<!-- /table:ru_live_evals -->

<!-- table:ru_live_tries -->
| sample | seed steps | ended by a try | Wilson 95% | first fit / second fit | failed tries | try evals | indexing calls | try s | indexing s |
|---|---|---|---|---|---|---|---|---|---|
| mid | 1,719 | 1,403 (81.6%) | 79.7% to 83.4% | 1208 / 195 | 315 | 447,703 | 321 | 37 | 78 |
| held | 3,736 | 2,721 (72.8%) | 71.4% to 74.2% | 2122 / 599 | 1014 | 1,081,713 | 1,059 | 90 | 296 |
<!-- /table:ru_live_tries -->

<!-- table:ru_live_wall -->
| sample | PXF wall s | RU wall s | RU label | RU load mean / max | PXF label |
|---|---|---|---|---|---|
| mid | 1479 | 1001 | contended | 2.3 / 3.1 | contended |
| held | 3066 | 2103 | contended | 2.4 / 3.1 | contended |
<!-- /table:ru_live_wall -->
