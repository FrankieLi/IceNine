<!-- table:pxf_gates -->
| gate | measured | limit | verdict |
|---|---|---|---|
| c1 | pairs 72.0% fewer, votes 71.3% fewer (PXF live runs, second tries included, against the stored seed voxels indexed again) | at least 50% fewer | met |
| c2 | 23,402,740 against 24,497,957 stored (-4.47%) | at most 2% more | met |
| c3 mid | broken minus fixed -7, median error change +0.00010 deg, grid searches -2 | at most 3; at most 0.001 deg; at most 2 more | met |
| c3 held | broken minus fixed -5, median error change -0.00020 deg, grid searches -2 | at most 3; at most 0.001 deg; at most 2 more | met |
| t2 | pooled wall 44.7% lower (8219 s to 4545 s; preflight verdicts in the wall table) | at least 25% lower | met |
| h1 | held7: 74.2% fewer evaluations than S8, broken minus fixed -54, median error change -0.00545 deg | at least 20% fewer; at most 5; at most 0.005 deg | met |
<!-- /table:pxf_gates -->

<!-- table:pxf_h1 -->
| arm | evaluations | fewer than S8 | wrong voxels | fixed vs S8 | broken vs S8 | median error change vs S8 (deg) | grid searches | wall s | preflight verdict | mean load |
|---|---|---|---|---|---|---|---|---|---|---|
| S8 | 77783446 |  | 96 |  |  |  |  | 11953 | contended | 2.4 |
| new J5SD (PXF) | 20075769 | 74.2% | 42 | 65 | 11 | -0.00545 | 6 | 3791 | contended | 2.4 |
| plain J5SD | 21788476 | 72.0% | 55 | 50 | 9 | -0.00523 | 12 | 9181 | contended | 2.4 |
<!-- /table:pxf_h1 -->

<!-- table:pxf_live -->
| sample | evaluations, stored | evaluations, PXF | PXF vs stored | wrong, stored | wrong, PXF | fixed | broken | median error change (deg) | grid searches, stored | grid searches, PXF |
|---|---|---|---|---|---|---|---|---|---|---|
| mid | 7383368 | 6822192 | -7.60% | 23 | 16 | 11 | 4 | +0.00010 | 2 | 0 |
| held | 17114589 | 16580548 | -3.12% | 45 | 40 | 20 | 15 | -0.00020 | 4 | 2 |
<!-- /table:pxf_live -->

<!-- table:pxf_steps -->
| sample | seed steps | failed the seed test at 50,000 | second tries | second try passed | sent to the rescue | evaluations of second tries | indexing seconds of second tries | second tries equal to the stored old step |
|---|---|---|---|---|---|---|---|---|
| mid | 1667 | 9 | 9 | 5 | 4 | 47406 | 7 | 9 of 9 |
| held | 3669 | 63 | 63 | 42 | 21 | 289213 | 112 | 60 of 60 |
<!-- /table:pxf_steps -->

<!-- table:pxf_wall -->
| sample | pairs, old | pairs, PXF | votes, old | votes, PXF | indexing s, old | indexing s, PXF | wall s, old | wall s, PXF | wall less | preflight verdict old / PXF | mean load old / PXF |
|---|---|---|---|---|---|---|---|---|---|---|---|
| mid | 163107708 | 83535544 | 85534715 | 44023434 | 953 | 441 | 2020 | 1479 | 26.8% | contended / contended | 2.4 / 2.5 |
| held | 855657758 | 201741534 | 362140868 | 84571595 | 4009 | 902 | 6199 | 3066 | 50.5% | contended / contended | 2.4 / 2.3 |
<!-- /table:pxf_wall -->
