<!-- table:px_gates -->
| gate | measured | limit | verdict |
|---|---|---|---|
| q1 | 6 of 2396 voxels lost from the top 20, 0 gained | at most 0.5% lost | met |
| c1 | pairs 73.9% fewer, votes 73.3% fewer (new live runs against the stored seed voxels indexed again) | at least 60% fewer | met |
| c2 | 25,033,884 against 24,497,957 stored (+2.19%) | at most 2% more | not met |
| c3 mid | broken minus fixed -6, median error change +0.00012 deg, grid searches -2 | at most 3; at most 0.001 deg; at most 2 more | met |
| c3 held | broken minus fixed -3, median error change -0.00027 deg, grid searches +3 | at most 3; at most 0.001 deg; at most 2 more | not met |
| t1 | indexing seconds 75.8% lower (contended) | at least 50% lower | met |
| t2 | pooled wall 44.9% lower (8219 s to 4531 s; contended) | at least 25% lower | met |
| h1 | held7 not built (c2 and c3 not met) |  | not evaluated |
<!-- /table:px_gates -->

<!-- table:px_live -->
| sample | evaluations, stored | old rerun | new | new vs stored | wrong, stored | wrong, old rerun | wrong, new | fixed | broken | median error change (deg) | grid searches vs stored |
|---|---|---|---|---|---|---|---|---|---|---|---|
| mid | 7383368 | 7383368 | 6801802 | -7.88% | 23 | 23 | 17 | 10 | 4 | +0.00012 | -2 |
| held | 17114589 | 17114589 | 18232082 | +6.53% | 45 | 45 | 42 | 19 | 16 | -0.00027 | 3 |
<!-- /table:px_live -->

<!-- table:px_q1 -->
| sample | voxels | in top 20, base | in top 20, new | lost | gained | rank 0, base | rank 0, new | pairs fewer |
|---|---|---|---|---|---|---|---|---|
| mid | 300 | 300 | 300 | 0 | 0 | 236 | 240 | 49.8% |
| held | 300 | 300 | 299 | 1 | 0 | 250 | 254 | 78.7% |
| held3 | 300 | 300 | 300 | 0 | 0 | 256 | 255 | 80.0% |
| held4 | 300 | 300 | 300 | 0 | 0 | 256 | 250 | 82.1% |
| held5 | 300 | 299 | 296 | 3 | 0 | 236 | 245 | 81.5% |
| held6 | 300 | 299 | 298 | 1 | 0 | 252 | 251 | 81.0% |
| grad | 300 | 300 | 299 | 1 | 0 | 227 | 226 | 69.6% |
| big | 300 | 298 | 298 | 0 | 0 | 223 | 223 | 14.9% |
| all | 2400 | 2396 | 2390 | 6 | 0 | 1936 | 1944 |  |
<!-- /table:px_q1 -->

<!-- table:px_t1 -->
| sample | voxels | old s | new s | seconds less | pairs less | votes less | machine |
|---|---|---|---|---|---|---|---|
| mid | 1626 | 951 | 417 | 56.1% | 50.8% | 50.7% | contended |
| held | 3685 | 4.01e+03 | 784 | 80.4% | 78.5% | 78.9% | contended |
| pooled | 5311 | 4.96e+03 | 1.2e+03 | 75.8% | 74.0% | 73.5% |  |
<!-- /table:px_t1 -->

<!-- table:px_wall -->
| sample | pairs, old | pairs, new | votes, old | votes, new | indexing s, old | indexing s, new | wall s, old | wall s, new | wall less | preflight verdict old / new | mean load old / new |
|---|---|---|---|---|---|---|---|---|---|---|---|
| mid | 163107708 | 82222776 | 85534715 | 43311536 | 953 | 432 | 2020 | 1470 | 27.2% | contended / contended | 2.4 / 2.4 |
| held | 855657758 | 183923622 | 362140868 | 76413880 | 4009 | 794 | 6199 | 3061 | 50.6% | contended / contended | 2.4 / 2.4 |
<!-- /table:px_wall -->
