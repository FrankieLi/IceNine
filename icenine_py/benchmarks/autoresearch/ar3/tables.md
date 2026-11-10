<!-- table:ar3_accuracy -->
| arm | median | 90th | 99th | largest |
|---|---|---|---|---|
| paint (no fit) | 0.0169 | 0.0362 | 0.0594 | 0.0594 |
| key on (fitted) | 0.0159 | 0.0384 | 0.0577 | 0.1085 |
<!-- /table:ar3_accuracy -->

<!-- table:ar3_direct -->
| set | n | true grain passes | lowest cost is the true grain | highest hit is the true grain | lowest cost and highest hit differ |
|---|---|---|---|---|---|
| all: all kinds | 126 | 126 of 126 (97 to 100%) | 126 of 126 (97 to 100%) | 126 of 126 (97 to 100%) | 0 of 126 (0 to 3%) |
| all: edge | 63 | 63 of 63 (94 to 100%) | 63 of 63 (94 to 100%) | 63 of 63 (94 to 100%) | 0 of 63 (0 to 6%) |
| all: corner | 32 | 32 of 32 (89 to 100%) | 32 of 32 (89 to 100%) | 32 of 32 (89 to 100%) | 0 of 32 (0 to 11%) |
| all: interior | 31 | 31 of 31 (89 to 100%) | 31 of 31 (89 to 100%) | 31 of 31 (89 to 100%) | 0 of 31 (0 to 11%) |
| lowest cost and highest hit different orientations: all kinds | 44 | 44 of 44 (92 to 100%) | 44 of 44 (92 to 100%) | 44 of 44 (92 to 100%) | 0 of 44 (0 to 8%) |
| lowest cost and highest hit different orientations: edge | 44 | 44 of 44 (92 to 100%) | 44 of 44 (92 to 100%) | 44 of 44 (92 to 100%) | 0 of 44 (0 to 8%) |
<!-- /table:ar3_direct -->

<!-- table:ar3_evals -->
| arm | total | seed paths | scoring | local fits | grid searches | wrong |
|---|---|---|---|---|---|---|
| key off (stored rerun) | 25,547,935 | 0 | 0 | 2,665,070 | 22,882,865 | 37 of 2000 (1.3 to 2.5%) |
| key on (stored) | 2,819,906 | 122,198 | 0 | 2,697,708 | 0 | 22 of 2000 (0.7 to 1.7%) |
| key on, FusedCost (stored) | 2,816,897 | 122,154 | 0 | 2,694,743 | 0 | 22 of 2000 (0.7 to 1.7%) |
| paint (FusedCost) | 156,766 | 119,033 | 6,697 | 31,036 | 0 | 0 of 2000 (0.0 to 0.2%) |
| paint + fit the contested (side arm) | 324,092 | 119,033 | 6,697 | 198,362 | 0 | 0 of 2000 (0.0 to 0.2%) |
<!-- /table:ar3_evals -->

<!-- table:ar3_flood -->
| distance from the grain's true region (voxel sides) | n |
|---|---|
| 0 to 1 | 44 |
| 1 to 2 | 0 |
| 2 to 3 | 0 |
| 3 to 4 | 0 |
| 4 or more | 0 |
<!-- /table:ar3_flood -->

<!-- table:ar3_margin -->
| fit if the two lowest costs differ by less than | voxels | fits | evals | fixed | broken | wrong |
|---|---|---|---|---|---|---|
| 0.001 | 8 | 16 | 11,094 | 0 | 0 | 0 |
| 0.003 | 12 | 24 | 16,390 | 0 | 0 | 0 |
| 0.01 | 20 | 40 | 26,982 | 0 | 0 | 0 |
| 0.03 | 66 | 132 | 87,886 | 0 | 0 | 0 |
| 0.1 | 120 | 240 | 159,382 | 0 | 0 | 0 |
| 0.3 | 126 | 252 | 167,326 | 0 | 0 | 0 |
| 1 | 126 | 252 | 167,326 | 0 | 0 | 0 |
| all | 126 | 252 | 167,326 | 0 | 0 | 0 |
<!-- /table:ar3_margin -->

<!-- table:ar3_paired -->
| vs | wrong there | wrong in paint | fixed | broken | wrong in both |
|---|---|---|---|---|---|
| key on (stored) | 22 | 0 | 22 | 0 | 0 |
| key on, FusedCost (stored) | 22 | 0 | 22 | 0 | 0 |
| key off (stored) | 37 | 0 | 37 | 0 | 0 |
<!-- /table:ar3_paired -->

<!-- table:ar3_premise -->
| scope | n_grains | n_voxels | max_spread_deg | n_grains_over_1e3_deg | n_grains_over_1_deg |
|---|---|---|---|---|---|
| region | 60 | 2000 | 0.00000 | 0 | 0 |
| layer | 497 | 24570 | 0.00000 | 0 | 0 |
<!-- /table:ar3_premise -->

<!-- table:ar3_profile -->
| part | s | calls or evaluations | share |
|---|---|---|---|
| local fits, cost evaluations | 121.2 | 2,694,743 | 54.1% |
| local fits, optimizer's own work | 65.5 | 4,095 | 29.2% |
| seed path: indexing | 11.9 | 64 | 5.3% |
| seed path: finisher, cost evaluations | 4.3 | 121,834 | 1.9% |
| seed path: finisher, optimizer's own work | 3.1 | 128 | 1.4% |
| seed path: strict costs on candidates and the rest of the seed path | 0.0 | 64 | 0.0% |
| overlap evaluations (the BFS's own) | 0.0 | 64 | 0.0% |
| grid search | 0.0 | 0 | 0.0% |
| one-time pairing (inside the first seed path) | 17.6 | n/a | 7.8% |
| everything else (BFS bookkeeping, neighbours, records) | 0.5 | n/a | 0.2% |
<!-- /table:ar3_profile -->

<!-- table:ar3_wall -->
| what | s |
|---|---|
| paint run (FusedCost, unbatched loop) | 37.9 |
|   of which the one-time pairing | 17.2 |
| key on, FusedCost, timers only, this session | 224.3 |
| key on, FusedCost, stored (AR2 step 2) | 213.1 |
| key on under cProfile (not a timing) | 235.6 |
| side arm (fits on the contested voxels) | 11.0 |
<!-- /table:ar3_wall -->

<!-- table:ar3_wrong -->
| what | n |
|---|---|
| wrong voxels of paint | 0 |
|   of which crossing a boundary, sharing an edge | 0 |
|   of which crossing a boundary, corner only | 0 |
|   of which interior voxels | 0 |
<!-- /table:ar3_wrong -->
