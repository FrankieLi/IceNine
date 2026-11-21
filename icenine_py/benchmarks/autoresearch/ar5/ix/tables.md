<!-- table:ix_e1 -->
| sample | voxels | identical | differ | old_s | new_s | less_pct |
|---|---|---|---|---|---|---|
| mid | 1626 | 1626 | 0 | 963.0 | 606.0 | 37.1 |
| held | 3685 | 3685 | 0 | 4036.4 | 2624.8 | 35.0 |
| held3 | 300 | 300 | 0 | 380.9 | 241.2 | 36.7 |
| held4 | 300 | 300 | 0 | 472.8 | 290.8 | 38.5 |
| held5 | 300 | 300 | 0 | 465.5 | 283.3 | 39.1 |
| held6 | 300 | 300 | 0 | 393.3 | 252.8 | 35.7 |
| grad | 300 | 300 | 0 | 135.0 | 93.4 | 30.8 |
| big | 300 | 300 | 0 | 147.4 | 81.0 | 45.1 |
<!-- /table:ix_e1 -->

<!-- table:ix_parts -->
| sample | part | old_s | new_s | less_pct |
|---|---|---|---|---|
| mid | vote | 50.3 | 50.0 | 0.6 |
| mid | cluster | 783.1 | 493.0 | 37.0 |
| mid | refine | 123.2 | 56.7 | 54.0 |
| mid | merge | 3.8 | 3.7 | 1.9 |
| held | vote | 215.6 | 215.9 | -0.1 |
| held | cluster | 3421.5 | 2202.1 | 35.6 |
| held | refine | 381.2 | 189.1 | 50.4 |
| held | merge | 8.7 | 8.6 | 1.0 |
<!-- /table:ix_parts -->

<!-- table:ix_w1 -->
| sample | voxels | old_s | new_s | less_pct | label |
|---|---|---|---|---|---|
| mid | 1626 | 963.0 | 606.0 | 37.1 | contended |
| held | 3685 | 4036.4 | 2624.8 | 35.0 | contended |
<!-- /table:ix_w1 -->

<!-- table:ix_w2 -->
| sample | evals_old | evals_new | map_new_equals_stored | wall_old | wall_new | less_pct | label |
|---|---|---|---|---|---|---|---|
| mid | 7,383,368 | 7,383,368 | True | 1674 | 1345 | 19.7 | contended / contended |
| held | 17,114,589 | 17,114,589 | True | 5632 | 4282 | 24.0 | contended / contended |
<!-- /table:ix_w2 -->
