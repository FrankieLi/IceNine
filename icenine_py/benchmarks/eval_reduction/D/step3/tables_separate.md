<!-- table:step3_separation -->
| wrong orientations | score | median best wrong / truth | max over voxels | voxels >= 0.5 truth | voxels >= truth |
|---|---|---|---|---|---|
| random axis, 1 degree and more | pixel overlap (C++ rule) | 0.110 | 0.173 | 0 | 0 |
| random axis, 1 degree and more | attributed, A r=15.7 um | 0.178 | 1.035 | 5 | 1 |
| random axis, 1 degree and more | attributed, A r=30 um | 0.138 | 0.983 | 2 | 0 |
| random axis, 1 degree and more | attributed, B footprint | 0.139 | 0.260 | 0 | 0 |
| random axis, 1 degree and more | attributed, oracle (true pairs) | 0.139 | 0.260 | 0 | 0 |
| random axis, 0.25 to 0.5 degree | pixel overlap (C++ rule) | 0.640 | 0.818 | 28 | 0 |
| random axis, 0.25 to 0.5 degree | attributed, A r=15.7 um | 0.962 | 1.035 | 23 | 11 |
| random axis, 0.25 to 0.5 degree | attributed, A r=30 um | 0.849 | 1.526 | 30 | 3 |
| random axis, 0.25 to 0.5 degree | attributed, B footprint | 0.830 | 0.921 | 30 | 0 |
| random axis, 0.25 to 0.5 degree | attributed, oracle (true pairs) | 0.830 | 0.921 | 30 | 0 |
| rotations about low-index crystal axes | pixel overlap (C++ rule) | 0.310 | 0.350 | 0 | 0 |
| rotations about low-index crystal axes | attributed, A r=15.7 um | 0.316 | 1.000 | 4 | 1 |
| rotations about low-index crystal axes | attributed, A r=30 um | 0.309 | 0.667 | 2 | 0 |
| rotations about low-index crystal axes | attributed, B footprint | 0.310 | 0.435 | 0 | 0 |
| rotations about low-index crystal axes | attributed, oracle (true pairs) | 0.307 | 0.435 | 0 | 0 |
<!-- /table:step3_separation -->

<!-- table:step3_sharpness -->
| misorientation (deg) | pixel overlap (C++ rule) | attributed, A r=15.7 um | attributed, A r=30 um | attributed, B footprint | attributed, oracle (true pairs) |
|---|---|---|---|---|---|
| 0.05 | 0.927 | 0.976 | 0.967 | 0.958 | 0.958 |
| 0.1 | 0.819 | 0.944 | 0.922 | 0.921 | 0.920 |
| 0.25 | 0.512 | 0.764 | 0.732 | 0.731 | 0.733 |
| 0.5 | 0.235 | 0.338 | 0.338 | 0.342 | 0.342 |
| 1 | 0.061 | 0.072 | 0.085 | 0.087 | 0.088 |
| 2 | 0.017 | 0.012 | 0.017 | 0.020 | 0.018 |
| 5 | 0.004 | 0.000 | 0.003 | 0.003 | 0.003 |
| 10 | 0.001 | 0.000 | 0.000 | 0.000 | 0.000 |
| 20 | 0.003 | 0.000 | 0.000 | 0.000 | 0.000 |
<!-- /table:step3_sharpness -->

<!-- table:step3_tolerance -->
| frame window | direction tolerance (deg) | B completeness at truth (median) | oracle purity at truth (median) | B half-score radius (deg) |
|---|---|---|---|---|
| 0 | 0.02 | 0.736 | 0.885 | 0.184 |
| 0 | 0.05 | 0.792 | 0.954 | 0.377 |
| 0 | 0.1 | 0.812 | 0.979 | 0.593 |
| 0 | 0.2 | 0.814 | 0.991 | 0.717 |
| 1 | 0.02 | 0.736 | 0.885 | 0.217 |
| 1 | 0.05 | 0.792 | 0.954 | 0.555 |
| 1 | 0.1 | 0.812 | 0.979 | 1.003 |
| 1 | 0.2 | 0.821 | 0.991 | 1.528 |
<!-- /table:step3_tolerance -->

<!-- table:step3_truth_scores -->
| score | median at the truth | minimum at the truth | half-score radius (deg) |
|---|---|---|---|
| pixel overlap (C++ rule) | 0.936 | 0.894 | 0.258 |
| attributed, A r=15.7 um | 0.078 | 0.000 | 0.384 |
| attributed, A r=30 um | 0.719 | 0.018 | 0.376 |
| attributed, B footprint | 0.792 | 0.696 | 0.377 |
| attributed, oracle (true pairs) | 0.789 | 0.696 | 0.378 |
<!-- /table:step3_truth_scores -->
