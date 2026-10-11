<!-- table:e0_fairness_broken -->
| setting | variant | n | off_wrong | on_wrong | broken | fixed | both_wrong | cpp_s0_wrong_on_broken | cpp_s1_wrong_on_broken | cpp_s2_wrong_on_broken |
|---|---|---|---|---|---|---|---|---|---|---|
| min_sin_eta 0.3 | clean | 200 | 40 | 61 | 38 | 17 | 23 | 24 | 13 | 17 |
| min_sin_eta 0.3 | realistic | 200 | 25 | 42 | 25 | 8 | 17 | 11 | 10 | 16 |
| min_sin_eta 0 | clean | 200 | 39 | 56 | 34 | 17 | 22 | 18 | 13 | 19 |
| min_sin_eta 0 | realistic | 200 | 18 | 56 | 44 | 6 | 12 | 22 | 19 | 22 |
<!-- /table:e0_fairness_broken -->

<!-- table:e0_fairness_counts -->
| variant | n | py_off | py_on | py_off_ms0 | py_on_ms0 | cpp_s0 | cpp_s1 | cpp_s2 | cpp_mean | cpp_wrong_all_seeds | cpp_wrong_any_seed |
|---|---|---|---|---|---|---|---|---|---|---|---|
| clean | 200 | 40 | 61 | 39 | 56 | 71 | 56 | 66 | 64.3 | 28 | 105 |
| realistic | 200 | 25 | 42 | 18 | 56 | 52 | 46 | 52 | 50.0 | 25 | 83 |
<!-- /table:e0_fairness_counts -->
