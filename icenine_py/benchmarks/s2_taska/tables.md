<!-- table:s2_taska_cpp -->
| arm | path | n | peak_counts_equal | pixel_counts_equal | quality_within_1e6 |
|---|---|---|---|---|---|
| screen, CppEdgeRules | torch (Task A, stored) | 2230 | 2229 | n/a | n/a |
| screen, CppEdgeRules | torch (fb42d8c) | 2230 | 2229 | n/a | n/a |
| screen, CppEdgeRules | fused portable | 2230 | 2229 | n/a | n/a |
| screen, CppEdgeRules | fused accelerate | 2230 | 2229 | n/a | n/a |
| complete cost, CppEdgeRules | torch (Task A, stored) | 2230 | 2230 | 1985 | 2228 |
| complete cost, CppEdgeRules | torch (fb42d8c) | 2230 | 2230 | 1985 | 2228 |
| complete cost, CppEdgeRules | fused portable | 2230 | 2230 | 1984 | 2229 |
| complete cost, CppEdgeRules | fused accelerate | 2230 | 2230 | 1985 | 2228 |
<!-- /table:s2_taska_cpp -->

<!-- table:s2_taska_inputs -->
| check | n | equal |
|---|---|---|
| records: label, kind, FZ and grid indices equal | 2230 | 2230 |
| pair product, order plain: bit-equal flag and largest difference | 2000 | 1446 |
| pair product, order clang_contract: bit-equal flag and largest difference | 2000 | 2000 |
| pair product, order fma_chain: bit-equal flag and largest difference | 2000 | 1286 |
| pair product, order fma_rev: bit-equal flag and largest difference | 2000 | 1081 |
| torch at fb42d8c = stored Python, screen | 2230 | 2230 |
| torch at fb42d8c = stored Python, complete cost | 2230 | 2230 |
| torch at fb42d8c = stored Python, screen, CppEdgeRules | 2230 | 2230 |
| torch at fb42d8c = stored Python, complete cost, CppEdgeRules | 2230 | 2230 |
<!-- /table:s2_taska_inputs -->

<!-- table:s2_taska_outlier -->
| mode | orientation | row | coords | torch_max_err_px | fused_max_err_px | torch_eq_fused | int_flips_vs_f64 | asin_arg | omega_err_rad | coord0_f32_px | coord0_f64_px |
|---|---|---|---|---|---|---|---|---|---|---|---|
| edge | 735 | 83 | 6 | 0.0669 | 0.0669 | True | 0 | 0.99969 | 1.54e-04 | n/a | n/a |
| edge | 735 | 195 | 6 | 0.0658 | 0.0658 | True | 2 | 0.99969 | 1.53e-04 | n/a | n/a |
| edge | 947 | 219 | 1 | 1.5005 | 1.5005 | True | 1 | 0.11969 | 6.87e-08 | -5.25e-05 | 5.37e-04 |
<!-- /table:s2_taska_outlier -->

<!-- table:s2_taska_paired -->
| arm | fused | counts | both | torch_only | fused_only | neither |
|---|---|---|---|---|---|---|
| screen, CppEdgeRules | portable | peak counts | 2229 | 0 | 0 | 1 |
| screen, CppEdgeRules | accelerate | peak counts | 2229 | 0 | 0 | 1 |
| complete cost, CppEdgeRules | portable | peak counts | 2230 | 0 | 0 | 0 |
| complete cost, CppEdgeRules | portable | peak and pixel counts | 1934 | 51 | 50 | 195 |
| complete cost, CppEdgeRules | accelerate | peak counts | 2230 | 0 | 0 | 0 |
| complete cost, CppEdgeRules | accelerate | peak and pixel counts | 1966 | 19 | 19 | 226 |
<!-- /table:s2_taska_paired -->

<!-- table:s2_taska_vs_torch -->
| arm | fused | n | differing | within_rule | causes | worst_px | pixon | pixov | pkov | pkon | quality | f64_side_torch | f64_side_fused | silent | silent_within_rule | silent_causes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| screen, CppEdgeRules | portable | 2230 | 0 | 0 | - | n/a | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 13 | 13 | pixel_integer |
| screen, CppEdgeRules | accelerate | 2230 | 0 | 0 | - | n/a | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 10 | 10 | pixel_integer |
| complete cost, CppEdgeRules | portable | 2230 | 108 | 108 | pixel_integer | 0.0019 | 107 | 4 | 0 | 0 | 1 | 59 | 53 | 99 | 99 | pixel_integer |
| complete cost, CppEdgeRules | accelerate | 2230 | 43 | 43 | pixel_integer | 0.0009 | 42 | 3 | 0 | 0 | 0 | 25 | 18 | 37 | 37 | pixel_integer |
<!-- /table:s2_taska_vs_torch -->
