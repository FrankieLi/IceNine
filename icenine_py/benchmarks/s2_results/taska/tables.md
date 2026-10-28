<!-- table:s2_lib_taska_cpp -->
| arm | path | n | peak_counts_equal | pixel_counts_equal | quality_within_1e6 |
|---|---|---|---|---|---|
| screen, CppEdgeRules | torch (Task A, stored) | 2230 | 2229 | n/a | n/a |
| screen, CppEdgeRules | torch (41efb45) | 2230 | 2229 | n/a | n/a |
| screen, CppEdgeRules | fused library | 2230 | 2229 | n/a | n/a |
| complete cost, CppEdgeRules | torch (Task A, stored) | 2230 | 2230 | 1985 | 2228 |
| complete cost, CppEdgeRules | torch (41efb45) | 2230 | 2230 | 1985 | 2228 |
| complete cost, CppEdgeRules | fused library | 2230 | 2230 | 1983 | 2228 |
<!-- /table:s2_lib_taska_cpp -->

<!-- table:s2_lib_taska_inputs -->
| check | n | equal |
|---|---|---|
| records: label, kind, FZ and grid indices equal | 2230 | 2230 |
| pair product, order plain: bit-equal flag and largest difference | 2000 | 1446 |
| pair product, order clang_contract: bit-equal flag and largest difference | 2000 | 2000 |
| pair product, order fma_chain: bit-equal flag and largest difference | 2000 | 1286 |
| pair product, order fma_rev: bit-equal flag and largest difference | 2000 | 1081 |
| torch at 41efb45 = stored Python, screen | 2230 | 2230 |
| torch at 41efb45 = stored Python, complete cost | 2230 | 2230 |
| torch at 41efb45 = stored Python, screen, CppEdgeRules | 2230 | 2230 |
| torch at 41efb45 = stored Python, complete cost, CppEdgeRules | 2230 | 2230 |
<!-- /table:s2_lib_taska_inputs -->

<!-- table:s2_lib_taska_paired -->
| arm | fused | counts | both | torch_only | fused_only | neither |
|---|---|---|---|---|---|---|
| screen, CppEdgeRules | library | peak counts | 2229 | 0 | 0 | 1 |
| complete cost, CppEdgeRules | library | peak counts | 2230 | 0 | 0 | 0 |
| complete cost, CppEdgeRules | library | peak and pixel counts | 1950 | 35 | 33 | 212 |
<!-- /table:s2_lib_taska_paired -->

<!-- table:s2_lib_taska_vs_torch -->
| arm | fused | n | differing | within_rule | causes | worst_px | pixon | pixov | pkov | pkon | quality | f64_side_torch | f64_side_fused | silent | silent_within_rule | silent_causes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| screen, CppEdgeRules | library | 2230 | 0 | 0 | - | n/a | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 14 | 14 | pixel_integer |
| complete cost, CppEdgeRules | library | 2230 | 73 | 73 | pixel_integer | 0.0010 | 72 | 3 | 0 | 0 | 0 | 39 | 35 | 80 | 80 | pixel_integer |
<!-- /table:s2_lib_taska_vs_torch -->
