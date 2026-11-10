<!-- table:interp_by_j -->
| j | local_rotation_deg | n | same_pass_pct | same_counts_pct | omega_err_median_of_max_deg | omega_err_worst_deg | pct_with_exact_peak_missing | pct_with_extra_peak |
|---|---|---|---|---|---|---|---|---|
| 0 | 0.000 | 102 | 100.0 | 100.0 | 0.0000 | 0.0000 | 0.0 | 0.0 |
| 1 | 2.837 | 92 | 96.7 | 93.5 | 0.0896 | 6.5947 | 2.2 | 1.1 |
| 2 | 2.837 | 98 | 100.0 | 98.0 | 0.0813 | 174.5162 | 3.1 | 2.0 |
| 3 | 2.837 | 109 | 99.1 | 98.2 | 0.1278 | 169.9071 | 3.7 | 2.8 |
| 4 | 2.837 | 94 | 96.8 | 93.6 | 0.1174 | 170.1070 | 2.1 | 2.1 |
| 5 | 2.837 | 134 | 99.3 | 96.3 | 0.0945 | 174.5162 | 2.2 | 2.2 |
| 6 | 2.837 | 104 | 100.0 | 95.2 | 0.1089 | 133.9467 | 4.8 | 3.8 |
| 7 | 2.837 | 73 | 97.3 | 97.3 | 0.1066 | 18.5853 | 1.4 | 0.0 |
| 8 | 2.837 | 94 | 100.0 | 94.7 | 0.1040 | 170.1794 | 6.4 | 3.2 |
<!-- /table:interp_by_j -->

<!-- table:interp_by_theta -->
| theta_deg | n | same_pass_pct | same_counts_pct | same_counts_lo | same_counts_hi | omega_err_median_of_max_deg | omega_err_p90_of_max_deg | omega_err_worst_deg | pct_with_exact_peak_missing | pct_with_extra_peak |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.05 | 75 | 100.0 | 100.0 | 95.1 | 100.0 | 0.0000 | 0.0005 | 0.602 | 0.0 | 0.0 |
| 0.10 | 75 | 100.0 | 100.0 | 95.1 | 100.0 | 0.0001 | 0.0011 | 0.631 | 0.0 | 0.0 |
| 0.20 | 75 | 100.0 | 100.0 | 95.1 | 100.0 | 0.0004 | 0.0033 | 0.074 | 0.0 | 0.0 |
| 0.50 | 75 | 100.0 | 98.7 | 92.8 | 99.8 | 0.0021 | 0.0818 | 174.464 | 1.3 | 1.3 |
| 1.00 | 75 | 100.0 | 100.0 | 95.1 | 100.0 | 0.0113 | 0.1912 | 7.303 | 0.0 | 0.0 |
| 2.00 | 75 | 100.0 | 98.7 | 92.8 | 99.8 | 0.0428 | 1.3502 | 170.192 | 1.3 | 1.3 |
| 3.00 | 75 | 98.7 | 96.0 | 88.9 | 98.6 | 0.1182 | 0.4056 | 173.865 | 2.7 | 1.3 |
| 5.00 | 75 | 100.0 | 97.3 | 90.8 | 99.3 | 0.3495 | 1.6521 | 8.550 | 4.0 | 4.0 |
<!-- /table:interp_by_theta -->

<!-- table:parity_by_group -->
| group | n | EXS | EXS_edge | EXF | EXF_edge |
|---|---|---|---|---|---|
| truth | 3 | 2 | 3 | 3 | 3 |
| off_0.1 | 9 | 6 | 9 | 5 | 9 |
| off_1.0 | 9 | 6 | 9 | 3 | 9 |
| off_3.0 | 9 | 9 | 9 | 4 | 9 |
| random | 200 | 178 | 200 | 128 | 200 |
| pair_random | 900 | 815 | 900 | 553 | 900 |
| sweep | 600 | 547 | 600 | 382 | 600 |
| pair_differ_pass | 442 | 356 | 441 | 261 | 442 |
| pair_differ_count | 58 | 46 | 58 | 32 | 58 |
<!-- /table:parity_by_group -->

<!-- table:parity_causes -->
| arm | peaks | frame_below_range | frame_boundary | pixel | counts | omega_only |
|---|---|---|---|---|---|---|
| screen, Python default | 108196 | 283 | 1 | 0 | 0 | 32487 |
| screen, Python with CppEdgeRules | 108196 | 0 | 1 | 0 | 0 | 32586 |
| complete cost, Python default | 463012 | 856 | 1 | 269 | 311 | 137918 |
| complete cost, Python with CppEdgeRules | 463012 | 0 | 1 | 259 | 0 | 138368 |
<!-- /table:parity_causes -->

<!-- table:parity_exact -->
| arm | n | cost | counts | frames | qual |
|---|---|---|---|---|---|
| screen, Python default | 2230 | 1965 (88.1%) | 1965 (88.1%) | 1846 (82.8%) | n/a |
| screen, Python with CppEdgeRules | 2230 | 2229 (99.96%) | 2229 (99.96%) | 2187 (98.1%) | n/a |
| complete cost, Python default | 2230 | 1371 (61.5%) | 1124 (50.4%) | 887 (39.8%) | 2020 (90.6%) |
| complete cost, Python with CppEdgeRules | 2230 | 2230 (100.0%) | 1985 (89.0%) | 1690 (75.8%) | 2228 (99.9%) |
<!-- /table:parity_exact -->

<!-- table:pass_patterns -->
| group | pattern | n |
|---|---|---|
| pair_random | ---- | 664 |
| pair_random | PPPP | 226 |
| pair_random | P--- | 6 |
| pair_random | -PPP | 4 |
| pair_differ_pass | P--- | 214 |
| pair_differ_pass | -PPP | 187 |
| pair_differ_pass | --P- | 40 |
| pair_differ_pass | PP-- | 1 |
| pair_differ_count | PPPP | 58 |
<!-- /table:pass_patterns -->
