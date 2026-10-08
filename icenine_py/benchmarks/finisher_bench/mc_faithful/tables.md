<!-- table:mcf_SW_clean -->
| method | before | after | ev_before | ev_after |
|---|---|---|---|---|
| (i) MC deployed (200 steps) | 0.1623 / 13% / 6.4% | 0.0271 / 44% / 0.7% | 201 | 208 |
| (i) MC April sweep (3500 steps) @ 250 | 0.1165 / 10% / 4.5% | 0.0267 / 43% / 0.2% | 250 | 137 |
| (i) MC April sweep (3500 steps) @ 1000 | 0.1165 / 10% / 4.4% | 0.0267 / 43% / 0.2% | 1000 | 137 |
| (i) MC April sweep (3500 steps) @ 2600 | 0.1164 / 10% / 4.4% | 0.0267 / 43% / 0.2% | 2600 | 137 |
| (vi) VarianceMinimizing, quarter box @ 250 | 0.0757 / 23% / 1.5% | 0.1001 / 20% / 1.2% | 250 | 250 |
| (vi) VarianceMinimizing, quarter box @ 1000 | 0.0136 / 65% / 0.5% | 0.0313 / 38% / 0.0% | 1000 | 1000 |
| (vi) VarianceMinimizing, quarter box @ 2600 | 0.0102 / 80% / 0.4% | 0.0230 / 46% / 0.0% | 2600 | 2600 |
| (vi) VarianceMinimizing, quarter box @ 10000 | 0.0075 / 90% / 0.1% | 0.0160 / 57% / 0.0% | 10000 | 10000 |
| default finisher (MC + VM) | 0.0382 / 29% / 0.6% | 0.0197 / 50% / 0.0% | 1677 | 569 |
| (iv) Nelder-Mead @ 250 (unchanged) | 0.0031 / 85% / 4.5% | 0.0031 / 85% / 4.5% | 250 | 250 |
| (v) CMA-ES, sigma0 0.2 deg @ 250 (unchanged) | 0.0032 / 95% / 0.4% | 0.0032 / 95% / 0.4% | 250 | 250 |
| (v) CMA-ES, sigma0 0.2 deg @ 1000 (unchanged) | 0.0020 / 99% / 0.4% | 0.0020 / 99% / 0.4% | 1000 | 1000 |
<!-- /table:mcf_SW_clean -->

<!-- table:mcf_SW_realistic -->
| method | before | after | ev_before | ev_after |
|---|---|---|---|---|
| (i) MC deployed (200 steps) | 0.1609 / 11% / 5.8% | 0.0249 / 45% / 0.5% | 201 | 208 |
| (i) MC April sweep (3500 steps) @ 250 | 0.1154 / 9% / 4.8% | 0.0239 / 45% / 0.2% | 250 | 137 |
| (i) MC April sweep (3500 steps) @ 1000 | 0.1154 / 9% / 4.7% | 0.0239 / 45% / 0.2% | 1000 | 137 |
| (i) MC April sweep (3500 steps) @ 2600 | 0.1154 / 9% / 4.7% | 0.0239 / 45% / 0.2% | 2600 | 137 |
| (vi) VarianceMinimizing, quarter box @ 250 | 0.0731 / 23% / 1.3% | 0.0866 / 19% / 1.3% | 250 | 250 |
| (vi) VarianceMinimizing, quarter box @ 1000 | 0.0148 / 64% / 0.3% | 0.0311 / 37% / 0.1% | 1000 | 1000 |
| (vi) VarianceMinimizing, quarter box @ 2600 | 0.0113 / 78% / 0.3% | 0.0229 / 46% / 0.0% | 2600 | 2600 |
| (vi) VarianceMinimizing, quarter box @ 10000 | 0.0083 / 88% / 0.3% | 0.0158 / 57% / 0.0% | 10000 | 10000 |
| default finisher (MC + VM) | 0.0357 / 30% / 0.0% | 0.0184 / 52% / 0.0% | 1776 | 584 |
| (iv) Nelder-Mead @ 250 (unchanged) | 0.0034 / 86% / 3.1% | 0.0034 / 86% / 3.1% | 250 | 250 |
| (v) CMA-ES, sigma0 0.2 deg @ 250 (unchanged) | 0.0038 / 92% / 0.2% | 0.0038 / 92% / 0.2% | 250 | 250 |
| (v) CMA-ES, sigma0 0.2 deg @ 1000 (unchanged) | 0.0025 / 95% / 0.2% | 0.0025 / 95% / 0.2% | 1000 | 1000 |
<!-- /table:mcf_SW_realistic -->

<!-- table:mcf_T5-H0_clean -->
| method | before | after | ev_before | ev_after |
|---|---|---|---|---|
| (i) MC deployed (200 steps) | 0.6892 / 4% / 37.0% | 0.1538 / 19% / 33.0% | 201 | 208 |
| (i) MC April sweep (3500 steps) @ 250 | 0.3493 / 2% / 21.0% | 0.0730 / 32% / 8.0% | 250 | 137 |
| (i) MC April sweep (3500 steps) @ 1000 | 0.3493 / 2% / 20.0% | 0.0730 / 32% / 8.0% | 1000 | 137 |
| (i) MC April sweep (3500 steps) @ 2600 | 0.3493 / 2% / 20.0% | 0.0730 / 32% / 8.0% | 2600 | 137 |
| (vi) VarianceMinimizing, quarter box @ 250 | 0.3376 / 11% / 33.0% | 0.4994 / 14% / 33.0% | 250 | 250 |
| (vi) VarianceMinimizing, quarter box @ 1000 | 0.0261 / 41% / 25.0% | 0.2331 / 24% / 32.0% | 1000 | 1000 |
| (vi) VarianceMinimizing, quarter box @ 2600 | 0.0157 / 60% / 22.0% | 0.1557 / 29% / 31.0% | 2600 | 2600 |
| (vi) VarianceMinimizing, quarter box @ 10000 | 0.0110 / 65% / 21.0% | 0.1520 / 35% / 30.0% | 10000 | 10000 |
| default finisher (MC + VM) | 0.0582 / 16% / 17.0% | 0.0560 / 27% / 16.0% | 929 | 538 |
| (iv) Nelder-Mead @ 250 (unchanged) | 0.0053 / 57% / 34.0% | 0.0053 / 57% / 34.0% | 250 | 250 |
| (v) CMA-ES, sigma0 0.2 deg @ 250 (unchanged) | 0.0043 / 71% / 18.0% | 0.0043 / 71% / 18.0% | 250 | 250 |
| (v) CMA-ES, sigma0 0.2 deg @ 1000 (unchanged) | 0.0026 / 82% / 17.0% | 0.0026 / 82% / 17.0% | 1000 | 1000 |
<!-- /table:mcf_T5-H0_clean -->

<!-- table:mcf_T5-H0_realistic -->
| method | before | after | ev_before | ev_after |
|---|---|---|---|---|
| (i) MC deployed (200 steps) | 0.6635 / 5% / 35.0% | 0.1199 / 31% / 34.0% | 201 | 208 |
| (i) MC April sweep (3500 steps) @ 250 | 0.3451 / 3% / 21.0% | 0.0538 / 29% / 9.0% | 250 | 137 |
| (i) MC April sweep (3500 steps) @ 1000 | 0.3451 / 3% / 21.0% | 0.0538 / 29% / 9.0% | 1000 | 137 |
| (i) MC April sweep (3500 steps) @ 2600 | 0.3451 / 3% / 21.0% | 0.0538 / 29% / 9.0% | 2600 | 137 |
| (vi) VarianceMinimizing, quarter box @ 250 | 0.4143 / 9% / 32.0% | 0.4739 / 9% / 33.0% | 250 | 250 |
| (vi) VarianceMinimizing, quarter box @ 1000 | 0.0244 / 45% / 24.0% | 0.2579 / 20% / 32.0% | 1000 | 1000 |
| (vi) VarianceMinimizing, quarter box @ 2600 | 0.0170 / 63% / 20.0% | 0.1806 / 28% / 29.0% | 2600 | 2600 |
| (vi) VarianceMinimizing, quarter box @ 10000 | 0.0096 / 72% / 20.0% | 0.1378 / 32% / 29.0% | 10000 | 10000 |
| default finisher (MC + VM) | 0.0888 / 20% / 18.0% | 0.0724 / 28% / 15.0% | 616 | 538 |
| (iv) Nelder-Mead @ 250 (unchanged) | 0.0046 / 60% / 34.0% | 0.0046 / 60% / 34.0% | 250 | 250 |
| (v) CMA-ES, sigma0 0.2 deg @ 250 (unchanged) | 0.0052 / 74% / 15.0% | 0.0052 / 74% / 15.0% | 250 | 250 |
| (v) CMA-ES, sigma0 0.2 deg @ 1000 (unchanged) | 0.0026 / 83% / 15.0% | 0.0026 / 83% / 15.0% | 1000 | 1000 |
<!-- /table:mcf_T5-H0_realistic -->

<!-- table:mcf_T5-H3_clean -->
| method | before | after | ev_before | ev_after |
|---|---|---|---|---|
| (i) MC deployed (200 steps) | 0.0385 / 30% / 2.5% | 0.0166 / 58% / 2.5% | 201 | 193 |
| (i) MC April sweep (3500 steps) @ 250 | 0.0657 / 19% / 1.5% | 0.0390 / 31% / 1.5% | 58 | 52 |
| (i) MC April sweep (3500 steps) @ 1000 | 0.0657 / 20% / 1.5% | 0.0390 / 31% / 1.5% | 58 | 52 |
| (i) MC April sweep (3500 steps) @ 2600 | 0.0657 / 20% / 1.5% | 0.0390 / 31% / 1.5% | 58 | 52 |
| (vi) VarianceMinimizing, quarter box @ 250 | 0.0243 / 43% / 2.5% | 0.0268 / 40% / 3.0% | 250 | 250 |
| (vi) VarianceMinimizing, quarter box @ 1000 | 0.0144 / 65% / 1.0% | 0.0180 / 54% / 2.5% | 1000 | 1000 |
| (vi) VarianceMinimizing, quarter box @ 2600 | 0.0122 / 73% / 0.5% | 0.0156 / 62% / 2.0% | 2600 | 2600 |
| (vi) VarianceMinimizing, quarter box @ 10000 | 0.0092 / 88% / 0.5% | 0.0120 / 75% / 2.0% | 10000 | 10000 |
| default finisher (MC + VM) | 0.0224 / 46% / 1.0% | 0.0151 / 59% / 0.0% | 2167 | 540 |
| (iv) Nelder-Mead @ 250 (unchanged) | 0.0023 / 98% / 2.0% | 0.0023 / 98% / 2.0% | 250 | 250 |
| (v) CMA-ES, sigma0 0.2 deg @ 250 (unchanged) | 0.0026 / 96% / 1.5% | 0.0026 / 96% / 1.5% | 250 | 250 |
| (v) CMA-ES, sigma0 0.2 deg @ 1000 (unchanged) | 0.0019 / 98% / 1.5% | 0.0019 / 98% / 1.5% | 1000 | 1000 |
<!-- /table:mcf_T5-H3_clean -->

<!-- table:mcf_T5-H3_realistic -->
| method | before | after | ev_before | ev_after |
|---|---|---|---|---|
| (i) MC deployed (200 steps) | 0.0382 / 26% / 3.0% | 0.0175 / 55% / 2.5% | 201 | 193 |
| (i) MC April sweep (3500 steps) @ 250 | 0.0635 / 18% / 1.0% | 0.0365 / 32% / 0.0% | 59 | 52 |
| (i) MC April sweep (3500 steps) @ 1000 | 0.0635 / 19% / 1.0% | 0.0365 / 32% / 0.0% | 59 | 52 |
| (i) MC April sweep (3500 steps) @ 2600 | 0.0635 / 19% / 1.0% | 0.0365 / 32% / 0.0% | 59 | 52 |
| (vi) VarianceMinimizing, quarter box @ 250 | 0.0245 / 42% / 2.5% | 0.0279 / 40% / 2.5% | 250 | 250 |
| (vi) VarianceMinimizing, quarter box @ 1000 | 0.0150 / 66% / 1.0% | 0.0176 / 54% / 2.5% | 1000 | 1000 |
| (vi) VarianceMinimizing, quarter box @ 2600 | 0.0122 / 76% / 1.0% | 0.0153 / 64% / 2.0% | 2600 | 2600 |
| (vi) VarianceMinimizing, quarter box @ 10000 | 0.0096 / 86% / 0.5% | 0.0130 / 73% / 2.0% | 10000 | 10000 |
| default finisher (MC + VM) | 0.0229 / 44% / 1.5% | 0.0168 / 56% / 0.5% | 2629 | 560 |
| (iv) Nelder-Mead @ 250 (unchanged) | 0.0031 / 94% / 2.5% | 0.0031 / 94% / 2.5% | 250 | 250 |
| (v) CMA-ES, sigma0 0.2 deg @ 250 (unchanged) | 0.0031 / 94% / 2.0% | 0.0031 / 94% / 2.0% | 250 | 250 |
| (v) CMA-ES, sigma0 0.2 deg @ 1000 (unchanged) | 0.0022 / 96% / 2.0% | 0.0022 / 96% / 2.0% | 1000 | 1000 |
<!-- /table:mcf_T5-H3_realistic -->

<!-- table:mcf_deployed -->
| set | n | noimp_before | noimp_after | early_before | early_after | both_before | both_after | evals_med_before | evals_med_after |
|---|---|---|---|---|---|---|---|---|---|
| T5-H3 clean | 200 | 36.5% (30.1-43.4) | 31.5% (25.5-38.2) | 36.5% (30.1-43.4) | 54.5% (47.6-61.3) | 73 | 63 | 201 | 193 |
| T5-H0 clean | 100 | 4.0% (1.6-9.8) | 6.0% (2.8-12.5) | 4.0% (1.6-9.8) | 45.0% (35.6-54.8) | 4 | 6 | 201 | 208 |
| SW clean | 1000 | 3.8% (2.8-5.2) | 2.2% (1.5-3.3) | 3.8% (2.8-5.2) | 28.3% (25.6-31.2) | 38 | 22 | 201 | 208 |
| T5-H3 realistic | 200 | 37.0% (30.6-43.9) | 33.0% (26.9-39.8) | 37.0% (30.6-43.9) | 56.0% (49.1-62.7) | 74 | 66 | 201 | 193 |
| T5-H0 realistic | 100 | 2.0% (0.6-7.0) | 4.0% (1.6-9.8) | 2.0% (0.6-7.0) | 39.0% (30.0-48.8) | 2 | 4 | 201 | 208 |
| SW realistic | 1000 | 4.2% (3.1-5.6) | 2.6% (1.8-3.8) | 4.3% (3.2-5.7) | 29.9% (27.1-32.8) | 42 | 26 | 201 | 208 |
<!-- /table:mcf_deployed -->
