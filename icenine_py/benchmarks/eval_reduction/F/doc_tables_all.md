<!-- table:evalF_e0_acc -->
| method | 0.1 deg | 0.5 deg | 1 deg | 2 deg | 3 deg | 5 deg |
|---|---|---|---|---|---|---|
| F4 Gauss-Newton on E4 rays | 480 accepted, 0 of them wrong; 0 of 0 rejected were right | 480 accepted, 0 of them wrong; 0 of 0 rejected were right | 480 accepted, 0 of them wrong; 0 of 0 rejected were right | 480 accepted, 0 of them wrong; 0 of 0 rejected were right | 478 accepted, 0 of them wrong; 0 of 2 rejected were right | 456 accepted, 0 of them wrong; 0 of 24 rejected were right |
| F5 expanding shells | 480 accepted, 0 of them wrong; 0 of 0 rejected were right | 480 accepted, 0 of them wrong; 0 of 0 rejected were right | 480 accepted, 0 of them wrong; 0 of 0 rejected were right | 480 accepted, 0 of them wrong; 0 of 0 rejected were right | 480 accepted, 0 of them wrong; 0 of 0 rejected were right | 480 accepted, 0 of them wrong; 0 of 0 rejected were right |
| F5 then F4 | 480 accepted, 0 of them wrong; 0 of 0 rejected were right | 480 accepted, 0 of them wrong; 0 of 0 rejected were right | 480 accepted, 0 of them wrong; 0 of 0 rejected were right | 480 accepted, 0 of them wrong; 0 of 0 rejected were right | 480 accepted, 0 of them wrong; 0 of 0 rejected were right | 480 accepted, 0 of them wrong; 0 of 0 rejected were right |
<!-- /table:evalF_e0_acc -->

<!-- table:evalF_e0_acc_margin -->
| method | runs | accepted but wrong | rejected but right | largest accepted strict cost | smallest rejected strict cost | largest error accepted (deg) | smallest error rejected (deg) |
|---|---|---|---|---|---|---|---|
| F4 Gauss-Newton on E4 rays | all offsets | 0 of 2854 (0.0-0.1%) | 0 of 26 (0.0-12.9%) | 0.602 | 0.951 | 0.168 | 3.14 |
| F4 Gauss-Newton on E4 rays | 5 deg | 0 of 456 (0.0-0.8%) | 0 of 24 (0.0-13.8%) | 0.472 | 0.951 | 0.168 | 3.71 |
| F4 with a 5-frame matching window | all offsets | 0 of 2880 (0.0-0.1%) | none | 0.797 | n/a | 0.245 | n/a |
| F4 with a 5-frame matching window | 5 deg | 0 of 480 (0.0-0.8%) | none | 0.469 | n/a | 0.169 | n/a |
<!-- /table:evalF_e0_acc_margin -->

<!-- table:evalF_e0_cap005 -->
| method | 0.1 deg | 0.5 deg | 1 deg | 2 deg | 3 deg | 5 deg |
|---|---|---|---|---|---|---|
| F1 FindOptimal, box = offset | 86% (83-89) | 66% (62-70) | 58% (53-62) | 32% (28-36) | 21% (18-25) | 7% (5-9) |
| F2 CMA-ES, sigma0 0.5 deg | 100% (99-100) | 99% (98-100) | 99% (98-100) | 81% (77-84) | 46% (42-51) | 11% (8-14) |
| F2 CMA-ES, sigma0 1 deg | 99% (98-100) | 99% (98-100) | 98% (96-99) | 96% (94-98) | 79% (75-82) | 32% (28-36) |
| F2 CMA-ES, sigma0 2 deg | 94% (92-96) | 90% (87-93) | 92% (89-94) | 86% (83-89) | 80% (76-83) | 43% (39-48) |
| F3 centroid Gauss-Newton (windows) | 96% (94-97) | 96% (94-97) | 96% (94-97) | 96% (94-97) | 96% (94-97) | 87% (83-89) |
| F4 Gauss-Newton on E4 rays | 91% (88-93) | 90% (87-93) | 91% (88-93) | 91% (89-94) | 90% (87-92) | 87% (84-90) |
| F4 with a 5-frame matching window | 92% (89-94) | 90% (87-93) | 93% (90-95) | 93% (91-95) | 93% (91-95) | 93% (90-95) |
| F5 expanding shells | 0% (0-1) | 0% (0-1) | 0% (0-1) | 0% (0-1) | 6% (4-8) | 0% (0-1) |
| F5 then F4 | 91% (88-93) | 91% (88-93) | 90% (87-92) | 91% (88-93) | 91% (88-93) | 90% (87-92) |
<!-- /table:evalF_e0_cap005 -->

<!-- table:evalF_e0_cap01 -->
| method | 0.1 deg | 0.5 deg | 1 deg | 2 deg | 3 deg | 5 deg |
|---|---|---|---|---|---|---|
| F1 FindOptimal, box = offset | 95% (93-97) | 79% (75-82) | 72% (67-76) | 47% (42-51) | 33% (29-37) | 12% (9-15) |
| F2 CMA-ES, sigma0 0.5 deg | 100% (99-100) | 100% (99-100) | 99% (98-100) | 81% (77-84) | 47% (42-51) | 11% (8-14) |
| F2 CMA-ES, sigma0 1 deg | 99% (98-100) | 99% (98-100) | 99% (97-99) | 96% (94-98) | 79% (75-82) | 32% (28-36) |
| F2 CMA-ES, sigma0 2 deg | 96% (94-98) | 90% (87-93) | 92% (89-94) | 87% (84-90) | 80% (76-83) | 43% (39-48) |
| F3 centroid Gauss-Newton (windows) | 100% (99-100) | 100% (99-100) | 100% (99-100) | 100% (99-100) | 100% (99-100) | 88% (85-91) |
| F4 Gauss-Newton on E4 rays | 99% (98-100) | 99% (97-99) | 99% (98-100) | 98% (97-99) | 98% (96-99) | 94% (91-96) |
| F4 with a 5-frame matching window | 100% (98-100) | 99% (97-99) | 99% (98-100) | 99% (97-99) | 99% (98-100) | 99% (97-99) |
| F5 expanding shells | 29% (25-33) | 0% (0-1) | 0% (0-1) | 0% (0-1) | 24% (20-28) | 0% (0-1) |
| F5 then F4 | 99% (98-100) | 98% (97-99) | 98% (96-99) | 99% (98-100) | 99% (97-99) | 99% (98-100) |
<!-- /table:evalF_e0_cap01 -->

<!-- table:evalF_e0_cap1 -->
| method | 0.1 deg | 0.5 deg | 1 deg | 2 deg | 3 deg | 5 deg |
|---|---|---|---|---|---|---|
| F1 FindOptimal, box = offset | 100% (99-100) | 100% (99-100) | 99% (98-100) | 85% (81-88) | 66% (62-70) | 42% (38-47) |
| F2 CMA-ES, sigma0 0.5 deg | 100% (99-100) | 100% (99-100) | 100% (99-100) | 82% (78-85) | 47% (43-52) | 11% (9-14) |
| F2 CMA-ES, sigma0 1 deg | 100% (99-100) | 100% (99-100) | 99% (98-100) | 98% (96-99) | 80% (76-84) | 32% (28-37) |
| F2 CMA-ES, sigma0 2 deg | 100% (99-100) | 99% (98-100) | 93% (91-95) | 89% (86-92) | 81% (77-84) | 46% (42-51) |
| F3 centroid Gauss-Newton (windows) | 100% (99-100) | 100% (99-100) | 100% (99-100) | 100% (99-100) | 100% (99-100) | 88% (85-91) |
| F4 Gauss-Newton on E4 rays | 100% (99-100) | 100% (99-100) | 100% (99-100) | 100% (99-100) | 100% (98-100) | 95% (93-97) |
| F4 with a 5-frame matching window | 100% (99-100) | 100% (99-100) | 100% (99-100) | 100% (99-100) | 100% (99-100) | 100% (99-100) |
| F5 expanding shells | 100% (99-100) | 100% (99-100) | 100% (99-100) | 100% (99-100) | 100% (99-100) | 100% (99-100) |
| F5 then F4 | 100% (99-100) | 100% (99-100) | 100% (99-100) | 100% (99-100) | 100% (99-100) | 100% (99-100) |
<!-- /table:evalF_e0_cap1 -->

<!-- table:evalF_e0_choice -->
| method | score 2 deg | score 5 deg | rank sum | mean score (tie-break) | carried |
|---|---|---|---|---|---|
| F4 Gauss-Newton on E4 rays | 3.29 | 2.65 | 3 | 2.97 | yes |
| F5 expanding shells | 4.56 | 0.733 | 3 | 2.65 | yes |
| F2 CMA-ES, sigma0 2 deg | 2.6 | 0.677 | 7 | 1.64 |  |
| F2 CMA-ES, sigma0 1 deg | 3.16 | 0.493 | 8 | 1.82 |  |
| F5 then F4 | 2.58 | 0.652 | 9 | 1.61 |  |
| F2 CMA-ES, sigma0 0.5 deg | 2.53 | 0.154 | 13 | 1.34 |  |
| F1 FindOptimal, box = offset | 1.81 | 0.489 | 13 | 1.15 |  |
<!-- /table:evalF_e0_choice -->

<!-- table:evalF_e0_choice_v1 -->
| method | score 2 deg | score 5 deg | rank sum | mean score (tie-break) | carried |
|---|---|---|---|---|---|
| F4 Gauss-Newton on E4 rays | 3.29 | 2.65 | 2 | 2.97 | yes |
| F2 CMA-ES, sigma0 1 deg | 3.16 | 0.493 | 6 | 1.82 | yes |
| F5 expanding shells | 2.8 | 0.503 | 6 | 1.65 |  |
| F2 CMA-ES, sigma0 2 deg | 2.6 | 0.677 | 6 | 1.64 |  |
| F2 CMA-ES, sigma0 0.5 deg | 2.53 | 0.154 | 12 | 1.34 |  |
| F5 then F4 | 1.9 | 0.464 | 12 | 1.18 |  |
| F1 FindOptimal, box = offset | 1.81 | 0.489 | 12 | 1.15 |  |
<!-- /table:evalF_e0_choice_v1 -->

<!-- table:evalF_e0_choice_wall -->
| method | score 2 deg | score 5 deg | rank sum | mean score (tie-break) | carried |
|---|---|---|---|---|---|
| F4 Gauss-Newton on E4 rays | 3.56 | 2.95 | 3 | 3.25 | yes |
| F5 expanding shells | 4.51 | 0.746 | 3 | 2.63 | yes |
| F5 then F4 | 2.62 | 0.665 | 8 | 1.64 |  |
| F2 CMA-ES, sigma0 2 deg | 2.59 | 0.672 | 8 | 1.63 |  |
| F2 CMA-ES, sigma0 1 deg | 3.13 | 0.492 | 9 | 1.81 |  |
| F1 FindOptimal, box = offset | 1.87 | 0.501 | 12 | 1.19 |  |
| F2 CMA-ES, sigma0 0.5 deg | 2.51 | 0.154 | 13 | 1.33 |  |
<!-- /table:evalF_e0_choice_wall -->

<!-- table:evalF_e0_choice_with_f4w -->
| method | score 2 deg | score 5 deg | rank sum | mean score (tie-break) | carried |
|---|---|---|---|---|---|
| F4 with a 5-frame matching window | 3.29 | 2.79 | 3 | 3.04 |  |
| F4 Gauss-Newton on E4 rays | 3.29 | 2.65 | 4 | 2.97 |  |
| F5 expanding shells | 4.56 | 0.733 | 4 | 2.65 |  |
| F2 CMA-ES, sigma0 2 deg | 2.6 | 0.677 | 9 | 1.64 |  |
| F2 CMA-ES, sigma0 1 deg | 3.16 | 0.493 | 10 | 1.82 |  |
| F5 then F4 | 2.58 | 0.652 | 11 | 1.61 |  |
| F2 CMA-ES, sigma0 0.5 deg | 2.53 | 0.154 | 15 | 1.34 |  |
| F1 FindOptimal, box = offset | 1.81 | 0.489 | 15 | 1.15 |  |
<!-- /table:evalF_e0_choice_with_f4w -->

<!-- table:evalF_e0_err -->
| method | 0.1 deg | 0.5 deg | 1 deg | 2 deg | 3 deg | 5 deg |
|---|---|---|---|---|---|---|
| F1 FindOptimal, box = offset | 0.0129 | 0.0212 | 0.0360 | 0.1118 | 0.2374 | 2.4255 |
| F2 CMA-ES, sigma0 0.5 deg | 0.0000 | 0.0000 | 0.0000 | 0.0076 | 1.9106 | 4.9791 |
| F2 CMA-ES, sigma0 1 deg | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0085 | 4.3222 |
| F2 CMA-ES, sigma0 2 deg | 0.0012 | 0.0014 | 0.0020 | 0.0040 | 0.0084 | 2.3745 |
| F3 centroid Gauss-Newton (windows) | 0.0110 | 0.0110 | 0.0110 | 0.0110 | 0.0110 | 0.0147 |
| F4 Gauss-Newton on E4 rays | 0.0232 | 0.0238 | 0.0241 | 0.0235 | 0.0241 | 0.0244 |
| F4 with a 5-frame matching window | 0.0220 | 0.0228 | 0.0220 | 0.0220 | 0.0228 | 0.0220 |
| F5 expanding shells | 0.0996 | 0.2249 | 0.1608 | 0.2266 | 0.1460 | 0.2238 |
| F5 then F4 | 0.0232 | 0.0244 | 0.0244 | 0.0241 | 0.0232 | 0.0232 |
<!-- /table:evalF_e0_err -->

<!-- table:evalF_e0_nearfar_cap005 -->
| method | near-axis, 2 deg | near-axis, 3 deg | near-axis, 5 deg | far-from-axis, 2 deg | far-from-axis, 3 deg | far-from-axis, 5 deg |
|---|---|---|---|---|---|---|
| F1 FindOptimal, box = offset | 22% (17-30) | 15% (10-21) | 4% (2-8) | 37% (31-42) | 24% (20-29) | 8% (6-12) |
| F2 CMA-ES, sigma0 0.5 deg | 77% (70-83) | 47% (39-55) | 12% (8-18) | 83% (79-87) | 46% (41-52) | 10% (7-14) |
| F2 CMA-ES, sigma0 1 deg | 94% (90-97) | 76% (68-82) | 31% (25-39) | 98% (95-99) | 80% (75-84) | 32% (27-37) |
| F2 CMA-ES, sigma0 2 deg | 84% (78-89) | 82% (76-88) | 51% (43-58) | 88% (83-91) | 78% (73-82) | 39% (34-45) |
| F3 centroid Gauss-Newton (windows) | 88% (81-92) | 88% (81-92) | 85% (79-90) | 100% (99-100) | 100% (99-100) | 88% (83-91) |
| F4 Gauss-Newton on E4 rays | 96% (92-98) | 97% (93-99) | 92% (87-95) | 89% (85-92) | 87% (82-90) | 84% (80-88) |
| F4 with a 5-frame matching window | 97% (93-99) | 99% (96-100) | 96% (92-98) | 92% (88-94) | 91% (87-93) | 91% (88-94) |
| F5 expanding shells | 0% (0-2) | 6% (3-11) | 0% (0-2) | 0% (0-1) | 6% (4-9) | 0% (0-1) |
| F5 then F4 | 94% (89-97) | 93% (88-96) | 92% (87-96) | 89% (85-92) | 90% (86-93) | 89% (85-92) |
<!-- /table:evalF_e0_nearfar_cap005 -->

<!-- table:evalF_e0_nearfar_cap1 -->
| method | near-axis, 2 deg | near-axis, 3 deg | near-axis, 5 deg | far-from-axis, 2 deg | far-from-axis, 3 deg | far-from-axis, 5 deg |
|---|---|---|---|---|---|---|
| F1 FindOptimal, box = offset | 84% (77-89) | 66% (59-73) | 44% (36-51) | 85% (81-88) | 66% (61-71) | 42% (37-47) |
| F2 CMA-ES, sigma0 0.5 deg | 77% (70-83) | 49% (42-57) | 12% (8-18) | 84% (80-88) | 46% (41-52) | 11% (8-14) |
| F2 CMA-ES, sigma0 1 deg | 97% (93-99) | 77% (70-83) | 34% (27-41) | 98% (96-99) | 82% (77-86) | 32% (27-37) |
| F2 CMA-ES, sigma0 2 deg | 89% (83-93) | 84% (78-89) | 53% (45-61) | 90% (86-93) | 79% (75-83) | 43% (38-48) |
| F3 centroid Gauss-Newton (windows) | 100% (98-100) | 100% (98-100) | 89% (84-93) | 100% (99-100) | 100% (99-100) | 88% (83-91) |
| F4 Gauss-Newton on E4 rays | 100% (98-100) | 99% (97-100) | 95% (90-97) | 100% (99-100) | 100% (98-100) | 95% (92-97) |
| F4 with a 5-frame matching window | 100% (98-100) | 100% (98-100) | 100% (98-100) | 100% (99-100) | 100% (99-100) | 100% (99-100) |
| F5 expanding shells | 100% (98-100) | 100% (98-100) | 100% (98-100) | 100% (99-100) | 100% (99-100) | 100% (99-100) |
| F5 then F4 | 100% (98-100) | 100% (98-100) | 100% (98-100) | 100% (99-100) | 100% (99-100) | 100% (99-100) |
<!-- /table:evalF_e0_nearfar_cap1 -->

<!-- table:evalF_e0_pervoxel_cap1 -->
| method | 2 deg: voxels with a lost run | 2 deg: fewest within 1 deg | 3 deg: voxels with a lost run | 3 deg: fewest within 1 deg | 5 deg: voxels with a lost run | 5 deg: fewest within 1 deg |
|---|---|---|---|---|---|---|
| F1 FindOptimal, box = offset | 23 of 24 | 12 of 20 | 24 of 24 | 7 of 20 | 24 of 24 | 6 of 20 |
| F2 CMA-ES, sigma0 1 deg | 10 of 24 | 18 of 20 | 24 of 24 | 13 of 20 | 24 of 24 | 3 of 20 |
| F2 CMA-ES, sigma0 2 deg | 22 of 24 | 16 of 20 | 24 of 24 | 13 of 20 | 24 of 24 | 3 of 20 |
| F4 Gauss-Newton on E4 rays | 0 of 24 | 20 of 20 | 2 of 24 | 19 of 20 | 14 of 24 | 17 of 20 |
| F4 with a 5-frame matching window | 0 of 24 | 20 of 20 | 0 of 24 | 20 of 20 | 0 of 24 | 20 of 20 |
| F5 expanding shells | 0 of 24 | 20 of 20 | 0 of 24 | 20 of 20 | 0 of 24 | 20 of 20 |
| F5 then F4 | 0 of 24 | 20 of 20 | 0 of 24 | 20 of 20 | 0 of 24 | 20 of 20 |
<!-- /table:evalF_e0_pervoxel_cap1 -->

<!-- table:evalF_e0_policy -->
| offset | runs | F4 rejected, F5+F4 run | within 1 deg | within 0.1 deg | within 0.05 deg | median evals (strict-equivalent) | 90th percentile |
|---|---|---|---|---|---|---|---|
| 0.1 deg | 480 | 0 | 100% (99-100) | 99% (98-100) | 91% (88-93) | 141 | 205 |
| 0.5 deg | 480 | 0 | 100% (99-100) | 99% (97-99) | 90% (87-93) | 214 | 288 |
| 1 deg | 480 | 0 | 100% (99-100) | 99% (98-100) | 91% (88-93) | 268 | 340 |
| 2 deg | 480 | 0 | 100% (99-100) | 98% (97-99) | 91% (89-94) | 304 | 377 |
| 3 deg | 480 | 2 | 100% (99-100) | 98% (97-99) | 90% (87-93) | 331 | 410 |
| 5 deg | 480 | 24 | 100% (99-100) | 99% (97-99) | 91% (89-94) | 358 | 468 |
<!-- /table:evalF_e0_policy -->

<!-- table:evalF_e0_policy_f4w -->
| offset | runs | F4 w5 rejected, F5+F4 run | within 1 deg | within 0.1 deg | within 0.05 deg | median evals (strict-equivalent) | 90th percentile |
|---|---|---|---|---|---|---|---|
| 0.1 deg | 480 | 0 | 100% (99-100) | 100% (98-100) | 92% (89-94) | 141 | 216 |
| 0.5 deg | 480 | 0 | 100% (99-100) | 99% (97-99) | 90% (87-93) | 214 | 288 |
| 1 deg | 480 | 0 | 100% (99-100) | 99% (98-100) | 93% (90-95) | 268 | 359 |
| 2 deg | 480 | 0 | 100% (99-100) | 99% (97-99) | 93% (91-95) | 304 | 377 |
| 3 deg | 480 | 0 | 100% (99-100) | 99% (98-100) | 93% (91-95) | 331 | 395 |
| 5 deg | 480 | 0 | 100% (99-100) | 99% (97-99) | 93% (90-95) | 358 | 425 |
<!-- /table:evalF_e0_policy_f4w -->

<!-- table:evalF_e0_se_vs_wall -->
| method | 1 deg | 2 deg | 5 deg |
|---|---|---|---|
| F1 FindOptimal, box = offset | 1.04 -> 1.04 | 1.03 -> 1.03 | 1.03 -> 1.03 |
| F2 CMA-ES, sigma0 0.5 deg | 0.99 -> 0.99 | 0.99 -> 0.99 | 0.99 -> 0.99 |
| F2 CMA-ES, sigma0 1 deg | 0.99 -> 0.99 | 0.99 -> 0.99 | 1.00 -> 1.00 |
| F2 CMA-ES, sigma0 2 deg | 1.00 -> 1.00 | 1.00 -> 1.00 | 0.99 -> 0.99 |
| F4 Gauss-Newton on E4 rays | 1.06 -> 1.06 | 1.08 -> 1.08 | 1.11 -> 1.11 |
| F4 with a 5-frame matching window | 1.05 -> 1.05 | 1.05 -> 1.05 | 1.05 -> 1.05 |
| F5 expanding shells | 1.62 -> 0.96 | 1.61 -> 0.99 | 1.48 -> 1.02 |
| F5 then F4 | 1.25 -> 1.01 | 1.41 -> 1.02 | 1.43 -> 1.02 |
<!-- /table:evalF_e0_se_vs_wall -->

<!-- table:evalF_e0_secs -->
| method | 0.1 deg | 0.5 deg | 1 deg | 2 deg | 3 deg | 5 deg |
|---|---|---|---|---|---|---|
| F1 FindOptimal, box = offset | 0.250 | 0.216 | 0.198 | 0.204 | 0.229 | 0.370 |
| F2 CMA-ES, sigma0 0.5 deg | 0.105 | 0.110 | 0.119 | 0.150 | 0.204 | 0.325 |
| F2 CMA-ES, sigma0 1 deg | 0.123 | 0.123 | 0.127 | 0.141 | 0.166 | 0.296 |
| F2 CMA-ES, sigma0 2 deg | 0.145 | 0.151 | 0.151 | 0.158 | 0.173 | 0.308 |
| F3 centroid Gauss-Newton (windows) | 0.176 | 0.180 | 0.186 | 0.189 | 0.191 | 0.183 |
| F4 Gauss-Newton on E4 rays | 0.057 | 0.087 | 0.107 | 0.122 | 0.131 | 0.138 |
| F4 with a 5-frame matching window | 0.060 | 0.090 | 0.113 | 0.128 | 0.136 | 0.150 |
| F5 expanding shells | 0.000 | 0.008 | 0.037 | 0.098 | 0.218 | 0.595 |
| F5 then F4 | 0.057 | 0.079 | 0.099 | 0.168 | 0.276 | 0.666 |
<!-- /table:evalF_e0_secs -->

<!-- table:evalF_e0_secs_single -->
| method | 0.1 deg | 0.5 deg | 1 deg | 2 deg | 3 deg | 5 deg |
|---|---|---|---|---|---|---|
| F1 FindOptimal, box = offset | n/a | n/a | 0.183 | 0.190 | n/a | 0.355 |
| F2 CMA-ES, sigma0 0.5 deg | n/a | n/a | 0.142 | 0.146 | n/a | 0.333 |
| F2 CMA-ES, sigma0 1 deg | n/a | n/a | 0.125 | 0.127 | n/a | 0.323 |
| F2 CMA-ES, sigma0 2 deg | n/a | n/a | 0.151 | 0.143 | n/a | 0.453 |
| F3 centroid Gauss-Newton (windows) | n/a | n/a | 0.177 | 0.182 | n/a | 0.185 |
| F4 Gauss-Newton on E4 rays | n/a | n/a | 0.091 | 0.112 | n/a | 0.124 |
| F4 with a 5-frame matching window | n/a | n/a | 0.091 | 0.113 | n/a | 0.133 |
| F5 expanding shells | n/a | n/a | 0.037 | 0.094 | n/a | 0.579 |
| F5 then F4 | n/a | n/a | 0.095 | 0.157 | n/a | 0.656 |
<!-- /table:evalF_e0_secs_single -->

<!-- table:evalF_e0_stop_raw -->
| method | 0.1 deg | 0.5 deg | 1 deg | 2 deg | 3 deg | 5 deg |
|---|---|---|---|---|---|---|
| F1 FindOptimal, box = offset | 573 (650) | 496 (540) | 452 (474) | 467 (505) | 529 (551) | 870 (870) |
| F2 CMA-ES, sigma0 0.5 deg | 225 (281) | 239 (302) | 260 (323) | 323 (562) | 456 (1030) | 715 (1422) |
| F2 CMA-ES, sigma0 1 deg | 267 (330) | 267 (330) | 274 (345) | 309 (407) | 365 (786) | 659 (1919) |
| F2 CMA-ES, sigma0 2 deg | 316 (457) | 330 (582) | 330 (541) | 344 (758) | 379 (1054) | 684 (2544) |
| F3 centroid Gauss-Newton (windows) | n/a (n/a) | n/a (n/a) | n/a (n/a) | n/a (n/a) | n/a (n/a) | n/a (n/a) |
| F4 Gauss-Newton on E4 rays | 14 (35) | 20 (40) | 26 (45) | 28 (43) | 30 (47) | 33 (48) |
| F4 with a 5-frame matching window | 14 (34) | 20 (40) | 26 (46) | 28 (45) | 30 (46) | 32 (50) |
| F5 expanding shells | 1 (1) | 32 (125) | 290 (290) | 1511 (2257) | 5916 (5916) | 22244 (26472) |
| F5 then F4 | 15 (36) | 52.5 (141) | 304 (326) | 1529 (2271) | 5928 (5943) | 22261 (26486) |
<!-- /table:evalF_e0_stop_raw -->

<!-- table:evalF_e0_stop_se -->
| method | 0.1 deg | 0.5 deg | 1 deg | 2 deg | 3 deg | 5 deg |
|---|---|---|---|---|---|---|
| F1 FindOptimal, box = offset | 573 (650) | 496 (540) | 452 (474) | 467 (505) | 529 (551) | 870 (870) |
| F2 CMA-ES, sigma0 0.5 deg | 225 (281) | 239 (302) | 260 (323) | 323 (562) | 456 (1030) | 715 (1422) |
| F2 CMA-ES, sigma0 1 deg | 267 (330) | 267 (330) | 274 (345) | 309 (407) | 365 (786) | 659 (1919) |
| F2 CMA-ES, sigma0 2 deg | 316 (457) | 330 (582) | 330 (541) | 344 (758) | 379 (1054) | 684 (2544) |
| F3 centroid Gauss-Newton (windows) | 426 (426) | 426 (426) | 426 (426) | 426 (426) | 426 (426) | 426 (426) |
| F4 Gauss-Newton on E4 rays | 141 (205) | 214 (288) | 268 (340) | 304 (377) | 331 (404) | 358 (431) |
| F4 with a 5-frame matching window | 141 (216) | 214 (288) | 268 (359) | 304 (377) | 331 (395) | 358 (425) |
| F5 expanding shells | 1 (1) | 15.7 (48) | 83.2 (83.2) | 219 (279) | 516 (516) | 1364 (1568) |
| F5 then F4 | 142 (206) | 184 (255) | 224 (314) | 388 (452) | 647 (713) | 1533 (1709) |
<!-- /table:evalF_e0_stop_se -->

<!-- table:evalF_e0_timing -->
| kind of evaluation | used for | milliseconds | times a strict evaluation |
|---|---|---|---|
| strict cost | all methods | 0.4219 | 1.000 |
| screen (pixel radius 3), per orientation, batches of 64 | v1 weights only | 0.0230 | 0.054 |
| E4 value, one rotation per call | F4 | 0.6584 | 1.561 |
| E4 value and Jacobian, one rotation per call | F4 | 10.7970 | 25.592 |
| screen, per orientation, a shell at a time (chunks of up to 256) | F5 | 0.0175 | 0.042 |
| E4 value, per orientation, 32 rotations per call | F5 | 0.3670 | 0.870 |
<!-- /table:evalF_e0_timing -->

<!-- table:evalF_e0_timing_voxels -->
| voxel | weights from it | F4 E4 value (ms) | F4 E4 value and Jacobian (ms) | F5 screen (ms per orientation) | F5 E4 (ms per orientation) | strict inside F5 (ms) |
|---|---|---|---|---|---|---|
| 10531 | yes | 0.649 | 11.04 | 0.0175 | 0.367 | 0.443 |
| 3108 |  | 0.630 | 8.76 | 0.0168 | 0.352 | 0.464 |
<!-- /table:evalF_e0_timing_voxels -->

<!-- table:evalF_e0_tocap005 -->
| method | 0.1 deg | 0.5 deg | 1 deg | 2 deg | 3 deg | 5 deg |
|---|---|---|---|---|---|---|
| F1 FindOptimal, box = offset | 32 (84.4) n=434 | 61 (138) n=327 | 87 (178) n=287 | 136 (204) n=167 | 140 (221) n=109 | 164 (243) n=33 |
| F2 CMA-ES, sigma0 0.5 deg | 83.5 (130) n=480 | 94 (145) n=479 | 115 (168) n=477 | 169 (253) n=389 | 215 (342) n=223 | 306 (654) n=52 |
| F2 CMA-ES, sigma0 1 deg | 122 (178) n=478 | 120 (179) n=475 | 131 (191) n=472 | 157 (234) n=463 | 200 (298) n=378 | 274 (429) n=153 |
| F2 CMA-ES, sigma0 2 deg | 167 (255) n=453 | 178 (275) n=434 | 176 (270) n=442 | 184 (306) n=416 | 208 (328) n=383 | 248 (464) n=208 |
| F3 centroid Gauss-Newton (windows) | 3 (3) n=460 | 3 (3) n=460 | 3 (3) n=460 | 3 (3) n=460 | 3 (3) n=460 | 3 (3) n=416 |
| F4 Gauss-Newton on E4 rays | 5 (9) n=441 | 13 (15) n=437 | 17 (19.6) n=445 | 19 (23) n=447 | 21 (25) n=445 | 23 (28.7) n=424 |
| F4 with a 5-frame matching window | 5 (9) n=442 | 13 (15) n=435 | 15 (19.5) n=446 | 19 (23) n=448 | 21 (25) n=449 | 23 (27) n=446 |
| F5 expanding shells | n/a | n/a | n/a | n/a | 5916 (5916) n=28 | n/a |
| F5 then F4 | 6 (10) n=441 | 43 (134) n=440 | 297 (301) n=434 | 1521 (2264) n=439 | 5921 (5925) n=438 | 22253 (26479) n=439 |
<!-- /table:evalF_e0_tocap005 -->

<!-- table:evalF_e0_tocap1 -->
| method | 0.1 deg | 0.5 deg | 1 deg | 2 deg | 3 deg | 5 deg |
|---|---|---|---|---|---|---|
| F1 FindOptimal, box = offset | 0 (0) n=480 | 0 (0) n=480 | 5 (30.8) n=479 | 32 (193) n=409 | 81 (223) n=325 | 126 (433) n=211 |
| F2 CMA-ES, sigma0 0.5 deg | 0 (0) n=480 | 0 (0) n=480 | 4 (13) n=479 | 32 (110) n=395 | 68 (177) n=227 | 133 (416) n=54 |
| F2 CMA-ES, sigma0 1 deg | 0 (0) n=480 | 0 (0) n=480 | 9 (27.4) n=477 | 23 (74) n=471 | 44.5 (126) n=386 | 109 (236) n=157 |
| F2 CMA-ES, sigma0 2 deg | 0 (0) n=480 | 0 (0) n=480 | 25 (82) n=451 | 30.5 (101) n=434 | 48 (135) n=392 | 82.5 (268) n=226 |
| F3 centroid Gauss-Newton (windows) | 0 (0) n=480 | 0 (0) n=480 | 3 (3) n=480 | 3 (3) n=480 | 3 (3) n=480 | 3 (3) n=423 |
| F4 Gauss-Newton on E4 rays | 0 (0) n=480 | 0 (0) n=480 | 5 (5) n=480 | 7 (7) n=480 | 9 (9) n=478 | 11 (11) n=456 |
| F4 with a 5-frame matching window | 0 (0) n=480 | 0 (0) n=480 | 5 (5) n=480 | 7 (7) n=480 | 9 (9) n=480 | 11 (11) n=480 |
| F5 expanding shells | 0 (0) n=480 | 0 (0) n=480 | 33 (33) n=480 | 558 (955) n=480 | 2258 (3221) n=480 | 15197 (18496) n=480 |
| F5 then F4 | 0 (0) n=480 | 0 (0) n=480 | 33 (33) n=480 | 558 (955) n=480 | 2258 (3221) n=480 | 15197 (18496) n=480 |
<!-- /table:evalF_e0_tocap1 -->

<!-- table:evalF_full_acc -->
| method | 1 deg | 2 deg | 5 deg |
|---|---|---|---|
| F4 Gauss-Newton on E4 rays | 120 accepted, 0 of them wrong; 0 of 0 rejected were right | 120 accepted, 0 of them wrong; 0 of 0 rejected were right | 81 accepted, 0 of them wrong; 0 of 39 rejected were right |
| F5 expanding shells | 120 accepted, 0 of them wrong; 0 of 0 rejected were right | 120 accepted, 0 of them wrong; 0 of 0 rejected were right | 120 accepted, 0 of them wrong; 0 of 0 rejected were right |
| F5 then F4 | 120 accepted, 0 of them wrong; 0 of 0 rejected were right | 120 accepted, 0 of them wrong; 0 of 0 rejected were right | 120 accepted, 0 of them wrong; 0 of 0 rejected were right |
<!-- /table:evalF_full_acc -->

<!-- table:evalF_full_acc_margin -->
| method | runs | accepted but wrong | rejected but right | largest accepted strict cost | smallest rejected strict cost | largest error accepted (deg) | smallest error rejected (deg) |
|---|---|---|---|---|---|---|---|
| F4 Gauss-Newton on E4 rays | all offsets | 0 of 321 (0.0-1.2%) | 0 of 39 (0.0-9.0%) | 0.437 | 0.953 | 0.295 | 3.37 |
| F4 Gauss-Newton on E4 rays | 5 deg | 0 of 81 (0.0-4.5%) | 0 of 39 (0.0-9.0%) | 0.159 | 0.953 | 0.086 | 3.37 |
| F4 with a 5-frame matching window | all offsets | 0 of 329 (0.0-1.2%) | 0 of 31 (0.0-11.0%) | 0.437 | 0.953 | 0.295 | 3.42 |
| F4 with a 5-frame matching window | 5 deg | 0 of 89 (0.0-4.1%) | 0 of 31 (0.0-11.0%) | 0.264 | 0.953 | 0.244 | 3.42 |
<!-- /table:evalF_full_acc_margin -->

<!-- table:evalF_full_cap005 -->
| method | 1 deg | 2 deg | 5 deg |
|---|---|---|---|
| F1 FindOptimal, box = offset | 81% (73-87) | 63% (54-71) | 22% (15-30) |
| F2 CMA-ES, sigma0 0.5 deg | 89% (82-94) | 82% (74-88) | 8% (5-15) |
| F2 CMA-ES, sigma0 1 deg | 92% (85-95) | 89% (82-94) | 41% (32-50) |
| F2 CMA-ES, sigma0 2 deg | 91% (84-95) | 92% (86-96) | 66% (57-74) |
| F4 Gauss-Newton on E4 rays | 33% (26-42) | 33% (26-42) | 25% (18-33) |
| F4 with a 5-frame matching window | 50% (41-59) | 50% (41-59) | 37% (29-46) |
| F5 expanding shells | 0% (0-3) | 0% (0-3) | 0% (0-3) |
| F5 then F4 | 33% (26-42) | 33% (26-42) | 33% (26-42) |
<!-- /table:evalF_full_cap005 -->

<!-- table:evalF_full_cap01 -->
| method | 1 deg | 2 deg | 5 deg |
|---|---|---|---|
| F1 FindOptimal, box = offset | 100% (97-100) | 87% (79-92) | 32% (25-41) |
| F2 CMA-ES, sigma0 0.5 deg | 100% (97-100) | 92% (86-96) | 10% (6-17) |
| F2 CMA-ES, sigma0 1 deg | 100% (97-100) | 100% (97-100) | 42% (33-51) |
| F2 CMA-ES, sigma0 2 deg | 97% (92-99) | 97% (92-99) | 73% (65-80) |
| F4 Gauss-Newton on E4 rays | 98% (94-100) | 99% (95-100) | 68% (59-75) |
| F4 with a 5-frame matching window | 98% (94-100) | 100% (97-100) | 72% (63-79) |
| F5 expanding shells | 0% (0-3) | 0% (0-3) | 0% (0-3) |
| F5 then F4 | 98% (94-100) | 98% (94-100) | 99% (95-100) |
<!-- /table:evalF_full_cap01 -->

<!-- table:evalF_full_cap1 -->
| method | 1 deg | 2 deg | 5 deg |
|---|---|---|---|
| F1 FindOptimal, box = offset | 100% (97-100) | 96% (91-98) | 56% (47-64) |
| F2 CMA-ES, sigma0 0.5 deg | 100% (97-100) | 92% (86-96) | 11% (6-18) |
| F2 CMA-ES, sigma0 1 deg | 100% (97-100) | 100% (97-100) | 42% (33-51) |
| F2 CMA-ES, sigma0 2 deg | 97% (92-99) | 98% (93-99) | 76% (67-83) |
| F4 Gauss-Newton on E4 rays | 100% (97-100) | 100% (97-100) | 68% (59-75) |
| F4 with a 5-frame matching window | 100% (97-100) | 100% (97-100) | 74% (66-81) |
| F5 expanding shells | 100% (97-100) | 100% (97-100) | 100% (97-100) |
| F5 then F4 | 100% (97-100) | 100% (97-100) | 100% (97-100) |
<!-- /table:evalF_full_cap1 -->

<!-- table:evalF_full_choice -->
| method | score 2 deg | score 5 deg | rank sum | mean score (tie-break) | carried |
|---|---|---|---|---|---|
| F2 CMA-ES, sigma0 2 deg | 2.99 | 1.82 | 3 | 2.4 | yes |
| F2 CMA-ES, sigma0 1 deg | 3.2 | 0.95 | 4 | 2.08 | yes |
| F4 Gauss-Newton on E4 rays | 2.69 | 1.52 | 6 | 2.1 |  |
| F2 CMA-ES, sigma0 0.5 deg | 2.86 | 0.276 | 10 | 1.57 |  |
| F5 expanding shells | 2.29 | 0.494 | 10 | 1.39 |  |
| F1 FindOptimal, box = offset | 2.05 | 0.642 | 10 | 1.35 |  |
| F5 then F4 | 1.41 | 0.436 | 13 | 0.923 |  |
<!-- /table:evalF_full_choice -->

<!-- table:evalF_full_choice_v1 -->
| method | score 2 deg | score 5 deg | rank sum | mean score (tie-break) | carried |
|---|---|---|---|---|---|
| F2 CMA-ES, sigma0 2 deg | 2.99 | 1.82 | 3 | 2.4 | yes |
| F2 CMA-ES, sigma0 1 deg | 3.2 | 0.95 | 4 | 2.08 | yes |
| F4 Gauss-Newton on E4 rays | 2.69 | 1.52 | 6 | 2.1 |  |
| F1 FindOptimal, box = offset | 2.05 | 0.642 | 9 | 1.35 |  |
| F2 CMA-ES, sigma0 0.5 deg | 2.86 | 0.276 | 10 | 1.57 |  |
| F5 expanding shells | 1.86 | 0.396 | 11 | 1.13 |  |
| F5 then F4 | 1.24 | 0.359 | 13 | 0.797 |  |
<!-- /table:evalF_full_choice_v1 -->

<!-- table:evalF_full_choice_wall -->
| method | score 2 deg | score 5 deg | rank sum | mean score (tie-break) | carried |
|---|---|---|---|---|---|
| F2 CMA-ES, sigma0 2 deg | 3.04 | 1.83 | 3 | 2.44 | yes |
| F2 CMA-ES, sigma0 1 deg | 3.23 | 0.966 | 4 | 2.1 | yes |
| F4 Gauss-Newton on E4 rays | 2.86 | 1.6 | 6 | 2.23 |  |
| F2 CMA-ES, sigma0 0.5 deg | 2.87 | 0.28 | 10 | 1.58 |  |
| F1 FindOptimal, box = offset | 2.14 | 0.669 | 10 | 1.4 |  |
| F5 expanding shells | 2.21 | 0.484 | 10 | 1.35 |  |
| F5 then F4 | 1.4 | 0.431 | 13 | 0.917 |  |
<!-- /table:evalF_full_choice_wall -->

<!-- table:evalF_full_choice_with_f4w -->
| method | score 2 deg | score 5 deg | rank sum | mean score (tie-break) | carried |
|---|---|---|---|---|---|
| F2 CMA-ES, sigma0 2 deg | 2.99 | 1.82 | 3 | 2.4 |  |
| F2 CMA-ES, sigma0 1 deg | 3.2 | 0.95 | 5 | 2.08 |  |
| F4 with a 5-frame matching window | 2.72 | 1.61 | 6 | 2.17 |  |
| F4 Gauss-Newton on E4 rays | 2.69 | 1.52 | 8 | 2.1 |  |
| F2 CMA-ES, sigma0 0.5 deg | 2.86 | 0.276 | 11 | 1.57 |  |
| F5 expanding shells | 2.29 | 0.494 | 12 | 1.39 |  |
| F1 FindOptimal, box = offset | 2.05 | 0.642 | 12 | 1.35 |  |
| F5 then F4 | 1.41 | 0.436 | 15 | 0.923 |  |
<!-- /table:evalF_full_choice_with_f4w -->

<!-- table:evalF_full_err -->
| method | 1 deg | 2 deg | 5 deg |
|---|---|---|---|
| F1 FindOptimal, box = offset | 0.0297 | 0.0371 | 0.4314 |
| F2 CMA-ES, sigma0 0.5 deg | 0.0274 | 0.0275 | 4.8527 |
| F2 CMA-ES, sigma0 1 deg | 0.0242 | 0.0276 | 3.8824 |
| F2 CMA-ES, sigma0 2 deg | 0.0255 | 0.0278 | 0.0353 |
| F4 Gauss-Newton on E4 rays | 0.0520 | 0.0520 | 0.0538 |
| F4 with a 5-frame matching window | 0.0486 | 0.0486 | 0.0519 |
| F5 expanding shells | 0.4391 | 0.5281 | 0.5176 |
| F5 then F4 | 0.0520 | 0.0514 | 0.0519 |
<!-- /table:evalF_full_err -->

<!-- table:evalF_full_pervoxel_cap1_5 -->
| method | voxel 5242 | voxel 15901 | voxel 17034 | voxel 18558 | voxel 19889 | voxel 22867 |
|---|---|---|---|---|---|---|
| F1 FindOptimal, box = offset | 11 of 20 | 12 of 20 | 10 of 20 | 15 of 20 | 7 of 20 | 12 of 20 |
| F2 CMA-ES, sigma0 1 deg | 9 of 20 | 9 of 20 | 6 of 20 | 7 of 20 | 7 of 20 | 12 of 20 |
| F2 CMA-ES, sigma0 2 deg | 15 of 20 | 14 of 20 | 18 of 20 | 15 of 20 | 15 of 20 | 14 of 20 |
| F4 Gauss-Newton on E4 rays | 14 of 20 | 18 of 20 | 14 of 20 | 12 of 20 | 13 of 20 | 10 of 20 |
| F4 with a 5-frame matching window | 14 of 20 | 19 of 20 | 17 of 20 | 14 of 20 | 14 of 20 | 11 of 20 |
| F5 expanding shells | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 |
| F5 then F4 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 | 20 of 20 |
<!-- /table:evalF_full_pervoxel_cap1_5 -->

<!-- table:evalF_full_policy -->
| offset | runs | F4 rejected, F5+F4 run | within 1 deg | within 0.1 deg | within 0.05 deg | median evals (strict-equivalent) | 90th percentile |
|---|---|---|---|---|---|---|---|
| 1 deg | 120 | 0 | 100% (97-100) | 98% (94-100) | 33% (26-42) | 317 | 431 |
| 2 deg | 120 | 0 | 100% (97-100) | 99% (95-100) | 33% (26-42) | 372 | 499 |
| 5 deg | 120 | 39 | 100% (97-100) | 99% (95-100) | 33% (26-42) | 547 | 2806 |
<!-- /table:evalF_full_policy -->

<!-- table:evalF_full_policy_f4w -->
| offset | runs | F4 w5 rejected, F5+F4 run | within 1 deg | within 0.1 deg | within 0.05 deg | median evals (strict-equivalent) | 90th percentile |
|---|---|---|---|---|---|---|---|
| 1 deg | 120 | 0 | 100% (97-100) | 98% (94-100) | 50% (41-59) | 319 | 403 |
| 2 deg | 120 | 0 | 100% (97-100) | 100% (97-100) | 50% (41-59) | 367 | 465 |
| 5 deg | 120 | 31 | 100% (97-100) | 98% (93-99) | 42% (34-51) | 502 | 2816 |
<!-- /table:evalF_full_policy_f4w -->

<!-- table:evalF_full_se_vs_wall -->
| method | 1 deg | 2 deg | 5 deg |
|---|---|---|---|
| F1 FindOptimal, box = offset | 1.05 -> 1.05 | 1.04 -> 1.04 | 1.04 -> 1.04 |
| F2 CMA-ES, sigma0 0.5 deg | 1.02 -> 1.02 | 1.00 -> 1.00 | 1.01 -> 1.01 |
| F2 CMA-ES, sigma0 1 deg | 1.01 -> 1.01 | 1.01 -> 1.01 | 1.02 -> 1.02 |
| F2 CMA-ES, sigma0 2 deg | 1.02 -> 1.02 | 1.02 -> 1.02 | 1.01 -> 1.01 |
| F4 Gauss-Newton on E4 rays | 1.07 -> 1.07 | 1.06 -> 1.06 | 1.05 -> 1.05 |
| F4 with a 5-frame matching window | 1.06 -> 1.06 | 1.06 -> 1.06 | 1.05 -> 1.05 |
| F5 expanding shells | 1.18 -> 0.96 | 1.19 -> 0.97 | 1.22 -> 0.98 |
| F5 then F4 | 1.10 -> 1.02 | 1.15 -> 0.99 | 1.20 -> 0.99 |
<!-- /table:evalF_full_se_vs_wall -->

<!-- table:evalF_full_secs -->
| method | 1 deg | 2 deg | 5 deg |
|---|---|---|---|
| F1 FindOptimal, box = offset | 0.223 | 0.219 | 0.390 |
| F2 CMA-ES, sigma0 0.5 deg | 0.134 | 0.153 | 0.187 |
| F2 CMA-ES, sigma0 1 deg | 0.132 | 0.146 | 0.202 |
| F2 CMA-ES, sigma0 2 deg | 0.157 | 0.154 | 0.201 |
| F4 Gauss-Newton on E4 rays | 0.142 | 0.160 | 0.200 |
| F4 with a 5-frame matching window | 0.147 | 0.160 | 0.200 |
| F5 expanding shells | 0.065 | 0.206 | 0.949 |
| F5 then F4 | 0.181 | 0.325 | 1.072 |
<!-- /table:evalF_full_secs -->

<!-- table:evalF_full_secs_single -->
| method | 1 deg | 2 deg | 5 deg |
|---|---|---|---|
| F1 FindOptimal, box = offset | 0.202 | 0.215 | 0.367 |
| F2 CMA-ES, sigma0 0.5 deg | 0.138 | 0.145 | 0.146 |
| F2 CMA-ES, sigma0 1 deg | 0.132 | 0.152 | 0.207 |
| F2 CMA-ES, sigma0 2 deg | 0.155 | 0.154 | 0.221 |
| F4 Gauss-Newton on E4 rays | 0.143 | 0.150 | 0.200 |
| F4 with a 5-frame matching window | 0.132 | 0.144 | 0.185 |
| F5 expanding shells | 0.062 | 0.199 | 0.909 |
| F5 then F4 | 0.160 | 0.293 | 1.013 |
<!-- /table:evalF_full_secs_single -->

<!-- table:evalF_full_stop_raw -->
| method | 1 deg | 2 deg | 5 deg |
|---|---|---|---|
| F1 FindOptimal, box = offset | 496 (540) | 467 (553) | 870 (902) |
| F2 CMA-ES, sigma0 0.5 deg | 281 (422) | 323 (527) | 393 (569) |
| F2 CMA-ES, sigma0 1 deg | 288 (415) | 312 (478) | 438 (625) |
| F2 CMA-ES, sigma0 2 deg | 334 (499) | 326 (464) | 418 (682) |
| F4 Gauss-Newton on E4 rays | 29 (39.1) | 31.5 (43) | 38 (57) |
| F4 with a 5-frame matching window | 29 (36) | 31 (39.1) | 38 (53.1) |
| F5 expanding shells | 125 (125) | 954 (955) | 18495 (18496) |
| F5 then F4 | 145 (156) | 979 (1001) | 18520 (18539) |
<!-- /table:evalF_full_stop_raw -->

<!-- table:evalF_full_stop_se -->
| method | 1 deg | 2 deg | 5 deg |
|---|---|---|---|
| F1 FindOptimal, box = offset | 496 (540) | 467 (553) | 870 (902) |
| F2 CMA-ES, sigma0 0.5 deg | 281 (422) | 323 (527) | 393 (569) |
| F2 CMA-ES, sigma0 1 deg | 288 (415) | 312 (478) | 438 (625) |
| F2 CMA-ES, sigma0 2 deg | 334 (499) | 326 (464) | 418 (682) |
| F4 Gauss-Newton on E4 rays | 317 (431) | 372 (499) | 444 (658) |
| F4 with a 5-frame matching window | 319 (403) | 367 (465) | 461 (657) |
| F5 expanding shells | 137 (137) | 437 (438) | 2026 (2027) |
| F5 then F4 | 392 (474) | 709 (813) | 2292 (2436) |
<!-- /table:evalF_full_stop_se -->

<!-- table:evalF_full_timing -->
| kind of evaluation | used for | milliseconds | times a strict evaluation |
|---|---|---|---|
| strict cost | all methods | 0.4398 | 1.000 |
| screen (pixel radius 3), per orientation, batches of 64 | v1 weights only | 0.0224 | 0.051 |
| E4 value, one rotation per call | F4 | 1.4978 | 3.406 |
| E4 value and Jacobian, one rotation per call | F4 | 12.0732 | 27.452 |
| screen, per orientation, a shell at a time (chunks of up to 256) | F5 | 0.0175 | 0.040 |
| E4 value, per orientation, 32 rotations per call | F5 | 1.2144 | 2.761 |
<!-- /table:evalF_full_timing -->

<!-- table:evalF_full_timing_voxels -->
| voxel | weights from it | F4 E4 value (ms) | F4 E4 value and Jacobian (ms) | F5 screen (ms per orientation) | F5 E4 (ms per orientation) | strict inside F5 (ms) |
|---|---|---|---|---|---|---|
| 15901 | yes | 1.456 | 11.63 | 0.0175 | 1.214 | 0.506 |
<!-- /table:evalF_full_timing_voxels -->

<!-- table:evalF_full_tocap005 -->
| method | 1 deg | 2 deg | 5 deg |
|---|---|---|---|
| F1 FindOptimal, box = offset | 75 (102) n=107 | 112 (184) n=80 | 168 (251) n=28 |
| F2 CMA-ES, sigma0 0.5 deg | 105 (151) n=114 | 146 (194) n=103 | 265 (352) n=11 |
| F2 CMA-ES, sigma0 1 deg | 118 (156) n=114 | 147 (190) n=113 | 239 (446) n=49 |
| F2 CMA-ES, sigma0 2 deg | 153 (208) n=111 | 158 (220) n=111 | 220 (296) n=82 |
| F4 Gauss-Newton on E4 rays | 15 (21) n=73 | 19 (25) n=77 | 25 (36.4) n=54 |
| F4 with a 5-frame matching window | 15 (25.4) n=84 | 19 (28) n=90 | 27 (37) n=66 |
| F5 expanding shells | n/a | n/a | n/a |
| F5 then F4 | 136 (141) n=70 | 967 (1518) n=70 | 18508 (18516) n=75 |
<!-- /table:evalF_full_tocap005 -->

<!-- table:evalF_full_tocap1 -->
| method | 1 deg | 2 deg | 5 deg |
|---|---|---|---|
| F1 FindOptimal, box = offset | 5 (15.1) n=120 | 23 (166) n=115 | 107 (244) n=67 |
| F2 CMA-ES, sigma0 0.5 deg | 4 (11) n=120 | 26 (73) n=111 | 111 (206) n=13 |
| F2 CMA-ES, sigma0 1 deg | 10 (29) n=120 | 18 (54) n=120 | 95 (256) n=50 |
| F2 CMA-ES, sigma0 2 deg | 26 (64.5) n=116 | 27 (65.8) n=117 | 81 (176) n=91 |
| F4 Gauss-Newton on E4 rays | 5 (5) n=120 | 7 (9) n=120 | 13 (21) n=81 |
| F4 with a 5-frame matching window | 5 (5) n=120 | 7 (7) n=120 | 13 (21) n=89 |
| F5 expanding shells | 33 (33) n=120 | 558 (558) n=120 | 15197 (15197) n=120 |
| F5 then F4 | 33 (33) n=120 | 558 (558) n=120 | 15197 (15197) n=120 |
<!-- /table:evalF_full_tocap1 -->
