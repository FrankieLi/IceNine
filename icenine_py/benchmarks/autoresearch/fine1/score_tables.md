<!-- table:fine1_projection_vs_measured -->
| arm | low | high | measured |
|---|---|---|---|
| CS | 22.7 | 34.9 | 25.3 |
| CF | 18.2 | 27.6 | 20.0 |
| PS | 33.7 | 51.0 | 50.9 |
| J5S | 0.6 | 0.8 | 2.1 |
<!-- /table:fine1_projection_vs_measured -->

<!-- table:fine1_ratios -->
| quantity | value |
|---|---|
| CS wall / CSP wall (oversubscribed) | 8.88 |
| CF wall / CFP wall (oversubscribed) | 7.05 |
| PS wall / J5S wall | 24.23 |
| PS evaluations / J5S evaluations | 86.95 |
| CS evaluations / CSP evaluations | 1.21 |
<!-- /table:fine1_ratios -->

<!-- table:fine1_whole_all -->
| arm | wrong | no_output | wrong_with_output | share | median |
|---|---|---|---|---|---|
| CS | 15,604 | 10,845 | 4,759 | 15.877% (15.650, 16.107) | 0.0852 |
| CF | 12,777 | 9,926 | 2,851 | 13.001% (12.792, 13.212) | 0.0753 |
| CSP | 14,914 | 10,411 | 4,503 | 15.175% (14.952, 15.401) | 0.0687 |
| CFP | 13,920 | 10,590 | 3,330 | 14.164% (13.947, 14.383) | 0.0708 |
| J5S | 561 | 0 | 561 | 0.571% (0.526, 0.620) | 0.0295 |
| PS | 714 | 0 | 714 | 0.726% (0.675, 0.782) | 0.0617 |
<!-- /table:fine1_whole_all -->

<!-- table:fine1_whole_common -->
| arm | workers | common_wrong | share | evals | wall_s | wall_h |
|---|---|---|---|---|---|---|
| CS | 1 worker | 4,681 | 5.464% (5.313, 5.618) | 1,670,595,806 | 91,206 | 25.3 |
| CF | 1 worker | 2,734 | 3.191% (3.075, 3.311) | 1,362,495,605 | 71,996 | 20.0 |
| CSP | 12 workers, 13 ranks on 12 cores | 4,416 | 5.154% (5.008, 5.304) | 1,380,035,950 | 10,267 | 2.9 |
| CFP | 12 workers, 13 ranks on 12 cores | 3,241 | 3.783% (3.657, 3.913) | 1,370,501,354 | 10,213 | 2.8 |
| J5S | 1 worker | 453 | 0.529% (0.482, 0.580) | 32,579,814 | 7,566 | 2.1 |
| PS | 1 worker | 616 | 0.719% (0.665, 0.778) | 2,832,907,439 | 183,312 | 50.9 |
<!-- /table:fine1_whole_common -->

<!-- table:fine1_whole_cpu -->
| arm | wall_h | cpu_h | cores | evals_m |
|---|---|---|---|---|
| CS | 25.3 | 50.6 | 2.0 | 1,670.6 |
| CF | 20.0 | 40.0 | 2.0 | 1,362.5 |
| CSP | 2.9 | 31.3 | 11.0 | 1,380.0 |
| CFP | 2.8 | 31.1 | 11.0 | 1,370.5 |
| J5S | 2.1 | 2.1 | 1.0 | 32.6 |
| PS | 50.9 | 50.9 | 1.0 | 2,832.9 |
<!-- /table:fine1_whole_cpu -->

<!-- table:fine1_whole_pairs -->
| pair | shared | same | over_005 | over_1 | median | both | only_a | only_b |
|---|---|---|---|---|---|---|---|---|
| CS vs CSP | 86,816 | 49 | 48,838 | 3,910 | 0.063 | 2,670 | 2,048 | 1,775 |
| CF vs CFP | 86,237 | 50 | 44,874 | 2,903 | 0.055 | 1,588 | 1,168 | 1,683 |
| CS vs CF | 86,367 | 2,083 | 43,522 | 3,858 | 0.051 | 1,832 | 2,875 | 927 |
| CS vs PS | 87,435 | 5 | 66,041 | 4,609 | 0.101 | 392 | 4,367 | 231 |
| CS vs J5S | 87,435 | 2 | 61,393 | 4,610 | 0.088 | 309 | 4,450 | 155 |
| PS vs J5S | 98,280 | 23 | 61,855 | 744 | 0.069 | 266 | 448 | 295 |
| CSP vs CFP | 87,073 | 4,158 | 43,727 | 3,788 | 0.051 | 2,007 | 2,463 | 1,278 |
<!-- /table:fine1_whole_pairs -->
