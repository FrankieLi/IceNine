<!-- table:es_predictor -->
| sample | AUC logistic | AUC tree | AUC boosted | saving (logistic) | rescued stopped |
|---|---|---|---|---|---|
| held3 | 0.588 | 0.570 | not reproduced | 2.34% | 6 of 85 |
| held4 | 0.574 | 0.575 | not reproduced | 1.59% | 11 of 129 |
| held5 | 0.588 | 0.559 | not reproduced | 1.39% | 7 of 91 |
| held6 | 0.560 | 0.566 | not reproduced | 2.33% | 11 of 101 |
| held | 0.547 | 0.563 | not reproduced | 2.47% | 13 of 120 |
| mid | 0.573 | 0.528 | not reproduced | 2.10% | 4 of 19 |
| grad | 0.352 | 0.594 | not reproduced | 1.92% | 7 of 10 |
<!-- /table:es_predictor -->

<!-- table:es_r0 -->
| sample | R0 saving | lost 0.95 ends | changed ends |
|---|---|---|---|
| held3 | 12.42% | 62 (1.43% to 2.34%) | 173 |
| held4 | 15.79% | 99 (1.98% to 2.92%) | 220 |
| held5 | 10.28% | 71 (1.59% to 2.51%) | 164 |
| held6 | 12.98% | 77 (1.50% to 2.34%) | 190 |
| held | 15.19% | 93 (2.06% to 3.08%) | 207 |
| mid | 9.86% | 19 (0.75% to 1.82%) | 48 |
<!-- /table:es_r0 -->

<!-- table:es_sizing -->
| sample | steps | reach D | rescued | never | after D, never | after D, rescued | seed test if all stop |
|---|---|---|---|---|---|---|---|
| held3 | 3,391 | 686 | 85 | 601 | 14.48% | 1.65% | 100 (90) |
| held4 | 4,109 | 729 | 129 | 600 | 15.92% | 3.06% | 101 (70) |
| held5 | 3,556 | 638 | 91 | 547 | 11.49% | 1.66% | 86 (104) |
| held6 | 4,109 | 688 | 101 | 587 | 13.52% | 1.86% | 91 (63) |
| held | 3,685 | 625 | 120 | 505 | 15.05% | 2.71% | 77 (66) |
| mid | 1,626 | 149 | 19 | 130 | 9.55% | 1.22% | 8 (12) |
| grad | 1,281 | 108 | 10 | 98 | 2.61% | 0.20% | 6 (11) |
| big | 330 | 5 | 1 | 4 | 1.52% | 0.30% | 0 (0) |
<!-- /table:es_sizing -->
