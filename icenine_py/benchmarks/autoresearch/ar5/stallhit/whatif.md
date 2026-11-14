What-if tables of AR5 step 12 (arithmetic on the stored arm-G traces of step 9 and step
11's arm S; a cost estimate, not a measurement). Written by `stallhit_whatif.py`.

Reading line, `stallhit_whatif_rules`: one row per stop rule; each cell is the share of
arm G's run (with the pass) saved WITHOUT the estimated extra candidates of a step that
loses its 0.95 stop (the plan's basis; `whatif.json` also holds the share with them),
then the seed steps whose taken end changes (a range where some stall hit ratios are not
on record: from the plain reading to the reading where every unknown stall at a G end of
0.95 or more stops). Step 11's row is the measured whole run
(`patience/stats_*.json`).

Reading line, `stallhit_whatif_cases`: refinements, then the evaluations they save, then
that as a share of G's run. Cases are exclusive and cover the refinements whose G end is
at 0.95 or more and that have a stall at k = 3000 (see the script text).

Reading line, `stallhit_whatif_k`: other k, gated high; each cell is the exact saving
(refinements cut after their last gain, G end at 0.95 or more), then the refinements cut
before a later gain with a G end at 0.95 or more (the count of ends that could change).

<!-- table:stallhit_whatif_cases -->
| case | mid | grad | big |
|---|---|---|---|
| A: cut after the last gain, G end at 0.95 or more (end unchanged) | 178; 893,957 (6.87%) | 12; 16,504 (0.16%) | 24; 145,105 (3.81%) |
| B: first stall known under 0.95, stops at a later stall | 7; 52,879 (0.41%) | 0; 0 (0.00%) | 2; 10,110 (0.27%) |
| C: first stall known at 0.95, stops before a later gain | 7; 32,874 (0.25%) | 0; 0 (0.00%) | 1; 17,292 (0.45%) |
| D: first stall unknown, G end at 0.95 or more | 4; 92,765 (0.71%) | 0; 0 (0.00%) | 0; 0 (0.00%) |
| E: first stall known under 0.95, no stop after it | 4; 0 (0.00%) | 0; 0 (0.00%) | 1; 0 (0.00%) |
<!-- /table:stallhit_whatif_cases -->

<!-- table:stallhit_whatif_k -->
| k | mid | grad | big |
|---|---|---|---|
| 2000 | 8.40%; 32 | 0.33%; 0 | 3.97%; 5 |
| 3000 | 6.87%; 22 | 0.16%; 0 | 3.81%; 4 |
| 5000 | 4.78%; 14 | 0.04%; 0 | 2.99%; 3 |
<!-- /table:stallhit_whatif_k -->

<!-- table:stallhit_whatif_rules -->
| rule (k = 3000) | mid | grad | big |
|---|---|---|---|
| `P(3000)`, step 11, measured whole run (changed) | 10.60%; 20 | 0.14%; 1 | 5.01%; 4 |
| gated high: cut only if the best is at 0.95 or more | 7.7% to 8.2%; 7 to 11 | 0.2%; 0 | 4.5%; 1 |
| gated low: cut only if the best is under 0.95 | 3.6%; 14 | 0.0%; 1 | 1.9%; 3 |
<!-- /table:stallhit_whatif_rules -->
