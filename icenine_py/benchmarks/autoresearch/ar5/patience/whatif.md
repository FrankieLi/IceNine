What-if tables of AR5 step 11 (arithmetic on the stored arm-G traces of step 9; a cost
estimate, not a measurement). Written by `patience_whatif.py`.

Reading line, `patience_whatif_roles`: one row per role of an MC refinement in its seed
step; each cell is refinements, then their evaluations as a share of G's run with the
pass, then (in brackets) the share spent after the refinement's last improvement.

Reading line, `patience_whatif_rules`: one row per stop rule; each cell is the saving
as a share of G's run with the pass, first WITH the estimated extra candidates (a step
that loses its 0.95 stop is assumed to run its remaining candidates, each priced at
that stored step's mean evaluations per candidate), then WITHOUT them (the plan's
basis), then the seed steps whose taken end changes. The extra-candidate price includes
the long refinements the rule cuts, so it is likely to overprice them; rule (iii)'s gap to the
plan comes from this pricing (its share without extras matches the plan's). Neither
share is a measurement: the gate's live replay of the whole changed path (`gate.json`,
every extra evaluation counted) supersedes both. The plan's table (ii) (best of 120
cost-gated rules, the in-sample checkpoint rule) is not rebuilt: its definitions are not
given to the digit.

<!-- table:patience_whatif_roles -->
| role | mid | grad | big |
|---|---|---|---|
| reached 0.95, taken | 410; 17.6% (14.5%) | 137; 2.5% (2.1%) | 45; 8.7% (6.9%) |
| taken in a step with no 0.95 end | 102; 2.2% (1.3%) | 66; 1.1% (0.7%) | 3; 0.3% (0.2%) |
| not taken (late gains idle) | 1,093; 17.5% (12.4%) | 788; 7.1% (5.9%) | 105; 6.6% (4.9%) |
<!-- /table:patience_whatif_roles -->

<!-- table:patience_whatif_rules -->
| rule | mid | grad | big |
|---|---|---|---|
| P(1000) | 19.6%; 22.5%; 90 | 1.7%; 1.9%; 14 | 6.2%; 10.9%; 8 |
| P(2000) | 14.2%; 15.1%; 39 | 0.4%; 0.4%; 2 | 7.1%; 8.0%; 6 |
| P(3000) | 10.9%; 11.4%; 25 | 0.2%; 0.2%; 1 | 5.2%; 6.1%; 4 |
| P(5000) | 6.7%; 7.2%; 15 | 0.0%; 0.0%; 0 | 4.5%; 4.5%; 3 |
| (iii) 500 while worse | 3.4%; 9.0%; 61 | 2.0%; 2.2%; 11 | -0.4%; 4.4%; 8 |
<!-- /table:patience_whatif_rules -->
