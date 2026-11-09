# AR0: what the stored results say about stopping a search early

Step AR0 of `docs/autoresearch/PLAN.md` (2026-11-09). No new search was run. The script
`scripts/autoresearch/ar0/stopping.py` reads stored results and writes
`benchmarks/autoresearch/ar0/` (`summary.json`, `tables.md`, two plots). The analysis and the
numbers are the AR0 implementer's (commit f4a1654), revised after review (commit 3f3757d). This
page was written by the main session from those files, because the implementer's write of this
file was refused by the tool, and revised with the review fixes.

## Terms

- **Evaluation**: one call of the cost function for one voxel and one orientation. All counts on
  this page are raw counts. The level-0 screen counts interpolated evaluations (one per
  orientation) and adds them to strict ones, so the counts are not the strict-equivalent units
  the plan asks for.
- **Cost**: the strict cost of an orientation. Lower is better. 0 is a perfect match.
- **Right / wrong**: the returned orientation is within 1 degree of the truth, or not.
- **Level**: one pass of the adaptive grid search. Level 0 is the coarse grid; levels 1 to 3 refine
  the candidates kept from the level before. The finisher runs last.
- **Screened list**: at each level, the grid orientations that pass the screen and the spacing
  filter, in order of screen score. Every one of them gets a short Monte-Carlo refinement (the
  quick MC).
- **Kept**: the first few of the screened list after the quick MC, re-ordered by strict cost. They
  seed the next level. The **best candidate** is the first kept one.
- **True basin**: a level holds the true basin when its screened list has an orientation within the
  level's grid diameter of the truth. The diameters are 5.0, 3.33, 2.22 and 1.48 degrees at levels
  0 to 3. At level 0 that is wider than 2 degrees, so this count is looser than the "within
  2 degrees" counts below.
- **Hit ratio**: the share of an orientation's predicted spots that are lit in the data. Higher is
  better.
- **E0**: 200 hard voxels, each searched once with no starting guess. Stored per voxel in
  `benchmarks/eval_reduction/E/cache/rank/e0_old_base/` (the default overlap rule) and
  `e0_cpp_base/` (the C++ overlap rule).
- **Pilot**: the 2,000-voxel BFS run (key-off rerun), stored in
  `benchmarks/eval_reduction/G/g2/bfs/runs/keyoff_rerun/`.

## Scope: E0 is not a full-sample search

E0 voxels cost about 48,700 evaluations per search (9,739,905 over 200 voxels). A search with no
starting guess on a full-sample voxel costs 332,046 on average. That is the mean of the six
full-sample voxels of the Task G2 ledger (`benchmarks/eval_reduction/G/g2/ledger/runs/`,
`g2_base_v*.json`); the pilot's 69 full searches average 331,636. The two kinds of search spread
their evaluations very differently, and most of this page is measured on E0 only.

Reading the table: one row per stage of a search with no starting guess. `e0_evals_pct` is the
stage's share of E0's raw evaluations. `full_evals_pct` is its share for the six full-sample
voxels. `full_time_pct` is its share of their wall time (single-worker, contended: the timing
check before the run did not pass, so the time shares are a rough guide). E0's wall time per stage
was not recorded.

<!-- table:split_stages -->
| stage | e0_evals_pct | full_evals_pct | full_time_pct |
|---|---|---|---|
| level 0 screen and re-cost | 90.7 | 23.7 | 15.8 |
| level 0 quick MC | 4.5 | 37.1 | 46.0 |
| level 1 | 2.3 | 26.3 | 25.4 |
| level 2 | 0.7 | 8.8 | 8.4 |
| level 3 | 0.2 | 2.8 | 2.7 |
| finisher | 1.7 | 1.3 | 1.6 |
| not in a stage (time only) | n/a | n/a | 0.0 |
<!-- /table:split_stages -->

On E0 the level-0 screen and re-cost take 90.7% of the evaluations. On full-sample voxels they take
23.7%, and the level-0 quick MC is the largest stage. The screen alone is 13.2% of the full-sample
evaluations but 0.9% of their wall time, because its interpolated evaluations are cheap. So raw
counts make the screen look more costly than it is.

Reading the table: one row per level k. Each column is the share spent after level k (the later
levels and the finisher). That is the most a stop after level k could save.

<!-- table:split_after -->
| after_level | e0_evals_pct | full_evals_pct | full_time_pct |
|---|---|---|---|
| 0 | 4.8 | 39.2 | 38.2 |
| 1 | 2.6 | 12.9 | 12.7 |
| 2 | 1.9 | 4.1 | 4.3 |
| 3 | 1.7 | 1.3 | 1.6 |
<!-- /table:split_after -->

A stop after level 0 could save at most 4.8% on E0, but up to 39.2% on full-sample voxels. After
level 1 the ceilings are 2.6% and 12.9%; after level 2, 1.9% and 4.1%. Whether a stop rule is safe
and worth it on full-sample voxels was **not tested**. No per-level record exists for full-sample
searches: the pilot's `records.json` holds one line per voxel (`voxel_idx`, `source`, `n_evals`,
`wall_s`, `hit_ratio` and reject and revisit flags), and the G2 ledger holds counts and times per
level but no costs or distances. A per-level record on full-sample voxels needs a new run. It
belongs to step AR5.

## The answer in short

1. **The final cost tells right from wrong with no mistakes on E0.** All 159 right voxels end at a
   cost of 0.390 or less. All 41 wrong voxels end at 0.673 or more: 0 of 200 misclassified (Wilson
   95% interval 0.0-1.9%). The separating range, any cut from 0.40 to 0.65, was read on the same
   200 voxels; it was not checked on new ones. So after a search we can tell that the answer is
   good. That is not the same as knowing, during the search, that more search will not help.
2. **The best cost after a level does not tell right from wrong.** After every level the two groups
   overlap.
3. **On E0, a rule "stop after level k when the best cost is at most c" saves very little.** With
   no voxel broken (0 of 159 right voxels; 0.0-2.4%) it saves at most 1.26% of E0's raw
   evaluations if a stopped voxel skips the finisher (gross), and at most 0.37% if the finisher
   still runs (net). On E0 the reason is where the evaluations are: 90.7% are the level-0 screen
   and re-cost, which are over before the rule can act. On full-sample voxels only 23.7% are, and
   39.2% come after level 0, so the same kind of rule could save far more there. That is untested.
4. **On E0, the true basin is near the top of the level-0 screened list.** For the 159 right voxels
   it sits at position 3 (median), in the first tenth of the list for 9 voxels in 10. The list is
   182 long (median over all 200 voxels). On E0 a stop rule can only pay inside the level-0
   screen. Better still is to replace the screen with a proposer that names a few candidates
   (steps AR1 and AR4).
5. **On the pilot, the hit ratio flags every wrong answer at a cut of 0.93.** All 37 wrong voxels
   fall below it (0 of 37 above; 0.0-9.4%), with 15 right ones. All 1,948 answers at or above it
   are right (0 wrong; 0.0-0.2%). The cut was picked and read on the same run, so this describes
   the run; it is not a tested threshold.

## E0: cost by level, and what a stop rule saves

Reading the plot. Left: the best cost after each level and the final cost, one dot per voxel; blue
is right, red is wrong. The groups overlap until the final column, where a gap opens. Middle: the
best candidate's distance to the truth against its cost, by level. Below a cost of about 0.69
every best candidate is within 2 degrees of the truth; from there up, the cost alone does not say.
Best candidates more than 2 degrees off appear from a cost of 0.777 at level 0 (0.773 at level 1)
and from 0.690 at levels 2 and 3 (`far_cost_min` in the second table). Right: for the stop rule,
the raw evaluations saved against the voxels broken (more than 1 degree off), as the threshold `c`
rises; solid lines are the gross saving, dashed lines the net saving (both defined below).

![E0 stopping analysis](../../benchmarks/autoresearch/ar0/fig_stopping.png)

Reading the table: one row per threshold `c` on the final cost. `right_le` and `wrong_le` count the
right and wrong voxels at or below `c`; `right_gt` and `wrong_gt` count those above. Any `c` from
0.40 to 0.65 splits the groups with no mistake on these 200 voxels. The range was read on the same
voxels it is applied to.

<!-- table:final_sep_old -->
| c | right_le | wrong_le | right_gt | wrong_gt |
|---|---|---|---|---|
| 0.20 | 149 | 0 | 10 | 41 |
| 0.30 | 153 | 0 | 6 | 41 |
| 0.40 | 159 | 0 | 0 | 41 |
| 0.50 | 159 | 0 | 0 | 41 |
| 0.60 | 159 | 0 | 0 | 41 |
| 0.65 | 159 | 0 | 0 | 41 |
| 0.70 | 159 | 4 | 0 | 37 |
| 0.75 | 159 | 24 | 0 | 17 |
| 0.80 | 159 | 34 | 0 | 7 |
| 0.85 | 159 | 36 | 0 | 5 |
| 0.90 | 159 | 39 | 0 | 2 |
| 0.95 | 159 | 41 | 0 | 0 |
<!-- /table:final_sep_old -->

Reading the table: one row per level, default rule. `best_le1_right` counts the right voxels whose
best candidate after the level is within 1 degree of the truth (of 159). `best_le2_right` and
`best_le2_wrong` count voxels whose best candidate is within 2 degrees, for the voxels that end
right and wrong. `kept_le2_right` and `kept_le2_wrong` count voxels with any kept candidate within
2 degrees. `far_cost_min` is the lowest cost of a best candidate more than 2 degrees off (all 200
voxels). `evals_after_pct` is the share of all raw evaluations spent after the level.
`evals_after_right_mean` is the mean number of evaluations a right voxel spends after the level.

<!-- table:levels_old -->
| level | best_le1_right | best_le2_right | best_le2_wrong | kept_le2_right | kept_le2_wrong | far_cost_min | evals_after_pct | evals_after_right_mean |
|---|---|---|---|---|---|---|---|---|
| 0 | 71 | 74 | 1 | 156 | 14 | 0.777 | 4.8 | 2411 |
| 1 | 141 | 141 | 1 | 159 | 7 | 0.773 | 2.6 | 1292 |
| 2 | 148 | 148 | 0 | 159 | 1 | 0.690 | 1.9 | 967 |
| 3 | 155 | 155 | 0 | 155 | 0 | 0.690 | 1.7 | 886 |
<!-- /table:levels_old -->

Two things to read from it. First, after level 0 the best candidate is within 1 degree for only 71
of the 159 right voxels; it takes level 1 to reach 141. Second, the search loses voxels it could
have solved: 14 wrong voxels still hold a kept candidate within 2 degrees after level 0, 7 after
level 1, 1 after level 2 and none after level 3. For those voxels the search later dropped a
candidate within 2 degrees. No stop rule helps them. A better ranking might: the smooth ranking key
of Task E fixed 34 of the 41 wrong voxels, which is consistent with this, but we did not check that
the 34 include these 14.

Reading the table: one row per level, default rule. `diameter_deg` is the level's grid diameter.
`screen_pct` and `quick_pct` are the shares of all raw evaluations spent in the level's screen and
re-cost, and in its quick Monte-Carlo. `n_disc_median` is the median length of the level's screened
list (all 200 voxels). `right_with_basin` and `wrong_with_basin` count voxels whose screened list
holds the true basin (an orientation within the grid diameter of the truth). `basin_pos_median` is
the true basin's position in the screened list (0 is the top), for the right voxels that have one;
`wrong_basin_pos_median` is the same for the wrong voxels. `basin_frac_median` and `basin_frac_p90`
are the right voxels' position as a share of the list: the median and the value that 9 voxels in
10 stay under.

<!-- table:stages_old -->
| level | diameter_deg | screen_pct | quick_pct | n_disc_median | right_with_basin | right_no_basin | wrong_with_basin | basin_pos_median | wrong_basin_pos_median | basin_frac_median | basin_frac_p90 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 5.00 | 90.7 | 4.5 | 182 | 159 | 0 | 41 | 3 | 6 | 0.02 | 0.10 |
| 1 | 3.33 | 1.1 | 1.2 | 48 | 159 | 0 | 28 | 0 | 8 | 0.00 | 0.05 |
| 2 | 2.22 | 0.4 | 0.3 | 12 | 159 | 0 | 7 | 0 | 3 | 0.00 | 0.00 |
| 3 | 1.48 | 0.1 | 0.1 | 3 | 159 | 0 | 1 | 0 | 2 | 0.00 | 0.00 |
<!-- /table:stages_old -->

The finisher takes the remaining 1.7%. All 41 wrong voxels hold the true basin in their level-0
screened list (median position 6). That is consistent with a ranking failure rather than a screen
that missed the basin. The level-0 diameter (5.0 degrees) is wider than 2 degrees, which is why
this count (41) is not the 14 wrong voxels with a kept candidate within 2 degrees in the table
above: the two counts use different lists and different radii.

Reading the table: one row per stop rule "stop after `level` when the best cost is at most `c`".
`n_stop` is the number of voxels the rule stops. `broken_1deg` and `broken_2deg` count voxels that
end right today and whose best candidate at the stop is more than 1 or 2 degrees from the truth.
`stop_wrong_base` counts stopped voxels that end wrong today anyway. The saving is in raw
evaluations, as a share of E0's total, counted two ways. `saved_gross_pct` assumes a stopped voxel
returns its best candidate with no finisher, so that candidate must already be within 1 degree:
read it with `broken_1deg`. `saved_net_pct` assumes the finisher still runs at its usual cost and
only the later levels are saved: read it with `broken_2deg`. That the Monte-Carlo finisher reaches
the truth from within 2 degrees is untested. In this grid the two broken counts are equal on every
row.

<!-- table:stop_grid_old -->
| level | c | n_stop | broken_1deg | broken_2deg | stop_wrong_base | saved_gross_pct | saved_net_pct |
|---|---|---|---|---|---|---|---|
| 0 | 0.30 | 2 | 0 | 0 | 0 | 0.05 | 0.03 |
| 0 | 0.40 | 2 | 0 | 0 | 0 | 0.05 | 0.03 |
| 0 | 0.50 | 2 | 0 | 0 | 0 | 0.05 | 0.03 |
| 0 | 0.60 | 5 | 0 | 0 | 0 | 0.12 | 0.08 |
| 0 | 0.70 | 9 | 0 | 0 | 0 | 0.21 | 0.13 |
| 0 | 0.80 | 28 | 1 | 1 | 0 | 0.70 | 0.44 |
| 0 | 0.90 | 83 | 27 | 27 | 3 | 2.06 | 1.30 |
| 1 | 0.30 | 6 | 0 | 0 | 0 | 0.08 | 0.02 |
| 1 | 0.40 | 7 | 0 | 0 | 0 | 0.09 | 0.03 |
| 1 | 0.50 | 19 | 0 | 0 | 0 | 0.25 | 0.08 |
| 1 | 0.60 | 38 | 0 | 0 | 0 | 0.51 | 0.16 |
| 1 | 0.70 | 61 | 0 | 0 | 0 | 0.80 | 0.25 |
| 1 | 0.80 | 107 | 2 | 2 | 0 | 1.42 | 0.45 |
| 1 | 0.90 | 170 | 16 | 16 | 19 | 2.22 | 0.70 |
| 2 | 0.30 | 16 | 0 | 0 | 0 | 0.15 | 0.01 |
| 2 | 0.40 | 30 | 0 | 0 | 0 | 0.29 | 0.02 |
| 2 | 0.50 | 55 | 0 | 0 | 0 | 0.54 | 0.05 |
| 2 | 0.60 | 85 | 0 | 0 | 0 | 0.84 | 0.07 |
| 2 | 0.70 | 127 | 0 | 0 | 1 | 1.26 | 0.11 |
| 2 | 0.80 | 155 | 7 | 7 | 6 | 1.54 | 0.13 |
| 2 | 0.90 | 192 | 11 | 11 | 33 | 1.85 | 0.16 |
<!-- /table:stop_grid_old -->

Reading the table: the best rule per level under a cap on broken voxels (`max_broken` 0 or 1,
counted at 1 degree): the largest gross saving; on a tie, the smaller `c`. The cap and `c` were
chosen on the same 200 voxels, so "0 broken" is an in-sample result (0 of 159; Wilson 95%
interval 0.0-2.4%).

<!-- table:safe_old -->
| max_broken | level | c | n_stop | broken_1deg | broken_2deg | saved_gross_pct | saved_net_pct |
|---|---|---|---|---|---|---|---|
| 0 | 0 | 0.75 | 17 | 0 | 0 | 0.39 | 0.25 |
| 0 | 1 | 0.75 | 88 | 0 | 0 | 1.16 | 0.37 |
| 0 | 2 | 0.70 | 127 | 0 | 0 | 1.26 | 0.11 |
| 1 | 0 | 0.80 | 28 | 1 | 1 | 0.70 | 0.44 |
| 1 | 1 | 0.75 | 88 | 0 | 0 | 1.16 | 0.37 |
| 1 | 2 | 0.70 | 127 | 0 | 0 | 1.26 | 0.11 |
<!-- /table:safe_old -->

Reading the table: the thresholds chosen on the default-rule results, applied unchanged to the
C++-rule results. These are the same 200 voxels and images, searched with a different overlap
rule; they are not new voxels, so this is only a weak check against fitting. The thresholds break
no voxel there either (0 of 144; 0.0-2.6%), and save as little.

<!-- table:cross_old_to_cpp -->
| level | c | n_stop | broken_1deg | broken_2deg | saved_gross_pct | saved_net_pct |
|---|---|---|---|---|---|---|
| 0 | 0.75 | 14 | 0 | 0 | 0.27 | 0.16 |
| 1 | 0.75 | 64 | 0 | 0 | 0.62 | 0.17 |
| 2 | 0.70 | 108 | 0 | 0 | 0.83 | 0.06 |
<!-- /table:cross_old_to_cpp -->

The C++-rule results give the same picture: 144 right voxels all at a final cost of 0.338 or less,
56 wrong voxels all at 0.678 or more, and a largest safe saving of 0.90% gross (0.17% net at
most). Their tables are in
`benchmarks/autoresearch/ar0/tables.md` (`final_sep_cpp`, `levels_cpp`, `stages_cpp`,
`stop_grid_cpp`, `safe_cpp`).

## Pilot: the hit ratio of the returned orientation

Reading the plot: the hit ratio of the returned orientation in the pilot's key-off rerun, for full
searches (left) and local fits (right); blue is right, red is wrong. The counts are on a log
scale.

![Pilot hit ratio](../../benchmarks/autoresearch/ar0/fig_pilot_hit_ratio.png)

Reading the table: one row per pilot run and kind of search. A `full_search` voxel is one that was
a seed or a reseed (from the records: source "seed", a rejected seed, or reseeded); every other
voxel is a `local_fit`. In the OutwardSearch runs (`outward_only`, `integrated`) a reseed can be
answered by the policy without a full search, so their `full_search` rows include a few voxels with
no full search. An earlier version of this table picked full searches by their evaluation count,
which miscounted those two runs; this is corrected. `n` and `n_wrong` are the voxel counts. The
next four columns are the lowest and the median hit ratio of the right voxels, and the median and
the highest hit ratio of the wrong voxels. `evals_in_wrong_pct` is the share of the group's
evaluations spent on voxels that end wrong.

<!-- table:pilot_hit -->
| run | group | n | n_wrong | hit_right_min | hit_right_median | hit_wrong_median | hit_wrong_max | evals_in_wrong_pct |
|---|---|---|---|---|---|---|---|---|
| keyoff_rerun | full_search | 69 | 4 | 0.936 | 1.000 | 0.777 | 0.904 | 5.8 |
| keyoff_rerun | local_fit | 1931 | 33 | 0.888 | 0.994 | 0.904 | 0.918 | 1.2 |
| outward_only | full_search | 70 | 6 | 0.951 | 1.000 | 0.875 | 0.892 | 2.1 |
| outward_only | local_fit | 1930 | 30 | 0.900 | 0.992 | 0.902 | 0.918 | 0.9 |
| integrated | full_search | 70 | 4 | 0.964 | 1.000 | 0.879 | 0.899 | 0.4 |
| integrated | local_fit | 1930 | 38 | 0.900 | 0.994 | 0.902 | 0.918 | 1.4 |
| smooth_only | full_search | 65 | 0 | 0.921 | 1.000 | n/a | n/a | 0.0 |
| smooth_only | local_fit | 1935 | 34 | 0.859 | 0.992 | 0.902 | 0.918 | 1.2 |
<!-- /table:pilot_hit -->

Reading the tables: the key-off rerun, one row per hit-ratio cut `h`. `right_ge` and `wrong_ge`
count right and wrong voxels at or above the cut; `right_lt` and `wrong_lt` count those below.
First the 69 full searches, then the 1,931 local fits.

<!-- table:pilot_thr_full_search -->
| h | right_ge | wrong_ge | right_lt | wrong_lt |
|---|---|---|---|---|
| 0.90 | 65 | 1 | 0 | 3 |
| 0.93 | 65 | 0 | 0 | 4 |
| 0.95 | 64 | 0 | 1 | 4 |
| 0.96 | 62 | 0 | 3 | 4 |
| 0.97 | 59 | 0 | 6 | 4 |
| 0.98 | 57 | 0 | 8 | 4 |
| 0.99 | 44 | 0 | 21 | 4 |
<!-- /table:pilot_thr_full_search -->

<!-- table:pilot_thr_local_fit -->
| h | right_ge | wrong_ge | right_lt | wrong_lt |
|---|---|---|---|---|
| 0.90 | 1896 | 28 | 2 | 5 |
| 0.93 | 1883 | 0 | 15 | 33 |
| 0.95 | 1852 | 0 | 46 | 33 |
| 0.96 | 1789 | 0 | 109 | 33 |
| 0.97 | 1590 | 0 | 308 | 33 |
| 0.98 | 1393 | 0 | 505 | 33 |
| 0.99 | 1076 | 0 | 822 | 33 |
<!-- /table:pilot_thr_local_fit -->

At a cut of 0.93 no wrong voxel passes in either group. Below it are all 37 wrong voxels (4 full
searches, 33 local fits) and 15 right local fits. The wrong local fits sit in a narrow band: 28 of
the 33 have a hit ratio from 0.90 to 0.918.

The BFS acceptance rule is relative, not a fixed 0.9. A neighbour's local fit is accepted when its
hit ratio divided by the best hit ratio so far in the BFS wave (starting from the seed's) is above
0.9 (C++ `Src/BreadthFirstReconstructor.tmpl.cpp`, lines 156 and 162-163; Python
`icenine/reconstructor.py`, lines 2121-2125). The best hit ratio so far is usually close to 1
(right voxels have a median of 0.994 or more), so local fits much below 0.9 are rejected. The narrow band is consistent with
selection by this rule, and that is the likely reason a cut just above it works on local fits. Not
checked: the best-so-far value at each accepted wrong fit.

The cut of 0.93 was picked and read on this run and was not tested on a second region. Read it as
a description of this run, not as a tuned threshold. In the other three pilot runs no wrong voxel
has a hit ratio above 0.918 either, but they cover the same region, so that is not an independent
check.

## Intervals on the zero counts

Reading the table: one row per count of zero that the page reports. `k` of `n`, and the Wilson 95%
interval in percent (`scripts/common/stats.py`). All of them are in-sample: the thresholds were
chosen on the same voxels.

<!-- table:zero_counts -->
| count | k | n | lo_pct | hi_pct |
|---|---|---|---|---|
| E0 old rule: voxels misclassified by the final cost (best c) | 0 | 200 | 0.0 | 1.9 |
| E0 C++ rule: voxels misclassified by the final cost (best c) | 0 | 200 | 0.0 | 1.9 |
| E0 old rule: broken at 1 degree, safe rule after level 0 | 0 | 159 | 0.0 | 2.4 |
| E0 old rule: broken at 1 degree, safe rule after level 1 | 0 | 159 | 0.0 | 2.4 |
| E0 old rule: broken at 1 degree, safe rule after level 2 | 0 | 159 | 0.0 | 2.4 |
| E0 C++ rule: broken at 1 degree, old-rule c after level 0 | 0 | 144 | 0.0 | 2.6 |
| E0 C++ rule: broken at 1 degree, old-rule c after level 1 | 0 | 144 | 0.0 | 2.6 |
| E0 C++ rule: broken at 1 degree, old-rule c after level 2 | 0 | 144 | 0.0 | 2.6 |
| pilot key-off: wrong voxels at or above the cut | 0 | 37 | 0.0 | 9.4 |
| pilot key-off: answers at or above the cut that are wrong | 0 | 1948 | 0.0 | 0.2 |
<!-- /table:zero_counts -->

## What this means for the plan

- On E0, a cost threshold between levels is not the early stop the owner asked for. It breaks no
  voxel on these 200, and it saves at most 1.26% of the raw evaluations gross, 0.37% net. On
  full-sample voxels 39.2% of the raw evaluations come after level 0, so a stop after level 0 could
  save far more there. Whether it is safe there is untested. Testing it needs the best cost and
  distance after each level on full-sample voxels, which is a new run and belongs to step AR5.
- On E0 the evaluations are in the level-0 screen. On full-sample voxels they are spread over the
  level-0 quick MC (37.1%), the level-0 screen and re-cost (23.7%) and levels 1 to 3. The level-0
  quick MC runs once per screened candidate, so a shorter screened list would cut it too; that is
  arithmetic, not a measured result. Two ways to get at the level-0 work: stop the screen once a
  candidate in a clear basin has shown up (on E0 the true basin is near the top of the list, but
  the list is ranked only after the whole grid is screened, so this needs an order to screen in);
  or do not screen the grid at all and let a proposer name the candidates (AR1, AR4).
- "Is this answer good" has a clean test after refinement: the final cost on E0 (measured after
  the finisher) and the hit ratio on the pilot (measured after the local fit or full search).
  Using the same test on a proposer's unrefined candidates, as the stop rule of AR5 (stop as soon
  as a candidate passes, search further only when none does), is a hypothesis. The untested part:
  before the finisher, candidate costs overlap at every level (headline 2), so it is not known how
  much refinement a candidate needs before the test applies.
- Not tested here: whether these thresholds hold on full-sample voxels with many more spots, and
  on images with noise.
