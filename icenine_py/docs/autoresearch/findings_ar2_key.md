# AR2 step 2: does indexing work as the BFS seed search?

Step 2 of "AR2" in `docs/autoresearch/PLAN.md` (2026-11-09), under D-20261109-N5 and the finisher arm that step 1 fixed (`findings_ar2_holdout.md`, D-20261109-N7 and N8). Library: the opt-in config key `SeedByIndexing` in `icenine/config_file.py`, `icenine/orientation_search.py` and `icenine/reconstructor.py`, and `seed_by_indexing` in `icenine/ray_indexing.py`. Scripts: `scripts/autoresearch/ar2/key_*.py`. Results: `benchmarks/autoresearch/ar2/key/`. Branch `feature/eval-reduction-ar2b`. The code, the runs and the text of this page are the implementer's. The main session placed the page from that text, because the tool refused the implementer's own write of the page; it filled the tables with `scripts/dev/sync_doc_tables.py` and ran `scripts/dev/audit_numbers.py` (929 of 930 numbers matched; the one left over is the date inside the decision name D-20261027-I2). The page was revised after review (review report 20261109-review-ar2b, decisions D-20261109-N9 and N10): the trace now covers both sides of each pair, the FusedCost differences are measured, and the claims are narrowed. All timings on this page are contended: the timing preflight refused (an editor process used about 200% of a CPU core), so the runs went ahead and are labelled.

## The result in one picture

Reading the picture. Left: cost evaluations of the 2,000-voxel BFS pilot, in millions, for the stored key-off run (grid search for every seed and reseed) and for the key-on run. Blue is the grid searches, orange the indexing path, green the local fits. Right: three paired comparisons, voxel by voxel. A voxel is wrong when its result is more than 1 degree from the true orientation. "Fixed" is wrong in the first run of the pair and right in the second. "Broken" is the reverse. The first group is the key-on run against the stored key-off run. The second group is the same comparison with `FusedCost` on in both runs. The third group compares two key-off runs that differ only in `FusedCost`. It shows how much a key-off run moves without the key. Read the key's split beside that third group.

![Evaluations and wrong voxels, key on against key off, and the key-off noise pair](../../benchmarks/autoresearch/ar2/key/ar2b_key.png)

- With the key on, the pilot costs 2,819,906 evaluations. The stored key-off run costs 25,547,935. That is about 9.1 times fewer. The grid searches were 22,882,865 of the key-off total. The 64 indexing paths cost 122,198. The local fits cost about the same in both runs (2,697,708 on, 2,665,070 off).
- Accuracy of the indexed seeds themselves (the key's own path, D-20261109-N9: the chain on the seed voxel itself). None of the 64 indexed seeds and reseeds is wrong: 0 of 64 (0.0 to 5.7%). The grid search had 2 wrong of the 62 voxels whose result is their own full search (0.9 to 11.0%). The intervals overlap, so the two are not shown different.
- Flips on the key's own path (D-20261109-N10). Fixes: 2 in the stored pair (voxels 5876 and 24334) and 1 in the FusedCost-on pair (voxel 16476). Breaks: 0 in both pairs. Each of these fixed voxels was a seed in both runs. The seed test rejected its grid result, and its indexed result is right. Under N10 these fixes count for the key. They are counts, not a rate.
- Whole-run wrong voxels: 22 of 2,000 with the key on, 37 for the stored key-off run, 24 for the key-off run with `FusedCost` on. Paired: 22 fixed, 7 broken, 15 wrong in both against the stored run; 14 fixed, 12 broken, 10 wrong in both with `FusedCost` on in both runs. The two key-off runs alone differ by 20 fixed and 7 broken. So the whole-run wrong count is not shown different from the grid search. No whole-run gain in wrong voxels is claimed (N9).
- Every fixed and broken voxel was traced on its wrong side, from the stored run files: the breaks in the key-on run, the fixes in the key-off run. The trace is an orientation match (see "The flips, traced"). All 7 breaks (12 with `FusedCost` on) and all 22 fixes (14) take the orientation of a grain that touches the voxel. Apart from the 3 fixes above, these are boundary crossings. Under N6 and N9 they count neither for nor against the key, on both sides. No break has a wrong indexed seed behind it.
- The key changes which voxels become seeds, so the BFS takes a different route over the layer. That effect is part of the whole-run counts above. It was not tested on its own.
- No seed or reseed fell back to the grid search: 0 of 64 (0.0 to 5.7%). So the fallback path ran on no real voxel. It is tested only with fakes.
- Costs outside the evaluation count: pairing the rays once took 17.6 s. Indexing took 0.186 s per voxel on average (11.9 s for all 64).
- Wall time, a first contended reading only (N9): 213 s with the key on and 1633 s with it off, about 7.7 times shorter, with `FusedCost` on, one worker and one run per arm. The preflight had refused. The quiet timing of review item 14 is still open. With `FusedCost` off the key-on run took 1272 s, and it overlapped another run.
- The stored baseline numbers were recomputed from the stored run files before anything else (no rerun): 37 wrong, 25,547,935 evaluations, 69 full searches.

## Terms

- **BFS**: breadth-first reconstruction. A seed voxel is searched in full. Its result is passed to its neighbours, which only get a local fit.
- **Seed**: a voxel that gets a full search because no fitted neighbour has reached it. **Reseed**: a full search on one voxel of a cluster of voxels the BFS could not fit.
- **Full search**: the grid search of the library (coarse grid, quick Monte Carlo, the finisher). It costs 331,636 evaluations on average here. **Local fit**: a short search from the neighbour's orientation, 658 evaluations on average in the key-off run.
- **Evaluation**: one call of the strict cost. Attaching rays and indexing are not evaluations.
- **Strict cost**: one minus the share of predicted spots that land on lit pixels. Lower is better.
- **Indexing**: proposing orientations from the observed rays (`findings_ar1_ray_indexing.md`). A **ray** is a spot seen on both detectors. **Pairing** matches the spots of the two detectors into rays. It is done once per run.
- **Candidate**: an orientation the indexer proposes for a voxel, ranked by the rays it explains.
- **Finisher**: the usual local refinement from one candidate (FindOptimal, then VarianceMinimizing).
- **The key's own path** (D-20261109-N9): the chain on the seed voxel itself. Index the voxel, take the strict cost of the first 5 candidates, run the finisher on the 2 with the lowest strict cost, keep the lower final cost, apply the seed test. This is the arm of step 1 ("best 2"). A break on this path counts against the key (N7, N8). A fix on this path counts for it (N10).
- **Seed test**: a result is accepted when its hit ratio is at least `MinAccelerationThreshold` (0.80 here). It is the same test as for a grid-search seed. A **rejected seed** is a seed or reseed whose result failed it.
- **Fallback**: no candidate, or the result fails the seed test, so the grid search runs.
- **Wrong**: more than 1 degree from the true orientation (symmetry taken into account).
- **Fixed, broken, wrong in both**: the paired split by voxel of the key-on run against the key-off run.
- **Boundary crossing**: a voxel whose wrong result is the orientation of a grain that touches it (D-20261109-N6). If local fits from a seed that is right for its own voxel carried that orientation there, it counts neither for nor against the key. This holds on both sides, for broken and for fixed voxels (N9). It does not excuse a break on the key's own path (N7, N8).
- **Match** (in the trace): a voxel of the same run whose stored result is its own full search and lies within 1 degree of the traced voxel's result. "Nearest" means nearest in orientation, not in space.
- **Edge-sharing, corner-only**: the voxels are triangles. A voxel is edge-sharing if a voxel of another grain shares a whole edge with it. It is corner-only if the other grain meets it only at a corner.
- **Contended**: other work was running on the machine. **FusedCost**: an opt-in key that evaluates the strict cost in one C call. It is not bit-for-bit: it rounds in float32 at pixel edges (D-20261027-I2).

## What was built

- The key `SeedByIndexing` (0 or 1, default 0). With it on, every seed and reseed of the BFS is indexed first. The pairing is built once per run (`AdaptiveVoxelReconstructor.index_seed`).
- The shared logic is `ray_indexing.seed_by_indexing`. Its seed logic is the same as in `holdout_run.py` and `ar1/e2e.py`: the same candidate counts, the same stable sort, `default_rng(0)` for each finisher run, and the same tie rule. The numbers of candidates costed (5) and refined (2) are function parameters, not config keys. The finisher uses its own `default_rng(0)`, so it does not use the BFS random stream. The ray set differs from the scripts'. The library pairs the rays with the hull of the run's sample, the scripts with the hull of the full sample (see "How the indexing here differs from AR1's").
- The result goes through the seed test. If there is no candidate, or the test fails, the grid search runs. All cost evaluations of the indexing path count in the seed's evaluations (the reseed's, for a reseed), used or not. `stats["index_*"]` and `BFSReconstruction.index_log` split them: used, wasted, grid fallback, indexing time, and the one-time pairing.
- With the key off the code path is not entered. Tests (`tests/test_seed_by_indexing.py`, 12 tests):
  - the key is parsed and is off by default;
  - with the key absent and with it written off, a small scripted BFS run gives the same orientations, evaluation counts, records and stats;
  - the scripted 14 by 14 grid still matches its golden;
  - with the key on, the indexing result is taken and the grid search is skipped;
  - with no candidate, and with a failed seed test, the grid search runs and its evaluations are counted;
  - reseeds are indexed and counted in the reseed bucket;
  - on the ThreeVoxels toy, the library path equals the script path written out by hand: candidates, strict costs, finisher costs, evaluations and orientation are all exactly equal.

  Every existing test is untouched. The key-off test uses a fake reconstructor, not real images. On real images the key-off path is guarded by the same `if` on the key. A real-image key-off run was not compared with an older commit.
- On real voxels (`key_path_check.py`, in the BFS pilot's setting) the library path equals the script path written out by hand on 6 of 6 voxels. It also equals the entry stored in the key-on run's log. Both this check and the toy test build their "script path" from the library's own context (`_index_context`: the paired rays and the crystal). So they check the steps after the context is built, not the context itself (the hull and the crystal symmetry).

## What was run

Clean full-sample images of `Examples/Example2.ManyGrains`, Q-max 8, the 2,000-voxel pilot region, BFS random seed 0, one thread per run. The config of every run is the stored key-off pilot's (`ReconstructPhaseD_mc_clean.config` with `BFSReseed 1`) plus the key (and, where stated, `FusedCost 1`). Runs, one after the other: key on (FusedCost off), key on (FusedCost on), key off (FusedCost on). The last is the key-off arm of the wall-time pair that review item 14 asks for. Item 14 asks for timing on a quiet machine. This run went ahead after the timing preflight refused, so it does not close item 14 (N9). The key-on FusedCost-off run was not a single run. It overlapped with the voxel 17361 run and with short checks whose times were not logged. The key-off FusedCost-off arm is the stored rerun `benchmarks/eval_reduction/G/g2/bfs/runs/keyoff_rerun`.

The first table is the plain comparison. Reading the table: `wrong` is the count out of 2,000 with a Wilson 95% interval. `fixed`, `broken` and `both_wrong` are paired by voxel against the key-off run of the same FusedCost setting (the first two rows are the stored rerun and its key-on partner, the last two are the FusedCost-on pair). `unresolved` is the voxels no fit accepted.

<!-- table:ar2b_accuracy -->
| run | wrong | fixed | broken | both_wrong | unresolved |
|---|---|---|---|---|---|
| key off (stored rerun) | 37 of 2000 (1.3 to 2.5%) |  |  |  | 2 |
| key on | 22 of 2000 (0.7 to 1.7%) | 22 | 7 | 15 | 0 |
| key off, FusedCost on | 24 of 2000 (0.8 to 1.8%) |  |  |  | 1 |
| key on, FusedCost on | 22 of 2000 (0.7 to 1.7%) | 14 | 12 | 10 | 0 |
<!-- /table:ar2b_accuracy -->

Reading the table: how much the key-off run moves when only `FusedCost` is switched on, paired by voxel. `wrong_a` is the stored key-off rerun, `wrong_b` the FusedCost-on key-off run.

<!-- table:ar2b_noise -->
| pair | wrong_a | wrong_b | fixed | broken | both_wrong |
|---|---|---|---|---|---|
| key off with FusedCost on, against key off without | 37 | 24 | 20 | 7 | 17 |
<!-- /table:ar2b_noise -->

Switching on `FusedCost` alone moves 27 voxels (20 fixed, 7 broken) of a key-off run. That is the size of movement to keep in mind when reading 22 fixed and 7 broken. `FusedCost` does change cost values: it rounds in float32 at pixel edges. The next table shows what that did to each pair of runs. Reading the table: `orientation_differs` counts the voxels whose stored orientation differs at all. The `.mic` file writes the first Euler angle with 4 decimals and the other two with 6, so smaller differences are not seen. `median_deg` and `max_deg` are the angles of those voxels. `mic_lines_differ` counts voxel lines of the `.mic` file that differ in any column. `evals_differ` counts voxels whose evaluation count differs.

<!-- table:ar2b_fused_diff -->
| pair | orientation_differs | median_deg | above_0_001_deg | above_0_01_deg | max_deg | mic_lines_differ | evals_differ |
|---|---|---|---|---|---|---|---|
| key on: FusedCost on against off | 2 of 2000 | 1.9e-02 | 2 | 2 | 0.019 | 45 of 2000 | 431 of 2000 |
| key off: FusedCost on against off | 2000 of 2000 | 1.3e-02 | 2000 | 1303 | 54.511 | 2000 of 2000 | 829 of 2000 |
<!-- /table:ar2b_fused_diff -->

- With the key on, the stored orientations differ on only 2 of 2,000 voxels, both by about 0.019 degrees. The evaluation counts differ on 431 voxels, so the searches did run differently. The wrong voxels are the same 22 with `FusedCost` off and on (paired: 0 fixed, 0 broken, 22 wrong in both).
- How this was measured: a voxel counts as differing when its rotation matrices, read from the two `.mic` files, are not exactly equal. Its angle is the robust angle `ray_indexing.misorientation_deg` (atan2 of the sine and cosine parts). The review's 637 came from the `arccos`-based angle of `score_bfs`, which is above zero also for identical orientations (up to about 2.4e-6 degrees here). 635 of those 637 voxels have identical stored orientations.
- With the key off, every voxel's orientation differs (median 0.013 degrees, 1303 voxels above 0.01 degrees), and the evaluation counts differ on 829 voxels. The cause of this larger movement was not tested. It is consistent with the key-off run choosing different seeds once the costs change (5 and 1 rejected seeds instead of 6 and 1), which changes its route. This was not checked voxel by voxel.

## The flips, traced

The trace works on the stored run files only. It works the same way on both sides of a pair (N9). A broken voxel is traced in the key-on run, where it is wrong. A fixed voxel is traced in the key-off run, where it is wrong. For each such voxel the trace asks three things:
- Which true grain does its result match?
- Does that grain touch the voxel?
- Which full searches of the same run does its result match, within 1 degree?

What the method can show:
- whether a wrong result is the orientation of a touching grain;
- whether that orientation is also the result of a full search of the run that is right for its own voxel ("right seed behind"), or of a full search that is wrong ("wrong seed behind");
- whether the voxel was itself a seed in either run.

What it cannot show:
- How the orientation got to the voxel. `records.json` stores no parent, so "it came from the match through local fits" is inferred, not traced.
- Where the source lies. A match is any full search within 1 degree, and the nearest match is nearest in orientation, not in space. Here it lies 0.6 to 8.1 voxel sides away.
- Rejected seeds. A seed that failed the seed test keeps no result of its own in `recon.mic`, so it cannot be a match.
- A wrong indexed seed behind a break. With 0 of 64 indexed seeds wrong, that column is "no" by construction.

So the trace adds only one thing to "0 of 64 indexed seeds wrong": every flipped voxel takes the orientation of a touching grain.

Reading the table: one row per pair and flip. `traced_in` is the run where the voxel is wrong. `crossing` counts voxels whose result is the orientation of a touching grain. `right_seed_behind` counts voxels whose result matches a full search of that run that is right for its own voxel. `wrong_grid_seed_behind` and `wrong_indexed_seed_behind` count voxels whose result is, or matches, a wrong full search. `no_match` counts voxels with no match. `grid_seed_off` counts voxels that were themselves grid seeds or reseeds of the key-off run (in brackets, how many of them failed the seed test). `indexed_seed_on` counts voxels that were themselves indexed seeds or reseeds of the key-on run.

<!-- table:ar2b_flips -->
| pair | flip | traced_in | n | crossing | right_seed_behind | wrong_grid_seed_behind | wrong_indexed_seed_behind | no_match | grid_seed_off | indexed_seed_on |
|---|---|---|---|---|---|---|---|---|---|---|
| stored pair | fixed | key off | 22 | 22 | 18 | 4 | 0 | 0 | 3 (2 rejected) | 2 |
| stored pair | broken | key on | 7 | 7 | 7 | 0 | 0 | 0 | 0 (0 rejected) | 0 |
| FusedCost-on pair | fixed | key off | 14 | 14 | 14 | 0 | 0 | 0 | 1 (1 rejected) | 1 |
| FusedCost-on pair | broken | key on | 12 | 12 | 12 | 0 | 0 | 0 | 1 (0 rejected) | 0 |
<!-- /table:ar2b_flips -->

- Every fixed and every broken voxel, in both pairs, takes the orientation of a grain that touches it.
- On the key's own path (N10) the counts are: fixes, 2 in the stored pair (5876 and 24334) and 1 in the FusedCost-on pair (16476); breaks, 0 in both pairs.
  - Each of these 3 fixed voxels was a seed in both runs. In the key-off run the seed test rejected its grid result, and the voxel ended unresolved. In the key-on run its indexed result passed the seed test and is right.
  - These fixes count for the key, the mirror of N7 and N8. They are counts, not a rate.
- Every other fixed or broken voxel is a boundary crossing that counts neither way (N6, N9). The fixes are treated like the breaks. The whole-run wrong count is still not shown different (N9).
- No broken voxel has a wrong indexed seed behind it. None of the broken voxels is an indexed seed of the key-on run. The trace did not test the path itself; this follows from 0 of 64 indexed seeds being wrong.
- 4 of the 22 fixes of the stored pair have a wrong grid seed behind them in the key-off run:
  - Voxel 15749 was itself a wrong grid seed.
  - Voxels 5876, 8880 and 21188 match grid seed 8866, which is wrong.

  With the key on, neither 15749 nor 8866 was a seed: both were revisits. So, apart from 5876, these fixes come with the changed route, not with indexing those voxels.
- The route change is an effect of the key. It shows in the whole-run counts and in the paired splits of both pairs. One question was not tested: would a key-off run with the key-on run's seeds break or fix the same voxels?

Reading the table: every voxel that is right with the key off and wrong with the key on, against the stored rerun, traced in the key-on run.
- `kind` says whether another grain shares an edge or only a corner.
- `source_on` and `source_off` are how each run got the voxel: a local fit from a neighbour, a revisit, its own seed, or unresolved.
- `taken` is the grain whose true orientation the result matches.
- `crossing` says whether that grain touches the voxel.
- `match` is what the voxel's result matches among the full searches of the key-on run.
- `wrong_seed_behind` is "indexed" or "grid" when the voxel's result is, or matches, a wrong full search.
- `seed_off` says whether the voxel was a grid seed of the key-off run, accepted or rejected by the seed test. `seed_on` says whether it was an indexed seed of the key-on run.

<!-- table:ar2b_breaks -->
| voxel | kind | source_on | source_off | err_on | taken | crossing | match | wrong_seed_behind | seed_off | seed_on |
|---|---|---|---|---|---|---|---|---|---|---|
| 7194 | edge | neighbor | revisit | 51.86 | grain 154 | yes | matches a right full search | no | no | no |
| 8890 | edge | neighbor | neighbor | 42.27 | grain 148 | yes | matches a right full search | no | no | no |
| 14299 | edge | neighbor | revisit | 44.52 | grain 24 | yes | matches a right full search | no | no | no |
| 16508 | edge | revisit | revisit | 33.49 | grain 400 | yes | matches a right full search | no | no | no |
| 16534 | edge | neighbor | revisit | 29.65 | grain 24 | yes | matches a right full search | no | no | no |
| 17286 | edge | neighbor | revisit | 31.44 | grain 247 | yes | matches a right full search | no | no | no |
| 21180 | edge | neighbor | neighbor | 33.73 | grain 396 | yes | matches a right full search | no | no | no |
<!-- /table:ar2b_breaks -->

Reading the table: the same 7 breaks in full. `result_cost` is the strict cost of the orientation the voxel ended with. `cost_at_truth` is the strict cost at the voxel's exact true orientation. `match` is the full-sample index of the nearest match (here an indexed seed or reseed). `match_err` is that seed's error in degrees. `match_dist` is its distance from the voxel in voxel sides. `match_rank` is the candidate it took (1 is the first one the indexer proposed). `match_cost` is its strict cost before and after the finisher.

<!-- table:ar2b_break_trace -->
| voxel | source | err_deg | result_cost | cost_at_truth | match | match_err | match_dist | match_rank | match_cost |
|---|---|---|---|---|---|---|---|---|---|
| 7194 | neighbor | 51.9 | 0.123 | 0.046 | 983 | 0.005 | 1.7 | 1 | 0.069 to 0.054 |
| 8890 | neighbor | 42.3 | 0.157 | 0.052 | 5876 | 0.023 | 1.7 | 1 | 0.123 to 0.078 |
| 14299 | neighbor | 44.5 | 0.137 | 0.053 | 9015 | 0.017 | 1.5 | 1 | 0.113 to 0.056 |
| 16508 | revisit | 33.5 | 0.192 | 0.057 | 8948 | 0.009 | 2.5 | 2 | 0.151 to 0.119 |
| 16534 | neighbor | 29.6 | 0.135 | 0.045 | 9015 | 0.017 | 2.5 | 1 | 0.113 to 0.056 |
| 17286 | neighbor | 31.4 | 0.121 | 0.024 | 6684 | 0.015 | 3.1 | 1 | 0.046 to 0.046 |
| 21180 | neighbor | 33.7 | 0.143 | 0.052 | 21170 | 0.004 | 2.6 | 1 | 0.086 to 0.066 |
<!-- /table:ar2b_break_trace -->

- All 7 breaks are edge-sharing boundary voxels that take the orientation of a grain that touches them. Each matches an indexed seed or reseed that is right for its own voxel (error 0.023 degrees or less), 1.5 to 3.1 voxel sides away. The nearest match took the first candidate in 6 cases and the second in 1.
- In all 7 the strict cost at the exact truth (0.024 to 0.057) is lower than the cost of the result (0.121 to 0.192). So the local fit stayed at the touching grain's orientation, although the cost prefers the truth.
- None is traced to a wrong indexing result. The route change was not tested.

Reading the table: every voxel that is wrong with the key off and right with the key on, against the stored rerun, traced in the key-off run. The columns are those of the breaks table, with the sides swapped. `source_off` and `err_off` belong to the run where the voxel is wrong. `match` is what its result matches among the key-off run's own full searches (grid seeds and reseeds).

<!-- table:ar2b_fixes -->
| voxel | kind | source_off | source_on | err_off | taken | crossing | match | wrong_seed_behind | seed_off | seed_on |
|---|---|---|---|---|---|---|---|---|---|---|
| 1384 | edge | neighbor | revisit | 20.54 | grain 75 | yes | matches a right full search | no | no | no |
| 1398 | edge | neighbor | revisit | 49.49 | grain 23 | yes | matches a right full search | no | no | no |
| 1981 | edge | neighbor | revisit | 34.20 | grain 23 | yes | matches a right full search | no | no | no |
| 3204 | edge | revisit | revisit | 34.20 | grain 23 | yes | matches a right full search | no | no | no |
| 5008 | edge | revisit | neighbor | 31.45 | grain 399 | yes | matches a right full search | no | no | no |
| 5876 | edge | unresolved | seed | 42.25 | grain 398 | yes | matches a wrong grid seed | grid | rejected | indexed |
| 8880 | edge | neighbor | neighbor | 33.73 | grain 398 | yes | matches a wrong grid seed | grid | no | no |
| 8899 | edge | revisit | revisit | 34.20 | grain 23 | yes | matches a right full search | no | no | no |
| 8909 | edge | neighbor | neighbor | 44.11 | grain 332 | yes | matches a right full search | no | no | no |
| 8940 | edge | neighbor | neighbor | 31.45 | grain 399 | yes | matches a right full search | no | no | no |
| 8978 | edge | revisit | revisit | 57.34 | grain 401 | yes | matches a right full search | no | no | no |
| 8988 | edge | neighbor | revisit | 52.01 | grain 153 | yes | matches a right full search | no | no | no |
| 13672 | edge | neighbor | revisit | 36.69 | grain 75 | yes | matches a right full search | no | no | no |
| 15749 | edge | seed | revisit | 45.75 | grain 331 | yes | own grid seed (wrong) | grid | accepted | no |
| 15907 | edge | neighbor | revisit | 46.02 | grain 75 | yes | matches a right full search | no | no | no |
| 16470 | edge | neighbor | revisit | 45.18 | grain 145 | yes | matches a right full search | no | no | no |
| 18179 | edge | neighbor | neighbor | 29.78 | grain 300 | yes | matches a right full search | no | no | no |
| 19418 | edge | neighbor | revisit | 36.69 | grain 75 | yes | matches a right full search | no | no | no |
| 19436 | edge | neighbor | revisit | 53.33 | grain 279 | yes | matches a right full search | no | no | no |
| 21188 | edge | neighbor | revisit | 47.14 | grain 398 | yes | matches a wrong grid seed | grid | no | no |
| 21215 | edge | neighbor | neighbor | 54.51 | grain 332 | yes | matches a right full search | no | no | no |
| 24334 | edge | unresolved | seed | 50.85 | grain 56 | yes | matches a right full search | no | rejected | indexed |
<!-- /table:ar2b_fixes -->

Reading the table: the same 22 fixes in full, in the key-off run. The columns are those of the break trace. `match` is "itself" when the voxel's result is its own grid seed. `match_rank` and `match_cost` are "n/a" because a grid search has no candidate list.

<!-- table:ar2b_fix_trace -->
| voxel | source | err_deg | result_cost | cost_at_truth | match | match_err | match_dist | match_rank | match_cost |
|---|---|---|---|---|---|---|---|---|---|
| 1384 | neighbor | 20.5 | 0.132 | 0.071 | 5872 | 0.001 | 5.3 | n/a | n/a |
| 1398 | neighbor | 49.5 | 0.124 | 0.060 | 19440 | 0.001 | 2.5 | n/a | n/a |
| 1981 | neighbor | 34.2 | 0.133 | 0.074 | 19440 | 0.001 | 2.5 | n/a | n/a |
| 3204 | revisit | 34.2 | 0.125 | 0.074 | 19440 | 0.001 | 5.8 | n/a | n/a |
| 5008 | revisit | 31.5 | 0.110 | 0.043 | 21243 | 0.005 | 2.1 | n/a | n/a |
| 5876 | unresolved | 42.3 | 0.342 | 0.054 | 8866 | 50.850 | 2.6 | n/a | n/a |
| 8880 | neighbor | 33.7 | 0.124 | 0.066 | 8866 | 50.850 | 3.0 | n/a | n/a |
| 8899 | revisit | 34.2 | 0.120 | 0.074 | 19440 | 0.001 | 4.9 | n/a | n/a |
| 8909 | neighbor | 44.1 | 0.127 | 0.048 | 21210 | 0.015 | 2.3 | n/a | n/a |
| 8940 | neighbor | 31.5 | 0.109 | 0.043 | 21243 | 0.005 | 1.2 | n/a | n/a |
| 8978 | revisit | 57.3 | 0.139 | 0.057 | 14283 | 0.006 | 2.1 | n/a | n/a |
| 8988 | neighbor | 52.0 | 0.152 | 0.071 | 21289 | 0.004 | 0.6 | n/a | n/a |
| 13672 | neighbor | 36.7 | 0.134 | 0.060 | 5872 | 0.001 | 4.0 | n/a | n/a |
| 15749 | seed | 45.8 | 0.150 | 0.011 | itself | 45.753 | 0.0 | n/a | n/a |
| 15907 | neighbor | 46.0 | 0.137 | 0.049 | 5872 | 0.001 | 8.1 | n/a | n/a |
| 16470 | neighbor | 45.2 | 0.157 | 0.055 | 24281 | 0.018 | 1.0 | n/a | n/a |
| 18179 | neighbor | 29.8 | 0.113 | 0.045 | 8932 | 0.009 | 4.0 | n/a | n/a |
| 19418 | neighbor | 36.7 | 0.133 | 0.060 | 5872 | 0.001 | 3.1 | n/a | n/a |
| 19436 | neighbor | 53.3 | 0.128 | 0.051 | 777 | 0.006 | 3.2 | n/a | n/a |
| 21188 | neighbor | 47.1 | 0.134 | 0.074 | 8866 | 50.850 | 4.0 | n/a | n/a |
| 21215 | neighbor | 54.5 | 0.130 | 0.044 | 21210 | 0.015 | 1.0 | n/a | n/a |
| 24334 | unresolved | 50.8 | 0.423 | 0.048 | 7944 | 0.009 | 1.5 | n/a | n/a |
<!-- /table:ar2b_fix_trace -->

- All 22 fixes are edge-sharing boundary voxels that take the orientation of a touching grain in the key-off run.
  - 18 match a grid seed or reseed that is right for its own voxel (error 0.018 degrees or less).
  - 3 match the wrong grid seed 8866 (error 50.850 degrees).
  - 1 (15749) is itself a wrong grid seed.
- In all 22 the strict cost at the exact truth (0.011 to 0.074) is lower than the cost of the key-off result (0.109 to 0.423). The highest result costs belong to the 2 unresolved voxels (5876 and 24334). Their stored result is the fit that failed.

Reading the table: the breaks of the FusedCost-on pair (key on against key off, both with FusedCost on), with the columns of the first breaks table.

<!-- table:ar2b_breaks_fused -->
| voxel | kind | source_on | source_off | err_on | taken | crossing | match | wrong_seed_behind | seed_off | seed_on |
|---|---|---|---|---|---|---|---|---|---|---|
| 1402 | edge | neighbor | revisit | 46.00 | grain 350 | yes | matches a right full search | no | no | no |
| 7179 | edge | neighbor | revisit | 30.74 | grain 350 | yes | matches a right full search | no | no | no |
| 7194 | edge | neighbor | revisit | 51.86 | grain 154 | yes | matches a right full search | no | no | no |
| 8866 | edge | revisit | revisit | 50.86 | grain 398 | yes | matches a right full search | no | no | no |
| 8890 | edge | neighbor | revisit | 42.27 | grain 148 | yes | matches a right full search | no | no | no |
| 14299 | edge | neighbor | revisit | 44.52 | grain 24 | yes | matches a right full search | no | no | no |
| 16508 | edge | revisit | seed | 33.49 | grain 400 | yes | matches a right full search | no | accepted | no |
| 16534 | edge | neighbor | revisit | 29.65 | grain 24 | yes | matches a right full search | no | no | no |
| 17286 | edge | neighbor | revisit | 31.44 | grain 247 | yes | matches a right full search | no | no | no |
| 21180 | edge | neighbor | revisit | 33.73 | grain 396 | yes | matches a right full search | no | no | no |
| 21234 | edge | revisit | revisit | 46.81 | grain 300 | yes | matches a right full search | no | no | no |
| 21328 | edge | revisit | revisit | 47.23 | grain 151 | yes | matches a right full search | no | no | no |
<!-- /table:ar2b_breaks_fused -->

Reading the table: the same 12 breaks of the FusedCost-on pair in full, with the columns of the first break trace.

<!-- table:ar2b_break_trace_fused -->
| voxel | source | err_deg | result_cost | cost_at_truth | match | match_err | match_dist | match_rank | match_cost |
|---|---|---|---|---|---|---|---|---|---|
| 1402 | neighbor | 46.0 | 0.111 | 0.051 | 19480 | 0.003 | 4.7 | 1 | 0.088 to 0.029 |
| 7179 | neighbor | 30.7 | 0.115 | 0.051 | 19480 | 0.003 | 1.2 | 1 | 0.088 to 0.029 |
| 7194 | neighbor | 51.9 | 0.123 | 0.046 | 983 | 0.005 | 1.7 | 1 | 0.069 to 0.054 |
| 8866 | revisit | 50.9 | 0.134 | 0.053 | 8889 | 0.012 | 3.6 | 1 | 0.063 to 0.063 |
| 8890 | neighbor | 42.3 | 0.157 | 0.052 | 5876 | 0.023 | 1.7 | 1 | 0.123 to 0.078 |
| 14299 | neighbor | 44.5 | 0.137 | 0.053 | 9015 | 0.017 | 1.5 | 1 | 0.113 to 0.056 |
| 16508 | revisit | 33.5 | 0.192 | 0.057 | 8948 | 0.009 | 2.5 | 2 | 0.151 to 0.119 |
| 16534 | neighbor | 29.6 | 0.135 | 0.045 | 9015 | 0.017 | 2.5 | 1 | 0.113 to 0.056 |
| 17286 | neighbor | 31.4 | 0.121 | 0.024 | 6684 | 0.015 | 3.1 | 1 | 0.046 to 0.046 |
| 21180 | neighbor | 33.7 | 0.143 | 0.052 | 21170 | 0.004 | 2.6 | 1 | 0.086 to 0.066 |
| 21234 | revisit | 46.8 | 0.118 | 0.057 | 3205 | 0.003 | 7.0 | 2 | 0.181 to 0.036 |
| 21328 | revisit | 47.2 | 0.131 | 0.045 | 9027 | 0.021 | 1.2 | 1 | 0.111 to 0.053 |
<!-- /table:ar2b_break_trace_fused -->

The 12 breaks of this pair are also all edge-sharing boundary crossings. Each matches an indexed seed or reseed that is right for its own voxel (error 0.023 degrees or less), 1.2 to 7.0 voxel sides away. None is traced to a wrong indexing result; the route change was not tested. All 7 breaks of the first pair are in this list too. Voxel 16508 was a right grid seed of the key-off run and a revisit of the key-on run, so it shows the changed route.

Reading the table: the fixes of the FusedCost-on pair, traced in the FusedCost-on key-off run, with the columns of the first fixes table.

<!-- table:ar2b_fixes_fused -->
| voxel | kind | source_off | source_on | err_off | taken | crossing | match | wrong_seed_behind | seed_off | seed_on |
|---|---|---|---|---|---|---|---|---|---|---|
| 1384 | edge | neighbor | revisit | 20.54 | grain 75 | yes | matches a right full search | no | no | no |
| 1964 | edge | neighbor | neighbor | 43.79 | grain 395 | yes | matches a right full search | no | no | no |
| 3213 | edge | neighbor | revisit | 54.04 | grain 24 | yes | matches a right full search | no | no | no |
| 8879 | edge | revisit | revisit | 47.13 | grain 147 | yes | matches a right full search | no | no | no |
| 8978 | edge | revisit | revisit | 57.34 | grain 401 | yes | matches a right full search | no | no | no |
| 8988 | edge | neighbor | revisit | 52.01 | grain 153 | yes | matches a right full search | no | no | no |
| 13250 | edge | neighbor | revisit | 46.25 | grain 148 | yes | matches a right full search | no | no | no |
| 13672 | edge | neighbor | revisit | 36.69 | grain 75 | yes | matches a right full search | no | no | no |
| 15907 | edge | neighbor | revisit | 46.02 | grain 75 | yes | matches a right full search | no | no | no |
| 16476 | edge | unresolved | seed | 49.08 | grain 144 | yes | matches a right full search | no | rejected | indexed |
| 16543 | edge | neighbor | revisit | 22.10 | grain 155 | yes | matches a right full search | no | no | no |
| 19418 | edge | neighbor | revisit | 36.69 | grain 75 | yes | matches a right full search | no | no | no |
| 19436 | edge | neighbor | revisit | 53.33 | grain 279 | yes | matches a right full search | no | no | no |
| 21164 | edge | revisit | revisit | 50.83 | grain 397 | yes | matches a right full search | no | no | no |
<!-- /table:ar2b_fixes_fused -->

Reading the table: the same 14 fixes in full, with the columns of the first fix trace.

<!-- table:ar2b_fix_trace_fused -->
| voxel | source | err_deg | result_cost | cost_at_truth | match | match_err | match_dist | match_rank | match_cost |
|---|---|---|---|---|---|---|---|---|---|
| 1384 | neighbor | 20.5 | 0.131 | 0.071 | 5872 | 0.003 | 5.3 | n/a | n/a |
| 1964 | neighbor | 43.8 | 0.124 | 0.028 | 21138 | 0.017 | 2.1 | n/a | n/a |
| 3213 | neighbor | 54.0 | 0.137 | 0.036 | 9015 | 0.022 | 4.6 | n/a | n/a |
| 8879 | revisit | 47.1 | 0.156 | 0.052 | 8899 | 0.004 | 5.3 | n/a | n/a |
| 8978 | revisit | 57.3 | 0.142 | 0.057 | 14283 | 0.012 | 2.1 | n/a | n/a |
| 8988 | neighbor | 52.0 | 0.152 | 0.071 | 21289 | 0.006 | 0.6 | n/a | n/a |
| 13250 | neighbor | 46.3 | 0.137 | 0.036 | 5876 | 0.006 | 2.5 | n/a | n/a |
| 13672 | neighbor | 36.7 | 0.136 | 0.060 | 5872 | 0.003 | 4.0 | n/a | n/a |
| 15907 | neighbor | 46.0 | 0.139 | 0.049 | 5872 | 0.003 | 8.1 | n/a | n/a |
| 16476 | unresolved | 49.1 | 0.387 | 0.028 | 399 | 0.012 | 3.1 | n/a | n/a |
| 16543 | neighbor | 22.1 | 0.115 | 0.086 | 16539 | 0.003 | 2.0 | n/a | n/a |
| 19418 | neighbor | 36.7 | 0.136 | 0.060 | 5872 | 0.003 | 3.1 | n/a | n/a |
| 19436 | neighbor | 53.3 | 0.131 | 0.051 | 777 | 0.013 | 3.2 | n/a | n/a |
| 21164 | revisit | 50.8 | 0.139 | 0.052 | 21153 | 0.029 | 3.0 | n/a | n/a |
<!-- /table:ar2b_fix_trace_fused -->

All 14 fixes of this pair are edge-sharing voxels that take a touching grain's orientation and match a right grid seed or reseed of the key-off run (error 0.029 degrees or less). Voxel 16476 was a rejected grid seed in the key-off run and a right indexed seed in the key-on run. It is the fix on the key's own path in this pair (N10). The other 13 are boundary crossings.

Reading the table: where the wrong voxels lie (edge-sharing, corner-only or interior), for each run, and the fixed and broken voxels of the stored pair.

<!-- table:ar2b_kinds -->
| run | edge | corner | interior |
|---|---|---|---|
| key off (stored rerun) | 37 | 0 | 0 |
| key on | 22 | 0 | 0 |
<!-- /table:ar2b_kinds -->

<!-- table:ar2b_flip_kinds -->
| flip | edge | corner | interior |
|---|---|---|---|
| fixed | 22 | 0 | 0 |
| broken | 7 | 0 | 0 |
<!-- /table:ar2b_flip_kinds -->

Every wrong voxel, in both runs, and every fixed and broken voxel is edge-sharing. No corner-only or interior voxel is wrong. This matches the holdout of step 1.

## Evaluations

Reading the table: cost evaluations of the whole pilot. `full` is the seed and reseed bucket (grid searches for key off; the indexing path, and any grid fallback, for key on). `local` is the neighbour and revisit fits. `full_searches` counts seeds plus reseeds. `local_fits` counts the local fits. The rows are the stored key-off rerun, key on, and the two FusedCost-on runs.

<!-- table:ar2b_evals -->
| run | total | full | local | full_searches | local_fits |
|---|---|---|---|---|---|
| key off (stored rerun) | 25,547,935 | 22,882,865 | 2,665,070 | 69 | 4049 |
| key on | 2,819,906 | 122,198 | 2,697,708 | 64 | 4095 |
| key off, FusedCost on | 25,602,383 | 22,904,126 | 2,698,257 | 69 | 4092 |
| key on, FusedCost on | 2,816,897 | 122,154 | 2,694,743 | 64 | 4095 |
<!-- /table:ar2b_evals -->

The stored baseline: 22,882,865 evaluations in 69 full searches and 2,665,070 in 4,049 local fits, 25,547,935 in all. These are the numbers of the task. They were recomputed from the stored files (no rerun) and match. A grid search costs 331,636 evaluations on average. An indexing path costs 1,909.

The number of full searches is 64 with the key on and 69 with it off.
- The key-off run had 40 seeds and 29 reseeds. Six of its seeds and one of its reseeds failed the seed test, so it accepted 34 seeds and 28 reseeds.
- The key-on run had 33 seeds and 31 reseeds, and accepted all of them.
- So the accepted seeds are about the same (34 against 33), and the 6 rejected seeds account for most of the 40 against 33.
- The key-on run has 2 more reseeds (31 against 29) and 3 more accepted reseeds (31 against 28). Why the reseed count rises was not tested. It is part of the route change.
- The FusedCost-on key-off run had 37 seeds and 32 reseeds (5 and 1 failed the seed test), so the counts also move with the other setting.

Reading the table: how every seed and reseed of the key-on run was handled, how often seeds failed the seed test, and how often the results of the grid search and of the indexing are wrong.

<!-- table:ar2b_index -->
| item | value |
|---|---|
| seeds and reseeds | 64 |
| indexed (every seed and reseed) | 64 |
| used (passed the seed test) | 64 |
| fallbacks to the grid search | 0 of 64 (0.0 to 5.7%) |
|   no candidate | 0 |
|   failed the seed test | 0 |
| evaluations of used indexing paths | 122,198 |
| evaluations of indexing paths that fell back | 0 |
| evaluations of the grid fallbacks | 0 |
| key off: grid seeds and reseeds that failed the seed test | 7 of 69 (5.0 to 19.5%) |
| key on: indexed seeds and reseeds that failed the seed test | 0 of 64 (0.0 to 5.7%) |
| key off: voxels whose result is their own grid search, wrong | 2 of 62 (0.9 to 11.0%) |
| key on: voxels whose result is their own indexed seed, wrong | 0 of 64 (0.0 to 5.7%) |
<!-- /table:ar2b_index -->

## The indexing path

The cost of the key's path, beside its result:
- Evaluations: 1,909 per seed, instead of 331,636 for a grid search.
- Wrong results: 0 among 64 indexed seeds against 2 among 62 grid seeds. The intervals overlap, so this is not shown different.
- Seed-test failures: 0 of 64 against 7 of 69 for the grid search (0.0 to 5.7% against 5.0 to 19.5%). The intervals touch, so this is not shown different either.
- Time outside the evaluation count: the path adds the one-off pairing and the indexing time of each voxel (next section).

Best 2 against best 1, as seen in this run: the second finisher run changed the candidate taken in 2 of 64 seeds. In 9 of 64 the taken result's finisher cost equals its strict cost before the finisher, so the finisher did not move it. This bears on voxel 17361 below.

## Times outside the evaluation count

Reading the table: the cost of the rays and of indexing, in seconds, from the key-on run with FusedCost off. Contended. The pairing happens inside the first seed. It is not an evaluation. The indexing time is the voting, clustering and refining of one voxel, with the rays already attached. "Whole path" adds the strict costs and the finisher runs.

<!-- table:ar2b_times -->
| item | value |
|---|---|
| one-time pairing of the rays (s) | 17.6 |
| indexing, total (s) | 11.9 |
| indexing per voxel, mean (s) | 0.186 |
| indexing per voxel, median (s) | 0.181 |
| indexing per voxel, max (s) | 0.333 |
| whole path per voxel, mean (s) | 1.01 |
<!-- /table:ar2b_times -->

Pairing the rays took 17.1 s in the FusedCost-on run. Both times are one-off costs of about 17 s. Indexing is 11.9 s of the 1272 s of the key-on run with FusedCost off.

Reading the table: wall time of whole runs, with the conditions. All of this task's runs are contended: the preflight refused, and an editor process used about 200% of one core. The key-on run with FusedCost off also overlapped another run. The stored key-off rerun ran under different conditions (three BFS runs at once on 2026-10-29), so its time cannot be compared with the others.

<!-- table:ar2b_wall -->
| run | wall_s | conditions |
|---|---|---|
| key off, stored rerun (FusedCost off) | 8997 | contended: 3 BFS runs at once, 1 thread each, 2026-10-29 |
| key on (FusedCost off) | 1272 | contended: timing preflight refused (an editor process at about 200% CPU); not a single run: overlapped with the voxel 17361 run (its log spans 23 s) and with short checks whose times were not logged |
| key on (FusedCost on) | 213 | contended: timing preflight refused (an editor process at about 200% CPU), single run |
| key off (FusedCost on) | 1633 | contended: timing preflight refused (an editor process at about 200% CPU), single run |
<!-- /table:ar2b_wall -->

The two FusedCost-on runs ran one after the other on the same machine, with the same editor process present in both (the load was not logged during the runs). They give 1633 s with the key off and 213 s with it on, about 7.7 times shorter. This is a first contended reading, one run per arm, taken after the preflight refused (N9). Review item 14 asks for this timing on a quiet machine, and it stays open. With the key on, the seeds and reseeds take about 36 s of the 213 s (17 s of it the pairing), and the local fits about 177 s. So most of the remaining wall time is the local fits. Their evaluation count is about the same as with the key off (2,694,743 against 2,698,257).

## Voxel 17361

This voxel (full-sample index, outside the pilot region) was wrong in all three arms of step 1. The library path was run on it. Reading the table: the truth candidate is the first candidate the indexer proposes. "Strict cost n deg from the truth" is the cost after turning the true orientation by n degrees about one axis. The ray Gauss-Newton step is `RayGaussNewton` (Task F's F4, one frame window, strict-cost acceptance 0.8).

<!-- table:ar2b_17361 -->
| item | value |
|---|---|
| truth candidate: rank in the indexer's list | 1 |
| truth candidate: distance from the truth (deg) | 0.0449 |
| truth candidate: strict cost | 0.1623 |
| strict cost at the exact truth | 0.0949 |
| strict cost 0.02 deg from the truth | 0.1556 |
| strict cost 0.045 deg from the truth | 0.2367 |
| strict cost 0.1 deg from the truth | 0.4356 |
| finisher, total evaluations (truth candidate) | 924 |
| stage 1 (FindOptimal, MC blocks): evaluations | 97 |
| stage 1: blocks run, blocks accepted | 3, 0 |
| stage 1: first MC step and box (deg) | 0.132, 0.329 |
| stage 1: cost before, after | 0.1623, 0.1623 |
| stage 2 (VarianceMinimizing): evaluations | 826 |
| stage 2: cost before, after | 0.1623, 0.1623 |
| restart of the finisher from its end point: cost, moved (deg) | 0.1623, 0.0000 |
| ray Gauss-Newton from the end point: cost, accepted | 0.1283, True |
| ray Gauss-Newton: moved (deg), distance from the truth after (deg) | 0.0627, 0.0190 |
| ray Gauss-Newton: strict, smooth evaluations | 2, 14 |
<!-- /table:ar2b_17361 -->

- The library path reproduces the step-1 numbers: the truth candidate starts 0.0449 degrees from the truth with strict cost 0.1623, the finisher uses 924 evaluations, and its end point equals its start.
- Both finisher stages ran. Stage 1 (FindOptimal, a Monte Carlo search in blocks of 31 steps) ran 3 blocks, accepted none, and stopped on its limit of successive failed blocks. Stage 2 (VarianceMinimizing) ran 826 evaluations and did not lower the cost.
- Restarting the finisher from its own end point gives the same end point. That restart used the same random stream (`default_rng(0)`), as the library does, so it repeats the same draws from the same place. By construction it cannot show anything new. No other random stream was tried.
- One `RayGaussNewton` step from the end point reaches strict cost 0.1283, accepted, 0.019 degrees from the truth. That is lower than 0.1623 and still above the 0.0949 of the exact truth. So one call does not reach a cost near 0.095.
- The strict cost rises from 0.0949 at the exact truth to 0.1556 at 0.02 degrees and 0.2367 at 0.045 degrees. The first Monte Carlo step is 0.132 degrees. This is consistent with the cost having a minimum narrower than the finisher's first step. It was not tested with a smaller step, and the probe turns about one axis only. The arm was not changed because of this.
- With the key on, this voxel would take the touching grain's orientation (strict cost 0.1329, 57.3 degrees from the truth), as in step 1, if the seed test passes it. The seed test was not run on this voxel.

## How the indexing here differs from AR1's

AR1 and step 1 paired the rays with the hull of the whole sample. The BFS pilot's sample is the 2,000-voxel region, so its hull is smaller and the paired rays differ.
- 47 voxels are both an indexed seed here and one of AR1's 69 full-search voxels. 9 of them have the same number of attached rays, and 8 of those 9 give the identical result (strict costs, finisher costs, evaluations). 24 of the other 38 also do (`summary.json`, `ar1_comparison`). The candidates in the tail of the list differ in some voxels.
- Indexing eight pilot voxels again in AR1's setting reproduces AR1's candidates and strict costs exactly (`thread_check_default.json`, 8 of 8).

So the difference is consistent with the setting (the hull), not with the code path. This was inferred from the code (`SmoothCost.observed_rays` uses `sample_hull_points`). It was not tested by swapping the hull. The library's context build was not compared with the script's (see "What was built").

## Criteria of the task

- Off path bit-identical: shown with fakes and the golden grid. The real-image key-off run was not compared with an older commit.
- Stored baseline: its numbers were recomputed from the stored run files (no rerun) and match (37, 25,547,935, 69).
- Wrong voxels, paired: 22 fixed, 7 broken, 15 wrong in both (stored pair); 14 fixed, 12 broken, 10 wrong in both (FusedCost-on pair).
  - On the key's own path (N10): fixes 5876 and 24334 (stored pair) and 16476 (FusedCost-on pair), which count for the key; 0 breaks. None traced to a wrong indexing result.
  - Every other fixed or broken voxel is a boundary crossing (all are edge-sharing) that counts neither way. The trace treats both sides the same way (N9). The route change was not tested.
  - Accuracy: 0 of 64 indexed seeds wrong against 2 of 62 grid seeds, not shown different. The whole-run wrong count is not shown different (N9).
- Evaluations: 2,819,906 against 25,547,935.
- Fallbacks: 0 of 64 (0.0 to 5.7%). The fallback path has no real-image count.
- Wall time: 213 s against 1633 s (FusedCost on), a first contended reading after a refused preflight. Review item 14's quiet timing stays open (N9). The key-on run with FusedCost off took 1272 s and overlapped another run.
- Voxel 17361: diagnosed; the arm is unchanged.

## What was not tested

- The fallback path on real images. No seed or reseed fell back, so its evaluation cost on this sample is not measured.
- Realistic images (neighbour overlap and detector noise). Clean images only.
- More than one sample: a copper sample, one layer, one 2,000-voxel region, one BFS random seed, one run per arm. The key-off movement under `FusedCost` shows that one run is not a stable baseline.
- The whole layer. In a larger sample the hull is larger and the pairing differs from the pilot's.
- The route change on its own: whether a key-off run with the key-on run's seeds would break or fix the same voxels. The trace shows only orientation matches; how an orientation reached a voxel is inferred.
- The cause of the key-off movement under `FusedCost`, and of the change in the number of seeds and reseeds.
- A real-image key-off run against an older commit (bit-identity on real images). The off path is tested with fakes, the golden grid and the unchanged existing tests.
- The library's context build (hull, crystal symmetry) against the scripts'.
- Quiet timing (review item 14). The preflight refused, so every time is contended. The load during each run was not logged. The FusedCost-on pair ran in sequence, not interleaved.
- Voxels with few attached rays. Step 1 had none below 30 explained rays.
- Why the finisher misses the narrow minimum of voxel 17361 (see above for what was observed).
- Any change of the finisher, the candidate counts (5 and 2), or the seed test. None was tried.
- C++ comparison. The key is Python only.
