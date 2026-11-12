# AR5 step 5: two opt-in keys for run S11

Task `ar5e`, 2026-11-11 (plan: `PLAN.md`, "AR5 step 5"; decision D-20261111-N21).

## Summary

Run `S11` used to exist only in scripts. It is `S8` with two changes. Both are now keys of the
library, both off by default. This page shows that a library run gives the same result as the
stored script run, on two samples.

"Same result" (parity) means: the same evaluation counts, the same seed steps and reseeds, the
same per-voxel records and the same map file, byte for byte. All four checks of the BFS part
were exact. The boundary pass was exact when started from the stored maps.

Inside a library run, the pass for `S11` made as many tries and evaluations as the stored pass, on
both samples. But the tries were not the same. Costs differ in 4,372 of 18,617 tries on `big` and
in 15,813 of 60,423 on the first gradient sample. The pass also changed a slightly different set
of voxels. For the other two runs on `big`, even the number of tries differed: by 1 and 4 tries
out of about 18,600.

No setting is recommended. No timing was made. See "What was not tested".

## Terms

- **`big`**: the sample with 32 grains. **First gradient sample**: the sample with 497 grains in
  which the orientation turns linearly with position inside each grain (AR3 step 3). Both have
  24,570 voxels and clean images, reconstructed at Q-max 8.
- **Evaluation**: one call of the cost function. All counts here are evaluations.
- **Hit ratio**: for one voxel, the share of an orientation's predicted detector pixels that fall
  on observed intensity (0 to 1). Computing it takes one overlap evaluation.
- **CMA fit** and **MC optimizer**: two local optimizations of one voxel's orientation. CMA is an
  evolution strategy (CMA-ES); MC is the stock Monte Carlo optimizer of the library.
- **Seed step**: ray indexing proposes candidate orientations for a new seed voxel. The first
  ones get the strict cost, and the cheapest are refined by a local optimization. **Reseed**: a
  new seed placed in a leftover cluster of voxels the BFS rejected.
- **Stop ratio**: the hit ratio set by `SeedIndexStopHitRatio` (0.95 here). The seed step stops
  at the first refined candidate whose end reaches it, and keeps that end.
- **Grid fallback**: the full grid search that runs when the seed step gives no candidate, or when
  its result fails the seed test (a hit ratio of at least 0.80).
- **Boundary pass**: a step after the BFS. For each voxel with a neighbour more than 0.5 degree
  away, it makes one try from each distinct neighbour orientation. The voxel keeps the best end:
  the highest hit ratio, then the lowest cost. **Try**: one CMA fit of one voxel, started from one
  neighbour's orientation (the try's candidate).
- **Mixed seed refinement**: for each indexed seed candidate, a CMA fit first. Its end is kept when
  its hit ratio reaches the stop ratio. If not, the stock MC optimizer runs from the same
  candidate, and the end with the strictly lower cost is kept.
- **No MC fallback**: in a neighbour fit or revisit that starts with a CMA fit, the CMA end is
  returned. The MC optimizer does not run after it, even if the CMA end did not lower the cost.
- **Wrong voxel**: an error above 1 degree against the truth. The library never sees the truth;
  only the check script does.

## The two keys

```
SeedByIndexing 1
BFSNeighborCMAFirst 1
CMASigma0 0.05
CMANeighborMaxEvals 60
SeedIndexCandidates 20 5
SeedIndexStopHitRatio 0.95
SeedIndexCMAFirst 200          # new: mixed seed refinement, at most 200 CMA evaluations
BFSNeighborCMANoFallback 1     # new: return the CMA end of a neighbour fit
```

- `SeedIndexCMAFirst n`: default 0 (off). A value of 2 or more switches the mixed seed refinement
  on. The first step of the fit is `CMASigma0`; the stop ratio is `SeedIndexStopHitRatio`. A CMA
  fit of 200 counts 201 evaluations. The hit ratio of the CMA end is one more overlap evaluation
  per refined candidate, and it is not counted: 436 candidates on `big` and 2,158 on the first
  gradient sample. When the CMA end is kept, the stop check of the seed step computes the same hit
  ratio again. The script runs counted it the same way, so parity is not affected. Refused unless `SeedByIndexing 1` is set and the stop ratio
  is above 0. Refused with `PaintGrains 1`, `LocalOptimizer cma` or the hybrid optimizer. A value
  of 1 is refused.
- `BFSNeighborCMANoFallback 0|1`: default 0. Refused unless `BFSNeighborCMAFirst 1` is set.

Counters that exist only when a key is on: `cma_not_fell` (CMA fits whose end cost was not below
the start cost), `index_mixed_cands`, `index_mixed_cma_reached`, `index_mixed_mc_ran`. Each
`index_log` entry also lists, per candidate, the CMA hit ratio and whether MC ran.

## Results

Reading line: each row is one library run compared with the stored script run of the same name.
"(=)" means equal to the stored value. All runs are single runs with one worker, on a machine with
other jobs running (contended); no wall time is a result.

| Check | Sample | Keys on | Evaluations before the pass | Seeds / reseeds | Grid fallbacks | Wrong voxels before the pass |
|---|---|---|---|---|---|---|
| 1 | big | S8 keys only | 7,201,863 (=) | 37 / 292 (=) | 0 | 67 (=) |
| 2 | big | + no MC fallback (vs `S10`) | 4,767,135 (=) | 36 / 297 (=) | 0 | 77 (=) |
| 3 | big | + `SeedIndexCMAFirst 200` (vs `S11`) | 2,879,705 (=) | 35 / 295 (=) | 0 | 70 (=) |
| 4 | first gradient | both keys (vs `S11`) | 8,622,191 (=) | 348 / 986 (=) | 11 (=) | 207 (=) |

In every row the per-voxel records (24,570 voxels), every numeric counter of the stored run, the
map before the pass (byte for byte) and the log of the seed steps were equal. Against the stored
map, 0 voxels were fixed and 0 broken before the pass in all four rows.

Reading line: counts of CMA fits, seed steps and mixed candidates, library against stored.

| Check | CMA fits in neighbour fits | `cma_not_fell` | Seed steps (stopped early) | Mixed candidates (CMA reached / MC ran) |
|---|---|---|---|---|
| 2 (big) | 32,398 (=) | 3,152 (=) | 333 (323) (=) | key not on |
| 3 (big) | 32,645 (=) | 3,212 (=) | 330 (327) (=) | 436 (266 / 170) (=) |
| 4 (gradient) | 53,345 (=) | 4,770 (=) | 1,334 (1,233) (=) | 2,158 (1,105 / 1,053) (=) |

On the gradient sample the MC optimizer ran for 1,053 of 2,158 candidates. The sequence of
(CMA hit ratio, MC ran) matched the stored log candidate by candidate. This run has 11 grid
fallbacks, so the mixed refinement was tested next to grid fallbacks, which `big` cannot do.

### The boundary pass (check 5)

Reading line: the "from stored map" columns start the library pass from the stored map of `S11`.
The "inside the run" columns are the pass inside the library run.

| Sample | From stored map: tries / evaluations | Total with pass | Planned total | Inside the run: tries / evaluations | Wrong after the pass (library / stored) |
|---|---|---|---|---|---|
| big | 18,617 / 1,135,637 (=) | 4,015,342 | 4,015,342 | 18,617 / 1,135,637 | 19 / 18 |
| first gradient | 60,423 / 3,685,803 (=) | 12,307,994 | 12,307,994 | 60,423 / 3,685,803 | 62 / 62 |

From the stored maps, every try (voxel, candidate, end, cost, hit ratio, confidence, evaluations)
and the set of changed voxels (983 on `big`, 681 on the gradient sample) were equal.

Inside the run, the number of tries and of evaluations equals the stored one for `S11` on both
samples. Equal counts do not mean equal contents. The next table compares the two passes try by
try. It is computed from the stored files by `keys_pass_contents.py`
(`benchmarks/autoresearch/ar5/keys/pass_contents.json`).

Reading line: one row per library run, against the stored pass of the same run. "Another
candidate" means the try started from another neighbour's orientation. The angle between two ends
is the rotation between them, with no symmetry reduction. "n/a" means the two lists of tries have
different lengths, so they cannot be paired try by try.

<!-- table:pass_contents -->
| Run | Sample | Tries (library / stored) | Tries with another candidate | Costs differ | Largest angle between ends (deg) | Ends over 0.01 deg apart | Changed voxels (library / stored) | Changed in one run only (library / stored) |
|---|---|---|---|---|---|---|---|---|
| S8L | big | 18,664 / 18,665 | n/a | n/a | n/a | n/a | 1,084 / 1,056 | 279 / 251 |
| S10L | big | 18,547 / 18,543 | n/a | n/a | n/a | n/a | 1,135 / 1,087 | 397 / 349 |
| S11L | big | 18,617 / 18,617 | 1 | 4,372 | 1.22 | 2,884 | 986 / 983 | 32 / 29 |
| S11L | first gradient | 60,423 / 60,423 | 0 | 15,813 | 1.11 | 10,541 | 688 / 681 | 28 / 21 |
<!-- /table:pass_contents -->

For `S11` on `big`, one try started from another candidate: at voxel 8415, candidate 1819 in the
library run and 8417 in the stored pass. Costs differ in 4,372 of 18,617 tries on `big` and in
15,813 of 60,423 tries on the first gradient sample. The ends are up to 1.22 degrees apart on
`big` and 1.11 degrees on the first gradient sample. Most are much closer: 2,884 and 10,541 tries
have ends more than 0.01 degree apart. The pass changed 986 voxels in the library run against 983
in the stored pass on `big` (32 only in the library run, 29 only in the stored pass). On the first
gradient sample it changed 688 against 681 (28 and 21).

After the pass, `big` has 19 wrong voxels in the library run and 18 in the stored map (0 fixed, 1
broken against the stored map). The first gradient sample has 62 in both (0 fixed, 0 broken).

The other two `big` runs differ in their counts as well. `S8L` had 18,664 tries against 18,665
(evaluations 1,138,504 against 1,138,565); after the pass it had 24 wrong voxels, as stored (5
fixed, 5 broken). `S10L` had 18,547 tries against 18,543 (evaluations 1,131,367 against
1,131,123); after the pass it had 26 wrong voxels against 22 stored (4 fixed, 8 broken).

These differences are consistent with one cause: the pass inside a run reads the unrounded map,
while the stored pass read the map of the stored `.mic` file, whose angles are rounded. That cause
was shown only for `S8L` (and `D4L`), in the AR3 step 10 fixes round (D-20261111-N9,
`findings_ar3_keys.md` part (3)). This step's `S8L` gives the counts of that round again (18,664
tries and 1,138,504 evaluations). The cause was not tested for `S10L` and `S11L`. The differences
are reported, not repaired.

## Predictions

Reading line: the four predictions of the plan, written before any code.

| # | Prediction | Result |
|---|---|---|
| 1 | Checks 1, 2 and 3 are exact at the first full run of each. | Met. |
| 2 | Check 4 is exact; if one check fails it is this one. | Met. Exact, with 11 grid fallbacks. |
| 3 | Check 5 is exact from the stored maps; the pass inside the run differs by under 0.1% of the pass evaluations. | Met. Exact from both stored maps. Inside the run, `S11` differed by 0 evaluations on both samples; the largest difference in the other two runs was 244 evaluations (`S10L`), about 0.02%. The prediction is about evaluations only. Equal counts do not mean equal contents: inside the run, the costs of 4,372 tries on `big` and 15,813 on the first gradient sample differ from the stored ones. |
| 4 | The step needs 25 to 50 new tests and no existing test changes. | Not met on the count: 62 new tests in the first round, 68 after the fixes round added 6 refusal rows (many are parametrized refusals). No existing test changed. |

The full suite with the keys in place, after the fixes round: 1,724 passed, 47 skipped, 1
deselected, of which 68 tests are in `tests/test_bfs_ar5_keys.py`
(`benchmarks/autoresearch/ar5/keys/suite_counts.json`). The first round counted 1,718 passed with
62 new tests; the 6 added refusal rows make the difference. No golden file changed. With both keys off, the library run of `S8` on `big` gave the stored 7,201,863
evaluations and 67 wrong voxels (check 1).

## What was not tested

- No timing was made. Evaluation counts are costs, not times.
- The keys were not run on `mid` or `held`, nor on a real sample or on realistic images.
- Parity is shown on two samples with clean images, one voxel order each.
- No default was changed. The documentation recommends no setting; that is the owner's choice.
- C++ was not touched.
- The stored runs were not changed and nothing was tuned to make a count match.

## Files

Code: `icenine/ray_indexing.py` (`mixed_refinement`, `seed_by_indexing(refine_fn=)`),
`icenine/reconstructor.py`, `icenine/config_file.py`, `icenine/orientation_search.py`. Tests:
`tests/test_bfs_ar5_keys.py`. Scripts: `scripts/autoresearch/ar5/keys_*.py`, `keys_run_grad.sh`.
Results: `benchmarks/autoresearch/ar5/keys/check_*.json` and `pass_check_*.json`. The contents of
the pass inside the runs (fixes round): `keys_pass_contents.py`, which reads stored files only,
writes `pass_contents.json` and `tables.md`.
