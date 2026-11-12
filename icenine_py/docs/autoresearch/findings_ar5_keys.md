# AR5 step 5: two opt-in keys for run S11

Task `ar5e`, 2026-11-11 (plan: `PLAN.md`, "AR5 step 5"; decision D-20261111-N21).

## Summary

Run `S11` used to exist only in scripts. It is `S8` with two changes. Both are now keys of the
library, both off by default. This page shows that a library run gives the same result as the
stored script run, on two samples.

"Same result" (parity) means: the same evaluation counts, the same seed steps and reseeds, the
same per-voxel records and the same map file, byte for byte. All four checks of the BFS part
were exact. The boundary pass was exact when started from the stored maps. Inside a library run
the pass gave the stored tries and evaluations for `S11` on both samples; for the other two runs
on `big` it differed by 1 and 4 tries out of about 18,600.

No setting is recommended. No timing was made. See "What was not tested".

## Terms

- **Evaluation**: one call of the cost function. All counts here are evaluations.
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
  fit of 200 counts 201 evaluations. Refused unless `SeedByIndexing 1` is set and the stop ratio
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
The "inside the run" columns are the pass inside the library run, which reads the unrounded map.

| Sample | From stored map: tries / evaluations | Total with pass | Planned total | Inside the run: tries / evaluations | Wrong after the pass (library / stored) |
|---|---|---|---|---|---|
| big | 18,617 / 1,135,637 (=) | 4,015,342 | 4,015,342 | 18,617 / 1,135,637 | 19 / 18 |
| first gradient | 60,423 / 3,685,803 (=) | 12,307,994 | 12,307,994 | 60,423 / 3,685,803 | 62 / 62 |

From the stored maps, every try (voxel, candidate, end, cost, hit ratio, confidence, evaluations)
and the set of changed voxels (983 on `big`, 681 on the gradient sample) were equal.

Inside the run the tries and evaluations equal the stored ones for `S11` on both samples. The
ends differ in the last digits because the run reads the unrounded map; the largest difference in
one entry of an end matrix is 0.021 on `big` and 0.018 on the gradient sample. On `big` the
library has 19 wrong voxels after the pass and the stored map 18 (0 fixed, 1 broken against the
stored map). On the gradient sample: 0 fixed, 0 broken.

The other two `big` runs differ inside the run. `S8L` had 18,664 tries against 18,665 (evaluations
1,138,504 against 1,138,565); after the pass 24 wrong voxels, as stored (5 fixed, 5 broken).
`S10L` had 18,547 tries against 18,543 (evaluations 1,131,367 against 1,131,123); after the pass
26 wrong voxels against 22 stored (4 fixed, 8 broken). Step 10 found the same kind of difference.
It is reported, not repaired.

## Predictions

Reading line: the four predictions of the plan, written before any code.

| # | Prediction | Result |
|---|---|---|
| 1 | Checks 1, 2 and 3 are exact at the first full run of each. | Met. |
| 2 | Check 4 is exact; if one check fails it is this one. | Met. Exact, with 11 grid fallbacks. |
| 3 | Check 5 is exact from the stored maps; the pass inside the run differs by under 0.1% of the pass evaluations. | Met. Exact from both stored maps. Inside the run, `S11` differed by 0 evaluations on both samples; the largest difference in the other two runs was 244 evaluations (`S10L`), about 0.02%. |
| 4 | The step needs 25 to 50 new tests and no existing test changes. | Not met on the count: 62 new tests (many are parametrized refusals). No existing test changed. |

The full suite with the keys in place: 1,718 passed, 47 skipped, 1 deselected. No golden file
changed. With both keys off, the library run of `S8` on `big` gave the stored 7,201,863
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
Results: `benchmarks/autoresearch/ar5/keys/check_*.json` and `pass_check_*.json`.
