---
title: "The batched low-Q feature pass"
subtitle: "How FeatureExtractor.features_batch works, what it must keep identical, and when it helps"
date: "2026-10-07"
geometry: margin=1in
fontsize: 11pt
---

# 1. What the low-Q features are

`FeatureExtractor` (`scripts/findoptimal_robustness/features.py`) turns one candidate orientation into
a feature vector from **one forward pass** against the detector images. For each eligible
(peak, detector) pair it records the reflection, its |q| family, the detector and whether the spot lands
on lit pixels: at the centre pixel (`hit0`), at the centre or any of the three triangle vertices
(`hit_any`), and anywhere in a +-1 or +-3 pixel box (`hit1`, `hit3`). The features are aggregates of
these hits (overall, by |q| family, by detector, "CSL-aware" aggregates over reflections shared with
Sigma relatives) plus two cost columns, `cost_local` (pixel radius 0) and `cost_global3` (radius 3).

`q_max` (constructor argument, Angstrom^-1, default `None`) restricts the extractor to the reflections
with |q| <= `q_max`. The peak table, the |q| families, the aggregates and the two cost columns
(from `VoxelCostFunction`s built with `max_q = q_max`) then involve only those reflections, so a
low-Q vector needs no high-|q| work. The study uses Q5 (26 of 112 reflections). The proxy of the
pruning rerank (`scripts/coarse_proxy/`) is a model on these features plus the free Q_max-8 cost of the
candidate.

# 2. The per-candidate path (the reference)

`FeatureExtractor.features(R, vertices, phase)`:

1. `peak_table(R, vertices)`: rotate the reflection vectors, `g_lab = (R @ g_hkl.T).T`; get the omega
   solutions and the observable mask (`get_scattering_omegas_torch`); keep omega1 and omega2 of
   each observable reflection.
2. Per peak, build the sample-to-lab rotation for that omega, the reflected-beam direction, the
   eta cut (`eta_limit`, `min_sin_eta`) and the omega -> frame lookup through `range_map.index_list`.
3. Per detector: intersect the rays from the three voxel vertices with the detector plane, convert to
   fractional column/row, take the centroid, discard pairs outside the detector.
4. Look the centre pixel, vertex pixels and the +-1 / +-3 boxes up in the sorted set of lit-pixel keys
   (`set_image(keys)`, `searchsorted`).
5. Reduce: means of the hit columns over all pairs, per family, per detector and per Sigma
   relative; then two further cost evaluations (`evaluate` at pixel radius 0 and 3) that redo
   the geometry internally.

Each call is many small tensor operations, so most of its roughly 1.5 ms is Python and
per-call overhead, not arithmetic (geometry only 0.742 ms against the whole pass 1.52 ms in the Task 2
timing).

# 3. The batched path

`features_batch(Rs, vertices, phase, with_cost=True, chunk=256)` returns the same matrix as stacking
`features(R, ...)` for every row of `Rs` (shape `(B, 3, 3)`).

**Vectorised across candidates and peaks.** `_batch_geometry` runs the arithmetic of `peak_table` once
for all B candidates: one omega solve on `B*K` reflection vectors, one batched rotation, ray-plane
intersection and lit-pixel lookup on all eligible (candidate, peak, detector) rows, each row carrying
its candidate index. `_fill_aggregates` replaces the per-candidate means by integer counts with
`np.bincount` (exact in float64, then one division), including the per-family, per-detector and
per-Sigma aggregates. The same geometry also feeds the two cost columns, which the per-candidate path
recomputes twice.

**Still done per candidate.** (i) The rotation `R @ g.T` (one 2-D matmul per candidate, see Section 5).
(ii) The cost columns: `_batch_costs` calls the C stage-D routine `_c_stage_d_overlap` once per candidate
and per pixel radius (0 and 3), on the rows of that candidate in the cost function's order, because the
quality is a running mean over peaks and its last bits depend on the order. The geometry that routine
needs is passed in from the batched pass, but the call itself is not batched. (iii) The loop over
relatives inside `_fill_aggregates` (a few Sigma relatives).

**Fallbacks and bookkeeping.** With `with_cost=True` and no C stage-D extension, it falls back to the
per-candidate loop. With `with_cost=False` the two cost columns are NaN. The cost-function evaluation
counters are advanced as the per-candidate path advances them (`+B` on each of the two low-Q cost
functions, or `+2B` on the local one when `q_max` is `None`). An empty batch returns `(0, n_features)`.

**Chunking.** `chunk` (default 256) bounds the working size: a larger batch is split and the pieces are
concatenated (recursion on `features_batch`). Rows are independent, so chunking is exact (tested with
`chunk=7`). The memory of a chunk is dominated by the `(B*K, 3)` reflection array and the per-pair
arrays for the +-3 boxes (49 lookups per pair).

## Sketch

```
features_batch(Rs, chunk):
    if no C stage-D and with_cost: return stack(features(R) for R in Rs)
    if len(Rs) > chunk: return concat(features_batch(piece) for piece in chunks(Rs))
    g_lab   = cat([R_i @ g.T for i in candidates], dim=1).T      # one 2-D matmul each, strided
    omegas  = scattering_omegas(g_lab)                            # all candidates, all reflections
    rows    = eligible (candidate, peak, detector) rows           # eta cut, frame lookup, detector
              ordered per candidate: omega1 rows, then omega2 rows
    hits    = lit-pixel lookups of centre, vertices, +-1 box, +-3 box     # vectorised over rows
    X       = bincount aggregates per candidate (overall, family, detector, Sigma)
    for each candidate b:                                         # still sequential
        for radius in (0, 3):
            X[b, cost] = 1 - stage_d_overlap(rows of b, radius)
    return X
```

# 4. The two bit-identity pitfalls

The batched features must equal the per-candidate ones exactly, not only to a tolerance: a one-ulp
difference in an omega can move a spot across a pixel or frame boundary, which flips a hit
(and so a hit rate by one peak count) or a cost. Two details were found necessary.

1. **One 2-D matmul per candidate, not a batched `torch.matmul`.** A batched matmul selects another
   kernel and differs from the 2-D one in the last bit (52346 of 271152 elements in one test). The
   rotation is therefore built as `torch.cat([Rt[i] @ gT for i in range(B)], dim=1).T`. This is a loop
   over candidates, but a cheap one.
2. **Keep the strided `(R @ g.T).T` layout.** The per-candidate path computes `(R @ g.T).T`, a transposed
   view of a contiguous result; the batched code reproduces that layout. A contiguous copy changes the
   last bit of the omegas.

Intermediate states recorded during development (5500 ad hoc candidates): with a batched matmul 5 rows
differed (up to 0.037 in a hit rate, 0.003 in a cost); with the per-candidate matmul but a contiguous
layout, 2 rows (cost columns). With both fixes the results are exactly equal. If a future torch version
changes the kernel selection, `test_features_batch_equals_per_candidate` is the place that will show it.

# 5. How equality is tested

`tests/test_coarse_proxy.py::test_features_batch_equals_per_candidate` builds the candidates from a
voxel's truth, random rotations and a sample of its stored E2 candidates, computes the per-candidate
matrix as the reference, and asserts `np.array_equal` (not `allclose`) for batch sizes 1, 7 and the whole
set, with Q_max 5 and `None`. It also checks the evaluation counters, `chunk=7` recursion,
`with_cost=False` (NaN cost columns, equal aggregates) and the empty batch. The timing script
(`scripts/coarse_proxy/timing_batch.py`) re-checks equality on its own 5000 candidates (25 cases x 200,
Q_max 5): 5000 of 5000 rows identical. The end-to-end check: a batched re-run of the keep-1/4 proxy
(`p_k4b`) reproduced all 400 stored runs exactly (`R_final`, evaluation counts). Not tested: other
detector geometries or crystal phases than the study's.

# 6. Timing

Single-worker (all thread counts 1, `preflight.require_quiet()` record saved as
`benchmarks/coarse_proxy/timing_batch_preflight.json`), Q5 pass with both cost columns, 25 random
cases x 200 candidates, 3 repetitions, methods interleaved, median over cases
(`benchmarks/coarse_proxy/timing_batch_tables.md`). `equiv` = evaluation equivalents of one Q_max-8
local evaluation (0.529 ms in the same run).

| path | batch | ms per candidate | speedup | equiv | proxy equiv |
|---|---|---|---|---|---|
| per-candidate (reference) | 1 | 1.555 | 1.0 | 2.94 | 2.95 |
| batched | 1 | 1.074 | 1.4 | 2.03 | 2.04 |
| batched | 50 | 0.151 | 10.3 | 0.28 | 0.30 |
| batched | 200 | 0.133 | 11.6 | 0.25 | 0.26 |

"Proxy" adds the prediction of the gradient-boosted model (0.0064 ms per candidate at batch 200).
Batch 1 is the batched code called per candidate, so its 1.07 ms is the fixed per-call overhead;
geometry alone costs 0.734 ms per candidate per-candidate and 0.124 ms batched at 200. What does not
speed up is the per-candidate stage-D call.

**When batching helps.** Only large batches gain. In a proxy run the scoring is called once per level:
level 0 has about 220 candidates per call (mean 220.9, range 95-562), later levels are small (means at keep 1/4:
63.3, 20.7 and 7.1 candidates at levels 1, 2 and 3; at keep 1/8: 32.0, 5.1 and 1.3), so only the first
call is in the fast regime. Even so the proxy took about 0.05 s per run in the single-worker timing, about 0.2% of a
run, against 0.37 s / 0.59 s per run per-candidate (Task 2 accounting). The remaining cost of a run is
the reconstructor's own evaluations, so batching removes the proxy's overhead but does not by itself
make a proxy run faster than the baseline (the keep-1/4 arm was 3.3% slower than baseline).

# 7. Using it, and why the per-candidate path remains the default

```
uv run python scripts/coarse_proxy/endtoend.py run --batched ...
```

`--batched` sets `cfg["batched"]`, and `proxy_features(..., batched=True)` calls
`fe.features_batch(np.stack([c.orientation ...]), vertices, phase)` once per `rank_key` call, then appends
the free Q8 cost column. Results are identical to the default, so use it for production proxy runs.

The per-candidate `features` stays the default and the **reference**: it is the readable definition
of the features, it does not need the C stage-D extension, and every equality test compares against it. If you change the features, change
`features` and `peak_table` first and make `features_batch` agree with them, then rerun the equality test.
A new geometry, phase or `q_max` should get an equality check before the batched pass is trusted there.
