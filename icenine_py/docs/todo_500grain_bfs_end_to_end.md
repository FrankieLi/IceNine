---
title: "TODO: full-sample BFS end-to-end test of the 500-grain sample with new orientations"
subtitle: "Classic BFS vs BFS with the NN and hybrid finisher: timing and accuracy (planned 2026-10-07)"
date: "2026-10-07"
geometry: margin=1in
fontsize: 11pt
---

# Status

**Data step done (2026-10-07); pilot done (2026-10-08, 2,000-voxel region, 5 runs; MIGRATION_HISTORY, "Phase D pilot results"): no-go for the full runs as configured, 28 of 60 grains were never seeded.** Requested by the project owner on 2026-10-07 as the next
step after the "Finisher and MC study". It is Phase D of that plan (MIGRATION_HISTORY, "Finisher and MC study
(2026-10-07)"). Done: the new-orientation sample, the full-sample clean and realistic images and the four
Python-only reconstruction configs (MIGRATION_HISTORY, "Phase D data"; `scripts/phase_d/`). No full-sample images
existed before. The full render takes about 1 min on 10 workers (contended). The BFS runs still wait for:

- the library work on the Phase D branch (neighbour budget, retry sigma, BFS refit keys);
- the owner's go-ahead for the long runs.

# What

Reconstruct the whole 500-grain Cu sample end to end with BFS reconstruction (`BFSReconstruction`), on full-sample
forward-simulated images. Compare two pipelines on timing and accuracy:

- **Classic BFS:** the C++-parity optimizers.
- **New BFS:** the network and hybrid finisher.

**Sample.**

- `Examples/Example2.ManyGrains/SimInput/rand_500grains_1mm_inFZ.mic`: one layer, 24,570 voxels, 497 distinct grain
  orientations.
- No second 500-grain sample exists. So keep the geometry and the voxel → grain map, and draw **new orientations**:
  one per grain, uniform on SO(3) with a fixed seed, reduced to the fundamental zone.
- Write `rand_500grains_1mm_neworient_s<seed>.mic` (same columns) plus the grain map.
- Check that no new orientation is within 1° of its old one, or of a neighbouring grain's new one.
- This makes the test orientation-novel for the networks: they were trained on the old orientations. It is not
  geometrically novel.

**Images.**

- Run the Python forward model (pixel-exact against C++) over the full sample, at both detector distances, as in
  `Example2.Simulation.config`.
- Make two variants:
  - clean;
  - realistic: detector noise added. Overlap is now physical, from neighbouring grains.
- This removes the per-voxel-image caveat of all earlier studies.
- Earlier notes wrongly claimed ManyGrains images exist. Check first.

**Pipelines.** All use BFS. A full no-start search runs only where BFS already runs one: seed voxels of new grains,
and voxels whose inherited start fails the cost check.

| Pipeline | Seed voxels (no start) | Neighbour voxels (inherited start) |
|---|---|---|
| Classic BFS (`mc`) | Coarse search → quick MC → FindOptimal (MC + VarianceMinimizing), C++ parity | Classic local MC refinement (`_fit_from_seed`) |
| CMA BFS (`cma`) | Coarse search with the proxy rerank → CMA-ES finisher | CMA-ES from the inherited start (`CMANeighborMaxEvals 250`), retry with `CMARetrySigma0 1.5` |
| CMA without retry (`cma_noretry`) | As `cma` | As `cma` with `CMARetrySigma0 0`: separates the retry from the optimizer |
| New BFS (hybrid) | Coarse search with the proxy rerank → the study's new finisher | Hybrid: net ×3 → finisher (H3); HG when the start may be ≥ 1.5° off |
| Ablation A | New | Classic |
| Ablation B | Classic | New |
| C++ IceNine BFS | (external reference on the same images) | |

`BFSRevisitRefit 1` (REFIT voxels revisited inside the BFS, as the C++ multi-client run does, instead of in a
post-pass) is planned in the library branch and will be enabled in both the `mc` and `cma` arms. The Phase D configs
(`scripts/phase_d/configs/`) have the `mc`, `cma` and `cma_noretry` arms for each image variant: `clean`,
`realistic` and `realistic_q16` (rendered at Q-max 16, reconstructed at Q-max 8). All read the grid-only mic.

# Why it might help

Every result so far comes from per-voxel images with at most 3 synthetic distractor sources. On those images:

- the reranks cut wrong answers 4–5x for no-start voxels;
- the hybrids are 3–13x faster than the alternatives when a start is 1.5–3° off;
- for starts within 0.1°, FindOptimal alone is best.

BFS gives every non-seed voxel a start, usually within the same grain (U1-like). Across a boundary it gives a wrong
start that must be detected. So the real-world mix of U0/U1/U2 cases, and the physical overlap, decide the net effect.
Only a full-sample run measures it.

# Framing

- **Question:** at matched or better accuracy, is new BFS faster than classic BFS, and where do the time and the
  errors come from?
- **Symmetry-agnostic:** report near-degenerate alternative solutions generally. CSL relatives are the cubic
  instance.

# Experiments (sketch)

1. **Pilot:** a 2,000-voxel region, both pipelines, clean images. Check:
   - the image pipeline;
   - that the net's ROI windows work on full-sample images zero-shot;
   - the per-voxel neighbour cost.
2. **Full runs:** both pipelines and both variants on 10 workers, plus the two ablations.
3. **Speed ratio:** single-worker, preflight-gated timing of classic vs new on a fixed region.
4. **Metrics:**
   - per-voxel misorientation map;
   - wrong (> 1°) rate with Wilson CI;
   - median error;
   - completeness;
   - grain-boundary voxels (within 1 voxel of a boundary) vs interior;
   - near vs far from the rotation axis;
   - grains found and lost;
   - wall time and stage split: seeds vs neighbours, and coarse search vs finisher vs network.
5. **Success:** new BFS has a wrong rate ≤ classic's Wilson upper bound, and median error within +0.005°, at lower
   wall time. Or it has a clearly lower wrong rate at ≤ 10% extra time.

**Metrics notes.**

- **Grain fragments:** with edge adjacency the 497 grains form 612 connected pieces (115 extra pieces of 1 to 3
  voxels). Count "grains found / lost" and the seed numbers with this in mind: an extra piece may need its own seed,
  or be reached only from a vertex neighbour.
- **Truth vs reconstruction mic:** score against `rand_500grains_1mm_neworient_s0.mic`; reconstruct from `..._grid.mic`
  (orientations zeroed).

**Compute estimate.** Measured by `scripts/phase_d/bfs_timing_probe.py` (clean images, classic `mc`, single process,
contended, 4 seeds and 20 neighbours; MIGRATION_HISTORY, "Phase D data"):

- A no-start seed took 177 to 711 s (mean 367 s). With roughly 600 to 1,000 seeds (497 grains, 612 pieces, plus
  failed inherited starts) that is about 60 to 100 h of seed work per serial run. The review's assumed 20 s per
  seed is not supported by this probe; the C++ adaptive step takes about 12.5 s per voxel.
- A classic neighbour step took about 0.4 s: 24,000 neighbours is about 2.7 h.
- So a serial `mc` BFS is dominated by the seeds, not the neighbours. The review's 6 to 9 h (MC) and 4.5 to 7 h (CMA)
  per serial run assume about 20 s per seed; I could not verify that, and I have no CMA timing on the full-sample
  images. The seeds could also be cheaper in the real run if the warm caches or less contention matter (the probe was
  contended and thin), so re-measure with a single-worker preflight-gated pilot before committing the budget.
- Plan: independent concurrent runs (one process per arm and variant) overnight, each on the memmap loader
  (0.37 GB measured after loading vs 7.9 GB dense). `BFSReconstruction` is single-process.

# Risks

- **Zero-shot networks:** they were trained on per-voxel windows. If they underperform on full-sample images, report
  the zero-shot result and a retrained one.
- **Error propagation:** BFS can carry errors across grain boundaries. Detection relies on the final cost.
- **Memory:** solved by the opt-in memmap uint8 loader (`ExperimentalData.from_binary_memmap`); the forward simulation takes about 1 min (Q-max 8) or 4 min (Q-max 16) on 10 workers.
- **C++ reference:** matching it needs an aligned config (eta limit, Q_max, data dir).

# Relation to other work

- `docs/findings_orientation_search_2026-10.md`: the findings this test checks at full scale.
- MIGRATION_HISTORY, "Finisher and MC study (2026-10-07)": Phases A–C choose the finisher used here.
- `docs/todo_hybrid_nn_refinement.md`: its one open item (full-sample renders) is addressed here.
- `docs/todo_future_ideas_nn_active_fourier.md`: Idea 1 (no-start network) would replace the seed-voxel coarse search.
