---
title: "TODO: full-sample BFS end-to-end test of the 500-grain sample with new orientations"
subtitle: "Classic BFS vs BFS with the NN and hybrid finisher: timing and accuracy (planned 2026-10-07)"
date: "2026-10-07"
geometry: margin=1in
fontsize: 11pt
---

# Status

**Planned, not started.** Requested by the project owner on 2026-10-07 as the next step after the "Finisher and MC
study". It is Phase D of that plan (MIGRATION_HISTORY, "Finisher and MC study (2026-10-07)"). It waits for:

- the study's Phase C, which chooses the new finisher;
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
| Classic BFS | Coarse search → quick MC → FindOptimal (MC + VarianceMinimizing), C++ parity | Classic local MC refinement (`_fit_from_seed`) |
| New BFS | Coarse search with the proxy rerank → the study's new finisher | Hybrid: net ×3 → finisher (H3); HG when the start may be ≥ 1.5° off |
| Ablation A | New | Classic |
| Ablation B | Classic | New |
| C++ IceNine BFS | (external reference on the same images) | |

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

**Compute estimate.** About (500 seeds × 25 s + 24,000 neighbours × 1.5 s) / 10 workers ≈ 1.4 h per pipeline and
variant. The forward simulation is still to be estimated. The pilot first.

# Risks

- **Zero-shot networks:** they were trained on per-voxel windows. If they underperform on full-sample images, report
  the zero-shot result and a retrained one.
- **Error propagation:** BFS can carry errors across grain boundaries. Detection relies on the final cost.
- **Memory:** full-sample image stacks need memory, and the forward-simulation time is unknown.
- **C++ reference:** matching it needs an aligned config (eta limit, Q_max, data dir).

# Relation to other work

- `docs/findings_orientation_search_2026-10.md`: the findings this test checks at full scale.
- MIGRATION_HISTORY, "Finisher and MC study (2026-10-07)": Phases A–C choose the finisher used here.
- `docs/todo_hybrid_nn_refinement.md`: its one open item (full-sample renders) is addressed here.
- `docs/todo_future_ideas_nn_active_fourier.md`: Idea 1 (no-start network) would replace the seed-voxel coarse search.
