---
title: "TODO: FindOptimal sometimes returns a wrong orientation"
subtitle: "Findings from the Stage 1 coarse-residual measurement (2026-09-29); cause not yet identified"
date: "2026-09-29"
geometry: margin=1in
fontsize: 11pt
---

# Status

**Investigated (2026-10-06); see `docs/findoptimal_robustness_report.md`.** Cause: at the pruning
step of each level (keep the best quarter by the post-quick-MC local cost) the candidate nearest the
truth, 1-2 degrees off, has a cost near 1 (the cost is sharp to 0.2-0.5 degrees) and loses a ranking
among near-1 costs to a wrong candidate, usually a CSL relative (Sigma3 62-66% of the wrong answers).
In 200 voxels x 3 seeds, 34.0% (clean) and 23.8% (realistic) of the full reconstructions are wrong;
in 83-88% of the wrong runs the truth basin is present and pruned (S2), 4-6% it is never a level-0
candidate (S1), 4-13% it reaches FindOptimal and is lost there (S3). Hypothesis 1 (coarse-cost
settings) is only partly right: coarse Q_max 8 from level 0 helps a little, a smaller pixel
tolerance is much worse; hypothesis 2 (FindOptimal logic) explains only the S3 minority. Fixes: a
CSL-relative check after the search (F1) takes the wrong rate to 3.5% / 3.3% for +11% evaluations; the
CSL relatives added at every level (F1b) gives 0% / 1.5%; keeping more candidates (F2a) 11% / 5.5%; a
learned one-pass classifier as the pruning rank 4.0% / 4.5%, not better than F1. Follow-up ideas
(a Q_max-8 cost proxy from low-Q_max evaluations; cheap de-duplication of candidates):
`docs/todo_coarse_cost_proxy_and_dedup.md`. The original
investigation notes below are kept as written.

The fixes exist as opt-in knobs (F1b: `extra_candidates`; F1: `scripts/findoptimal_robustness/f1_run.py`), default off by owner decision (2026-10-06) to match the C++ behaviour; see `MIGRATION_HISTORY.md`, "FindOptimal robustness".

# Observation

`scripts/checks/measure_coarse_residual.py` runs `AdaptiveVoxelReconstructor` with
`Examples/Example2.ThreeVoxels/ConfigFiles/ReconstructQ8.config` (Q_max = 8) against
`ScatteringData_Python`, 5 seeds per voxel (15 runs; data in
`benchmarks/toy_orientation_stage1/coarse_residual.npz`).

- The best candidate after the coarse levels is a wrong solution in 4 of 10 runs on
  voxels 0 and 1, even though a correct candidate (0.07–1.09° off) is usually among
  the ones passed on.
- FindOptimal then sometimes keeps a wrong solution despite receiving a good one:
  voxel 1 seed 0 (good hand-off 0.35°, final 54°) and seed 1 (good hand-off 1.09°,
  final 42.5°).
- Voxel 2: the coarse search fails on every seed (best hand-off ~60° off; final
  error ~60° in all 5 runs).

# What the wrong answers are (checked)

- **Voxels 0 and 1** are the two halves of one 1.5 µm rhombus (one triangle points
  up, one down; they share an edge at the same position). Their true orientations
  differ by 53.94°. Every ~54° wrong hand-off for voxel 0 is **voxel 1's
  orientation** to within 0.04–0.51°. Both voxels' spots are in the simulated data
  and, at 1.48 µm pixels, land within about a pixel of each other, so the
  neighbour's orientation is a plausible match.
- **Voxel 2's** wrong solutions are the copper **Σ3 twin** of its true orientation:
  59.6–60.0° about a ⟨111⟩ axis in the crystal frame (within 0.04–0.58° of ⟨111⟩) in
  all 5 seeds. A twin shares a subset of reflections with the true orientation.
- Voxel 1 seed 1's 42.5° final is neither a neighbour nor a twin; unexplained.

# What is NOT explained: why the search prefers them

The hard cost function evaluated at each orientation (default settings, all
reflections up to the config's Q_max, no pixel-radius tolerance) **prefers the truth**:

| Voxel | Cost at true orientation | Cost at the wrong answer |
|---|---|---|
| 0 | 0.190 | 0.331 (voxel 1's orientation) |
| 1 | 0.186 | 0.362 (voxel 0's orientation) |
| 2 | 0.333 | 0.820 (its twin, seed 2) |

So the data can distinguish them, and the failure is in the search, not in an
inherently ambiguous cost. Hypotheses to test:

1. **Coarse-stage cost settings.** The coarse levels use a smaller Q_max
   (`n_q_max` starts at 5 + `min_local_resolution` and increases per level) and
   `pixel_radius=3`. With a 3-pixel tolerance and few reflections, a neighbour's or
   a twin's spots may match nearly as well as the truth. My check above did not
   reproduce those settings.
2. **FindOptimal logic.** How it picks among candidates and what it compares (costs
   computed with different settings? early stop on a convergence threshold, e.g.
   `max_convergence_cost`?), and whether it can return a candidate it just refined
   away from the good one.
3. **Voxel geometry.** Whether the cost for voxel 1 credits spots that belong to voxel
   0 (the two triangles overlap in their projected footprints), so that voxel 0's
   orientation earns a better score than it deserves.
4. **Voxel 2 specifically.** Its 0.75 µm side is about half a pixel, so its spots are
   tiny and its true cost is already 0.33; check whether the coarse grid ever
   samples close enough to the truth (the closest hand-off was 2.9–3.3° in two runs).

# Suggested investigation

1. Pick one failing case (voxel 1, seed 0 or 1) and log, per level, each candidate's
   orientation, cost and error to ground truth, plus FindOptimal's per-candidate
   result and its selection rule.
2. Re-evaluate the final candidates with the same cost settings FindOptimal used,
   to see whether the ranking reverses between stages.
3. Repeat with `pixel_radius` and the level-wise Q_max varied, and with the
   neighbouring voxel's spots removed from the data, to isolate hypotheses 1 and 3.
4. Check the C++ reconstruction on the same case: the MIGRATION history records that
   voxel 2 fails in C++ too, which suggests part of this is algorithmic rather than
   a port bug.

# Relevance to the toy NN

A network refiner should be evaluated per candidate, as FindOptimal is, and should
be able to reject a wrong candidate that the cost function accepts. Generalizing
across voxels (Stage 3) makes twin and neighbour confusions a realistic test case.
