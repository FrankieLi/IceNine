---
title: "TODO: multiple detectors with many grains and voxels far from the rotation axis"
subtitle: "Deferred investigation, raised after the Step 3b pairing fix (2026-10-01); not started"
date: "2026-10-01"
geometry: margin=1in
fontsize: 11pt
---

# Status

**Deferred, not started.** Raised by the project owner after detector pairing showed no
benefit in the toy study (`MIGRATION_HISTORY.md`, Architecture Step 3b;
`docs/orientation_nn_design.md` §3.4, §7). The expectation is that multiple detectors
matter most where the toy study has not looked yet: **many grains lighting the detector**
and **voxels far from the rotation axis**.

# Why the toy study is not the right test

- **Clean data:** `GNLayerNet` already puts every peak's Jacobian from both detectors into
  one set of normal equations, so the linear solve already combines the two views. Gauss–
  Newton with both detectors improves ⊥ only 1.0–1.8× and δ_z by 0–37% over one detector
  (`benchmarks/toy_orientation_arch/diag_multi_summary.txt`): with a single scatterer the
  detectors are largely redundant.
- **Distractors:** the toy neighbours are the 2 nearest voxels (9.4 µm away, same grain for
  23 of 30 targets) plus a Σ3 twin of the nearest neighbour. Their diffracted rays start
  almost where the target's do, so a two-detector ray-consistency check cannot separate
  them. Encoder-level pairing (fixed in Step 3b) and the inert version both performed
  within seed spread of the unpaired network.
- **Geometry:** most test voxels sit within r⊥ ≤ 0.5 mm of the axis in a 1 mm sample, and
  only a few voxels contribute signal to each window.

# Hypothesis to test

With many grains, a target window also collects spots from **distant grains**, whose rays
originate far from the target voxel. Spot positions on L₁ and L₂ are related by
x_⊥ + (L − x_∥)·k̂′_⊥/k̂′_∥: a spot is from the target only if its L₁–L₂ pair points back
to the target's position. This geometric filter should reject distant-grain spots that
single-detector data cannot, and the effect should grow with the number of grains and
with r⊥ (larger parallax lever arm; the z column of the spot Jacobian is pure parallax,
∝ r⊥, TN §3.6.3).

# Suggested investigation

1. **Data:** render target windows with spots from *all* grains of a large sample
   (e.g. ManyGrains `rand_500grains_1mm_inFZ.mic`, or a denser/larger synthetic sample),
   not only nearest neighbours; sweep the number of contributing grains and the target's r⊥
   (including r⊥ > 0.5 mm with a larger sample).
2. **Baselines:** Gauss–Newton and Huber Gauss–Newton with detector 1 only, detector 2 only,
   and both (detector mask exists: `CentroidGaussNewton`), plus an explicit ray-consistency
   pre-filter (reject spots whose L₁–L₂ line misses the voxel by more than a tolerance).
3. **Networks:** `GNLayerNet` unpaired vs `--pairing` (fixed, Step 3b), trained with the
   many-grain distractors; check whether the pair path's learned weights down-weight
   inconsistent pairs.
4. **Metrics:** error vs number of grains and vs r⊥, per axis (z vs ⊥), held-out voxels,
   mean Mahalanobis² and coverage; ≥ 2 seeds.
5. **Real data (later):** the same comparison on experimental data with both detector
   distances, where far-field contamination is real.

# Related

- `docs/orientation_nn_design.md` §2.7 (distractor model), §3.4 (pairing), §8 (limitations).
- `docs/todo_findoptimal_wrong_candidates.md` (neighbour/twin confusions in the search).
- `docs/research_ideas_nerf.md` (render-and-compare with full-sample rendering).
