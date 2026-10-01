---
title: "TODO: intensity inputs and distribution-based losses for the orientation NN"
subtitle: "Discussion notes from Stage 3 (2026-09-29); deferred, not started"
date: "2026-09-29"
geometry: margin=1in
fontsize: 11pt
---

# Status

**Deferred.** Raised while reviewing the Stage 3 set-network results
(`MIGRATION_HISTORY.md`, "Toy Orientation NN — Stage 3"). Two questions: (1) does
binarising the windows make learning harder, and should the network see real
intensities instead? (2) since a peak is a distribution over pixels and frames,
should the measure of fit compare distributions? Neither is the current bottleneck
(see "Order" below); this note records the reasoning so it can be picked up later.

# 1. Binarised vs. intensity inputs

## Is binarisation limiting us now? No.

- On the same binarised data, the exact-Bayes posterior mean reaches ~0.004° and
  Gauss–Newton ~0.03° (voxel 0, r⊥ = 12 µm) / ~0.01° (ManyGrains voxel 77,
  r⊥ = 399 µm). The pooled set networks were 5–10× above the Bayes floor; `GNLayerNet`
  (Step 2) is 2–7× above the exact-Bayes median at voxel 0 and 9–16× at voxel 77
  (per bin, `docs/orientation_nn_design.md` Section 8), so the information is in the
  binary data; the networks do not extract all of it.
- For supervised training, binary inputs do not make the loss surface discontinuous:
  gradients are with respect to the weights, not δ. Binarisation only makes several δ
  map to the same input, which is already priced into the Bayes floor.
- The Stage 3 stage-axis failure was a training problem, not missing information:
  the same windows gave the fc net and the sum-pooled set net a good stage axis.

## What binarisation does cost

- **Perpendicular axes:** sub-pixel information. An intensity-weighted centroid is
  more precise than the lit-pixel centroid. Small next to the current errors.
- **Stage axis:** nothing while α = 0. With a perfect crystal (current renderer) a
  spot sits in exactly one frame, so intensity carries no sub-frame ω information.
  With α > 0 (rocking-curve spread, Stage 1 item 4, deferred) a spot spreads across
  adjacent frames and the **intensity ratio between frames locates ω inside a
  frame**. This is the strongest argument for intensity inputs, and it is tied to the
  α > 0 renderer.

## Why not feed the simulator's intensities directly

- The forward model's intensity is geometric (partial-volume overlap): no structure
  factor, strain, detector response, background or noise. Real intensities differ.
- IceNine reconstructs from reduced, thresholded data (`.d` / `.bin` peak files), so
  the deployed input is binary or near-binary.
- A network trained on exact simulated intensities would learn to rely on detail real
  data do not reproduce (sim-to-real gap). This is why the dataset generator stores
  thresholded windows.

## Proposed middle ground

- Soft but physically defensible inputs: the **fraction of each pixel covered by the
  spot**, and with α > 0 the **fraction of the rocking curve falling in each frame**.
  Both are geometric quantities we trust, and they recover the sub-pixel and sub-frame
  information.
- Train with noise and threshold jitter so the network does not rely on exact values.
- Cheap first step that stays binary: give each peak its **second moments** (spot
  covariance: width, elongation, orientation) as measurement features next to the
  centroid (`toy_orientation_model.measurement_features`). This describes each peak as
  a distribution without a new input format or loss.

# 2. Comparing distributions

Two distinct places a distribution comparison fits.

## 2a. Data space: predicted spots vs. observed spots

Render the spots at the predicted orientation and compare each peak's pixel (and
frame) distribution with the observed one, e.g. a per-peak Wasserstein / Sinkhorn
distance.

- **Helps optimisers.** IceNine's hard cost counts pixel overlap, which has zero
  gradient once predicted and observed spots stop touching; that is part of why
  Riemannian Adam and MC stall at larger δ (Stage 2). A transport distance still says
  how far the mass must move, so its gradient stays informative for disjoint spots.
- **Enables label-free training / fine-tuning on real data**, where the true δ is
  unknown: render-and-compare. This is the iNeRF-style refinement parked in
  `docs/research_ideas_nerf.md`.
- **A fit check on real data**, where there is no ground truth to score against.

## 2b. Orientation space: predicted vs. true posterior

The networks already output a distribution (Gaussian mean + full covariance, trained
with NLL). With the exact-Bayes samples (`scripts/exact_bayes_baseline.py`) the output
can also be scored against the true posterior (KL divergence, calibration / coverage),
not only by mean error. Stage 3 showed why this matters: the 30-voxel set net's mean
Mahalanobis² was 3.2–4.1 on training voxels but 17–34 on held-out voxels (overconfident;
superseded: `GNLayerNet` has 3.5–4.5 on held-out clean data, see
`docs/orientation_nn_design.md` Section 7; coverage was still not checked).

# Related training note (done: hypothesis confirmed, 2026-09-30)

The Stage 3 stage-axis failure was a Gaussian-NLL pathology, not an input problem: the NLL
gradient on the mean is scaled by 1/σ², so once σ_⊥ ≈ 0.01° and σ_z ≈ 0.5° the z-mean gradient is
~2500x weaker (the variance ratio (0.5/0.01)²) and z stays at "don't know" (Seitzer et al. 2022). Tested on voxel 0
(`MIGRATION_HISTORY.md`, "Stage 3 — NLL diagnostic"; runs in
`benchmarks/toy_orientation_stage3/nll/`):

1. Per-sample `--beta-nll 0.5 / 1.0`: **no effect** on z (final-epoch z RMS 0.42-0.43°). It
   scales each sample by one scalar and cannot rebalance z against ⊥ inside a sample.
2. `--loss decoupled` (MSE on the mean + NLL of the covariance at stopgrad(mean)): **fixes z**
   with meanmax pooling (z RMS 0.024-0.040°, median 0.023-0.034° at all magnitudes; ⊥ 0.010-0.014°).
   `--loss mse-then-cov` gives the same. `--pool all` + decoupled is best (0.020-0.030°).
3. Frame-only probe (`--arch probe`): recovers z to 0.023-0.046° with mean pooling, so the
   information survives mean pooling.
4. On the far voxel the loss lets meanmax reach 0.019-0.024° (it was at the prior), but does
   not beat plain NLL with `--pool all` (0.013-0.018°).

The earlier explanations (max pooling discards counts; mean pooling divides by the peak count;
sum pooling rescales the signal) are superseded: sum pooling only compensated for the loss.
Use `--loss decoupled` for the offset head from now on. [Update 2026-10-01: the 30-voxel
network has since been retrained with it (decoupled set net, then `GNLayerNet`, which replaced it);
see `docs/orientation_nn_design.md`.]

# Order

1. Fix training: done, decoupled loss above (the 30-voxel network was retrained with it; `GNLayerNet` replaced the set net).
2. Per-peak second moments as extra measurement features (cheap, stays binary).
3. α > 0 renderer with geometric per-pixel and per-frame fractions as soft inputs;
   noise and threshold jitter.
4. Per-peak Sinkhorn render-and-compare loss: refinement stage and label-free
   training on real data (with `research_ideas_nerf.md`).
5. Posterior scoring against exact Bayes (KL, coverage) in the evaluation script.
