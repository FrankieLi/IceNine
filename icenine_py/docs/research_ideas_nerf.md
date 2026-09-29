---
title: "Research Idea: NeRF-Style Methods for Orientation Refinement"
subtitle: "Side-quest literature scan for the toy orientation NN (parked, to revisit)"
date: "2026-09-29"
geometry: margin=1in
fontsize: 11pt
---

# Status

A quick, not in-depth literature scan done on 2026-09-29 during the toy
orientation-NN work (after Stage 0; see `MIGRATION_HISTORY.md`). Parked as a
research idea to come back to. Nothing here has been implemented or tested.
References were checked by web search on that date.

# NeRF basics

- **Idea** (Mildenhall et al., ECCV 2020). A small MLP maps a 3D point and viewing
  direction to density and colour; the scene lives entirely in the network's
  weights.
- **Training.** A differentiable volume renderer turns the network into images
  from known camera poses, and the weights are optimised per scene to match the
  photos. There is no training set; it is analysis-by-synthesis.
- **Why it works.** The renderer is smooth, so gradients reach the weights; and a
  Fourier-feature positional encoding of the input coordinates (Tancik et al.,
  NeurIPS 2020) overcomes the MLP's bias toward smooth functions so it can fit
  sharp detail.

# Follow-ups closest to our problem

- **iNeRF** (Lin et al., IROS 2021): inverts a fixed NeRF to estimate a camera
  pose by gradient descent on the pixel residual from an initial guess, sampling
  pixels near interest points. This is the analog of the single-voxel toy
  problem: fixed physics, unknown orientation, loss restricted to regions around
  the peaks.
- **BARF** (Lin et al., ICCV 2021): high-frequency encodings shrink the basin of
  convergence for pose; coarse-to-fine annealing widens it. Same idea as IceNine's
  `MultiScaleImageStack` blurring and the coarse Sukharev grid.
- **Soft rasterisation** (SoftRas, Liu et al., ICCV 2019; 3D Gaussian Splatting,
  Kerbl et al., SIGGRAPH 2023): replace hard rasterisation with smooth splats so
  gradients flow. The direct fix for the zero-gradient problem of the thresholded
  forward model (`nn_inverse_problem_formulation.md` §2.2).
- **Scientific neural fields**:
  - tomography: NAF (Zha et al., MICCAI 2022) and IntraTomo (Zang et al.,
    ICCV 2021) fit a neural field to CT projections with no training labels;
  - cryo-EM: cryoDRGN (Zhong et al., Nature Methods 2021) pairs a neural-field
    volume with an encoder that infers per-image pose and heterogeneity.
- **HEDM specifically**: PARA-X (Cocke et al., arXiv 2609.13619, Sept 2026) refines
  near-field HEDM intragranular orientation and strain fields by gradient descent
  through a differentiable forward model, comparing simulated and measured data
  with an optimal-transport objective because it copes with misaligned peaks. The
  closest published work to where this project could go.

# How this maps onto the toy NN

Near-field HEDM already has NeRF's shape: an unknown field over the sample
(orientation), known "cameras" (the ω angles and detector geometry), a physics
renderer (Bragg condition and ray tracing), and measured images. The single-voxel
toy problem, where everything but the offset δ is known, is iNeRF rather than
NeRF.

Ideas, in rough priority order:

1. **Soft ROI renderer.** Splat each ROI peak as a Gaussian in pixels and spread
   it across frames with an ω-profile, instead of hard-rasterising it.
   `BatchedObserver` (`icenine/orientation_eval.py`) already computes each peak's
   ω\* and spot vertices in batched torch, and the closed-form spot-motion
   Jacobian is in the docs, so this is a modest change. It renders only ROI peaks,
   not full frames, so it is fast.
2. **iNeRF-style refinement after the network.** Stage 0 networks reach about
   0.3–0.5° RMS error; the noise-free floor is about 1e-3°. Use the network's
   prediction to start gradient refinement through the soft renderer ("the
   encoder proposes, the renderer refines"). Loss options: soft-image L2, or
   optimal transport / centroid matching as in PARA-X, which keeps a gradient when
   simulated and measured spots do not overlap (centroid matching is essentially
   the planned Gauss–Newton baseline).
3. **BARF-style annealing.** Start with splat widths that cover the network's
   error (about 8 px along the ring and a few frames for near-axis peaks, from the
   spot-motion Jacobian and 1/|sin η|), then shrink to pixel scale.
4. **Self-supervised training for real data.** With a differentiable renderer the
   network could be trained on unlabeled measured data by penalising re-rendering
   error, as NAF and IntraTomo do. A route to Stage 4, where real data has no
   ground-truth orientations.
5. **Longer term: a true neural field.** Sample position → orientation, fitted to
   all frames at once (whole-sample reconstruction, not the toy). Known weakness:
   grain boundaries are discontinuities, and coordinate MLPs smooth them.

# Caveats

- NeRF relies on a smooth renderer and dense, graded pixel data. HEDM data is
  sparse, near-binary spots, so a soft renderer is an approximation; results must
  still be checked against the hard cost function.
- IceNine's differentiable cost with Riemannian Adam is already a blurred-image
  version of idea 2. The new parts would be rendering only ROI peaks analytically,
  the annealing schedule, and the network initialisation.

# Suggested first experiment

"iNeRF for one voxel": build the soft ROI renderer, then refine from both
predict-nominal and the Stage 0 network's prediction on the same 120 Stage 0 test
cases, and report per-axis error against the exact-Bayes floor. Directly
comparable to the Stage 0 table; roughly a day of work.

# References

- Cocke, C. K., Camacho, E., Gorske, S. F., Faber, K. T. & Bhattacharya, K. (2026).
  Physics-aware global Rietveld refinement for high-energy X-ray diffraction
  microscopy with application to reconstructing intragranular orientation and
  strain fields. arXiv:2609.13619.
- Kerbl, B., Kopanas, G., Leimkühler, T. & Drettakis, G. (2023). 3D Gaussian
  splatting for real-time radiance field rendering. *ACM Trans. Graph.* 42(4)
  (SIGGRAPH).
- Lin, C.-H., Ma, W.-C., Torralba, A. & Lucey, S. (2021). BARF: Bundle-adjusting
  neural radiance fields. *ICCV*, 5741–5751.
- Lin, Y.-C., Florence, P., Barron, J. T., Rodriguez, A., Isola, P. & Lin, T.-Y.
  (2021). iNeRF: Inverting neural radiance fields for pose estimation. *IROS*.
- Liu, S., Li, T., Chen, W. & Li, H. (2019). Soft Rasterizer: A differentiable
  renderer for image-based 3D reasoning. *ICCV*.
- Mildenhall, B., Srinivasan, P. P., Tancik, M., Barron, J. T., Ramamoorthi, R. &
  Ng, R. (2020). NeRF: Representing scenes as neural radiance fields for view
  synthesis. *ECCV*.
- Tancik, M. et al. (2020). Fourier features let networks learn high frequency
  functions in low dimensional domains. *NeurIPS*.
- Xie, Y. et al. (2022). Neural fields in visual computing and beyond. *Computer
  Graphics Forum* 41, 641–676.
- Zang, G., Idoughi, R., Li, R., Wonka, P. & Heidrich, W. (2021). IntraTomo:
  Self-supervised learning-based tomography via sinogram synthesis and prediction.
  *ICCV*.
- Zha, R., Zhang, Y. & Li, H. (2022). NAF: Neural attenuation fields for
  sparse-view CBCT reconstruction. *MICCAI*.
- Zhong, E. D., Bepler, T., Berger, B. & Davis, J. H. (2021). CryoDRGN:
  reconstruction of heterogeneous cryo-EM structures using neural networks.
  *Nature Methods* 18, 176–185.
