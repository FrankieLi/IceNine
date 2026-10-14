---
title: "Spot Translation Under a Change of Voxel Position"
subtitle: "For a fixed orientation, every diffraction spot moves rigidly by an affine function of the voxel's position"
date: "2026-10-13"
geometry: margin=1in
fontsize: 11pt
---

# Question

A voxel centred at sample-frame position $\mathbf{x}$ with crystal orientation $O$ produces, at
rotation angle $\omega$, a set of lit detector pixels $(u, v)$. Is there a transform that gives the
pixels the same orientation would produce from a voxel at $\mathbf{x}'$?

**Answer.** Yes, per diffraction peak $p$: the spot translates rigidly by
$$
\Delta\mathbf{u}_p \;=\; A_p\,(\mathbf{x}' - \mathbf{x}),
\qquad
A_p \;=\; \Pi\, Q_p\, R_z(\omega_p)\, B \;\in\; \mathbb{R}^{2\times 3},
$$
exactly, in continuous pixel coordinates. The translation differs from peak to peak, so there is no
single transform of a whole frame. The symbols are defined below.

# Notation

All directions are unit vectors. Quantities follow `VoxelCostFunction.evaluate` and
`calculate_diffraction_overlap_batched` in `icenine/cost_functions.py`.

| Symbol | Meaning |
|---|---|
| $O$ | crystal orientation (crystal $\to$ sample frame) |
| $\mathbf{g}_{hkl}$ | reciprocal lattice vector of reflection $p$ in the crystal frame |
| $B,\ \mathbf{t}$ | sample-to-lab rotation and translation (`sample_to_lab_matrix`) |
| $R_z(\omega)$ | rotation by $\omega$ about the stage axis $\hat{\mathbf{z}}$ |
| $\hat{\mathbf{b}}$ | incident beam direction |
| $\mathbf{x}$ | voxel position (centroid or a vertex) in the sample frame |
| $\hat{\mathbf{n}},\ D$ | detector plane $\{\mathbf{X} : \hat{\mathbf{n}}\cdot\mathbf{X} + D = 0\}$ |
| $\hat{\mathbf{e}}_c,\ \hat{\mathbf{e}}_r$ | lab directions of the detector's column and row axes |
| $\mathbf{o}$ | lab position of the detector's pixel origin |
| $a_c,\ a_r$ | pixel width and height |

# Step 1: everything that sets a peak's direction is independent of position

The Bragg condition is solved for $\mathbf{g}_s = O\,\mathbf{g}_{hkl}$ and the beam alone
(`get_scattering_omegas_torch`), giving the two rotation angles $\omega_p$ at which the reflection
fires. The lab-frame plane normal and the diffracted direction at that angle are
$$
\hat{\mathbf{g}}_p = R_z(\omega_p)\,B\,\frac{\mathbf{g}_s}{|\mathbf{g}_s|},
\qquad
\hat{\mathbf{k}}'_p = \hat{\mathbf{b}} - 2\,(\hat{\mathbf{b}}\cdot\hat{\mathbf{g}}_p)\,\hat{\mathbf{g}}_p .
$$
None of $\omega_p$, $\hat{\mathbf{k}}'_p$, the frame index
$f_p = \lfloor(\omega_p - \omega_{\text{low}})/\Delta\omega\rfloor$, the eta angle or the
observability of the reflection involves $\mathbf{x}$.

# Step 2: the voxel enters only as the ray's starting point

At $\omega_p$ the voxel point $\mathbf{x}$ sits at the lab position
$$
\mathbf{x}_{\text{lab}} = R_z(\omega_p)\,B\,\mathbf{x} + \mathbf{t},
$$
and the diffracted ray is $\mathbf{X}(\lambda) = \mathbf{x}_{\text{lab}} + \lambda\,\hat{\mathbf{k}}'_p$.

# Step 3: the detector hit is a parallel projection, hence affine

Solving $\hat{\mathbf{n}}\cdot\mathbf{X}(\lambda) + D = 0$ gives
$$
\lambda^\ast = -\,\frac{\hat{\mathbf{n}}\cdot\mathbf{x}_{\text{lab}} + D}{\hat{\mathbf{n}}\cdot\hat{\mathbf{k}}'_p},
\qquad
\mathbf{s}_p = \mathbf{x}_{\text{lab}} + \lambda^\ast\,\hat{\mathbf{k}}'_p
= Q_p\,\mathbf{x}_{\text{lab}} - \frac{D}{\hat{\mathbf{n}}\cdot\hat{\mathbf{k}}'_p}\,\hat{\mathbf{k}}'_p ,
$$
with the oblique projector (projection onto the detector plane *along the ray*)
$$
Q_p = I - \frac{\hat{\mathbf{k}}'_p\,\hat{\mathbf{n}}^{\top}}{\hat{\mathbf{n}}\cdot\hat{\mathbf{k}}'_p} .
$$
$Q_p$ depends only on $\hat{\mathbf{k}}'_p$ and $\hat{\mathbf{n}}$, so $\mathbf{s}_p$ is an affine function of
$\mathbf{x}_{\text{lab}}$.

# Step 4: pixel coordinates are affine in the hit point

The continuous pixel coordinates (before flooring) are
$$
\mathbf{u}_p = \Pi\,(\mathbf{s}_p - \mathbf{o}) + \tfrac12\,\mathbf{1},
\qquad
\Pi = \begin{pmatrix} \hat{\mathbf{e}}_c^{\top}/a_c \\[2pt] \hat{\mathbf{e}}_r^{\top}/a_r \end{pmatrix}
\in \mathbb{R}^{2\times 3},
$$
and the pixel index is $\lfloor\mathbf{u}_p\rfloor$ (with the C++ rule that a negative millimetre
coordinate is pixel $-1$).

# Step 5: compose

Substituting Steps 2 to 4,
$$
\boxed{\;
\mathbf{u}_p(\mathbf{x}) = A_p\,\mathbf{x} + \mathbf{c}_p,
\qquad
A_p = \Pi\,Q_p\,R_z(\omega_p)\,B,
\qquad
\mathbf{c}_p = \Pi\Big(Q_p\,\mathbf{t} - \frac{D}{\hat{\mathbf{n}}\cdot\hat{\mathbf{k}}'_p}\,\hat{\mathbf{k}}'_p - \mathbf{o}\Big) + \tfrac12\,\mathbf{1}
\;}
$$
Both $A_p$ and $\mathbf{c}_p$ depend on the orientation and the peak, never on the position. Therefore,
for any displacement $\Delta\mathbf{x} = \mathbf{x}' - \mathbf{x}$ (not only a small one),
$$
\mathbf{u}_p(\mathbf{x}') - \mathbf{u}_p(\mathbf{x}) = A_p\,\Delta\mathbf{x} .
$$
Applied to the three vertices $\mathbf{x} + \boldsymbol{\delta}_i$ of a triangular voxel,
$\mathbf{u}_p(\mathbf{x} + \boldsymbol{\delta}_i) = \mathbf{u}_p(\mathbf{x}) + A_p\boldsymbol{\delta}_i$: the
spot's shape $\{A_p\boldsymbol{\delta}_i\}$ is fixed by the orientation and the voxel size, and moving
the voxel translates that shape rigidly.

This is the exact, finite-displacement form of the parallax term of the spot-motion Jacobian in
`docs/nn_inverse_problem_formulation.md` §3.6.3,
$\Gamma_p\hat{\mathbf{z}} = -\Pi\,Q_p\,(\hat{\mathbf{z}}\times\mathbf{x}_{v,p})$.

# Geometric reading

Write the lab displacement $\mathbf{d} = R_z(\omega_p)\,B\,\Delta\mathbf{x}$. For a detector facing
the beam ($\hat{\mathbf{n}} \approx \hat{\mathbf{x}}$) write
$\hat{\mathbf{k}}'_p = (\cos 2\theta_p,\ \sin 2\theta_p\,\hat{\boldsymbol{\rho}}_p)$, where
$\hat{\boldsymbol{\rho}}_p$ is the spot's radial direction in the detector plane. Then
$$
Q_p\,\mathbf{d} \;=\; \mathbf{d}_\perp \;-\; \tan 2\theta_p\; d_x\;\hat{\boldsymbol{\rho}}_p ,
$$
where $\mathbf{d}_\perp$ is the part of $\mathbf{d}$ transverse to the beam and $d_x$ the part along it.

- Motion across the beam moves the spot one-for-one: $\mathbf{d}_\perp / a$ pixels.
- Motion along the beam slides the spot radially by $\tan 2\theta_p\, d_x / a$ pixels (parallax).
- For an in-layer move ($\Delta x_3 = 0$, $B \approx I$), $\mathbf{d}$ is horizontal, so the vertical
  shift comes only from parallax:
  $\Delta u_{\text{row}} \approx -\tan 2\theta_p\,\rho_{p,z}\, d_x / a_r$ (taking $\hat{\mathbf{e}}_r \approx \hat{\mathbf{z}}$), with
  $d_x = \Delta x_1\cos\omega_p - \Delta x_2\sin\omega_p$.

# Several detectors

$A_p = \Pi\,Q_p\,R_z(\omega_p)\,B$ contains the detector only through $\hat{\mathbf{n}}$,
$\hat{\mathbf{e}}_c$, $\hat{\mathbf{e}}_r$ and the pixel size, not through the distance $D$. Detectors
with the same normal, axes and pixel size therefore see the **same** shift in pixels at every distance.
The calibrated IceNine detectors are close to parallel, so their shifts are close but not equal.

# Pixelisation

The translation law is exact for continuous coordinates. Lit pixels come from flooring and
rasterising, so split the shift into an integer and a fractional part,
$$
\Delta\mathbf{u}_p = \mathbf{m}_p + \boldsymbol{\varphi}_p,
\qquad
\mathbf{m}_p = \lfloor\Delta\mathbf{u}_p\rfloor \in \mathbb{Z}^2,\quad
\boldsymbol{\varphi}_p \in [0, 1)^2 .
$$
Since $\lfloor u + m\rfloor = \lfloor u\rfloor + m$ for integer $m$, the rasterised pixel set $\mathcal{R}$
of a triangle $T$ satisfies
$$
\mathcal{R}(T + \Delta\mathbf{u}_p) = \mathcal{R}(T + \boldsymbol{\varphi}_p) + \mathbf{m}_p ,
$$
so only the sub-pixel phase $\boldsymbol{\varphi}_p$ needs re-rasterising (in real arithmetic; float32
rounding can flip a pixel lying within about $10^{-3}$ px of a boundary).

What can change besides the pixel set:

- A shifted spot can leave the detector, which changes the "detector lit" flags.
- $\lambda^\ast > 0$ (the ray points towards the detector) cannot change, since the whole sample lies on the beam side of the detector plane.

# Numerical check

`VoxelCostFunction.evaluate` (final cost, triangle mode, Q-max 8) on Example2.ThreeVoxels voxel 0,
20 orientations within 3 degrees of the truth, two random in-layer moves $\mathbf{a}$, $\mathbf{b}$ of
up to 150 $\mu$m each, projected vertices read at the input of the stage-D C call
(`scripts/speed_parity/translation_check.py`):

| Check | Largest deviation (px) |
|---|---|
| Rigid: all three vertices of a peak shift equally | $9.8\times10^{-4}$ |
| Linear: shift($2\mathbf{a}$) = 2 shift($\mathbf{a}$) | $1.2\times10^{-3}$ |
| Additive: shift($\mathbf{a}+\mathbf{b}$) = shift($\mathbf{a}$) + shift($\mathbf{b}$) | $1.2\times10^{-3}$ |
| Detector 0 shift = detector 1 shift | $0.37$ (of shifts up to 54 px; the detectors are not exactly parallel) |

Frame indices were unchanged in every case. The deviations of the first three rows are float32
rounding of coordinates near 1000 px.

# Consequence for the cost of an evaluation

The orientation-dependent part of an evaluation (Bragg solve, eta filter, frame lookup, $A_p$ and
$\mathbf{c}_p$) does not depend on the voxel. Per voxel, projecting a peak costs one $2\times3$
matrix-vector product, $\mathbf{u}_p = A_p\mathbf{x} + \mathbf{c}_p$, followed by the pixel test.

- The level-0 orientation set (4,886 FZ points $\times$ 9 local-grid points = 43,974 orientations;
  the local grid is generated deterministically, `reconstructor.py:513`) is identical for every
  voxel. A table of $(f_p, A_p, \mathbf{c}_p)$ per orientation and peak can be built once per sample
  (up to about 50 peaks per orientation at Q-max 5, so one to two million peaks; with two detectors,
  of order 100-250 MB, to be measured). Plan: `docs/todo_translation_screen_prototype.md`; after
  `BatchedScreen` the level-0 screen is only about 1% of a no-start search, so the payoff depends on
  ROADMAP S2 and S3.
- A BFS neighbour's first evaluation (the inherited orientation at a new position) is a translation
  of the parent's spots.
- The quick MC, FindOptimal and the higher levels propose new orientations per voxel and gain
  nothing from this; they need a faster `evaluate`.

The speed-up of a table-driven level-0 screen has not been measured.
