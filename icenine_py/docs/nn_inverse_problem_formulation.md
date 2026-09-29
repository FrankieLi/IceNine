---
title: "Inverting an Integrated Forward Model"
subtitle: "What the IceNine orientation network should approximate, and why"
date: "2026-09-28"
geometry: margin=1in
fontsize: 11pt
---

# Summary

The toy orientation network maps binned detector data to a voxel's orientation.
This note makes precise what that map should be.

1. The experiment integrates the diffraction signal over rotation frames and
   detector pixels. The resulting forward operator is many-to-one and
   discontinuous, so it has no inverse in Hadamard's sense (Section 2).
2. For each diffraction peak, the orientations at which the peak changes frame
   form smooth surfaces in orientation space. To first order these are parallel,
   equally spaced planes, and one direction (rotation about the rotation axis)
   behaves exactly linearly. The orientations consistent with a measurement form
   a cell bounded by these surfaces, and its size sets an error floor no method
   can beat. The plane picture holds well for most peaks but fails for the
   near-axis peaks, which are the most informative (Section 3). The same
   geometry gives closed-form angular-resolution estimates in terms of pixel
   size, frame width, number of frames and the voxel's distance from the
   rotation axis. Perpendicular to the axis the pixel size sets the resolution;
   about the axis the frame width sets it for voxels near the axis, and parallax
   (pixels again) for voxels far from it (Section 3.6).
3. The principled replacement for "the inverse" is the Bayesian posterior over
   orientation given the data. In the noise-free limit it is the prior
   restricted to that cell (Section 4).
4. A network trained with squared loss on simulated pairs converges to the
   posterior mean; with a likelihood loss it can also learn the posterior
   covariance. This is amortised, simulation-based inference (Section 5).
5. The simulator contains information real data will not (exact intensities,
   all-or-nothing pixel fill), which the network could exploit (Section 6).

Physics results used here are derived and numerically verified in
`omega_peak_width_derivation.md` (the "derivation note"). The frame-boundary
geometry in Section 3.1 and the spot-motion and resolution results in Section 3.6
have been checked numerically (Sections 3.1.6, 3.6.3 and 3.6.5). Section
numbers are written into the headings, so render without pandoc's `-N`, e.g.
`pandoc nn_inverse_problem_formulation.md -o nn_inverse_problem_formulation.pdf`. The pixel level-set
structure in Section 3.2 has not.

# 1. Setup and notation

## 1.1 Geometry and running example

- **Lab frame.** The incident X-ray beam travels along $+\hat{\mathbf{x}}$ with
  wavevector $\mathbf{k}$, $|\mathbf{k}| = 2\pi/\lambda$. The sample rotates about
  $\hat{\mathbf{z}}$, the **rotation axis**.
- **Rotation angle.** $\omega$ is the rotation-stage angle and
  $R_z(\omega)$ the rotation by $\omega$ about $\hat{\mathbf{z}}$.
- **Sample frame.** Coordinates fixed to the sample. It coincides with the lab
  frame at $\omega = 0$ (true for the running example, whose base sample
  rotation is the identity), and a sample-frame vector $\mathbf{v}$ appears in the
  lab at $R_z(\omega)\mathbf{v}$.
- **Frames.** The detector records one image, a **frame**, per rotation
  interval of width $\Delta\omega_f$. Frame $j$ is centred on $\omega_j$ and
  integrates over $I_j = [\omega_j - \Delta\omega_f/2,\ \omega_j + \Delta\omega_f/2]$.
  Its **frame edges** are the endpoints of these intervals.
- **Pixels.** Detector position $\mathbf{u} = (x, y)$ is measured in pixel
  units (column, row); pixel $m$ covers the unit square $P_m$.
- **Voxel.** A small triangular element of the sample, assumed to be a single
  perfect crystal.
- **Running example.** `Example2.ThreeVoxels`: copper, 64.351 keV, two
  detectors, 180 frames with $\Delta\omega_f = 1^\circ$, 1.48 $\mu$m pixels.

## 1.2 The unknown

Local refinement starts from a known **nominal orientation** $R_{\text{nom}}$,
for example the output of IceNine's coarse Sukharev orientation-grid search. The
unknown is a small rotation vector $\delta \in \mathbb{R}^3$ applied in the sample
frame:
$$
R(\delta) = \exp([\delta]_\times)\,R_{\text{nom}},
$$
where $[\delta]_\times$ is the $3\times3$ skew-symmetric matrix with
$[\delta]_\times\mathbf{v} = \delta\times\mathbf{v}$, and $\exp$ is the matrix
exponential. Equivalently, $\delta$ is a rotation by angle $|\delta|$ about the
axis $\delta/|\delta|$ (an element of the Lie algebra $\mathfrak{so}(3)$). This
is how the toy-NN prototype composes perturbations (a left multiplication). The
prototype draws its perturbations uniformly in a box of the near-identity
quaternion parameterisation used by IceNine's Monte Carlo search, not uniformly
in $\delta$, so the prior it implies is not exactly uniform in $\delta$. The true offset
is $\delta_{\text{true}}$. The training perturbation distribution is a prior
$\pi(\delta)$. $\nabla$ denotes the gradient with respect to $\delta$, and
$\mathbf{e}_1, \mathbf{e}_2, \mathbf{e}_3$ the standard basis of $\mathbb{R}^3$.

## 1.3 Peaks

- **Reflection.** Each crystal reflection has a reciprocal-lattice vector
  $\mathbf{g}_{hkl}$ in the crystal frame; in the sample frame it is
  $\mathbf{g}_s = R_{\text{nom}}\,\mathbf{g}_{hkl}$ at the nominal orientation,
  and $\exp([\delta]_\times)\,\mathbf{g}_s$ after the offset. Its Bragg angle $\theta$ is set
  by $|\mathbf{g}| = 2|\mathbf{k}|\sin\theta$.
- **Branch and peak.** As the sample rotates, a reflection satisfies the Bragg
  condition at up to two angles (two **branches**). A (reflection, branch) pair
  is **Bragg-observable** if that angle lies in the scanned range. A **peak** $p$
  is a Bragg-observable pair whose spot also lands on a detector; the prototype
  code calls the set of peaks of one voxel its **ROI set**. There are $P$ peaks.
- **Bragg crossing.** $\omega^*_p(\delta)$ is the rotation angle at which peak
  $p$ satisfies the Bragg condition.
- **Angles.** $\chi_p$ is the angle between $\mathbf{g}_s$ and the rotation axis,
  and $\varphi_{0,p}$ the azimuth of $\mathbf{g}_s$ about $\hat{\mathbf{z}}$,
  measured from $\hat{\mathbf{x}}$. $\eta_p$ is the azimuth of the diffracted
  beam around the incident beam, measured from the rotation-axis direction; it
  satisfies $\cos\eta_p = \cos\chi_p/\cos\theta_p$ (derivation note, §5.3).
  "Near-axis" peaks are those with small $|\sin\eta_p|$, i.e. $\chi_p$ close to
  $\theta_p$.
- **Detector position.** $\mathbf{u}_p(\delta) = (x_p, y_p)$ is where peak $p$
  lands, in pixel units. The spot's footprint is the projected voxel triangle,
  with vertices $\mathbf{u}_{p,i} = (x_{p,i}, y_{p,i})$, $i = 0, 1, 2$.
- **Intensity.** $L_p(\delta)$ is the peak's integrated intensity, including the
  **Lorentz factor** $1/(\sin 2\theta_p\,|\sin\eta_p|)$ applied by
  `XDMEtaAcceptFn` and its batched form `batch_eta_filter`. Its
  $1/|\sin\eta_p|$ part accounts for how long the reflection stays in the
  diffracting condition (derivation note, §6.1); the $1/\sin 2\theta_p$ part is
  not derived in these notes.
- **Rocking width.** Because of the beam's energy bandwidth and divergence, the
  Bragg condition is satisfied over a small range $\alpha$ of glancing angles
  (the **tolerance**; there is no mosaic contribution under the
  perfect-crystal assumption). Peak $p$ is then excited over a range of
  rotation angle of width $\Delta\omega_p \approx \alpha/|\sin\eta_p|$
  (derivation note, §5). $\alpha = 0$ is the **single-frame case**: each peak is
  excited at a single instant and lands in exactly one frame, which is what the
  current simulator assumes.

## 1.4 Instantaneous signal

The signal from peak $p$ arriving at rotation angle $\omega$ and detector point
$\mathbf{u}$ is modelled as
$$
\rho_p(\omega, \mathbf{u};\delta) = L_p(\delta)\;\Psi_p\big(\omega - \omega^*_p(\delta)\big)\;s_p\big(\mathbf{u} - \mathbf{u}_p(\delta)\big),
$$
where $\Psi_p$ is the **rocking profile** (unit area, width $\Delta\omega_p$) and
$s_p$ is the spot footprint.

## 1.5 Integrated forward operator

The detector records, for each peak, frame and pixel,
$$
D_p[j, m] = \int_{I_j}\!\int_{P_m} \rho_p(\omega, \mathbf{u};\delta)\,d\mathbf{u}\,d\omega .
$$
Collecting all entries gives the data $D = A(\delta)$; $A$ is the **integrated
forward operator**. It helps to split $A = B \circ G$:

- $G$: $\delta \mapsto \{\omega^*_p, \mathbf{u}_p, L_p\}_{p=1}^{P}$ is smooth
  (rigid-body geometry and the Bragg condition).
- $B$: integration over frames and pixels, i.e. binning.

All of the difficulty below comes from $B$.

## 1.6 Observations and thresholding

Real data are $D_{\text{obs}} = A(\delta_{\text{true}}) + \varepsilon$, with
$\varepsilon$ covering noise and model error. IceNine's reconstruction then
binarises: experimental frames are thresholded (`get_binary_numpy`,
`SparseImageStack(binary=True)`), and the hard cost counts overlapping lit pixels
(`cost_functions.py`, described there as "inherently non-differentiable (binary
pixel test, integer counting)").

Write $\mathbf{1}[\cdot]$ for the indicator (1 if the condition holds, else 0),
$\tau$ for the detection threshold, and $T(D) = \mathbf{1}[D > \tau]$
(elementwise). The **measurement map** is
$$
M = T \circ A ,
$$
from orientation offset to thresholded data. **In Sections 2–5, "data" and $D$
mean thresholded data, $D = M(\delta)$.** Unthresholded data are discussed in
Section 6.

# 2. Why the integrated operator has no inverse

Hadamard (1902) called a problem *well-posed* if a solution exists, is unique,
and depends continuously on the data. Problems that violate any of these are
*ill-posed*; see Engl, Hanke & Neubauer (1996) for the classical theory.
Inverting $M$ fails all three conditions.

## 2.1 Uniqueness fails: the measurement map is piecewise constant

Take the single-frame case ($\alpha = 0$). Peak $p$ lands in frame $j$ exactly
when $\omega^*_p(\delta) \in I_j$, which is the indicator of a region of
$\delta$-space. Likewise, whether a pixel is lit is the indicator of the region
of $\delta$ where the spot covers it (the simulator's fill is all-or-nothing per
pixel). So $M$ is **constant on the cells of a partition of $\delta$-space** and
jumps between cells.

Consequently, each **fibre** $M^{-1}(D) = \{\delta : M(\delta) = D\}$ has positive
volume: a whole region of orientations produces identical data. No function of
$D$ can recover $\delta_{\text{true}}$ exactly.

(The unthresholded $A$ is not piecewise constant, because $L_p(\delta)$ varies
smoothly. That residual dependence is the intensity leak discussed in
Section 6.)

## 2.2 Stability fails: jumps, and zero gradient almost everywhere

At cell boundaries $M$ is discontinuous, so an arbitrarily small change in
orientation can change the data by a whole frame or pixel. Inside a cell its
derivative is zero. This is why IceNine's hard cost cannot be optimised by
gradient descent, and why the differentiable cost blurs the data first
(`MultiScaleImageStack` max-pooling and omega blending). That blurring is a
mollification of the forward model, a standard regularisation device, rather
than a change to the underlying measurement.

## 2.3 Existence fails for real data

With noise and model error (for example, the simulator placing each peak in a
single frame when real peaks spread), $T(D_{\text{obs}})$ is generally not in
the range of $M$, so $M(\delta) = T(D_{\text{obs}})$ has no exact solution at
all.

## 2.4 How this differs from the textbook case

Classical ill-posedness usually comes from a compact integral operator whose
smoothing makes the inverse unbounded (Engl et al. 1996). Here the unknown is
only three-dimensional, and the problem is not ill-conditioning but
**quantisation**: sampling on a grid of frames and pixels. The relevant
mathematics is that of quantisation and bounded-error estimation (Sections 3–4)
as much as regularisation.

# 3. Geometry of the consistent set

## 3.1 Frame boundaries: level sets that are nearly parallel planes

This subsection states what is exact and what is first-order, and reports a
numerical check. Throughout it, $\alpha = 0$.

### 3.1.1 The Bragg condition as a function of $\omega$ and $\delta$

Fix one peak and drop the subscript $p$. With $\mathbf{g}_s$ the reflection's
sample-frame vector at nominal orientation, the Bragg condition including the
perturbation is
$$
F(\omega, \delta) = \mathbf{k}\cdot R_z(\omega)\,\exp([\delta]_\times)\,\mathbf{g}_s + \tfrac12|\mathbf{g}|^2 = 0,
$$
and $\omega^*(\delta)$ is the solution for the peak's branch. Define
$f'(\omega) = \partial F/\partial\omega$. From the derivation note (§4 and §5.3),
$|f'(\omega^*)| = |\mathbf{k}||\mathbf{g}|\cos\theta\,|\sin\eta|$, which is
nonzero whenever $|\sin\eta| > 0$.

### 3.1.2 Exact statement: boundaries are level sets of $\omega^*$

With $\alpha = 0$, peak $p$'s contribution to the data depends on $\delta$ only
through the index $j_p(\delta)$ of the frame containing $\omega^*_p(\delta)$.
The region of $\delta$-space where $j_p = j$ is bounded by the **level sets**
$$
\mathcal{S}_{p,b} = \{\delta : \omega^*_p(\delta) = b\}, \qquad b \in \{\text{frame edges}\}.
$$
Since $f'(\omega^*) \neq 0$, the implicit function theorem makes $\omega^*_p$ a
smooth function of $\delta$ and each $\mathcal{S}_{p,b}$ a smooth surface. So
far nothing is approximate. Surfaces of the same peak never intersect, because
they are level sets of one function at different values.

### 3.1.3 The gradient of $\omega^*$

At $\delta = 0$, with lab-frame vector $\mathbf{g} = R_z(\omega^*)\mathbf{g}_s$:

- $\partial F/\partial\omega = \mathbf{k}\cdot(\hat{\mathbf{z}}\times\mathbf{g}) = f'(\omega^*)$,
  since $\tfrac{d}{d\omega}R_z(\omega)\mathbf{v} = \hat{\mathbf{z}}\times R_z(\omega)\mathbf{v}$;
- $\partial F/\partial\delta_i = \mathbf{k}\cdot R_z(\omega^*)(\mathbf{e}_i\times\mathbf{g}_s)
  = \big(R_z(\omega^*)\mathbf{e}_i\big)\cdot(\mathbf{g}\times\mathbf{k})$, so
  $\nabla_\delta F = R_z(\omega^*)^{\!\top}(\mathbf{g}\times\mathbf{k})$.

By implicit differentiation,
$$
\nabla\omega^*_p = -\frac{R_z(\omega^*_p)^{\!\top}(\mathbf{g}_p\times\mathbf{k})}{f'(\omega^*_p)},
\qquad
|\nabla\omega^*_p| = \frac{|\mathbf{g}_p||\mathbf{k}|\cos\theta_p}{|\mathbf{k}||\mathbf{g}_p|\cos\theta_p\,|\sin\eta_p|} = \frac{1}{|\sin\eta_p|},
$$
using $|\mathbf{g}\times\mathbf{k}| = |\mathbf{g}||\mathbf{k}|\cos\theta$ (the
Bragg condition fixes the angle between $\mathbf{k}$ and $\mathbf{g}$ at
$90^\circ + \theta$).

### 3.1.4 First order: parallel, equally spaced planes

Linearising, $\omega^*_p(\delta) \approx \omega^*_p(0) + \nabla\omega^*_p\cdot\delta$.
The level sets of a linear function are **parallel planes**, all with the unit
normal $\mathbf{n}_p = \nabla\omega^*_p/|\nabla\omega^*_p|$. Frame edges are
equally spaced by $\Delta\omega_f$, so the planes are also **equally spaced**,
by
$$
h_p = \frac{\Delta\omega_f}{|\nabla\omega^*_p|} = \Delta\omega_f\,|\sin\eta_p| .
$$
Each peak therefore confines $\delta$ to a **slab** of thickness $h_p$ with
normal $\mathbf{n}_p$. Different peaks (including the two branches of one
reflection) have different normals, which is what makes the intersection of
slabs a bounded cell (Section 3.3), provided the normals span $\mathbb{R}^3$.

In the running example $|\sin\eta_p| \in [0.033, 1]$, so slab thickness ranges
from $0.033^\circ$ to $1^\circ$. Near-axis peaks give the thinnest slabs.

### 3.1.5 One direction is exact: rotation about the rotation axis

Let $\beta$ be a rotation angle. Rotations about $\hat{\mathbf{z}}$ commute with
$R_z(\omega)$: $R_z(\omega)\exp(\beta[\hat{\mathbf{z}}]_\times) = R_z(\omega + \beta)$.
So rotating *any* orientation by $\beta$ about the rotation axis (in the sample
frame) shifts *every* peak's crossing by exactly $-\beta$, to all orders:
$$
\omega^*_p\big(\exp(\beta[\hat{\mathbf{z}}]_\times)\,R\big) = \omega^*_p(R) - \beta .
$$
Correspondingly, $\nabla\omega^*_p\cdot\hat{\mathbf{z}} = -1$ for every peak,
since $\hat{\mathbf{z}}\cdot(\mathbf{g}\times\mathbf{k}) = \mathbf{k}\cdot(\hat{\mathbf{z}}\times\mathbf{g}) = f'$.
Three consequences:

1. **Every slab is exactly one frame wide along $\hat{\mathbf{z}}$.** Since
   $\mathbf{n}_p\cdot\hat{\mathbf{z}} = -|\sin\eta_p|$, the slab's extent along
   $\hat{\mathbf{z}}$ is $h_p/|\mathbf{n}_p\cdot\hat{\mathbf{z}}| = \Delta\omega_f$.
2. **Pixels see rotation about $\hat{\mathbf{z}}$ only through parallax.** The
   rotated crystal at $\omega^* - \beta$ has the same lab orientation as before,
   so every diffracted ray has the same *direction*. The spot moves only because
   the ray now starts from a different *place*: the voxel's lab position is
   rotated by $-\beta$ about the axis, a displacement of $\beta\,r_\perp$, where
   $r_\perp$ is the voxel's distance from the axis. Section 3.6.2 shows that the
   spot then shifts by about $\beta\,r_\perp|\cos\phi_v|/a$ pixels, where
   $\phi_v$ is the voxel's lab azimuth when the peak fires and $a$ is the pixel
   size. Near the axis this is negligible (running example: $r_\perp = 12\,\mu$m,
   at most 0.14 pixel for a $1^\circ$ rotation), and the $\hat{\mathbf{z}}$
   component of $\delta$ is determined almost entirely by frame membership. Far
   from the axis it is not: beyond roughly $r_\perp \approx \sqrt{2}\,a/\Delta\omega_f$
   (about $120\,\mu$m for the running example's pixels and frames), parallax
   becomes the main source of $\hat{\mathbf{z}}$ information (Section 3.6.4).
3. **Staggered slabs give a Vernier effect along $\hat{\mathbf{z}}$.** Each slab
   is a full frame wide along $\hat{\mathbf{z}}$, but peaks sit at different
   positions within their frames, so the slabs are offset from one another. If
   those offsets are spread uniformly, the intersection of $P$ one-frame-wide
   intervals that all contain the true point has expected length about
   $2\Delta\omega_f/(P+1)$. Resolution along the rotation axis therefore improves
   with the number of peaks, down to the accuracy with which the frame edges are
   known (rotation-stage calibration).

### 3.1.6 Where the planes bend: numerical check

The planes are a first-order picture, and their curvature grows sharply as
$\eta \to 0$. For one peak, the Bragg condition reads
$\cos(\omega^* + \varphi_0) = -\sin\theta/\sin\chi$ (derivation note, §3).
Differentiating gives $|d\omega^*/d\chi| = \tan\theta\cot\chi/|\sin\eta|$, and
differentiating once more gives a second derivative of order
$\tan\theta/|\sin\eta|^3$ for near-axis peaks. Over the ball $|\delta| \le \beta$,
the linearisation error in $\omega^*$ is therefore roughly
$$
e_{\text{lin}} \sim \frac{\beta^2\tan\theta}{|\sin\eta|^3},
$$
and the plane picture holds (error well below a frame) when
$\beta \ll \sqrt{\Delta\omega_f\,|\sin\eta|^3/\tan\theta}$, with angles in
radians. In the running example $\tan\theta$ ranges from 0.05 to 0.25; for
$\tan\theta = 0.1$ the condition gives about $4^\circ$ at $|\sin\eta| = 0.3$ but
only about $0.14^\circ$ at $|\sin\eta| = 0.033$. This is a heuristic scaling,
not a bound.

**Check.** All 5,100 Bragg-observable (reflection, branch) pairs of the running
example's three voxels, whether or not they hit a detector (the geometry of
$\omega^*$ does not depend on the detector).
For each peak and radius $\beta$, 200 random directions on the sphere
$|\delta| = \beta$: the exact $\omega^*$ from `get_scattering_omegas_torch` was
compared with the linear prediction, and the gradient (by central differences)
was re-evaluated at 12 of those points to measure how much the plane normal
tilts and the spacing changes.

- **Gradient formula:** $|\nabla\omega^*|\,|\sin\eta| = 1$ to $6\times10^{-14}$;
  $\nabla\omega^*\cdot\hat{\mathbf{z}} = -1$ exactly; the $\hat{\mathbf{z}}$ shift
  is exact to $4\times10^{-15}$ rad for $\beta$ up to $10^\circ$.
- **Linearisation error** (max over directions, in frames; median / max over
  peaks):

| $\beta$ | $\lvert\sin\eta\rvert < 0.1$ (20 peaks) | $0.1$–$0.3$ (180 peaks) | $\ge 0.3$ (4,900 peaks) |
|---|---|---|---|
| $0.1^\circ$ | 0.14 / 0.70 | 0.003 / 0.018 | 0.0001 / 0.001 |
| $0.5^\circ$ | 4.3 / 5.9 (13% stop diffracting) | 0.07 / 0.53 | 0.002 / 0.036 |
| $1^\circ$ | 12 / 15 (20% stop diffracting) | 0.31 / 2.8 | 0.008 / 0.15 |
| $2^\circ$ | 23 / 36 (27% stop diffracting) | 1.4 / 15 (1% stop) | 0.03 / 0.63 |

- **Parallelism** (largest tilt of the plane normal; largest factor by which the
  local plane spacing changes, up or down; median / max over peaks). Gradients
  are only evaluated where the peak still diffracts, so for $|\sin\eta| < 0.1$ at
  larger $\beta$ the most extreme points drop out; that is why the maximum
  spacing factor falls from ×7.5 at $0.5^\circ$ to ×2.7 at $1^\circ$:

| $\beta$ | $\lvert\sin\eta\rvert < 0.1$ | $0.1$–$0.3$ | $\ge 0.3$ |
|---|---|---|---|
| $0.1^\circ$ | $2.5^\circ$ / $3.0^\circ$; ×1.09 / ×1.39 | $0.45^\circ$ / $1.0^\circ$; ×1.01 / ×1.04 | $0.07^\circ$ / $0.33^\circ$; ×1.001 / ×1.007 |
| $0.5^\circ$ | $8^\circ$ / $18^\circ$; ×1.8 / ×7.5 | $2.3^\circ$ / $5.3^\circ$; ×1.05 / ×1.24 | $0.37^\circ$ / $1.7^\circ$; ×1.005 / ×1.04 |
| $1^\circ$ | $15^\circ$ / $21^\circ$; ×1.9 / ×2.7 | $4.5^\circ$ / $12^\circ$; ×1.11 / ×1.77 | $0.74^\circ$ / $3.4^\circ$; ×1.009 / ×1.08 |

- **Heuristic scaling:** a log–log fit of the linearisation error at
  $\beta = 1^\circ$ against $\beta^2\tan\theta/|\sin\eta|^3$ gives slope 0.84
  (1 would be exact agreement), so the heuristic captures the trend but not the
  exact exponent.

**Reading the results.** For 96% of peaks ($|\sin\eta| \ge 0.3$), the parallel
plane picture is accurate to a few hundredths of a frame out to $0.5^\circ$, with
normals tilting under $2^\circ$ and spacing changing under 4%. For the 0.4% of
peaks with $|\sin\eta| < 0.1$, the boundaries are strongly curved even at
$0.1^\circ$, and beyond about $0.5^\circ$ these peaks often stop diffracting
altogether, which adds a different kind of boundary (the surface $\chi = \theta$
where the peak appears or disappears). These are the same peaks that give the
thinnest slabs, so **the most informative peaks are exactly where the plane
picture fails**. Any exact computation should use the full nonlinear
$\omega^*_p(\delta)$.

## 3.2 Pixel boundaries

In the simulator, pixel boundaries have the same structure. `add_triangle_scanline`
first truncates each of the spot triangle's three vertices to integer pixel
coordinates (`image_data.py`, matching the C++ rasteriser), then fills between
them. The lit-pixel set is therefore a function only of the six integers
$\lfloor x_{p,i}(\delta)\rfloor$, $\lfloor y_{p,i}(\delta)\rfloor$, $i = 0,1,2$.
Its boundaries in $\delta$-space are the level sets $\{x_{p,i}(\delta) = \ell\}$
and $\{y_{p,i}(\delta) = \ell\}$ for integers $\ell$, i.e. six families of
surfaces per spot. To first order each family is a set of parallel planes
spaced one pixel divided by $|\nabla x_{p,i}|$ or $|\nabla y_{p,i}|$.

Three further points:

- **Vernier again.** The three vertices of a small triangle move almost
  together, so the three column families are nearly parallel to one another
  but offset, and likewise the rows.
- **Spots are about a pixel wide.** The running example's voxel sides (1.5 and
  0.75 $\mu$m) are comparable to its 1.48 $\mu$m pixels, so each spot lights only
  one or a few pixels.
- **Pixel planes carry $\hat{\mathbf{z}}$ information only through parallax**
  (Section 3.1.5, point 2): negligible for voxels near the rotation axis,
  important far from it.

The spatial sensitivity is derived in Section 3.6.2 and given exactly, with a
numerical check, as the spot-motion Jacobian in Section 3.6.3. For a voxel on
the rotation axis, a strain-free spot can only slide along its instantaneous
Debye–Scherrer ring, so to first order all six families share one normal
direction, perpendicular to $\hat{\mathbf{z}}$. For an off-axis voxel, parallax
adds a second direction. The level-set structure itself has not been checked
numerically. Real detectors integrate partial
pixel coverage and a point-spread function, so real pixel boundaries are soft,
like the frame boundaries of spread peaks (Section 3.5), whereas the simulator's
are hard.

## 3.3 The consistent set

A **hyperplane arrangement** is a finite collection of planes; its **cells** are
the connected regions left when the planes are removed (Stanley 2007). In the
regime where the plane picture holds, the frame and pixel planes of all peaks
form such an arrangement, and its cells are the fibres of Section 2.1. Given
data $D$, let $j_p^{\text{obs}}$ be the observed frame of peak $p$. The
**consistent set** is (in this linear regime) the convex polytope
$$
C(D) = \bigcap_{p=1}^{P} \{\delta : \omega^*_p(\delta) \in I_{j_p^{\text{obs}}} \text{ and peak } p \text{ lights its observed pixels}\},
$$
an intersection of slabs. Computing it is **set-membership (bounded-error)
estimation**: each measurement confines the parameter to a set, and the
estimate is their intersection (Milanese & Vicino 1991).

With one peak, $C$ is a single slab. With several peaks whose normals
$\mathbf{n}_p$ span $\mathbb{R}^3$, the slabs intersect in a small polytope; more
peaks with diverse normals give smaller cells.

Outside that regime (near-axis peaks, larger perturbations) $C(D)$ is defined by
the same conditions, but its boundaries are curved surfaces, so it need not be
convex.

## 3.4 One-dimensional intuition: quantisation error

If a single measurement only says "$\delta\cdot\mathbf{n}_p$ lies in an interval
of width $h$" and the prior is locally flat, the posterior along $\mathbf{n}_p$
is uniform on that interval, with standard deviation $h/\sqrt{12}$. This is the
classical quantisation-noise result (Widrow & Kollár 2008; Gray & Neuhoff 1998).
In the running example that gives a root-mean-square (RMS) error of roughly
$0.01^\circ$ to $0.29^\circ$ per peak, along that peak's normal, from frame
information alone.

## 3.5 Spread peaks carry sub-bin information

When $\alpha > 0$ and a peak's rocking width $\Delta\omega_p$ is comparable to
$\Delta\omega_f$, the fraction of its intensity in each frame varies
*continuously* with $\omega^*_p(\delta)$. Along $\mathbf{n}_p$ the data are then
no longer piecewise constant, and the slab collapses toward a plane: the peak
reports $\omega^*_p$ more finely than one frame. This is the same mechanism as
**dither** in quantisation, where spreading a signal before quantising lets
averaged outputs resolve below one step (Widrow & Kollár 2008), and as
**centroid localisation** of a blurred spot, which reaches precision far below
the pixel or resolution limit (Thompson, Larson & Webb 2002).

The same holds spatially. Real detectors record partial pixel coverage and a
point-spread function, so real data carry sub-pixel information. The
simulator's all-or-nothing fill does not.

## 3.6 Angular resolution from pixel size, frame width and number of frames

This section turns the geometry above into an estimate of how precisely the
data determine $\delta$, as a function of the pixel size, the frame width and
the number of frames. The main result is that the resolution is anisotropic.
For a voxel near the rotation axis, the frame width sets it about the axis and
the pixel size sets it perpendicular to the axis. Farther from the axis,
parallax lets the pixels constrain rotation about the axis too, and the
anisotropy shrinks.

### 3.6.1 Assumptions and new notation

- Single-frame case ($\alpha = 0$) and a strain-free crystal (lattice spacings
  fixed).
- Each detector is flat, perpendicular to the beam, at distance $d$ from the
  rotation axis, with square pixels of side $a$ (both in the same length unit).
  $a/d$ is then roughly the angle one pixel subtends at the sample. Each peak is
  measured on one detector.
- The voxel sits at distance $r_\perp$ from the rotation axis. When peak $p$
  fires, the voxel's lab position is $\mathbf{x}_{v,p} = R_z(\omega^*_p)\mathbf{x}_v$
  (with $\mathbf{x}_v$ its sample-frame position), and $\phi_{v,p}$ is its
  azimuth about $\hat{\mathbf{z}}$ measured from $\hat{\mathbf{x}}$. Results are
  first order in $r_\perp/d$ and assume $2\theta$ is small enough that rays hit
  the detector nearly head-on.
- The scan has $N$ frames of width $\Delta\omega_f$. The number of peaks is
  $P = \nu\,N\,\Delta\omega_f$, where $\nu$ is the number of observable peaks per
  radian of rotation (it depends on the crystal, beam energy, detector coverage,
  and $Q_{\max}$, the largest $|\mathbf{g}|$ included).
- Quantisation is treated as independent measurement error: a frame index has
  variance $\Delta\omega_f^2/12$ and a pixel coordinate $1/12$ px$^2$ (the
  variance of a uniform error over one bin, Section 3.4), independent across
  peaks and coordinates. The three vertices of a spot are not counted
  separately.

For peak $p$, work in the lab frame at the moment it diffracts. Let
$\hat{\mathbf{g}}_p = \mathbf{g}_p/|\mathbf{g}_p|$ with $\mathbf{g}_p = R_z(\omega^*_p)\mathbf{g}_{s,p}$,
$$
\hat{\mathbf{c}}_p = \frac{\mathbf{g}_p\times\mathbf{k}}{|\mathbf{g}_p\times\mathbf{k}|},
\qquad
\hat{\mathbf{t}}_p = \hat{\mathbf{g}}_p\times\hat{\mathbf{c}}_p ,
$$
so $(\hat{\mathbf{g}}_p, \hat{\mathbf{c}}_p, \hat{\mathbf{t}}_p)$ is an orthonormal
basis, and write $\delta_L = R_z(\omega^*_p)\,\delta$ for the orientation offset
expressed in that lab frame. Because $R_z$ preserves $\hat{\mathbf{z}}$ and
lengths, statements about the $\hat{\mathbf{z}}$ component and about magnitudes
carry over to $\delta$ unchanged. $\nabla_L$ is the gradient with respect to
$\delta_L$.

### 3.6.2 What one peak measures

**Blind direction.** Rotating the crystal about $\hat{\mathbf{g}}_p$ leaves
$\mathbf{g}_p$, and hence the whole peak, unchanged. Each peak measures at most
two components of $\delta$.

**Frame.** From Section 3.1.3, in the lab frame,
$\nabla_L\omega^*_p = -\hat{\mathbf{c}}_p/(\hat{\mathbf{z}}\cdot\hat{\mathbf{c}}_p)$,
with $|\hat{\mathbf{z}}\cdot\hat{\mathbf{c}}_p| = |\sin\eta_p|$.

**What "the ring" means here.** With no strain, $|\mathbf{g}_p|$ is fixed, and
elastic scattering fixes $|\mathbf{k}'| = |\mathbf{k}|$ for the diffracted
wavevector $\mathbf{k}' = \mathbf{k} + \mathbf{g}$. So the *direction* of the
diffracted ray always lies on the cone of half-angle $2\theta_p$ about the beam
(the Debye–Scherrer cone). The ray starts at the voxel, so at the instant peak
$p$ fires, the possible spot positions form the intersection of the detector
plane with that cone placed with its apex at $\mathbf{x}_{v,p}$: for a detector
perpendicular to the beam, a circle of radius $d_p\tan 2\theta_p$ centred on the
voxel's projection along the beam, where
$d_p = d - \mathbf{x}_{v,p}\cdot\hat{\mathbf{x}}$ is the voxel-to-detector
distance. This **instantaneous ring** is a conic section.

The ring seen in data accumulated over a scan is not. Different peaks of the
same reflection family fire at different $\omega^*$, when the voxel is at a
different lab position, so their instantaneous circles have different centres
(offset by up to $r_\perp$) and radii. The accumulated locus is a union of
shifted circles rather than a single conic section. The shift is negligible when
$r_\perp \ll d$ and pixels are coarse (far field), but not in near-field
geometry, where $r_\perp$ can be hundreds of micrometres and pixels are about
$1.5\,\mu$m. The derivation below uses only the instantaneous ring, locally
around one spot, and never assumes a global ring.

A small orientation offset moves a spot in two ways: the ray's direction
changes, sliding the spot along its instantaneous ring, and the ray's starting
point changes, because $\omega^*$ shifts and the voxel moves with the stage
(parallax).

**Position, part 1: the direction slides along the instantaneous ring.**
After a small offset, the lab-frame vector
changes by $d\mathbf{g} = \Omega\times\mathbf{g}$ with
$\Omega = \delta_L + d\omega^*\,\hat{\mathbf{z}}$: the crystal rotation plus the
extra stage rotation needed to stay in the Bragg condition. Staying in the
condition ($\mathbf{k}\cdot d\mathbf{g} = 0$) means
$\Omega\cdot(\mathbf{g}\times\mathbf{k}) = 0$, i.e. $\Omega$ has no
$\hat{\mathbf{c}}$ component. Writing
$\Omega = \Omega_g\hat{\mathbf{g}} + \Omega_t\hat{\mathbf{t}}$ and using
$\hat{\mathbf{t}}\times\hat{\mathbf{g}} = \hat{\mathbf{c}}$,
$$
d\mathbf{k}' = d\mathbf{g} = |\mathbf{g}|\,\Omega_t\,\hat{\mathbf{c}} .
$$
On the cone, $|d\mathbf{k}'| = |\mathbf{k}|\sin 2\theta\,|d\eta|$, so
$|d\eta| = 2\sin\theta\,|\Omega_t|/\sin 2\theta = |\Omega_t|/\cos\theta$. The
instantaneous ring has radius $d_p\tan 2\theta$, so the spot moves along it by
$d_p\tan 2\theta\,|d\eta|$. Let $\xi_p$ be the spot's position along its
instantaneous ring, in pixels. Since
$\Omega_t = \delta_L\cdot\hat{\mathbf{t}} + d\omega^*\,(\hat{\mathbf{z}}\cdot\hat{\mathbf{t}})$
and $d\omega^* = \nabla_L\omega^*\cdot\delta_L$,
$$
\nabla_L\xi_p = \frac{d_p\,\tan 2\theta_p}{a\cos\theta_p}\,
\Big[\hat{\mathbf{t}}_p + (\hat{\mathbf{z}}\cdot\hat{\mathbf{t}}_p)\,\nabla_L\omega^*_p\Big].
$$
Because $\hat{\mathbf{z}}\cdot\nabla_L\omega^* = -1$,
$\hat{\mathbf{z}}\cdot\nabla_L\xi_p = \hat{\mathbf{z}}\cdot\hat{\mathbf{t}}_p - \hat{\mathbf{z}}\cdot\hat{\mathbf{t}}_p = 0$:
the direction term is exactly blind to rotation about the axis.

**Position, part 2: parallax.** When $\omega^*$ shifts by $d\omega^*$, the voxel
moves with the stage by $d\omega^*\,(\hat{\mathbf{z}}\times\mathbf{x}_{v,p})
= r_\perp\,d\omega^*\,(-\sin\phi_{v,p},\ \cos\phi_{v,p},\ 0)$, a horizontal
displacement. For a ray hitting the detector nearly head-on, the spot moves by
this displacement's component in the detector plane, which is along the
detector's horizontal axis $\hat{\mathbf{e}}_h$ with magnitude
$r_\perp\cos\phi_{v,p}\,d\omega^*$. In pixels,
$$
\Delta\mathbf{u}_p^{\text{par}} \approx \frac{r_\perp\cos\phi_{v,p}}{a}\,(\nabla_L\omega^*_p\cdot\delta_L)\;\hat{\mathbf{e}}_h .
$$

These two parts are the intuition. Section 3.6.3 assembles them into a single
matrix, exact to first order in $\delta$ for any flat detector, without the
head-on approximation used for part 2.

### 3.6.3 The spot-motion Jacobian

**Definition.** The spot-motion Jacobian $\Gamma_p$ of peak $p$ is the
$2\times3$ matrix giving how far its spot moves on the detector, in pixels, for
a small orientation offset:
$$
\Delta\mathbf{u}_p \approx \Gamma_p\,\delta .
$$
Its rows are the detector column and row; its columns are rotation about the
sample $\hat{\mathbf{x}}$, $\hat{\mathbf{y}}$ and $\hat{\mathbf{z}}$ axes; its
entries are in pixels per radian. "Full" means it includes both ways a spot
moves: the change of the diffracted ray's direction and the change of the ray's
starting point (parallax).

**Ingredients.** All evaluated at the moment peak $p$ fires, in the lab frame:

- $\mathbf{g}_p = R_z(\omega^*_p)\,\mathbf{g}_{s,p}$ and the diffracted
  direction $\hat{\mathbf{k}}'_p = (\mathbf{k} + \mathbf{g}_p)/|\mathbf{k}|$;
- $\nabla_L\omega^*_p = -\hat{\mathbf{c}}_p/(\hat{\mathbf{z}}\cdot\hat{\mathbf{c}}_p)$
  (Section 3.6.2);
- $\mathbf{x}_{v,p}$, the lab position of the voxel's centroid, which is the
  ray's starting point;
- the detector: unit normal $\hat{\mathbf{n}}$, lab directions of its column
  and row axes $\hat{\mathbf{e}}_{\text{col}}, \hat{\mathbf{e}}_{\text{row}}$,
  and pixel size $a$. Write $\Pi$ for the $2\times3$ matrix with rows
  $\hat{\mathbf{e}}_{\text{col}}^{\top}/a$ and $\hat{\mathbf{e}}_{\text{row}}^{\top}/a$
  (converts a lab displacement in the detector plane to pixels);
- $\Lambda_p$, the length of the ray from the voxel to the detector;
- $[\mathbf{v}]_\times$, the skew-symmetric matrix with
  $[\mathbf{v}]_\times\mathbf{w} = \mathbf{v}\times\mathbf{w}$ (Section 1.2).

**Step 1: where a ray hits the detector.** A ray from $\mathbf{x}_o$ along
$\hat{\mathbf{k}}'$ meets the detector plane at
$\mathbf{s} = \mathbf{x}_o + \Lambda\,\hat{\mathbf{k}}'$, with $\Lambda$ fixed by
$\hat{\mathbf{n}}\cdot\mathbf{s}$ equalling the plane's offset. Differentiating,
and eliminating $d\Lambda$ with the plane condition
$\hat{\mathbf{n}}\cdot d\mathbf{s} = 0$,
$$
d\mathbf{s} = Q_p\,\big(d\mathbf{x}_o + \Lambda_p\,d\hat{\mathbf{k}}'\big),
\qquad
Q_p = I - \frac{\hat{\mathbf{k}}'_p\,\hat{\mathbf{n}}^{\top}}{\hat{\mathbf{n}}\cdot\hat{\mathbf{k}}'_p} .
$$
$Q_p$ projects a displacement onto the detector plane *along the ray*. It
accounts for rays hitting the detector at an angle, which Section 3.6.2's
head-on approximation ignored.

**Step 2: the direction changes.** The incident beam is fixed, so
$d\hat{\mathbf{k}}' = d\mathbf{g}/|\mathbf{k}|$. From Section 3.6.2,
$d\mathbf{g} = \Omega\times\mathbf{g}$ with
$\Omega = \delta_L + (\nabla_L\omega^*\cdot\delta_L)\,\hat{\mathbf{z}}$, so
$$
d\hat{\mathbf{k}}' = -\frac{1}{|\mathbf{k}|}\,[\mathbf{g}_p]_\times\big(I + \hat{\mathbf{z}}\,\nabla_L\omega^{*\top}_p\big)\,\delta_L .
$$

**Step 3: the starting point moves.** The voxel rotates with the stage through
the extra angle $d\omega^* = \nabla_L\omega^*\cdot\delta_L$:
$$
d\mathbf{x}_o = (\hat{\mathbf{z}}\times\mathbf{x}_{v,p})\,\nabla_L\omega^{*\top}_p\,\delta_L .
$$

**Step 4: assemble.** Convert to pixels with $\Pi$ and to the sample frame with
$\delta_L = R_z(\omega^*_p)\,\delta$:
$$
\Gamma_p = \Pi\,Q_p\Big[\underbrace{-\frac{\Lambda_p}{|\mathbf{k}|}\,[\mathbf{g}_p]_\times\big(I + \hat{\mathbf{z}}\,\nabla_L\omega^{*\top}_p\big)}_{\text{ring (direction) term}}
\;+\;\underbrace{(\hat{\mathbf{z}}\times\mathbf{x}_{v,p})\,\nabla_L\omega^{*\top}_p}_{\text{parallax term}}\Big]\,R_z(\omega^*_p) .
$$
This is first order in $\delta$ and makes no further approximation: it holds for
any flat detector orientation and any voxel position. It applies to each
vertex of the spot's footprint, using that vertex's own starting point and ray
length; because projection along parallel rays is affine, the Jacobian of the
spot's centroid is the average of the three vertex Jacobians, which is what
using the voxel's centroid as $\mathbf{x}_{v,p}$ gives.

**Special case: detector perpendicular to the beam.** With
$\hat{\mathbf{n}} = \hat{\mathbf{x}}$, $\hat{\mathbf{c}}_p$ lies in the detector
plane, $Q_p\hat{\mathbf{c}}_p = \hat{\mathbf{c}}_p$, and
$\Lambda_p = d_p/\cos 2\theta_p$. The ring term becomes
$(d_p\tan 2\theta_p/\cos\theta_p)\;\Pi\,\hat{\mathbf{c}}_p\,[\hat{\mathbf{t}}_p + (\hat{\mathbf{z}}\cdot\hat{\mathbf{t}}_p)\nabla_L\omega^*_p]^{\top}$,
Section 3.6.2's result ($\Pi$ supplies the $1/a$), with the spot moving along
$\hat{\mathbf{c}}_p$, the ring tangent. The parallax term's in-plane part is $r_\perp\cos\phi_{v,p}/a$ along
the horizontal axis, as in Section 3.6.2, plus an obliquity correction
proportional to $r_\perp\sin\phi_{v,p}\,\tan 2\theta_p$ that the head-on
approximation dropped.

**Structural properties.**

1. **Blind axis.** With $\hat{\mathbf{g}}_{s,p}$ the unit vector along
   $\mathbf{g}_{s,p}$, $\Gamma_p\,\hat{\mathbf{g}}_{s,p} = 0$ exactly, because
   $[\mathbf{g}]_\times\mathbf{g} = 0$ and $\nabla\omega^*\cdot\hat{\mathbf{g}} = 0$
   (since $\hat{\mathbf{c}} \perp \mathbf{g}$). Rotation about the peak's own
   reciprocal vector changes nothing, so $\Gamma_p$ has rank at most two.
2. **The $\hat{\mathbf{z}}$ column is pure parallax.** Because
   $(I + \hat{\mathbf{z}}\nabla_L\omega^{*\top})\hat{\mathbf{z}} = \hat{\mathbf{z}} - \hat{\mathbf{z}} = 0$
   and $\nabla_L\omega^*\cdot\hat{\mathbf{z}} = -1$,
   $$
   \Gamma_p\,\hat{\mathbf{z}} = -\Pi\,Q_p\,(\hat{\mathbf{z}}\times\mathbf{x}_{v,p}) :
   $$
   minus the voxel's velocity per radian of stage rotation, as seen on the
   detector along the ray. It vanishes for a voxel on the axis. This is the exact
   form of Section 3.1.5, point 2.
3. **Rank.** The ring term alone has rank one (it always moves the spot along
   $Q_p\hat{\mathbf{c}}_p$). Near the axis the parallax term is small, so the
   spot moves essentially along a line; off the axis $\Gamma_p$ has rank two.
4. **Full per-peak observation Jacobian.** Stacking the frame gradient
   $\nabla\omega^{*\top}_p = \nabla_L\omega^{*\top}_p R_z(\omega^*_p)$ on top of
   $\Gamma_p$ gives the $3\times3$ map from $\delta$ to (frame angle, column,
   row). It also annihilates $\hat{\mathbf{g}}_{s,p}$, so each peak constrains
   at most the two components of $\delta$ perpendicular to its reciprocal
   vector.
5. **Geometric reading.** By the singular value decomposition, $\Gamma_p$ maps a
   ball of orientation offsets of radius $\beta$ to an ellipse of spot positions
   with semi-axes $\beta$ times its singular values. That ellipse is the region
   a window around the spot must cover.

**Worked example.** Reflection 362, branch 2, of the running example's voxel 0
($\omega^* = 7.7^\circ$, spot near column 819, row 1575, voxel lab azimuth
$-172^\circ$, i.e. upstream of the axis). Entries in px/rad; columns are
rotation about sample $x$, $y$, $z$:

| $r_\perp$ | | ring term | parallax term | $\Gamma_p$ (closed form) | finite differences |
|---|---|---|---|---|---|
| $12\,\mu$m | col | $(-339.5,\ 108.9,\ 0)$ | $(1.7,\ 12.9,\ 7.5)$ | $(-337.8,\ 121.9,\ 7.5)$ | $(-337.7,\ 121.8,\ 7.5)$ |
| | row | $(197.3,\ -63.3,\ 0)$ | $(0.0,\ -0.3,\ -0.2)$ | $(197.2,\ -63.6,\ -0.2)$ | $(197.3,\ -63.6,\ -0.1)$ |
| $500\,\mu$m | col | $(-388.1,\ 124.5,\ 0)$ | $(76.7,\ 569.0,\ 329.8)$ | $(-311.4,\ 693.5,\ 329.8)$ | $(-311.4,\ 693.5,\ 329.8)$ |
| | row | $(225.5,\ -72.3,\ 0)$ | $(-2.1,\ -15.9,\ -9.2)$ | $(223.3,\ -88.2,\ -9.2)$ | $(223.3,\ -88.3,\ -9.2)$ |

Reading it:

- Near the axis the parallax term is tiny and all three columns point the same
  way on the detector: singular values 415 and 7 px/rad, so the spot moves along
  a line (the ring tangent), about $7\,$px per degree.
- The blind axis, $(-0.155, -0.483, 0.862)$, equals the direction of
  $\mathbf{g}_{s,p}$, $(-0.154, -0.480, 0.864)$.
- At $500\,\mu$m the voxel sits $0.5$ mm upstream, so the ray is 15% longer and
  the ring term 15% larger. The parallax term, $(\hat{\mathbf{z}}\times\mathbf{x}_v)$
  times this peak's frame gradient $\nabla\omega^* = (-0.23, -1.73, -1)$, adds
  mainly to the column entries and makes the $z$ column large. Singular values
  become 845 and 175 px/rad: the spot now moves over a 2D region.
- The small row entries of the parallax term are the obliquity correction the
  head-on approximation missed.

**Check.** For all 790 ROI peaks of voxel 0, at $r_\perp$ = 12, 100, 250 and
$500\,\mu$m, the closed form was compared with central finite differences
(step $10^{-3}$ rad) of the simulator's projected spot centroid. Relative
Frobenius error: median $1.5\times10^{-4}$, 90th percentile $3\times10^{-4}$,
maximum 1.1% (near-axis peaks, where the finite-difference step itself sees
curvature). $|\Gamma_p\,\hat{\mathbf{g}}_{s,p}|$ is zero to machine precision.

### 3.6.4 Combining peaks

With the error model of Section 3.6.1, the information matrix is
$$
J = \frac{12}{\Delta\omega_f^2}\sum_{p=1}^{P}\nabla\omega^*_p\,\nabla\omega^{*\top}_p
\;+\; 12\sum_{p=1}^{P}\Gamma_p^{\top}\Gamma_p ,
$$
with $\Gamma_p$ the spot-motion Jacobian of Section 3.6.3, and the covariance of
a least-squares estimate of $\delta$ is approximately $J^{-1}$. Using the
simplified ring and parallax terms of Section 3.6.2, and treating the two
displacements as independent pixel measurements, gives closed forms.

- **About the rotation axis.** Every peak contributes exactly
  $12/\Delta\omega_f^2$ through its frame and
  $12\,(r_\perp\cos\phi_{v,p}/a)^2$ through parallax. With peaks firing at
  well-spread $\omega^*$, $\cos^2\phi_{v,p}$ averages to $1/2$, so
  $J_{zz} \approx 12P\,[\,1/\Delta\omega_f^2 + r_\perp^2/(2a^2)\,]$. When pixel
  information dominates the perpendicular block, the coupling between blocks is
  negligible and
  $$
  \sigma_z \approx \frac{1}{\sqrt{12\,P\left(\dfrac{1}{\Delta\omega_f^2} + \dfrac{r_\perp^2}{2a^2}\right)}} .
  $$
  Two limits:
  - *near the axis* ($r_\perp \ll \sqrt{2}\,a/\Delta\omega_f$): frames dominate,
    $\sigma_z \approx \Delta\omega_f/\sqrt{12P} = \sqrt{\Delta\omega_f/(12\nu N)}$;
  - *far from the axis* ($r_\perp \gg \sqrt{2}\,a/\Delta\omega_f$): parallax
    dominates, $\sigma_z \approx a/(r_\perp\sqrt{6P})$, independent of the frame
    width.

  The crossover is at $r_\perp^* = \sqrt{2}\,a/\Delta\omega_f$, about $120\,\mu$m
  for $1.48\,\mu$m pixels and $1^\circ$ frames.
- **Perpendicular to the axis.** The $\nabla\xi_p$ lie in the plane
  perpendicular to $\hat{\mathbf{z}}$; if their directions are spread over it,
  each perpendicular axis receives half of $12\sum_p|\nabla\xi_p|^2$.
  Approximating the bracket in $\nabla\xi_p$ by a unit vector and $d_p$ by $d$
  (the bracket's squared length is exactly $1 + \sin^2\theta_p\cot^2\eta_p$, so
  this underestimates the information from near-axis peaks, by up to a factor of
  about 10 in squared length for the running example's nearest-axis peaks),
  $|\nabla\xi_p| \approx \kappa_p\,d/a$ with $\kappa_p = \tan 2\theta_p/\cos\theta_p \approx 2\theta_p$,
  $$
  \sigma_\perp \approx \frac{a/d}{\kappa\,\sqrt{6\,P}} = \frac{a/d}{\kappa\,\sqrt{6\,\nu\,N\,\Delta\omega_f}},
  $$
  where $\kappa$ is the root-mean-square of $\kappa_p$ over peaks. (Parallax
  adds a little perpendicular information too, which this neglects.)

The RMS misorientation is then $\sqrt{2\sigma_\perp^2 + \sigma_z^2}$. The two
formulas share a form: each is $a$ divided by a lever arm and by $\sqrt{6P}$.
Perpendicular to the axis the lever arm is $\kappa d \approx$ the ring radius;
about the axis it is $r_\perp$ (far from the axis) or the frame-equivalent
$\sqrt{2}\,a/\Delta\omega_f$ (near it). The anisotropy is therefore
$$
\frac{\sigma_z}{\sigma_\perp} \approx \frac{\kappa\,d}{\sqrt{2a^2/\Delta\omega_f^2 + r_\perp^2}}
\;\longrightarrow\;
\begin{cases}
\kappa\,\Delta\omega_f/(\sqrt{2}\,a/d) & \text{near the axis,}\\[2pt]
\kappa\,d/r_\perp & \text{far from it.}
\end{cases}
$$

### 3.6.5 Numerical check

Running example, voxel 0: 790 ROI peaks, all assigned to the first detector
($d = 3.36$ mm, $a = 1.48\,\mu$m, so $a/d = 4.4\times10^{-4}$ rad), 788 with
usable finite differences. Spot positions were taken from the simulator's
projected centroid and differentiated numerically.

- **Ring motion:** voxel 0 is only about $12\,\mu$m from the axis, so parallax is
  negligible and the spot-position Jacobian is rank one, as predicted: its
  second singular value is 0.6% of the first (median; at most 3%).
- **Magnitude of $\nabla\xi_p$:** matches the formula with median relative error
  0.75% (90th percentile 2.8%). A typical spot moves about 1,010 px/rad
  (17.6 px/deg) for rotations perpendicular to the axis, and 5.9 px/rad for
  rotation about the axis (parallax only, consistent with the 8 px/rad scale).
- **Resolution** ($\kappa = 0.43$):

| | $\sigma_x$, $\sigma_y$ | $\sigma_z$ |
|---|---|---|
| $J^{-1}$, frames only | $0.0098^\circ$, $0.0062^\circ$ | $0.010^\circ$ |
| $J^{-1}$, pixels only | $0.00079^\circ$ | $0.10^\circ$ (parallax only) |
| $J^{-1}$, frames + pixels | $0.00079^\circ$ | $0.0102^\circ$ |
| Closed forms (3.6.4) | $0.00086^\circ$ | $0.0103^\circ$ |

The predicted anisotropy is 12; $J^{-1}$ gives 13.

**Off-axis voxels.** The same voxel was moved to larger $r_\perp$ (keeping its
azimuth and orientation) and the analysis repeated on a random subset of 250 of
its peaks, so the $\sigma$ values here are larger than in the table above, which
used all 788. The prediction for the $\hat{\mathbf{z}}$ sensitivity is
$r_\perp|\cos\phi_v|/a$, whose median over uniformly spread $\phi_v$ is
$0.71\,r_\perp/a$. Here $r_\perp$ labels the distance of the voxel's reference
vertex (its position in the `.mic` file) from the axis; the triangle's centroid,
which is what enters the formulas, differs by less than $1\,\mu$m (about
$11.3\,\mu$m for the $12\,\mu$m row).

| $r_\perp$ | Jacobian 2nd/1st singular value (median / max) | $\hat{\mathbf{z}}$ sensitivity, measured / $(r_\perp/a)$ | $\sigma_z$ from $J^{-1}$ | $\sigma_z$ closed form |
|---|---|---|---|---|
| $12\,\mu$m | 0.006 / 0.03 | 0.73 | $0.018^\circ$ | $0.018^\circ$ |
| $50\,\mu$m | 0.026 / 0.14 | 0.72 | $0.017^\circ$ | $0.017^\circ$ |
| $100\,\mu$m | 0.052 / 0.25 | 0.73 | $0.014^\circ$ | $0.014^\circ$ |
| $250\,\mu$m | 0.13 / 0.46 | 0.73 | $0.0075^\circ$ | $0.0079^\circ$ |
| $500\,\mu$m | 0.23 / 0.82 | 0.76 | $0.0039^\circ$ | $0.0043^\circ$ |

Over the same range the perpendicular errors stayed between $0.0008^\circ$ and
$0.0015^\circ$, so the anisotropy fell from about 12 to about 4. The closed
forms predict about 2.8 at $500\,\mu$m; the difference is that parallax also
adds perpendicular information, which the $\sigma_\perp$ closed form neglects,
so the measured $\sigma_\perp$ ($\approx 0.001^\circ$) is below its closed form
($0.0015^\circ$ for $P = 250$). The Jacobian
becomes clearly rank two as the voxel moves off the axis, confirming that spots
no longer move along a single line.

### 3.6.6 Scaling and caveats

- **Fixed total sweep** ($N\Delta\omega_f$ constant): for voxels near the axis,
  $\sigma_z \propto \Delta\omega_f$, so halving the frame width halves the error
  about the axis; far from the axis, $\sigma_z$ no longer depends on the frame
  width. $\sigma_\perp$ is unchanged either way.
- **Fixed frame width:** both errors fall as $1/\sqrt{N}$.
- **Pixel size** enters $\sigma_\perp$ through $a/d$ and, far from the axis,
  $\sigma_z$ through $a/r_\perp$.
- **Voxel position matters.** Resolution about the rotation axis is not uniform
  across the sample: it improves with distance from the axis once
  $r_\perp \gtrsim \sqrt{2}\,a/\Delta\omega_f$.
- **Validity.** The parallax result is first order in $r_\perp/d$ and assumes
  near-normal incidence on the detector; it has been checked only up to
  $r_\perp = 0.5$ mm with $d = 3.36$ mm.

These are idealised floors, well below the roughly $0.1^\circ$ typical of real
high-energy X-ray diffraction microscopy (HEDM) reconstructions. The check uses a noise-free, perfectly calibrated
simulator and all 790 of voxel 0's peaks (reflections up to the config's
$Q_{\max} = 16$ Å$^{-1}$, whereas reconstruction typically uses
$Q_{\max} = 8$ Å$^{-1}$). Real resolution
is degraded by detector point-spread, spot footprint, calibration errors in $d$,
beam centre and detector tilt, strain (which lets spots move radially),
intensity noise and thresholding, and overlapping peaks. It would be improved by
the second detector, which this estimate ignores.

The independent-error model gives the $1/\sqrt{P}$ scaling. The noise-free
set-membership picture (Section 3.1.5, point 3) scales as $1/P$ instead, but
only if frame edges are known to better than $\Delta\omega_f/P$, which is not
realistic.

# 4. The Bayesian inverse

## 4.1 Posterior

Write $\mathcal{P}(D\mid\delta)$ for the likelihood, the probability of data $D$
given offset $\delta$. The Bayesian formulation replaces "solve $M(\delta) = D$"
with the **posterior**
$$
\pi(\delta \mid D) \propto \mathcal{P}(D \mid \delta)\,\pi(\delta),
$$
which exists and is well-behaved even when the inverse is not (Stuart 2010;
Kaipio & Somersalo 2005; Tarantola 2005). It describes every orientation
consistent with the data and how plausible each is.

## 4.2 Noise-free limit: the prior restricted to the cell

Without noise the likelihood is an indicator, $\mathcal{P}(D\mid\delta) = \mathbf{1}[M(\delta) = D]$, so
$$
\pi(\delta \mid D) = \frac{\pi(\delta)\,\mathbf{1}[\delta \in C(D)]}{\pi\big(C(D)\big)} ,
$$
where $\pi(C(D))$ is the prior probability of the consistent set. The posterior
is the prior, cut down to the consistent set and renormalised.

## 4.3 Estimators

$E[\cdot]$, $\mathrm{Cov}[\cdot]$ and $\operatorname{tr}$ denote expectation,
covariance and matrix trace.

- **Posterior mean** $E[\delta \mid D]$: the prior-weighted centroid of $C(D)$.
- **Posterior covariance** $\mathrm{Cov}[\delta \mid D]$: the prior-weighted
  second moment of $C(D)$ about its centroid. It is anisotropic, because the
  slabs have different orientations and widths.
- **Maximum a posteriori (MAP) / minimum cost:** with a flat prior, every point
  of $C(D)$ is equally optimal, so the maximiser is not unique. An argmin-cost
  search (IceNine's Monte Carlo (MC) search or Riemannian Adam) returns *some*
  point of the set, plus search error.

## 4.4 The error floor (Bayes risk)

The smallest achievable mean squared error, over all estimators, is the
**minimum mean squared error**
$$
\text{MMSE} = E_D\big[\operatorname{tr}\mathrm{Cov}[\delta \mid D]\big],
$$
where $E_D$ averages over data generated by drawing $\delta \sim \pi$ and
setting $D = M(\delta)$. It is attained by the posterior mean (Lehmann &
Casella 1998). For comparison, an estimate drawn at random from the posterior
has expected squared error $2\operatorname{tr}\mathrm{Cov}[\delta\mid D]$, twice
the floor. A search that lands on an arbitrary point of the consistent set is in
this spirit, though its answer is not literally a posterior draw.

# 5. What supervised training learns

## 5.1 Squared loss learns the posterior mean

Generate training pairs by sampling $\delta \sim \pi$ and simulating
$D = M(\delta)$ (plus noise, if modelled). For any estimator $f$, with
expectations over these pairs,
$$
E\,\|f(D) - \delta\|^2
= E\,\|f(D) - E[\delta\mid D]\|^2 + E\,\|E[\delta\mid D] - \delta\|^2 ,
$$
because the cross term vanishes by the tower property of conditional
expectation. The second term does not depend on $f$, so the minimiser is
$f^*(D) = E[\delta \mid D]$, and the minimum value is the MMSE of Section 4.4.
This is the standard result that the squared-loss optimal predictor is the
conditional mean (Bishop 2006, §1.5.5; Lehmann & Casella 1998). Adler & Öktem
(2018, §4.2 and App. C.3) use exactly this to compute posterior means for CT
reconstruction with a directly trained network.

Two consequences:

- **The prior is part of the answer.** $f^*$ depends on the training
  distribution $\pi$. It must match how the network will be used, for example
  the residual error left by the coarse Sukharev search.
- **The best possible test error is the MMSE, not zero.** Network error should be
  reported relative to it.

## 5.2 Learning the uncertainty

Two established routes:

1. **Residual regression.** Train a second network on the squared residuals of
   the first; it converges to the conditional variance (Adler & Öktem 2018).
2. **Mean–variance estimation.** Output a mean and a covariance and train with
   the Gaussian negative log-likelihood (NLL) (Nix & Weigend 1994). Kendall &
   Gal (2017) call this heteroscedastic *aleatoric* uncertainty: uncertainty that
   is in the data and does not go away with more training data. Binning
   uncertainty is exactly of this kind. Plain Gaussian NLL can train poorly;
   Seitzer et al. (2022) diagnose the problem and propose a reweighted loss they
   call $\beta$-NLL (their $\beta$ is a weighting exponent, unrelated to the
   rotation angle $\beta$ above).

Because the consistent set is anisotropic (Section 4.3), the network should
output a full $3\times3$ covariance (for example via its Cholesky factor), not
three independent variances.

## 5.3 Amortised, simulation-based inference

Solving each voxel by search pays the full cost of inference every time. A
network trained once on simulations and then applied in a single forward pass
*amortises* that cost (Gershman & Goodman 2014). Learning posteriors from a
simulator this way is simulation-based inference (Cranmer, Brehmer & Louppe
2020).

A Gaussian is only an approximation here: with a locally flat prior, the
noise-free posterior is uniform on the consistent set (Section 4.2), which in
the linear regime is a polytope. If the full posterior shape matters, neural
posterior estimation learns a flexible conditional density (Papamakarios &
Murray 2016). On the rotation group, parametric options are the matrix Fisher
(Mohlin, Sullivan & Bianchi 2020) and Bingham (Gilitschenski et al. 2020)
distributions, and Implicit-PDF (Murphy et al. 2021) represents arbitrary,
multimodal densities on SO(3).

## 5.4 Rotation-specific points

- **Output a tangent vector, not an absolute quaternion.** For small $\delta$,
  $\exp$ is nearly an isometry near the identity, so Euclidean means and
  covariances of $\delta$ are meaningful. Quaternion and Euler-angle outputs are
  discontinuous representations that networks learn poorly (Zhou et al. 2019).
- **Larger spreads need rotation means.** If perturbations grow, the posterior
  mean should be a proper mean on SO(3), such as the Riemannian (Karcher) mean
  (Moakher 2002).

# 6. Simulator artefacts the network could exploit

The network learns $E[\delta\mid D]$ under the *simulator's* likelihood. Anything
the simulator encodes more precisely than a real experiment is a shortcut that
will not transfer.

- **Exact intensities.** $L_p(\delta)$ varies smoothly with orientation (through
  $\eta_p$ and the Lorentz factor). If the network is given unthresholded,
  noise-free data $A(\delta)$ rather than $M(\delta)$, the exact intensity values
  leak orientation information that real, noisy data do not reliably provide,
  and that IceNine's own binarised cost ignores. Mitigation: train on
  thresholded data, or add realistic intensity noise and scaling; at minimum,
  evaluate on thresholded inputs.
- **All-or-nothing pixel fill.** The simulator lacks partial coverage and
  point-spread, so it has *less* sub-pixel information than real data. A network
  trained on it will not learn to use that information.
- **Single-frame peaks** ($\alpha = 0$). Real near-axis peaks spread across
  frames (Section 3.5); the simulator currently discards that information.

# 7. Relation to IceNine's existing methods

| Method | What it computes | Relation to the posterior |
|---|---|---|
| Hard cost + MC search | Minimises binary pixel overlap by random search | Returns some point of $C(D)$, plus search error |
| Differentiable cost + Riemannian Adam | Minimises a blurred (mollified) cost by gradient descent | Point estimate of a smoothed problem; blur trades resolution for gradients |
| Orientation network | One forward pass | Approximates $E[\delta\mid D]$ and, with NLL, $\mathrm{Cov}[\delta\mid D]$ |

All three are bounded by the same floor, the MMSE of Section 4.4. For the
physics of the forward model itself, see Suter et al. (2006) and Li & Suter
(2013). BraggNN (Liu et al. 2022) is a related use of deep learning in
HEDM, for localising individual
Bragg peaks rather than orientations.

# 8. Implications for the next stages

The staged plan these implications feed into, with its open decisions, is
recorded in `icenine_py/MIGRATION_HISTORY.md` (toy orientation NN section).

1. **Output and loss.** Predict $\delta$ and a full covariance; train with
   Gaussian NLL ($\beta$-NLL if training is unstable).
2. **Exact Bayes baseline (new).** Sample $\delta \sim \pi$ and keep the samples
   whose *exact* frame assignments (computed with `get_scattering_omegas_torch`)
   match the observed ones; their mean and covariance are the posterior mean and
   covariance, and averaging over test cases gives the MMSE floor. Do not use the
   linear polytope for this: Section 3.1.6 shows it is wrong for exactly the
   near-axis peaks that matter most. The polytope is still useful for intuition
   and as a fast approximation when only peaks with $|\sin\eta| \gtrsim 0.3$ are
   used. Start with frame constraints, then add pixel constraints. Compare the
   network against this, not only against predict-nominal.
3. **Prior.** Choose $\pi$ to match deployment: the residual error typically left
   by the coarse Sukharev search (an open decision).
4. **Avoid intensity leakage.** Train on thresholded data $M(\delta)$, or add
   noise to $A(\delta)$; at least evaluate on thresholded inputs.
5. **Frame spread.** Model $\alpha > 0$ in the renderer so the network can use
   sub-frame information from near-axis peaks.
6. **Report errors per axis.** Because resolution about the rotation axis
   (frame-limited) and perpendicular to it (pixel-limited) differ by an order of
   magnitude (Section 3.6), report network and baseline errors separately for
   the $\hat{\mathbf{z}}$ and perpendicular components, and compare each with
   its Section 3.6 estimate.

# References

- Adler, J. & Öktem, O. (2018). Deep Bayesian inversion. arXiv:1811.05910.
  <https://arxiv.org/abs/1811.05910>
- Arridge, S., Maass, P., Öktem, O. & Schönlieb, C.-B. (2019). Solving inverse
  problems using data-driven models. *Acta Numerica* 28, 1–174.
  doi:10.1017/S0962492919000059
- Bishop, C. M. (2006). *Pattern Recognition and Machine Learning*. Springer.
  §1.5.5, Loss functions for regression.
- Cranmer, K., Brehmer, J. & Louppe, G. (2020). The frontier of simulation-based
  inference. *PNAS* 117(48), 30055–30062. doi:10.1073/pnas.1912789117
- Engl, H. W., Hanke, M. & Neubauer, A. (1996). *Regularization of Inverse
  Problems*. Mathematics and Its Applications 375. Kluwer.
- Gershman, S. J. & Goodman, N. D. (2014). Amortized inference in probabilistic
  reasoning. *Proc. 36th Annual Conference of the Cognitive Science Society*,
  517–522.
- Gilitschenski, I., Sahoo, R., Schwarting, W., Amini, A., Karaman, S. & Rus, D.
  (2020). Deep orientation uncertainty learning based on a Bingham loss. *ICLR*.
- Gray, R. M. & Neuhoff, D. L. (1998). Quantization. *IEEE Trans. Information
  Theory* 44(6), 2325–2383.
- Hadamard, J. (1902). Sur les problèmes aux dérivées partielles et leur
  signification physique. *Princeton University Bulletin* 13, 49–52.
- Kaipio, J. & Somersalo, E. (2005). *Statistical and Computational Inverse
  Problems*. Applied Mathematical Sciences 160. Springer.
- Kendall, A. & Gal, Y. (2017). What uncertainties do we need in Bayesian deep
  learning for computer vision? *NeurIPS*.
- Lehmann, E. L. & Casella, G. (1998). *Theory of Point Estimation*, 2nd ed.
  Springer.
- Li, S. F. & Suter, R. M. (2013). Adaptive reconstruction method for
  three-dimensional orientation imaging. *J. Appl. Cryst.* 46, 512–524.
  doi:10.1107/S0021889813005268
- Liu, Z., Sharma, H., Park, J.-S., Kenesei, P., Miceli, A., Almer, J.,
  Kettimuthu, R. & Foster, I. (2022). BraggNN: fast X-ray Bragg peak analysis
  using deep learning. *IUCrJ* 9, 104. doi:10.1107/S2052252521011258
- Milanese, M. & Vicino, A. (1991). Optimal estimation theory for dynamic systems
  with set membership uncertainty: an overview. *Automatica* 27(6), 997–1009.
- Moakher, M. (2002). Means and averaging in the group of rotations. *SIAM J.
  Matrix Anal. Appl.* 24(1), 1–16.
- Mohlin, D., Sullivan, J. & Bianchi, G. (2020). Probabilistic orientation
  estimation with matrix Fisher distributions. *NeurIPS*. arXiv:2006.09740.
- Murphy, K., Esteves, C., Jampani, V., Ramalingam, S. & Makadia, A. (2021).
  Implicit-PDF: non-parametric representation of probability distributions on
  the rotation manifold. *ICML*, PMLR 139, 7882–7893. arXiv:2106.05965.
- Nix, D. A. & Weigend, A. S. (1994). Estimating the mean and variance of the
  target probability distribution. *Proc. IEEE Int. Conf. Neural Networks*,
  vol. 1, 55–60. doi:10.1109/ICNN.1994.374138
- Papamakarios, G. & Murray, I. (2016). Fast $\varepsilon$-free inference of simulation
  models with Bayesian conditional density estimation. *NeurIPS*.
  arXiv:1605.06376.
- Seitzer, M., Tavakoli, A., Antic, D. & Martius, G. (2022). On the pitfalls of
  heteroscedastic uncertainty estimation with probabilistic neural networks.
  *ICLR*. arXiv:2203.09168.
- Stanley, R. P. (2007). An introduction to hyperplane arrangements. In
  *Geometric Combinatorics*, IAS/Park City Mathematics Series 13, 389–496. AMS.
- Stuart, A. M. (2010). Inverse problems: a Bayesian perspective. *Acta
  Numerica* 19, 451–559. doi:10.1017/S0962492910000061
- Suter, R. M., Hennessy, D., Xiao, C. & Lienert, U. (2006). Forward modeling
  method for microstructure reconstruction using x-ray diffraction microscopy:
  single-crystal verification. *Rev. Sci. Instrum.* 77, 123905.
  doi:10.1063/1.2400017
- Tarantola, A. (2005). *Inverse Problem Theory and Methods for Model Parameter
  Estimation*. SIAM.
- Thompson, R. E., Larson, D. R. & Webb, W. W. (2002). Precise nanometer
  localization analysis for individual fluorescent probes. *Biophys. J.* 82(5),
  2775–2783. doi:10.1016/S0006-3495(02)75618-X
- Widrow, B. & Kollár, I. (2008). *Quantization Noise: Roundoff Error in Digital
  Computation, Signal Processing, Control, and Communications*. Cambridge
  University Press.
- Zhou, Y., Barnes, C., Lu, J., Yang, J. & Li, H. (2019). On the continuity of
  rotation representations in neural networks. *CVPR*, 5738–5746.
