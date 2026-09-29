---
title: "Angular Width of a Diffraction Peak in a Rotation Scan"
subtitle: "Why peaks near the rotation axis spread across multiple omega frames"
date: "2026-09-28"
geometry: margin=1in
fontsize: 11pt
---

# Status

This is a first-order derivation written for the IceNine toy orientation-NN work.
Its relations have been checked numerically against the simulator's Bragg
solver (`get_scattering_omegas_torch`) on Example2.ThreeVoxels and agree to
floating-point precision (Section 7.1). It has **not** yet been checked against
real data or cross-checked against a crystallography reference (Section 7.2).

# 1. Problem

In a rotation scan the sample turns about a fixed axis while the detector
integrates each frame over an interval $[\omega_j - \Delta\omega_f/2,\;
\omega_j + \Delta\omega_f/2]$. A reflection $\mathbf{g}_{hkl}$ produces
intensity only while the Bragg condition is satisfied.

For a perfect crystal illuminated by a perfectly parallel, monochromatic beam,
the Bragg condition holds at a single instant $\omega^*$, so every peak falls
in exactly one frame. Real data show peaks spread over several adjacent
frames, much more often for reflections whose $\mathbf{g}$ lies close to the
rotation axis. This note derives how long (in $\omega$) a reflection stays in
the diffraction condition, given a small intrinsic tolerance on the Bragg
condition.

# 2. Geometry and assumptions

- The rotation axis is $\hat{\mathbf{z}}$. The sample orientation at rotation
  angle $\omega$ is $R_z(\omega)$ applied to the voxel's lattice.
- The incident wavevector $\mathbf{k}$ is perpendicular to the rotation axis
  ($\mathbf{k}\cdot\hat{\mathbf{z}} = 0$). For Example2 the beam direction is
  $(1, 0, 0)$, so this holds.
- $|\mathbf{k}| = 2\pi/\lambda$ and $|\mathbf{g}| = 2|\mathbf{k}|\sin\theta$,
  where $\theta$ is the Bragg angle.
- $\chi$ is the angle between $\mathbf{g}$ and the rotation axis, and
  $\varphi_0$ is the azimuth of $\mathbf{g}$ about the axis at $\omega = 0$:
  $$
  \mathbf{g}(\omega) = |\mathbf{g}|
  \begin{pmatrix}
  \sin\chi\cos(\omega+\varphi_0) \\
  \sin\chi\sin(\omega+\varphi_0) \\
  \cos\chi
  \end{pmatrix}.
  $$
- Each voxel is treated as a perfect crystal. Any finite width comes from the
  beam (energy bandwidth, divergence), not from mosaic spread inside the voxel.

# 3. The Bragg condition as a function of omega

Elastic scattering requires $|\mathbf{k} + \mathbf{g}| = |\mathbf{k}|$, which is
equivalent to
$$
f(\omega) \equiv \mathbf{k}\cdot\mathbf{g}(\omega) + \tfrac{1}{2}|\mathbf{g}|^2 = 0 .
$$
Only the component of $\mathbf{g}$ perpendicular to the axis rotates, and
$\mathbf{k}$ lies in that plane, so
$$
\mathbf{k}\cdot\mathbf{g}(\omega) = |\mathbf{k}||\mathbf{g}|\sin\chi\,\cos(\omega + \varphi_0).
$$
Using $|\mathbf{g}| = 2|\mathbf{k}|\sin\theta$, the condition $f(\omega^*) = 0$
becomes
$$
\cos(\omega^* + \varphi_0) = -\frac{\sin\theta}{\sin\chi}.
$$
A solution exists only when $\sin\chi \ge \sin\theta$. When $\chi < \theta$ the
reflection never diffracts during the rotation. When $\chi > \theta$ there are
two solutions (the two branches `omega1`/`omega2` returned by
`get_scattering_omegas_torch`), and they merge as $\chi \to \theta$.

# 4. Rate at which the condition is crossed

Differentiating,
$$
f'(\omega) = -|\mathbf{k}||\mathbf{g}|\sin\chi\,\sin(\omega + \varphi_0).
$$
At a solution, $\sin(\omega^* + \varphi_0) = \pm\sqrt{1 - \sin^2\theta / \sin^2\chi}$,
so
$$
\left|f'(\omega^*)\right|
= |\mathbf{k}||\mathbf{g}|\sin\chi\sqrt{1 - \frac{\sin^2\theta}{\sin^2\chi}}
= |\mathbf{k}||\mathbf{g}|\sqrt{\sin^2\chi - \sin^2\theta}.
$$

# 5. Peak width in omega

## 5.1 Tolerance on the Bragg condition

Let $\psi$ be the glancing angle between $\mathbf{k}$ and the reflecting lattice
planes, so that $\mathbf{k}\cdot\hat{\mathbf{g}} = -|\mathbf{k}|\sin\psi$. Then
$$
f = |\mathbf{g}|\left(\mathbf{k}\cdot\hat{\mathbf{g}} + \tfrac{1}{2}|\mathbf{g}|\right)
  = |\mathbf{g}||\mathbf{k}|\,(\sin\theta - \sin\psi)
  \approx -|\mathbf{g}||\mathbf{k}|\cos\theta\,(\psi - \theta).
$$
Suppose the reflection is excited whenever $\psi$ lies within a small total
angular window $\alpha$ around $\theta$. That corresponds to a window in $f$ of
total width
$$
\varepsilon = |\mathbf{g}||\mathbf{k}|\cos\theta\;\alpha .
$$

## 5.2 First-order width

To first order in $\omega - \omega^*$, the time spent inside the window is
$\Delta\omega \approx \varepsilon / |f'(\omega^*)|$:
$$
\boxed{\;\Delta\omega \approx \frac{\alpha\cos\theta}{\sqrt{\sin^2\chi - \sin^2\theta}}\;}
$$

## 5.3 The same result in terms of the detector azimuth

Let $\eta$ be the azimuth of the diffracted beam on the detector, measured
around the incident beam from the rotation axis. The diffracted wavevector
$\mathbf{k}' = \mathbf{k} + \mathbf{g}$ has
$$
k'_z = g_z = |\mathbf{g}|\cos\chi
\qquad\text{and}\qquad
k'_z = |\mathbf{k}|\sin 2\theta\cos\eta ,
$$
so
$$
\cos\eta = \frac{|\mathbf{g}|\cos\chi}{|\mathbf{k}|\sin 2\theta} = \frac{\cos\chi}{\cos\theta}
\quad\Longrightarrow\quad
\sin^2\eta = \frac{\sin^2\chi - \sin^2\theta}{\cos^2\theta}.
$$
Substituting $\sqrt{\sin^2\chi - \sin^2\theta} = \cos\theta\,|\sin\eta|$ into
the boxed result:
$$
\boxed{\;\Delta\omega \approx \frac{\alpha}{|\sin\eta|}\;}
$$

Peaks near $\eta = 0$ (diffracted toward the rotation-axis direction, which is
the same as $\chi \to \theta$) stay in the diffraction condition longest and
spread across the most frames.

## 5.4 Sources of the tolerance alpha

- **Energy bandwidth.** From $\sin\theta = \lambda|\mathbf{g}|/4\pi$,
  $\delta\theta = \tan\theta\,(\delta\lambda/\lambda)$, so
  $\alpha_\lambda \approx \tan\theta\,(\Delta\lambda/\lambda)$.
- **Beam divergence.** Roughly $\alpha_{\text{div}} \approx$ the divergence
  component that changes $\psi$.
- **Mosaic spread.** Zero under the perfect-crystal-per-voxel assumption.

Independent contributions combine roughly in quadrature. The Example2 config
contains a `BeamEnergyWidth` parameter; its units and whether the simulator
uses it have not been checked.

# 6. Consequences

## 6.1 Consistency with the existing Lorentz factor

`XDMEtaAcceptFn` (Python) and `PeakFilters.h` (C++) scale each peak's intensity
by $1/(\sin 2\theta\,|\sin\eta|)$, with $\eta$ computed as
`atan2(|y|, |z|)` of the diffracted direction. With the beam along $x$ and the
axis along $z$, that is the same $\eta$ as in Section 5.3. The $1/|\sin\eta|$
factor matches the dwell time derived here. The $1/\sin 2\theta$ factor is not
derived in this note.

The current forward simulator therefore credits a near-axis peak with the extra
integrated intensity from its long dwell, but deposits all of it in the single
frame containing $\omega^*$.

## 6.2 Distributing a peak across frames

Model the excitation as a box of width $\Delta\omega$ centred on $\omega^*$.
The fraction of the peak's integrated intensity recorded in frame $j$ is
$$
w_j = \frac{\left|\,[\omega^* - \tfrac{\Delta\omega}{2},\,\omega^* + \tfrac{\Delta\omega}{2}]
\;\cap\; [\omega_j - \tfrac{\Delta\omega_f}{2},\,\omega_j + \tfrac{\Delta\omega_f}{2}]\,\right|}{\Delta\omega},
\qquad \sum_j w_j = 1 .
$$
A Gaussian profile is a smoother alternative. The split $\{w_j\}$ carries
information about $\omega^*$ finer than one frame.

## 6.3 Sensitivity of omega* to orientation

A small rotation $\beta$ about a unit axis $\mathbf{n}$ changes $\mathbf{g}$ by
$\delta\mathbf{g} = \beta\,\mathbf{n}\times\mathbf{g}$ and $f$ by
$\delta f = \beta\,\mathbf{k}\cdot(\mathbf{n}\times\mathbf{g})$. The crossing
moves by
$$
\delta\omega^* = -\frac{\delta f}{f'(\omega^*)} ,
$$
which has the same $1/|\sin\eta|$ denominator.

**Worst-case sensitivity.** Writing
$\delta f = \beta\,\mathbf{n}\cdot(\mathbf{g}\times\mathbf{k})$, the rotation
axis that moves $\omega^*$ the most is $\mathbf{n}^* \propto \mathbf{g}\times\mathbf{k}$
(evaluated at $\omega^*$), giving $\max_{\mathbf{n}} |\delta f| = \beta\,|\mathbf{g}\times\mathbf{k}|$.
The Bragg condition $\mathbf{k}\cdot\mathbf{g} = -\tfrac{1}{2}|\mathbf{g}|^2$
fixes the angle between $\mathbf{k}$ and $\mathbf{g}$ at $90^\circ + \theta$,
so $|\mathbf{g}\times\mathbf{k}| = |\mathbf{g}||\mathbf{k}|\cos\theta$.
Dividing by $|f'(\omega^*)| = |\mathbf{k}||\mathbf{g}|\cos\theta\,|\sin\eta|$:
$$
\boxed{\;\max_{\mathbf{n}} \left|\frac{\delta\omega^*}{\beta}\right| = \frac{1}{|\sin\eta|}\;}
$$
radians of $\omega$ per radian of orientation error. For Example2 this ranges
from 1 to about 30, so a $0.1^\circ$ orientation error can move a near-axis
peak by up to about $3^\circ$ (three $1^\circ$ frames).

Near-axis peaks are therefore both the most spread out in $\omega$ and the most
sensitive to orientation, so they carry the most $\omega$ information.

## 6.4 Where the first-order result breaks down

At $\chi = \theta$ ($\eta = 0$), $f'(\omega^*) = 0$ and the first-order width
diverges. Expanding to second order around the tangent point, with
$|f''| \approx |\mathbf{k}||\mathbf{g}|\sin\theta$, gives a finite maximum width
of order
$$
\Delta\omega_{\max} \sim 2\sqrt{\alpha / \tan\theta}.
$$
This is an order-of-magnitude estimate. The first-order formula should be
reliable when $|\sin\eta| \gg \sqrt{\alpha\tan\theta}$.

In Example2 no observable peak comes close to this regime: the smallest
$|\sin\eta|$ is 0.033 ($\eta \approx 1.9^\circ$), because reflections with
$\chi < \theta$ never diffract. At $\Delta E/E = 10^{-3}$ the first-order
width is accurate to $6\times10^{-4}$ even for the closest peaks (Section 7.1).

# 7. Verification

## 7.1 Simulator check (done)

**Setup.** All three voxels of Example2.ThreeVoxels (Cu, 64.351 keV, beam along
$+x$, rotation about $z$), all observable reflections and both branches:
5,100 (reflection, branch) pairs, of which 2,352 are in the voxels' ROI sets.
Bragg angles span $2.65^\circ$–$14.16^\circ$ and $|\sin\eta|$ spans 0.033–1.
Everything was computed in float64. The sample's base rotation
(`sample_to_lab_matrix`) is the identity, so the solver frame and the
simulator's lab frame coincide.

**Method.**

- *Width (Section 5).* An energy change $E \to E(1 \pm \epsilon)$ shifts the
  Bragg angle by $\alpha = \tan\theta\,\epsilon$, so the resulting shift in
  $\omega^*$ (central difference from the solver) should equal
  $\alpha/|\sin\eta|$.
- *Sensitivity (Section 6.3).* Rotate the orientation by a small angle
  $\beta = 10^{-7}$ rad, re-solve, and compare with $-\delta f/f'$ (random
  axis) and with $1/|\sin\eta|$ (worst-case axis $\mathbf{g}\times\mathbf{k}$).
- $\eta$ is computed from the diffracted ray $\mathbf{k}+\mathbf{g}$ as
  `atan2(|y|, |z|)`, the same expression as `XDMEtaAcceptFn`.

**Results (relative error, over all 5,100 pairs).**

| Check | Median | Max |
|---|---|---|
| $f(\omega^*) = 0$ at the solver's $\omega^*$ (normalized by $\lvert\mathbf{k}\rvert\lvert\mathbf{g}\rvert$) | $9\times10^{-17}$ | $1.5\times10^{-15}$ |
| Section 4: $\lvert f'(\omega^*)\rvert = \lvert\mathbf{k}\rvert\lvert\mathbf{g}\rvert\sqrt{\sin^2\chi-\sin^2\theta}$ | $1\times10^{-16}$ | $6\times10^{-14}$ |
| Section 5.3: $\lvert\sin\eta\rvert = \sqrt{\sin^2\chi-\sin^2\theta}/\cos\theta$ | $1\times10^{-16}$ | $6\times10^{-14}$ |
| Section 5.3: $\Delta\omega = \alpha/\lvert\sin\eta\rvert$, $\Delta E/E = 10^{-6}$ | $3\times10^{-10}$ | $1\times10^{-8}$ |
| Section 5.3: $\Delta\omega = \alpha/\lvert\sin\eta\rvert$, $\Delta E/E = 10^{-3}$ | $1\times10^{-6}$ | $6\times10^{-4}$ |
| Section 6.3: random axis, measured vs $-\delta f/f'$ | $4\times10^{-8}$ | $3\times10^{-5}$ |
| Section 6.3: worst-case axis, measured vs $1/\lvert\sin\eta\rvert$ | $7\times10^{-10}$ | $7\times10^{-9}$ |

The largest errors at $\Delta E/E = 10^{-3}$ come from the 8 peaks with
$|\sin\eta| < 10\sqrt{\alpha\tan\theta}$, consistent with the second-order
breakdown in Section 6.4.

**What this does and does not show.** It confirms that the formulas are a
correct first-order description of the Bragg geometry the simulator
implements. It does not show that the formulas describe the real experiment:
that depends on the beam-perpendicular-to-axis geometry holding in practice
and on the true value and sources of $\alpha$.

## 7.2 Still to do

1. **Real-data check.** For peaks with known orientation, plot the number of
   frames each peak spans against $1/|\sin\eta|$. A linear trend confirms the
   model, and its slope estimates $\alpha$. Note that a peak narrower than one
   frame still spans two frames when it straddles a frame boundary, with
   probability roughly equal to its width divided by the frame width.
2. **Literature check.** Compare with the rotation-method Lorentz factor in a
   standard crystallography text.
3. **Config check.** Confirm the units of `BeamEnergyWidth` and whether any
   code path uses it. For scale: for the most spread-out Example2 peak
   ($|\sin\eta| \approx 0.033$) to span a full $1^\circ$ frame requires
   $\alpha \gtrsim 0.033^\circ$; from energy bandwidth alone at
   $\theta \approx 5^\circ$ that is $\Delta E/E \approx 0.7\%$.
