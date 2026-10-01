---
title: "Toy Orientation Network: Design Document"
subtitle: "Data pipeline, learned Gauss-Newton architecture, training, evaluation, tests and current results (branch feature/nn-orientation-arch, through commit 4d48dd1)"
date: "2026-10-01"
geometry: margin=1in
fontsize: 11pt
---

# Summary

The toy orientation network refines the orientation of one voxel given its nominal
orientation $R_{\text{nom}}$ and the thresholded detector data around the spots that
orientation predicts. It outputs an offset $\hat\delta$ (a rotation vector, degrees) and a full
$3\times3$ covariance $\hat\Sigma$. The current architecture, `GNLayerNet`, is a learned
Gauss-Newton layer: the physics (spot-motion Jacobian, exact nominal prediction) is built in,
and a small shared per-peak network only decides how much to trust each measurement.

History and per-step decisions are in `MIGRATION_HISTORY.md` (sections "Toy Orientation NN --
Theory Phase" through "Architecture (parallax)"); the physics the network approximates is in
`docs/nn_inverse_problem_formulation.md` (the "theory note"). Numbers come from saved files in
`benchmarks/toy_orientation_arch/` and `benchmarks/toy_orientation_stage*/`; anything not
checked against a saved file is marked.

Headline results (Section 7; 30 ManyGrains voxels, offsets up to $1^\circ$, simulated data):

- Clean data: `GNLayerNet` is close to a per-voxel Gauss-Newton (GN) fit on the 24 voxels it
  trained on (0.9--1.2$\times$ GN's median error) and 1.1--1.7$\times$ GN on 6 unseen voxels, with GN's
  fall in error with distance from the rotation axis (parallax) and a calibrated covariance. The
  earlier pooled set network was 3$\times$ worse than GN.
- With neighbour/twin spots and pixel noise, plain GN degrades about $20\times$ and Huber-robust GN
  recovers a third of it. The network trained on corrupted windows is 2--3$\times$ better than robust
  GN on distractor data, but not on pixel noise alone; trained on clean data only it is no better
  than robust GN and confidently wrong. There is no real-data result.

# 1. Problem and scope

## 1.1 The task

Local refinement. A voxel has a known nominal orientation $R_{\text{nom}}$ (in practice the coarse-search
hand-off; here the voxel's true orientation in the `.mic`). The unknown is a small rotation vector
$\delta\in\mathbb R^3$, in degrees, in the sample frame (`orientation_eval.offsets_to_matrices`):
$$
R(\delta)=\exp([\delta]_\times)\,R_{\text{nom}},\qquad [\delta]_\times v=\delta\times v .
$$
The network maps the data of the nominal peak set to $(\hat\delta,\hat\Sigma)$, $\hat\Sigma=LL^\top$.
$z$ is the stage (rotation) axis (Example2's base sample rotation is the identity, so sample $z$ = lab $z$;
`error_summary` assumes this); $\perp$ means the $x,y$ components. Pixels constrain $\perp$; frames and
parallax constrain $z$. Symbols: $M$ peak entries per voxel (padded to $M_{\max}=130$), $W=32$ window side (px),
$K=4$ window half-width (frames), $\Delta\omega_f=1^\circ$ frame width, $r_\perp$ distance of the voxel from the
rotation axis, $\omega^*_p(\delta)$ rotation angle at which peak $p$ satisfies the Bragg condition, $B$ batch size.

## 1.2 What the network approximates

The thresholded measurement map is piecewise constant and has no inverse. Squared loss on simulated pairs
drives a network to the posterior mean $E[\delta\mid D]$ under the training prior (uniform in the ball
$\lvert\delta\rvert\le1^\circ$, `sample_prior_offsets`); a Gaussian likelihood loss adds the posterior covariance.
The Bayes risk floors every method (theory note Sections 2, 4, 5). Three results shape the design:

1. Rotating by $\beta$ about $\hat z$ shifts every $\omega^*_p$ by exactly $-\beta$ (3.1.5; tested), so frames carry
   $\delta_z$ as a *global* sum over peaks.
2. A spot's detector position moves through the spot-motion Jacobian $\Gamma_p$ (ring sliding plus parallax, 3.6.3).
   The parallax part of its $z$ column grows with $r_\perp$, so $\delta_z$ is better determined far from the axis.
3. In the independent-quantisation model (frame variance $\Delta\omega_f^2/12$, pixel variance $1/12$ px$^2$) the covariance is
   $(\sum_p J_p^\top J_p)^{-1}$. The noise-free exact-Bayes floor is 13--22$\times$ below it (3.6.7), not a practical target.

## 1.3 Toy assumptions and geometry

| Assumption | Value |
|---------|---------------------------|
| Data | Thresholded lit/unlit pixels from the observer's forward model; rocking width $\alpha=0$ (each peak in one frame; $\alpha>0$ not implemented); no noise except the corruption layer (2.7) |
| Peaks | $Q_{\max}=8$ Angstrom$^{-1}$ (`--max-q 8`); both detectors; only $\lvert\sin\eta\rvert\ge0.3$ at nominal (`--min-sin-eta 0.3`, near-axis peaks drift many frames per degree); identity fixed at nominal (new peaks ignored, vanished peaks give blank windows) |
| Windows | $32\times32$ px about the nominal centroid, $\pm4$ frames about the nominal frame; anything outside is blanked |
| Prior, test | train: uniform ball $\lvert\delta\rvert\le1^\circ$; test: $\lvert\delta\rvert\in\{0.1,0.25,0.5,1.0\}^\circ$, random directions |
| Sample | copper; triangular voxels, ManyGrains side 9.375 $\mu$m (24,570 voxels), ThreeVoxels 0.75--1.5 $\mu$m |

Geometry (`Examples/Example2.ThreeVoxels/ConfigFiles/StandardGeometry.2Det`, `Example2.Simulation.config`; ManyGrains is
identical except `MaxQ`). Near-field HEDM, distances in mm: beam along $+\hat x$ at 64.351 keV; rotation about $\hat z$ in 180 frames
of $1^\circ$ from $-90^\circ$ to $+90^\circ$ (`omega_180_2L.dat`), `EtaLimit 86`; two $2048\times2048$ detectors with $1.48\ \mu$m
pixels at lab $x$ = L1 3.36157 mm (index 0; beam centre J,K = 1022.77, 1925.81 px) and L2 5.38718 mm (index 1; 1020.30, 1921.31 px).
$\Delta L\approx2$ mm is small next to the rotation-axis lever arm, so the detectors are largely redundant for $\delta_z$.
Voxels used: ThreeVoxels voxel 0 ($r_\perp=12\ \mu$m, 113 peaks); ManyGrains voxel 77 ($r_\perp=398.7\ \mu$m, 118 peaks); 30 ManyGrains
voxels over $r_\perp=0$--$501\ \mu$m (105--130 peaks, mean 119).

Data flow (code: `icenine/orientation_{nn,eval}.py`, `toy_orientation_model.py`, `orientation_baselines.py`, `scripts/`, `tests/test_orientation_*.py`):

```
 .mic voxel (R_nom, vertices) + config
   | define_roi_set(detectors="all"); keep |sin eta| >= 0.3
 ROI list: M entries (reflection, branch, detector)
   | BatchedObserver (float64)
   |- peak_context()          -> context (M,16)
   |- nominal_offsets(), pair_index() -> aux sidecar
   |- WindowSpec.from_nominal -> window origin, frame0
   |- render_windows(delta)   -> uint8 windows (N,M,W,W)
   |- render_distractor_windows -> dis_windows (optional)
 dataset .pt -> corrupt_windows -> decode_windows -> x (B,M,2,W,W)
 GNLayerNet(x, context[voxel_id], aux[voxel_id])
   -> delta_hat (B,3) deg, L (B,3,3)
 training: decoupled loss; evaluation: error_summary vs truth,
           GN, Huber GN, exact Bayes
```

# 2. Data pipeline

## 2.1 ROI peak set and the "recorded peak" rule

`define_roi_set(..., detectors="all")` ray-traces the voxel once at nominal. An entry is `ROIPeak(reflection_index,
omega_branch in {1,2}, detector_index, g_hkl, form_intensity, sin_2theta, nominal_omega, nominal_row, nominal_col)`; its
identity and window centre are fixed for the dataset. A (reflection, branch) pair gives entries only if, as in
`ForwardSimulation._simulate_peaks`: (1) it is Bragg-observable with $\omega$ in the scan and passes the eta filter; (2) all three
spot vertices hit *every* detector plane (else it is dropped on all detectors); (3) per detector, the spot overlaps the pixel grid
(`spot_overlaps_grid`: the bounding box of the rasteriser-truncated vertices meets the grid): one entry per such detector.
Rule 3 is the "recorded peak" rule: of 790 peaks hitting voxel 0's detector plane only 364 are ever recorded (the rasteriser clips
the rest); all numbers here use the corrected set. The generator then drops entries below `--min-sin-eta` and rejects voxels with
fewer than `--min-peaks` (40) entries. `detectors="first"` (first detector only) is the Stage 0 behaviour.

## 2.2 `BatchedObserver`

A float64, vectorised re-implementation of the simulator's per-peak maths: `observe(delta_deg (B,3))` returns `present (B,M)`, `frame (B,M)` ($-1$ out of range), `omega (B,M)` (rad) and `verts (B,M,3,2)` ((col,row) on the peak's own detector); also `observe_frames`, `sin_eta`, `vertex_keys`, `lit_pixel_set` (the rasteriser's clip, round and scanline fill). It agrees with the simulator on presence, frame and centroids to $10^{-3}$ px (tested), at about 1 ms per candidate.

## 2.3 Frame-coded windows

`WindowSpec.from_nominal(observer, W, K)` fixes per entry the integer origin
$(\text{col}_0,\text{row}_0)=(\lfloor c_{\text{nom}}\rfloor-W/2,\ \lfloor r_{\text{nom}}\rfloor-W/2)$ and the nominal frame $\text{frame}_0$ (all
ROI spots must be present at nominal). `render_windows(observer, spec, deltas_deg)` returns `windows uint8 (N,M,W,W)` and `status uint8
(N,M)`. Window pixel $(y,x)$ is detector pixel (row $\text{row}_0+y$, col $\text{col}_0+x$); a lit pixel of a spot recorded in frame $f$
holds the code $1+(f-\text{frame}_0)+K\in\{1,\dots,9\}$, unlit is 0. The lit set is the rasteriser's own, so with $\alpha=0$ the encoding
is lossless inside the window. `status`: 0 inside, 1 absent, 2 frame outside $\pm K$, 3 lit pixels partly outside (1 and 2 give an
all-zero window). `decode_windows(w, K)` gives float channels `(..., 2, W, W)`: lit (0/1) and $(w-1-K)/K\in[-1,1]$ at lit pixels, else 0.
A peak is "present" in the network iff its window has a lit pixel, so blank and padded entries coincide (96.8--98.1% of voxel 77's spots lie inside their windows at the $1^\circ$ prior; MIGRATION_HISTORY).

## 2.4 Per-peak context

`BatchedObserver.peak_context(h_deg=0.01)` returns `(M, 14 + n_det)` float32 (16 for two detectors) by central differences about
nominal; spots missing at a probe get zeros in columns 0--8. Column order (read from the code):

| Columns | Content | Units |
|------|------------------|------------|
| 0:6 | $\Gamma$: $\partial(\text{col})/\partial\delta_{x,y,z}$, then $\partial(\text{row})/\partial\delta_{x,y,z}$ | px/deg, over 20 |
| 6:9 | $\partial\omega^*/\partial\delta_{x,y,z}$ (column 8 is $-1$ exactly; tested) | rad/rad |
| 9, 10 | $\lvert\sin\eta\rvert$, $\sin\theta$ | |
| 11:13 | detector one-hot | |
| 13:15 | nominal centroid (col,row) over (ncols,nrows) | fraction |
| 15 | nominal frame index | $2\,\text{frame}/n_{\text{frames}}-1$ |

## 2.5 Exact nominal measurement and pairing

`measurement_features(x, K, nom_off)` maps decoded windows to `(B,M,5)`: `[present, lit count/20, mean frame offset (frames), centroid
col (px), centroid row (px)]`, all zero for absent peaks. Without `nom_off` the centroid is relative to the window centre and the
frame to `frame0`, which carry a fixed per-peak offset (0--1 px, up to $\pm0.5$ frame) because windows are cut on integer pixels.
`nominal_offsets(observer)` `(M,3)` is the exact nominal value: the fractional part of the nominal centroid and
$(\omega_{\text{nom}}-\omega_{\text{centre}}(\text{frame}_0))/\Delta\omega_f\in[-0.5,0.5)$. Subtracting it gives *measurement minus exact nominal
prediction*, the pure quantisation error at $\delta=0$ (tested); the features match `extract_measurements` (the GN baseline's measurements)
up to this offset (tested). `pair_index(roi_list)` `(M,)` links the two detector entries of one diffracted ray (same `reflection_index` and
`omega_branch`; with more detectors, the next detector cyclically), $-1$ if unpaired; about 94% of entries are paired.

## 2.6 Dataset files

`torch.save` dicts, gitignored (`scripts/*.pt`), written as `toy_orientation_<tag>_{train,test}.pt`. Single-voxel files (no `--n-voxels`) hold `windows (N,M,W,W)`, `offsets_deg (N,3)`, meta (`window_size, n_peaks, voxel_index, R_nom, prior_radius_deg, max_q, detectors, min_sin_eta, renderer, frame_half_width, seed, example, r_perp_um, side_um`) and, with `--renderer observer`, `status (N,M)` and `context (M,16)`; test files add `magnitudes_deg`. (`--renderer simulator`, the default there, is the legacy float renderer.) The current multi-voxel format (`--n-voxels V`; shapes checked on `toy_orientation_arch_dis_*.pt`):

| Field | Shape / dtype | Meaning |
|------------|---------|------------------|
| `windows`, `dis_windows` | `(N,130,32,32)` uint8 | frame-coded windows zero-padded to $M_{\max}$; optional distractor layer |
| `offsets_deg`, `voxel_id`, `magnitudes_deg` | `(N,3)` float32, `(N,)` int64, `(N,)` (test only) | targets $\delta$; voxel $0..V-1$ by increasing $r_\perp$; $\lvert\delta\rvert$ bin |
| `context`, `R_nom` | `(V,130,16)`, `(V,3,3)` float32 | per-voxel context table (zero rows for padding); nominal orientations |
| `n_peaks_per_voxel`, `voxel_indices`, `r_perp_um`, `held_out` | `(V,)` | peak counts, mic indices, $r_\perp$ ($\mu$m), held-out flag |
| meta | scalars | `window_size=32, n_peaks=130, prior_radius_deg=1.0, max_q=8.0, detectors="all", min_sin_eta=0.3, renderer="observer", frame_half_width=4, seed=42, example="manygrains", multi_voxel=True` |

Train $N=12{,}000$ (24 voxels $\times$ 500 offsets, 3.2 GB); test $N=1200$ (30 voxels $\times$ 4 magnitudes $\times$ 10). The trainer gathers
`context[voxel_id]` per batch; padding rows are all zero (inert; tested). **Aux sidecar** (`make_dataset_aux.py`): `nom_off (V,M,3)`, `pair_index
(V,M)` ($-1$ none), `det_idx (V,M)`, `frame_width_rad`, `multi_voxel`, `n_peaks`, `voxel_indices` (the trainer asserts the last two match);
the Stage 3 aux file is reused for the distractor data. `toy_orientation_arch_dis_val.pt` (200 samples, 4 training voxels) is used only to choose
the Huber threshold; its generation command is not recorded (unverified).

Voxel selection. `select_voxels` picks, for each of 30 radii evenly spaced in $[0,500]\ \mu$m, up to 60 random candidates (`--voxel-seed 0`); `accept_voxels` takes the first usable one per radius and never reuses a grain (tested with fakes). Sorted by $r_\perp$, every 5th voxel from index 2 is **held out** (indices 2, 7, 12, 17, 22, 27; $r_\perp=32,126,215,301,372,457\ \mu$m): no training samples, 240 of the 1200 test cases.

## 2.7 Distractors and corruption

Both change only the windows, never $\delta$. **Distractor layer** (`--neighbors 2 --twin`; defaults `--neighbor-radius-um 30`, `--neighbor-p 0.5`,
`--neighbor-sigma-deg 0.3`): `build_distractor_sources` takes up to `--neighbors` mic voxels nearest the target in the sample plane (not closer than
1 $\mu$m), each with its own orientation, vertices and ROI set, plus a $\Sigma3$ twin ($60^\circ$ about $[111]$, `R @ T`) of the nearest. Per sample each
source follows the *target's* $\delta$ plus an independent Gaussian offset ($0.3^\circ/\sqrt3$ per component, $0.3^\circ$ rms) and is present with probability
0.5. `render_distractor_windows` draws each present source spot that overlaps a target window (same detector, within $\pm K$ frames) with the same frame
coding; `combine_windows` overlays it *under* the target, whose pixels always win (tested).

Neighbour-model caveat (checked from the `.mic`). In the 30-voxel data all 60 neighbours (2 per voxel) lie at 9.375 $\mu$m (one voxel pitch), and in 45 of 60
the neighbour's orientation equals the target's (same grain, tolerance $10^{-6}$). For those the distractor is the target's own grain seen from an adjacent
voxel, following the target's $\delta$ plus $0.3^\circ$ rms; the other 15 cross a grain boundary; the twin is unrelated by construction. Results hold for this
mixture only.

**Pixel corruption** (`CorruptionConfig`; per sample and entry independently, in this order): `neighbours=True` (overlay `dis_windows`); `p_flip=0.05` (each lit
pixel dropped, each 4-neighbour of a lit pixel lit with its code: edge jitter); `p_hot=0.05` (one isolated hot pixel, random frame code); `p_blob=0.1` (one
$2$--$4\times2$--$4$ px blob, random frame code); `p_miss=0.1` (whole window zeroed). Variants (`CorruptionConfig.named`): `clean`/`none`, `neighbours` (layer only),
`noise` (pixel terms only), `all`. `corrupt_dataset(windows, distractors, name, K, seed=12345, chunk=100)` gives a deterministic copy so all methods see identical
test inputs. Quirk (from the code, not tested): corruption also
hits zero-padded entries, which can then count as "present"; their $J=0$ so they add nothing to $A$ or $b$ and can only perturb the pooled covariance features. The
GN baselines slice windows to the true peak count and do not see them.

# 3. Model

## 3.1 `GNLayerNet`: data flow

```
 x (B,M,2,W,W) -> measurement_features(nom_off) -> y (B,M,3) [sigma units]
 context (B,M,16): Gamma, d omega*/d delta        -> J (B,M,3,3)
 (x, ctx, meas) -> PeakEncoder -> f (B,M,64)
   [pairing: f <- f + MLP([f, f_partner, has_partner])]
 repeat K times (shared head weights), starting at delta' = 0:
   rho = y - J delta'
   head([f, asinh(rho), delta' (, asinh(rho_partner), has)])
        -> w (B,M,3), dy (B,M,3)
   A  = sum_i J_i^T diag(w_i) J_i + ridge I
   delta' = A^-1 sum_i J_i^T diag(w_i) (y_i + dy_i)
 pooled = mean over present peaks of f -> D = diag(exp(g(pooled)))
 delta_hat = s delta' (deg);  Sigma = s^2 D A^-1 D;  L = chol3(Sigma)
```

## 3.2 Physics inputs

For present peaks (absent: $y=0$, $J=0$), with $(\bar c,\bar r)$ the lit centroid (px), $f$ the mean frame offset (frames) and $(c_0,r_0,f_0)$ the exact
nominal values (`nom_off`), the measurement in units of its quantisation sigma ($1/\sqrt{12}$ px, $1/\sqrt{12}$ frame) is
$$
y_i=\sqrt{12}\,\big(\bar c_i-c_{0,i},\ \bar r_i-r_{0,i},\ f_i-f_{0,i}\big)\in\mathbb R^3 .
$$
Its Jacobian with respect to $\delta'=\delta/s$, with unit scale $s=$ `delta_scale` $=0.05^\circ$, is
$$
J_i=\sqrt{12}\,s\begin{bmatrix}\Gamma_i\\ \gamma_i^\top\end{bmatrix}\in\mathbb R^{3\times3},
$$
where $\Gamma_i\in\mathbb R^{2\times3}$ (px/deg) is context columns 0:6 times 20 and $\gamma_i\in\mathbb R^3$ is columns 6:9 in frames per degree (times
$(\pi/180)/\Delta\omega_f$, `frame_width_rad` from the aux file). This is where the $r_\perp$ parallax enters.

## 3.3 Learned weights, pooled normal equations, covariance

For each of the three rows of peak $i$ the encoder and head emit a weight $w_{ir}=\mathrm{softplus}(\cdot)\ge0$ (zero for absent peaks) and a correction
$\Delta y_{ir}$. The estimate solves the weighted normal equations pooled over all peaks and rows:
$$
A=\sum_i J_i^\top\,\mathrm{diag}(w_i)\,J_i+\lambda I,\qquad
\hat\delta'=A^{-1}\sum_i J_i^\top\,\mathrm{diag}(w_i)\,(y_i+\Delta y_i),\qquad \lambda=10^{-3}.
$$
The $3\times3$ inverse (`inv3`, adjugate) and Cholesky factor (`chol3`) are closed form in plain tensor ops, so the layer is differentiable and runs on Apple
MPS (which lacks `linalg.solve`). At initialisation the head's last layer has zero weights and bias $(\mathrm{softplus}^{-1}(1),0)$, so $w=1$, $\Delta y=0$ and the layer *is* one
undamped Gauss-Newton step from nominal, `CentroidGaussNewton.solve_linear` (tested: $2\times10^{-3}$ deg, covariance diagonal within 3%). That step equals converged GN to
within quantisation noise up to $1^\circ$ (median error by $r_\perp$ bin 0.0194/0.0138/0.0120/0.0076 vs 0.0191/0.0127/0.0110/0.0076, `diag_multi_summary.txt`), so learning
only has to down-weight bad peaks and correct biased ones.

Covariance. $A^{-1}$ in $\sigma$ units is the independent-quantisation covariance. The output is
$$
\hat\Sigma=s^2\,D\,A^{-1}D,\qquad D=\mathrm{diag}\big(\exp g(\bar f)\big),
$$
with $\bar f$ the mean of the peak encodings over present peaks and $g$ an MLP ($64\to64\to3$, output layer zero-initialised so $D=I$). $D$ is a learned per-sample, per-axis
calibration: it rescales the standard deviations but keeps the correlation structure that $A$ takes from the geometry. It is needed because quantisation errors are not independent
across peaks (GN's measured $z$ error is $1.5\times$ its predicted $\sigma$, MIGRATION_HISTORY Stage 2) and because $w,\Delta y$ change what $A$ means.

## 3.4 Iterations, pairing, shapes

`n_iter` $=K$ (`--gn-iters`) unrolls IRLS-style rounds with *shared* head weights. Each round the head also sees the residual at the current estimate, $\mathrm{asinh}(y-J\delta')$,
and $\delta'$, so it can down-weight outliers; each round re-solves from the nominal linearisation (the Jacobian is not recomputed). On clean data iterations add nothing; their
purpose is outlier rejection. `pairing=True` (`--pairing`): each entry also receives its other-detector partner's encoding (`aux["pair_index"]`),
$f\leftarrow f+\mathrm{MLP}([f,f_{\text{partner}},\text{has}])$ with the MLP output layer zero-initialised (the network starts unpaired; a first attempt without that was unstable and lr 1e-4 was also needed, MIGRATION_HISTORY), and the head sees
$\mathrm{asinh}$ of the partner's residual and the flag, a pair-consistency check. Both detectors' rows still enter $A$ through their own $J$ rows.

| Tensor | Shape | Notes |
|---------|------------|------------|
| windows $x$; context | `(B,M,2,32,32)` float32; `(M,16)` or `(B,M,16)` | trainer passes `context[voxel_id]` |
| aux `nom_off`; `pair_index` | `(M or B, M, 3)`; `(M or B, M)` long | `pair_index` only with pairing |
| `y`, `w`, `dy`; `J` | `(B,M,3)`; `(B,M,3,3)` | rows col/row/frame, columns $\delta_{x,y,z}$ |
| encoding $f$; head | `(B,M,64)`; in $64+6(+4)$, hidden 64, out 6 | `feat_dim=64` |
| $A$, $A^{-1}$; output | `(B,3,3)`; `delta_hat (B,3)` deg, `L (B,3,3)` | $\hat\Sigma=LL^\top$ |

`PeakEncoder`: windows $\mathrm{conv}(4\to8,3,\text{stride }2)\to\mathrm{conv}(8\to16,3,\text{stride }2)\to\mathrm{fc}(1024\to64)$ (4 channels: lit, frame, two coordinate channels in
$[-1,1]$); context MLP $16\to64\to64$; measurement MLP $5\to64\to64$; concatenation (192) $\to64\to64$; ReLU. `GNLayerNet` has 102,657 parameters (115,393 with pairing; computed),
independent of $M$ and $K$, and is invariant to peak order and zero padding (tested, with and without pairing).

## 3.5 Other architectures still in the code

- `ToyOrientationNet` (`--head quat`): flatten, FC 512/256/128, unit quaternion (191M parameters at 364 peaks); the v0 baseline, kept to re-evaluate under the Stage 0 protocol.
- `ToyOffsetNet` (`--arch fc`): same trunk with $\delta$ and a Cholesky factor, fixed peak count (136M parameters at 130 two-channel peaks); Stage 0--1 reference, strongest single-voxel net at small offsets on voxel 77 (median 0.008--0.009, MIGRATION_HISTORY); unusable on multi-voxel data.
- `PeakSetNet` (`--arch set`): shared conv encoder, context, optional measurement features (`--no-meas`), masked pooling `--pool {meanmax,meansum,all}` (127,681 / 135,873 parameters); the pre-GN comparison point. Mean+max pooling loses $\delta_z$ under NLL (fixed by `--pool all` or the decoupled loss), but on 30 voxels it is $3\times$ worse than GN with no parallax trend. Has its own copy of the `PeakEncoder` layers.
- `FrameProbeNet` (`--arch probe`): per-peak MLP of the mean frame offset and $\partial\omega^*/\partial\delta$ only, mean pooled, identity covariance, `--loss mse` (8,835 parameters); a diagnostic that $\delta_z$ survives mean pooling.

# 4. Training

`scripts/train_toy_orientation_nn.py`. `--arch gn` needs `--aux` and `--head offset`; `--arch set` needs `--aux --subpixel` to use the exact nominal measurement.

**Losses** (`--loss`; `icenine/orientation_nn.py`): `nll` (default; `gaussian_nll_loss` $=\tfrac12\lVert L^{-1}(\delta-\hat\delta)\rVert^2+\sum\log L_{ii}$, `--beta-nll b>0` weights each sample by
$\mathrm{stopgrad}(\det\Sigma^{b/3})$); `mse` (`mse_deg_loss` $=\tfrac12\sum_{\text{axes}}(\delta-\hat\delta)^2/s_m^2$, $s_m=$ `--mse-scale` $=0.1^\circ$); `decoupled` (`decoupled_nll_loss`: `mse` on $\hat\delta$ plus the NLL
at $\mathrm{stopgrad}(\hat\delta)$, so the mean gets exactly the MSE gradient and $\hat\Sigma$ exactly the NLL gradient; tested); `mse-then-cov`/`mse-then-nll` (`mse` for the first half of the epochs, then `decoupled` or `nll`;
only second-half epochs can be the best-validation checkpoint). *The NLL pathology:* the NLL gradient on the mean is $\Sigma^{-1}(\hat\delta-\delta)$, so the $z$ gradient is smaller than the $\perp$ gradient by
$(\sigma_\perp/\sigma_z)^2$ ($1/625$ for $0.02^\circ$ vs $0.5^\circ$; tested); once $\sigma_z$ is large the network stays at a calibrated "I do not know $z$". This stalled the pooled set network on voxel 0; changing only the loss to
`decoupled` fixed it, while per-sample $\beta$-NLL did not (one scalar cannot rebalance axes within a sample). All current `GNLayerNet` runs use `decoupled`.

**Optimiser.** Adam, `--lr` (3e-4 for the first GN runs, 1e-4 for Step 3/4), `--batch-size 64`, `--epochs 60`, `--clip 1.0` (gradient norm), `--cosine` (decay to 0). `--ema d` keeps a per-step exponential average of
the weights ($d=0.998$; buffers copied) and validation is evaluated on it. `--checkpoint`: `best` (default; best-validation epoch), `ema` (EMA weights at the last epoch), `last`. Current runs use `--ema 0.998 --checkpoint ema`:
validation loss on four voxels is noisy and best-epoch selection picked early epochs (12--20 of 60) that were worse everywhere (MIGRATION_HISTORY, Stage 3 retrain).

**Validation split.** Default: a random `--val-frac` (10%) of training *samples*, in-distribution, blind to degradation on new voxels. `--val-voxels N` (multi-voxel): `split_by_voxel` sorts the training voxels by $r_\perp$,
cuts $N$ equal strata and draws one voxel per stratum with `--seed`; their samples leave training and are used only for selection and monitoring. With $N=4$: seed 0 gives voxels 6, 11, 19, 24 ($r_\perp=99,181,326,410\ \mu$m; 10,000
train and 2,000 validation samples), seed 1 gives 3, 11, 20, 29 (computed with `split_by_voxel`; MIGRATION_HISTORY quotes only seed 0). The six held-out test voxels are never used for selection; the validation voxels are among the 24
"in-dist" test voxels, so that group contains unseen-voxel cases that differ by seed.

**Devices.** `--device {cpu,mps}`: network training and inference in float32 (headline runs on the Apple GPU: about 5 s per epoch clean, 15 s with corruption); data stay on the CPU, each batch is corrupted and decoded there
then moved. Always float64 on the CPU: the observer, `ExactBayes`, the GN baselines and all error statistics. A one-epoch CPU/MPS check (set net, seed 0) gave losses $-0.87531$ vs $-0.87529$.

**Corruption-aware training.** `--corrupt-train {none,neighbours,noise,all}`: each training batch is corrupted on the fly (`corrupt_windows`, generator `--seed + 7`); validation windows are corrupted once (`--seed + 99`) so the
validation loss is comparable across epochs; `neighbours` and `all` need `dis_windows`. Test variants: `--eval-variants clean,neighbours,noise,all` (deterministic, seed 12345).

Headline runs use `--head offset --loss decoupled --device mps --clip 1 --cosine --batch-size 64 --epochs 60 --ema 0.998 --checkpoint ema --val-voxels 4 --arch gn`, seeds 0 and 1 (Section 9; the clean runs' flags are as stated in MIGRATION_HISTORY, the result files do not store them).

# 5. Evaluation protocol

**Test sets.** Fixed $\lvert\delta\rvert\in\{0.1,0.25,0.5,1.0\}^\circ$, random directions (`sample_fixed_magnitude_offsets`), 10 per voxel and magnitude: 1200 cases over 30 voxels, 960 on the 24 training-set voxels ("in-dist": seen *voxel*,
new offsets) and 240 on the 6 held-out voxels. Single-voxel sets have 30 cases per magnitude. (The float32 value of $0.1^\circ$ prints as `0.10000000149011612`, the key in the result jsons.)

**Metrics.** `error_summary(delta_hat, delta_true)` returns, in degrees, `rms_z`, `rms_perp` ($\sqrt{\text{mean}((e_x^2+e_y^2)/2)}$), `rms_total`, `median_angle` (median rotation angle between $\exp([\hat\delta]_\times)$ and $\exp([\delta]_\times)$),
`success_0p5` (misorientation $<0.5^\circ$, the `bench_hp_sweep.py` criterion) and `success_0p1` ($<0.1^\circ$, used from Stage 1); errors are rotation-vector differences (accurate to $1^\circ$). *Calibration:* Mahalanobis$^2$
$=(\delta-\hat\delta)^\top\hat\Sigma^{-1}(\delta-\hat\delta)$ averaged over cases; 3 for a calibrated Gaussian, $\gg3$ overconfident (1500 = confidently wrong), $<3$ underconfident. *Parallax test:* per voxel the median angle over its 40 test
cases; the Pearson correlation with $r_\perp$ over the 30 voxels (`corr_err_rperp`) and the median per-voxel error in-dist and held-out. GN on clean data has corr $-0.63$ (error falls as parallax grows); a method not using parallax has corr $\ge0$.

**Baselines.**

- *predict-nominal*: $\hat\delta=0$ (median angle $=\lvert\delta\rvert$).
- *Exact Bayes* (`ExactBayes`, `exact_bayes_baseline.py`, seed 7): posterior of $\delta$ given exact noise-free data (presence, frame, rasteriser pixel sets) by adaptive importance sampling on the consistent set; its mean is the best squared-error estimate and
  $\sqrt{\mathrm{tr}\,\mathrm{cov}}$ the floor. Limits: single-voxel datasets only (7 min per 120 cases, not run on the 30-voxel data); ignores newly appearing peaks and spot overlaps; uses exact pixels, not the network's windows (a generous reference).
- *Plain GN* (`CentroidGaussNewton`, `gauss_newton_baseline.py`): Levenberg-Marquardt on three residuals per spot (frame-centre $\omega$ with $\sigma=\Delta\omega_f/\sqrt{12}$; centroid col, row with $1/\sqrt{12}$ px), the observer as exact model,
  central-difference Jacobians (`fd_step_deg=0.002`), start at nominal, `max_iter=30`, `tol_deg=1e-7`; returns `delta`, `cov` $=(J^\top WJ)^{-1}$, `n_used`, `status`, $\chi^2$/dof. Limits: it is given the peak
  identities and the basin and sees only per-spot centroid and frame; about 7 ms per case.
- *Huber GN* (`huber_c`, `--huber c`): the same with Huber loss on normalised residuals via IRLS (weight $\min(1,c/\lvert r\rvert)$); one start, one robust loss. $c$ was chosen on a validation split, not the test set: 200 prior-ball samples from the four validation voxels
  (6, 11, 19, 24), windows corrupted with `all` (`step4_huber_val.txt`). Median angle: plain GN 0.0945; Huber $c=0.5/1/2/3/5$: 0.0336 / **0.0325** / 0.0328 / 0.0369 / 0.0426 (on clean validation windows plain 0.0135, Huber 0.0179 / 0.0158 / 0.0125 / 0.0134 /
  0.0135). $c=1$ is used for all variants, including clean data (about 20% worse held-out median: 0.0144 vs 0.0119). With `--max-iter 200` all fits converge and the median is unchanged (MIGRATION_HISTORY).
- *MC and Riemannian Adam* (`optimizer_baselines.py`): the `bench_hp_sweep.py` protocol on the Python-simulated ThreeVoxels images generated at voxel 0's truth; each optimiser starts at $\exp([\delta]_\times)R_{\text{nom}}$ and searches for the truth. MC: 3500 steps, 2
  restarts, step fraction 0.5, search box $1.5\lvert\delta\rvert$ (it is told $\lvert\delta\rvert$). Adam (geoopt): scale 2, $\omega$ window 1, lr $10^{-4}$, 100 steps. $\lvert g\rvert\le8$. Limits: voxel 0 only (ManyGrains has no images); data fixed and start displaced (equivalent to
  first order to the networks' setting, not identical); a harder spot-association problem; one run, 30 cases per bin; never run on corrupted data.

**Seeds.** Dataset: generator `--seed 42` (one stream), `--voxel-seed 0`, distractor draws `seed+1000`. Trainer `--seed s` sets initialisation, batch order and the validation-voxel draw (headline nets $s\in\{0,1\}$); training/validation corruption `s+7`/`s+99`; test
corruption fixed 12345 for all methods; exact Bayes `--seed 7`; MC `seed + case index`. Network results are 2-seed means (per-seed values where tabulated); one dataset draw, no confidence intervals; GN is deterministic. Seed spread of the held-out median reaches
$0.003^\circ$ clean and $0.04^\circ$ for the clean-trained net on `all`.

# 6. Test matrix

86 tests pass in 23 s (`uv run pytest tests/test_orientation_*.py`, 2026-10-01). **nn** = `test_orientation_nn.py`, **ev** = `test_orientation_eval.py`, **bl** = `test_orientation_baselines.py`; the `test_` prefix of method names is dropped; "phys" = needs
ThreeVoxels voxel 0 (skipped if absent; every physics test uses that voxel).

| Test (file::class::name) | Property guaranteed |
|---------------|---------------|
| nn::TestDefineROISet::{roi_set_nonempty, roi_peaks_have_valid_branch_and_detector, roi_set_deterministic} | non-empty set, branch in {1,2}, valid detector; redefining reproduces identities and restores sample state (phys) |
| nn::TestRenderLocalWindows::{window_shape, nominal_orientation_has_no_missing_peaks, sample_rotation_state_restored, roi_stability_across_perturbation_range} | legacy renderer: shapes, all peaks at nominal, state restored, $<50\%$ dropout at $1^\circ$ |
| nn::TestSampleLocalPerturbations::{returns_valid_rotation_matrices, perturbations_stay_close_to_nominal} | legacy sampler: proper rotations, within $10^\circ$ |
| nn::TestLossAndMetric (4), TestToyOrientationNet::output_is_unit_quaternion | quaternion loss zero at identity, sign-invariant, finite gradient; misorientation matches `orientation_search`; unit output |
| nn::split_by_voxel_disjoint_deterministic_and_spans_r_perp | validation voxels disjoint, seed-deterministic, one per $r_\perp$ stratum; bad `n_val` raises |
| ev::TestRotations (6) | Rodrigues equals scipy; offset-quaternion round trip; $R=\exp([\delta]_\times)R_{\text{nom}}$ acts in the sample frame; ball prior; fixed-magnitude sampler |
| ev::TestErrorSummary (3) | zero error; success fractions; $z$/$\perp$ separation |
| ev::TestPixelGrid::{spot_overlaps_grid_rule, lit_pixel_set_matches_the_rasteriser, roi_peaks_are_recorded_peaks} | recorded-peak rule; `lit_pixel_set` equals the rasteriser on 400 random triangles incl. clipped; every ROI peak on the grid at nominal (phys) |
| ev::TestOffsetHead (7) | Cholesky lower-triangular, positive diagonal; NLL equals `MultivariateNormal`; $\beta$-NLL finite; decoupled mean gradient = MSE gradient (independent of $\Sigma$), covariance gradient = NLL gradient; NLL mean gradient $\propto\Sigma^{-1}$; `ToyOffsetNet` shapes |
| ev::TestBatchedObserver::matches_simulator[first,all] | presence, frame, centroids ($10^{-3}$ px) equal the serial simulator at 4 offsets (phys) |
| ev::TestBatchedObserver::{all_detectors_mode, invalid_detectors_mode, missing_any_detector_plane_drops_peak_everywhere} | `"all"` extends `"first"`, all present at nominal; invalid mode raises; drop-everywhere rule |
| ev::TestBatchedObserver::{max_q_override_limits_reflections, nominal_offset_reproduces_roi_set, sample_state_untouched} | `max_q` limits reflections; ROI present at nominal; sample not mutated |
| ev::TestBatchedObserver::rotation_about_stage_axis_shifts_every_omega_exactly | $\omega^*\to\omega^*-\beta$ for a $z$ rotation, to $10^{-9}$ |
| ev::TestExactBayes::{posterior_contains_truth_and_is_tight, pixels_tighten_the_posterior} | truth inside, $\sigma<0.005^\circ$, ESS $>100$; pixels tighten it vs frames only |
| ev::TestRenderWindows::{decode_windows, windows_match_lit_pixel_sets, sin_eta_near_axis} | decoding; window pixels = rasteriser pixels with code $1+f-f_0+K$ and correct status; $\lvert\sin\eta\rvert$ range |
| ev::TestEndToEndVsForwardSimulation | windows at nominal equal `ForwardSimulation` thresholded images pixel for pixel; at an offset nothing spurious, $<2\%$ missing |
| bl::TestPeakSetNet (7) | shapes, positive-definite covariance, parameter count independent of $M$, permutation invariance, absent peaks ignored, per-sample context gather, finite with no peaks, gradients |
| bl::TestMeasurementFeatures (5) | feature values; zero frame feature without frame channel; padding leaves `PeakSetNet` unchanged (3 pools) |
| bl::TestMeasurements (2), TestMeasurementFeaturesMatchExtraction | `frame_center_omega`; extraction consistent with the observer; features equal `extract_measurements` (phys) |
| bl::TestGaussNewton (3), TestGaussNewtonStatus | GN far better than nominal at $0.6^\circ$ ($\perp<0.02$, $z<0.12$); small at nominal; $\sigma_z>3\sigma_\perp$; status consistent with `converged` (phys) |
| bl::TestPeakContext, TestFrameProbeNet | context shape (M,16), column 8 $=-1$, $\lvert\sin\eta\rvert\ge0.3$, one-hot (phys); probe uses only frame and $\partial\omega^*/\partial\delta$ |
| bl::TestVoxelSelection (3) | radius selection; a grain is never reused; a rejected voxel's grain stays available (synthetic mics) |
| bl::TestGNLayer::{normal_equations_match_weighted_lstsq, unit_weights_reproduce_linear_gauss_newton_step} | `gn_normal_equations`/`inv3` equal weighted least squares; $w=1,\Delta y=0$ gives `solve_linear` and covariance diagonal within 3% (phys) |
| bl::TestGNLayer::{permutation_and_padding_invariance[False,True], gradients_flow_and_are_finite} | order and padding invariance with and without pairing; finite gradients ($K=2$, pairing) |
| bl::TestPairingAndNominalOffsets (4) | `pair_index` symmetric, same reflection/branch/$\omega$, other detector, unpaired and 3-way cases; features minus `nominal_offsets` = measurement minus exact nominal; detector mask partitions entries (phys) |
| bl::TestDistractors (3), TestRobustGaussNewton | identical source reproduces target; displaced neighbour adds pixels, never alters target pixels or codes; corruptions deterministic and bounded; Huber with huge $c$ equals plain GN and resists 20 gross outliers (phys) |

**Not tested:**

- The trainer end to end (loss selection, EMA, checkpoint choice, two-phase schedule, `--val-voxels`, corruption-aware training, MPS vs CPU); only `split_by_voxel` and the loss functions are unit tested.
- The generator's multi-voxel path (`main_multi`), `build_distractor_sources`, held-out assignment, dataset fields and shapes, `make_dataset_aux.py` and aux/dataset consistency; `gauss_newton_baseline.py`, `exact_bayes_baseline.py`, `optimizer_baselines.py`, the `summarize_*` scripts.
- Physics on ManyGrains voxels or at large $r_\perp$: all physics tests use ThreeVoxels voxel 0 ($r_\perp=12\ \mu$m), so the parallax columns of the context are untested far from the axis.
- Accuracy and calibration of any trained network (no regression test on result numbers).
- `GNLayerNet` beyond equivalence at initialisation: outlier rejection for $K>1$, that the zero-initialised pair MLP equals the unpaired network, off-diagonals of `chol3`, the learned $D$.
- Corruption statistics (rates, blob sizes, code ranges), corruption of padded entries, multi-source/twin/`source_active` distractors with real neighbours.
- `ExactBayes` calibration (measured ad hoc in Stage 0) and three or more detectors; other geometries, $\alpha>0$, other noise models.

# 7. Current results

Setup unless stated: 30 ManyGrains voxels, 24 train (4 used for validation) and 6 held out; train offsets in the $1^\circ$ ball; test $\lvert\delta\rvert=0.1/0.25/0.5/1.0^\circ$, 10 offsets per voxel and magnitude; `GNLayerNet` trained 60 epochs, decoupled loss, EMA
weights, two seeds averaged; medians in degrees; one dataset draw, no confidence intervals.

## 7.1 Clean data

Median angle at the four bins, from `step2_summary.txt`, `step3_summary.txt`, `summarize_results.py` on `multi_res_set_subpx_s*.json`, and `step4_summary.txt` (plain GN).

| Method (2 seeds) | In-dist | Held-out | corr(err, $r_\perp$) | Voxel median err in / held | Mah$^2$ in / held |
|------------|---------|---------|------|------|------|
| `PeakSetNet`, `--subpixel` | .033/.035/.033/.043 | .029/.034/.039/.053 | +0.29 | .0354 / .0362 | 3.1--5.1 / 4.0--9.3 |
| `GNLayerNet` $K=1$, lr 3e-4 | .013/.015/.012/.013 | .017/.016/.015/.013 | $-0.72$ | .0136 / .0188 | 3.1--3.2 / 3.5--4.0 |
| `GNLayerNet` $K=3$, lr 3e-4 | .013/.015/.012/.015 | .018/.016/.017/.016 | $-0.64$ | .0131 / .0169 | 3.1--3.4 / 3.8--4.2 |
| `GNLayerNet` $K=3$, lr 1e-4 | .013/.016/.013/.013 | .016/.016/.015/.015 | $-0.61$ | .0138 / .0163 | 2.8--3.2 / 3.5--3.7 |
| `GNLayerNet` $K=3$ paired, lr 1e-4 | .011/.014/.012/.014 | .013/.015/.016/.016 | $-0.60$ | .0143 / .0145 | 2.8--3.2 / 3.9--4.5 |
| plain GN | .012/.015/.012/.014 | .014/.012/.012/.010 | $-0.63$ | .0127 / .0143 | -- |

**Corrected parity statement.** The Step 2 commit subject (e6f4072) says the net "reaches Gauss-Newton parity"; that overstates the held-out result. Ratio of the net's median angle to GN's (mean of the four bins), per seed, from the saved jsons:

| Run | In-dist, seed 0 / 1 | Held-out, seed 0 / 1 | corr(err, $r_\perp$), seed 0 / 1 |
|---------|---------|---------|---------|
| $K=1$, lr 3e-4 | 0.93 / 1.12 | 1.15 / 1.43 | $-0.64$ / $-0.79$ |
| $K=3$, lr 3e-4 | 1.17 / 0.91 | **1.74** / 1.10 | $-0.43$ / $-0.85$ |
| $K=3$, lr 1e-4 | 0.92 / 1.19 | 1.24 / 1.38 | $-0.75$ / $-0.47$ |
| $K=3$ paired, lr 1e-4 | 0.93 / 1.04 | 1.13 / 1.40 | $-0.70$ / $-0.50$ |

So: close to GN on seen voxels (0.9--1.2$\times$), 1.1--1.7$\times$ GN on unseen ones, with GN's parallax trend (negative correlation in every run, strength seed dependent); the set network's $3\times$ gap is removed in every run. Pairing neither helps nor hurts on clean data beyond seed spread. Single-voxel checks (30 cases per bin, one run): on voxel 0 ($r_\perp=12\ \mu$m) `GNLayerNet` $K=1$ has median .0142/.0176/.0184/.0249 at the four bins vs GN .0401/.0242/.0348/.0362 and exact Bayes .0069/.0040/.0028/.0040; on voxel 77 ($r_\perp=399\ \mu$m) .0098/.0078/.0086/.0095 vs GN .0118/.0106/.0094/.0079 and Bayes .0006--.0009 (`single_v0_gn.json`, `single_v77_gn.json`, `far_res_bayes.json`). MC and Riemannian Adam on voxel 0 (from `benchmarks/toy_orientation_stage2/pred_*.npz`) reach medians .048/.145/.360/.672 and .056/.131/.185/.523 and leave $\delta_z$ largely uncorrected on these near-axis voxels (their 96%/92% success at $1^\circ$ in the ManyGrains sweep was a different, far-from-axis sample; the $r_\perp$ dependence was not measured).

## 7.2 Corrupted data (Step 4)

Dataset `toy_orientation_arch_dis_*` (seed 42; 12,000/1,200 samples; clean windows identical to the Stage 3 data). Median angle pooled over the four bins (it does not depend on $\lvert\delta\rvert$ in any row; per-bin values are in `step4_summary.txt`). "Clean-trained":
the same architecture trained on the clean windows of this file; "corr-trained": `--corrupt-train all`; Huber $c=1$; nets are 2-seed means, held-out per-seed values in brackets; $<0.1^\circ$ is the held-out fraction; Mah$^2$ is in-dist / held-out.

| Test set | Method | Median in-dist | Median held-out [seeds] | $<0.1^\circ$ | Mah$^2$ |
|------|---------|---------|---------------|------|---------|
| clean | plain GN | .0131 | .0119 | 1.00 | |
| | Huber GN | .0125 | .0144 | 1.00 | |
| | net, clean-trained | .0138 | .0152 [.0143, .0161] | 1.00 | 3.0 / 3.6 |
| | net, corr-trained | .0185 | .0193 [.0176, .0210] | 1.00 | 2.1 / 2.2 |
| | paired, corr-trained | .0185 | .0186 [.0185, .0186] | 1.00 | 2.3 / 2.3 |
| neighbours | plain GN | .2253 | .3428 | 0.20 | |
| | Huber GN | .1375 | .2223 | 0.33 | |
| | net, clean-trained | .1934 | .2706 [.2916, .2496] | 0.23 | 1556 / 1951 |
| | net, corr-trained | .0576 | .0719 [.0723, .0715] | 0.66 | 3.0 / 2.9 |
| | paired, corr-trained | .0594 | .0710 [.0702, .0718] | 0.65 | 3.0 / 2.7 |
| noise | plain GN | .0545 | .0573 | 0.85 | |
| | Huber GN | **.0177** | **.0215** | 1.00 | |
| | net, clean-trained | .0704 | .0741 [.0746, .0736] | 0.71 | 223 / 272 |
| | net, corr-trained | .0241 | .0242 [.0227, .0258] | 0.99 | 3.0 / 3.0 |
| | paired, corr-trained | .0233 | .0231 [.0228, .0233] | 1.00 | 3.0 / 3.0 |
| all | plain GN | .2228 | .3311 | 0.20 | |
| | Huber GN | .1405 | .2239 | 0.32 | |
| | net, clean-trained | .2033 | .2757 [.3016, .2498] | 0.16 | 1462 / 1796 |
| | net, corr-trained | **.0615** | **.0821** [.0827, .0815] | 0.64 | 2.9 / 2.6 |
| | paired, corr-trained | **.0616** | **.0762** [.0734, .0789] | 0.62 | 2.8 / 2.5 |

1. Distractor spots break GN (about $20\times$ worse than clean, 20--34% of cases within $0.1^\circ$); Huber GN removes about a third of the error and loses the parallax trend (corr $+0.13/+0.14$ vs $-0.63$ clean). Pixel noise alone costs plain GN $4\times$; Huber GN recovers nearly all of it.
2. The corruption-trained net is 2--3$\times$ better than Huber GN on `neighbours` and `all`, with calibrated covariance (Mah$^2$ 2.5--3.0, held-out too); on `noise` alone Huber GN is as good or better (0.0215 vs 0.0242; every net seed is worse). Errors are heavy-tailed: held-out $\perp$ RMS 0.05
   vs 0.004 clean, 34--38% of held-out cases above $0.1^\circ$.
3. Corruption-aware training, not the iterative head, is what matters: the clean-trained net is no better than plain GN on distractors, worse on `noise` (0.074 vs 0.057) and wildly overconfident; corruption-aware training costs clean-data accuracy (held-out 0.0176--0.0210 vs 0.0143--0.0161 vs 0.0119 GN).
4. Held-out voxels are harder under distractors for every method (held-out/in-dist: net 1.3$\times$, plain GN 1.5$\times$, Huber GN 1.6$\times$); the $r_\perp$ trend is lost for all (corr $-0.1$ to $-0.18$ for plain GN and nets).
5. Pairing has no clear benefit (differences $\le0.001^\circ$ except held-out `all`: paired .0734/.0789 vs unpaired .0827/.0815, 2 seeds, 6 voxels: suggestive). Untested hypothesis: a neighbour spot lands consistently on both detectors, so the partner check cannot reject it.

Caveats. Two seeds, six held-out voxels; seed spread of the held-out median up to $0.003^\circ$ clean, $0.0055^\circ$ (corr-trained nets on `all`), $0.04^\circ$ (clean-trained net on `all`). Train and test corruptions are the same family with the same neighbour sets per voxel (new random draws, not new kinds of nuisance), so the
advantage over robust GN may shrink out of family; the neighbour caveat of 2.7 applies (75% same-grain). Robust GN is one fixed construction. Learning rate and $K$ were not re-tuned for corrupted data. What the weight head learned has not been analysed.

# 8. Known limitations, open questions, roadmap

- **Simulation only.** No real data; $\alpha=0$, known peak identities, nominal orientation within $1^\circ$. Real near-field HEDM resolution is about $0.1^\circ$ (theory note); the toy's clean errors ($0.01^\circ$) are noise-free quantisation only and the Bayes floor ($10^{-3}$ deg) is not reachable.
- **Fixed peaks and windows.** Peaks absent at nominal are invisible; about 2--3% of spots at the $1^\circ$ prior are not fully inside their windows (voxel 77); near-axis peaks are dropped; windows are not sized from the Jacobian. The layer linearises at nominal (negligible up to $1^\circ$, no re-linearisation for larger priors).
- **Generalisation and calibration.** Six held-out voxels from one sample, geometry and structure; held-out error is 1.1--1.7$\times$ in-dist on clean data. $\hat\Sigma$ is calibrated on this distribution only: the older set network was 2--3$\times$ overconfident on unseen voxels (Mah$^2$ 17--34) and the clean-trained GN net is
  catastrophically overconfident under corruption.
- **Corruption realism.** Same-family train/test, at most 3 sources, neighbours mostly same-grain (2.7), no intensity effects.
- **No integration with `FindOptimal`.** The coarse-search hand-off sometimes selects a wrong ~54$^\circ$ solution (neighbour voxel or $\Sigma3$ twin); a refiner should be evaluated per candidate.

Open questions: what the weight head learns on distractor data and why pairing does not reject neighbour spots; out-of-family corruption sweeps and per-sample severity as augmentation; more seeds and confidence intervals; recalibration of $\hat\Sigma$ on held-out voxels; the $\sim1.2\times$ held-out gap to GN on clean data; the heavy $\perp$ tail
on distractor data.

Roadmap (all unimplemented): per-peak second moments, an $\alpha>0$ renderer with soft per-pixel/per-frame inputs, a Sinkhorn render-and-compare loss and posterior scoring against exact Bayes (KL, coverage) -- `docs/todo_intensity_and_distribution_losses.md`; a soft ROI renderer, iNeRF-style refinement after the network, BARF-style annealing
and self-supervised training on real data (parked) -- `docs/research_ideas_nerf.md`; why the coarse search and `FindOptimal` return neighbour or twin orientations, and per-candidate evaluation -- `docs/todo_findoptimal_wrong_candidates.md`; rocking width $\alpha$ for real data (decisions D1, D4) -- `MIGRATION_HISTORY.md`, Stage 1 decisions.

# 9. How to reproduce

From `icenine_py/`, always `uv run python`. Scripts change directory (the physics setup `chdir`s to the example), so use paths relative to `icenine_py/` or absolute. Flags below exist in the argparse definitions; where a result file does not store its command that is stated. Datasets are gitignored and large (3.2 GB).

```bash
cd icenine_py
G="--example manygrains --max-q 8 --detectors all --min-sin-eta 0.3 \
   --frame-half-width 4 --prior-radius 1 --n-voxels 30 --voxel-seed 0 \
   --test-magnitudes 0.1 0.25 0.5 1.0 --per-voxel-train 500 \
   --per-voxel-test 10"
D=benchmarks/toy_orientation_arch
S=scripts/toy_orientation

# 1. datasets (clean + aux sidecar; with distractor layer, seed 42 default)
uv run python scripts/generate_toy_orientation_dataset.py $G \
    --tag stage3_multi
uv run python scripts/make_dataset_aux.py \
    --data ${S}_stage3_multi_test.pt --out ${S}_stage3_multi_aux.pt
uv run python scripts/generate_toy_orientation_dataset.py $G \
    --neighbors 2 --twin --tag arch_dis

# 2. GN and Huber-GN on each test variant (clean = --corrupt none)
for V in clean neighbours noise all; do
  C=$V; [ "$V" = clean ] && C=none
  uv run python scripts/gauss_newton_baseline.py \
      --test ${S}_arch_dis_test.pt --corrupt $C \
      --out $D/dis_pred_gn_$V.npz
  uv run python scripts/gauss_newton_baseline.py \
      --test ${S}_arch_dis_test.pt --corrupt $C --huber 1 \
      --out $D/dis_pred_huber1_$V.npz
done

# 3. nets: this is the corruption-trained net. Drop --corrupt-train for
#    the clean-trained control; add --pairing for the paired net; seeds 0
#    and 1. The clean multi-voxel nets use stage3_multi_{train,test}.pt,
#    --lr 3e-4 or 1e-4, --gn-iters 1 or 3 (flags as in MIGRATION_HISTORY)
uv run python scripts/train_toy_orientation_nn.py --head offset \
    --arch gn --gn-iters 3 --loss decoupled --device mps --lr 1e-4 \
    --clip 1 --cosine --batch-size 64 --epochs 60 --ema 0.998 \
    --checkpoint ema --val-voxels 4 --corrupt-train all \
    --eval-variants clean,neighbours,noise,all --seed 0 \
    --train ${S}_arch_dis_train.pt --test ${S}_arch_dis_test.pt \
    --aux ${S}_stage3_multi_aux.pt \
    --results-json $D/dis_res_gn_k3_corr_s0.json \
    --save-predictions $D/dis_res_gn_k3_corr_s0.npz

# 4. tables
uv run python scripts/summarize_results.py "K3=$D/multi_res_gn_k3_s?.json"
uv run python scripts/summarize_arch_step4.py \
    --test ${S}_arch_dis_test.pt --huber 1 --seeds 0 1 \
    --runs "corr-trained=dis_res_gn_k3_corr" \
    "paired=dis_res_gnpair_k3_corr" \
    "clean-trained=dis_res_gn_k3_cleantrain" \
    --out $D/step4_summary.json
```

Single-voxel baselines (exact Bayes; MC and Adam also need the Python-simulated ThreeVoxels images; the dataset command is as in MIGRATION_HISTORY, not re-run):

```bash
uv run python scripts/generate_toy_orientation_dataset.py \
    --renderer observer --max-q 8 --detectors all --min-sin-eta 0.3 \
    --frame-half-width 4 --prior-radius 1 --n-train 10000 \
    --test-magnitudes 0.1 0.25 0.5 1.0 --tag stage1
uv run python scripts/exact_bayes_baseline.py --test ${S}_stage1_test.pt
uv run python scripts/optimizer_baselines.py --test ${S}_stage1_test.pt \
    --out-dir benchmarks/toy_orientation_stage2
```

Unverified: the exact commands for `toy_orientation_arch_dis_val.pt`, for the Huber-validation GN predictions under `benchmarks/toy_orientation_arch/val/`, and for the Stage 1 single-voxel set (tag and sample count follow MIGRATION_HISTORY; not re-run). Tests: `uv run pytest tests/test_orientation_nn.py tests/test_orientation_eval.py tests/test_orientation_baselines.py`.
