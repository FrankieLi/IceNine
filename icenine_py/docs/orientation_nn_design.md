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
`docs/nn_inverse_problem_formulation.md` (the "theory note", cited as "TN §"). Numbers come from saved files in
`benchmarks/toy_orientation_arch/` and `benchmarks/toy_orientation_stage*/`; anything not
checked against a saved file is marked.

Headline results (Section 7; 30 ManyGrains voxels, offsets up to $1^\circ$, simulated data):

- Clean data: `GNLayerNet` is close to a per-voxel Gauss-Newton (GN) fit on the 24 training-set
  voxels (20 trained on, 4 used only for validation; which 4 depends on the seed): 0.9--1.2$\times$
  GN's median error, and 1.1--1.7$\times$ GN on 6 unseen voxels. It shows GN's fall in error with
  distance from the rotation axis (parallax), and its mean Mahalanobis$^2$ is near 3 (coverage not
  checked). The earlier pooled set network was 3$\times$ worse than GN.
- With neighbour/twin spots and pixel noise, plain GN degrades about $20\times$ and Huber-robust GN
  recovers a third of it. The network trained on realistic windows (deliberately perturbed synthetic data, Section 2.7) is 2--3$\times$ better than robust
  GN on distractor data, but not on pixel noise alone; trained on clean data only, it lies between
  plain and robust GN on distractor data and is confidently wrong. There is no real-data result.

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
$K=4$ window half-width (frames), $T$ number of unrolled IRLS iterations in `GNLayerNet` (the run labels "K=1", "K=3"
in `MIGRATION_HISTORY.md` and result file names mean $T=1$, $T=3$), $\Delta\omega_f=1^\circ$ frame width, $r_\perp$ distance of the voxel from the
rotation axis, $\omega^*_p(\delta)$ rotation angle at which peak $p$ satisfies the Bragg condition, $B$ batch size,
$\varphi$ a frame offset (Section 3.2), $\Gamma_p$ the $2\times3$ spot-motion Jacobian (px per degree here, px per radian in
the theory note) and $J_p$ the $3\times3$ per-peak Jacobian in $\sigma$ units (Section 3.2; the theory note's $J$ is the
information matrix).

## 1.2 What the network approximates

The thresholded measurement map is piecewise constant and has no inverse. Squared loss on simulated pairs
drives a network to the posterior mean $E[\delta\mid\text{data}]$ under the training prior (uniform in the ball
$\lvert\delta\rvert\le1^\circ$, `sample_prior_offsets`); a Gaussian likelihood loss adds the posterior covariance.
The Bayes risk floors every method (TN §2, §4, §5). Three results shape the design:

1. Rotating by $\beta$ about $\hat z$ shifts every $\omega^*_p$ by exactly $-\beta$ (TN §3.1.5; tested), so frames carry
   $\delta_z$ as a *global* sum over peaks.
2. A spot's detector position moves through the spot-motion Jacobian $\Gamma_p$ (ring sliding plus parallax, TN §3.6.3).
   Its $z$ column is pure parallax, $-\Pi Q_p(\hat z\times\mathbf x_v)$ (TN §3.6.3), proportional to $r_\perp$ and zero on the
   axis, so pixels constrain $\delta_z$ only away from the axis.
3. In the independent-quantisation model (frame variance $\Delta\omega_f^2/12$, pixel variance $1/12$ px$^2$) the covariance is
   $(\sum_p J_p^\top J_p)^{-1}$ ($J_p$ the $\sigma$-unit Jacobian of Section 3.2). The noise-free exact-Bayes floor is
   13--22$\times$ below it (TN §3.6.7), not a practical target.

## 1.3 Toy assumptions and geometry

| Assumption | Value |
|---------|---------------------------|
| Data | Thresholded lit/unlit pixels from the observer's forward model; rocking width $\alpha=0$ (each peak in one frame; $\alpha>0$ not implemented); no noise except the realism layer (2.7) |
| Peaks | $Q_{\max}=8$ Angstrom$^{-1}$ (`--max-q 8`); both detectors; only $\lvert\sin\eta\rvert\ge0.3$ at nominal (`--min-sin-eta 0.3`, near-axis peaks drift many frames per degree); identity fixed at nominal (new peaks ignored, vanished peaks give blank windows) |
| Windows | $32\times32$ px about the nominal centroid, $\pm4$ frames about the nominal frame; anything outside is blanked |
| Prior, test | train: uniform ball $\lvert\delta\rvert\le1^\circ$; test: $\lvert\delta\rvert\in\{0.1,0.25,0.5,1.0\}^\circ$, random directions |
| Sample | copper; triangular voxels, ManyGrains side 9.375 $\mu$m (24,570 voxels), ThreeVoxels 0.75--1.5 $\mu$m |

Geometry (`Examples/Example2.ThreeVoxels/ConfigFiles/StandardGeometry.2Det`, `Example2.Simulation.config`; ManyGrains is
identical except `MaxQ`). Near-field HEDM, distances in mm: beam along $+\hat x$ at 64.351 keV; rotation about $\hat z$ in 180 frames
of $1^\circ$ from $-90^\circ$ to $+90^\circ$ (`omega_180_2L.dat`), `EtaLimit 86`; two $2048\times2048$ detectors with $1.48\ \mu$m
pixels at lab $x$ = L1 3.36157 mm (index 0; beam centre J,K = 1022.77, 1925.81 px) and L2 5.38718 mm (index 1; 1020.30, 1921.31 px).
Measured: Gauss-Newton with both detectors improves the $\perp$ error by 1.0--1.8$\times$ and the $z$ error by 0--37%
over either detector alone (`diag_multi_summary.txt`, medians by $r_\perp$ bin), so the detectors are partly redundant for $\delta_z$.
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
 dataset .pt -> make_realistic_windows -> decode_windows -> x (B,M,2,W,W)
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
Rule 3 is the "recorded peak" rule: of 790 peaks hitting voxel 0's detector plane only 364 are recorded at nominal (the Stage 0 set: first detector, config `MaxQ`, no
$\eta$ cut; the rasteriser clips the rest); all numbers here use the corrected set. The generator then drops entries below
`--min-sin-eta` and, in multi-voxel mode only, rejects voxels with fewer than `--min-peaks` (40) entries. `detectors="first"` (first detector only) is the Stage 0 behaviour.

## 2.2 `BatchedObserver`

A float64, vectorised re-implementation of the simulator's per-peak maths: `observe(delta_deg (B,3))` returns `present (B,M)`, `frame (B,M)` ($-1$ out of range), `omega (B,M)` (rad) and `verts (B,M,3,2)` ((col,row) on the peak's own detector); also `observe_frames`, `sin_eta`, `vertex_keys`, `lit_pixel_set` (the rasteriser's clip, round and scanline fill). It agrees with the simulator on presence, frame and centroids to $10^{-3}$ px (tested); cost per candidate about 1 ms (unverified: no saved log).

## 2.3 Frame-coded windows

`WindowSpec.from_nominal(observer, W, K)` fixes per entry the integer origin
$(\text{col}_0,\text{row}_0)=(\lfloor c_{\text{nom}}\rfloor-W/2,\ \lfloor r_{\text{nom}}\rfloor-W/2)$ and the nominal frame $\text{frame}_0$ (all
ROI spots must be present at nominal). `render_windows(observer, spec, deltas_deg)` returns `windows uint8 (N,M,W,W)` and `status uint8
(N,M)`. Window pixel $(y,x)$ is detector pixel (row $\text{row}_0+y$, col $\text{col}_0+x$); a lit pixel of a spot recorded in frame $n$
holds the code $1+(n-\text{frame}_0)+K\in\{1,\dots,9\}$, unlit is 0. The lit set is the rasteriser's own, so with $\alpha=0$ the encoding
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

`torch.save` dicts, gitignored (`scripts/*.pt`), written as `toy_orientation_<tag>_{train,test}.pt`. Single-voxel files (no `--n-voxels`) hold `windows (N,M,W,W)`, `offsets_deg (N,3)`, meta (`window_size, n_peaks, voxel_index, R_nom, prior_radius_deg, max_q, detectors, min_sin_eta, renderer, frame_half_width, seed, example, r_perp_um, side_um`) and, with `--renderer observer`, `status (N,M)` and `context (M,16)`; test files add `magnitudes_deg`. (`--renderer simulator`, the default there, is the legacy renderer: it writes binarised ($w>0$) uint8 single-channel windows
from `render_local_windows`, without frame coding.) The current multi-voxel format (`--n-voxels V`; shapes checked on `toy_orientation_arch_dis_*.pt`):

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
the Stage 3 aux file is reused for the distractor data. `toy_orientation_arch_dis_val.pt` (200 samples = 50 per voxel for training voxels 6, 11, 19, 24, `magnitudes_deg=0`, numpy seed 0, drawn
from `toy_orientation_arch_dis_train.pt`) is used only to choose the Huber threshold; `scripts/make_dis_val.py` recreates it (verified sample-for-sample against the file used).

Voxel selection. `select_voxels` picks, for each of 30 radii evenly spaced in $[0,500]\ \mu$m, up to 60 random candidates (`--voxel-seed 0`); `accept_voxels` takes the first usable one per radius and never reuses a grain (tested with fakes). Sorted by $r_\perp$, every 5th voxel from index 2 is **held out** (indices 2, 7, 12, 17, 22, 27; $r_\perp=32,126,215,301,372,457\ \mu$m): no training samples, 240 of the 1200 test cases.

## 2.7 Realistic data: distractors and detector noise

**Terminology: realistic data.** *Realistic* windows (test sets, training) are the exact synthetic windows with two kinds of deliberately simulated complexity added, to approach what real detector data contain: (1) *overlap* (`neighbours` variant): genuine diffraction spots of neighbouring voxels and a Σ3 twin that land in the target's windows — real signal from other grains, making the data more complex rather than noisier; (2) *detector noise* (`noise` variant): missing spots, threshold jitter at spot edges, hot pixels and spurious blobs. `all` combines both. *Clean* (ideal) windows have neither. Orientation labels are always exact; only the windows a method sees change. *Realism-trained* = trained on windows made realistic on the fly (`--realistic-train all`). The model is still a toy: same-family train/test nuisances, at most 3 sources, no intensity effects (see the caveats). Earlier versions of this work (and committed logs, result filenames with `corr`, and the old aliases `--corrupt-train`, `--corrupt`, `CorruptionConfig`, `corrupt_windows`, `corrupt_dataset`) called this *corruption*; it never meant faulty data.

Both change only the windows, never $\delta$. **Distractor layer** (`--neighbors 2 --twin`; defaults `--neighbor-radius-um 30`, `--neighbor-p 0.5`,
`--neighbor-sigma-deg 0.3`): `build_distractor_sources` takes up to `--neighbors` mic voxels nearest the target in the sample plane (not closer than
1 $\mu$m), each with its own orientation, vertices and ROI set, plus a $\Sigma3$ twin ($60^\circ$ about $[111]$, `R @ T`) of the nearest. Per sample each
source follows the *target's* $\delta$ plus an independent Gaussian offset ($0.3^\circ/\sqrt3$ per component, $0.3^\circ$ rms) and is present with probability
0.5. `render_distractor_windows` draws each present source spot that overlaps a target window (same detector, within $\pm K$ frames) with the same frame
coding; `combine_windows` overlays it *under* the target, whose pixels always win (tested).

Neighbour-model caveat (checked from the `.mic`). In the 30-voxel data all 60 neighbours (2 per voxel) lie at 9.375 $\mu$m (one voxel pitch), and in 45 of 60
the neighbour's orientation equals the target's (same grain, misorientation $<10^{-3}$ degrees). For those the distractor is the target's own grain seen from an adjacent
voxel, following the target's $\delta$ plus $0.3^\circ$ rms; the other 15 cross a grain boundary. The twin is a $\Sigma3$ ($60^\circ$ about $[111]$) twin of the
nearest neighbour; since that neighbour is in the target's own grain for 23 of 30 voxels (computed from the `.mic`), the twin is usually $\Sigma3$-related to the target, and its
twin-invariant reflections land on the target's own spots. Results hold for this mixture only. Distractor sources render only spots in their *own* filtered ROI set (`min_sin_eta`, `max_q`), so the distractor model understates contamination.

**Detector noise** (`RealismConfig`; per sample and entry independently, in this order): `neighbours=True` (overlay `dis_windows`); `p_flip=0.05` (each lit
pixel dropped, each 4-neighbour of a lit pixel lit with its code: edge jitter); `p_hot=0.05` (one isolated hot pixel, random frame code); `p_blob=0.1` (one
$2$--$4\times2$--$4$ px blob, random frame code); `p_miss=0.1` (whole window zeroed). Variants (`RealismConfig.named`): `clean`/`none`, `neighbours` (layer only),
`noise` (pixel terms only), `all`. `make_realistic_dataset(windows, distractors, name, K, seed=12345, chunk=100)` gives a deterministic copy so all methods see identical test inputs.
Quirk (from the code, not tested): the realism layer also hits zero-padded entries, which can then count as "present"; their $J=0$ so they add nothing to $A$ or $b$ and can only perturb the pooled covariance features.
The GN baselines slice windows to the true peak count and do not see them. `make_realistic_windows` / `make_realistic_dataset` take an optional `valid` mask that zeroes padded entries (tested); it was not used in the reported runs.

# 3. Model

## 3.1 `GNLayerNet`: data flow

Here $y$ is the measurement in $\sigma$ units, $J$ the per-peak Jacobian, $s=0.05^\circ$ the unit of $\delta'=\delta/s$ and $f$ the
64-dimensional peak encoding (all defined in Sections 3.2--3.3).

```
 x (B,M,2,W,W) -> measurement_features(nom_off) -> y (B,M,3) [sigma units]
 context (B,M,16): Gamma, d omega*/d delta        -> J (B,M,3,3)
 (x, ctx, meas) -> PeakEncoder -> f (B,M,64)
   [pairing: f <- f + has * pair2(relu(pair1([f, f_partner, has])))  (pair2 zero-init, linear; see 3.4)]
 repeat T times (shared head weights), starting at delta' = 0:
   rho = y - J delta'
   head([f, asinh(rho), delta' (, asinh(rho_partner), has)])
        -> w (B,M,3), dy (B,M,3)
   A  = sum_i J_i^T diag(w_i) J_i + ridge I
   delta' = A^-1 sum_i J_i^T diag(w_i) (y_i + dy_i)
 pooled = mean over present peaks of f -> D = diag(exp(g(pooled)))
 delta_hat = s delta' (deg);  Sigma = s^2 D A^-1 D;  L = chol3(Sigma)
```

## 3.2 Physics inputs

For present peaks (absent: $y=0$, $J=0$), with $(\bar c,\bar r)$ the lit centroid (px), $\varphi$ the mean frame offset (frames) and $(c_0,r_0,\varphi_0)$ the exact
nominal values (`nom_off`), the measurement in units of its quantisation sigma ($1/\sqrt{12}$ px, $1/\sqrt{12}$ frame) is
$$
y_i=\sqrt{12}\,\big(\bar c_i-c_{0,i},\ \bar r_i-r_{0,i},\ \varphi_i-\varphi_{0,i}\big)\in\mathbb R^3 .
$$
Its Jacobian with respect to $\delta'=\delta/s$, with unit scale $s=$ `delta_scale` $=0.05^\circ$, is
$$
J_i=\sqrt{12}\,s\begin{bmatrix}\Gamma_i\\ \gamma_i^\top\end{bmatrix}\in\mathbb R^{3\times3},
$$
where $\Gamma_i\in\mathbb R^{2\times3}$ (px/deg) is context columns 0:6 times 20 and $\gamma_i\in\mathbb R^3$ is columns 6:9 in frames per degree (times
$(\pi/180)/\Delta\omega_f$, `frame_width_rad` from the aux file). $\Gamma_i$ here is per degree (the theory note's $\Gamma_p$ is per radian) and
$J_i$ is the per-peak Jacobian in $\sigma$ units, not the theory note's information matrix. This is where the $r_\perp$ parallax enters.

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
calibration: it rescales the standard deviations but keeps the correlation structure that $A$ takes from the geometry. Hypothesised need (not isolated by an experiment): quantisation errors
are not independent across peaks (the clean GN run in `dis_log_gn_clean.txt` gives rms$_z$ 0.019--0.024$^\circ$ against a predicted $\sigma_z$ of 0.0154$^\circ$, $1.3$--$1.5\times$), and $w,\Delta y$ change what $A$ means.

## 3.4 Iterations, pairing, shapes

`n_iter` $=T$ (`--gn-iters`) unrolls IRLS-style rounds with *shared* head weights. Each round the head also sees the residual at the current estimate, $\mathrm{asinh}(y-J\delta')$,
and $\delta'$, so it can down-weight outliers; each round re-solves from the nominal linearisation (the Jacobian is not recomputed). On clean data iterations add nothing. Realism-aware
training is what matters on realistic data; whether the iterations help there was not tested (no $T=1$ realism-trained run).

`pairing=True` (`--pairing`) adds two things, both using the other-detector entry of the same diffracted ray (`aux["pair_index"]`; its encoding $f_{\text{partner}}$, zero if absent):

1. **Encoder-level mixing:** $f\leftarrow f+\text{has}\cdot W_2\,\mathrm{ReLU}(W_1[f,f_{\text{partner}},\text{has}])$, `pair1` $129\to64$, `pair2` $64\to64$ (12,480 parameters), masked by the has-partner flag so unpaired entries are unchanged.
   `pair2` is zero-initialised and has *no* ReLU after it, so the network starts as the unpaired one and the update trains (`pair2` gets gradient immediately, `pair1` once `pair2` has left zero, as `head2`/`cov2` are also zero-initialised).
   This was a bug until 2026-10-01: the update used to be $\mathrm{ReLU}(W_2\ldots)$; with $W_2=0$ the pre-activation is exactly 0, $\mathrm{ReLU}'(0)=0$ in PyTorch, and `pair1`/`pair2` received exactly zero gradient forever, so all earlier paired runs ("inert mixing" in Section 7) used only item 2. An even earlier non-zero-initialised version (`gnpairv1`) trained this path and was unstable; the zero-initialised linear version was stable at lr 1e-4 in all four fixed runs, no gate or normalisation was needed.
2. **Partner residual:** the head also sees $\mathrm{asinh}$ of the partner entry's residual and the has-flag, a pair-consistency check. Both detectors' rows still enter $A$ through their own $J$ rows.

| Tensor | Shape | Notes |
|---------|------------|------------|
| windows $x$; context | `(B,M,2,32,32)` float32; `(M,16)` or `(B,M,16)` | trainer passes `context[voxel_id]` |
| aux `nom_off`; `pair_index` | `(M or B, M, 3)`; `(M or B, M)` long | `pair_index` only with pairing |
| `y`, `w`, `dy`; `J` | `(B,M,3)`; `(B,M,3,3)` | rows col/row/frame, columns $\delta_{x,y,z}$ |
| encoding $f$; head | `(B,M,64)`; in $64+6(+4)$, hidden 64, out 6 | `feat_dim=64` |
| $A$, $A^{-1}$; output | `(B,3,3)`; `delta_hat (B,3)` deg, `L (B,3,3)` | $\hat\Sigma=LL^\top$ |

`PeakEncoder`: windows $\mathrm{conv}(4\to8,3,\text{stride }2)\to\mathrm{conv}(8\to16,3,\text{stride }2)\to\mathrm{fc}(1024\to64)$ (4 channels: lit, frame, two coordinate channels in
$[-1,1]$); context MLP $16\to64\to64$; measurement MLP $5\to64\to64$; concatenation (192) $\to64\to64$; ReLU. `GNLayerNet` has 102,657 parameters (115,393 with pairing, of which 12,480 are in the pair MLP; computed),
independent of $M$ and $T$, and is invariant to peak order and zero padding (tested, with and without pairing).

## 3.5 Other architectures still in the code

- `ToyOrientationNet` (`--head quat`): flatten, FC trunk, unit quaternion; the v0 baseline, kept to re-evaluate under the Stage 0 protocol.
- `ToyOffsetNet` (`--arch fc`): same trunk with $\delta$ and a Cholesky factor, fixed peak count; Stage 0--1 reference, single-voxel only.
- `PeakSetNet` (`--arch set`): shared conv encoder, context, optional measurement features (`--no-meas`; `--subpixel` for the exact nominal measurement), masked pooling `--pool {meanmax,meansum,all}`; the pre-GN comparison point.
- `FrameProbeNet` (`--arch probe`): per-peak MLP of the mean frame offset and $\partial\omega^*/\partial\delta$ only, mean pooled, `--loss mse`; a diagnostic that $\delta_z$ survives mean pooling.

# 4. Training

`scripts/train_toy_orientation_nn.py`. `--arch gn` needs `--aux` and `--head offset`; `--arch set` needs `--aux --subpixel` to use the exact nominal measurement.

**Losses** (`--loss`; `icenine/orientation_nn.py`): `nll` (default; `gaussian_nll_loss` $=\tfrac12\lVert L^{-1}(\delta-\hat\delta)\rVert^2+\sum\log L_{ii}$, `--beta-nll b>0` weights each sample by
$\mathrm{stopgrad}(\det\Sigma^{b/3})$); `mse` (`mse_deg_loss` $=\tfrac12\sum_{\text{axes}}(\delta-\hat\delta)^2/s_m^2$, $s_m=$ `--mse-scale` $=0.1^\circ$); `decoupled` (`decoupled_nll_loss`: `mse` on $\hat\delta$ plus the NLL
at $\mathrm{stopgrad}(\hat\delta)$, so the $\hat\delta$ output receives exactly the MSE gradient and $\hat\Sigma$ exactly the NLL gradient; tested). In `GNLayerNet` this is not a full separation: the NLL
term still trains the shared encoder and the weights $w$ (through $A^{-1}$), and $w$ also sets $\hat\delta$. `mse-then-cov`/`mse-then-nll` (`mse` for the first half of the epochs, then `decoupled` or `nll`;
only second-half epochs can be the best-validation checkpoint). *Why decoupled:* the NLL gradient on the mean is $\Sigma^{-1}(\hat\delta-\delta)$, so the $z$ gradient is smaller than the $\perp$ gradient by
$(\sigma_\perp/\sigma_z)^2$ and a network that is unsure of $z$ stops learning it; this stalled the pooled set network on voxel 0, `decoupled` fixed it and per-sample $\beta$-NLL did not
(details in MIGRATION_HISTORY). All current `GNLayerNet` runs use `decoupled`.

**Optimiser.** Adam, `--lr` (3e-4 for the first GN runs, 1e-4 for Step 3/4), `--batch-size 64`, `--epochs 60`, `--clip 1.0` (gradient norm), `--cosine` (decay to 0) (values used; trainer defaults are lr 1e-3, batch 32, 30 epochs, clip 0, no cosine, ema 0, checkpoint `best`, device `cpu`, arch `fc`, loss `nll`). `--ema d` keeps a per-step exponential average of
the weights and validation is evaluated on it; `--checkpoint`: `best` (best-validation epoch), `ema` (EMA weights at the last epoch), `last`. Current runs use `--ema 0.998 --checkpoint ema`:
best-epoch selection on four validation voxels picked early epochs (12--20 of 60) that were worse everywhere (MIGRATION_HISTORY, Stage 3 retrain).

**Validation split.** Default: a random `--val-frac` (10%) of training *samples*, in-distribution, blind to degradation on new voxels. `--val-voxels N` (multi-voxel): `split_by_voxel` sorts the training voxels by $r_\perp$,
cuts $N$ equal strata and draws one voxel per stratum with `--seed`; their samples leave training and are used only for selection and monitoring. With $N=4$: seed 0 gives voxels 6, 11, 19, 24 ($r_\perp=99,181,326,410\ \mu$m; 10,000
train and 2,000 validation samples), seed 1 gives 3, 11, 20, 29 (computed with `split_by_voxel`; MIGRATION_HISTORY quotes seed 0 only). The six held-out test voxels are never used for selection; the validation voxels are among the 24
"in-dist" test voxels, so that group contains unseen-voxel cases that differ by seed.

**Devices.** `--device {cpu,mps}`: network training and inference in float32 (headline runs on the Apple GPU: 60 epochs took 400 s ($T=1$) to 540 s ($T=3$) on clean data, 950--1,050 s with the realism layer; one clean $T=3$ run, lr 1e-4, logged 4,441 s). Data stay on the CPU; the observer, `ExactBayes`, the GN baselines (about 10 ms per case, Huber GN 24--34 ms) and all error statistics are float64 on the CPU.

**Realism-aware training.** `--realistic-train {none,neighbours,noise,all}`: each training batch is made realistic on the fly (`make_realistic_windows`, generator `--seed + 7`); validation windows are made realistic once (`--seed + 99`) so the
validation loss is comparable across epochs; `neighbours` and `all` need `dis_windows`. Test variants: `--eval-variants clean,neighbours,noise,all` (deterministic, seed 12345).

Headline runs: `--arch gn --head offset --loss decoupled --device mps` with the optimiser settings above and `--val-voxels 4`, seeds 0 and 1 (Section 9; the clean runs' flags are as in MIGRATION_HISTORY, the result files do not store them).

# 5. Evaluation protocol

**Test sets.** Fixed $\lvert\delta\rvert\in\{0.1,0.25,0.5,1.0\}^\circ$, random directions (`sample_fixed_magnitude_offsets`), 10 per voxel and magnitude: 1200 cases over 30 voxels, 960 on the 24 training-set voxels ("in-dist": 20 trained on plus 4 used only for validation, with no training samples; 160 of the 960 cases; which 4 depends on the seed)
and 240 on the 6 held-out voxels. Single-voxel sets have 30 cases per magnitude.

**Metrics.** `error_summary(delta_hat, delta_true)` returns, in degrees, `rms_z`, `rms_perp` ($\sqrt{\text{mean}((e_x^2+e_y^2)/2)}$), `rms_total`, `median_angle` (median rotation angle between $\exp([\hat\delta]_\times)$ and $\exp([\delta]_\times)$),
`success_0p5` (misorientation $<0.5^\circ$, the `bench_hp_sweep.py` criterion) and `success_0p1` ($<0.1^\circ$, used from Stage 1); errors are rotation-vector differences (accurate to $1^\circ$). *Calibration:* Mahalanobis$^2$
$=(\delta-\hat\delta)^\top\hat\Sigma^{-1}(\delta-\hat\delta)$ averaged over cases; 3 for a calibrated Gaussian, $\gg3$ overconfident (1500 = confidently wrong), $<3$ underconfident (a necessary, not sufficient, check: coverage and tails are not examined). *Parallax test:* per voxel the median angle over its 40 test
cases; the Pearson correlation with $r_\perp$ over the 30 voxels (`corr_err_rperp`) and the median per-voxel error in-dist and held-out. GN on clean data has corr $-0.63$ (error falls as parallax grows); a method not using parallax has corr $\ge0$.

**Baselines.**

- *predict-nominal*: $\hat\delta=0$ (median angle $=\lvert\delta\rvert$).
- *Exact Bayes* (`ExactBayes`, `exact_bayes_baseline.py`, seed 7): posterior of $\delta$ given exact noise-free data (presence, frame, rasteriser pixel sets) by adaptive importance sampling on the consistent set; its mean is the best squared-error estimate and
  $\sqrt{\mathrm{tr}\,\mathrm{cov}}$ the floor. Limits: single-voxel datasets only (not run on the 30-voxel data); ignores newly appearing peaks and spot overlaps; uses exact pixels (a generous reference).
- *Plain GN* (`CentroidGaussNewton`, `gauss_newton_baseline.py`): Levenberg-Marquardt on three residuals per spot (frame-centre $\omega$ with $\sigma=\Delta\omega_f/\sqrt{12}$; centroid col, row with $1/\sqrt{12}$ px), the observer as exact model,
  central-difference Jacobians (`fd_step_deg=0.002`), start at nominal, `max_iter=30`, `tol_deg=1e-7`; returns `delta`, `cov` $=(J^\top\mathrm{diag}(w)J)^{-1}$ with $J$ the Jacobian of the normalised residuals and $w$ the IRLS weights (all 1 for plain GN), `n_used`, `status`, $\chi^2$/dof. Limits: it is given the peak
  identities and the basin and sees only per-spot centroid and frame; about 7 ms per case.
- *Huber GN* (`huber_c`, `--huber c`): the same with Huber loss on normalised residuals via IRLS (weight $\min(1,c/\lvert r\rvert)$); one start, one robust loss. $c$ was chosen on a validation split, not the test set: 200 prior-ball samples from the four validation voxels
  (6, 11, 19, 24), windows made realistic with `all`; $c=1$ minimised the realistic-validation median (`step4_huber_val.txt`). $c=1$ is used for all variants, including clean data (about 20% worse held-out median there: 0.0144 vs 0.0119). With `--max-iter 200` all fits converge and the median is unchanged (MIGRATION_HISTORY).
- *MC and Riemannian Adam* (`optimizer_baselines.py`): the `bench_hp_sweep.py` protocol on the Python-simulated ThreeVoxels images at voxel 0's truth; each optimiser starts at $\exp([\delta]_\times)R_{\text{nom}}$ and searches for the truth. MC: 3500 steps, 2
  restarts, step fraction 0.5, search box $1.5\lvert\delta\rvert$ (it is told $\lvert\delta\rvert$). Adam (geoopt): scale 2, $\omega$ window 1, lr $10^{-4}$, 100 steps. Limits: voxel 0 only (ManyGrains has no images); a harder spot-association problem; one run, 30 cases per bin; never run on realistic data.

**Seeds.** Dataset: generator `--seed 42` (one stream), `--voxel-seed 0`, distractor draws `seed+1000`. Trainer `--seed s` sets initialisation, batch order and the validation-voxel draw (headline nets $s\in\{0,1\}$); training/validation realism draws `s+7`/`s+99`; test
realism draws fixed 12345 for all methods; exact Bayes `--seed 7`; MC `seed + case index`. Network results are 2-seed means (per-seed values where tabulated); one dataset draw, no confidence intervals; GN is deterministic. Seed spread of the held-out median (two seeds)
reaches 0.006--0.008$^\circ$ for Step 2/3 nets on clean data (mean of the four bin medians; $T=3$, lr 3e-4: 0.0206 vs 0.0130), 0.003$^\circ$ for Step 4 nets on clean data, 0.0055$^\circ$ for paired nets on `all`
and 0.05$^\circ$ for the clean-trained net on `all` (0.3016 vs 0.2498).

# 6. Test matrix

86 tests pass in 23 s (`uv run pytest tests/test_orientation_nn.py tests/test_orientation_eval.py tests/test_orientation_baselines.py`, 2026-10-01; the glob `tests/test_orientation_*.py` also matches `test_orientation_search.py`, 107 tests in all). **nn** = `test_orientation_nn.py`, **ev** = `test_orientation_eval.py`, **bl** = `test_orientation_baselines.py`; the `test_` prefix of method names is dropped; "phys" = needs
ThreeVoxels voxel 0 (skipped if absent; every physics test uses that voxel).

| Test (file::class::name) | Property guaranteed |
|---------------|---------------|
| nn::{TestDefineROISet (3), TestRenderLocalWindows (4), TestSampleLocalPerturbations (2), TestLossAndMetric (4), TestToyOrientationNet (1)} | legacy Stage-0 pipeline: ROI set deterministic, state restored (phys); legacy renderer and sampler; quaternion loss and unit-norm output |
| nn::split_by_voxel_disjoint_deterministic_and_spans_r_perp | validation voxels disjoint, seed-deterministic, one per $r_\perp$ stratum; bad `n_val` raises |
| ev::TestRotations (6), TestErrorSummary (3) | Rodrigues equals scipy; $R=\exp([\delta]_\times)R_{\text{nom}}$ acts in the sample frame; ball prior and fixed-magnitude sampler; error-summary values, $z$/$\perp$ separation |
| ev::TestPixelGrid (3) | recorded-peak rule; `lit_pixel_set` equals the rasteriser on 400 random triangles; every ROI peak on the grid at nominal (phys) |
| ev::TestOffsetHead (7) | Cholesky lower-triangular, positive diagonal; NLL equals `MultivariateNormal`; $\beta$-NLL finite; decoupled mean gradient = MSE gradient, covariance gradient = NLL gradient; NLL mean gradient $\propto\Sigma^{-1}$ |
| ev::TestBatchedObserver::{matches_simulator[first,all], all_detectors_mode, invalid_detectors_mode, missing_any_detector_plane_drops_peak_everywhere, max_q_override_limits_reflections, nominal_offset_reproduces_roi_set, sample_state_untouched} | presence, frame, centroids ($10^{-3}$ px) equal the serial simulator at 4 offsets (phys); detector-mode, `max_q` and drop-everywhere rules; sample not mutated |
| ev::TestBatchedObserver::rotation_about_stage_axis_shifts_every_omega_exactly | $\omega^*\to\omega^*-\beta$ for a $z$ rotation, to $10^{-9}$ |
| ev::TestExactBayes (2), TestRenderWindows (3), TestEndToEndVsForwardSimulation | truth inside the posterior, $\sigma<0.005^\circ$, pixels tighten it; window decoding and codes $1+(n-\text{frame}_0)+K$ with correct status; windows at nominal equal `ForwardSimulation` images pixel for pixel |
| bl::TestPeakSetNet (7), TestMeasurementFeatures (5) | shapes, positive-definite covariance, parameter count independent of $M$, permutation invariance, padding ignored, finite with no peaks; feature values |
| bl::TestMeasurements (2), TestMeasurementFeaturesMatchExtraction, TestGaussNewton (3), TestGaussNewtonStatus | extraction consistent with the observer; GN far better than nominal at $0.6^\circ$, small at nominal, $\sigma_z>3\sigma_\perp$, status consistent (phys) |
| bl::TestPeakContext, TestFrameProbeNet, TestVoxelSelection (3) | context shape (M,16), column 8 $=-1$, one-hot (phys); probe uses only frame and $\partial\omega^*/\partial\delta$; radius selection, a grain never reused (synthetic mics) |
| bl::TestGNLayer::{normal_equations_match_weighted_lstsq, unit_weights_reproduce_linear_gauss_newton_step} | `gn_normal_equations`/`inv3` equal weighted least squares; $w=1,\Delta y=0$ gives `solve_linear` and covariance diagonal within 3% (phys) |
| bl::TestGNLayer::{permutation_and_padding_invariance[False,True], gradients_flow_and_are_finite[gn,gn_paired,set], pair_mlp_trains_from_the_zero_init, paired_net_equals_unpaired_net_at_initialisation} | order and padding invariance with and without pairing (the test re-initialises the zero-initialised pair parameters, so it exercises the pair path with non-zero weights); every parameter of `GNLayerNet` (paired, unpaired, $T=3$) and `PeakSetNet` gets a non-zero finite gradient after a few small steps off the zero-initialised layers (fails on the pre-fix code); `pair2` leaves zero under SGD and then `pair1` trains; at initialisation the paired net equals the unpaired one (extra head inputs zeroed), `pair1` is irrelevant while `pair2`=0, and the pair path is live once `pair2`$\neq0$ |
| bl::TestPairingAndNominalOffsets (4), TestDistractors (3), TestRobustGaussNewton | `pair_index` symmetric and correct incl. unpaired/3-way; features minus `nominal_offsets` = measurement minus exact nominal; distractors never alter target pixels, realism variants deterministic and bounded; Huber with huge $c$ equals plain GN and resists 20 gross outliers (phys) |

**Not tested:**

- The trainer end to end (loss selection, EMA, checkpoints, `--val-voxels`, realism-aware training, MPS vs CPU) and the generator's multi-voxel path, `build_distractor_sources`, `make_dataset_aux.py`, the GN/exact-Bayes/optimiser scripts and `summarize_*` (only `split_by_voxel` and the loss functions are unit tested).
- Physics on ManyGrains voxels or at large $r_\perp$: all physics tests use ThreeVoxels voxel 0 ($r_\perp=12\ \mu$m), so the parallax columns of the context are untested far from the axis.
- Accuracy and calibration of any trained network (no regression test on result numbers).
- `GNLayerNet` beyond equivalence at initialisation: outlier rejection for $T>1$, off-diagonals of `chol3`, the learned $D$; whether the trained pair MLP learns anything useful (only that it receives gradient; weights are not saved).
- Nuisance statistics (rates, blob sizes, code ranges), realism of padded entries, multi-source/twin distractors with real neighbours; `ExactBayes` calibration and three or more detectors; other geometries, $\alpha>0$, other noise models.

# 7. Current results

Setup unless stated: 30 ManyGrains voxels, 24 train (4 used for validation) and 6 held out; train offsets in the $1^\circ$ ball; test $\lvert\delta\rvert=0.1/0.25/0.5/1.0^\circ$, 10 offsets per voxel and magnitude; `GNLayerNet` trained 60 epochs, decoupled loss, EMA
weights, two seeds averaged; medians in degrees; one dataset draw, no confidence intervals.

## 7.1 Clean data

Median angle at the four $\lvert\delta\rvert$ bins, from `step2_summary.txt`, `step3_summary.txt`, `step1_subpixel_summary.txt`, `summarize_results.py` on `multi_res_*_s*.json`, and `step4_summary.txt` (plain GN). "In-dist" means the 24 training-set voxels,
which includes the 4 validation voxels (no training samples; seed dependent, Section 4). Mah$^2$ ranges are min--max over the four bins; corr and voxel medians are 2-seed means (seed 0 / seed 1 in brackets for corr). The last column is the net's median angle divided by GN's
(mean of the four bin medians; a ratio of the four-bin means), in-dist and held-out, seed 0 / seed 1, from the saved jsons.

| Method (2 seeds) | In-dist | Held-out | corr(err, $r_\perp$) [s0 / s1] | Voxel median err in / held | Mah$^2$ in / held | Net/GN in; held |
|------------|---------|---------|------|------|------|------|
| `PeakSetNet`, `--pool all --loss decoupled --subpixel` | .033/.035/.033/.043 | .029/.034/.039/.053 | +0.29 | .0354 / .0362 | 3.1--5.1 / 4.0--9.3 | -- |
| `GNLayerNet` $T=1$, lr 3e-4 | .013/.015/.012/.013 | .017/.016/.015/.013 | $-0.72$ [$-0.64$/$-0.79$] | .0136 / .0188 | 3.1--3.2 / 3.5--4.0 | .93/1.12; 1.15/1.43 |
| `GNLayerNet` $T=3$, lr 3e-4 | .013/.015/.012/.015 | .018/.016/.017/.016 | $-0.64$ [$-0.43$/$-0.85$] | .0131 / .0169 | 3.1--3.4 / 3.8--4.2 | 1.17/.91; **1.74**/1.10 |
| `GNLayerNet` $T=3$, lr 1e-4 | .013/.016/.013/.013 | .016/.016/.015/.015 | $-0.61$ [$-0.75$/$-0.47$] | .0138 / .0163 | 2.8--3.2 / 3.5--3.7 | .92/1.19; 1.24/1.38 |
| `GNLayerNet` $T=3$ paired, inert mixing$^*$, lr 1e-4 | .011/.014/.012/.014 | .013/.015/.016/.016 | $-0.60$ [$-0.70$/$-0.50$] | .0143 / .0145 | 2.8--3.2 / 3.9--4.5 | .93/1.04; 1.13/1.40 |
| `GNLayerNet` $T=3$ paired, fixed mixing$^\dagger$, lr 1e-4 | .013/.015/.013/.014 | .015/.016/.015/.013 | $-0.63$ [$-0.49$/$-0.77$] | .0145 / .0166 | 2.7--3.4 / 3.1--3.7 | -- |
| plain GN | .012/.015/.012/.014 | .014/.012/.012/.010 | $-0.63$ | .0127 / .0143 | -- | -- |

$^*$ "Inert mixing": the runs made before 2026-10-01, in which only the partner-residual input was active; the pair-encoding MLP got zero gradient (Section 3.4).
$^\dagger$ After the fix (Section 3.4): encoder-level mixing trains; two seeds, `multi_res_gnpairfix_k3_lr1e-4_s{0,1}`, `pairfix_clean_summary.txt`. Held-out per seed (s0 / s1), medians at the four bins: .014/.015/.015/.013 / .016/.017/.016/.014; z RMS .019--.023 and $\perp$ RMS .003--.005 (inert .017--.028 / .003--.007, unpaired .021--.033 / .003--.004; GN .016--.020 / .003--.004).

The Step 2 commit subject (e6f4072) says the net "reaches Gauss-Newton parity"; that overstates the held-out result. Reading the table: the net is close to GN on the 24 training-set voxels (0.9--1.2$\times$) and 1.1--1.7$\times$ GN on unseen ones, with GN's parallax trend (negative correlation in every run, strength seed dependent); on clean data the net's held-out error is
1.05--1.35$\times$ its own in-dist error. The set network's $3\times$ gap is removed in every run. Pairing neither helps nor hurts on clean data beyond seed spread, with inert mixing (partner residual only) and with live mixing alike (fixed pairing: same medians, $\perp$ and $z$ RMS inside the spread; held-out Mah$^2$ 3.3--3.5 vs 3.9--4.5 inert and 3.5--3.7 unpaired).

Single-voxel checks (30 cases per bin, one run each; `single_v0_gn.json`, `single_v77_gn.json`, `benchmarks/toy_orientation_stage3/far_res_bayes.json`):

| Voxel | Method | Median at the four bins |
|-------|--------|-------------------------|
| 0 ($r_\perp=12\ \mu$m) | `GNLayerNet` $T=1$ / GN / exact Bayes | .0142/.0176/.0184/.0249 / .0401/.0242/.0348/.0362 / .0069/.0040/.0028/.0040 |
| 77 ($r_\perp=399\ \mu$m) | `GNLayerNet` $T=1$ / GN / exact Bayes | .0098/.0078/.0086/.0095 / .0118/.0106/.0094/.0079 / .0006--.0009 |

MC and Riemannian Adam on voxel 0 (from `benchmarks/toy_orientation_stage2/pred_*.npz`) reach medians .048/.145/.360/.672 and .056/.131/.185/.523 and leave $\delta_z$ largely uncorrected on these near-axis voxels (their 96%/92% success at $1^\circ$ in the ManyGrains sweep was a different, far-from-axis sample; the $r_\perp$ dependence was not measured).

## 7.2 Realistic data (Step 4)

"Realistic" means the deliberately simulated complexity of Section 2.7 (overlapping neighbour/twin spots and detector noise on exact synthetic windows), not faulty data; labels stay exact.

Dataset `toy_orientation_arch_dis_*` (seed 42; 12,000/1,200 samples; clean windows identical to the Stage 3 data). Median angle pooled over the four bins (it does not depend on $\lvert\delta\rvert$ in any row; per-bin values are in `step4_summary.txt`). "Clean-trained":
the same architecture trained on the clean windows of this file; "realism-trained": `--realistic-train all`; "paired, inert": realism-trained with `--pairing` before the fix (partner-residual input only, Section 3.4); "paired, fixed": the same after the fix (live encoder mixing; `dis_res_gnpairfix_k3_corr`, `step4_pairfix_summary.txt`); all nets are $T=3$ `GNLayerNet`; Huber $c=1$; nets are 2-seed means, held-out per-seed values in brackets; $<0.1^\circ$ is the held-out fraction; Mah$^2$ is in-dist / held-out. "In-dist" includes the 4 validation voxels (Section 4).

| Test set | Method | Median in-dist | Median held-out [seeds] | $<0.1^\circ$ | Mah$^2$ |
|------|---------|---------|---------------|------|---------|
| clean | plain GN | .0131 | .0119 | 1.00 | |
| | Huber GN | .0125 | .0144 | 1.00 | |
| | net, clean-trained | .0138 | .0152 [.0143, .0161] | 1.00 | 3.0 / 3.6 |
| | net, realism-trained | .0185 | .0193 [.0176, .0210] | 1.00 | 2.1 / 2.2 |
| | paired, inert mixing | .0185 | .0186 [.0185, .0186] | 1.00 | 2.3 / 2.3 |
| | paired, fixed | .0188 | .0199 [.0214, .0183] | 1.00 | 2.3 / 2.5 |
| neighbours | plain GN | .2253 | .3428 | 0.20 | |
| | Huber GN | .1375 | .2223 | 0.33 | |
| | net, clean-trained | .1934 | .2706 [.2916, .2496] | 0.23 | 1556 / 1951 |
| | net, realism-trained | .0576 | .0719 [.0723, .0715] | 0.66 | 3.0 / 2.9 |
| | paired, inert mixing | .0594 | .0710 [.0702, .0718] | 0.65 | 3.0 / 2.7 |
| | paired, fixed | .0585 | .0706 [.0739, .0673] | 0.68 | 3.0 / 2.9 |
| noise | plain GN | .0545 | .0573 | 0.85 | |
| | Huber GN | **.0177** | **.0215** | 1.00 | |
| | net, clean-trained | .0704 | .0741 [.0746, .0736] | 0.71 | 223 / 272 |
| | net, realism-trained | .0241 | .0242 [.0227, .0258] | 0.99 | 3.0 / 3.0 |
| | paired, inert mixing | .0233 | .0231 [.0228, .0233] | 1.00 | 3.0 / 3.0 |
| | paired, fixed | .0250 | .0257 [.0262, .0253] | 0.99 | 3.0 / 3.3 |
| all | plain GN | .2228 | .3311 | 0.20 | |
| | Huber GN | .1405 | .2239 | 0.32 | |
| | net, clean-trained | .2033 | .2757 [.3016, .2498] | 0.16 | 1462 / 1796 |
| | net, realism-trained | **.0615** | **.0821** [.0827, .0815] | 0.64 | 2.9 / 2.6 |
| | paired, inert mixing | **.0616** | **.0762** [.0734, .0789] | 0.62 | 2.8 / 2.5 |
| | paired, fixed | .0629 | .0775 [.0780, .0770] | 0.64 | 2.8 / 2.8 |

1. Distractor spots break GN (about $20\times$ worse than clean, 20--34% of cases within $0.1^\circ$); Huber GN removes about a third of the error and loses the parallax trend. Pixel noise alone costs plain GN $4\times$; Huber GN recovers nearly all of it.
2. The realism-trained net is 2--3$\times$ better than Huber GN on `neighbours` and `all`, with mean Mahalanobis$^2$ near 3 (2.5--3.0 held-out too; coverage not checked); on `noise` alone Huber GN is as good or better (0.0215 vs 0.0242; every net seed is worse). Errors are heavy-tailed: held-out $\perp$ RMS 0.05
   vs 0.004 clean, 34--38% of held-out cases above $0.1^\circ$.
3. Realism-aware training is what matters; whether the iterations help on realistic data was not tested (no $T=1$ realism-trained run). The clean-trained net lies between plain and Huber GN on distractors (9--21% better than plain GN: in-dist .1934/.2033 vs .2253/.2228, held-out .2706/.2757 vs .3428/.3311), is worse on `noise` (0.074 vs 0.057) and wildly overconfident; realism-aware training costs clean-data accuracy (held-out 0.0176--0.0210 vs 0.0143--0.0161 vs 0.0119 GN) and makes the covariance underconfident on clean data (Mah$^2$ 2.1--2.3).
4. Held-out voxels are harder under distractors for every method (held-out/in-dist: net 1.3$\times$, plain GN 1.5$\times$, Huber GN 1.6$\times$); the $r_\perp$ trend is lost for all (corr $-0.09$ to $-0.18$ for plain GN and the realism-trained nets, $-0.17$ to $+0.04$ for the clean-trained net, $+0.13$/$+0.14$ for Huber GN).
5. Pairing has no clear benefit, with or without live encoder mixing. Inert mixing (partner residual only): differences $\le0.002^\circ$ except held-out `all` (.0734/.0789 vs unpaired .0827/.0815). Fixed mixing (2 seeds, 6 voxels): `neighbours` .0706 vs .0719 unpaired, `noise` .0257 vs .0242 (slightly worse), `all` .0775 [.0780, .0770] vs .0821 [.0827, .0815], i.e. the small held-out `all` gain is reproduced at the same size, not enlarged, so it is suggestive at best; error-vs-$r_\perp$ correlations are unchanged (clean $-0.58$/$-0.74$, `all` $-0.08$/$-0.10$). Untested hypothesis: a neighbour spot lands consistently on both detectors, so the partner check cannot reject it.

Caveats. Two seeds, six held-out voxels; seed spreads of the held-out median are given in Section 5. Train and test nuisances are the same family with the same neighbour sets per voxel (new random draws, not new kinds of nuisance), so the
advantage over robust GN may shrink out of family; the neighbour-model caveat of Section 2.7 applies. Robust GN is one fixed construction. Learning rate and $T$ were not re-tuned for realistic data. What the weight head learned has not been analysed.

# 8. Known limitations, open questions, roadmap

- **Simulation only.** No real data; $\alpha=0$, known peak identities, nominal orientation within $1^\circ$. Real near-field HEDM resolution is about $0.1^\circ$ (theory note); the toy's clean errors ($0.01^\circ$) are noise-free quantisation only and the Bayes floor (0.001$^\circ$ at voxel 77, 0.004--0.008$^\circ$ at voxel 0; `GNLayerNet` median errors are 2--7$\times$ the exact-Bayes median at voxel 0 and 9--16$\times$ at voxel 77, per bin) is not reachable.
- **Fixed peaks and windows.** Peaks absent at nominal are invisible; about 2--3% of spots at the $1^\circ$ prior are not fully inside their windows (voxel 77); near-axis peaks are dropped; windows are not sized from the Jacobian. The layer linearises at nominal (negligible up to $1^\circ$, no re-linearisation for larger priors).
- **Generalisation and calibration.** Six held-out voxels from one sample, geometry and structure; on clean data held-out error is 1.1--1.7$\times$ GN's and 1.05--1.35$\times$ the same net's in-dist error. $\hat\Sigma$ is checked on this distribution only, and only through the mean Mahalanobis$^2$ (3.5--4.5 held-out clean, up to 1.5$\times$ overconfident; realism-trained nets underconfident on clean data): the older set network was overconfident on unseen voxels (superseded; Mah$^2$ 17--34) and the clean-trained GN net is
  catastrophically overconfident on realistic data.
- **Realism caveats.** Same-family train/test, at most 3 sources, no intensity effects; neighbours are mostly same-grain and the twin is usually $\Sigma3$-related to the target (2.7). The realism layer also hits padded entries of the multi-voxel arrays in all reported runs (`--mask-padding` was not used): they have J = 0, so they do not move the estimate delta, but hot pixels/blobs can make them look "present", so they enter the pooled covariance features and the count n; the Gauss-Newton baseline slices `[:n_pk]`, so the comparison is slightly asymmetric against the net. Use `--mask-padding` for future runs.
- **No integration with `FindOptimal`.** The coarse-search hand-off sometimes selects a wrong ~54$^\circ$ solution (neighbour voxel or $\Sigma3$ twin); a refiner should be evaluated per candidate.

Open questions: what the weight head learns on distractor data and why the partner check does not reject neighbour spots; out-of-family nuisance sweeps; more seeds and confidence intervals; recalibration of $\hat\Sigma$ on held-out voxels; the $\sim1.2\times$ held-out gap to GN on clean data; the heavy $\perp$ tail on distractor data.

Roadmap (all unimplemented): per-peak second moments, an $\alpha>0$ renderer, a Sinkhorn render-and-compare loss and posterior scoring against exact Bayes (KL, coverage) -- `docs/todo_intensity_and_distribution_losses.md`; iNeRF-/BARF-style refinement and self-supervised training on real data (parked) -- `docs/research_ideas_nerf.md`; why the coarse search and `FindOptimal` return neighbour or twin orientations -- `docs/todo_findoptimal_wrong_candidates.md`; multiple detectors with many grains and voxels far from the rotation axis, where detector pairing is expected to matter most -- `docs/todo_multidetector_many_grains.md`; rocking width $\alpha$ for real data -- `MIGRATION_HISTORY.md`, Stage 1 decisions D1, D4.

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

# 1b. Huber-c validation set (200 samples; numpy seed 0, voxels 6 11 19 24)
uv run python scripts/make_dis_val.py

# 2. GN and Huber-GN on each test variant (clean = --realistic none)
for V in clean neighbours noise all; do
  C=$V; [ "$V" = clean ] && C=none
  uv run python scripts/gauss_newton_baseline.py \
      --test ${S}_arch_dis_test.pt --realistic $C \
      --out $D/dis_pred_gn_$V.npz
  uv run python scripts/gauss_newton_baseline.py \
      --test ${S}_arch_dis_test.pt --realistic $C --huber 1 \
      --out $D/dis_pred_huber1_$V.npz
done

# 3. nets: this is the realism-trained net. Drop --realistic-train for
#    the clean-trained control; add --pairing for the paired net (run the fixed code: names with a
#    `_pairfix` suffix, e.g. dis_res_gnpairfix_k3_corr, multi_res_gnpairfix_k3_lr1e-4); seeds 0
#    and 1. The clean multi-voxel nets use stage3_multi_{train,test}.pt,
#    --lr 3e-4 or 1e-4, --gn-iters 1 or 3 (flags as in MIGRATION_HISTORY)
#    and --extra gn=benchmarks/toy_orientation_stage3/\
#    multi_pred_gauss_newton.npz (without it the results have no gn key,
#    GN row or GN corr; the file equals dis_pred_gn_clean.npz to 2e-7)
uv run python scripts/train_toy_orientation_nn.py --head offset \
    --arch gn --gn-iters 3 --loss decoupled --device mps --lr 1e-4 \
    --clip 1 --cosine --batch-size 64 --epochs 60 --ema 0.998 \
    --checkpoint ema --val-voxels 4 --realistic-train all \
    --eval-variants clean,neighbours,noise,all --seed 0 \
    --train ${S}_arch_dis_train.pt --test ${S}_arch_dis_test.pt \
    --aux ${S}_stage3_multi_aux.pt \
    --results-json $D/dis_res_gn_k3_corr_s0.json \
    --save-predictions $D/dis_res_gn_k3_corr_s0.npz

# 4. tables (the run label "corr-trained" = realism-trained; kept so the committed summaries reproduce)
uv run python scripts/summarize_results.py "K3=$D/multi_res_gn_k3_s?.json"
uv run python scripts/summarize_arch_step4.py \
    --test ${S}_arch_dis_test.pt --huber 1 --seeds 0 1 \
    --runs "corr-trained=dis_res_gn_k3_corr" \
    "paired-inert=dis_res_gnpair_k3_corr" \
    "paired-fixed=dis_res_gnpairfix_k3_corr" \
    "clean-trained=dis_res_gn_k3_cleantrain" \
    --out $D/step4_summary.json
```

The clean-variant predictions of `dis_res_gn_k3_cleantrain_s{0,1}.npz` are byte-identical to `multi_res_gn_k3_lr1e-4_s{0,1}.npz` (same deterministic run); the realistic-variant files exist only under the `cleantrain` name, so both are kept.

Saving weights and the perturbation sweep (2026-10-04). `--save-model` was added to the trainer afterwards; the four unpaired T=3 nets were retrained one at a time with the commands above plus `--save-model scripts/toy_orientation_sweep_model_{clean,realistic}_s{0,1}.pt` (clean-trained: the `multi_res_gn_k3_lr1e-4` flags on `stage3_multi`, no `--realistic-train`; realism-trained: as in item 3) and reproduce the committed predictions bit-for-bit. Two trainings sharing the MPS GPU are not bit-reproducible, so run them sequentially. Then, from `icenine_py/`:

```bash
M=scripts/toy_orientation_sweep_model
uv run python scripts/perturbation_sweep.py run --workers 10 --out-dir benchmarks/toy_orientation_sweep \
    --models clean_s0=${M}_clean_s0.pt clean_s1=${M}_clean_s1.pt \
             realistic_s0=${M}_realistic_s0.pt realistic_s1=${M}_realistic_s1.pt
uv run python scripts/perturbation_sweep.py summarize --out-dir benchmarks/toy_orientation_sweep
uv run pytest tests/test_perturbation_sweep.py
```

(about 22 minutes of wall time with 10 CPU workers; per-voxel results are cached in `scripts/perturbation_sweep_cache/` so an interrupted run resumes.) Protocol and results: `MIGRATION_HISTORY.md`, "Perturbation sweep".

Existing-optimizer baselines on the same sweep cases (MC, Riemannian Adam, plain and Huber Gauss-Newton; per-case images built from pixel sets, no full-sample render; about 2 h with 10 CPU workers, resumable via `scripts/optimizer_sweep_cache/`):

```bash
uv run python scripts/optimizer_sweep.py run --workers 10 --out-dir benchmarks/toy_orientation_sweep
uv run python scripts/optimizer_sweep.py summarize --out-dir benchmarks/toy_orientation_sweep   # txt, json, perturbation_sweep_vs_optimizers.png
uv run pytest tests/test_optimizer_sweep.py
```

Protocol and results: `MIGRATION_HISTORY.md`, "Comparison with existing optimizers".

Multi-level reconstruction (FindOptimal) on the same cases: experiment A runs `AdaptiveVoxelReconstructor.reconstruct_voxel` from scratch once per voxel and variant (no starting guess), experiment B runs only its final stage (`refine_from_candidates`: FindOptimal + VarianceMinimizing + final evaluation) from each case's perturbed start; errors are reduced by cubic symmetry. About 5 min (A) + 32 min (B) (the 1 h was the pilot projection) with 10 CPU workers, resumable via `scripts/findoptimal_sweep_cache/`:

```bash
uv run python scripts/findoptimal_sweep.py pilot --workers 10      # timing pilot (2 voxels)
uv run python scripts/findoptimal_sweep.py run-a --workers 10
uv run python scripts/findoptimal_sweep.py run-b --workers 10 --n-dirs 20
uv run python scripts/findoptimal_sweep.py summarize   # findoptimal_sweep_summary.{txt,json}, perturbation_sweep_vs_findoptimal.png
uv run pytest tests/test_findoptimal_sweep.py tests/test_findoptimal_refactor.py
```

Protocol and results: `MIGRATION_HISTORY.md`, "Comparison with multi-level reconstruction (FindOptimal)".

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

Gaps in the record. (1) The single-voxel `GNLayerNet` runs are not in the commands above: `single_v77_gn` used `--arch gn` on `stage3_far_{train,test}.pt` with `--aux stage3_far_aux.pt --extra gn=benchmarks/toy_orientation_stage3/far_pred_gauss_newton.npz`, lr 3e-4; `single_v0_gn` is the analogous run on the voxel-0 files (its exact command is unverified). (2) The Table 7.1 `PeakSetNet --subpixel` row is `--arch set --pool all --loss decoupled --subpixel` (flags as reported, not stored in the result files). (3) The original run scripts (`run_s2.sh`, `run_s3*.sh`, `run_s4.sh`) lived in a session scratchpad and are not in the repo. (4) Unverified: for the Huber-validation GN predictions under `benchmarks/toy_orientation_arch/val/`, and for the Stage 1 single-voxel set (tag and sample count follow MIGRATION_HISTORY; not re-run).

Tests: `uv run pytest tests/test_orientation_nn.py tests/test_orientation_eval.py tests/test_orientation_baselines.py`.
