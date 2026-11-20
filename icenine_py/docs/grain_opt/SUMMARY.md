# AR6 (`grainopt`) summary

Program closed 2026-11-20 (owner, decision card "Wrap up with a summary"). Branch
`feature/grain-opt`, not merged. Details per step: `findings.md`; ideas and verdicts: `IDEAS.md`;
judging rules written before each run: `PLAN.md`. All numbers are on emulated (rendered) images of
the `grad`, `rand`, `corr3` and `strain*` samples; errors are against the ground truth. Timings are
contended single runs unless stated.

## At a glance

- **Orientation (G1a):** fit each grain as one smooth orientation field on the strict cost
  (`Astrict/100`): 3 of 853 C voxels above 0.1 deg against 18 for AR3's per-voxel CMA, in half its
  wall time (3 ms per voxel); a polished boundary pass leaves 1 of 853.
- **Deformed grains (G1b):** ladder + gated seed path + neighbour pass + wider seed path + grid +
  a second pass at gate 0.11 + per-voxel polish: above 0.1 deg rand 4.5 / 0.4 / 3.3 / 3.2% and corr3
  1 / 1 / 0.4 / 7% at 0.5 / 1 / 2 / 5 deg voxel-to-voxel misorientation (ladder alone: up to 72%).
- **Strain (G1c):** fast Gauss-Newton centroid fit (`GP4@P`, Tukey weights, ratio filter 2.5,
  region from the cleaned start): strain error medians 1.3 / 1.2 / 1.4e-4 at 1e-4 / 5e-4 / 1e-3, no
  stalls, 29 s per grain; the CMA hybrid stays the reference. An intragranular strain gradient is
  only partly observable (internal-field relative error 0.72 at best).
- **Intensity (G2a):** only the Poisson deviance comes near the strict cost (optimum 0.010 / 0.006
  deg from the ground truth); MMD, moments, least squares and Sinkhorn are biased or noisy.
- **Simulation (G2b):** model-side PSF and rocking (library, opt-in) did not help; a strict-or-Poisson
  acceptance gate rescues rocked images (20 -> 5 of 853 above 0.1 deg; clean 3) and is now the
  scripts' default.
- **Learned models (G2c):** GNLayerNet trained on single-voxel windows does not transfer to observed
  grain windows.

## Answers per goal

**G1a, cost against quality for a fixed objective.** The grain field on the strict cost is the
frontier on G20 (GO2): 100 evaluations per voxel give 3 of 853 C voxels above 0.1 deg (AR3 C: 18;
intervals apart), p95 0.035 deg (C 0.055); 400 evaluations lower the median from 0.0146 to 0.0062
deg. Smooth objectives first (least squares, Sinkhorn) cost 20-250 x more time for no gain. Paint errors
on grain boundaries are recovered only when neighbour candidates are polished before comparison
(GO2c: 1 of 853, 4.7 ms per voxel). Near the ground truth the smooth objectives' offsets lie within
their evaluation noise (GO1b), so they are basin finders, not finishers.

**G1b, robustness to the kind of gradient.** Far starts need a global step: GO4d cached the paired
observed rays per sample so the gated seed path runs on `rand` / `corr3`; GO4e-GO4i added a
neighbour pass, a wider seed path, a grid search and a lower second gate. Result above. The
remaining errors on `corr3` at 5 deg (7%) are voxels whose right orientation no step proposes.

**G1c, elastic strain.** Near field resolves about 0.2-0.5 px per 1e-3 strain (GO5), so strain is
recoverable but needs sub-pixel centroids. Coverage costs alone fail; the hybrid (coverage for
orientation, centroids for strain) works (GO5g/h); Gauss-Newton with frozen Tukey weights, a CMA
fallback for stalls and a polish above 6e-4 matches it in half the time (GO5l); a 2.5 ratio filter
and windows from the outlier-cleaned start remove the stalls (GO5v). On noisy images the fast fit
nearly matches the CMA hybrid (misses it only at 5e-4: 2.6e-4 against 2.4e-4, GO5p). A 1 px detector
PSF breaks the fit unless the ratio filter is widened (GO7e: 2.1 / 3.5e-4). Intragranular gradients:
linear, quadratic and GP-prior fields recover the internal field only partly (relative error 0.72
at 1e-3 at best, GO5q-GO5u).

**G2a, intensity.** Poisson deviance with per-reflection scales: optimum 0.0101 / 0.0060 deg (GO9),
missing the 0.01 rule by a hair; structure-factor scales, per-spot MMD and moment matching are worse
(GO9b-d). As an acceptance check it earns its place (GO7g).

**G2b, simulation.** Rocking (alpha 1.0) and a 1 px PSF on the images cost one grain and the strain
fit respectively (GO7, GO7b). The opt-in model PSF and rocking curve (library change D-20261119-I1,
three reviews, defaults bit-identical) did not fix either (GO7f). The rocking damage is threshold
dilution: weak reflections fall under the binarisation threshold (unlit predicted peaks 0.3% ->
12%). Lower thresholds restore the data but not the gate; the strict-or-Poisson gate does (GO7g).
A soft "pixel touched" surrogate is unbiased at the ground truth but too wide to finish (GO1c/d).

**G2c, learned models.** GNLayerNet retrained to its documented accuracy (GO8a) does not transfer to
grain windows with intragranular spread (GO8b); retraining on grain-like windows (I66) is open.

## Recommended pipelines (scripts in `scripts/grain_opt/`)

- **Orientation, mildly deformed grains:** `run_pipeline.py --pipeline Astrict/100` (grain field on
  the strict cost; strict-or-Poisson gate by default, `--gate strict` for the old rule); add
  `P2/60/0.12` for boundary paint errors.
- **Strongly deformed / far starts:** `L/366/3,G/15/0.15,N/30/0.3,W/20/0.15,X/1/0.15,W/20/0.11,
  X/1/0.11,V/61/0.02` with the per-sample ray cache (`ray_cache.py`).
- **Strain:** `fit_strain.py --arm GP4@P --cost hyb --lam 8 --ratio_hi 2.5 --clean_start`; the CMA
  hybrid (`--arm S@P --cost hyb --lam 8`) as the reference when time allows.

## What did not work, and why

- Sinkhorn / least squares as finishers: flat or noisy bottoms at the 0.01 deg scale (GO1, GO1b).
- Per-spot MMD, moment matching: whole-window clouds include merged spots and neighbour pixels
  (GO9c/d). Moment residuals work only with robust weights inside the strain fit.
- Model PSF and strict-cost rocking window: the damage is not in the model's footprint but in
  binarisation (GO7f). Lower thresholds fix the data, but the strict gate on grain 197 is a knife
  edge (GO7g).
- GNLayerNet inside grains: trained on single-voxel windows (GO8b).
- Splitting merged spots, PSF-aware coverage: did not undo the 1 px PSF damage (GO7c/d).

## Open questions that need real data

- Rocking and noise thresholds: the emulation has no noise or background, so lower thresholds and
  the Poisson gate look better here than they will on real counts.
- The Poisson gate on real counts below the binarisation threshold (noise, background, gain).
- The strain-gradient signal limit: how much of an intragranular field real spot shapes carry.
- PSF width and rocking width of the actual instrument (the GO7f options are ready for both).

## Open owner decisions

- GO7f defaults: renderer options as attributes and flags (built) or config keys (plan); PSF lit
  threshold 1e-3 of the frame maximum; |sin eta| floor 0.01; how the high500 frozen-file test treats
  approved library changes.
- Moving the strict-or-Poisson gate into the library (needs approval).
- I35: G1b from a real painted start on `rand` / `corr3` (the AR3 paint likely exceeds 2 h).
- I66: retrain GNLayerNet on grain-like windows.
- Review and merge of `feature/grain-opt` (nothing was merged to develop or master).
