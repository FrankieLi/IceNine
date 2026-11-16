# AR6 cost functions tried

One row per cost; add a row for every new cost when its step closes (the tracking page reads this
table, section "Cost functions tried"). Columns: what the cost compares; what it was used to fit;
steps that tested it; orientation result (the optimum's distance from the ground truth, deg; GO1
grains 381 / 302 unless noted); strain result (median strain error at 1e-4 / 5e-4 / 1e-3, strain
only at the ground-truth orientation, "@P" = from the painted start where run); bias and
sharpness; verdict (keep, finisher, basin finder, dropped, running). Details: `findings.md`.

| cost | compares | used for | steps | orientation result (deg) | strain result | bias / sharpness | verdict |
|---|---|---|---|---|---|---|---|
| strict cost (strict binary overlap) | per voxel, predicted voxel spots against lit pixels (`VoxelCostFunction`), binary | orientation (per voxel, grain field) | GO1, GO2, GO4, GO4b, GO4c, GO5d | 0.006 / 0.008 (only qualifier in GO1); grain field 100 eval: median 0.015 on G20 | no strain model; as a field leg on `strain` (GO5d) the field takes up strain as rotation (0.015 / 0.021 / 0.039 deg), strain after it @P 4.7 / 5.4 / 6.3e-4 | sharpest orientation cost; noisy per voxel (few pixels), pooled by the grain field | keep: orientation finisher and referee |
| coverage | per spot, share of predicted coverage mass on lit pixels | strain (with orientation field in strain fits) | GO5, GO5b, GO5c, GO5d (reference in GO5e, GO5f) | staged joint strain fit @P: 0.012 voxel median | 1.7 / 1.6 / 2.7e-4 (mean of 2 runs); @P staged 7.1 / 3.4 / 3.0e-4; @P long joint 6.5 / 5.7 / 4.5e-4 | flat band 1-3e-4 (saturates once all mass is lit); trades 0.01 deg against 2-4e-4 strain | keep: strain reference, orientation field in strain fits |
| Sinkhorn unbalanced (reach 3 px; blur 1 or 0.5 px; with or without pixel-box footprint blur) | transport distance between predicted coverage and lit pixels | orientation (basin finder) | GO1, GO2 | 0.031-0.041 for every variant | GO5 sensitivity: a 10% rise needs about 1e-3 per component (flat) | flat bottom below 0.04 deg; blur and footprint dilation change nothing | dropped (basin finder adds nothing from the painted start, about 700 x wall time) |
| Sinkhorn balanced (blur 1) | as above, balanced transport | orientation | GO1 | 0.095 / 0.013 | not run | worse than unbalanced (contamination) | dropped |
| intensity-weighted Sinkhorn | unbalanced Sinkhorn with lit pixels weighted by clean intensity | orientation (basin finder) | GO1, GO2, GO5 | 0.0074 / 0.026 | GO5 sensitivity: 10% rise needs 0.9-2.2e-3 per component | sharper than binary Sinkhorn on one grain only; 1,391 ms per voxel in GO2 | dropped (too slow, no gain) |
| least squares with per-spot scale | predicted intensity against observed, per-spot scale | orientation (basin finder) | GO1, GO2, GO4 | 0.012 / 0.042 | not run | between Sinkhorn and strict; helps far starts (GO4) | basin finder |
| signed distance to lit pixels (`sd`) | mass-weighted signed distance of predicted coverage to the lit set (sub-pixel) | strain | GO5e | not run | 1.8e-3 / 1.2e-3 / 8.9e-4 | sharp (scatter 6e-5) but biased: rewards shrinking the footprint into the spot (one-sided) | dropped |
| clean-intensity cost (`int`) | 1 - mass-weighted clean intensity under the predicted coverage | strain | GO5e | not run | 2.5e-3 / 1.9e-3 / 1.5e-3 | sharp (scatter 1e-4) but biased like `sd` | dropped |
| centroid match (`cen`) | per spot, predicted mass centroid against observed intensity centroid (truncated at 2 px) | strain | GO5f | @P staged: 0.06-0.09 voxel median (field gradient drifts) | 1.3 / 2.2 / 1.7e-4 (scatter 5e-6); @P staged 1.4 / 2.2 / 3.7e-4 | most consistent strain cost; small bias 1-2e-4; blind to the orientation gradient | keep for strain (with coverage for orientation: GO5g) |
| two-sided distance (`ch`) | predicted mass off the lit set plus lit pixels' distance to the nearest predicted point | strain | GO5f | not run | 1.7e-3 / 9.0e-4 / 7.4e-4 | biased like `sd` | dropped |
| hybrid (`hyb`) | `cen` + lambda x coverage | strain and orientation together | GO5g (running) | pending | pending | pending | running |
