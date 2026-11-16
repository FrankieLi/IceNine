# AR6 ideas

One row per idea. Status: proposed, running, done, dropped. Result in one line with its number;
details in the findings page named in the last column. New ideas get the next free id (I32).
`build_page.py` reads this table: keep the columns and the `|` layout.

| id | goal | idea | status | result | findings |
|---|---|---|---|---|---|
| I1 | G1a | Bias study: where each objective's minimum sits against the ground truth (6 variants + LS) | done | GO1: Sinkhorn variants flat within 0.031-0.041 deg of ground truth (objective <1% lower); strict cost 0.006-0.008 deg | findings.md |
| I2 | G2a | Pixel-box dilation of the prediction (lit pixels are those a spot touches) | done | GO1: no effect (0.034 / 0.041 deg, same as without) | findings.md |
| I3 | G2b | Soft "pixel touched" model: sigmoid of signed distance pixel-square to spot union | proposed | | |
| I4 | G1a | Sub-pixel noise floor: objective curvature vs strict cost curvature, per spot | proposed | | |
| I5 | G1a | Stage A CMA: budget sweep 50 / 100 / 200 / 400 evaluations, sigma0 0.02 / 0.05 deg | running | GO2: 50/100/200/400 evals: >0.1 deg 20/3/3/3 of 853, median 0.0148/0.0146/0.0111/0.0062 | findings.md |
| I6 | G1a | Stage B Adam / L-BFGS on per-voxel offsets with smoothness weight sweep | proposed | | |
| I7 | G1a | Gauss-Newton / LM on per-spot centroid (and second-moment) residuals, autograd Jacobian | proposed | | |
| I8 | G1a | Coarse-to-fine blur (eps annealing) and warm-started Sinkhorn potentials | proposed | | |
| I9 | G1a | Per-voxel CMA on the strict cost (AR3 C) and on the grain objective, same voxels | running | GO2: per-voxel CMA (sigma 0.02) 61 evals: 34 of 853 >0.1; AR3 C (sigma 0.05): 18 | findings.md |
| I10 | G1a | Hybrid: grain field (stage A) then per-voxel strict-cost polish only where needed | running | GO2: field 100 + per-voxel 10/30/61: median 0.0151/0.0130/0.0114, >0.1 stays 3 | findings.md |
| I11 | G1b | Sample generator: uncorrelated random voxel rotations about the grain mean, 0.5-5 deg | done | `rand` built and rendered (22.4M lit px); 69%/97% of neighbour pairs over 2 deg at 2/5 deg | PLAN.md GO3 |
| I12 | G1b | Sample generator: spatially correlated random fields (GND-like), correlation 2-10 voxels | done | `corr3` (length 3 voxel sides) built and rendered (19.2M lit px) | PLAN.md GO3 |
| I13 | G1b | Robustness of the best G1a methods across linear, quadratic, correlated, random fields | proposed | | |
| I14 | G1c | Strained reciprocal vectors in SpotModel (opt-in): g' = R (I + eps)^-T g, Bragg from abs g' | proposed | | |
| I15 | G1c | Strain in the stock forward model render path (opt-in, from the .mic deformation tensor) | proposed | | |
| I16 | G1c | Strain sensitivity: objective curvature in the 6 strain components; recovery at 1e-4 to 1e-3 (resolution scales with pixel size / distance) | proposed | | |
| I17 | G2a | Intensity-weighted Sinkhorn (observed mass = clean intensity) | done | GO1: sharpens on grain 381 (0.0074 deg) but not 302 (0.026 deg) | findings.md |
| I18 | G2a | Poisson NLL on pixel counts with a per-reflection scale and background | proposed | | |
| I19 | G2a | Per-spot moment matching (centroid + covariance), as objective and as LM residual | proposed | | |
| I20 | G2a | MMD with a Gaussian kernel per spot | proposed | | |
| I21 | G2a | Least squares with per-reflection scale on intensity (baseline arm) | done | GO1: 0.012 / 0.042 deg; basin-finder arm in GO2 | findings.md |
| I22 | G2a | Use structure factor x Lorentz to predict relative spot intensity across reflections | proposed | | |
| I23 | G2b | Rocking curve alpha > 0: spots spread over frames; sub-frame omega from frame ratios | proposed | | |
| I24 | G2b | Partial-volume coverage per pixel and a detector point-spread function in the render | proposed | | |
| I25 | G2c | GNLayerNet (orientation_nn_design.md) as a per-voxel start inside a grain | proposed | | |
| I26 | G2c | Learned residual / amortised Gauss-Newton step for the grain field (todo_hybrid_nn) | proposed | | |
| I27 | G2c | Neural field for orientation (and strain) over the grain, render-and-compare (nerf note) | proposed | | |
| I28 | G1a | Grain-field CMA on the mean strict cost of the group's voxels (sharp objective, cheap) | running | GO2: field CMA on strict, 100 evals: 3 of 853 voxels >0.1 deg vs 18 for per-voxel CMA (C) | findings.md |
| I29 | G1a | Two-phase: Sinkhorn (or LS) as basin finder, strict cost as finisher | running | GO2: Sinkhorn stage rejected by the gate on early grains | findings.md |
| I30 | G1c | Borrow HEXRD far-field strain fitting: per-grain lattice-strain fit from spot centroids; sensitivity vs pixel/distance ratio (e.g. a far-field-like detector distance in our config) | proposed | Reference saved (owner, 2026-11-16): hexrd.github.io, local ~/Research/github/hexrd | |
| I31 | G1c | Survey strain-fitting code in older IceNine C++ branches (GrainFitting, IceNine_v3.0.x, fix_pmc_opt, archive snapshot, master) and ~/Research/IceNine.others/IceNine.Strain | done | Strain only in the C++ image simulator (g = O M g0, M ~ I - eps, Bragg from abs g); no fitter exists; branches identical; reuse convention + .mic cols 13-18 + config keys | findings_cpp_strain.md |
