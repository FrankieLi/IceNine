# AR6 ideas

One row per idea. Status: proposed, running, done, dropped. Result in one line with its number;
details in the findings page named in the last column. New ideas get the next free id (I28).
`build_page.py` reads this table: keep the columns and the `|` layout.

| id | goal | idea | status | result | findings |
|---|---|---|---|---|---|
| I1 | G1a | Bias study: where each objective's minimum sits against the truth (6 variants + LS) | running | Smoke test: stage A moved voxels 0.016 to 0.050 deg away while lowering Sinkhorn | findings.md |
| I2 | G2a | Pixel-box dilation of the prediction (lit pixels are those a spot touches) | running | in GO1 | findings.md |
| I3 | G2b | Soft "pixel touched" model: sigmoid of signed distance pixel-square to spot union | proposed | | |
| I4 | G1a | Sub-pixel noise floor: objective curvature vs strict cost curvature, per spot | proposed | | |
| I5 | G1a | Stage A CMA: budget sweep 50 / 100 / 200 / 400 evaluations, sigma0 0.02 / 0.05 deg | proposed | | |
| I6 | G1a | Stage B Adam / L-BFGS on per-voxel offsets with smoothness weight sweep | proposed | | |
| I7 | G1a | Gauss-Newton / LM on per-spot centroid (and second-moment) residuals, autograd Jacobian | proposed | | |
| I8 | G1a | Coarse-to-fine blur (eps annealing) and warm-started Sinkhorn potentials | proposed | | |
| I9 | G1a | Per-voxel CMA on the strict cost (AR3 C) and on the grain objective, same voxels | proposed | | |
| I10 | G1a | Hybrid: grain field (stage A) then per-voxel strict-cost polish only where needed | proposed | | |
| I11 | G1b | Sample generator: uncorrelated random voxel rotations about the grain mean, 0.5-5 deg | proposed | | |
| I12 | G1b | Sample generator: spatially correlated random fields (GND-like), correlation 2-10 voxels | proposed | | |
| I13 | G1b | Robustness of the best G1a methods across linear, quadratic, correlated, random fields | proposed | | |
| I14 | G1c | Strained reciprocal vectors in SpotModel (opt-in): g' = R (I + eps)^-T g, Bragg from abs g' | proposed | | |
| I15 | G1c | Strain in the stock forward model render path (opt-in, from the .mic deformation tensor) | proposed | | |
| I16 | G1c | Strain sensitivity: objective curvature in the 6 strain components; recovery at 1e-4 to 1e-3 | proposed | | |
| I17 | G2a | Intensity-weighted Sinkhorn (observed mass = clean intensity) | running | in GO1 | findings.md |
| I18 | G2a | Poisson NLL on pixel counts with a per-reflection scale and background | proposed | | |
| I19 | G2a | Per-spot moment matching (centroid + covariance), as objective and as LM residual | proposed | | |
| I20 | G2a | MMD with a Gaussian kernel per spot | proposed | | |
| I21 | G2a | Least squares with per-reflection scale on intensity (baseline arm) | running | in GO1 | findings.md |
| I22 | G2a | Use structure factor x Lorentz to predict relative spot intensity across reflections | proposed | | |
| I23 | G2b | Rocking curve alpha > 0: spots spread over frames; sub-frame omega from frame ratios | proposed | | |
| I24 | G2b | Partial-volume coverage per pixel and a detector point-spread function in the render | proposed | | |
| I25 | G2c | GNLayerNet (orientation_nn_design.md) as a per-voxel start inside a grain | proposed | | |
| I26 | G2c | Learned residual / amortised Gauss-Newton step for the grain field (todo_hybrid_nn) | proposed | | |
| I27 | G2c | Neural field for orientation (and strain) over the grain, render-and-compare (nerf note) | proposed | | |
