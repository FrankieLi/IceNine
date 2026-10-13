# What flipping the default of CppOverlapRule would change (measured list; the default is NOT flipped)

Method: temporarily set the defaults of `VoxelCostFunction(cpp_overlap_rule, center_screen)`,
`SearchParameters.cpp_overlap_rule` and `ConfigFile.cpp_overlap_rule` to True, run
`uv run pytest tests/` (305 s), restore. Result: 13 failed, 843 passed, 34 skipped, 1 deselected.
Two of the 13 are this task's own tests that assert the default (False). The other 11:

| Test file | Failing tests | What pins the old value |
|---|---|---|
| `test_findoptimal_refactor.py` | `test_reconstruct_voxel_unchanged_by_refactor`, `test_recorder_hook_leaves_reconstruct_voxel_bit_identical`, `test_rank_key_with_post_mc_costs_reproduces_golden`, `test_extra_candidates_returning_nothing_reproduces_golden` | `GOLDEN` toy reconstruct_voxel results |
| `test_cma_optimizer.py` | `test_default_local_optimization_unchanged`, `test_default_refine_from_candidates_unchanged`, `test_local_optimization_cma_uses_cma` | `GOLDEN_REFINE` and default-path goldens |
| `test_coarse_proxy.py` | `test_q_max_none_leaves_features_unchanged`, `test_features_batch_equals_per_candidate[5.0]`, `[None]` | proxy features computed against the cost function |
| `test_nn_hybrid.py` | `test_findoptimal_alone_reproduces_stored_experiment_b` | stored experiment B |

Not run in this checkout: the tests that need `cpp_outputs/` (the directory is absent here, 34 skipped);
they compare diffraction geometry (omegas, symmetry, ranges), not cost values, and are not expected to
change, but this was not verified.

Stored studies that use `VoxelCostFunction` on the 2-detector data and so would be re-run or re-labelled
(all are expected to shift; none was re-run here): the Phase A to D finisher/MC study results
(`benchmarks/cost_sensitivity`, `finisher_bench`, `mc_mechanism`, `cma_finisher` E0, `findoptimal_robustness`),
the seed-cost diagnosis, parity and pilot tables in `benchmarks/phase_d_seed_diag` (candidate counts, evaluations,
BFS seed hit ratios against the 0.8 threshold, reseed pilot), and the NN/proxy work that harvests candidates or
labels by cost (`coarse_proxy`, harvested candidates, hybrid). The rule changes hit ratio and quality as well as
the screen, so every threshold on them (BFS seed test, `MaxDeepeningHitRatio`) moves.
