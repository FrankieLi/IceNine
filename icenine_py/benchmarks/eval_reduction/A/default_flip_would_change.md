# Tests that fail if the Task A switches became the default

Method (as `benchmarks/phase_d_seed_diag/level0/default_flip_would_change.md`): the defaults were set to True
temporarily, `uv run pytest tests/` was run (no extra `-q`), and the defaults were restored. With the switches off
(the committed state) the suite is 889 passed, 37 skipped, 1 deselected (303 s), and no golden was regenerated or
changed.

## `CppEdgeRules` default on (cost function, config, `SearchParameters`, reconstructor): 10 failed, 879 passed

| Test | What it pins |
|---|---|
| `test_bfs_refit.py::test_default_bfs_matches_recorded_old_result` | recorded BFS result of the default search |
| `test_coarse_proxy.py::test_q_max_none_leaves_features_unchanged` | coarse-proxy features from a default search |
| `test_coarse_proxy.py::test_features_batch_equals_per_candidate[5.0]` | same |
| `test_coarse_proxy.py::test_features_batch_equals_per_candidate[None]` | same |
| `test_findoptimal_refactor.py::test_reconstruct_voxel_unchanged_by_refactor` | golden `reconstruct_voxel` result |
| `test_findoptimal_refactor.py::test_recorder_hook_leaves_reconstruct_voxel_bit_identical` | same golden |
| `test_findoptimal_refactor.py::test_rank_key_with_post_mc_costs_reproduces_golden` | same golden |
| `test_findoptimal_refactor.py::test_extra_candidates_returning_nothing_reproduces_golden` | same golden |
| `test_ray_overlap.py::test_rule_only_in_coarse_levels` | counts of a default search |
| `test_cpp_edge_rules.py::test_default_flags_unchanged_and_key_reaches_the_flags` | a new test that pins the default flags (would be edited with the flip) |

## `CppEdgeRules`, `CppOmegaInterp` and `BatchedScreen` all default on: 12 failed, 877 passed

The ten above, plus `test_profiling.py::test_instrumented_reconstruct_voxel_is_bit_identical_and_uninstalled` and
`test_batched_screen.py::test_config_keys_and_default_off` (a new test that pins the default).

Not tested here, and stored results that were computed through the cost function with the old frame lookup (about
0.3% of the peaks of a 180-frame, one-degree scan) would also move: the Phase D, level-0, ray-overlap and centroid
tables, the E0 baselines and the BFS pilot results. None were regenerated.
