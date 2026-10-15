# S4: translation correctness and loop ordering (results)

Plan: `docs/todo_translation_screen_prototype.md`. Write-up: MIGRATION_HISTORY "S4: translation correctness and
loop ordering". Scripts: `scripts/speed_parity/s4/` (run from `icenine_py/`).

| File | Content | Made by |
|---|---|---|
| `summary.json`, `tables.md` | every headline number and the doc tables | `summary.py` |
| `step0_timing.json`, `preflight_step0*.json` | BatchedScreen split, passers, re-cost, no-start search split, neighbour step (contended) | `step0_timing.py` |
| `step1_law.json` | 1a (affine law, peak sets) and 1b (lit pixels: integer shift, phase re-raster, pure translation) | `translation_correctness.py law` |
| `step1_bitmap.json` | 1e (dilated bitmap against the C pixel-radius test) | `translation_correctness.py bitmap` |
| `step1_cost.json` | 1c (cost-level effect of a pure pixel-set translation), both parts merged | `translation_correctness.py cost`, `summary.py` |
| `step2_table_phase_d.json`, `step2_equality_*.json` | the peak table and its equality with BatchedScreen (every flip listed) | `screen_kernel.py` |
| `step3_check.json` | every ordering against the numpy reference; streamed O2 and MPS flips | `orderings.py check` |
| `step3_timing.json`, `preflight_step3*.json` | ordering timings on the 2,000-voxel patch (contended) | `orderings.py time` |
| `step3_full_layer.json` | the best ordering on the whole layer | `orderings.py full` |
| `step3_model.json`, `step3_o4.json` | O4 cost model inputs and the O4 timing | `orderings.py model`, `time_o4` |
| `step3_rule_rates.json` | screen pass counts under the default and the C++ qualified-peak rule (not timed) | `orderings.py rules` |
| `step3_o3_threads.json`, `preflight_step3_o3*.json` | O3 tile 2,000 and 256 against O2 at 1 and 10 threads (contended) | `orderings.py time_o3` |
| `step4_recost_equal.json`, `step4_recost_timing.json` | orientation-major re-cost: equality and timing | `recost.py` |
| `step4_recost_memory.json` | bytes per orientation of a stored local-cost table (not timed) | `recost.py memory` |

`cache/` (ignored) holds the Step 0 reference (`ref/`), the peak tables, the bit-packed dilated bitmaps, the
compiled kernel and the logs. Order of a full rerun: `step0_reference.py --set three`, `--set phase_d --workers 4`,
`translation_correctness.py law | bitmap | cost --part 0/1 --parts 2`, `screen_kernel.py table/equal`, then alone:
`step0_timing.py`, `orderings.py time`, `orderings.py full`, `recost.py time`, `orderings.py time_o4`,
`orderings.py time_o3`; then `orderings.py check`, `orderings.py rules`, `orderings.py model`, `recost.py equal`,
`recost.py memory`, `summary.py`.
