---
name: mic-maps
description: Publish an Artifact page of orientation maps and misorientation maps (against a reference) for one or more .mic reconstructions, so spatial patterns of reconstruction differences can be seen. Use when comparing reconstructions of 500 or more voxels per layer (BFS variants, optimizers, image sets, Python vs C++), or when asked for orientation or misorientation maps.
user-invocable: true
allowed-tools: Bash(uv run *), Bash(ls *), Bash(find *), Bash(head *), Read, Write, Artifact
argument-hint: "<reference.mic> <recon.mic>... | description of the runs to compare"
---

Tools: `icenine_py/scripts/plot_mic.py` (one figure) and `icenine_py/scripts/mic_maps_report.py`
(spec → `index.html` + `maps/*.png` + `stats.json`). Both run from `icenine_py/` with `uv run`.

1. **Find the files.** Pick the reference and the reconstructions to compare.
   - The reference is the truth if one exists (for Phase D: `Examples/Example2.ManyGrains/SimInput/rand_500grains_1mm_neworient_s0.mic`
     with the region's `region_index.npy`). Otherwise use the reconstruction the user names, and say on the page
     that the reference is a reconstruction, not the truth.
   - Check the size with `head -3` and a voxel count. If a layer has fewer than 500 voxels, maps show little:
     run `plot_mic.py --stats` and give the user the numbers in chat instead, unless they asked for the page.
   - Check that all files are on the same grid (same side length and generation). If the reconstruction covers a
     region of a larger reference, pass the region's index file as `reference.index`, or let position matching
     do it (reference voxels without a match are drawn grey and the page says how many).
   - Multi-layer files: give `z` per figure (one figure per layer). The tools refuse multi-layer input without it.
   - Crystal symmetry is cubic only. For another symmetry, stop and tell the user the tools need extending.
2. **Group the figures.** One figure per condition (for example per optimizer and image set). Each figure is the
   reference plus at most 4 reconstructions, in the same panel order in every figure. Labels in plain words, no
   run-folder shorthand (write "Reseed, voxel rule", not "v2"). Per-figure `reference` overrides the default.
3. **Write the spec** to the scratchpad (format in the docstring of `mic_maps_report.py`): a 2–4 word `title`, a
   one-paragraph `summary` saying what is compared and on what, figure `title` and `caption`, and `notes` with
   caveats (simulation Q-max and reconstruction Q-max, arm mismatches, one region or seed).
4. **Build:** `cd icenine_py && uv run python scripts/mic_maps_report.py <spec> --out-dir <scratchpad>/<name>`.
   The page it writes already follows the Artifact page contract (title, colour tokens, dark mode, phone width).
5. **Check once.** Compare the printed wrong counts with any scored numbers for the same runs (e.g. a run's
   `voxels.npz` or `score.json`); a mismatch means a wrong reference, index or layer, so fix that first. Look at
   one PNG with Read. Do not loop on cosmetics.
6. **Publish** with the Artifact tool: `file_path` = `<out-dir>/index.html`, `files` = `{"maps/<name>.png":
   "<out-dir>/maps/<name>.png", ...}` for every figure, `icon: "map"` on the first publish, and a one-sentence
   `description`. To update the page later, rebuild into the same out-dir and publish the same `file_path`.
7. **Reply** with the link, one line per figure giving the wrong count of each panel, and any caveat (unmatched
   voxels, reference not the truth). If the maps belong with a study, offer to copy the PNGs into the run's
   `benchmarks/.../maps/` folder and commit them.
