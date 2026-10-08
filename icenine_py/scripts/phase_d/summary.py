"""Write the Phase D data doc tables (results/phase_d_tables.md) from the stored JSONs.

Usage (from icenine_py/): uv run python scripts/phase_d/summary.py
Then copy the tables into the "Phase D data" section of MIGRATION_HISTORY.md. sync_doc_tables.py
fails on the whole file (other sections have markers without generated blocks), so cut that
section into a scratch file, run
  uv run python ../scripts/dev/sync_doc_tables.py --doc SCRATCH.md \
      --tables scripts/phase_d/results/phase_d_tables.md
on it and splice it back.
"""

import json
import sys
from pathlib import Path
from typing import Union

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))
import doc_tables as DT  # noqa: E402

R = HERE / "results"
EX = HERE.parents[2] / "Examples" / "Example2.ManyGrains" / "SimInput"


def _s(v: Union[int, float]) -> str:
    """int -> plain digits, float -> 4 significant digits without exponent where sensible."""
    if isinstance(v, int):
        return str(v)
    return f"{v:.1f}" if abs(v) >= 100 else f"{v:.3g}"


def main() -> None:
    st = json.loads((EX / "rand_500grains_1mm_neworient_s0_stats.json").read_text())
    full = json.loads((R / "render_full.json").read_text())
    pilot = json.loads((R / "render_pilot.json").read_text())
    q16 = json.loads((R / "render_full_q16.json").read_text())
    clean = json.loads((R / "sanity_full_clean.json").read_text())
    real = json.loads((R / "sanity_full_realistic.json").read_text())
    real16 = json.loads((R / "sanity_full_realistic_q16.json").read_text())

    sample_rows = [
        {"quantity": "voxels", "value": _s(st["n_voxels"])},
        {"quantity": "grains (distinct orientations)", "value": _s(st["n_grains"])},
        {
            "quantity": "voxels per grain min / median / max",
            "value": "{} / {:g} / {}".format(*st["grain_size_min_median_max"]),
        },
        {"quantity": "adjacent grain pairs", "value": _s(st["n_adjacent_pairs"])},
        {"quantity": "redrawn grains (of 497)", "value": _s(st["n_redrawn"])},
        {"quantity": "min new vs any old (deg)", "value": _s(st["min_new_vs_old_deg"])},
        {"quantity": "min new vs adjacent new (deg)", "value": _s(st["min_neighbour_new_deg"])},
        {"quantity": "min old vs old (deg)", "value": _s(st["min_old_vs_old_deg"])},
        {
            "quantity": "mic round-trip max error (deg)",
            "value": _s(st["roundtrip_max_misorientation_deg"]),
        },
    ]
    render_rows = []
    for name, d in (("pilot (500 voxels)", pilot), ("full", full), ("full Q-max 16", q16)):
        render_rows.append(
            {
                "run": name,
                "voxels": d["n_voxels"],
                "workers": d["workers"],
                "wall_s": d["wall_s"],
                "max_worker_sim_s": d["worker_t_sim_s_max"],
                "max_worker_write_noise_s": d["worker_t_write_noise_s_max"],
                "max_worker_rss_gb": d["worker_peak_rss_gb_max"],
            }
        )
    frames_rows = [
        {"quantity": "frames (180 omega x 2 detectors)", "value": _s(full["n_frames"])},
        {
            "quantity": "lit pixels per frame, median",
            "value": _s(full["lit_clean_per_frame_median"]),
        },
        {
            "quantity": "lit pixels per frame, min",
            "value": _s(full["lit_clean_per_frame_min_max"][0]),
        },
        {
            "quantity": "lit pixels per frame, max",
            "value": _s(full["lit_clean_per_frame_min_max"][1]),
        },
        {"quantity": "lit pixels, clean total", "value": _s(full["lit_clean_total"])},
        {"quantity": "lit pixels, realistic total", "value": _s(full["lit_real_total"])},
        {"quantity": "spots (connected components), all frames", "value": _s(full["spots_total"])},
        {"quantity": "spots missed", "value": _s(full["missed_total"])},
        {"quantity": "hot pixels added", "value": _s(full["hot_total"])},
        {"quantity": "blobs added", "value": _s(full["blob_total"])},
    ]
    cost_rows = []
    for name, d in (("clean", clean), ("realistic", real), ("realistic_q16", real16)):
        rows = d["cost"]["rows"]
        row = {"images": name, "voxels": len(rows)}
        for k, lab in (
            ("q_true", "q_true"),
            ("q_0.25deg_mean", "q_0.25"),
            ("q_0.5deg_mean", "q_0.5"),
            ("q_1.0deg_max", "q_1.0_max"),
        ):
            v = np.array([x[k] for x in rows])
            row[lab] = "{:.2f} - {:.2f}".format(v.min(), v.max())
        row["min_gap_true_minus_1deg"] = min(x["q_true"] - x["q_1.0deg_max"] for x in rows)
        cost_rows.append(row)

    mm_rows = []
    for v in ("clean", "realistic", "realistic_q16"):
        m = json.loads((R / f"memmap_check_{v}.json").read_text())
        mm_rows.append(
            {
                "images": v,
                "evaluations": m["n_evaluations"],
                "identical": str(m["all_overlap_fields_identical"]),
                "dense_rss_gb": m["dense_peak_rss_gb"],
                "memmap_rss_gb": m["memmap_peak_rss_gb"],
                "dense_load_s": m["dense_load_s"],
                "memmap_load_s": m["memmap_load_s"],
            }
        )
    q16_rows = [
        {"quantity": "render wall time, 10 workers (s)", "value": _s(q16["wall_s"])},
        {"quantity": "max worker sim time (s)", "value": _s(q16["worker_t_sim_s_max"])},
        {"quantity": "max worker RSS (GB)", "value": _s(q16["worker_peak_rss_gb_max"])},
        {"quantity": "lit pixels, clean total", "value": _s(q16["lit_clean_total"])},
        {"quantity": "lit pixels, realistic total", "value": _s(q16["lit_real_total"])},
        {
            "quantity": "lit fraction of all pixels (Q-max 8: %s)"
            % _s(full["lit_clean_fraction_of_pixels"]),
            "value": _s(q16["lit_clean_fraction_of_pixels"]),
        },
        {"quantity": "spots, all frames", "value": _s(q16["spots_total"])},
        {
            "quantity": "spot size px: median / p95 / max",
            "value": "{:g} / {:g} / {}".format(
                q16["spot_px_pooled_median"], q16["spot_px_pooled_p95"], q16["spot_px_max"]
            ),
        },
        {"quantity": "spots larger than 2000 px", "value": _s(q16["spot_px_n_gt_2000"])},
    ]
    size_rows = [
        {
            "quantity": "spot size px: median / p95 / max (Q-max 8)",
            "value": "{:g} / {:g} / {}".format(
                full["spot_px_pooled_median"], full["spot_px_pooled_p95"], full["spot_px_max"]
            ),
        },
        {
            "quantity": "lit fraction of all pixels (Q-max 8)",
            "value": _s(full["lit_clean_fraction_of_pixels"]),
        },
    ]
    pr = json.loads((R / "bfs_timing_probe_clean.json").read_text())
    probe_rows = [
        {"quantity": "no-start seed voxel (reconstruct_voxel), n", "value": str(len(pr["seed_s"]))},
        {
            "quantity": "  wall time per seed (s)",
            "value": ", ".join("%.0f" % x for x in pr["seed_s"]),
        },
        {
            "quantity": "  final error (deg)",
            "value": ", ".join("%.3f" % x for x in pr["seed_err_deg"]),
        },
        {
            "quantity": "neighbour step (local_optimization from 0.3 deg off), n",
            "value": str(pr["n_neighbours"]),
        },
        {
            "quantity": "  median / mean wall time (s)",
            "value": "%.1f / %.1f" % (pr["neighbour_median_s"], pr["neighbour_mean_s"]),
        },
        {
            "quantity": "  median final error (deg)",
            "value": "%.3f" % pr["neighbour_err_deg_median"],
        },
    ]
    tables = {
        "phase_d_probe": DT.markdown_table(probe_rows, ["quantity", "value"]),
        "phase_d_memmap": DT.markdown_table(
            mm_rows,
            [
                "images",
                "evaluations",
                "identical",
                "dense_rss_gb",
                "memmap_rss_gb",
                "dense_load_s",
                "memmap_load_s",
            ],
            {
                "dense_rss_gb": ".2f",
                "memmap_rss_gb": ".2f",
                "dense_load_s": ".1f",
                "memmap_load_s": ".1f",
            },
        ),
        "phase_d_q16": DT.markdown_table(q16_rows, ["quantity", "value"]),
        "phase_d_sizes": DT.markdown_table(size_rows, ["quantity", "value"]),
        "phase_d_sample": DT.markdown_table(sample_rows, ["quantity", "value"], {}),
        "phase_d_render": DT.markdown_table(
            render_rows,
            [
                "run",
                "voxels",
                "workers",
                "wall_s",
                "max_worker_sim_s",
                "max_worker_write_noise_s",
                "max_worker_rss_gb",
            ],
            {
                "wall_s": ".1f",
                "max_worker_sim_s": ".1f",
                "max_worker_write_noise_s": ".1f",
                "max_worker_rss_gb": ".2f",
            },
        ),
        "phase_d_frames": DT.markdown_table(frames_rows, ["quantity", "value"], {}),
        "phase_d_cost_check": DT.markdown_table(
            cost_rows,
            [
                "images",
                "voxels",
                "q_true",
                "q_0.25",
                "q_0.5",
                "q_1.0_max",
                "min_gap_true_minus_1deg",
            ],
            {"min_gap_true_minus_1deg": ".2f"},
        ),
    }
    DT.write_tables(R / "phase_d_tables.md", tables)
    print((R / "phase_d_tables.md").read_text())


if __name__ == "__main__":
    main()
