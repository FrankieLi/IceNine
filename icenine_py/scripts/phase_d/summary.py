"""Write the Phase D data doc tables (results/phase_d_tables.md) from the stored JSONs.

Usage (from icenine_py/): uv run python scripts/phase_d/summary.py
Then: uv run python ../scripts/dev/sync_doc_tables.py --doc MIGRATION_HISTORY.md \
          --tables scripts/phase_d/results/phase_d_tables.md
"""

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))
import doc_tables as DT  # noqa: E402

R = HERE / "results"
EX = HERE.parents[2] / "Examples" / "Example2.ManyGrains" / "SimInput"


def _s(v):
    """int -> plain digits, float -> 4 significant digits without exponent where sensible."""
    if isinstance(v, int):
        return str(v)
    return f"{v:.1f}" if abs(v) >= 100 else f"{v:.3g}"


def main() -> None:
    st = json.loads((EX / "rand_500grains_1mm_neworient_s0_stats.json").read_text())
    full = json.loads((R / "render_full.json").read_text())
    pilot = json.loads((R / "render_pilot.json").read_text())
    clean = json.loads((R / "sanity_full_clean.json").read_text())
    real = json.loads((R / "sanity_full_realistic.json").read_text())

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
    for name, d in (("pilot (500 voxels)", pilot), ("full", full)):
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
    for name, d in (("clean", clean), ("realistic", real)):
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

    tables = {
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
