"""S4 summary: summary.json and tables.md (doc tables) from the Step 0-4 result files.

Usage (from icenine_py/): uv run python scripts/speed_parity/s4/summary.py
Then: uv run python ../scripts/dev/sync_doc_tables.py --doc MIGRATION_HISTORY.md \
          --tables benchmarks/s4_translation/tables.md
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "common"))
import common as C  # noqa: E402
import doc_tables  # noqa: E402
import stats  # noqa: E402

OUT = C.OUT
SEEDS = (350, 395, 504)  # BFS seed counts of the Phase D seed oracle (MIGRATION_HISTORY Phase D)


def load(name: str) -> Optional[Dict[str, Any]]:
    p = OUT / name
    return json.loads(p.read_text()) if p.exists() else None


def wilson_str(k: int, n: int) -> str:
    lo, hi = stats.wilson(k, n)
    return f"{lo:.3f}-{hi:.3f}"


def main() -> None:
    S: Dict[str, Any] = {}
    T: Dict[str, str] = {}
    s0 = load("step0_timing.json")
    law = load("step1_law.json")
    bit = load("step1_bitmap.json")
    eq_d = load("step2_equality_phase_d.json")
    eq_3 = load("step2_equality_three.json")
    tab = load("step2_table_phase_d.json")
    chk = load("step3_check.json")
    tim = load("step3_timing.json")
    full = load("step3_full_layer.json")
    model = load("step3_model.json")
    rce = load("step4_recost_equal.json")
    rct = load("step4_recost_timing.json")
    cost_rows: List[Dict[str, Any]] = []
    for p in sorted((OUT / "cache").glob("step1_cost_part*.json")):
        cost_rows += json.loads(p.read_text())["rows"]
    ref_idx = json.loads((C.CACHE / "ref" / "ref_index_phase_d.json").read_text())
    if cost_rows:  # the 1c rows of both parts, kept with the results (the parts live in cache/)
        (OUT / "step1_cost.json").write_text(json.dumps(dict(rows=cost_rows), indent=1) + "\n")

    # ---------------- Step 0
    if s0:
        scr = s0["screen"]
        med = lambda k: float(np.median([r[k] for r in scr]))  # noqa: E731
        sea = s0["search"]
        smed = lambda k: float(np.median([r[k] for r in sea]))  # noqa: E731
        S["step0"] = dict(
            label=s0["label"],
            preflight_load=s0["preflight"],
            n_screen_voxels=len(scr),
            screen_s=med("t_library_s"),
            screen_geometry_s=med("t_geometry_s") + med("t_bragg_s"),
            screen_stage_d_loop_s=med("t_stage_d_loop_s"),
            passers_median=med("passers"),
            passers_min=int(min(r["passers"] for r in scr)),
            passers_max=int(max(r["passers"] for r in scr)),
            recost_us_per_eval=float(np.median([r["us_per_eval"] for r in s0["recost"]])),
            recost_s_per_voxel=float(np.median([r["s"] for r in s0["recost"]])),
            search_wall_s=smed("wall_s"),
            search_screen_s=smed("screen_s"),
            search_recost_s=smed("recost_s"),
            search_recost_calls=smed("recost_calls"),
            search_discrete_other_s=smed("discrete_other_s"),
            search_quick_mc_s=smed("quick_mc_s"),
            search_quick_mc_level0_s=float(np.median([r["quick_mc_levels_s"][0] for r in sea])),
            search_rest_s=smed("rest_s"),
            search_err_deg_max=float(max(r["err_deg"] for r in sea)),
            n_search=len(sea),
            neighbour_median_s=s0["neighbour_median_s"],
            n_neighbour=len(s0["neighbour_s"]),
            images=s0["images"],
            peaks_per_orientation=s0["peaks_per_orientation"],
        )
        T["s4_step0"] = doc_tables.markdown_table(
            [
                dict(
                    part="BatchedScreen, whole level-0 screen of one voxel (chunk 32)",
                    value=S["step0"]["screen_s"],
                    unit="s",
                    n=len(scr),
                ),
                dict(
                    part="  of which torch geometry (Bragg solve, eta, stages B-C)",
                    value=S["step0"]["screen_geometry_s"],
                    unit="s",
                    n=len(scr),
                ),
                dict(
                    part="  of which the per-orientation stage-D loop",
                    value=S["step0"]["screen_stage_d_loop_s"],
                    unit="s",
                    n=len(scr),
                ),
                dict(
                    part="screen passers per voxel (of 43,974), median",
                    value=S["step0"]["passers_median"],
                    unit="orientations",
                    n=len(scr),
                ),
                dict(
                    part="re-cost of the passers (local evaluate), per voxel",
                    value=S["step0"]["recost_s_per_voxel"],
                    unit="s",
                    n=len(s0["recost"]),
                ),
                dict(
                    part="re-cost, per evaluate",
                    value=S["step0"]["recost_us_per_eval"],
                    unit="us",
                    n=len(s0["recost"]),
                ),
                dict(
                    part="no-start search (reconstruct_voxel), wall",
                    value=S["step0"]["search_wall_s"],
                    unit="s",
                    n=len(sea),
                ),
                dict(
                    part="  screens of all levels",
                    value=S["step0"]["search_screen_s"],
                    unit="s",
                    n=len(sea),
                ),
                dict(
                    part="  re-costs of all levels",
                    value=S["step0"]["search_recost_s"],
                    unit="s",
                    n=len(sea),
                ),
                dict(
                    part="  quick MC, all levels (level 0)",
                    value=f"{S['step0']['search_quick_mc_s']:.1f} ({S['step0']['search_quick_mc_level0_s']:.1f})",
                    unit="s",
                    n=len(sea),
                ),
                dict(
                    part="  rest of the discrete stage (spacing filter, loop)",
                    value=S["step0"]["search_discrete_other_s"],
                    unit="s",
                    n=len(sea),
                ),
                dict(
                    part="  FindOptimal, VarianceMinimizing, final",
                    value=S["step0"]["search_rest_s"],
                    unit="s",
                    n=len(sea),
                ),
                dict(
                    part="BFS neighbour step (CMA, 250 evaluations)",
                    value=S["step0"]["neighbour_median_s"],
                    unit="s",
                    n=len(s0["neighbour_s"]),
                ),
            ],
            ["part", "value", "unit", "n"],
            {"value": ".3g"},
        )
    # ---------------- Step 1
    if law and bit:
        tl = law["tallies"]
        iii = law["iii"]
        S["step1"] = dict(
            law={k: v for k, v in law.items() if k not in ("fails",)},
            bitmap=dict(
                all_equal=bit["all_equal"], centres=bit["centres_total"], images=len(bit["images"])
            ),
        )
        rows = [
            dict(
                check="1a float64: u(x + dx) - (u(x) + A dx), largest",
                n=law["peakset_compared_cf"],
                result=f"{law['dev64_max']:.1e} px",
                exceptions="none",
            ),
            dict(
                check="1a float32 (evaluate's stage-D input) against the float64 law, largest",
                n=law["peakset_compared_cf"],
                result=f"{law['dev32_max']:.1e} px",
                exceptions=f"float32 alone at x: {law['base32_max']:.1e} px",
            ),
            dict(
                check="1a peak sets (omegas, normals, frames, plane hits) identical",
                n=law["peakset_compared_cf"],
                result=f"{law['peakset_identical_cf']} of {law['peakset_compared_cf']} moves",
                exceptions="none",
            ),
            dict(
                check="1b simulator peak sets identical (count and frames)",
                n=law["peakset_compared_sim"],
                result=f"{law['peakset_identical_sim']} of {law['peakset_compared_sim']} moves",
                exceptions="none",
            ),
            dict(
                check="1b(i) integer shift: S(x') = S(x) + m (float32 simulator)",
                n=sum(tl["i32"][k] for k in ("equal", "edge", "float", "fail")),
                result=f"{tl['i32']['equal']} equal, {tl['i32']['fail']} failures",
                exceptions=f"{tl['i32']['edge']} at the detector edge ({tl['i32']['edge_off']} off it), {tl['i32']['float']} float flips",
            ),
            dict(
                check="1b(i) the same in float64 vertices",
                n=sum(tl["i64"][k] for k in ("equal", "edge", "float", "fail")),
                result=f"{tl['i64']['equal']} equal, {tl['i64']['fail']} failures",
                exceptions=f"{tl['i64']['edge']} at the detector edge, {tl['i64']['float']} float flips",
            ),
            dict(
                check="1b(ii) phase re-rasterised: S(x') = R(T + phi) + m (float32 simulator)",
                n=sum(tl["ii32"][k] for k in ("equal", "edge", "float", "fail")),
                result=f"{tl['ii32']['equal']} equal, {tl['ii32']['fail']} failures",
                exceptions=f"{tl['ii32']['edge']} at the detector edge ({tl['ii32']['edge_off']} off it), {tl['ii32']['float']} float flips (largest distance {law['float_flip_dist_max']['ii32']:.1e} px)",
            ),
            dict(
                check="1b(ii) the same in float64 vertices",
                n=sum(tl["ii64"][k] for k in ("equal", "edge", "float", "fail")),
                result=f"{tl['ii64']['equal']} equal, {tl['ii64']['fail']} failures",
                exceptions=f"{tl['ii64']['edge']} at the detector edge, {tl['ii64']['float']} float flips",
            ),
            dict(
                check="1b(iii) pure pixel-set translation S(x) + round(A dx) = S(x')",
                n=iii["n_spots"],
                result=f"{iii['n_equal']} equal ({100 * iii['n_equal'] / iii['n_spots']:.1f}%, Wilson {wilson_str(iii['n_equal'], iii['n_spots'])})",
                exceptions=f"{iii['n_edge_excluded']} edge spots excluded",
            ),
            dict(
                check="1e dilated bitmap and box bounds = pixel_radius_overlap_uint8",
                n=bit["centres_total"],
                result="equal" if bit["all_equal"] else "differs",
                exceptions="none (4 frames x 2 detectors, border band included)",
            ),
        ]
        T["s4_q1"] = doc_tables.markdown_table(rows, ["check", "n", "result", "exceptions"])
        T["s4_q1_pure"] = doc_tables.markdown_table(
            [
                dict(stat="spots compared (both on the detector, unclipped)", value=iii["n_spots"]),
                dict(stat="exactly equal", value=iii["n_equal"]),
                dict(
                    stat="symmetric difference, median / 95th pct / max (px)",
                    value=f"{iii['symdiff_px_median']:.0f} / {iii['symdiff_px_p95']:.0f} / {iii['symdiff_px_max']}",
                ),
                dict(
                    stat="symmetric difference / spot size, median / 95th pct",
                    value=f"{iii['rel_median']:.2f} / {iii['rel_p95']:.2f}",
                ),
                dict(
                    stat="spots with every differing pixel within 1 px (8-neighbour) of the other set",
                    value=iii["all_adjacent_to_other_set"],
                ),
                dict(
                    stat="spots with every differing pixel on its set's 4-boundary",
                    value=iii["all_on_boundary"],
                ),
                dict(
                    stat="spots leaving / entering the detector over the moves",
                    value=f"{law['spots_leaving_detector']} / {law['spots_entering_detector']}",
                ),
            ],
            ["stat", "value"],
        )
    if cost_rows:
        n = len(cost_rows)
        agg = dict(
            voxels=n,
            passers=int(sum(r["n_passers"] for r in cost_rows)),
            passers_recost_changed=int(sum(r["passers_recost_changed"] for r in cost_rows)),
            passers_any_field_changed=int(sum(r["passers_any_field_changed"] for r in cost_rows)),
            recost_abs_diff_max=float(max(r["recost_abs_diff_max"] for r in cost_rows)),
            quality_abs_diff_max=float(max(r["quality_abs_diff_max"] for r in cost_rows)),
            quality_abs_diff_median_of_medians=float(
                np.median([r["quality_abs_diff_median"] for r in cost_rows])
            ),
            cliques=int(sum(r["cliques"] for r in cost_rows)),
            cliques_same=int(sum(r["cliques_same_filtered_set"] for r in cost_rows)),
            best_same=int(sum(r["best_same"] for r in cost_rows)),
            top30_overlap_median=float(np.median([r["top30_overlap"] for r in cost_rows])),
            peaks_compared=int(sum(r["peaks_compared"] for r in cost_rows)),
            peaks_integer_vertices_changed=int(
                sum(r["peaks_integer_vertices_changed"] for r in cost_rows)
            ),
            disp_um_min=float(min(r["disp_um"] for r in cost_rows)),
            disp_um_max=float(max(r["disp_um"] for r in cost_rows)),
            kept_exact_median=float(np.median([r["kept_total_exact"] for r in cost_rows])),
            kept_exact_min=int(min(r["kept_total_exact"] for r in cost_rows)),
            kept_exact_max=int(max(r["kept_total_exact"] for r in cost_rows)),
        )
        S["step1c"] = agg
        T["s4_q1_cost"] = doc_tables.markdown_table(
            [
                dict(
                    measure="passers re-costed, exact and translated (20 voxels)",
                    value=agg["passers"],
                ),
                dict(
                    measure="passers with any OverlapInfo field changed",
                    value=f"{agg['passers_any_field_changed']} ({100 * agg['passers_any_field_changed'] / agg['passers']:.1f}%)",
                ),
                dict(
                    measure="passers with the re-cost (1 - confidence) changed",
                    value=f"{agg['passers_recost_changed']} ({100 * agg['passers_recost_changed'] / agg['passers']:.1f}%)",
                ),
                dict(measure="largest re-cost change", value=f"{agg['recost_abs_diff_max']:.3f}"),
                dict(
                    measure="largest quality change (final cost = 1 - quality)",
                    value=f"{agg['quality_abs_diff_max']:.3f}",
                ),
                dict(
                    measure="cliques whose spacing-filtered set is unchanged",
                    value=f"{agg['cliques_same']} of {agg['cliques']} ({100 * agg['cliques_same'] / agg['cliques']:.1f}%)",
                ),
                dict(
                    measure="voxels whose best filtered candidate is unchanged",
                    value=f"{agg['best_same']} of {n} (Wilson {wilson_str(agg['best_same'], n)})",
                ),
                dict(
                    measure="top-30 filtered candidates in common, median over voxels",
                    value=f"{agg['top30_overlap_median']:.1f} of 30",
                ),
            ],
            ["measure", "value"],
        )
    # ---------------- Step 2
    if eq_d and eq_3:
        rows = []
        for name, e in (("Phase D (20 voxels)", eq_d), ("ThreeVoxels (3 voxels)", eq_3)):
            t = e["totals"]
            pairs = e["n_voxels"] * e["n_orient"]
            for prec in ("f64", "f32"):
                lst = [x for x in e["listed"] if x["table"] == prec]
                dmax = max(
                    (min(i["dist_to_boundary"] for i in x["flipped"]) for x in lst if x["flipped"]),
                    default=None,
                )
                rows.append(
                    dict(
                        set=name,
                        table="float64" if prec == "f64" else "float32",
                        pairs=pairs,
                        peak_overlap_differs=t[f"{prec}_peak_overlap_diff"],
                        pass_differs=t[f"{prec}_pass_diff"],
                        peak_coordinate_flips=t[f"{prec}_peak_coord_flips"],
                        unexplained=sum(1 for x in lst if not x["flipped"]),
                        largest_boundary_distance=dmax,
                    )
                )
        T["s4_equality"] = doc_tables.markdown_table(
            rows,
            [
                "set",
                "table",
                "pairs",
                "peak_overlap_differs",
                "pass_differs",
                "peak_coordinate_flips",
                "unexplained",
                "largest_boundary_distance",
            ],
            {"largest_boundary_distance": ".1e"},
        )
        S["step2"] = dict(rows=rows, table=tab)
    # ---------------- Step 3
    if tim:
        arms = tim["arms"]
        S["step3"] = dict(
            label=tim["label"],
            preflight_load=tim["preflight_load"],
            table_build_s=tim["table_build_s"],
            frame_sort_s=tim["frame_sort_s"],
            table_bytes=tim["table_bytes"],
            bitmap_bytes=tim["bitmap_bytes"],
            out_bytes_2000=tim["out_bytes_2000"],
            peak_rss_gb=tim["peak_rss_gb"],
            all_equal=all(a["equal"] for a in arms.values()),
            arms={
                k: dict(
                    median_s=v["median_s"],
                    ns_per_test=v["ns_per_test"],
                    projected_full_layer_s=v["projected_full_layer_s"],
                    equal=v["equal"],
                )
                for k, v in arms.items()
            },
            stream=tim["stream"],
            mps=tim["mps"],
            check_all_equal=chk["all_equal"] if chk else None,
        )
        main_rows = []
        for k, v in arms.items():
            if v["n"] != 2000 or "th=" in k:
                continue
            main_rows.append(
                dict(
                    arm=k,
                    threads=v["threads"],
                    median_s=v["median_s"],
                    ns_per_test=v["ns_per_test"],
                    pairs_per_s=v["pairs_per_s"],
                    full_layer_s=v["projected_full_layer_s"],
                    equal=v["equal"],
                )
            )
        T["s4_orderings"] = doc_tables.markdown_table(
            main_rows,
            ["arm", "threads", "median_s", "ns_per_test", "pairs_per_s", "full_layer_s", "equal"],
            {"median_s": ".3g", "ns_per_test": ".3g", "pairs_per_s": ".3g", "full_layer_s": ".3g"},
        )
        thr_rows = []
        for base in (
            "O2 counts u8 b=32FZ",
            "O2 pass u8 b=32FZ",
            "O1 pass packed",
            "O3 counts u8 tile=256",
            "O5 counts u8 tile=256",
        ):
            row: Dict[str, Any] = dict(arm=base)
            one = arms.get(base) or arms.get(base.replace("u8", "u8"))
            for th in (1, 2, 4, 6, 8, 10):
                key = base if th == 1 else f"{base} th={th}"
                if key in arms:
                    row[f"{th}"] = arms[key]["median_s"]
            if one and f"{base} th=10" in arms:
                row["speedup_10"] = one["median_s"] / arms[f"{base} th=10"]["median_s"]
            thr_rows.append(row)
        T["s4_threads"] = doc_tables.markdown_table(
            thr_rows,
            ["arm", "1", "2", "4", "6", "8", "10", "speedup_10"],
            {c: ".3g" for c in ("1", "2", "4", "6", "8", "10", "speedup_10")},
        )
        lin_rows = []
        for base in (
            "O1 pass packed",
            "O2 counts u8 b=32FZ",
            "O3 counts u8 tile=256",
            "O5 counts u8 tile=256",
            "O2 counts u8 b=32FZ th=10",
        ):
            r = dict(arm=base)
            for n in (500, 1000, 2000):
                key = base if n == 2000 else f"{base} n={n}"
                if key in arms:
                    r[f"s_{n}"] = arms[key]["median_s"]
                    r[f"us_per_voxel_{n}"] = 1e6 * arms[key]["median_s"] / n
            lin_rows.append(r)
        T["s4_linearity"] = doc_tables.markdown_table(
            lin_rows,
            [
                "arm",
                "s_500",
                "s_1000",
                "s_2000",
                "us_per_voxel_500",
                "us_per_voxel_1000",
                "us_per_voxel_2000",
            ],
            {
                c: ".3g"
                for c in (
                    "s_500",
                    "s_1000",
                    "s_2000",
                    "us_per_voxel_500",
                    "us_per_voxel_1000",
                    "us_per_voxel_2000",
                )
            },
        )
        var_rows = []
        for k, v in tim["stream"].items():
            var_rows.append(
                dict(
                    arm=f"O2 streamed (no table), {k}",
                    median_s=float(np.median(v["times"])),
                    peaks_s=v["t_peaks_s"],
                    kernel_s=v["t_kernel_s"],
                    frames_changed=v.get("frames_changed_vs_chunk32"),
                    counts_differ=0 if v["equal_table_kernel"] else "some (patch)",
                )
            )
        for k, v in tim["mps"].items():
            var_rows.append(
                dict(
                    arm=f"O2 on MPS (float32 offsets), {k}",
                    median_s=float(np.median(v["times"])),
                    counts_differ=v["po_diff"],
                    pass_differ=v["pass_diff"],
                )
            )
        T["s4_variants"] = doc_tables.markdown_table(
            var_rows,
            [
                "arm",
                "median_s",
                "peaks_s",
                "kernel_s",
                "frames_changed",
                "counts_differ",
                "pass_differ",
            ],
            {"median_s": ".3g", "peaks_s": ".3g", "kernel_s": ".3g"},
        )
    if full:
        S["step3_full"] = full
        T["s4_full_layer"] = doc_tables.markdown_table(
            [
                dict(
                    item="ordering",
                    value=f"{full['kind'].upper()} counts, {full['layout']}, {full['threads']} threads, voxel tiles of {full['voxel_tile']}",
                ),
                dict(
                    item="voxels x orientations", value=f"{full['n_voxels']} x {full['n_orient']}"
                ),
                dict(
                    item="wall time (s), kernel + pass reduction",
                    value=f"{full['wall_s']:.2f} ({full['kernel_s']:.2f} + {full['reduce_s']:.2f})",
                ),
                dict(
                    item="BatchedScreen for the same work, single worker (h)",
                    value=(
                        f"{full['n_voxels'] * S['step0']['screen_s'] / 3600:.2f}" if s0 else "n/a"
                    ),
                ),
                dict(
                    item="passers (voxel, orientation) pairs",
                    value=f"{full['passers_total']} ({100 * full['passers_total'] / (full['n_voxels'] * full['n_orient']):.1f}%)",
                ),
                dict(
                    item="passers per voxel, median (min-max)",
                    value=f"{full['passers_per_voxel_median']:.0f} ({full['passers_per_voxel_min']}-{full['passers_per_voxel_max']})",
                ),
                dict(
                    item="passer list as int32 pairs / as a bitset / counts as uint8 (GB)",
                    value=f"{full['sparse_list_bytes_int32_pairs'] / 1e9:.2f} / {full['pass_bitset_bytes'] / 1e9:.3f} / {full['counts_matrix_bytes_uint8'] / 1e9:.2f}",
                ),
                dict(item="peak RSS of the run (GB)", value=f"{full['peak_rss_gb']:.2f}"),
                dict(item="label", value=full["label"]),
            ],
            ["item", "value"],
        )
    if model:
        S["step3_model"] = model
    o4 = load("step3_o4.json")
    if model and o4:
        S["step3_o4"] = {k: v for k, v in o4.items() if k != "arms"}
        S["step3_o4"]["arms"] = {
            k: dict(median_s=v["median_s"], equal=v["equal"]) for k, v in o4["arms"].items()
        }
        T["s4_o4"] = doc_tables.markdown_table(
            [
                dict(
                    quantity="bright rate per O2 test (patch, sampled peaks)",
                    model=model["bright_rate_per_test"],
                    measured=None,
                ),
                dict(
                    quantity="voxels per distinct pixel hit",
                    model=model["voxels_per_distinct_pixel"],
                    measured=None,
                ),
                dict(
                    quantity="bright pixels visited per peak",
                    model=model["bright_distinct_pixels_per_peak_mean"],
                    measured=o4["bright_pixels_per_peak"],
                ),
                dict(
                    quantity="candidate voxels tested per bright pixel",
                    model=model["voxels_per_distinct_pixel"],
                    measured=o4["candidates_per_bright_pixel"],
                ),
                dict(
                    quantity="cost per bright pixel (ns)",
                    model="5 / 10 / 20 (assumed)",
                    measured=o4["ns_per_bright_pixel_1thread"],
                ),
                dict(
                    quantity="O2 time / O4 time, 1 thread",
                    model=f"{model['o2_over_o4_tbp5']:.2f} / {model['o2_over_o4_tbp10']:.2f} / {model['o2_over_o4_tbp20']:.2f}",
                    measured=1.0 / o4["o4_over_o2_1thread"],
                ),
                dict(
                    quantity="O4 / O2 wall, 1 thread (s)",
                    model=None,
                    measured=f"{o4['arms']['O4 counts packed th=1']['median_s']:.2f} / {o4['arms']['O2 counts packed th=1']['median_s']:.2f}",
                ),
                dict(
                    quantity="O4 / O2 wall, 10 threads (s)",
                    model=None,
                    measured=f"{o4['arms']['O4 counts packed th=10']['median_s']:.2f} / {o4['arms']['O2 counts packed th=10']['median_s']:.2f}",
                ),
            ],
            ["quantity", "model", "measured"],
            {"model": ".3g", "measured": ".3g"},
        )
    # ---------------- Step 4
    if rce:
        S["step4_equal"] = {k: v for k, v in rce.items() if k != "rows"}
    if rct:
        S["step4_timing"] = rct
    if rce and rct:
        T["s4_recost"] = doc_tables.markdown_table(
            [
                dict(
                    item="today: evaluate per (orientation, voxel) (us)",
                    value=f"{rct['evaluate_us']:.0f}",
                    n=rct["n_evaluate"],
                ),
                dict(
                    item="orientation-major: per pair, float64 vertices + stage-D C call (us)",
                    value=f"{rct['pair_us']:.1f}",
                    n=rct["pairs"],
                ),
                dict(
                    item="  the same with one vertex product per orientation (us)",
                    value=f"{rct['pair_us_batched_vertices']:.1f}",
                    n=rct["pairs"],
                ),
                dict(
                    item="  set-up per orientation: peaks, A and c (us)",
                    value=f"{rct['setup_us_per_orientation']:.0f}",
                    n=rct["n_orient"],
                ),
                dict(
                    item="speed-up per pair",
                    value=f"{rct['speedup_pair_vs_evaluate']:.1f}x",
                    n=rct["pairs"],
                ),
                dict(
                    item="equality: pairs compared with the Step 0 per-voxel evaluate",
                    value=rce["pairs"],
                    n=20,
                ),
                dict(
                    item="  pairs with any OverlapInfo field different",
                    value=f"{rce['any_field_diff']} ({100 * rce['any_field_diff'] / rce['pairs']:.1f}%)",
                    n=rce["pairs"],
                ),
                dict(
                    item="  pairs with the re-cost (1 - confidence) different",
                    value=rce["recost_diff"],
                    n=rce["pairs"],
                ),
                dict(
                    item="  largest quality difference",
                    value=f"{rce['quality_abs_diff_max']:.4f}",
                    n=rce["pairs"],
                ),
                dict(
                    item="  sampled differing pairs explained by an integer-vertex flip",
                    value=f"{rce['explain_integer_vertex_differs']} of {rce['explain_sampled']} (largest distance {rce['explain_flip_dist_max']:.1e} px)",
                    n=rce["explain_sampled"],
                ),
            ],
            ["item", "value", "n"],
        )
    # ---------------- projections
    if s0 and full and rct and tim:
        st0 = S["step0"]
        V = full["n_voxels"]
        T_ns = st0["search_wall_s"]
        L0 = st0["screen_s"] + st0["recost_s_per_voxel"]
        rest = T_ns - L0
        T_nb = st0["neighbour_median_s"]
        t_pair = rct["pair_us"] * 1e-6
        T_G = full["wall_s"]
        # single-worker global screen: the measured 1-thread O2 time of the 2,000-voxel patch,
        # scaled linearly to the layer (the 10-thread full-layer wall is T_G)
        T_G1 = tim["arms"]["O2 counts u8 b=32FZ"]["projected_full_layer_s"]
        pairs = full["passers_total"]
        pv = full["passers_per_voxel_median"]
        prj = []
        for s in SEEDS:
            a_ = dict(
                scenario="(a) today's BFS",
                seeds=s,
                screen=s * st0["screen_s"],
                recost=s * st0["recost_s_per_voxel"],
                refine=s * rest,
                neighbours=(V - s) * T_nb,
            )
            b_ = dict(
                scenario="(b) global screen and re-cost, per-voxel refinement",
                seeds=s,
                screen=T_G1,
                recost=pairs * t_pair,
                refine=V * rest,
                neighbours=0.0,
            )
            c_ = dict(
                scenario="(c) hybrid: BFS, seeds take the global level-0 list",
                seeds=s,
                screen=T_G1,
                recost=s * pv * t_pair,
                refine=s * rest,
                neighbours=(V - s) * T_nb,
            )
            for r in (a_, b_, c_):
                r["total_h"] = (r["screen"] + r["recost"] + r["refine"] + r["neighbours"]) / 3600
                tot = r["total_h"] * 3600
                r["share_screen"] = r["screen"] / tot
                r["share_recost"] = r["recost"] / tot
                r["share_refine"] = r["refine"] / tot
                r["share_neighbours"] = r["neighbours"] / tot
                prj.append(r)
        S["projection"] = dict(
            inputs=dict(
                V=V,
                nostart_s=T_ns,
                level0_screen_s=st0["screen_s"],
                level0_recost_s=st0["recost_s_per_voxel"],
                rest_of_search_s=rest,
                neighbour_s=T_nb,
                recost_pair_us=rct["pair_us"],
                global_screen_wall_s=T_G,
                global_screen_threads=full["threads"],
                global_screen_single_worker_equiv_s=T_G1,
                passers_total=pairs,
                passers_per_voxel_median=pv,
            ),
            rows=prj,
        )
        T["s4_projection"] = doc_tables.markdown_table(
            prj,
            [
                "scenario",
                "seeds",
                "total_h",
                "share_screen",
                "share_recost",
                "share_refine",
                "share_neighbours",
            ],
            {
                "total_h": ".3g",
                "share_screen": ".3f",
                "share_recost": ".3f",
                "share_refine": ".3f",
                "share_neighbours": ".3f",
            },
        )
    if "projection" in S and tim:
        P = S["projection"]["inputs"]
        arms = tim["arms"]
        o2 = arms["O2 counts u8 b=32FZ"]["median_s"]
        rows = {r["scenario"][:3] + str(r["seeds"]): r for r in S["projection"]["rows"]}
        S["derived"] = dict(
            batchedscreen_layer_h=P["V"] * P["level0_screen_s"] / 3600,
            screen_speedup_single_thread=P["V"]
            * P["level0_screen_s"]
            / P["global_screen_single_worker_equiv_s"],
            screen_speedup_10_threads_wall=P["V"]
            * P["level0_screen_s"]
            / P["global_screen_wall_s"],
            o1_counts_over_o2=arms["O1 counts u8"]["median_s"] / o2,
            o1_pass_packed_over_o2=arms["O1 pass packed"]["median_s"] / o2,
            o3_tile2000_over_o2=arms["O3 counts u8 tile=2000"]["median_s"] / o2,
            o5_tile256_over_o2=arms["O5 counts u8 tile=256"]["median_s"] / o2,
            o2_pass_over_counts=arms["O2 pass u8 b=32FZ"]["median_s"] / o2,
            mps_pairs_patch=int(tim["n_orient"] * 2000),
            bitmap_u8_det0_mb=tim["bitmap_bytes"]["u8_det0"] / 1e6,
            bitmap_packed_det_mb=tim["bitmap_bytes"]["packed_both"] / 2e6,
            level0_screen_recost_s=P["level0_screen_s"] + P["level0_recost_s"],
            mps_over_o2_1thread=min(float(np.median(v["times"])) for v in tim["mps"].values()) / o2,
            stream_b32_over_table=float(np.median(tim["stream"]["b=32FZ"]["times"])) / o2,
            level0_share_of_nostart=(P["level0_screen_s"] + P["level0_recost_s"]) / P["nostart_s"],
            quick_mc_level0_share_of_nostart=S["step0"]["search_quick_mc_level0_s"]
            / P["nostart_s"],
            seed_saving_s=P["level0_screen_s"]
            + P["level0_recost_s"]
            - P["passers_per_voxel_median"] * P["recost_pair_us"] * 1e-6,
            c_saving_h={s: rows[f"(a){s}"]["total_h"] - rows[f"(c){s}"]["total_h"] for s in SEEDS},
            c_saving_fraction={
                s: 1 - rows[f"(c){s}"]["total_h"] / rows[f"(a){s}"]["total_h"] for s in SEEDS
            },
            b_over_a={s: rows[f"(b){s}"]["total_h"] / rows[f"(a){s}"]["total_h"] for s in SEEDS},
            recost_layer_h_orientation_major=P["passers_total"] * P["recost_pair_us"] * 1e-6 / 3600,
            recost_layer_h_evaluate=P["passers_total"]
            * S["step0"]["recost_us_per_eval"]
            * 1e-6
            / 3600,
            pass_rate=P["passers_total"] / (P["V"] * full["n_orient"]),
        )
    (OUT / "summary.json").write_text(json.dumps(S, indent=1, default=float) + "\n")
    doc_tables.write_tables(OUT / "tables.md", T)
    print("tables:", sorted(T))


if __name__ == "__main__":
    main()
