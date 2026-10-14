"""Task A step 0: summary.json and tables.md from results.json.gz of parity_compare.py.

  uv run python scripts/eval_reduction/parity_summary.py --dir DIR \
      --out benchmarks/eval_reduction/A/step0

``DIR`` holds in_<voxel>.json (labels and the local matrices). Every number in the report comes
from the summary.json written here.
"""

import argparse
import gzip
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(ICE / "scripts" / "common"))
import doc_tables as DT  # noqa: E402
from stats import wilson  # noqa: E402

GROUPS = [
    "truth",
    "off_0.1",
    "off_1.0",
    "off_3.0",
    "random",
    "pair_random",
    "sweep",
    "pair_differ_pass",
    "pair_differ_count",
]


def share(k: int, n: int) -> Dict[str, Any]:
    lo, hi = wilson(k, n)
    return dict(k=k, n=n, pct=100.0 * k / n if n else float("nan"), lo=100 * lo, hi=100 * hi)


def fmt_share(s: Dict[str, Any]) -> str:
    """'k (pct%)' with one decimal, or two when one decimal would round a partial share to
    100.0% or 0.0% (2,229 of 2,230 prints 99.96%, not 100.0%)."""
    pct = s["pct"]
    one = f"{pct:.1f}"
    if s["k"] != s["n"] and s["k"] != 0 and one in ("100.0", "0.0"):
        return f"{s['k']} ({pct:.2f}%)"
    return f"{s['k']} ({one}%)"


def rot_angle_deg(m: np.ndarray) -> float:
    return float(np.degrees(np.arccos(np.clip((np.trace(m) - 1.0) / 2.0, -1.0, 1.0))))


def arm_block(rows: List[Dict[str, Any]], arm: str) -> Dict[str, Any]:
    n = len(rows)
    A = [r[arm] for r in rows]
    no_pix = lambda a: not any(
        k in a["stages"] for k in ("frame_below_range", "frame_boundary", "hit", "pixel")
    )  # noqa: E731
    out = dict(
        n=n,
        conf_equal=share(sum(a["conf_equal"] for a in A), n),
        counts_equal=share(sum(a["counts_equal"] for a in A), n),
        quality_close=share(sum(a["quality_close"] for a in A), n),
        quality_bits_equal=share(sum(a["quality_bits_equal"] for a in A), n),
        same_frame_and_pixels=share(sum(no_pix(a) for a in A), n),
        no_peak_changes_counts=share(sum(not a["affecting"] for a in A), n),
        peaks_compared=int(sum(a["n_matched"] for a in A)),
        peaks_omega_bits_equal=(
            int(sum(a["n_omega_bits_equal"] for a in A)) if arm.startswith("EXF") else None
        ),
        peaks_unmatched_python=int(sum(a["py_only"] for a in A)),
        peaks_cpp_only_with_pixels=int(sum(a["cpp_only_live"] for a in A)),
        peaks_cpp_only_no_pixel=int(sum(a["cpp_only_dead"] for a in A)),
    )
    st: Counter = Counter()
    af: Counter = Counter()
    for a in A:
        st.update(a["stages"])
        af.update(a["affecting"])
    out["peak_stages"] = dict(st)
    out["peaks_changing_counts_by_cause"] = dict(af)
    out["orientations_with_cause"] = {
        k: sum(1 for a in A if a["affecting"].get(k)) for k in sorted(af)
    }
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    d, out = Path(a.dir), Path(a.out)
    res = json.load(gzip.open(out / "results.json.gz"))
    rows: List[Dict[str, Any]] = res["rows"]
    meta: Dict[Any, Dict[str, Any]] = {}
    for v in sorted({r["voxel"] for r in rows}):
        for m in json.load(open(d / f"in_{v}.json")):
            meta[(v, m["id"])] = m
    S: Dict[str, Any] = dict(n_orientations=len(rows), voxel_check=res["voxel_check"])
    S["by_label_n"] = dict(Counter(r["label"] for r in rows))

    # ---- exact paths: all orientations, by arm and by group ----
    S["exact"] = {arm: arm_block(rows, arm) for arm in ("EXS", "EXS_edge", "EXF", "EXF_edge")}
    S["exact_by_group"] = {}
    for g in GROUPS:
        gr = [r for r in rows if r["label"] == g]
        if gr:
            S["exact_by_group"][g] = {
                arm: arm_block(gr, arm) for arm in ("EXS", "EXS_edge", "EXF", "EXF_edge")
            }
    # omega bits and absolute differences (complete-cost arm, all matched peaks)
    om = np.array([x for r in rows for x in r["EXF_edge"]["omega_abs_list"]])
    S["omega"] = dict(
        n=int(len(om)),
        bits_equal_pct=100.0 * float(np.mean(om == 0.0)),
        abs_median_rad=float(np.median(om)),
        abs_p99_rad=float(np.percentile(om, 99)),
        abs_max_rad=float(om.max()),
        abs_max_deg=float(np.degrees(om.max())),
    )
    # product
    pe = [r["product_bit_equal"] for r in rows if "product_bit_equal" in r]
    S["product"] = dict(
        pairs=len(pe),
        bit_equal=int(sum(pe)),
        max_abs=max(r["product_max_abs_diff"] for r in rows if "product_max_abs_diff" in r),
    )
    # cost values (re-cost 1 - confidence) difference, python default vs edge against C++
    cost_diff = np.array([abs(r["cost_cpp_exf"] - r["cost_py_exf"]) for r in rows])
    S["recost_default_abs_diff"] = dict(
        eq=int((cost_diff == 0).sum()),
        lt1e6=int((cost_diff < 1e-6).sum()),
        lt1e3=int((cost_diff < 1e-3).sum()),
        n=len(rows),
        max=float(cost_diff.max()),
    )
    ce = np.array(
        [
            abs((1.0 - r["EXF_edge"]["cpp"]["conf"]) - (1.0 - r["EXF_edge"]["py"]["conf"]))
            for r in rows
        ]
    )
    S["recost_edge_abs_diff"] = dict(
        eq=int((ce == 0).sum()),
        lt1e6=int((ce < 1e-6).sum()),
        lt1e3=int((ce < 1e-3).sum()),
        n=len(rows),
        max=float(ce.max()),
    )
    qd = np.array(
        [abs(r["EXF_edge"]["cpp"]["quality"] - r["EXF_edge"]["py"]["quality"]) for r in rows]
    )
    S["quality_edge_abs_diff"] = dict(
        eq=int((qd == 0).sum()),
        lt1e6=int((qd < 1e-6).sum()),
        lt1e3=int((qd < 1e-3).sum()),
        n=len(rows),
        max=float(qd.max()),
    )
    # pixel mismatches against the double-precision position
    pxs: Dict[str, Any] = {}
    for arm in ("EXS_edge", "EXF_edge"):
        xs = [x for r in rows for pk in r[arm]["pixel_peaks"] for x in pk.get("pix64", [])]
        dist = np.array([x["dist"] for x in xs]) if xs else np.array([0.0])
        pxs[arm] = dict(
            coordinates=len(xs),
            peaks=int(sum(len(r[arm]["pixel_peaks"]) for r in rows)),
            dist_to_integer_max=float(dist.max()),
            dist_to_integer_median=float(np.median(dist)),
            cpp_matches_double=int(sum(x["cpp_right"] for x in xs)),
            python_matches_double=int(sum(x["py_right"] for x in xs)),
            neither=int(sum((not x["cpp_right"]) and (not x["py_right"]) for x in xs)),
        )
        ul = np.array([u for r in rows for u in r[arm]["pixel_omega_ulp"]])
        pxs[arm]["omega_identical_bits"] = int((ul == 0).sum()) if len(ul) else 0
        pxs[arm]["omega_differs"] = int((ul > 0).sum()) if len(ul) else 0
    S["pixel_boundary"] = pxs
    # degenerate-lit
    S["pixel_less_triangle"] = dict(
        peaks=int(sum(r["EXF"]["stages"].get("counts", 0) for r in rows)),
        with_lit_no_overlap_detector=int(sum(r["EXF"]["counts_lit_no_ov"] for r in rows)),
        cpp_pixon_zero=int(sum(r["EXF"]["counts_cpp_pixon_zero"] for r in rows)),
        orientations=int(sum(1 for r in rows if r["EXF"]["stages"].get("counts"))),
    )
    # frame
    S["frame_below_range"] = dict(
        peaks=int(sum(r["EXF"]["stages"].get("frame_below_range", 0) for r in rows)),
        peaks_matched_total=int(sum(r["EXF"]["n_matched"] for r in rows)),
        orientations_screen=int(
            sum(1 for r in rows if r["EXS"]["stages"].get("frame_below_range"))
        ),
        orientations_screen_cost_changed=int(
            sum(1 for r in rows if r["EXS"]["affecting"].get("frame_below_range"))
        ),
        orientations_full_cost_changed=int(
            sum(1 for r in rows if r["EXF"]["affecting"].get("frame_below_range"))
        ),
    )
    # pass flips from the frame bug (screen): python default pass vs C++
    S["screen_pass"] = dict(
        cpp_pass=int(sum(r["EXS"]["cpp"]["pkov"] > 0 for r in rows)),
        py_default_pass=int(sum(r["EXS"]["py"]["pkov"] > 0 for r in rows)),
        py_edge_pass=int(sum(r["EXS_edge"]["py"]["pkov"] > 0 for r in rows)),
        default_differs_from_cpp=int(
            sum((r["EXS"]["cpp"]["pkov"] > 0) != (r["EXS"]["py"]["pkov"] > 0) for r in rows)
        ),
        edge_differs_from_cpp=int(
            sum(
                (r["EXS_edge"]["cpp"]["pkov"] > 0) != (r["EXS_edge"]["py"]["pkov"] > 0)
                for r in rows
            )
        ),
    )

    # ---- interpolated screen (pairs only) ----
    pairs = [r for r in rows if r["kind"] == 1]
    kind0 = [r for r in rows if r["kind"] == 0]
    S["interp_control_identity_local"] = dict(
        n=len(kind0),
        same_counts=int(sum(r["INT_vs_EXS_cpp"]["same_counts"] for r in kind0)),
        omega_err_max_deg=float(max(r["INT_vs_EXS_cpp"]["omega_err_deg_max"] for r in kind0)),
    )
    byj: Dict[int, List[Dict[str, Any]]] = defaultdict(list)
    ang: Dict[int, float] = {}
    for r in pairs:
        if (
            r["label"] != "pair_random"
        ):  # the unselected pairs only (the others were picked for a difference)
            continue
        m = meta[(r["voxel"], r["id"])]
        j = int(m["j"])
        byj[j].append(r)
        ang[j] = rot_angle_deg(np.asarray(m["local"], dtype=np.float64))
    interp_rows = []
    for j in sorted(byj):
        rr = byj[j]
        mx = np.array([r["INT_vs_EXS_cpp"]["omega_err_deg_max"] for r in rr])
        md = np.array([r["INT_vs_EXS_cpp"]["omega_err_deg_med"] for r in rr])
        interp_rows.append(
            dict(
                j=j,
                local_rotation_deg=ang[j],
                n=len(rr),
                same_pass_pct=100.0 * np.mean([r["INT_vs_EXS_cpp"]["same_pass"] for r in rr]),
                same_counts_pct=100.0 * np.mean([r["INT_vs_EXS_cpp"]["same_counts"] for r in rr]),
                omega_err_median_of_median_deg=float(np.median(md)),
                omega_err_median_of_max_deg=float(np.median(mx)),
                omega_err_worst_deg=float(mx.max()),
                pct_with_exact_peak_missing=100.0
                * np.mean([r["INT_vs_EXS_cpp"]["exact_live_missing_in_int"] > 0 for r in rr]),
                pct_with_extra_peak=100.0
                * np.mean([r["INT_vs_EXS_cpp"]["int_live_not_in_exact"] > 0 for r in rr]),
            )
        )
    S["interp_by_j"] = interp_rows
    byt: Dict[float, List[Dict[str, Any]]] = defaultdict(list)
    for r in rows:
        if r["label"] == "sweep":
            byt[float(meta[(r["voxel"], r["id"])]["theta"])].append(r)
    sweep_rows = []
    for th in sorted(byt):
        rr = byt[th]
        mx = np.array([r["INT_vs_EXS_cpp"]["omega_err_deg_max"] for r in rr])
        n = len(rr)
        k_pass = sum(r["INT_vs_EXS_cpp"]["same_pass"] for r in rr)
        k_cnt = sum(r["INT_vs_EXS_cpp"]["same_counts"] for r in rr)
        lo, hi = wilson(k_cnt, n)
        sweep_rows.append(
            dict(
                theta_deg=th,
                n=n,
                same_pass_pct=100.0 * k_pass / n,
                same_counts_pct=100.0 * k_cnt / n,
                same_counts_lo=100 * lo,
                same_counts_hi=100 * hi,
                omega_err_median_of_max_deg=float(np.median(mx)),
                omega_err_p90_of_max_deg=float(np.percentile(mx, 90)),
                omega_err_worst_deg=float(mx.max()),
                pct_with_exact_peak_missing=100.0
                * np.mean([r["INT_vs_EXS_cpp"]["exact_live_missing_in_int"] > 0 for r in rr]),
                pct_with_extra_peak=100.0
                * np.mean([r["INT_vs_EXS_cpp"]["int_live_not_in_exact"] > 0 for r in rr]),
            )
        )
    S["interp_by_theta"] = sweep_rows
    allmx = np.array([r["INT_vs_EXS_cpp"]["omega_err_deg_max"] for r in pairs])
    S["interp_all_pairs"] = dict(
        n=len(pairs),
        same_pass=share(sum(r["INT_vs_EXS_cpp"]["same_pass"] for r in pairs), len(pairs)),
        same_counts=share(sum(r["INT_vs_EXS_cpp"]["same_counts"] for r in pairs), len(pairs)),
        omega_err_max_deg_median=float(np.median(allmx)),
        omega_err_max_deg_worst=float(allmx.max()),
    )

    # ---- pass patterns for the pairs ----
    def pattern(r: Dict[str, Any]) -> str:
        bits = (
            r["INT"]["pkov"] > 0,  # C++ interpolated (production)
            r["EXS"]["cpp"]["pkov"] > 0,  # C++ exact
            r["EXS"]["py"]["pkov"] > 0,  # Python exact, default
            r["EXS_edge"]["py"]["pkov"] > 0,  # Python exact, edge rules
        )
        return "".join("P" if b else "-" for b in bits)

    pat: Dict[str, Counter] = {}
    for g in ("pair_random", "pair_differ_pass", "pair_differ_count"):
        pat[g] = Counter(pattern(r) for r in pairs if r["label"] == g)
    S["pass_patterns"] = {g: dict(c) for g, c in pat.items()}
    # stored C++ trace agreement: the production C++ pass of the driver vs the stored trace
    dp = [r for r in pairs if r["label"] == "pair_differ_pass"]
    S["stored_trace_check"] = dict(
        n=len(dp),
        driver_cpp_pass_matches_stored=int(
            sum((r["INT"]["pkov"] > 0) == (r["meta_cpp"] > 0) for r in dp)
        ),
        driver_cpp_overlap_equals_stored=int(sum(r["INT"]["pkov"] == r["meta_cpp"] for r in dp)),
        stored_python_pass_matches_py_default=int(
            sum((r["EXS"]["py"]["pkov"] > 0) == (r["meta_py"] > 0) for r in dp)
        ),
    )
    # attribution of the pass difference between C++ (production, interpolated) and Python
    # (default exact)
    attr = Counter()
    for r in dp:
        ci, ce_, pd, pe_ = (
            r["INT"]["pkov"] > 0,
            r["EXS"]["cpp"]["pkov"] > 0,
            r["EXS"]["py"]["pkov"] > 0,
            r["EXS_edge"]["py"]["pkov"] > 0,
        )
        if ci == pd:
            attr["agree_now"] += 1
        elif ce_ == pd and ci != ce_:
            attr["interpolation (C++ interpolated differs from C++ exact; exact paths agree)"] += 1
        elif pd != pe_ and pe_ == ci:
            attr[
                "frame lookup below the range (Python default differs from edge rules;"
                " edge rules = C++)"
            ] += 1
        elif ci == ce_ and ce_ != pe_:
            attr["exact-path float32 boundary (C++ exact differs from Python edge exact)"] += 1
        else:
            attr["combination / other"] += 1
    S["differ_pass_attribution"] = dict(attr)
    # same for all pairs
    attr2 = Counter()
    for r in pairs:
        ci, ce_, pd, pe_ = (
            r["INT"]["pkov"] > 0,
            r["EXS"]["cpp"]["pkov"] > 0,
            r["EXS"]["py"]["pkov"] > 0,
            r["EXS_edge"]["py"]["pkov"] > 0,
        )
        attr2["C++ interpolated pass != C++ exact pass"] += int(ci != ce_)
        attr2["C++ exact pass != Python default exact pass"] += int(ce_ != pd)
        attr2["C++ exact pass != Python edge-rule exact pass"] += int(ce_ != pe_)
        attr2["C++ interpolated pass != Python default exact pass"] += int(ci != pd)
        attr2["C++ interpolated pass != Python edge-rule exact pass"] += int(ci != pe_)
    S["pairs_pass_disagreements"] = dict(attr2, n=len(pairs))

    (out / "summary.json").write_text(json.dumps(S, indent=1, default=float))

    # ---- tables ----
    T: Dict[str, str] = {}
    labels = {
        "EXS": "screen, Python default",
        "EXS_edge": "screen, Python with CppEdgeRules",
        "EXF": "complete cost, Python default",
        "EXF_edge": "complete cost, Python with CppEdgeRules",
    }
    trows = []
    for arm, lab in labels.items():
        b = S["exact"][arm]
        trows.append(
            dict(
                arm=lab,
                n=b["n"],
                cost=fmt_share(b["conf_equal"]),
                counts=fmt_share(b["counts_equal"]),
                frames=fmt_share(b["same_frame_and_pixels"]),
                qual=fmt_share(b["quality_close"]) if arm.startswith("EXF") else "n/a",
            )
        )
    T["parity_exact"] = DT.markdown_table(
        trows,
        ["arm", "n", "cost", "counts", "frames", "qual"],
        labels=None,
    )
    T["interp_by_theta"] = DT.markdown_table(
        sweep_rows,
        [
            "theta_deg",
            "n",
            "same_pass_pct",
            "same_counts_pct",
            "same_counts_lo",
            "same_counts_hi",
            "omega_err_median_of_max_deg",
            "omega_err_p90_of_max_deg",
            "omega_err_worst_deg",
            "pct_with_exact_peak_missing",
            "pct_with_extra_peak",
        ],
        formats={
            "theta_deg": ".2f",
            "same_pass_pct": ".1f",
            "same_counts_pct": ".1f",
            "same_counts_lo": ".1f",
            "same_counts_hi": ".1f",
            "omega_err_median_of_max_deg": ".4f",
            "omega_err_p90_of_max_deg": ".4f",
            "omega_err_worst_deg": ".3f",
            "pct_with_exact_peak_missing": ".1f",
            "pct_with_extra_peak": ".1f",
        },
    )
    prow = []
    for arm, lab in labels.items():
        b = S["exact"][arm]
        row = dict(arm=lab)
        for k in ("frame_below_range", "frame_boundary", "pixel", "counts"):
            row[k] = b["peaks_changing_counts_by_cause"].get(k, 0)
        row["omega_only"] = b["peak_stages"].get("omega_benign", 0)
        row["peaks"] = b["peaks_compared"]
        prow.append(row)
    T["parity_causes"] = DT.markdown_table(
        prow,
        ["arm", "peaks", "frame_below_range", "frame_boundary", "pixel", "counts", "omega_only"],
    )
    grow = []
    for g, blk in S["exact_by_group"].items():
        row = dict(group=g, n=blk["EXS"]["n"])
        for arm in ("EXS", "EXS_edge", "EXF", "EXF_edge"):
            row[arm] = blk[arm]["conf_equal"]["k"]
        grow.append(row)
    T["parity_by_group"] = DT.markdown_table(
        grow, ["group", "n", "EXS", "EXS_edge", "EXF", "EXF_edge"]
    )
    T["interp_by_j"] = DT.markdown_table(
        interp_rows,
        [
            "j",
            "local_rotation_deg",
            "n",
            "same_pass_pct",
            "same_counts_pct",
            "omega_err_median_of_max_deg",
            "omega_err_worst_deg",
            "pct_with_exact_peak_missing",
            "pct_with_extra_peak",
        ],
        formats={
            "local_rotation_deg": ".3f",
            "same_pass_pct": ".1f",
            "same_counts_pct": ".1f",
            "omega_err_median_of_max_deg": ".4f",
            "omega_err_worst_deg": ".4f",
            "pct_with_exact_peak_missing": ".1f",
            "pct_with_extra_peak": ".1f",
        },
    )
    prow2 = []
    for g, c in S["pass_patterns"].items():
        for pt, k in sorted(c.items(), key=lambda kv: -kv[1]):
            prow2.append(dict(group=g, pattern=pt, n=k))
    T["pass_patterns"] = DT.markdown_table(prow2, ["group", "pattern", "n"])
    (out / "tables.md").write_text("")
    DT.write_tables(out / "tables.md", T)
    print("wrote summary.json and tables.md")


if __name__ == "__main__":
    main()
