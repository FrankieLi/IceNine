"""Task A step 3: tables of the omega-interpolation and batched-screen results.

  uv run python scripts/eval_reduction/step3_summary.py --a benchmarks/eval_reduction/A

Reads step3/screen_v15901*.json (pass sets against the stored C++ trace, screen_check.py),
step3/screen_timing.json (microsecond split, screen_timing.py) and the ledger runs (per-voxel
wall time and per-level screen time of the arms with the switches off and on). Writes
step3/summary.json and step3/tables.md.
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(ICE / "scripts" / "common"))
import doc_tables as DT  # noqa: E402

SIX = [15901, 17034, 18558, 19889, 22867, 5242]
ARM_LABEL = {
    "unbatched": "exact, one orientation per call (today)",
    "batched": "exact, 32 FZ points per call (BatchedScreen)",
    "interp1": "interpolated, one FZ point per call (CppOmegaInterp)",
    "interp32": "interpolated, 32 FZ points per call (both)",
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", required=True)
    a = ap.parse_args()
    root = Path(a.a).resolve()
    st3 = root / "step3"
    S: Dict[str, Any] = {}
    T: Dict[str, str] = {}

    # ---- level-0 pass sets of voxel 15901 against the stored C++ trace ----
    rows = []
    for tag, edge in (("v15901", "off"), ("v15901_edge", "on")):
        p = st3 / f"screen_{tag}.json"
        if not p.exists():
            continue
        r = json.loads(p.read_text())
        for arm in ("unbatched", "batched", "interp1", "interp32"):
            x = r[arm]
            rows.append(
                dict(
                    arm=ARM_LABEL[arm],
                    edge_rules=edge,
                    pass_pairs=x["pass_pairs"],
                    cpp_pass_pairs=r["cpp_pass_pairs"],
                    fz_identical=x["fz_identical_pass_set"],
                    pairs_pass_differs=x["pairs_pass_differs"],
                    overlap_equal=x["overlap_equal_pairs"],
                    overlap_equal_pct=x["overlap_equal_pct"],
                    screen_s_contended=r["wall_s"][arm],
                )
            )
        S[f"screen_{tag}"] = r
    cols = [
        "arm",
        "edge_rules",
        "pass_pairs",
        "cpp_pass_pairs",
        "fz_identical",
        "pairs_pass_differs",
        "overlap_equal",
        "overlap_equal_pct",
    ]
    T["passsets_l0"] = DT.markdown_table(
        rows,
        cols,
        formats={
            "pass_pairs": ",d",
            "cpp_pass_pairs": ",d",
            "fz_identical": ",d",
            "pairs_pass_differs": ",d",
            "overlap_equal": ",d",
            "overlap_equal_pct": ".3f",
        },
    )
    # exact-vs-batched identity
    if (st3 / "screen_v15901_edge.json").exists():
        r = json.loads((st3 / "screen_v15901_edge.json").read_text())
        S["batched_vs_unbatched_overlap_differs"] = r["between_arms_overlap_differs"][
            "unbatched|batched"
        ]
        S["interp1_vs_interp32_overlap_differs"] = r["between_arms_overlap_differs"][
            "interp1|interp32"
        ]

    # ---- microsecond split (quiet, single worker) ----
    p = st3 / "screen_timing.json"
    if p.exists():
        r = json.loads(p.read_text())
        S["screen_timing"] = r
        rows = []
        base = r["unbatched"]["us_per_orientation"]
        for arm in ("unbatched", "batched", "interp1", "interp32"):
            x = r[arm]
            rows.append(
                dict(
                    arm=ARM_LABEL[arm],
                    us_total=x["us_per_orientation"],
                    us_solve=x["us_solve"],
                    us_overlap_c=x["us_stage_d"],
                    us_rest=x["us_other"],
                    speedup=base / x["us_per_orientation"],
                )
            )
        T["speed_screen"] = DT.markdown_table(
            rows,
            ["arm", "us_total", "us_solve", "us_overlap_c", "us_rest", "speedup"],
            formats={k: ".1f" for k in ("us_total", "us_solve", "us_overlap_c", "us_rest")}
            | {"speedup": ".1f"},
        )

    # ---- per-voxel: ledger runs ----
    def run(tag: str, v: int) -> Any:
        p = root / "ledger" / "runs" / f"{tag}_v{v}.json"
        return json.loads(p.read_text()) if p.exists() else None

    arms = [
        ("cpp_edge", "off"),
        ("cpp_edge_batched", "BatchedScreen"),
        ("cpp_edge_interp", "CppOmegaInterp"),
        ("cpp_edge_interp_batched", "both"),
    ]
    rows = []
    for tag, lab in arms:
        rr = [run(tag, v) for v in SIX]
        rr = [x for x in rr if x]
        if not rr:
            continue
        l0 = [x["levels"][0] for x in rr]
        scr_us = [1e6 * c["screen_eval_s"] / c["screen_evals"] for c in l0]
        wall = [x["wall_s"] for x in rr]
        l02 = [sum(c["discrete_s"] + c["quick_s"] for c in x["levels"][:3]) for x in rr]
        rows.append(
            dict(
                switches=lab,
                n=len(rr),
                l0_screen_us=float(np.mean(scr_us)),
                l0_screen_s=float(np.mean([c["screen_eval_s"] for c in l0])),
                levels012_s=float(np.mean(l02)),
                wall_s=float(np.mean(wall)),
                evals=float(np.mean([x["evals_total"] for x in rr])),
                contended=sum("contended" in x["label"] for x in rr),
            )
        )
    S["voxel_speed"] = rows
    T["speed_voxel"] = DT.markdown_table(
        rows,
        [
            "switches",
            "n",
            "l0_screen_us",
            "l0_screen_s",
            "levels012_s",
            "wall_s",
            "evals",
            "contended",
        ],
        formats={
            "l0_screen_us": ".1f",
            "l0_screen_s": ".1f",
            "levels012_s": ".1f",
            "wall_s": ".1f",
            "evals": ",.0f",
        },
    )
    # quiet repeat
    rows = []
    for tag in ("timing_cpp_edge", "timing_cpp_edge_interp_batched"):
        for v in (15901, 5242):
            x = run(tag, v)
            if x:
                l0 = x["levels"][0]
                rows.append(
                    dict(
                        arm=tag.replace("timing_", ""),
                        voxel=v,
                        label=x["label"],
                        wall_s=x["wall_s"],
                        evals=x["evals_total"],
                        l0_screen_us=1e6 * l0["screen_eval_s"] / l0["screen_evals"],
                        us_per_eval=x["mean_us_per_eval"],
                        load=x["preflight"]["loadavg"][0],
                    )
                )
    if rows:
        S["quiet_repeat"] = rows
        T["speed_voxel_quiet"] = DT.markdown_table(
            rows,
            ["arm", "voxel", "label", "wall_s", "evals", "l0_screen_us", "us_per_eval", "load"],
            formats={
                "voxel": "d",
                "wall_s": ".1f",
                "evals": ",d",
                "l0_screen_us": ".1f",
                "us_per_eval": ".0f",
                "load": ".2f",
            },
        )
    (st3 / "summary.json").write_text(json.dumps(S, indent=1, default=float) + "\n")
    DT.write_tables(st3 / "tables.md", T)
    print("wrote step3/summary.json and tables.md")


if __name__ == "__main__":
    main()
