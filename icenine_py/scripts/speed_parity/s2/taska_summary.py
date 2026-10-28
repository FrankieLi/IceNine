"""S2 Task A check summary: summary.json and tables.md (doc_tables) from ``taska_check.py``.

Inputs (``--dir``, default benchmarks/s2_taska): inputs.json, run.json.gz, outlier.json; and the
stored Task A results (``benchmarks/eval_reduction/A/step0/results.json.gz``) for the C++ and the
Task A Python aggregates. Every number in the tables comes from these files.

C++ agreement is counted per orientation, as MIGRATION_HISTORY "A-exact" counts it: the peak counts
(peak overlap, peak on detector) equal to C++'s; for the complete cost also the pixel counts.

  uv run python scripts/speed_parity/s2/taska_summary.py [--dir benchmarks/s2_taska]
"""

import argparse
import gzip
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Tuple

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
sys.path.insert(0, str(ICE / "scripts" / "common"))

import doc_tables as DT  # noqa: E402
import stats  # noqa: E402

STORED = ICE / "benchmarks" / "eval_reduction" / "A" / "step0" / "results.json.gz"
ARMS = ("EXS", "EXF", "EXS_edge", "EXF_edge")
FUSED_ARMS = ("EXS_edge", "EXF_edge")
ARM_LABEL = {
    "EXS": "screen",
    "EXF": "complete cost",
    "EXS_edge": "screen, CppEdgeRules",
    "EXF_edge": "complete cost, CppEdgeRules",
}
VARIANTS = ("portable", "accelerate")
PEAK = ("pkov", "pkon")
PIXEL = ("pixov", "pixon")
ALL5 = ("pixov", "pixon", "pkov", "pkon", "quality")
QUALITY_TOL = 1e-6  # the C++ ground-truth tolerance of the port (CLAUDE.md)


def eq(a: Dict[str, Any], b: Dict[str, Any], keys: Tuple[str, ...]) -> bool:
    return all(a[k] == b[k] for k in keys)


def frac(k: int, n: int) -> Dict[str, Any]:
    lo, hi = stats.wilson(k, n)
    return dict(k=k, n=n, frac=k / n if n else None, lo=lo, hi=hi)


def paired(a: List[bool], b: List[bool]) -> Dict[str, int]:
    """a = torch matches C++, b = fused matches C++ (per orientation)."""
    return dict(
        both=sum(x and y for x, y in zip(a, b)),
        torch_only=sum(x and not y for x, y in zip(a, b)),
        fused_only=sum(y and not x for x, y in zip(a, b)),
        neither=sum(not x and not y for x, y in zip(a, b)),
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dir", default=str(ICE / "benchmarks" / "s2_taska"))
    a = ap.parse_args()
    d = Path(a.dir).resolve()
    inputs = json.loads((d / "inputs.json").read_text())
    run = json.loads(gzip.open(d / "run.json.gz").read())
    outlier = json.loads((d / "outlier.json").read_text())
    stored = {
        (int(r["voxel"]), int(r["id"])): r for r in json.loads(gzip.open(STORED).read())["rows"]
    }
    rows = run["rows"]
    n = len(rows)
    S: Dict[str, Any] = dict(meta=dict(run=run["meta"], inputs=inputs["meta"]), n=n)

    # ---- regeneration ------------------------------------------------------------------------
    regen: Dict[str, Any] = {}
    for arm in ARMS:
        regen[arm] = sum(
            eq(r[arm]["torch"], stored[(r["voxel"], r["id"])][arm]["py"], ALL5) for r in rows
        )
    S["regeneration"] = dict(
        records=inputs["n_records"],
        stored=inputs["n_stored"],
        labels_and_indices_differ=len(inputs["labels_and_indices_differ"]),
        records_by_kind=dict(Counter(str(r["kind"]) for r in rows)),
        single_records_identity_local=inputs.get("single_records_identity_local"),
        pair_records_two_or_more_elements_differing=inputs.get(
            "pair_records_two_or_more_elements_differing"
        ),
        orders=inputs["orders"],
        selected_order=inputs["selected_order"],
        torch_equals_stored_python=regen,
        mirror_equal={arm: sum(r[arm]["mirror_equal"] for r in rows) for arm in FUSED_ARMS},
    )
    t_inputs = [
        dict(
            check="records: label, kind, FZ and grid indices equal",
            n=inputs["n_stored"],
            equal=inputs["n_records"] - len(inputs["labels_and_indices_differ"]),
        )
    ]
    for o, x in inputs["orders"].items():
        t_inputs.append(
            dict(
                check=f"pair product, order {o}: bit-equal flag and largest difference",
                n=x["n_pairs"],
                equal=x["both_match"],
            )
        )
    for arm in ARMS:
        t_inputs.append(
            dict(
                check=f"torch at {run['meta']['head']} = stored Python, {ARM_LABEL[arm]}",
                n=n,
                equal=regen[arm],
            )
        )

    # ---- against C++ -------------------------------------------------------------------------
    cpp_rows: List[Dict[str, Any]] = []
    pair_rows: List[Dict[str, Any]] = []
    S["cpp"] = {}
    for arm in FUSED_ARMS:
        cpp = [stored[(r["voxel"], r["id"])][arm]["cpp"] for r in rows]
        stp = [stored[(r["voxel"], r["id"])][arm]["py"] for r in rows]
        paths: Dict[str, List[Dict[str, Any]]] = {
            "torch (Task A, stored)": stp,
            f"torch ({run['meta']['head']})": [r[arm]["torch"] for r in rows],
        }
        for v in VARIANTS:
            paths[f"fused {v}"] = [r[arm][v]["agg"] for r in rows]
        S["cpp"][arm] = {}
        # the stored C++ screen values have pixel overlap and quality 0 (counted here; the C++ side
        # was not traced), so only its peak counts compare; pixel counts and quality are compared on
        # the complete cost only
        S["cpp"][arm]["cpp_pixov_and_quality_zero"] = sum(
            c["pixov"] == 0 and c["quality"] == 0.0 for c in cpp
        )
        peak_match: Dict[str, List[bool]] = {}
        pix_match: Dict[str, List[bool]] = {}
        for p, aggs in paths.items():
            pk = [eq(x, c, PEAK) for x, c in zip(aggs, cpp)]
            px = [eq(x, c, PEAK + PIXEL) for x, c in zip(aggs, cpp)]
            qu = [abs(x["quality"] - c["quality"]) <= QUALITY_TOL for x, c in zip(aggs, cpp)]
            peak_match[p], pix_match[p] = pk, px
            ent: Dict[str, Any] = dict(peak_counts=frac(sum(pk), n))
            if arm == "EXF_edge":
                ent["peak_and_pixel_counts"] = frac(sum(px), n)
                ent["quality_within_tol"] = sum(qu)
            S["cpp"][arm][p] = ent
            cpp_rows.append(
                dict(
                    arm=ARM_LABEL[arm],
                    path=p,
                    n=n,
                    peak_counts_equal=sum(pk),
                    pixel_counts_equal=sum(px) if arm == "EXF_edge" else None,
                    quality_within_1e6=sum(qu) if arm == "EXF_edge" else None,
                )
            )
        tkey = f"torch ({run['meta']['head']})"
        for v in VARIANTS:
            for what, mm in (("peak counts", peak_match), ("peak and pixel counts", pix_match)):
                if what == "peak and pixel counts" and arm != "EXF_edge":
                    continue
                sp = paired(mm[tkey], mm[f"fused {v}"])
                S["cpp"][arm][f"paired_{v}_{what.replace(' ', '_')}"] = sp
                pair_rows.append(dict(arm=ARM_LABEL[arm], fused=v, counts=what, **sp))
        # the orientations where some path disagrees with C++ on the peak counts
        S["cpp"][arm]["peak_mismatch_records"] = [
            dict(
                voxel=r["voxel"],
                id=r["id"],
                label=r["label"],
                cpp=[c[k] for k in PEAK],
                **{p: [paths[p][i][k] for k in PEAK] for p in paths},
            )
            for i, (r, c) in enumerate(zip(rows, cpp))
            if not all(peak_match[p][i] for p in paths)
        ]

    # ---- fused against torch on this set -----------------------------------------------------
    vt_rows: List[Dict[str, Any]] = []
    S["fused_vs_torch"] = {}
    for arm in FUSED_ARMS:
        for v in VARIANTS:
            diff = [r for r in rows if r[arm][v]["differs_from_torch"]]
            causes: Counter = Counter()
            fields: Counter = Counter()
            for r in diff:
                causes.update(r[arm][v]["verdict"]["causes"])
                for k in ALL5:
                    fields[k] += int(r[arm][v]["agg"][k] != r[arm]["torch"][k])
            worst = max((r[arm][v]["verdict"]["worst_near_px"] for r in diff), default=None)
            n_rule = sum(r[arm][v]["verdict"]["within_rule"] for r in diff)
            flipped = [x for r in diff for x in r[arm][v]["trace"]]
            side = Counter(x.get("ref_side") for x in flipped if "ref_side" in x)
            sil = [r[arm][v]["silent_verdict"] for r in rows if r[arm][v].get("silent_rows")]
            silent = len(sil)
            silent_rule = sum(x["within_rule"] for x in sil)
            silent_causes = sorted({c for x in sil for c in x["causes"]})
            silent_worst = max((x["worst_near_px"] for x in sil), default=None)
            ent = dict(
                n=n,
                differing=frac(len(diff), n),
                within_rule=n_rule,
                causes=dict(causes),
                fields=dict(fields),
                worst_near_px=worst,
                flipped_rows=len(flipped),
                # the flipped rows are pooled over these orientations (not independent rows)
                flipped_rows_orientations=sum(bool(r[arm][v]["trace"]) for r in diff),
                ref_side=dict(side),
                ref_side_torch=frac(side.get("torch", 0), sum(side.values())),
                silent=silent,
                silent_within_rule=silent_rule,
                silent_causes=silent_causes,
                silent_worst_near_px=silent_worst,
                differing_records=[
                    dict(
                        voxel=r["voxel"], id=r["id"], label=r["label"], verdict=r[arm][v]["verdict"]
                    )
                    for r in diff
                ],
            )
            S["fused_vs_torch"][f"{arm}/{v}"] = ent
            vt_rows.append(
                dict(
                    arm=ARM_LABEL[arm],
                    fused=v,
                    n=n,
                    differing=len(diff),
                    within_rule=n_rule,
                    causes=", ".join(sorted(causes)) or "-",
                    worst_px=worst,
                    pixon=fields["pixon"],
                    pixov=fields["pixov"],
                    pkov=fields["pkov"],
                    pkon=fields["pkon"],
                    quality=fields["quality"],
                    f64_side_torch=side.get("torch", 0),
                    f64_side_fused=side.get("fused", 0),
                    silent=silent,
                    silent_within_rule=silent_rule,
                    silent_causes=", ".join(silent_causes) or "-",
                )
            )

    # ---- the Step 0 outlier ------------------------------------------------------------------
    # one table row per (mode, orientation, Bragg row); the detector coordinate (in pixels, from
    # the coordinate 0 of the -1 rule) only where the -1 rule decides the pixel
    groups: Dict[Tuple[str, int, int], List[Dict[str, Any]]] = {}
    for mode, x in outlier["modes"].items():
        for e in x["entries"]:
            groups.setdefault((mode, e["orientation"], e["row"]), []).append(e)
    out_rows: List[Dict[str, Any]] = []
    for (mode, o, r), es in groups.items():
        neg = [e for e in es if e["torch_stage_d_int"] == -1 or e["ref_stage_d_int"] == -1]
        out_rows.append(
            dict(
                mode=mode,
                orientation=o,
                row=r,
                coords=len(es),
                torch_max_err_px=max(e["torch_err"] for e in es),
                fused_max_err_px=max(e["fused_err"] for e in es),
                torch_eq_fused=all(e["torch_pix"] == e["fused_pix"] for e in es),
                int_flips_vs_f64=sum(e["torch_stage_d_int"] != e["ref_stage_d_int"] for e in es),
                asin_arg=es[0]["asin_arg_torch"],
                omega_err_rad=max(abs(e["omega_torch"] - e["omega_ref"]) for e in es),
                coord0_f32_px=neg[0]["coord_torch_over_pixel"] if neg else None,
                coord0_f64_px=neg[0]["coord_f64_over_pixel"] if neg else None,
            )
        )
    S["outlier"] = dict(
        meta=outlier["meta"],
        n_stage_d_rows={m: x["n_stage_d_rows"] for m, x in outlier["modes"].items()},
        n_entries={m: len(x["entries"]) for m, x in outlier["modes"].items()},
        groups=out_rows,
        # coordinates beyond min_px, and their integer pixel flips against float64, per orientation
        per_orientation=[
            dict(
                mode=mode,
                orientation=o,
                rows=sorted(
                    g["row"] for g in out_rows if (g["mode"], g["orientation"]) == (mode, o)
                ),
                coords=sum(
                    g["coords"] for g in out_rows if (g["mode"], g["orientation"]) == (mode, o)
                ),
                int_flips_vs_f64=sum(
                    g["int_flips_vs_f64"]
                    for g in out_rows
                    if (g["mode"], g["orientation"]) == (mode, o)
                ),
            )
            for mode, o in sorted({(g["mode"], g["orientation"]) for g in out_rows})
        ],
    )

    tables = {
        "s2_taska_inputs": DT.markdown_table(t_inputs, ["check", "n", "equal"]),
        "s2_taska_cpp": DT.markdown_table(
            cpp_rows,
            ["arm", "path", "n", "peak_counts_equal", "pixel_counts_equal", "quality_within_1e6"],
        ),
        "s2_taska_paired": DT.markdown_table(
            pair_rows, ["arm", "fused", "counts", "both", "torch_only", "fused_only", "neither"]
        ),
        "s2_taska_vs_torch": DT.markdown_table(
            vt_rows,
            [
                "arm",
                "fused",
                "n",
                "differing",
                "within_rule",
                "causes",
                "worst_px",
                "pixon",
                "pixov",
                "pkov",
                "pkon",
                "quality",
                "f64_side_torch",
                "f64_side_fused",
                "silent",
                "silent_within_rule",
                "silent_causes",
            ],
            formats=dict(worst_px=".4f"),
        ),
        "s2_taska_outlier": DT.markdown_table(
            out_rows,
            [
                "mode",
                "orientation",
                "row",
                "coords",
                "torch_max_err_px",
                "fused_max_err_px",
                "torch_eq_fused",
                "int_flips_vs_f64",
                "asin_arg",
                "omega_err_rad",
                "coord0_f32_px",
                "coord0_f64_px",
            ],
            formats=dict(
                torch_max_err_px=".4f",
                fused_max_err_px=".4f",
                asin_arg=".5f",
                omega_err_rad=".2e",
                coord0_f32_px=".2e",
                coord0_f64_px=".2e",
            ),
        ),
    }
    (d / "summary.json").write_text(json.dumps(S, indent=1, default=float) + "\n")
    DT.write_tables(d / "tables.md", tables)
    for k, t in tables.items():
        print(f"## {k}\n{t}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
