#!/usr/bin/env python3
"""Q8: the E0 per-voxel images of the three variants at the truth (cost, rejected real peaks).

The E0 counterpart of q8_sanity.py (which reads the Phase D full-sample images): for each of the
200 E0 voxels and each image variant (clean, all = old realistic, realistic_nodrop) it evaluates
the local hard cost (pixel radius 0, Q-max 8) at the TRUE orientation under four overlap rules
(old = lenient, cpp, rayR2, rayA) and two min_sin_eta settings (0.3 = E0 and the Q3 rows, 0 = the
ray-overlap rows of step 4). At the truth every detector-0 overlap is a real peak, so
rejected = detector-0 overlaps (case det0_bright of the ray-rule arms) - peaks kept. Cost = 1 -
quality.

  uv run python scripts/phase_d/q8_e0_truth.py run --workers 10
  uv run python scripts/phase_d/q8_e0_truth.py summary

Writes benchmarks/q8_nodrop/e0_truth_summary.json and e0_truth_tables.md. Not a timing run.
"""

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(ICE / "scripts" / "findoptimal_robustness"))
sys.path.insert(0, str(ICE / "scripts" / "common"))
import common as C  # noqa: E402
import doc_tables  # noqa: E402
import stats as S  # noqa: E402

CACHE = HERE / "cache" / "q8_e0_truth"
OUT = ICE / "benchmarks" / "q8_nodrop"
VARIANTS = ("clean", "all", "realistic_nodrop")
VLABEL = {"clean": "clean", "all": "realistic (10% drop)", "realistic_nodrop": "realistic, no drop"}
ARMS = ("old", "cpp", "rayR2", "rayA")
LABEL = {"old": "lenient", "cpp": "C++ rule", "rayR2": "ray 2px", "rayA": "ray geometric"}


def pretty(key: str) -> str:
    """Table header of a row key: arm names spelled out, `<arm>_ci` -> `<label> interval`."""
    if key.endswith("_ci") and key[:-3] in LABEL:
        return LABEL[key[:-3]] + " 95% interval"
    return LABEL.get(key, key)


def pretty_rows(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    return [{pretty(k): v for k, v in r.items()} for r in rows]


MSINS = (0.3, 0.0)
N_CASES = 8
DET0_BRIGHT = 5 + 7  # column of the case count det0_bright in a result row


def task(item: Tuple[int, int, str, str]) -> str:
    vidx, vpos, variant, path = item
    from icenine.ray_overlap import RayRule
    from optimizer_sweep import voxel_context

    W = C.get_worker()
    keys = np.load(C.CACHE_DIR / "images" / f"v{vidx}_{variant}.npz")["keys"]
    C.attach(keys)
    vctx = voxel_context(W.ctx, vidx)
    rec = W.rec
    sp = rec.setup.search_params
    reg = float(rec.setup.exp_setup.config_file.detector_dist_deviation)
    R = np.asarray(vctx.R_true, dtype=np.float32)
    res: Dict[str, Any] = dict(voxel=vidx, variant=variant)
    try:
        for msin in MSINS:
            rec.min_sin_eta = msin
            for arm in ARMS:
                sp.cpp_overlap_rule = arm != "old"
                ray = None
                if arm == "rayR2":
                    ray = RayRule(arm="R", radius_px=2, registration_mm=reg)
                elif arm == "rayA":
                    ray = RayRule(arm="A", registration_mm=reg)
                info = rec._make_local_cost_fn(ray).evaluate(R, vctx.vertices, vctx.voxel.phase)
                res[f"{msin}|{arm}"] = [
                    info.peak_overlap,
                    info.peak_on_detector,
                    info.pixel_overlap,
                    info.pixel_on_detector,
                    info.quality,
                ] + list(info.ray_cases if info.ray_cases is not None else [0] * N_CASES)
    finally:
        sp.cpp_overlap_rule = False
        rec.min_sin_eta = 0.3
    Path(path).write_text(json.dumps(res, default=float) + "\n")
    return f"voxel {vidx} {variant}"


def voxels() -> List[int]:
    info = dict(np.load(C.OUT_DIR / "voxels.npz"))
    return [int(v) for v in info["voxel_indices"]][: C.N_VOXELS]


def run(workers: int) -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    items = [
        (v, vpos, var, str(CACHE / f"v{v}_{var}.json"))
        for vpos, v in enumerate(voxels())
        for var in VARIANTS
        if (C.CACHE_DIR / "images" / f"v{v}_{var}.npz").exists()
    ]
    todo = [it for it in items if not Path(it[3]).exists()]
    print(f"{len(items)} cases, {len(todo)} to run", flush=True)
    t0 = time.time()
    C.run_pool(task, todo, workers, "q8_e0_truth")
    print(f"wall {time.time() - t0:.0f}s (not a timing claim)")


def summary() -> None:
    vox = voxels()
    recs: Dict[str, Dict[int, Dict[str, Any]]] = {v: {} for v in VARIANTS}
    for var in VARIANTS:
        for v in vox:
            p = CACHE / f"v{v}_{var}.json"
            if p.exists():
                recs[var][v] = json.loads(p.read_text())
    common = [v for v in vox if all(v in recs[var] for var in VARIANTS)]
    summ: Dict[str, Any] = {"n_voxels": len(common)}
    rej_rows: List[Dict[str, Any]] = []
    cost_rows: List[Dict[str, Any]] = []
    for msin in MSINS:
        for var in VARIANTS:
            rs = [recs[var][v] for v in common]
            ref = int(sum(r[f"{msin}|rayR2"][DET0_BRIGHT] for r in rs))
            assert ref == int(sum(r[f"{msin}|rayA"][DET0_BRIGHT] for r in rs))
            row: Dict[str, Any] = {"msin": msin, "images": VLABEL[var], "det0_overlaps": ref}
            crow: Dict[str, Any] = {"msin": msin, "images": VLABEL[var]}
            for arm in ARMS:
                kept = int(sum(r[f"{msin}|{arm}"][0] for r in rs))
                rej = ref - kept
                lo, hi = S.wilson(max(rej, 0), ref)
                q = np.array([r[f"{msin}|{arm}"][4] for r in rs], dtype=float)
                summ[f"{var}|{msin}|{arm}"] = dict(
                    kept=kept, det0=ref, rejected=rej, lo=lo, hi=hi,
                    cost_median=float(np.median(1 - q)), cost_mean=float((1 - q).mean()),
                )  # fmt: skip
                if arm != "old":
                    row[arm] = f"{rej} of {ref} ({100.0 * rej / ref:.1f}%)"
                    row[arm + "_ci"] = f"{100 * lo:.2f} to {100 * hi:.2f}%"
                crow[arm] = float(np.median(1 - q))
            rej_rows.append(row)
            cost_rows.append(crow)
    # paired per voxel: cost at the truth against clean, C++ rule, min_sin_eta 0.3
    for var in ("all", "realistic_nodrop"):
        d = np.array(
            [
                (1 - recs[var][v]["0.3|cpp"][4]) - (1 - recs["clean"][v]["0.3|cpp"][4])
                for v in common
            ]
        )
        summ[f"{var}|cost_minus_clean|0.3|cpp"] = dict(
            median=float(np.median(d)), n=len(d), n_higher=int((d > 1e-12).sum()),
            n_equal=int((np.abs(d) <= 1e-12).sum()), n_lower=int((d < -1e-12).sum()),
        )  # fmt: skip
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "e0_truth_summary.json").write_text(json.dumps(summ, indent=1, default=float) + "\n")
    cols = ["msin", "images", "det0_overlaps"] + [
        pretty(c) for a in ARMS if a != "old" for c in (a, a + "_ci")
    ]
    t_rej = doc_tables.markdown_table(pretty_rows(rej_rows), cols)
    arm_cols = [LABEL[a] for a in ARMS]
    t_cost = doc_tables.markdown_table(
        pretty_rows(cost_rows),
        ["msin", "images"] + arm_cols,
        formats={"msin": ".1f", **{a: ".4f" for a in arm_cols}},
    )
    doc_tables.write_tables(
        OUT / "e0_truth_tables.md", {"q8_e0_truth_cost": t_cost, "q8_e0_truth_rejected": t_rej}
    )
    print(t_cost, "\n")
    print(t_rej)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("cmd", choices=["run", "summary"])
    ap.add_argument("--workers", type=int, default=10)
    a = ap.parse_args()
    run(a.workers) if a.cmd == "run" else summary()
