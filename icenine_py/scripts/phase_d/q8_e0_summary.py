#!/usr/bin/env python3
"""Q8 step 4 summary: E0 (200 voxels, seed 0) on clean, old realistic (all) and realistic_nodrop.

Reads the per-voxel results of scripts/phase_d/run_q8_step4.sh (scripts/phase_d/cache/
<arm>[_ms0]_q8/v<voxel>_<variant>.npz, written by e0_cpprule.py) and the run_meta files, and writes
benchmarks/q8_nodrop/e0_summary.json, e0_per_voxel.npz (error, evaluations, cost and runtime of every
voxel, arm and variant) and e0_tables.md (doc_tables blocks).

  uv run python scripts/phase_d/q8_e0_summary.py

A voxel is wrong when its reduced misorientation to the truth exceeds 1 degree (C.WRONG_DEG).
Paired splits are counted per voxel: for an arm A against a reference B, "broken" = right under B
and wrong under A, "fixed" = wrong under B and right under A, "both wrong" = wrong under both.
Intervals are 95% Wilson intervals on the wrong fraction (voxels treated as independent). No
p-values. Wall times are those of one 10-worker pool per (arm, variant): contended multi-worker,
not a timing claim; the per-voxel runtime is the same contended time.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(ICE / "scripts" / "findoptimal_robustness"))
sys.path.insert(0, str(ICE / "scripts" / "common"))
import common as C  # noqa: E402
import doc_tables  # noqa: E402
import stats as S  # noqa: E402

CACHE = HERE / "cache"
OUT = ICE / "benchmarks" / "q8_nodrop"
META = OUT / "runs"
EARLIER_Q3 = ICE / "benchmarks" / "phase_d_seed_diag" / "level0" / "e0" / "summary_on.json"
EARLIER_RAY = ICE / "benchmarks" / "phase_d_seed_diag" / "ray_overlap" / "e2e_summary.json"

# arm id -> (e0_cpprule.py --arm, --msin, label)
ARMS: Dict[str, Tuple[str, float, str]] = {
    "off": ("off", 0.3, "lenient rule, min_sin_eta 0.3"),
    "on": ("on", 0.3, "C++ rule, min_sin_eta 0.3"),
    "off_ms0": ("off", 0.0, "lenient rule, min_sin_eta 0"),
    "on_ms0": ("on", 0.0, "C++ rule, min_sin_eta 0"),
    "rayR2_ms0": ("rayR2", 0.0, "ray rule 2-pixel window, min_sin_eta 0"),
    "rayA_ms0": ("rayA", 0.0, "ray rule geometric window, min_sin_eta 0"),
}
VARIANTS = ("clean", "all", "realistic_nodrop")
VLABEL = {"clean": "clean", "all": "realistic (10% drop)", "realistic_nodrop": "realistic, no drop"}
# rule comparisons within a min_sin_eta setting: (arm, reference arm)
RULE_PAIRS = [
    ("on", "off"),
    ("on_ms0", "off_ms0"),
    ("rayR2_ms0", "off_ms0"),
    ("rayA_ms0", "off_ms0"),
    ("rayR2_ms0", "on_ms0"),
    ("rayA_ms0", "on_ms0"),
    ("rayA_ms0", "rayR2_ms0"),
]
# stored earlier counts (other machine and code revision): arm id -> key of e2e_summary.json
EARLIER_RAY_KEY = {"off_ms0": "old", "on_ms0": "cpp", "rayR2_ms0": "R2", "rayA_ms0": "A"}


def voxels() -> List[int]:
    info = dict(np.load(C.OUT_DIR / "voxels.npz"))
    return [int(v) for v in info["voxel_indices"]][: C.N_VOXELS]


def load(arm_id: str, var: str, vox: List[int]) -> Optional[Dict[str, np.ndarray]]:
    root = CACHE / f"{arm_id}_q8"
    files = [root / f"v{v}_{var}.npz" for v in vox]
    if not all(f.exists() for f in files):
        return None
    z = [np.load(f) for f in files]
    Rt = np.stack([x["R_true"] for x in z])
    R = np.stack([x["R"] for x in z])
    err = C.err_deg(R, Rt)
    return dict(
        err=err,
        wrong=err > C.WRONG_DEG,
        evals_global=np.array([int(x["evals_global"]) for x in z]),
        evals_local=np.array([int(x["evals_local"]) for x in z]),
        runtime=np.array([float(x["runtime"]) for x in z]),
        cost=np.array([float(x["cost"]) for x in z]),
    )


def image_check(vox: List[int]) -> Dict[str, Any]:
    """Lit pixels of the cached E0 case images of the three variants (keys, sorted unique)."""
    out: Dict[str, Any] = dict(n_voxels=len(vox), lit=dict.fromkeys(VARIANTS, 0))
    all_not_in_nodrop = nodrop_not_in_all = clean_not_in_nodrop = 0
    for v in vox:
        k = {var: np.load(C.CACHE_DIR / "images" / f"v{v}_{var}.npz")["keys"] for var in VARIANTS}
        for var in VARIANTS:
            out["lit"][var] += len(k[var])
        all_not_in_nodrop += len(np.setdiff1d(k["all"], k["realistic_nodrop"]))
        nodrop_not_in_all += len(np.setdiff1d(k["realistic_nodrop"], k["all"]))
        clean_not_in_nodrop += len(np.setdiff1d(k["clean"], k["realistic_nodrop"]))
    out.update(
        old_realistic_pixels_not_in_nodrop=all_not_in_nodrop,
        nodrop_pixels_not_in_old_realistic=nodrop_not_in_all,
        clean_pixels_not_in_nodrop=clean_not_in_nodrop,
    )
    return out


def split(a: np.ndarray, ref: np.ndarray) -> Dict[str, int]:
    """Paired split of arm wrong flags `a` against reference wrong flags `ref`, per voxel."""
    broken, fixed = S.paired_discordant(a, ref)
    return dict(
        fixed=fixed,
        broken=broken,
        both_wrong=int((a & ref).sum()),
        both_right=int((~a & ~ref).sum()),
    )


def wall(arm_id: str, var: str) -> Optional[float]:
    arm, msin, _ = ARMS[arm_id]
    suffix = "_ms0" if msin == 0 else ""
    p = META / f"run_meta_run_{arm}{suffix}_q8_{var}.json"
    return float(json.loads(p.read_text())["wall_s"]) if p.exists() else None


def main() -> None:
    vox = voxels()
    data: Dict[str, Dict[str, Dict[str, np.ndarray]]] = {}
    for arm_id in ARMS:
        for var in VARIANTS:
            d = load(arm_id, var, vox)
            if d is not None:
                data.setdefault(arm_id, {})[var] = d
    earlier_q3 = json.loads(EARLIER_Q3.read_text()) if EARLIER_Q3.exists() else {}
    earlier_ray = json.loads(EARLIER_RAY.read_text()) if EARLIER_RAY.exists() else {}
    summ: Dict[str, Any] = dict(n_voxels=len(vox), wrong_deg=C.WRONG_DEG, arms={}, rules={})
    summ["images"] = image_check(vox)
    wrong_rows: List[Dict[str, Any]] = []
    var_rows: List[Dict[str, Any]] = []
    rule_rows: List[Dict[str, Any]] = []
    for arm_id, (arm, msin, label) in ARMS.items():
        for var in VARIANTS:
            d = data.get(arm_id, {}).get(var)
            if d is None:
                continue
            n, k = len(d["wrong"]), int(d["wrong"].sum())
            lo, hi = S.wilson(k, n)
            right = ~d["wrong"]
            ev = d["evals_global"] + d["evals_local"]
            earlier = ""
            if var in ("clean", "all"):
                if arm_id == "off" and earlier_q3:
                    earlier = f"{earlier_q3[var]['wrong_before']}"
                elif arm_id == "on" and earlier_q3:
                    earlier = f"{earlier_q3[var]['wrong_after']}"
                elif arm_id in EARLIER_RAY_KEY and earlier_ray:
                    key = f"{'clean' if var == 'clean' else 'realistic'}|{EARLIER_RAY_KEY[arm_id]}"
                    earlier = f"{earlier_ray[key]['wrong']}" if key in earlier_ray else ""
            rec = dict(
                n=n, wrong=k, lo=lo, hi=hi, mean_evals=float(ev.mean()),
                mean_evals_global=float(d["evals_global"].mean()),
                mean_evals_local=float(d["evals_local"].mean()),
                mean_runtime_s=float(d["runtime"].mean()), wall_s=wall(arm_id, var),
                median_err_right=float(np.median(d["err"][right])), n_right=int(right.sum()),
                mean_cost=float(d["cost"].mean()), earlier_wrong=earlier,
            )  # fmt: skip
            summ["arms"].setdefault(arm_id, {})[var] = rec
            wrong_rows.append(
                dict(
                    arm=label, images=VLABEL[var], wrong=f"{k} of {n}",
                    interval=f"{lo:.3f} to {hi:.3f}", evals=rec["mean_evals"],
                    runtime=rec["mean_runtime_s"], wall_min=(rec["wall_s"] or float("nan")) / 60,
                    med_err=rec["median_err_right"], stored_earlier=earlier or "n/a",
                )
            )  # fmt: skip
        dd = data.get(arm_id, {})
        if all(v in dd for v in VARIANTS):
            nd, old, cl = (dd[v]["wrong"] for v in ("realistic_nodrop", "all", "clean"))
            e_nd = dd["realistic_nodrop"]["evals_global"] + dd["realistic_nodrop"]["evals_local"]
            e_old = dd["all"]["evals_global"] + dd["all"]["evals_local"]
            s_clean, s_old, s_oc = split(nd, cl), split(nd, old), split(old, cl)
            summ["arms"][arm_id]["paired_variants"] = dict(
                nodrop_vs_clean=s_clean, nodrop_vs_old_realistic=s_old,
                old_realistic_vs_clean=s_oc,
                evals_nodrop_minus_old_mean=float((e_nd - e_old).mean()),
            )  # fmt: skip
            for name, s_ in (
                ("no drop vs clean", s_clean),
                ("no drop vs old realistic", s_old),
                ("old realistic vs clean", s_oc),
            ):
                var_rows.append(dict(arm=label, comparison=name, **s_))
    for arm_id, ref_id in RULE_PAIRS:
        for var in VARIANTS:
            a, r = data.get(arm_id, {}).get(var), data.get(ref_id, {}).get(var)
            if a is None or r is None:
                continue
            s_ = split(a["wrong"], r["wrong"])
            summ["rules"][f"{var}|{arm_id}|vs|{ref_id}"] = s_
            rule_rows.append(
                dict(
                    images=VLABEL[var], rule=ARMS[arm_id][2], reference=ARMS[ref_id][2],
                    wrong_rule=int(a["wrong"].sum()), wrong_ref=int(r["wrong"].sum()), **s_,
                )
            )  # fmt: skip
    summ["total_wall_min"] = float(
        sum(
            (rec["wall_s"] or 0.0)
            for by_var in summ["arms"].values()
            for rec in by_var.values()
            if isinstance(rec, dict) and "wall_s" in rec
        )
        / 60
    )
    OUT.mkdir(parents=True, exist_ok=True)
    # the per-voxel numbers behind every paired split (the result caches are gitignored)
    per_voxel: Dict[str, np.ndarray] = {"voxel": np.array(vox)}
    for arm_id, by_var in data.items():
        for var, d in by_var.items():
            for key in ("err", "evals_global", "evals_local", "cost", "runtime"):
                per_voxel[f"{arm_id}|{var}|{key}"] = d[key]
    np.savez_compressed(OUT / "e0_per_voxel.npz", **per_voxel)
    (OUT / "e0_summary.json").write_text(json.dumps(summ, indent=1, default=float) + "\n")
    tables = {
        "q8_e0_wrong": doc_tables.markdown_table(
            wrong_rows,
            [
                "arm",
                "images",
                "wrong",
                "interval",
                "evals",
                "runtime",
                "wall_min",
                "med_err",
                "stored_earlier",
            ],
            formats={"evals": ".0f", "runtime": ".1f", "wall_min": ".1f", "med_err": ".4f"},
        ),
        "q8_e0_paired_variants": doc_tables.markdown_table(
            var_rows, ["arm", "comparison", "fixed", "broken", "both_wrong", "both_right"]
        ),
        "q8_e0_paired_rules": doc_tables.markdown_table(
            rule_rows,
            [
                "images",
                "rule",
                "reference",
                "wrong_rule",
                "wrong_ref",
                "fixed",
                "broken",
                "both_wrong",
            ],
        ),
    }
    doc_tables.write_tables(OUT / "e0_tables.md", tables)
    for t in tables.values():
        print(t, "\n")


if __name__ == "__main__":
    main()
