"""Task E measurements 2 and 3, summary.

Measurement 2 (ranking only, on the baseline runs ``<kind>_<rule>_base``): per level 0-2, the
post-quick-MC list of the strict search is re-ranked by each arm. The truth basin at a level is the
set of candidates within ``diameter / 1.5`` (the next level's grid radius: a kept candidate there
brings the truth into the next level's grid) or within 1 deg. For every case where the basin is
present after the quick MC: is it kept by the quarter cut under the strict ranking and under the
arm's ranking (fixed = only the arm keeps it, broken = only the strict ranking keeps it), and the
rank of the best basin candidate.

Measurement 3 (end to end, ``<kind>_<rule>_<arm tag>`` against ``<kind>_<rule>_base``): wrong
counts (error above 1 deg) with Wilson intervals, the paired split by voxel (fixed / broken / both
wrong), evaluations of the strict cost and of the smooth cost. Every broken voxel is traced
(owner rule D-20261014-N13):
  traced to the arm   at a level the arm ranks, the basin was present after the quick MC, the
                      arm's key cut it and the strict cost would have kept it (same list)
  traced to the draw  the arm's key kept the basin at every level it ranks, or the strict cost
                      would have cut it too, or the basin was already absent after a quick MC
                      (lost in the screen or the random quick MC, where the arm does not act)
  untraced            records missing
Fixed voxels are traced the same way (the arm kept a basin the strict cost would have cut).

  uv run python scripts/eval_reduction/smooth_rank_summary.py rank
  uv run python scripts/eval_reduction/smooth_rank_summary.py e2e --arms e123a,e12a
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICE / "scripts" / "common"))
import doc_tables as DT  # noqa: E402
import stats as ST  # noqa: E402
from smooth_rank_probe import CACHE, RANK_ARMS, RANK_LEVELS, SIX  # noqa: E402

OUT = ICE / "benchmarks" / "eval_reduction" / "E"
WRONG_DEG = 1.0
SETS = [("e0", "old"), ("e0", "cpp"), ("full", "old"), ("full", "cpp")]
SET_LABEL = {
    ("e0", "old"): "E0 200, lenient rule",
    ("e0", "cpp"): "E0 200, C++ rule",
    ("full", "old"): "six full-sample voxels, lenient rule",
    ("full", "cpp"): "six full-sample voxels, C++ rule",
    ("full30", "old"): "30 full-sample voxels, lenient rule (supplementary arms)",
}


def load_runs(d: Path) -> Dict[int, Dict[str, Any]]:
    out = {}
    for f in sorted(d.glob("v*.npz")):
        z = np.load(f)
        out[int(z["voxel"])] = {k: z[k] for k in z.files}
    return out


def basin_radius(z: Dict[str, Any], lv: int, kind: str) -> float:
    if kind == "1deg":
        return 1.0
    return float(np.degrees(z[f"L{lv}_diameter"])) / 1.5


def kept(dist: np.ndarray, order: np.ndarray, n_keep: int, B: float) -> Tuple[bool, int]:
    """(basin kept by the cut, rank of the best basin candidate) under ``order``."""
    inb = dist[order] <= B
    if not inb.any():
        return False, -1
    r = int(np.argmax(inb))
    return r < n_keep, r


# --------------------------------------------------------------------------------------------
# measurement 2
# --------------------------------------------------------------------------------------------
def rank_summary(
    kind: str, rule: str, bdef: str, tags: Tuple[str, ...] = ("base",), arm_list: Any = None
) -> Optional[Dict[str, Any]]:
    runs: Dict[int, Dict[str, Any]] = {}
    for t in tags:
        runs.update(load_runs(CACHE / f"{kind}_{rule}_{t}"))
    if not runs:
        return None
    arms = ["strict"] + list(arm_list if arm_list is not None else RANK_ARMS)
    res: Dict[str, Any] = {"n_voxels": len(runs), "levels": {}}
    for lv in range(RANK_LEVELS):
        L: Dict[str, Any] = {a: dict(kept=0, fixed=0, broken=0, ranks=[], rel=[]) for a in arms}
        present = 0
        for v, z in runs.items():
            if f"L{lv}_dist" not in z:
                continue
            dist = z[f"L{lv}_dist"].astype(np.float64)
            nk = int(z[f"L{lv}_n_keep"])
            B = basin_radius(z, lv, bdef)
            if not (dist <= B).any():
                continue
            present += 1
            ks, rs = kept(dist, np.arange(len(dist)), nk, B)
            for a in arms:
                if a == "strict":
                    ka, ra = ks, rs
                else:
                    order = np.argsort(z[f"L{lv}_key_{a}"], kind="stable")
                    ka, ra = kept(dist, order, nk, B)
                L[a]["kept"] += int(ka)
                L[a]["fixed"] += int(ka and not ks)
                L[a]["broken"] += int(ks and not ka)
                L[a]["ranks"].append(ra)
                L[a]["rel"].append(ra / len(dist))
        res["levels"][str(lv)] = dict(present=present)
        for a in arms:
            x = L[a]
            res["levels"][str(lv)][a] = dict(
                kept=x["kept"],
                fixed=x["fixed"],
                broken=x["broken"],
                kept_wilson95=list(ST.wilson(x["kept"], present)) if present else [None, None],
                median_best_rank=float(np.median(x["ranks"])) if x["ranks"] else None,
                median_best_rank_fraction=float(np.median(x["rel"])) if x["rel"] else None,
            )
    # pooled over levels 0-2
    pool: Dict[str, Any] = {}
    tot_present = sum(res["levels"][str(lv)]["present"] for lv in range(RANK_LEVELS))
    for a in arms:
        k = sum(res["levels"][str(lv)][a]["kept"] for lv in range(RANK_LEVELS))
        pool[a] = dict(
            kept=k,
            fixed=sum(res["levels"][str(lv)][a]["fixed"] for lv in range(RANK_LEVELS)),
            broken=sum(res["levels"][str(lv)][a]["broken"] for lv in range(RANK_LEVELS)),
            kept_wilson95=list(ST.wilson(k, tot_present)) if tot_present else [None, None],
        )
    res["pooled"] = dict(present=tot_present, **pool)
    return res


SUPP_ARMS = ["E4a", "E12a", "E123a", "E34a", "E1234a", "E34t25"]


def rank_tables(S: Dict[str, Any]) -> Dict[str, str]:
    T: Dict[str, str] = {}
    for bdef in ("next", "1deg"):
        for kind, rule in SETS + [("full30", "old")]:
            R = S.get(f"{kind}_{rule}_{bdef}")
            if R is None:
                continue
            rows = []
            for a in ["strict"] + (SUPP_ARMS if kind == "full30" else list(RANK_ARMS)):
                row: Dict[str, Any] = {"ranking key": a if a != "strict" else "strict (baseline)"}
                for lv in range(RANK_LEVELS):
                    x = R["levels"][str(lv)]
                    y = x[a]
                    cell = f"{y['kept']} of {x['present']}"
                    if a != "strict":
                        cell += f" (+{y['fixed']} / -{y['broken']})"
                    row[f"level {lv}: basin kept"] = cell
                    row[f"level {lv}: best rank (median)"] = y["median_best_rank"]
                p = R["pooled"]
                lo, hi = p[a]["kept_wilson95"]
                row["levels 0-2: kept (95% Wilson)"] = (
                    f"{p[a]['kept']} of {p['present']} ({100 * lo:.1f}-{100 * hi:.1f}%)"
                )
                row["levels 0-2: fixed / broken"] = (
                    "" if a == "strict" else f"{p[a]['fixed']} / {p[a]['broken']}"
                )
                rows.append(row)
            T[f"evalE_rank_{kind}_{rule}_{bdef}"] = DT.markdown_table(
                rows,
                list(rows[0].keys()),
                {k: ".0f" for k in rows[0] if "best rank" in k},
            )
    return T


# --------------------------------------------------------------------------------------------
# measurement 3
# --------------------------------------------------------------------------------------------
def trace(z: Dict[str, Any], levels: int) -> Tuple[str, str]:
    """(class, detail) of the loss of the truth basin in an arm run (see the module docstring)."""
    for lv in range(levels):
        if f"L{lv}_dist" not in z:
            return "untraced", f"no record of level {lv}"
        dist = z[f"L{lv}_dist"].astype(np.float64)
        cost = z[f"L{lv}_cost"].astype(np.float64)
        nk = int(z[f"L{lv}_n_keep"])
        B = basin_radius(z, lv, "next")
        if not (dist <= B).any():
            return "draw", f"basin absent after the quick MC of level {lv}"
        arm_kept = bool((dist[:nk] <= B).any())
        strict_order = np.argsort(cost, kind="stable")
        strict_kept = bool((dist[strict_order[:nk]] <= B).any())
        if not arm_kept and strict_kept:
            return "arm", f"level {lv}: the arm cut the basin, the strict cost would have kept it"
        if not arm_kept:
            return "draw", f"level {lv}: the strict cost would have cut the basin too"
    return "draw", "basin kept through the smoothed levels; lost later (level 3 or FindOptimal)"


def fixed_trace(z: Dict[str, Any], levels: int) -> str:
    for lv in range(levels):
        if f"L{lv}_dist" not in z:
            return "untraced"
        dist = z[f"L{lv}_dist"].astype(np.float64)
        cost = z[f"L{lv}_cost"].astype(np.float64)
        nk = int(z[f"L{lv}_n_keep"])
        B = basin_radius(z, lv, "next")
        if not (dist <= B).any():
            return "draw"
        arm_kept = bool((dist[:nk] <= B).any())
        strict_kept = bool((dist[np.argsort(cost, kind="stable")[:nk]] <= B).any())
        if arm_kept and not strict_kept:
            return "arm"
    return "draw"


def e2e(kind: str, rule: str, tag: str, levels: int, subset: Optional[List[int]] = None):
    base = load_runs(CACHE / f"{kind}_{rule}_base")
    arm = load_runs(CACHE / f"{kind}_{rule}_{tag}")
    if not arm or not base:
        return None
    vox = sorted(set(base) & set(arm))
    if subset is not None:
        vox = [v for v in vox if v in subset]
    if not vox:
        return None
    spec = (CACHE / f"{kind}_{rule}_{tag}" / "spec.txt").read_text().strip()
    eb = np.array([float(base[v]["err_deg"]) for v in vox])
    ea = np.array([float(arm[v]["err_deg"]) for v in vox])
    wb, wa = eb > WRONG_DEG, ea > WRONG_DEG
    only_a, only_b = ST.paired_discordant(wa, wb)  # broken (arm wrong only), fixed
    n = len(vox)
    broken = [v for v, x, y in zip(vox, wa, wb) if x and not y]
    fixed = [v for v, x, y in zip(vox, wa, wb) if y and not x]
    tr = {v: trace(arm[v], levels) for v in broken}
    ft = {v: fixed_trace(arm[v], levels) for v in fixed}
    both_right = ~wa & ~wb
    out = dict(
        spec=spec,
        n=n,
        wrong_base=int(wb.sum()),
        wrong_arm=int(wa.sum()),
        wrong_base_wilson95=list(ST.wilson(int(wb.sum()), n)),
        wrong_arm_wilson95=list(ST.wilson(int(wa.sum()), n)),
        broken=int(only_a),
        fixed=int(only_b),
        both_wrong=int((wa & wb).sum()),
        broken_voxels=broken,
        fixed_voxels=fixed,
        broken_traced_arm=sum(1 for c, _ in tr.values() if c == "arm"),
        broken_traced_draw=sum(1 for c, _ in tr.values() if c == "draw"),
        broken_untraced=sum(1 for c, _ in tr.values() if c == "untraced"),
        broken_trace={str(v): list(t) for v, t in tr.items()},
        fixed_traced_arm=sum(1 for c in ft.values() if c == "arm"),
        fixed_traced_draw=sum(1 for c in ft.values() if c == "draw"),
        mean_evals_base=float(np.mean([int(base[v]["evals"]) for v in vox])),
        mean_evals_arm=float(np.mean([int(arm[v]["evals"]) for v in vox])),
        mean_smooth_evals_arm=float(np.mean([int(arm[v]["smooth_evals"]) for v in vox])),
        median_err_both_right_base=float(np.median(eb[both_right])) if both_right.any() else None,
        median_err_both_right_arm=float(np.median(ea[both_right])) if both_right.any() else None,
        n_both_right=int(both_right.sum()),
        strict_evals_change=float(
            np.mean([int(arm[v]["evals"]) for v in vox])
            / np.mean([int(base[v]["evals"]) for v in vox])
            - 1.0
        ),
        # where the wrong voxels end (deg from the truth, sorted), and how many end within 1 deg
        # of 60 deg (a near-degenerate solution)
        err_wrong_base=sorted(float(e) for e in eb[wb]),
        err_wrong_arm=sorted(float(e) for e in ea[wa]),
        n_wrong_base_near60=int((wb & (np.abs(eb - 60.0) < 1.0)).sum()),
        n_wrong_arm_near60=int((wa & (np.abs(ea - 60.0) < 1.0)).sum()),
        n_fixed_not_near60=int((wb & ~wa & (np.abs(eb - 60.0) >= 1.0)).sum()),
    )
    if kind == "e0" and rule == "old":
        out.update(base_vs_task_c(base, vox))
    return out


def base_vs_task_c(base: Dict[int, Dict[str, Any]], vox: List[int]) -> Dict[str, Any]:
    """The E0 lenient baseline (``SmoothCostRank`` off on the Task E branch) against Task C's
    baseline cache (``C/cache/e0/base``, run before Task E): voxels with equal error and equal
    evaluation count."""
    root = OUT.parent / "C" / "cache" / "e0" / "base"
    n_cmp = n_same = 0
    for v in vox:
        f = root / f"v{v}.json"
        if not f.exists():
            continue
        c = json.loads(f.read_text())
        n_cmp += 1
        n_same += int(
            float(c["err_deg"]) == float(base[v]["err_deg"])
            and int(c["evals"]) == int(base[v]["evals"])
        )
    return dict(task_c_base_compared=n_cmp, task_c_base_identical=n_same)


def e2e_tables(E: Dict[str, Any]) -> Dict[str, str]:
    rows = []
    for key, x in E.items():
        if x is None:
            continue
        kind, rule, tag, sub = key.split("|")
        lo_b, hi_b = x["wrong_base_wilson95"]
        lo_a, hi_a = x["wrong_arm_wilson95"]
        rows.append(
            {
                "voxel set": sub,
                "rule": "lenient" if rule == "old" else "C++",
                "arm (SmoothCostRank)": x["spec"],
                "n": x["n"],
                "wrong, baseline": f"{x['wrong_base']} ({100 * lo_b:.1f}-{100 * hi_b:.1f}%)",
                "wrong, arm": f"{x['wrong_arm']} ({100 * lo_a:.1f}-{100 * hi_a:.1f}%)",
                "fixed": x["fixed"],
                "broken": x["broken"],
                "both wrong": x["both_wrong"],
                "broken: traced to arm / draw / untraced": " / ".join(
                    str(x[k])
                    for k in ("broken_traced_arm", "broken_traced_draw", "broken_untraced")
                ),
                "fixed: traced to arm / draw": " / ".join(
                    str(x[k]) for k in ("fixed_traced_arm", "fixed_traced_draw")
                ),
                "median error of voxels right in both (deg), baseline / arm": (
                    f"{x['median_err_both_right_base']:.4f} / {x['median_err_both_right_arm']:.4f}"
                ),
                "strict evaluations per voxel, baseline": x["mean_evals_base"],
                "strict evaluations per voxel, arm": x["mean_evals_arm"],
                "smooth evaluations per voxel": x["mean_smooth_evals_arm"],
            }
        )
    if not rows:
        return {}
    fm = {k: ".0f" for k in rows[0] if "evaluations" in k}
    return {"evalE_e2e": DT.markdown_table(rows, list(rows[0].keys()), fm)}


def timing_tables() -> Dict[str, str]:
    """Measurement 4.

    * ``evalE_timing_ms``: ms per evaluation (``smooth_timing.py``, single-worker, preflight-gated).
    * ``evalE_timing_voxel_wall``: measured wall time of whole voxels, key off against on
      (``smooth_voxel_timing.py``; labelled single-worker or contended from its preflight).
    * ``evalE_timing_per_voxel``: a MODEL of the full-sample voxel (evaluations per voxel from the
      end-to-end runs x ms per evaluation at the level-0 widths, plus the build), as shipped (a
      ``SmoothCost`` build in every ``reconstruct_voxel``) and as a projection for a build cache per
      image set that is not built (Task G). Full-sample rows only: on E0 most counted evaluations
      are screen and re-cost evaluations, far cheaper than one local call, so the model does not
      hold there (Task E review); every E0 voxel has its own images, so a cache would not help."""
    d = OUT / "timing"
    T: Dict[str, str] = {}
    tim = {}
    for which in ("full", "e0"):
        f = d / f"timing_{which}.json"
        if f.exists():
            tim[which] = json.loads(f.read_text())
    if not tim:
        return T
    rows = []
    for name in tim[next(iter(tim))]["arms"]:
        row: Dict[str, Any] = {"cost": name}
        for which, t in tim.items():
            lab = "full-sample" if which == "full" else "E0"
            a = t["arms"][name]
            row[f"{lab}: ms per evaluation, batch of 256"] = a["ms_per_eval_batch256"]
            row[f"{lab}: ms per evaluation, one per call"] = a["ms_per_eval_single"]
            row[f"{lab}: x strict"] = a["ratio_batch256_to_strict"]
        rows.append(row)
    strict = {"cost": "strict cost (one call per orientation)"}
    for which, t in tim.items():
        lab = "full-sample" if which == "full" else "E0"
        strict[f"{lab}: ms per evaluation, batch of 256"] = t["strict_ms_per_eval"]
        strict[f"{lab}: ms per evaluation, one per call"] = t["strict_ms_per_eval"]
        strict[f"{lab}: x strict"] = 1.0
    rows.insert(0, strict)
    fm = {k: ".3f" for k in rows[0] if k != "cost"}
    T["evalE_timing_ms"] = DT.markdown_table(rows, list(rows[0].keys()), fm)
    wall_rows = []
    for which in ("e0", "full"):
        f = d / f"voxel_wall_{which}.json"
        if not f.exists():
            continue
        W = json.loads(f.read_text())
        by: Dict[int, Dict[str, List[float]]] = {}
        same: Dict[int, bool] = {}
        for r in W["runs"]:
            by.setdefault(r["voxel"], {}).setdefault(r["spec"], []).append(r["wall_s"])
            same[r["voxel"]] = same.get(r["voxel"], True) and bool(r["same_as_cache"])
        for v, x in by.items():
            off = float(np.median(x["off"]))
            row = {
                "voxel set": "E0" if which == "e0" else "full-sample",
                "voxel": v,
                "timing": f"{W['label']}, median of {len(x['off'])}",
                "key off (s)": off,
            }
            for spec in ("E3+E4,tau=0.25", "E1+E2+E3+E4"):
                on = float(np.median(x[spec]))
                row[f"{spec} (s)"] = on
                row[f"{spec}: extra time"] = f"{on / off - 1:+.0%}"
            row["same error as the end-to-end run"] = "yes" if same[v] else "no"
            wall_rows.append(row)
    if wall_rows:
        T["evalE_timing_voxel_wall"] = DT.markdown_table(
            wall_rows,
            list(wall_rows[0].keys()),
            {k: ".1f" for k in wall_rows[0] if k.endswith("(s)")},
        )
    e2e_path = OUT / "e2e" / "summary.json"
    if e2e_path.exists():
        E = json.loads(e2e_path.read_text())
        rows = []
        arm_of = {"E3+E4,tau=0.25": "E3+E4", "E1+E2+E3+E4": "E1+E2+E3+E4"}
        for key, x in E.items():
            kind, rule, tag, sub = key.split("|")
            if rule != "old" or kind != "full" or kind not in tim:
                continue
            t = tim[kind]
            ms = t["arms"][arm_of[x["spec"]]]["ms_per_eval_batch256"]
            smooth_s = x["mean_smooth_evals_arm"] * ms / 1000.0
            build_s = t["build_trees_s"] + t["build_rays_s"]
            strict_s = x["mean_evals_arm"] * t["strict_ms_per_eval"] / 1000.0
            base_s = x["mean_evals_base"] * t["strict_ms_per_eval"] / 1000.0
            rows.append(
                {
                    "voxel set": sub,
                    "arm": x["spec"],
                    "modelled strict cost time per voxel, baseline (s)": base_s,
                    "modelled strict cost time per voxel, arm (s)": strict_s,
                    "smooth evaluations per voxel": x["mean_smooth_evals_arm"],
                    "modelled smooth cost time per voxel (s)": smooth_s,
                    "SmoothCost build (s)": build_s,
                    "extra time as shipped (a build per voxel)": (
                        f"{(strict_s + smooth_s + build_s) / base_s - 1:+.1%}"
                    ),
                    "projection: one build per image set (cache not built)": (
                        f"{(strict_s + smooth_s) / base_s - 1:+.1%}"
                    ),
                }
            )
        if rows:
            T["evalE_timing_per_voxel"] = DT.markdown_table(
                rows,
                list(rows[0].keys()),
                {k: ".1f" for k in rows[0] if "(s)" in k} | {"smooth evaluations per voxel": ".0f"},
            )
    DT.write_tables(d / "tables.md", T)
    return T


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["rank", "e2e", "timing"])
    ap.add_argument("--arms", default="")
    ap.add_argument("--levels", type=int, default=RANK_LEVELS)
    a = ap.parse_args()
    if a.cmd == "rank":
        S: Dict[str, Any] = {}
        for kind, rule in SETS:
            for bdef in ("next", "1deg"):
                S[f"{kind}_{rule}_{bdef}"] = rank_summary(kind, rule, bdef)
        # supplementary: the 30-voxel lenient set (six with all arms + 24 with SUPP_ARMS)
        if (CACHE / "full_old_base30").exists():
            for bdef in ("next", "1deg"):
                S[f"full30_old_{bdef}"] = rank_summary(
                    "full", "old", bdef, ("base", "base30"), SUPP_ARMS
                )
        d = OUT / "rank"
        d.mkdir(parents=True, exist_ok=True)
        (d / "summary.json").write_text(json.dumps(S, indent=1, default=float) + "\n")
        T = rank_tables(S)
        DT.write_tables(d / "tables.md", T)
        for k, t in T.items():
            print(k)
            print(t)
            print()
        return
    if a.cmd == "timing":
        for k, t in timing_tables().items():
            print(k)
            print(t)
        return
    E: Dict[str, Any] = {}
    for tag in [t for t in a.arms.split(",") if t]:
        for kind, rule in SETS:
            if kind == "e0":
                E[f"{kind}|{rule}|{tag}|E0 200"] = e2e(kind, rule, tag, a.levels)
            else:
                E[f"{kind}|{rule}|{tag}|six full-sample"] = e2e(kind, rule, tag, a.levels, SIX)
                if len(load_runs(CACHE / f"{kind}_{rule}_{tag}")) > len(SIX):
                    E[f"{kind}|{rule}|{tag}|30 full-sample"] = e2e(kind, rule, tag, a.levels)
    E = {k: v for k, v in E.items() if v is not None}
    d = OUT / "e2e"
    d.mkdir(parents=True, exist_ok=True)
    (d / "summary.json").write_text(json.dumps(E, indent=1, default=float) + "\n")
    T = e2e_tables(E)
    DT.write_tables(d / "tables.md", T)
    for k, t in T.items():
        print(k)
        print(t)


if __name__ == "__main__":
    main()
