"""Task G2 summary: summary.json and tables.md under benchmarks/eval_reduction/G/g2/.

Parts (each skipped when its inputs are missing):
  ledger  step 2, six interior full-sample voxels (clean, Q-max 8, seed 1): evaluations by stage
          of the key-off search (``g2_base``, against Task A's ``base``) and of the integrated
          no-start configuration (``g2_integrated``: SmoothCostRank E3+E4 tau 0.25; OutwardSearch
          does not act on a no-start search)
  e0      step 2, the E0 200 voxels (clean, seed 0, mc arm, min_sin_eta 0): the integrated
          configuration paired against the stored lenient baseline (Task E's ``e0_old_base``),
          breaks traced with Task E's rule (D-20261014-N13), and equality with Task E's own run of
          the same key (``e0_old_e34t25``)
  pilot   the 150-voxel sub-region pilot (key-off check, OutwardSearch only, integrated)
  bfs     step 3, the 2,000-voxel pilot region: the stored pilot
          (phase_d_pilot_reseed_v2/runs/mc_clean) against the integrated configuration and
          against OutwardSearch alone

  uv run python scripts/eval_reduction/g2_summary.py
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICE / "scripts" / "common"))
import doc_tables as DT  # noqa: E402
import g2_bfs  # noqa: E402
from stats import wilson  # noqa: E402

G = ICE / "benchmarks" / "eval_reduction" / "G"
OUT = G / "g2"
CACHE = G / "cache" / "g2"
SIX = [15901, 17034, 18558, 19889, 22867, 5242]
STAGES = ["screen", "recost", "quick", "find", "variance", "final"]
STORED = ICE / "benchmarks" / "phase_d_pilot_reseed_v2" / "runs" / "mc_clean"
# Task F's single-worker weights on the full sample (F/timing/full.json, full_f5.json), ms
F_TIMING = ICE / "benchmarks" / "eval_reduction" / "F" / "timing"


# --------------------------------------------------------------------------------------------
def ledger() -> Optional[Dict[str, Any]]:
    d = OUT / "ledger" / "runs"
    a = ICE / "benchmarks" / "eval_reduction" / "A" / "ledger" / "runs"
    if not (d / "g2_integrated_v5242.json").exists():
        return None
    rows: List[Dict[str, Any]] = []
    stage_sum: Dict[str, Dict[str, int]] = {}
    out: Dict[str, Any] = dict(voxels={})
    for v in SIX:
        ta = json.loads((a / f"base_v{v}.json").read_text())
        b = json.loads((d / f"g2_base_v{v}.json").read_text())
        g = json.loads((d / f"g2_integrated_v{v}.json").read_text())
        same = ta["evals_total"] == b["evals_total"] and ta["err_deg"] == b["err_deg"]
        for tag, r in (("task_a_base", ta), ("g2_integrated", g)):
            s = stage_sum.setdefault(tag, {k: 0 for k in STAGES})
            for k in STAGES:
                s[k] += int(r["totals"][k])
        row = dict(
            voxel=v,
            evals_task_a=ta["evals_total"],
            evals_key_off=b["evals_total"],
            key_off_equals_task_a=bool(same),
            evals_integrated=g["evals_total"],
            change=g["evals_total"] / ta["evals_total"] - 1.0,
            smooth=int(g.get("smooth_evals", 0)),
            err_task_a=ta["err_deg"],
            err_integrated=g["err_deg"],
            wall_key_off=b["wall_s"],
            wall_integrated=g["wall_s"],
        )
        rows.append(row)
        out["voxels"][str(v)] = row
    out["stages"] = stage_sum
    tot_a = sum(r["evals_task_a"] for r in rows)
    tot_g = sum(r["evals_integrated"] for r in rows)
    out["totals"] = dict(
        task_a=tot_a,
        integrated=tot_g,
        change=tot_g / tot_a - 1.0,
        smooth=sum(r["smooth"] for r in rows),
        wrong_task_a=sum(1 for r in rows if r["err_task_a"] > 1.0),
        wrong_integrated=sum(1 for r in rows if r["err_integrated"] > 1.0),
        all_key_off_equal_task_a=all(r["key_off_equals_task_a"] for r in rows),
        median_wall_key_off_s=float(np.median([r["wall_key_off"] for r in rows])),
        median_wall_integrated_s=float(np.median([r["wall_integrated"] for r in rows])),
    )
    T = {}
    T["evalG2_ledger_voxels"] = DT.markdown_table(
        rows,
        ["voxel", "evals_task_a", "evals_key_off", "evals_integrated", "change", "smooth",
         "err_task_a", "err_integrated", "wall_key_off", "wall_integrated"],
        formats=dict(change="+.1%", err_task_a=".4f", err_integrated=".4f", wall_key_off=".0f",
                     wall_integrated=".0f"),
    )  # fmt: skip
    srows = []
    for k in STAGES:
        x, y = stage_sum["task_a_base"][k], stage_sum["g2_integrated"][k]
        srows.append(dict(stage=k, task_a=x, integrated=y, change=(y / x - 1.0) if x else None))
    srows.append(dict(stage="all", task_a=tot_a, integrated=tot_g, change=tot_g / tot_a - 1.0))
    srows.append(
        dict(
            stage="smooth (E3+E4, outside the strict count)",
            task_a=0,
            integrated=out["totals"]["smooth"],
        )
    )
    T["evalG2_ledger_stages"] = DT.markdown_table(
        srows, ["stage", "task_a", "integrated", "change"], formats=dict(change="+.1%")
    )
    out["tables"] = T
    return out


# --------------------------------------------------------------------------------------------
def e0() -> Optional[Dict[str, Any]]:
    import smooth_rank_summary as SRS

    root = ICE / "benchmarks" / "eval_reduction" / "E" / "cache" / "rank"
    if not (root / "e0_old_g2").exists() or len(list((root / "e0_old_g2").glob("v*.npz"))) < 200:
        return None
    res = SRS.e2e("e0", "old", "g2", 3)
    # equality with Task E's run of the same key
    n_same = n_cmp = 0
    diff: List[int] = []
    for p in sorted((root / "e0_old_g2").glob("v*.npz")):
        q = root / "e0_old_e34t25" / p.name
        if not q.exists():
            continue
        a, b = np.load(p), np.load(q)
        n_cmp += 1
        same = (
            float(a["err_deg"]) == float(b["err_deg"])
            and int(a["evals"]) == int(b["evals"])
            and int(a["smooth_evals"]) == int(b["smooth_evals"])
            and np.array_equal(a["R_final"], b["R_final"])
        )
        n_same += int(same)
        if not same:
            diff.append(int(a["voxel"]))
    res["task_e_arm_compared"] = n_cmp
    res["task_e_arm_identical"] = n_same
    res["task_e_arm_differing_voxels"] = diff
    row = dict(
        n=res["n"],
        wrong_base=f"{res['wrong_base']} ({100 * res['wrong_base_wilson95'][0]:.1f}-"
        f"{100 * res['wrong_base_wilson95'][1]:.1f}%)",
        wrong_integrated=f"{res['wrong_arm']} ({100 * res['wrong_arm_wilson95'][0]:.1f}-"
        f"{100 * res['wrong_arm_wilson95'][1]:.1f}%)",
        fixed=res["fixed"],
        broken=res["broken"],
        both_wrong=res["both_wrong"],
        broken_traced=f"{res['broken_traced_arm']} / {res['broken_traced_draw']} / "
        f"{res['broken_untraced']}",
        evals_base=res["mean_evals_base"],
        evals_integrated=res["mean_evals_arm"],
        smooth=res["mean_smooth_evals_arm"],
        identical_to_task_e=f"{n_same} of {n_cmp}",
    )
    res["tables"] = {
        "evalG2_e0": DT.markdown_table(
            [row],
            ["n", "wrong_base", "wrong_integrated", "fixed", "broken", "both_wrong",
             "broken_traced", "evals_base", "evals_integrated", "smooth", "identical_to_task_e"],
            formats=dict(evals_base=".0f", evals_integrated=".0f", smooth=".0f"),
        )  # fmt: skip
    }
    return res


# --------------------------------------------------------------------------------------------
def _weights() -> Dict[str, float]:
    f = json.loads((F_TIMING / "full.json").read_text())["ms_per_evaluation"]
    return dict(
        strict=f["strict"],
        screen=f["screen"],
        e4_value=f["e4_value"],
        e4_jac=f["e4_value_and_jacobian"],
    )


def strict_equiv(o: Dict[str, Any], w: Dict[str, float]) -> float:
    """Policy evaluations in strict-evaluation units with Task F's single-worker full-sample
    weights (every non-Jacobian smooth evaluation at F4's one-at-a-time E4 weight: an upper
    estimate for F5's batched ones)."""
    return (
        o["strict"]
        + o["screen"] * w["screen"] / w["strict"]
        + (o["smooth"] - o["jac"]) * w["e4_value"] / w["strict"]
        + o["jac"] * w["e4_jac"] / w["strict"]
    )


def bfs_rows(runs: List[Dict[str, Any]], w: Dict[str, float]) -> Dict[str, str]:
    acc, ev, wall, pol = [], [], [], []
    for r in runs:
        acc.append(
            dict(
                run=r["run"],
                n=r["n"],
                wrong=f"{r['wrong']} ({100 * r['wrong_wilson95'][0]:.1f}-"
                f"{100 * r['wrong_wilson95'][1]:.1f}%)",
                unresolved=r["unresolved"],
                grains_lost=len(r["grains_lost"]),
                median_err_right=r["median_err_right_deg"],
                seeds=r["seeds"],
                reseeds=f"{r['reseeds']} ({r['reseed_accepted']} accepted)",
                neighbor_fits=r["neighbor_fits"],
                accepted_first=r["neighbor_accepted_first"],
                revisits=f"{r['revisit_attempted']} ({r['revisit_accepted']} accepted)",
            )
        )
        e = r["evals"]
        ev.append(
            dict(
                run=r["run"],
                seed=e["seed"],
                reseed=e["reseed"],
                neighbor=e["neighbor"],
                revisit=e["revisit"],
                total=e["total"],
                smooth_full_searches=r.get("smooth_full_searches", 0),
            )
        )
        wl = r["wall"]
        wall.append(dict(run=r["run"], **{k: wl[k] for k in wl}))
        for ph, o in sorted(r.get("outward_log", {}).items()):
            pol.append(
                dict(
                    run=r["run"],
                    phase=ph,
                    runs=o["runs"],
                    accepted=o["accepted"],
                    fallback=o["fallback"],
                    fallback_accepted=o["fallback_accepted"],
                    strict=o["strict"],
                    screen=o["screen"],
                    smooth=o["smooth"],
                    jac=o["jac"],
                    strict_equiv=strict_equiv(o, w),
                    wall_s=o["total_wall_s"],
                    median_wall_s=o["median_wall_s"],
                    max_wall_s=o["max_wall_s"],
                )
            )
    return dict(
        accuracy=DT.markdown_table(
            acc,
            ["run", "n", "wrong", "unresolved", "grains_lost", "median_err_right", "seeds",
             "reseeds", "neighbor_fits", "accepted_first", "revisits"],
            formats=dict(median_err_right=".4f"),
        ),  # fmt: skip
        evals=DT.markdown_table(
            ev, ["run", "seed", "reseed", "neighbor", "revisit", "total", "smooth_full_searches"]
        ),
        wall=DT.markdown_table(
            wall,
            ["run", "seed", "reseed", "neighbor", "revisit", "total"],
            formats={k: ".0f" for k in ("seed", "reseed", "neighbor", "revisit", "total")},
        ),
        policy=DT.markdown_table(
            pol,
            ["run", "phase", "runs", "accepted", "fallback", "fallback_accepted", "strict",
             "screen", "smooth", "jac", "strict_equiv", "wall_s", "median_wall_s", "max_wall_s"],
            formats=dict(strict_equiv=".0f", wall_s=".0f", median_wall_s=".2f", max_wall_s=".1f"),
        ),  # fmt: skip
    )


def _with_smooth(r: Dict[str, Any], run_dir: Path) -> Dict[str, Any]:
    run = json.loads((run_dir / "run.json").read_text())
    r["smooth_full_searches"] = int(sum(run.get("smooth_evals_per_full_search", [])))
    return r


def paired_rows(prs: List[Dict[str, Any]]) -> str:
    return DT.markdown_table(prs, ["a", "b", "wrong_a", "wrong_b", "fixed", "broken", "both_wrong"])


def bfs_part(
    base: Path, others: List[Path], prefix: str, also_base: Optional[Path] = None
) -> Optional[Dict[str, Any]]:
    """Runs ``base`` + ``others`` (those that exist), each paired against ``base``; with
    ``also_base`` (a run in ``others``) the later runs are paired against it too."""
    dirs = [d for d in [base] + others if (d / "run.json").exists()]
    if len(dirs) < 2 or dirs[0] != base:
        return None
    w = _weights()
    runs = [_with_smooth(g2_bfs.run_summary(d), d) for d in dirs]
    prs = [g2_bfs.paired(base, d) for d in dirs[1:]]
    if also_base is not None and also_base in dirs:
        prs += [g2_bfs.paired(also_base, d) for d in dirs[dirs.index(also_base) + 1 :]]
    T = {f"{prefix}_{k}": v for k, v in bfs_rows(runs, w).items()}
    T[f"{prefix}_paired"] = paired_rows(prs)
    return dict(runs=runs, paired=prs, weights_ms=w, tables=T)


# --------------------------------------------------------------------------------------------
def main() -> None:
    summary: Dict[str, Any] = {}
    tables: Dict[str, str] = {}
    for name, fn in (("ledger", ledger), ("e0", e0)):
        res = fn()
        if res is not None:
            tables.update(res.pop("tables"))
            summary[name] = res
    pilot = CACHE / "pilot"
    res = bfs_part(
        pilot / "keyoff_branch", [pilot / "outward_only", pilot / "integrated"], "evalG2_pilot150"
    )
    if res is not None:
        tables.update(res.pop("tables"))
        res["keyoff_check"] = json.loads((OUT / "keyoff_pilot_compare.json").read_text())
        summary["pilot150"] = res
    bfs = CACHE / "bfs"
    # the stored pilot predates the centroid-screen default (2026-10-12), so the key-off rerun on
    # this branch is the like-for-like baseline; both are reported
    res = bfs_part(
        STORED,
        [bfs / "keyoff_rerun", bfs / "integrated", bfs / "outward_only"],
        "evalG2_bfs",
        also_base=bfs / "keyoff_rerun",
    )
    if res is not None:
        tables.update(res.pop("tables"))
        summary["bfs"] = res
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, default=float) + "\n")
    DT.write_tables(OUT / "tables.md", tables)
    print(f"parts: {sorted(summary)}; tables: {sorted(tables)}")


if __name__ == "__main__":
    main()
