"""Task G1 step 2b: Task E's deferred full-sample test on the other 24 voxels of the 30-voxel set.

The 30-voxel ray-study set (``benchmarks/phase_d_seed_diag/ray_overlap/voxels.json``) starts with
the six interior voxels that Task E ran; this script runs the other 24 (``VOXELS_30[6:]``). One
no-start search per voxel through Task E's probe (``smooth_rank_probe.task``: the Task C probe,
clean full-sample images, Q-max 8 configuration, seed 1, ``min_sin_eta`` 0), for

* the strict search (``--spec off --tag base``; the re-ranking of measurement 2 is skipped),
* the carried arms E3+E4 (tau 0.25 deg, ``--tag e34t25``) and E1+E2+E3+E4 (auto widths,
  ``--tag e1234a``) as the ``SmoothCostRank`` key of levels 0-2,

each under the lenient rule (``old``, library default) and the C++ rule (``cpp``).

  uv run python scripts/eval_reduction/g1_full24.py run --rule old --spec off --tag base
  uv run python scripts/eval_reduction/g1_full24.py run --rule cpp --spec "E3+E4,tau=0.25" \
      --tag e34t25 --workers 3
  uv run python scripts/eval_reduction/g1_full24.py summary

One npz per voxel (and a json with its wall time, contended) under
benchmarks/eval_reduction/G/cache/g1/full_<rule>_<tag>/ (gitignored); the summary writes
benchmarks/eval_reduction/G/g1/summary.json and tables.md. Wrong = error above 1 deg; every
broken voxel is traced with Task E's rule (``smooth_rank_summary.trace``, D-20261014-N13).
"""

import argparse
import json
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICE / "scripts" / "common"))

import smooth_rank_probe as SRP  # noqa: E402

G = ICE / "benchmarks" / "eval_reduction" / "G"
CACHE = G / "cache" / "g1"
OUT = G / "g1"
LOGS = G / "cache" / "logs"
VOXELS_24: List[int] = [int(v) for v in SRP.VOXELS_30[6:]]
ARMS: Dict[str, str] = {"e34t25": "E3+E4,tau=0.25", "e1234a": "E1+E2+E3+E4"}
RULES = ("old", "cpp")
RULE_LABEL = {"old": "lenient", "cpp": "C++"}
WRONG_DEG = 1.0


def run_dir(rule: str, tag: str) -> Path:
    return CACHE / f"full_{rule}_{tag}"


def _init() -> None:
    SRP._init("full")
    SRP._W["no_rerank"] = True


def task(item: Tuple[Any, ...]) -> str:
    t0 = time.time()
    msg = SRP.task(item)
    path = Path(item[-1])
    side = path.with_suffix(".json")
    if path.exists() and not side.exists():
        side.write_text(json.dumps(dict(wall_s=time.time() - t0, label="contended")) + "\n")
    return msg


def run(rule: str, spec: str, tag: str, workers: int, limit: int, reverse: bool) -> None:
    root = run_dir(rule, tag)
    root.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    (root / "spec.txt").write_text(spec + "\n")
    items = [
        ("full", v, p, rule, spec, str(root / f"v{v}.npz"))
        for p, v in enumerate(VOXELS_24)
        if not (root / f"v{v}.npz").exists()
    ][: limit or None]
    if reverse:
        items = items[::-1]
    print(f"{root.name} ({spec}): {len(items)} tasks, {workers} workers", flush=True)
    t0 = time.time()
    if workers <= 1:
        _init()
        for i, it in enumerate(items):
            print(f"[{i + 1}/{len(items)} {time.time() - t0:.0f}s] {task(it)}", flush=True)
    else:
        with mp.get_context("spawn").Pool(workers, initializer=_init) as pool:
            for i, msg in enumerate(pool.imap_unordered(task, items)):
                print(f"[{i + 1}/{len(items)} {time.time() - t0:.0f}s] {msg}", flush=True)
    (root / "run_meta.json").write_text(
        json.dumps(
            dict(
                workers=workers,
                wall_s=time.time() - t0,
                tasks=len(items),
                spec=spec,
                rule=rule,
                note="contended multi-worker wall, not a timing claim",
            )
        )
        + "\n"
    )
    sfx = "_rev" if reverse else ""
    (LOGS / f"g1_{root.name}{sfx}.done").write_text(f"{time.time() - t0:.0f}s\n")


# --------------------------------------------------------------------------------------------
# summary
# --------------------------------------------------------------------------------------------
def paired(rule: str, tag: str) -> Optional[Dict[str, Any]]:
    import stats as ST
    from smooth_rank_summary import fixed_trace, load_runs, trace

    base = load_runs(run_dir(rule, "base"))
    arm = load_runs(run_dir(rule, tag))
    vox = [v for v in VOXELS_24 if v in base and v in arm]
    if not vox:
        return None
    eb = np.array([float(base[v]["err_deg"]) for v in vox])
    ea = np.array([float(arm[v]["err_deg"]) for v in vox])
    wb, wa = eb > WRONG_DEG, ea > WRONG_DEG
    n_broken, n_fixed = ST.paired_discordant(wa, wb)  # (arm wrong only, baseline wrong only)
    broken = [v for v, x, y in zip(vox, wa, wb) if x and not y]
    fixed = [v for v, x, y in zip(vox, wa, wb) if y and not x]
    tr = {v: trace(arm[v], SRP.RANK_LEVELS) for v in broken}
    ft = {v: fixed_trace(arm[v], SRP.RANK_LEVELS) for v in fixed}
    both_right = ~wa & ~wb
    n = len(vox)
    return dict(
        rule=rule,
        tag=tag,
        spec=ARMS[tag],
        n=n,
        voxels=vox,
        wrong_base=int(wb.sum()),
        wrong_arm=int(wa.sum()),
        wrong_base_wilson95=list(ST.wilson(int(wb.sum()), n)),
        wrong_arm_wilson95=list(ST.wilson(int(wa.sum()), n)),
        fixed=int(n_fixed),
        broken=int(n_broken),
        both_wrong=int((wa & wb).sum()),
        net=int(n_fixed) - int(n_broken),
        fixed_voxels=fixed,
        broken_voxels=broken,
        wrong_base_voxels=[v for v, x in zip(vox, wb) if x],
        wrong_arm_voxels=[v for v, x in zip(vox, wa) if x],
        broken_traced_arm=sum(1 for c, _ in tr.values() if c == "arm"),
        broken_traced_draw=sum(1 for c, _ in tr.values() if c == "draw"),
        broken_untraced=sum(1 for c, _ in tr.values() if c == "untraced"),
        broken_trace={str(v): list(t) for v, t in tr.items()},
        fixed_traced_arm=sum(1 for c in ft.values() if c == "arm"),
        fixed_traced_draw=sum(1 for c in ft.values() if c == "draw"),
        fixed_trace={str(v): c for v, c in ft.items()},
        err_base={str(v): float(e) for v, e in zip(vox, eb)},
        err_arm={str(v): float(e) for v, e in zip(vox, ea)},
        n_both_right=int(both_right.sum()),
        strict_evals_change_pct=round(
            100.0
            * (
                np.mean([int(arm[v]["evals"]) for v in vox])
                / np.mean([int(base[v]["evals"]) for v in vox])
                - 1.0
            ),
            1,
        ),
        median_err_both_right_base=float(np.median(eb[both_right])) if both_right.any() else None,
        median_err_both_right_arm=float(np.median(ea[both_right])) if both_right.any() else None,
        mean_evals_base=float(np.mean([int(base[v]["evals"]) for v in vox])),
        mean_evals_arm=float(np.mean([int(arm[v]["evals"]) for v in vox])),
        mean_smooth_evals_arm=float(np.mean([int(arm[v]["smooth_evals"]) for v in vox])),
    )


def decision(P: Dict[str, Any]) -> Dict[str, Any]:
    """The default of the plan's decision point: the arm with more fixed minus broken voxels across
    both rules; ``off`` if neither arm fixes more than it breaks."""
    net: Dict[str, Optional[int]] = {}
    complete: Dict[str, bool] = {}
    for tag in ARMS:
        xs = [P.get(f"{rule}|{tag}") for rule in RULES]
        complete[tag] = all(x is not None and x["n"] == len(VOXELS_24) for x in xs)
        net[tag] = sum(x["net"] for x in xs if x is not None) if any(xs) else None
    cands = {t: n for t, n in net.items() if n is not None and n > 0}
    ms = ms_per_eval()
    if not cands:
        choice = "off"
        reason = "neither arm fixes more voxels than it breaks across both rules"
    else:
        best = max(cands.values())
        tops = [t for t, n in cands.items() if n == best]
        # tie-break (as Task E's decision point): the cheaper arm per smooth evaluation
        tops.sort(key=lambda t: ms.get(t, float("inf")))
        choice = ARMS[tops[0]]
        reason = f"net fixed minus broken across both rules: {best}" + (
            f"; tie between {', '.join(ARMS[t] for t in tops)} broken by ms per evaluation"
            if len(tops) > 1
            else ""
        )
    return dict(
        net_across_rules=net,
        complete=complete,
        choice=choice,
        reason=reason,
        tie_break_ms_per_eval_batch256_full=ms,
    )


def ms_per_eval() -> Dict[str, float]:
    """Task E measurement 4 (single-worker, full-sample images, batch of 256): ms per smooth
    evaluation of each carried arm (tau enters only the final ratio, so E3+E4 stands for tau 0.25).
    """
    f = ICE / "benchmarks" / "eval_reduction" / "E" / "timing" / "timing_full.json"
    arms = json.loads(f.read_text())["arms"]
    names = {"e34t25": "E3+E4", "e1234a": "E1+E2+E3+E4"}
    return {t: float(arms[n]["ms_per_eval_batch256"]) for t, n in names.items() if n in arms}


def tables(P: Dict[str, Any], base_wrong: Dict[str, Any]) -> Dict[str, str]:
    import doc_tables as DT

    rows = []
    for key, x in P.items():
        lo_b, hi_b = x["wrong_base_wilson95"]
        lo_a, hi_a = x["wrong_arm_wilson95"]
        rows.append(
            {
                "rule": RULE_LABEL[x["rule"]],
                "arm (SmoothCostRank)": x["spec"],
                "n": x["n"],
                "wrong, strict": f"{x['wrong_base']} ({100 * lo_b:.1f}-{100 * hi_b:.1f}%)",
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
                "median error, right in both (deg), strict / arm": (
                    "-"
                    if x["median_err_both_right_base"] is None
                    else f"{x['median_err_both_right_base']:.4f} / "
                    f"{x['median_err_both_right_arm']:.4f}"
                ),
                "strict evaluations per voxel, strict search": x["mean_evals_base"],
                "strict evaluations per voxel, arm": x["mean_evals_arm"],
                "smooth evaluations per voxel": x["mean_smooth_evals_arm"],
            }
        )
    T: Dict[str, str] = {}
    if rows:
        fm = {k: ".0f" for k in rows[0] if "evaluations" in k}
        T["evalG1_full24"] = DT.markdown_table(rows, list(rows[0].keys()), fm)
    vrows = []
    for v in VOXELS_24:
        row: Dict[str, Any] = {"voxel": v}
        for rule in RULES:
            b = base_wrong.get(rule, {}).get(str(v))
            row[f"strict, {RULE_LABEL[rule]} (deg)"] = "-" if b is None else f"{b:.3f}"
            for tag in ARMS:
                x = P.get(f"{rule}|{tag}")
                e = None if x is None else x["err_arm"].get(str(v))
                row[f"{tag}, {RULE_LABEL[rule]} (deg)"] = "-" if e is None else f"{e:.3f}"
        vrows.append(row)
    T["evalG1_full24_voxels"] = DT.markdown_table(vrows, list(vrows[0].keys()), {})
    return T


def basin_ranks(v: int) -> Dict[str, Any]:
    """Per run and level of a voxel: list length after the quick MC, the number kept, the basin
    radius (deg) and the 0-based rank of the first truth-basin candidate in the recorded (run's own)
    order and under the strict cost's order of the same list (-1: no basin candidate)."""
    from smooth_rank_summary import basin_radius

    out: Dict[str, Any] = {}
    for rule in RULES:
        for tag in ("base",) + tuple(ARMS):
            f = run_dir(rule, tag) / f"v{v}.npz"
            if not f.exists():
                continue
            z = dict(np.load(f))
            levels = []
            for lv in range(int(z["n_levels"])):
                d = z[f"L{lv}_dist"].astype(np.float64)
                c = z[f"L{lv}_cost"].astype(np.float64)
                B = basin_radius(z, lv, "next")
                own = np.nonzero(d <= B)[0]
                strict = np.nonzero(d[np.argsort(c, kind="stable")] <= B)[0]
                levels.append(
                    dict(
                        level=lv,
                        n=int(len(d)),
                        n_keep=int(z[f"L{lv}_n_keep"]),
                        basin_deg=round(B, 2),
                        first_basin_rank=int(own[0]) if len(own) else -1,
                        first_basin_rank_strict_order=int(strict[0]) if len(strict) else -1,
                    )
                )
            out[f"{rule}|{tag}"] = dict(err_deg=float(z["err_deg"]), levels=levels)
    return out


def summary() -> None:
    import doc_tables as DT
    from smooth_rank_summary import load_runs

    P: Dict[str, Any] = {}
    for rule in RULES:
        for tag in ARMS:
            x = paired(rule, tag)
            if x is not None:
                P[f"{rule}|{tag}"] = x
    base_err = {
        rule: {str(v): float(z["err_deg"]) for v, z in load_runs(run_dir(rule, "base")).items()}
        for rule in RULES
    }
    wrong_any = sorted(
        {int(v) for x in P.values() for v in x["wrong_base_voxels"] + x["wrong_arm_voxels"]}
    )
    S = dict(
        voxels=VOXELS_24,
        wrong_deg=WRONG_DEG,
        paired=P,
        strict_err=base_err,
        basin_ranks_of_wrong_voxels={str(v): basin_ranks(v) for v in wrong_any},
        decision=decision(P),
        note="no-start search, clean full-sample images, Q-max 8 config, seed 1, min_sin_eta 0",
    )
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(S, indent=1, default=float) + "\n")
    T = tables(P, base_err)
    DT.write_tables(OUT / "tables.md", T)
    for k, t in T.items():
        print(k)
        print(t)
        print()
    print(json.dumps(S["decision"], indent=1))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("cmd", choices=["run", "summary"])
    ap.add_argument("--rule", default="old", choices=list(RULES))
    ap.add_argument("--spec", default="off")
    ap.add_argument("--tag", default="base")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--reverse", action="store_true")
    a = ap.parse_args()
    if a.cmd == "summary":
        summary()
        return
    if a.tag != "base" and ARMS.get(a.tag) != a.spec:
        raise SystemExit(f"tag {a.tag} expects spec {ARMS.get(a.tag)!r}, got {a.spec!r}")
    run(a.rule, a.spec, a.tag, a.workers, a.limit, a.reverse)


if __name__ == "__main__":
    main()
