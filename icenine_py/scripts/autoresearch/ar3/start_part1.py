"""AR3 step 9, part 1: what the stored runs say about the start of a neighbour fit, and class A.

  scripts/autoresearch/ar3/start_bg.sh big part1_big start_part1.py     # logs to start/cache/logs
  scripts/autoresearch/ar3/start_bg.sh grad part1_grad start_part1.py
  ... --limit N     # pilot: the first N fit rows and no class A step

Reads the logged run of ``S5`` (``fits.npz``: one row per local optimization of the BFS, with the
start and end orientation, the end's hit ratio, the flood's best hit ratio before the fit and the
evaluations): step 8's for ``big`` (the worktree of step 8, read only), this task's rerun for
``grad`` (``<RESEED_CACHE>/S5_grad_logged``). Nothing in a fit is run again except where said.

Part 1.1. The hit ratio of each fit's start is one strict evaluation at the stored start
orientation (analysis only; it is not in any evaluation total). With the truth, the error at the
start and at the end. Written to ``start/feat_<sample>.npz`` (read by the part 2 scripts) and
summarized in ``start/part1_<sample>.json`` by phase (neighbour, revisit), by gradient level of the
voxel's true grain: the fits, the share whose start has a hit ratio of at least 0.90, 0.95 and
0.98, the evaluations, and for those fits how far they moved (angle start to end) and how their
error changed.

Part 1.2, class A (``reseed_slim.py``'s definition: the last fit of a voxel of a reseed cluster
before the reseed, ended right, error at most 1 degree, and was rejected by the relative test; a
voxel is counted at its first cluster). For each: the hit ratio at the true orientation against
0.9 times the flood's best at that moment (the relative test is hit / max(best, hit) > 0.9, so
the truth passes when its hit ratio is above 0.9 times the best), and then the same after one more
local optimization from the rejected fit's end (CMA, first step 0.1 degrees, at most 60
evaluations, ``generator default_rng(0)``). Whether the stored fit includes the BFS's own retry is
read from the run's configuration (``retry_on`` needs ``local_optimizer == "cma"``).
"""

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import reseed_common as RC

import boundary_common as B  # noqa: E402
import ar3_common as C  # noqa: E402

START = C.ICE / "benchmarks" / "autoresearch" / "ar3" / "start"
STEP8 = RC.WORKTREES / "agent-a635b2f904279bdc9/icenine_py/benchmarks/autoresearch/ar3/reseed/cache"
LOGGED = {"big": STEP8 / "S5_big_logged", "grad": RC.CACHE / "S5_grad_logged"}
THRESH = (0.90, 0.95, 0.98)
SIGMA_A = 0.1  # degrees: the extra local optimization of class A (the plan)
MAX_EVALS = 60
PHASES = {0: "neighbour", 1: "revisit", 2: "restart"}


def extra_fit(env: Any, v: int, R0: np.ndarray, sigma: float) -> Dict[str, Any]:
    """One more local optimization: CMA, first step ``sigma`` degrees, at most 60 evaluations (the
    start's own included) plus the final one; returns the end, its hit ratio and the evaluations."""
    rec = env.hooks.rec
    vox = env.hooks.mic.voxels[v]
    old = rec.params.local_optimizer
    rec.params.local_optimizer = "cma"
    try:
        r = rec.local_optimization(
            voxel_vertices=env.hooks.vertices(v),
            phase_index=vox.phase,
            initial_orientation=np.asarray(R0, np.float64),
            rng=np.random.default_rng(0),
            cma_sigma0_deg=sigma,
            cma_max_evals=MAX_EVALS,
        )
    finally:
        rec.params.local_optimizer = old
    hit, conf = env.bfs._ratios(r.overlap_info)
    return dict(
        end=np.asarray(r.orientation, np.float64),
        hit=float(hit),
        evals=int(rec.last_local_optimization_evals),
    )


def features(env: Any, d: Path, limit: int) -> Dict[str, np.ndarray]:
    z = np.load(d / "fits.npz")
    meta, R = z["meta"], z["R"].astype(np.float64)
    n = len(meta) if not limit else min(limit, len(meta))
    v = meta[:n, 0].astype(int)
    Rs = B.orthonormal(R[:n, :9].reshape(-1, 3, 3))
    Re = B.orthonormal(R[:n, 9:].reshape(-1, 3, 3))
    start_hit = np.zeros(n)
    end_hit_check = np.zeros(n)
    t0 = time.time()
    for i in range(n):
        start_hit[i] = env.score3(Rs[i], int(v[i]))[1]
        end_hit_check[i] = env.score3(Re[i], int(v[i]))[1]
        if i % 2000 == 0:
            print(f"  features {i}/{n} {time.time() - t0:.0f}s", flush=True)
    err_s, err_e = env.err_rows(Rs, v), env.err_rows(Re, v)
    move = env.mis(Rs, Re)
    return dict(
        voxel=v.astype(np.int32),
        flood=meta[:n, 1].astype(np.int32),
        expander=meta[:n, 2].astype(np.int32),
        phase=meta[:n, 3].astype(np.int8),
        hit_end=meta[:n, 4].astype(np.float32),
        best_before=meta[:n, 5].astype(np.float32),
        evals=meta[:n, 6].astype(np.int32),
        start_hit=start_hit.astype(np.float32),
        end_hit_check=end_hit_check.astype(np.float32),
        err_start=err_s.astype(np.float32),
        err_end=err_e.astype(np.float32),
        move=move.astype(np.float32),
        grain=np.asarray(env.gmap)[v].astype(np.int32),
        level=np.asarray(env.level)[v].astype(np.float32),
    )


def pct(a: np.ndarray, q: float) -> float:
    return float(np.percentile(a, q)) if len(a) else float("nan")


def by_level(f: Dict[str, np.ndarray], level_of_voxel: np.ndarray) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    lv = level_of_voxel[f["voxel"]]
    for ph in (0, 1):
        for L in ["all"] + sorted(set(lv.tolist())):
            sel = (f["phase"] == ph) & (np.ones(len(lv), bool) if L == "all" else lv == L)
            n = int(sel.sum())
            row: Dict[str, Any] = dict(fits=n, evals=int(f["evals"][sel].sum()))
            for t in THRESH:
                s2 = sel & (f["start_hit"] >= t)
                row[f"start_ge_{t:.2f}"] = int(s2.sum())
                row[f"evals_of_those_{t:.2f}"] = int(f["evals"][s2].sum())
                row[f"move_median_{t:.2f}"] = pct(f["move"][s2], 50)
                row[f"move_p90_{t:.2f}"] = pct(f["move"][s2], 90)
                row[f"hit_gain_median_{t:.2f}"] = pct((f["hit_end"] - f["start_hit"])[s2], 50)
                row[f"wrong_start_{t:.2f}"] = int((s2 & (f["err_start"] > RC.RIGHT)).sum())
                row[f"right_end_wrong_start_{t:.2f}"] = int(
                    (s2 & (f["err_start"] > RC.RIGHT) & (f["err_end"] <= RC.RIGHT)).sum()
                )
                row[f"right_start_wrong_end_{t:.2f}"] = int(
                    (s2 & (f["err_start"] <= RC.RIGHT) & (f["err_end"] > RC.RIGHT)).sum()
                )
            row["err_start_median"] = pct(f["err_start"][sel], 50)
            row["err_end_median"] = pct(f["err_end"][sel], 50)
            row["start_hit_median"] = pct(f["start_hit"][sel], 50)
            out[f"{PHASES[ph]}/{L}"] = row
    return out


def class_a(env: Any, d: Path, f: Dict[str, np.ndarray]) -> Dict[str, Any]:
    """Part 1.2 on the voxels of class A (first cluster of each voxel)."""
    clusters = RC.read_gz(d / "clusters.json.gz")
    floods = RC.read_gz(d / "floods.json.gz")
    reseed_floods = [x for x in floods if x["reseed"]]
    assert len(reseed_floods) == len(clusters)
    flood = f["flood"]
    first_row: Dict[int, int] = {}
    for i, fl in enumerate(flood):
        first_row.setdefault(int(fl), i)
    n = len(flood)
    seen: set = set()
    rows: List[Dict[str, Any]] = []
    n_members = 0
    for k, (c, fl) in enumerate(zip(clusters, reseed_floods)):
        lim = first_row.get(fl["flood"], n)
        last: Dict[int, int] = {}
        for i in range(lim):
            last[int(f["voxel"][i])] = i
        for m in c["members"]:
            m = int(m)
            if m in seen:
                continue
            seen.add(m)
            n_members += 1
            i = last.get(m)
            if i is None or f["err_end"][i] > RC.RIGHT:
                continue  # class N, B or C
            hit, best = float(f["hit_end"][i]), float(f["best_before"][i])
            best_a = max(hit, best)
            if not (best_a > 0 and hit / best_a <= 0.9):
                continue  # the last fit was accepted: not a rejected fit
            rows.append(dict(k=k, voxel=m, row=int(i), hit=hit, best=best_a))
    Rv = env.Rtrue
    Re_all = np.load(d / "fits.npz")["R"][:, 9:].astype(np.float64)
    out_rows: List[Dict[str, Any]] = []
    t0 = time.time()
    for j, r in enumerate(rows):
        v, i = r["voxel"], r["row"]
        hit_truth = env.score3(Rv[v], v)[1]
        # the end of the stored fit is the start of the extra fit
        out_rows.append(r)
        r["hit_truth"] = float(hit_truth)
        R_end = B.orthonormal(Re_all[i].reshape(1, 3, 3))[0]
        ef = extra_fit(env, v, R_end, SIGMA_A)
        r["hit_after"] = ef["hit"]
        r["err_after"] = float(env.err_rows(ef["end"][None], np.array([v]))[0])
        r["evals_extra"] = ef["evals"]
        r["err_end"] = float(f["err_end"][i])
        r["level"] = float(env.level[v])
        if j % 200 == 0:
            print(f"  class A {j}/{len(rows)} {time.time() - t0:.0f}s", flush=True)
    a = out_rows
    best = np.array([x["best"] for x in a])
    hit = np.array([x["hit"] for x in a])
    ht = np.array([x["hit_truth"] for x in a])
    ha = np.array([x["hit_after"] for x in a])
    truth_pass = ht > 0.9 * best
    after_pass = ha / np.maximum(ha, best) > 0.9
    summary = dict(
        voxels_in_clusters_first_appearance=n_members,
        class_a=len(a),
        truth_passes_relative_test=int(truth_pass.sum()),
        truth_fails_relative_test=int((~truth_pass).sum()),
        after_extra_fit_passes=int(after_pass.sum()),
        after_pass_of_truth_pass=int((after_pass & truth_pass).sum()),
        after_pass_of_truth_fail=int((after_pass & ~truth_pass).sum()),
        truth_fail_and_after_fail=int((~after_pass & ~truth_pass).sum()),
        truth_pass_and_after_fail=int((~after_pass & truth_pass).sum()),
        median_hit_stored=float(np.median(hit)) if len(a) else None,
        median_best=float(np.median(best)) if len(a) else None,
        median_hit_truth=float(np.median(ht)) if len(a) else None,
        median_hit_after=float(np.median(ha)) if len(a) else None,
        median_err_end=float(np.median([x["err_end"] for x in a])) if len(a) else None,
        median_err_after=float(np.median([x["err_after"] for x in a])) if len(a) else None,
        evals_extra_total=int(sum(x["evals_extra"] for x in a)),
        by_level={
            str(L): dict(
                n=int(sum(1 for x in a if x["level"] == L)),
                truth_pass=int(sum(1 for x, p in zip(a, truth_pass) if x["level"] == L and p)),
                after_pass=int(sum(1 for x, p in zip(a, after_pass) if x["level"] == L and p)),
            )
            for L in sorted(set(x["level"] for x in a))
        },
    )
    return dict(summary=summary, rows=out_rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    START.mkdir(parents=True, exist_ok=True)
    d = LOGGED[RC.SAMPLE]
    env = RC.make_env()
    env.err_rows = lambda Rm, v: np.asarray(env.T.err(Rm, v), dtype=float)
    f = features(env, d, a.limit)
    tag = RC.SAMPLE + ("_pilot" if a.limit else "")
    np.savez_compressed(START / f"feat_{tag}.npz", **f)
    chk = dict(
        rows=len(f["voxel"]),
        max_abs_end_hit_diff=float(np.max(np.abs(f["end_hit_check"] - f["hit_end"]))),
        n_end_hit_diff_gt_1e3=int((np.abs(f["end_hit_check"] - f["hit_end"]) > 1e-3).sum()),
    )
    out: Dict[str, Any] = dict(sample=RC.SAMPLE, check_end_hit=chk)
    out["by_level"] = by_level(f, np.asarray(env.level))
    if not a.limit:
        ca = class_a(env, d, f)
        out["class_a"] = ca["summary"]
        (START / f"classA_{RC.SAMPLE}.json").write_text(json.dumps(ca["rows"]) + "\n")
    (START / f"part1_{tag}.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(json.dumps(chk))
    if "class_a" in out:
        print(json.dumps(out["class_a"], indent=1)[:2500])
    sys.stdout.flush()


if __name__ == "__main__":
    main()
