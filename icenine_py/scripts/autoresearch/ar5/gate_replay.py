"""AR5 step 7, Phase A: replay the skip on the stored passes of the disclosed samples.

  uv run python scripts/autoresearch/ar5/gate_replay.py

For each stored pass (``gate_common.DISCLOSED``) and each skip cut ``h`` (the frozen 0.95, and a
what-if curve 0.90 to 0.99): remove the tries of the voxels whose own hit ratio is at least ``h``
(a skipped voxel keeps its end from before the pass), then score the map against the truth.
Stop rule checks, written to the output: the stored pass evaluations equal the plan's totals; the
replay at ``h`` above 1 (nothing skipped) reproduces the stored map after the pass; the replay at
0.95 meets the plan's decision rule on the pass-level parts. Truth is used for scoring only.
Writes ``benchmarks/autoresearch/ar5/gate/replay.json`` (slim, tracked).
"""

import json
from typing import Any, Dict, List, Sequence

import gate_common as GC
import gate_rules as GR
import numpy as np

import stats as ST  # noqa: E402  (scripts/common is put on the path by gate_common)

H_CURVE = [0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99]


def side(err: np.ndarray) -> Dict[str, Any]:
    k = int((err > GC.WRONG).sum())
    lo, hi = ST.wilson(k, len(err))
    return dict(wrong=k, wilson_lo=lo, wilson_hi=hi, median_err_deg=float(np.median(err)))


def replay_map(P: Dict[str, Any], rep: Dict[str, Any]) -> np.ndarray:
    """The map after the replayed pass: the stored map after the full pass, with every skipped
    voxel reset to its end before the pass. (A kept voxel has all its tries, so its end is the
    stored one.)"""
    R = P["after"].copy()
    sk = rep["skipped"]
    R[sk] = P["before"][sk]
    return R


def check_map(P: Dict[str, Any], rep0: Dict[str, Any]) -> Dict[str, Any]:
    """The replay with nothing skipped: the winner of every boundary voxel is the end in the
    stored map after the pass (float32 end of the try, or the old end)."""
    bad = 0
    n_changed = 0
    ends = P["end"].astype(np.float64)
    for v in np.flatnonzero(np.isin(np.arange(len(P["after"])), P["voxel"])):
        w = int(rep0["winner"][v])
        want = P["before"][v] if w < 0 else ends[w]
        if not np.allclose(P["after"][v], want, atol=1e-6, rtol=0):
            bad += 1
        n_changed += int(w >= 0)
    return dict(boundary_voxels=int(rep0["boundary_voxels"]), mismatched=bad, changed=n_changed)


def run_rows(
    sample: str,
    run: str,
    P: Dict[str, Any],
    err: Any,
    s8: Dict[str, Any],
    hs: Sequence[float],
) -> Dict[str, Any]:
    """All numbers of one stored pass. ``err(R)`` scores a map; ``s8`` is the S8 baseline."""
    old_hit = P["old_score"][:, 1]
    old_cost = P["old_score"][:, 0]
    e_before, e_full = err(P["before"]), err(P["after"])
    w_before, w_full = e_before > GC.WRONG, e_full > GC.WRONG
    w_s8 = s8["err_after"] > GC.WRONG
    rep0 = GR.replay(old_hit, old_cost, P["voxel"], P["hit"], P["cost"], P["evals"], 2.0)
    fixed_by_pass = np.flatnonzero(w_before & ~w_full)
    bnd = rep0["boundary_voxels"]
    isb = np.zeros(len(old_hit), dtype=bool)
    isb[P["voxel"]] = True
    out: Dict[str, Any] = dict(
        sample=sample,
        run=run,
        n=len(old_hit),
        pass_evals_full=rep0["evals_full"],
        pass_json_evals=int(P["pass_json"]["evals"]),
        n_tries=int(len(P["voxel"])),
        boundary_voxels=bnd,
        bfs_evals=P["bfs_evals"],
        s8_total=int(s8["stages"]["total"]),
        full=dict(before=side(e_before), after=side(e_full)),
        s8=dict(after=side(s8["err_after"])),
        map_check=check_map(P, rep0),
        own_hit_of_fixed=[round(float(x), 4) for x in sorted(old_hit[fixed_by_pass])],
        max_own_hit_wrong_before=float(old_hit[w_before & isb].max()),
        max_own_hit_wrong_after=float(old_hit[w_full & isb].max()),
        by_h=[],
    )
    fx8, br8 = GR.fixed_broken(w_s8, w_full)
    out["full_vs_s8"] = dict(fixed=fx8, broken=br8)
    for h in hs:
        rep = GR.replay(old_hit, old_cost, P["voxel"], P["hit"], P["cost"], P["evals"], h)
        e_skip = err(replay_map(P, rep))
        w_skip = e_skip > GC.WRONG
        fixed, broken = GR.fixed_broken(w_full, w_skip)
        fx8h, br8h = GR.fixed_broken(w_s8, w_skip)
        total_skip = P["bfs_evals"] + rep["evals"]
        s = side(e_skip)
        rule = GR.gate_rule(
            cut=rep["cut"],
            wrong_before=int(w_before.sum()),
            wrong_full=int(w_full.sum()),
            fixed=fixed,
            broken=broken,
            med_full=out["full"]["after"]["median_err_deg"],
            med_skip=s["median_err_deg"],
            total_skip=total_skip,
            total_s8=out["s8_total"],
            fixed_s8=fx8h,
            broken_s8=br8h,
            med_s8=out["s8"]["after"]["median_err_deg"],
        )
        out["by_h"].append(
            dict(
                h=h,
                pass_evals=rep["evals"],
                cut=rep["cut"],
                tries_skipped=rep["tries_skipped"],
                voxels_skipped=rep["voxels_skipped"],
                wrong=s["wrong"],
                wilson_lo=s["wilson_lo"],
                wilson_hi=s["wilson_hi"],
                median_err_deg=s["median_err_deg"],
                fixed_vs_full=fixed,
                broken_vs_full=broken,
                fixed_vs_s8=fx8h,
                broken_vs_s8=br8h,
                total_with_skip=total_skip,
                total_full=P["bfs_evals"] + rep["evals_full"],
                net_gain_full=GR.net_gain(int(w_before.sum()), int(w_full.sum())),
                net_gain_lost=GR.net_gain(int(w_before.sum()), int(w_full.sum()))
                - GR.net_gain(int(w_before.sum()), s["wrong"]),
                rule=rule,
            )
        )
    return out


def main() -> None:
    rows: List[Dict[str, Any]] = []
    stop: List[str] = []
    for sample, run in GC.DISCLOSED:
        P = GC.load_pass(sample, run)
        T = GC.Truth(sample)
        s8 = GC.s8_baseline(sample)
        r = run_rows(sample, run, P, T.err, s8, H_CURVE)
        want = GC.PLAN_PASS_EVALS[(sample, run)]
        if r["pass_evals_full"] != want or r["pass_json_evals"] != want:
            stop.append(f"{sample} {run}: stored total {r['pass_evals_full']} != plan {want}")
        at95 = [x for x in r["by_h"] if abs(x["h"] - GR.H_S) < 1e-9][0]
        if round(100 * at95["cut"], 1) != GC.PLAN_CUT_PCT[(sample, run)]:
            stop.append(f"{sample} {run}: cut {100 * at95['cut']:.2f} != plan")
        if r["map_check"]["mismatched"]:
            stop.append(f"{sample} {run}: replay of the full pass != stored map")
        rows.append(r)
        print(
            f"{sample:5s} {run:5s} pass {r['pass_evals_full']:>9,} "
            f"cut@0.95 {100 * at95['cut']:.1f}% "
            f"wrong {r['full']['after']['wrong']} -> {at95['wrong']} fixed/broken "
            f"{at95['fixed_vs_full']}/{at95['broken_vs_full']} median "
            f"{r['full']['after']['median_err_deg']:.5f} -> {at95['median_err_deg']:.5f} "
            f"rule {at95['rule']['passes']}",
            flush=True,
        )
    GC.GATE.mkdir(parents=True, exist_ok=True)
    (GC.GATE / "replay.json").write_text(
        json.dumps(dict(h_s=GR.H_S, rows=rows, stop_rule_flags=stop), indent=1) + "\n"
    )
    print("STOP RULE FLAGS:", stop or "none")


if __name__ == "__main__":
    main()
