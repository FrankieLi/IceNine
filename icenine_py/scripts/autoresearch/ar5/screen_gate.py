"""AR5 step 21, phase A step 4 (gates g0, g1, g2): score the gate runs of the six disclosed samples.

  uv run python scripts/autoresearch/ar5/screen_gate.py

Reads, per sample (``screen_gate_runs.sh`` made them): ``<cache>/full_<s>/`` (J5's pass re-run with
its own budget, 60, and the gain trace), ``<cache>/J5C_<s>/`` (step 20's repair on the re-run's
tries: the full-tries replay), ``<cache>/J5S_<s>/`` (the screen pass and the repair on the wide
flag), and the stored files (links): J5's tries and maps, S8's map, step 20's J5C. The truth is
read for scoring only. Checks:

* g0: the re-run reproduces J5's stored tries (count, voxel, candidate, evaluations; cost and hit
  ratio to 1e-12; ends to 1e-6, the stored ones are float32) and the stored after-pass map; the
  replay reproduces J5C: repair evaluations (137,214 / 146,894 / 124,146 / 76,835 / 116,765 /
  19,602), the wrong voxels and the map.
* validity: the screen tries equal the re-run's trace cut at one evaluation (first strict gain),
  try by try: end, cost, hit ratio exactly; 2 evaluations a try; end equal to the candidate's start.
* g1: prediction 1 exactly (pass evaluations, flagged, repair evaluations, totals, fixed / broken
  against J5C).
* g2: per sample at least 10% fewer evaluations than J5C, broken minus fixed against J5C at most 0,
  median error not higher; on mid broken minus fixed against S8 at most +1.

Writes ``benchmarks/autoresearch/ar5/screen/gate.json`` when it all passes, else
``gate_fail.json`` (and exits 1): then stop, no tuning, no other form, no held5. Also writes
``gate_tables.md``. Wall times are contended and are not timing results.
"""

import argparse
import json
import sys
from typing import Any, Dict, List

import screen_common as SC  # noqa: E402  (first)

SC.SF.check_or_exit()
import numpy as np  # noqa: E402
import screen_data as SD  # noqa: E402
import screen_rules as SR  # noqa: E402

sys.path.insert(0, str(SC.ICE / "scripts" / "common"))
import stats as ST  # noqa: E402

TOL_COST = 1e-12
TOL_END = 1e-6


def wilson(k: int, n: int) -> List[float]:
    lo, hi = ST.wilson(k, n)
    return [float(lo), float(hi)]


def load_run(name: str) -> Dict[str, Any]:
    d = SC.RCACHE / name
    out: Dict[str, Any] = dict(rep=json.loads((d / "repair.json").read_text()))
    out["after"] = np.load(d / "after.npy").astype(np.float64)
    out["after_pass"] = np.load(d / "after_pass.npy").astype(np.float64)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--samples", default=",".join(SR.SAMPLES), help="a subset: a test; no gate.json"
    )
    a = ap.parse_args()
    use = tuple(a.samples.split(","))
    partial = use != tuple(SR.SAMPLES)
    y = SD.sym()
    rows: Dict[str, Any] = {}
    g0: Dict[str, Any] = {}
    g1: Dict[str, Any] = {}
    g2: List[Dict[str, Any]] = []
    valid: Dict[str, Any] = {}
    for s in use:
        c = SC.RCACHE
        d = SD.load_inputs(s, truth=True)
        Rt = d["R_t"]
        st = SD.stored_stats(s)
        j5, j5c_tot, s8_tot = SR.TOTALS[s]
        stored = d["stored"]
        passes = c / f"passes_{s}"
        R_j5 = np.load(passes / "after_J5.npy").astype(np.float64)
        R_s8 = np.load(passes / "after_S8.npy").astype(np.float64)
        R_j5c = np.load(c / f"stored_J5C_{s}" / "after.npy").astype(np.float64)
        rep_j5c_stored = json.loads((c / f"stored_J5C_{s}" / "repair.json").read_text())

        def err(R: np.ndarray) -> np.ndarray:
            return np.asarray(SR.misorientation_deg(R, Rt, y), dtype=float)

        e_j5, e_s8, e_j5c = err(R_j5), err(R_s8), err(R_j5c)
        w_j5, w_s8, w_j5c = (e > SR.WRONG_DEG for e in (e_j5, e_s8, e_j5c))

        # ---- g0 (a): the re-run of J5's tries
        full = np.load(c / f"full_{s}" / "tries.npz")
        n = len(stored["voxel"])
        same_n = bool(len(full["voxel"]) == n)
        same_ids = same_n and bool(
            (full["voxel"] == stored["voxel"]).all()
            and (full["cand"] == stored["cand"]).all()
            and (full["evals"] == stored["evals"]).all()
        )
        d_cost = float(np.abs(full["cost"] - stored["cost"]).max()) if same_n else np.inf
        d_hit = float(np.abs(full["hit"] - stored["hit"]).max()) if same_n else np.inf
        d_end = float(np.abs(full["end"] - stored["end"]).max()) if same_n else np.inf
        R_after_full, _c, _h = SR.after_pass_state(
            d["R_before"], d["old_score"], full["voxel"], full["end"], full["cost"], full["hit"]
        )
        d_after = float(np.abs(R_after_full - R_j5).max())
        rerun_ok = bool(
            same_ids
            and d_cost <= TOL_COST
            and d_hit <= TOL_COST
            and d_end <= TOL_END
            and d_after <= TOL_END
        )
        # ---- g0 (b): the replay of J5C
        rep_c = json.loads((c / f"J5C_{s}" / "repair.json").read_text())
        R_replay = np.load(c / f"J5C_{s}" / "after.npy").astype(np.float64)
        w_replay = err(R_replay) > SR.WRONG_DEG
        replay_ok = bool(
            rep_c["n_evals"] == SR.G0_J5C_EVALS[s]
            and rep_j5c_stored["n_evals"] == SR.G0_J5C_EVALS[s]
            and rep_c["n_flagged"] == rep_j5c_stored["n_flagged"]
            and rep_c["n_fits"] == rep_j5c_stored["n_fits"]
            and (w_replay == w_j5c).all()
            and j5 + SR.G0_J5C_EVALS[s] == j5c_tot
            and st["total"] == j5
            and int(w_j5.sum()) == st["wrong"]
        )
        g0[s] = dict(
            ok=bool(rerun_ok and replay_ok),
            rerun_ok=rerun_ok,
            replay_ok=replay_ok,
            n_tries=int(n),
            same_ids_and_evals=bool(same_ids),
            max_abs_cost=d_cost,
            max_abs_hit=d_hit,
            max_abs_end=d_end,
            max_abs_after_pass_map=d_after,
            replay_evals=int(rep_c["n_evals"]),
            replay_flagged=int(rep_c["n_flagged"]),
            replay_fits=int(rep_c["n_fits"]),
            replay_wrong=int(w_replay.sum()),
            stored_J5C_wrong=int(w_j5c.sum()),
            replay_map_max_abs_diff=float(np.abs(R_replay - R_j5c).max()),
            j5_wrong=int(w_j5.sum()),
            stored_j5_wrong=st["wrong"],
        )

        # ---- the screen run (J5S)
        sc_run = load_run(f"J5S_{s}")
        rep = sc_run["rep"]
        R_s = sc_run["after"]
        e_s = err(R_s)
        w_s = e_s > SR.WRONG_DEG
        # validity: screen tries equal the trace cut at one evaluation
        tr = np.load(c / f"full_{s}" / "trace.npz")
        scr = np.load(c / f"J5S_{s}" / "tries.npz")
        off, g_n = tr["off"], tr["g_n"]
        first = off[:-1]
        has_first = bool((g_n[first] == 1).all()) and len(first) == len(scr["voxel"])
        ok_ids = bool((scr["voxel"] == full["voxel"]).all() and (scr["cand"] == full["cand"]).all())
        start = d["R_before"][scr["cand"]]
        n_bad_cost = int((tr["g_cost"][first] != scr["cost"]).sum()) if has_first else -1
        n_bad_hit = int((tr["g_hit"][first] != scr["hit"]).sum()) if has_first else -1
        n_bad_end = (
            int((np.abs(tr["g_R"][first] - scr["end"]).reshape(len(first), -1).max(1) > 0).sum())
            if has_first
            else -1
        )
        n_bad_start = int((np.abs(scr["end"] - start).reshape(len(start), -1).max(1) > 0).sum())
        valid[s] = dict(
            ok=bool(
                has_first
                and ok_ids
                and n_bad_cost == 0
                and n_bad_hit == 0
                and n_bad_end == 0
                and n_bad_start == 0
                and bool((scr["evals"] == SR.EVALS_PER_SCREEN).all())
            ),
            n_tries=int(len(scr["voxel"])),
            ids_equal=ok_ids,
            first_gain_is_evaluation_1=has_first,
            n_cost_differ=n_bad_cost,
            n_hit_differ=n_bad_hit,
            n_end_differ_from_trace=n_bad_end,
            n_end_differ_from_start=n_bad_start,
            evals_all_two=bool((scr["evals"] == SR.EVALS_PER_SCREEN).all()),
        )

        # ---- g1 and g2
        total = SR.total_evals(j5, rep["stored_pass_evals"], rep["pass_evals"], rep["n_evals"])
        fixed, broken = SR.fixed_broken(w_j5c, w_s)
        f_j5, b_j5 = SR.fixed_broken(w_j5, w_s)
        f_s8_c, b_s8_c = SR.fixed_broken(w_s8, w_j5c)
        f_s8, b_s8 = SR.fixed_broken(w_s8, w_s)
        net_c, net_s = b_s8_c - f_s8_c, b_s8 - f_s8
        want = (
            SR.P1_PASS_EVALS[s],
            SR.P1_FLAGGED[s],
            SR.P1_REPAIR_EVALS[s],
            SR.P1_TOTAL[s],
            SR.P1_FIXED[s],
            SR.P1_BROKEN[s],
        )
        got = (rep["pass_evals"], rep["n_flagged"], rep["n_evals"], total, fixed, broken)
        g1[s] = dict(ok=bool(got == want), got=list(got), want=list(want))
        med, med_c, med_s8, med_j5 = (float(np.median(e)) for e in (e_s, e_j5c, e_s8, e_j5))
        g = SR.gate_g2(s, total, j5c_tot, fixed, broken, med, med_c, net_s)
        g2.append(g)
        rows[s] = dict(
            sample=s,
            tries=rep["n_tries"],
            pass_evals=rep["pass_evals"],
            stored_pass_evals=rep["stored_pass_evals"],
            flagged=rep["n_flagged"],
            contested=rep["n_contested"],
            extra=rep["n_extra"],
            fits=rep["n_fits"],
            changed=rep["n_changed"],
            repair_evals=rep["n_evals"],
            total=total,
            j5_total=j5,
            j5c_total=j5c_tot,
            s8_total=s8_tot,
            saving_vs_j5c_pct=SR.saving_pct(j5c_tot, total),
            cut_vs_s8_pct=SR.saving_pct(s8_tot, total),
            repair_share_of_j5_pct=100.0 * rep["n_evals"] / j5,
            pass_share_of_j5_pct=100.0 * rep["pass_evals"] / j5,
            wrong_J5=int(w_j5.sum()),
            wrong_J5C=int(w_j5c.sum()),
            wrong_J5S=int(w_s.sum()),
            wrong_S8=int(w_s8.sum()),
            wilson_J5S=wilson(int(w_s.sum()), len(w_s)),
            fixed_vs_J5C=fixed,
            broken_vs_J5C=broken,
            fixed_vs_J5=f_j5,
            broken_vs_J5=b_j5,
            net_J5C_vs_S8=net_c,
            net_J5S_vs_S8=net_s,
            median_J5S=med,
            median_J5C=med_c,
            median_J5=med_j5,
            median_S8=med_s8,
            judge=SR.judge(total, j5c_tot, s8_tot, broken - fixed, net_s, med, med_c, med_s8),
            gate=g,
            wall_pass_s=rep["wall_pass_s"],
            wall_repair_s=rep["wall_repair_s"],
            timing_label=rep["timing_label"],
        )

    ok0 = all(v["ok"] for v in g0.values())
    ok_valid = all(v["ok"] for v in valid.values())
    ok1 = all(v["ok"] for v in g1.values())
    ok2 = all(v["ok"] for v in g2)
    ok = ok0 and ok_valid and ok1 and ok2
    out = dict(
        gate_ok=ok,
        g0_ok=ok0,
        validity_ok=ok_valid,
        g1_ok=ok1,
        g2_ok=ok2,
        g0=g0,
        validity=valid,
        g1=g1,
        g2=g2,
        rows=rows,
        note="wall times are contended, not timing results",
    )
    name = "gate_partial.json" if partial else ("gate.json" if ok else "gate_fail.json")
    (SC.RD / name).write_text(json.dumps(out, indent=1) + "\n")
    lines = [
        "| sample | pass evals | flagged (extra) | repair evals | total | vs J5C | fixed / broken "
        "vs J5C | vs S8 net (J5C, J5S) | median J5C -> J5S | g0 | valid | g1 | g2 |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for s, r in rows.items():
        lines.append(
            f"| {s} | {r['pass_evals']:,} | {r['flagged']:,} ({r['extra']:,}) | "
            f"{r['repair_evals']:,} | {r['total']:,} | {r['saving_vs_j5c_pct']:.1f}% | "
            f"{r['fixed_vs_J5C']} / {r['broken_vs_J5C']} | {r['net_J5C_vs_S8']:+d}, "
            f"{r['net_J5S_vs_S8']:+d} | {r['median_J5C']:.5f} -> {r['median_J5S']:.5f} | "
            f"{'pass' if g0[s]['ok'] else 'FAIL'} | {'pass' if valid[s]['ok'] else 'FAIL'} | "
            f"{'pass' if g1[s]['ok'] else 'FAIL'} | {'pass' if r['gate']['ok'] else 'FAIL'} |"
        )
    (SC.RD / "gate_tables.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print("GATE", "PASSED" if ok else "FAILED", "->", name)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
