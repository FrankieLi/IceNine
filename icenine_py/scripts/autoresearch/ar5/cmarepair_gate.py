"""AR5 step 20, phase A step 4 (gates g1 and g2): score the CMA repair of the five stored maps and
of J5L.

  uv run python scripts/autoresearch/ar5/cmarepair_gate.py

Reads ``<cache>/<arm>C_<sample>/after.npy`` and ``repair.json`` (``cmarepair_run.py``), the stored
J5 maps after the pass and the truth (scoring only). Per sample: wrong voxels, fixed and broken
against J5 and against S8 (by voxel), median error, repair evaluations as a share of J5's run with
the pass, per-fit evaluations. The gate (``cmarepair_rules.GATE``): g0 (``flags.json``), g1
(prediction 1, exactly: flagged and fits, repair evaluations, fixed / broken, and J5L), g2 (per
disclosed sample: cost at most 1.5% of J5, 2.0% on grad; fixed minus broken against J5 at least 3;
on mid broken minus fixed against S8 at most +1). Writes ``gate.json`` when it passes, else
``gate_fail.json`` (and exits 1): then stop, no tuning, no other form, no held4. Also writes
``gate_tables.md``. J5L is a measurement only for g2 (its count is part of g1).
"""

import json
import sys
from typing import Any, Dict

import cmarepair_common as RC  # noqa: E402  (first)
import cmarepair_data as RD  # noqa: E402
import cmarepair_rules as CR  # noqa: E402
import numpy as np  # noqa: E402


def main() -> int:
    RC.CF.check_or_exit()
    y = RD.sym()
    flags = json.loads((RC.RD / "flags.json").read_text())
    rows: Dict[str, Any] = {}
    g1: Dict[str, Any] = {}
    g2 = []
    for s in CR.SAMPLES + ("mid_J5L",):
        arm = "J5L" if s == "mid_J5L" else "J5"
        smp = "mid" if s == "mid_J5L" else s
        d = RD.load(smp, arm)
        st = RD.stored_stats(smp, arm)
        run = RC.RCACHE / f"{arm}C_{smp}"
        rep = json.loads((run / "repair.json").read_text())
        Rn = np.load(run / "after.npy")
        Rt = d["R_t"]
        w_old = CR.wrong_mask(d["R"], Rt, y)
        w_new = CR.wrong_mask(Rn, Rt, y)
        fixed, broken = CR.fixed_broken(w_old, w_new)
        err_old = np.asarray(CR.misorientation_deg(d["R"], Rt, y), float)
        err_new = np.asarray(CR.misorientation_deg(Rn, Rt, y), float)
        share = 100.0 * rep["n_evals"] / st["total"]
        row: Dict[str, Any] = dict(
            sample=smp,
            arm=arm,
            wrong_J5=int(w_old.sum()),
            wrong_J5C=int(w_new.sum()),
            fixed_vs_J5=fixed,
            broken_vs_J5=broken,
            flagged=rep["n_flagged"],
            fits=rep["n_fits"],
            changed=rep["n_changed"],
            repair_evals=rep["n_evals"],
            j5_total=st["total"],
            share_pct=share,
            median_evals_per_fit=rep["median_evals_per_fit"],
            median_err_J5=float(np.median(err_old)),
            median_err_J5C=float(np.median(err_new)),
            timing_label=rep["timing_label"],
        )
        if "R_s8" in d:
            w8 = CR.wrong_mask(d["R_s8"], Rt, y)
            f8_old, b8_old = CR.fixed_broken(w8, w_old)
            f8, b8 = CR.fixed_broken(w8, w_new)
            row.update(
                wrong_S8=int(w8.sum()),
                J5_vs_S8=dict(fixed=f8_old, broken=b8_old, net=b8_old - f8_old),
                J5C_vs_S8=dict(fixed=f8, broken=b8, net=b8 - f8),
            )
        rows[s] = row
        if s == "mid_J5L":
            want = CR.P1_J5L
            g1[s] = dict(
                ok=bool(
                    (rep["n_flagged"], rep["n_fits"], rep["n_evals"], fixed, broken)
                    == (want["flagged"], want["fits"], want["evals"], want["fixed"], want["broken"])
                ),
                got=[rep["n_flagged"], rep["n_fits"], rep["n_evals"], fixed, broken],
                want=[want["flagged"], want["fits"], want["evals"], want["fixed"], want["broken"]],
            )
        else:
            want_t = (CR.FLAGGED[s], CR.FITS[s], CR.P1_EVALS[s], CR.P1_FIXED[s], CR.P1_BROKEN[s])
            got_t = (rep["n_flagged"], rep["n_fits"], rep["n_evals"], fixed, broken)
            g1[s] = dict(ok=bool(got_t == want_t), got=list(got_t), want=list(want_t))
            g = CR.gate_g2(smp, fixed, broken, share)
            if smp == "mid":
                g["s8_net"] = row["J5C_vs_S8"]["net"]
                g["ok_s8"] = bool(g["s8_net"] <= int(CR.GATE["mid_net_vs_s8_max"]))
                g["ok"] = bool(g["ok"] and g["ok_s8"])
            g2.append(g)
            row["gate"] = g
    g0 = bool(flags["g0_ok"])
    g1_ok = all(v["ok"] for v in g1.values())
    g2_ok = all(g["ok"] for g in g2)
    ok = g0 and g1_ok and g2_ok
    out = dict(
        gate_ok=ok,
        g0_ok=g0,
        g1_ok=g1_ok,
        g2_ok=g2_ok,
        g1=g1,
        g2=g2,
        rows=rows,
        note="J5L is a measurement only for g2; wall times are contended, not timing results",
    )
    name = "gate.json" if ok else "gate_fail.json"
    (RC.RD / name).write_text(json.dumps(out, indent=1) + "\n")
    lines = [
        "| sample | flagged | fits | changed | evals | share of J5 | median evals/fit "
        "| wrong J5 -> J5C | fixed / broken vs J5 | vs S8 net (J5, J5C) | g1 | g2 |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for s, r in rows.items():
        g = r.get("gate")
        s8 = f"{r['J5_vs_S8']['net']:+d}, {r['J5C_vs_S8']['net']:+d}" if "J5C_vs_S8" in r else ""
        lines.append(
            f"| {s} | {r['flagged']} | {r['fits']} | {r['changed']} | {r['repair_evals']:,} | "
            f"{r['share_pct']:.2f}% | {r['median_evals_per_fit']:.0f} | {r['wrong_J5']} -> "
            f"{r['wrong_J5C']} | {r['fixed_vs_J5']} / {r['broken_vs_J5']} | {s8} | "
            f"{'pass' if g1[s]['ok'] else 'FAIL'} | "
            f"{('pass' if g['ok'] else 'FAIL') if g else 'measurement'} |"
        )
    (RC.RD / "gate_tables.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print("GATE", "PASSED" if ok else "FAILED", "->", name)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
