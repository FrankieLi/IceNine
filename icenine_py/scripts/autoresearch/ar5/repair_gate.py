"""AR5 step 19, phase A step 5 (gate g1, g2): score the repair of the five stored maps and of J5L.

  uv run python scripts/autoresearch/ar5/repair_gate.py

Reads ``<cache>/<arm>R_<sample>/after.npy`` and ``repair.json`` (``repair_run.py``), the stored
J5 maps after the pass and the truth (scoring only). Per sample: wrong voxels, fixed and broken
against J5 and against S8 (by voxel), median error, repair evaluations as a share of J5's run with
the pass, per-fit evaluations. The gate (``repair_rules.GATE``): g0 (``flags.json``), g1 on mid, g2
on held3, held, grad, big. Writes ``gate.json`` when it passes, else ``gate_fail.json`` (and
exits 1): then stop, no tuning, no other form, no held4. Also writes ``gate_tables.md``.
J5L is a measurement only (against J5, J5R and S8 on mid).
"""

import json
import sys
from typing import Any, Dict

import repair_common as RC  # noqa: E402  (first)
import repair_data as RD
import repair_rules as RR

import numpy as np


def main() -> int:
    RC.RF.check_or_exit()
    y = RD.sym()
    flags = json.loads((RC.RD / "flags.json").read_text())
    rows: Dict[str, Any] = {}
    gate = []
    for s in RR.SAMPLES + ("mid_J5L",):
        arm = "J5L" if s == "mid_J5L" else "J5"
        smp = "mid" if s == "mid_J5L" else s
        d = RD.load(smp, arm)
        st = RD.stored_stats(smp, arm)
        run = RC.RCACHE / f"{arm}R_{smp}"
        rep = json.loads((run / "repair.json").read_text())
        Rn = np.load(run / "after.npy")
        Rt = d["R_t"]
        w_old = RR.wrong_mask(d["R"], Rt, y)
        w_new = RR.wrong_mask(Rn, Rt, y)
        fixed, broken = RR.fixed_broken(w_old, w_new)
        err_old = np.asarray(RR.misorientation_deg(d["R"], Rt, y), float)
        err_new = np.asarray(RR.misorientation_deg(Rn, Rt, y), float)
        share = 100.0 * rep["n_evals"] / st["total"]
        row: Dict[str, Any] = dict(
            sample=smp,
            arm=arm,
            wrong_J5=int(w_old.sum()),
            wrong_J5R=int(w_new.sum()),
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
            median_err_J5R=float(np.median(err_new)),
            timing_label=rep["timing_label"],
        )
        if "R_s8" in d:
            w8 = RR.wrong_mask(d["R_s8"], Rt, y)
            f8_old, b8_old = RR.fixed_broken(w8, w_old)  # J5 against S8: S8 wrong, J5 right ...
            f8, b8 = RR.fixed_broken(w8, w_new)
            row.update(
                wrong_S8=int(w8.sum()),
                J5_vs_S8=dict(fixed=f8_old, broken=b8_old, net=b8_old - f8_old),
                J5R_vs_S8=dict(fixed=f8, broken=b8, net=b8 - f8),
            )
        rows[s] = row
        if s != "mid_J5L":
            g = RR.gate_g12(smp, fixed, broken, share)
            gate.append(g)
            row["gate"] = g
    g0 = bool(flags["g0_ok"])
    ok = g0 and all(g["ok"] for g in gate)
    out = dict(
        gate_ok=ok,
        g0_ok=g0,
        tests=gate,
        rows=rows,
        note="J5L is a measurement only; wall times are contended and are not timing results",
    )
    name = "gate.json" if ok else "gate_fail.json"
    (RC.RD / name).write_text(json.dumps(out, indent=1) + "\n")
    lines = [
        "| sample | flagged | fits | changed | evals | share of J5 | median evals/fit "
        "| wrong J5 -> J5R | fixed / broken vs J5 | gate |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for s, r in rows.items():
        g = r.get("gate")
        lines.append(
            f"| {s} | {r['flagged']} | {r['fits']} | {r['changed']} | {r['repair_evals']:,} | "
            f"{r['share_pct']:.2f}% | {r['median_evals_per_fit']:.0f} | {r['wrong_J5']} -> "
            f"{r['wrong_J5R']} | {r['fixed_vs_J5']} / {r['broken_vs_J5']} | "
            f"{('pass' if g['ok'] else 'FAIL') if g else 'measurement'} |"
        )
    (RC.RD / "gate_tables.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print("GATE", "PASSED" if ok else "FAILED", "->", name)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
