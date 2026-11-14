"""AR5 step 9, phase A step 4: the gate, from the two replay files (disclosed samples only).

  uv run python scripts/autoresearch/ar5/rescue_gate.py

Reads ``rescue/replay_grad.json`` and ``rescue/replay_mid.json`` and writes ``rescue/gate.json``:
rescued of 23 (Wilson 95%), wrong rescued ends, tries per case, evaluations (all tries), the
errors of the rescued ends, and the plan's predictions 2. Gate (plan, "Decision rules"): rescued
in at least 12 of 23 and no wrong rescued end; else stop, tune nothing. Exit 1 when it fails.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent.parent / "common"))
RESCUE = HERE.parents[2] / "benchmarks" / "autoresearch" / "ar5" / "rescue"

import stats as ST  # noqa: E402


def gate(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    n = len(rows)
    resc = [r for r in rows if r["rescued"]]
    k = len(resc)
    lo, hi = ST.wilson(k, n)
    errs = [r["taken_err"] for r in resc]
    tries = [r["n_tries"] for r in rows]
    ev = [r["evals"] for r in rows]
    ok = k >= 12 and not any(r["wrong_rescued"] for r in rows)
    return dict(
        cases=n,
        rescued=k,
        wilson95=[lo, hi],
        wrong_rescued=int(sum(1 for r in rows if r["wrong_rescued"])),
        tries_mean=float(np.mean(tries)),
        tries_total=int(sum(tries)),
        tries_min_max=[int(min(tries)), int(max(tries))],
        evals_total=int(sum(ev)),
        evals_per_case_mean=float(np.mean(ev)),
        rescued_err_median=float(np.median(errs)) if errs else None,
        rescued_err_max=float(max(errs)) if errs else None,
        rescued_hit_min=float(min(r["taken_hit"] for r in resc)) if resc else None,
        grid_end_reached_095=int(sum(1 for r in rows if r["grid_end_reached_095"])),
        rescued_where_grid_end_below_095=int(sum(1 for r in resc if not r["grid_end_reached_095"])),
        prediction2=dict(
            rescued_12_to_21=bool(12 <= k <= 21),
            every_rescued_right=bool(all(not r["wrong_rescued"] for r in resc)),
            median_err_le_0_05=bool(errs and float(np.median(errs)) <= 0.05),
            tries_mean_le_2_5=bool(np.mean(tries) <= 2.5),
        ),
        gate_passed=bool(ok),
    )


def main() -> None:
    rows: List[Dict[str, Any]] = []
    per: Dict[str, Any] = {}
    for s in ("grad", "mid"):
        d = json.loads((RESCUE / f"replay_{s}.json").read_text())
        per[s] = gate(d["rows"])
        rows += d["rows"]
    out = dict(all=gate(rows), per_sample=per)
    (RESCUE / "gate.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out["all"], indent=1))
    sys.exit(0 if out["all"]["gate_passed"] else 1)


if __name__ == "__main__":
    main()
