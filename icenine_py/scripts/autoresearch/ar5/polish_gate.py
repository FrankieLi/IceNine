"""AR5 step 14, phase A steps 3 and 5: the gates (pure rules; the CLI reads the probe files).

  uv run python scripts/autoresearch/ar5/polish_gate.py --gate 1     # writes polish/gate1.json
  uv run python scripts/autoresearch/ar5/polish_gate.py --gate 2     # writes polish/gate2.json

Gate 1, disclosed probe (``mid``, ``grad``, ``held``: 28 cases), for each form: (g0) prediction 1
(every stored rescue end re-scores to its stored cost and hit ratio, 1e-9); (g1) at least 13 of 28
accepted; (g2) no accepted end wrong (over 1 degree) where the grid end is right; (g3) median error
of accepted ends at most 0.05 degree, none over 0.3; (g4) the grid evaluations of accepted cases
minus all ``G2`` evaluations on the 28 at least 40% of 9,595,587. If both forms pass take the
larger (g4); a tie goes to ``G2b``. If neither passes: stop, no tuning, no third form.
Gate 2, ``held2`` probe of the chosen form (21 cases; after the freeze commit): (t1) at least 7 of
21 accepted; (t2) as (g2); (t3) net saving at least 25% of 7,252,187. Else stop.
The counts are checked against the stored ones (probe files carry the grid evaluations of each
case). Truth errors are for scoring only.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
PO = ICE / "benchmarks" / "autoresearch" / "ar5" / "polish"

GRID_G1 = 9_595_587
GRID_G2 = 7_252_187
N_G1 = 28
N_G2 = 21
G1_SAMPLES = ("mid", "grad", "held")
FORMS = ("G2a", "G2b")


def metrics(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Counts of a probe's rows (one per failed rescue)."""
    acc = [r for r in rows if r["accepted"]]
    errs = [float(r["taken_err"]) for r in acc]
    saved = int(sum(int(r["saved"]) for r in rows))
    spent = int(sum(int(r["probe_evals"]) for r in rows))
    bad = [r for r in acc if r["wrong_accepted"] and r["grid_end_err"] <= 1.0]
    return dict(
        n=len(rows),
        accepted=len(acc),
        wrong_accepted=int(sum(1 for r in acc if r["wrong_accepted"])),
        wrong_where_grid_right=len(bad),
        median_err=float(np.median(errs)) if errs else None,
        max_err=float(max(errs)) if errs else None,
        saved_grid_evals=saved,
        probe_evals=spent,
        net_saving=saved - spent,
        grid_evals=int(sum(int(r["grid_evals"]) for r in rows)),
        accepted_angle_to_grid_end=[float(r["taken_angle_to_grid_end"]) for r in acc],
    )


def gate1_form(m: Dict[str, Any], pred1_ok: bool) -> Dict[str, Any]:
    g0 = bool(pred1_ok)
    g1 = m["accepted"] >= 13
    g2 = m["wrong_where_grid_right"] == 0
    g3 = m["accepted"] > 0 and m["median_err"] <= 0.05 and m["max_err"] <= 0.3
    g4 = m["net_saving"] * 10 >= 4 * GRID_G1  # at least 40% of 9,595,587
    return dict(
        g0=g0,
        g1=bool(g1),
        g2=bool(g2),
        g3=bool(g3),
        g4=bool(g4),
        passes=bool(g0 and g1 and g2 and g3 and g4),
    )


def choose(res: Dict[str, Dict[str, Any]]) -> Optional[str]:
    """The form that passes with the larger net saving (g4); a tie goes to G2b; None if none."""
    ok = [f for f in FORMS if res[f]["gate"]["passes"]]
    if not ok:
        return None
    if len(ok) == 1:
        return ok[0]
    a, b = res["G2a"]["metrics"]["net_saving"], res["G2b"]["metrics"]["net_saving"]
    return "G2a" if a > b else "G2b"


def gate2(m: Dict[str, Any]) -> Dict[str, Any]:
    t1 = m["accepted"] >= 7
    t2 = m["wrong_where_grid_right"] == 0
    t3 = m["net_saving"] * 100 >= 25 * GRID_G2
    return dict(t1=bool(t1), t2=bool(t2), t3=bool(t3), passes=bool(t1 and t2 and t3))


def _rows(form: str, samples: Any) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for s in samples:
        out += json.loads((PO / f"probe_{form}_{s}.json").read_text())["rows"]
    return out


def main() -> None:
    which = sys.argv[sys.argv.index("--gate") + 1]
    if which == "1":
        pred1 = all(json.loads((PO / f"pred1_{s}.json").read_text())["ok"] for s in G1_SAMPLES)
        res: Dict[str, Dict[str, Any]] = {}
        for f in FORMS:
            rows = _rows(f, G1_SAMPLES)
            m = metrics(rows)
            assert m["n"] == N_G1 and m["grid_evals"] == GRID_G1, (m["n"], m["grid_evals"])
            res[f] = dict(metrics=m, gate=gate1_form(m, pred1))
        ch = choose(res)
        out = dict(gate=1, pred1_ok=pred1, forms=res, chosen=ch, passes=ch is not None)
        (PO / "gate1.json").write_text(json.dumps(out, indent=1) + "\n")
        for f in FORMS:
            m = res[f]["metrics"]
            print(
                f,
                {
                    k: m[k]
                    for k in (
                        "accepted",
                        "wrong_accepted",
                        "median_err",
                        "max_err",
                        "net_saving",
                        "probe_evals",
                    )
                },
                res[f]["gate"],
            )
        print("chosen:", ch)
        sys.exit(0 if ch else 3)
    g1 = json.loads((PO / "gate1.json").read_text())
    form = g1["chosen"]
    rows = _rows(form, ("held2",))
    m = metrics(rows)
    assert m["n"] == N_G2 and m["grid_evals"] == GRID_G2, (m["n"], m["grid_evals"])
    p1 = json.loads((PO / "pred1_held2.json").read_text())["ok"]
    g = gate2(m)
    g["pred1"] = bool(p1)
    g["passes"] = bool(g["passes"] and p1)
    out = dict(gate=2, form=form, metrics=m, gate_result=g, passes=g["passes"])
    (PO / "gate2.json").write_text(json.dumps(out, indent=1) + "\n")
    print(
        form,
        {
            k: m[k]
            for k in (
                "accepted",
                "wrong_accepted",
                "median_err",
                "max_err",
                "net_saving",
                "probe_evals",
            )
        },
        g,
    )
    sys.exit(0 if g["passes"] else 3)


if __name__ == "__main__":
    main()
