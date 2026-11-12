"""AR5 step 8, Phase A step 1: the stored ``S11`` runs of step 3 (big, grad, mid) are found and
their totals agree with step 3's tracked tables. Reads ``run.json`` of each stored run (kept
worktree) and ``benchmarks/autoresearch/ar5/mix/part2_S11_<sample>.json``. Writes
``order/stored.json``.
``held`` is not read."""

import json
import sys
from typing import Any, Dict

import order_common as OC  # noqa: E402  (first)

STAGES = ("seed_step", "grid_fallback", "neighbour", "revisit", "before_pass", "pass_", "total")


def main() -> None:
    out: Dict[str, Any] = {}
    bad = []
    for s in ("big", "grad", "mid"):
        run = json.loads((OC.STORED_S11[s] / "run.json").read_text())
        tab = json.loads(
            (OC.ICE / "benchmarks/autoresearch/ar5/mix" / f"part2_S11_{s}.json").read_text()
        )
        after = int(tab["stages"]["before_pass"])
        total_with_pass = int(tab["stages"]["total"])
        passj = json.loads((OC.STORED_S11_PASS[s] / "D4_S11" / "pass.json").read_text())
        pass_evals = int(passj["evals"] + passj["score_evals"]) if "score_evals" in passj else -1
        ok = int(run["evals_total"]) == after
        out[s] = dict(
            run_evals_total=int(run["evals_total"]),
            table_before_pass=after,
            table_total_with_pass=total_with_pass,
            table_pass=int(tab["stages"]["pass_"]),
            pass_json_evals=int(passj["evals"]),
            run_matches_table=ok,
            seed_steps=int(tab["seed_steps"]),
            grid_fallbacks=int(tab["grid_fallbacks"]),
            pass_total_matches=int(passj["evals"]) == int(tab["stages"]["pass_"]),
        )
        del pass_evals
        if not (ok and out[s]["pass_total_matches"]):
            bad.append(s)
    out["mismatch"] = bad
    (OC.ORDER / "stored.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
