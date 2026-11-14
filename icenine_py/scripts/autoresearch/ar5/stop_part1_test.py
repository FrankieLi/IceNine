"""AR5 step 1, part 1c: every setting of the try test on the three samples; the pass rule on mid.

  RESEED_SAMPLE=mid uv run python scripts/autoresearch/ar5/stop_part1_starts.py   # start scores
  uv run python scripts/autoresearch/ar5/stop_part1_test.py

Reads ``part1_pick.json`` (committed before ``mid`` was read) and ``starts_<sample>.npz``.
Nothing is changed after the pick: the settings are the plan's, P3 uses the picked ``d``.
The pass rule (plan): on ``mid`` the picked setting cuts the evaluations of the pass by at least
half and loses at most 5% of its net gain. Also writes the prediction 1 check: some setting
removes at least 70% of the pass's evaluations on the pick samples, loses at most 5% of its net
gain there, and passes on mid.
"""

import json
from typing import Any, Dict

import stop_part1_lib as P
import stop_rules as SR


def main() -> None:
    pick = json.loads((P.AR5 / "part1_pick.json").read_text())
    d_pick = pick["p1_pick_d"]
    data = {s: P.load(s) for s in ("big", "grad", "mid")}
    table: Dict[str, Dict[str, Any]] = {}
    for name, kind, d, m in P.settings(d_pick):
        table[name] = {s: P.evaluate(data[s], kind, d, m) for s in data}
    full = {s: table["full"][s] for s in data}
    pooled_full = full["big"]["net_gain"] + full["grad"]["net_gain"]
    evals_full = full["big"]["evals"] + full["grad"]["evals"]
    out: Dict[str, Any] = dict(table=table, d_pick=d_pick, pick=pick["pick"])
    # the pass rule for the picked setting, on mid
    if pick["pick"] is not None:
        r = table[pick["pick"]["name"]]["mid"]
        f = full["mid"]
        out["pass_rule_mid"] = dict(
            setting=pick["pick"]["name"],
            evals=r["evals"],
            evals_full=f["evals"],
            cut=1 - r["evals"] / f["evals"],
            net_gain=r["net_gain"],
            net_gain_full=f["net_gain"],
            passes=SR.try_test_passes(r["evals"], f["evals"], r["net_gain"], f["net_gain"]),
        )
    else:
        out["pass_rule_mid"] = dict(setting=None, passes=False)
    # prediction 1: any setting with >= 70% cut and <= 5% loss on the pick samples, passing on mid
    pred = []
    for name, per in table.items():
        if name == "full":
            continue
        ev = per["big"]["evals"] + per["grad"]["evals"]
        gain = per["big"]["net_gain"] + per["grad"]["net_gain"]
        cut = 1 - ev / evals_full
        ok_pick = cut >= 0.70 and gain >= (1 - SR.MAX_LOSS) * pooled_full - SR.TOL
        mid = per["mid"]
        ok_mid = SR.try_test_passes(
            mid["evals"], full["mid"]["evals"], mid["net_gain"], full["mid"]["net_gain"]
        )
        pred.append(
            dict(setting=name, cut_pick=cut, gain_pick=gain, ok_pick=ok_pick, ok_mid=ok_mid)
        )
    out["prediction1"] = dict(
        rows=pred,
        holds=any(p["ok_pick"] and p["ok_mid"] for p in pred),
        holds_on_pick_samples=any(p["ok_pick"] for p in pred),
    )
    P.AR5.joinpath("part1_results.json").write_text(json.dumps(out, indent=1) + "\n")
    for name, per in table.items():
        print(
            f"{name:<16}",
            " | ".join(
                f"{s}: ev {per[s]['evals']:>8} kept {per[s]['tries_kept']:>6} "
                f"wrong {per[s]['wrong_after']:>3} fix {per[s]['fixed']} brk {per[s]['broken']}"
                for s in data
            ),
        )
    print("pass rule on mid:", out["pass_rule_mid"])
    print(
        "prediction 1 holds:",
        out["prediction1"]["holds"],
        "on pick samples:",
        out["prediction1"]["holds_on_pick_samples"],
    )


if __name__ == "__main__":
    main()
