"""AR5 step 8, Phase A step 3: the what-if table of the plan, reproduced from the stored replay.

  uv run python scripts/autoresearch/ar5/order_whatif.py --sample big|grad|mid ... (default: all)

For each of ``big``, ``grad``, ``mid`` (``held`` is not read here) it applies both orders to every
stored seed step of step 2's replay (``mix_rules.mixed_step`` for the interleaved order,
``order_rules.sweep_step`` for the sweep order) and writes ``benchmarks/autoresearch/ar5/order/
whatif.json`` and ``whatif.md`` (the table of the plan). The numbers of the plan are in
``PLAN``; a difference is a stop (the script exits 1 and says which).
"""

import argparse
import json
import sys
from typing import Any, Dict, List

import order_common as OC  # noqa: E402  (first)

import fit_part1_lib as L  # noqa: E402
import fit_rules as FR  # noqa: E402
import mix_rules as MR  # noqa: E402
import numpy as np  # noqa: E402
import order_rules as OR  # noqa: E402

PLAN: Dict[str, Dict[str, Any]] = {
    "big": dict(
        steps=329, x200=938875, sweep=230604, cut=75.4, passes=329, right=329, split=(0, 0),
        differs=37, med_a=0.01232, med_y=0.01232,
    ),
    "grad": dict(
        steps=1279, x200=1563116, sweep=1301274, cut=16.8, passes=1267, right=1260, split=(0, 0),
        differs=7, med_a=0.00898, med_y=0.00906,
    ),
    "mid": dict(
        steps=1611, x200=5643575, sweep=3268675, cut=42.1, passes=1600, right=1599, split=(0, 0),
        differs=188, med_a=0.01135, med_y=0.01121,
    ),
}  # fmt: skip


def one(sample: str) -> Dict[str, Any]:
    rows = L.load(sample)
    a_steps: List[Dict[str, Any]] = []
    y_steps: List[Dict[str, Any]] = []
    for r in rows:
        c, m = r["arms"].get("C200"), r["arms"].get("M")
        a_steps.append(MR.mixed_step(r["costed"], c, m))
        y_steps.append(OR.sweep_step(r["costed"], c, m))
    ea = int(sum(s["evals"] for s in a_steps))
    ey = int(sum(s["evals"] for s in y_steps))
    pa = [s["passed"] for s in a_steps]
    py = [s["passed"] for s in y_steps]
    ra = [s["right"] for s in a_steps]
    ry = [s["right"] for s in y_steps]
    sp = FR.paired_split(pa, py)
    sr = FR.paired_split(ra, ry)
    # "taken end differs": the end taken by the two orders is not the same end (same candidate
    # and same end kind C/M, or both without a taken end)
    differs = 0
    for r, a, y in zip(rows, a_steps, y_steps):
        if a["taken"] is None and y["taken"] is None:
            continue
        ka = (a["taken"], _kind(r, a))
        ky = (y["taken"], y["end"])
        differs += int(ka != ky)
    ea_err = np.array([s["err"] for s in a_steps])
    ey_err = np.array([s["err"] for s in y_steps])
    return dict(
        sample=sample,
        steps=len(rows),
        x200=ea,
        sweep=ey,
        cut=round(100.0 * (1 - ey / ea), 1),
        passes_a=int(sum(pa)),
        passes_y=int(sum(py)),
        right_a=int(sum(ra)),
        right_y=int(sum(ry)),
        passes_split=[sp["only_a"], sp["only_b"]],
        right_split=[sr["only_a"], sr["only_b"]],
        differs=differs,
        median_a=float(np.median(ea_err[np.isfinite(ea_err)])),
        median_y=float(np.median(ey_err[np.isfinite(ey_err)])),
        missing_y=int(sum(s["missing"] for s in y_steps)),
        missing_a=int(sum(s["missing"] for s in a_steps)),
    )


def _kind(row: Dict[str, Any], a: Dict[str, Any]) -> str:
    """The end kind the interleaved order took ('C' or 'M'); ``mixed_step`` does not store it."""
    p = a["taken"]
    if p is None:
        return ""
    c, m = row["arms"]["C200"], row["arms"]["M"]
    if c["fin_hit"][p] < FR.TAU and m["fin_cost"][p] < c["fin_cost"][p]:
        return "M"
    return "C"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", nargs="*", default=["big", "grad", "mid"])
    a = ap.parse_args()
    assert "held" not in a.sample, "held is not read by the what-if"
    res = {s: one(s) for s in a.sample}
    bad: List[str] = []
    for s, r in res.items():
        p = PLAN[s]
        chk = dict(
            steps=r["steps"] == p["steps"],
            x200=r["x200"] == p["x200"],
            sweep=r["sweep"] == p["sweep"],
            cut=r["cut"] == p["cut"],
            passes=r["passes_a"] == p["passes"] and r["passes_y"] == p["passes"],
            right=r["right_a"] == p["right"] and r["right_y"] == p["right"],
            split=tuple(r["passes_split"]) == p["split"] and tuple(r["right_split"]) == p["split"],
            differs=r["differs"] == p["differs"],
            med=round(r["median_a"], 5) == p["med_a"] and round(r["median_y"], 5) == p["med_y"],
        )
        r["plan_match"] = chk
        bad += [f"{s}:{k}" for k, v in chk.items() if not v]
    out = dict(samples=res, mismatches=bad)
    OC.ORDER.mkdir(parents=True, exist_ok=True)
    (OC.ORDER / "whatif.json").write_text(json.dumps(out, indent=1) + "\n")
    lines = [
        "| sample | seed steps | `X200` | sweep order | cut | passes, both | right, both "
        "| paired split | taken end differs | median error |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for s, r in res.items():
        lines.append(
            f"| {s} | {r['steps']:,} | {r['x200']:,} | {r['sweep']:,} | {r['cut']}% "
            f"| {r['passes_y']:,} | {r['right_y']:,} "
            f"| {r['passes_split'][0]} / {r['passes_split'][1]} | {r['differs']} "
            f"| {r['median_a']:.5f} to {r['median_y']:.5f} |"
        )
    (OC.ORDER / "whatif.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print("mismatches against the plan:", bad or "none")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
