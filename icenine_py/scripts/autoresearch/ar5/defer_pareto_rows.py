"""AR5 step 22, phase B: the pareto rows of this step, in the form of
``ar5/screen/pareto_rows.json``: ``J5SD`` on the eight samples and, on the clean test sample
``held6``, ``S8``, ``J5`` and ``J5S`` as well.

  uv run python scripts/autoresearch/ar5/defer_pareto_rows.py

Reads ``ar5/defer/summary.json`` only (no sample is opened). ``in_sample`` is true for the seven
disclosed samples (the rule was chosen on them) and false for ``held6``. Writes
``ar5/defer/pareto_rows.json``.
"""

import json
from pathlib import Path
from typing import Any, Dict, List

HERE = Path(__file__).resolve().parent
RD = HERE.parents[2] / "benchmarks" / "autoresearch" / "ar5" / "defer"
N = 24570
ORDER = ("held6", "held5", "held4", "held3", "held", "mid", "grad", "big")


def row(sample: str, arm: str, total: int, wrong: Dict[str, Any]) -> Dict[str, Any]:
    return dict(
        sample=sample,
        in_sample=sample != "held6",
        arm=arm,
        evals_total=int(total),
        evals_per_voxel=round(total / N, 2),
        wrong=int(wrong["n"]),
        n=N,
        fraction_right=round(1 - wrong["n"] / N, 6),
        wrong_wilson95=[round(x, 6) for x in wrong["wilson"]],
    )


def main() -> None:
    d = json.loads((RD / "summary.json").read_text())["rows"]
    rows: List[Dict[str, Any]] = []
    for s in ORDER:
        r = d[s]
        arms = ("s8", "j5", "j5s", "j5sd") if s == "held6" else ("j5sd",)
        for a in arms:
            rows.append(row(s, a.upper(), r["totals"][a], r["wrong"][a]))
    (RD / "pareto_rows.json").write_text(json.dumps(rows, indent=1) + "\n")
    print(len(rows), "rows")


if __name__ == "__main__":
    main()
