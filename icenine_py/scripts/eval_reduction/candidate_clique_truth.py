"""Task C step 1 addendum: is the truth-nearest candidate of a clique also its best-scored one?

For every stored step-1 run and level: the discrete-stage candidates within 3 deg (and 1 deg) of the
truth, and how many of them are NOT the best-scored candidate of their clique (the one a
"best 1 per clique" filter would keep). Also, per run and level, whether the candidate nearest to
the truth (within 3 deg) survives such a filter. Counts only.

  uv run python scripts/eval_reduction/candidate_clique_truth.py
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICE / "scripts" / "findoptimal_robustness"))
sys.path.insert(0, str(ICE / "scripts" / "common"))
RAW = ICE / "benchmarks" / "eval_reduction" / "C" / "cache" / "step1_raw"
OUT = ICE / "benchmarks" / "eval_reduction" / "C" / "step1"


def main() -> None:
    import doc_tables as DT
    from candidate_overlap import load
    from csl import reduced_misorientation_deg

    rows: List[Dict[str, Any]] = []
    summ: Dict[str, Any] = {}
    for which in ("full", "e0"):
        for rule in ("cpp", "lenient"):
            acc: Dict[int, Dict[str, int]] = {}
            for fn in sorted(RAW.glob(f"{which}_{rule}_v*.npz")):
                run = load(fn)
                Rt = np.array(run["meta"]["R_true"])
                for l, L in run["levels"].items():
                    R = L["R_disc"]
                    cl = L["clique"].astype(int)
                    d = reduced_misorientation_deg(R, np.broadcast_to(Rt, R.shape).copy())
                    first = np.zeros(len(R), bool)  # best-scored of its clique (list = rank order)
                    seen: set = set()
                    for i, c in enumerate(cl):
                        if c not in seen:
                            first[i] = True
                            seen.add(c)
                    a = acc.setdefault(
                        l,
                        dict(runs=0, near3=0, near3_not_first=0, near1=0, near1_not_first=0,
                             runs_with_near3=0, nearest_not_first=0, nearest_in_multi=0,
                             nearest_not_first_multi=0),
                    )  # fmt: skip
                    a["runs"] += 1
                    n3, n1 = d < 3.0, d < 1.0
                    a["near3"] += int(n3.sum())
                    a["near3_not_first"] += int((n3 & ~first).sum())
                    a["near1"] += int(n1.sum())
                    a["near1_not_first"] += int((n1 & ~first).sum())
                    if n3.any():
                        a["runs_with_near3"] += 1
                        j = int(d.argmin())
                        a["nearest_not_first"] += int(not first[j])
                        multi = int((cl == cl[j]).sum()) >= 2
                        a["nearest_in_multi"] += int(multi)
                        a["nearest_not_first_multi"] += int(multi and not first[j])
            summ[f"{which}_{rule}"] = acc
            for l, a in sorted(acc.items()):
                rows.append(dict(set=which, rule=rule, level=l, **a))
    (OUT / "clique_truth.json").write_text(json.dumps(summ, indent=1) + "\n")
    cols = ["set", "rule", "level", "runs", "near3", "near3_not_first", "near1", "near1_not_first",
            "runs_with_near3", "nearest_not_first", "nearest_in_multi",
            "nearest_not_first_multi"]  # fmt: skip
    DT.write_tables(
        OUT / "tables_clique_truth.md", {"c_clique_truth": DT.markdown_table(rows, cols)}
    )
    print(DT.markdown_table(rows, cols))


if __name__ == "__main__":
    main()
