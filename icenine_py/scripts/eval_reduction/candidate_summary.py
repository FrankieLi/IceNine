"""Task C step 3 summary: the six-voxel ledger runs and the E0 arms side by side, the arm chosen for
Task G, and the doc tables.

  uv run python scripts/eval_reduction/candidate_summary.py

Reads ``benchmarks/eval_reduction/C/ledger/runs/<arm>_v<voxel>.json`` (eval_ledger.py runs with
``--opt candidate_dedup=<spec>``) and ``step3/summary_e0.json`` (candidate_arms.py summary); writes
``step3/summary.json``, ``step3/tables.md`` and ``step3/arms.png``.

Choice rule (plan, Task C decision point): the arm with the fewest evaluations among those with no
broken voxel on the E0 set (paired against the baseline) and no more wrong voxels than the baseline
on the six full-sample voxels; the baseline when there is none. Evaluations are compared on the six
full-sample voxels (E0 voxels have few candidates, so little overlap) and quoted for E0 too.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(ICE / "scripts" / "common"))
BASE = ICE / "benchmarks" / "eval_reduction" / "C"
SIX = [15901, 17034, 18558, 19889, 22867, 5242]
WRONG_DEG = 1.0  # as common.WRONG_DEG (E0)


def arm_order() -> List[Dict[str, str]]:
    out = []
    for line in (BASE / "arms.txt").read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            a, s = line.split()
            out.append(dict(arm=a, spec=s))
    return out


def main() -> None:
    import doc_tables as DT
    import stats as S

    runs: Dict[str, Dict[int, Dict[str, Any]]] = {}
    for fn in sorted((BASE / "ledger" / "runs").glob("*_v*.json")):
        r = json.loads(fn.read_text())
        runs.setdefault(r["tag"], {})[r["voxel"]] = r
    e0 = json.loads((BASE / "step3" / "summary_e0.json").read_text())
    base6 = runs["base"]
    summ: Dict[str, Any] = {}
    rows_six, rows_all, rows_stage = [], [], []
    for a in arm_order():
        arm = a["arm"]
        R = runs.get(arm, {})
        if len(R) < len(SIX):
            continue
        ev = np.array([R[v]["evals_total"] for v in SIX], dtype=float)
        evb = np.array([base6[v]["evals_total"] for v in SIX], dtype=float)
        err = np.array([R[v]["err_deg"] for v in SIX])
        errb = np.array([base6[v]["err_deg"] for v in SIX])
        wrong, wrongb = err > WRONG_DEG, errb > WRONG_DEG
        broken, fixed = S.paired_discordant(wrong, wrongb)
        lo, hi = S.wilson(int(wrong.sum()), len(SIX))
        stages = {k: int(sum(R[v]["totals"][k] for v in SIX)) for k in R[SIX[0]]["totals"]}
        labels = [R[v]["label"] for v in SIX]
        wall = float(sum(R[v]["wall_s"] for v in SIX))
        s6: Dict[str, Any] = dict(
            spec=a["spec"], evals_total=int(ev.sum()), evals_ratio=float(ev.sum() / evb.sum()),
            evals_per_voxel=[int(x) for x in ev], wrong=int(wrong.sum()), wilson=[lo, hi],
            broken=broken, fixed=fixed, err_per_voxel=[float(x) for x in err],
            max_err=float(err.max()), median_err=float(np.median(err)), stages=stages,
            wall_s=wall, n_single_worker_quiet=int(sum(1 for x in labels if x == "single-worker")),
        )  # fmt: skip
        d: Dict[str, Any] = dict(six=s6)
        if arm in e0:
            d["e0"] = e0[arm]
        summ[arm] = d
        rows_six.append(
            dict(
                arm=arm, spec=a["spec"], evals=int(ev.sum()),
                vs_base=f"{100 * ev.sum() / evb.sum():.1f}%",
                wrong=f"{int(wrong.sum())}/6", broken=broken, fixed=fixed,
                median_err=float(np.median(err)), max_err=float(err.max()),
                wall_min=wall / 60, quiet=f"{s6['n_single_worker_quiet']}/6",
            )
        )  # fmt: skip
        rows_stage.append(dict(arm=arm, **{k: stages[k] for k in stages}))
        if arm in e0:
            x = e0[arm]
            rows_all.append(
                dict(
                    arm=arm, spec=a["spec"],
                    six_evals=int(ev.sum()), six_vs_base=f"{100 * ev.sum() / evb.sum():.1f}%",
                    six_wrong=f"{int(wrong.sum())}/6",
                    e0_mean_evals=x["mean_evals"], e0_vs_base=f"{100 * x['evals_ratio']:.1f}%",
                    e0_wrong=f"{x['wrong']}/{x['n']}", e0_broken=x["broken"], e0_fixed=x["fixed"],
                    e0_both_wrong=x["both_wrong"],
                )
            )  # fmt: skip
    # choice
    eligible = [
        a for a, d in summ.items()
        if a != "base" and "e0" in d and d["e0"]["broken"] == 0
        and d["six"]["wrong"] <= summ["base"]["six"]["wrong"]
    ]  # fmt: skip
    choice = min(eligible, key=lambda a: summ[a]["six"]["evals_total"]) if eligible else "base"
    summ["_choice"] = dict(
        arm=choice, eligible=eligible,
        rule="fewest six-voxel evaluations among arms with no broken E0 voxel and no more wrong "
        "voxels than the baseline on the six; the baseline if none",
    )  # fmt: skip
    (BASE / "step3").mkdir(exist_ok=True)
    (BASE / "step3" / "summary.json").write_text(json.dumps(summ, indent=1, default=float) + "\n")
    tabs = {
        "c_six": DT.markdown_table(
            rows_six,
            ["arm", "spec", "evals", "vs_base", "wrong", "broken", "fixed", "median_err", "max_err",
             "wall_min", "quiet"],
            formats={"evals": ",d", "median_err": ".4f", "max_err": ".4f", "wall_min": ".1f"},
        ),  # fmt: skip
        "c_stage": DT.markdown_table(
            rows_stage,
            list(rows_stage[0].keys()),
            formats={k: ",d" for k in rows_stage[0] if k != "arm"},
        ),
        "c_all": DT.markdown_table(
            rows_all, list(rows_all[0].keys()), formats={"six_evals": ",d", "e0_mean_evals": ",.0f"}
        ),
    }
    DT.write_tables(BASE / "step3" / "tables.md", tabs)
    print(tabs["c_all"])
    print("choice:", summ["_choice"])


if __name__ == "__main__":
    main()
