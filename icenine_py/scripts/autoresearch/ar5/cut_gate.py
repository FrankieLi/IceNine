"""AR5 step 10, phase A steps 1 and 4: the reproduction check and the gate (disclosed cases only).

  uv run python scripts/autoresearch/ar5/cut_gate.py

Reads ``cut/replay_<sample>_n0_repro.json`` (the 2 full searches rerun), ``cut/replay_<sample>_n1000
.json`` (the cut search on the 23 disclosed cases) and step 9's stored record (``cache/levels_
<sample>.json``). Writes ``cut/gate.json`` and prints the verdict. Stop rule: a stored total that
does not reproduce (335,077 for mid 7167; 341,304 for grad 19058; the record's totals 4,007,241
and 3,726,723). Gate (plan "Decision rules"): the cut search is right (error at most 1 degree) in
at least 22 of 23, reaches hit ratio 0.95 in at least 20, and the 23 need at least 50% fewer
evaluations than the full searches (7,733,964). Also prints prediction 1's items. No ``held`` file.
No tests of significance; counts, Wilson intervals (``scripts/common/stats.wilson``) and the paired
split against the full search, by voxel (``stats.paired_discordant``). Also writes
``cut/tables.md`` (``doc_tables.write_tables``) for the MIGRATION_HISTORY section.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Sequence

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent / "common"))
CUT = HERE.parents[2] / "benchmarks" / "autoresearch" / "ar5" / "cut"
STORED_TOTAL = {"mid": 4_007_241, "grad": 3_726_723}
STORED_ONE = {("mid", 7167): 335_077, ("grad", 19058): 341_304}
FULL_ALL = 7_733_964

import doc_tables as DT  # noqa: E402
import stats as ST  # noqa: E402


def wilson(k: int, n: int) -> List[float]:
    """``stats.wilson`` (95%) as a list, the form stored in ``gate.json``."""
    return list(ST.wilson(k, n))


def paired(ok_cut: Sequence[bool], ok_full: Sequence[bool]) -> Dict[str, int]:
    """The paired split by voxel: both, fixed (cut only), broken (full search only), neither."""
    broken, fixed = ST.paired_discordant([not x for x in ok_cut], [not x for x in ok_full])
    both = sum(1 for a, b in zip(ok_cut, ok_full) if a and b)
    return dict(both=both, fixed=fixed, broken=broken, neither=len(ok_cut) - both - fixed - broken)


def yn(ok: bool) -> str:
    return "yes" if ok else "no"


def fmt_wilson(w: Sequence[float]) -> str:
    return f"{w[0]:.2f} to {w[1]:.2f}"


def tables(rec: Dict[str, Any], rows: List[Dict[str, Any]]) -> Dict[str, str]:
    """The doc tables of the section (generated; copy with ``sync_doc_tables.py``)."""
    rp = rec["replay"]
    n = rp["cases"]
    gate_rows = [
        {
            "item": "right (at most 1 degree)",
            "cut search": f"{rp['right']} of {n}",
            "Wilson 95%": fmt_wilson(rp["right_wilson95"]),
            "gate": "at least 22",
            "full search": f"{rp['full_search_right']} of {n}",
            "met": yn(rp["right"] >= 22),
        },
        {
            "item": "hit ratio at least 0.95",
            "cut search": f"{rp['reach_095']} of {n}",
            "Wilson 95%": fmt_wilson(rp["reach_095_wilson95"]),
            "gate": "at least 20",
            "full search": f"{rp['full_search_reach_095']} of {n} "
            f"({fmt_wilson(rp['full_search_reach_095_wilson95'])})",
            "met": yn(rp["reach_095"] >= 20),
        },
        {
            "item": "evaluations, 23 searches",
            "cut search": f"{rp['evals_total']:,}",
            "Wilson 95%": "n/a",
            "gate": "at least 50% fewer",
            "full search": f"{rp['evals_full_total']:,}",
            "met": f"{yn(rp['fewer_share'] >= 0.5)} ({100 * rp['fewer_share']:.1f}% fewer)",
        },
    ]
    cols = ["item", "cut search", "Wilson 95%", "gate", "full search", "met"]
    pr = rp["paired_vs_full"]
    pair_rows = [
        dict(count=k, both=v["both"], fixed=v["fixed"], broken=v["broken"], neither=v["neither"])
        for k, v in (("right", pr["right"]), ("hit ratio at least 0.95", pr["reach_095"]))
    ]
    case_rows = [
        {
            "case": f"{r['sample']} {int(r['voxel'])}",
            "cut end, degrees": f"{r['err']:.2f}",
            "cut hit": f"{r['hit']:.2f}",
            "full end, degrees": f"{r['full_end_err']:.3f}",
            "full hit": f"{r['full_end_hit']:.2f}",
            "cut end to full end, degrees": f"{r['mis_to_full_end']:.2f}",
        }
        for r in rows
        if r["err"] > 1.0 or r["hit"] < 0.95
    ]
    case_cols = list(case_rows[0]) if case_rows else ["case"]
    p1 = rec["prediction_1"]
    pred_rows = [
        dict(prediction="right in 23", result=f"{rp['right']}", met=yn(p1["right_23"])),
        dict(
            prediction="hit ratio at least 0.95 in 20 or 21",
            result=f"{rp['reach_095']}",
            met=yn(p1["reach_20_or_21"]),
        ),
        dict(
            prediction="within 0.1 degree of the full end in at least 18",
            result=f"{rp['within_0p1_of_full_end']}",
            met=yn(p1["within_0p1_at_least_18"]),
        ),
        dict(
            prediction="85,000 to 135,000 evaluations each",
            result=f"{rp['evals_per_search_min']:,} to {rp['evals_per_search_max']:,}",
            met=yn(p1["evals_85k_135k_each"]),
        ),
    ]
    return {
        "cut_gate": DT.markdown_table(gate_rows, cols),
        "cut_paired": DT.markdown_table(pair_rows, ["count", "both", "fixed", "broken", "neither"]),
        "cut_cases": DT.markdown_table(case_rows, case_cols),
        "cut_pred1": DT.markdown_table(pred_rows, ["prediction", "result", "met"]),
    }


def main() -> None:
    rec: Dict[str, Any] = {}
    # reproduction
    for s in ("mid", "grad"):
        stored = json.loads((CUT / "cache" / f"levels_{s}.json").read_text())
        rep = json.loads((CUT / f"replay_{s}_n0_repro.json").read_text())["rows"]
        assert len(rep) == 1
        v = int(rep[0]["voxel"])
        srow = [r for r in stored["rows"] if int(r["voxel"]) == v][0]
        rec[f"repro_{s}"] = dict(
            voxel=v,
            rerun_evals=int(rep[0]["evals"]),
            expected=STORED_ONE[(s, v)],
            stored_row_evals=int(srow["grid_evals"]),
            reproduces=bool(int(rep[0]["evals"]) == STORED_ONE[(s, v)] == int(srow["grid_evals"])),
            record_total=int(stored["grid_evals_total"]),
            record_total_expected=STORED_TOTAL[s],
            record_total_ok=bool(int(stored["grid_evals_total"]) == STORED_TOTAL[s]),
            record_cases=int(stored["cases"]),
        )
    # replay
    rows: List[Dict[str, Any]] = []
    for s in ("mid", "grad"):
        rows += json.loads((CUT / f"replay_{s}_n1000.json").read_text())["rows"]
    n = len(rows)
    right = sum(1 for r in rows if r["err"] <= 1.0)
    reach = sum(1 for r in rows if r["hit"] >= 0.95)
    full_reach = sum(1 for r in rows if r["full_end_hit"] >= 0.95)
    full_right = sum(1 for r in rows if r["full_end_err"] <= 1.0)
    near = sum(1 for r in rows if r["mis_to_full_end"] is not None and r["mis_to_full_end"] <= 0.1)
    n_mis = sum(1 for r in rows if r["mis_to_full_end"] is not None)
    ev = [int(r["evals"]) for r in rows]
    total = int(sum(ev))
    saved = 1.0 - total / FULL_ALL
    rec["replay"] = dict(
        cases=n,
        right=right,
        right_wilson95=wilson(right, n),
        reach_095=reach,
        full_search_reach_095=full_reach,
        within_0p1_of_full_end=near,
        n_with_full_end=n_mis,
        evals_total=total,
        evals_full_total=FULL_ALL,
        fewer_share=saved,
        evals_per_search_min=int(min(ev)),
        evals_per_search_max=int(max(ev)),
        evals_per_search_median=float(np.median(ev)),
        wrong_cases=[
            dict(
                sample=r["sample"], voxel=int(r["voxel"]), err=float(r["err"]), hit=float(r["hit"])
            )
            for r in rows
            if r["err"] > 1.0
        ],
        not_reaching_095=[
            dict(
                sample=r["sample"],
                voxel=int(r["voxel"]),
                hit=float(r["hit"]),
                full_end_hit=float(r["full_end_hit"]),
            )
            for r in rows
            if r["hit"] < 0.95
        ],
        level0_run_values=sorted({int(r["level0_run"]) for r in rows}),
        reach_095_wilson95=wilson(reach, n),
        full_search_right=full_right,
        full_search_reach_095_wilson95=wilson(full_reach, n),
        wrong_cases_full_end=[
            dict(
                sample=r["sample"],
                voxel=int(r["voxel"]),
                full_end_err=float(r["full_end_err"]),
                full_end_hit=float(r["full_end_hit"]),
            )
            for r in rows
            if r["err"] > 1.0
        ],
        paired_vs_full=dict(
            right=paired([r["err"] <= 1.0 for r in rows], [r["full_end_err"] <= 1.0 for r in rows]),
            reach_095=paired(
                [r["hit"] >= 0.95 for r in rows], [r["full_end_hit"] >= 0.95 for r in rows]
            ),
        ),
    )
    repro_ok = all(
        rec[f"repro_{s}"]["reproduces"] and rec[f"repro_{s}"]["record_total_ok"]
        for s in ("mid", "grad")
    )
    gate_ok = bool(right >= 22 and reach >= 20 and saved >= 0.5)
    pred = dict(
        right_23=bool(right == 23),
        reach_20_or_21=bool(reach in (20, 21)),
        within_0p1_at_least_18=bool(near >= 18),
        evals_85k_135k_each=bool(min(ev) >= 85_000 and max(ev) <= 135_000),
    )
    rec["reproduction_ok"] = bool(repro_ok)
    rec["gate_pass"] = gate_ok
    rec["prediction_1"] = pred
    rec["verdict"] = "READY TO RUN" if (repro_ok and gate_ok) else "STOP"
    (CUT / "gate.json").write_text(json.dumps(rec, indent=1) + "\n")
    DT.write_tables(CUT / "tables.md", tables(rec, rows))
    print(json.dumps({k: rec[k] for k in ("reproduction_ok", "gate_pass", "verdict")}))
    print(json.dumps(rec["replay"], indent=1))
    print(json.dumps(pred))
    for s in ("mid", "grad"):
        print(json.dumps(rec[f"repro_{s}"]))


if __name__ == "__main__":
    main()
