"""AR5 step 23, phase A step 1: do step 22's stored records reproduce the two tables of the plan?

  uv run python scripts/autoresearch/ar5/ball_stored.py [--out ar5/ball/stored_check.json]

Reads (read only, through the links of ``ball_link_cache.sh``) step 22's ``J5D_<s>`` runs
(``run.json``, ``polish.json``, ``e1_log``, ``mc_traces``, ``defer_steps``), the J5SD arm's
``repair.json`` and writes ``stored_check.json``. Arithmetic only; no run, no image, no truth.
Nothing of ``held2``, ``held7`` or ``fine7`` (``ball_guard.guarded``).

Definitions, as in the plan: the seed step's share is ``index + rescue + polish + grid`` of the
J5SD total (J5D search + screen pass + wide repair); MC refinements: the evaluations of the rows of
``mc_traces`` (the seed step's MC refinements; a refinement DEF resumed is one row); BFS fits: the
neighbour and revisit evaluations; tail: the evaluations of a row after its last gain (its last
trace entry); *in resumed*: the part of the tail in rows flagged ``resumed``; stops at 0.95: the
rows that ``PH(500)`` stopped at 0.95 or more (``stop_kind == "high"``), and in brackets those
whose last gain is after evaluation 208 (after the first stage); ceiling at cap c: the evaluations
saved if every row ended at its own count or at ``208 + c + 1``, whichever is smaller.
"""

import gzip
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
sys.path.insert(0, str(HERE))
import ball_guard as BG  # noqa: E402

CACHE = BG.base_cache()
SAMPLES = ("held5", "held4", "held3", "held", "mid", "grad", "big", "held6")
FIRST_STAGE = 208

# the plan's tables (first: shares in %, then the second table)
T1 = {
    "held5": (20_076_940, 69.3, 34.9, 22.8, 24.0, 14, 7.9),
    "held4": (20_685_850, 67.2, 43.4, 24.0, 13.3, 8, 8.8),
    "held3": (18_363_331, 66.7, 39.3, 24.1, 16.9, 9, 9.2),
    "held": (17_114_589, 62.7, 43.4, 27.4, 8.1, 4, 9.9),
    "mid": (7_383_368, 45.4, 27.7, 43.0, 9.0, 2, 11.5),
    "grad": (6_134_231, 32.3, 13.3, 53.1, 11.0, 2, 14.6),
    "big": (2_452_596, 10.7, 6.7, 81.2, 0.0, 0, 8.1),
    "held6": (22_690_096, 70.0, 36.1, 22.0, 24.4, 16, 8.0),
}
# refinements, mean evals, tail %, in resumed %, stops at 0.95, (late last gain), c120, c300
T2 = {
    "held5": (5_926, 1_183, 22.5, 12.8, 835, 450, 25.2, 19.9),
    "held4": (6_695, 1_341, 26.8, 16.2, 1_027, 568, 32.8, 27.0),
    "held3": (5_797, 1_245, 25.3, 15.6, 774, 391, 28.9, 23.3),
    "held": (5_919, 1_256, 27.0, 15.4, 944, 506, 32.1, 25.9),
    "mid": (1_724, 1_184, 18.2, 9.9, 361, 165, 20.0, 15.8),
    "grad": (1_026, 793, 10.0, 4.5, 127, 29, 7.8, 4.8),
    "big": (172, 953, 4.5, 1.0, 44, 15, 4.4, 3.1),
    "held6": (6_382, 1_283, 23.5, 14.2, 950, 479, 26.8, 21.8),
}


def gz(p: Path) -> Any:
    with gzip.open(BG.guarded(p), "rt") as f:
        return json.load(f)


def js(p: Path) -> Any:
    return json.loads(BG.guarded(p).read_text())


def one(s: str) -> Dict[str, Any]:
    run = CACHE / f"J5D_{s}"
    rep = js(CACHE / f"stored_J5SD_{s}" / "repair.json")
    st = js(run / "run.json")["stats"]
    pol = js(run / "polish.json")
    search = int(
        st["n_evals_seed"] + st["n_evals_reseed"] + st["n_evals_neighbor"] + st["n_evals_revisit"]
    )
    total = search + int(rep["pass_evals"]) + int(rep["n_evals"])
    idx, resc, pols, grid = (
        int(st["index_evals"]),
        int(pol["rescue_evals"]),
        int(pol["polish_evals"]),
        int(pol["grid_evals"]),
    )
    bfs = int(st["n_evals_neighbor"]) + int(st["n_evals_revisit"])
    rows: List[Dict[str, Any]] = gz(run / "mc_traces.json.gz")
    ev = [int(r["evals"]) for r in rows]
    tail = [
        int(r["evals"]) - int(r["trace"][-1][0]) if r["trace"] else int(r["evals"]) for r in rows
    ]
    last = [int(r["trace"][-1][0]) if r["trace"] else 0 for r in rows]
    resumed = [bool(r.get("resumed")) for r in rows]
    high = [r for r in rows if r.get("stop_kind") == "high"]  # PH(500) stops at 0.95 or more
    out: Dict[str, Any] = dict(
        sample=s,
        total=total,
        seed_step_pct=100 * (idx + resc + pols + grid) / total,
        mc_pct=100 * sum(ev) / total,
        bfs_pct=100 * bfs / total,
        grid_pct=100 * grid / total,
        grid_count=int(pol["grid_fallbacks_left"]),
        screen_repair_pct=100 * (int(rep["pass_evals"]) + int(rep["n_evals"])) / total,
        refinements=len(rows),
        mean_evals=sum(ev) / max(1, len(ev)),
        tail_pct=100 * sum(tail) / total,
        resumed_pct=100 * sum(t for t, r in zip(tail, resumed) if r) / total,
        stops_095=len(high),
        stops_095_after_stage1=sum(1 for r in high if r["trace"][-1][0] > FIRST_STAGE),
    )
    for c in (120, 300):
        cap = FIRST_STAGE + c + 1
        out[f"ceiling_{c}_pct"] = 100 * sum(max(0, e - min(e, cap)) for e in ev) / total
    return out


def main() -> int:
    out = {s: one(s) for s in SAMPLES}
    res: Dict[str, Any] = dict(rows=out, table1={}, table2={})
    for s in SAMPLES:
        o = out[s]
        t = T1[s]
        res["table1"][s] = dict(
            total=(o["total"], t[0]),
            seed=(round(o["seed_step_pct"], 1), t[1]),
            mc=(round(o["mc_pct"], 1), t[2]),
            bfs=(round(o["bfs_pct"], 1), t[3]),
            grid=(round(o["grid_pct"], 1), t[4]),
            grid_count=(o["grid_count"], t[5]),
            screen_repair=(round(o["screen_repair_pct"], 1), t[6]),
        )
        u = T2[s]
        res["table2"][s] = dict(
            refinements=(o["refinements"], u[0]),
            mean_evals=(round(o["mean_evals"]), u[1]),
            tail=(round(o["tail_pct"], 1), u[2]),
            resumed=(round(o["resumed_pct"], 1), u[3]),
            stops=(o["stops_095"], u[4]),
            stops_late=(o["stops_095_after_stage1"], u[5]),
            ceil120=(round(o["ceiling_120_pct"], 1), u[6]),
            ceil300=(round(o["ceiling_300_pct"], 1), u[7]),
        )
    bad: List[str] = []
    for tab in ("table1", "table2"):
        for s, d in res[tab].items():
            for k, (mine, plan) in d.items():
                if mine != plan:
                    bad.append(f"{tab} {s} {k}: stored {mine} plan {plan}")
    res["differences"] = bad
    res["reproduces"] = not bad
    outp = Path(sys.argv[sys.argv.index("--out") + 1]) if "--out" in sys.argv else None
    if outp:
        outp.write_text(json.dumps(res, indent=1) + "\n")
    print(
        json.dumps({k: v for k, v in res.items() if k in ("differences", "reproduces")}, indent=1)
    )
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
