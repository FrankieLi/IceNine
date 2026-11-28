"""AR5 step 24, phase A step 1 and gate g1: the sizing of rule ``BF`` from stored rows (arithmetic
only; no run, no image, no truth).

  uv run python scripts/autoresearch/ar5/bf_size.py [--write]

Reads (read only, through the links of ``bf_link_cache.sh``) step 23's replay rows
(``replay_on_<s>/rows.jsonl``: the ``CB`` step for every stored seed step of step 22's live ``J5D``
run), step 22's ``J5D_<s>`` records (``e1_log``, ``mc_traces``, ``run.json``) and the ``J5SD``
arm's ``repair.json``. A seed step under ``BF`` either ends in the ball pass (step 23's row) or
falls back to ``J5D``'s own stage 2 (step 22's trace), so its cost is exact: the ball pass's
evaluations plus, for every refinement ``J5D`` ran, its evaluations beyond the first stage.

Checks (g1, frozen in ``bf_rules.GATE``): the three tables of the plan are reproduced exactly (the
plan's strings are in ``PLAN_T1``, ``PLAN_T2``, ``PLAN_T3`` below); the accounting identity
(costing + CMA fits + MC refinements = ``J5D``'s step) holds in every step; lost 0.95 ends (J5D has
one, BF has none) are 0 (at most 3% of the steps; a Wilson upper bound is written); the index stage
is at least 10% fewer than ``J5D``'s on the held types and ``mid``. With ``--write`` it writes
``ar5/bf/size.json``, ``size_tables.md`` (tables through ``doc_tables``) and exits 1 on any
difference. Nothing of ``held2``, ``held7`` or ``fine7`` is opened (``bf_guard.guarded``).
"""

import gzip
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICE / "scripts" / "common"))
import bf_guard as BG  # noqa: E402
import bf_rules as BRU  # noqa: E402
import doc_tables as DT  # noqa: E402
import stats as ST  # noqa: E402

CACHE = BG.base_cache()
RD = ICE / "benchmarks" / "autoresearch" / "ar5" / "bf"
SAMPLES = ("held5", "held4", "held3", "held", "mid", "grad", "big", "held6")

# the plan's first table: J5D index, CB index, fall back, BF index, BF fewer %, of J5SD %, lost,
# seed test J5D / CB / BF
PLAN_T1 = {
    "held5": (
        "8,690,640",
        "3,398,238",
        "773",
        "7,376,267",
        "15.1%",
        "6.55%",
        "0",
        "104 / 540 / 96",
    ),
    "held4": (
        "10,882,299",
        "3,884,631",
        "899",
        "9,076,725",
        "16.6%",
        "8.73%",
        "0",
        "70 / 558 / 68",
    ),
    "held3": ("8,840,378", "3,287,015", "813", "7,856,813", "11.1%", "5.36%", "0", "90 / 528 / 85"),
    "held": ("9,129,537", "3,433,723", "777", "7,666,627", "16.0%", "8.55%", "0", "66 / 480 / 64"),
    "mid": ("2,644,672", "1,099,560", "173", "2,099,764", "20.6%", "7.38%", "0", "12 / 67 / 12"),
    "grad": ("1,255,626", "690,503", "97", "1,023,858", "18.5%", "3.78%", "0", "11 / 35 / 10"),
    "big": ("261,653", "143,064", "6", "192,919", "26.3%", "2.80%", "0", "0 / 1 / 0"),
    "held6": (
        "10,054,329",
        "3,684,745",
        "804",
        "8,368,681",
        "16.8%",
        "7.43%",
        "0",
        "63 / 504 / 62",
    ),
}
# second table: steps, first seeds, first-seed saving % of J5SD, first lost, gained, moved over 0.05
PLAN_T2 = {
    "held5": ("3,556", "282", "1.13%", "7", "83", "90"),
    "held4": ("4,109", "338", "1.37%", "19", "107", "124"),
    "held3": ("3,391", "293", "1.38%", "10", "72", "79"),
    "held": ("3,685", "297", "1.52%", "11", "80", "92"),
    "mid": ("1,626", "162", "0.74%", "2", "27", "51"),
    "grad": ("1,281", "348", "0.52%", "0", "20", "19"),
    "big": ("330", "39", "0.07%", "0", "0", "5"),
    "held6": ("4,109", "339", "1.19%", "19", "113", "105"),
}
# third table: lost, never, stock stage 2 in never %, never steps %, ball hit lost / never,
# J5D never, best cand.
PLAN_ROUNDING = {
    "held6 t3": "the plan prints 16.4%; the stored count is 3,709,750 of 22,690,096 = 16.3496%, "
    "which rounds to 16.3%; a display slip of the plan, the count reproduces"
}
PLAN_T3 = {
    "held5": ("309", "464", "15.7%", "20.3%", "0.761 / 0.736", "0.886", "92"),
    "held4": ("406", "493", "19.1%", "23.8%", "0.771 / 0.759", "0.896", "94"),
    "held3": ("284", "529", "20.6%", "26.3%", "0.765 / 0.744", "0.893", "70"),
    "held": ("352", "425", "18.8%", "23.7%", "0.762 / 0.760", "0.899", "112"),
    "mid": ("70", "103", "11.1%", "13.9%", "0.847 / 0.817", "0.888", "18"),
    "grad": ("19", "78", "5.0%", "7.5%", "0.829 / 0.858", "0.886", "15"),
    "big": ("2", "4", "1.9%", "2.3%", "0.849 / 0.894", "0.945", "1"),
    "held6": ("330", "474", "16.4%", "20.5%", "0.765 / 0.759", "0.901", "98"),
}
T1_COLS = (
    "J5D index",
    "CB index",
    "fall back",
    "BF index",
    "BF fewer",
    "of J5SD",
    "lost",
    "seed test J5D/CB/BF",
)
T1_PLAN_COLS = T1_COLS[:-1] + ("seed test J5D/CB/BF",)


def gz(p: Path) -> Any:
    with gzip.open(BG.guarded(p), "rt") as f:
        return json.load(f)


def js(p: Path) -> Any:
    return json.loads(BG.guarded(p).read_text())


def rows_of(s: str) -> List[Dict[str, Any]]:
    p = BG.guarded(CACHE / f"replay_on_{s}" / "rows.jsonl")
    return [json.loads(x) for x in p.read_text().splitlines() if x]


def j5sd_total(s: str) -> int:
    st = js(CACHE / f"J5D_{s}" / "run.json")["stats"]
    rep = js(CACHE / f"stored_J5SD_{s}" / "repair.json")
    search = int(
        st["n_evals_seed"] + st["n_evals_reseed"] + st["n_evals_neighbor"] + st["n_evals_revisit"]
    )
    return search + int(rep["pass_evals"]) + int(rep["n_evals"])


def med(a: List[float]) -> float:
    return float(np.median(np.asarray(a, float))) if len(a) else float("nan")


def one(s: str) -> Dict[str, Any]:
    rows = rows_of(s)
    e1 = gz(CACHE / f"J5D_{s}" / "e1_log.json.gz")
    total = j5sd_total(s)
    j5d = sum(r["j5d_evals"] for r in rows)
    cb = sum(r["evals"] for r in rows)
    id_bad = 0
    for r in rows:
        e = e1[r["index"]]
        c_ev = sum(f[5] for f in e["fits"] if f[2] == "C")
        m_ev = sum(f[5] for f in e["fits"] if f[2] == "M")
        if e["costed"] + c_ev + m_ev != r["j5d_evals"]:
            id_bad += 1
    bf = 0
    fb_steps = 0
    first = dict(j5d=0, cb=0, n=0, lost=0)
    lost = gained = moved_found = 0
    st2_never = never_step_j5d = 0
    n_lost = n_never = best_cand = 0
    hit_lost: List[float] = []
    hit_never: List[float] = []
    j5d_never: List[float] = []
    for r in rows:
        if not r["reseed"]:
            first["j5d"] += r["j5d_evals"]
            first["cb"] += r["evals"]
            first["n"] += 1
            first["lost"] += int(r["lost_end"])
        gained += int(r["gained_end"])
        # BF's end: the ball pass's when it has one, else J5D's own step (the fallback)
        bf_stop = r["stop"] if r["stop"] >= 0 else r["j5d_stop"]
        lost += int(BRU.lost_end(r["j5d_stop"], bf_stop))
        if r["stop"] >= 0:
            bf += r["evals"]
            moved_found += int(r["taken_moved_deg"] > 0.05)
            continue
        fb_steps += 1
        extra = 0
        st2 = 0
        for rf in r.get("refinements", []):
            if rf["j5d_evals"] is not None:
                extra += rf["j5d_evals"] - rf["stage1"]
                st2 += rf["j5d_evals"] - rf["stage1"]
        bf += r["evals"] + extra
        if r["lost_end"]:
            n_lost += 1
            hit_lost.append(r["taken"][2])
            if r["j5d_stop"] == r["taken"][0]:
                best_cand += 1
        else:
            n_never += 1
            hit_never.append(r["taken"][2])
            j5d_never.append(r["j5d_taken"][2])
            st2_never += st2
            never_step_j5d += r["j5d_evals"]
    seed_j5d = sum(r["j5d_taken"][2] < 0.80 for r in rows)
    seed_cb = sum(r["taken"][2] < 0.80 for r in rows)
    seed_bf = sum(r["j5d_taken"][2] < 0.80 for r in rows if r["stop"] < 0)
    lo, hi = ST.wilson(lost, len(rows))
    return dict(
        sample=s,
        steps=len(rows),
        id_bad=id_bad,
        j5sd=total,
        j5d=j5d,
        cb=cb,
        bf=bf,
        fb_steps=fb_steps,
        fewer=100 * (1 - bf / j5d),
        of_j5sd=100 * (j5d - bf) / total,
        lost=lost,
        lost_wilson_pct=[100 * lo, 100 * hi],
        seed_j5d=int(seed_j5d),
        seed_cb=int(seed_cb),
        seed_bf=int(seed_bf),
        first_seeds=first["n"],
        first_saving=100 * (first["j5d"] - first["cb"]) / total,
        first_lost=first["lost"],
        gained=gained,
        moved_found=moved_found,
        n_lost_cb=n_lost,
        n_never=n_never,
        never_st2=100 * st2_never / total,
        never_steps=100 * never_step_j5d / total,
        hit_lost=med(hit_lost),
        hit_never=med(hit_never),
        j5d_never=med(j5d_never),
        best_cand=best_cand,
    )


def f1(x: float) -> str:
    return f"{x:.1f}"


def strings(o: Dict[str, Any]) -> Dict[str, Any]:
    t1 = (
        f"{o['j5d']:,}",
        f"{o['cb']:,}",
        f"{o['fb_steps']:,}",
        f"{o['bf']:,}",
        f"{o['fewer']:.1f}%",
        f"{o['of_j5sd']:.2f}%",
        f"{o['lost']:,}",
        f"{o['seed_j5d']} / {o['seed_cb']} / {o['seed_bf']}",
    )
    t2 = (
        f"{o['steps']:,}",
        f"{o['first_seeds']:,}",
        f"{o['first_saving']:.2f}%",
        f"{o['first_lost']}",
        f"{o['gained']}",
        f"{o['moved_found']}",
    )
    t3 = (
        f"{o['n_lost_cb']}",
        f"{o['n_never']}",
        f"{o['never_st2']:.1f}%",
        f"{o['never_steps']:.1f}%",
        f"{o['hit_lost']:.3f} / {o['hit_never']:.3f}",
        f"{o['j5d_never']:.3f}",
        f"{o['best_cand']}",
    )
    return dict(t1=t1, t2=t2, t3=t3)


def main() -> int:
    write = "--write" in sys.argv
    res = {s: one(s) for s in SAMPLES}
    diffs: List[str] = []
    notes: List[str] = []
    for s, o in res.items():
        st = strings(o)
        for tab, plan in (("t1", PLAN_T1), ("t2", PLAN_T2), ("t3", PLAN_T3)):
            if st[tab] != plan[s]:
                msg = f"{s} {tab}: stored {st[tab]} plan {plan[s]}"
                if f"{s} {tab}" in PLAN_ROUNDING:  # a display slip of the plan, not a count
                    notes.append(msg + f" ({PLAN_ROUNDING[f'{s} {tab}']})")
                else:
                    diffs.append(msg)
        if o["id_bad"]:
            diffs.append(f"{s}: accounting identity fails in {o['id_bad']} steps")
        if o["lost"] > BRU.GATE["g1_lost_max"] * o["steps"]:
            diffs.append(f"{s}: lost ends {o['lost']} over the limit")
        if s in BRU.GATE["g1_fewer_samples"] and o["fewer"] < 100 * BRU.GATE["g1_fewer_min"]:
            diffs.append(f"{s}: index stage {o['fewer']:.1f}% fewer, under the gate")
    out = dict(
        rows=res,
        differences=diffs,
        plan_rounding_notes=notes,
        reproduces=not diffs,
        g1_passes=not diffs,
    )
    if write:
        RD.mkdir(parents=True, exist_ok=True)
        (RD / "size.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
        tabs = {}
        t1_rows = [dict(sample=s, **dict(zip(T1_COLS, strings(res[s])["t1"]))) for s in SAMPLES]
        tabs["ar5_bf_size_index"] = DT.markdown_table(t1_rows, ("sample",) + T1_COLS)
        c2 = (
            "steps",
            "first seeds",
            "first-seed saving",
            "first lost",
            "gained",
            "moved over 0.05",
        )
        t2_rows = [dict(sample=s, **dict(zip(c2, strings(res[s])["t2"]))) for s in SAMPLES]
        tabs["ar5_bf_size_first"] = DT.markdown_table(t2_rows, ("sample",) + c2)
        c3 = (
            "lost",
            "never",
            "stock stage 2 in never",
            "never steps",
            "ball hit lost / never",
            "J5D never",
            "best cand.",
        )
        t3_rows = [dict(sample=s, **dict(zip(c3, strings(res[s])["t3"]))) for s in SAMPLES]
        tabs["ar5_bf_size_fallback"] = DT.markdown_table(t3_rows, ("sample",) + c3)
        DT.write_tables(RD / "size_tables.md", tabs)
    print(
        json.dumps(
            {
                "reproduces": out["reproduces"],
                "differences": diffs,
                "lost": {s: res[s]["lost"] for s in SAMPLES},
                "lost_wilson_upper_pct": {
                    s: round(res[s]["lost_wilson_pct"][1], 3) for s in SAMPLES
                },
            },
            indent=1,
        )
    )
    return 0 if not diffs else 1


if __name__ == "__main__":
    sys.exit(main())
