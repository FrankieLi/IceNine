"""AR3 step 8, part 2, remedy E1: refine in order of cost; stop at the first that passes a test.

  uv run python scripts/autoresearch/ar3/reseed_e1.py

Reads only the slim step files (``reseed/slim/steps_<sample>_<run>.json.gz``; for ``big`` run S5 the
rows of the logged run, ``S5L``). Nothing is run again: a seed step costs the strict cost on the
costed candidates and then one refinement per refined candidate, in the order of the strict cost.
The stored step runs all of them and keeps the one with the lowest final cost. E1 stops after the
first refinement whose end has hit ratio >= tau (the stop test) and keeps it; if none passes it
is the stored step. So E1's evaluations at a step are the costing plus the refinements up to and
including that one, and its result is that refinement's end. The pieces are those logged by the
replay (identical to the stored run's counts, see ``check_*.json``).

The plan's E1 stops at the seed test (tau 0.80). It failed the plan's rule (it broke more cases
than it fixed), so a tuned stop test was added after those results had been seen: the grid
``TAUS`` (0.80 to 0.98) is tried on the pick set of run S5, pooled over both samples. The rule:
the tau with the largest cut among those with broken <= fixed; if none, 0.80. The pick table and
the chosen tau are written to ``remedies_pick.json`` before the test rows are summed; then
``remedies_e1.json`` holds every set at tau 0.80 and at the chosen tau.

The pick set is the seed steps whose voxel lies in a true grain with an even index, the test set
those with an odd index.

Per set and per sample and run: steps; evaluations of the stored step and of E1; the cut; steps
where the first passing refinement is not the stored one's choice (``switched``); paired with the
stored step by whether the kept end is right (error at most 1 degree): fixed (stored wrong, E1
right), broken (stored right, E1 wrong); steps that pass the seed test at the end (the test the
caller applies); and the median error change.
"""

import json
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE.parents[2] / "benchmarks" / "autoresearch" / "ar3" / "reseed"
SLIM = OUT / "slim"
TEST = 0.80
RIGHT = 1.0
TAUS = (0.80, 0.85, 0.90, 0.95, 0.98)  # stop tests tried on the pick set (0.80 is the plan's)


def load(sample: str, run: str) -> List[Dict[str, Any]]:
    import gzip

    with gzip.open(SLIM / f"steps_{sample}_{run}.json.gz", "rt") as f:
        return json.load(f)


def e1_step(r: Dict[str, Any], tau: float = TEST) -> Dict[str, Any]:
    """Stored and E1 outcome of one seed step."""
    c = int(np.argmin(r["fin_costs"]))  # stored: lowest final cost (first of equal)
    s_hit, s_err = r["fin_hit"][c], r["fin_err"][c]
    stop = next((i for i, h in enumerate(r["fin_hit"]) if h >= tau), None)
    k = c if stop is None else stop
    n_ref = len(r["fin_evals"])
    ev_e1 = r["costed"] + sum(r["fin_evals"][: (n_ref if stop is None else stop + 1)])
    return dict(
        grain=r["grain"],
        level=r["level"],
        evals_stored=r["evals"],
        evals_e1=int(ev_e1),
        refine_stored=int(sum(r["fin_evals"])),
        refine_e1=int(sum(r["fin_evals"][: (n_ref if stop is None else stop + 1)])),
        costing=r["costed"],
        stopped=stop is not None,
        stop_at=-1 if stop is None else int(stop),
        switched=k != c,
        right_stored=bool(s_err <= RIGHT),
        right_e1=bool(r["fin_err"][k] <= RIGHT),
        pass_stored=bool(s_hit >= TEST),
        pass_e1=bool(r["fin_hit"][k] >= TEST),
        err_stored=s_err,
        err_e1=r["fin_err"][k],
    )


def summarize(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    if not rows:
        return dict(n=0)
    es = np.array([x["evals_stored"] for x in rows], float)
    e1 = np.array([x["evals_e1"] for x in rows], float)
    rs = np.array([x["right_stored"] for x in rows])
    r1 = np.array([x["right_e1"] for x in rows])
    return dict(
        n=len(rows),
        evals_stored=int(es.sum()),
        evals_e1=int(e1.sum()),
        cut=float(1 - e1.sum() / es.sum()),
        stopped=int(sum(x["stopped"] for x in rows)),
        switched=int(sum(x["switched"] for x in rows)),
        right_stored=int(rs.sum()),
        right_e1=int(r1.sum()),
        fixed=int((~rs & r1).sum()),
        broken=int((rs & ~r1).sum()),
        both_wrong=int((~rs & ~r1).sum()),
        pass_stored=int(sum(x["pass_stored"] for x in rows)),
        pass_e1=int(sum(x["pass_e1"] for x in rows)),
        median_err_change_deg=float(np.median([x["err_e1"] - x["err_stored"] for x in rows])),
        refine_share_stored=float(sum(x["refine_stored"] for x in rows) / es.sum()),
    )


def main() -> None:
    cases = [("big", "S5L"), ("big", "D4"), ("grad", "S5"), ("grad", "D4")]
    data = {c: load(*c) for c in cases}
    # ---- the pick: the stop test hit >= tau, tau chosen on the pick sets (even grains) of run S5
    # on both samples pooled, before any test row is summed. Rule: the tau with the largest cut of
    # evaluations among those with broken <= fixed on the pooled pick set; none: E1 is reported
    # with the plain seed test (0.80) and fails the part 2 rule.
    pick_rows = {
        t: [
            e1_step(r, t)
            for c in (("big", "S5L"), ("grad", "S5"))
            for r in data[c]
            if r["grain"] % 2 == 0
        ]
        for t in TAUS
    }
    pick_tab = {str(t): summarize(rows) for t, rows in pick_rows.items()}
    ok = [t for t in TAUS if pick_tab[str(t)]["broken"] <= pick_tab[str(t)]["fixed"]]
    chosen = max(ok, key=lambda t: pick_tab[str(t)]["cut"]) if ok else TEST
    (OUT / "remedies_pick.json").write_text(
        json.dumps(
            dict(
                E1="stop at the first refinement (strict-cost order) with hit ratio >= tau",
                sets="pick = even true-grain index, test = odd",
                taus=list(TAUS),
                pick_pooled_S5_both_samples=pick_tab,
                chosen_tau=chosen,
                chosen_by="largest cut among tau with broken <= fixed on the pooled pick set; "
                "none: tau 0.80 (the seed test)",
            ),
            indent=1,
        )
        + "\n"
    )
    out: Dict[str, Any] = {}
    for sample, run in cases:
        for tag, tau in (("tau0.80", TEST), ("tau_chosen", chosen)):
            rows = [e1_step(r, tau) for r in data[(sample, run)]]
            for name, sel in (
                ("pick", [x for x in rows if x["grain"] % 2 == 0]),
                ("test", [x for x in rows if x["grain"] % 2 == 1]),
                ("all", rows),
            ):
                out[f"{sample}/{run}/{tag}/{name}"] = summarize(sel)
            for L in sorted({x["level"] for x in rows}):
                out[f"{sample}/{run}/{tag}/level{L}"] = summarize(
                    [x for x in rows if x["level"] == L]
                )
    (OUT / "remedies_e1.json").write_text(
        json.dumps(dict(chosen_tau=chosen, taus=list(TAUS), results=out), indent=1) + "\n"
    )
    print("pick pooled", json.dumps(pick_tab, indent=0)[:1500])
    print("chosen tau", chosen)
    for k in out:
        if k.endswith(("test", "all")):
            print(k, out[k])


if __name__ == "__main__":
    main()
