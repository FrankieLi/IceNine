"""AR5 step 22 (task ar5v): rule ``DEF``, defer a stalled seed refinement while untried candidates
remain (pure: no image, no reconstructor, no truth; numpy not needed).

Arm ``J5D`` is arm ``J5`` (step 16: rule ``PH(500)`` in the seed step's MC refinement) with rule
``DEF`` added; ``J5SD`` is ``J5D`` with step 21's screen pass and wide repair (imported unchanged).

**Terms.** *Check*: ``J5``'s ``PH(500)`` check, one counted evaluation of the running best's OWN hit
ratio (never a ratio to the neighbours) after 500 evaluations with no gain. At 0.95 or more the
refinement stops, else it goes on. *Pause* / *resume*: keep a refinement's state and generator /
continue it exactly where it stopped.

**Rule ``DEF``.** At a check that goes on, if the seed step has an untried candidate, the
refinement is paused and the next candidate is tried. The step stops at the first end at 0.95 or
more, as now. If no later candidate reaches 0.95, the paused refinements are resumed in candidate
order under ``J5``'s rules until one reaches 0.95; else the lowest-cost end is taken, as now. The
last candidate is never paused. No new check, no threshold, no neighbour, no model.

This file holds the constants, the decision rules and predictions of the plan (frozen), the
walk of one stored seed step (``walk_step``: the replay's arithmetic, from ``J5``'s stored
``e1_log`` row and ``mc_traces``) and the seed step's state machine (``run_seed_step``: shared by
the live arm and the tests; the fits are callbacks).
"""

from typing import Any, Callable, Dict, List, Optional, Tuple

TAU = 0.95  # the seed step's stop value (unchanged)
K = 500  # the stall size of PH(500) (unchanged)
N_CAND = 5  # candidates refined in the seed step (unchanged)
ARM_D = "J5D"
ARM_SD = "J5SD"
NEW_SAMPLE = "held6"
GENERATOR_SEED = 6  # held6 = held with generator seed 6: a new sample, not an extra random seed
CMA_EVALS = 201  # one CMA fit of the seed step (C200: 200 plus the final evaluation)
DISCLOSED = ("held5", "held4", "held3", "held", "mid", "grad", "big")
LIVE_REST = ("held5", "held4", "held3", "grad", "big")  # the gate's mid and held come first

RULE = (
    "Rule DEF: in the seed step (interleaved order, 20 costed, the best 5 refined in turn: a CMA "
    "fit and, under 0.95, an MC refinement), at J5's PH(500) check that goes on (the running "
    "best's own hit ratio under 0.95 after 500 evaluations with no gain), if the step has an "
    "untried candidate the refinement is paused (its state and generator kept) and the next "
    "candidate is tried. The step stops at the first end at 0.95 or more, as now. If no later "
    "candidate reaches 0.95 the paused refinements are resumed in candidate order, exactly where "
    "they paused and under J5's rules, until one reaches 0.95; else the lowest-cost end is taken, "
    "as now. The last candidate is never paused. CMA fits, 5 candidates, 0.80, 0.95, PH(500), G1, "
    "G2b, the grid search, BFS, the screen pass and the wide repair are unchanged. Steps of type A "
    "(found at a later candidate than every paused refinement) and N (no end at 0.95) are exactly "
    "J5's; only C steps (the found refinement would itself be paused) change. Arm J5D = arm J5 "
    "with DEF; arm J5SD = J5D with step 21's screen pass and wide repair."
)

# ---- the plan's table (prediction 1): reproduced exactly by the tracked replay (g1)
# total = J5S total (J5 search + screen pass + wide repair); a_pct: exact saving of type A steps as
# a share of it (2 decimals); net bounds (1 decimal): C's extra cost from 201 evaluations a step up
# to (4 - k) x the mean cost of a candidate run up to its pause.
PLAN_TABLE: Dict[str, Dict[str, Any]] = {
    "held5": dict(
        total=22_706_889, a_steps=998, a_pct=10.54, c_steps=191, n_steps=546, lo=9.1, hi=11.7
    ),
    "held4": dict(
        total=23_612_288, a_steps=1_106, a_pct=13.14, c_steps=251, n_steps=593, lo=11.7, hi=15.1
    ),
    "held3": dict(
        total=20_240_655, a_steps=869, a_pct=9.75, c_steps=172, n_steps=599, lo=8.4, hi=11.0
    ),
    "held": dict(
        total=19_457_817, a_steps=1_010, a_pct=12.28, c_steps=238, n_steps=512, lo=10.2, hi=14.1
    ),
    "mid": dict(
        total=8_397_057, a_steps=358, a_pct=10.65, c_steps=79, n_steps=127, lo=8.8, hi=11.9
    ),
    "grad": dict(total=6_257_077, a_steps=269, a_pct=2.40, c_steps=10, n_steps=98, lo=2.0, hi=2.4),
    "big": dict(total=2_631_446, a_steps=62, a_pct=4.95, c_steps=9, n_steps=4, lo=4.5, hi=5.6),
}

PREDICTIONS: Dict[str, str] = {
    "1": "The tracked replay reproduces the plan's table exactly on the seven disclosed samples: A "
    "steps, exact savings, C steps and N steps.",
    "2": "Live J5SD against J5S, disclosed samples: savings inside the net bounds widened by 1.5 "
    "points; 30% to 70% of C steps end at another candidate; broken minus fixed -3 to +3; grid "
    "searches within 2 of J5's.",
    "3": "held6: 8% to 14% fewer evaluations than J5S; broken minus fixed against J5S -3 to +3; "
    "median error within 0.0005 degree; the house rule against S8 passes.",
}

GATE: Dict[str, Any] = dict(
    g0="with DEF off the scripts reproduce J5's seed steps on mid (ends, evaluations); a resumed "
    "refinement equals J5's unpaused one, evaluation by evaluation",
    g1="prediction 1, exactly",
    g2="the live J5D run on mid equals the replay in every seed step up to the first C step in "
    "run order",
    g3_samples=("mid", "held"),
    g3_fewer_min=0.07,  # J5SD against J5S, share of evaluations
    g3_broken_minus_fixed_max=3,
    g3_median_error_higher_max_deg=0.001,
    g3_grid_more_max=2,  # grid searches more than J5's
    fails="stop and report: no tuning, no threshold, no held6",
)

JUDGE: Dict[str, Any] = dict(
    r1_fewer_min=dict(
        held5=0.07, held4=0.07, held3=0.07, held=0.07, mid=0.07, grad=0.015, big=0.03
    ),
    r2_broken_minus_fixed_max=3,
    r3_median_error_higher_max_deg=0.001,
    r4_vs_s8=dict(fewer=0.20, broken_minus_fixed_max=5, median_error_higher_max_deg=0.005),
    passes="J5SD passes on a sample if (1) to (4) hold against J5S (and S8 for 4)",
    step_passes="J5SD passes on held6 and on the seven disclosed samples",
    fails_on_held6="no key value; study merge",
)

CAPS_MIN: Dict[str, int] = dict(S8=480, rest=300)


# ---- the walk of one stored seed step (the replay's arithmetic) ---------------------------------
def first_went_on(trace: Dict[str, Any]) -> Optional[List[Any]]:
    """The first check of a stored MC refinement that went on (hit ratio under 0.95), or None.
    Checks are ``[evaluation, hit ratio, stopped]``; every check before a first one that goes on
    would have stopped the refinement, so this is the check where DEF pauses it."""
    return next((c for c in trace["checks"] if c[1] < TAU), None)


def walk_step(
    entry: Dict[str, Any], traces: Dict[int, Dict[str, Any]], n_pos: int = N_CAND
) -> Dict[str, Any]:
    """DEF on one stored J5 seed step (``entry``: the ``e1_log`` row; ``traces``: its MC
    refinements' ``mc_traces`` rows by candidate position).

    Returns ``type`` (``"A"``, ``"N"``, ``"C"`` or ``"-"`` when no refinement is paused before the
    step's end), ``stop`` (J5's stop position, -1 none), ``paused`` (positions and the evaluations
    each would keep paying after the pause: ``tail``), ``saving`` (A: the exact saving; C: the
    earlier paused tails, which are saved whatever happens later), ``expected_evals`` (A and N:
    exactly the step's evaluations under DEF; None for C), and, for N steps, the cost of each
    candidate run up to its pause (``cand_cost``, the sizing of C's extra cost)."""
    stop = int(entry["stop"])
    paused: Dict[int, int] = {}
    for p, t in traces.items():
        c = first_went_on(t)
        if c is not None and p < n_pos - 1:
            paused[int(p)] = int(t["evals"]) - (int(c[0]) + 1)
    out: Dict[str, Any] = dict(
        voxel=int(entry["voxel"]),
        stop=stop,
        paused={int(p): int(v) for p, v in paused.items()},
        j5_evals=int(entry["evals"]),
        cand_cost=[],
        saving=0,
    )
    if stop == -1:
        out["type"] = "N"
        out["expected_evals"] = int(entry["evals"])
        for p in range(n_pos):
            if p in traces:
                c = first_went_on(traces[p])
                run = (
                    (int(c[0]) + 1)
                    if (c is not None and p < n_pos - 1)
                    else int(traces[p]["evals"])
                )
                out["cand_cost"].append(CMA_EVALS + run)
            else:
                out["cand_cost"].append(CMA_EVALS)
        return out
    if stop in paused:
        out["type"] = "C"
        out["saving"] = int(sum(v for p, v in paused.items() if p < stop))
        out["expected_evals"] = None
        return out
    earlier = {p: v for p, v in paused.items() if p < stop}
    out["type"] = "A" if earlier else "-"
    out["saving"] = int(sum(earlier.values()))
    out["expected_evals"] = int(entry["evals"]) - out["saving"]
    return out


def size_sample(steps: List[Dict[str, Any]], total: int) -> Dict[str, Any]:
    """The plan's table row of one sample from the walks of its seed steps (``walk_step`` rows)."""
    a_sav = c_sav = 0
    n = dict(A=0, C=0, N=0)
    c_k: List[int] = []
    cand: List[int] = []
    for w in steps:
        if w["type"] in n:
            n[w["type"]] += 1
        if w["type"] == "A":
            a_sav += w["saving"]
        elif w["type"] == "C":
            c_sav += w["saving"]
            c_k.append(w["stop"])
        elif w["type"] == "N":
            cand.extend(w["cand_cost"])
    mean_c = sum(cand) / max(1, len(cand))
    extra_hi = sum((N_CAND - 1 - k) * mean_c for k in c_k)
    extra_lo = CMA_EVALS * len(c_k)
    return dict(
        total=int(total),
        a_steps=n["A"],
        a_saving=int(a_sav),
        a_pct=100 * a_sav / total,
        c_steps=n["C"],
        c_earlier_saving=int(c_sav),
        n_steps=n["N"],
        mean_cand=mean_c,
        extra_lo=int(extra_lo),
        extra_hi=float(extra_hi),
        lo=100 * (a_sav + c_sav - extra_hi) / total,
        hi=100 * (a_sav + c_sav - extra_lo) / total,
    )


def matches_plan(sample: str, r: Dict[str, Any]) -> List[str]:
    """The differences between a replay row and the plan's table (empty: reproduced exactly)."""
    t = PLAN_TABLE[sample]
    bad: List[str] = []
    for k in ("total", "a_steps", "c_steps", "n_steps"):
        if int(r[k]) != int(t[k]):
            bad.append(f"{sample} {k}: replay {r[k]} plan {t[k]}")
    for k, nd in (("a_pct", 2), ("lo", 1), ("hi", 1)):
        if f"{r[k]:.{nd}f}" != f"{t[k]:.{nd}f}":
            bad.append(f"{sample} {k}: replay {r[k]:.{nd}f} plan {t[k]:.{nd}f}")
    return bad


# ---- the seed step's state machine (live arm and tests) ------------------------------------------
Fit = Dict[str, Any]


def run_seed_step(
    n_pos: int,
    run_c: Callable[[int], Fit],
    run_m: Callable[[int, bool], Fit],
    resume_m: Callable[[int, Fit], Fit],
    tau: float = TAU,
    defer: bool = True,
) -> Dict[str, Any]:
    """One seed step (interleaved order) with DEF.

    ``run_c(p)``: the CMA fit of candidate position ``p``: a dict with ``hit``, ``cost``,
    ``evals``, ``R``. ``run_m(p, allow_pause)``: its MC refinement; with ``allow_pause`` the first
    check that goes on may pause it, then the dict has ``paused`` True, ``hit`` and ``cost`` of the
    running best, ``evals`` paid so far (the check included). ``resume_m(p, fit)``: continues a
    paused refinement under J5's rules and returns the complete fit (``evals`` as J5 counts it).

    A paused refinement is never allowed at the last candidate (``p == n_pos - 1``), and with
    ``defer`` False nothing is paused (the step is J5's). Returns ``kept`` (per tried position:
    ``e`` the end, ``evals`` of the position's fits, ``paused`` True if never resumed, ``c``,
    ``m``), ``stop`` (position, or None), ``calls`` (the order the fits ran), ``paused`` (positions
    paused when their check went on) and ``resumed`` (positions resumed), ``found`` (``tried``,
    ``resume`` or ``none``)."""
    c_of: Dict[int, Fit] = {}
    m_of: Dict[int, Fit] = {}
    calls: List[Tuple[Any, ...]] = []
    paused: List[int] = []
    resumed: List[int] = []
    stop: Optional[int] = None
    found = "none"
    end: Dict[int, Fit] = {}
    for p in range(n_pos):
        c = run_c(p)
        calls.append(("C", p))
        c_of[p] = c
        e = c
        if c["hit"] < tau:
            allow = defer and p < n_pos - 1
            m = run_m(p, allow)
            calls.append(("M", p, "pause" if m.get("paused") else "run"))
            m_of[p] = m
            if m.get("paused"):
                assert allow, "a paused refinement at a candidate that may not pause"
                paused.append(p)
            if m["cost"] < c["cost"]:
                e = m
        end[p] = e
        if not (m_of.get(p) or {}).get("paused") and e["hit"] >= tau:
            stop, found = p, "tried"
            break
    if stop is None:
        for p in paused:
            m = resume_m(p, m_of[p])
            calls.append(("R", p))
            m_of[p] = m
            resumed.append(p)
            c = c_of[p]
            end[p] = m if m["cost"] < c["cost"] else c
            if end[p]["hit"] >= tau:
                stop, found = p, "resume"
                break
    kept = []
    for p in sorted(c_of):
        m = m_of.get(p)
        kept.append(
            dict(
                pos=p,
                e=end[p],
                evals=int(c_of[p]["evals"] + (m["evals"] if m else 0)),
                paused=bool(m and m.get("paused")),
                c=c_of[p],
                m=m,
            )
        )
    return dict(kept=kept, stop=stop, calls=calls, paused=paused, resumed=resumed, found=found)
