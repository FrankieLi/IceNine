"""AR5 step 8 (task ar5h): the pure rules of the sweep order of the mixed seed refinement.

Interleaved order (arm ``A``, today's ``SeedIndexCMAFirst 200``, ``mix_rules.mixed_step``): per
candidate a CMA fit, then the MC refinement if the CMA end hit ratio is under ``tau``.

Sweep order (arm ``Y``): CMA fits on the candidates in order until one CMA end reaches ``tau``
(that end is taken and the step stops). If none does, the MC refinement runs on the candidates in
order; the end of a candidate is the lower-cost one of its CMA end and its MC end (a cost tie keeps
the CMA end); the step stops at the first such kept end with hit ratio >= ``tau`` (and takes it);
if none reaches ``tau``, the lowest cost of all ends is taken (the first on a tie).

``sweep_step`` works on a stored replay row (the CMA arm ``c`` and the MC arm ``m``: dicts with
lists ``fin_hit``, ``fin_cost``, ``fin_evals``, ``fin_err``); ``sweep_choose`` is the same rule on
ends produced one at a time (used by the running script layer). No truth, no images, numpy only.
Tests: ``tests/test_ar5_order.py``.
"""

from typing import Any, Callable, Dict, List, Optional, Tuple

import fit_rules as FR
import mix_rules as MR
import numpy as np

TAU = FR.TAU


class Missing(Exception):
    """The rule needed an MC end that the replay did not compute."""


def sweep_choose(
    n_cand: int,
    cma_fit: Callable[[int], Tuple[float, float, int]],
    mc_fit: Callable[[int], Tuple[float, float, int]],
    tau: float = TAU,
) -> Dict[str, Any]:
    """Run the sweep rule over ``n_cand`` candidates (in refinement order). ``cma_fit(p)`` and
    ``mc_fit(p)`` return (hit ratio, cost, evaluations) of candidate ``p``'s CMA end and MC end;
    ``mc_fit`` is called only when no CMA end reached ``tau``. Returns ``taken`` (position),
    ``end`` ('C' or 'M'), ``hit``, ``cost``, ``evals`` (all fits that ran), ``mc_ran`` (number of
    MC refinements), ``n_cma`` (CMA fits that ran), ``stopped`` and ``phase`` ('cma', 'mc' or
    'none')."""
    evals = 0
    c_hit: List[float] = []
    c_cost: List[float] = []
    for p in range(n_cand):
        h, c, e = cma_fit(p)
        evals += int(e)
        c_hit.append(float(h))
        c_cost.append(float(c))
        if h >= tau:
            return dict(
                taken=p,
                end="C",
                hit=float(h),
                cost=float(c),
                evals=evals,
                mc_ran=0,
                n_cma=p + 1,
                stopped=True,
                phase="cma",
            )
    k_hit: List[float] = []
    k_cost: List[float] = []
    k_end: List[str] = []
    for p in range(n_cand):
        h, c, e = mc_fit(p)
        evals += int(e)
        if c < c_cost[p]:
            k_hit.append(float(h))
            k_cost.append(float(c))
            k_end.append("M")
        else:
            k_hit.append(c_hit[p])
            k_cost.append(c_cost[p])
            k_end.append("C")
        if k_hit[-1] >= tau:
            return dict(
                taken=p,
                end=k_end[-1],
                hit=k_hit[-1],
                cost=k_cost[-1],
                evals=evals,
                mc_ran=p + 1,
                n_cma=n_cand,
                stopped=True,
                phase="mc",
            )
    if not k_hit:
        return dict(
            taken=None,
            end="",
            hit=0.0,
            cost=float("inf"),
            evals=evals,
            mc_ran=0,
            n_cma=n_cand,
            stopped=False,
            phase="none",
        )
    j = int(np.argmin(np.asarray(k_cost, dtype=float)))
    return dict(
        taken=j,
        end=k_end[j],
        hit=k_hit[j],
        cost=k_cost[j],
        evals=evals,
        mc_ran=n_cand,
        n_cma=n_cand,
        stopped=False,
        phase="mc",
    )


def sweep_step(
    costed: int,
    c: Optional[Dict[str, Any]],
    m: Optional[Dict[str, Any]],
    tau: float = TAU,
) -> Dict[str, Any]:
    """One seed step in the sweep order from the replay row's arms (``c`` the CMA arm, ``m`` the MC
    arm; None for a step without candidates). Same keys as ``mix_rules.mixed_step``, plus ``end``.
    ``missing`` is True when the rule needed an MC end the replay did not compute; the step is
    then ended there and must be reported."""
    out: Dict[str, Any] = dict(
        evals=0,
        taken=None,
        end="",
        hit=0.0,
        err=float("inf"),
        stopped=False,
        mc_ran=0,
        n_refined=0,
        passed=False,
        right=False,
        missing=False,
    )
    if c is None or m is None:
        return out
    n = len(c["fin_hit"])

    def cma_fit(p: int) -> Tuple[float, float, int]:
        return c["fin_hit"][p], c["fin_cost"][p], int(c["fin_evals"][p])

    def mc_fit(p: int) -> Tuple[float, float, int]:
        if p >= len(m["fin_hit"]):
            raise Missing()
        return m["fin_hit"][p], m["fin_cost"][p], int(m["fin_evals"][p])

    try:
        r = sweep_choose(n, cma_fit, mc_fit, tau)
    except Missing:
        out["evals"] = int(costed)
        out["missing"] = True
        return out
    out["evals"] = int(costed) + r["evals"]
    if r["taken"] is None:
        return out
    src = m if r["end"] == "M" else c
    out.update(
        taken=int(r["taken"]),
        end=r["end"],
        hit=float(r["hit"]),
        err=float(src["fin_err"][r["taken"]]),
        stopped=bool(r["stopped"]),
        mc_ran=int(r["mc_ran"]),
        n_refined=n,
    )
    out["passed"] = FR.seed_test(out["hit"])
    out["right"] = bool(out["err"] <= MR.RIGHT_DEG)
    return out
