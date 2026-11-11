"""AR5 step 3: the pure rules of the mixed seed refinement (numpy only; no truth, no images).

Mixed seed refinement (``X60``, ``X200``): for each indexed candidate, in S8's order, a CMA fit
(first step 0.05 degree, at most 60 or 200 evaluations). If its end hit ratio is at least 0.95 that
end is kept and the seed step stops (S8's early stop). If not, the stock MC optimizer runs from the
same candidate and the candidate's end is the one of the two with the lower cost (a tie keeps the
CMA end). The early-stop test is then made on the hit ratio of that kept end. If no end reaches
0.95, the lowest cost of up to 5 ends is kept, as in S8. Evaluations count both refinements.

Part 1 builds this from the replay cache (the CMA arm and the MC arm of the same candidates);
the script layer (``mix_common.mixed_refine``) runs it. Tests: ``tests/test_ar5_mix_rules.py``.
"""

from typing import Any, Dict, Optional

import fit_rules as FR

RIGHT_DEG = 1.0  # the plan: "within 1 degree of the truth"


def mixed_candidate(
    c_hit: float,
    c_cost: float,
    c_evals: int,
    m_hit: float,
    m_cost: float,
    m_evals: int,
    tau: float = FR.TAU,
) -> Dict[str, Any]:
    """The mixed refinement of one candidate from its CMA end (``c_*``) and its MC end (``m_*``;
    only read when the CMA end hit ratio is below ``tau``). Returns the kept end's ``end``
    ('C' or 'M'), ``hit``, ``cost``, the evaluations of what ran and whether MC ran."""
    if c_hit >= tau:
        return dict(end="C", hit=c_hit, cost=c_cost, evals=int(c_evals), mc_ran=False)
    if m_cost < c_cost:
        return dict(end="M", hit=m_hit, cost=m_cost, evals=int(c_evals + m_evals), mc_ran=True)
    return dict(end="C", hit=c_hit, cost=c_cost, evals=int(c_evals + m_evals), mc_ran=True)


def mixed_step(
    costed: int,
    c: Optional[Dict[str, Any]],
    m: Optional[Dict[str, Any]],
    tau: float = FR.TAU,
) -> Dict[str, Any]:
    """One seed step with the mixed refinement, from the replay row's arms ``c`` (CMA) and ``m``
    (MC): dicts with lists ``fin_hit``, ``fin_cost``, ``fin_evals`` and ``fin_err`` (one entry per
    refined candidate, in refinement order). ``costed`` is the number of candidates costed (one
    evaluation each). ``c`` and ``m`` are None for a step without candidates.

    Returns ``evals``, ``taken`` (position of the kept end), ``hit``, ``err`` of it, ``stopped``
    (stopped early), ``mc_ran`` (candidates on which MC ran), ``n_refined``, the flags ``passed``
    (seed test) and ``right`` (within 1 degree), and ``missing``: True when the rule needed the MC
    end of a candidate that the replay did not compute (MC had stopped earlier in the replay). In
    that case the CMA fit of that candidate is counted and the step ends there; the caller must
    report such steps."""
    out: Dict[str, Any] = dict(
        evals=0,
        taken=None,
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
    evals = int(costed)
    hits, costs, errs = [], [], []
    stopped = False
    for p in range(len(c["fin_hit"])):
        evals += int(c["fin_evals"][p])
        if c["fin_hit"][p] < tau and p >= len(m["fin_hit"]):
            out["missing"] = True
            break
        mc_p = c["fin_hit"][p] < tau
        k = mixed_candidate(
            c["fin_hit"][p],
            c["fin_cost"][p],
            int(c["fin_evals"][p]),
            m["fin_hit"][p] if mc_p else 0.0,
            m["fin_cost"][p] if mc_p else 0.0,
            int(m["fin_evals"][p]) if mc_p else 0,
            tau,
        )
        evals += k["evals"] - int(c["fin_evals"][p])
        out["mc_ran"] += int(k["mc_ran"])
        src = m if k["end"] == "M" else c
        hits.append(k["hit"])
        costs.append(k["cost"])
        errs.append(src["fin_err"][p])
        if k["hit"] >= tau:
            stopped = True
            break
    if not hits:
        out["evals"] = evals
        return out
    taken, stopped2 = FR.early_stop(hits, costs, tau)
    assert stopped2 == stopped
    out.update(
        evals=evals,
        taken=int(taken),
        hit=float(hits[taken]),
        err=float(errs[taken]),
        stopped=bool(stopped),
        n_refined=len(hits),
    )
    out["passed"] = FR.seed_test(out["hit"])
    out["right"] = bool(out["err"] <= RIGHT_DEG)
    return out
