"""AR5 step 17 (task ar5q): the rule family for an early skip of the BFS local optimizations of arm
``J5``, its replay on a recording, its selection and its decision rules (pure; numpy and the
standard library only).

Terms (plan "AR5 step 17"). A **fit** is one BFS local optimization (CMA, 60 evaluations plus one
final hit evaluation = 61) from an inherited start: a **neighbour fit** (a voxel's first) or a
**revisit** (a refit of a rejected voxel from a later neighbour). **Best**: the expansion's running
best hit ratio. A fit is accepted if its hit ratio is over 0.9 times the best. A **wrong skip** is a
skipped or stopped fit that ``J5`` accepted.

A skipped or stopped fit is a REJECTED fit (the voxel becomes or stays a refit voxel; a skipped
neighbour fit holds its start, a skipped revisit keeps the old end). Every check is one counted
evaluation, except rule ``D`` (none).

* ``S(t, scope)``: before the fit, check the start; skip if its hit ratio is under ``t`` times the
  best. One evaluation per fit in scope (skipped or not). A skip saves the fit's evaluations.
* ``C(t, scope)``: after 20 evaluations, check the running best; stop if under ``t`` times the
  best. One evaluation per fit in scope. A stop saves the fit's evaluations minus 20.
* ``D(d)``: skip a revisit whose start lies within ``d`` degrees of a start already rejected at
  that voxel since its last reset (no evaluation). A reset is a revisit counter that did not go up
  (a reseed reset the voxel's revisit budget).

``t`` in {0.3, 0.4, 0.5, 0.6}; scope revisit, neighbour or both; ``d`` in {0.1, 0.25, 0.5}: 27
rules. Nothing here is tuned after a recording is read.

A fit row of the recording (``fitskip_record``) holds at least: ``voxel``, ``kind`` (``neighbour``
or ``revisit``), ``best`` (the running best hit ratio before the fit), ``start_R`` (9 numbers),
``start_hit``, ``hit20`` (hit ratio of the running best after 20 evaluations), ``evals`` (counted),
``accepted`` (bool), ``n_revisits`` (the voxel's revisit count at the fit, 0 for a neighbour fit)
and, after ``fitskip_final``, ``final_before`` (the end is the voxel's final answer before the
pass).
"""

from fractions import Fraction
from typing import Any, Dict, Iterable, List, NamedTuple, Optional, Sequence, Tuple

import numpy as np

TS = (0.3, 0.4, 0.5, 0.6)
SCOPES = ("revisit", "neighbour", "both")
DS = (0.1, 0.25, 0.5)
CHECK_AT = 20  # evaluations before rule C's check
FORM_RANK = {"D": 0, "S": 1, "C": 2}  # the tie rule: D before S before C
SCOPE_RANK = {"revisit": 0, "neighbour": 1, "both": 2}
WRONG_SKIP_MAX_SHARE = Fraction(1, 1000)  # selection (a): at most 0.1% of the sample's fits
FINAL_CHANGED_MAX = 2  # selection (b): at most 2 final answers before the pass change
ACCEPT = 0.9

RULE = (
    "Early skip of the BFS local optimizations of arm J5. A skipped or stopped fit is a rejected "
    "fit (the voxel becomes or stays a refit voxel; a skipped neighbour fit holds its start, a "
    "skipped revisit keeps the old end). Every check is one counted evaluation. S(t, scope): "
    "before the fit, check the start; skip if its hit ratio is under t times the best. C(t, "
    "scope): after 20 evaluations, check the running best; stop if under t times the best. D(d): "
    "skip a revisit whose start lies within d degrees of a start already rejected at that voxel "
    "since its last reset (no evaluation). t in {0.3, 0.4, 0.5, 0.6}; scope revisits, neighbour "
    "fits or both; d in {0.1, 0.25, 0.5}: 27 rules."
)
SELECTION = (
    "A rule is admissible if on each of the five samples (a) wrong skips are at most 0.1% of that "
    "sample's fits and (b) at most 2 final answers before the pass change. Among admissible rules "
    "take the largest sum over the five samples of the net saving as a share of J5 (saved "
    "evaluations of skipped fits that J5 rejected, minus every check; wrong skips earn nothing). "
    "Ties: fewer wrong skips in all, then D before S before C, then the smaller t or d."
)
TIE_COMPLETION = (
    "Not stated by the plan, added here before any recording is read: when the three stated tie "
    "rules leave rules equal, the scope order revisit, neighbour, both decides (an equal net "
    "saving to the last evaluation is already unlikely)."
)

# stored counts of arm J5 (plan table; stop rule): total with the pass; fits and rejected fits
EXPECT: Dict[str, Dict[str, int]] = {
    "held3": dict(
        total=23_503_385, revisits=48_181, rej_revisits=38_314, nbr=24_277, nbr_rej=12_953
    ),
    "held": dict(
        total=22_939_812, revisits=52_638, rej_revisits=42_756, nbr=24_273, nbr_rej=13_295
    ),
    "mid": dict(total=10_589_493, revisits=27_618, rej_revisits=18_778, nbr=24_409, nbr_rej=10_297),
    "grad": dict(total=9_041_686, revisits=29_176, rej_revisits=15_416, nbr=24_222, nbr_rej=14_693),
    "big": dict(total=3_567_098, revisits=8_010, rej_revisits=2_873, nbr=24_533, nbr_rej=5_424),
}
SAMPLES = tuple(EXPECT)  # the five disclosed samples, in the plan's table order
FIT_EVALS = 61

GATE: Dict[str, Any] = dict(
    g0="prediction 1: J5r gives J5's map before the pass at every voxel, its fit counts (table) "
    "and evaluations",
    g1="at least one admissible rule",
    g2_net_saving_min=dict(held=0.04, mid=0.04, held3=0.03, grad=0.03, big=0.01),
    g3="exact replay: on mid the live J5C run, up to its first wrong skip in run order, equals "
    "the replay fit by fit (skips, ends, evaluations), and over the whole run if it has none",
    fails="stop, report, no tuning, no second pick, no whole layer",
)
JUDGE: Dict[str, Any] = dict(
    r1_total_fewer=dict(held4=0.03, held3=0.03, held=0.03, mid=0.03, grad=0.03, big=0.01),
    r2_broken_minus_fixed_max=5,
    r3_median_error_higher_max_deg=0.001,
    r4_wrong_seeds_more_max=1,
    r4_grid_within=1,
    r4_reseeds_more_max_share=0.02,
    r5_vs_s8=dict(fewer=0.20, broken_minus_fixed_max=5, median_error_higher_max_deg=0.005),
    r5_note="if J5 fails (5) on a sample, J5C is judged there by (1) to (4)",
    pass_samples=("held4", "held3", "held", "mid", "grad", "big"),
    step_passes="J5C passes on held4 and on the five disclosed samples; fails on held4: no key",
)
PREDICTIONS = (
    "1. J5r gives J5's map before the pass at every voxel, its fit counts (table) and "
    "evaluations. 2. The chosen rule is S(0.3, both) or a D rule (about even odds); no C rule is "
    "chosen. 3. Net saving of the chosen rule on the recording: held 5% to 8%, held3 4.5% to "
    "7.5%, mid 5% to 9%, grad 7% to 10%, big 2% to 4% of J5; wrong skips 5 to 40 per sample. 4. "
    "J5C against J5, with the pass: held4 4% to 9% fewer; the disclosed samples within 1.5 points "
    "of their recording saving; broken minus fixed -2 to 3; median error up at most 0.0005. 5. "
    "mid passes rule 5 (about 3 in 5)."
)


class Rule(NamedTuple):
    form: str  # "S", "C" or "D"
    t: Optional[float]  # S, C
    scope: Optional[str]  # S, C
    d: Optional[float]  # D

    @property
    def name(self) -> str:
        if self.form == "D":
            return f"D({self.d})"
        return f"{self.form}({self.t},{self.scope})"

    def sort_key(self) -> Tuple[Any, ...]:
        """The tie order after 'fewer wrong skips in all': D, S, C; smaller t or d; scope."""
        return (
            FORM_RANK[self.form],
            self.t if self.t is not None else self.d,
            SCOPE_RANK.get(self.scope or "", 0),
        )


def family() -> List[Rule]:
    """The 27 rules: D(d) x3, S(t, scope) x12, C(t, scope) x12."""
    out = [Rule("D", None, None, d) for d in DS]
    for form in ("S", "C"):
        out += [Rule(form, t, sc, None) for t in TS for sc in SCOPES]
    assert len(out) == 27
    return out


def in_scope(kind: str, scope: Optional[str]) -> bool:
    return scope == "both" or scope == kind


def accepted_test(hit: float, best_before: float) -> bool:
    """The caller's test (``_expand``): ``best = max(hit, best_before)``; accepted if best > 0 and
    hit / best > 0.9."""
    best = max(hit, best_before)
    return bool(best > 0 and (hit / best) > ACCEPT)


def below_line(hit: float, best: float, t: float) -> bool:
    """True when ``hit`` is under ``t`` times ``best`` (skip / stop only strictly below)."""
    return bool(hit < t * best)


def mis_deg(R1: Any, R2: Any, sym: np.ndarray) -> float:
    """Smallest angle in degrees between ``R1`` and ``R2 @ s`` over the proper symmetry rotations
    ``sym`` (S, 3, 3)."""
    a = np.asarray(R1, np.float64).reshape(3, 3)
    b = np.asarray(R2, np.float64).reshape(3, 3)
    m = a.T @ b  # (3, 3)
    tr = np.einsum("ij,sji->s", m, np.asarray(sym, np.float64))  # trace(m @ s)
    c = np.clip((tr - 1.0) / 2.0, -1.0, 1.0)
    return float(np.degrees(np.arccos(c.max())))


class VoxelHistory:
    """Per voxel: the starts of fits rejected since the last reset, and the last revisit count.
    The same object serves the replay and the live rule ``D``."""

    def __init__(self) -> None:
        self.rejected: Dict[int, List[np.ndarray]] = {}
        self.last_rev: Dict[int, int] = {}

    def begin(self, voxel: int, kind: str, n_revisits: int) -> None:
        """Called before a fit: a revisit counter that did not go up is a reset."""
        if kind == "revisit":
            if n_revisits <= self.last_rev.get(voxel, 0):
                self.rejected[voxel] = []
            self.last_rev[voxel] = n_revisits

    def near(self, voxel: int, start: Any, d: float, sym: np.ndarray) -> bool:
        return any(mis_deg(start, r, sym) < d for r in self.rejected.get(voxel, []))

    def reject(self, voxel: int, start: Any) -> None:
        self.rejected.setdefault(voxel, []).append(np.asarray(start, np.float64).reshape(3, 3))


def replay(
    rule: Rule, fits: Iterable[Dict[str, Any]], sym: Optional[np.ndarray] = None
) -> Dict[str, Any]:
    """The rule on a recording (fits in run order). Returns counts: ``n_fits``, ``skipped``,
    ``wrong`` (skipped fits that J5 accepted), ``saved`` (evaluations of skipped fits that J5
    rejected; for C the evaluations beyond 20), ``checks``, ``net`` (= saved - checks),
    ``final_changed`` (skipped fits whose end is the voxel's final answer before the pass).
    The replay is exact for skips of fits J5 rejected (the plan); a wrong skip is counted and
    earns nothing, later fits are taken as J5 made them."""
    hist = VoxelHistory()
    n = skipped = wrong = saved = checks = final_changed = 0
    for f in fits:
        n += 1
        kind = str(f["kind"])
        v = int(f["voxel"])
        skip = False
        if rule.form == "D":
            assert sym is not None, "rule D needs the symmetry rotations"
            hist.begin(v, kind, int(f.get("n_revisits", 0)))
            skip = kind == "revisit" and hist.near(v, f["start_R"], float(rule.d or 0), sym)
            if skip or not bool(f["accepted"]):
                hist.reject(v, f["start_R"])
        else:
            assert rule.t is not None
            if in_scope(kind, rule.scope):
                checks += 1
                hit = float(f["start_hit"] if rule.form == "S" else f["hit20"])
                skip = below_line(hit, float(f["best"]), rule.t)
        if skip:
            skipped += 1
            if bool(f["accepted"]):
                wrong += 1
            else:
                saved += int(f["evals"]) - (CHECK_AT if rule.form == "C" else 0)
                if bool(f.get("final_before", False)):
                    final_changed += 1
    return dict(
        rule=rule.name,
        n_fits=n,
        skipped=skipped,
        wrong=wrong,
        saved=saved,
        checks=checks,
        net=saved - checks,
        final_changed=final_changed,
    )


def admissible(per_sample: Dict[str, Dict[str, Any]]) -> bool:
    """(a) wrong skips at most 0.1% of that sample's fits and (b) at most 2 final answers
    before the pass change, on every sample."""
    for r in per_sample.values():
        if Fraction(int(r["wrong"]), 1) > WRONG_SKIP_MAX_SHARE * int(r["n_fits"]):
            return False
        if int(r["final_changed"]) > FINAL_CHANGED_MAX:
            return False
    return True


def net_share_sum(per_sample: Dict[str, Dict[str, Any]]) -> Fraction:
    """The sum over the samples of net / (J5 total with the pass), exactly."""
    return sum(
        (Fraction(int(r["net"]), EXPECT[s]["total"]) for s, r in per_sample.items()), Fraction(0)
    )


def select(results: Dict[str, Dict[str, Dict[str, Any]]]) -> Dict[str, Any]:
    """``results[rule name][sample]`` = ``replay`` output. Returns the pick (or ``None`` when no
    rule is admissible) and the full ranking of the admissible rules."""
    by_name = {r.name: r for r in family()}
    ok = []
    for name, per in results.items():
        if set(per) != set(SAMPLES):
            raise ValueError(f"{name}: need all five samples, got {sorted(per)}")
        if admissible(per):
            wrong_all = sum(int(r["wrong"]) for r in per.values())
            ok.append((name, net_share_sum(per), wrong_all))
    ok.sort(key=lambda x: (-x[1], x[2], by_name[x[0]].sort_key()))
    return dict(
        pick=ok[0][0] if ok else None,
        admissible=[dict(rule=n, net_share_sum=float(s), wrong_all=w) for n, s, w in ok],
        n_admissible=len(ok),
    )


def gate_g2(pick: Dict[str, Dict[str, Any]]) -> Dict[str, bool]:
    """g2 per sample: the chosen rule's net saving as a share of J5 at least the plan's number."""
    return {
        s: bool(
            Fraction(int(pick[s]["net"]), EXPECT[s]["total"])
            >= Fraction(str(GATE["g2_net_saving_min"][s]))
        )
        for s in SAMPLES
    }


def counts_ok(s: str, stats: Dict[str, Any], total_with_pass: int) -> Dict[str, Any]:
    """Prediction 1, the counts part: the run's ``stats`` (``run.json``) and the total with the
    pass against the table."""
    e = EXPECT[s]
    got = dict(
        total=int(total_with_pass),
        revisits=int(stats["n_revisit_attempted"]),
        rej_revisits=int(stats["n_revisit_rejected"]),
        nbr=int(stats["n_neighbor_fits"]),
        nbr_rej=int(stats["n_neighbor_fits"]) - int(stats["n_neighbor_accepted_first"]),
    )
    return dict(sample=s, expected=e, got=got, ok=bool(got == e))


def rule_by_name(name: str) -> Rule:
    for r in family():
        if r.name == name:
            return r
    raise ValueError(f"unknown rule {name!r}")
