"""AR5 step 8, phase A step 5: the numbers of the finished runs of one sample.

  RESEED_SAMPLE=<held|mid|grad|big> uv run python scripts/autoresearch/ar5/order_stats.py

For each of the arms ``A`` and ``Y`` (and ``Yl`` on ``big``) it runs step 2's ``fit_part2_stats``
(evaluations by stage with the pass, grid fallbacks, wrong voxels before and after the pass with
Wilson 95% intervals, median error, fixed and broken against the stored ``S8``, the house rule),
then reads the seed-step logs (``e1_log.json.gz``; the MC traces are NOT read) and writes
``benchmarks/autoresearch/ar5/order/stats_<sample>.json``:

* per arm: seed-step evaluations, seed steps, CMA fits, MC refinements and their evaluations, steps
  that ran MC, the validity counts (a step of ``Y`` that ran MC although a CMA end had reached 0.95;
  a CMA fit of ``Y`` after one that reached 0.95; MC in ``A`` after a CMA end at 0.95);
* ``Y`` against ``A`` (and, on ``big``, ``Yl`` against ``Y`` and ``A``): fixed and broken by voxel
  before and after the pass (fixed = wrong in the first, right in the second), the seed steps
  matched by (voxel, occurrence), those whose taken end differs and how many of the differing ends
  are wrong (more than 1 degree) in each arm, and prediction 5: the fits both arms ran at the same
  voxel and candidate (kind, hit ratio, cost, evaluations, end orientation compared exactly).
Truth is used for scoring only. The error arrays stay in the cache (untracked).
"""

import gzip
import json
import sys
from typing import Any, Dict, List, Tuple

import order_common as OC  # noqa: E402  (first)

import fit_common as FC  # noqa: E402
import fit_part2_stats as P2  # noqa: E402
import fit_rules as FR  # noqa: E402
import numpy as np  # noqa: E402

import gradient_common as G  # noqa: E402

import stats as ST  # noqa: E402

_HELD_AR3 = OC.WT / "agent-af8961874e131bda2/icenine_py/benchmarks/autoresearch/ar3/held/cache"
FC.STORED_S8.setdefault("held", _HELD_AR3 / "S8_held")
FC.STORED_PASS.setdefault("held", _HELD_AR3 / "passes_held")

SAMPLE = FC.SAMPLE
STATS_DIR = OC.ORDER_CACHE / "stats"
ARMS = ["A", "Y"] + (["Yl"] if SAMPLE == "big" else [])


def read_log(arm: str) -> List[Dict[str, Any]]:
    with gzip.open(OC.ORDER_CACHE / f"{arm}_{SAMPLE}" / "e1_log.json.gz", "rt") as f:
        return json.load(f)


def kept_ends(e: Dict[str, Any]) -> List[Dict[str, Any]]:
    """The kept end of each refined candidate position (lower cost of its fits; CMA on a tie)."""
    by: Dict[int, Dict[str, Any]] = {}
    for pos, j, kind, hit, cost, ev, R in e["fits"]:
        d = dict(pos=pos, j=j, kind=kind, hit=hit, cost=cost, evals=ev, R=R)
        if pos not in by or (kind == "M" and cost < by[pos]["cost"]):
            by[pos] = d
    return [by[p] for p in sorted(by)]


def taken_end(e: Dict[str, Any]) -> Dict[str, Any]:
    ends = kept_ends(e)
    if e["stop"] >= 0:
        return ends[e["stop"]]
    return ends[int(np.argmin([x["cost"] for x in ends]))]


def log_counts(log: List[Dict[str, Any]], order: str) -> Dict[str, Any]:
    tau = FR.TAU
    c = dict(
        seed_steps=len(log),
        steps_with_fits=0,
        cma_fits=0,
        mc_fits=0,
        cma_evals=0,
        mc_evals=0,
        steps_with_mc=0,
        bad_mc_after_cma_pass=0,
        bad_cma_after_cma_pass=0,
        steps_taken_cma_end=0,
        steps_taken_mc_end=0,
    )
    for e in log:
        if not e["fits"]:
            continue
        c["steps_with_fits"] += 1
        cs = [f for f in e["fits"] if f[2] == "C"]
        ms = [f for f in e["fits"] if f[2] == "M"]
        c["cma_fits"] += len(cs)
        c["mc_fits"] += len(ms)
        c["cma_evals"] += sum(f[5] for f in cs)
        c["mc_evals"] += sum(f[5] for f in ms)
        c["steps_with_mc"] += int(bool(ms))
        passed = [f[0] for f in cs if f[3] >= tau]
        if passed:  # a CMA end reached 0.95 at position p
            p = min(passed)
            if order == "Y":
                c["bad_mc_after_cma_pass"] += len(ms)  # Y: no MC at all
                c["bad_cma_after_cma_pass"] += sum(1 for f in cs if f[0] > p)
            else:
                c["bad_mc_after_cma_pass"] += sum(1 for f in ms if f[0] >= p)
                c["bad_cma_after_cma_pass"] += sum(1 for f in cs if f[0] > p)
        t = taken_end(e)
        c["steps_taken_mc_end" if t["kind"] == "M" else "steps_taken_cma_end"] += 1
    return c


def index(log: List[Dict[str, Any]]) -> Dict[Tuple[int, int], Dict[str, Any]]:
    seen: Dict[int, int] = {}
    out: Dict[Tuple[int, int], Dict[str, Any]] = {}
    for e in log:
        v = int(e["voxel"])
        k = seen.get(v, 0)
        seen[v] = k + 1
        out[(v, k)] = e
    return out


def seed_err(T: Any, v: int, R: List[float]) -> float:
    Rm = np.asarray(R, np.float64).reshape(3, 3)
    return float(T.err(Rm[None], [v])[0])


def compare_logs(la: List[Dict[str, Any]], lb: List[Dict[str, Any]], T: Any) -> Dict[str, Any]:
    ia, ib = index(la), index(lb)
    keys = sorted(set(ia) & set(ib))
    differs = wrong_a = wrong_b = 0
    fits_cmp = fits_diff = 0
    for k in keys:
        ea, eb = ia[k], ib[k]
        if not ea["fits"] or not eb["fits"]:
            continue
        fa = {(f[1], f[2]): f for f in ea["fits"]}
        fb = {(f[1], f[2]): f for f in eb["fits"]}
        for key in set(fa) & set(fb):
            fits_cmp += 1
            x, y = fa[key], fb[key]
            if (x[3], x[4], x[5], x[6]) != (y[3], y[4], y[5], y[6]):
                fits_diff += 1
        ta, tb = taken_end(ea), taken_end(eb)
        if (ta["j"], ta["kind"]) != (tb["j"], tb["kind"]):
            differs += 1
            wrong_a += int(seed_err(T, k[0], ta["R"]) > G.WRONG)
            wrong_b += int(seed_err(T, k[0], tb["R"]) > G.WRONG)
    return dict(
        steps_a=len(la),
        steps_b=len(lb),
        matched=len(keys),
        only_a=len(set(ia) - set(ib)),
        only_b=len(set(ib) - set(ia)),
        taken_differs=differs,
        differing_wrong_first=wrong_a,
        differing_wrong_second=wrong_b,
        shared_fits_compared=fits_cmp,
        shared_fits_different=fits_diff,
    )


def pair_voxels(ea: Dict[str, np.ndarray], eb: Dict[str, np.ndarray]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for tag in ("before", "after"):
        wa, wb = ea[tag] > G.WRONG, eb[tag] > G.WRONG
        fixed, broken = ST.paired_discordant(wa, wb)  # fixed: wrong in a, right in b
        out[tag] = dict(fixed=int(fixed), broken=int(broken), both_wrong=int((wa & wb).sum()))
    return out


def main() -> None:
    STATS_DIR.mkdir(parents=True, exist_ok=True)
    FC.FIT = STATS_DIR  # fit_part2_stats writes its json and error arrays here (untracked)
    per: Dict[str, Any] = {}
    for arm in ARMS:
        sys.argv = ["fit_part2_stats.py", "--name", arm]
        P2.main()
        d = json.loads((STATS_DIR / f"part2_{arm}_{SAMPLE}.json").read_text())
        assert not d["unfinished"], arm
        order = "Y" if arm.startswith("Y") else "A"
        d["log"] = log_counts(read_log(arm), order)
        d["order_json"] = json.loads(
            (OC.ORDER_CACHE / f"{arm}_{SAMPLE}" / "order.json").read_text()
        )
        d["order_json"].pop("state", None)
        per[arm] = d
    T = G.Truth()
    errs = {a: dict(np.load(STATS_DIR / f"part2_{a}_{SAMPLE}_err.npz")) for a in ARMS}
    logs = {a: read_log(a) for a in ARMS}
    pairs: Dict[str, Any] = {}
    pairs["Y_vs_A"] = dict(
        voxels=pair_voxels(errs["A"], errs["Y"]), seed_steps=compare_logs(logs["A"], logs["Y"], T)
    )
    if "Yl" in ARMS:
        pairs["Yl_vs_Y"] = dict(
            voxels=pair_voxels(errs["Y"], errs["Yl"]),
            seed_steps=compare_logs(logs["Y"], logs["Yl"], T),
        )
        pairs["Yl_vs_A"] = dict(voxels=pair_voxels(errs["A"], errs["Yl"]))
    out = dict(sample=SAMPLE, arms=per, pairs=pairs)
    (OC.ORDER / f"stats_{SAMPLE}.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({"sample": SAMPLE, "pairs": pairs["Y_vs_A"]}, indent=1))


if __name__ == "__main__":
    main()
