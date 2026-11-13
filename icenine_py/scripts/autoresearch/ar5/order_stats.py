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
  The differing-end counts cover the MATCHED steps only; the steps run in one arm only are
  counted (``only_a``, ``only_b``) and their ends are scored by ``order_trace.py``.
Truth is used for scoring only. The error arrays stay in the cache (untracked).

The log readers (``kept_ends``, ``taken_end``, ``log_counts``, ``index``, ``compare_logs``) live in
``order_logs.py`` (pure, unit-tested); they were moved there unchanged after the runs.
"""

import gzip
import json
import sys
from typing import Any, Dict, List

import order_common as OC  # noqa: E402  (first)

import fit_common as FC  # noqa: E402
import fit_part2_stats as P2  # noqa: E402
import fit_rules as FR  # noqa: E402
import numpy as np  # noqa: E402
import order_logs as OL  # noqa: E402

import gradient_common as G  # noqa: E402

import stats as ST  # noqa: E402

_HELD_AR3 = OC.WT / "agent-af8961874e131bda2/icenine_py/benchmarks/autoresearch/ar3/held/cache"
FC.STORED_S8.setdefault("held", _HELD_AR3 / "S8_held")
FC.STORED_PASS.setdefault("held", _HELD_AR3 / "passes_held")

SAMPLE = FC.SAMPLE
STATS_DIR = OC.ORDER_CACHE / "stats"
ARMS = ["A", "Y"] + (["Yl"] if SAMPLE == "big" else [])
assert OL.TAU == FR.TAU


def read_log(arm: str) -> List[Dict[str, Any]]:
    with gzip.open(OC.ORDER_CACHE / f"{arm}_{SAMPLE}" / "e1_log.json.gz", "rt") as f:
        return json.load(f)


def seed_err(T: Any, v: int, R: List[float]) -> float:
    Rm = np.asarray(R, np.float64).reshape(3, 3)
    return float(T.err(Rm[None], [v])[0])


def compare_logs(la: List[Dict[str, Any]], lb: List[Dict[str, Any]], T: Any) -> Dict[str, Any]:
    return OL.compare_logs(la, lb, lambda v, R: seed_err(T, v, R) > G.WRONG)


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
        d["log"] = OL.log_counts(read_log(arm), order, FR.TAU)
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
