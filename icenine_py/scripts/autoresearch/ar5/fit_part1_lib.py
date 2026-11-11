"""AR5 step 2, Part 1: read the replay rows and count (shared by the pick and the summary).

``load(sample)`` reads every ``cache/replay_<sample>_*of*.json.gz`` shard (sorted by step) and
returns the rows; ``arrays(rows, arm)`` gives, per step, the evaluations, the end hit ratio, the
end error, and the two flags of the plan's rules: ``passed`` (seed test: end hit ratio >= 0.80)
and ``right`` (end within 1 degree of the truth). A step without candidates has no arms in its
row: it fails both flags for every arm and costs 0 (as in the log of S8).
"""

import glob
import gzip
import json
from typing import Any, Dict, List

import fit_common as FC
import fit_rules as FR
import numpy as np

RIGHT_DEG = 1.0  # the plan: "within 1 degree of the truth"


def load(sample: str) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    files = sorted(glob.glob(str(FC.CACHE / f"replay_{sample}_*of*.json.gz")))
    files = [f for f in files if "_limit" not in f]
    for f in files:
        with gzip.open(f, "rt") as fh:
            rows += json.load(fh)
    rows.sort(key=lambda r: r["step"])
    return rows


def arrays(rows: List[Dict[str, Any]], arm: str) -> Dict[str, np.ndarray]:
    n = len(rows)
    evals = np.zeros(n, np.int64)
    hit = np.zeros(n)
    err = np.full(n, np.inf)
    for i, r in enumerate(rows):
        d = r["arms"].get(arm)
        if d is None:
            continue
        evals[i], hit[i], err[i] = d["evals"], d["hit"], d["err"]
    passed = np.array([FR.seed_test(h) for h in hit], dtype=bool)
    passed &= np.array([bool(r["arms"]) for r in rows], dtype=bool)
    right = err <= RIGHT_DEG
    return dict(evals=evals, hit=hit, err=err, passed=passed, right=right)


def counts(rows: List[Dict[str, Any]]) -> Dict[str, Dict[str, int]]:
    out: Dict[str, Dict[str, int]] = {}
    for arm in FC.ARMS:
        a = arrays(rows, arm)
        out[arm] = dict(
            n_pass=int(a["passed"].sum()),
            n_right=int(a["right"].sum()),
            evals=int(a["evals"].sum()),
        )
    return out


def err_quantiles(rows: List[Dict[str, Any]]) -> Dict[str, Dict[str, float]]:
    """Median and 90th percentile of the end error (degrees) of each arm over the steps that have
    candidates; ``n`` is the number of those steps. ``C0`` is the error of the candidate as costed
    (the start of every refinement)."""
    out: Dict[str, Dict[str, float]] = {}
    for arm in FC.ARMS:
        e = arrays(rows, arm)["err"]
        e = e[np.isfinite(e)]
        out[arm] = dict(
            n=int(len(e)),
            median=float(np.median(e)),
            p90=float(np.percentile(e, 90)),
            within_0_2=float((e <= 0.2).mean()),
        )
    return out


def paired(rows: List[Dict[str, Any]], arm: str, key: str) -> Dict[str, int]:
    """Paired split of a flag (``passed`` or ``right``): a = M, b = the arm."""
    m = arrays(rows, "M")[key]
    x = arrays(rows, arm)[key]
    d = FR.paired_split(m, x)
    return dict(both=d["both"], m_only=d["only_a"], arm_only=d["only_b"], neither=d["neither"])
