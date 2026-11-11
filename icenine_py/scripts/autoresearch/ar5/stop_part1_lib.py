"""AR5 step 1, part 1: a rule of the try test applied to the stored tries of one sample.

Reads ``benchmarks/autoresearch/ar5/starts_<sample>.npz`` (written by ``stop_part1_starts.py``).
The rule functions see only hit ratios, costs and voxel numbers; the errors in the file are used
by ``evaluate`` to count wrong voxels (scoring), never by a rule.
"""

from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import stop_rules as SR

WRONG = 1.0
AR5 = Path(__file__).resolve().parents[3] / "benchmarks" / "autoresearch" / "ar5"


def load(sample: str) -> Dict[str, np.ndarray]:
    z = np.load(AR5 / f"starts_{sample}.npz")
    return {k: z[k] for k in z.files}


def settings(d_pick: float) -> List[Tuple[str, str, Any, Any]]:
    """(name, kind, d, m) of every setting in the plan's table; P3 uses ``d_pick``."""
    out: List[Tuple[str, str, Any, Any]] = [("full", "none", None, None)]
    out += [(f"P1 d={d:g}", "P1", d, None) for d in SR.D_VALUES]
    out += [(f"P2 m={m}", "P2", None, m) for m in SR.M_VALUES]
    out += [(f"P3 d={d_pick:g} m={m}", "P3", d_pick, m) for m in SR.M_VALUES]
    return out


def evaluate(S: Dict[str, np.ndarray], kind: str, d: Any, m: Any) -> Dict[str, Any]:
    """Tries kept, evaluations, wrong voxels after the pass and the split against the full pass."""
    voxel = S["voxel"]
    old_hit = S["old_hit"].astype(np.float64)
    old_cost = S["old_cost"].astype(np.float64)
    s_hit, s_cost = S["start_hit"].astype(np.float64), S["start_cost"].astype(np.float64)
    e_hit, e_cost = S["end_hit"].astype(np.float64), S["end_cost"].astype(np.float64)
    if kind == "none":
        keep = np.ones(len(voxel), bool)
    else:
        keep = SR.rule_keep(kind, voxel, s_hit, s_cost, old_hit[voxel], d=d, m=m)
    ch = SR.choose_try(voxel, old_hit, old_cost, e_hit, e_cost, keep)
    wrong_before = S["err_before"] > WRONG
    err = S["err_before"].astype(np.float64).copy()
    took = ch >= 0
    err[took] = S["end_err"][ch[took]]
    wrong_after = err > WRONG
    full_ch = SR.choose_try(voxel, old_hit, old_cost, e_hit, e_cost, np.ones(len(voxel), bool))
    err_f = S["err_before"].astype(np.float64).copy()
    err_f[full_ch >= 0] = S["end_err"][full_ch[full_ch >= 0]]
    wrong_full = err_f > WRONG
    pv = SR.paired_vs_full(wrong_before, wrong_full, wrong_after)
    ev = SR.pass_evals(S["evals"], keep) if kind != "none" else int(S["evals"].sum())
    return dict(
        tries=int(len(voxel)),
        tries_kept=int(keep.sum()),
        evals=int(ev),
        wrong_before=int(wrong_before.sum()),
        wrong_after=int(wrong_after.sum()),
        net_gain=SR.net_gain(wrong_before, wrong_after),
        voxels_changed=int(took.sum()),
        **pv,
    )
