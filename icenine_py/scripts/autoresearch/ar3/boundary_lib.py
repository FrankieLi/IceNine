"""AR3 step 5: keep rules applied to stored tries, and the skip-rule features.

Nothing here reads the truth: a keep rule sees the strict score (cost, hit ratio, confidence) of
the voxel's old end and of the ends of its tries; the skip features see the old scores of the voxel
and of its neighbours. Ties go to the first option, the old end first.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import boundary_common as B
import numpy as np

RULES = ("cost", "hit", "conf")


class Tries:
    """The stored tries of one pass (``tries.npz``)."""

    def __init__(self, run: Path) -> None:
        z = np.load(run / "tries.npz")
        self.voxel = z["voxel"]
        self.cand = z["cand"]
        self.end = B.orthonormal(z["end"].astype(np.float64))
        self.cost, self.hit, self.conf = z["cost"], z["hit"], z["conf"]
        self.evals = z["evals"]
        self.old_score = z["old_score"]  # (n, 3): cost, hit, conf of every voxel's old end
        self.old_R = z["old_R"].astype(np.float64)
        self.info = __import__("json").loads((run / "pass.json").read_text())
        u, first = np.unique(self.voxel, return_index=True)
        self.vox = u
        self.lo = first
        self.hi = np.append(first[1:], len(self.voxel))


def apply_rule(
    t: Tries, R0: np.ndarray, rule: str, skip: Optional[np.ndarray] = None
) -> Tuple[np.ndarray, np.ndarray]:
    """(new map, took_a_try (n,) bool). ``skip`` (n,) bool: those voxels keep the old end."""
    R = R0.copy()
    took = np.zeros(len(R0), dtype=bool)
    for v, lo, hi in zip(t.vox, t.lo, t.hi):
        if skip is not None and skip[v]:
            continue
        oc, oh, of = t.old_score[v]
        cost = np.concatenate([[oc], t.cost[lo:hi]])
        hit = np.concatenate([[oh], t.hit[lo:hi]])
        conf = np.concatenate([[of], t.conf[lo:hi]])
        if rule == "cost":
            order = np.lexsort((np.arange(len(cost)), cost))
        elif rule == "hit":
            order = np.lexsort((np.arange(len(cost)), cost, -hit))
        elif rule == "conf":
            order = np.lexsort((np.arange(len(cost)), cost, -conf))
        else:
            raise KeyError(rule)
        k = int(order[0])
        if k > 0:
            R[v] = t.end[lo + k - 1]
            took[v] = True
    return R, took


def rel_hit(env: B.Env, hit: np.ndarray) -> np.ndarray:
    """Hit ratio of each voxel over the highest hit ratio among its neighbours (map only)."""
    best = np.zeros(env.n)
    np.maximum.at(best, env.pi, hit[env.pj])
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where(best > 0, hit / best, np.inf)
    return r


def tries_per_voxel(t: Tries, n: int) -> np.ndarray:
    c = np.zeros(n, dtype=np.int64)
    c[t.vox] = t.hi - t.lo
    return c


def counts_after(e_before: np.ndarray, e_after: np.ndarray, sel: np.ndarray) -> Dict[str, int]:
    wa, wb = (e_before > B.WRONG)[sel], (e_after > B.WRONG)[sel]
    return dict(
        n=int(sel.sum()),
        wrong_before=int(wa.sum()),
        wrong_after=int(wb.sum()),
        fixed=int((wa & ~wb).sum()),
        broken=int((~wa & wb).sum()),
        both_wrong=int((wa & wb).sum()),
    )


def groups(env: B.Env) -> List[Tuple[str, np.ndarray]]:
    g: List[Tuple[str, np.ndarray]] = [("all", np.ones(env.n, bool))]
    if env.level is not None:
        g += [(f"level {L:g}", env.level == L) for L in B.G.LEVELS]
    return g
