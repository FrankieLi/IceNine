"""AR5 step 18 (task ar5r): the 16 features of a BFS local optimization, known before the fit
(pure; numpy and the standard library only). One class serves the replay on the recording ``J5r``
and the live hook, so the two cannot differ in what a feature means.

A **fit** is one BFS local optimization (CMA, 60 evaluations plus one final hit evaluation): a
**neighbour fit** (a voxel's first) or a **revisit** (a refit of a rejected voxel from a later
neighbour). **Best**: the expansion's running best hit ratio. A fit is accepted if its hit ratio is
over 0.9 times the best. Features 4 to 6 come from the check of the start (one counted evaluation).

The 16 features (plan "AR5 step 18", frozen):
 1 kind (1 revisit, 0 neighbour fit)        2 revisit count of the voxel (0 for a neighbour fit)
 3 best                                      4 start hit ratio
 5 start cost                                6 start confidence
 7 start hit ratio over best (0 if best 0)   8 distance to rejected starts (degrees, crystal
 symmetry, capped at 10, 10 if none; rejected starts at this voxel since its last reset, a revisit
 count that did not go up)                   9 parent unknown (1 if the expanding voxel has no kept
 fit end yet)                               10 parent hit (best if unknown)
11 parent cost (-1 if unknown)              12 start hit over parent hit (1 if unknown)
13 own hit: hit ratio of this voxel's kept fit end (-1 if none)
14 fits so far at the voxel                 15 rejected fits so far at the voxel
16 place in the expansion (the fit's index in it, from 0)

Only fit ends count as kept ends: seed, reseed and pass ends do not, and a SKIPPED fit has no fit
end (it never updates a kept end). The recording's rows are fits that were not skipped.
"""

import itertools
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

FEATURES: Tuple[str, ...] = (
    "kind",
    "n_revisits",
    "best",
    "start_hit",
    "start_cost",
    "start_conf",
    "start_over_best",
    "dist_rejected",
    "parent_unknown",
    "parent_hit",
    "parent_cost",
    "start_over_parent",
    "own_hit",
    "fits_so_far",
    "rejected_so_far",
    "place",
)
N_FEATURES = len(FEATURES)
DIST_CAP = 10.0


def cubic_rotations() -> np.ndarray:
    """The 24 proper rotations of the cubic point group (the replay's symmetry; the live hook reads
    the group from the config, and a test checks the two agree for this sample's crystal)."""
    out = []
    for p in itertools.permutations(range(3)):
        for s in itertools.product((1.0, -1.0), repeat=3):
            m = np.zeros((3, 3))
            for i in range(3):
                m[i, p[i]] = s[i]
            if np.linalg.det(m) > 0:
                out.append(m)
    return np.array(out)


def min_misorientation_deg(start: Any, others: Sequence[np.ndarray], sym: np.ndarray) -> float:
    """Smallest angle in degrees from ``start`` to any of ``others`` over the proper rotations
    ``sym`` (same quantity as ``fitskip_rules.mis_deg``); 99 if there is no other start (the
    caller caps at ``DIST_CAP``)."""
    if not len(others):
        return 99.0
    a = np.asarray(start, np.float64).reshape(3, 3)
    b = np.asarray(others, np.float64).reshape(-1, 3, 3)
    m = np.einsum("ji,kjl->kil", a, b)  # a.T @ b_k
    tr = np.einsum("kij,sji->ks", m, np.asarray(sym, np.float64))
    c = np.clip((tr - 1.0) / 2.0, -1.0, 1.0).max()
    return float(np.degrees(np.arccos(c)))


class FeatureState:
    """Per-voxel history for the features. Call ``features`` before a fit (once per fit) and
    ``update`` after it (also for a skipped fit, with ``kept`` False: it has no fit end)."""

    def __init__(self, sym: np.ndarray) -> None:
        self.sym = np.asarray(sym, np.float64)
        self.kept: Dict[int, Tuple[float, float]] = {}  # voxel -> (end hit, end cost)
        self.rejected: Dict[int, List[np.ndarray]] = {}
        self.last_rev: Dict[int, int] = {}
        self.nfits: Dict[int, int] = {}
        self.nrej: Dict[int, int] = {}

    def features(
        self,
        voxel: int,
        kind: str,
        n_revisits: int,
        best: float,
        start_R: Any,
        start_hit: float,
        start_cost: float,
        start_conf: float,
        parent: int,
        pos: int,
    ) -> List[float]:
        is_rev = kind == "revisit"
        if is_rev:
            if n_revisits <= self.last_rev.get(voxel, 0):
                self.rejected[voxel] = []  # a reset: the revisit counter did not go up
            self.last_rev[voxel] = n_revisits
        R = np.asarray(start_R, np.float64).reshape(3, 3)
        mis = min_misorientation_deg(R, self.rejected.get(voxel, []), self.sym)
        par = self.kept.get(int(parent))
        own = self.kept.get(voxel)
        return [
            float(is_rev),
            float(n_revisits),
            float(best),
            float(start_hit),
            float(start_cost),
            float(start_conf),
            float(start_hit) / best if best > 0 else 0.0,
            min(mis, DIST_CAP),
            float(par is None),
            par[0] if par else float(best),
            par[1] if par else -1.0,
            (float(start_hit) / par[0]) if par and par[0] > 0 else 1.0,
            own[0] if own else -1.0,
            float(self.nfits.get(voxel, 0)),
            float(self.nrej.get(voxel, 0)),
            float(pos),
        ]

    def update(
        self,
        voxel: int,
        start_R: Any,
        accepted: bool,
        kept: bool,
        end_hit: Optional[float],
        end_cost: Optional[float],
    ) -> None:
        self.nfits[voxel] = self.nfits.get(voxel, 0) + 1
        if not accepted:
            self.nrej[voxel] = self.nrej.get(voxel, 0) + 1
            self.rejected.setdefault(voxel, []).append(
                np.asarray(start_R, np.float64).reshape(3, 3)
            )
        if kept and end_hit is not None and end_cost is not None:
            self.kept[voxel] = (float(end_hit), float(end_cost))


def rows_to_arrays(
    rows: Sequence[Dict[str, Any]], sym: Optional[np.ndarray] = None
) -> Dict[str, Any]:
    """The features and labels of a recording (rows in run order; the fields of ``J5r``). Returns
    ``X`` (n, 16), ``accepted``, ``evals``, ``revisit`` (arrays)."""
    st = FeatureState(cubic_rotations() if sym is None else sym)
    X = np.zeros((len(rows), N_FEATURES))
    for j, r in enumerate(rows):
        v = int(r["voxel"])
        X[j] = st.features(
            v,
            str(r["kind"]),
            int(r["n_revisits"]),
            float(r["best"]),
            r["start_R"],
            float(r["start_hit"]),
            float(r["start_cost"]),
            float(r["start_conf"]),
            int(r["parent"]),
            int(r["pos"]),
        )
        st.update(
            v,
            r["start_R"],
            bool(r["accepted"]),
            bool(r["kept"]),
            float(r["end_hit"]),
            float(r["end_cost"]),
        )
    return dict(
        X=X,
        accepted=np.array([bool(r["accepted"]) for r in rows]),
        evals=np.array([int(r["evals"]) for r in rows]),
        revisit=np.array([str(r["kind"]) == "revisit" for r in rows]),
    )
