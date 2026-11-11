"""Boundary pass (AR3 step 10, opt-in key ``BoundaryPass``): one fix-up pass over a finished map.

A *boundary voxel* has a neighbour (radius 2 side lengths) whose orientation differs from its own
by more than 0.5 degree (symmetry-reduced). Its *candidates* are those neighbours' orientations,
with candidates within 0.5 degree of one another counted once (the first in neighbour-index order
is kept). A *try* is one CMA fit (first step 0.05 degree, at most 60 evaluations, the start's own
evaluation included, plus the final overlap evaluation) started from one candidate. The rule
``hit`` keeps, among the old end and the ends of the tries, the highest hit ratio, then the lowest
cost, then the earliest (the old end first).

The pass is a port of the step 5 study code (``scripts/autoresearch/ar3/boundary_common.py`` and
``boundary_lib.py``, rule ``hit``). It reads no truth. The candidates of every boundary voxel are
read from the map as it was before the pass (a snapshot), so the order of the voxels does not
matter for them. Its settings are constants of this module (the measured ones), not keys, and it
draws from its own ``default_rng(0)`` (one generator for the whole pass, voxels in index order,
candidates in neighbour-index order), not from the run's generator.

Voxels that were never fitted (``reconstruction_id`` NOT_VISITED, e.g. with ``max_voxels``) are
skipped as targets and as sources of candidates.
"""

import time
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from .mic_file import ReconstructionState
from .ray_indexing import misorientation_deg, proper_rotations
from .reconstructor import _get_voxel_vertices

BOUNDARY_THRESHOLD_DEG = 0.5  # boundary test and candidate merge distance
PASS_SIGMA0_DEG = 0.05  # first step of a try (degrees)
PASS_MAX_EVALS = 60  # evaluation budget of a try, the start's own evaluation included
PASS_SEED = 0  # the pass's own generator


@dataclass
class BoundaryPassResult:
    """What one pass did. ``tries`` has one dict per try (voxel, cand, end, cost, hit, conf,
    evals) in the order made; ``old_score`` maps each boundary voxel to the (cost, hit, conf) of
    its old end; ``changed`` lists the voxels that took a try's end."""

    boundary_voxels: List[int] = field(default_factory=list)
    tries: List[Dict[str, Any]] = field(default_factory=list)
    old_score: Dict[int, Tuple[float, float, float]] = field(default_factory=dict)
    changed: List[int] = field(default_factory=list)
    n_evals_tries: int = 0
    n_evals_score: int = 0
    wall_s: float = 0.0


def boundary_candidates(
    R: np.ndarray,
    neighbours: Sequence[Sequence[int]],
    sym: np.ndarray,
    active: Optional[np.ndarray] = None,
    threshold_deg: float = BOUNDARY_THRESHOLD_DEG,
) -> Tuple[np.ndarray, Dict[int, List[int]]]:
    """(is_boundary (n,), {boundary voxel: candidate voxels}) of a map ``R`` (n, 3, 3).

    ``neighbours[v]`` are the voxels near ``v`` (``v`` excluded); they are taken in index order.
    A candidate is a neighbour whose orientation is more than ``threshold_deg`` from ``v``'s;
    candidates within ``threshold_deg`` of an earlier kept candidate are dropped. ``active``
    (n,) bool: voxels that are False are neither boundary voxels nor candidates. No truth is used.
    """
    n = len(R)
    nbr = [np.sort(np.asarray(x, dtype=np.int64)) for x in neighbours]
    isb = np.zeros(n, dtype=bool)
    cands: Dict[int, List[int]] = {}
    if n == 0 or not any(len(x) for x in nbr):
        return isb, cands
    pi = np.repeat(np.arange(n), [len(x) for x in nbr])
    pj = np.concatenate(nbr)
    far = np.asarray(misorientation_deg(R[pi], R[pj], sym), dtype=float) > threshold_deg
    if active is not None:
        far &= active[pi] & active[pj]
    np.logical_or.at(isb, pi, far)
    off = np.concatenate([[0], np.cumsum([len(x) for x in nbr])])
    for v in np.where(isb)[0]:
        nb = nbr[v][far[off[v] : off[v + 1]]]
        keep: List[int] = []
        for j in nb:
            if keep:
                a = misorientation_deg(np.repeat(R[j][None], len(keep), 0), R[keep], sym)
                if (np.asarray(a) <= threshold_deg).any():
                    continue
            keep.append(int(j))
        cands[int(v)] = keep
    return isb, cands


def keep_hit(old: Tuple[float, float], tries: Sequence[Tuple[float, float]]) -> int:
    """Rule ``hit``: index into ``[old] + tries`` of the end to keep (0: the old end). ``old`` and
    each try are (cost, hit ratio). Highest hit ratio, then lowest cost, then the earliest."""
    cost = np.array([old[0]] + [t[0] for t in tries], dtype=float)
    hit = np.array([old[1]] + [t[1] for t in tries], dtype=float)
    return int(np.lexsort((np.arange(len(cost)), cost, -hit))[0])


def run_boundary_pass(bfs: Any, mic: Any) -> BoundaryPassResult:
    """Run the pass once on ``mic`` (the BFS's map, in place) with rule ``hit``.

    Every voxel that has a fitted orientation is a target and a source. The strict score of a
    boundary voxel's old end is one evaluation (counted in ``n_evals_boundary_score``); the tries'
    evaluations are ``n_evals_boundary_pass``. A voxel that changes takes the chosen end's
    orientation, cost, confidence and hit ratio. Writes ``bfs.stats`` (boundary_voxels,
    boundary_tries, boundary_changed, n_evals_boundary_pass, n_evals_boundary_score,
    wall_boundary_pass_s) and ``bfs.boundary_log`` (one entry per changed voxel with its old and
    new orientation, cost, confidence and hit ratio, so the map before the pass can be rebuilt)."""
    t0 = time.time()
    rec = bfs.reconstructor
    voxels = mic.voxels
    n = len(voxels)
    active = np.array(
        [v.reconstruction_id != ReconstructionState.NOT_VISITED for v in voxels], dtype=bool
    )
    R = np.stack([np.asarray(v.orientation, dtype=np.float64).copy() for v in voxels])
    nbrs = [
        sorted(
            int(j) for j in mic.get_neighbors(v, radius=2.0 * voxels[v].side_length) if int(j) != v
        )
        for v in range(n)
    ]
    sym = proper_rotations(bfs.setup.exp_setup.get_sample_symmetry())
    _isb, cands = boundary_candidates(R, nbrs, sym, active=active)
    cost_fn = rec._make_local_cost_fn()
    rng = np.random.default_rng(PASS_SEED)
    res = BoundaryPassResult(boundary_voxels=sorted(cands))
    log: List[Dict[str, Any]] = []
    for v in res.boundary_voxels:
        vox = voxels[v]
        verts = _get_voxel_vertices(vox)
        info = cost_fn.evaluate(R[v], verts, vox.phase)
        hit0, conf0 = bfs._ratios(info)
        old = (float(info.cost), float(hit0), float(conf0))
        res.old_score[v] = old
        res.n_evals_score += 1
        rows: List[Dict[str, Any]] = []
        for c in cands[v]:
            r = rec.local_optimization(
                voxel_vertices=verts,
                phase_index=vox.phase,
                initial_orientation=R[c],
                rng=rng,
                cma_sigma0_deg=PASS_SIGMA0_DEG,
                cma_max_evals=PASS_MAX_EVALS,
                optimizer="cma",
            )
            hit, conf = bfs._ratios(r.overlap_info)
            rows.append(
                dict(
                    voxel=v,
                    cand=c,
                    end=np.asarray(r.orientation, dtype=np.float64),
                    cost=float(r.cost),
                    hit=float(hit),
                    conf=float(conf),
                    evals=int(rec.last_local_optimization_evals),
                )
            )
        res.tries += rows
        res.n_evals_tries += sum(t["evals"] for t in rows)
        k = keep_hit((old[0], old[1]), [(t["cost"], t["hit"]) for t in rows])
        if k > 0:
            new = rows[k - 1]
            log.append(
                dict(
                    voxel=v,
                    cand=new["cand"],
                    old=R[v].copy(),
                    old_cost=float(vox.cost),
                    old_confidence=float(vox.confidence),
                    old_overlap_ratio=float(vox.overlap_ratio),
                    new=new["end"].copy(),
                    new_cost=new["cost"],
                    new_confidence=new["conf"],
                    new_overlap_ratio=new["hit"],
                )
            )
            vox.orientation = new["end"].copy()
            vox.cost = new["cost"]
            vox.confidence = new["conf"]
            vox.overlap_ratio = new["hit"]
            res.changed.append(v)
    res.wall_s = time.time() - t0
    st = bfs.stats
    st["boundary_voxels"] = len(res.boundary_voxels)
    st["boundary_tries"] = len(res.tries)
    st["boundary_changed"] = len(res.changed)
    st["n_evals_boundary_pass"] = res.n_evals_tries
    st["n_evals_boundary_score"] = res.n_evals_score
    st["wall_boundary_pass_s"] = res.wall_s
    bfs.boundary_log = log
    print(
        f"Boundary pass: {len(res.boundary_voxels)} boundary voxels, {len(res.tries)} tries, "
        f"{len(res.changed)} changed, {res.n_evals_tries} evaluations in {res.wall_s:.1f}s",
        flush=True,
    )
    return res
