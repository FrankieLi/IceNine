"""Boundary-pass repair (AR5 step 20, opt-in key ``BoundaryPassRepair``): a second look at the
contested voxels after the boundary pass, with a cheap local optimization.

After the pass, a voxel is **contested** when its hit ratio (the share of simulated peaks that land
on lit pixels) is under 0.95 of the median hit ratio of its neighbours (radius 2 side lengths) AND
its **margin** is under 0.15. The margin is its hit ratio minus the best hit ratio among its *other*
pass ends: the old end and the try ends more than 1 degree from the kept end. A voxel with no other
pass end is never contested. No truth is used.

Each contested voxel, in index order, gets one **CMA local optimization** (the library's local CMA
run, first step 0.05 degree, at most 120 evaluations with the start's own, plus the final
evaluation: 121 in all) from its own end and one from each **candidate** (a neighbour's orientation
more than 0.5 degree from its own, merged within 0.5 degree, read from a snapshot of the map after
the pass). Each one has its own generator ``default_rng([0, v, c + 1])`` (``c`` = -1 for the voxel's
own end), so no draw depends on any other fit. The kept end is the one with the highest hit ratio,
then the lowest cost, then the current end (``boundary_pass.keep_hit``). A voxel whose stored hit
ratio is missing costs one counted evaluation. Every evaluation is counted.

This is the library form of the script rule ``scripts/autoresearch/ar5/cmarepair_rules.py`` (study
step 20, arm ``J5C``), which imports step 19's frozen flag; the functions below are ports of that
flag, and ``tests/test_boundary_repair.py`` checks them against it. The settings are constants of
this module, not keys. The state after the pass is the voxel's own ``cost`` and ``overlap_ratio``,
except for boundary voxels, whose old-end score is the pass's strict re-score (as the script's
``old_score``).
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from .boundary_pass import boundary_candidates, keep_hit
from .ray_indexing import misorientation_deg

RATIO_MAX = 0.95  # contested: hit ratio under this share of the neighbours' median
MARGIN_MAX = 0.15  # ... and margin under this
FAR_DEG = 1.0  # a pass end is "other" when more than this from the kept end
SIGMA0_DEG = 0.05  # first step of a repair local optimization (degrees)
MAX_EVALS = 120  # its budget, the start's own evaluation included (121 with the final one)
SEED = 0  # first entry of every per-fit generator
NO_OTHER = 1.0  # margin of a voxel with no other pass end (never contested)


@dataclass
class RepairInfo:
    """What one repair did. ``fits`` has one dict per local optimization (voxel, start voxel or
    None for the own end, end, cost, hit, conf, evals) in the order made; ``changed`` lists the
    voxels that took a fit's end. ``R_pre``, ``cost_pre``, ``hit_pre``: the state after the pass,
    before the repair (copies)."""

    flagged: List[int] = field(default_factory=list)
    fits: List[Dict[str, Any]] = field(default_factory=list)
    changed: List[int] = field(default_factory=list)
    n_evals_fits: int = 0
    n_evals_score: int = 0
    wall_s: float = 0.0
    R_pre: Optional[np.ndarray] = None
    cost_pre: Optional[np.ndarray] = None
    hit_pre: Optional[np.ndarray] = None

    @property
    def n_evals(self) -> int:
        return self.n_evals_fits + self.n_evals_score


def fit_rng(v: int, c: Optional[int]) -> np.random.Generator:
    """The generator of the fit at voxel ``v`` from start ``c`` (None: the voxel's own end)."""
    cc = -1 if c is None else int(c)
    return np.random.default_rng([SEED, int(v), cc + 1])


def other_end_hit(
    R_after: np.ndarray,
    R_before: np.ndarray,
    old_hit: np.ndarray,
    t_voxel: np.ndarray,
    t_end: np.ndarray,
    t_hit: np.ndarray,
    sym: np.ndarray,
) -> np.ndarray:
    """(n,) best hit among each voxel's other pass ends: the try ends more than ``FAR_DEG`` from the
    kept end, and the old end when it is more than ``FAR_DEG`` from the kept end; -1 for none."""
    n = len(R_after)
    other = np.full(n, -1.0)
    if len(t_voxel):
        tv = np.asarray(t_voxel)
        far = np.asarray(misorientation_deg(t_end, R_after[tv], sym), dtype=float) > FAR_DEG
        np.maximum.at(other, tv[far], np.asarray(t_hit, dtype=float)[far])
    oldfar = np.asarray(misorientation_deg(R_before, R_after, sym), dtype=float) > FAR_DEG
    return np.where(oldfar, np.maximum(other, np.asarray(old_hit, dtype=float)), other)


def margin_of(hit: np.ndarray, other: np.ndarray) -> np.ndarray:
    """Margin: hit minus the best other end's hit; ``NO_OTHER`` where there is none."""
    return np.where(other < 0, NO_OTHER, hit - other)


def ratio_of(hit: np.ndarray, nbrs: Sequence[Sequence[int]]) -> np.ndarray:
    """Hit over the median hit of the neighbours (inf when a voxel has no neighbour)."""
    med = np.array(
        [np.median(hit[np.asarray(x, dtype=np.int64)]) if len(x) else np.nan for x in nbrs]
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        r = hit / med
    return np.where(np.isnan(med), np.inf, r)


def contested(
    hit: np.ndarray, other: np.ndarray, nbrs: Sequence[Sequence[int]]
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(flag (n,) bool, ratio, margin): ratio under ``RATIO_MAX``, margin under ``MARGIN_MAX``."""
    ratio = ratio_of(hit, nbrs)
    margin = margin_of(hit, other)
    return (ratio < RATIO_MAX) & (margin < MARGIN_MAX), ratio, margin


def run_repair(
    bfs: Any,
    mic: Any,
    res: Any,
    R_before: np.ndarray,
    nbrs: Sequence[Sequence[int]],
    sym: np.ndarray,
    active: Optional[np.ndarray],
    cost_fn: Any,
) -> RepairInfo:
    """Repair the contested voxels of ``mic`` (the BFS's map after the pass, in place).

    ``res``: the pass's ``BoundaryPassResult`` (its tries and old-end scores); ``R_before``: the
    map before the pass (n, 3, 3); ``nbrs``: the pass's neighbour lists; ``cost_fn``: the strict
    cost function the pass scored with. Returns the ``RepairInfo``; also writes ``bfs.stats``
    (repair_flagged, repair_fits, repair_changed, n_evals_boundary_repair,
    wall_boundary_repair_s) and ``bfs.repair_log`` (one entry per changed voxel, with its state
    just before the repair, so the map before the repair can be rebuilt)."""
    import time

    from .reconstructor import _get_voxel_vertices

    t0 = time.time()
    rec = bfs.reconstructor
    voxels = mic.voxels
    n = len(voxels)
    R = np.stack([np.asarray(v.orientation, dtype=np.float64).copy() for v in voxels])
    cost = np.array([float(v.cost) for v in voxels], dtype=float)
    hit = np.array([float(v.overlap_ratio) for v in voxels], dtype=float)
    changed_pass = set(res.changed)
    for v, (c0, h0, _f0) in res.old_score.items():  # the pass's re-score of the old end
        if v not in changed_pass:
            cost[v], hit[v] = c0, h0
    info = RepairInfo(R_pre=R.copy(), cost_pre=cost.copy(), hit_pre=hit.copy())
    t_voxel = np.array([t["voxel"] for t in res.tries], dtype=np.int64)
    t_end = np.array([t["end"] for t in res.tries], dtype=np.float64).reshape(-1, 3, 3)
    t_hit = np.array([t["hit"] for t in res.tries], dtype=float)
    old_hit = np.array(
        [res.old_score[v][1] if v in res.old_score else hit[v] for v in range(n)], dtype=float
    )
    other = other_end_hit(R, R_before, old_hit, t_voxel, t_end, t_hit, sym)
    flag, _ratio, _margin = contested(hit, other, nbrs)
    flagged = [int(v) for v in np.where(flag)[0]]
    _isb, cands = boundary_candidates(R, nbrs, sym, active=active)
    info.flagged = flagged
    log: List[Dict[str, Any]] = []
    for v in flagged:
        vox = voxels[v]
        verts = _get_voxel_vertices(vox)
        if np.isnan(hit[v]) or np.isnan(cost[v]):
            e = cost_fn.evaluate(R[v], verts, vox.phase)
            h, _c = bfs._ratios(e)
            cost[v], hit[v] = float(e.cost), float(h)
            info.n_evals_score += 1
        starts: List[Tuple[Optional[int], np.ndarray]] = [(None, R[v])]
        starts += [(int(c), R[int(c)]) for c in cands.get(v, [])]
        rows: List[Dict[str, Any]] = []
        for c, R0 in starts:
            r = rec.local_optimization(
                voxel_vertices=verts,
                phase_index=vox.phase,
                initial_orientation=np.asarray(R0, dtype=np.float64),
                rng=fit_rng(v, c),
                cma_sigma0_deg=SIGMA0_DEG,
                cma_max_evals=MAX_EVALS,
                optimizer="cma",
            )
            h, conf = bfs._ratios(r.overlap_info)
            rows.append(
                dict(
                    voxel=v,
                    start=c,
                    end=np.asarray(r.orientation, dtype=np.float64),
                    cost=float(r.cost),
                    hit=float(h),
                    conf=float(conf),
                    evals=int(rec.last_local_optimization_evals),
                )
            )
        info.fits += rows
        info.n_evals_fits += sum(x["evals"] for x in rows)
        k = keep_hit((float(cost[v]), float(hit[v])), [(x["cost"], x["hit"]) for x in rows])
        if k > 0:
            new = rows[k - 1]
            log.append(
                dict(
                    voxel=v,
                    start=new["start"],
                    old=np.asarray(vox.orientation, dtype=np.float64).copy(),
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
            info.changed.append(v)
    info.wall_s = time.time() - t0
    st = bfs.stats
    st["repair_flagged"] = len(info.flagged)
    st["repair_fits"] = len(info.fits)
    st["repair_changed"] = len(info.changed)
    st["n_evals_boundary_repair"] = info.n_evals
    st["wall_boundary_repair_s"] = info.wall_s
    bfs.repair_log = log
    print(
        f"Boundary repair: {len(info.flagged)} contested voxels, {len(info.fits)} local "
        f"optimizations, {len(info.changed)} changed, {info.n_evals} evaluations "
        f"in {info.wall_s:.1f}s",
        flush=True,
    )
    return info
