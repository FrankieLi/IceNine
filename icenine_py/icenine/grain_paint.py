"""Painting grains (autoresearch AR3, step 1): list the grains by indexing, flood each grain with
single strict scores, give every voxel the passing grain with the lowest strict cost.

Terms (see ``docs/autoresearch/PLAN.md``, AR3):

* A **grain entry** is one accepted seed: its voxel, its orientation and its hit ratio.
* To **score** a grain on a voxel is one strict evaluation of the grain's orientation there (no
  fit); it gives a cost and a hit ratio.
* A grain **passes** on a voxel by the BFS neighbour test, unchanged: hit ratio above 0.9 times the
  best hit ratio seen so far in that grain's flood (the best starts at the seed's hit ratio and is
  raised by every scored voxel, passing or not, as in ``BFSReconstruction._expand``).
* To **flood** a grain is to score it on the neighbours of its seed, then on the neighbours of every
  voxel where it passes (the BFS neighbourhood). A flood does not look at any other grain, so it
  does not depend on them.
* A voxel is **contested** when two or more grains pass on it, **unpainted** when none does.

Two design choices of the loop (neither is in the plan's text; both matter once fits are accepted):

* A voxel settled "by a fit" keeps its fit, even if a grain seeded later passes there with a lower
  strict cost. The painted map takes the fit first (``PaintResult.final_orientation``).
* The fit test compares the fit's hit ratio with the grain's best hit ratio at the end of that
  grain's whole flood (``GrainEntry.best_hit``), not with a running best at the moment the voxel is
  reached, as ``BFSReconstruction._expand`` does. The flood is finished before any fit runs, so the
  whole-flood best is the only one there is.

Evaluations are counted by stage (``STAGES``): the seed path's evaluations go to ``seed_index``
when its result is accepted and to ``index_rejected`` when not; the grid search's go to
``seed_grid`` when accepted and to ``grid_rejected`` when not.

``paint_grains`` is independent of the physics: it works through a ``PaintHooks`` object (a fake
cost in the unit tests, ``ReconstructorHooks`` for a real run). The opt-in config key
``PaintGrains`` (step 2) runs it through ``BFSReconstruction.reconstruct_sample``; with the key off
nothing calls this module.
"""

import time
from collections import deque
from dataclasses import dataclass, field
from typing import Any, Deque, Dict, List, Optional, Protocol, Sequence, Set, Tuple

import numpy as np

PASS_FRACTION = 0.9  # the BFS neighbour test: hit / best > 0.9
STAGES = (
    "seed_index",
    "seed_grid",
    "score",
    "fit_accepted",
    "fit_rejected",
    "index_rejected",
    "grid_rejected",
)


class PaintHooks(Protocol):
    """What the loop needs from the cost function and the reconstructor."""

    min_seed_hit: float  # the seed test (``MinAccelerationThreshold``)

    def neighbours(self, v: int) -> Sequence[int]:
        """Voxels the BFS expansion reaches from ``v`` (radius 2 side lengths), ``v`` excluded."""

    def score(self, R: np.ndarray, v: int) -> Tuple[float, float]:
        """One strict evaluation of orientation ``R`` on voxel ``v``: (cost, hit ratio)."""

    def local_fit(self, v: int, R0: np.ndarray) -> Tuple[np.ndarray, float, float, int]:
        """One local fit on ``v`` from ``R0``: (orientation, cost, hit ratio, evaluations)."""

    def index_seed(self, v: int) -> Tuple[Optional[np.ndarray], float, float, int]:
        """The seed path on ``v`` (index, strict cost on 5, finisher on 2): (orientation or None,
        cost, hit ratio, evaluations). Orientation None: no candidate."""

    def grid_seed(self, v: int) -> Tuple[np.ndarray, float, float, int]:
        """The grid search on ``v`` (the fallback): (orientation, cost, hit ratio, evaluations)."""

    def misorientation_deg(self, Ra: np.ndarray, Rb: np.ndarray) -> float:
        """Symmetry-reduced misorientation angle in degrees."""


@dataclass
class GrainEntry:
    """One accepted seed. ``kind``: "index" (the seed path) or "grid" (the fallback)."""

    voxel: int
    R: np.ndarray
    cost: float
    hit: float
    kind: str
    best_hit: float = 0.0  # best hit ratio seen in the flood
    near_earlier: Optional[int] = None  # index of an earlier entry within ``near_deg``
    near_deg: float = float("nan")  # misorientation to that entry
    flood_scores: int = 0
    flood_passes: int = 0


@dataclass
class VoxelFit:
    """A voxel settled "by a fit" (or a rejected seed's result for a failed voxel)."""

    R: np.ndarray
    cost: float
    hit: float
    from_grain: int = -1


@dataclass
class PaintResult:
    entries: List[GrainEntry] = field(default_factory=list)
    # per grain: voxel -> (cost, hit, passed)
    scores: List[Dict[int, Tuple[float, float, bool]]] = field(default_factory=list)
    stage: Dict[int, str] = field(default_factory=dict)  # "seed", "fit", "failed" (loop outcome)
    fits: Dict[int, VoxelFit] = field(default_factory=dict)  # fit-settled voxels
    failed: Dict[int, VoxelFit] = field(default_factory=dict)  # failed voxels: rejected result
    evals: Dict[str, int] = field(default_factory=dict)
    wall_s: Dict[str, float] = field(default_factory=dict)
    counts: Dict[str, int] = field(default_factory=dict)
    log: List[Dict[str, Any]] = field(default_factory=list)  # one row per loop decision

    def _add(self, key: str, n: int, dt: float) -> None:
        self.evals[key] = self.evals.get(key, 0) + int(n)
        self.wall_s[key] = self.wall_s.get(key, 0.0) + float(dt)

    def passing(self, v: int) -> List[Tuple[int, float, float]]:
        """(grain, cost, hit) of the grains that pass on ``v``, in grain order."""
        return [(g, s[v][0], s[v][1]) for g, s in enumerate(self.scores) if v in s and s[v][2]]

    def ranked(self, v: int) -> List[Tuple[int, float, float]]:
        """The passing grains of ``v`` by strict cost (ties: lowest grain index first)."""
        return sorted(self.passing(v), key=lambda t: (t[1], t[0]))

    def painted_grain(self, v: int) -> Optional[int]:
        """The grain that paints ``v``: the passing grain with the lowest strict cost."""
        r = self.ranked(v)
        return r[0][0] if r else None

    def final_orientation(self, v: int) -> Tuple[Optional[np.ndarray], str]:
        """(orientation, source) of voxel ``v`` in the painted map. A voxel settled by a fit keeps
        its fit; otherwise the painting grain; a failed voxel nobody painted keeps the rejected
        search's result. Source: "fit", "grain", "failed" or "none"."""
        if v in self.fits:
            return self.fits[v].R, "fit"
        g = self.painted_grain(v)
        if g is not None:
            return self.entries[g].R, "grain"
        if v in self.failed:
            return self.failed[v].R, "failed"
        return None, "none"


def paint_grains(
    order: Sequence[int], n_voxels: int, hooks: PaintHooks, near_deg: float = 1.0
) -> PaintResult:
    """The loop of AR3 step 1 (see the module docstring and the plan).

    ``order`` is the shuffled voxel order of the BFS runs (no new random draw here). The next voxel
    that is unpainted, not settled and not failed becomes a seed. If at least one grain was scored
    on it, it first gets one local fit from the lowest-cost of those grains; if that fit passes the
    neighbour test (against that grain's best hit ratio over its whole flood) the voxel is settled
    "by a fit", keeps that fit whatever grain later passes there, and the loop goes on. Otherwise
    the seed path runs; if its result passes the seed test it becomes a grain entry and is flooded;
    if not, the grid search runs and the same test applies; a voxel whose search is rejected is
    failed (it may still be painted by a later grain).
    """
    res = PaintResult()
    res.counts.update({k: 0 for k in ("seeds", "fits_accepted", "fits_rejected", "failed")})
    has_pass: Set[int] = set()
    settled: Set[int] = set()  # seeds, fit-settled and failed voxels
    scored_on: Dict[int, List[int]] = {}  # voxel -> grains scored on it (any outcome)

    def flood(g: int) -> None:
        e = res.entries[g]
        sc = res.scores[g]
        best = e.hit
        queue: Deque[int] = deque(hooks.neighbours(e.voxel))
        seen = {e.voxel} | set(queue)
        while queue:
            u = queue.popleft()
            t0 = time.perf_counter()
            cost, hit = hooks.score(e.R, u)
            res._add("score", 1, time.perf_counter() - t0)
            best = max(best, hit)
            ok = best > 0 and hit / best > PASS_FRACTION
            sc[u] = (float(cost), float(hit), bool(ok))
            scored_on.setdefault(u, []).append(g)
            e.flood_scores += 1
            if ok:
                e.flood_passes += 1
                has_pass.add(u)
                for w in hooks.neighbours(u):
                    if w not in seen:
                        seen.add(w)
                        queue.append(w)
        e.best_hit = best

    def add_entry(v: int, R: np.ndarray, cost: float, hit: float, kind: str) -> None:
        g = len(res.entries)
        e = GrainEntry(voxel=v, R=np.asarray(R, dtype=np.float64), cost=cost, hit=hit, kind=kind)
        for j, o in enumerate(res.entries):
            d = hooks.misorientation_deg(o.R, e.R)
            if d < near_deg and (e.near_earlier is None or d < e.near_deg):
                e.near_earlier, e.near_deg = j, float(d)
        res.entries.append(e)
        res.scores.append({v: (float(cost), float(hit), True)})
        scored_on.setdefault(v, []).append(g)
        has_pass.add(v)
        settled.add(v)
        res.stage[v] = "seed"
        res.counts["seeds"] += 1
        flood(g)

    for v in order:
        if v in settled or v in has_pass:
            continue
        row: Dict[str, Any] = dict(voxel=int(v))
        grains = scored_on.get(v, [])
        if grains:
            g = min(grains, key=lambda j: (res.scores[j][v][0], j))
            t0 = time.perf_counter()
            R, cost, hit, n = hooks.local_fit(v, res.entries[g].R)
            dt = time.perf_counter() - t0
            best = max(res.entries[g].best_hit, hit)
            ok = best > 0 and hit / best > PASS_FRACTION
            res._add("fit_accepted" if ok else "fit_rejected", n, dt)
            row.update(fit_from=g, fit_hit=float(hit), fit_cost=float(cost), fit_ok=bool(ok))
            if ok:
                res.fits[v] = VoxelFit(np.asarray(R, dtype=np.float64), cost, hit, g)
                settled.add(v)
                res.stage[v] = "fit"
                res.counts["fits_accepted"] += 1
                res.log.append(row)
                continue
            res.counts["fits_rejected"] += 1
        t0 = time.perf_counter()
        R_i, cost_i, hit_i, n_i = hooks.index_seed(v)
        dt = time.perf_counter() - t0
        row.update(index_evals=int(n_i), index_hit=float(hit_i) if R_i is not None else None)
        if R_i is not None and hit_i >= hooks.min_seed_hit:
            res._add("seed_index", n_i, dt)
            row["seed"] = "index"
            res.log.append(row)
            add_entry(v, R_i, float(cost_i), float(hit_i), "index")
            continue
        res._add("index_rejected", n_i, dt)  # the seed path's evaluations that were not used
        t0 = time.perf_counter()
        R_g, cost_g, hit_g, n_g = hooks.grid_seed(v)
        dt = time.perf_counter() - t0
        row.update(grid_evals=int(n_g), grid_hit=float(hit_g))
        if hit_g >= hooks.min_seed_hit:
            res._add("seed_grid", n_g, dt)
            row["seed"] = "grid"
            res.log.append(row)
            add_entry(v, R_g, float(cost_g), float(hit_g), "grid")
            continue
        res._add("grid_rejected", n_g, dt)  # a rejected grid search: not a seed path
        row["seed"] = "failed"
        res.log.append(row)
        res.failed[v] = VoxelFit(np.asarray(R_g, dtype=np.float64), cost_g, hit_g)
        settled.add(v)
        res.stage[v] = "failed"
        res.counts["failed"] += 1
    for k in STAGES:
        res.evals.setdefault(k, 0)
        res.wall_s.setdefault(k, 0.0)
    return res


# --------------------------------------------------------------------------------------------
# the adapter for the real reconstructor
# --------------------------------------------------------------------------------------------
class ReconstructorHooks:
    """``PaintHooks`` on a ``BFSReconstruction`` (its setup, strict cost, local fit and seed paths).

    One strict cost function is built once and used for every score, so one score is one evaluation
    (``cost_fn.eval_count`` rises by one). The seed path and the grid search are the reconstructor's
    own (``AdaptiveVoxelReconstructor.index_seed`` and ``reconstruct_voxel``); their acceptance
    hit ratio comes from ``evaluate_overlap``, which the BFS does not count either."""

    def __init__(self, bfs: Any, mic: Any, rng: np.random.Generator) -> None:
        from .ray_indexing import proper_rotations
        from .reconstructor import _get_voxel_vertices

        self.bfs = bfs
        self.rec = bfs.reconstructor
        self.mic = mic
        self.rng = rng
        self.min_seed_hit = float(bfs.setup.config.min_acceleration_threshold)
        self.cost_fn = self.rec._make_local_cost_fn()
        self.sym = proper_rotations(bfs.setup.exp_setup.get_sample_symmetry())
        self._vertices = _get_voxel_vertices
        self._nb: Dict[int, List[int]] = {}
        self._vcache: Dict[int, Any] = {}
        # study record (AR3 step 2): one row per seed path run, from ``SeedResult``
        self.seed_records: List[Dict[str, Any]] = []

    def vertices(self, v: int) -> Any:
        if v not in self._vcache:
            self._vcache[v] = self._vertices(self.mic.voxels[v])
        return self._vcache[v]

    def neighbours(self, v: int) -> List[int]:
        if v not in self._nb:
            r = 2.0 * self.mic.voxels[v].side_length
            self._nb[v] = [int(n) for n in self.mic.get_neighbors(v, radius=r)]
        return self._nb[v]

    def score(self, R: np.ndarray, v: int) -> Tuple[float, float]:
        info = self.cost_fn.evaluate(
            np.asarray(R, np.float64), self.vertices(v), self.mic.voxels[v].phase
        )
        hit, _ = self.bfs._ratios(info)
        return float(info.cost), float(hit)

    def local_fit(self, v: int, R0: np.ndarray) -> Tuple[np.ndarray, float, float, int]:
        vox = self.mic.voxels[v]
        r = self.rec.local_optimization(
            voxel_vertices=self.vertices(v),
            phase_index=vox.phase,
            initial_orientation=np.asarray(R0, np.float64),
            rng=self.rng,
        )
        hit, _ = self.bfs._ratios(r.overlap_info)
        return (
            np.asarray(r.orientation, np.float64),
            float(r.cost),
            float(hit),
            int(self.rec.last_local_optimization_evals),
        )

    def index_seed(self, v: int) -> Tuple[Optional[np.ndarray], float, float, int]:
        vox = self.mic.voxels[v]
        s = self.rec.index_seed(self.vertices(v), vox.phase)
        self.seed_records.append(
            dict(
                voxel=int(v),
                n_cand=int(s.n_cand),
                evals=int(s.evals),
                taken=int(s.taken),
                strict_costs=[float(c) for c in s.strict_costs],
                refined=[int(j) for j in s.refined],
                finisher_costs=[float(c) for c in s.finisher_costs],
                finisher_evals=[int(e) for e in s.finisher_evals],
                orientation=None if s.orientation is None else np.asarray(s.orientation).tolist(),
                cost=float(s.cost) if s.orientation is not None else None,
            )
        )
        if s.orientation is None:
            return None, 1.0, 0.0, int(s.evals)
        info = self.rec.evaluate_overlap(s.orientation, self.vertices(v), vox.phase)
        hit, _ = self.bfs._ratios(info)
        return np.asarray(s.orientation, np.float64), float(s.cost), float(hit), int(s.evals)

    def grid_seed(self, v: int) -> Tuple[np.ndarray, float, float, int]:
        vox = self.mic.voxels[v]
        r = self.rec.reconstruct_voxel(
            voxel_vertices=self.vertices(v), phase_index=vox.phase, rng=self.rng
        )
        n = int(self.rec.last_eval_counts[2])
        info = self.rec.evaluate_overlap(r.orientation, self.vertices(v), vox.phase)
        hit, _ = self.bfs._ratios(info)
        return np.asarray(r.orientation, np.float64), float(r.cost), float(hit), n

    def misorientation_deg(self, Ra: np.ndarray, Rb: np.ndarray) -> float:
        from .ray_indexing import misorientation_deg

        return float(misorientation_deg(Ra[None], Rb[None], self.sym)[0])
