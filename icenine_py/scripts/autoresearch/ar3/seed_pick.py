"""AR3 step 6, part 2: the picked seed remedy as a script wrapper (no library change).

The pick of part 1 (see ``seed/summary_part01.json``): ``neighbours_only then more20_5``. At a seed
or reseed voxel, in this order:

  1. one try (a CMA fit, first step 0.05 degrees, at most 60 evaluations plus the final one) from
     each distinct orientation (0.5 degree merge, symmetry-reduced, first in neighbour order)
     among the voxel's neighbours that are FITTED at that moment; keep rule ``hit`` (highest hit
     ratio, then lowest cost, then neighbour order). If the kept end has hit ratio >= 0.80 it is
     the seed result and nothing else is spent;
  2. otherwise the seed step with the strict cost on the first 20 candidates (in place of 5) and
     the local optimization on the best 5 (in place of 2); the lowest final cost is the result.

The caller (``BFSReconstruction._index_seed``) applies the seed test as before; a failed result
goes to the grid search as before. Two instance methods are replaced and nothing else:
``bfs._index_seed`` (only to tell the wrapper which voxel it is for) and
``bfs.reconstructor.index_seed`` (which returns a ``SeedResult`` with the evaluations of both
steps). Unlike part 1, the neighbours are the ones FITTED at that moment of the run, not the final
map.
"""

import time
from typing import Any, Dict, List, Optional

import numpy as np

TEST = 0.80
SIGMA0_DEG = 0.05
MAX_EVALS = 60
MERGE_DEG = 0.5


class SeedPick:
    def __init__(self, bfs: Any, mic: Any, orig_local_optimization: Any, sym: np.ndarray) -> None:
        self.bfs = bfs
        self.mic = mic
        self.rec = bfs.reconstructor
        self.orig_lo = orig_local_optimization  # the unwrapped method, called with cma keywords
        self.sym = sym
        self.cur: Optional[int] = None
        self.nb: Dict[int, List[int]] = {}
        self.log: List[Dict[str, Any]] = []
        self.orig_index_seed = self.rec.index_seed
        self.orig_bfs_index_seed = bfs._index_seed

    def install(self) -> None:
        self.bfs._index_seed = self._bfs_index_seed
        self.rec.index_seed = self._index_seed

    # ---- helpers
    def _neighbours(self, v: int) -> List[int]:
        if v not in self.nb:
            r = 2.0 * self.mic.voxels[v].side_length
            self.nb[v] = [int(n) for n in self.mic.get_neighbors(v, radius=r) if int(n) != v]
        return self.nb[v]

    def _distinct(self, starts: List[int]) -> List[int]:
        from icenine.ray_indexing import misorientation_deg

        keep: List[int] = []
        for j in starts:
            Rj = np.asarray(self.mic.voxels[j].orientation, np.float64)
            if keep:
                Rk = np.stack(
                    [np.asarray(self.mic.voxels[k].orientation, np.float64) for k in keep]
                )
                if (
                    misorientation_deg(np.repeat(Rj[None], len(keep), 0), Rk, self.sym) <= MERGE_DEG
                ).any():
                    continue
            keep.append(j)
        return keep

    def _bfs_index_seed(self, voxel: Any, vertices: Any, rec: Any, reseed: bool) -> Any:
        self.cur = int(rec.voxel_idx)
        self.reseed = bool(reseed)
        try:
            return self.orig_bfs_index_seed(voxel, vertices, rec, reseed)
        finally:
            self.cur = None

    def _try(self, vertices: Any, phase: int, R0: np.ndarray) -> Dict[str, Any]:
        p = self.rec.params
        old = p.local_optimizer
        p.local_optimizer = "cma"
        try:
            r = self.orig_lo(
                voxel_vertices=vertices,
                phase_index=phase,
                initial_orientation=np.asarray(R0, np.float64),
                rng=np.random.default_rng(0),
                cma_sigma0_deg=SIGMA0_DEG,
                cma_max_evals=MAX_EVALS,
            )
        finally:
            p.local_optimizer = old
        hit, _ = self.bfs._ratios(r.overlap_info)
        return dict(
            R=np.asarray(r.orientation, np.float64),
            cost=float(r.cost),
            hit=float(hit),
            evals=int(self.rec.last_local_optimization_evals),
        )

    # ---- the replacement of ``AdaptiveVoxelReconstructor.index_seed``
    def _index_seed(
        self, vertices: Any, phase_index: int, n_cost: int = 5, n_refine: int = 2
    ) -> Any:
        from icenine import ray_indexing as RI
        from icenine.orientation_search import SearchCandidate
        from icenine.reconstructor import ReconstructionState

        t0 = time.perf_counter()
        v = self.cur
        assert v is not None
        entry: Dict[str, Any] = dict(voxel=v, reseed=self.reseed)
        starts = [
            j
            for j in self._neighbours(v)
            if self.mic.voxels[j].reconstruction_id == ReconstructionState.FITTED
        ]
        starts = self._distinct(starts)
        entry["n_starts"] = len(starts)
        tries: List[Dict[str, Any]] = []
        for j in starts:
            tries.append(self._try(vertices, phase_index, self.mic.voxels[j].orientation))
        ev_tries = int(sum(t["evals"] for t in tries))
        entry["tries_evals"] = ev_tries
        if tries:
            k = int(
                np.lexsort(
                    (np.arange(len(tries)), [t["cost"] for t in tries], [-t["hit"] for t in tries])
                )[0]
            )
            best = tries[k]
            entry["tries_hit"] = best["hit"]
            if best["hit"] >= TEST:
                entry["by"] = "tries"
                entry["evals"] = ev_tries
                self.log.append(entry)
                dt = time.perf_counter() - t0
                return RI.SeedResult(
                    orientation=best["R"],
                    cost=best["cost"],
                    taken=-1,
                    n_cand=0,
                    n_att=0,
                    evals=ev_tries,
                    index_s=0.0,
                    total_s=dt,
                )
        # the seed step with the strict cost on 20 and the local optimization on the best 5
        res = self.orig_index_seed(vertices, phase_index, n_cost=20, n_refine=5)
        res.evals = int(res.evals) + ev_tries
        res.total_s = time.perf_counter() - t0
        entry["by"] = "seed_step"
        entry["evals"] = int(res.evals)
        entry["found"] = res.orientation is not None
        self.log.append(entry)
        return res
