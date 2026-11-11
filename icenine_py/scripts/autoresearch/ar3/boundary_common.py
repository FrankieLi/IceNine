"""AR3 step 5 (task ar3e): shared set-up of the boundary study.

The six stored maps, the two environments they were made in (the gradient sample with its images;
step 2's layer with its images), the reconstructor's own neighbour relation, the boundary voxels
and candidates of a map, the strict score of one orientation, and one try (a CMA fit started from
a candidate). Truth is exposed for scoring and for part 0 only; nothing here feeds a rule.

Terms (plan, step 5): a boundary voxel has a neighbour whose orientation differs from its own by
more than 0.5 degrees (symmetry-reduced); its candidates are those neighbours' orientations, with
candidates within 0.5 degrees of each other counted once (the first in neighbour-index order is
kept); a try is one CMA fit (first step 0.05 degrees, at most 60 evaluations plus the final one).
"""

from typing import Any, Dict, List, Optional, Tuple

import ar3_common as C
import gradient_common as G
import numpy as np
import run_fine as RF
import run_layer as RL

BND = C.OUT / "boundary"
CACHE = BND / "cache"
RUNS = BND / "runs"
FINE_RUNS = C.OUT / "fine" / "runs"
L2 = C.OUT / "layer" / "runs"
THR = 0.5  # degrees: boundary threshold and candidate merge distance (set, not tuned)
WRONG = G.WRONG

# name -> (environment, what it is, name of the compared map)
MAPS: Dict[str, Tuple[str, str, Optional[str]]] = {
    "D4": ("grad", "step 4 D4 (BFS, indexed seeds, picked fit)", "D"),
    "P4R4": ("grad", "step 4 P4 plus R4 (paint, picked fit)", "PR"),
    "E4": ("layer", "step 4 E4 (BFS, step 2's sample)", "D2"),
    "D": ("grad", "step 3 D (BFS, stock fit)", None),
    "PR": ("grad", "step 3 P plus R (paint, stock refine)", None),
    "D2": ("layer", "step 2 D (BFS, stock fit)", None),
}
STEP4 = ("D4", "P4R4", "E4")
# part 2 maps (made by run_keys.py); compared with the step 4 map of the same sample
MAPS_KEYS: Dict[str, Tuple[str, str, Optional[str]]] = {
    "K": ("grad", "run K (BFS, indexed seeds, CMA keys, gradient sample)", "D4"),
    "KE": ("layer", "run KE (BFS, indexed seeds, CMA keys, step 2's sample)", "E4"),
}
ALL_MAPS = {**MAPS, **MAPS_KEYS}


def orthonormal(R: np.ndarray) -> np.ndarray:
    """Nearest rotation of each matrix (the step 3 map is stored as float32)."""
    U, _s, Vt = np.linalg.svd(R)
    return U @ Vt


def load_map(name: str, n: int) -> np.ndarray:
    """The stored orientations (n, 3, 3), float64, of one of the six maps."""
    if name == "D4":
        return G.read_map(FINE_RUNS / "D4" / "recon.mic.gz", n)[0]
    if name == "P4R4":
        base = G.read_map(FINE_RUNS / "P4" / "recon.mic.gz", n)[0].copy()
        for r in G.load_gz(FINE_RUNS / "P4" / "refine_R4.json.gz"):
            base[r["voxel"]] = np.array(r["end"])
        return base
    if name == "E4":
        return G.read_map(FINE_RUNS / "E4" / "recon.mic.gz", n)[0]
    if name == "D":
        return G.read_map(G.GRAD / "runs" / "D" / "recon.mic.gz", n)[0]
    if name == "PR":
        R = np.load(G.GRAD / "runs" / "P" / "maps.npz")["R"].astype(np.float64)
        return orthonormal(R)
    if name == "D2":
        return G.read_map(L2 / "D" / "recon.mic", n)[0]
    if name in MAPS_KEYS:
        return G.read_map(RUNS / name / "recon.mic.gz", n)[0]
    raise KeyError(name)


class Env:
    """One reconstruction environment: reconstructor, strict cost, neighbour relation, truth."""

    def __init__(self, kind: str) -> None:
        assert kind in ("grad", "layer")
        self.kind = kind
        out = CACHE / f"env_{kind}"
        out.mkdir(parents=True, exist_ok=True)
        if kind == "grad":
            setup, bfs, _t = C.build(C.CONFIG, out, sample=G.GRID_MIC, stack_name=G.STACK)
        else:
            setup, bfs, _t = C.build(C.CONFIG, out, sample=RL.LAYER_SAMPLE)
        from icenine import grain_paint as GP

        self.setup, self.bfs = setup, bfs
        self.mic = setup.sample.get_mic()
        self.n = len(self.mic.voxels)
        self.rng = np.random.default_rng(0)
        self.hooks = GP.ReconstructorHooks(bfs, self.mic, self.rng)
        self.sym = self.hooks.sym
        if kind == "grad":
            T = G.Truth()
            assert T.n == self.n
            self.T: Any = T
            self.Rtrue = T.R
            self.gmap = np.asarray(T.gmap)
            self.level: Optional[np.ndarray] = T.level
            self._err = T.err
        else:
            import layer_summary as LS

            W = LS.World()
            assert W.n == self.n
            self.T = W
            self.Rtrue = W.Rtrue
            self.gmap = np.asarray(W.L.gmap)
            self.level = None
            self._err = W.err
        self.isA = (self.gmap % 2) == 0  # even true-grain id: pick set A; odd: test set B
        self.nbr: List[np.ndarray] = [
            np.asarray(sorted(self.hooks.neighbours(v)), dtype=np.int64) for v in range(self.n)
        ]
        self.pi = np.repeat(np.arange(self.n), [len(x) for x in self.nbr])
        self.pj = np.concatenate(self.nbr)

    # ---- orientation helpers
    def mis(self, Ra: np.ndarray, Rb: np.ndarray) -> np.ndarray:
        from icenine.ray_indexing import misorientation_deg

        return np.asarray(misorientation_deg(Ra, Rb, self.sym), dtype=float)

    def err(self, R: np.ndarray) -> np.ndarray:
        """Error against the truth (analysis and scoring only)."""
        return np.asarray(self._err(R), dtype=float)

    # ---- the boundary and the candidates (the map only; no truth)
    def boundary(self, R: np.ndarray) -> Tuple[np.ndarray, Dict[int, List[int]]]:
        """(is_boundary (n,), {voxel: [neighbour voxels whose orientation is a candidate]})."""
        ang = self.mis(R[self.pi], R[self.pj])
        far = ang > THR
        isb = np.zeros(self.n, dtype=bool)
        np.logical_or.at(isb, self.pi, far)
        cands: Dict[int, List[int]] = {}
        for v in np.where(isb)[0]:
            nb = [int(j) for j in self.nbr[v] if far[self._slot(v, j)]]
            cands[int(v)] = self._merge(R, nb)
        return isb, cands

    def _slot(self, v: int, j: int) -> int:
        if not hasattr(self, "_off"):
            self._off = np.concatenate([[0], np.cumsum([len(x) for x in self.nbr])])
        k = int(np.searchsorted(self.nbr[v], j))
        return int(self._off[v]) + k

    def _merge(self, R: np.ndarray, nb: List[int]) -> List[int]:
        """Candidates within THR of each other are counted once; the first one is kept."""
        keep: List[int] = []
        for j in nb:
            if keep:
                a = self.mis(np.repeat(R[j][None], len(keep), 0), R[keep])
                if (a <= THR).any():
                    continue
            keep.append(j)
        return keep

    # ---- scores and tries
    def score3(self, Rv: np.ndarray, v: int) -> Tuple[float, float, float]:
        """(strict cost, hit ratio, confidence) of one orientation at voxel ``v``."""
        info = self.hooks.cost_fn.evaluate(
            np.asarray(Rv, np.float64), self.hooks.vertices(v), self.mic.voxels[v].phase
        )
        hit, conf = self.bfs._ratios(info)
        return float(info.cost), float(hit), float(conf)


def try_fit(hooks: Any, bfs: Any, v: int, R0: np.ndarray, rng: Any) -> Dict[str, Any]:
    """One try: the CMA fit of step 4 (first step 0.05 degrees, at most 60 evaluations, the start's
    own included, plus the final evaluation), as ``run_fine.Fitter.one`` calls it for arm C."""
    rec = hooks.rec
    vox = hooks.mic.voxels[v]
    old = rec.params.local_optimizer
    rec.params.local_optimizer = "cma"
    try:
        r = rec.local_optimization(
            voxel_vertices=hooks.vertices(v),
            phase_index=vox.phase,
            initial_orientation=np.asarray(R0, np.float64),
            rng=rng,
            cma_sigma0_deg=RF.CMA_SIGMA0_DEG,
            cma_max_evals=RF.CMA_MAX_EVALS,
        )
    finally:
        rec.params.local_optimizer = old
    hit, conf = bfs._ratios(r.overlap_info)
    return dict(
        end=np.asarray(r.orientation, np.float64),
        cost=float(r.cost),
        hit=float(hit),
        conf=float(conf),
        evals=int(rec.last_local_optimization_evals),
    )
