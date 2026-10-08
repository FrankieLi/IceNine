"""Scripted 14 x 14 voxel, 4 grain grid for BFSReconstruction tests (no images, no physics).

The fake reconstructor draws from the BFS generator in every full and local fit, so any change in
the number or order of fits changes the final generator state: the golden in
``data/bfs_grid_golden.npz`` (recorded from the commit BEFORE the Phase D library options) pins
the whole fit sequence, not just the end result. ``run_grid`` takes the package modules as
arguments so the golden can be recorded with the parent commit's package.
"""

from types import SimpleNamespace
from typing import Any, Dict, List, Tuple

import numpy as np
from scipy.spatial.transform import Rotation

N = 14
POS = [(float(i), float(j)) for i in range(N) for j in range(N)]
KEY = {p: k for k, p in enumerate(POS)}


def _grain(p: Tuple[float, float]) -> int:
    x, y = p
    if x + 0.3 * y < 5:
        return 0
    if y < 7:
        return 1
    return 2 if x < 10 else 3


GRAIN = [_grain(p) for p in POS]
_GR = [
    Rotation.from_rotvec(np.radians(v)).as_matrix()
    for v in ([0, 0, 0], [0, 0, 6], [3, 0, 12], [0, 5, 0])
]
_jit = np.random.default_rng(99)
TRUTH = [
    Rotation.from_rotvec(np.radians(_jit.normal(0, 0.2, 3))).as_matrix() @ _GR[g] for g in GRAIN
]


def ang(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.degrees(np.linalg.norm(Rotation.from_matrix(a @ b.T).as_rotvec())))


def info(hit: float) -> Any:
    return SimpleNamespace(
        pixel_overlap=int(round(hit * 10000)),
        pixel_on_detector=10000,
        peak_overlap=int(round(hit * 1000)),
        peak_on_detector=1000,
        cost=1 - hit,
    )


def make_fake(candidate_cls: Any, log: List[Tuple[str, int, float]]) -> Any:
    """A reconstructor stand-in: full fits are right 85% of the time (else 0.5 deg-scale noise),
    local fits converge to the truth within 2 deg of it, else stay near the start."""

    class Fake:
        def __init__(self, setup: Any) -> None:
            self.last_local_optimization_evals = 0
            self._c = (0, 0, 0)

        @property
        def last_eval_counts(self) -> Tuple[int, int, int]:
            return self._c

        def _i(self, v: Any) -> int:
            return KEY[(round(float(v[0][0]), 6), round(float(v[0][1]), 6))]

        def reconstruct_voxel(self, voxel_vertices: Any, phase_index: int, rng: Any) -> Any:
            i = self._i(voxel_vertices)
            u = rng.random()
            R = TRUTH[i] if u > 0.15 else Rotation.from_rotvec(rng.normal(0, 0.5, 3)).as_matrix()
            self._c = (10, 5, 15)
            log.append(("full", i, float(u)))
            return candidate_cls(R.astype(np.float64), 0.1, info(0.9))

        def evaluate_overlap(self, R: np.ndarray, v: Any, phase: int) -> Any:
            return info(max(0.0, 0.97 - 0.05 * ang(R, TRUTH[self._i(v)])))

        def local_optimization(
            self,
            voxel_vertices: Any,
            phase_index: int,
            initial_orientation: np.ndarray,
            rng: Any = None,
            **kw: Any,
        ) -> Any:
            i = self._i(voxel_vertices)
            n = rng.normal(0, 0.05, 3)
            base = TRUTH[i] if ang(initial_orientation, TRUTH[i]) < 2.0 else initial_orientation
            R = Rotation.from_rotvec(np.radians(n)).as_matrix() @ base
            hit = max(0.0, 0.97 - 0.05 * ang(R, TRUTH[i]) - abs(n[0]))
            self.last_local_optimization_evals = 7
            log.append(("local", i, float(n[0])))
            return candidate_cls(R.astype(np.float64), 1 - hit, info(hit))

    return Fake


def run_grid(
    reconstructor_module: Any,
    mic_module: Any,
    candidate_cls: Any,
    params_cls: Any,
    params: Dict[str, Any],
    seed: int = 5,
) -> Dict[str, Any]:
    """Run a BFS on the scripted grid with ``params`` and return the outcome (and the BFS)."""
    log: List[Tuple[str, int, float]] = []
    original = reconstructor_module.AdaptiveVoxelReconstructor
    reconstructor_module.AdaptiveVoxelReconstructor = make_fake(candidate_cls, log)
    voxels = [
        mic_module.Voxel(position=np.array([x, y, 0.0]), orientation=np.eye(3), side_length=0.6)
        for x, y in POS
    ]
    mic = mic_module.MicFile(voxels)
    setup = SimpleNamespace(
        sample=SimpleNamespace(get_mic=lambda: mic),
        config=SimpleNamespace(min_acceleration_threshold=0.8),
        search_params=params_cls(**params),
    )
    try:
        bfs = reconstructor_module.BFSReconstruction(setup)
    finally:
        reconstructor_module.AdaptiveVoxelReconstructor = original
    rng = np.random.default_rng(seed)
    done = bfs.reconstruct_sample(rng=rng)
    return dict(
        done=np.array(done),
        states=np.array([int(v.reconstruction_id) for v in mic.voxels]),
        R=np.stack([v.orientation for v in mic.voxels]).astype(np.float64),
        hit=np.array([float(v.overlap_ratio) for v in mic.voxels]),
        log_kind=np.array([e[0] for e in log]),
        log_idx=np.array([e[1] for e in log]),
        log_val=np.array([e[2] for e in log]),
        rng_after=np.array(rng.random()),
        bfs=bfs,
        mic=mic,
    )


def n_wrong(out: Dict[str, Any], deg: float = 1.0) -> int:
    return sum(1 for i in range(len(POS)) if ang(out["R"][i], TRUTH[i]) > deg)
