"""scripts/finisher_diagnosis: LoggedMC reproduces MCOptimizer, geometry helpers, case selection."""

import os
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).parent.parent
for sub in ("scripts/finisher_diagnosis", "scripts/nn_hybrid", "scripts", "benchmarks"):
    sys.path.insert(0, str(ROOT / sub))


def _diag():
    keep = dict(os.environ)
    try:
        return __import__("diagnose")
    finally:
        for k in set(os.environ) - set(keep):
            del os.environ[k]


class _Toy:
    """Cost = rough function of the orientation's distance to a target (no pixels needed)."""

    def __init__(self, target: np.ndarray) -> None:
        self.target = target
        self.n = 0

    def evaluate(self, R, vertices, phase):
        self.n += 1
        ang = np.linalg.norm(Rotation.from_matrix(np.asarray(R, float) @ self.target.T).as_rotvec())
        return SimpleNamespace(cost=float(np.floor(ang * 200) / 200))


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_logged_mc_matches_mcoptimizer(seed):
    D = _diag()
    from icenine.orientation_search import MCOptimizer

    target = Rotation.from_rotvec([0.01, -0.02, 0.015]).as_matrix()
    start = np.eye(3)
    kw = dict(angular_box_side=0.006, angular_step=0.006)
    outs = []
    for cls in (MCOptimizer, D.LoggedMC):
        mc = cls(_Toy(target), None, 0, np.random.default_rng(seed))
        r1 = mc.optimize(start, max_mc_steps=300, max_restarts=2, **kw)
        r2 = mc.variance_minimizing_optimize(r1.orientation, 0.006, 300, 2, 0.0, 0.02**2)
        outs.append((r1.orientation, r1.cost, r2.orientation, r2.cost))
        last = mc
    for a, b in zip(*outs):
        assert np.array_equal(a, b)
    log = last.mc_logs[0]
    assert log["stop"] in (0, 1, 2) and log["steps_run"] <= 300
    assert last.vm_logs[0]["steps_taken"] >= 300 or last.vm_logs[0]["capped"] == 0


def test_vm_step_cap_stops_runaway_budget():
    D = _diag()

    class Noisy:
        def evaluate(self, R, vertices, phase):
            return SimpleNamespace(cost=float(np.random.default_rng().uniform(0.0, 1.0)))

    mc = D.LoggedMC(Noisy(), None, 0, np.random.default_rng(0))
    mc.vm_step_cap = 500
    mc.variance_minimizing_optimize(np.eye(3), 0.006, 100, 2, 0.0, 0.02**2)
    assert mc.vm_logs[0]["capped"] == 1


def test_geodesic_and_tiny_rotations():
    D = _diag()
    rng = np.random.default_rng(0)
    R0 = Rotation.random(random_state=1).as_matrix()
    R1 = Rotation.random(random_state=2).as_matrix()
    ts = np.array([0.0, 0.25, 0.5, 1.0])
    pts = D.geodesic_points(R0, R1, ts)
    assert np.allclose(pts[0], R0) and np.allclose(pts[-1], R1)
    total = D.angle_deg(R1, R0)
    assert np.isclose(D.angle_deg(pts[1], R0), 0.25 * total)
    assert np.isclose(D.angle_deg(pts[2], R1), 0.5 * total)
    ang = np.array([0.001, 0.05])
    rots = D.tiny_rotations(R0, ang, 5, rng)
    assert rots.shape == (2, 5, 3, 3)
    for i, a in enumerate(ang):
        for k in range(5):
            assert np.isclose(D.angle_deg(rots[i, k], R0), a, rtol=1e-6)


def test_select_cases_stratified():
    D = _diag()
    nv, nr, nd = 20, 10, 20
    raw = dict(
        ran=np.ones((nv, nr, nd, 2), bool),
        fallback=np.zeros((nv, nr, nd, 2), bool),
        fail_pass1=np.zeros((nv, nr, nd, 2), np.int8),
        voxel_indices=np.arange(nv)[::-1].copy(),
    )
    sweep_vox = np.arange(nv)
    for pipe, n in (("H3", 200), ("H0", 100)):
        cases = D.select_cases(pipe, raw, sweep_vox)
        assert sum(len(d) for _, _, d in cases) == n
        per_r = np.bincount([r for _, r, d in cases for _ in d], minlength=nr)
        assert per_r[9] == 0 and per_r[:9].min() >= n // 9
        assert all(len(set(d)) == len(d) for _, _, d in cases)
