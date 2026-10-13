"""scripts/mc_mechanism/mc_trace.py: TracedMC reproduces MCOptimizer; its trace obeys the rules."""

import os
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).parent.parent
for sub in (
    "scripts/mc_mechanism",
    "scripts/finisher_diagnosis",
    "scripts/nn_hybrid",
    "scripts",
    "scripts/common",
    "benchmarks",
):
    sys.path.insert(0, str(ROOT / sub))


def _trace():
    keep = dict(os.environ)
    try:
        return __import__("mc_trace")
    finally:
        for k in set(os.environ) - set(keep):
            del os.environ[k]


class _Toy:
    def __init__(self, target: np.ndarray) -> None:
        self.target = target

    def evaluate(self, R, vertices, phase):
        ang = np.linalg.norm(Rotation.from_matrix(np.asarray(R, float) @ self.target.T).as_rotvec())
        return SimpleNamespace(cost=float(np.floor(ang * 200) / 200))


@pytest.mark.parametrize("seed", [0, 1, 2, 3])
def test_traced_mc_matches_mcoptimizer_and_obeys_rules(seed):
    T = _trace()
    from icenine.orientation_search import MCOptimizer

    target = Rotation.from_rotvec([0.01, -0.02, 0.015]).as_matrix()
    start = np.eye(3)
    kw = dict(angular_box_side=0.006, angular_step=0.0024, max_mc_steps=300, max_restarts=2)
    ref = MCOptimizer(_Toy(target), None, 0, np.random.default_rng(seed)).optimize(start, **kw)
    mc = T.TracedMC(_Toy(target), None, 0, np.random.default_rng(seed))
    mc.R_true = target
    res = mc.optimize(start, **kw)
    assert np.array_equal(ref.orientation, res.orientation) and ref.cost == res.cost
    tr = mc.traces[0]
    ev, step, merg, ns = tr["event"], tr["step_deg"], tr["min_erg"], tr["n_since"]
    ran = ev != -1
    log = mc.mc_logs[0]
    assert int(ran.sum()) == log["steps_run"]
    # block-end events sit on the last step of a block; blocks have the fixed length min_erg
    ends = np.nonzero((ev == 2) | (ev == 3) | (ev == 4))[0]
    assert len(ends) == log["n_blocks"]
    lens = np.diff(np.concatenate([[-1], ends]))
    assert np.all(lens[:-1] == merg[0])
    assert lens[-1] <= merg[0]
    # a global improvement (successful block) halves the next block's step; a failed block resets
    # it to the initial step; the successive-failure count is 0 after a success
    for t in ends:
        if t + 1 < len(ev) and ran[t + 1]:
            if ev[t] == 2:
                assert step[t + 1] == 0.5 * step[t] and ns[t] == 0
            else:
                assert step[t + 1] == step[0]
    assert int((ev == 2).sum()) == log["n_accept"]
    # exhaustion fires only when the successive failures exceed max_restarts (2)
    for t in np.nonzero(ev == 4)[0]:
        assert ns[t] == kw["max_restarts"] + 1
    assert log["stop"] != 1 or ev[ends[-1]] == 4
