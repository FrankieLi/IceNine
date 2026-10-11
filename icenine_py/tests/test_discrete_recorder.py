"""run_discrete_search_spaced: the optional `recorder` hook leaves the result unchanged and
records, per FZ point, the screen (peak_overlap per local-grid index), the re-cost and the
indices kept by the spacing filter (Phase D level-0 trace)."""

import math
from typing import Any, Dict, List, Tuple

import numpy as np

from icenine.cost_functions import OverlapInfo
from icenine.orientation_search import run_discrete_search_spaced


class _StubCost:
    """Deterministic stand-in: the overlap depends on the orientation's z-rotation angle."""

    def __init__(self, pixel_radius: int) -> None:
        self.pixel_radius = pixel_radius
        self.max_q = 5.0

    def evaluate(self, orientation: np.ndarray, voxel_vertices: Any, phase_index: int = 0):
        ang = math.atan2(orientation[1, 0], orientation[0, 0])
        k = int(round(abs(ang) * 1000)) % 4
        info = OverlapInfo()
        info.peak_overlap = k if self.pixel_radius > 0 else 0
        info.peak_on_detector = 3
        info.quality = 0.1 * k
        info.n_quality_points = 1
        return info


def _rz(a: float) -> np.ndarray:
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def _args() -> Dict[str, Any]:
    fz = np.stack([_rz(0.01 * i) for i in range(6)])
    lg = np.stack([_rz(0.002 * j) for j in range(5)])
    return dict(
        global_cost_fn=_StubCost(3),
        local_cost_fn=_StubCost(0),
        fz_orientations=fz,
        local_grid=lg,
        voxel_vertices=None,
        angular_radius=math.radians(0.3),
        symmetry_quats=np.array([[1.0, 0.0, 0.0, 0.0]]),
    )


def _run(recorder=None) -> Tuple[List[Any], List[Tuple[str, Dict[str, Any]]]]:
    events: List[Tuple[str, Dict[str, Any]]] = []
    kw = _args()
    if recorder:
        kw["recorder"] = lambda e, d: events.append((e, d))
    return run_discrete_search_spaced(**kw), events


def test_recorder_default_path_unchanged():
    plain, _ = _run()
    traced, events = _run(recorder=True)
    assert len(plain) == len(traced)
    for p, t in zip(plain, traced):
        assert np.array_equal(p.orientation, t.orientation)
        assert p.cost == t.cost
    assert events[0][0] == "level" and events[-1][0] == "level_end"


def test_recorder_contents_consistent():
    out, events = _run(recorder=True)
    fz = [d for e, d in events if e == "fz"]
    assert [d["fz"] for d in fz] == list(range(6))
    for d in fz:
        assert len(d["peak_overlap"]) == 5
        assert d["j"] == [j for j, o in enumerate(d["peak_overlap"]) if o > 0]
        assert len(d["cost"]) == len(d["j"])
        assert set(d["kept"]) <= set(d["j"])
    assert events[0][1]["n_fz"] == 6 and events[0][1]["n_local"] == 5
    assert events[-1][1]["n_candidates"] == len(out) == sum(len(d["kept"]) for d in fz)
