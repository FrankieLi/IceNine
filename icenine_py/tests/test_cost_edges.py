"""Task B tooling: event classification and the attribution of a change of the cost (synthetic
traces), and the invariants of a stored ray walk when its cache is present."""

import sys
from pathlib import Path
from typing import Any, Dict

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "eval_reduction"))
import cost_edges as CE  # noqa: E402

DIMS = (2048, 2048)


def trace(P: int = 4) -> Dict[str, Any]:
    """P peaks, all evaluated on 2 detectors, one vertex-free (screen) signature each, q = 1."""
    return dict(
        status=np.full(P, CE.PRESENT, np.int8),
        frame=np.arange(P, dtype=np.int32),
        hit=np.ones((P, 2), np.uint8),
        sig=np.full((P, 2, 1, 2), 100, np.int32),
        cnt=np.tile(np.array([[5, 5, 1, 1, 2]], np.int32), (P, 1)),
        q=np.ones(P),
        rays=np.zeros((P, 8), np.int32),
        q_sum=float(P),
        n_q=P,
    )


def test_no_change_is_no_event() -> None:
    a, b = trace(), trace()
    ev = CE.classify(a, b, DIMS, 3)
    assert (ev == -1).all()
    d, parts, flags = CE.attribute(a, b, ev)
    assert d == 0 and parts == {} and flags["n_den_only"] == 0


def test_event_types_and_priority() -> None:
    a, b = trace(5), trace(5)
    b["frame"][0] += 1  # frame
    b["sig"][1, 0, 0, 0] += 1  # pixel
    b["sig"][2, 1, 0, 1] = 5000  # leaves the detector area: detector (wins over the pixel change)
    b["status"][3] = 1  # eta cut: peak not evaluated any more
    b["cnt"][4, 0] += 1  # counts only: subpixel
    ev = CE.classify(a, b, DIMS, 3)
    names = [CE.EVENT_TYPES[i] if i >= 0 else None for i in ev]
    assert names == ["frame", "pixel", "detector", "eta_cut", "subpixel"]


def test_attribution_sums_to_the_change_of_the_cost() -> None:
    a, b = trace(4), trace(4)
    # peak 0: frame crossing drops its quality to 0; peak 1: pixel event lowers it to 0.5
    b["frame"][0] += 1
    b["q"][0] = 0.0
    b["sig"][1, 0, 0, 0] += 1
    b["q"][1] = 0.5
    b["q_sum"], b["n_q"] = 2.5, 4
    ev = CE.classify(a, b, DIMS, 3)
    d, parts, _ = CE.attribute(a, b, ev)
    assert d == pytest.approx(1 - 2.5 / 4 - 0.0)  # cost 0 -> 0.375
    assert sum(parts.values()) == pytest.approx(d)
    assert parts["frame"] == pytest.approx(0.25) and parts["pixel"] == pytest.approx(0.125)


def test_denominator_change_with_no_numerator_change() -> None:
    a, b = trace(5), trace(5)
    a["status"][4] = 1  # the fifth peak is eta-cut in a ...
    a["q"][4] = 0.0
    a["q_sum"], a["n_q"] = 4.0, 4
    b["q_sum"], b["n_q"] = (
        5.0,
        5,
    )  # ... and evaluated in b with quality 1: the mean stays 1, N changes
    ev = CE.classify(a, b, DIMS, 3)
    d, parts, flags = CE.attribute(a, b, ev)
    assert d == pytest.approx(0.0)
    assert flags["dN"] == 1
    assert sum(parts.values()) == pytest.approx(0.0, abs=1e-12)
    assert set(parts) == {"eta_cut", "denominator"}


def test_ray_radii_windows() -> None:
    w = CE.ray_radii()
    assert set(w) == {"fin", "mc", "grid"}
    assert w["fin"].shape == (1, 101) and w["grid"].shape == (4, 101)
    assert w["grid"][-1, -1] == pytest.approx(5.0)


STORED = ROOT / "benchmarks" / "eval_reduction" / "B" / "cache" / "full" / "v15901_d0.npz"


@pytest.mark.skipif(not STORED.exists(), reason="stored ray walk (gitignored cache) missing")
def test_stored_walk_parts_sum_to_cost_change() -> None:
    z = np.load(STORED)
    for combo in ("full|cpp", "screen|old"):
        c = z[f"{combo}|cost|dense"].astype(np.float64)
        parts = z[f"{combo}|parts|dense"].astype(np.float64)
        assert np.abs(parts.sum(1) - np.diff(c)).max() < 1e-4
    assert float(z["maxdiff"]) < 1e-5
