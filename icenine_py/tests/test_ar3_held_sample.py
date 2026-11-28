"""Tests of the pure parts of the sample builder for ``held`` (AR3 step 11)."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "autoresearch" / "ar3"))
import held_sample_lib as H  # noqa: E402


def _grid(n: int = 20) -> np.ndarray:
    xs, ys = np.meshgrid(np.arange(n, dtype=float), np.arange(n, dtype=float))
    return np.stack([xs.ravel(), ys.ravel()], axis=1)


def test_cell_seeds_sorted_unique_and_depend_on_seed() -> None:
    a = H.draw_cell_seeds(400, 25, 1)
    b = H.draw_cell_seeds(400, 25, 1)
    c = H.draw_cell_seeds(400, 25, 0)
    assert len(a) == 25 and len(set(a.tolist())) == 25
    assert np.all(np.diff(a) > 0)
    assert np.array_equal(a, b)
    assert not np.array_equal(a, c)


def test_voronoi_assigns_nearest_seed_and_seeds_to_themselves() -> None:
    xy = _grid()
    seeds = H.draw_cell_seeds(len(xy), 9, 1)
    cell = H.voronoi_cells(xy, seeds)
    assert np.array_equal(cell[seeds], np.arange(9))
    d = np.linalg.norm(xy[:, None, :] - xy[seeds][None], axis=2)
    assert np.array_equal(cell, d.argmin(axis=1))


def test_connected_grains_split_a_disconnected_cell() -> None:
    xy = np.array([[0.0, 0], [1, 0], [2, 0], [10, 0], [11, 0]])
    cell = np.array([0, 0, 0, 0, 0])  # one cell, two separate pieces
    g = H.connected_grains(xy, cell, 1.2)
    assert len(np.unique(g)) == 2
    assert g[0] == g[1] == g[2] and g[3] == g[4] and g[0] != g[3]
    cell2 = np.array([0, 1, 1, 1, 1])  # a cell boundary also separates neighbours
    assert len(np.unique(H.connected_grains(xy, cell2, 1.2))) == 3


def test_every_voxel_has_a_grain() -> None:
    xy = _grid()
    seeds = H.draw_cell_seeds(len(xy), 12, 1)
    g = H.connected_grains(xy, H.voronoi_cells(xy, seeds), 1.2)
    assert len(g) == len(xy) and g.min() == 0
    assert set(np.unique(g).tolist()) == set(range(g.max() + 1))
    assert g.max() + 1 >= 12


def test_levels_same_seed_same_draw_and_only_allowed_values() -> None:
    levels = (0.0, 0.2, 0.5, 1.0, 2.0, 4.0)
    _g1, a = H.draw_levels(1, 250, levels)
    _g2, b = H.draw_levels(1, 250, levels)
    _g3, c = H.draw_levels(0, 250, levels)
    assert a == b and a != c and len(a) == 250
    assert set(a) <= set(levels) and len(set(a)) == len(levels)


def test_levels_generator_state_continues() -> None:
    g, _lv = H.draw_levels(1, 10, (0.0, 1.0))
    ref = np.random.default_rng(1)
    ref.integers(0, 2, size=10)
    assert g.random() == ref.random()


def test_grains_to_redraw_takes_the_later_grain_of_close_pairs() -> None:
    a, b = [0, 1, 4, 2], [1, 2, 3, 5]
    ang = [3.0, 6.0, 4.9, 5.0]
    assert H.grains_to_redraw(a, b, ang, 5.0) == [1, 4]
    assert H.grains_to_redraw(a, b, [9.0] * 4, 5.0) == []
