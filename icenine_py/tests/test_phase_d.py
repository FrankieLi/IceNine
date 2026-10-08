"""scripts/phase_d helpers: grain grouping, FZ reduction and separation checks, mic round trip,
detector-noise model and ASCII image I/O."""

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT / "scripts" / "phase_d"))

import grains as G  # noqa: E402
import noise as N  # noqa: E402
from icenine.mic_file import MicFile  # noqa: E402
from icenine.symmetry import create_cubic_symmetry  # noqa: E402
from icenine.orientation_search import get_symmetry_quaternions  # noqa: E402

SYM = get_symmetry_quaternions(create_cubic_symmetry(1.0))


def test_group_grains_first_appearance_order() -> None:
    e = np.array([[10, 20, 30], [1, 2, 3], [10, 20, 30], [1, 2, 3], [5, 5, 5]], dtype=float)
    gid, uniq = G.group_grains(e)
    assert gid.tolist() == [0, 1, 0, 1, 2]
    assert uniq.tolist() == [[10, 20, 30], [1, 2, 3], [5, 5, 5]]


def test_misorientation_symmetry_reduced() -> None:
    q = G.euler_to_quat(np.array([[30.0, 40.0, 50.0]]))
    # a cubic symmetry copy of the same orientation has zero misorientation
    q2 = np.array([G._qmul(q[0], s) for s in SYM])
    d = G.misorientation_matrix_deg(q, q2, SYM)
    assert d.max() < 1e-4
    # 5 deg about z, no symmetry help
    r = G.euler_to_quat(np.array([[35.0, 40.0, 50.0]]))
    assert G.misorientation_matrix_deg(q, r, SYM)[0, 0] == pytest.approx(5.0, abs=1e-3)


def test_uniform_draws_reduce_into_fz() -> None:
    q = G.reduce_quats(G.draw_uniform_quats(200, np.random.default_rng(0)), SYM)
    # reduced: |w| is already the maximum over the group
    for x in q:
        assert max(abs(G._qmul(x, s)[0]) for s in SYM) <= abs(x[0]) + 1e-9
    assert np.allclose(np.linalg.norm(q, axis=1), 1.0)


def test_draw_separated_respects_min_separation_and_is_deterministic() -> None:
    old = G.reduce_quats(G.draw_uniform_quats(30, np.random.default_rng(1)), SYM)
    adj = [(i, i + 1) for i in range(29)]
    # a huge exclusion forces redraws; 15 deg is still satisfiable
    a, st = G.draw_separated(30, old, adj, SYM, seed=3, min_sep_deg=15.0)
    b, _ = G.draw_separated(30, old, adj, SYM, seed=3, min_sep_deg=15.0)
    assert np.array_equal(a, b)
    assert st["min_new_vs_old_deg"] >= 15.0 and st["min_neighbour_new_deg"] >= 15.0
    assert st["n_redrawn"] > 0
    c, _ = G.draw_separated(30, old, adj, SYM, seed=4, min_sep_deg=15.0)
    assert not np.array_equal(a, c)


def test_grain_adjacency() -> None:
    pos = np.array([[0.0, 0], [1, 0], [2, 0], [10, 0]])
    grain = np.array([0, 0, 1, 2])
    assert G.grain_adjacency(pos, grain, 1.5) == [(0, 1)]


def test_mic_roundtrip(tmp_path: Path) -> None:
    src = tmp_path / "a.mic"
    src.write_text(
        "0.600000\n"
        "-0.295312\t0.008119\t0.000000\t1\t6\t1\t197.642000\t48.197400\t134.216000\t1.000000\n"
        "-0.290625\t-0.519615\t0.000000\t2\t6\t1\t197.642000\t48.197400\t134.216000\t1.000000\n"
    )
    new = np.array([[12.5, 33.25, 301.125]])
    dst = tmp_path / "b.mic"
    G.write_mic_with_euler(str(src), str(dst), np.array([0, 0]), new)
    mic = MicFile.read(str(dst))
    assert len(mic.voxels) == 2
    q_back = G.euler_to_quat(np.array([G.matrix_to_euler(v.orientation) for v in mic.voxels]))
    q_exp = G.euler_to_quat(np.repeat(new, 2, axis=0))
    ang = np.degrees(2 * np.arccos(np.clip(np.abs((q_back * q_exp).sum(1)), 0, 1)))
    assert ang.max() < 1e-3
    assert [v.points_up for v in mic.voxels] == [True, False]


def test_noise_is_deterministic_and_only_touches_expected_pixels() -> None:
    img = np.zeros((128, 128), dtype=np.float32)
    for cy, cx in ((20, 20), (60, 90), (100, 40)):
        img[cy : cy + 4, cx : cx + 4] = 5.0
    p = N.NoiseParams(p_miss=0.0, p_flip=0.0, p_hot=0.0, p_blob=0.0)
    out, c = N.add_detector_noise(img, np.random.default_rng(0), p)
    assert np.array_equal(out, img) and c["n_spots"] == 3
    allmiss = N.NoiseParams(p_miss=1.0, p_flip=0.0, p_hot=1.0, p_blob=1.0)
    out, c = N.add_detector_noise(img, np.random.default_rng(0), allmiss)
    assert out.sum() == 0 and c["n_missed"] == 3 and c["n_hot"] == 0
    full = N.NoiseParams()
    a, _ = N.add_detector_noise(img, np.random.default_rng(7), full)
    b, _ = N.add_detector_noise(img, np.random.default_rng(7), full)
    assert np.array_equal(a, b)
    hot = N.NoiseParams(p_miss=0.0, p_flip=0.0, p_hot=1.0, p_blob=0.0)
    out, c = N.add_detector_noise(img, np.random.default_rng(0), hot)
    assert c["n_hot"] >= 1 and (out > 0).sum() == (img > 0).sum() + c["n_hot"]
    assert (img > 0).sum() == 48  # input untouched


def test_noise_flip_rates() -> None:
    img = np.zeros((400, 400), dtype=np.float32)
    img[100:300, 100:300] = 1.0
    p = N.NoiseParams(p_miss=0.0, p_flip=0.05, p_hot=0.0, p_blob=0.0)
    _, c = N.add_detector_noise(img, np.random.default_rng(0), p)
    assert c["n_dropped_px"] == pytest.approx(0.05 * 200 * 200, rel=0.1)
    assert c["n_grown_px"] == pytest.approx(0.05 * 800, rel=0.3)


def test_ascii_image_roundtrip(tmp_path: Path) -> None:
    img = np.zeros((16, 20), dtype=np.float32)
    img[3, 5], img[15, 19], img[0, 0] = 1.5, 33793.3, 2.0
    f = tmp_path / "x.d0"
    assert N.write_ascii_image(str(f), img) == 3
    assert f.read_text().splitlines()[0].count(",") == 2 and "#" not in f.read_text()
    back = N.read_ascii_image(str(f), 16, 20)
    assert np.allclose(back, img, rtol=1e-5)
    empty = tmp_path / "e.d0"
    N.write_ascii_image(str(empty), np.zeros((4, 4), dtype=np.float32))
    assert N.read_ascii_image(str(empty), 4, 4).sum() == 0


def test_python_loader_reads_written_image(tmp_path: Path) -> None:
    from icenine.image_data import ImageData

    img = np.zeros((16, 20), dtype=np.float32)
    img[3, 5] = 7.0
    f = tmp_path / "y.d0"
    N.write_ascii_image(str(f), img)
    im = ImageData(16, 20)
    im.load_ascii(str(f))
    assert float(im._pixels_dense[3, 5]) == 7.0 and int(im.num_nonzero) == 1
