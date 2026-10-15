"""Tests of icenine.spot_pairing on synthetic two-detector patterns, and of the opt-in
ForwardSimulation.peak_hook used by the provenance writer."""

import os
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import torch

from icenine import spot_pairing as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "eval_reduction"))
import pairing_common as pc  # noqa: E402

PW = 0.00148


def make_geoms(x0: float = 3.0, x1: float = 5.0):
    """Two ideal detectors: normal along x, columns along y, rows along -z, 2048 pixels."""
    out = []
    for x in (x0, x1):
        half = 1024 * PW
        out.append(
            sp.DetectorGeometry(
                pos=np.array([x, 0.0, 0.0]),
                org=np.array([0.0, -half, half]),
                bj=np.array([0.0, 1.0, 0.0]),
                bk=np.array([0.0, 0.0, -1.0]),
                pw=PW,
                ph=PW,
                hw=half,
                hh=half,
                n=np.array([1.0, 0.0, 0.0]),
                d=-x,
            )
        )
    return out


def ray_pixels(geoms, origin, two_theta, eta):
    """Centroid pixels (col0, row0, col1, row1) of the ray from ``origin`` (lab, z = 0)."""
    d = np.array(
        [np.cos(two_theta), np.sin(two_theta) * np.sin(eta), np.sin(two_theta) * np.cos(eta)]
    )
    pix = []
    for g in geoms:
        t = (g.distance_x() - origin[0]) / d[0]
        col, row = sp.lab_to_pixel(g, origin + t * d)
        pix += [float(col), float(row)]
    return pix


def test_pixel_round_trip():
    g = make_geoms()[0]
    col, row = np.array([10.5, 1500.25]), np.array([700.0, 33.5])
    c2, r2 = sp.lab_to_pixel(g, sp.pixel_to_lab(g, col, row))
    assert np.allclose(c2, col) and np.allclose(r2, row)


def test_line_origin_recovers_the_origin_and_sigma_grows_for_flat_rays():
    gs = make_geoms()
    origin = np.array([0.12, -0.2, 0.0])
    steep = ray_pixels(gs, origin, np.radians(15), np.radians(10))
    flat = ray_pixels(gs, origin, np.radians(15), np.radians(85))
    res = []
    for p in (steep, flat):
        lo = sp.line_origin(gs[0], gs[1], p[0], p[1], p[2], p[3], 1.0, 1.0)
        assert np.allclose(lo["origin"][:2], origin[:2], atol=1e-9)
        assert abs(lo["origin"][2]) < 1e-12 and lo["s"] < 0
        assert np.isclose(np.linalg.norm(lo["dirn"]), 1.0)
        res.append(lo)
    assert res[1]["sigma_mm"] > 5 * res[0]["sigma_mm"]
    assert res[1]["cond"] < res[0]["cond"]


def test_distance_x_is_the_beam_axis_hit_not_the_pixel_corner():
    # a tilted detector whose pixel-(0, 0) corner sits 0.1 mm closer than the beam-axis hit
    n = np.array([0.999, -0.02, 0.03])
    n = n / np.linalg.norm(n)
    pos = np.array([3.36, 0.0, 0.0])
    g = sp.DetectorGeometry(
        pos=pos,
        org=np.array([-0.1, -1.5, 2.85]),
        bj=np.array([0.0, 1.0, 0.0]),
        bk=np.array([0.0, 0.0, -1.0]),
        pw=PW,
        ph=PW,
        hw=1024 * PW,
        hh=1024 * PW,
        n=n,
        d=-float(n @ pos),
    )
    assert np.isclose(g.distance_x(), 3.36)
    assert not np.isclose(g.distance_x(), (g.pos + g.org)[0])


def test_detector1_tolerance_is_its_own_pixel_not_scaled_by_l1_over_l0():
    gs, origins, tth, pix = frame_pattern(n=1)
    hull = sp.SampleHull(np.array([[-0.5, -0.5], [0.5, -0.5], [0.5, 0.5], [-0.5, 0.5]]))
    c0 = {"c_col": pix[:, 0] - 0.5, "c_row": pix[:, 1] - 0.5, "npix": np.array([40.0])}
    c1 = {"c_col": pix[:, 2] - 0.5, "c_row": pix[:, 3] - 0.5, "npix": np.array([40.0])}
    p = sp.PairingParams(voxel_size_mm=0.009, omega_step=np.radians(1.0), centroid_tol_px=1.0)
    r = sp.pair_frame(gs, hull, 0.0, c0, c1, p)
    lo = sp.line_origin(gs[0], gs[1], pix[0, 0], pix[0, 1], pix[0, 2], pix[0, 3], 1.0, 1.0)
    assert len(r["i0"]) == 1
    assert np.isclose(r["sigma_mm"][0], lo["sigma_mm"])


def test_frame_check_accepts_identity_and_rejects_other_setups():
    good_cfg = SimpleNamespace(beam_direction=np.array([1.0, 0, 0]), beam_deflection_chi_laue=0.0)
    good = SimpleNamespace(sample_to_lab_matrix=torch.eye(4))
    pc.check_frame_assumption(good_cfg, good)
    moved = torch.eye(4)
    moved[0, 3] = 0.01
    with pytest.raises(ValueError):
        pc.check_frame_assumption(good_cfg, SimpleNamespace(sample_to_lab_matrix=moved))
    tilted = SimpleNamespace(beam_direction=np.array([1.0, 0.01, 0]), beam_deflection_chi_laue=0.0)
    with pytest.raises(ValueError):
        pc.check_frame_assumption(tilted, good)
    bent = SimpleNamespace(beam_direction=np.array([1.0, 0, 0]), beam_deflection_chi_laue=0.001)
    with pytest.raises(ValueError):
        pc.check_frame_assumption(bent, good)


def test_sample_hull_distance():
    hull = sp.SampleHull(np.array([[-1, -1], [1, -1], [1, 1], [-1, 1.0]]))
    assert hull.outside_distance(np.array([0.2, 0.3])) == 0.0
    assert np.isclose(hull.outside_distance(np.array([1.5, 0.0])), 0.5)
    assert np.isclose(hull.radius, np.sqrt(2.0))


def test_rotation_helpers_are_inverse():
    q = np.array([[0.3, -0.1], [0.0, 0.5]])
    assert np.allclose(sp.lab_to_sample_xy(sp.sample_to_lab_xy(q, 0.7), 0.7), q)
    assert np.allclose(
        sp.rot_z(np.array(0.4))[:2, :2] @ np.array([1.0, 0.0]), [np.cos(0.4), np.sin(0.4)]
    )


def frame_pattern(seed=0, n=6):
    """n rays with distinct origins inside a 1 mm hull; detector-1 spots in a shuffled order."""
    rng = np.random.default_rng(seed)
    gs = make_geoms()
    origins = np.zeros((n, 3))
    origins[:, :2] = rng.uniform(-0.4, 0.4, size=(n, 2))
    tth = np.radians(rng.uniform(6, 14, n))
    eta = np.radians(rng.uniform(30, 80, n)) * rng.choice([-1, 1], n)
    pix = np.array([ray_pixels(gs, origins[i], tth[i], eta[i]) for i in range(n)])
    return gs, origins, tth, pix


def test_pair_frame_recovers_shuffled_pairs():
    gs, origins, tth, pix = frame_pattern()
    hull = sp.SampleHull(np.array([[-0.5, -0.5], [0.5, -0.5], [0.5, 0.5], [-0.5, 0.5]]))
    n = len(pix)
    perm = np.random.default_rng(3).permutation(n)
    c0 = {"c_col": pix[:, 0] - 0.5, "c_row": pix[:, 1] - 0.5, "npix": np.full(n, 40.0)}
    c1 = {"c_col": pix[perm, 2] - 0.5, "c_row": pix[perm, 3] - 0.5, "npix": np.full(n, 40.0)}
    params = sp.PairingParams(
        voxel_size_mm=0.009, omega_step=np.radians(1.0), ring_two_theta=np.sort(tth)
    )
    r = sp.pair_frame(gs, hull, 0.0, c0, c1, params)
    ch = r["chosen"].astype(bool)
    got = {(int(a), int(b)) for a, b in zip(r["i0"][ch], r["i1"][ch])}
    inv = np.argsort(perm)
    assert got == {(i, int(inv[i])) for i in range(n)}
    for a, b, o in zip(r["i0"][ch], r["i1"][ch], r["origin_lab"][ch]):
        assert np.allclose(o[:2], origins[a, :2], atol=1e-6)


def test_ring_gate_and_hull_gate_reject_wrong_candidates():
    gs, origins, tth, pix = frame_pattern(n=4)
    hull = sp.SampleHull(np.array([[-0.5, -0.5], [0.5, -0.5], [0.5, 0.5], [-0.5, 0.5]]))
    c0 = {"c_col": pix[:1, 0] - 0.5, "c_row": pix[:1, 1] - 0.5, "npix": np.array([40.0])}
    # a detector-1 spot shifted far along the columns: the line leaves the hull
    far = {"c_col": pix[:1, 2] - 0.5 + 900.0, "c_row": pix[:1, 3] - 0.5, "npix": np.array([40.0])}
    base = dict(voxel_size_mm=0.009, omega_step=np.radians(1.0))
    r = sp.pair_frame(gs, hull, 0.0, c0, far, sp.PairingParams(**base), return_gates=True)
    assert len(r["i0"]) == 0 and not r["gate_hull"][0, 0]
    # right place but the wrong ring list: ring gate rejects
    c1 = {"c_col": pix[:1, 2] - 0.5, "c_row": pix[:1, 3] - 0.5, "npix": np.array([40.0])}
    wrong = sp.PairingParams(ring_two_theta=np.array([tth[0] + 0.1]), **base)
    r = sp.pair_frame(gs, hull, 0.0, c0, c1, wrong, return_gates=True)
    assert len(r["i0"]) == 0 and not r["gate_ring"][0, 0]
    ok = sp.PairingParams(ring_two_theta=np.array([tth[0]]), **base)
    assert len(sp.pair_frame(gs, hull, 0.0, c0, c1, ok)["i0"]) == 1


def test_size_gate_skips_edge_spots():
    gs, origins, tth, pix = frame_pattern(n=1)
    hull = sp.SampleHull(np.array([[-0.5, -0.5], [0.5, -0.5], [0.5, 0.5], [-0.5, 0.5]]))
    c0 = {"c_col": pix[:, 0] - 0.5, "c_row": pix[:, 1] - 0.5, "npix": np.array([100.0])}
    c1 = {"c_col": pix[:, 2] - 0.5, "c_row": pix[:, 3] - 0.5, "npix": np.array([20.0])}
    p = sp.PairingParams(voxel_size_mm=0.009, omega_step=0.017, shape_weight=1.0)
    assert len(sp.pair_frame(gs, hull, 0.0, c0, c1, p)["i0"]) == 0
    c1["edge"] = np.array([True])
    assert len(sp.pair_frame(gs, hull, 0.0, c0, c1, p)["i0"]) == 1


def test_assign_pairs_caps_and_tie_rule():
    i0 = np.array([0, 0, 1])
    i1 = np.array([0, 1, 1])
    cost = np.array([0.0, 0.0, 0.5])
    p = sp.PairingParams(max_partners=1)
    ch = sp.assign_pairs(i0, i1, cost, np.array([10.0, 10.0]), np.array([10.0, 10.0]), p)
    # equal costs: the first listed wins (stable); spot 0 is then used up
    assert ch.tolist() == [True, False, True]
    p2 = sp.PairingParams(max_partners=1, multi_npix=50, multi_partners=2)
    ch2 = sp.assign_pairs(i0, i1, cost, np.array([100.0, 10.0]), np.array([10.0, 10.0]), p2)
    # spot 0 of detector 0 may take two partners; detector-1 spot 1 is then used up for spot 1
    assert ch2.tolist() == [True, True, False]


def test_attribute_peaks_radius_and_per_peak_radius():
    origins = np.array([[0.0, 0.0], [0.02, 0.0], [0.5, 0.5]])
    vox = np.array([[0.0, 0.0], [0.5, 0.5]])
    out = sp.attribute_peaks(origins, vox, 0.03)
    assert sorted(out[0].tolist()) == [0, 1] and out[1].tolist() == [2]
    out = sp.attribute_peaks(origins, vox, np.array([0.01, 0.01, 0.01]))
    assert out[0].tolist() == [0]


EX = Path(__file__).resolve().parents[2] / "Examples" / "Example2.ManyGrains"


@pytest.mark.skipif(not (EX / "SimInput").exists(), reason="Example2.ManyGrains not present")
def test_peak_hook_is_off_by_default_and_matches_the_rasterized_peaks():
    cwd = os.getcwd()
    os.chdir(EX)
    try:
        from icenine.config_file import ConfigFile
        from icenine.forward_simulation import ForwardSimulation
        from icenine.sample import Sample
        from icenine.simulation import Simulation

        class Sink:
            def __init__(self):
                self.n = 0

            def add_triangle_scanline(self, *a, **k):
                self.n += 1

        cfg = ConfigFile.from_file(str(EX / "ConfigFiles" / "Example2.Simulation.config"))
        cfg.sample_filename = "SimInput/rand_500grains_1mm_neworient_s0.mic"
        sim = ForwardSimulation(cfg)
        assert sim.peak_hook is None
        sim.exp_setup.initialize_experiment()
        dets = sim.exp_setup.get_detector_list()
        rmap = sim.exp_setup.get_range_to_index_map()
        sim.simulator = Simulation(sim.exp_setup)
        sample = Sample()
        sim.exp_setup.initialize_sample(sample, dets[0])
        sample.get_mic().voxels = sample.get_mic().voxels[:3]
        n_om = len(sim.exp_setup.get_omega_range_list())

        def run():
            sinks = [[Sink() for _ in dets] for _ in range(n_om)]
            sim._simulate_peaks_batched(sinks, dets, sample, rmap)
            return sum(s.n for row in sinks for s in row)

        n_off = run()
        calls = []
        sim.peak_hook = lambda v, r, oi, om: calls.append((v, r, oi, om))
        n_on = run()
        assert n_on == n_off > 0
        assert len(calls) * len(dets) == n_on
        assert {c[0] for c in calls} <= {0, 1, 2}
        assert all(0 <= c[2] < n_om for c in calls)

        # several chunks (one voxel each): the hook's global voxel index must name the voxel of
        # the chunk, checked against a single-voxel run of each voxel
        mic = sample.get_mic()
        three = list(mic.voxels)

        def calls_of(batch_size, voxel_offset=0):
            got = []
            sim.peak_hook = lambda v, r, oi, om: got.append((v + voxel_offset, r, oi, om))
            sinks = [[Sink() for _ in dets] for _ in range(n_om)]
            sim._simulate_peaks_batched(sinks, dets, sample, rmap, batch_size=batch_size)
            return sorted(got)

        def same(a, b):
            # identical (voxel, reflection, image); omega equal up to float32 batching noise
            assert [c[:3] for c in a] == [c[:3] for c in b]
            assert np.allclose([c[3] for c in a], [c[3] for c in b], rtol=0.0, atol=1e-6)

        chunked = calls_of(1)
        same(chunked, calls_of(None))
        for k, vx in enumerate(three):
            mic.voxels = [vx]
            single = calls_of(None, voxel_offset=k)
            assert single
            same(single, [c for c in chunked if c[0] == k])
        mic.voxels = three
    finally:
        os.chdir(cwd)
