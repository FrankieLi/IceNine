"""Opt-in strain in the stock render path (``ForwardSimulation.apply_deformation``, GO5).

Default off is unchanged by a deformation in the .mic; on, the rendered triangle vertices of a
strained voxel equal ``SpotModel.predict(..., strain=)`` positions; the .mic writer keeps its
default format and writes more digits on request."""

import os
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import pytest
import torch

EX = Path(__file__).resolve().parents[2] / "Examples" / "Example2.ManyGrains"
pytestmark = pytest.mark.skipif(
    not (EX / "SimInput").exists(), reason="Example2.ManyGrains not present"
)


class _Sink:
    """Records the vertices of every triangle drawn into it, with the peak context."""

    def __init__(self, ctx: Dict[str, Any], det: int, log: List[Any]) -> None:
        self.ctx, self.det, self.log = ctx, det, log

    def add_triangle_scanline(self, v0: Any, v1: Any, v2: Any, intensity: float) -> None:
        self.log.append(
            (self.ctx["peak"], self.det, np.array([tuple(v0), tuple(v1), tuple(v2)], dtype=np.float64), intensity)
        )


def _render(sim: Any, sample: Any, dets: Any, rmap: Any, n_om: int) -> List[Any]:
    log: List[Any] = []
    ctx: Dict[str, Any] = {"peak": None}
    sim.peak_hook = lambda v, r, oi, om: ctx.__setitem__("peak", (v, r, oi, om))
    sinks = [[_Sink(ctx, d, log) for d in range(len(dets))] for _ in range(n_om)]
    sim._simulate_peaks_batched(sinks, dets, sample, rmap)
    sim.peak_hook = None
    return log


def _setup(n_vox: int = 3) -> Tuple[Any, Any, Any, Any, int]:
    from icenine.config_file import ConfigFile
    from icenine.forward_simulation import ForwardSimulation
    from icenine.sample import Sample
    from icenine.simulation import Simulation

    cfg = ConfigFile.from_file(str(EX / "ConfigFiles" / "Example2.Simulation.config"))
    cfg.sample_filename = "SimInput/rand_500grains_1mm_neworient_s0.mic"
    sim = ForwardSimulation(cfg)
    sim.exp_setup.initialize_experiment()
    dets = sim.exp_setup.get_detector_list()
    rmap = sim.exp_setup.get_range_to_index_map()
    sim.simulator = Simulation(sim.exp_setup)
    sample = Sample()
    sim.exp_setup.initialize_sample(sample, dets[0])
    sample.get_mic().voxels = sample.get_mic().voxels[:n_vox]
    return sim, sample, dets, rmap, len(sim.exp_setup.get_omega_range_list())


def _same(a: List[Any], b: List[Any]) -> bool:
    return len(a) == len(b) and all(
        x[0] == y[0] and x[1] == y[1] and np.array_equal(x[2], y[2]) and x[3] == y[3]
        for x, y in zip(a, b)
    )


def test_default_off_ignores_deformation_and_flag_is_false() -> None:
    cwd = os.getcwd()
    os.chdir(EX)
    try:
        sim, sample, dets, rmap, n_om = _setup()
        assert sim.apply_deformation is False
        base = _render(sim, sample, dets, rmap, n_om)
        assert len(base) > 0
        eps = np.array([[1e-3, 0, 0], [0, 1e-3, 0], [0, 0, 1e-3]])
        for v in sample.get_mic().voxels:
            v.deformation = np.linalg.inv(np.eye(3) + eps).astype(np.float32)
        assert _same(base, _render(sim, sample, dets, rmap, n_om))  # off: deformation ignored
        sim.apply_deformation = True
        assert not _same(base, _render(sim, sample, dets, rmap, n_om))  # on: spots move
    finally:
        os.chdir(cwd)


def test_identity_deformation_matches_unstrained_render() -> None:
    cwd = os.getcwd()
    os.chdir(EX)
    try:
        sim, sample, dets, rmap, n_om = _setup()
        base = _render(sim, sample, dets, rmap, n_om)
        sim.apply_deformation = True  # voxels carry no deformation: identity
        on = _render(sim, sample, dets, rmap, n_om)
        assert len(on) == len(base)
        d = max(float(np.abs(x[2] - y[2]).max()) for x, y in zip(base, on))
        assert d < 5e-3  # pixels; float32 rounding of g = O (I g0)
    finally:
        os.chdir(cwd)


def test_strained_render_matches_spot_model_positions() -> None:
    import icenine.grain_fit as GF
    from icenine.reconstructor import _get_voxel_vertices
    from icenine.smooth_cost import SmoothCost

    cwd = os.getcwd()
    os.chdir(EX)
    try:
        sim, sample, dets, rmap, n_om = _setup(2)
        rng = np.random.default_rng(1)
        eps6 = rng.uniform(-1.5e-3, 1.5e-3, size=(2, 6))
        mic = sample.get_mic()
        for v, e in zip(mic.voxels, eps6):
            v.deformation = GF.strain_to_m(torch.as_tensor(e[None])).numpy()[0].astype(np.float32)
        sim.apply_deformation = True
        log = _render(sim, sample, dets, rmap, n_om)
        assert len(log) > 20

        sc = SmoothCost(
            sim.simulator, dets, rmap, None, sample, sample.get_structure_list(),
            eta_limit=sim.exp_setup.get_eta_limit(),
        )
        model = GF.SpotModel(sc, phase=int(mic.voxels[0].phase))
        Rs = torch.as_tensor(np.stack([v.orientation for v in mic.voxels]), dtype=torch.float64)
        Vs = torch.as_tensor(
            np.stack([_get_voxel_vertices(v).numpy() for v in mic.voxels]), dtype=torch.float64
        )
        rows = model.predict(Rs, Vs, torch.as_tensor(eps6, dtype=torch.float64))
        om_pred = (rows.fpos * sc.width + sc.low).numpy()
        worst, n_cmp = 0.0, 0
        for (v, r, oi, om), det, tri, _ in log:
            sel = np.nonzero(
                (rows.vox.numpy() == v) & (rows.refl.numpy() == r) & (np.abs(om_pred - om) < 1e-4)
            )[0]
            assert len(sel) == 1
            i = sel[0]
            pred = np.stack([rows.col[i, det, :3].numpy(), rows.row[i, det, :3].numpy()], axis=1)
            worst = max(worst, float(np.abs(pred - tri).max()))
            n_cmp += 1
        assert n_cmp == len(log)
        assert worst < 0.05, worst  # pixels
        # the strain moved the spots well beyond that tolerance
        sim.apply_deformation = False
        plain = _render(sim, sample, dets, rmap, n_om)
        assert not _same(plain, log)
    finally:
        os.chdir(cwd)


def test_mic_writer_default_unchanged_and_digits_option(tmp_path: Path) -> None:
    from icenine.mic_file import MicFile

    mic = MicFile.read(str(EX / "SimInput" / "rand_500grains_1mm_neworient_s0.mic"))
    mic.voxels = mic.voxels[:3]
    M = np.linalg.inv(np.eye(3) + np.diag([1e-4, -2e-4, 3e-4]) + 1e-4 * (np.ones((3, 3)) - np.eye(3)))
    for v in mic.voxels:
        v.deformation = M
    a, b, c = tmp_path / "a.mic", tmp_path / "b.mic", tmp_path / "c.mic"
    mic.write(str(a))
    mic.write(str(b), deformation_digits=6)
    mic.write(str(c), deformation_digits=9)
    assert a.read_bytes() == b.read_bytes()
    row6 = np.loadtxt(a, skiprows=1)[0, 13:19]
    row9 = np.loadtxt(c, skiprows=1)[0, 13:19]
    # the default 6-digit rounding is 1e-6 absolute (1% of a 1e-4 strain); 9 digits is 1e-9
    f32 = np.array([M[0, 0], M[1, 1], M[2, 2], M[1, 0], M[1, 2], M[2, 0]])
    assert np.abs(row9 - f32).max() < 2e-9  # written as float64 here
    assert np.abs(row6 - f32).max() > 2e-9
    back = MicFile.read(str(c)).voxels[0].deformation
    assert np.abs(back - M).max() < 6e-8  # read back as float32


def test_serial_path_strain_matches_batched() -> None:
    import icenine.grain_fit as GF

    cwd = os.getcwd()
    os.chdir(EX)
    try:
        sim, sample, dets, rmap, n_om = _setup(2)
        for v, e in zip(sample.get_mic().voxels, ([1e-3, -5e-4, 2e-4, 4e-4, 0, 3e-4], [0] * 6)):
            v.deformation = GF.strain_to_m(torch.tensor([e], dtype=torch.float64))[0].numpy()
        sim.apply_deformation = True
        log_b = _render(sim, sample, dets, rmap, n_om)
        log_s: List[Any] = []
        ctx: Dict[str, Any] = {"peak": None}
        sinks = [[_Sink(ctx, d, log_s) for d in range(len(dets))] for _ in range(n_om)]
        sim._simulate_peaks(sinks, dets, sample, rmap)

        def key(log: List[Any]) -> np.ndarray:
            a = np.array([np.concatenate([[t[1]], t[2].ravel()]) for t in log])
            return a[np.lexsort(a.T[::-1])]

        a, b = key(log_b), key(log_s)
        assert a.shape == b.shape and a.shape[0] > 10
        assert np.abs(a - b).max() < 0.05
    finally:
        os.chdir(cwd)
