"""Per-evaluation timing of VoxelCostFunction.evaluate on Example2.ThreeVoxels.

Usage (from icenine_py/): uv run python scripts/speed_parity/prof_eval.py [--lp "<arm name>"]
(--lp needs: uv run --with line_profiler ...)
"""

import os
import sys
import time
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "icenine_py"))
torch.set_num_threads(1)

from icenine.config_file import ConfigFile  # noqa: E402
from icenine.experiment_setup import XDMExperimentSetup  # noqa: E402
from icenine.experimental_data import ExperimentalData  # noqa: E402
from icenine.mic_file import MicFile  # noqa: E402
from icenine.reconstructor import _get_voxel_vertices  # noqa: E402
from icenine.sample import Sample  # noqa: E402
from icenine.simulation import Simulation  # noqa: E402
from icenine.cost_functions import VoxelCostFunction  # noqa: E402
import icenine.cost_functions as cf  # noqa: E402

ex = ROOT / "Examples" / "Example2.ThreeVoxels"
os.chdir(ex)
config = ConfigFile.from_file(str(ex / "ConfigFiles" / "Example2.Simulation.config"))
exp_data = ExperimentalData.from_image_directory(
    directory=str(ex / "ScatteringData_Python"),
    basename="3Grains.sim",
    ext="d",
    serial_length=5,
    n_omega=180,
    n_detectors=2,
    num_rows=2048,
    num_cols=2048,
)
if hasattr(exp_data, "prepare_for_reconstruction"):
    exp_data.prepare_for_reconstruction()
setup = XDMExperimentSetup(config)
setup.initialize_experiment()
dets = setup.get_detector_list()
rmap = setup.get_range_to_index_map()
sample = Sample()
setup.initialize_sample(sample, dets[0])
sim = Simulation(setup)
structs = sample.get_structure_list()
mic = MicFile.read(str(ex / "SimInput" / "three_voxels.mic"))
vox = mic.voxels[0]
verts = _get_voxel_vertices(vox)


def make(pixel_radius, max_q, cpp_rule):
    return VoxelCostFunction(
        simulator=sim,
        detector_list=dets,
        range_map=rmap,
        exp_data=exp_data,
        sample=sample,
        structure_list=structs,
        mode="hard",
        pixel_radius=pixel_radius,
        max_q=max_q,
        cpp_overlap_rule=cpp_rule,
    )


def rot(axis, ang):
    axis = axis / np.linalg.norm(axis)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * K @ K


rng = np.random.default_rng(0)
truth = np.asarray(vox.orientation, dtype=np.float64)
orients = [rot(rng.normal(size=3), np.deg2rad(rng.uniform(0, 3))) @ truth for _ in range(400)]

ARMS = {
    "screen (r=3, Qmax 5, cpp rule)": (3, 5.0, True),
    "recost/quickMC (r=3, Qmax 5..8)": (3, 8.0, True),
    "final (triangle, Qmax 8, cpp rule)": (0, 8.0, True),
}


def bench(fn, reps=3):
    best = None
    for _ in range(reps):
        t0 = time.perf_counter()
        for o in orients:
            fn.evaluate(o, verts, vox.phase)
        dt = (time.perf_counter() - t0) / len(orients) * 1e6
        best = dt if best is None else min(best, dt)
    return best


def peak_counts(fn):
    """K reflections, N=2K omega solutions, M valid peaks (median over orients)."""
    g, m = fn._phase_recip_vecs[vox.phase]
    ks, ns = [], []
    from icenine.diffraction_core import get_scattering_omegas_torch

    for o in orients[:50]:
        gl = (torch.from_numpy(o).float() @ g.T).T
        r = get_scattering_omegas_torch(gl, m, sim.beam_energy, sim.beam_deflection_chi)
        ks.append(len(g))
        ns.append(int(r.observable.sum()) * 2)
    return int(np.median(ks)), int(np.median(ns))


def count_torch_ops(fn):
    from torch.profiler import profile, ProfilerActivity

    with profile(activities=[ProfilerActivity.CPU]) as p:
        for o in orients[:20]:
            fn.evaluate(o, verts, vox.phase)
    ev = p.key_averages()
    n_ops = sum(e.count for e in ev if e.key.startswith("aten::")) / 20
    return n_ops


if __name__ == "__main__":
    print(
        f"C ext: {cf._HAS_C_RASTERIZE}, torch {torch.__version__},"
        f" threads {torch.get_num_threads()}"
    )
    for name, args in ARMS.items():
        fn = make(*args)
        for o in orients[:20]:
            fn.evaluate(o, verts, vox.phase)  # warm caches
        us = bench(fn)
        K, N = peak_counts(fn)
        ops = count_torch_ops(fn)
        print(
            f"{name:40s} {us:7.1f} us/eval   reflections={K} omega-solutions={N}"
            f"  aten ops/eval={ops:.0f}"
        )

    if "--lp" in sys.argv:
        from line_profiler import LineProfiler
        from icenine import diffraction_core

        fn = make(
            *(
                ARMS[sys.argv[sys.argv.index("--lp") + 1]]
                if len(sys.argv) > sys.argv.index("--lp") + 1
                else ARMS["final (triangle, Qmax 8, cpp rule)"]
            )
        )
        lp = LineProfiler(
            VoxelCostFunction.evaluate,
            cf.calculate_diffraction_overlap_batched,
            diffraction_core.get_scattering_omegas_torch,
        )
        lp.runcall(lambda: [fn.evaluate(o, verts, vox.phase) for o in orients])
        lp.print_stats(output_unit=1e-6, stripzeros=True)
