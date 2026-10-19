"""Task E review fix 11: record the key-off fixture of ``tests/test_smooth_cost.py`` on a commit.

One no-start ``reconstruct_voxel`` on a ThreeVoxels toy problem (voxel 0, its real two-detector
images) with a candidate list long enough for the ranking to matter: 16 start orientations 0.5 to 4
deg from the truth plus 16 random ones (fixed seeds), two search levels (``max_local_resolution``
1), rng ``default_rng(7)``, one torch thread. Level 0 keeps about 20 candidates after the quick MC,
so a changed ranking changes the recorded order (the three-orientation toy of
``test_findoptimal_refactor`` keeps one, which no ranking can change).

Saved: the final orientation and cost, every recorder event's numeric arrays (in event order) and
``last_level_best``. The fixture is made on the BASE commit (before Task E), with the library of
that commit first on the path, so that the test on this branch compares the key-off search with the
search before the key existed:

  git archive <base> icenine_py/icenine | tar -x -C <scratch>
  cp icenine/_rasterize.cpython-*.so <scratch>/icenine_py/icenine/   # extension unchanged by E
  PYTHONPATH=<scratch>/icenine_py uv run python scripts/eval_reduction/smooth_keyoff_fixture.py \
      --out tests/data/smooth_keyoff_ref.npz --commit <base>

The script prints which ``icenine`` it imported; check that it is the scratch copy.
"""

import argparse
import contextlib
import io
import json
import os
import platform
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import torch

ICE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ICE / "tests"))
N_NEAR = 16
N_RANDOM = 16


def build_toy() -> Tuple[Any, Any]:
    """(reconstructor, voxel) of the toy problem (the library on the path, any commit)."""
    from scipy.spatial.transform import Rotation

    import test_findoptimal_refactor as T
    from icenine.config_file import ConfigFile
    from icenine.experimental_data import ExperimentalData
    from icenine.mic_file import MicFile
    from icenine.reconstructor import AdaptiveVoxelReconstructor, setup_reconstruction

    data_dir = T._require_data()
    cwd = os.getcwd()
    os.chdir(T.EXAMPLE)
    try:
        config = ConfigFile.from_file("ConfigFiles/ReconstructQ8.config")
        config.out_file_basename = "3Grains.sim"
        exp_data = ExperimentalData.from_image_directory(
            directory=str(data_dir),
            basename="3Grains.sim",
            ext="d",
            serial_length=5,
            n_omega=180,
            n_detectors=2,
            num_rows=2048,
            num_cols=2048,
        )
        voxel = MicFile.read("SimInput/three_voxels.mic").voxels[0]
        R_true = np.asarray(voxel.orientation, dtype=np.float64)
        rng = np.random.default_rng(3)
        rv = rng.normal(size=(N_NEAR, 3)) / np.sqrt(3)
        rv *= np.radians(rng.uniform(0.5, 4.0, size=(N_NEAR, 1)))
        near = Rotation.from_rotvec(rv).as_matrix() @ R_true
        fz = np.concatenate([near, Rotation.random(N_RANDOM, random_state=5).as_matrix()])
        setup = setup_reconstruction(config, exp_data=exp_data, fz_orientations=fz)
        setup.search_params.max_local_resolution = 1
        setup.search_params.max_mc_steps = 150
        rec = AdaptiveVoxelReconstructor(setup, min_sin_eta=0.0)
    finally:
        os.chdir(cwd)
    return rec, voxel


def toy_run(smooth_rank: Optional[str] = None) -> Tuple[Any, List[Tuple[str, Any]], Any]:
    """(result, recorder events, reconstructor) of one toy search; ``smooth_rank`` None leaves
    the attribute untouched (the base commit has no such key)."""
    from icenine.reconstructor import _get_voxel_vertices

    rec, voxel = build_toy()
    if smooth_rank is not None:
        rec.setup.search_params.smooth_rank = smooth_rank
    events: List[Tuple[str, Dict[str, Any]]] = []
    rec.recorder = lambda name, data: events.append((name, data))
    with contextlib.redirect_stdout(io.StringIO()):
        res = rec.reconstruct_voxel(
            _get_voxel_vertices(voxel), voxel.phase, rng=np.random.default_rng(7)
        )
    return res, events, rec


def arrays_of(
    res: Any, events: List[Tuple[str, Dict[str, Any]]], level_best: List[Tuple[int, Any, float]]
) -> Dict[str, np.ndarray]:
    out: Dict[str, np.ndarray] = {
        "final/orientation": np.asarray(res.orientation, dtype=np.float64),
        "final/cost": np.asarray(np.float64(res.cost)),
        "events/names": np.array([n for n, _ in events]),
    }
    for i, (name, data) in enumerate(events):
        for k, v in data.items():
            a = np.asarray(v)
            if a.dtype.kind in "biuf":
                out[f"events/{i:03d}_{name}/{k}"] = a
    for j, (lv, R, c) in enumerate(level_best):
        out[f"level_best/{j}/level"] = np.asarray(np.int64(lv))
        out[f"level_best/{j}/orientation"] = np.asarray(R, dtype=np.float64)
        out[f"level_best/{j}/cost"] = np.asarray(np.float64(c))
    return out


def main() -> None:
    import icenine

    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--commit", required=True, help="the commit whose library is on the path")
    a = ap.parse_args()
    torch.set_num_threads(1)
    print("icenine imported from", Path(icenine.__file__).parent, flush=True)
    res, events, rec = toy_run()
    arrays = arrays_of(res, events, rec.last_level_best)
    meta = dict(
        commit=a.commit,
        platform_system=platform.system(),
        platform_machine=platform.machine(),
        torch=torch.__version__,
        threads=torch.get_num_threads(),
        toy=f"ThreeVoxels voxel 0, {N_NEAR} near + {N_RANDOM} random starts, 2 levels, rng 7",
    )
    arrays["meta"] = np.array(json.dumps(meta))
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(a.out, **arrays)  # type: ignore[arg-type]
    n_l0 = [len(d["cost"]) for n, d in events if n == "quick_mc" and d["level"] == 0]
    print(f"{a.out}: {len(arrays)} arrays, level-0 list {n_l0}, meta {meta}")


if __name__ == "__main__":
    main()
