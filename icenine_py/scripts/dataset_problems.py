"""Rebuild the per-voxel physics objects (observer, window spec) behind a saved toy dataset.

Shared by the Gauss-Newton baseline, the architecture diagnostics and the aux-table builder.
Datasets store windows, offsets and (for multi-voxel data) per-voxel context, but not the
observers, so they are rebuilt from the generator's own helpers with the dataset's settings.
"""

import sys
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).parent))
from generate_toy_orientation_dataset import (  # noqa: E402
    build_problem,
    example_dir_for,
    setup_example,
)


def load_meta(path):
    """Dataset dict with the windows memory-mapped (cheap to open)."""
    return torch.load(Path(path).resolve(), mmap=True)


def iter_voxels(data):
    """Yield dict(v, r_perp_um, n_pk, sample_idx, obs, spec, problem) for each voxel of the dataset.

    Changes the working directory (the physics setup does). Open the dataset file with an
    absolute path before calling this."""
    from icenine.orientation_eval import BatchedObserver, WindowSpec

    max_q = data.get("max_q", float("nan"))
    max_q = None if max_q != max_q else max_q
    example_dir = example_dir_for(data.get("example"))
    multi = bool(data.get("multi_voxel", False))
    if multi:
        vid = data["voxel_id"].numpy()
        voxels = [(v, int(data["voxel_indices"][v])) for v in range(len(data["voxel_indices"]))]
        setup = setup_example(example_dir, max_q=max_q)
    else:
        vid = None
        voxels = [(0, data["voxel_index"])]
        setup = None
    for v, mic_index in voxels:
        problem = build_problem(
            example_dir,
            mic_index,
            max_q=max_q,
            detectors=data.get("detectors", "first"),
            min_sin_eta=float(data.get("min_sin_eta", 0.0)),
            setup=setup,
        )
        n_pk = int(data["n_peaks_per_voxel"][v]) if multi else int(data["n_peaks"])
        assert len(problem["roi_list"]) == n_pk, "ROI set differs from the dataset's"
        obs = BatchedObserver(
            problem["R_nom"],
            problem["vertices"],
            problem["sample"],
            problem["detector_list"],
            problem["range_map"],
            problem["exp_setup"],
            problem["roi_list"],
        )
        spec = WindowSpec.from_nominal(obs, data["window_size"], data["frame_half_width"])
        pos = np.asarray(problem["voxel"].position, dtype=float)
        r_perp = float(np.hypot(pos[0], pos[1]) * 1e3)
        sample_idx = np.nonzero(vid == v)[0] if multi else np.arange(len(data["offsets_deg"]))
        yield dict(
            v=v,
            r_perp_um=r_perp,
            n_pk=n_pk,
            sample_idx=sample_idx,
            obs=obs,
            spec=spec,
            problem=problem,
        )
