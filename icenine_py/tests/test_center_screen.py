"""Pin the default start point of the pixel-radius screen (the voxel centroid, C++ VoxelToCenter).

Commit d96c413 changed the pixel_radius > 0 screen from projecting triangle vertex 0 to projecting
the voxel centroid (C++ CostFunctions::VoxelToCenter, used by GlobalSearchCostFunctions in
DiscreteAdaptive.h and Reconstructor.h). No numeric golden moved with it, so these tests pin it:
on a constructed triangle whose vertex 0 and centroid project to different pixels, the default cost
equals the cost of the triangle collapsed to its centroid and differs from the vertex-0 opt-out.
"""

import inspect
import os
from pathlib import Path
from typing import Any, Dict

import pytest
import torch

from icenine.cost_functions import OverlapInfo, VoxelCostFunction

EXAMPLE_DIR = Path(__file__).parent.parent.parent / "Examples" / "Example2.ThreeVoxels"
_CACHE: Dict[str, Any] = {}


def test_center_screen_default_is_true() -> None:
    """The constructor default is the centroid screen (identical start point to C++)."""
    default = inspect.signature(VoxelCostFunction.__init__).parameters["center_screen"].default
    assert default is True


@pytest.fixture(scope="module")
def three_voxels() -> Dict[str, Any]:
    """ThreeVoxels forward-simulated images (sparse) and the setup to build cost functions."""
    if _CACHE:
        return _CACHE
    from icenine.config_file import ConfigFile
    from icenine.experiment_setup import XDMExperimentSetup
    from icenine.experimental_data import ExperimentalData
    from icenine.mic_file import MicFile
    from icenine.reconstructor import _get_voxel_vertices
    from icenine.sample import Sample
    from icenine.simulation import Simulation

    config_path = EXAMPLE_DIR / "ConfigFiles" / "Example2.Simulation.config"
    data_dir = EXAMPLE_DIR / "ScatteringData_Python"
    mic_path = EXAMPLE_DIR / "SimInput" / "three_voxels.mic"
    if not config_path.exists() or not mic_path.exists():
        pytest.skip("ThreeVoxels config or .mic not found")
    if not data_dir.exists() or len(list(data_dir.glob("*.d*"))) < 360:
        pytest.skip(f"Forward sim output not found or incomplete: {data_dir}")

    old_cwd = os.getcwd()
    os.chdir(EXAMPLE_DIR)
    try:
        config = ConfigFile.from_file(str(config_path))
        exp_data = ExperimentalData.from_image_directory(
            directory=str(data_dir),
            basename="3Grains.sim",
            ext="d",
            serial_length=5,
            n_omega=180,
            n_detectors=2,
            num_rows=2048,
            num_cols=2048,
            mode="sparse",
        )
        setup = XDMExperimentSetup(config)
        setup.initialize_experiment()
        detector_list = setup.get_detector_list()
        sample = Sample()
        setup.initialize_sample(sample, detector_list[0])
        mic = MicFile.read(str(mic_path))
        _CACHE.update(
            kwargs=dict(
                simulator=Simulation(setup),
                detector_list=detector_list,
                range_map=setup.get_range_to_index_map(),
                exp_data=exp_data,
                sample=sample,
                structure_list=sample.get_structure_list(),
            ),
            voxel=mic.voxels[0],
            vertices=_get_voxel_vertices(mic.voxels[0]),
        )
    finally:
        os.chdir(old_cwd)
    return _CACHE


def _collapse(tri: torch.Tensor) -> torch.Tensor:
    return tri.mean(dim=0, keepdim=True).repeat(3, 1)


def test_center_screen_projects_the_centroid(three_voxels: Dict[str, Any]) -> None:
    """Default (centroid) screen == vertex-0 screen on the collapsed triangle; != vertex 0."""
    c = three_voxels
    voxel = c["voxel"]
    # Constructed triangle: vertex 0 at the true voxel centre (bright at the true orientation),
    # the other two vertices 0.06 mm away, so the centroid sits about 0.028 mm off (tens of
    # pixels on the detector) and lands off the voxel's spots.
    centre = c["vertices"].mean(dim=0)
    d = 0.06
    tri = torch.stack(
        [centre, centre + torch.tensor([d, 0.0, 0.0]), centre + torch.tensor([0.0, d, 0.0])]
    )
    orient = voxel.orientation

    on = VoxelCostFunction(pixel_radius=3, **c["kwargs"])
    off = VoxelCostFunction(pixel_radius=3, center_screen=False, **c["kwargs"])
    assert on.center_screen is True

    r_on = on.evaluate(orient, tri, phase_index=voxel.phase)
    r_collapsed = off.evaluate(orient, _collapse(tri), phase_index=voxel.phase)
    r_vertex0 = off.evaluate(orient, tri, phase_index=voxel.phase)
    assert isinstance(r_on, OverlapInfo)
    assert r_on.peak_on_detector > 0
    assert r_on == r_collapsed
    assert r_vertex0 != r_on
    # vertex 0 is the true voxel centre, the centroid is not: the opt-out scores higher here
    assert r_vertex0.quality > r_on.quality


def test_center_screen_has_no_effect_without_pixel_radius(three_voxels: Dict[str, Any]) -> None:
    """pixel_radius 0 rasterises the full triangle, so center_screen does not change the cost."""
    c = three_voxels
    voxel = c["voxel"]
    tri = c["vertices"]
    on = VoxelCostFunction(pixel_radius=0, **c["kwargs"])
    off = VoxelCostFunction(pixel_radius=0, center_screen=False, **c["kwargs"])
    r_on = on.evaluate(voxel.orientation, tri, phase_index=voxel.phase)
    assert r_on.peak_on_detector > 0
    assert r_on == off.evaluate(voxel.orientation, tri, phase_index=voxel.phase)
