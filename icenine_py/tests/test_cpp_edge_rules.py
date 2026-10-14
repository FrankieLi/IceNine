"""Opt-in C++ low-edge conventions (CppEdgeRules): frame lookup below the range, negative detector
coordinate, pixel-less triangle lit through an in-bounds vertex; the switch off is bit-identical."""

from pathlib import Path

import numpy as np
import pytest
import torch

from icenine import cost_functions as CF

CONFIGS = Path(__file__).resolve().parents[1] / "scripts" / "phase_d" / "configs"


def test_omega_bins_default_truncates_toward_zero():
    low, width = -np.pi / 2, np.pi / 180
    om = np.array([low - 0.5 * width, low - 1.5 * width, low + 0.5 * width, low + 3.5 * width])
    assert CF._omega_bins(om, low, width, False).tolist() == [0, -1, 0, 3]
    # C++ AngleToIndex: f < 0 -> -1 (no frame); otherwise truncation
    assert CF._omega_bins(om, low, width, True).tolist() == [-1, -2, 0, 3]


def test_pixel_coordinate_negative_is_minus_one_only_with_the_switch():
    d = torch.tensor([-0.2, -0.0001, 0.0, 0.3])
    half, size = 0.1, 0.2
    default = CF._pixel_coordinate(d, half, size, False)
    assert torch.equal(default, (d + half) / size)
    edge = CF._pixel_coordinate(d, half, size, True)
    assert edge[0] == -1.0 and edge[1] == -1.0
    assert torch.equal(edge[2:], default[2:])


@pytest.mark.skipif(
    not CF._HAS_C_RASTERIZE or not CF._C_HAS_VERTEX_LIT,
    reason="C ext missing or built before cpp_edge_rules (rebuild _rasterize)",
)
def test_pixel_less_triangle_is_lit_by_an_in_bounds_vertex():
    rows = cols = 32
    img = [np.zeros((rows, cols), dtype=np.uint8)]
    wedge = np.zeros(1, dtype=np.int32)
    hit = np.ones((1, 1), dtype=np.uint8)

    def run(verts, bits):
        v = np.ascontiguousarray(np.array(verts, dtype=np.float32).reshape(1, 1, 3, 2))
        return CF._c_stage_d_overlap(img, wedge, hit, None, v, 1, 0, 1, rows, cols, bits)

    # a triangle on the last image row produces no pixel (the fill stops one row early, in C++ too)
    flat = [[10.0, 31.0], [12.0, 31.0], [14.0, 31.0]]
    base = run(flat, 1)
    assert base[1] == 0 and base[3] == 0  # the Python rule needs a pixel: not on the detector
    assert run(flat, 3)[3] == 1  # C++ lights the detector through the in-bounds vertex
    out = [[-50.0, 10.0], [-50.0, 10.0], [-50.0, 10.0]]  # no vertex in bounds: still unlit
    assert run(out, 3)[3] == 0
    tri = [[5.0, 5.0], [12.0, 5.0], [5.0, 12.0]]  # a real triangle: both rules agree
    assert run(tri, 1) == run(tri, 3)


def test_default_flags_unchanged_and_key_reaches_the_flags(tmp_path):
    from icenine.config_file import ConfigFile
    from icenine.orientation_search import SearchParameters
    from icenine.ray_overlap import RayRule
    from icenine.reconstructor import AdaptiveVoxelReconstructor

    text = (CONFIGS / "ReconstructPhaseD_mc_clean.config").read_text()
    off, on = tmp_path / "off.config", tmp_path / "on.config"
    off.write_text(text)
    on.write_text(text + "\nCppEdgeRules 1\n")
    assert ConfigFile.from_file(str(off)).cpp_edge_rules is False
    assert ConfigFile.from_file(str(on)).cpp_edge_rules is True
    assert SearchParameters().cpp_edge_rules is False
    assert SearchParameters.from_config(ConfigFile.from_file(str(on))).cpp_edge_rules is True

    class _Setup:
        search_params = SearchParameters()

    rec = AdaptiveVoxelReconstructor.__new__(AdaptiveVoxelReconstructor)
    rec.setup = _Setup()  # type: ignore[assignment]
    assert rec._cost_flags() == dict(center_screen=True)
    assert rec._cost_flags(RayRule()) == dict(cpp_overlap_rule=True, center_screen=True)
    rec.setup.search_params.cpp_edge_rules = True  # type: ignore[attr-defined]
    assert rec._cost_flags() == dict(center_screen=True, cpp_edge_rules=True)
    assert rec._cost_flags(RayRule()) == dict(
        cpp_overlap_rule=True, center_screen=True, cpp_edge_rules=True
    )
