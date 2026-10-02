"""Small helpers of the toy orientation scripts: prediction file naming, padding mask,
results-json layouts."""

import sys
from pathlib import Path

import pytest
import torch

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))
import summarize_results as sr  # noqa: E402
import train_toy_orientation_nn as tr  # noqa: E402


@pytest.mark.parametrize(
    "save, variant, expected",
    [
        ("out/p.npz", "clean", "out/p.npz"),
        ("out/p.npz", "all", "out/p_all.npz"),
        ("out/p", "all", "out/p_all.npz"),  # no suffix: variants must not collide with clean
        ("out/p.v2.npz", "noise", "out/p.v2_noise.npz"),
    ],
)
def test_variant_path(save, variant, expected):
    assert tr.variant_path(save, variant) == expected


def test_variant_paths_are_distinct_without_suffix():
    paths = {tr.variant_path("p", v) for v in ("clean", "neighbours", "noise", "all")}
    assert len(paths) == 4


def test_padding_mask():
    data = {
        "n_peaks": 4,
        "n_peaks_per_voxel": torch.tensor([2, 4]),
        "voxel_id": torch.tensor([0, 1, 0]),
    }
    m = tr.padding_mask(data)
    assert m.tolist() == [[1, 1, 0, 0], [1, 1, 1, 1], [1, 1, 0, 0]]
    assert tr.padding_mask({"n_peaks": 4}) is None


def test_summarize_results_reads_both_layouts():
    rows = {"in-dist/0.5": {}, "summary": {}}
    assert sr.rows_of(rows, "clean") is rows  # older single-variant layout
    wrapped = {"clean": rows, "all": {"x": 1}}
    assert sr.rows_of(wrapped, "clean") is rows
    assert sr.rows_of({"clean": rows}, "clean") is rows
    r = {"in-dist/0.25": {}, "in-dist/0.10000000149011612": {}, "in-dist/1.0": {}}
    assert sr.mag_keys(r, "in-dist") == ["0.10000000149011612", "0.25", "1.0"]


def test_gn_layer_net_warns_without_nominal_offsets():
    from icenine.toy_orientation_model import GNLayerNet

    net = GNLayerNet(window_size=8, in_channels=2, context_dim=9, frame_width_rad=0.01)
    x = torch.zeros(1, 3, 2, 8, 8)
    with pytest.warns(UserWarning, match="sub-pixel bias"):
        net.physics(x, torch.zeros(3, 9), None)
