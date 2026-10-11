"""Opt-in C++ CountQualifiedPeaks rule (CppOverlapRule): hand-built detector patterns, Python and
C paths identical, switch off bit-identical to the original rule."""

import itertools

import numpy as np
import pytest

from icenine import cost_functions as CF
from icenine.cost_functions import count_qualified_peaks, count_qualified_peaks_cpp


def cpp_ref(lit, ov):
    """Independent transcription of Src/OverlapInfo.tmpl.cpp:78-139."""
    n = len(lit)
    valid_sim, on_prev, n_ov = lit[0], lit[0], 0
    for i in range(1, n):
        if not valid_sim:
            break
        if lit[i]:
            valid_sim = True
        elif on_prev:
            valid_sim, on_prev = True, False
        else:
            valid_sim = False
    prev_ov, valid_ov = ov[0], ov[0]
    if ov[0]:
        n_ov += 1
    for i in range(1, n):
        if not valid_ov:
            break
        if ov[i]:
            n_ov += 1
            valid_ov = True
        elif prev_ov and not lit[i]:
            valid_ov, prev_ov = True, False
        else:
            valid_ov = False
    return int(valid_sim), int(valid_ov), n_ov if (valid_sim and valid_ov) else 0


def test_two_detector_patterns():
    # (lit, overlap) -> (on, overlap, n_det)
    cases = {
        ((True, True), (True, True)): (1, 1, 2),
        ((True, True), (True, False)): (1, 0, 0),  # the rule that differs from the original
        ((True, False), (True, False)): (1, 1, 1),
        ((True, True), (False, True)): (1, 0, 0),
        ((False, True), (False, True)): (0, 0, 0),
    }
    for (lit, ov), want in cases.items():
        assert count_qualified_peaks(list(lit), list(ov), 2, cpp_rule=True) == want
    # the original rule counts the second case as overlap
    assert count_qualified_peaks([True, True], [True, False], 2) == (1, 1, 1)


def test_three_detector_quirks():
    # lit, unlit, lit is a valid peak in C++ (bOnPreviousDetector), invalid in the original rule
    assert count_qualified_peaks([True, False, True], [False] * 3, 3, cpp_rule=True)[0] == 1
    assert count_qualified_peaks([True, False, True], [False] * 3, 3)[0] == 0
    # lit, unlit, unlit is invalid in C++, valid in the original rule
    assert count_qualified_peaks([True, False, False], [False] * 3, 3, cpp_rule=True)[0] == 0
    assert count_qualified_peaks([True, False, False], [False] * 3, 3)[0] == 1
    # bPreviousOverlap: one unlit non-overlapping detector tolerated, then a second invalidates
    assert count_qualified_peaks([True, False, True], [True, False, True], 3, cpp_rule=True)[1] == 1
    assert (
        count_qualified_peaks([True, False, False], [True, False, False], 3, cpp_rule=True)[1] == 0
    )


@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_exhaustive_matches_reference(n):
    for lit in itertools.product([False, True], repeat=n):
        for ov in itertools.product([False, True], repeat=n):
            ov = tuple(o and l for o, l in zip(ov, lit))  # overlap needs a lit detector
            assert count_qualified_peaks_cpp(list(lit), list(ov), n) == cpp_ref(lit, ov)


def test_default_rule_unchanged():
    """cpp_rule=False reproduces the original rule on every 2- and 3-detector pattern."""
    for n in (2, 3):
        for lit in itertools.product([False, True], repeat=n):
            for ov in itertools.product([False, True], repeat=n):
                assert count_qualified_peaks(list(lit), list(ov), n) == count_qualified_peaks(
                    list(lit), list(ov), n, cpp_rule=False
                )


@pytest.mark.skipif(not CF._HAS_C_RASTERIZE or not CF._C_HAS_CPP_RULE, reason="C ext missing")
def test_c_path_equals_python_path():
    """stage_d_overlap with the flag equals the Python fallback on random two-detector input."""
    rng = np.random.default_rng(0)
    rows = cols = 64
    n_det, n_peaks, rad = 2, 200, 2
    images = [(rng.random((rows, cols)) < 0.25).astype(np.uint8) for _ in range(n_det * 3)]
    wedge = rng.integers(0, 3, n_peaks).astype(np.int32)
    hit = (rng.random((n_peaks, n_det)) < 0.85).astype(np.uint8)
    centers = rng.integers(-3, 67, (n_peaks, n_det, 2)).astype(np.int32)
    res = CF._c_stage_d_overlap(
        images, wedge, hit, np.ascontiguousarray(centers), None, n_det, rad, n_peaks, rows, cols, 1
    )
    po = pon = 0
    for p in range(n_peaks):
        lit, ov = [False] * n_det, [False] * n_det
        for d in range(n_det):
            if not hit[p, d]:
                continue
            c, r = centers[p, d]
            inb = bright = False
            for dx in range(-rad, rad + 1):
                for dy in range(-rad, rad + 1):
                    x, y = c + dx, r + dy
                    if 0 <= x < cols and 0 <= y < rows:
                        inb = True
                        bright = bright or bool(images[wedge[p] * n_det + d][y, x])
            lit[d], ov[d] = inb, bright
        on, o, _ = count_qualified_peaks(lit, ov, n_det, cpp_rule=True)
        pon += on
        po += o
    assert (res[2], res[3]) == (po, pon)
    off = CF._c_stage_d_overlap(
        images, wedge, hit, np.ascontiguousarray(centers), None, n_det, rad, n_peaks, rows, cols
    )
    assert off[2] >= res[2]  # the original rule never counts fewer overlapping peaks


def test_config_key(tmp_path):
    from pathlib import Path

    from icenine.config_file import ConfigFile

    base = Path(__file__).resolve().parents[1] / "scripts" / "phase_d" / "configs"
    text = (base / "ReconstructPhaseD_mc_clean.config").read_text()
    off = tmp_path / "off.config"
    on = tmp_path / "on.config"
    off.write_text(text)
    on.write_text(text + "\nCppOverlapRule 1\n")
    assert ConfigFile.from_file(str(off)).cpp_overlap_rule is False
    assert ConfigFile.from_file(str(on)).cpp_overlap_rule is True


def test_search_parameters_carry_the_key(tmp_path):
    from pathlib import Path

    from icenine.config_file import ConfigFile
    from icenine.orientation_search import SearchParameters

    base = Path(__file__).resolve().parents[1] / "scripts" / "phase_d" / "configs"
    text = (base / "ReconstructPhaseD_mc_clean.config").read_text()
    on = tmp_path / "on.config"
    on.write_text(text + "\nCppOverlapRule 1\n")
    assert SearchParameters.from_config(ConfigFile.from_file(str(on))).cpp_overlap_rule is True
    off = tmp_path / "off.config"
    off.write_text(text)
    assert SearchParameters().cpp_overlap_rule is False
    assert SearchParameters.from_config(ConfigFile.from_file(str(off))).cpp_overlap_rule is False
