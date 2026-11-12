"""Tests of the pure parts of the cpptime scripts (scripts/autoresearch/cpptime/).

C++ time on the pilot: one client. The tests pin the config conversion (make_config), the stop
rules of the guard, the log parsing (parse_logs), the voxel matching and scoring (score_cpp), the
identity check of two C++ runs (same_cpp_run) and the quiet-limit choice.
"""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ICE = Path(__file__).resolve().parents[1]
SCRIPTS = ICE / "scripts" / "autoresearch" / "cpptime"
sys.path.insert(0, str(SCRIPTS))

import guard as G  # noqa: E402
import make_config as MC  # noqa: E402
import parse_logs as PL  # noqa: E402
import same_cpp_run as SC  # noqa: E402
import score_cpp as SCORE  # noqa: E402

PY_CONFIG = """InfileBasename \t ScatteringData_PhaseD/full/clean/500Grains.sim
OutStructureBasename \t PhaseD_mc_clean
DetectorFilename  ConfigFiles/StandardGeometry.2Det
SampleFilename           SimInput/x_grid.mic
StructureFilename        DataFiles/copper.dat
RotationRangeFilename    DataFiles/omega_180_2L.dat
FundamentalZoneFilename  DataFiles/MyFZ.dat
SecondsBetweenSave            1800
OptimizationFilename            ConfigFiles/Optimization.2Det
DetectionLimitFilename          ConfigFiles/Sensitivity.2Det
OptimizationConstrainFilename ConfigFiles/Constrain.2L
MaxMCSteps                    200
LazyBFS
BFSRevisitRefit 1
BFSReseed 1
FusedCost 1
"""


def test_convert_drops_python_keys_and_absolutizes() -> None:
    out = MC.convert(PY_CONFIG, Path("/r/region_grid.mic"), "cpp_pilot", ex=Path("/ex"))
    keys = [ln.split()[0] for ln in out.splitlines() if ln.split()]
    assert not {"BFSRevisitRefit", "BFSReseed", "FusedCost"} & set(keys)
    assert "SampleFilename /r/region_grid.mic" in out
    assert "OutStructureBasename cpp_pilot" in out
    assert "SecondsBetweenSave 100000" in out
    assert "DetectorFilename /ex/ConfigFiles/StandardGeometry.2Det" in out
    assert "InfileBasename /ex/ScatteringData_PhaseD/full/clean/500Grains.sim" in out
    # the search lines are untouched
    assert "MaxMCSteps                    200" in out and "LazyBFS" in keys


def test_convert_needs_the_path_keys() -> None:
    with pytest.raises(AssertionError):
        MC.convert("LazyBFS\n", Path("/r/g.mic"), "x", ex=Path("/ex"))


def test_real_arm_a_config_converts() -> None:
    out = MC.convert(MC.SRC.read_text(), MC.REGION_GRID, "cpp_pilot")
    keys = [ln.split()[0] for ln in out.splitlines() if ln.split() and ln[0] != "#"]
    assert not {"BFSRevisitRefit", "BFSReseed", "FusedCost"} & set(keys)
    assert keys.count("SampleFilename") == 1


def test_guard_passes_the_exit_code(tmp_path: Path) -> None:
    rc = G.run([sys.executable, "-c", "import sys; sys.exit(3)"], 30, 30, tmp_path / "e", [])
    assert rc == 3


def test_guard_stamps_stderr_lines(tmp_path: Path) -> None:
    cmd = [sys.executable, "-c", "import sys; print('hello', file=sys.stderr)"]
    assert G.run(cmd, 30, 30, tmp_path / "e", []) == 0
    stamp, text = (tmp_path / "e").read_text().split(" ", 1)
    assert text.strip() == "hello" and float(stamp) > 1e9


def test_guard_kills_on_runtime_assert(tmp_path: Path) -> None:
    code = "import sys,time; print('RUNTIME_ASSERT x', file=sys.stderr, flush=True); time.sleep(60)"
    assert G.run([sys.executable, "-c", code], 30, 30, tmp_path / "e", [], poll_s=0.1) == G.FAILED


def test_guard_flags_a_failure_marker_after_exit(tmp_path: Path) -> None:
    code = "import sys; print('[ConfigFile] Error: x', file=sys.stderr)"
    assert G.run([sys.executable, "-c", code], 30, 30, tmp_path / "e", [], poll_s=0.1) == G.FAILED


def test_guard_time_limit(tmp_path: Path) -> None:
    rc = G.run(
        [sys.executable, "-c", "import time; time.sleep(60)"], 1, 30, tmp_path / "e", [], 0.1
    )
    assert rc == G.LIMIT


def test_guard_stall(tmp_path: Path) -> None:
    log = tmp_path / "log.txt"
    log.write_text("a")
    rc = G.run(
        [sys.executable, "-c", "import time; time.sleep(60)"],
        60,
        1,
        tmp_path / "e",
        [str(log)],
        poll_s=0.1,
    )
    assert rc == G.STALL


CLIENT = """StartTime: Thu Nov 12 00:17:37 2026
Beginning Lazy BFS: Triangle
  First Voxel Required: 15 Sec
    Orientation: 91.0 8.8 228.7
 [Conf, HitRatio ] 0.985294 1
Fitting neighbors via Breadth First Search
QueueSize 5 Fitted 2 continuous voxels
  [ Cost ] 0.1  1/2 3 4/5 0.9 |=|
QueueSize 4 Fitted 3 continuous voxels
  First Voxel Required: 16 Sec
 [Conf, HitRatio ] 0.5 0.5
Total Processing Time: (sec) 40... Client Exiting
Stop: Thu Nov 12 00:37:37 2026
Total Time elapsed (sec): 77
"""


def test_parse_client_log() -> None:
    r = PL.parse_client_log(CLIENT)
    assert r["global_searches"] == 2
    assert r["first_voxel_s"] == [15.0, 16.0]
    assert r["local_optimizations"] == 2
    assert r["total_processing_s"] == 40.0
    assert r["rejected_seeds"] == 1


def test_parse_server_log() -> None:
    r = PL.parse_server_log("x\nTotal Time elapsed: 77\n")
    assert r["total_time_s"] == 77.0


def test_parse_stderr_read_time() -> None:
    txt = "100.5 Before initialization\n101.0 Before reading experimental data\n"
    txt += "131.25 After reading experimental data\n140.0 Sending raw data\n"
    r = PL.parse_stderr(txt)
    assert r["read_s"] == pytest.approx(30.25)
    assert r["after_read_epoch"] == 131.25


def test_unit_costs() -> None:
    u = PL.unit_costs(global_s=30.0, bfs_s=130.0, n_global=3, n_local=100, n_evals=2_000_000)
    assert u["s_per_global"] == pytest.approx(10.0)
    assert u["s_per_local"] == pytest.approx(1.0)
    assert u["us_per_eval"] == pytest.approx(65.0)


def _grid(tmp_path: Path, pts: np.ndarray) -> Path:
    p = tmp_path / "g.mic"
    rows = ["0.600000"] + [
        f"{x:.6f}\t{y:.6f}\t0.000000\t{int(i % 2) + 1}\t6\t1\t0\t0\t0\t1"
        for i, (x, y) in enumerate(pts)
    ]
    p.write_text("\n".join(rows) + "\n")
    return p


def test_read_grid_positions(tmp_path: Path) -> None:
    pts = np.array([[0.1, 0.2], [0.3, -0.4]])
    g = SCORE.read_positions(_grid(tmp_path, pts))
    assert g.shape == (2, 3)
    assert np.allclose(g[:, :2], pts) and list(g[:, 2]) == [1, 2]


def test_match_voxels_by_centroid_and_updown() -> None:
    region = np.array([[0.1, 0.2, 1], [0.1, 0.2, 2], [0.5, 0.5, 1]], dtype=float)
    out = np.array([[0.5, 0.5, 1], [0.1, 0.2, 2]], dtype=float)
    idx, missing = SCORE.match_voxels(region, out, tol=1e-3)
    assert list(idx) == [2, 1]
    assert list(missing) == [0]


def test_match_voxels_rejects_unmatched_output() -> None:
    region = np.array([[0.1, 0.2, 1]], dtype=float)
    out = np.array([[0.9, 0.9, 1]], dtype=float)
    with pytest.raises(AssertionError):
        SCORE.match_voxels(region, out, tol=1e-3)


def test_subgrid_takes_nearest_rows_in_grid_order(tmp_path: Path) -> None:
    import make_subgrid as MS

    pts = np.array([[0.0, 0.0], [5.0, 5.0], [0.1, 0.0], [0.0, 0.2], [9.0, 9.0]])
    g = _grid(tmp_path, pts)
    text, rows = MS.make(g, centre=0, n=3)
    assert list(rows) == [0, 2, 3]
    lines = text.splitlines()
    assert lines[0] == "0.600000" and len(lines) == 4
    assert lines[2].split()[0] == "0.100000"


def test_stack_vs_ascii_pixels(tmp_path: Path) -> None:
    import stacks_check as SK

    f = tmp_path / "frame"
    f.write_text("3, 1, 100\n4, 1, 0\n5, 2, 7\n")
    frame = np.zeros((4, 8), dtype=np.uint8)
    frame[1, 3] = 1  # row 1, column 3: ASCII "x, y" = "3, 1"
    frame[2, 5] = 1
    assert SK.ascii_pixels(f) == {(3, 1), (5, 2)}
    assert SK.stack_pixels(frame, swap=False) == {(3, 1), (5, 2)}
    assert SK.stack_pixels(frame, swap=True) != SK.ascii_pixels(f)


def test_split_wrong_overlap() -> None:
    py = np.array([True, True, False, False])
    cpp = np.array([True, False, True, False])
    r = SCORE.paired_wrong(py, cpp)
    assert r == dict(both=1, only_cpp=1, only_python=1, neither=1)


def _cpp_run(d: Path, mic: bytes, log: str) -> Path:
    d.mkdir(parents=True)
    (d / "cpp_pilot.mic.LBFS").write_bytes(mic)
    (d / "ParallelReconLog_1.log").write_text(log)
    return d


def test_same_cpp_run_ignores_clock(tmp_path: Path) -> None:
    a = _cpp_run(tmp_path / "a", b"m", CLIENT)
    slower = (
        CLIENT.replace("15 Sec", "16 Sec")
        .replace("40...", "41...")
        .replace("00:37", "00:38")
        .replace("(sec): 77", "(sec): 78")
    )
    b = _cpp_run(tmp_path / "b", b"m", slower)
    assert SC.diff_runs(a, b) is None


def test_same_cpp_run_flags_mic_and_log(tmp_path: Path) -> None:
    a = _cpp_run(tmp_path / "a", b"m", CLIENT)
    b = _cpp_run(tmp_path / "b", b"n", CLIENT)
    assert "LBFS" in str(SC.diff_runs(a, b))
    c = _cpp_run(tmp_path / "c", b"m", CLIENT + "QueueSize 1 Fitted 9\n")
    assert "log" in str(SC.diff_runs(a, c))


def test_cli_same_cpp_run_exit_code(tmp_path: Path) -> None:
    a = _cpp_run(tmp_path / "a", b"m", CLIENT)
    b = _cpp_run(tmp_path / "b", b"m", CLIENT)
    r = subprocess.run([sys.executable, str(SCRIPTS / "same_cpp_run.py"), str(a), str(b)])
    assert r.returncode == 0


def test_quiet_verdict_load_and_busy() -> None:
    import quiet as Q

    ok = [dict(t=0, load=1.2), dict(t=30, load=2.4)]
    assert Q.verdict(ok, [], limit=2.5)["quiet"] is True
    hot = ok + [dict(t=60, load=2.6)]
    v = Q.verdict(hot, [], limit=2.5)
    assert v["quiet"] is False and v["max_load"] == pytest.approx(2.6)
    busy = [dict(t=0, procs=[dict(pid=9, cpu=35.0, command="/x/cpptools")])]
    v = Q.verdict(ok, busy, limit=2.5)
    assert v["quiet"] is False and v["busy"][0]["pid"] == 9


def test_quiet_limit_rule() -> None:
    import quiet as Q

    assert Q.choose_limit(idle=1.0, own=1.1) == 2.5
    assert Q.choose_limit(idle=1.4, own=1.3) == pytest.approx(1.5 + 1.3)
