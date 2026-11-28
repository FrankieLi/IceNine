"""Tests of the pure parts of the cpptime scripts (scripts/autoresearch/cpptime/).

C++ time on the pilot: one client. The tests pin the config conversion (make_config), the stop
rules of the guard, the log parsing (parse_logs), the voxel matching and scoring (score_cpp), the
identity check of two C++ runs (same_cpp_run) and the quiet-limit choice.
"""

import importlib.util
import json
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest

ICE = Path(__file__).resolve().parents[1]
SCRIPTS = ICE / "scripts" / "autoresearch" / "cpptime"
ITEM14 = ICE / "scripts" / "autoresearch" / "item14"


def _load(folder: Path, name: str, as_name: str) -> Any:
    """Load folder/name.py by path under a unique module name, so that another script folder's
    module of the same name (summary, quiet, ...) cannot be picked up in the full suite."""
    spec = importlib.util.spec_from_file_location(as_name, folder / f"{name}.py")
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


# the scripts import each other by bare name (score_cpp, make_region, ...): this folder first
sys.path.insert(0, str(SCRIPTS))
G = _load(SCRIPTS, "guard", "cpptime_guard")
MC = _load(SCRIPTS, "make_config", "cpptime_make_config")
PL = _load(SCRIPTS, "parse_logs", "cpptime_parse_logs")
SC = _load(SCRIPTS, "same_cpp_run", "cpptime_same_cpp_run")
SCORE = _load(SCRIPTS, "score_cpp", "cpptime_score_cpp")
MS = _load(SCRIPTS, "make_subgrid", "cpptime_make_subgrid")
SK = _load(SCRIPTS, "stacks_check", "cpptime_stacks_check")
TR = _load(SCRIPTS, "trace_cpp_log", "cpptime_trace_cpp_log")
SM = _load(SCRIPTS, "summary", "cpptime_summary")
Q = _load(ITEM14, "quiet", "cpptime_item14_quiet")  # the one shared quiet.py

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

    pts = np.array([[0.0, 0.0], [5.0, 5.0], [0.1, 0.0], [0.0, 0.2], [9.0, 9.0]])
    g = _grid(tmp_path, pts)
    text, rows = MS.make(g, centre=0, n=3)
    assert list(rows) == [0, 2, 3]
    lines = text.splitlines()
    assert lines[0] == "0.600000" and len(lines) == 4
    assert lines[2].split()[0] == "0.100000"


def test_stack_vs_ascii_pixels(tmp_path: Path) -> None:

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


def test_same_cpp_run_ignores_only_the_fitting_time_column(tmp_path: Path) -> None:
    head = "0.6\n"
    row = "0.1 0.2 0 1 6 1 10 20 30 1 0.125 1 {} 1 1 1 0 0 0\n"  # index 12 is fFittingTime
    a = _cpp_run(tmp_path / "a", (head + row.format("0")).encode(), CLIENT)
    b = _cpp_run(tmp_path / "b", (head + row.format("1.4e-45")).encode(), CLIENT)
    assert SC.diff_runs(a, b) is None
    c = _cpp_run(tmp_path / "c", (head + row.format("0").replace(" 30 ", " 31 ")).encode(), CLIENT)
    assert "LBFS" in str(SC.diff_runs(a, c))
    d = _cpp_run(
        tmp_path / "d", (head + row.format("0").replace(" 1 0 0 0", " 1 0 0 1")).encode(), CLIENT
    )
    assert "LBFS" in str(SC.diff_runs(a, d))  # a deformation column is still checked


LOG_RECORDS = """First Voxel Required: 13 Sec
Fitting neighbors via Breadth First Search
QueueSize 4 Fitted 2 continuous voxels
  [ Cost ] 0.05 |=| 1 1    Original : 1 2 3    Optimized: 1 2 3 VoxelLocation 0.18 -0.08 0
QueueSize 5 Fitted 3 continuous voxels
  [ Cost ] 0.29 |=| 1 0.7    Original : 1 2 3    Optimized: 1 2 4 VoxelLocation 0.19 -0.08 0
First Voxel Required: 14 Sec
QueueSize 4 Fitted 2 continuous voxels
  [ Cost ] 0.05 |=| 0.95 0.9    Original : 5 6 7    Optimized: 5 6 7 VoxelLocation 0.1 0.2 0
"""


def test_trace_parse_records() -> None:
    recs = TR.parse_records(LOG_RECORDS)
    assert [r["seed"] for r in recs] == [0, 0, 1]
    assert [r["refit"] for r in recs] == [False, True, False]  # 0.7/1 <= 0.9; 0.9/0.95 > 0.9
    assert (
        recs[1]["original"] != recs[1]["optimized"] and recs[0]["original"] == recs[0]["optimized"]
    )
    assert (recs[0]["x"], recs[0]["y"]) == (0.18, -0.08)
    s = TR.summarize(recs, {0: False, 1: True, 2: False}, [0, 1, 2])
    assert s["refit"] == 1 and s["refit_wrong"] == 1 and s["accepted_wrong"] == 0
    assert s["unchanged"] == 2


def test_summary_timing_label_wording() -> None:
    base = dict(
        busy=[],
        samples=dict(
            n_over_2_5=7,
            n_load=43,
            n_over_limit=3,
            max_load=6.2,
            proc_period_s=300,
            load_period_s=30,
        ),
    )
    failed = SM.timing_label(dict(base, strict_quiet=False))
    assert "no process over 20% CPU caught by the sampler" in failed
    assert "load samples above 2.5 (7 of 43" in failed and "sampled every 300 s" in failed
    assert "condition (ii) of N1 held only because nothing was caught" in failed
    strict = SM.timing_label(dict(base, strict_quiet=True))
    assert "assumed idle load of 1.5" in strict and "limit 3.5 only" in strict


def test_cli_same_cpp_run_exit_code(tmp_path: Path) -> None:
    a = _cpp_run(tmp_path / "a", b"m", CLIENT)
    b = _cpp_run(tmp_path / "b", b"m", CLIENT)
    r = subprocess.run([sys.executable, str(SCRIPTS / "same_cpp_run.py"), str(a), str(b)])
    assert r.returncode == 0


def test_quiet_verdict_load_and_busy() -> None:

    ok = [dict(t=0, load=1.2), dict(t=30, load=2.4)]
    assert Q.verdict(ok, [], limit=2.5)["quiet"] is True
    hot = ok + [dict(t=60, load=2.6)]
    v = Q.verdict(hot, [], limit=2.5)
    assert v["quiet"] is False and v["max_load"] == pytest.approx(2.6)
    busy = [dict(t=0, procs=[dict(pid=9, cpu=35.0, command="/x/cpptools")])]
    v = Q.verdict(ok, busy, limit=2.5)
    assert v["quiet"] is False and v["busy"][0]["pid"] == 9


def _proc(path: str, cpu: float = 99.0) -> dict:
    return dict(pid=7, cpu=cpu, command=path)


def test_r4b_verdict_strict_has_empty_label() -> None:

    v = Q.verdict([dict(t=0, load=1.0)], [], limit=2.5)
    assert v["quiet"] is True and v["relaxed_ok"] is True and v["label"] == ""
    assert Q.exit_code(v) == 0


def test_r4b_verdict_relaxed_daemon_and_load_spike() -> None:

    loads = [dict(t=0, load=1.2), dict(t=30, load=3.1), dict(t=60, load=1.3)]
    procs = [dict(t=0, procs=[_proc("/usr/libexec/spotlightknowledged")])]
    v = Q.verdict(loads, procs, limit=2.5)
    assert v["quiet"] is False and v["relaxed_ok"] is True
    assert (
        v["label"] == "strict quiet check failed (load max 3.10 > 2.5 in 1 of 3 samples; "
        "busy: spotlightknowledged)"
    )
    assert Q.exit_code(v) == 3


def test_r4b_verdict_neither() -> None:

    v = Q.verdict([dict(t=0, load=1.0)], [dict(procs=[_proc("/opt/app/cpptools")])], limit=2.5)
    assert v["relaxed_ok"] is False and Q.exit_code(v) == 1


def test_r4b_pair_check() -> None:

    assert Q.pair_report("C", 1000.0, 1010.0, "", "")["reported"] is True
    p = Q.pair_report("C", 1000.0, 1030.0, "", "strict quiet check failed (d)")
    assert p["reported"] is False and p["note"] == "no time reported for arm C"
    assert Q.pair_report("C", 1000.0, 1030.0, "", "")["reported"] is True


def test_run_script_parses_and_has_r4b() -> None:
    r = subprocess.run(["bash", "-n", str(SCRIPTS / "run_cpp_time.sh")], capture_output=True)
    assert r.returncode == 0
    text = (SCRIPTS / "run_cpp_time.sh").read_text()
    assert "quiet.py pair" in text and "LABEL" in text and "pair C C1 C2" in text


def test_quiet_limit_rule() -> None:

    assert Q.choose_limit(idle=1.0, own=1.1) == 2.5
    assert Q.choose_limit(idle=1.4, own=1.3) == pytest.approx(1.5 + 1.3)
