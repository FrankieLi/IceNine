"""Tests of the pure parts of the item 14 scripts (scripts/autoresearch/item14/).

Item 14 and 3: both speed keys together (FusedCost + SmoothCostRank), with a timing A/B. The
tests pin: the identity check of two runs (same_run), the E0 arm bookkeeping (e0_arms), the
quiet verdict of the load sampler (quiet), the section shares of split, and the config file.
"""

import importlib.util
import json
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict

import numpy as np
import pytest

ICE = Path(__file__).resolve().parents[1]
SCRIPTS = ICE / "scripts" / "autoresearch" / "item14"


def _load(name: str) -> Any:
    """Load scripts/autoresearch/item14/<name>.py by path under a unique module name, so that
    another script folder's module of the same name (summary, split, ...) cannot be used."""
    spec = importlib.util.spec_from_file_location(f"item14_{name}", SCRIPTS / f"{name}.py")
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


E = _load("e0_arms")
Q = _load("quiet")
S = _load("same_run")
SP = _load("split")
SM = _load("summary")

RECORD: Dict[str, Any] = dict(voxel_idx=3, n_evals=10, wall_s=1.5, hit_ratio=0.9)


def _make_run(d: Path, wall: float = 10.0, mic: bytes = b"mic", n_evals: int = 10) -> Path:
    d.mkdir(parents=True)
    (d / "recon.mic").write_bytes(mic)
    rec = {"3": dict(RECORD, n_evals=n_evals, wall_s=wall)}
    (d / "records.json").write_text(json.dumps(rec))
    stats = dict(
        n_seeds=2,
        wall_seed_s=wall,
        wall_total_s=wall,
        reseeds=[dict(voxel=1, seconds=wall, full_search_seconds=wall, evaluations=5)],
    )
    run = dict(wall_total_s=wall, setup_s=wall, timing_label="x" * int(wall), stats=stats)
    (d / "run.json").write_text(json.dumps(run))
    return d


def test_same_run_ignores_the_clock_only(tmp_path: Path) -> None:
    a = _make_run(tmp_path / "a", wall=10.0)
    b = _make_run(tmp_path / "b", wall=11.0)
    assert S.diff_runs(a, b) is None


def test_same_run_names_a_changed_count(tmp_path: Path) -> None:
    a = _make_run(tmp_path / "a")
    b = _make_run(tmp_path / "b", n_evals=11)
    msg = S.diff_runs(a, b)
    assert msg is not None and "records.json" in msg and "n_evals" in msg


def test_same_run_names_a_changed_mic(tmp_path: Path) -> None:
    a = _make_run(tmp_path / "a")
    b = _make_run(tmp_path / "b", mic=b"mid")
    msg = S.diff_runs(a, b)
    assert msg is not None and "recon.mic" in msg


def test_same_run_exit_code(tmp_path: Path) -> None:
    a = _make_run(tmp_path / "a")
    b = _make_run(tmp_path / "b", n_evals=12)
    c = _make_run(tmp_path / "c", wall=3.0)
    py = [sys.executable, str(SCRIPTS / "same_run.py")]
    assert subprocess.run([*py, str(a), str(c)], capture_output=True).returncode == 0
    r = subprocess.run([*py, str(a), str(b)], capture_output=True, text=True)
    assert r.returncode == 1 and "differ" in (r.stdout + r.stderr)


def test_strip_clock_drops_nested_wall_fields() -> None:
    obj = {"wall_x": 1, "a": [{"seconds": 2, "keep": 3}], "stats": {"wall_s": 4, "n": 5}}
    assert S.strip_clock(obj) == {"a": [{"keep": 3}], "stats": {"n": 5}}


def test_e0_arm_tags() -> None:
    assert E.arm_tag("off", True) == "old_off_fused"
    assert E.arm_tag("E3+E4,tau=0.25", True) == "old_e34t25_fused"
    assert E.arm_tag("E3+E4,tau=0.25", False) == "old_e34t25"
    assert E.stored_tag("off") == "old_base"
    assert E.stored_tag("E3+E4,tau=0.25") == "old_e34t25"


def _voxel(**kw: Any) -> Dict[str, Any]:
    d: Dict[str, Any] = dict(
        voxel=np.int64(5),
        err_deg=np.float64(0.2),
        evals=np.int64(100),
        smooth_evals=np.int64(0),
        R_final=np.eye(3),
        cost_final=np.float64(0.1),
    )
    d.update(kw)
    return d


def test_diff_voxel() -> None:
    assert E.diff_voxel(_voxel(), _voxel()) is None
    assert "evals" in (E.diff_voxel(_voxel(evals=np.int64(101)), _voxel()) or "")
    r = np.eye(3)
    r[0, 1] = 1e-9
    assert "R_final" in (E.diff_voxel(_voxel(R_final=r), _voxel()) or "")


def test_wrong_count() -> None:
    errs = np.array([0.1, 1.5, 1.0, 3.0])
    assert E.n_wrong(errs) == 2  # over 1 degree: strictly greater


def test_quiet_verdict_load_and_busy() -> None:
    ok = [dict(t=0, load=1.2), dict(t=30, load=2.4)]
    assert Q.verdict(ok, [], limit=2.5)["quiet"] is True
    hot = ok + [dict(t=60, load=2.6)]
    v = Q.verdict(hot, [], limit=2.5)
    assert v["quiet"] is False and v["max_load"] == pytest.approx(2.6)
    busy = [dict(t=0, procs=[dict(pid=9, cpu=35.0, command="/x/cpptools")])]
    v = Q.verdict(ok, busy, limit=2.5)
    assert v["quiet"] is False and v["busy"][0]["pid"] == 9


def test_quiet_limit_rule() -> None:
    # idle + own under 2.5: keep 2.5; over: preflight limit + own
    assert Q.choose_limit(idle=1.0, own=1.1) == 2.5
    assert Q.choose_limit(idle=1.4, own=1.3) == pytest.approx(1.5 + 1.3)


def test_shares_sum_to_one() -> None:
    s = SP.shares(wall=10.0, screen=3.0, smooth=4.0, strict=1.0)
    assert s["rest"] == pytest.approx(0.2)
    assert sum(s.values()) == pytest.approx(1.0)


def test_both_keys_config_is_arm_a_plus_one_line() -> None:
    a = (ICE / "scripts/autoresearch/ar2/key_off_fused.config").read_text().rstrip("\n").split("\n")
    b = (SCRIPTS / "both_keys.config").read_text().rstrip("\n").split("\n")
    assert [x for x in b if x not in a] == ["SmoothCostRank E3+E4,tau=0.25"]
    assert [x for x in a if x not in b] == []
    assert "FusedCost 1" in b


def test_run_script_parses() -> None:
    r = subprocess.run(["bash", "-n", str(SCRIPTS / "run_item14.sh")], capture_output=True)
    assert r.returncode == 0


def _proc(path: str, cpu: float = 99.0) -> Dict[str, Any]:
    return dict(pid=7, cpu=cpu, command=path)


def test_verdict_strict_has_empty_label() -> None:
    v = Q.verdict([dict(t=0, load=1.0)], [], limit=2.5)
    assert v["quiet"] is True and v["relaxed_ok"] is True and v["label"] == ""
    assert Q.exit_code(v) == 0


def test_verdict_relaxed_daemon_and_load_spike() -> None:
    loads = [dict(t=0, load=1.2), dict(t=30, load=3.1), dict(t=60, load=1.3)]
    procs = [dict(t=0, procs=[_proc("/usr/libexec/spotlightknowledged")])]
    v = Q.verdict(loads, procs, limit=2.5)
    assert v["quiet"] is False and v["relaxed_ok"] is True
    assert v["label"] == (
        "strict quiet check failed (load max 3.10 > 2.5 in 1 of 3 samples; "
        "busy: spotlightknowledged)"
    )
    assert Q.exit_code(v) == 3


def test_verdict_system_prefixes() -> None:
    for path in ("/usr/libexec/x", "/System/Library/y/z", "/sbin/w"):
        v = Q.verdict([dict(t=0, load=1.0)], [dict(procs=[_proc(path)])], limit=2.5)
        assert v["relaxed_ok"] is True, path


def test_verdict_neither_non_daemon_or_median_over() -> None:
    v = Q.verdict([dict(t=0, load=1.0)], [dict(procs=[_proc("/opt/app/cpptools")])], limit=2.5)
    assert v["relaxed_ok"] is False and Q.exit_code(v) == 1
    mixed = [dict(procs=[_proc("/usr/libexec/d"), _proc("/Applications/App")])]
    assert Q.verdict([dict(t=0, load=1.0)], mixed, limit=2.5)["relaxed_ok"] is False
    high = [dict(t=i, load=3.0) for i in range(3)]
    v = Q.verdict(high, [dict(procs=[_proc("/usr/libexec/d")])], limit=2.5)
    assert v["relaxed_ok"] is False and Q.exit_code(v) == 1
    # R4b (i) is a median of at most 2.5 even when the strict limit is larger (cpptime: 3.5)
    mid = [dict(t=i, load=3.0) for i in range(3)]
    v = Q.verdict(mid, [dict(procs=[_proc("/usr/libexec/d")])], limit=3.5)
    assert v["median_load"] == 3.0 and v["relaxed_ok"] is False
    v = Q.verdict([dict(t=i, load=2.4) for i in range(3)], [], limit=3.5)
    assert v["relaxed_ok"] is True


def test_pair_check() -> None:
    p = Q.pair_report("A", 1000.0, 1010.0, "", "")
    assert p["reported"] is True and p["labels"] == []
    p = Q.pair_report("B", 4494.29, 4480.81, "strict quiet check failed (d)", "")
    assert p["reported"] is True and p["labels"] == ["strict quiet check failed (d)"]
    p = Q.pair_report("B", 1000.0, 1030.0, "", "strict quiet check failed (d)")
    assert p["reported"] is False and p["note"] == "no time reported for arm B"
    p = Q.pair_report("A", 1000.0, 1030.0, "", "")  # unlabelled pair: no check needed
    assert p["reported"] is True


def test_label_when_nothing_busy_was_caught() -> None:
    loads = [dict(t=i, load=v) for i, v in enumerate([1.0, 4.46, 1.2, 3.0, 1.1, 1.3])]
    v = Q.verdict(loads, [dict(procs=[])], limit=2.5)
    assert v["label"] == (
        "strict quiet check failed (load max 4.46 > 2.5 in 2 of 6 samples; "
        "no process over 20% CPU caught)"
    )
    assert v["relaxed_ok"] is True and Q.exit_code(v) == 3


def test_system_applications_is_not_a_daemon() -> None:
    p = [dict(procs=[_proc("/System/Applications/Photos.app/Contents/MacOS/Photos")])]
    v = Q.verdict([dict(t=0, load=1.0)], p, limit=2.5)
    assert v["relaxed_ok"] is False


def test_summary_run_label() -> None:
    assert SM.run_label(dict(strict_quiet=True, label="")) == "quiet"
    assert SM.run_label(dict(strict_quiet=False, label="strict quiet check failed (x)")) == (
        "strict quiet check failed (x)"
    )


def test_summary_leaves_out_an_unreported_arm() -> None:
    def run(t: float) -> Dict[str, float]:
        return dict(real_s=t, wall_total_s=t, full_s=t - 1, local_s=1.0)

    runs = {"A1": run(100.0), "A2": run(101.0), "B1": run(300.0), "B2": run(400.0)}
    ok = dict(reported=True)
    arm, ratio, share = SM.arm_times(runs, dict(A=ok, B=ok))
    assert set(arm) == {"A", "B"} and ratio == pytest.approx(350.0 / 100.5)
    arm, ratio, share = SM.arm_times(runs, dict(A=ok, B=dict(reported=False)))
    assert set(arm) == {"A"} and np.isnan(ratio) and np.isnan(share)
