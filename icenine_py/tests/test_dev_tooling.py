"""Tests of the developer tooling: hooks, scripts/dev and scripts/common helpers."""

import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pytest
from scipy.spatial.transform import Rotation
from scipy.stats import binomtest

REPO = Path(__file__).resolve().parents[2]
DEV = REPO / "scripts" / "dev"
sys.path.insert(0, str(REPO / "icenine_py" / "scripts" / "common"))
sys.path.insert(0, str(DEV))

import audit_numbers  # noqa: E402
import doc_tables  # noqa: E402
import precommit_check  # noqa: E402
import preflight  # noqa: E402
import stats  # noqa: E402
import sync_doc_tables  # noqa: E402

HAVE_JQ_PERL = shutil.which("jq") is not None and shutil.which("perl") is not None
GIT_ENV = {
    "GIT_AUTHOR_NAME": "t",
    "GIT_AUTHOR_EMAIL": "t@example.com",
    "GIT_COMMITTER_NAME": "t",
    "GIT_COMMITTER_EMAIL": "t@example.com",
}


def _git(repo: Path, *args: str) -> str:
    env = {**os.environ, **GIT_ENV}
    out = subprocess.run(["git", *args], cwd=repo, env=env, capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    return out.stdout.strip()


def _init_repo(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    _git(path, "init", "-q", "-b", "develop")
    (path / "README").write_text("x\n")
    _git(path, "add", "README")
    _git(path, "commit", "-q", "-m", "init")
    return path


# ---------------------------------------------------------------------------------------------
# A. Bash PreToolUse hook
# ---------------------------------------------------------------------------------------------

HOOK = REPO / ".claude" / "hooks" / "check-bash-command.sh"

HOOK_CASES = [
    # (command, blocked)
    ("python3 script.py", True),
    ("pytest tests/", True),
    ("pip install foo", True),
    ("cd icenine_py && pytest -q", True),
    ("ls; python x.py", True),
    ("cat f | python3 -c 'print(1)'", True),
    ("echo hi\npytest tests/", True),
    ("git add -A", True),
    ("git add --all", True),
    ("git add .", True),
    ("git add -u", True),
    ("git commit -a -m 'x'", True),
    ("git commit -am 'x'", True),
    ("git commit --all -m x", True),
    ("git status && git add -A && git commit -m ok", True),
    ("git -C /tmp/repo add .", True),
    ("uv run python script.py", False),
    ("uv run pytest tests/", False),
    ("uv pip install -e .", False),
    ("git add path/with-A-in-name.py", False),
    ("git add icenine_py/scripts/ALL-data.py src/-A.txt", False),
    ("git add -u icenine_py/foo.py", False),
    ('git commit -m "fix the -a flag handling"', False),
    ('git commit -m "use git add -A never" && echo done', False),
    ("git commit --amend -m 'x'", False),
    ("git add ./foo/bar.py", False),
    ('echo "run pytest later"', False),
    ("ls -la && git status", False),
    ("git commit -m \"$(cat <<'EOF'\nmsg with git add . inside\nEOF\n)\"", False),
    # review round: continuations, heredocs, wrappers, subshells, versions, :/ and overrides
    ("uv run \\\n  pytest tests/", False),
    ("cat > run.sh <<EOF\npython x.py\npytest\nEOF\nls", False),
    ("cat <<'EOF' > x.md\ngit add -A\nEOF", False),
    ("python3 \\\n  script.py", True),
    ("nohup python3 x.py &", True),
    ("nice -n 19 python3 x.py", True),
    ("nice -n 19 uv run python x.py", False),
    ("env -u VIRTUAL_ENV python3 x.py", True),
    ("env -u VIRTUAL_ENV uv run pytest", False),
    ("FOO=1 BAR=2 pytest -q", True),
    ("time pip install x", True),
    ("sudo -n pip3 install x", True),
    ("(python3 x.py)", True),
    ("{ pytest; }", True),
    ("echo $(python3 -c 1)", True),
    ("pip3 install x", True),
    ("python3.11 x.py", True),
    ("git add :/", True),
    ("git add '*'", False),  # quoted glob is stripped; a bare * is blocked below
    ("git add *", True),
    ("git commit --no-verify -m x", True),
    ("git commit -nm x", True),
    ("git commit -m 'about --no-verify'", False),
    ("ALLOW_LARGE=1 git commit -m x", True),
    ("ALLOW_CLAUDE_CONFIG=1 git commit -m x", True),
    ("export ALLOW_ABS_PATHS=1", True),
    ("echo ALLOW_LARGE", False),
]


@pytest.mark.skipif(not HAVE_JQ_PERL, reason="needs jq and perl")
@pytest.mark.parametrize("command,blocked", HOOK_CASES)
def test_bash_hook(command: str, blocked: bool) -> None:
    payload = json.dumps({"tool_input": {"command": command}})
    out = subprocess.run(
        [str(HOOK)], input=payload, capture_output=True, text=True, check=True
    ).stdout.strip()
    if blocked:
        spec = json.loads(out)["hookSpecificOutput"]
        assert spec["hookEventName"] == "PreToolUse" and spec["permissionDecision"] == "deny"
        assert spec["permissionDecisionReason"], command
    else:
        assert out == "", command


@pytest.mark.skipif(not HAVE_JQ_PERL, reason="needs jq and perl")
def test_bash_hook_empty_input() -> None:
    out = subprocess.run([str(HOOK)], input="{}", capture_output=True, text=True, check=True)
    assert out.stdout.strip() == ""


def test_settings_uses_hook_script_and_tests_hook_uses_uv() -> None:
    settings = json.loads((REPO / ".claude" / "settings.json").read_text())
    pre = settings["hooks"]["PreToolUse"][0]["hooks"][0]
    assert pre["command"] == '"$CLAUDE_PROJECT_DIR"/.claude/hooks/check-bash-command.sh'
    assert pre["timeout"] == 5
    run_tests = (REPO / ".claude" / "hooks" / "run-tests-on-change.sh").read_text()
    assert "&& pytest" not in run_tests and "uv run pytest" in run_tests


# ---------------------------------------------------------------------------------------------
# B. pre-commit checks
# ---------------------------------------------------------------------------------------------


def _stage(repo: Path, rel: str, content: str, binary_size: Optional[int] = None) -> None:
    p = repo / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    if binary_size is not None:
        p.write_bytes(os.urandom(binary_size))
    else:
        p.write_text(content)
    _git(repo, "add", rel)


@pytest.fixture()
def repo(tmp_path: Path) -> Path:
    return _init_repo(tmp_path / "r")


def test_precommit_clean(repo: Path) -> None:
    _stage(repo, "notes.md", "all fine, realistic data\n")
    assert precommit_check.check(repo, {}) == []


def test_precommit_forbidden_paths(repo: Path) -> None:
    for rel in ["m.pt", "d/x.pkl", "d/y.joblib", "a/cache/z.txt"]:
        _stage(repo, rel, "x")
    msgs = precommit_check.check(repo, {})
    for rel in ["m.pt", "d/x.pkl", "d/y.joblib", "a/cache/z.txt"]:
        assert any(m.startswith(rel + ":") for m in msgs), rel


def test_precommit_claude_config_override(repo: Path) -> None:
    _stage(repo, "CLAUDE.md", "x\n")
    _stage(repo, ".claude/agents/a.md", "x\n")
    assert len(precommit_check.check(repo, {})) == 2
    assert precommit_check.check(repo, {"ALLOW_CLAUDE_CONFIG": "1"}) == []


def test_precommit_large_files(repo: Path) -> None:
    _stage(repo, "big.bin", "", binary_size=11 * 1024 * 1024)
    _stage(repo, "mid.npz", "", binary_size=6 * 1024 * 1024)
    msgs = precommit_check.check(repo, {})
    assert any(m.startswith("big.bin") for m in msgs)
    assert any(m.startswith("mid.npz") for m in msgs)
    assert precommit_check.check(repo, {"ALLOW_LARGE": "1"}) == []


def test_precommit_abs_paths(repo: Path) -> None:
    home = "/" + "Users" + "/someone/data"
    _stage(repo, "a.md", f"see {home}\nfine\n")
    msgs = precommit_check.check(repo, {})
    assert msgs == ["a.md:1: absolute home path (set ALLOW_ABS_PATHS=1)"]
    assert precommit_check.check(repo, {"ALLOW_ABS_PATHS": "1"}) == []


def test_precommit_python_format(repo: Path) -> None:
    _stage(repo, "icenine_py/new_bad.py", "x = {  'a':1 }\n")
    _stage(repo, "icenine_py/new_good.py", 'x = {"a": 1}\n')
    msgs = precommit_check.check(repo, {})
    assert msgs == ["icenine_py/new_bad.py:1: new file fails black --check -l 100"]


def test_precommit_long_added_line_in_modified_file(repo: Path) -> None:
    _stage(repo, "icenine_py/old.py", "# " + "x" * 120 + "\n")  # legacy content, committed
    _git(repo, "commit", "-q", "-m", "legacy")
    (repo / "icenine_py" / "old.py").write_text(
        "# " + "x" * 120 + "\n" + "y = '" + "z" * 100 + "'\n"
    )
    _git(repo, "add", "icenine_py/old.py")
    msgs = precommit_check.check(repo, {})
    assert msgs == ["icenine_py/old.py:2: line longer than 100 characters"]


def test_precommit_realistic_terminology(repo: Path) -> None:
    body = "\n".join(
        [
            "the corrupted data",  # 1 flagged  # noqa: realistic
            "Corruption of the beam",  # 2 flagged  # noqa: realistic
            "corrupt_windows(x)",  # 3 alias
            "--corrupt-train flag",  # 4 alias
            "a corrupt thing  noqa: realistic",  # 5 noqa
            "uncorrupted is fine",  # 6 not a word match
            "it corrupts and is corrupting",  # 7 flagged  # noqa: realistic
            "corrupt_windows and corrupted together",  # 8  # noqa: realistic
            "",
        ]
    )
    _stage(repo, "t.md", body)
    msgs = precommit_check.check(repo, {})
    assert [m.split(":")[1] for m in msgs] == ["1", "2", "7", "8"]


def test_precommit_cli_silent_when_clean_and_exit_code(repo: Path) -> None:
    script = DEV / "precommit_check.py"
    _stage(repo, "ok.md", "fine\n")
    res = subprocess.run([sys.executable, str(script)], cwd=repo, capture_output=True, text=True)
    assert res.returncode == 0 and res.stdout == ""
    _stage(repo, "bad.pt", "x")
    res = subprocess.run([sys.executable, str(script)], cwd=repo, capture_output=True, text=True)
    assert res.returncode == 1 and res.stdout.startswith("bad.pt:1:")


def test_githook_and_installer_exist_and_are_not_activated_by_default(tmp_path: Path) -> None:
    """The repo ships the hooks but never activates them: a fresh clone has no core.hooksPath.
    (A developer's own checkout may have run install_hooks.sh, so local config is not checked.)"""
    assert os.access(REPO / ".githooks" / "pre-commit", os.X_OK)
    assert os.access(DEV / "install_hooks.sh", os.X_OK)
    assert "scripts/dev/precommit_check.py" in (REPO / ".githooks" / "pre-commit").read_text()
    clone = tmp_path / "clone"
    subprocess.run(
        ["git", "clone", "-q", "--no-checkout", str(REPO), str(clone)],
        check=True,
        capture_output=True,
    )
    res = subprocess.run(
        ["git", "config", "--local", "--get", "core.hooksPath"],
        cwd=clone,
        capture_output=True,
        text=True,
    )
    assert res.stdout.strip() == "", "a fresh clone must not have the hooks activated"


# ---------------------------------------------------------------------------------------------
# C. task scripts
# ---------------------------------------------------------------------------------------------


def _script_env(
    repo: Path, extra: Dict[str, str], name: str, *args: str
) -> subprocess.CompletedProcess:
    env = {**os.environ, **GIT_ENV, **extra}
    return subprocess.run(
        [str(DEV / name), *args], cwd=repo, env=env, capture_output=True, text=True
    )


def _script(repo: Path, name: str, *args: str, stdin: str = "") -> subprocess.CompletedProcess:
    env = {**os.environ, **GIT_ENV}
    return subprocess.run(
        [str(DEV / name), *args],
        cwd=repo,
        env=env,
        capture_output=True,
        text=True,
        input=stdin,
    )


def test_start_and_finish_task(repo: Path) -> None:
    assert _script(repo, "start_task.sh", "t1").returncode != 0  # develop is not feature/*
    _git(repo, "checkout", "-q", "-b", "feature/x")
    res = _script(repo, "start_task.sh", "t1")
    assert res.returncode == 0, res.stderr
    assert _git(repo, "rev-parse", "--abbrev-ref", "HEAD") == "feature/x-t1"
    assert _git(repo, "config", "branch.feature/x-t1.parent") == "feature/x"
    (repo / "f.txt").write_text("hello\n")
    _git(repo, "add", "f.txt")
    _git(repo, "commit", "-q", "-m", "add f")
    res = _script(repo, "finish_task.sh", "--no-tests")
    assert res.returncode == 0, res.stderr
    assert _git(repo, "rev-parse", "--abbrev-ref", "HEAD") == "feature/x"
    assert "feature/x-t1" not in _git(repo, "branch", "--list")
    assert _git(repo, "log", "-1", "--format=%s") == "Merge task t1: add f"
    assert (repo / "f.txt").exists()


def test_finish_task_guesses_parent_with_yes(repo: Path) -> None:
    _git(repo, "checkout", "-q", "-b", "feature/y")
    _git(repo, "checkout", "-q", "-b", "feature/y-job")  # no recorded parent
    (repo / "g.txt").write_text("g\n")
    _git(repo, "add", "g.txt")
    _git(repo, "commit", "-q", "-m", "add g")
    res = _script(repo, "finish_task.sh", "--no-tests", "-m", "custom merge", stdin="y\n")
    assert res.returncode == 0, res.stderr
    assert _git(repo, "log", "-1", "--format=%s") == "custom merge"


def _task_branch(repo: Path, fname: str = "h.py") -> None:
    _git(repo, "checkout", "-q", "-b", "feature/z")
    _git(repo, "checkout", "-q", "-b", "feature/z-job")
    _git(repo, "config", "branch.feature/z-job.parent", "feature/z")
    (repo / fname).write_text("x = 1\n")
    _git(repo, "add", fname)
    _git(repo, "commit", "-q", "-m", "add " + fname)


def test_finish_task_refusals(repo: Path) -> None:
    res = _script(repo, "finish_task.sh", "--no-tests")  # develop is not a task branch
    assert res.returncode == 1 and "feature/*" in res.stderr
    _task_branch(repo)
    (repo / "README").write_text("dirty\n")  # tracked file modified
    res = _script(repo, "finish_task.sh", "--no-tests")
    assert res.returncode == 1 and "uncommitted" in res.stderr
    _git(repo, "checkout", "--", "README")
    _git(repo, "config", "branch.feature/z-job.parent", "develop")
    res = _script(repo, "finish_task.sh", "--no-tests")
    assert res.returncode == 1 and "refusing" in res.stderr
    assert _git(repo, "rev-parse", "--abbrev-ref", "HEAD") == "feature/z-job"


def test_finish_task_failing_tests_or_build_do_not_merge(repo: Path) -> None:
    _task_branch(repo)
    head = _git(repo, "rev-parse", "feature/z")
    env = {"FINISH_TASK_TEST_CMD": "exit 3"}
    res = _script_env(repo, env, "finish_task.sh")
    assert res.returncode != 0
    assert _git(repo, "rev-parse", "feature/z") == head
    assert _git(repo, "rev-parse", "--abbrev-ref", "HEAD") == "feature/z-job"
    # a C++ change triggers the build command; a failing build must stop the merge
    (repo / "src.cpp").write_text("int x;\n")
    _git(repo, "add", "src.cpp")
    _git(repo, "commit", "-q", "-m", "cpp")
    env = {"FINISH_TASK_TEST_CMD": "true", "FINISH_TASK_BUILD_CMD": "false"}
    res = _script_env(repo, env, "finish_task.sh")
    assert res.returncode != 0 and _git(repo, "rev-parse", "feature/z") == head
    env = {"FINISH_TASK_TEST_CMD": "true", "FINISH_TASK_BUILD_CMD": "true"}
    assert _script_env(repo, env, "finish_task.sh").returncode == 0
    assert _git(repo, "rev-parse", "--abbrev-ref", "HEAD") == "feature/z"


def test_finish_task_txt_does_not_trigger_build(repo: Path) -> None:
    _task_branch(repo, "notes.txt")
    env = {"FINISH_TASK_TEST_CMD": "true", "FINISH_TASK_BUILD_CMD": "false"}
    assert _script_env(repo, env, "finish_task.sh").returncode == 0


def test_finish_task_guess_rules(repo: Path) -> None:
    _git(repo, "checkout", "-q", "-b", "feature/q-job")  # parent feature/q does not exist
    res = _script(repo, "finish_task.sh", "--no-tests", stdin="y\n")
    assert res.returncode == 1 and "does not exist" in res.stderr
    _git(repo, "branch", "feature/q")
    res = _script(repo, "finish_task.sh", "--no-tests", "--yes")  # --yes never guesses
    assert res.returncode == 1 and "does not guess" in res.stderr
    res = _script(repo, "finish_task.sh", "--no-tests", stdin="n\n")
    assert res.returncode == 1


# ---------------------------------------------------------------------------------------------
# D. stats
# ---------------------------------------------------------------------------------------------


def test_wilson_known_values() -> None:
    lo, hi = stats.wilson(50, 100)
    assert lo == pytest.approx(0.4038, abs=1e-3) and hi == pytest.approx(0.5962, abs=1e-3)
    lo, hi = stats.wilson(0, 10)
    assert lo == 0.0 and hi == pytest.approx(0.2775, abs=1e-3)
    lo, hi = stats.wilson(204, 600)
    assert (lo, hi) == (pytest.approx(0.3028, abs=1e-3), pytest.approx(0.3786, abs=1e-3))
    assert all(math.isnan(v) for v in stats.wilson(0, 0))


@pytest.mark.parametrize("b,c", [(0, 0), (3, 3), (1, 9), (12, 4), (0, 7)])
def test_mcnemar_matches_scipy(b: int, c: int) -> None:
    expected = 1.0 if b + c == 0 else binomtest(min(b, c), b + c, 0.5).pvalue
    assert stats.mcnemar_exact(b, c) == pytest.approx(expected)


def test_paired_discordant() -> None:
    a = [True, True, False, False, True]
    b = [True, False, True, False, False]
    assert stats.paired_discordant(a, b) == (2, 1)
    with pytest.raises(ValueError):
        stats.paired_discordant([True], [True, False])


def test_win_rate_ties_and_nan() -> None:
    a = np.array([0.1, 0.5, 0.3, np.nan, 0.2])
    b = np.array([0.2, 0.4, 0.3005, 0.1, np.nan])
    # pairs kept: (0.1,0.2) win, (0.5,0.4) loss, (0.3,0.3005) tie -> (1 + 0.5) / 3
    rate, tie_frac, n = stats.win_rate(a, b)
    assert (rate, tie_frac, n) == (pytest.approx(0.5), pytest.approx(1 / 3), 3)
    assert stats.win_rate(a, b, tie=0.0) == (pytest.approx(2 / 3), 0.0, 3)
    rate, tie_frac, n = stats.win_rate([np.nan], [1.0])
    assert math.isnan(rate) and math.isnan(tie_frac) and n == 0
    # the tie is strict: a difference of exactly `tie` is not a tie; inf pairs are dropped
    assert stats.win_rate([0.0, np.inf], [0.5, 1.0], tie=0.5) == (1.0, 0.0, 1)
    assert stats.win_rate([0.0], [0.4999], tie=0.5) == (0.5, 1.0, 1)


def test_reorder() -> None:
    vals = np.array([[10.0], [20.0], [30.0]])
    have = [5, 7, 9]
    assert stats.reorder(vals, have, [9, 5, 7]).ravel().tolist() == [30.0, 10.0, 20.0]
    assert stats.reorder(vals, have, [7]).ravel().tolist() == [20.0]  # partial
    m = np.arange(6).reshape(2, 3)
    assert stats.reorder(m, [1, 2, 3], [3, 1], axis=1).tolist() == [[2, 0], [5, 3]]
    with pytest.raises(ValueError, match="missing"):
        stats.reorder(vals, have, [5, 8])
    with pytest.raises(ValueError):
        stats.reorder(vals, [1, 1, 2], [1])


def test_misorientation_wrapper_symmetry_equivalent() -> None:
    R = Rotation.random(random_state=1).as_matrix()
    S = Rotation.from_euler("z", 90, degrees=True).as_matrix()  # a cubic symmetry operator
    assert float(stats.misorientation_deg_cubic(R, R @ S)) == pytest.approx(0.0, abs=1e-6)
    d = Rotation.from_rotvec(np.radians(2.0) * np.array([0, 0, 1.0])).as_matrix()
    assert float(stats.misorientation_deg_cubic(R, R @ d)) == pytest.approx(2.0, abs=1e-6)
    batch = np.stack([R, R @ S])
    assert stats.misorientation_deg_cubic(batch, R).shape == (2,)


# ---------------------------------------------------------------------------------------------
# E. doc tables
# ---------------------------------------------------------------------------------------------


def test_markdown_table_formatting() -> None:
    rows: List[Optional[Dict[str, object]]] = [
        {"name": "a", "x": 0.123456, "n": 12},
        {"name": "b", "x": float("nan"), "n": None},
        None,
    ]
    t = doc_tables.markdown_table(rows, ["name", "x", "n"], {"n": "d"})
    assert t.splitlines() == [
        "| name | x | n |",
        "|---|---|---|",
        "| a | 0.123 | 12 |",
        "| b | n/a | n/a |",
        "| not evaluated | not evaluated | not evaluated |",
    ]
    assert doc_tables.markdown_table(rows, ["name", "x", "n"], {"n": "d"}) == t


def test_markdown_table_ints_floats_and_row_labels() -> None:
    rows: List[Optional[Dict[str, object]]] = [
        {"name": "a", "k": 123456, "x": 123456.0, "np": np.int64(7), "f": np.float32(0.5)},
        None,
    ]
    t = doc_tables.markdown_table(rows, ["name", "k", "x", "np", "f"], labels=["a", "missing row"])
    lines = t.splitlines()
    assert lines[2] == "| a | 123456 | 1.23e+05 | 7 | 0.5 |"  # ints exact, floats 3 sig. figs
    ne = "not evaluated"
    assert lines[3] == f"| missing row | {ne} | {ne} | {ne} | {ne} |"
    # a per-column format still wins
    assert "| 1.2e+05 |" in doc_tables.markdown_table([{"x": 123456.0}], ["x"], {"x": ".2g"})


def test_sync_round_trip_idempotent_and_check(tmp_path: Path) -> None:
    gen = tmp_path / "gen.md"
    doc_tables.write_tables(gen, {"one": "| a |\n|---|\n| 1 |", "two": "| b |\n|---|\n| 2 |"})
    assert doc_tables.read_tables(gen)["one"] == "| a |\n|---|\n| 1 |"
    doc = tmp_path / "doc.md"
    doc.write_text("intro\n<!-- table:one -->\nold\n<!-- /table:one -->\nend\n")
    script = DEV / "sync_doc_tables.py"
    cmd = [sys.executable, str(script), "--doc", str(doc), "--tables", str(gen)]
    check = subprocess.run([*cmd, "--check"], capture_output=True, text=True)
    assert check.returncode == 1 and doc.read_text().count("old") == 1  # --check never writes
    assert "unused generated blocks: two" in check.stderr
    assert subprocess.run(cmd, capture_output=True, text=True).returncode == 0
    first = doc.read_text()
    assert "| 1 |" in first and "old" not in first and first.startswith("intro\n")
    assert subprocess.run(cmd, capture_output=True, text=True).returncode == 0
    assert doc.read_text() == first  # idempotent
    assert subprocess.run([*cmd, "--check"], capture_output=True, text=True).returncode == 0


def test_sync_missing_marker_errors(tmp_path: Path) -> None:
    gen = tmp_path / "gen.md"
    doc_tables.write_tables(gen, {"one": "x"})
    doc = tmp_path / "doc.md"
    doc.write_text("<!-- table:nope -->\n<!-- /table:nope -->\n")
    res = subprocess.run(
        [sys.executable, str(DEV / "sync_doc_tables.py"), "--doc", str(doc), "--tables", str(gen)],
        capture_output=True,
        text=True,
    )
    assert res.returncode == 2 and "nope" in res.stderr
    with pytest.raises(KeyError):
        sync_doc_tables.sync(doc.read_text(), {"one": "x"})


# ---------------------------------------------------------------------------------------------
# F. audit_numbers
# ---------------------------------------------------------------------------------------------

DOC = """# Report

## Intro 2026-10-06
Not in scope 777.

## Results (2026-10-06)

1. Step one uses 8 seeds in `scripts/foo_1234.py`.
Wrong rate 34.0% (CI 30.3-37.9), 204/600 runs; the cost fell by 0.0197-0.0218 over 360k points.
Another 41.7% claim, and 55 voxels.

| a | b |
|---|---|
| 1.5 | 2.5 |

## Next
Out of scope 999.
"""


def test_audit_numbers_synthetic(tmp_path: Path, capsys: pytest.CaptureFixture) -> None:
    doc = tmp_path / "d.md"
    doc.write_text(DOC)
    src = tmp_path / "res.json"
    src.write_text(
        json.dumps(
            {"wrong": 204, "n": 600, "ci": [0.303, 0.379], "rng": [0.0197, 0.0218], "pts": 360000}
        )
    )
    code = audit_numbers.main(["--doc", str(doc), "--section", "Results", "--sources", str(src)])
    out = capsys.readouterr().out
    assert code == 0  # advisory
    unmatched = [line for line in out.splitlines() if "unmatched" in line]
    texts = " ".join(unmatched)
    assert "'41.7%'" in texts and "'55'" in texts
    for ok in ["'34.0%'", "'30.3'", "'37.9'", "'204/600'", "'0.0197'", "'0.0218'", "'360k'"]:
        assert ok not in texts, ok
    for skipped in ["777", "999", "2026", "1234", "'8'"]:
        assert skipped not in texts, skipped
    assert (
        audit_numbers.main(
            ["--doc", str(doc), "--section", "Results", "--sources", str(src), "--strict"]
        )
        == 1
    )


def test_audit_numbers_percent_fraction_and_missing_section(tmp_path: Path) -> None:
    nums = audit_numbers.extract_numbers(["Rate 98-99% and 0.34 ok, ratio 21/50."])
    assert [n["text"] for n in nums] == ["21/50", "98", "99%", "0.34"]
    st = audit_numbers.Store()
    for v in [0.99, 0.98, 21.0, 50.0, 0.34]:
        st.add(v, 0, "k")
    st.files.append("src.json")
    st.seal()
    assert audit_numbers.audit(nums, st) == []
    with pytest.raises(ValueError):
        audit_numbers.find_section(["# A"], "B")


def test_audit_scaling_only_for_percentages_and_fractions(tmp_path: Path) -> None:
    st = audit_numbers.Store()
    st.files.append("s.json")
    for v in [0.34, 3400.0, -12.0]:
        st.add(v, 0, "k")
    st.seal()

    def un(text: str) -> List[str]:
        nums = audit_numbers.extract_numbers([text])
        return [d["text"] for d in audit_numbers.audit(nums, st)]

    assert un("rate 34.0%") == []  # percentage <- source fraction
    assert un("fraction 0.34") == []
    assert un("count 34") == ["34"]  # plain integer: no /100
    assert un("count 3400") == []
    assert un("big 340000") == ["340000"]  # no x100 for plain integers
    assert un("minus -12") == []  # compared on |value|


def test_audit_verbose_and_per_file(tmp_path: Path, capsys: pytest.CaptureFixture) -> None:
    doc = tmp_path / "docs" / "d.md"
    doc.parent.mkdir()
    doc.write_text("# S\n\nRates 34.0% and 12.5% here.\n\n| a |\n|---|\n| 34.0% |\n")
    (tmp_path / "one.json").write_text(json.dumps({"a": {"rate": 0.34}}))
    (tmp_path / "two.txt").write_text("other 12.5 value\n")
    argv = ["--doc", str(doc), "--section", "S", "--sources", str(tmp_path / "*.*")]
    assert audit_numbers.main([*argv, "-v", "--per-file", "--strict"]) == 1
    out = capsys.readouterr().out
    assert "one.json:a.rate" in out and "two.txt:line 1" in out
    assert "no single source file matches all of '34.0%', '12.5%'" in out
    assert "chance-match rate" in out and "dp:" in out  # per precision bucket


# ---------------------------------------------------------------------------------------------
# G. preflight
# ---------------------------------------------------------------------------------------------


def test_preflight_structure_with_mocked_subprocess(monkeypatch: pytest.MonkeyPatch) -> None:
    def fake_run(cmd: List[str]) -> Optional[str]:
        if cmd[0] == "pmset":
            return "Now drawing from 'Battery Power'\n -InternalBattery-0 80%"
        return "  PID  %CPU COMM\n  1 0.5 init\n 4242 95.0 /bin/hog\n"

    monkeypatch.setattr(preflight, "_run", fake_run)
    monkeypatch.setattr(preflight, "SAMPLE_GAP_S", 0.0)
    info = preflight.preflight()
    for key in [
        "timestamp",
        "host",
        "loadavg",
        "cpu_count",
        "power",
        "omp_num_threads",
        "mkl_num_threads",
        "torch_threads",
        "busy_processes",
    ]:
        assert key in info
    assert info["power"] == "battery"
    assert [p["pid"] for p in info["busy_processes"]] == [4242]
    json.dumps(info)  # serialisable


def test_busy_processes_two_samples_and_ignore_list(monkeypatch: pytest.MonkeyPatch) -> None:
    samples = [
        {1: _p(1, 90, "/bin/steady"), 2: _p(2, 90, "/bin/spike"), 3: _p(3, 90, "/x/WindowServer")},
        {1: _p(1, 95, "/bin/steady"), 2: _p(2, 5, "/bin/spike"), 3: _p(3, 99, "/x/WindowServer")},
    ]
    monkeypatch.setattr(preflight, "_sample_cpu", lambda: samples.pop(0))
    monkeypatch.setattr(preflight, "_ancestor_pids", lambda: {os.getpid()})
    monkeypatch.setattr(preflight, "SAMPLE_GAP_S", 0.0)
    assert [p["pid"] for p in preflight._busy_processes()] == [1]
    # the ancestors of this process are never reported
    samples[:] = [{1: _p(1, 90, "uv")}, {1: _p(1, 90, "uv")}]
    monkeypatch.setattr(preflight, "_ancestor_pids", lambda: {1})
    assert preflight._busy_processes() == []


def _p(pid: int, cpu: float, command: str) -> Dict[str, object]:
    return {"pid": pid, "cpu": cpu, "command": command}


def test_require_quiet_logic() -> None:
    quiet = {"loadavg": [0.3, 0.5, 0.5], "power": "ac", "busy_processes": []}
    assert preflight.require_quiet(info=quiet) is quiet
    with pytest.raises(preflight.MachineBusyError, match="load"):
        preflight.require_quiet(info={**quiet, "loadavg": [3.0, 1, 1]})
    with pytest.raises(preflight.MachineBusyError, match="battery"):
        preflight.require_quiet(info={**quiet, "power": "battery"})
    preflight.require_quiet(allow_battery=True, info={**quiet, "power": "battery"})
    busy = [{"pid": 9, "cpu": 80.0, "command": "x"}]
    with pytest.raises(preflight.MachineBusyError, match="busy process"):
        preflight.require_quiet(info={**quiet, "busy_processes": busy})
    preflight.require_quiet(max_load=5.0, info={**quiet, "loadavg": [3.0, 1, 1]})


def test_timing_preflight_cli_writes_json(tmp_path: Path) -> None:
    out = tmp_path / "pf.json"
    res = subprocess.run(
        [sys.executable, str(DEV / "timing_preflight.py"), "--json", str(out)],
        capture_output=True,
        text=True,
    )
    assert res.returncode == 0, res.stderr
    assert "loadavg" in json.loads(out.read_text())


# ---------------------------------------------------------------------------------------------
# H, I. job status, checkpoint, new_todo
# ---------------------------------------------------------------------------------------------


def test_job_status_and_checkpoint(repo: Path) -> None:
    res = _script(repo, "job_status.sh")
    assert res.returncode == 0 and "== git ==" in res.stdout and "branch: develop" in res.stdout
    res = _script(repo, "checkpoint.sh", "a note")
    assert res.returncode == 0, res.stderr
    path = Path(res.stdout.strip())
    assert path.exists() and "Note: a note" in path.read_text()
    assert path.parent.name == "checkpoints"


def test_new_todo(repo: Path) -> None:
    res = _script(repo, "new_todo.sh", "my_idea", "My idea", "sub line")
    assert res.returncode == 0, res.stderr
    path = Path(res.stdout.strip())
    assert path.name == "todo_my_idea.md"
    text = path.read_text()
    assert text.startswith('---\ntitle: "TODO: My idea"\nsubtitle: "sub line"\ndate: "')
    for h in [
        "# Status",
        "# What",
        "# Why it might help",
        "# Framing",
        "# Experiments (sketch)",
        "# Risks",
        "# Relation to other work",
    ]:
        assert h in text
    assert _script(repo, "new_todo.sh", "my_idea", "My idea").returncode == 1  # no overwrite
    assert _script(repo, "new_todo.sh").returncode == 2


def test_agent_files_present() -> None:
    impl = (REPO / ".claude" / "agents" / "implementer.md").read_text()
    assert "model: sonnet" in impl and len(impl.splitlines()) <= 90
    assert "## Claims audit" in (REPO / ".claude" / "agents" / "code-reviewer.md").read_text()


# ---------------------------------------------------------------------------------------------
# G. session tooling: finish_task extensions, checkpoint, agent/skill files
# ---------------------------------------------------------------------------------------------


def test_finish_task_ignores_dirty_claude_md_but_never_stages_it(repo: Path) -> None:
    (repo / "CLAUDE.md").write_text("a\n")
    _git(repo, "add", "CLAUDE.md")
    _git(repo, "commit", "-q", "-m", "claude")
    _task_branch(repo)
    (repo / "CLAUDE.md").write_text("dirty\n")
    res = _script(repo, "finish_task.sh", "--no-tests")
    assert res.returncode == 0, res.stderr
    assert (repo / "CLAUDE.md").read_text() == "dirty\n"
    assert _git(repo, "diff", "--cached", "--name-only") == ""
    assert "CLAUDE.md" not in _git(repo, "show", "--stat", "--format=", "HEAD")


def test_finish_task_removes_task_worktree_pushes_and_starts_next(
    repo: Path, tmp_path: Path
) -> None:
    remote = tmp_path / "remote.git"
    _git(tmp_path, "init", "-q", "--bare", str(remote))
    _git(repo, "remote", "add", "origin", str(remote))
    _task_branch(repo)
    wt = tmp_path / "wt"
    _git(repo, "checkout", "-q", "feature/z")
    _git(repo, "worktree", "add", "-q", str(wt), "feature/z-job")  # task lives in a worktree
    res = _script(wt, "finish_task.sh", "--no-tests", "--yes", "--push", "--next", "t2")
    assert res.returncode == 0, res.stderr
    assert not wt.exists()
    assert "feature/z-job" not in _git(repo, "branch", "--list")
    assert _git(repo, "rev-parse", "--abbrev-ref", "HEAD") == "feature/z-t2"
    assert _git(repo, "config", "branch.feature/z-t2.parent") == "feature/z"
    assert _git(repo, "log", "-1", "--format=%s", "feature/z").startswith("Merge task job:")
    assert _git(remote, "rev-parse", "feature/z") == _git(repo, "rev-parse", "feature/z")


def test_checkpoint_content_memory_pointer_and_length(repo: Path, tmp_path: Path) -> None:
    (repo / "scripts" / "dev").mkdir(parents=True)
    shutil.copy(DEV / "job_status.sh", repo / "scripts" / "dev" / "job_status.sh")
    shutil.copy(DEV / "checkpoint.sh", repo / "scripts" / "dev" / "checkpoint.sh")
    _git(repo, "checkout", "-q", "-b", "feature/c")
    _git(repo, "checkout", "-q", "-b", "feature/c-job")
    _git(repo, "config", "branch.feature/c-job.parent", "feature/c")
    (repo / "k.txt").write_text("k\n")
    plans = tmp_path / "plans"
    plans.mkdir()
    (plans / "my-plan.md").write_text("plan\n")
    (repo / ".claude" / "reports").mkdir(parents=True)
    (repo / ".claude" / "reports" / "20260101-0000-t.md").write_text("r\n")
    mem = tmp_path / "mem.md"
    mem.write_text("# mem\n- Latest checkpoint: /old/path\n")
    env = {**os.environ, **GIT_ENV, "CHECKPOINT_PLANS_DIR": str(plans)}
    cmd = [str(repo / "scripts/dev/checkpoint.sh"), "--memory-file", str(mem), "a note"]
    res = subprocess.run(cmd, cwd=repo, env=env, capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    out = Path(res.stdout.strip())
    text = out.read_text()
    assert out.parent.resolve() == (repo / ".claude" / "checkpoints").resolve()
    for needle in [
        "- branch: feature/c-job",
        "- parent: feature/c",
        "k.txt",
        "my-plan.md",
        ".claude/reports/20260101-0000-t.md",
        "Note: a note",
        "## Done this session",
        "## Waiting on the owner",
        "## Next step",
    ]:
        assert needle in text, needle
    assert len(text.splitlines()) < 80
    assert mem.read_text().count("Latest checkpoint:") == 1
    assert str(out) in mem.read_text() and "/old/path" not in mem.read_text()


def test_agent_and_skill_files_are_consistent() -> None:
    agents = REPO / ".claude" / "agents"
    assert (agents / "study-conventions.md").exists()
    for name in ["implementer.md", "code-reviewer.md"]:
        body = (agents / name).read_text()
        assert "study-conventions.md" in body and ".claude/reports/" in body
    for skill in ["handoff", "resume", "review-task", "merge-task"]:
        text = (REPO / ".claude" / "skills" / skill / "SKILL.md").read_text()
        assert f"name: {skill}" in text and "user-invocable: true" in text
    assert ".claude/reports/" in (REPO / ".gitignore").read_text().splitlines()
