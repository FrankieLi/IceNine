"""Tests of scripts/dev/verify_handback.sh and scripts/dev/read_on_fetch.sh.

Every test works in temporary repositories (a bare repository stands in for origin; no network) and
never reads or writes the real clone's git config. shmem and pytest are replaced by stubs.
"""

import os
import stat
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional

import pytest

REPO = Path(__file__).resolve().parents[2]
DEV = REPO / "scripts" / "dev"
GIT_ENV = {
    "GIT_AUTHOR_NAME": "t",
    "GIT_AUTHOR_EMAIL": "t@example.com",
    "GIT_COMMITTER_NAME": "t",
    "GIT_COMMITTER_EMAIL": "t@example.com",
    "GIT_CONFIG_GLOBAL": os.devnull,
    "GIT_CONFIG_NOSYSTEM": "1",
}
ROADMAP = "icenine_py/docs/ROADMAP.md"
DECISIONS = "icenine_py/docs/DECISIONS.md"
ILIUM_BRANCH = "feature/python-speed-parity"
NEWT_BRANCH = "feature/eval-reduction"

ROADMAP_TEXT = """# Between agents

## From Newt Agent

- 2026-10-14 heads-up: newt first

## From Ilium Agent

- 2026-10-13 answer: ilium first

# Tracks

| E | cost | `feature/eval-reduction` | Newt Agent |

# Track E: eval

| E1 | newt row | open |

# Track S: cost

| S1 | ilium row | open |
| S2 | ilium row 2 | open |

# Track H: housekeeping

| H11 | other | open |
| H12 | tooling | in progress |
"""

DECISIONS_TEXT = """# Decision log

### D-20261014-I1 · S · Ilium decision one
- 2026-10-14, owner. Recorded by Ilium Agent. Supersedes: none.
- Ilium text. Where: S1.
"""


def _env(extra: Optional[Dict[str, str]] = None) -> Dict[str, str]:
    env = {k: v for k, v in os.environ.items() if not k.startswith(("VERIFY_", "READ_ON_"))}
    env.update(GIT_ENV)
    env.update(extra or {})
    return env


def _git(repo: Path, *args: str) -> str:
    out = subprocess.run(["git", *args], cwd=repo, env=_env(), capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    return out.stdout.strip()


def _cfg(repo: Path, pattern: str) -> str:
    """Local config entries matching a pattern ("" when there is none)."""
    out = subprocess.run(
        ["git", "config", "--local", "--get-regexp", pattern],
        cwd=repo,
        env=_env(),
        capture_output=True,
        text=True,
    )
    return out.stdout.strip()


def _commit(repo: Path, files: Dict[str, str], msg: str = "change") -> str:
    for rel, text in files.items():
        p = repo / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text)
        _git(repo, "add", rel)
    _git(repo, "commit", "-q", "-m", msg)
    return _git(repo, "rev-parse", "HEAD")


def _run(
    script: str, cwd: Path, *args: str, env: Optional[Dict[str, str]] = None
) -> subprocess.CompletedProcess:
    return subprocess.run(
        [str(DEV / script), *args], cwd=cwd, env=_env(env), capture_output=True, text=True
    )


def _executable(path: Path, text: str) -> Path:
    path.write_text(text)
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return path


# ---------------------------------------------------------------------------------------------
# read_on_fetch.sh
# ---------------------------------------------------------------------------------------------

SHMEM_STUB = """#!/bin/bash
# stub: STUB_KNOWN lists the IDs that "are saved"; every call is logged to STUB_LOG
printf '%s\\n' "$*" >> "$STUB_LOG"
case "$1" in
  query)
    [ "${STUB_QUERY_RC:-0}" = 0 ] || exit "$STUB_QUERY_RC"
    for id in ${STUB_KNOWN:-}; do
      if [ "$2" = "Owner decision $id" ]; then
        echo "[f1] (user, 2026-10-18 10:00)  Owner decision $id (2026-10-14, IceNine E): saved"
      fi
    done
    echo "[f2] (tool, 2026-10-18 10:00)  [Bash] noise that mentions D-20261014-N2 and D-20261014-N1"
    ;;
  add) exit "${STUB_ADD_RC:-0}" ;;
esac
"""


class World:
    """A bare origin, the caller's clone `own` and a peer clone that writes the other branch."""

    def __init__(self, root: Path, own_branch: str, peer_branch: str) -> None:
        self.root = root
        self.origin = root / "origin.git"
        subprocess.run(["git", "init", "-q", "--bare", str(self.origin)], check=True)
        self.own = root / "own"
        self.own.mkdir()
        _git(self.own, "init", "-q", "-b", own_branch)
        _commit(
            self.own,
            {ROADMAP: ROADMAP_TEXT, DECISIONS: DECISIONS_TEXT, "README": "x\n"},
            "base",
        )
        _git(self.own, "remote", "add", "origin", str(self.origin))
        _git(self.own, "branch", peer_branch)
        _git(self.own, "push", "-q", "origin", own_branch, peer_branch)
        self.peer = root / "peer"
        subprocess.run(
            ["git", "clone", "-q", str(self.origin), str(self.peer)], check=True, env=_env()
        )
        _git(self.peer, "checkout", "-q", peer_branch)
        self.peer_branch = peer_branch
        self.stub = _executable(root / "shmem_stub.sh", SHMEM_STUB)
        self.log = root / "shmem.log"

    def peer_push(self, files: Dict[str, str], msg: str) -> str:
        h = _commit(self.peer, files, msg)
        _git(self.peer, "push", "-q", "origin", self.peer_branch)
        return h

    def env(self, **extra: str) -> Dict[str, str]:
        env = {"READ_ON_FETCH_SHMEM_CMD": str(self.stub), "STUB_LOG": str(self.log)}
        env.update(extra)
        return env

    def shmem_calls(self) -> List[str]:
        return self.log.read_text().splitlines() if self.log.exists() else []


def _newt_files(mailbox_line: str, decisions: str) -> Dict[str, str]:
    return {
        ROADMAP: ROADMAP_TEXT.replace(
            "- 2026-10-14 heads-up: newt first\n",
            "- 2026-10-14 heads-up: newt first\n" + mailbox_line + "\n",
        ),
        DECISIONS: DECISIONS_TEXT + decisions,
        "icenine_py/docs/newt_note.md": "note\n",
    }


N_DECISIONS = """
### D-20261014-N1 · E · First Newt decision
- 2026-10-14, owner. Recorded by Newt Agent. Supersedes: none.
- Newt text one. Where: E1.

### D-20261014-N2 · E (O3) · Second Newt decision with a 'quote' and `code`
- 2026-10-15, owner. Recorded by Newt Agent. Supersedes: none.
- Newt text two, with 5 $HOME "dq" signs. Where: E2.
"""


@pytest.fixture()
def world(tmp_path: Path) -> World:
    w = World(tmp_path, ILIUM_BRANCH, NEWT_BRANCH)
    w.peer_push(_newt_files("- 2026-10-15 request: newt second", N_DECISIONS), "newt work")
    return w


def test_read_on_fetch_refuses_without_agent(world: World) -> None:
    res = _run("read_on_fetch.sh", world.own, env=world.env())
    assert res.returncode == 1 and "no agent letter" in res.stderr
    assert world.shmem_calls() == []
    res = _run("read_on_fetch.sh", world.own, "--agent", "X", env=world.env())
    assert res.returncode == 1 and "must be I or N" in res.stderr
    assert _cfg(world.own, "icenine") == ""


def test_read_on_fetch_ilium_reads_newt_branch(world: World) -> None:
    _git(world.own, "config", "--local", "icenine.agent", "I")
    res = _run("read_on_fetch.sh", world.own, env=world.env(STUB_KNOWN="D-20261014-N1"))
    assert res.returncode == 0, res.stderr
    out = res.stdout
    assert "agent I; other branch origin/feature/eval-reduction" in out
    assert "newt work" in out and "newt_note.md" in out
    assert "+- 2026-10-15 request: newt second" in out  # mailbox diff
    assert "D-20261014-N2" in out and "not in shmem: D-20261014-N2" in out
    assert "not in shmem: D-20261014-N1" not in out  # saved: found by the exact phrase
    assert "D-20261014-I1" not in out  # own letter is ignored
    assert "git merge origin/feature/eval-reduction" in out
    # printed, not run: only queries reached the stub
    assert all(c.startswith("query ") for c in world.shmem_calls())
    assert not any(c.startswith("add") for c in world.shmem_calls())
    # the ready line is a valid shell line that reproduces the fact exactly
    ready = [ln for ln in out.splitlines() if " add --source user " in ln]
    assert len(ready) == 1
    run_ready = subprocess.run(
        ["bash", "-c", ready[0]], env=_env(world.env()), capture_output=True, text=True
    )
    assert run_ready.returncode == 0
    added = world.shmem_calls()[-1]
    assert added == (
        "add --source user Owner decision D-20261014-N2 (2026-10-15, IceNine E (O3)): "
        "Second Newt decision with a 'quote' and `code`. Newt text two, with 5 $HOME \"dq\" "
        "signs. Where: E2."
    )


def test_read_on_fetch_save_runs_the_add_lines(world: World) -> None:
    res = _run(
        "read_on_fetch.sh",
        world.own,
        "--agent",
        "I",
        "--save",
        env=world.env(STUB_KNOWN="D-20261014-N1"),
    )
    assert res.returncode == 0, res.stderr
    assert "saved D-20261014-N2" in res.stdout and "ready lines" not in res.stdout
    adds = [c for c in world.shmem_calls() if c.startswith("add ")]
    assert len(adds) == 1 and "D-20261014-N2" in adds[0]
    # a failing add is reported and the read point does not move
    res = _run(
        "read_on_fetch.sh",
        world.own,
        "--agent",
        "I",
        "--save",
        "--mark",
        env=world.env(STUB_ADD_RC="1"),
    )
    assert res.returncode == 1 and "FAILED to save D-20261014-N1" in res.stdout
    assert _cfg(world.own, "readonfetch") == ""


def test_read_on_fetch_mark_moves_the_point(world: World) -> None:
    key = "icenine.feature/eval-reduction.readOnFetch"
    tip = _git(world.peer, "rev-parse", "HEAD")
    first = _run("read_on_fetch.sh", world.own, "--agent", "I", env=world.env())
    assert "newt work" in first.stdout
    assert _cfg(world.own, "readonfetch") == ""  # no --mark
    marked = _run("read_on_fetch.sh", world.own, "--agent", "I", "--mark", env=world.env())
    assert marked.returncode == 0 and "newt work" in marked.stdout
    assert _git(world.own, "config", "--local", "--get", key) == tip
    again = _run("read_on_fetch.sh", world.own, "--agent", "I", env=world.env())
    assert again.returncode == 0
    assert "newt work" not in again.stdout and "none" in again.stdout
    assert "unchanged" in again.stdout
    # a new peer commit shows alone; --since overrides the stored point
    world.peer_push({"icenine_py/docs/more.md": "m\n"}, "newt more")
    third = _run("read_on_fetch.sh", world.own, "--agent", "I", env=world.env())
    assert "newt more" in third.stdout and "newt work" not in third.stdout
    base = _git(world.own, "rev-list", "--max-parents=0", "HEAD")
    since = _run("read_on_fetch.sh", world.own, "--agent", "I", "--since", base, env=world.env())
    assert "newt more" in since.stdout and "newt work" in since.stdout


def test_read_on_fetch_leaves_out_own_commits_the_peer_merged(world: World) -> None:
    _commit(world.own, {"icenine_py/docs/ilium_note.md": "i\n"}, "ilium own work")
    _git(world.own, "push", "-q", "origin", ILIUM_BRANCH)
    _git(world.peer, "fetch", "-q", "origin")
    _git(world.peer, "merge", "-q", "--no-edit", f"origin/{ILIUM_BRANCH}")
    _git(world.peer, "push", "-q", "origin", NEWT_BRANCH)
    base = _git(world.own, "rev-list", "--max-parents=0", "HEAD")
    res = _run("read_on_fetch.sh", world.own, "--agent", "I", "--since", base, env=world.env())
    assert res.returncode == 0, res.stderr
    commits = res.stdout.split("== commits on", 1)[1].split("== mailbox", 1)[0]
    assert "newt work" in commits
    assert "ilium own work" not in commits
    assert "(1 commit(s) already on your branch left out)" in commits
    assert "newt_note.md" in commits and "ilium_note.md" not in commits
    # once the peer branch is merged here too, nothing is new
    _git(world.own, "merge", "-q", "--no-edit", f"origin/{NEWT_BRANCH}")
    again = _run("read_on_fetch.sh", world.own, "--agent", "I", "--since", base, env=world.env())
    commits = again.stdout.split("== commits on", 1)[1].split("== mailbox", 1)[0]
    assert "none" in commits and "already on your branch)" in commits


def test_read_on_fetch_unreachable_shmem_is_reported(world: World) -> None:
    res = _run(
        "read_on_fetch.sh",
        world.own,
        "--agent",
        "I",
        "--mark",
        env=world.env(STUB_QUERY_RC="7"),
    )
    assert res.returncode == 3
    assert "shmem unreachable" in res.stdout
    assert "not in shmem" not in res.stdout and "ready lines" not in res.stdout
    assert _cfg(world.own, "readonfetch") == ""
    missing = _run(
        "read_on_fetch.sh",
        world.own,
        "--agent",
        "I",
        env={"READ_ON_FETCH_SHMEM_CMD": str(world.root / "no_such_shmem")},
    )
    assert missing.returncode == 3 and "shmem unreachable" in missing.stdout


def test_read_on_fetch_newt_reads_ilium_branch(tmp_path: Path) -> None:
    w = World(tmp_path, NEWT_BRANCH, ILIUM_BRANCH)
    _git(w.own, "config", "--local", "icenine.agent", "N")
    files = {
        ROADMAP: ROADMAP_TEXT.replace(
            "- 2026-10-13 answer: ilium first\n",
            "- 2026-10-13 answer: ilium first\n- 2026-10-16 heads-up: ilium second\n",
        ),
        DECISIONS: DECISIONS_TEXT
        + "\n### D-20261016-I2 · S · Ilium two\n- 2026-10-16, owner. Recorded by Ilium Agent."
        " Supersedes: none.\n- Ilium text two. Where: S2.\n"
        + "\n### D-20261016-N1 · E · A Newt entry\n- 2026-10-16, owner. Recorded by Newt Agent."
        " Supersedes: none.\n- Own text. Where: E1.\n",
    }
    w.peer_push(files, "ilium work")
    res = _run("read_on_fetch.sh", w.own, env=w.env(STUB_KNOWN="D-20261014-I1"))
    assert res.returncode == 0, res.stderr
    out = res.stdout
    assert "agent N; other branch origin/feature/python-speed-parity" in out
    assert "ilium work" in out and "+- 2026-10-16 heads-up: ilium second" in out
    assert "not in shmem: D-20261016-I2" in out
    assert "D-20261014-I1" not in out  # saved in the stub, so not listed
    assert "D-20261016-N1" not in out  # Newt's own entry is not read as the other agent's
    assert "git merge origin/feature/python-speed-parity" in out
    assert not any(c.startswith("add") for c in w.shmem_calls())


def test_read_on_fetch_other_override(world: World) -> None:
    _git(world.own, "checkout", "-q", "-b", "feature/side", ILIUM_BRANCH)
    _git(world.own, "push", "-q", "origin", "feature/side")
    _git(world.peer, "fetch", "-q", "origin")
    _git(world.peer, "checkout", "-q", "-b", "feature/side", "origin/feature/side")
    _commit(world.peer, {"side.txt": "s\n"}, "side commit")
    _git(world.peer, "push", "-q", "origin", "feature/side")
    res = _run(
        "read_on_fetch.sh", world.own, "--agent", "I", "--other", "feature/side", env=world.env()
    )
    assert res.returncode == 0, res.stderr
    assert "origin/feature/side" in res.stdout and "side commit" in res.stdout
    assert "newt work" not in res.stdout
    bad = _run("read_on_fetch.sh", world.own, "--agent", "I", "--other", "nope", env=world.env())
    assert bad.returncode == 1 and "does not exist" in bad.stderr


# ---------------------------------------------------------------------------------------------
# verify_handback.sh
# ---------------------------------------------------------------------------------------------

REPORT_FIELDS = [
    "Branch / commits",
    "Files",
    "Tests",
    "Results",
    "Criteria",
    "Deviations from the plan",
    "Suspicious / unexpected",
    "Left undone",
]


def _report(
    path: Path,
    hashes: List[str],
    *,
    drop: Optional[str] = None,
    tests: str = "tests/test_x.py::test_a, tests/test_y.py; full suite deferred",
    results: str = "none",
    extra: str = "",
) -> Path:
    bodies = {
        "Branch / commits": "branch t, " + " ".join(h[:9] for h in hashes),
        "Files": "scripts/dev/x.sh",
        "Tests": tests,
        "Results": results,
        "Criteria": "met",
        "Deviations from the plan": "none",
        "Suspicious / unexpected": "none",
        "Left undone": "none",
    }
    lines = [f"{k}: {v}" for k, v in bodies.items() if k != drop]
    path.write_text("\n".join(lines) + "\n" + extra)
    return path


class Task:
    def __init__(self, root: Path) -> None:
        self.root = root
        self.repo = root / "r"
        self.repo.mkdir()
        _git(self.repo, "init", "-q", "-b", ILIUM_BRANCH)
        _commit(
            self.repo,
            {ROADMAP: ROADMAP_TEXT, DECISIONS: DECISIONS_TEXT, "icenine_py/keep.txt": "k\n"},
            "base",
        )
        self.parent = ILIUM_BRANCH
        self.branch = ILIUM_BRANCH + "-t"
        _git(self.repo, "checkout", "-q", "-b", self.branch)
        _git(self.repo, "config", "branch." + self.branch + ".parent", self.parent)
        self.args_file = root / "pytest_args.txt"
        self.hashes: List[str] = []

    def commit(self, files: Dict[str, str], msg: str = "work") -> str:
        h = _commit(self.repo, files, msg)
        self.hashes.append(h)
        return h

    def env(self, rc: int = 0, out: str = "3 passed in 0.1s", **extra: str) -> Dict[str, str]:
        stub = _executable(
            self.root / "pytest_stub.sh",
            f'#!/bin/bash\necho "$@" > "{self.args_file}"\npwd >> "{self.args_file}"\n'
            f"echo '{out}'\nexit {rc}\n",
        )
        return {"VERIFY_HANDBACK_PYTEST_CMD": str(stub), **extra}

    def report(self, **kw: Any) -> Path:
        return _report(self.root / "report.md", self.hashes, **kw)

    def run(
        self,
        *args: str,
        env: Optional[Dict[str, str]] = None,
        agent: Optional[str] = "I",
        cwd: Optional[Path] = None,
        report: Optional[str] = None,
    ) -> subprocess.CompletedProcess:
        if agent is not None:
            _git(self.repo, "config", "--local", "icenine.agent", agent)
        return _run(
            "verify_handback.sh",
            cwd or self.repo,
            self.parent,
            report if report is not None else str(self.root / "report.md"),
            *args,
            env=env if env is not None else self.env(),
        )


@pytest.fixture()
def task(tmp_path: Path) -> Task:
    return Task(tmp_path)


def test_verify_clean_handback_passes_and_runs_named_tests(task: Task) -> None:
    task.commit({"scripts/dev/x.sh": "echo\n"})
    task.report()
    res = task.run()
    assert res.returncode == 0, res.stdout
    assert "summary: 0 FAIL, 0 FLAG" in res.stdout
    assert "all eight report fields present" in res.stdout
    assert "commit hashes agree" in res.stdout
    assert "named tests passed (2 path(s))" in res.stdout
    assert len(res.stdout.splitlines()) <= 30
    called = task.args_file.read_text().splitlines()
    assert called[0] == "tests/test_x.py::test_a tests/test_y.py"
    assert called[1].endswith("/icenine_py")


def test_verify_failing_test_is_fail(task: Task) -> None:
    task.commit({"scripts/dev/x.sh": "echo\n"})
    task.report()
    res = task.run(env=task.env(rc=1, out="1 failed, 2 passed"))
    assert res.returncode == 1
    assert "FAIL  named tests failed" in res.stdout and "1 failed, 2 passed" in res.stdout
    assert "summary: 1 FAIL" in res.stdout


def test_verify_tests_flag_skip_and_no_tests(task: Task) -> None:
    task.commit({"scripts/dev/x.sh": "echo\n"})
    task.report(tests="none named")
    res = task.run()
    assert res.returncode == 0 and "FLAG  report names no tests/" in res.stdout
    task.report()
    res = task.run("--no-tests")
    assert "tests skipped" in res.stdout and not task.args_file.exists()


def test_verify_bare_test_names_resolve_under_icenine_py_tests(task: Task) -> None:
    task.commit({"icenine_py/tests/test_here.py": "x\n", "icenine_py/tests/test_two.py": "x\n"})
    task.report(tests="`test_here.py::test_a`, test_two.py; also tests/test_two.py. Full suite: 5")
    res = task.run()
    assert res.returncode == 0, res.stdout
    assert "named tests passed (2 path(s))" in res.stdout
    assert "summary: 0 FAIL, 0 FLAG" in res.stdout
    called = task.args_file.read_text().splitlines()[0]
    assert called == "tests/test_here.py::test_a tests/test_two.py"


def test_verify_missing_bare_test_name_is_flagged_and_not_run(task: Task) -> None:
    task.commit({"icenine_py/tests/test_here.py": "x\n"})
    task.report(tests="test_here.py and test_gone.py::test_b")
    res = task.run()
    assert res.returncode == 0, res.stdout
    assert "FLAG  test file(s) named on the Tests field not found" in res.stdout
    assert "test_gone.py::test_b" in res.stdout
    assert task.args_file.read_text().splitlines()[0] == "tests/test_here.py"
    # nothing left to run: flagged, no pytest call
    task.args_file.unlink()
    task.report(tests="test_gone.py")
    res = task.run()
    assert "FLAG  test file(s)" in res.stdout and "tests not run" in res.stdout
    assert not task.args_file.exists()


def test_verify_no_worktree_flags_and_skips_tests(task: Task) -> None:
    task.commit({"scripts/dev/x.sh": "echo\n"})
    task.report()
    other = task.root / "other_wt"
    _git(task.repo, "checkout", "-q", "--detach")  # nobody holds the task branch now
    _git(task.repo, "worktree", "add", "-q", "-f", str(other), "--detach")
    res = _run(
        "verify_handback.sh",
        other,
        task.parent,
        str(task.root / "report.md"),
        "--branch",
        task.branch,
        "--agent",
        "I",
        env=task.env(),
    )
    assert res.returncode == 0
    assert "no worktree holds" in res.stdout and not task.args_file.exists()
    # the same branch checked out in a worktree is found from another checkout
    _git(task.repo, "checkout", "-q", task.branch)
    _git(task.repo, "checkout", "-q", "--detach")
    wt = task.root / "task_wt"
    _git(task.repo, "worktree", "add", "-q", str(wt), task.branch)
    res = _run(
        "verify_handback.sh",
        other,
        task.parent,
        str(task.root / "report.md"),
        "--branch",
        task.branch,
        "--agent",
        "I",
        env=task.env(),
    )
    assert res.returncode == 0 and "named tests passed" in res.stdout
    assert task.args_file.read_text().splitlines()[1].endswith("task_wt/icenine_py")


def test_verify_no_commits_is_fail(task: Task) -> None:
    task.report()
    res = task.run("--no-tests")
    assert res.returncode == 1 and "FAIL  no commits on" in res.stdout


def test_verify_changed_file_classes(task: Task) -> None:
    task.commit(
        {
            "icenine_py/icenine/cost_functions.py": "a\n",
            "icenine_py/icenine/other_module.py": "b\n",
            "Src/x.cpp": "c\n",
            "icenine_py/pyproject.toml": "d\n",
            "icenine_py/tests/data/x_golden.npz": "e\n",
        }
    )
    task.report()
    res = task.run("--no-tests")
    out = res.stdout
    assert res.returncode == 0, out
    assert "FLAG  shared Track E/S file changed: icenine_py/icenine/cost_functions.py" in out
    assert "heads-up is due" in out
    assert "FLAG  library change: icenine_py/icenine/other_module.py" in out
    assert "cost_functions.py" not in out.split("FLAG  library change")[1].split("\n")[0]
    assert "FLAG  C++ change: Src/x.cpp" in out
    assert "FLAG  dependency files: icenine_py/pyproject.toml" in out
    assert "FLAG  golden file: icenine_py/tests/data/x_golden.npz" in out
    assert "summary: 0 FAIL, 5 FLAG" in out


def test_verify_committed_forbidden_files_are_fail(task: Task) -> None:
    task.commit(
        {
            ".claude/agents/a.md": "a\n",
            "CLAUDE.md": "c\n",
            DECISIONS: DECISIONS_TEXT + "\n### D-20261018-I9 · X · y\n",
        }
    )
    task.report()
    res = task.run("--no-tests")
    assert res.returncode == 1
    for needle in [
        "FAIL  DECISIONS.md committed",
        "FAIL  CLAUDE.md committed",
        "FAIL  .claude/** committed: .claude/agents/a.md",
    ]:
        assert needle in res.stdout
    assert "summary: 3 FAIL" in res.stdout


def test_verify_missing_report_field_is_fail(task: Task) -> None:
    task.commit({"scripts/dev/x.sh": "echo\n"})
    task.report(drop="Suspicious / unexpected")
    res = task.run("--no-tests")
    assert res.returncode == 1
    assert "FAIL  report fields missing: suspicious / unexpected" in res.stdout
    # the other seven are all accepted in heading and bold forms
    md = "\n".join(
        f"## {f}\ntext" if i % 2 else f"**{f}**: text" for i, f in enumerate(REPORT_FIELDS)
    )
    (task.root / "report.md").write_text(md + "\n")
    res = task.run("--no-tests")
    assert "all eight report fields present" in res.stdout


FIELD_FORMS = [
    # (form with {} for the field name, accepted)
    ("{}:", True),
    ("{}", True),
    ("**{}:**", True),
    ("**{}**:", True),
    ("## {}", True),
    ("- {}:", True),
    ("- **{}:** x", True),
    ("{} (what the tests cover):", True),
    ("{} (x)", True),
    ("**{}** (x):", True),
    ("### {} (what the tests cover)", True),
    ("4. {}:", True),
    ("4) {}:", True),
    ("{} (plan Steps 1 and 2, owner refinements):", True),
    ("## {} — x", True),
    ("{}foo:", False),
    ("{} — x", False),
    ("{} of the run", False),
]


@pytest.mark.parametrize("form,accepted", FIELD_FORMS)
@pytest.mark.parametrize("name", ["Results", "RESULTS", "results"])
def test_verify_report_field_forms(task: Task, name: str, form: str, accepted: bool) -> None:
    task.commit({"a.txt": "1\n"})
    task.report(drop="Results")
    with open(task.root / "report.md", "a") as fh:
        fh.write(form.format(name) + "\n")
    res = task.run("--no-tests")
    assert ("report fields missing: results" in res.stdout) is (not accepted), res.stdout


def test_verify_report_field_lookalikes_and_slashes(task: Task) -> None:
    task.commit({"a.txt": "1\n"})
    task.report(drop="Results")
    with open(task.root / "report.md", "a") as fh:
        fh.write("Result:\n")  # singular: not the field
    assert "report fields missing: results" in task.run("--no-tests").stdout
    lines = [
        "Branch/commits: b " + task.hashes[0][:9],
        "Files (list):",
        "### Tests (targeted)",
        "**Results** (x):",
        "3. Criteria (plan Steps 3 and 4):",
        "Deviations from the plan (none):",
        "Suspicious/unexpected: none",
        "Left undone (none):",
    ]
    (task.root / "report.md").write_text("\n".join(lines) + "\n")
    res = task.run("--no-tests")
    assert "all eight report fields present" in res.stdout, res.stdout
    assert "commit hashes agree" in res.stdout


def test_verify_relative_paths_from_a_subdirectory(task: Task) -> None:
    task.commit({"a.txt": "1\n"})
    task.report()
    review = task.root / "review.md"
    review.write_text("## Required fixes\n1. one\n")
    sub = task.repo / "icenine_py"
    sub.mkdir(exist_ok=True)
    rel = os.path.relpath(task.root / "report.md", sub)
    res = task.run("--no-tests", "--fixes", os.path.relpath(review, sub), cwd=sub, report=rel)
    assert res.returncode in (0, 1) and "not found" not in res.stderr, res.stderr
    assert "report fields missing" not in res.stdout and "summary:" in res.stdout
    assert "no numbered items" not in res.stdout and "required fixes not named" in res.stdout


def test_verify_hash_cross_check_both_ways(task: Task) -> None:
    h1 = task.commit({"a.txt": "1\n"})
    h2 = task.commit({"b.txt": "2\n"})
    parent_hash = _git(task.repo, "rev-parse", task.parent)
    _report(
        task.root / "report.md",
        [h1],
        extra="",
    )
    res = task.run("--no-tests")
    assert f"branch commits not named in the report: {h2[:9]}" in res.stdout
    text = (
        (task.root / "report.md")
        .read_text()
        .replace(h1[:9], h1[:9] + " abcdef1 " + parent_hash[:9])
    )
    (task.root / "report.md").write_text(text)
    res = task.run("--no-tests")
    assert "report names hashes not on the branch" in res.stdout
    assert "abcdef1" in res.stdout and parent_hash[:9] in res.stdout
    assert res.returncode == 0  # FLAG only
    # a report with no hashes at all
    (task.root / "report.md").write_text(
        "Branch / commits: none\n" + "\n".join(f"{f}: x" for f in REPORT_FIELDS[1:]) + "\n"
    )
    res = task.run("--no-tests")
    assert "names no commit hash" in res.stdout


def test_verify_worktree_changes(task: Task) -> None:
    task.commit({"a.txt": "1\n"})
    task.report()
    (task.repo / "stray.txt").write_text("s\n")
    (task.repo / ".claude").mkdir()
    (task.repo / ".claude" / "skill.md").write_text("k\n")
    (task.repo / "CLAUDE.md").write_text("c\n")
    res = task.run("--no-tests")
    assert "FLAG  uncommitted changes in the worktree: stray.txt" in res.stdout
    assert "skill.md" not in res.stdout.split("note  left uncommitted")[0]
    assert "note  left uncommitted for the owner: .claude/, CLAUDE.md" in res.stdout
    (task.repo / "stray.txt").unlink()
    res = task.run("--no-tests")
    assert "ok    worktree clean outside CLAUDE.md and .claude/**" in res.stdout


@pytest.mark.parametrize(
    "agent,old,new,flagged",
    [
        ("I", "| S1 | ilium row | open |", "| S1 | ilium row | **done** |", False),
        ("I", "| H12 | tooling | in progress |", "| H12 | tooling | merged |", False),
        ("I", "| H11 | other | open |", "| H11 | other | closed |", True),
        ("I", "| E1 | newt row | open |", "| E1 | newt row | closed |", True),
        (
            "I",
            "- 2026-10-13 answer: ilium first\n",
            "- 2026-10-13 answer: ilium first\n- new\n",
            False,
        ),
        (
            "I",
            "- 2026-10-14 heads-up: newt first\n",
            "- 2026-10-14 heads-up: newt first\n- x\n",
            True,
        ),
        ("N", "| E1 | newt row | open |", "| E1 | newt row | closed |", False),
        (
            "N",
            "- 2026-10-14 heads-up: newt first\n",
            "- 2026-10-14 heads-up: newt first\n- x\n",
            False,
        ),
        ("N", "| S1 | ilium row | open |", "| S1 | ilium row | **done** |", True),
        ("N", "| H12 | tooling | in progress |", "| H12 | tooling | merged |", True),
    ],
)
def test_verify_roadmap_territory(
    task: Task, agent: str, old: str, new: str, flagged: bool
) -> None:
    assert old in ROADMAP_TEXT
    task.commit({ROADMAP: ROADMAP_TEXT.replace(old, new)})
    task.report()
    res = task.run("--no-tests", agent=agent)
    assert res.returncode == 0, res.stdout
    assert ("FLAG  ROADMAP.md changed outside" in res.stdout) is flagged
    if not flagged:
        assert "ok    ROADMAP.md changes are inside" in res.stdout


def test_verify_refuses_without_agent_letter(task: Task) -> None:
    task.commit({ROADMAP: ROADMAP_TEXT.replace("| S1 | ilium row | open |", "| S1 | x | open |")})
    task.report()
    res = task.run("--no-tests", agent=None)
    assert res.returncode == 2 and "no agent letter" in res.stderr
    assert "summary" not in res.stdout and not task.args_file.exists()
    # --agent works instead of the config key
    res = task.run("--no-tests", "--agent", "I", agent=None)
    assert res.returncode == 0 and "summary: 0 FAIL, 0 FLAG" in res.stdout


def test_verify_roadmap_hunk_spanning_two_territories_is_flagged(task: Task) -> None:
    text = ROADMAP_TEXT.replace("| S2 | ilium row 2 | open |", "| S2 | x | open |").replace(
        "| H11 | other | open |", "| H11 | y | open |"
    )
    task.commit({ROADMAP: text})
    task.report()
    res = task.run("--no-tests")
    assert "FLAG  ROADMAP.md changed outside" in res.stdout


def test_verify_timing_label_heuristic(task: Task) -> None:
    task.commit({"a.txt": "1\n"})
    task.report(results="\n- evaluation 12 ms per voxel\n- screen 3.5 us\n")
    res = task.run("--no-tests")
    assert "FLAG  2 Results line(s) with a time unit and no" in res.stdout
    task.report(
        results="\n- evaluation 12 ms per voxel (single-worker)\n- screen 3.5 us contended\n"
    )
    res = task.run("--no-tests")
    assert "no unlabelled timing line" in res.stdout and "FLAG  2 Results" not in res.stdout
    task.report(results="single-worker timings\n- evaluation 12 ms per voxel\n")
    assert "no unlabelled timing line" in task.run("--no-tests").stdout
    task.report(results="\n- 5 samples, 3 settings, 12 voxels\n")
    assert "no unlabelled timing line" in task.run("--no-tests").stdout


def test_verify_long_report_is_flagged(task: Task) -> None:
    task.commit({"a.txt": "1\n"})
    task.report(extra="\n".join(f"line {i}" for i in range(80)) + "\n")
    res = task.run("--no-tests")
    assert "FLAG  report has 88 lines" in res.stdout


def test_verify_fixes_numbering(task: Task) -> None:
    task.commit({"a.txt": "1\n"})
    review = task.root / "review.md"
    review.write_text(
        "# Review\n\n1. not in the list\n\n## Required fixes\n\n"
        "1. first thing\n2. second thing\n3) third thing\n4. fourth thing\n\n## Other\n\n9. nope\n"
    )
    task.report(extra="Fix 1: done. Fixes 3-4 done (abc1234 fix 2 not)\n")
    res = task.run("--no-tests", "--fixes", str(review))
    assert "all 4 required fixes named" in res.stdout
    task.report(extra="Fix 1 and fix 4 done.\n")
    res = task.run("--no-tests", "--fixes", str(review))
    assert "FLAG  required fixes not named in the report: 2, 3" in res.stdout
    task.report(extra="")
    (task.root / "report.md").write_text(
        (task.root / "report.md").read_text().replace("Files: scripts/dev/x.sh", "Files: x")
        + "\n1. one\n2. two\n3. three\n4. four\n"
    )
    assert "all 4 required fixes named" in task.run("--no-tests", "--fixes", str(review)).stdout
    review.write_text("# Review\n\nnothing numbered\n")
    res = task.run("--no-tests", "--fixes", str(review))
    assert "no numbered items under" in res.stdout


def test_verify_argument_errors(task: Task) -> None:
    task.commit({"a.txt": "1\n"})
    task.report()
    res = _run("verify_handback.sh", task.repo)
    assert res.returncode == 2 and "usage" in res.stderr
    res = _run("verify_handback.sh", task.repo, "nope", str(task.root / "report.md"))
    assert res.returncode == 2 and "no such parent" in res.stderr
    res = _run("verify_handback.sh", task.repo, task.parent, str(task.root / "gone.md"))
    assert res.returncode == 2 and "not found" in res.stderr
    res = task.run("--no-tests", "--agent", "Q")
    assert res.returncode == 2 and "must be I or N" in res.stderr


def test_tooling_scripts_never_touch_the_real_clone_config(world: World) -> None:
    """The fixtures set icenine.agent in temporary repositories only."""
    _git(world.own, "config", "--local", "icenine.agent", "N")
    assert _git(world.own, "config", "--local", "--get", "icenine.agent") == "N"
    assert str(world.own).startswith(str(world.root))
    assert not str(world.own).startswith(str(REPO))
