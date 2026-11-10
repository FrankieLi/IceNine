"""Tests of scripts/dev/log_decision.sh and scripts/dev/roadmap_row.py.

Temporary repositories only: a bare repository stands in for `origin` (no network), a stub stands in
for shmem, and `icenine.agent` is set in the temporary repository, never in the real clone. The
user's and the system's git configuration are shut out (see `_env`). One test reads the real
ROADMAP.md, read-only (`--show`).
"""

import os
import stat
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional

import pytest

REPO = Path(__file__).resolve().parents[2]
DEV = REPO / "scripts" / "dev"
LOG_DECISION = DEV / "log_decision.sh"
ROADMAP_ROW = DEV / "roadmap_row.py"
DOC_REL = "icenine_py/docs/DECISIONS.md"
TRAILER = "Co-Authored-By: Test <test@example.com>"

GIT_ENV = {
    "GIT_AUTHOR_NAME": "t",
    "GIT_AUTHOR_EMAIL": "t@example.com",
    "GIT_COMMITTER_NAME": "t",
    "GIT_COMMITTER_EMAIL": "t@example.com",
}

HEADER = "# Decision log\n\n## Entries\n"

_HOME = tempfile.TemporaryDirectory(prefix="icenine-test-home-")  # an empty HOME for every command


def _env(extra: Optional[Dict[str, str]] = None) -> Dict[str, str]:
    """Environment without GIT_* and LOG_DECISION_* and without the user's or system git config."""
    env = {k: v for k, v in os.environ.items() if not k.startswith(("LOG_DECISION_", "GIT_"))}
    env.update(GIT_ENV)
    env.update(HOME=_HOME.name, GIT_CONFIG_GLOBAL=os.devnull, GIT_CONFIG_NOSYSTEM="1")
    env.update(extra or {})
    return env


def _git(repo: Path, *args: str) -> str:
    out = subprocess.run(["git", *args], cwd=repo, env=_env(), capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    return out.stdout.strip()


def _entry(id_: str, title: str = "t") -> str:
    return (
        f"\n### {id_} · row · {title}\n- 2026-01-01, owner. Recorded by X. Supersedes: none.\n"
        "- body. Where: here.\n"
    )


class Setup:
    """A temporary clone with `origin`, branch feature/z (has DECISIONS.md), task feature/z-job."""

    def __init__(self, tmp: Path, entries: str = "") -> None:
        assert tmp != REPO and REPO not in tmp.parents  # a temporary directory, never the clone
        self.tmp = tmp
        self.origin = tmp / "origin.git"
        subprocess.run(
            ["git", "init", "-q", "--bare", "-b", "develop", str(self.origin)],
            check=True, env=_env(),
        )  # fmt: skip
        self.repo = tmp / "r"
        self.repo.mkdir()
        _git(self.repo, "init", "-q", "-b", "develop")
        (self.repo / "README").write_text("x\n")
        _git(self.repo, "add", "README")
        _git(self.repo, "commit", "-q", "-m", "init")
        _git(self.repo, "remote", "add", "origin", str(self.origin))
        _git(self.repo, "checkout", "-q", "-b", "feature/z")
        (self.repo / DOC_REL).parent.mkdir(parents=True)
        (self.repo / DOC_REL).write_text(HEADER + entries)
        _git(self.repo, "add", DOC_REL)
        _git(self.repo, "commit", "-q", "-m", "decisions")
        _git(self.repo, "branch", "feature/z-job")
        _git(self.repo, "config", "branch.feature/z-job.parent", "feature/z")
        _git(self.repo, "push", "-q", "origin", "develop", "feature/z")
        self.shmem_log = tmp / "shmem.log"
        stub = tmp / "shmem_stub.sh"
        stub.write_text(f'#!/bin/bash\nprintf "%s\\n" "$*" >> "{self.shmem_log}"\n')
        stub.chmod(stub.stat().st_mode | stat.S_IXUSR)
        self.stub = stub

    def run(
        self,
        cwd: Path,
        *args: str,
        env: Optional[Dict[str, str]] = None,
        agent: Optional[str] = "I",
        base: bool = True,
    ) -> subprocess.CompletedProcess:
        full: List[str] = []
        if base:  # defaults first: a later option in `args` overrides them
            full += [
                "--row", "S9", "--title", "Tool title", "--decision", "The decision text",
                "--where", "ROADMAP S9",
            ]  # fmt: skip
        full += list(args)
        if agent is not None:
            full += ["--agent", agent]
        e = _env({"LOG_DECISION_SHMEM_CMD": str(self.stub), "LOG_DECISION_TRAILER": TRAILER})
        e.update(env or {})
        return subprocess.run(
            [str(LOG_DECISION), *full], cwd=cwd, env=e, capture_output=True, text=True
        )

    def doc(self, branch: str = "feature/z") -> str:
        return _git(self.repo, "show", f"{branch}:{DOC_REL}")

    def shmem_lines(self) -> List[str]:
        return self.shmem_log.read_text().splitlines() if self.shmem_log.exists() else []


@pytest.fixture()
def s(tmp_path: Path) -> Setup:
    return Setup(tmp_path)


def _task_worktree(s: Setup) -> Path:
    """Main checkout on feature/z (holds the target); the task branch lives in its own worktree."""
    wt = s.tmp / "wt"
    _git(s.repo, "worktree", "add", "-q", str(wt), "feature/z-job")
    return wt


def _heads(doc: str) -> List[str]:
    return [ln for ln in doc.splitlines() if ln.startswith("### ")]


# ---------------------------------------------------------------------------------------------
# log_decision.sh: IDs
# ---------------------------------------------------------------------------------------------


def test_id_on_a_new_date_is_n1(tmp_path: Path) -> None:
    s = Setup(tmp_path, _entry("D-20261014-I21"))
    res = s.run(s.repo, "--date", "2026-10-18")
    assert res.returncode == 0, res.stderr
    assert _heads(s.doc())[-1] == "### D-20261018-I1 · S9 · Tool title"


def test_id_counts_on_the_date_and_ignores_the_other_agent(tmp_path: Path) -> None:
    entries = (
        _entry("D-20261016-I22") + _entry("D-20261016-N7") + _entry("D-20261014-I21")
        + _entry("D-20261016-I2", "second")
    )  # fmt: skip
    s = Setup(tmp_path, entries)
    res = s.run(s.repo, "--date", "2026-10-16")
    assert res.returncode == 0, res.stderr
    assert _heads(s.doc())[-1].startswith("### D-20261016-I23 · ")
    res = s.run(s.repo, "--date", "2026-10-16", "--title", "Another", base=False,
                *["--row", "S9", "--decision", "d", "--where", "w"])  # fmt: skip
    assert res.returncode == 0, res.stderr
    assert _heads(s.doc())[-1].startswith("### D-20261016-I24 · ")


def test_agent_n_counts_its_own_entries_and_names_the_recorder(tmp_path: Path) -> None:
    s = Setup(tmp_path, _entry("D-20261016-I22") + _entry("D-20261016-N3"))
    res = s.run(s.repo, "--date", "2026-10-16", agent="N")
    assert res.returncode == 0, res.stderr
    doc = s.doc()
    assert _heads(doc)[-1] == "### D-20261016-N4 · S9 · Tool title"
    assert "Recorded by Newt Agent." in doc.splitlines()[-2]
    assert "Owner decision D-20261016-N4 (2026-10-16, IceNine S9)" in s.shmem_lines()[0]


def test_agent_i_recorder_and_recorder_override(s: Setup) -> None:
    assert s.run(s.repo, "--date", "2026-10-18").returncode == 0
    assert "Recorded by Ilium Agent." in s.doc()
    res = s.run(s.repo, "--date", "2026-10-19", "--recorder", "Someone Else")
    assert res.returncode == 0, res.stderr
    assert "Recorded by Someone Else." in s.doc()


def test_agent_from_git_config_and_refusal_when_unset(s: Setup) -> None:
    res = s.run(s.repo, "--date", "2026-10-18", agent=None)
    assert res.returncode == 2 and "agent unknown" in res.stderr
    assert s.doc().count("### ") == 0
    _git(s.repo, "config", "--local", "icenine.agent", "N")
    res = s.run(s.repo, "--date", "2026-10-18", agent=None)
    assert res.returncode == 0, res.stderr
    assert _heads(s.doc())[-1].startswith("### D-20261018-N1 · ")
    res = s.run(s.repo, "--date", "2026-10-18", "--title", "x2", agent="Q")
    assert res.returncode == 2  # an unknown letter is not accepted either


# ---------------------------------------------------------------------------------------------
# log_decision.sh: entry format
# ---------------------------------------------------------------------------------------------


def test_entry_format_and_wrapping(s: Setup) -> None:
    long = " ".join(f"word{i}" for i in range(60))
    res = s.run(s.repo, "--date", "2026-10-18", "--supersedes", "none", base=False,
                *["--row", "S4", "--title", "A title", "--decision", long + ".",
                  "--where", "`docs/x.md`"])  # fmt: skip
    assert res.returncode == 0, res.stderr
    doc = s.doc()
    tail = doc.split("### D-20261018-I1")[1].splitlines()
    assert tail[0] == " · S4 · A title"
    assert tail[1] == "- 2026-10-18, owner. Recorded by Ilium Agent. Supersedes: none."
    body = tail[2:]
    assert body[0].startswith("- word0 word1 ")
    assert len(body) > 2 and all(ln.startswith("  ") for ln in body[1:])
    assert all(len(ln) <= 98 for ln in body)
    assert body[-1].endswith("Where: `docs/x.md`.")
    assert " ".join(w for ln in body for w in ln.split())[2:].count("..") == 0
    assert (s.repo / DOC_REL).read_text().endswith("\n") and "\n\n### D-20261018-I1" in doc


def test_supersedes_and_commit_message(s: Setup) -> None:
    assert s.run(s.repo, "--date", "2026-10-18").returncode == 0
    res = s.run(s.repo, "--date", "2026-10-18", "--title", "Second", "--supersedes",
                "D-20261018-I1")  # fmt: skip
    assert res.returncode == 0, res.stderr
    assert "Supersedes: D-20261018-I1." in s.doc()
    msg = _git(s.repo, "log", "-1", "--format=%B", "feature/z")
    assert msg.startswith("DECISIONS D-20261018-I2: Second")
    assert msg.endswith(TRAILER)
    assert _git(s.repo, "show", "--stat", "--format=", "feature/z").count("|") == 1  # one file


# ---------------------------------------------------------------------------------------------
# log_decision.sh: where the commit lands
# ---------------------------------------------------------------------------------------------


def test_commit_in_the_worktree_that_holds_the_parent(s: Setup) -> None:
    wt = _task_worktree(s)
    before = _git(wt, "rev-parse", "HEAD")
    res = s.run(wt, "--date", "2026-10-18")
    assert res.returncode == 0, res.stderr
    assert _git(s.repo, "rev-parse", "--abbrev-ref", "HEAD") == "feature/z"
    assert "D-20261018-I1" in (s.repo / DOC_REL).read_text()  # the holder's working tree updated
    assert _git(s.repo, "status", "--porcelain") == ""
    assert _git(wt, "rev-parse", "HEAD") == before  # the task branch did not move
    assert _git(wt, "status", "--porcelain") == ""
    assert _git(s.repo, "log", "-1", "--format=%s").startswith("DECISIONS D-20261018-I1")


def test_commit_with_a_temporary_worktree_when_no_one_holds_the_parent(s: Setup) -> None:
    _git(s.repo, "checkout", "-q", "feature/z-job")
    z_before = _git(s.repo, "rev-parse", "feature/z")
    doc_before = (s.repo / DOC_REL).read_text()
    res = s.run(s.repo, "--date", "2026-10-18")
    assert res.returncode == 0, res.stderr
    assert _git(s.repo, "rev-parse", "--abbrev-ref", "HEAD") == "feature/z-job"
    assert _git(s.repo, "status", "--porcelain") == ""
    assert (s.repo / DOC_REL).read_text() == doc_before  # the caller's tree is unchanged
    assert _git(s.repo, "rev-parse", "feature/z") != z_before
    assert "D-20261018-I1" in s.doc("feature/z")
    assert "D-20261018-I1" not in s.doc("feature/z-job")
    assert len(_git(s.repo, "worktree", "list").splitlines()) == 1  # the temp worktree is gone
    assert _git(s.origin, "rev-parse", "feature/z") == _git(s.repo, "rev-parse", "feature/z")


def test_commit_from_the_target_branch_itself(s: Setup) -> None:
    res = s.run(s.repo, "--date", "2026-10-18")
    assert res.returncode == 0, res.stderr
    assert "D-20261018-I1" in (s.repo / DOC_REL).read_text()


def test_chain_of_parents_reaches_the_feature_branch(s: Setup) -> None:
    _git(s.repo, "branch", "feature/z-job-sub", "feature/z-job")
    _git(s.repo, "config", "branch.feature/z-job-sub.parent", "feature/z-job")
    _git(s.repo, "checkout", "-q", "feature/z-job-sub")
    res = s.run(s.repo, "--date", "2026-10-18")
    assert res.returncode == 0, res.stderr
    assert "D-20261018-I1" in s.doc("feature/z")


def test_dirty_decisions_in_the_holder_worktree_is_refused(s: Setup) -> None:
    wt = _task_worktree(s)
    (s.repo / DOC_REL).write_text((s.repo / DOC_REL).read_text() + "\nscribble\n")
    res = s.run(wt, "--date", "2026-10-18")
    assert res.returncode == 1 and "uncommitted changes" in res.stderr
    assert s.shmem_lines() == []
    assert not (s.repo.parent / "r" / ".git" / "log_decision.lock").exists()


# ---------------------------------------------------------------------------------------------
# log_decision.sh: refusals
# ---------------------------------------------------------------------------------------------


def test_refuses_develop_detached_and_branches_without_a_feature_parent(s: Setup) -> None:
    before = _git(s.repo, "rev-parse", "feature/z")
    _git(s.repo, "checkout", "-q", "develop")
    res = s.run(s.repo)
    assert res.returncode == 1 and "not allowed" in res.stderr
    _git(s.repo, "checkout", "-q", "-b", "task-without-parent")
    assert s.run(s.repo).returncode == 1
    _git(s.repo, "config", "branch.task-without-parent.parent", "develop")
    res = s.run(s.repo)
    assert res.returncode == 1  # the chain ends at develop
    _git(s.repo, "checkout", "-q", "--detach", "feature/z")
    res = s.run(s.repo)
    assert res.returncode == 1 and "detached" in res.stderr
    assert _git(s.repo, "rev-parse", "feature/z") == before
    assert s.shmem_lines() == []


def test_refuses_a_duplicate_title_on_the_same_date(s: Setup) -> None:
    assert s.run(s.repo, "--date", "2026-10-18").returncode == 0
    res = s.run(s.repo, "--date", "2026-10-18")
    assert res.returncode == 1 and "exists already" in res.stderr
    assert s.doc().count("### ") == 1 and len(s.shmem_lines()) == 1
    assert s.run(s.repo, "--date", "2026-10-19").returncode == 0  # another date is fine


def test_refuses_unknown_supersedes_and_newlines(s: Setup) -> None:
    res = s.run(s.repo, "--date", "2026-10-18", "--supersedes", "D-20990101-I1")
    assert res.returncode == 1 and "not an existing entry" in res.stderr
    res = s.run(s.repo, "--date", "2026-10-18", "--title", "two\nlines")
    assert res.returncode == 2 and "newline" in res.stderr
    res = s.run(s.repo, "--date", "18-10-2026")
    assert res.returncode == 2
    assert s.doc().count("### ") == 0


def test_supersedes_is_a_literal_existing_id_not_a_pattern(tmp_path: Path) -> None:
    s = Setup(tmp_path, _entry("D-20261014-I21"))
    before = _git(s.repo, "rev-parse", "feature/z")
    for pattern in ["D-.*", "D-20261014-I2.", "D-20261014-I2", "D-2026101[0-9]-I21"]:
        res = s.run(s.repo, "--date", "2026-10-18", "--supersedes", pattern)
        assert res.returncode in (1, 2), pattern  # 2: malformed, 1: well formed but unknown
    assert s.run(s.repo, "--date", "2026-10-18", "--supersedes", "D-.*").returncode == 2
    res = s.run(s.repo, "--date", "2026-10-18", "--supersedes", "D-20261014-I2")
    assert res.returncode == 1 and "not an existing entry" in res.stderr  # a prefix is no match
    assert _git(s.repo, "rev-parse", "feature/z") == before and s.shmem_lines() == []
    res = s.run(s.repo, "--date", "2026-10-18", "--supersedes", "D-20261014-I21")
    assert res.returncode == 0, res.stderr
    assert "Supersedes: D-20261014-I21." in s.doc()


def test_refuses_without_a_trailer_but_dry_run_does_not_need_one(s: Setup) -> None:
    res = s.run(s.repo, "--date", "2026-10-18", env={"LOG_DECISION_TRAILER": ""})
    assert res.returncode == 2 and "trailer" in res.stderr
    res = s.run(s.repo, "--date", "2026-10-18", "--trailer", "Co-Authored-By: Flag <f@x>",
                env={"LOG_DECISION_TRAILER": ""})  # fmt: skip
    assert res.returncode == 0, res.stderr
    assert _git(s.repo, "log", "-1", "--format=%B", "feature/z").endswith(
        "Co-Authored-By: Flag <f@x>"
    )
    res = s.run(s.repo, "--date", "2026-10-19", "--dry-run", env={"LOG_DECISION_TRAILER": ""})
    assert res.returncode == 0, res.stderr


def test_refuses_when_origin_is_ahead(s: Setup) -> None:
    other = s.tmp / "other"
    _git(s.tmp, "clone", "-q", str(s.origin), str(other))
    _git(other, "checkout", "-q", "feature/z")
    (other / "new.txt").write_text("n\n")
    _git(other, "add", "new.txt")
    _git(other, "commit", "-q", "-m", "ahead")
    _git(other, "push", "-q", "origin", "feature/z")
    before = _git(s.repo, "rev-parse", "feature/z")
    res = s.run(s.repo, "--date", "2026-10-18")
    assert res.returncode == 1 and "merge origin/feature/z" in res.stderr
    assert _git(s.repo, "rev-parse", "feature/z") == before
    assert s.shmem_lines() == []
    _git(s.repo, "merge", "-q", "origin/feature/z")
    assert s.run(s.repo, "--date", "2026-10-18").returncode == 0


def test_origin_ahead_without_a_holder_hints_fast_forward_not_a_merge_into_the_task(
    s: Setup,
) -> None:
    other = s.tmp / "other"
    _git(s.tmp, "clone", "-q", str(s.origin), str(other))
    _git(other, "checkout", "-q", "feature/z")
    (other / "new.txt").write_text("n\n")
    _git(other, "add", "new.txt")
    _git(other, "commit", "-q", "-m", "ahead")
    _git(other, "push", "-q", "origin", "feature/z")
    _git(s.repo, "checkout", "-q", "feature/z-job")  # the caller is on the task; no holder
    res = s.run(s.repo, "--date", "2026-10-18")
    assert res.returncode == 1
    assert "git fetch origin feature/z:feature/z" in res.stderr
    assert "merge origin/feature/z in a checkout of feature/z" in res.stderr  # only the fallback
    assert "-C" not in res.stderr
    _git(s.repo, "fetch", "-q", "origin", "feature/z:feature/z")  # the hinted command
    assert s.run(s.repo, "--date", "2026-10-18").returncode == 0


def test_lock_is_released_after_a_failure_and_a_held_lock_refuses(s: Setup) -> None:
    lock = s.repo / ".git" / "log_decision.lock"
    res = s.run(s.repo, "--date", "2026-10-18", "--supersedes", "D-20990101-I1")
    assert res.returncode == 1 and not lock.exists()
    lock.mkdir()
    res = s.run(s.repo, "--date", "2026-10-18")
    assert res.returncode == 1 and "lock" in res.stderr
    assert lock.exists()  # a refused run must not remove another run's lock
    assert "owner unknown" in res.stderr and "rmdir" in res.stderr
    (lock / "owner").write_text("pid 4242, started 2026-10-18 09:00:00\n")
    res = s.run(s.repo, "--date", "2026-10-18")
    assert res.returncode == 1 and "pid 4242, started 2026-10-18 09:00:00" in res.stderr
    assert (lock / "owner").exists()
    (lock / "owner").unlink()
    lock.rmdir()
    assert s.run(s.repo, "--date", "2026-10-18").returncode == 0
    assert not lock.exists()


# ---------------------------------------------------------------------------------------------
# log_decision.sh: shmem, push, dry run
# ---------------------------------------------------------------------------------------------


def test_shmem_stub_gets_the_exact_line_and_push_updates_origin(s: Setup) -> None:
    res = s.run(s.repo, "--date", "2026-10-18", "--row", "S4", "--title", "Park it.",
                "--decision", "Parked, see plan.", "--where", "ROADMAP S4",
                base=False)  # fmt: skip
    assert res.returncode == 0, res.stderr
    assert s.shmem_lines() == [
        "add --source user Owner decision D-20261018-I1 (2026-10-18, IceNine S4): Park it. "
        "Parked, see plan."
    ]
    assert _git(s.origin, "rev-parse", "feature/z") == _git(s.repo, "rev-parse", "feature/z")
    assert "D-20261018-I1" in res.stdout and "ok: pushed feature/z" in res.stdout


def test_no_push_leaves_origin_alone(s: Setup) -> None:
    before = _git(s.origin, "rev-parse", "feature/z")
    res = s.run(s.repo, "--date", "2026-10-18", "--no-push")
    assert res.returncode == 0, res.stderr
    assert _git(s.origin, "rev-parse", "feature/z") == before
    assert _git(s.repo, "rev-parse", "feature/z") != before
    assert len(s.shmem_lines()) == 1


def test_dry_run_prints_and_changes_nothing(s: Setup) -> None:
    before = _git(s.repo, "rev-parse", "feature/z")
    res = s.run(s.repo, "--date", "2026-10-18", "--dry-run")
    assert res.returncode == 0, res.stderr
    assert "target: feature/z" in res.stdout and "id: D-20261018-I1" in res.stdout
    assert (
        "### D-20261018-I1 · S9 · Tool title" in res.stdout
        and "shmem: Owner decision" in res.stdout
    )
    assert _git(s.repo, "rev-parse", "feature/z") == before
    assert s.shmem_lines() == [] and not (s.repo / ".git" / "log_decision.lock").exists()


def test_failed_shmem_keeps_the_commit_and_prints_the_rerun_command(s: Setup) -> None:
    res = s.run(s.repo, "--date", "2026-10-18", env={"LOG_DECISION_SHMEM_CMD": "false"})
    assert res.returncode == 3  # partial failure: committed, shmem failed
    assert "FAILED: false add --source user" in res.stdout
    assert "D-20261018-I1" in s.doc()
    assert "ok: pushed" in res.stdout  # the push still ran
    assert _git(s.origin, "rev-parse", "feature/z") == _git(s.repo, "rev-parse", "feature/z")


def test_failed_push_keeps_the_commit_and_exits_3(s: Setup) -> None:
    _git(s.repo, "remote", "set-url", "origin", str(s.tmp / "missing.git"))
    res = s.run(s.repo, "--date", "2026-10-18", "--no-fetch")
    assert res.returncode == 3  # partial failure: committed and in shmem, push failed
    assert "FAILED: git push origin feature/z" in res.stdout
    assert "D-20261018-I1" in s.doc()
    assert len(s.shmem_lines()) == 1


def test_unreachable_origin_refuses_and_no_fetch_skips_the_check_with_a_warning(s: Setup) -> None:
    _git(s.repo, "remote", "set-url", "origin", str(s.tmp / "missing.git"))
    before = _git(s.repo, "rev-parse", "feature/z")
    res = s.run(s.repo, "--date", "2026-10-18", "--no-push")
    assert res.returncode == 1 and "cannot reach origin" in res.stderr
    assert "--no-fetch" in res.stderr
    assert _git(s.repo, "rev-parse", "feature/z") == before and s.shmem_lines() == []
    res = s.run(s.repo, "--date", "2026-10-18", "--no-push", "--no-fetch")
    assert res.returncode == 0, res.stderr
    assert "warning: --no-fetch" in res.stderr
    assert "D-20261018-I1" in s.doc()


def test_origin_without_the_branch_is_skipped_with_a_note(s: Setup) -> None:
    _git(s.origin, "branch", "-D", "feature/z")
    res = s.run(s.repo, "--date", "2026-10-18", "--no-push")
    assert res.returncode == 0, res.stderr
    assert "origin has no feature/z" in res.stderr


def test_printed_shmem_rerun_command_is_safe_to_paste(s: Setup, tmp_path: Path) -> None:
    stub = s.tmp / "flaky_shmem.sh"  # fails while STUB_FAIL is set, else records its arguments
    stub.write_text(
        '#!/bin/bash\n[ -z "${STUB_FAIL:-}" ] || exit 1\n'
        f'printf "%s\\n" "$4" >> "{s.shmem_log}"\n'  # $1 add, $2 --source, $3 user, $4 the line
    )
    stub.chmod(stub.stat().st_mode | stat.S_IXUSR)
    text = "Use `touch PWNED` and $(touch PWNED2) here, \"quoted\" and 'single' ; done"
    res = s.run(s.repo, "--date", "2026-10-18", "--decision", text,
                env={"LOG_DECISION_SHMEM_CMD": str(stub), "STUB_FAIL": "1"})  # fmt: skip
    assert res.returncode == 3
    printed = [ln for ln in res.stdout.splitlines() if ln.startswith("FAILED: ")]
    assert len(printed) == 1
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    again = subprocess.run(
        ["bash", "-c", printed[0][len("FAILED: "):]], cwd=scratch, env=_env(),
        capture_output=True, text=True,
    )  # fmt: skip
    assert again.returncode == 0, again.stderr
    assert list(scratch.iterdir()) == []  # no backtick or $(...) span ran
    (line,) = s.shmem_lines()
    assert line.startswith("Owner decision D-20261018-I1 (2026-10-18, IceNine S9): Tool title. ")
    assert text in line


def test_unfinished_commit_removes_the_temporary_worktree_and_the_lock(s: Setup) -> None:
    _git(s.repo, "checkout", "-q", "feature/z-job")  # no worktree holds feature/z now
    hook = s.repo / ".git" / "hooks" / "pre-commit"
    hook.write_text("#!/bin/sh\nexit 1\n")
    hook.chmod(hook.stat().st_mode | stat.S_IXUSR)
    tmpdir = s.tmp / "tmpdir"
    tmpdir.mkdir()
    before = _git(s.repo, "rev-parse", "feature/z")
    res = s.run(s.repo, "--date", "2026-10-18", env={"TMPDIR": str(tmpdir)})
    assert res.returncode == 1 and "commit on feature/z failed" in res.stderr
    assert list(tmpdir.iterdir()) == []  # the temporary worktree is gone
    assert _git(s.repo, "worktree", "list", "--porcelain").count("worktree ") == 1
    assert not (s.repo / ".git" / "log_decision.lock").exists()
    assert _git(s.repo, "rev-parse", "feature/z") == before and s.shmem_lines() == []
    hook.unlink()
    assert s.run(s.repo, "--date", "2026-10-18", env={"TMPDIR": str(tmpdir)}).returncode == 0


# ---------------------------------------------------------------------------------------------
# roadmap_row.py
# ---------------------------------------------------------------------------------------------

S3_ROW = (
    "| S3 | batched quick MC: all candidates in lockstep (about 58% of evaluations) | per-call "
    "overhead shared | **done** (2026-10-14, merged as 31881fa; `docs/todo_batched_quick_mc.md`). "
    "Opt-in keys `BatchedQuickMC` and `SearchDevice` (`cpu`, `mps`); the default stays off "
    "until a later large-scale test (D-20261014-I2)"
)  # a copy of the real S3 row's shape: no closing pipe

TRAIL = " "  # a space after a closing pipe (kept as written, so it is spelled out here)
ROADMAP = f"""# Tracks

| Track | What |
|---|---|
| E | evaluation |

# Track E: evaluation reduction (Newt Agent)

| Task | Owner ask | Status |
|---|---|---|
| A | counts differ | **done** 2026-10-13: ledger |
| E | smooth the cost | **started** 2026-10-14 |

# Track S: per-evaluation cost (Ilium Agent)

Text with a | in prose is not a row.

| Item | What | Expected | Status |
|---|---|---|---|
| S1 | remove loops | 1.2x | **in progress** on the branch; detail kept |
| S2 | fused C | 15 us | not started; after S3 |
{S3_ROW}
| S4 | level-0 screen | at most 2x | **started** |{TRAIL}
| S5 | translation | small | idea |

# Track H: housekeeping

| Item | Where | Fix |
|---|---|---|
| H1 | `sample.py` | `.save` to `.write` |
| H12 | `scripts/dev/` | workflow tooling |
| A | dup of the Track E row A | x |

# Track Q: open questions

| Item | What | Status |
|---|---|---|
| Q1 | a question | open |
"""


@pytest.fixture()
def road(tmp_path: Path) -> Path:
    repo = tmp_path / "rr"
    repo.mkdir()
    _git(repo, "init", "-q", "-b", "develop")
    (repo / "ROADMAP.md").write_text(ROADMAP)
    return repo


def _row(repo: Path, *args: str, agent: Optional[str] = "I") -> subprocess.CompletedProcess:
    full = [sys.executable, str(ROADMAP_ROW), "--doc", str(repo / "ROADMAP.md"), *args]
    if agent is not None:
        full += ["--agent", agent]
    return subprocess.run(full, cwd=repo, env=_env(), capture_output=True, text=True)


def _lines(repo: Path) -> List[str]:
    return (repo / "ROADMAP.md").read_text().split("\n")


def test_row_status_replaces_the_last_cell(road: Path) -> None:
    res = _row(road, "--row", "S2", "--status", "**done** 2026-10-18")
    assert res.returncode == 0, res.stderr
    assert "- | S2 | fused C | 15 us | not started; after S3 |" in res.stdout
    assert "+ | S2 | fused C | 15 us | **done** 2026-10-18 |" in res.stdout
    new = _lines(road)
    old = ROADMAP.split("\n")
    assert [i for i, (a, b) in enumerate(zip(old, new)) if a != b] == [old.index(
        "| S2 | fused C | 15 us | not started; after S3 |")]  # fmt: skip


def test_row_state_keeps_the_detail(road: Path) -> None:
    res = _row(road, "--row", "S1", "--state", "parked")
    assert res.returncode == 0, res.stderr
    assert "| S1 | remove loops | 1.2x | **parked** on the branch; detail kept |" in _lines(road)
    res = _row(road, "--row", "S1", "--state", "**done**")  # asterisks are accepted
    assert "**done** on the branch; detail kept |" in "\n".join(_lines(road))


def test_row_state_without_a_bold_span_exits_2(road: Path) -> None:
    res = _row(road, "--row", "S5", "--state", "parked")
    assert res.returncode == 2 and "bold" in res.stderr
    assert (road / "ROADMAP.md").read_text() == ROADMAP


def test_row_without_a_closing_pipe_keeps_its_form(road: Path) -> None:
    res = _row(road, "--row", "S3", "--state", "parked")
    assert res.returncode == 0, res.stderr
    new = [ln for ln in _lines(road) if ln.startswith("| S3 |")]
    assert len(new) == 1 and not new[0].endswith("|")
    assert "| **parked** (2026-10-14, merged as 31881fa;" in new[0]
    assert new[0].endswith("until a later large-scale test (D-20261014-I2)")
    res = _row(road, "--row", "S3", "--status", "**done** again")
    assert [ln for ln in _lines(road) if ln.startswith("| S3 |")][0].endswith("| **done** again")
    assert len([ln for ln in _lines(road) if ln.startswith("| S3 |")][0].split(" | ")) == 4


def test_row_trailing_space_after_the_pipe_is_kept(road: Path) -> None:
    assert _row(road, "--row", "S4", "--status", "**done**").returncode == 0
    assert "| S4 | level-0 screen | at most 2x | **done** | " in _lines(road)


def test_row_missing_and_duplicate(road: Path) -> None:
    res = _row(road, "--row", "S99", "--status", "x")
    assert res.returncode == 2 and "S99" in res.stderr and "Track S: S1, S2, S3" in res.stderr
    res = _row(road, "--row", "A", "--status", "x")  # in Track E and Track H
    assert res.returncode == 2 and "ambiguous" in res.stderr
    assert (road / "ROADMAP.md").read_text() == ROADMAP


def test_row_header_and_prose_pipes_are_not_rows(road: Path) -> None:
    res = _row(road, "--row", "Item", "--show")
    assert res.returncode == 2
    res = _row(road, "--row", "Text", "--show")
    assert res.returncode == 2


@pytest.mark.parametrize(
    "agent,row,ok",
    [
        ("I", "S1", True),
        ("I", "H12", True),
        ("I", "H1", False),
        ("I", "E", False),
        ("I", "A", False),
        ("N", "E", True),
        ("N", "S1", False),
        ("N", "H12", False),
    ],
)
def test_row_ownership_both_agents(road: Path, agent: str, row: str, ok: bool) -> None:
    if row == "A":  # make the Track H copy of A unique, so ownership (not ambiguity) decides
        (road / "ROADMAP.md").write_text(
            ROADMAP.replace("| A | dup of the Track E row A | x |\n", "")
        )
    res = _row(road, "--row", row, "--status", "**x** new", agent=agent)
    assert res.returncode == (0 if ok else 2), res.stderr
    if not ok:
        assert "not an agent-" in res.stderr
        assert (road / "ROADMAP.md").read_text() in (ROADMAP, ROADMAP.replace(
            "| A | dup of the Track E row A | x |\n", ""))  # fmt: skip


def test_row_agent_from_git_config_or_refusal(road: Path) -> None:
    res = _row(road, "--row", "S1", "--status", "x", agent=None)
    assert res.returncode == 2 and "agent unknown" in res.stderr
    _git(road, "config", "--local", "icenine.agent", "N")
    assert _row(road, "--row", "S1", "--status", "x", agent=None).returncode == 2
    assert _row(road, "--row", "E", "--status", "**x**", agent=None).returncode == 0
    _git(road, "config", "--local", "icenine.agent", "I")
    assert _row(road, "--row", "S1", "--status", "**x**", agent=None).returncode == 0


def test_row_show_needs_no_agent_and_changes_nothing(road: Path) -> None:
    res = _row(road, "--row", "S2", "--show", agent=None)
    assert (
        res.returncode == 0 and res.stdout == "| S2 | fused C | 15 us | not started; after S3 |\n"
    )
    res = _row(road, "--row", "H1", "--show", agent=None)  # read-only: any row may be shown
    assert res.returncode == 0 and "`.save` to `.write`" in res.stdout
    assert (road / "ROADMAP.md").read_text() == ROADMAP


def test_row_pipe_and_newline_in_the_text_are_refused(road: Path) -> None:
    for text in ["a | b", "line\nbreak", "  "]:
        res = _row(road, "--row", "S2", "--status", text)
        assert res.returncode == 2, text
    assert (road / "ROADMAP.md").read_text() == ROADMAP
    res = _row(road, "--row", "S2", "--status", "a \\| b")  # an escaped bar is allowed
    assert res.returncode == 0 and "| a \\| b |" in "\n".join(_lines(road))
    res = _row(road, "--row", "S1", "--state", "a|b")
    assert res.returncode == 2


def test_row_second_run_is_idempotent_and_check_reports_a_change(road: Path) -> None:
    res = _row(road, "--row", "S2", "--status", "**done**", "--check")
    assert res.returncode == 1 and "would change" in res.stdout
    assert (road / "ROADMAP.md").read_text() == ROADMAP  # --check writes nothing
    assert _row(road, "--row", "S2", "--status", "**done**").returncode == 0
    after = (road / "ROADMAP.md").read_text()
    res = _row(road, "--row", "S2", "--status", "**done**")
    assert res.returncode == 0 and "unchanged" in res.stdout
    assert (road / "ROADMAP.md").read_text() == after
    assert _row(road, "--row", "S2", "--status", "**done**", "--check").returncode == 0


def test_row_cell_count_mismatch_is_refused(road: Path) -> None:
    text = ROADMAP.replace("| S5 | translation | small | idea |", "| S5 | a | b | c | d | idea |")
    (road / "ROADMAP.md").write_text(text)
    res = _row(road, "--row", "S5", "--status", "x")
    assert res.returncode == 2 and "cells" in res.stderr
    assert (road / "ROADMAP.md").read_text() == text


def test_row_crlf_line_endings_are_kept(road: Path) -> None:
    (road / "ROADMAP.md").write_bytes(ROADMAP.replace("\n", "\r\n").encode())
    assert _row(road, "--row", "S2", "--status", "**done**").returncode == 0
    data = (road / "ROADMAP.md").read_bytes()
    assert b"| S2 | fused C | 15 us | **done** |\r\n" in data
    assert data.count(b"\r\n") == data.count(b"\n")


def test_row_real_roadmap_s3_row_has_no_closing_pipe_and_parses() -> None:
    real = REPO / "icenine_py" / "docs" / "ROADMAP.md"
    rows = [ln for ln in real.read_text().splitlines() if ln.startswith("| S3 |")]
    res = subprocess.run(
        [sys.executable, str(ROADMAP_ROW), "--doc", str(real), "--row", "S3", "--show"],
        env=_env(), capture_output=True, text=True,
    )  # fmt: skip
    assert res.returncode == 0, res.stderr
    assert res.stdout.rstrip("\n") == rows[0]
