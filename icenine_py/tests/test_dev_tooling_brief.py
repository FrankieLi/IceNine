"""Tests of scripts/dev/make_brief.sh on fixture plans and a fixture DECISIONS.md.

Every test runs in a temporary repository; the real clone's git config is never read or written
(the global and system configs are switched off for the subprocess, and the script runs with the
temporary repository as its working directory).
"""

import os
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "dev" / "make_brief.sh"
GIT_ENV = {
    "GIT_AUTHOR_NAME": "t",
    "GIT_AUTHOR_EMAIL": "t@example.com",
    "GIT_COMMITTER_NAME": "t",
    "GIT_COMMITTER_EMAIL": "t@example.com",
    "GIT_CONFIG_GLOBAL": os.devnull,
    "GIT_CONFIG_NOSYSTEM": "1",
}
TRAILER = "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
NOW = "20261018-1530"

PLAN = """\
---
title: "TODO: demo"
---

# Status

Approved.

# Steps

## 1. First step

Do the first thing.

## 2. Second step

Do the second thing.
"""

STOP_ONE = """\
**Stop rule:** if step 1 fails, stop and
report. The premise is then wrong.
"""

DECISIONS = """\
# Decision log

Header text that mentions todo_demo.md and approved, before the entries.

## Entries

### D-20261001-I1 · H1 · Something unrelated
- 2026-10-01, owner. Recorded by Ilium Agent. Supersedes: none.
- Nothing about plans here.

### D-20261002-I1 · H2 · Demo plan discussed
- 2026-10-02, owner. Recorded by Ilium Agent. Supersedes: none.
- The plan `docs/todo_demo.md` needs a rewrite first. Where: the plan.

### D-20261003-I1 · H2 · Other plan approved
- 2026-10-03, owner. Recorded by Ilium Agent. Supersedes: none.
- The plan `docs/todo_other.md` is approved. Where: the plan.

### D-20261004-I1 · H2 · Demo plan approved
- 2026-10-04, owner. Recorded by Ilium Agent. Supersedes: none.
- The plan `docs/todo_demo.md` is Approved with its grouping. Where: the plan.
"""


GIT_LOCATION_VARS = ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE", "GIT_COMMON_DIR")


def _clean_environ() -> Dict[str, str]:
    """os.environ without MAKE_BRIEF_TRAILER and without git variables that pick a repository."""
    skip = ("MAKE_BRIEF_TRAILER", *GIT_LOCATION_VARS)
    return {k: v for k, v in os.environ.items() if k not in skip}


def _git(repo: Path, *args: str) -> str:
    env = {**_clean_environ(), **GIT_ENV}
    out = subprocess.run(["git", *args], cwd=repo, env=env, capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    return out.stdout.strip()


def _write(repo: Path, rel: str, text: str) -> None:
    path = repo / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


@pytest.fixture()
def repo(tmp_path: Path) -> Path:
    path = tmp_path / "clone"
    path.mkdir()
    _git(path, "init", "-q", "-b", "feature/par")
    _write(path, "icenine_py/docs/todo_demo.md", PLAN)
    _write(path, "icenine_py/docs/DECISIONS.md", DECISIONS)
    _git(path, "add", "icenine_py/docs/todo_demo.md", "icenine_py/docs/DECISIONS.md")
    _git(path, "commit", "-q", "-m", "init")
    _git(path, "config", "--local", "icenine.maxHeavyProcs", "6")
    return path


def _brief(
    repo: Path,
    *args: str,
    env: Optional[Dict[str, str]] = None,
    trailer: Optional[str] = TRAILER,
) -> subprocess.CompletedProcess:
    full_env = _clean_environ()
    full_env.update(GIT_ENV)
    full_env["MAKE_BRIEF_NOW"] = NOW
    full_env.update(env or {})
    cmd: List[str] = [str(SCRIPT), *args]
    if trailer is not None:
        cmd += ["--trailer", trailer]
    return subprocess.run(cmd, cwd=repo, env=full_env, capture_output=True, text=True)


def _lines(*parts: str) -> str:
    return "\n".join(parts) + "\n"


NO_DEP = (
    "No new dependency (`pyproject.toml`, `uv.lock` unchanged); no default flip; no golden change."
)
COMMITS = "Commits: stage by explicit path; end each commit message with:"
REFUSALS = (
    "Refusals: if a tool or hook refuses an action (Edit or Write refused for a path, a hook "
    "denial, a declined permission), stop and report it; never redo the action another way, "
    "such as writing files through Bash."
)
INTERRUPTIONS = (
    "Interruptions: commit each finished piece as you go (sleep or an app restart can stop you "
    "at any point); when resumed, start from `git status` and `git log`."
)
REPORT = (
    "Report: write the full report to `.claude/reports/{name}.md`; "
    "the final message is at most about 250 words."
)
APPROVAL = "Approval: D-20261004-I1 · H2 · Demo plan approved"


def test_full_brief_with_stop_rule_default_scope(repo: Path) -> None:
    plan = PLAN + "\n" + STOP_ONE
    _write(repo, "icenine_py/docs/todo_demo.md", plan)
    out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "demo-a")
    assert out.returncode == 0, out.stderr
    assert out.stdout == _lines(
        "Follow study-conventions.",
        "Implement the plan icenine_py/docs/todo_demo.md.",
        APPROVAL,
        "Branch: in your own worktree run `git checkout -b feature/par-demo-a feature/par` and "
        "`git config branch.feature/par-demo-a.parent feature/par`.",
        "Scope: no change under `icenine_py/icenine/`.",
        NO_DEP,
        "Stop rule of the plan:",
        "> **Stop rule:** if step 1 fails, stop and",
        "> report. The premise is then wrong.",
        REFUSALS,
        INTERRUPTIONS,
        "Compute: at most 6 heavy processes at once on this machine.",
        COMMITS,
        TRAILER,
        REPORT.format(name=f"{NOW}-demo-a"),
    )


def test_section_and_default_stop_rule(repo: Path) -> None:
    out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "t", "--section", "2. Second")
    assert out.returncode == 0, out.stderr
    assert out.stdout == _lines(
        "Follow study-conventions.",
        'Implement section "2. Second" of icenine_py/docs/todo_demo.md.',
        APPROVAL,
        "Branch: in your own worktree run `git checkout -b feature/par-t feature/par` and "
        "`git config branch.feature/par-t.parent feature/par`.",
        "Scope: no change under `icenine_py/icenine/`.",
        NO_DEP,
        "Stop rule: stop and report if the plan cannot be followed.",
        REFUSALS,
        INTERRUPTIONS,
        "Compute: at most 6 heavy processes at once on this machine.",
        COMMITS,
        TRAILER,
        REPORT.format(name=f"{NOW}-t"),
    )


def test_unknown_section_refused(repo: Path) -> None:
    out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "t", "--section", "9. Missing")
    assert out.returncode == 1
    assert out.stdout == ""
    assert "no heading containing '9. Missing'" in out.stderr


def test_library_changes_exception_and_two_stop_rules(repo: Path) -> None:
    plan = (
        PLAN
        + "\n**Library changes:** `icenine/cost_functions.py` may change (flag-gated).\n"
        + "\n**Process limit exception:** 8 workers for the sweep in step 2.\n"
        + "\n"
        + STOP_ONE
        + "\nSome text between.\n\n"
        + "  **Stop rule:** indented, so not a rule.\n\n"
        + "**Stop rule:** a second rule\n  on two lines.\n\nAfter.\n"
    )
    _write(repo, "icenine_py/docs/todo_demo.md", plan)
    out = _brief(repo, "--plan", "icenine_py/docs/todo_demo.md", "--task", "t")
    assert out.returncode == 0, out.stderr
    assert out.stdout == _lines(
        "Follow study-conventions.",
        "Implement the plan icenine_py/docs/todo_demo.md.",
        APPROVAL,
        "Branch: in your own worktree run `git checkout -b feature/par-t feature/par` and "
        "`git config branch.feature/par-t.parent feature/par`.",
        "Scope: **Library changes:** `icenine/cost_functions.py` may change (flag-gated).",
        NO_DEP,
        "Stop rule of the plan:",
        "> **Stop rule:** if step 1 fails, stop and",
        "> report. The premise is then wrong.",
        "> **Stop rule:** a second rule",
        ">   on two lines.",
        REFUSALS,
        INTERRUPTIONS,
        "Compute (plan exception to the machine limit of 6 heavy processes): "
        "**Process limit exception:** 8 workers for the sweep in step 2.",
        COMMITS,
        TRAILER,
        REPORT.format(name=f"{NOW}-t"),
    )


def test_process_limit_comes_from_config(repo: Path) -> None:
    _git(repo, "config", "--local", "icenine.maxHeavyProcs", "3")
    out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "t")
    assert out.returncode == 0, out.stderr
    assert "Compute: at most 3 heavy processes at once on this machine.\n" in out.stdout


def test_process_limit_unset_or_invalid_refused(repo: Path) -> None:
    for value in ("", "0", "six", "-1", "06"):
        _git(repo, "config", "--local", "icenine.maxHeavyProcs", value)
        out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "t")
        assert out.returncode == 1, value
        assert out.stdout == ""
        assert "icenine.maxHeavyProcs" in out.stderr
    _git(repo, "config", "--local", "--unset", "icenine.maxHeavyProcs")
    out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "t")
    assert out.returncode == 1
    assert "icenine.maxHeavyProcs is not set" in out.stderr
    assert out.stdout == ""


def test_exception_does_not_waive_the_machine_limit_check(repo: Path) -> None:
    _write(
        repo,
        "icenine_py/docs/todo_demo.md",
        PLAN + "\n**Process limit exception:** 8 workers.\n",
    )
    _git(repo, "config", "--local", "--unset", "icenine.maxHeavyProcs")
    out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "t")
    assert out.returncode == 1
    assert "icenine.maxHeavyProcs" in out.stderr


def test_trailer_from_environment_and_missing_trailer(repo: Path) -> None:
    out = _brief(
        repo,
        "--plan",
        "docs/todo_demo.md",
        "--task",
        "t",
        trailer=None,
        env={"MAKE_BRIEF_TRAILER": "Co-Authored-By: Env <e@example.com>"},
    )
    assert out.returncode == 0, out.stderr
    assert "\nCo-Authored-By: Env <e@example.com>\nReport:" in out.stdout
    out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "t", trailer=None)
    assert out.returncode == 1
    assert out.stdout == ""
    assert "no commit trailer" in out.stderr


def test_approval_refusals(repo: Path) -> None:
    # an entry that names the plan without "approved", and an approval of another plan
    _write(repo, "icenine_py/docs/todo_new.md", PLAN)
    _write(
        repo,
        "icenine_py/docs/DECISIONS.md",
        DECISIONS.replace("is Approved with", "is accepted with").replace(
            "Demo plan approved", "Demo plan accepted"
        )
        + "\n### D-20261005-I1 · H2 · Notes\n- 2026-10-05, owner.\n- todo_new.md exists.\n",
    )
    out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "t")
    assert out.returncode == 1
    assert out.stdout == ""
    assert "no entry" in out.stderr and "todo_demo.md" in out.stderr
    out = _brief(repo, "--plan", "docs/todo_new.md", "--task", "t")
    assert out.returncode == 1
    assert "todo_new.md" in out.stderr


def test_approval_is_case_insensitive_and_takes_the_latest(repo: Path) -> None:
    text = DECISIONS.replace("is Approved with", "is APPROVED with") + (
        "\n### D-20261006-I2 · H2 · Demo plan approved again\n"
        "- 2026-10-06, owner. Recorded by Ilium Agent. Supersedes: none.\n"
        "- Plan `docs/todo_demo.md` approved for the second group.\n"
    )
    _write(repo, "icenine_py/docs/DECISIONS.md", text)
    out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "t")
    assert out.returncode == 0, out.stderr
    assert "\nApproval: D-20261006-I2 · H2 · Demo plan approved again\n" in out.stdout


def test_approval_from_decisions_option(repo: Path, tmp_path: Path) -> None:
    other = tmp_path / "other_decisions.md"
    other.write_text("### D-20261007-N1 · E · Demo\n- x.\n- todo_demo.md approved.\n")
    out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "t", "--decisions", str(other))
    assert out.returncode == 0, out.stderr
    assert "\nApproval: D-20261007-N1 · E · Demo\n" in out.stdout


def test_fixes_on_existing_task_branch(repo: Path) -> None:
    _git(repo, "checkout", "-q", "-b", "feature/par-demo-a")
    _git(repo, "config", "branch.feature/par-demo-a.parent", "feature/par")
    _write(repo, ".claude/reports/20261018-1600-review-demo-a.md", "# Review\n")
    _git(repo, "checkout", "-q", "feature/par")
    review = ".claude/reports/20261018-1600-review-demo-a.md"
    out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "demo-a", "--fixes", review)
    assert out.returncode == 0, out.stderr
    assert out.stdout == _lines(
        "Follow study-conventions.",
        "Plan: icenine_py/docs/todo_demo.md.",
        f"Apply the required fixes in {repo.resolve() / review}, follow study-conventions, "
        "re-run the affected tests, and add the new commit hashes to your report. Name the fix "
        "numbers in your commit subjects.",
        APPROVAL,
        "Branch: in your own worktree, run `git checkout feature/par-demo-a` (parent "
        "`feature/par`; the main session frees it first) and add new commits there.",
        "Scope: no change under `icenine_py/icenine/`.",
        NO_DEP,
        "Stop rule: stop and report if the plan cannot be followed.",
        REFUSALS,
        INTERRUPTIONS,
        "Compute: at most 6 heavy processes at once on this machine.",
        COMMITS,
        TRAILER,
        REPORT.format(name=f"{NOW}-demo-a-fixes"),
    )


def test_fixes_from_the_task_branch_itself(repo: Path) -> None:
    _git(repo, "checkout", "-q", "-b", "feature/par-demo-a")
    _git(repo, "config", "branch.feature/par-demo-a.parent", "feature/par")
    _write(repo, "review.md", "# Review\n")
    out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "demo-a", "--fixes", "review.md")
    assert out.returncode == 0, out.stderr
    assert (
        "Branch: in your own worktree, run `git checkout feature/par-demo-a` (parent "
        "`feature/par`; the main session frees it first)" in out.stdout
    )


def test_fixes_after_the_task_was_merged_use_a_fixes_branch(repo: Path) -> None:
    _write(repo, "review.md", "# Review\n")
    out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "demo-a", "--fixes", "review.md")
    assert out.returncode == 0, out.stderr
    assert (
        "Branch: the task is already merged, so in your own worktree run "
        "`git checkout -b feature/par-demo-a-fixes feature/par` and "
        "`git config branch.feature/par-demo-a-fixes.parent feature/par`.\n"
    ) in out.stdout
    assert f"{NOW}-demo-a-fixes.md" in out.stdout


def test_fixes_report_missing_refused(repo: Path) -> None:
    out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "t", "--fixes", "nope.md")
    assert out.returncode == 1
    assert "review report not found" in out.stderr


def test_other_refusals_and_usage(repo: Path) -> None:
    out = _brief(repo, "--plan", "docs/missing.md", "--task", "t")
    assert out.returncode == 1 and "plan not found" in out.stderr
    out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "bad name")
    assert out.returncode == 1 and "task name" in out.stderr
    out = _brief(repo, "--plan", "docs/todo_demo.md")
    assert out.returncode == 2
    out = _brief(repo, "--bogus")
    assert out.returncode == 2 and "unknown option" in out.stderr
    _git(repo, "checkout", "-q", "-b", "scratch")
    out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "t")
    assert out.returncode == 1 and "not a feature/* branch" in out.stderr
    assert out.stdout == ""


def test_plan_given_as_absolute_path_prints_a_repo_relative_path(repo: Path) -> None:
    out = _brief(repo, "--plan", str(repo / "icenine_py/docs/todo_demo.md"), "--task", "t")
    assert out.returncode == 0, out.stderr
    assert "\nImplement the plan icenine_py/docs/todo_demo.md.\n" in out.stdout
    assert str(repo) not in out.stdout
    assert str(Path.home()) not in out.stdout


def test_real_plan_and_decisions_of_this_repository(tmp_path: Path) -> None:
    """The plan of this very tooling is approved in the committed DECISIONS.md.

    The test pins no entry ID: a later entry may re-approve the plan; it fails only if the log has
    no approval of the plan or a later entry withdraws it.
    """
    env = _clean_environ()
    env.update(GIT_ENV)
    plan = REPO / "icenine_py" / "docs" / "todo_ilium_workflow_tooling.md"
    decisions = REPO / "icenine_py" / "docs" / "DECISIONS.md"
    if not plan.exists() or not decisions.exists():
        pytest.skip("plan or decision log not present")
    work = tmp_path / "w"
    work.mkdir()
    _git(work, "init", "-q", "-b", "feature/par")
    _git(work, "config", "--local", "icenine.maxHeavyProcs", "6")
    _write(work, "icenine_py/docs/" + plan.name, plan.read_text())
    _write(work, "icenine_py/docs/DECISIONS.md", decisions.read_text())
    out = subprocess.run(
        [str(SCRIPT), "--plan", plan.name, "--task", "tooling-c", "--trailer", TRAILER],
        cwd=work,
        env=env,
        capture_output=True,
        text=True,
    )
    assert out.returncode == 0, out.stderr
    assert "\nApproval: D-" in out.stdout
    assert "Implement the plan icenine_py/docs/todo_ilium_workflow_tooling.md." in out.stdout


def _entry(entry_id: str, title: str, text: str, supersedes: str = "none") -> str:
    return (
        f"\n### {entry_id} · H2 · {title}\n"
        f"- 2026-10-09, owner. Recorded by Ilium Agent. Supersedes: {supersedes}.\n"
        f"- {text}\n"
    )


def _without_approval() -> str:
    """The fixture log without its approval of todo_demo.md (the last entry)."""
    return DECISIONS[: DECISIONS.index("### D-20261004-I1")]


def test_a_refusal_is_not_an_approval(repo: Path) -> None:
    for text in (
        "The plan `docs/todo_demo.md` is not approved; rewrite first.",
        "The plan `docs/todo_demo.md` is unapproved.",
        "The plan `docs/todo_demo.md` has not yet been approved.",
        "The plan `docs/todo_demo.md` was never approved.",
        "The plan `docs/todo_demo.md` is disapproved.",
    ):
        log = _without_approval() + _entry("D-20261005-I1", "Demo plan", text)
        _write(repo, "icenine_py/docs/DECISIONS.md", log)
        out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "t")
        assert out.returncode == 1, text
        assert out.stdout == "" and "no entry" in out.stderr, text
    # a negation next to a real approval in the same entry does not hide the approval
    text = "The plan `docs/todo_demo.md` is not approved as drafted but approved with the edits."
    _write(
        repo,
        "icenine_py/docs/DECISIONS.md",
        _without_approval() + _entry("D-20261005-I1", "Demo plan", text),
    )
    out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "t")
    assert out.returncode == 0, out.stderr
    assert "\nApproval: D-20261005-I1 · H2 · Demo plan\n" in out.stdout


def test_a_withdrawn_or_superseded_approval_is_refused(repo: Path) -> None:
    withdrawn = _entry(
        "D-20261005-I1",
        "Demo plan withdrawn",
        "The plan `docs/todo_demo.md` is withdrawn.",
        supersedes="D-20261004-I1",
    )
    superseded_elsewhere = _entry(
        "D-20261005-I2", "Approval replaced", "Pause; no file named.", supersedes="D-20261004-I1"
    )
    for later, shown in ((withdrawn, "D-20261005-I1"), (superseded_elsewhere, "D-20261005-I2")):
        _write(repo, "icenine_py/docs/DECISIONS.md", DECISIONS + later)
        out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "t")
        assert out.returncode == 1, shown
        assert out.stdout == ""
        assert "D-20261004-I1" in out.stderr and shown in out.stderr, out.stderr
    # withdrawn without a Supersedes line
    log = DECISIONS + _entry("D-20261005-I3", "Demo plan", "`docs/todo_demo.md` needs a rewrite.")
    _write(repo, "icenine_py/docs/DECISIONS.md", log)
    out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "t")
    assert out.returncode == 1 and "D-20261005-I3" in out.stderr
    # a later entry that only mentions the plan does not cancel the approval
    log = DECISIONS + _entry("D-20261005-I4", "Cost", "Run `docs/todo_demo.md` on 6 processes.")
    _write(repo, "icenine_py/docs/DECISIONS.md", log)
    out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "t")
    assert out.returncode == 0, out.stderr
    assert APPROVAL in out.stdout
    # a new approval after a withdrawal counts again, with or without superseding the old one
    log = (
        DECISIONS
        + withdrawn
        + _entry(
            "D-20261006-I1",
            "Demo plan approved again",
            "The plan `docs/todo_demo.md` is approved.",
            supersedes="D-20261005-I1",
        )
    )
    _write(repo, "icenine_py/docs/DECISIONS.md", log)
    out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "t")
    assert out.returncode == 0, out.stderr
    assert "\nApproval: D-20261006-I1 · H2 · Demo plan approved again\n" in out.stdout


def test_plan_name_is_matched_on_a_path_boundary(repo: Path) -> None:
    for text in (
        "The plan `docs/old_todo_demo.md` is approved.",
        "The plan `docs/todo_demo.md2` is approved.",
        "The plan `docs/todo_demo.mdx` is approved.",
    ):
        log = _without_approval() + _entry("D-20261005-I1", "Other plan", text)
        _write(repo, "icenine_py/docs/DECISIONS.md", log)
        out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "t")
        assert out.returncode == 1, text
        assert "no entry" in out.stderr, text
    for text in (
        "todo_demo.md approved.",
        "Plan (todo_demo.md) approved.",
        'The plan "todo_demo.md" is approved.',
        "Approved: icenine_py/docs/todo_demo.md.",
        "The plan todo_demo.md. Approved.",
    ):
        log = _without_approval() + _entry("D-20261005-I1", "Demo", text)
        _write(repo, "icenine_py/docs/DECISIONS.md", log)
        out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "t")
        assert out.returncode == 0, (text, out.stderr)


def test_wrapped_library_changes_and_exception_are_quoted_whole(repo: Path) -> None:
    plan = (
        PLAN
        + "\n**Library changes:** `icenine/x.py` may change, but only behind a new flag that is\n"
        + "off by default; nothing else under `icenine/`.\n"
        + "\n**Process limit exception:** 8 workers for the sweep in step 2,\n"
        + "and 10 for the final check.\n\nAfter.\n"
    )
    _write(repo, "icenine_py/docs/todo_demo.md", plan)
    out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "t")
    assert out.returncode == 0, out.stderr
    assert (
        "Scope: **Library changes:** `icenine/x.py` may change, but only behind a new flag "
        "that is\noff by default; nothing else under `icenine/`.\n"
    ) in out.stdout
    assert (
        "Compute (plan exception to the machine limit of 6 heavy processes): "
        "**Process limit exception:** 8 workers for the sweep in step 2,\n"
        "and 10 for the final check.\n"
        "Commits:"
    ) in out.stdout
    assert "After." not in out.stdout


def test_markers_as_bullets_and_not_inside_code_fences(repo: Path) -> None:
    plan = (
        PLAN
        + "\n```\n**Stop rule:** inside a fence, an example only.\n**Library changes:** same.\n"
        + "```\n"
        + "\n- **Stop rule:** a bullet rule\n  that wraps.\n- **Stop rule:** a second bullet.\n"
        + "\n* **Library changes:** `icenine/y.py` only.\n\n"
        + "Text.\n```text\n**Process limit exception:** fenced, ignored.\n```\n"
        + "- **Process limit exception:** 9 for step 2.\n"
        + "# Next heading\n"
    )
    _write(repo, "icenine_py/docs/todo_demo.md", plan)
    out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "t")
    assert out.returncode == 0, out.stderr
    assert out.stdout == _lines(
        "Follow study-conventions.",
        "Implement the plan icenine_py/docs/todo_demo.md.",
        APPROVAL,
        "Branch: in your own worktree run `git checkout -b feature/par-t feature/par` and "
        "`git config branch.feature/par-t.parent feature/par`.",
        "Scope: **Library changes:** `icenine/y.py` only.",
        NO_DEP,
        "Stop rule of the plan:",
        "> **Stop rule:** a bullet rule",
        ">   that wraps.",
        "> **Stop rule:** a second bullet.",
        REFUSALS,
        INTERRUPTIONS,
        "Compute (plan exception to the machine limit of 6 heavy processes): "
        "**Process limit exception:** 9 for step 2.",
        COMMITS,
        TRAILER,
        REPORT.format(name=f"{NOW}-t"),
    )


def test_fixes_report_outside_the_repository_is_printed_as_an_absolute_path(
    repo: Path, tmp_path: Path
) -> None:
    other = tmp_path / "other worktree" / ".claude" / "reports"
    other.mkdir(parents=True)
    (other / "review.md").write_text("# Review\n")
    out = _brief(
        repo, "--plan", "docs/todo_demo.md", "--task", "t", "--fixes", str(other / "review.md")
    )
    assert out.returncode == 0, out.stderr
    expected = (other / "review.md").resolve()
    assert f"Apply the required fixes in {expected}, follow" in out.stdout
    # a relative path (from the current directory) is resolved the same way
    out = _brief(
        repo,
        "--plan",
        "docs/todo_demo.md",
        "--task",
        "t",
        "--fixes",
        os.path.relpath(other / "review.md", repo),
    )
    assert out.returncode == 0, out.stderr
    assert f"Apply the required fixes in {expected}, follow" in out.stdout


def test_section_must_match_exactly_one_heading_outside_code_fences(repo: Path) -> None:
    plan = (
        PLAN
        + "\n## 11. Eleven\n\nText.\n\n```bash\n# not a heading\n```\n\n"
        + "## 1. First step again\n\nText.\n"
    )
    _write(repo, "icenine_py/docs/todo_demo.md", plan)

    def run(section: str) -> subprocess.CompletedProcess:
        return _brief(repo, "--plan", "docs/todo_demo.md", "--task", "t", "--section", section)

    out = run("not a heading")
    assert out.returncode == 1 and out.stdout == ""
    assert "no heading containing 'not a heading'" in out.stderr
    for ambiguous, count in (("1", 3), ("First step", 2)):
        out = run(ambiguous)
        assert out.returncode == 1 and out.stdout == ""
        assert f"'{ambiguous}' matches {count} headings" in out.stderr, out.stderr
    # an exact heading wins over longer headings that contain it; one containing heading is enough
    for ok in ("1. First step", "Eleven", "Second step"):
        out = run(ok)
        assert out.returncode == 0, (ok, out.stderr)
        assert f'Implement section "{ok}" of icenine_py/docs/todo_demo.md.\n' in out.stdout


def test_bad_clock_override_refused(repo: Path) -> None:
    out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "t", env={"MAKE_BRIEF_NOW": "now"})
    assert out.returncode == 1 and out.stdout == ""
    assert "MAKE_BRIEF_NOW" in out.stderr
