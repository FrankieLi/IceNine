"""Tests of the S2 follow-ups to make_brief.sh and verify_handback.sh (tooling task (d)).

- `make_brief.sh --section` reads the "**Library changes:**" and "**Process limit exception:**"
  markers only from the section or from the plan-wide text (outside every "## " section under
  "# Steps"), never from another task's section; with neither, "none stated in this section".
  In a plan with no "# Steps" heading every "## " section is its own task's, and a stderr note
  names each marker paragraph left out (review fix 1, D-20261028-I5).
- `verify_handback.sh` does not flag a named hash the branch contains (its base or parent tip).
- `verify_handback.sh` reports pytest's real counts and exit code, never a warning line; exit 0
  with no counts is a FLAG (review fix 2).
- `finish_task.sh` runs the suite without an extra `-q` (D-20261028-I5).

Every test works in temporary repositories; pytest is a stub or this environment's own pytest.
"""

import os
import subprocess
import sys
from pathlib import Path
from typing import Dict

import pytest
from test_dev_tooling_brief import GIT_ENV, PLAN, REPO, _brief, _clean_environ
from test_dev_tooling_brief import _git as _bgit
from test_dev_tooling_brief import _write
from test_dev_tooling_brief import repo  # noqa: F401  (fixture)
from test_dev_tooling_handback import Task, _executable, _git, _report
from test_dev_tooling_handback import task  # noqa: F401  (fixture)

# ---------------------------------------------------------------------------------------------
# make_brief.sh --section
# ---------------------------------------------------------------------------------------------

S2_LIKE = """\
---
title: "TODO: demo"
---

# Status

Approved.

# Steps

Text before the first step.

## Task A check (scripts only; in parallel with Steps 1-2)

**Library changes:** none. Scripts only.

**Process limit exception:** at most 3 heavy processes at once (3 each).

## Steps 1-2: the kernel

**Library changes:** new kernel file.

**Process limit exception:** at most 3 heavy processes at once (3 each).

## Step 3: the batched entry

Do step 3.

**Stop rule:** stop if step 3 is slower.

## Follow-ups

Intro of the follow-ups.

### Task (d): briefs

**Library changes:** none (task d).

**Process limit exception:** 3 for task d.

### Task (e): tables

**Process limit exception:** 3 for task e.

### Task (f): nothing stated

Do f.

# Rules

Plain rules.
"""


def _plan(repo: Path, text: str) -> None:
    _write(repo, "icenine_py/docs/todo_demo.md", text)


def test_section_without_markers_does_not_take_another_sections(repo: Path) -> None:
    """The S2 case: Step 3 got the parallel tasks' "3 each" from the whole plan."""
    _plan(repo, S2_LIKE)
    out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "t", "--section", "Step 3")
    assert out.returncode == 0, out.stderr
    assert "3 each" not in out.stdout and "Scripts only" not in out.stdout
    assert "new kernel file" not in out.stdout
    assert (
        "Scope: library changes none stated in this section, so no change under "
        "`icenine_py/icenine/`.\n"
    ) in out.stdout
    assert (
        "Compute: process limit exception none stated in this section; at most 6 heavy "
        "processes at once on this machine.\n"
    ) in out.stdout
    # one stderr note per missing marker
    assert out.stderr.count('note: section "Step 3" has no') == 2
    assert "**Library changes:**" in out.stderr and "**Process limit exception:**" in out.stderr
    # the section's own stop rule is quoted
    assert "> **Stop rule:** stop if step 3 is slower." in out.stdout


def test_section_markers_are_used_and_no_note(repo: Path) -> None:
    _plan(repo, S2_LIKE)
    out = _brief(
        repo, "--plan", "docs/todo_demo.md", "--task", "t", "--section", "Steps 1-2: the kernel"
    )
    assert out.returncode == 0, out.stderr
    assert "Scope: **Library changes:** new kernel file.\n" in out.stdout
    assert (
        "Compute (plan exception to the machine limit of 6 heavy processes): "
        "**Process limit exception:** at most 3 heavy processes at once (3 each).\n"
    ) in out.stdout
    assert "Scripts only" not in out.stdout
    assert out.stderr == ""


def test_nested_task_sections_keep_to_their_own_markers(repo: Path) -> None:
    _plan(repo, S2_LIKE)
    d = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "d", "--section", "Task (d)")
    assert d.returncode == 0, d.stderr
    assert "Scope: **Library changes:** none (task d).\n" in d.stdout
    assert "**Process limit exception:** 3 for task d.\n" in d.stdout
    assert "task e" not in d.stdout
    e = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "e", "--section", "Task (e)")
    assert e.returncode == 0, e.stderr
    assert "**Process limit exception:** 3 for task e.\n" in e.stdout
    # task (e) states no library changes: not task (d)'s, not Task A's
    assert "library changes none stated in this section" in e.stdout
    assert "task d" not in e.stdout and "Scripts only" not in e.stdout
    assert e.stderr.count("note:") == 1
    f = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "f", "--section", "Task (f)")
    assert f.returncode == 0, f.stderr
    assert "task d" not in f.stdout and "task e" not in f.stdout and "3 each" not in f.stdout
    assert f.stdout.count("none stated in this section") == 2
    # the stop rule still falls back to the whole plan (a rule written elsewhere still applies)
    assert "> **Stop rule:** stop if step 3 is slower." in f.stdout


@pytest.mark.parametrize(
    "where",
    ["intro", "rules"],
)
def test_plan_wide_markers_outside_the_step_sections_apply(repo: Path, where: str) -> None:
    lib = "**Library changes:** plan-wide: none anywhere.\n"
    exc = "**Process limit exception:** plan-wide: 2 at once.\n"
    if where == "intro":
        text = S2_LIKE.replace("Text before the first step.\n", "Text.\n\n" + lib + "\n" + exc)
    else:
        text = S2_LIKE.replace("Plain rules.\n", "Plain rules.\n\n" + lib + "\n" + exc)
    _plan(repo, text)
    out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "t", "--section", "Step 3")
    assert out.returncode == 0, out.stderr
    assert "Scope: **Library changes:** plan-wide: none anywhere.\n" in out.stdout
    assert "**Process limit exception:** plan-wide: 2 at once.\n" in out.stdout
    assert "3 each" not in out.stdout and "none stated" not in out.stdout
    assert out.stderr == ""
    # a section's own marker still wins over the plan-wide one
    own = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "t", "--section", "Task (d)")
    assert "Scope: **Library changes:** none (task d).\n" in own.stdout
    assert "plan-wide: none anywhere" not in own.stdout


def test_follow_ups_intro_is_not_plan_wide(repo: Path) -> None:
    """Text of a "## " section under "# Steps" belongs to that section, not to the plan."""
    text = S2_LIKE.replace(
        "Intro of the follow-ups.\n",
        "Intro of the follow-ups.\n\n**Library changes:** follow-ups intro only.\n",
    )
    _plan(repo, text)
    out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "t", "--section", "Step 3")
    assert "follow-ups intro only" not in out.stdout
    assert "library changes none stated in this section" in out.stdout
    # inside the follow-ups, task (f) is in that section's range, so the intro applies
    f = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "f", "--section", "Follow-ups")
    assert "Scope: **Library changes:** follow-ups intro only.\n" in f.stdout


def test_fenced_headings_do_not_open_step_sections(repo: Path) -> None:
    text = S2_LIKE.replace(
        "Plain rules.\n",
        "Plain rules.\n\n```\n# Steps\n## fake\n```\n\n**Library changes:** rules-wide.\n",
    )
    _plan(repo, text)
    out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "t", "--section", "Step 3")
    assert "Scope: **Library changes:** rules-wide.\n" in out.stdout


def test_without_section_behaviour_is_unchanged(repo: Path) -> None:
    _plan(repo, S2_LIKE)
    out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "t")
    assert out.returncode == 0, out.stderr
    # every marker paragraph of the whole plan, as before
    assert "Scope: **Library changes:** none. Scripts only." in out.stdout
    assert "new kernel file" in out.stdout and "3 for task e" in out.stdout
    assert "none stated" not in out.stdout and out.stderr == ""
    # a plan with no markers at all keeps the old default lines
    _plan(repo, PLAN)
    plain = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "t")
    assert "Scope: no change under `icenine_py/icenine/`.\n" in plain.stdout
    assert "Compute: at most 6 heavy processes at once on this machine.\n" in plain.stdout


# A plan from the new_todo.sh template: no "# Steps" heading (D-20261028-I5, review fix 1).
NO_STEPS = """\
---
title: "TODO: demo"
---

# Status

Approved.

# Experiments (sketch)

Plan-wide text.

## Task A

**Library changes:** A changes cost_functions.py.

**Process limit exception:** 3 each for A.

**Stop rule:** stop if A is slower.

## Task B

Do B.

# Risks

## Compute

**Process limit exception:** risks-section exception.
"""


def _line(text: str, needle: str) -> int:
    """The 1-based line number of the first line of text that contains needle."""
    return next(i for i, x in enumerate(text.splitlines(), 1) if needle in x)


def test_no_steps_heading_every_section_is_its_own_tasks(repo: Path) -> None:
    """The review's leak: Task B got Task A's permissions, with nothing on stderr."""
    _plan(repo, NO_STEPS)
    out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "b", "--section", "Task B")
    assert out.returncode == 0, out.stderr
    assert "cost_functions" not in out.stdout and "3 each" not in out.stdout
    assert "risks-section" not in out.stdout  # a "## " section outside Task B, not plan-wide
    assert out.stdout.count("none stated in this section") == 2
    # each skipped marker paragraph is named with its line and its section
    plan = "icenine_py/docs/todo_demo.md"
    for needle, marker, head in [
        ("A changes", "**Library changes:**", "## Task A"),
        ("3 each for A", "**Process limit exception:**", "## Task A"),
        ("risks-section", "**Process limit exception:**", "## Compute"),
    ]:
        note = (
            f"note: skipped the {marker} paragraph at line {_line(NO_STEPS, needle)} of {plan}: "
            f'it is in "{head}" and the plan has no "# Steps" heading, so every "## " section '
            "is its own task's\n"
        )
        assert note in out.stderr, out.stderr
    assert out.stderr.count("note: skipped") == 3
    assert out.stderr.count('(the plan has no "# Steps" heading)') == 2
    # the stop rule keeps its whole-plan fallback
    assert "> **Stop rule:** stop if A is slower." in out.stdout
    # the named section's own markers: used, no note
    a = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "a", "--section", "Task A")
    assert a.returncode == 0, a.stderr
    assert "Scope: **Library changes:** A changes cost_functions.py.\n" in a.stdout
    assert "**Process limit exception:** 3 each for A.\n" in a.stdout
    assert "risks-section" not in a.stdout and a.stderr == ""


def test_no_steps_heading_text_outside_every_section_is_plan_wide(repo: Path) -> None:
    text = NO_STEPS.replace(
        "Plan-wide text.\n", "Plan-wide text.\n\n**Library changes:** plan-wide: none.\n"
    )
    _plan(repo, text)
    out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "b", "--section", "Task B")
    assert out.returncode == 0, out.stderr
    assert "Scope: **Library changes:** plan-wide: none.\n" in out.stdout
    assert "cost_functions" not in out.stdout
    # Task A's library paragraph was still left out, and the note says so
    assert f"**Library changes:** paragraph at line {_line(text, 'A changes')} " in out.stderr
    assert "has no **Library changes:** paragraph" not in out.stderr


def test_with_steps_heading_other_hash_hash_sections_are_plan_wide(repo: Path) -> None:
    """With "# Steps", only the "## " sections under it are task-scoped (unchanged)."""
    text = NO_STEPS.replace("# Experiments (sketch)\n", "# Steps\n")
    _plan(repo, text)
    out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "b", "--section", "Task B")
    assert out.returncode == 0, out.stderr
    assert "**Process limit exception:** risks-section exception.\n" in out.stdout
    assert "3 each" not in out.stdout and "cost_functions" not in out.stdout
    assert "library changes none stated in this section" in out.stdout
    assert out.stderr == (
        'note: section "Task B" has no **Library changes:** paragraph and the plan has none '
        'outside the "## " sections under "# Steps"; the brief says "none stated in this '
        'section"\n'
    )


def test_no_steps_heading_without_section_is_unchanged(repo: Path) -> None:
    _plan(repo, NO_STEPS)
    out = _brief(repo, "--plan", "docs/todo_demo.md", "--task", "t")
    assert out.returncode == 0, out.stderr
    assert "Scope: **Library changes:** A changes cost_functions.py.\n" in out.stdout
    assert "3 each for A." in out.stdout and "risks-section exception." in out.stdout
    assert "none stated" not in out.stdout and out.stderr == ""


# ---------------------------------------------------------------------------------------------
# finish_task.sh: default test command
# ---------------------------------------------------------------------------------------------


def test_finish_task_default_pytest_has_no_extra_q(repo: Path, tmp_path: Path) -> None:
    """addopts has -q; a second -q (verbosity -2) drops pytest's count line from the output."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    log = tmp_path / "uv_args.txt"
    _executable(bin_dir / "uv", f'#!/bin/bash\necho "$(basename "$PWD") $*" > "{log}"\n')
    _bgit(repo, "checkout", "-q", "-b", "feature/par-job")
    _bgit(repo, "config", "branch.feature/par-job.parent", "feature/par")
    _write(repo, "icenine_py/x.py", "x = 1\n")
    _bgit(repo, "add", "icenine_py/x.py")
    _bgit(repo, "commit", "-q", "-m", "add x")
    env = {**_clean_environ(), **GIT_ENV, "PATH": f"{bin_dir}:{os.environ['PATH']}"}
    env.pop("FINISH_TASK_TEST_CMD", None)
    res = subprocess.run(
        [str(REPO / "scripts/dev/finish_task.sh")],
        cwd=repo,
        env=env,
        capture_output=True,
        text=True,
    )
    assert res.returncode == 0, res.stderr
    assert log.read_text() == "icenine_py run pytest tests/\n"
    assert _bgit(repo, "rev-parse", "--abbrev-ref", "HEAD") == "feature/par"


# ---------------------------------------------------------------------------------------------
# verify_handback.sh: hashes
# ---------------------------------------------------------------------------------------------


def _named(task: Task, *extra: str) -> None:
    """A report naming the branch's commits and the extra hashes."""
    _report(task.root / "report.md", task.hashes + list(extra))


def test_base_and_parent_tip_are_not_flagged(task: Task) -> None:
    base = _git(task.repo, "rev-parse", task.parent)
    task.commit({"a.txt": "1\n"})
    _named(task, base)
    res = task.run("--no-tests")
    assert res.returncode == 0, res.stdout
    assert "hashes not on the branch" not in res.stdout
    assert "ok    commit hashes agree between report and branch" in res.stdout
    assert any(x.startswith("note") and base[:9] in x for x in res.stdout.splitlines())
    assert "summary: 0 FAIL, 0 FLAG" in res.stdout


def test_hashes_the_branch_does_not_contain_are_flagged(task: Task) -> None:
    task.commit({"a.txt": "1\n"})
    # the parent moves on after the branch started: its new tip is not on the branch
    _git(task.repo, "checkout", "-q", task.parent)
    later = _commit_here(task, {"p.txt": "p\n"})
    _git(task.repo, "checkout", "-q", "-b", "other", task.parent + "~1")
    side = _commit_here(task, {"o.txt": "o\n"})
    _git(task.repo, "checkout", "-q", task.branch)
    _named(task, later, side)
    res = task.run("--no-tests")
    line = next(x for x in res.stdout.splitlines() if "hashes not on the branch" in x)
    assert line.startswith("FLAG") and later[:9] in line and side[:9] in line
    assert res.returncode == 0  # FLAG only


def _commit_here(task: Task, files: Dict[str, str]) -> str:
    for rel, text in files.items():
        (task.repo / rel).write_text(text)
        _git(task.repo, "add", rel)
    _git(task.repo, "commit", "-q", "-m", "elsewhere")
    return _git(task.repo, "rev-parse", "HEAD")


# ---------------------------------------------------------------------------------------------
# verify_handback.sh: test summary
# ---------------------------------------------------------------------------------------------

DOCS = "-- Docs: https://docs.pytest.org/en/stable/how-to/capture-warnings.html"


def _stub(task: Task, body: str) -> Dict[str, str]:
    stub = _executable(task.root / "pytest_stub2.sh", "#!/bin/bash\n" + body)
    return {"VERIFY_HANDBACK_PYTEST_CMD": str(stub)}


def _junit_writer(tests: int, failures: int, errors: int, skipped: int, rc: int) -> str:
    xml = (
        '<?xml version="1.0" encoding="utf-8"?><testsuites name="pytest tests">'
        f'<testsuite name="pytest" errors="{errors}" failures="{failures}" '
        f'skipped="{skipped}" tests="{tests}" time="0.1">'
        "</testsuite></testsuites>"
    )
    return (
        "j=$(printf '%s' \"$PYTEST_ADDOPTS\" | sed -E \"s/.*--junitxml='([^']*)'.*/\\1/\")\n"
        f"printf '%s\\n' '{xml}' > \"$j\"\n"
        f"echo '{DOCS}'\nexit {rc}\n"
    )


def _test_line(stdout: str) -> str:
    return next(x for x in stdout.splitlines() if "named tests" in x)


def test_counts_from_junit_never_the_warning_line(task: Task) -> None:
    task.commit({"a.txt": "1\n"})
    task.report()
    res = task.run(env=_stub(task, _junit_writer(7, 0, 0, 2, 0)))
    assert res.returncode == 0, res.stdout
    line = _test_line(res.stdout)
    assert (
        line == "ok    named tests passed (2 path(s)): 5 passed, 0 failed, 2 skipped; pytest exit 0"
    )
    assert "Docs" not in res.stdout


def test_junit_failures_and_errors_are_fail(task: Task) -> None:
    task.commit({"a.txt": "1\n"})
    task.report()
    res = task.run(env=_stub(task, _junit_writer(6, 1, 1, 0, 1)))
    assert res.returncode == 1
    assert "FAIL  named tests failed (2 path(s)): 4 passed, 1 failed, 0 skipped, 1 errors; " in (
        res.stdout
    )
    assert "pytest exit 1" in res.stdout
    # failures in the counts stay FAIL even when the command exits 0
    res = task.run(env=_stub(task, _junit_writer(3, 1, 0, 0, 0)))
    assert res.returncode == 1 and "FAIL  named tests failed" in res.stdout
    assert "pytest exit 0" in res.stdout


def test_count_line_fallback_and_missing_count_line(task: Task) -> None:
    task.commit({"a.txt": "1\n"})
    task.report()
    out = "==== 4 passed, 1 skipped, 3 warnings in 0.52s (0:00:00) ===="
    res = task.run(env=_stub(task, f"echo '{out}'\necho\nexit 0\n"))
    assert res.returncode == 0, res.stdout
    assert _test_line(res.stdout).endswith(": 4 passed, 0 failed, 1 skipped; pytest exit 0")
    # verbosity -2: no count line, the last line is the warnings' Docs link; with exit 0 the
    # gate cannot show any count, so it is a FLAG (not a FAIL: exit code 0)
    res = task.run(env=_stub(task, f"echo '.. [100%]'\necho '{DOCS}'\nexit 0\n"))
    assert res.returncode == 0, res.stdout
    line = _test_line(res.stdout)
    assert line.startswith("FLAG  named tests (2 path(s)): no counts"), line
    assert line.endswith("; pytest exit 0")
    assert "summary: 0 FAIL, 1 FLAG" in res.stdout
    assert "Docs" not in res.stdout
    res = task.run(env=_stub(task, "echo '1 failed, 2 passed in 0.1s'\nexit 1\n"))
    assert res.returncode == 1
    assert "2 passed, 1 failed, 0 skipped; pytest exit 1" in res.stdout


PYPROJECT = """\
[tool.pytest.ini_options]
addopts = "-ra -q"
"""

REAL_TESTS = """\
import warnings

import pytest


def test_warns():
    warnings.warn("a warning", UserWarning)


def test_ok():
    pass


@pytest.mark.skip(reason="skipped")
def test_skipped():
    pass
"""


@pytest.mark.parametrize("extra_q", ["", " -q"])
def test_real_pytest_with_quiet_addopts_and_a_warning(task: Task, extra_q: str) -> None:
    """End to end with this environment's pytest; " -q" is the old default (verbosity -2)."""
    task.commit(
        {"icenine_py/pyproject.toml": PYPROJECT, "icenine_py/tests/test_real.py": REAL_TESTS}
    )
    task.report(tests="tests/test_real.py")
    env = {"VERIFY_HANDBACK_PYTEST_CMD": f"{sys.executable} -m pytest -p no:cacheprovider{extra_q}"}
    res = task.run(env=env)
    assert res.returncode == 0, res.stdout
    line = _test_line(res.stdout)
    assert (
        line == "ok    named tests passed (1 path(s)): 2 passed, 0 failed, 1 skipped; pytest exit 0"
    )
    # a failing test: FAIL with the counts and pytest's exit code
    (task.repo / "icenine_py/tests/test_real.py").write_text(
        REAL_TESTS + "\n\ndef test_bad():\n" "    assert False\n"
    )
    res = task.run(env=env)
    assert res.returncode == 1
    assert "2 passed, 1 failed, 1 skipped; pytest exit 1" in _test_line(res.stdout)
