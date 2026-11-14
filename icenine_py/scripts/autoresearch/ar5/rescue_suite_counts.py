"""AR5 step 9, phase B (second review fix 5): store the key-off full suite's counts as a tracked
source.

  uv run python scripts/autoresearch/ar5/rescue_suite_counts.py <junit.xml> <commit> [<note>]

``<junit.xml>`` is the file of the full suite run with the key off
(``uv run pytest tests/ --junitxml=<junit.xml>``); ``<commit>`` the commit whose library and tests
it ran; ``<note>`` how that commit is known. Reads the junit counts (tests, failures, errors,
skipped; passed = tests - failures - errors - skipped) and the test cases of
``tests/test_bfs_seed_rescue.py``. The junit file does not
hold the deselected count: it is the number of tests marked ``slow``, which ``pyproject.toml``
deselects (``addopts -m 'not slow'``), collected here with ``pytest --collect-only -m slow``.
Writes ``benchmarks/autoresearch/ar5/rescue/suite_counts.json``. No wall time is kept (the suite
ran with other jobs on the machine).
"""

import json
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict

ICE = Path(__file__).resolve().parents[3]
OUT = ICE / "benchmarks" / "autoresearch" / "ar5" / "rescue"


def slow_count() -> Dict[str, Any]:
    """``pytest --collect-only -m slow``: the tests ``addopts -m 'not slow'`` deselects."""
    r = subprocess.run(
        ["uv", "run", "pytest", "tests/", "--collect-only", "-m", "slow", "-p", "no:warnings"],
        cwd=ICE,
        capture_output=True,
        text=True,
    )
    line = [ln for ln in r.stdout.splitlines() if "collected" in ln][-1].strip()
    m = re.match(r"(\d+)/(\d+) tests collected", line)
    assert m, line
    line = re.sub(r" in [0-9.]+s.*$", "", line)
    return dict(line=line, slow=int(m.group(1)), all_tests=int(m.group(2)))


def main() -> None:
    junit, commit = Path(sys.argv[1]), sys.argv[2]
    root = ET.parse(junit).getroot()
    suite = root if root.tag == "testsuite" else root.find("testsuite")
    assert suite is not None
    n = {k: int(suite.get(k, 0)) for k in ("tests", "failures", "errors", "skipped")}
    cases = suite.findall("testcase")
    sc = slow_count()
    out: Dict[str, Any] = dict(
        command="uv run pytest tests/ --junitxml=junit.xml",
        commit=commit,
        commit_note=sys.argv[3] if len(sys.argv) > 3 else "",
        junit_timestamp=suite.get("timestamp", ""),
        key="SeedRescueCMA off (absent: the default)",
        junit=n,
        passed=n["tests"] - n["failures"] - n["errors"] - n["skipped"],
        skipped=n["skipped"],
        deselected=sc["slow"],
        deselected_source=f"pytest --collect-only -m slow: '{sc['line']}'",
        test_bfs_seed_rescue_cases=sum(
            1 for c in cases if c.get("classname", "").endswith("test_bfs_seed_rescue")
        ),
        timing_label="contended: other jobs ran on the machine; no time is a result",
    )
    assert out["junit"]["tests"] + out["deselected"] == sc["all_tests"], (out, sc)
    (OUT / "suite_counts.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
