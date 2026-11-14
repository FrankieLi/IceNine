"""AR5 step 5 (fixes round): store the full test suite's counts as a tracked source.

  uv run pytest tests/ --junitxml=<scratch>/junit.xml > <scratch>/suite.log 2>&1
  uv run python scripts/autoresearch/ar5/keys_suite_counts.py <scratch>/junit.xml \\
      <scratch>/suite.log [<commit the suite ran at>]

Reads the junit file (tests, failures, errors, skipped, and the test cases of
``tests/test_bfs_ar5_keys.py``) and the last line of the pytest log (the line with the deselected
count, which the junit file does not hold). Writes
``benchmarks/autoresearch/ar5/keys/suite_counts.json``. The suite time is contended (another long
job ran on the machine) and is not a result.
"""

import json
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict

OUT = Path(__file__).resolve().parents[3] / "benchmarks" / "autoresearch" / "ar5" / "keys"


def main() -> None:
    junit, log = Path(sys.argv[1]), Path(sys.argv[2])
    root = ET.parse(junit).getroot()
    suite = root if root.tag == "testsuite" else root.find("testsuite")
    assert suite is not None
    cases = suite.findall("testcase")
    lines = [ln for ln in log.read_text().splitlines() if ln.strip()]
    last = re.sub(r"^=+ | =+$", "", lines[-1]).strip()
    last = re.sub(r" in [0-9.]+s.*$", "", last)  # the wall time is contended, not a result
    # the commit the suite ran at (third argument), else HEAD
    commit = (
        sys.argv[3]
        if len(sys.argv) > 3
        else subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True, check=True
        ).stdout.strip()
    )
    out: Dict[str, Any] = dict(
        command="uv run pytest tests/ --junitxml=junit.xml",
        commit=commit,
        pytest_count_line=last,
        junit=dict(
            tests=int(suite.get("tests", 0)),
            failures=int(suite.get("failures", 0)),
            errors=int(suite.get("errors", 0)),
            skipped=int(suite.get("skipped", 0)),
        ),
        test_bfs_ar5_keys_cases=sum(
            1 for c in cases if c.get("classname", "").endswith("test_bfs_ar5_keys")
        ),
        timing_label="contended: another long job ran on the machine; no time is a result",
    )
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "suite_counts.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
