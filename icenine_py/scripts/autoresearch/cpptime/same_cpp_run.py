"""Identity check of two C++ pilot runs (cpptime).

  uv run python scripts/autoresearch/cpptime/same_cpp_run.py RUN_A RUN_B

Equal means: ``cpp_pilot.mic.LBFS`` byte for byte, and the client log
(``ParallelReconLog_1.log``) equal once its time fields are dropped. Exit 0 if equal, 1 with the
first difference named.
"""

import re
import sys
from pathlib import Path
from typing import List, Optional

MIC = "cpp_pilot.mic.LBFS"
CLIENT_LOG = "ParallelReconLog_1.log"
# time fields: whole-second timers and the start-time stamp
NUM = r"[-+]?(?:\d+(?:\.\d+)?(?:e[-+]?\d+)?|nan|inf)"
TIME_PATTERNS = [
    (re.compile(r"^StartTime:.*$", re.M), "StartTime: <t>"),
    (re.compile(r"(First Voxel Required:)\s*" + NUM), r"\1 <t>"),
    (re.compile(r"(Time elapsed for this set of voxels \(sec\):)\s*" + NUM), r"\1 <t>"),
    (re.compile(r"(Time per voxel \(sec\):)\s*" + NUM), r"\1 <t>"),
    (re.compile(r"(Time spent waiting for command \(sec\):)\s*" + NUM), r"\1 <t>"),
    (re.compile(r"(Total Processing Time: \(sec\))\s*" + NUM), r"\1 <t>"),
]


def strip_time(text: str) -> List[str]:
    """The log lines with every time field replaced by ``<t>``."""
    for pat, rep in TIME_PATTERNS:
        text = pat.sub(rep, text)
    return text.splitlines()


def diff_runs(a: Path, b: Path) -> Optional[str]:
    """First difference between two run folders, or None if they are identical."""
    ma, mb = (a / MIC).read_bytes(), (b / MIC).read_bytes()
    if ma != mb:
        return f"{MIC} differs ({len(ma)} vs {len(mb)} bytes)"
    la, lb = strip_time((a / CLIENT_LOG).read_text()), strip_time((b / CLIENT_LOG).read_text())
    if len(la) != len(lb):
        return f"client log length {len(la)} vs {len(lb)} lines"
    for i, (x, y) in enumerate(zip(la, lb)):
        if x != y:
            return f"client log line {i + 1}: {x!r} vs {y!r}"
    return None


def main() -> int:
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    d = diff_runs(Path(sys.argv[1]), Path(sys.argv[2]))
    print("identical" if d is None else f"DIFFERENT: {d}")
    return 0 if d is None else 1


if __name__ == "__main__":
    sys.exit(main())
