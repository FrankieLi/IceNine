"""Read the logs of a C++ pilot run (cpptime).

  uv run python scripts/autoresearch/cpptime/parse_logs.py RUN_DIR [--time-file FILE] [--out FILE]

RUN_DIR holds ``ParallelReconLog_0.log`` (server), ``ParallelReconLog_1.log`` (client) and
``stderr.log`` (epoch-stamped by guard.py). Clock fields are whole seconds (``time_t``): a
single value below 1 s is not meaningful. Counts: a global search is one ``First Voxel
Required`` line; a rejected seed is a global search whose ``[Conf, HitRatio ]`` line is below the
acceptance and is not followed by the BFS line; a local optimization is one ``QueueSize`` line.
"""

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

NUM = r"(\d+(?:\.\d+)?(?:e[-+]?\d+)?)"
FIRST_VOXEL = re.compile(r"First Voxel Required:\s*" + NUM + r"\s*Sec")
TOTAL_PROC = re.compile(r"Total Processing Time: \(sec\)\s*" + NUM)
SERVER_TOTAL = re.compile(r"Total Time elapsed[^:]*:\s*" + NUM)
EPOCH_LINE = re.compile(r"^([0-9]+\.[0-9]+) (.*)$")
BFS_LINE = "Fitting neighbors via Breadth First Search"


def parse_client_log(text: str) -> Dict[str, Any]:
    """Global searches, local optimizations and the client timers of a client log."""
    lines = text.splitlines()
    first = [float(m.group(1)) for ln in lines if (m := FIRST_VOXEL.search(ln))]
    n_local = sum(1 for ln in lines if ln.startswith("QueueSize"))
    total = [float(m.group(1)) for ln in lines if (m := TOTAL_PROC.search(ln))]
    # a global search without the BFS line after it is a rejected seed
    rejected, pending = 0, False
    for ln in lines:
        if FIRST_VOXEL.search(ln):
            if pending:
                rejected += 1
            pending = True
        elif BFS_LINE in ln:
            pending = False
    rejected += 1 if pending else 0
    return dict(
        global_searches=len(first),
        first_voxel_s=first,
        global_s_sum=float(sum(first)),
        local_optimizations=n_local,
        rejected_seeds=rejected,
        total_processing_s=total[-1] if total else None,
    )


def parse_server_log(text: str) -> Dict[str, Any]:
    m = SERVER_TOTAL.findall(text)
    return dict(total_time_s=float(m[-1]) if m else None)


def parse_stderr(text: str) -> Dict[str, Any]:
    """Epoch of the image-read markers and the read time between them."""
    marks: Dict[str, float] = {}
    first_epoch: Optional[float] = None
    for ln in text.splitlines():
        m = EPOCH_LINE.match(ln)
        if not m:
            continue
        t, rest = float(m.group(1)), m.group(2)
        first_epoch = t if first_epoch is None else first_epoch
        for key, tag in (
            ("before_read", "Before reading experimental data"),
            ("after_read", "After reading experimental data"),
        ):
            if tag in rest and key not in marks:
                marks[key] = t
    out: Dict[str, Any] = dict(first_epoch=first_epoch)
    out["before_read_epoch"] = marks.get("before_read")
    out["after_read_epoch"] = marks.get("after_read")
    if "before_read" in marks and "after_read" in marks:
        out["read_s"] = marks["after_read"] - marks["before_read"]
    return out


def unit_costs(
    global_s: float, bfs_s: float, n_global: int, n_local: int, n_evals: Optional[int]
) -> Dict[str, Optional[float]]:
    """Seconds per global search, per local optimization, microseconds per evaluation.

    The local-optimization time is the BFS time minus the global-search time."""
    local_s = bfs_s - global_s
    return dict(
        s_per_global=global_s / n_global if n_global else None,
        s_per_local=local_s / n_local if n_local else None,
        us_per_eval=1e6 * bfs_s / n_evals if n_evals else None,
    )


def parse_time_file(text: str) -> Dict[str, float]:
    """``/usr/bin/time -p`` output: real, user, sys."""
    out: Dict[str, float] = {}
    for ln in text.splitlines():
        p = ln.split()
        if len(p) == 2 and p[0] in ("real", "user", "sys"):
            out[p[0]] = float(p[1])
    return out


def parse_run(run: Path, time_file: Optional[Path] = None) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    out.update(parse_client_log((run / "ParallelReconLog_1.log").read_text()))
    out.update(
        {
            "server_" + k: v
            for k, v in parse_server_log((run / "ParallelReconLog_0.log").read_text()).items()
        }
    )
    se = run / "stderr.log"
    if se.exists():
        out.update(parse_stderr(se.read_text()))
    if time_file is not None and time_file.exists():
        out.update({"time_" + k: v for k, v in parse_time_file(time_file.read_text()).items()})
    return out


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run", type=Path)
    ap.add_argument("--time-file", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args(argv)
    r = parse_run(a.run, a.time_file)
    text = json.dumps(r, indent=1) + "\n"
    if a.out:
        a.out.write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
