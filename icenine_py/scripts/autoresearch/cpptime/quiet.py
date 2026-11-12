"""Load sampler and quiet verdict of the timed runs (item 14, rule R4).

quiet.py sample --out FILE.jsonl --watch-pid PID [--period 30]
    logs the 1-minute load every 30 s and the top CPU processes every 300 s, until PID is gone
quiet.py verdict --samples FILE.jsonl --limit 2.5 --out FILE.json
    exit 0 if quiet, 1 if not: no load sample above the limit and no other process above 20% CPU
"""

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
import preflight as PF  # noqa: E402

BUSY_PCT = 20.0
LOAD_PERIOD_S = 30.0
PROC_PERIOD_S = 300.0
PRESET_LIMIT = 2.5
PREFLIGHT_LIMIT = 1.5


def choose_limit(idle: float, own: float) -> float:
    """R4's limit: 2.5 unless idle load plus the run's own load exceeds it, then the preflight's
    own limit plus the run's own load."""
    return PRESET_LIMIT if idle + own <= PRESET_LIMIT else PREFLIGHT_LIMIT + own


def verdict(
    load_samples: List[Dict[str, Any]], proc_samples: List[Dict[str, Any]], limit: float
) -> Dict[str, Any]:
    loads = [float(s["load"]) for s in load_samples]
    busy: List[Dict[str, Any]] = []
    for s in proc_samples:
        busy += [p for p in s.get("procs", []) if float(p["cpu"]) > BUSY_PCT]
    max_load = max(loads) if loads else float("nan")
    return dict(
        quiet=bool(loads) and max_load <= limit and not busy,
        limit=limit,
        max_load=max_load,
        median_load=sorted(loads)[len(loads) // 2] if loads else float("nan"),
        n_load=len(loads),
        busy=busy,
    )


def _descendants(pid: int) -> Set[int]:
    out = subprocess.run(["ps", "-Ao", "pid,ppid"], capture_output=True, text=True).stdout
    kids: Dict[int, List[int]] = {}
    for line in out.splitlines()[1:]:
        p = line.split()
        if len(p) == 2 and p[0].isdigit() and p[1].isdigit():
            kids.setdefault(int(p[1]), []).append(int(p[0]))
    seen, todo = {pid}, [pid]
    while todo:
        for c in kids.get(todo.pop(), []):
            if c not in seen:
                seen.add(c)
                todo.append(c)
    return seen


def _alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
        return True
    except OSError:
        return False


def sample(out: Path, watch_pid: int, period: float) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    t0, next_proc = time.time(), 0.0
    with out.open("w") as f:
        while True:
            t = time.time() - t0
            f.write(json.dumps(dict(kind="load", t=t, load=os.getloadavg()[0])) + "\n")
            if t >= next_proc:
                skip = _descendants(watch_pid) | PF._ancestor_pids()
                top = sorted(PF._sample_cpu().items(), key=lambda kv: -kv[1]["cpu"])[:8]
                procs = [
                    p
                    for pid, p in top
                    if pid not in skip and os.path.basename(p["command"]) not in PF.IGNORED_COMMANDS
                ]
                f.write(json.dumps(dict(kind="procs", t=t, procs=procs)) + "\n")
                next_proc = t + PROC_PERIOD_S
            f.flush()
            if not _alive(watch_pid):
                return
            time.sleep(period)


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("cmd", choices=["sample", "verdict"])
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--watch-pid", type=int, default=0)
    ap.add_argument("--period", type=float, default=LOAD_PERIOD_S)
    ap.add_argument("--samples", type=Path, default=None)
    ap.add_argument("--limit", type=float, default=PRESET_LIMIT)
    a = ap.parse_args(argv)
    if a.cmd == "sample":
        sample(a.out, a.watch_pid, a.period)
        return 0
    assert a.samples is not None
    rows = [json.loads(x) for x in a.samples.read_text().splitlines() if x.strip()]
    v = verdict(
        [r for r in rows if r["kind"] == "load"], [r for r in rows if r["kind"] == "procs"], a.limit
    )
    a.out.write_text(json.dumps(v, indent=1) + "\n")
    print(json.dumps(v))
    return 0 if v["quiet"] else 1


if __name__ == "__main__":
    sys.exit(main())
