"""Load sampler and quiet verdict of the timed runs (item 14, rule R4).

quiet.py sample --out FILE.jsonl --watch-pid PID [--period 30]
    logs the 1-minute load every 30 s and the top CPU processes every 300 s, until PID is gone
quiet.py verdict --samples FILE.jsonl --limit 2.5 --out FILE.json
    strict quiet: no load sample above the limit and no other process above 20% CPU.
    Exit 0 strictly quiet; 3 not strict but relaxed_ok (R4b, D-20261114-N1: median load within
    the limit and every busy process a macOS daemon under /usr/libexec/, /System/Library/ or
    /sbin/);
    1 neither. The JSON also holds relaxed_ok and label.
quiet.py pair --arm A --runs A1 A2 --cache DIR
    R4b pair check: if a run of the pair carried a label, the real times must agree within 2%;
    writes DIR/pair_<arm>.json; exit 0 reported, 6 no time reported for the arm.
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
    median = sorted(loads)[len(loads) // 2] if loads else float("nan")  # upper median, n even
    quiet = bool(loads) and max_load <= limit and not busy
    # R4b (D-20261114-N1): median load at most 2.5 (PRESET_LIMIT, whatever the strict limit is)
    # and every busy process a macOS daemon
    relaxed_ok = bool(loads) and median <= PRESET_LIMIT and all(_is_daemon(p) for p in busy)
    names = sorted({os.path.basename(p["command"]) for p in busy})
    n_over = sum(x > limit for x in loads)
    parts = []
    if n_over:
        parts.append(f"load max {max_load:.2f} > {limit} in {n_over} of {len(loads)} samples")
    parts.append(f"busy: {', '.join(names)}" if names else "no process over 20% CPU caught")
    label = "" if quiet else f"strict quiet check failed ({'; '.join(parts)})"
    return dict(
        quiet=quiet,
        limit=limit,
        max_load=max_load,
        median_load=median,
        n_load=len(loads),
        busy=busy,
        relaxed_ok=relaxed_ok,
        label=label,
    )


DAEMON_PREFIXES = ("/usr/libexec/", "/System/Library/", "/sbin/")  # not /System/Applications/
PAIR_TOLERANCE = 0.02


def _is_daemon(p: Dict[str, Any]) -> bool:
    return str(p["command"]).startswith(DAEMON_PREFIXES)


def exit_code(v: Dict[str, Any]) -> int:
    """0 strictly quiet, 3 not strict but relaxed_ok (R4b), 1 neither."""
    return 0 if v["quiet"] else 3 if v["relaxed_ok"] else 1


def pair_report(arm: str, t1: float, t2: float, label1: str, label2: str) -> Dict[str, Any]:
    """R4b (iii): if either run of a repeat pair carried a label, the real times must agree
    within 2% (difference over the smaller time), else no time is reported for the arm."""
    labels = [x for x in (label1, label2) if x]
    rel = abs(t1 - t2) / min(t1, t2)
    ok = not labels or rel <= PAIR_TOLERANCE
    return dict(
        arm=arm,
        real_s=[t1, t2],
        rel_diff=rel,
        labels=labels,
        reported=ok,
        note="" if ok else f"no time reported for arm {arm}",
    )


def _real_s(path: Path) -> float:
    for line in path.read_text().splitlines():
        if line.startswith("real"):
            return float(line.split()[1])
    raise ValueError(f"no real time in {path}")


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
    ap.add_argument("cmd", choices=["sample", "verdict", "pair"])
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--arm", default="")
    ap.add_argument("--runs", nargs=2, default=None, help="pair: the two run names (A1 A2)")
    ap.add_argument("--cache", type=Path, default=None, help="pair: the item14 cache folder")
    ap.add_argument("--watch-pid", type=int, default=0)
    ap.add_argument("--period", type=float, default=LOAD_PERIOD_S)
    ap.add_argument("--samples", type=Path, default=None)
    ap.add_argument("--limit", type=float, default=PRESET_LIMIT)
    a = ap.parse_args(argv)
    if a.cmd == "pair":
        assert a.cache is not None and a.runs is not None
        n1, n2 = a.runs
        lab = [json.loads((a.cache / "quiet" / f"{n}.json").read_text())["label"] for n in a.runs]
        t = [_real_s(a.cache / "logs" / f"{n}.time") for n in a.runs]
        p = pair_report(a.arm, t[0], t[1], lab[0], lab[1])
        p["runs"] = [n1, n2]
        (a.cache / f"pair_{a.arm}.json").write_text(json.dumps(p, indent=1) + "\n")
        print(json.dumps(p))
        return 0 if p["reported"] else 6
    assert a.out is not None
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
    return exit_code(v)


if __name__ == "__main__":
    sys.exit(main())
