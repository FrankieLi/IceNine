"""Machine-state preflight for timing runs: load, power source, thread settings, busy processes."""

import datetime
import os
import platform
import subprocess
import time
from typing import Any, Dict, List, Optional, Set


class MachineBusyError(RuntimeError):
    """Raised by ``require_quiet`` when the machine is not quiet enough for a timing run."""


def _run(cmd: List[str]) -> Optional[str]:
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=10, check=False)
        return out.stdout if out.returncode == 0 else None
    except (OSError, subprocess.SubprocessError):
        return None


def _power_source() -> str:
    """'ac', 'battery' or 'unknown' (macOS pmset; tolerates failure)."""
    out = _run(["pmset", "-g", "batt"])
    if not out:
        return "unknown"
    if "AC Power" in out:
        return "ac"
    if "Battery Power" in out:
        return "battery"
    return "unknown"


def _torch_threads() -> Optional[int]:
    try:
        import torch  # noqa: PLC0415

        return int(torch.get_num_threads())
    except Exception:  # noqa: BLE001 - torch absent or broken
        return None


IGNORED_COMMANDS = {
    "WindowServer",
    "mds",
    "mds_stores",
    "mdworker",
    "mdworker_shared",
    "kernel_task",
}
SAMPLE_GAP_S = 1.0


def _ancestor_pids() -> Set[int]:
    """This process and its ancestors (the uv / python wrappers that launched it)."""
    pids = {os.getpid()}
    out = _run(["ps", "-Ao", "pid,ppid"])
    if not out:
        return pids | {os.getppid()}
    parent = {}
    for line in out.splitlines()[1:]:
        parts = line.split()
        if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
            parent[int(parts[0])] = int(parts[1])
    pid = os.getpid()
    while pid in parent and parent[pid] > 1 and parent[pid] not in pids:
        pid = parent[pid]
        pids.add(pid)
    return pids


def _sample_cpu() -> Dict[int, Dict[str, Any]]:
    out = _run(["ps", "-Ao", "pid,pcpu,comm"])
    procs: Dict[int, Dict[str, Any]] = {}
    for line in (out or "").splitlines()[1:]:
        parts = line.split(None, 2)
        if len(parts) < 3:
            continue
        try:
            procs[int(parts[0])] = {
                "pid": int(parts[0]),
                "cpu": float(parts[1]),
                "command": parts[2],
            }
        except ValueError:
            continue
    return procs


def _busy_processes(threshold: float = 20.0) -> List[Dict[str, Any]]:
    """Other processes above ``threshold`` % CPU in two samples SAMPLE_GAP_S apart, excluding
    this process's ancestors and a default ignore list of system daemons."""
    first = _sample_cpu()
    time.sleep(SAMPLE_GAP_S)
    second = _sample_cpu()
    skip = _ancestor_pids()
    busy = []
    for pid, p in second.items():
        if pid in skip or os.path.basename(p["command"]) in IGNORED_COMMANDS:
            continue
        if p["cpu"] > threshold and first.get(pid, {"cpu": 0.0})["cpu"] > threshold:
            busy.append(p)
    return busy


def preflight() -> Dict[str, Any]:
    """Record the machine state relevant to a timing run (JSON-serialisable)."""
    try:
        load = list(os.getloadavg())
    except OSError:
        load = []
    return {
        "timestamp": datetime.datetime.now().isoformat(timespec="seconds"),
        "host": platform.node(),
        "loadavg": load,
        "cpu_count": os.cpu_count(),
        "power": _power_source(),
        "omp_num_threads": os.environ.get("OMP_NUM_THREADS"),
        "mkl_num_threads": os.environ.get("MKL_NUM_THREADS"),
        "torch_threads": _torch_threads(),
        "busy_processes": _busy_processes(),
    }


def require_quiet(
    max_load: float = 1.5, allow_battery: bool = False, info: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """Raise ``MachineBusyError`` if the 1-minute load exceeds ``max_load``, the machine is on
    battery (unless allowed), or another process uses more than 20% CPU. Returns the preflight
    dict (pass ``info`` to reuse one)."""
    info = info if info is not None else preflight()
    problems = []
    if info.get("loadavg") and info["loadavg"][0] > max_load:
        problems.append(f"1-min load {info['loadavg'][0]:.2f} > {max_load}")
    if info.get("power") == "battery" and not allow_battery:
        problems.append("running on battery power")
    for p in info.get("busy_processes", []):
        problems.append(f"busy process pid {p['pid']} ({p['cpu']:.0f}% CPU): {p['command']}")
    if problems:
        raise MachineBusyError("machine not quiet for timing: " + "; ".join(problems))
    return info
