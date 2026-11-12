"""Run a command under the stop rules of the cpptime plan (IceNine exits 0 on every failure).

  guard.py --limit-s 3600 --stall-s 600 --stderr FILE --watch GLOB [--watch GLOB] -- CMD ...

stdout of CMD goes to the caller's stdout; stderr goes to --stderr, each line prefixed with the
epoch seconds it arrived (the image-read time is read from these). Kills the whole
process group and exits with a code of its own when:
  124  the run passed --limit-s;
  125  no watched log grew for --stall-s (a blocked rank);
  126  stderr (or a watched log) holds RUNTIME_ASSERT or "[ConfigFile] Error".
Otherwise exits with CMD's code.
"""

import argparse
import glob
import os
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import IO, List, Optional

BAD = ("RUNTIME_ASSERT", "[ConfigFile] Error")
LIMIT, STALL, FAILED = 124, 125, 126


def log_size(patterns: List[str]) -> int:
    """Total size in bytes of the files matching the patterns."""
    return sum(os.path.getsize(f) for p in patterns for f in glob.glob(p) if os.path.isfile(f))


def has_bad(path: Path) -> Optional[str]:
    """The first failure marker found in the file, or None."""
    if not path.exists():
        return None
    text = path.read_text(errors="replace")
    return next((b for b in BAD if b in text), None)


def stamp_lines(src: IO[bytes], dst: IO[bytes]) -> None:
    """Copy src to dst line by line, each line prefixed with the epoch seconds it arrived."""
    for line in iter(src.readline, b""):
        dst.write(b"%.3f " % time.time() + line)
        dst.flush()


def kill_group(proc: "subprocess.Popen[bytes]") -> None:
    try:
        os.killpg(proc.pid, signal.SIGTERM)
        time.sleep(2)
        os.killpg(proc.pid, signal.SIGKILL)
    except (ProcessLookupError, PermissionError):
        pass


def run(
    cmd: List[str],
    limit_s: float,
    stall_s: float,
    stderr_path: Path,
    watch: List[str],
    poll_s: float = 1.0,
    cwd: Optional[Path] = None,
) -> int:
    stderr_path.parent.mkdir(parents=True, exist_ok=True)
    with stderr_path.open("wb") as err:
        proc = subprocess.Popen(cmd, stderr=subprocess.PIPE, start_new_session=True, cwd=cwd)
        assert proc.stderr is not None
        reader = threading.Thread(target=stamp_lines, args=(proc.stderr, err), daemon=True)
        reader.start()
        t0 = time.time()
        last_size, last_change = log_size(watch), time.time()
        while proc.poll() is None:
            time.sleep(poll_s)
            now = time.time()
            if now - t0 > limit_s:
                print(f"guard: passed the {limit_s:.0f} s limit", file=sys.stderr)
                kill_group(proc)
                return LIMIT
            size = log_size(watch)
            if size != last_size:
                last_size, last_change = size, now
            elif watch and now - last_change > stall_s:
                print(f"guard: no log grew for {stall_s:.0f} s", file=sys.stderr)
                kill_group(proc)
                return STALL
            bad = has_bad(stderr_path)
            if bad:
                print(f"guard: stderr holds {bad}", file=sys.stderr)
                kill_group(proc)
                return FAILED
        rc = proc.returncode
        reader.join(timeout=5)
    bad = has_bad(stderr_path)
    if bad:
        print(f"guard: stderr holds {bad}", file=sys.stderr)
        return FAILED
    return rc


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--limit-s", type=float, required=True)
    ap.add_argument("--stall-s", type=float, default=600.0)
    ap.add_argument("--stderr", type=Path, required=True)
    ap.add_argument("--watch", action="append", default=[])
    ap.add_argument("--cwd", type=Path, default=None, help="working directory of CMD")
    ap.add_argument("cmd", nargs=argparse.REMAINDER)
    a = ap.parse_args()
    cmd = a.cmd[1:] if a.cmd and a.cmd[0] == "--" else a.cmd
    stderr = a.stderr.resolve()
    watch = [str((a.cwd / w) if a.cwd else w) for w in a.watch]
    return run(cmd, a.limit_s, a.stall_s, stderr, watch, cwd=a.cwd)


if __name__ == "__main__":
    sys.exit(main())
