"""AR3 step 4, part 3: start D4, P4 and E4 side by side, each killed at its cap.

  uv run python scripts/autoresearch/ar3/wrapped_bg.py [D4 P4 E4]

Caps (minutes): D4 150, P4 (with R4) 120, E4 90. A killed run is reported as not finished: a
``<kind>.killed`` marker is written beside the log and ``<kind>.done`` is written in both cases.
Logs: ``fine/cache/logs/<kind>.log``. Each child is started under ``caffeinate -i``.
"""

import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
LOGS = HERE.parents[2] / "benchmarks" / "autoresearch" / "ar3" / "fine" / "cache" / "logs"
CAPS_MIN = {"D4": 150, "P4": 120, "E4": 90}


def main() -> None:
    kinds = sys.argv[1:] or list(CAPS_MIN)
    LOGS.mkdir(parents=True, exist_ok=True)
    procs = {}
    t0 = time.time()
    for k in kinds:
        log = open(LOGS / f"{k}.log", "w")
        procs[k] = subprocess.Popen(
            [
                "caffeinate",
                "-i",
                "uv",
                "run",
                "python",
                "-u",
                str(HERE / "run_wrapped.py"),
                "--kind",
                k,
            ],
            stdout=log,
            stderr=subprocess.STDOUT,
            cwd=str(HERE),
            start_new_session=True,
        )
    left = dict(procs)
    while left:
        time.sleep(15)
        for k, p in list(left.items()):
            if p.poll() is not None:
                (LOGS / f"{k}.done").write_text(f"exit {p.returncode}\n")
                del left[k]
            elif time.time() - t0 > CAPS_MIN[k] * 60:
                import os
                import signal

                os.killpg(os.getpgid(p.pid), signal.SIGKILL)
                (LOGS / f"{k}.killed").write_text(f"killed at the cap of {CAPS_MIN[k]} min\n")
                (LOGS / f"{k}.done").write_text("killed\n")
                del left[k]


if __name__ == "__main__":
    main()
