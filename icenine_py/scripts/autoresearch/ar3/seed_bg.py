"""AR3 step 6, part 2: start S and K2 side by side, each killed at its cap (60 minutes).

  uv run python scripts/autoresearch/ar3/seed_bg.py [S K2]

A killed run is reported as not finished: ``<kind>.killed`` is written beside the log, and
``<kind>.done`` in both cases. Logs: ``seed/cache/logs/<kind>.log``. Each child runs under
``caffeinate -dimsu``. The cap counts awake time only: a gap of more than 60 s between two checks
(the loop sleeps 15 s) is taken as sleep and not counted (as in ``keys_bg.py``).
"""

import os
import signal
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
LOGS = HERE.parents[2] / "benchmarks" / "autoresearch" / "ar3" / "seed" / "cache" / "logs"
CAPS_MIN = {"S": 60, "K2": 60, "S5": 60}


def main() -> None:
    kinds = sys.argv[1:] or list(CAPS_MIN)
    LOGS.mkdir(parents=True, exist_ok=True)
    procs = {}
    for k in kinds:
        log = open(LOGS / f"{k}.log", "w")
        procs[k] = subprocess.Popen(
            [
                "caffeinate",
                "-dimsu",
                "uv",
                "run",
                "python",
                "-u",
                str(HERE / "run_seed.py"),
                "--kind",
                k,
            ],
            stdout=log,
            stderr=subprocess.STDOUT,
            cwd=str(HERE),
            start_new_session=True,
        )
    left = dict(procs)
    awake = 0.0
    last = time.time()
    while left:
        time.sleep(15)
        now = time.time()
        gap = now - last
        last = now
        if gap <= 60:
            awake += gap
        for k, p in list(left.items()):
            if p.poll() is not None:
                (LOGS / f"{k}.done").write_text(f"exit {p.returncode}\n")
                del left[k]
            elif awake > CAPS_MIN[k] * 60:
                os.killpg(os.getpgid(p.pid), signal.SIGKILL)
                (LOGS / f"{k}.killed").write_text(f"killed at the cap of {CAPS_MIN[k]} min\n")
                (LOGS / f"{k}.done").write_text("killed\n")
                del left[k]


if __name__ == "__main__":
    main()
