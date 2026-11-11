"""AR3 step 7: the four whole-layer runs on the sample ``big``, each with a coded time cap.

  big_run.py --kind D      # stock BFS, grid search (cap 180 min)
  big_run.py --kind D4     # step 5's run (cap 90 min)
  big_run.py --kind S5     # D4 with step 6's seed step (cap 90 min)
  big_run.py --kind P4     # paint, then R4 (cap 90 min); the pass scripts make the map P4R4
  ... --smoke N            # pilot: first N voxels (a BFS may still flood the whole layer)

This script only sets ``AR3_SAMPLE=big`` and calls the existing run functions of steps 3, 4 and 6
(``run_gradient.run_d``, ``run_wrapped.run_bfs`` / ``run_p4``, ``run_seed.run``) unchanged. The
cap is a watchdog thread: when it fires, it writes ``capped.json`` (elapsed time, the run's
counters and the number of voxels with a record at that moment) beside the run's other files and
ends the process; the run is not repeated. A smoke run writes to ``big/smoke/`` instead.
"""

import os

os.environ["AR3_SAMPLE"] = "big"

import argparse  # noqa: E402
import json  # noqa: E402
import threading  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any, Dict  # noqa: E402

import ar3_common as C  # noqa: E402
import run_gradient as RG  # noqa: E402
import run_seed as RS  # noqa: E402
import run_wrapped as RW  # noqa: E402

CAP_MIN = {"D": 180, "D4": 90, "S5": 90, "P4": 90}
state: Dict[str, Any] = {}


def patch_build(smoke: int) -> None:
    """Keep the run's reconstructor (for the watchdog); in a smoke run limit the voxel count."""
    orig_build = C.build

    def build(*a: Any, **k: Any) -> Any:
        setup, bfs, t = orig_build(*a, **k)
        state["bfs"] = bfs
        if smoke:
            rs = bfs.reconstruct_sample

            def limited(*aa: Any, **kk: Any) -> Any:
                kk["max_voxels"] = smoke
                return rs(*aa, **kk)

            bfs.reconstruct_sample = limited
        return setup, bfs, t

    C.build = build
    for m in (RG, RW, RS):  # modules that call ``C.build`` through the module object
        assert m.C is C


def watchdog(out: Path, kind: str, cap_s: float, t0: float) -> None:
    def loop() -> None:
        time.sleep(cap_s)
        bfs = state.get("bfs")
        st = {}
        n_rec = 0
        if bfs is not None:
            st = {k: v for k, v in bfs.stats.items() if isinstance(v, (int, float))}
            n_rec = len(bfs.records)
        (out / "capped.json").write_text(
            json.dumps(
                dict(
                    kind=kind,
                    cap_s=cap_s,
                    elapsed_s=time.time() - t0,
                    voxels_with_record=n_rec,
                    stats=st,
                ),
                default=float,
            )
            + "\n"
        )
        print(f"{kind}: CAP of {cap_s} s reached; stopped", flush=True)
        os._exit(3)

    threading.Thread(target=loop, daemon=True).start()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kind", required=True, choices=sorted(CAP_MIN))
    ap.add_argument("--smoke", type=int, default=0)
    ap.add_argument("--cap-min", type=float, default=0.0, help="override the cap (tests only)")
    a = ap.parse_args()
    base = C.OUT / "smoke" if a.smoke else C.OUT
    if a.smoke:
        RG.RUNS = base / "gradient_runs"
        RW.FRUNS = base / "fine_runs"
        RS.RUNS = base / "seed_runs"
    patch_build(a.smoke)
    cap_s = 60.0 * (a.cap_min or CAP_MIN[a.kind])
    t0 = time.time()
    if a.kind == "D":
        out = RG.RUNS / "D"
        out.mkdir(parents=True, exist_ok=True)
        watchdog(out, a.kind, cap_s, t0)
        RG.run_d(argparse.Namespace(no_timing_gate=False))
    elif a.kind == "D4":
        out = RW.FRUNS / "D4"
        out.mkdir(parents=True, exist_ok=True)
        watchdog(out, a.kind, cap_s, t0)
        RW.run_bfs("D4", out, a.smoke or None)
    elif a.kind == "S5":
        out = RS.RUNS / "S5"
        out.mkdir(parents=True, exist_ok=True)
        watchdog(out, a.kind, cap_s, t0)
        RS.run("S5", out, a.smoke or None)
    else:
        out = RW.FRUNS / "P4"
        out.mkdir(parents=True, exist_ok=True)
        watchdog(out, a.kind, cap_s, t0)
        RW.run_p4(out)


if __name__ == "__main__":
    main()
