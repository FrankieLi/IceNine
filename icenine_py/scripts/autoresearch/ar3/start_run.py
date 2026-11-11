"""AR3 step 9, part 3: the whole-layer runs D4 and S5 on a sample, each with a coded time cap.

  RESEED_SAMPLE=mid RESEED_CACHE=<start/cache> \\
      uv run python scripts/autoresearch/ar3/start_run.py --kind D4|S5    # cap 90 minutes
  ... --smoke N       # pilot: first N voxels (a BFS may still flood the whole layer)

The run functions are step 7's (``big_run.py``): ``run_wrapped.run_bfs("D4", ...)`` and
``run_seed.run("S5", ...)``, unchanged. ``S8`` is ``reseed_s8.py``. The cap is a watchdog thread
that writes ``capped.json`` and ends the process; a run that hits it is not repeated. Writes to
``<RESEED_CACHE>/<kind>_<sample>/`` (``recon.mic.gz``, ``records.json.gz``, ``run.json`` ...).
"""

import argparse
import json
import os
import threading
import time
from pathlib import Path
from typing import Any, Dict

import reseed_common as RC

import ar3_common as C  # noqa: E402
import run_seed as RS  # noqa: E402
import run_wrapped as RW  # noqa: E402

CAP_MIN = 90
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
    for m in (RW, RS):
        assert m.C is C


def watchdog(out: Path, kind: str, cap_s: float, t0: float) -> None:
    def loop() -> None:
        time.sleep(cap_s)
        bfs = state.get("bfs")
        st: Dict[str, Any] = {}
        n_rec = 0
        if bfs is not None:
            st = {k: v for k, v in bfs.stats.items() if isinstance(v, (int, float))}
            n_rec = len(bfs.records)
        (out / "capped.json").write_text(
            json.dumps(
                dict(kind=kind, cap_s=cap_s, elapsed_s=time.time() - t0, voxels=n_rec, stats=st),
                default=float,
            )
            + "\n"
        )
        print(f"{kind}: CAP of {cap_s} s reached; stopped", flush=True)
        os._exit(3)

    threading.Thread(target=loop, daemon=True).start()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kind", required=True, choices=["D4", "S5"])
    ap.add_argument("--smoke", type=int, default=0)
    a = ap.parse_args()
    name = f"{a.kind}_{RC.SAMPLE}" + ("_smoke" if a.smoke else "")
    out = RC.CACHE / name
    out.mkdir(parents=True, exist_ok=True)
    patch_build(a.smoke)
    t0 = time.time()
    watchdog(out, a.kind, CAP_MIN * 60.0, t0)
    if a.kind == "D4":
        RW.run_bfs("D4", out, a.smoke or None)
    else:
        RS.run("S5", out, a.smoke or None)
    (out / f"{a.kind}.done").write_text(f"{time.time() - t0:.1f}\n")


if __name__ == "__main__":
    main()
