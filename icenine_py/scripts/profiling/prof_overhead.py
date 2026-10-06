#!/usr/bin/env python3
"""Overhead of the stage wrappers (Instrument): (1) per-call cost of the `evaluate` wrapper on a
dummy method (so the cost of the wrapper alone, without the evaluation), as a share of a U0 run
(number of evaluate calls per run x per-call overhead / wall time); (2) one paired check on a
real run: plain vs wrapped baseline reconstruct_voxel of the same voxel, alternating, 2 + 2 runs.

  uv run python scripts/profiling/prof_overhead.py
"""

import sys
import time
from pathlib import Path
from typing import Any, Dict

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import prof_common as PC  # noqa: E402


class _Dummy:
    pixel_radius = 0

    def evaluate(self, x: int) -> int:
        return x + 1


def per_call_overhead(n: int = 400_000) -> Dict[str, float]:
    T = PC.ST.StageTimer()
    inst = PC.Instrument(T, 3500, net=False)
    d = _Dummy()

    def loop() -> float:
        t0 = time.perf_counter()
        for i in range(n):
            d.evaluate(i)
        return (time.perf_counter() - t0) / n

    plain = min(loop() for _ in range(3))
    T.patch(_Dummy, "evaluate", inst._eval_name)
    try:
        wrapped = min(loop() for _ in range(3))
    finally:
        T.uninstall()
    return dict(plain_s=plain, wrapped_s=wrapped, overhead_s=wrapped - plain, n=n)


def paired_run() -> Dict[str, Any]:
    import prof_u0 as U

    C = U.C
    U.init_worker(C.worker_args())  # installs the wrappers permanently
    from optimizer_sweep import voxel_context

    W = C.get_worker()
    T, I = U._S["T"], U._S["I"]  # noqa: E741
    v, vpos, var = U.B.voxels()[0], 0, "clean"
    C.attach(np.load(C.CACHE_DIR / "images" / f"v{v}_{var}.npz")["keys"])
    vctx = voxel_context(W.ctx, v)
    walls: Dict[str, list] = {"plain": [], "wrapped": []}
    ref = None
    for mode in ("wrapped", "plain", "wrapped", "plain", "wrapped"):
        if mode == "plain" and I._installed:
            I.uninstall()
        elif mode == "wrapped" and not I._installed:
            I.install()
        t0 = time.perf_counter()
        with C.quiet():
            res = W.rec.reconstruct_voxel(vctx.vertices, vctx.voxel.phase, rng=C.run_seed(vpos, 0))
        walls[mode].append(time.perf_counter() - t0)
        R = np.asarray(res.orientation)
        ref = R if ref is None else ref
        assert np.array_equal(R, ref)
        n_ev = sum(W.rec.last_eval_counts[:2])
    return dict(walls=walls, evaluate_calls=int(n_ev), voxel=int(v), variant=var)


def main() -> None:
    PC.OUT.mkdir(parents=True, exist_ok=True)
    out: Dict[str, Any] = dict(isolation=PC.isolation_record("overhead"))
    out["per_call"] = per_call_overhead()
    out["paired"] = paired_run()
    w = out["paired"]["walls"]
    out["paired"]["median_plain_s"] = float(np.median(w["plain"]))
    out["paired"]["median_wrapped_s"] = float(np.median(w["wrapped"]))
    out["paired"]["wrapped_over_plain"] = (
        out["paired"]["median_wrapped_s"] / out["paired"]["median_plain_s"]
    )
    calls = out["paired"]["evaluate_calls"]
    out["share_of_run_from_per_call"] = (
        calls * out["per_call"]["overhead_s"] / out["paired"]["median_plain_s"]
    )
    PC.write_json(PC.OUT / "overhead.json", out)
    print(
        f"wrapper per call {out['per_call']['overhead_s'] * 1e6:.2f} us; {calls} evaluate calls "
        f"-> {100 * out['share_of_run_from_per_call']:.2f}% of a run; paired run wrapped/plain "
        f"{out['paired']['wrapped_over_plain']:.4f} "
        f"({walls_str(w)})"
    )


def walls_str(w: Dict[str, list]) -> str:
    return (
        f"plain {[round(x, 2) for x in w['plain']]}, wrapped {[round(x, 2) for x in w['wrapped']]}"
    )


if __name__ == "__main__":
    main()
