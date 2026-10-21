"""S2 Step 0c: us per call of the fused prototype against the torch ``evaluate``.

S1's per-call set (``prof_eval``: ThreeVoxels voxel 0, ``Example2.Simulation.config``, 400
orientations within 3 degrees of the truth, seed 0). Arms: the Q-max 8 triangle cost and the radius-3
Q-max 5 screen (centroid), each with the default rule and with ``CppOverlapRule``. Implementations:
``torch`` (the library ``evaluate``), ``fused`` (``s2fused.FusedEval``, portable C, building the
``OverlapInfo``), ``fused_acc`` (the Accelerate-sgemm build), ``c_only`` (the bare C call of the
portable build, no ``OverlapInfo``: the Python share of ``fused``). Each implementation has its own
cost function object; one untimed warm-up pass fills the frames.

Single worker (torch threads 1), preflight-gated (``preflight.require_quiet``; ``--allow-busy`` runs
anyway and labels the timings "contended"). ``--rounds`` interleaved rounds (the implementation
order rotates each round); the per-call time of an (arm, implementation) is the median over rounds
of one pass over the 400 orientations.

  uv run python scripts/speed_parity/s2/s2_timing.py --out benchmarks/s2_step0
"""

import argparse
import json
import statistics
import sys
import time
from pathlib import Path
from typing import Any, Callable, Dict, List

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
sys.path.insert(0, str(ICE / "scripts" / "common"))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))

import preflight  # noqa: E402

ARMS = {
    "triangle_q8_default": (0, 8.0, False),
    "triangle_q8_cpp": (0, 8.0, True),
    "screen_r3_q5_default": (3, 5.0, False),
    "screen_r3_q5_cpp": (3, 5.0, True),
}
IMPLS = ("torch", "fused", "fused_acc", "c_only")


def gate(allow_busy: bool) -> Dict[str, Any]:
    info = preflight.preflight()
    try:
        preflight.require_quiet(info=info)
        info["quiet"] = True
    except preflight.MachineBusyError as exc:
        if not allow_busy:
            raise
        info["quiet"] = False
        info["not_quiet_reason"] = str(exc)
        print("WARNING: machine not quiet, timings are 'contended':", exc, flush=True)
    info["label"] = "single-worker" if info["quiet"] else "contended"
    return info


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--rounds", type=int, default=5)
    ap.add_argument("--n", type=int, default=400, help="orientations per pass (pilot: fewer)")
    ap.add_argument("--allow-busy", action="store_true")
    a = ap.parse_args()
    out = Path(a.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    info = gate(a.allow_busy)
    (out / "preflight_timing.json").write_text(json.dumps(info, indent=1) + "\n")

    import prof_eval as P  # chdirs into the example, builds the physics, torch threads 1
    import s2fused as SF

    torch.set_num_threads(1)
    orients = P.orients[: a.n]
    verts, phase = P.verts, P.vox.phase
    calls: Dict[str, Dict[str, Callable[[], None]]] = {}
    for arm, args in ARMS.items():
        cf_t = P.make(*args)
        fe = SF.FusedEval(P.make(*args), "portable")
        fa = SF.FusedEval(P.make(*args), "accelerate")
        fr = SF.FusedEval(P.make(*args), "portable")
        raw = SF.ext("portable").evaluate

        def run_torch(cf: Any = cf_t) -> None:
            for o in orients:
                cf.evaluate(o, verts, phase)

        def run_fused(f: Any = fe) -> None:
            for o in orients:
                f.evaluate(o, verts, phase)

        def run_acc(f: Any = fa) -> None:
            for o in orients:
                f.evaluate(o, verts, phase)

        for o in orients:  # fill the frames of the bare-call context
            fr.evaluate(o, verts, phase)
        cap = fr.context(phase).capsule
        vn = fr._verts(verts)

        def run_raw(c: Any = cap, v: np.ndarray = vn) -> None:
            for o in orients:
                raw(c, o, v)

        calls[arm] = dict(torch=run_torch, fused=run_fused, fused_acc=run_acc, c_only=run_raw)
        for fn in calls[arm].values():
            fn()  # warm-up: frames filled, contexts built
    per: Dict[str, Dict[str, List[float]]] = {arm: {i: [] for i in IMPLS} for arm in ARMS}
    for rnd in range(a.rounds):
        order = list(IMPLS[rnd % len(IMPLS) :] + IMPLS[: rnd % len(IMPLS)])
        for arm in ARMS:
            for impl in order:
                t0 = time.perf_counter()
                calls[arm][impl]()
                per[arm][impl].append(1e6 * (time.perf_counter() - t0) / len(orients))
        print(
            f"round {rnd}: "
            + "  ".join(
                f"{arm} " + "/".join(f"{per[arm][i][-1]:.1f}" for i in IMPLS) for arm in ARMS
            ),
            flush=True,
        )
    end = preflight.preflight()
    res = dict(
        label=info["label"],
        preflight_start=info,
        preflight_end_loadavg=end.get("loadavg"),
        torch=torch.__version__,
        torch_threads=torch.get_num_threads(),
        n_orientations=len(orients),
        rounds=a.rounds,
        impls=list(IMPLS),
        arms={
            k: dict(pixel_radius=v[0], max_q=v[1], cpp_overlap_rule=v[2]) for k, v in ARMS.items()
        },
        us_per_call={
            arm: {
                i: dict(
                    median=statistics.median(per[arm][i]),
                    min=min(per[arm][i]),
                    max=max(per[arm][i]),
                    rounds=per[arm][i],
                )
                for i in IMPLS
            }
            for arm in ARMS
        },
    )
    (out / "timing.json").write_text(json.dumps(res, indent=1) + "\n")
    for arm in ARMS:
        print(
            f"{arm:22s} "
            + "  ".join(f"{i} {res['us_per_call'][arm][i]['median']:.1f}" for i in IMPLS)
            + " us/call (median)"
        )
    print("label:", info["label"])


if __name__ == "__main__":
    main()
