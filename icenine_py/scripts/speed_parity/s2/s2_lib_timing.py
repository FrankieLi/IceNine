"""S2 Step 5: us per call of the library fused evaluate (key FusedCost) against the torch
``evaluate``, single orientation and batched, and the GPU_MIN_BATCH check.

Plan ``docs/todo_s2_fused_c_eval.md``, "Steps 4-6", item 5 (Step 0c repeated on the library path).

``single``: S1's per-call set (``prof_eval``: ThreeVoxels voxel 0, ``Example2.Simulation.config``,
400 orientations within 3 degrees of the truth, seed 0) and Step 0c's four arms (the Q-max 8
triangle cost and the radius-3 Q-max 5 centroid screen, each with the default rule and with
``CppOverlapRule``). Implementations: ``torch`` (``evaluate``, key off), ``fused`` (``evaluate``,
key on: the library path including the ``OverlapInfo``), ``c_only`` (the bare
``_rasterize.fused_evaluate`` call on the cost function's context: the Python share of ``fused``
is the difference). One untimed warm-up pass fills the frames and builds the contexts.

``batched``: ``BatchedCostEvaluator.evaluate_many`` on the local (triangle) cost of ThreeVoxels
voxel 0 (Q-max 8, and the search config's local cost) and of the first E0 clean voxel
(``ReconstructQ8`` physics; the reflection count K is in the output), 512
orientations within 3 degrees of the truth, scored in chunks of B; us per orientation for the
variants ``compact`` (torch, key off), ``fused`` (key on, ``fused_evaluate_many``) and, when a GPU
is available, ``fixed`` on it (the variant ``evaluate_many`` picks at B >= ``GPU_MIN_BATCH``). With
the key on and a GPU device, B >= ``GPU_MIN_BATCH`` (64) still goes to the GPU: the ``fixed``
against ``fused`` columns measure whether it should (reported only; the constant is not changed,
D-20261028-I1). ``evaluate`` (one orientation per call, key off and on) is the B = 1 reference.

Single worker (torch threads 1), preflight-gated (``--allow-busy`` labels the run "contended");
``--rounds`` interleaved rounds, the implementation order rotating; the median over rounds.

  uv run python scripts/speed_parity/s2/s2_lib_timing.py single --out benchmarks/s2_results
  uv run python scripts/speed_parity/s2/s2_lib_timing.py batched --out benchmarks/s2_results
"""

import argparse
import contextlib
import io
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
sys.path.insert(0, str(ICE))
sys.path.insert(0, str(ICE / "scripts" / "common"))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))

import preflight  # noqa: E402
from s2_timing import gate  # noqa: E402

SINGLE_ARMS = {
    "triangle_q8_default": (0, 8.0, False),
    "triangle_q8_cpp": (0, 8.0, True),
    "screen_r3_q5_default": (3, 5.0, False),
    "screen_r3_q5_cpp": (3, 5.0, True),
}
SINGLE_IMPLS = ("torch", "fused", "c_only")
BATCHES = (1, 8, 16, 32, 64, 128, 256, 512)


def summarize(per: Dict[str, Dict[str, List[float]]]) -> Dict[str, Any]:
    return {
        arm: {
            i: dict(median=statistics.median(v), min=min(v), max=max(v), rounds=v)
            for i, v in d.items()
        }
        for arm, d in per.items()
    }


def cmd_single(a: argparse.Namespace, info: Dict[str, Any]) -> Dict[str, Any]:
    import prof_eval as P  # chdirs into the example, builds the physics, torch threads 1

    from icenine import cost_functions as CF

    torch.set_num_threads(1)
    orients = P.orients[: a.n]
    verts, phase = P.verts, P.vox.phase
    raw = CF._c_fused_evaluate
    calls: Dict[str, Dict[str, Callable[[], None]]] = {}
    for arm, args in SINGLE_ARMS.items():
        cf_t = P.make(*args)
        cf_f = P.make(*args)
        cf_f.fused_cost = True
        cf_r = P.make(*args)
        cf_r.fused_cost = True

        def run_torch(cf: Any = cf_t) -> None:
            for o in orients:
                cf.evaluate(o, verts, phase)

        def run_fused(cf: Any = cf_f) -> None:
            for o in orients:
                cf.evaluate(o, verts, phase)

        for o in orients:  # builds the bare call's context and fills its frames
            cf_r.evaluate(o, verts, phase)
        cap = cf_r._fused_context(phase).capsule
        vn = cf_r._fused_vertices(verts)

        def run_raw(c: Any = cap, v: np.ndarray = vn) -> None:
            for o in orients:
                raw(c, o, v)

        calls[arm] = dict(torch=run_torch, fused=run_fused, c_only=run_raw)
        for fn in calls[arm].values():
            fn()  # warm-up
        if cf_f.eval_count != len(orients) or cf_f._fused_ctx.get(phase) is None:
            raise SystemExit(f"{arm}: the key-on cost function did not run the fused path")
    per: Dict[str, Dict[str, List[float]]] = {arm: {i: [] for i in SINGLE_IMPLS} for arm in calls}
    for rnd in range(a.rounds):
        k = rnd % len(SINGLE_IMPLS)
        order = list(SINGLE_IMPLS[k:] + SINGLE_IMPLS[:k])
        for arm in calls:
            for impl in order:
                t0 = time.perf_counter()
                calls[arm][impl]()
                per[arm][impl].append(1e6 * (time.perf_counter() - t0) / len(orients))
        print(
            f"round {rnd}: "
            + "  ".join(
                f"{arm} " + "/".join(f"{per[arm][i][-1]:.1f}" for i in SINGLE_IMPLS)
                for arm in calls
            ),
            flush=True,
        )
    return dict(
        set="ThreeVoxels voxel 0 (prof_eval)",
        n_orientations=len(orients),
        arms={
            k: dict(pixel_radius=v[0], max_q=v[1], cpp_overlap_rule=v[2])
            for k, v in SINGLE_ARMS.items()
        },
        impls=list(SINGLE_IMPLS),
        us_per_call=summarize(per),
    )


def _near(truth: np.ndarray, n: int, seed: int) -> np.ndarray:
    import s1_parity as S1

    rng = np.random.default_rng(seed)
    return np.stack(
        [S1._rot(rng.normal(size=3), np.deg2rad(rng.uniform(0.0, 3.0))) @ truth for _ in range(n)]
    )


def _cases() -> List[Dict[str, Any]]:
    """(name, local cost function factory, vertices, phase, truth) of the batched sets: ThreeVoxels
    voxel 0 at Q-max 8 (``three0``, Step 0c's triangle arm) and with the local cost of its search
    config (``three0_search``, ``Example2.Simulation.config``, the reconstructor's local cost), and
    the first E0 clean voxel (``e0_0``)."""
    import prof_eval as P
    from s3_step0 import E0Set, ThreeVoxels

    with contextlib.redirect_stdout(io.StringIO()):
        rec3 = ThreeVoxels().reconstructor("Example2.Simulation.config")
    out = [
        dict(
            name="three0",
            make=lambda: P.make(0, 8.0, False),
            verts=P.verts,
            phase=P.vox.phase,
            truth=P.truth,
        ),
        dict(
            name="three0_search",
            make=rec3._make_local_cost_fn,
            verts=P.verts,
            phase=P.vox.phase,
            truth=P.truth,
        ),
    ]
    with contextlib.redirect_stdout(io.StringIO()):
        SE = E0Set()
        verts, truth, phase, rec = SE.load(0)
    out.append(
        dict(
            name="e0_0", make=rec._make_local_cost_fn, verts=verts, phase=phase, truth=truth
        )  # fmt: skip
    )
    return out


def cmd_batched(a: argparse.Namespace, info: Dict[str, Any]) -> Dict[str, Any]:
    from icenine import batched_cost as BC

    torch.set_num_threads(1)
    gpu = "mps" if BC.device_available("mps") else ("cuda" if BC.device_available("cuda") else None)
    res: Dict[str, Any] = dict(gpu=gpu, gpu_min_batch=BC.GPU_MIN_BATCH, batches=list(BATCHES))
    for case in _cases():
        R = _near(case["truth"], a.n_batch, 7)
        verts, phase = case["verts"], case["phase"]
        cf_off = case["make"]()
        cf_on = case["make"]()
        cf_on.fused_cost = True
        K = int(cf_off._phase_recip_vecs[phase][0].shape[0])
        ev_off = BC.BatchedCostEvaluator(cf_off)
        ev_on = BC.BatchedCostEvaluator(cf_on)
        ev_gpu = BC.BatchedCostEvaluator(case["make"](), device=gpu) if gpu else None
        impls: Dict[str, Callable[[int], None]] = {}

        def chunks(fn: Callable[[np.ndarray], Any], B: int) -> None:
            for s in range(0, len(R), B):
                fn(R[s : s + B])

        impls["compact"] = lambda B: chunks(
            lambda x: ev_off.evaluate_many(x, verts, phase, variant="compact"), B
        )
        impls["fused"] = lambda B: chunks(
            lambda x: ev_on.evaluate_many(x, verts, phase, variant="fused"), B
        )
        if ev_gpu is not None:
            impls["fixed_" + gpu] = lambda B: chunks(
                lambda x: ev_gpu.evaluate_many(x, verts, phase, variant="fixed"), B
            )

        def single(cf: Any) -> Callable[[int], None]:
            def run(B: int) -> None:
                for r in R:
                    cf.evaluate(r, verts, phase)

            return run

        cf_e_off, cf_e_on = case["make"](), case["make"]()
        cf_e_on.fused_cost = True
        names = list(impls)
        for B in BATCHES:  # warm-up: frames, contexts, GPU kernels
            for nm in names:
                impls[nm](B)
        single(cf_e_off)(1)
        single(cf_e_on)(1)
        if ev_on.counts["orient_fused"] == 0 or cf_e_on._fused_ctx.get(phase) is None:
            raise SystemExit(f"{case['name']}: the key-on path did not run the fused entry")
        per: Dict[str, Dict[str, List[float]]] = {
            f"B={B}": {nm: [] for nm in names} for B in BATCHES
        }
        per["evaluate"] = {"torch": [], "fused": []}
        for rnd in range(a.rounds):
            k = rnd % len(names)
            order = names[k:] + names[:k]
            for B in BATCHES:
                for nm in order:
                    t0 = time.perf_counter()
                    impls[nm](B)
                    per[f"B={B}"][nm].append(1e6 * (time.perf_counter() - t0) / len(R))
            for nm, cf in (("torch", cf_e_off), ("fused", cf_e_on)):
                t0 = time.perf_counter()
                single(cf)(1)
                per["evaluate"][nm].append(1e6 * (time.perf_counter() - t0) / len(R))
            print(
                f"{case['name']} round {rnd}: "
                + "  ".join(
                    f"B={B} " + "/".join(f"{per[f'B={B}'][nm][-1]:.1f}" for nm in names)
                    for B in BATCHES
                )
                + "  evaluate "
                + "/".join(f"{per['evaluate'][nm][-1]:.1f}" for nm in ("torch", "fused")),
                flush=True,
            )
        res[case["name"]] = dict(
            K=K, n_orientations=len(R), impls=names, us_per_orientation=summarize(per)
        )
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("cmd", choices=["single", "batched"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--rounds", type=int, default=5)
    ap.add_argument("--n", type=int, default=400, help="single: orientations per pass")
    ap.add_argument("--n-batch", type=int, default=512, help="batched: orientations per pass")
    ap.add_argument("--allow-busy", action="store_true")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    out = Path(a.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    info = gate(a.allow_busy)
    (out / f"preflight_{a.cmd}{a.tag}.json").write_text(json.dumps(info, indent=1) + "\n")
    import s1_parity as S1

    meta = dict(
        head=S1.git_head(),
        dirty=S1.git_dirty(),
        label=info["label"],
        torch=torch.__version__,
        torch_threads=1,
        rounds=a.rounds,
    )
    body = cmd_single(a, info) if a.cmd == "single" else cmd_batched(a, info)
    end = preflight.preflight()
    res = dict(meta=meta, preflight_start=info, preflight_end_loadavg=end.get("loadavg"), **body)
    (out / f"timing_{a.cmd}{a.tag}.json").write_text(
        json.dumps(res, indent=1, default=float) + "\n"
    )
    print("label:", info["label"])


if __name__ == "__main__":
    main()
