"""AR5 step 21: the boundary pass of J5's stored map before the pass, with try budget 1 (a screen
try, arm ``J5S``) or 60 (J5's own tries; the gate's re-run, ``--trace`` keeps the gain trace).

  RESEED_SAMPLE=<sample> SCREEN_CACHE=<dir> uv run python \
      scripts/autoresearch/ar5/screen_pass.py --mode full|screen [--limit N] [--no-trace]

Reads ``<cache>/passes_<sample>/before_J5.npy`` and the old ends' scores of ``D4_J5/tries.npz``
(``screen_data``). The boundary voxels and candidates are those of the pass (``Env.boundary``:
voxels in index order, candidates in neighbour order), each try draws from ``order_draws.make_rng(v,
'pass', k)`` (k the number of tries made at the voxel so far: J5's pass used the same, ``--aligned
1``), and a try is ``screen_rules.pass_try``. ``--mode full`` writes ``<cache>/full_<sample>/``
(``tries.npz`` with float64 ends, ``trace.npz``: every strict new low of each try's running best
cost, ``pass.json``); ``--mode screen`` writes ``<cache>/screen_<sample>/`` (``tries.npz``,
``pass.json``). ``--limit N`` takes the first N boundary voxels (a pilot; the folder gets
``_pilot``). The freeze check runs first for every sample.
"""

import argparse
import json
import os
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

import screen_common as SC  # noqa: E402  (first: guards, freeze for the held samples, cache)

for _k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_k, "1")
import numpy as np  # noqa: E402
import screen_data as SD  # noqa: E402
import screen_rules as SR  # noqa: E402
import torch  # noqa: E402

torch.set_num_threads(1)


class Tracer:
    """Records, per try, every strict new low of the running best cost (evaluation number, cost,
    hit ratio, orientation), through a subclass of the reconstructor's cost function."""

    def __init__(self, bfs: Any) -> None:
        self.bfs = bfs
        self.on = False
        self.n = 0
        self.best = np.inf
        self.gains: List[Any] = []
        import icenine.reconstructor as REC
        from icenine.cost_functions import VoxelCostFunction

        outer = self

        class Traced(VoxelCostFunction):  # type: ignore[misc]
            def evaluate(self, orientation: Any, voxel_vertices: Any, phase_index: int = 0) -> Any:
                info = super().evaluate(orientation, voxel_vertices, phase_index)
                if outer.on:
                    outer.n += 1
                    if info.cost < outer.best:
                        outer.best = info.cost
                        hit, _c = outer.bfs._ratios(info)
                        outer.gains.append(
                            (outer.n, float(info.cost), float(hit), np.array(orientation))
                        )
                return info

        self._REC = REC
        self._orig = REC.VoxelCostFunction
        REC.VoxelCostFunction = Traced

    def start(self) -> None:
        self.on, self.n, self.best, self.gains = True, 0, np.inf, []

    def stop(self) -> List[Any]:
        self.on = False
        return self.gains

    def close(self) -> None:
        self._REC.VoxelCostFunction = self._orig


def run_pass(
    env: Any,
    Rb: np.ndarray,
    budget: int,
    trace: bool = False,
    limit: int = 0,
    log: Optional[Any] = None,
) -> Dict[str, Any]:
    """Every try of every boundary voxel of ``Rb`` with the given budget (see the module text)."""
    import order_draws as OD

    isb, cands = env.boundary(Rb)
    voxels = sorted(cands)
    if limit:
        voxels = voxels[:limit]
    tracer = Tracer(env.bfs) if trace else None
    cols: Dict[str, List[Any]] = {k: [] for k in ("voxel", "cand", "cost", "hit", "conf", "evals")}
    ends: List[np.ndarray] = []
    g_n: List[int] = []
    g_cost: List[float] = []
    g_hit: List[float] = []
    g_R: List[np.ndarray] = []
    off = [0]
    t0 = time.time()
    last = t0
    try:
        for i, v in enumerate(voxels):
            for k, c in enumerate(cands[v]):
                rng = OD.make_rng(int(v), "pass", k)
                if tracer is not None:
                    tracer.start()
                t = SR.pass_try(env.hooks, env.bfs, int(v), Rb[c], rng, budget)
                if tracer is not None:
                    for n_, c_, h_, R_ in tracer.stop():
                        g_n.append(n_)
                        g_cost.append(c_)
                        g_hit.append(h_)
                        g_R.append(R_)
                    off.append(len(g_n))
                cols["voxel"].append(int(v))
                cols["cand"].append(int(c))
                for key in ("cost", "hit", "conf", "evals"):
                    cols[key].append(t[key])
                ends.append(t["end"])
            if log is not None and time.time() - last > 120:
                last = time.time()
                log(f"boundary voxels done {i + 1} of {len(voxels)}, {last - t0:.0f}s")
    finally:
        if tracer is not None:
            tracer.close()
    out: Dict[str, Any] = dict(
        voxel=np.array(cols["voxel"], np.int32),
        cand=np.array(cols["cand"], np.int32),
        end=np.array(ends, np.float64).reshape(-1, 3, 3),
        cost=np.array(cols["cost"]),
        hit=np.array(cols["hit"]),
        conf=np.array(cols["conf"]),
        evals=np.array(cols["evals"], np.int32),
        n_boundary_voxels=len(voxels),
        n_boundary_all=int(isb.sum()),
        wall_s=time.time() - t0,
    )
    if tracer is not None:
        out["trace"] = dict(
            off=np.array(off, np.int64),
            g_n=np.array(g_n, np.int16),
            g_cost=np.array(g_cost),
            g_hit=np.array(g_hit),
            g_R=np.array(g_R, np.float64).reshape(-1, 3, 3),
        )
    return out


def save_pass(out: Path, res: Dict[str, Any], info: Dict[str, Any]) -> None:
    out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out / "tries.npz",
        voxel=res["voxel"],
        cand=res["cand"],
        end=res["end"],
        cost=res["cost"],
        hit=res["hit"],
        conf=res["conf"],
        evals=res["evals"],
    )
    if "trace" in res:
        np.savez_compressed(out / "trace.npz", evals=res["evals"], **res["trace"])
    (out / "pass.json").write_text(json.dumps(info, indent=2) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True, choices=("full", "screen"))
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--no-trace", action="store_true")
    a = ap.parse_args()
    SC.SF.check_or_exit()  # the stop rule: nothing runs before frozen.json is committed
    s = SC.SAMPLE
    if not s:
        raise SystemExit("screen_pass.py: set RESEED_SAMPLE")
    name = ("full_" if a.mode == "full" else "screen_") + s + ("_pilot" if a.limit else "")
    out = SC.RCACHE / name
    if (out / "pass.json").exists():
        raise SystemExit(f"refusing: {out / 'pass.json'} exists")
    d = SD.load_inputs(s, truth=False)
    import boundary_common as B  # noqa: E402  (after screen_common: the sample is set)

    env = B.Env("grad")
    assert env.n == len(d["R_before"]), (env.n, len(d["R_before"]))
    budget = SR.FULL_MAX_EVALS if a.mode == "full" else SR.SCREEN_MAX_EVALS
    trace = a.mode == "full" and not a.no_trace
    res = run_pass(env, d["R_before"], budget, trace=trace, limit=a.limit, log=print)
    info = dict(
        sample=s,
        mode=a.mode,
        budget=budget,
        trace=trace,
        limit=a.limit,
        n_boundary_voxels=res["n_boundary_voxels"],
        n_boundary_all=res["n_boundary_all"],
        n_tries=int(len(res["voxel"])),
        evals=int(res["evals"].sum()),
        wall_s=res["wall_s"],
        timing_label="contended (not a timing result)",
        rng="order_draws.make_rng(v, 'pass', k)",
    )
    save_pass(out, res, info)
    print("done", json.dumps(info), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
