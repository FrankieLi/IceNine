"""AR3 step 8, part 3: run S8 = S5 with remedy E1 (stop the refinements at the first that passes).

  RESEED_SAMPLE=big  uv run python scripts/autoresearch/ar3/reseed_s8.py          # cap 90 minutes
  RESEED_SAMPLE=grad uv run python scripts/autoresearch/ar3/reseed_s8.py
  ... --smoke N     # pilot: first N voxels (a BFS may still flood the whole layer)

S8 is step 6's run S5 (``run_seed.run("S5", ...)``, unchanged: indexed seeds, the step 4 fit
wrapper, rule-free BFS with reseed) with one change in a script, none in the library: the seed
step's refinements stop at the first one (in the order of the strict cost of the 20 costed
candidates) whose end has hit ratio >= ``TAU`` (set on the pick set in ``reseed_e1.py``,
``remedies_pick.json``), and that end is the seed step's result; if none passes, all 5 are
refined and the lowest final cost is kept, as S5 does. The seed step's costing, indexing and
refinement code are the library's own (``ray_indexing.index_voxel``, ``refine_from_candidates``);
only the loop of ``ray_indexing.seed_by_indexing`` is repeated here with the early stop.

The hit ratio of a refinement's end is one overlap evaluation (``evaluate_overlap``, the call the
library makes after a seed step); like the library's, it is not counted in the evaluation totals
(at most 5 per seed step against several thousand counted).

Writes ``reseed/cache/S8_<sample>/`` with the files of ``run_seed.run`` and ``e1_log.json.gz`` (one
row per seed step: voxel, reseed, evaluations, refinement evaluations and hit ratios, the index
of the refinement that stopped it).
"""

import json
import os
import sys
import threading
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import reseed_common as RC

import ar3_common as C  # noqa: E402
import run_seed as RS  # noqa: E402
from seed_pick import SeedPick  # noqa: E402

TAU = float(json.loads((RC.OUT / "remedies_pick.json").read_text())["chosen_tau"])
CAP_MIN = 90
OUTD = RC.CACHE / f"S8_{RC.SAMPLE}"


class E1Pick(SeedPick):
    """``SeedPick`` (step 6) with its seed-step call replaced by the E1 loop; no tries."""

    def __init__(self, *a: Any, **k: Any) -> None:
        super().__init__(*a, **k)
        assert not self.use_tries
        self.e1_log: List[Dict[str, Any]] = []
        self.tau = TAU
        self.orig_index_seed = self.e1_step  # SeedPick calls this with n_cost=20, n_refine=5

    def e1_step(
        self, vertices: Any, phase_index: int, n_cost: int = 5, n_refine: int = 2, tau: float = -1
    ) -> Any:
        from icenine import ray_indexing as RI
        from icenine.orientation_search import SearchCandidate

        tau = self.tau if tau < 0 else tau
        rec = self.rec
        t0 = time.perf_counter()
        sc, crystal, params, strict_cf = rec._index_context(phase_index)

        def strict(R: np.ndarray) -> float:
            return float(strict_cf.evaluate(np.asarray(R, np.float64), vertices, phase_index).cost)

        res, att = RI.index_voxel(sc, vertices, crystal, params, RI.SEED_TOP_K)
        index_s = time.perf_counter() - t0
        out = RI.SeedResult(
            n_cand=len(res.R),
            n_att=int(len(att)),
            index_s=index_s,
            total_s=index_s,
            candidates=res.R,
        )
        entry: Dict[str, Any] = dict(
            voxel=self.cur,
            reseed=self.reseed,
            n_cand=out.n_cand,
            stop=-1,
            fin_evals=[],
            fin_hit=[],
            costed=0,
        )
        if len(res.R) == 0:
            self.e1_log.append(entry)
            return out
        R = res.R
        k = min(n_cost, len(R))
        costs = np.array([strict(R[j]) for j in range(k)], dtype=np.float64)
        out.strict_costs = [float(c) for c in costs]
        out.evals = k
        entry["costed"] = k
        order = [int(j) for j in np.argsort(costs, kind="stable")[:n_refine]]
        stop: Optional[int] = None
        for pos, j in enumerate(order):
            r = rec.refine_from_candidates(
                [SearchCandidate(orientation=R[j].copy(), cost=float(costs[j]))],
                vertices,
                phase_index,
                rng=np.random.default_rng(0),
            )
            out.refined.append(j)
            out.finisher_costs.append(float(r.cost))
            out.finisher_evals.append(int(rec.last_eval_counts[1]))
            out.finisher_R.append(np.asarray(r.orientation, dtype=np.float64))
            info = rec.evaluate_overlap(
                np.asarray(r.orientation, np.float64), vertices, phase_index
            )
            hit = float(self.bfs._ratios(info)[0])
            entry["fin_evals"].append(out.finisher_evals[-1])
            entry["fin_hit"].append(hit)
            if hit >= tau:
                stop = pos
                break
        out.evals += int(sum(out.finisher_evals))
        best = int(np.argmin(out.finisher_costs)) if stop is None else stop
        out.taken = out.refined[best]
        out.orientation = out.finisher_R[best]
        out.cost = out.finisher_costs[best]
        out.total_s = time.perf_counter() - t0
        entry["stop"] = -1 if stop is None else stop
        entry["evals"] = int(out.evals)
        self.e1_log.append(entry)
        return out


state: Dict[str, Any] = {}


def patch_pick() -> None:
    class Keep(E1Pick):
        def __init__(self, *a: Any, **k: Any) -> None:
            super().__init__(*a, **k)
            state["pick"] = self

    RS.SeedPick = Keep


def patch_smoke(smoke: int) -> None:
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


def watchdog(out: Path, cap_s: float, t0: float) -> None:
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
                dict(cap_s=cap_s, elapsed_s=time.time() - t0, voxels_with_record=n_rec, stats=st),
                default=float,
            )
            + "\n"
        )
        print(f"S8: CAP of {cap_s} s reached; stopped", flush=True)
        os._exit(3)

    threading.Thread(target=loop, daemon=True).start()


def main() -> None:
    smoke = 0
    if "--smoke" in sys.argv:
        smoke = int(sys.argv[sys.argv.index("--smoke") + 1])
    out = OUTD if not smoke else RC.CACHE / f"S8_{RC.SAMPLE}_smoke"
    out.mkdir(parents=True, exist_ok=True)
    patch_pick()
    patch_smoke(smoke)
    t0 = time.time()
    watchdog(out, CAP_MIN * 60.0, t0)
    RS.run("S5", out, smoke or None)
    pick = state["pick"]
    RC.write_gz(pick.e1_log, out / "e1_log.json.gz")
    (out / "s8.json").write_text(
        json.dumps(
            dict(
                run="S8 = S5 + E1",
                sample=RC.SAMPLE,
                tau=TAU,
                seed_steps=len(pick.e1_log),
                stopped_early=sum(1 for e in pick.e1_log if e["stop"] >= 0),
                wall_s=time.time() - t0,
                cap_min=CAP_MIN,
            ),
            indent=1,
        )
        + "\n"
    )
    (out / "s8.done").write_text("done\n")


if __name__ == "__main__":
    main()
