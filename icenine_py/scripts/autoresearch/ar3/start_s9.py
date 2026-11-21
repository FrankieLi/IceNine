"""AR3 step 9, part 3: run S9 = S8 plus every remedy that passed Part 2 (cap 90 minutes).

  scripts/autoresearch/ar3/start_bg.sh <big|grad|mid> S9_<sample> start_s9.py
  ... --smoke N          # pilot: first N voxels (a BFS may still flood the whole layer)
  ... --f1 0.95 --f2 0.1 # settings by hand (a pilot); the default reads part2_verdict.json

S8 is ``reseed_s8.py`` (``S5`` with the seed step's refinements stopped at the first with hit
ratio >= 0.95). S9 is that class (``reseed_s8.E1Pick``) with one more instance method replaced,
``bfs._local_fit``; the library is not changed.

F1, the start test: before a neighbour fit or revisit, one evaluation of the voxel's start (the
orientation inherited from the expander; one strict cost evaluation, counted in the voxel's and
the phase's evaluations). If its hit ratio is at least the threshold the start is kept as the
fit's result (the cost and overlap of that evaluation) and the fit is skipped; the BFS then
applies its own relative test to the kept start as it does to any fit's result. Otherwise the fit
runs as before (the evaluation is not reused: it is counted).

F2, one more fit from the end of a rejected fit: after a fit whose result the BFS will reject
(hit ratio / max(flood best before, hit ratio) <= 0.9, the flood's best read from the caller's
frame as ``reseed_logged_run.py`` does), one more local optimization from the fit's end (CMA, the
setting's first step, at most 60 evaluations plus the final one, own generator ``default_rng(0)``
so the run's generator is not consumed). If its hit ratio is higher it replaces the fit's result.
The BFS then judges the result by the same relative test. Its evaluations are counted in the
voxel's and the phase's.

Writes ``<RESEED_CACHE>/S9_<sample>/`` with the files of ``run_seed.run`` and ``s9.json``
(counters of F1 and F2), ``e1_log.json.gz`` (the seed steps, as S8's) and ``f1_log.json.gz``
(one row per start test: voxel, phase, hit ratio of the start, kept).
"""

import argparse
import json
import sys
import time
from typing import Any, Dict, List, Optional

import numpy as np
import reseed_common as RC
import reseed_s8 as S8

import ar3_common as C  # noqa: E402
import run_seed as RS  # noqa: E402

OUTD = RC.CACHE / f"S9_{RC.SAMPLE}"
CFG: Dict[str, Optional[float]] = {"f1": None, "f2": None}
PHASES = ("neighbor", "revisit")


class S9Pick(S8.E1Pick):
    def install(self) -> None:
        super().install()
        b = self.bfs
        rec = self.rec
        o_lf = b._local_fit
        orig_lo = self.orig_lo
        cf = rec._make_local_cost_fn()
        st: Dict[str, Any] = dict(
            f1_tests=0,
            f1_kept=0,
            f1_kept_by_phase={p: 0 for p in PHASES},
            f1_tests_by_phase={p: 0 for p in PHASES},
            f2_runs=0,
            f2_evals=0,
            f2_replaced=0,
            f2_would_pass=0,
        )
        self.s9 = st
        self.f1_log: List[List[float]] = []

        def local_fit(
            voxel: Any, vertices: Any, rng: Any, rec_: Any, phase: str, *a: Any, **k: Any
        ) -> Any:
            from icenine.orientation_search import SearchCandidate

            plain = phase in PHASES and not a and not k
            if CFG["f1"] is not None and plain:
                R0 = np.asarray(voxel.orientation, np.float64)
                info = cf.evaluate(R0, vertices, voxel.phase)
                rec_.n_evals += 1
                b.stats[f"n_evals_{phase}"] += 1
                hit, conf = b._ratios(info)
                keep = bool(hit >= CFG["f1"])
                st["f1_tests"] += 1
                st["f1_tests_by_phase"][phase] += 1
                self.f1_log.append([int(rec_.voxel_idx), PHASES.index(phase), float(hit), keep])
                if keep:
                    st["f1_kept"] += 1
                    st["f1_kept_by_phase"][phase] += 1
                    cand = SearchCandidate(orientation=R0.copy(), cost=float(info.cost))
                    cand.overlap_info = info
                    return cand, hit, conf
            res, hit, conf = o_lf(voxel, vertices, rng, rec_, phase, *a, **k)
            if CFG["f2"] is not None and plain:
                best = float(sys._getframe(1).f_locals.get("best_conf", -1.0))
                ba = max(hit, best)
                if ba > 0 and hit / ba <= 0.9:  # the BFS will reject this fit
                    p = rec.params
                    old = p.local_optimizer
                    p.local_optimizer = "cma"
                    try:
                        r2 = orig_lo(
                            voxel_vertices=vertices,
                            phase_index=voxel.phase,
                            initial_orientation=np.asarray(res.orientation, np.float64),
                            rng=np.random.default_rng(0),
                            cma_sigma0_deg=float(CFG["f2"]),
                            cma_max_evals=60,
                        )
                    finally:
                        p.local_optimizer = old
                    n2 = int(rec.last_local_optimization_evals)
                    rec_.n_evals += n2
                    b.stats[f"n_evals_{phase}"] += n2
                    st["f2_runs"] += 1
                    st["f2_evals"] += n2
                    h2, c2 = b._ratios(r2.overlap_info)
                    if h2 > hit:
                        res, hit, conf = r2, h2, c2
                        st["f2_replaced"] += 1
                    ba2 = max(hit, best)
                    st["f2_would_pass"] += int(ba2 > 0 and hit / ba2 > 0.9)
            return res, hit, conf

        b._local_fit = local_fit


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", type=int, default=0)
    ap.add_argument("--f1", type=float, default=None)
    ap.add_argument("--f2", type=float, default=None)
    ap.add_argument("--pilot-s", type=float, default=0.0, help="pilot: stop after this many s")
    a = ap.parse_args()
    if a.f1 is None and a.f2 is None:
        v = json.loads((RC.OUT.parent / "start" / "part2_verdict.json").read_text())
        if v["F1"]["verdict"]["passes"]:
            CFG["f1"] = float(v["F1"]["verdict"]["threshold"])
        if v.get("F2", {}).get("verdict", {}).get("passes"):
            CFG["f2"] = float(v["F2"]["verdict"]["sigma_deg"])
    else:
        CFG["f1"], CFG["f2"] = a.f1, a.f2
    out = OUTD if not a.smoke else RC.CACHE / f"S9_{RC.SAMPLE}_smoke"
    out.mkdir(parents=True, exist_ok=True)

    class Keep(S9Pick):
        def __init__(self, *aa: Any, **kk: Any) -> None:
            super().__init__(*aa, **kk)
            S8.state["pick"] = self

    RS.SeedPick = Keep
    S8.patch_smoke(a.smoke)
    t0 = time.time()
    S8.watchdog(out, S8.CAP_MIN * 60.0, t0)
    if a.pilot_s:  # a pilot: after this many seconds write the counters so far and stop

        def stop() -> None:
            import os
            import threading

            def loop() -> None:
                time.sleep(a.pilot_s)
                p = S8.state["pick"]
                st = {k: v for k, v in S8.state["bfs"].stats.items() if isinstance(v, (int, float))}
                (out / "pilot.json").write_text(
                    json.dumps(dict(s9=p.s9, stats=st, f1_rows=len(p.f1_log)), default=float)
                )
                os._exit(0)

            threading.Thread(target=loop, daemon=True).start()

        stop()
    RS.run("S5", out, a.smoke or None)
    pick = S8.state["pick"]
    RC.write_gz(pick.e1_log, out / "e1_log.json.gz")
    RC.write_gz(pick.f1_log, out / "f1_log.json.gz")
    (out / "s9.json").write_text(
        json.dumps(
            dict(
                run="S9 = S8 + F1/F2",
                sample=RC.SAMPLE,
                settings=CFG,
                tau=S8.TAU,
                seed_steps=len(pick.e1_log),
                stopped_early=sum(1 for e in pick.e1_log if e["stop"] >= 0),
                counters=pick.s9,
                wall_s=time.time() - t0,
                cap_min=S8.CAP_MIN,
            ),
            indent=1,
        )
        + "\n"
    )
    (out / "s9.done").write_text("done\n")


if __name__ == "__main__":
    main()
