"""AR5 step 2, Part 2: whole-layer runs of S8 and S10 (the script layer; the library is unchanged).

  scripts/autoresearch/ar5/fit_bg.sh <big|grad|mid> <name>_<sample> fit_run.py --name S10 \
      --seed-arm C60 --no-mc-fallback [--smoke N]

``S8`` (``--seed-arm M`` and no ``--no-mc-fallback``) is step 8's ``reseed_s8.py``: run ``S5``
(indexed seeds, the step 4 fit wrapper, rule-free BFS with reseed) with the seed step's refinements
stopped at the first end with hit ratio >= 0.95. This script reaches the same code through the same
classes (``reseed_s8.E1Pick``, ``run_seed.run("S5", ...)``, ``run_wrapped.Wrap``); the base check
(Part 2) is that it reproduces the stored S8 exactly when both changes are off.

The two changes of ``S10`` (plan, AR5 step 2), each off by default:

  (i)  ``--seed-arm C60|C200|C0``: the seed step refines each indexed candidate with that arm
       (``fit_common.refine_arm``) in place of the stock MC refinement (arm ``M``);
  (ii) ``--no-mc-fallback``: a neighbour fit or revisit returns what the CMA fit returned. The
       wrapper's stock MC fit is not run when the CMA end cost is not below the start cost. The
       number of such fits is counted as ``cma_not_fell`` (the fallback would have run).

The boundary pass is made after by ``fit_pass.py``. Writes ``<cache>/<name>_<sample>/`` with the
files of ``run_seed.run`` and ``e1_log.json.gz`` (one row per seed step), ``fit.json``,
``fit.done``; a run that reaches the 90 minute cap writes ``capped.json`` and stops (reported as
unfinished).
"""

import argparse
import json
import sys
import time
from typing import Any, Dict, List

import fit_common as FC
import fit_rules as FR
import numpy as np
import reseed_common as RC
import reseed_s8 as S8

import run_seed as RS  # noqa: E402
import run_wrapped as RW  # noqa: E402

CFG: Dict[str, Any] = {"arm": "M", "no_fb": False}


class FitPick(S8.E1Pick):
    """``E1Pick`` (step 8) with the refinement of a candidate chosen by ``CFG['arm']``."""

    def e1_step(
        self, vertices: Any, phase_index: int, n_cost: int = 5, n_refine: int = 2, tau: float = -1
    ) -> Any:
        from icenine import ray_indexing as RI

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
        stop = None
        for pos, j in enumerate(order):
            Rend, c_end, ev = FC.refine_arm(
                CFG["arm"], rec, self.orig_lo, R[j], float(costs[j]), vertices, phase_index
            )
            out.refined.append(j)
            out.finisher_costs.append(float(c_end))
            out.finisher_evals.append(int(ev))
            out.finisher_R.append(Rend)
            info = rec.evaluate_overlap(Rend, vertices, phase_index)
            hit = float(self.bfs._ratios(info)[0])
            entry["fin_evals"].append(out.finisher_evals[-1])
            entry["fin_hit"].append(hit)
            if hit >= tau:
                stop = pos
                break
        out.evals += int(sum(out.finisher_evals))
        best, _stopped = FR.early_stop(entry["fin_hit"], out.finisher_costs, tau)
        assert (stop is None) == (not _stopped) and (stop is None or stop == best)
        out.taken = out.refined[best]
        out.orientation = out.finisher_R[best]
        out.cost = out.finisher_costs[best]
        out.total_s = time.perf_counter() - t0
        entry["stop"] = -1 if stop is None else stop
        entry["evals"] = int(out.evals)
        self.e1_log.append(entry)
        return out


class FitWrap(RW.Wrap):
    """Step 4's wrapper; with ``no_fb`` the stock MC fit never runs (the CMA end is returned)."""

    def __init__(self, bfs: Any) -> None:
        self.cma_not_fell = 0
        super().__init__(bfs)
        if CFG["no_fb"]:
            self.rec.local_optimization = self.wrapped_no_fallback

    def wrapped_no_fallback(
        self,
        voxel_vertices: Any,
        phase_index: int,
        initial_orientation: np.ndarray,
        rng: Any = None,
        **kw: Any,
    ) -> Any:
        assert not kw, "the wrapper expects the classic-mode call (no cma keywords)"
        self.calls[sys._getframe(1).f_code.co_name] += 1
        R0 = np.asarray(initial_orientation, np.float64)
        r1, n1, c0 = self.cma_fit(voxel_vertices, phase_index, R0, rng)
        self.evals_picked += n1
        if r1.cost < c0:
            self.n_picked_fell += 1
        else:
            self.cma_not_fell += 1  # the stock MC fit would have run here
        self.rec.last_local_optimization_evals = n1
        return r1

    def counts(self) -> Dict[str, Any]:
        d = super().counts()
        d["cma_not_fell"] = self.cma_not_fell
        d["mc_fallback_off"] = bool(CFG["no_fb"])
        return d


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True, help="S8x (base check), S10, S10a, S10b")
    ap.add_argument("--seed-arm", default="M", choices=list(FC.ARMS))
    ap.add_argument("--no-mc-fallback", action="store_true")
    ap.add_argument("--smoke", type=int, default=0)
    ap.add_argument("--cap-min", type=float, default=float(S8.CAP_MIN))
    a = ap.parse_args()
    CFG["arm"], CFG["no_fb"] = a.seed_arm, bool(a.no_mc_fallback)
    out = FC.CACHE / f"{a.name}_{FC.SAMPLE}{'_smoke' if a.smoke else ''}"
    out.mkdir(parents=True, exist_ok=True)

    class Keep(FitPick):
        def __init__(self, *aa: Any, **kk: Any) -> None:
            super().__init__(*aa, **kk)
            S8.state["pick"] = self

    RS.SeedPick = Keep
    RW.Wrap = FitWrap
    S8.patch_smoke(a.smoke)
    t0 = time.time()
    S8.watchdog(out, a.cap_min * 60.0, t0)
    RS.run("S5", out, a.smoke or None)
    pick: FitPick = S8.state["pick"]
    RC.write_gz(pick.e1_log, out / "e1_log.json.gz")
    log: List[Dict[str, Any]] = pick.e1_log
    (out / "fit.json").write_text(
        json.dumps(
            dict(
                run=a.name,
                sample=FC.SAMPLE,
                seed_arm=a.seed_arm,
                no_mc_fallback=bool(a.no_mc_fallback),
                tau=S8.TAU,
                seed_steps=len(log),
                stopped_early=sum(1 for e in log if e["stop"] >= 0),
                wall_s=time.time() - t0,
                wall_label="contended (other jobs on the machine); not a timing result",
                cap_min=a.cap_min,
            ),
            indent=1,
        )
        + "\n"
    )
    (out / "fit.done").write_text("done\n")


if __name__ == "__main__":
    main()
