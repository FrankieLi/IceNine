"""AR5 step 22 (task ar5v): the seed step with rule ``DEF`` (script layer; library unchanged).

``DeferPick`` is step 16's ``High500Pick`` (arm ``J5``: interleaved order, 20 costed, 5 refined,
rule ``PH(500)`` in the MC refinement) with its seed step run by ``defer_rules.run_seed_step``:

* ``DEFER['on']`` True: arm ``J5D`` (rule ``DEF``). At the first ``PH(500)`` check of an MC
  refinement that goes on, and only at a candidate that is not the last, the refinement is
  PAUSED: its state (evaluation count, running best, its trace, the check's hit ratio) is kept and
  the next candidate is tried. If no later candidate reaches 0.95, the paused refinements are
  resumed in candidate order, each exactly where it stopped, under ``J5``'s rules.
* ``DEFER['on']`` False: the same state machine with nothing paused (the step is ``J5``'s; this is
  gate g0).

**A resume is an exact re-run of the paused part (phase A; the library's MC refinement has no
pause).** The refinement has its own generator (``default_rng(0)`` in the stock MC local
optimization, not keyed by a call count), so re-running it from the candidate reproduces the first
run evaluation by evaluation. The re-run's first part (the evaluations up to the pause and the
check) is CHECKED identical to the first run (check evaluation number, hit ratio, best cost, the
trace of gains) and is NOT counted; it is listed apart (``rerun_evals``). The evaluations a
refinement is charged are those of ``J5``'s unpaused refinement: the part before the pause (the
check included) when paused, the whole refinement's ``J5`` count when resumed. A resume is not a
new call of anything keyed: the aligned draws (voxel, kind of call, count) are not touched by the
seed step's refinements.

Draws no random number; the stall check's hit ratio is the voxel's OWN hit ratio.
"""

from typing import Any, Dict, List, Optional

import defer_common as DC  # noqa: E402,F401  (first: guards, freeze, cache, held6)

import high500_pick as HP  # noqa: E402
import defer_rules as DR  # noqa: E402
import fit_rules as FR  # noqa: E402
import mix_common as MC  # noqa: E402
import numpy as np  # noqa: E402
import order_draws as OD  # noqa: E402
import order_run as ORUN  # noqa: E402
import stallhit_rules as SR  # noqa: E402

import time  # noqa: E402

_BASE = HP.High500Pick  # step 16's class (its ``_fit`` keeps the one tracer on this class)
DEFER: Dict[str, Any] = dict(on=True)
STATS: Dict[str, int] = dict(
    steps=0,
    steps_with_pause=0,
    pauses=0,
    resumes=0,
    found_tried_after_pause=0,
    found_in_resume=0,
    none_all_resumed=0,
    rerun_evals=0,
    rerun_checked=0,
)
LOG: List[Dict[str, Any]] = []  # one row per seed step (gz-written by defer_run)


class DeferPause(Exception):
    """Raised by the stall check's hit function when a refinement is paused."""


class DeferRerunMismatch(AssertionError):
    """The re-run of a paused part is not identical to the first run (stop: the plan's rule)."""


class DeferPick(HP.High500Pick):
    """``High500Pick`` whose seed step runs ``defer_rules.run_seed_step`` (see the module text)."""

    def _tracer(self) -> Any:
        if _BASE.ptracer is None:
            _BASE.ptracer = HP.High500Tracer()
            ORUN.OrderPick.tracer = _BASE.ptracer
            _BASE.tracer = _BASE.ptracer
        return _BASE.ptracer

    def _fit_m(
        self,
        rec: Any,
        j: int,
        Rj: np.ndarray,
        cj: float,
        vertices: Any,
        phase: int,
        allow_pause: bool,
        prior: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        """One MC refinement of the seed step under ``PH(500)``; paused at its first check that
        goes on when ``allow_pause``; ``prior`` (the state of a paused first run) makes it a
        resume: the first part is re-run, checked and not counted."""
        tr = self._tracer()
        tr.begin(int(HP.H5["k"]), None, tuple(HP.H5["aux"]))
        state: Dict[str, Any] = {}
        seen = {"first": False}

        def hit_fn(R: np.ndarray) -> float:
            info = rec.evaluate_overlap(np.asarray(R, np.float64), vertices, phase)
            h = float(self.bfs._ratios(info)[0])
            if prior is not None and not seen["first"]:
                seen["first"] = True
                same = (
                    int(tr.n) == int(prior["n"])
                    and h == prior["hit"]
                    and float(tr.best) == prior["best"]
                    and [list(x) for x in tr.trace] == prior["trace"]
                )
                if not same:
                    raise DeferRerunMismatch(
                        f"re-run differs at the pause: n {tr.n} vs {prior['n']}, hit {h} vs "
                        f"{prior['hit']}, best {tr.best} vs {prior['best']}"
                    )
                STATS["rerun_checked"] += 1
                STATS["rerun_evals"] += int(prior["n"]) + 1
            if allow_pause and h < DR.TAU:
                state.update(
                    n=int(tr.n),
                    hit=h,
                    best=float(tr.best),
                    best_R=np.array(tr.best_R, dtype=np.float64, copy=True),
                    trace=[list(x) for x in tr.trace],
                )
                raise DeferPause()
            return h

        tr.hit_fn = hit_fn
        cut = False
        try:
            Rend, c, ev = MC._orig_refine_arm("M", rec, self.orig_lo, Rj, cj, vertices, phase)
            ev = int(ev) + tr.continued
        except DeferPause:
            tr.on = False
            STATS["pauses"] += 1
            if prior is None:
                ORUN.COUNTS["mc_fits"] += 1
            return dict(
                j=int(j),
                kind="M",
                hit=float(state["hit"]),
                cost=float(state["best"]),
                evals=int(state["n"]) + 1,
                R=np.asarray(state["best_R"]),
                paused=True,
                state=state,
            )
        except SR.StallHitStop:
            cut = True
            assert tr.best_R is not None
            Rend = tr.best_R
            final = rec._make_local_cost_fn().evaluate(
                orientation=Rend, voxel_vertices=vertices, phase_index=phase
            )
            c = float(final.cost)
            ev = tr.n + 1 + tr.continued
        tr.on = False
        info = rec.evaluate_overlap(np.asarray(Rend, np.float64), vertices, phase)
        hit = float(self.bfs._ratios(info)[0])
        if prior is None:
            ORUN.COUNTS["mc_fits"] += 1
        HP.STATS["mc_refinements"] += 1
        HP.STATS["checks_total"] += len(tr.checks)
        HP.STATS["checks_went_on"] += tr.continued
        if cut:
            HP.STATS["stops_high" if tr.stop_kind == "high" else "stops_low"] += 1
            HP.STATS["evals_in_stopped"] += int(ev)
        ORUN.TRACES.append(
            dict(
                voxel=self.cur,
                step=ORUN.COUNTS["seed_steps"],
                j=int(j),
                evals=int(ev),
                n_calls=int(ev),
                trace=tr.trace,
                cut=cut,
                stop_kind=tr.stop_kind,
                checks=[list(x) for x in tr.checks],
                continued=int(tr.continued),
                events={str(a): [list(x) for x in v] for a, v in tr.events.items()},
                resumed=prior is not None,
            )
        )
        return dict(j=int(j), kind="M", hit=hit, cost=float(c), evals=int(ev), R=np.asarray(Rend))

    def e1_step(
        self, vertices: Any, phase_index: int, n_cost: int = 5, n_refine: int = 2, tau: float = -1
    ) -> Any:
        from icenine import ray_indexing as RI

        assert ORUN.ORD["order"] == "A", "rule DEF is defined for the interleaved order"
        tau = self.tau if tau < 0 else tau
        rec = self.rec
        t0 = time.perf_counter()
        ORUN.COUNTS["seed_steps"] += 1
        STATS["steps"] += 1
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
            fits=[],
            order=ORUN.ORD["order"],
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

        def run_c(p: int) -> Dict[str, Any]:
            j = order[p]
            f = self._fit("C", rec, j, R[j], float(costs[j]), vertices, phase_index)
            f["pos"] = p
            return f

        def run_m(p: int, allow: bool) -> Dict[str, Any]:
            j = order[p]
            f = self._fit_m(rec, j, R[j], float(costs[j]), vertices, phase_index, allow)
            f["pos"] = p
            return f

        rer = {"n": 0}

        def resume_m(p: int, prior: Dict[str, Any]) -> Dict[str, Any]:
            j = order[p]
            rer["n"] += int(prior["state"]["n"]) + 1
            f = self._fit_m(
                rec, j, R[j], float(costs[j]), vertices, phase_index, False, prior["state"]
            )
            f["pos"] = p
            return f

        st = DR.run_seed_step(len(order), run_c, run_m, resume_m, tau, defer=bool(DEFER["on"]))
        stop = st["stop"]
        fits: List[Dict[str, Any]] = []
        for kp in st["kept"]:
            p = kp["pos"]
            fits.append(kp["c"])
            if kp["m"] is not None:
                fits.append(kp["m"])
            e = kp["e"]
            out.refined.append(order[p])
            out.finisher_costs.append(float(e["cost"]))
            out.finisher_evals.append(int(kp["evals"]))
            out.finisher_R.append(np.asarray(e["R"], dtype=np.float64))
            entry["fin_evals"].append(int(kp["evals"]))
            entry["fin_hit"].append(float(e["hit"]))
            if kp["paused"]:  # never resumed: its row in the traces is the paused part
                s = kp["m"]["state"]
                ORUN.TRACES.append(
                    dict(
                        voxel=self.cur,
                        step=ORUN.COUNTS["seed_steps"],
                        j=int(kp["m"]["j"]),
                        evals=int(kp["m"]["evals"]),
                        n_calls=int(kp["m"]["evals"]),
                        trace=s["trace"],
                        cut=False,
                        stop_kind="pause",
                        checks=[[s["n"], s["hit"], False]],
                        continued=1,
                        events={},
                        resumed=False,
                    )
                )
        out.evals += int(sum(out.finisher_evals))
        best, stopped = FR.early_stop(entry["fin_hit"], out.finisher_costs, tau)
        assert (stop is None) == (not stopped) and (stop is None or stop == best)
        out.taken = out.refined[best]
        out.orientation = out.finisher_R[best]
        out.cost = out.finisher_costs[best]
        out.total_s = time.perf_counter() - t0
        entry["stop"] = -1 if stop is None else int(stop)
        entry["evals"] = int(out.evals)
        entry["fits"] = [
            [
                f["pos"],
                f["j"],
                f["kind"],
                f["hit"],
                f["cost"],
                f["evals"],
                np.round(f["R"], 12).ravel().tolist(),
            ]
            for f in fits
        ]
        n_p, n_r = len(st["paused"]), len(st["resumed"])
        rerun = int(rer["n"])
        entry["defer"] = dict(
            on=bool(DEFER["on"]),
            calls=[list(x) for x in st["calls"]],
            paused=[int(p) for p in st["paused"]],
            resumed=[int(p) for p in st["resumed"]],
            found=st["found"],
            rerun_evals=rerun,
        )
        if n_p:
            STATS["steps_with_pause"] += 1
            STATS["resumes"] += n_r
            if st["found"] == "tried":
                STATS["found_tried_after_pause"] += 1
            elif st["found"] == "resume":
                STATS["found_in_resume"] += 1
            else:
                STATS["none_all_resumed"] += 1
        LOG.append(
            dict(
                voxel=int(self.cur),
                reseed=bool(self.reseed),
                **entry["defer"],
                stop=entry["stop"],
                evals=entry["evals"],
            )
        )
        self.e1_log.append(entry)
        kept = [
            dict(R=Rr, cost=float(c), hit=float(h))
            for Rr, c, h in zip(out.finisher_R, out.finisher_costs, entry["fin_hit"])
        ]
        self.last_step = dict(
            key=OD.voxel_key(vertices), voxel=int(self.cur), kept=kept, used=False
        )
        return out


def summary() -> Dict[str, Any]:
    return dict(rule=DR.RULE, defer_on=bool(DEFER["on"]), **STATS)
