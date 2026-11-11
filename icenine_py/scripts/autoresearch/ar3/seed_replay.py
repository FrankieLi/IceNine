"""AR3 step 6, parts 0 and 1: replay of the seed step and of the six remedies at one voxel.

Nothing here changes the library. The seed step does not read the BFS state (it reads the voxel's
observed rays and uses ``default_rng(0)``), so it can be run again at any voxel. The truth is used
for analysis only (classes, errors); no remedy reads it.

Remedies (plan, step 6), each returns ``dict(passed, err, hit, cost, evals, seed_in)`` where
``evals`` counts everything the remedy spent at the voxel from scratch (the stock seed step
included when the remedy runs it or follows it, and ``seed_in`` says so) and ``passed`` is the seed
test (hit ratio >= 0.80) on the remedy's end:

  more20_2, more20_5  strict cost on the first 20 candidates; local optimization (stock) from the
                      best 2, or from the best 5
  cma_seed            the seed step with ``LocalOptimizer cma`` (CMASigma0 0.05, 1000 evaluations
                      per run), as run KE used it: cost on 5, one CMA run from each of the best 2
  wide                after a failed stock seed step: one CMA run (first step 1.5 degrees, 1000
                      evaluations) from its best end
  neighbours          one try (CMA, 0.05 degrees, 60 evaluations plus the final one) from each
                      distinct orientation among the voxel's neighbours in the final map; keep
                      rule ``hit`` (highest hit ratio, then lowest cost); the stock seed step
                      runs after them only if no try passes (``neighbours_only`` is the tries alone,
                      which is also remedy 5, ``defer``)
  outward             after a failed stock seed step: the library's outward search from its best
                      end, then (as ``_outward_reseed``) a stock local optimization from its end
"""

import time
from typing import Any, Dict, List, Optional

import boundary_common as B
import numpy as np
import seed_common as S

N_INDEX = 50  # candidates the indexer returns in the replay
N_LOG = 20  # candidates whose strict cost and hit ratio are logged


class Replayer:
    def __init__(self, env: B.Env) -> None:
        self.env = env
        self.rec = env.hooks.rec
        self.bfs = env.bfs
        thr = float(env.setup.config.min_acceleration_threshold)
        assert abs(thr - S.TEST) < 1e-12, thr
        self.maps: Dict[str, np.ndarray] = {}

    # ---- pieces
    def vertices(self, v: int) -> Any:
        return self.env.hooks.vertices(v)

    def phase(self, v: int) -> int:
        return int(self.env.mic.voxels[v].phase)

    def hit_of(self, R: np.ndarray, v: int) -> float:
        info = self.rec.evaluate_overlap(np.asarray(R, np.float64), self.vertices(v), self.phase(v))
        return float(self.bfs._ratios(info)[0])

    def score(self, R: np.ndarray, v: int) -> Dict[str, float]:
        c, h, _f = self.env.score3(R, v)
        return dict(cost=c, hit=h)

    def end_row(
        self, R: Optional[np.ndarray], v: int, cost: Optional[float] = None
    ) -> Dict[str, Any]:
        if R is None:
            return dict(err=None, hit=0.0, cost=1.0, passed=False)
        s = self.score(R, v)
        return dict(
            err=S.err1(self.env, R, v),
            hit=s["hit"],
            cost=s["cost"] if cost is None else float(cost),
            passed=bool(s["hit"] >= S.TEST),
        )

    # ---- the stock seed step (the library's own call)
    def stock(self, v: int) -> Any:
        return self.rec.index_seed(self.vertices(v), self.phase(v))

    def candidates(self, v: int) -> Any:
        from icenine import ray_indexing as RI

        sc, crystal, params, strict_cf = self.rec._index_context(self.phase(v))
        res, _att = RI.index_voxel(sc, self.vertices(v), crystal, params, N_INDEX)
        return res.R, strict_cf

    def seed_from(
        self, v: int, R: np.ndarray, strict_cf: Any, n_cost: int, n_refine: int
    ) -> Dict[str, Any]:
        """The tail of ``seed_by_indexing`` on a given candidate list: strict cost on the first
        ``n_cost``, local optimization (``refine_from_candidates``) on the ``n_refine`` lowest."""
        from icenine.orientation_search import SearchCandidate

        vert, ph = self.vertices(v), self.phase(v)
        out: Dict[str, Any] = dict(evals=0, R=None, cost=1.0, refined=[], fin=[])
        if len(R) == 0:
            return out
        k = min(n_cost, len(R))
        costs = np.array(
            [
                float(strict_cf.evaluate(np.asarray(R[j], np.float64), vert, ph).cost)
                for j in range(k)
            ]
        )
        out["evals"] = k
        out["strict_costs"] = [float(c) for c in costs]
        order = [int(j) for j in np.argsort(costs, kind="stable")[:n_refine]]
        best = None
        for j in order:
            r = self.rec.refine_from_candidates(
                [SearchCandidate(orientation=R[j].copy(), cost=float(costs[j]))],
                vert,
                ph,
                rng=np.random.default_rng(0),
            )
            ne = int(self.rec.last_eval_counts[1])
            out["evals"] += ne
            out["refined"].append(j)
            Rj = np.asarray(r.orientation, np.float64)
            out["fin"].append(dict(cand=j, evals=ne, cost=float(r.cost), R=Rj))
            if best is None or float(r.cost) < best[0]:
                best = (float(r.cost), Rj)
        assert best is not None
        out["cost"], out["R"] = best
        return out

    # ---- tries (remedy 4 and 5)
    def try_one(self, v: int, R0: np.ndarray) -> Dict[str, Any]:
        t = B.try_fit(self.env.hooks, self.bfs, v, R0, np.random.default_rng(0))
        return t

    def neighbour_starts(self, v: int, Rmap: np.ndarray) -> List[int]:
        """Distinct orientations among the voxel's neighbours in the map (merge: 0.5 degrees, the
        first in neighbour-index order is kept); the voxel's own entry is not used."""
        nb = [int(j) for j in self.env.nbr[v] if int(j) != v]
        return self.env._merge(Rmap, nb)

    def tries(self, v: int, Rmap: np.ndarray, truth: bool = True) -> Dict[str, Any]:
        starts = self.neighbour_starts(v, Rmap)
        rows: List[Dict[str, Any]] = []
        for j in starts:
            t = self.try_one(v, Rmap[j])
            rows.append(
                dict(
                    nbr=j,
                    evals=t["evals"],
                    cost=t["cost"],
                    hit=t["hit"],
                    err=S.err1(self.env, t["end"], v) if truth else None,
                    start_err=S.err1(self.env, Rmap[j], v) if truth else None,
                    R=t["end"],
                )
            )
        ev = int(sum(r["evals"] for r in rows))
        if rows:  # keep rule ``hit``: highest hit ratio, then lowest cost, then neighbour order
            k = int(
                np.lexsort(
                    (np.arange(len(rows)), [r["cost"] for r in rows], [-r["hit"] for r in rows])
                )[0]
            )
            best = rows[k]
            return dict(
                n_starts=len(rows),
                evals=ev,
                err=best["err"],
                hit=best["hit"],
                cost=best["cost"],
                passed=bool(best["hit"] >= S.TEST),
                R=best["R"],
                start_errs=[r["start_err"] for r in rows],
                end_errs=[r["err"] for r in rows],
                hits=[r["hit"] for r in rows],
            )
        return dict(
            n_starts=0,
            evals=0,
            err=None,
            hit=0.0,
            cost=1.0,
            passed=False,
            R=None,
            start_errs=[],
            end_errs=[],
            hits=[],
        )

    # ---- remedies that follow the seed step
    def wide(self, v: int, R0: np.ndarray) -> Dict[str, Any]:
        rec = self.rec
        old = rec.params.local_optimizer
        rec.params.local_optimizer = "cma"
        try:
            r = rec.local_optimization(
                voxel_vertices=self.vertices(v),
                phase_index=self.phase(v),
                initial_orientation=np.asarray(R0, np.float64),
                rng=np.random.default_rng(0),
                cma_sigma0_deg=1.5,
                cma_max_evals=1000,
            )
        finally:
            rec.params.local_optimizer = old
        return dict(
            R=np.asarray(r.orientation, np.float64),
            cost=float(r.cost),
            evals=int(rec.last_local_optimization_evals),
        )

    def outward(self, v: int, R0: np.ndarray) -> Dict[str, Any]:
        res = self.rec.outward_search(self.vertices(v), self.phase(v), np.asarray(R0, np.float64))
        ev = int(res.n_strict + res.n_screen)
        if not res.accepted:
            return dict(R=None, cost=1.0, evals=ev, accepted=False)
        r = self.rec.local_optimization(
            voxel_vertices=self.vertices(v),
            phase_index=self.phase(v),
            initial_orientation=np.asarray(res.orientation, np.float64),
            rng=np.random.default_rng(0),
        )
        ev += int(self.rec.last_local_optimization_evals)
        return dict(
            R=np.asarray(r.orientation, np.float64), cost=float(r.cost), evals=ev, accepted=True
        )

    def cma_seed(self, v: int, R: np.ndarray, strict_cf: Any) -> Dict[str, Any]:
        p = self.rec.params
        old = (p.local_optimizer, p.cma_sigma0_deg, p.cma_max_evals)
        p.local_optimizer, p.cma_sigma0_deg, p.cma_max_evals = "cma", 0.05, 1000
        try:
            return self.seed_from(v, R[:20], strict_cf, 5, 2)
        finally:
            p.local_optimizer, p.cma_sigma0_deg, p.cma_max_evals = old

    # ---- the whole voxel
    def voxel(
        self,
        v: int,
        maps: Dict[str, np.ndarray],
        stored: Dict[str, Dict[str, Any]],
        remedies: bool = True,
    ) -> Dict[str, Any]:
        """``maps``: run name -> its final map (n, 3, 3); ``stored``: run name -> dict(R, hit) the
        stored end of the voxel in that run (the grid search's end at a fallback voxel)."""
        env = self.env
        t0 = time.time()
        out: Dict[str, Any] = dict(
            voxel=int(v),
            grain=int(env.gmap[v]),
            isA=bool(env.isA[v]),
            level=None if env.level is None else float(env.level[v]),
        )
        Rt = env.Rtrue[v]
        out["truth"] = self.score(Rt, v)
        # the stock seed step, replayed
        s = self.stock(v)
        st_hit = self.hit_of(s.orientation, v) if s.orientation is not None else 0.0
        out["stock"] = dict(
            n_cand=int(s.n_cand),
            evals=int(s.evals),
            n_att=int(s.n_att),
            strict_costs=[float(c) for c in s.strict_costs],
            refined=[int(j) for j in s.refined],
            fin_costs=[float(c) for c in s.finisher_costs],
            fin_evals=[int(e) for e in s.finisher_evals],
            fin_err=[S.err1(env, r, v) for r in s.finisher_R],
            fin_hit=[self.hit_of(r, v) for r in s.finisher_R],
            hit=st_hit,
            passed=bool(s.orientation is not None and st_hit >= S.TEST),
            err=None if s.orientation is None else S.err1(env, s.orientation, v),
            cost=float(s.cost),
        )
        # the 50-candidate list
        R50, strict_cf = self.candidates(v)
        n = len(R50)
        errs = env.mis(R50, np.repeat(Rt[None], n, 0)).tolist() if n else []
        rank1 = next((i for i, e in enumerate(errs) if e <= B.WRONG), None)
        first = []
        for j in range(min(N_LOG, n)):
            info = strict_cf.evaluate(
                np.asarray(R50[j], np.float64), self.vertices(v), self.phase(v)
            )
            first.append(
                dict(cost=float(info.cost), hit=float(self.bfs._ratios(info)[0]), err=errs[j])
            )
        out["cand"] = dict(n=int(n), errs=[float(e) for e in errs], rank1=rank1, first=first)
        if n >= min(N_LOG, s.n_cand) and s.n_cand:
            k = min(20, n, s.n_cand)
            out["cand"]["same_first20"] = bool(
                np.allclose(R50[:k], np.asarray(s.candidates)[:k], atol=1e-12)
            )
        # class (plan, part 0 item 4)
        refined = out["stock"]["refined"]
        if out["stock"]["passed"]:
            cls = "f"
        elif out["truth"]["hit"] < S.TEST:
            cls = "a"
        elif rank1 is None:
            cls = "b"
        elif rank1 >= 5:
            cls = "c"
        elif rank1 not in refined:
            cls = "d"
        else:
            cls = "e"
        out["class"] = cls
        # the stored end(s)
        out["stored"] = {}
        for run, d in stored.items():
            Rm = maps.get(run)
            nb = [int(j) for j in env.nbr[v] if int(j) != v]
            nerr = (
                env.mis(Rm[nb], np.repeat(Rt[None], len(nb), 0)) if (Rm is not None and nb) else []
            )
            out["stored"][run] = dict(
                err=S.err1(env, d["R"], v),
                hit=d["hit"],
                n_nbr=len(nb),
                n_nbr_within1=int(np.sum(np.asarray(nerr) <= B.WRONG)),
                n_distinct=len(self.neighbour_starts(v, Rm)) if Rm is not None else 0,
            )
        if not remedies:
            out["wall_s"] = time.time() - t0
            return out
        rem: Dict[str, Dict[str, Any]] = {}
        S_ev = int(s.evals)
        # remedy 1
        for nm, nr in (("more20_2", 2), ("more20_5", 5)):
            r = self.seed_from(v, R50[:20], strict_cf, 20, nr)
            rem[nm] = dict(self.end_row(r["R"], v, r["cost"]), evals=int(r["evals"]), seed_in=False)
        # remedy 2
        r = self.cma_seed(v, R50, strict_cf)
        rem["cma_seed"] = dict(
            self.end_row(r["R"], v, r["cost"]), evals=int(r["evals"]), seed_in=False
        )
        # remedies 3 and 6 follow the stock seed step when it fails
        best_R = s.orientation
        for nm in ("wide", "outward"):
            if out["stock"]["passed"]:
                rem[nm] = dict(self.end_row(best_R, v), evals=S_ev, seed_in=True, ran=False)
            elif best_R is None:
                rem[nm] = dict(
                    err=None, hit=0.0, cost=1.0, passed=False, evals=S_ev, seed_in=True, ran=False
                )
            else:
                r = self.wide(v, best_R) if nm == "wide" else self.outward(v, best_R)
                row = self.end_row(r["R"], v, r["cost"])
                rem[nm] = dict(row, evals=S_ev + int(r["evals"]), seed_in=True, ran=True)
        # remedies 4 and 5 (per final map)
        rem["neighbours"] = {}
        rem["neighbours_only"] = {}
        for run, Rm in maps.items():
            t = self.tries(v, Rm)
            base = dict(err=t["err"], hit=t["hit"], cost=t["cost"], passed=t["passed"])
            rem["neighbours_only"][run] = dict(
                base,
                evals=t["evals"],
                n_starts=t["n_starts"],
                start_errs=t["start_errs"],
                end_errs=t["end_errs"],
                hits=t["hits"],
                seed_in=False,
            )
            if t["passed"]:
                rem["neighbours"][run] = dict(base, evals=t["evals"], seed_in=False, seed_ran=False)
            else:  # the stock seed step runs after the tries (its result is the one replayed above)
                row = (
                    self.end_row(best_R, v)
                    if best_R is not None
                    else dict(err=None, hit=0.0, cost=1.0, passed=False)
                )
                rem["neighbours"][run] = dict(
                    row, evals=t["evals"] + S_ev, seed_in=True, seed_ran=True
                )
        out["remedies"] = rem
        out["wall_s"] = time.time() - t0
        return out
