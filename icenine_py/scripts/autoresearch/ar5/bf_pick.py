"""AR5 step 24 (task ar5x): the seed step with rule ``BF`` (ball first; script layer; library
unchanged).

``BFPick`` is step 23's ``BallPick`` (arm ``J5DC``: step 22's ``DeferPick`` whose seed-step MC
refinement runs rule ``CB`` when ``BallPick``'s switch is on). The seed step runs in two passes
(``bf_rules`` has the rule text):

1. Ball pass: ``BallPick.e1_step`` with the ball on (exactly step 23's arm ``J5DC`` for the step).
   If it ends at 0.95 or more, that is the step. Its refinements' kept states are the first stages
   (``BP.LOG`` rows: the stage-1 evaluations, the trace of gains).
2. Fallback, only if the ball pass found no end at 0.95: ``DeferPick.e1_step`` with the ball off
   (step 22's arm ``J5D``), with the CMA fits taken from the ball pass (not run again) and each MC
   refinement run again from its candidate (the refinement has its own generator, so its first
   stage is the same: CHECKED equal to the kept one, ``BR.first_stage_equal``, and its evaluations
   are counted once, in the ball pass). The step's evaluations are the ball pass's plus the
   fallback's refinements beyond their first stages (``bf_rules.fallback_extra``).

The bookkeeping of the ball pass is undone before the fallback (one seed step in ``e1_log``, in
``DP.LOG``, in the step counter) so that the run's files have one row per seed step. The ball
pass's refinements stay in ``BP.LOG`` (``ball_steps.json.gz``); the MC traces of a step that falls
back are the fallback's (``J5D``'s). ``BF['on']`` False: the seed step is ``J5D``'s.

Draws no random number of its own. Reads no truth.
"""

from typing import Any, Dict, List

import bf_common as BC  # noqa: E402,F401  (first: guards, freeze, cache, held7)

import ball_pick as BP  # noqa: E402
import ball_rules as BR  # noqa: E402
import bf_rules as BRU  # noqa: E402
import defer_pick as DP  # noqa: E402
import numpy as np  # noqa: E402
import order_run as ORUN  # noqa: E402

BF: Dict[str, Any] = dict(on=True)
STATS: Dict[str, int] = dict(
    steps=0,
    ball_found=0,
    fallback_steps=0,
    fallback_found=0,
    fallback_none=0,
    no_candidates=0,
    stage1_checked=0,
    stage1_mismatch=0,
    ball_waste_evals=0,
    fallback_extra_evals=0,
)
LOG: List[Dict[str, Any]] = []  # one row per seed step (gz-written by bf_run)
FEATURES: List[Dict[str, Any]] = []  # one row per fallback step (written by bf_run)
_SYM: Dict[str, Any] = {}


class BFStageMismatch(AssertionError):
    """The fallback's first stage differs from the ball pass's kept one (stop rule)."""


def _sym() -> Any:
    if "s" not in _SYM:
        from icenine.ray_indexing import proper_rotations
        from icenine.symmetry import create_cubic_symmetry

        _SYM["s"] = proper_rotations(create_cubic_symmetry(3.61))
    return _SYM["s"]


def _mis(a: Any, b: Any) -> float:
    x = np.asarray(a, np.float64).reshape(3, 3)
    y = np.asarray(b, np.float64).reshape(3, 3)
    best = 180.0
    for o in _sym():
        m = x @ (np.asarray(o, np.float64) @ y).T
        c = np.clip((np.trace(m) - 1.0) / 2.0, -1.0, 1.0)
        best = min(best, float(np.degrees(np.arccos(c))))
    return best


class BFPick(BP.BallPick):
    """``BallPick`` whose seed step is rule ``BF`` when ``BF['on']`` (else ``J5D``'s)."""

    _cma_use: Any = None
    _cma_rec: Any = None

    def _fit(
        self, kind: str, rec: Any, j: int, Rj: np.ndarray, cj: float, vertices: Any, phase: int
    ) -> Dict[str, Any]:
        if kind == "C":
            if self._cma_use is not None:  # the fallback: the CMA fit of the ball pass
                f = self._cma_use[int(j)]
                return dict(f, R=np.array(f["R"], dtype=np.float64, copy=True))
            f = super()._fit(kind, rec, j, Rj, cj, vertices, phase)
            if self._cma_rec is not None:
                self._cma_rec[int(j)] = dict(f, R=np.array(f["R"], dtype=np.float64, copy=True))
            return f
        return super()._fit(kind, rec, j, Rj, cj, vertices, phase)

    def e1_step(
        self, vertices: Any, phase_index: int, n_cost: int = 5, n_refine: int = 2, tau: float = -1
    ) -> Any:
        if not BF["on"]:
            BP.BALL["on"] = False
            return super().e1_step(vertices, phase_index, n_cost, n_refine, tau)
        BP.BALL["on"] = True
        self._cma_use = None
        self._cma_rec = {}
        n_tr, n_ball, n_dlog = len(ORUN.TRACES), len(BP.LOG), len(DP.LOG)
        STATS["steps"] += 1
        out_b = super().e1_step(vertices, phase_index, n_cost, n_refine, tau)
        eb = self.e1_log[-1]
        refs = [dict(r) for r in BP.LOG[n_ball:]]
        if not eb["fits"]:
            STATS["no_candidates"] += 1
            eb["bf"] = dict(fallback=False, ball_evals=int(eb["evals"]) if "evals" in eb else 0)
            LOG.append(dict(voxel=int(self.cur), reseed=bool(self.reseed), fallback=False, evals=0))
            return out_b
        if eb["stop"] >= 0:
            STATS["ball_found"] += 1
            eb["bf"] = dict(fallback=False, ball_evals=int(eb["evals"]))
            LOG.append(
                dict(
                    voxel=int(self.cur),
                    reseed=bool(self.reseed),
                    fallback=False,
                    stop=int(eb["stop"]),
                    evals=int(eb["evals"]),
                )
            )
            return out_b
        # ---- the fallback: J5D's seed step from the kept states
        STATS["fallback_steps"] += 1
        cma = dict(self._cma_rec)
        ball_traces = ORUN.TRACES[n_tr:]
        del ORUN.TRACES[n_tr:]
        self.e1_log.pop()
        del DP.LOG[n_dlog:]
        ORUN.COUNTS["seed_steps"] -= 1
        DP.STATS["steps"] -= 1
        stage1_by_j = {int(r["j"]): int(r["stage1"]) for r in refs}
        waste = sum(int(r["proposals"]) + 1 for r in refs)
        self._cma_use = cma
        BP.BALL["on"] = False
        try:
            out_f = super().e1_step(vertices, phase_index, n_cost, n_refine, tau)
        finally:
            BP.BALL["on"] = True
            self._cma_use = None
        ef = self.e1_log[-1]
        fb_traces = ORUN.TRACES[n_tr:]
        by_ball = {int(r["j"]): r for r in refs}
        n_chk = 0
        for t in fb_traces:
            b = by_ball.get(int(t["j"]))
            assert b is not None, "the fallback ran a refinement the ball pass never ran"
            n_chk += 1
            if not BR.first_stage_equal(int(b["stage1"]), b["trace"], t["trace"]):
                STATS["stage1_mismatch"] += 1
                raise BFStageMismatch(
                    f"voxel {self.cur}: the fallback's first stage differs from the kept one "
                    f"(candidate {t['j']})"
                )
        STATS["stage1_checked"] += n_chk
        fx = BRU.fallback_extra(ef["fits"], stage1_by_j)
        total = BRU.step_evals(int(eb["evals"]), int(eb["stop"]), fx["extra"])
        # per position: the CMA fit and the ball refinement of the ball pass, and the fallback's
        # refinement beyond its first stage (the lists feed the run's per-position logs)
        c_ev = {int(f["j"]): int(f["evals"]) for f in cma.values()}
        m_ev = {int(f[1]): int(f[5]) for f in ef["fits"] if f[2] == "M"}
        new_fin: List[int] = []
        for jpos, jj in enumerate(out_f.refined):
            jj = int(jj)
            ev = c_ev.get(jj, 0)
            if jj in by_ball:
                ev += int(by_ball[jj]["evals"])
            if jj in m_ev:
                ev += m_ev[jj] - stage1_by_j[jj]
            new_fin.append(int(ev))
        out_f.finisher_evals = new_fin
        out_f.evals = int(total)
        ef["evals"] = int(total)
        found = bool(ef["stop"] >= 0)
        STATS["fallback_found" if found else "fallback_none"] += 1
        STATS["ball_waste_evals"] += int(waste)
        STATS["fallback_extra_evals"] += int(fx["extra"])
        ball_ends = [f for f in eb["fits"] if f[2] == "M"]
        miso = []
        for a in range(len(ball_ends)):
            for b2 in range(a + 1, len(ball_ends)):
                miso.append(_mis(ball_ends[a][6], ball_ends[b2][6]))
        ball_rows = [dict(hit=f[3], cost=f[4]) for f in ball_ends]
        feat = BRU.feature_row(
            len(self.e1_log) - 1,
            ball_rows,
            [dict(hit=f[3], cost=f[4]) for f in eb["fits"] if f[2] == "C"],
            found,
            miso,
        )
        feat.update(voxel=int(self.cur), reseed=bool(self.reseed))
        FEATURES.append(feat)
        ef["bf"] = dict(
            fallback=True,
            ball_evals=int(eb["evals"]),
            ball_waste=int(waste),
            extra=int(fx["extra"]),
            found=found,
            n_ball_traces=len(ball_traces),
        )
        LOG.append(
            dict(
                voxel=int(self.cur),
                reseed=bool(self.reseed),
                fallback=True,
                stop=int(ef["stop"]),
                ball_evals=int(eb["evals"]),
                extra=int(fx["extra"]),
                evals=int(total),
                found=found,
                ball_waste=int(waste),
            )
        )
        return out_f


def summary() -> Dict[str, Any]:
    return dict(rule=BRU.RULE, bf_on=bool(BF["on"]), **STATS)
