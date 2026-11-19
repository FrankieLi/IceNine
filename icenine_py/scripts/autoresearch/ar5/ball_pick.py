"""AR5 step 23 (task ar5w): the seed step with rule ``CB`` (script layer; library unchanged).

``BallPick`` is step 22's ``DeferPick`` (arm ``J5D``: interleaved order, 20 costed, 5 refined,
rules ``PH(500)`` and ``DEF`` in the MC refinement) whose MC refinement of the seed step has its
second stage (the library's ``MCOptimizer.variance_minimizing_optimize``, the variance stage)
replaced by ``ball_rules.ball_stage``:

* ``BALL['on']`` True: arm ``J5DC``. While a seed step's MC refinement runs (``BALL['active']``,
  set only inside ``_fit_m``) the method is the converging ball; at any other time (the grid
  search's own refinement, a rescue, BFS) it is the library's, untouched.
* ``BALL['on']`` False: ``DeferPick`` exactly (arm ``J5D``; gate g0's first half).

Stage 1 is the library's, so it runs exactly as in ``J5D``; the generator is the refinement's own
(``default_rng(0)`` in ``fit_common.refine_arm``) and the ball continues it. Every evaluation of
the ball goes through the cost function the tracer and the stall meter wrap, so every evaluation
is counted. Under ``CB`` the stall check and the pause must never fire (the stage has at most
about 400 evaluations with no gain): ``BallCheckFired`` is raised if one does (stop rule).

Draws no random number of its own beyond the stage's four per proposal.
"""

from typing import Any, Dict, List, Optional

import ball_common as BC  # noqa: E402,F401  (first: guards, freeze, cache, held7)

import ball_rules as BR  # noqa: E402
import defer_pick as DP  # noqa: E402
import high500_pick as HP  # noqa: E402
import numpy as np  # noqa: E402
import order_run as ORUN  # noqa: E402

BALL: Dict[str, Any] = dict(on=True, active=False)
STATS: Dict[str, int] = dict(
    refinements=0,
    ends_floor=0,
    ends_cap=0,
    proposals=0,
    gains=0,
    halvings=0,
    stage1_evals=0,
    stage2_evals=0,
    checks=0,
    pauses=0,
)
RECORD: List[Dict[str, Any]] = []  # filled by the stage (one row per call)
LOG: List[Dict[str, Any]] = []  # one row per seed-step MC refinement (gz-written by ball_run)


class BallCheckFired(AssertionError):
    """A stall check or a pause fired under rule CB (stop rule of the plan)."""


_ORIG: Any = None


def install() -> None:
    """Replace ``MCOptimizer.variance_minimizing_optimize`` by the guarded ball (once)."""
    global _ORIG
    from icenine.orientation_search import MCOptimizer, SearchCandidate

    if _ORIG is not None:
        return
    _ORIG = MCOptimizer.variance_minimizing_optimize

    def patched(
        self: Any, initial_orientation: Any, search_box_side: float, *args: Any, **kw: Any
    ) -> Any:
        if not BALL["active"]:
            return _ORIG(self, initial_orientation, search_box_side, *args, **kw)
        tr = HP.High500Pick.ptracer
        assert tr is not None and tr.on, "the ball runs only inside a traced seed-step refinement"
        stage1 = int(tr.n)

        def ev(R: np.ndarray) -> Any:
            info = self.cost_fn.evaluate(R, self.voxel_vertices, self.phase_index)
            return float(info.cost), info

        out = BR.ball_stage(
            ev, np.asarray(initial_orientation, np.float64), search_box_side, self._rng
        )
        RECORD.append(
            dict(
                stage1=stage1,
                proposals=int(out["proposals"]),
                gains=int(out["gains"]),
                halvings=int(out["halvings"]),
                end=out["end"],
                r_final_deg=float(np.degrees(out["r_final"])),
                start_cost=float(out["start_cost"]),
                end_cost=float(out["cost"]),
                box_deg=float(np.degrees(search_box_side)),
            )
        )
        return SearchCandidate(orientation=out["R"], cost=out["cost"], overlap_info=out["info"])

    MCOptimizer.variance_minimizing_optimize = patched  # type: ignore[assignment,method-assign]


class BallPick(DP.DeferPick):
    """``DeferPick`` whose seed-step MC refinement runs rule ``CB`` when ``BALL['on']``."""

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
        if not BALL["on"]:
            return super()._fit_m(rec, j, Rj, cj, vertices, phase, allow_pause, prior)
        assert prior is None, "rule CB never resumes a refinement"
        install()
        del RECORD[:]
        BALL["active"] = True
        try:
            f = super()._fit_m(rec, j, Rj, cj, vertices, phase, allow_pause, prior)
        finally:
            BALL["active"] = False
        tr = HP.High500Pick.ptracer
        assert tr is not None
        n_checks = len(tr.checks)
        if f.get("paused") or n_checks:
            STATS["checks"] += n_checks
            STATS["pauses"] += int(bool(f.get("paused")))
            raise BallCheckFired(
                f"voxel {self.cur}: a stall check ({n_checks}) or a pause fired under CB"
            )
        assert len(RECORD) == 1, "the ball stage ran once per refinement"
        rc = dict(RECORD[0])
        sp = BR.split_refinement(int(f["evals"]), rc["stage1"], rc["proposals"])
        assert sp["final"] == 1, f"the final evaluation is one evaluation, got {sp}"
        row = dict(
            voxel=int(self.cur),
            step=int(ORUN.COUNTS["seed_steps"]),
            j=int(j),
            evals=int(f["evals"]),
            trace=[list(x) for x in ORUN.TRACES[-1]["trace"]],
            **rc,
        )
        LOG.append(row)
        STATS["refinements"] += 1
        STATS["ends_" + rc["end"]] += 1
        STATS["proposals"] += rc["proposals"]
        STATS["gains"] += rc["gains"]
        STATS["halvings"] += rc["halvings"]
        STATS["stage1_evals"] += rc["stage1"]
        STATS["stage2_evals"] += sp["stage2"]
        return f


def summary() -> Dict[str, Any]:
    return dict(rule=BR.RULE, ball_on=bool(BALL["on"]), **STATS)
