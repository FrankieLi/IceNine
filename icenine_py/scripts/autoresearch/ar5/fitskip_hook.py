"""AR5 step 17 (task ar5q): the hook on the BFS local optimizations (script layer; library
unchanged). One object patches three methods of the ``BFSReconstruction`` INSTANCE (no class, no
library file is touched): ``_expand`` and ``_insert_seed`` (to know the expansion's running best and
which voxel passed its orientation on) and ``_local_fit`` (the one call that makes a fit). Two jobs,
independent of each other:

* **recording** (arm ``J5r``): after every fit it stores a row (see ``fitskip_rules``'s module text
  and ``ROW_KEYS``). It draws no random number and changes no counted value. Its two overlap
  evaluations per fit (the start's hit ratio and the running best's hit ratio after 20 evaluations)
  are made outside the counted run, each on a fresh cost function that is not the fit's (as
  ``evaluate_overlap`` does for the seed step), and are listed apart (``rec_evals``).
* **acting** (arm ``J5C``, rules ``S`` and ``D``): a fit the rule skips is a REJECTED fit that
  still takes and discards its aligned draw (the book's generator of (voxel, ``_local_fit``)), so
  no other fit's draw moves. A skipped neighbour fit holds its start (orientation, cost, hit ratio
  and confidence of the start); a skipped revisit keeps the old end (its confidence is returned as
  -1, so the push rule of ``_expand`` restores the old fit). Every check of ``S`` is one counted
  evaluation (added to the voxel record and to the stats bucket of the phase). Rule ``C`` has no
  live form here: it would stop the CMA run after 20 evaluations, which the library cannot do
  without a change; it is replayed only (``fitskip_rules.replay``).
"""

import json
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

import fitskip_rules as FR

KIND = {"neighbor": "neighbour", "revisit": "revisit", "restart": "restart"}
ROW_KEYS = (
    "i exp pos voxel kind parent best start_R start_hit start_conf start_cost runbest n_calls "
    "hit20 evals end_R end_hit end_cost end_conf accepted kept n_revisits rec_evals"
).split()


class FitHook:
    def __init__(
        self,
        bfs: Any,
        wrap: Any = None,
        rule: Optional[FR.Rule] = None,
        sym: Optional[np.ndarray] = None,
        out_path: Optional[Path] = None,
        book: Any = None,
        voxel_key: Any = None,
    ) -> None:
        self.bfs = bfs
        self.wrap = wrap  # the run's ``Wrap`` (its CostRecorder holds the calls of the last fit)
        self.rule = rule
        if rule is not None and rule.form == "C":
            raise ValueError("rule C has no live form (replay only)")
        self.sym = sym
        self.book = book
        self.voxel_key = voxel_key
        self.hist = FR.VoxelHistory()
        self.out = open(out_path, "a") if out_path is not None else None
        self.best = 0.0
        self.parent: Dict[int, int] = {}
        self.exp = 0
        self.pos = 0
        self.i = 0
        self.counts: Dict[str, int] = dict(
            fits=0, neighbour=0, revisit=0, restart=0, rec_evals=0, checks=0, skipped=0
        )
        self.orig_local_fit = bfs._local_fit
        self.orig_expand = bfs._expand
        self.orig_insert = bfs._insert_seed
        bfs._local_fit = self.local_fit
        bfs._expand = self.expand
        bfs._insert_seed = self.insert_seed

    # ---- the three patched methods
    def expand(self, mic: Any, start_idx: int, best_conf: float, *a: Any, **k: Any) -> Any:
        self.best = float(best_conf)
        self.parent = {}
        self.exp += 1
        self.pos = 0
        return self.orig_expand(mic, start_idx, best_conf, *a, **k)

    def insert_seed(self, mic: Any, voxel_idx: int, queue: Any, *a: Any, **k: Any) -> Any:
        before = len(queue)
        out = self.orig_insert(mic, voxel_idx, queue, *a, **k)
        grew = len(queue) - before
        if grew > 0:
            for n in list(queue)[-grew:]:
                self.parent[int(n)] = int(voxel_idx)
        return out

    def _ratios(self, info: Any) -> Tuple[float, float]:
        hit, conf = self.bfs._ratios(info)
        return float(hit), float(conf)

    def _add_checks(self, rec: Any, phase: str, n: int) -> None:
        rec.n_evals += n
        self.bfs.stats[f"n_evals_{phase}"] += n
        self.counts["checks"] += n

    def local_fit(
        self,
        voxel: Any,
        vertices: Any,
        rng: Any,
        rec: Any,
        phase: str,
        sigma0_deg: Optional[float] = None,
    ) -> Tuple[Any, float, float]:
        kind = KIND.get(phase, phase)
        v = int(rec.voxel_idx)
        best = self.best
        start = np.array(voxel.orientation, np.float64)
        old_conf = float(voxel.confidence)
        n_rev = int(rec.n_revisits)
        if kind == "revisit":
            self.hist.begin(v, kind, n_rev)
        rule = self.rule
        skipped = False
        check_info: Any = None
        result: Tuple[Any, float, float]
        if rule is not None and sigma0_deg is None and kind in ("neighbour", "revisit"):
            if rule.form == "S" and FR.in_scope(kind, rule.scope):
                assert rule.t is not None
                check_info = self.bfs.reconstructor.evaluate_overlap(start, vertices, voxel.phase)
                s_hit, s_conf = self._ratios(check_info)
                self._add_checks(rec, phase, 1)
                skipped = FR.below_line(s_hit, best, rule.t)
                if skipped:
                    cost = float(check_info.cost)
                    res = SimpleNamespace(orientation=start.copy(), cost=cost, overlap_info=None)
                    result = (res, s_hit, s_conf if kind == "neighbour" else -1.0)
            elif rule.form == "D" and kind == "revisit":
                assert rule.d is not None and self.sym is not None
                skipped = self.hist.near(v, start, rule.d, self.sym)
                if skipped:
                    res = SimpleNamespace(
                        orientation=np.array(voxel.orientation, np.float64).copy(),
                        cost=float(voxel.cost),
                        overlap_info=None,
                    )
                    result = (res, 0.0, -1.0)
        if skipped:
            # the fit still takes and discards its aligned draw: no other fit's draw moves
            if self.book is not None:
                self.book.rng(self.voxel_key(vertices), "_local_fit")
            self.counts["skipped"] += 1
        else:
            result = self.orig_local_fit(voxel, vertices, rng, rec, phase, sigma0_deg)
        res, hit, conf = result
        best_after = max(hit, best)
        accepted = bool(best_after > 0 and (hit / best_after) > FR.ACCEPT)
        if self.rule is not None and self.rule.form == "D" and kind in ("neighbour", "revisit"):
            if skipped or not accepted:
                self.hist.reject(v, start)
        kept = bool(kind != "revisit" or accepted or conf >= old_conf)
        if self.out is not None and not skipped:
            self._record(
                kind, v, best, start, voxel, vertices, rec, res, hit, conf, accepted, kept, n_rev
            )
        self.best = best_after
        self.pos += 1
        return result

    # ---- recording
    def _record(
        self,
        kind: str,
        v: int,
        best: float,
        start: np.ndarray,
        voxel: Any,
        vertices: Any,
        rec: Any,
        res: Any,
        hit: float,
        conf: float,
        accepted: bool,
        kept: bool,
        n_rev: int,
    ) -> None:
        calls = list(self.wrap.cr.calls) if self.wrap is not None else []
        costs = [float(c) for c, _ in calls]
        runbest = np.minimum.accumulate(np.array(costs)).tolist() if costs else []
        evals = int(self.bfs.reconstructor.last_local_optimization_evals)
        ev = self.bfs.reconstructor.evaluate_overlap
        s_hit, s_conf = self._ratios(ev(start, vertices, voxel.phase))
        hit20 = None
        if len(calls) >= FR.CHECK_AT:
            j = int(np.argmin(np.array(costs[: FR.CHECK_AT])))
            hit20 = self._ratios(ev(np.asarray(calls[j][1], np.float64), vertices, voxel.phase))[0]
        self.counts["rec_evals"] += 2 if hit20 is not None else 1
        self.counts["fits"] += 1
        self.counts[kind] = self.counts.get(kind, 0) + 1
        row = dict(
            i=self.i,
            exp=self.exp,
            pos=self.pos,
            voxel=v,
            kind=kind,
            parent=self.parent.get(v, -1),
            best=best,
            start_R=np.round(start, 12).ravel().tolist(),
            start_hit=s_hit,
            start_conf=s_conf,
            start_cost=costs[0] if costs else None,
            runbest=[round(x, 9) for x in runbest],
            n_calls=len(calls),
            hit20=hit20,
            evals=evals,
            end_R=np.round(np.asarray(res.orientation, np.float64), 12).ravel().tolist(),
            end_hit=hit,
            end_cost=float(res.cost),
            end_conf=conf,
            accepted=accepted,
            kept=kept,
            n_revisits=n_rev,
            rec_evals=2 if hit20 is not None else 1,
        )
        self.i += 1
        assert self.out is not None
        self.out.write(json.dumps(row, separators=(",", ":")) + "\n")
        self.out.flush()

    def summary(self) -> Dict[str, Any]:
        return dict(
            rule=self.rule.name if self.rule else None,
            recording=self.out is not None,
            counts=self.counts,
            expansions=self.exp,
            time=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        )

    def close(self) -> None:
        if self.out is not None:
            self.out.close()


def read_rows(path: Path) -> List[Dict[str, Any]]:
    rows = []
    with open(path) as f:
        for line in f:
            if line.strip():
                rows.append(json.loads(line))
    return rows
