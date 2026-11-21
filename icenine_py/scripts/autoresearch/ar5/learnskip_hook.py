"""AR5 step 18 (task ar5r): the learned skip of BFS local optimizations, live (script layer; the
library is unchanged). ``LearnHook`` is step 17's ``fitskip_hook.FitHook`` (it patches the
``_expand``, ``_insert_seed`` and ``_local_fit`` methods of the ``BFSReconstruction`` INSTANCE) with
the rule ``L``:

* every neighbour fit and revisit is CHECKED: one counted evaluation of the start gives its hit
  ratio, confidence and cost (as rule ``S`` of step 17); the 16 features of
  ``learnskip_features`` follow; the model (``learnskip_model.PlainModel``) scores the fit (the
  predicted chance that it is accepted); a fit whose score is under the threshold is SKIPPED;
* a skipped fit costs no further evaluation. It is a REJECTED fit: it still takes and discards its
  aligned draw (the book's generator of (voxel, ``_local_fit``)), so no other fit's draw moves; a
  skipped neighbour fit holds its start (orientation, cost, hit ratio and confidence of the
  start); a skipped revisit keeps the old end (its confidence is returned as -1, so the push rule
  of ``_expand`` restores the old fit). A skipped fit has no fit end: it never updates a kept end;
* a fit that is not skipped runs unchanged after its check;
* a model call is NOT a cost-function evaluation (counted zero). Its wall time and the feature
  time are kept apart (``learnskip.json``).

One row per fit goes to ``<out_path>`` (jsonl): index, voxel, kind, the 16 features, the score,
whether it was skipped, and for the fits that ran the accepted flag, the end hit ratio and cost and
the evaluations of the fit; ``learnskip_gate`` compares them with the replay.
"""

import json
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, Optional, Tuple

import numpy as np

import fitskip_hook as FH
import fitskip_rules as FR
import learnskip_features as LFE
import learnskip_model as LM


class LearnHook(FH.FitHook):
    def __init__(
        self,
        bfs: Any,
        model: LM.PlainModel,
        sym: np.ndarray,
        wrap: Any = None,
        out_path: Optional[Path] = None,
        book: Any = None,
        voxel_key: Any = None,
    ) -> None:
        # the base class patches the three methods; its own rule and recording stay off
        super().__init__(
            bfs, wrap=wrap, rule=None, sym=sym, out_path=None, book=book, voxel_key=voxel_key
        )
        self.model = model
        self.fs = LFE.FeatureState(sym)
        self.rows = open(out_path, "a") if out_path is not None else None
        self.counts.update(scored=0, kept_fit_ends=0, skipped_start_passes_test=0)
        self.t_model = 0.0
        self.t_feat = 0.0
        self.i = 0

    def local_fit(
        self,
        voxel: Any,
        vertices: Any,
        rng: Any,
        rec: Any,
        phase: str,
        sigma0_deg: Optional[float] = None,
    ) -> Tuple[Any, float, float]:
        kind = FH.KIND.get(phase, phase)
        if sigma0_deg is not None or kind not in ("neighbour", "revisit"):
            # not a fit of this rule (restarts): unchanged, and not counted as a position
            return self.orig_local_fit(voxel, vertices, rng, rec, phase, sigma0_deg)
        v = int(rec.voxel_idx)
        best = self.best
        start = np.array(voxel.orientation, np.float64)
        old_conf = float(voxel.confidence)
        n_rev = int(rec.n_revisits)
        info = self.bfs.reconstructor.evaluate_overlap(start, vertices, voxel.phase)
        s_hit, s_conf = self._ratios(info)
        s_cost = float(info.cost)
        self._add_checks(rec, phase, 1)
        t0 = time.perf_counter()
        x = self.fs.features(
            v, kind, n_rev, best, start, s_hit, s_cost, s_conf, self.parent.get(v, -1), self.pos
        )
        t1 = time.perf_counter()
        score = self.model.score(x)
        skipped = self.model.skips(score)
        t2 = time.perf_counter()
        self.t_feat += t1 - t0
        self.t_model += t2 - t1
        self.counts["scored"] += 1
        result: Tuple[Any, float, float]
        if skipped:
            res = SimpleNamespace(orientation=start.copy(), cost=s_cost, overlap_info=None)
            result = (res, s_hit, s_conf if kind == "neighbour" else -1.0)
            if self.book is not None:  # the fit still takes and discards its aligned draw
                self.book.rng(self.voxel_key(vertices), "_local_fit")
            self.counts["skipped"] += 1
            evals = 0
        else:
            before = int(rec.n_evals)
            result = self.orig_local_fit(voxel, vertices, rng, rec, phase, sigma0_deg)
            evals = int(rec.n_evals) - before
        res, hit, conf = result
        best_after = max(hit, best)
        accepted = bool(best_after > 0 and (hit / best_after) > FR.ACCEPT)
        if skipped:
            # a skipped fit has no fit end (the voxel holds its start). It is a rejected fit
            # unless its start already passes the caller's acceptance test (then ``_expand``
            # accepts it as returned; counted apart, so such a case is visible)
            self.counts["skipped_start_passes_test"] += int(accepted)
            self.fs.update(v, start, accepted, False, None, None)
        else:
            kept = bool(kind != "revisit" or accepted or conf >= old_conf)
            self.fs.update(v, start, accepted, kept, hit, float(res.cost))
            self.counts["kept_fit_ends"] += int(kept)
        if self.rows is not None:
            self.rows.write(
                json.dumps(
                    dict(
                        i=self.i,
                        voxel=v,
                        kind=kind,
                        x=x,
                        score=score,
                        skipped=bool(skipped),
                        accepted=accepted,
                        end_hit=None if skipped else float(hit),
                        end_cost=None if skipped else float(res.cost),
                        evals=evals,
                    ),
                    separators=(",", ":"),
                )
                + "\n"
            )
        self.i += 1
        self.counts["fits"] += 1
        self.counts[kind] = self.counts.get(kind, 0) + 1
        self.best = best_after
        self.pos += 1
        return result

    def summary(self) -> Dict[str, Any]:
        return dict(
            rule="L",
            threshold=self.model.threshold,
            n_rounds=self.model.d.get("n_rounds"),
            counts=self.counts,
            expansions=self.exp,
            model_wall_s=self.t_model,
            feature_wall_s=self.t_feat,
            time=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        )

    def close(self) -> None:
        if self.rows is not None:
            self.rows.close()
