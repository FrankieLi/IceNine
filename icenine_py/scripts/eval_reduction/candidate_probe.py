"""Probe for Task C (candidate overlap and de-duplication): runs one no-start search with the
existing recorder hooks and returns, per level, the candidates, their cliques, the evaluations each
quick-MC call used and the post-MC orientations.

No library change: the events are those of ``AdaptiveVoxelReconstructor.discrete_recorder`` and
``.recorder``; the FZ points and local grid of each level are captured by wrapping
``reconstructor.run_discrete_search_spaced``; evaluations are counted by wrapping
``VoxelCostFunction.evaluate`` and the batched screen / cost entry points (counts only, no timing).

The clique of a discrete-stage candidate is recovered from the "fz" events: the candidates of one
FZ point are the survivors of its spacing filter (``kept`` = local-grid indices j), and the level's
list is their concatenation sorted by cost with a stable sort. The reconstruction is checked
against the recorded orientations and scores (``local_grid[j] @ fz`` equals the recorded R).
"""

import contextlib
import io
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(ICE / "scripts" / "findoptimal_robustness"))
sys.path.insert(0, str(ICE / "scripts" / "common"))


class Probe:
    """Install around one ``reconstruct_voxel`` call; ``levels`` and ``totals`` hold the result."""

    def __init__(self, rec: Any, keep_raw: bool = True):
        self.rec = rec
        self.keep_raw = keep_raw
        self.levels: List[Dict[str, Any]] = []
        self.cnt = dict(n3=0, n0=0)
        self._cur: Dict[str, Any] = {}
        self._fz: List[Dict[str, Any]] = []
        self._qcalls: List[int] = []
        self._level_evals0 = 0

    # ------------------------------------------------------------------ install
    def install(self) -> None:
        from icenine import batched_cost as BCost
        from icenine import batched_screen as BS
        from icenine import cost_functions as CF
        from icenine import orientation_search as OS
        from icenine import reconstructor as RC

        self._mods = (BCost, BS, CF, OS, RC)
        cnt = self.cnt
        self._orig = dict(
            evaluate=CF.VoxelCostFunction.evaluate,
            exact=BS.ScreenEvaluator.overlaps_exact,
            interp=BS.ScreenEvaluator.overlaps_interp,
            many=BCost.BatchedCostEvaluator.evaluate_many,
            mc=OS.MCOptimizer.optimize,
            disc=RC.run_discrete_search_spaced,
        )
        o = self._orig

        def eval_w(this: Any, *a: Any, **k: Any) -> Any:
            cnt["n3" if this.pixel_radius > 0 else "n0"] += 1
            return o["evaluate"](this, *a, **k)

        def exact_w(this: Any, orientations: Any, *a: Any, **k: Any) -> Any:
            cnt["n3"] += len(orientations)
            return o["exact"](this, orientations, *a, **k)

        def interp_w(this: Any, centres: Any, local_grid: Any, *a: Any, **k: Any) -> Any:
            cnt["n3"] += len(centres) * len(local_grid)
            return o["interp"](this, centres, local_grid, *a, **k)

        def many_w(this: Any, orientations: Any, *a: Any, **k: Any) -> Any:
            before = cnt["n0"]
            r = o["many"](this, orientations, *a, **k)
            cnt["n0"] = before + len(orientations)  # replaces what evaluate() counted inside
            return r

        def mc_w(this: Any, **k: Any) -> Any:
            before = cnt["n0"] + cnt["n3"]
            r = o["mc"](this, **k)
            if k.get("max_mc_steps") == 10:
                self._qcalls.append(cnt["n0"] + cnt["n3"] - before)
            return r

        def disc_w(**kw: Any) -> Any:
            self._cur = dict(
                fz=np.asarray(kw["fz_orientations"]).copy(),
                local_grid=np.asarray(kw["local_grid"]).copy(),
                radius=float(kw["angular_radius"]),
                evals0=cnt["n3"] + cnt["n0"],
            )
            self._fz = []
            return o["disc"](**kw)

        CF.VoxelCostFunction.evaluate = eval_w  # type: ignore[method-assign]
        BS.ScreenEvaluator.overlaps_exact = exact_w  # type: ignore[method-assign]
        BS.ScreenEvaluator.overlaps_interp = interp_w  # type: ignore[method-assign]
        BCost.BatchedCostEvaluator.evaluate_many = many_w  # type: ignore[method-assign]
        OS.MCOptimizer.optimize = mc_w  # type: ignore[method-assign]
        RC.run_discrete_search_spaced = disc_w
        self.rec.discrete_recorder = self._disc_hook
        self.rec.recorder = self._hook

    def uninstall(self) -> None:
        BCost, BS, CF, OS, RC = self._mods
        o = self._orig
        CF.VoxelCostFunction.evaluate = o["evaluate"]  # type: ignore[method-assign]
        BS.ScreenEvaluator.overlaps_exact = o["exact"]  # type: ignore[method-assign]
        BS.ScreenEvaluator.overlaps_interp = o["interp"]  # type: ignore[method-assign]
        BCost.BatchedCostEvaluator.evaluate_many = o["many"]  # type: ignore[method-assign]
        OS.MCOptimizer.optimize = o["mc"]  # type: ignore[method-assign]
        RC.run_discrete_search_spaced = o["disc"]
        self.rec.discrete_recorder = None
        self.rec.recorder = None

    # ------------------------------------------------------------------ hooks
    def _disc_hook(self, ev: str, d: Dict[str, Any]) -> None:
        if ev == "fz":
            self._fz.append(d)

    def _hook(self, name: str, d: Dict[str, Any]) -> None:
        if name == "discrete":
            self._on_discrete(d)
        elif name == "quick_mc":
            self._on_quick(d)

    def _on_discrete(self, d: Dict[str, Any]) -> None:
        cur = self._cur
        fz_of, j_of, cost_of = [], [], []
        for ev in self._fz:
            cmap = dict(zip(ev["j"], ev["cost"]))
            for j in ev["kept"]:
                fz_of.append(ev["fz"])
                j_of.append(j)
                cost_of.append(cmap[j])
        fz_of_a, j_of_a, cost_a = np.array(fz_of), np.array(j_of), np.array(cost_of)
        order = np.argsort(cost_a, kind="stable")
        R = np.asarray(d["R"])
        # the reconstruction of the clique labels must reproduce the recorded list
        if len(order) == len(R):
            rebuilt = cur["local_grid"][j_of_a[order]] @ cur["fz"][fz_of_a[order]]
            ok = bool(np.allclose(rebuilt, R, atol=1e-12)) and bool(
                np.allclose(cost_a[order], np.asarray(d["score"]), atol=1e-12)
            )
            clique = fz_of_a[order]
        else:  # extra candidates were appended: not used by this study
            ok = False
            clique = np.full(len(R), -1)
        self._level = dict(
            level=int(d["level"]),
            n_q_max=float(d["n_q_max"]),
            diameter=float(d["diameter"]),
            n_fz=int(len(cur["fz"])),
            n_local=int(len(cur["local_grid"])),
            n_cliques=int(len({e["fz"] for e in self._fz if len(e["kept"]) > 0})),
            n_disc=int(len(R)),
            clique_ok=ok,
            screen_recost_evals=int(self.cnt["n3"] + self.cnt["n0"] - cur["evals0"]),
            clique=clique,
            R_disc=R.copy(),
            score_disc=np.asarray(d["score"]).copy(),
        )
        if self.keep_raw:
            self._level["fz"] = cur["fz"]
            self._level["local_grid"] = cur["local_grid"]
        self._qcalls = []
        self._level["_evals_at_mc_start"] = self.cnt["n3"] + self.cnt["n0"]

    def _on_quick(self, d: Dict[str, Any]) -> None:
        L = self._level
        perm = np.asarray(d["perm"])  # sorted position -> index into this level's discrete list
        n = L["n_disc"]
        Rs, cs = np.asarray(d["R"]), np.asarray(d["cost"])
        L["mc_idx"] = perm.copy()  # discrete indices that went through the quick MC
        L["n_mc"] = int(len(perm))
        L["R_after_ranked"] = Rs.copy()  # post-MC, ranked best first
        L["cost_after_ranked"] = cs.copy()
        L["n_keep"] = int(d["n_keep"])
        L["quick_evals_total"] = int(self.cnt["n3"] + self.cnt["n0"] - L.pop("_evals_at_mc_start"))
        # per-call evaluations in the order the candidates entered the quick MC
        calls = np.array(self._qcalls, dtype=np.int64)
        if len(calls) == len(perm):
            qe = np.zeros(n, dtype=np.int64)
            in_order = np.sort(perm)  # sequential loop = discrete (list) order of the survivors
            qe[in_order] = calls
            L["quick_evals"] = qe  # 0 for candidates that did not go through the quick MC
            L["quick_calls_ok"] = True
        else:
            L["quick_evals"] = np.zeros(n, dtype=np.int64)
            L["quick_calls_ok"] = False
        R_after = np.zeros((n, 3, 3))
        c_after = np.full(n, np.nan)
        R_after[perm] = Rs
        c_after[perm] = cs
        L["R_after"] = R_after
        L["cost_after"] = c_after
        for key, val in d.items():  # additive keys of the dedup study (prefix "dd_")
            if key.startswith("dd_"):
                L[key] = np.asarray(val).copy()
        self.levels.append(L)

    # ------------------------------------------------------------------ run
    def run(self, vertices: Any, phase: int, rng: Any) -> Any:
        self.install()
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                res = self.rec.reconstruct_voxel(
                    voxel_vertices=vertices, phase_index=phase, rng=rng
                )
        finally:
            self.uninstall()
        return res

    @property
    def total_evals(self) -> int:
        return int(self.cnt["n3"] + self.cnt["n0"])
