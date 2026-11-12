"""AR5 step 8, Phase A: one whole-layer run of the interleaved or the sweep order (script layer).

  RESEED_SAMPLE=<big|grad|mid|held> uv run python scripts/autoresearch/ar5/order_run.py \
      --arm A|Y|Yl|Al [--cap-min 240] [--smoke N] [--guards-only]

Arms (the run ``S11`` of step 3, the merged keys ``SeedIndexCMAFirst 200`` and
``BFSNeighborCMANoFallback 1``, without the ``N1`` rule; D-20261112-N8):

* ``A``   interleaved order (CMA on a candidate, then MC if its hit ratio is under 0.95),
  aligned draws;
* ``Y``   sweep order (CMA on every candidate first, MC only if no CMA end reaches 0.95),
  aligned draws;
* ``Yl``  the sweep order with the library's own draws (one shared generator);
* ``Al``  the interleaved order with the library's own draws (a check: it must equal step 3's
  ``S11``).

It is step 2's ``fit_run.py`` with the seed step's refinement loop replaced (``OrderPick``) and, for
the aligned arms, the generator of every neighbour/revisit fit and every grid fallback replaced by
``order_draws`` (voxel, kind of call, count at that voxel). Nothing under ``icenine/`` changes.
Besides the files of ``fit_run.py`` it writes, in the run folder: ``order.json`` (counts of the
run), and ``mc_traces.json.gz`` (one row per MC refinement: the best cost after each improvement
with the evaluation count; recorded, not read in this step). The ``e1_log.json.gz`` rows carry the
fits of the step (candidate position, kind C or M, hit ratio, cost, evaluations, end orientation).
A run that reaches the cap writes ``capped.json`` and exits with code 3. The boundary pass is made
after by ``order_pass.py``. Refuses to start when the run folder exists; ``held`` is refused until
``ar5/order/frozen.json`` is committed (``order_common``) and ``order_freeze.py --check`` agrees.
"""

import argparse
import json
import subprocess
import sys
import time
from typing import Any, Dict, List, Optional

import order_common as OC  # noqa: E402  (first: sets the cache, installs the arms)

import fit_common as FC  # noqa: E402
import fit_rules as FR  # noqa: E402
import fit_run as FRUN  # noqa: E402
import mix_common as MC  # noqa: E402
import numpy as np  # noqa: E402
import order_draws as OD  # noqa: E402
import order_rules as OR  # noqa: E402
import reseed_common as RC  # noqa: E402

ORD: Dict[str, Any] = {"order": "A", "aligned": True, "book": None}
ARMS: Dict[str, Dict[str, Any]] = {
    "A": dict(order="A", aligned=True),
    "Y": dict(order="Y", aligned=True),
    "Yl": dict(order="Y", aligned=False),
    "Al": dict(order="A", aligned=False),
}
COUNTS: Dict[str, int] = dict(
    seed_steps=0, cma_fits=0, mc_fits=0, steps_cma_stop=0, steps_mc_phase=0, steps_none=0
)
TRACES: List[Dict[str, Any]] = []


class Tracer:
    """Records the running best cost of the strict cost calls while ``on`` (one MC refinement).
    It draws no random number and changes no value."""

    def __init__(self) -> None:
        from icenine.cost_functions import VoxelCostFunction

        self.on = False
        self.n = 0
        self.best = float("inf")
        self.trace: List[List[float]] = []
        orig = VoxelCostFunction.evaluate
        me = self

        def patched(self_: Any, orientation: Any, *args: Any, **kw: Any) -> Any:
            out = orig(self_, orientation, *args, **kw)
            if me.on:
                me.n += 1
                c = float(out.cost)
                if c < me.best:
                    me.best = c
                    me.trace.append([me.n, c])
            return out

        VoxelCostFunction.evaluate = patched

    def start(self) -> None:
        self.on, self.n, self.best, self.trace = True, 0, float("inf"), []

    def stop(self) -> Dict[str, Any]:
        self.on = False
        return dict(n_calls=self.n, trace=self.trace)


class OrderPick(FRUN.FitPick):
    """``FitPick`` with the refinement loop of the seed step in the order ``ORD['order']``."""

    tracer: Optional[Tracer] = None

    def _fit(
        self, kind: str, rec: Any, j: int, Rj: np.ndarray, cj: float, vertices: Any, phase: int
    ) -> Dict[str, Any]:
        if OrderPick.tracer is None:
            OrderPick.tracer = Tracer()
        tr = OrderPick.tracer
        if kind == "M":
            tr.start()
        Rend, c, ev = MC._orig_refine_arm(
            "C200" if kind == "C" else "M", rec, self.orig_lo, Rj, cj, vertices, phase
        )
        t = tr.stop() if kind == "M" else None
        info = rec.evaluate_overlap(np.asarray(Rend, np.float64), vertices, phase)
        hit = float(self.bfs._ratios(info)[0])
        COUNTS["cma_fits" if kind == "C" else "mc_fits"] += 1
        f = dict(j=int(j), kind=kind, hit=hit, cost=float(c), evals=int(ev), R=np.asarray(Rend))
        if t is not None:
            TRACES.append(
                dict(
                    voxel=self.cur,
                    step=COUNTS["seed_steps"],
                    j=int(j),
                    evals=int(ev),
                    n_calls=t["n_calls"],
                    trace=t["trace"],
                )
            )
        return f

    def e1_step(
        self, vertices: Any, phase_index: int, n_cost: int = 5, n_refine: int = 2, tau: float = -1
    ) -> Any:
        from icenine import ray_indexing as RI

        tau = self.tau if tau < 0 else tau
        rec = self.rec
        t0 = time.perf_counter()
        COUNTS["seed_steps"] += 1
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
            order=ORD["order"],
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
        fits: List[Dict[str, Any]] = []

        def run(kind: str, p: int) -> Dict[str, Any]:
            j = order[p]
            f = self._fit(kind, rec, j, R[j], float(costs[j]), vertices, phase_index)
            f["pos"] = p
            fits.append(f)
            return f

        kept: List[Dict[str, Any]] = []
        stop: Optional[int] = None
        if ORD["order"] == "A":
            for p in range(len(order)):
                c = run("C", p)
                e = c
                if c["hit"] < tau:
                    m = run("M", p)
                    if m["cost"] < c["cost"]:
                        e = m
                kept.append(dict(e, evals=sum(f["evals"] for f in fits if f["pos"] == p)))
                if e["hit"] >= tau:
                    stop = p
                    break
        else:
            cmas: Dict[int, Dict[str, Any]] = {}
            mcs: Dict[int, Dict[str, Any]] = {}

            def cma_fit(p: int) -> Any:
                cmas[p] = run("C", p)
                return cmas[p]["hit"], cmas[p]["cost"], cmas[p]["evals"]

            def mc_fit(p: int) -> Any:
                mcs[p] = run("M", p)
                return mcs[p]["hit"], mcs[p]["cost"], mcs[p]["evals"]

            r = OR.sweep_choose(len(order), cma_fit, mc_fit, tau)
            for p in range(len(cmas)):
                e = cmas[p]
                if p in mcs and mcs[p]["cost"] < e["cost"]:
                    e = mcs[p]
                ev = cmas[p]["evals"] + (mcs[p]["evals"] if p in mcs else 0)
                kept.append(dict(e, evals=ev))
            stop = r["taken"] if r["stopped"] else None
            COUNTS[
                {"cma": "steps_cma_stop", "mc": "steps_mc_phase"}.get(r["phase"], "steps_none")
            ] += 1
            assert len(cmas) == r["n_cma"] and len(mcs) == r["mc_ran"] or r["phase"] == "none"
        for p, e in enumerate(kept):
            out.refined.append(order[p])
            out.finisher_costs.append(float(e["cost"]))
            out.finisher_evals.append(int(e["evals"]))
            out.finisher_R.append(np.asarray(e["R"], dtype=np.float64))
            entry["fin_evals"].append(int(e["evals"]))
            entry["fin_hit"].append(float(e["hit"]))
        out.evals += int(sum(out.finisher_evals))
        best, stopped = FR.early_stop(entry["fin_hit"], out.finisher_costs, tau)
        assert (stop is None) == (not stopped) and (stop is None or stop == best)
        if ORD["order"] == "Y":
            assert best == r["taken"]
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
        self.e1_log.append(entry)
        return out


class OrderWrap(FRUN.FitWrap):
    """``FitWrap`` (the neighbour/revisit CMA fit without the MC fallback); for the aligned arms
    the generator of each call is the aligned one, and so is the grid fallback's."""

    def __init__(self, bfs: Any) -> None:
        super().__init__(bfs)
        self.book = OD.DrawBook()
        ORD["book"] = self.book
        if ORD["aligned"]:
            self.rec.reconstruct_voxel = OD.aligned_reconstruct_voxel(
                self.book, self.rec.reconstruct_voxel
            )

    def wrapped_no_fallback(
        self,
        voxel_vertices: Any,
        phase_index: int,
        initial_orientation: np.ndarray,
        rng: Any = None,
        **kw: Any,
    ) -> Any:
        if ORD["aligned"]:
            rng = self.book.rng(OD.voxel_key(voxel_vertices), sys._getframe(1).f_code.co_name)
        return self._no_fallback(voxel_vertices, phase_index, initial_orientation, rng, **kw)

    def _no_fallback(
        self, voxel_vertices: Any, phase_index: int, R0_in: np.ndarray, rng: Any, **kw: Any
    ) -> Any:
        """``FitWrap.wrapped_no_fallback`` with the caller taken two frames up (the counter of
        wrapper calls by caller stays the one of the unchanged run)."""
        assert not kw, "the wrapper expects the classic-mode call (no cma keywords)"
        self.calls[sys._getframe(2).f_code.co_name] += 1
        R0 = np.asarray(R0_in, np.float64)
        r1, n1, c0 = self.cma_fit(voxel_vertices, phase_index, R0, rng)
        self.evals_picked += n1
        if r1.cost < c0:
            self.n_picked_fell += 1
        else:
            self.cma_not_fell += 1
        self.rec.last_local_optimization_evals = n1
        return r1


def library_state() -> Dict[str, Any]:
    h = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=str(OC.ICE), capture_output=True, text=True
    ).stdout.strip()
    d = subprocess.run(
        ["git", "status", "--porcelain", "--", "icenine", "scripts/autoresearch"],
        cwd=str(OC.ICE),
        capture_output=True,
        text=True,
    ).stdout.strip()
    return dict(head=h, dirty_library_or_scripts=d.splitlines()[:20])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", required=True, choices=list(ARMS))
    ap.add_argument("--cap-min", type=float, default=240.0)
    ap.add_argument("--smoke", type=int, default=0)
    ap.add_argument("--guards-only", action="store_true")
    ap.add_argument("--name", default="", help="run folder name (default: the arm)")
    a = ap.parse_args()
    if FC.SAMPLE == "held":
        import order_freeze as OF

        OF.check_or_exit()
    if a.guards_only:
        print("guards passed")
        return
    name = a.name or a.arm
    ORD.update(ARMS[a.arm])
    out = FC.CACHE / f"{name}_{FC.SAMPLE}{'_smoke' if a.smoke else ''}"
    if out.exists() and any(out.iterdir()):
        raise SystemExit(f"refusing to start: {out} already exists")
    FRUN.FitPick = OrderPick
    FRUN.FitWrap = OrderWrap
    sys.argv = [
        "fit_run.py",
        "--name",
        name,
        "--seed-arm",
        "X200",
        "--no-mc-fallback",
        "--cap-min",
        str(a.cap_min),
    ] + (["--smoke", str(a.smoke)] if a.smoke else [])
    t0 = time.time()
    FRUN.main()
    book: OD.DrawBook = ORD["book"]
    by_kind: Dict[str, int] = {}
    for (_k, kind), n in book.counts.items():
        by_kind[kind] = by_kind.get(kind, 0) + n
    RC.write_gz(TRACES, out / "mc_traces.json.gz")
    (out / "order.json").write_text(
        json.dumps(
            dict(
                arm=a.arm,
                name=name,
                sample=FC.SAMPLE,
                order=ORD["order"],
                aligned_draws=bool(ORD["aligned"]),
                aligned_calls=book.n_issued,
                aligned_calls_by_kind=by_kind,
                counts=COUNTS,
                mc_traces=len(TRACES),
                mc_traces_with_trace_eq_evals=int(
                    sum(1 for t in TRACES if t["n_calls"] == t["evals"])
                ),
                wall_s=time.time() - t0,
                wall_label="contended (other jobs on the machine); not a timing result",
                state=library_state(),
            ),
            indent=1,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
