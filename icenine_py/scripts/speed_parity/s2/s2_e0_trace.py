"""S2 Step 6: trace the E0 wrong flips of the key FusedCost (D-20261014-N13: a break counts
against the key only when traced to it).

For each voxel whose wrong flag differs between a key-off arm and its key-on arm (``off`` vs ``on``,
or ``bqmc`` vs ``bqmc_on``), the key-on search is rerun in one process with a shadow: every fused
result (``VoxelCostFunction.evaluate`` with the key on, and ``BatchedCostEvaluator.evaluate_many``
on its fused entry) is recomputed on the torch path (``evaluate`` with the key off;
``evaluate_many(variant="compact")``) without counting it, and the search goes on with the fused
result. Up to the first evaluation where the two differ, the key-on search makes exactly the key-off
search's calls (the key changes nothing else, and the shadow consumes no random number), so the
first differing evaluation is where the two searches can first part. It is traced with
``s2_parity.trace(every=True)`` (``fused_debug`` intermediates against ``torch_mirror``) and
checked with the parity rule of ``s2_lib_parity.check``. The rerun must reproduce the key-on arm's
stored record (final orientation and evaluation counts), else the trace does not describe that run.

``--changed`` also traces every voxel whose final orientation differs between the two arms while
its wrong flag does not (change "both right" / "both wrong"): the search took another path.

Per voxel: the pair, the flag change (fixed / broken), the number of shadowed and differing
evaluations, the first differing one (call number, path, pixel radius, which ``OverlapInfo`` fields
differ, the traced entries and the rule verdict), and whether the rerun reproduced the record.

  uv run python scripts/speed_parity/s2/s2_e0_trace.py --cache scripts/findoptimal_robustness/cache \
      --out benchmarks/s2_results/e0_trace.json
"""

import argparse
import contextlib
import io
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
sys.path.insert(0, str(ICE))
sys.path.insert(0, str(ICE / "scripts" / "findoptimal_robustness"))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))

import common as C  # noqa: E402
import e0_run as E  # noqa: E402
import s1_parity as S1  # noqa: E402
import s2_lib_parity as LP  # noqa: E402
import s2_parity as SP  # noqa: E402
import s2fused as SF  # noqa: E402
import torch_mirror as TM  # noqa: E402

from icenine import batched_cost as BC  # noqa: E402
from icenine import cost_functions as CF  # noqa: E402

#: arm -> (cache tag, BatchedQuickMC, FusedCost)
ARMS = {
    "off": ("_s2off", False, False),
    "on": ("_s2on", False, True),
    "bqmc": ("_s2bq", True, False),
    "bqmc_on": ("_s2bqon", True, True),
}
PAIRS = (("off", "on"), ("bqmc", "bqmc_on"))


class _Ds:
    """The ``ds`` of ``s2_lib_parity.near_px`` (detector widths)."""

    def __init__(self, setup: Any) -> None:
        self.setup = setup


class Shadow:
    def __init__(self) -> None:
        self.active = False
        self.n = 0
        self.n_diff = 0
        self.first: Optional[Dict[str, Any]] = None
        self.ds: Any = None

    def record(
        self, path: str, cf: Any, R: np.ndarray, verts: Any, phase: int, f: Any, t: Any
    ) -> None:
        self.n += 1
        if f == t:
            return
        self.n_diff += 1
        if self.first is not None:
            return
        fields = [SF.FIELDS[j] for j in range(len(f)) if f[j] != t[j]]
        with contextlib.redirect_stdout(io.StringIO()):
            m = TM.mirror(cf, R, verts, phase)
            d = cf.fused_debug(R, verts, phase)
            ref = TM.reference64(cf, R, verts, phase)
        if tuple(m["result"]) != t or tuple(d["result"]) != f:
            raise SystemExit("shadow: mirror or fused_debug does not reproduce the results")
        entries = SP.trace(m, d, ref, cf, len(cf.detector_list), every=True)
        bad = LP.check(entries, self.ds, 0.01) if entries else ["no traced entry"]
        self.first = dict(
            call=self.n,
            path=path,
            pixel_radius=int(cf.pixel_radius),
            fields=fields,
            torch=list(t),
            fused=list(f),
            entries=entries,
            within_rule=not bad,
            failures=bad,
            max_near_px=max((LP.near_px(e, self.ds) for e in entries), default=None),
        )


SH = Shadow()
_eval = CF.VoxelCostFunction.evaluate
_many = BC.BatchedCostEvaluator.evaluate_many


def evaluate(self: Any, orientation: Any, voxel_vertices: Any, phase_index: int = 0) -> Any:
    r = _eval(self, orientation, voxel_vertices, phase_index)
    if SH.active and self.fused_cost and CF.fused_unsupported_reason(self) is None:
        n = self.eval_count
        self.fused_cost = False
        try:
            t = _eval(self, orientation, voxel_vertices, phase_index)
        finally:
            self.fused_cost = True
            self.eval_count = n
        R = np.asarray(orientation, dtype=np.float64)
        SH.record(
            "evaluate",
            self,
            R,
            voxel_vertices,
            phase_index,
            SF.result_tuple(r),
            SF.result_tuple(t),
        )
    return r


def evaluate_many(
    self: Any, orientations: Any, voxel_vertices: Any, phase_index: int = 0, variant: Any = None
) -> Any:
    n_fused = self.counts["orient_fused"]
    r = _many(self, orientations, voxel_vertices, phase_index, variant)
    if SH.active and self.counts["orient_fused"] > n_fused:
        n = self.cf.eval_count
        t = _many(self, orientations, voxel_vertices, phase_index, "compact")
        self.cf.eval_count = n
        R = np.asarray(orientations, dtype=np.float64)
        for b in range(len(R)):
            SH.record(
                "evaluate_many",
                self.cf,
                R[b],
                voxel_vertices,
                phase_index,
                SF.result_tuple(r[b]),
                SF.result_tuple(t[b]),
            )
    return r


def wrong(path: Path) -> Optional[bool]:
    d = np.load(path)
    if "unbuildable" in d.files:
        return None
    return bool(float(C.err_deg(d["s0_R_final"], d["R_true"])) > C.WRONG_DEG)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cache", default=str(C.CACHE_DIR))
    ap.add_argument("--out", required=True)
    ap.add_argument("--variant", default="clean")
    ap.add_argument("--limit", type=int, default=0, help="first N flips per pair (pilot)")
    ap.add_argument(
        "--changed",
        action="store_true",
        help="also every voxel whose final orientation differs while its wrong flag does not",
    )
    a = ap.parse_args()
    out_path = Path(a.out).resolve()
    cache = Path(a.cache).resolve()
    CF.VoxelCostFunction.evaluate = evaluate  # type: ignore[method-assign]
    BC.BatchedCostEvaluator.evaluate_many = evaluate_many  # type: ignore[method-assign]
    info = E.load_info()
    vox = [int(v) for v in info["voxel_indices"]]
    flips: List[Dict[str, Any]] = []
    for off, on in PAIRS:
        n = 0
        for vpos, v in enumerate(vox):
            name = f"v{v}_{a.variant}.npz"
            pa, pb = cache / f"e0{ARMS[off][0]}" / name, cache / f"e0{ARMS[on][0]}" / name
            if not (pa.exists() and pb.exists()):
                continue
            wa, wb = wrong(pa), wrong(pb)
            if wa is None or wb is None:
                continue
            moved = not np.array_equal(np.load(pa)["s0_R_final"], np.load(pb)["s0_R_final"])
            if wa == wb and not (a.changed and moved):
                continue
            flips.append(
                dict(pair=f"{off} vs {on}", arm=on, vidx=v, vpos=vpos, wrong_off=wa, wrong_on=wb)
            )
            n += 1
            if a.limit and n >= a.limit:
                break
    print(f"{len(flips)} voxels to trace", flush=True)
    C.init_worker(C.worker_args())
    W = C.get_worker()
    SH.ds = _Ds(W.rec.setup)
    res: List[Dict[str, Any]] = []
    for f in flips:
        t0 = time.perf_counter()
        tag, bqmc, fused = ARMS[f["arm"]]
        case = E.get_case(f["vidx"], f["vpos"], a.variant)
        assert case is not None
        keys, R_true, vctx = case
        C.attach(keys)
        rec = W.rec
        sp = rec.setup.search_params
        sp.batched_quick_mc, sp.search_device = bqmc, "cpu"
        sp.cpp_center_screen, sp.batched_screen, sp.fused_cost = True, True, fused
        SH.active, SH.n, SH.n_diff, SH.first = True, 0, 0, None
        try:
            with C.quiet():
                r = rec.reconstruct_voxel(
                    vctx.vertices, vctx.voxel.phase, rng=C.run_seed(f["vpos"], 0)
                )
        finally:
            SH.active = False
        sp.fused_cost = False
        stored = np.load(cache / f"e0{tag}" / f"v{f['vidx']}_{a.variant}.npz")
        g, loc, _ = rec.last_eval_counts
        reproduced = bool(
            np.array_equal(np.asarray(r.orientation, np.float64), stored["s0_R_final"])
            and int(g) == int(stored["s0_evals_global"])
            and int(loc) == int(stored["s0_evals_local"])
        )
        ent = dict(
            f,
            change=(
                ("fixed" if f["wrong_off"] else "broken")
                if f["wrong_off"] != f["wrong_on"]
                else ("both wrong" if f["wrong_off"] else "both right")
            ),
            reproduced=reproduced,
            shadowed=SH.n,
            differing=SH.n_diff,
            first=SH.first,
            seconds=time.perf_counter() - t0,
        )
        res.append(ent)
        fs = SH.first
        print(
            f"{f['pair']:16s} v{f['vidx']:<6d} {ent['change']:7s} reproduced {reproduced}  "
            f"shadowed {SH.n} differing {SH.n_diff}  first "
            + (
                f"call {fs['call']} {fs['path']} r={fs['pixel_radius']} {fs['fields']} "
                f"rule {fs['within_rule']} {fs['max_near_px']}"
                if fs
                else "none"
            ),
            flush=True,
        )
    out = dict(
        meta=dict(
            head=S1.git_head(),
            dirty=S1.git_dirty(),
            variant=a.variant,
            arms=ARMS,
            changed=a.changed,
        ),
        flips=res,
        totals=dict(
            n=len(res),
            reproduced=sum(r["reproduced"] for r in res),
            first_within_rule=sum(bool(r["first"] and r["first"]["within_rule"]) for r in res),
            no_difference=sum(r["first"] is None for r in res),
        ),
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(json.dumps(out["totals"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
