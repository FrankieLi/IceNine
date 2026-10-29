"""S2 Step 3: the batched fused entry (``_rasterize.fused_evaluate_many``, used by
``BatchedCostEvaluator`` when the key FusedCost is on) against the single-orientation fused path.

Two subcommands, on S1's datasets and orientation set (``s1_parity.orientation_set``):

* ``dump-fused --out X.npz``: every ``OverlapInfo`` field of the single-orientation fused
  ``evaluate`` (``VoxelCostFunction(fused_cost=True)``) for every covered mode
  (``s2_parity.MODES``). Run on the base commit and on the new code, then ``compare`` (from
  ``s1_parity.py``) shows that moving the per-orientation code of ``_fused_eval.c`` into shared
  functions left the single-orientation path bit for bit unchanged.
* ``batched --ref X.npz``: ``BatchedCostEvaluator(fused_cost=True).evaluate_many`` in chunks of
  each ``--batch`` size, for the modes the batched cost supports (pixel radius 0), against the
  single-orientation fused results of ``--ref`` (normally the dump of the new code), every field
  with ``==``; also the evaluation count and that every orientation went through the fused entry.

  uv run python scripts/speed_parity/s2/s2_batched_check.py dump-fused \
      --out scripts/speed_parity/cache/s2_batched/fused_base.npz
  uv run python scripts/speed_parity/s2/s2_batched_check.py batched \
      --ref scripts/speed_parity/cache/s2_batched/fused_new.npz --json benchmarks/s2_batched/x.json

Single-threaded (``torch.set_num_threads(1)``), as ``s1_parity.py``.
"""

import argparse
import contextlib
import io
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
sys.path.insert(0, str(ICE))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))

import s1_parity as S1  # noqa: E402
import s2_parity as SP  # noqa: E402

DATASETS = ["three0"] + [f"full{v}" for v in S1.FULL_VOXELS]
BATCHED_MODES = [m for m, kw in SP.MODES.items() if kw.get("pixel_radius", 0) == 0]


def _fields(infos: List[Any]) -> Dict[str, np.ndarray]:
    return {f: np.array([getattr(o, f) for o in infos]) for f in S1.FIELDS}


def cmd_dump_fused(a: argparse.Namespace) -> int:
    out: Dict[str, np.ndarray] = {}
    for name in a.datasets:
        with contextlib.redirect_stdout(io.StringIO()):
            ds = S1.load_dataset(name)
        R = S1.orientation_set(ds, a.n_near, a.n_fz)
        out[f"{name}/R"] = R
        for mode in a.modes:
            cf = ds.cost_fn(dict(SP.MODES[mode], fused_cost=True))
            t0 = time.perf_counter()
            with contextlib.redirect_stdout(io.StringIO()):
                infos = [cf.evaluate(r, ds.verts, ds.phase) for r in R]
            dt = time.perf_counter() - t0
            if not cf._fused_ctx or any(c.reason for c in cf._fused_ctx.values()):
                raise RuntimeError(f"{name}/{mode}: the fused path was not used")
            for f, v in _fields(infos).items():
                out[f"{name}/{mode}/{f}"] = v
            print(f"{name:10s} {mode:15s} n={len(R)}  {dt * 1e3:.0f} ms", flush=True)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(a.out, **out)  # type: ignore[arg-type]
    print("wrote", a.out)
    return 0


def cmd_batched(a: argparse.Namespace) -> int:
    from icenine.batched_cost import BatchedCostEvaluator

    ref = np.load(a.ref, allow_pickle=False)
    rows: List[Dict[str, Any]] = []
    bad = 0
    print(f"{'dataset':10s} {'mode':14s} {'B':>5s} {'n':>5s}  differing orientations per field")
    for name in a.datasets:
        with contextlib.redirect_stdout(io.StringIO()):
            ds = S1.load_dataset(name)
        R = ref[f"{name}/R"]
        for mode in a.modes:
            for B in a.batch:
                cf = ds.cost_fn(dict(SP.MODES[mode], fused_cost=True))
                ev = BatchedCostEvaluator(cf, "cpu")
                got: List[Any] = []
                with contextlib.redirect_stdout(io.StringIO()):
                    for s0 in range(0, len(R), B):
                        got += ev.evaluate_many(R[s0 : s0 + B], ds.verts, ds.phase)
                diff = {}
                for f, v in _fields(got).items():
                    d = int((v != ref[f"{name}/{mode}/{f}"]).sum())
                    if d:
                        diff[f] = d
                n_fused = ev.counts["orient_fused"]
                ok = not diff and cf.eval_count == len(R) and n_fused == len(R)
                bad += not ok
                rows.append(
                    dict(
                        dataset=name,
                        mode=mode,
                        batch=B,
                        n=len(R),
                        differing=diff,
                        eval_count=cf.eval_count,
                        orient_fused=n_fused,
                        identical=ok,
                    )
                )
                msg = "identical" if ok else f"DIFFER {diff} evals {cf.eval_count} fused {n_fused}"
                print(f"{name:10s} {mode:14s} {B:5d} {len(R):5d}  {msg}", flush=True)
    verdict = "batched fused equals single fused (every field)" if bad == 0 else "DIFFERENCES"
    print("RESULT:", verdict)
    if a.json:
        Path(a.json).parent.mkdir(parents=True, exist_ok=True)
        Path(a.json).write_text(
            json.dumps(dict(ref=Path(a.ref).name, rows=rows, failing=bad), indent=1) + "\n"
        )
    return 0 if bad == 0 else 1


def main() -> int:
    torch.set_num_threads(1)
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("dump-fused")
    d.add_argument("--out", required=True)
    d.add_argument("--datasets", nargs="+", default=DATASETS)
    d.add_argument("--modes", nargs="+", default=list(SP.MODES), choices=list(SP.MODES))
    d.add_argument("--n-near", type=int, default=500)
    d.add_argument("--n-fz", type=int, default=500)
    b = sub.add_parser("batched")
    b.add_argument("--ref", required=True)
    b.add_argument("--datasets", nargs="+", default=DATASETS)
    b.add_argument("--modes", nargs="+", default=BATCHED_MODES, choices=BATCHED_MODES)
    b.add_argument("--batch", nargs="+", type=int, default=[1, 7, 64, 1000])
    b.add_argument("--json", default=None)
    a = ap.parse_args()
    for k in ("out", "ref", "json"):  # the dataset loaders change the working directory
        if getattr(a, k, None):
            setattr(a, k, str(Path(getattr(a, k)).resolve()))
    return cmd_dump_fused(a) if a.cmd == "dump-fused" else cmd_batched(a)


if __name__ == "__main__":
    sys.exit(main())
