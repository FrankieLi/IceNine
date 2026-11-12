"""Item 14: one full search (global search, no start) on full-sample voxels, arms F and FS, timed
by section: batched screen, smooth ranking key, strict evaluations, rest.

Uses the section timers of ``scripts/speed_parity/s2/s2_search.py`` (imported, not edited) plus one
timer on ``SmoothRanker.__call__``. Setup as the Task E ``full`` set (clean full-sample images,
Q-max 8 configuration, ``min_sin_eta`` 0, rng seed 1). Arm F = ``FusedCost`` on, ``SmoothCostRank``
off; arm FS = both on. Single worker; the caller gates the preflight (``run_item14.sh``); the label
is written to the result.

  uv run python scripts/autoresearch/item14/split.py --voxels 15901,17034 --out DIR [--label L]

Shares: ``screen`` = ScreenEvaluator and pixel-radius screen evaluations, ``smooth`` = the ranker's
calls, ``strict`` = strict ``evaluate`` calls plus the batched local cost (neither of the first
two), ``rest`` = wall minus these three.
"""

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
sys.path.insert(0, str(ICE))
sys.path.insert(0, str(ICE / "scripts" / "common"))
sys.path.insert(0, str(ICE / "scripts" / "eval_reduction"))
sys.path.insert(0, str(ICE / "scripts" / "speed_parity"))
sys.path.insert(0, str(ICE / "scripts" / "speed_parity" / "s2"))

SPEC_S = "E3+E4,tau=0.25"
SMOOTH: Dict[str, float] = {"s": 0.0, "n": 0.0}


def shares(wall: float, screen: float, smooth: float, strict: float) -> Dict[str, float]:
    """Fractions of the wall time: screen, smooth, strict and the rest (they sum to 1)."""
    return dict(
        screen=screen / wall,
        smooth=smooth / wall,
        strict=strict / wall,
        rest=(wall - screen - smooth - strict) / wall,
    )


def _patch_ranker() -> None:
    from icenine import smooth_cost as SC

    orig = SC.SmoothRanker.__call__

    def timed(self: Any, level: int, candidates: Any, diameter_rad: float) -> Any:
        t0 = time.perf_counter()
        try:
            return orig(self, level, candidates, diameter_rad)
        finally:
            SMOOTH["s"] += time.perf_counter() - t0
            SMOOTH["n"] += len(candidates)

    SC.SmoothRanker.__call__ = timed  # type: ignore[method-assign]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--voxels", default="15901,17034")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--label", default="single-worker, preflight gated by the caller")
    ap.add_argument("--no-warmup", action="store_true")
    a = ap.parse_args()

    import numpy as np
    import s2_search as S2
    import parity_common as PC

    S2.install()
    _patch_ranker()
    ctx = PC.PyContext()
    rec, sp = ctx.rec, ctx.rec.setup.search_params
    vox = [int(v) for v in a.voxels.split(",")]
    arms = {"F": "off", "FS": SPEC_S}

    def one(v: int, arm: str) -> Dict[str, Any]:
        sp.fused_cost = True
        sp.smooth_rank = arms[arm]
        SMOOTH.update(s=0.0, n=0.0)
        r = S2.timed_run(rec, ctx.verts(v), int(ctx.mic.voxels[v].phase),
                         np.random.default_rng(1), ctx.truth(v))  # fmt: skip
        wall = r["wall_s"]
        strict = r["evaluate_s"] + r["batched_quick_mc"]["s"]
        r.update(
            voxel=v, arm=arm, smooth_s=SMOOTH["s"], smooth_evals=int(SMOOTH["n"]),
            strict_s=strict, shares=shares(wall, r["screen_s"], SMOOTH["s"], strict),
        )  # fmt: skip
        r["fused_calls_all_evaluate"] = r["evaluate_fused_n"] == r["evaluate_n"]
        return r

    if not a.no_warmup:  # untimed: fills the lit-pixel trees and observed rays of the smooth key
        for arm in arms:
            one(vox[0], arm)
    runs: List[Dict[str, Any]] = []
    for v in vox:
        for arm in arms:
            r = one(v, arm)
            runs.append(r)
            s = r["shares"]
            print(
                f"v{v} {arm:2s} wall {r['wall_s']:6.2f}s screen {s['screen']:.2f} "
                f"smooth {s['smooth']:.2f} strict {s['strict']:.2f} rest {s['rest']:.2f}  "
                f"evals {r['evals']['total']} smooth_evals {r['smooth_evals']} "
                f"fused {r['evaluate_fused_n']}/{r['evaluate_n']}  err {r['err_deg']:.3f}",
                flush=True,
            )
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "split.json").write_text(
        json.dumps(dict(label=a.label, voxels=vox, spec=SPEC_S, runs=runs), indent=1, default=float)
        + "\n"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
