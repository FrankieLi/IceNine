#!/usr/bin/env python3
"""Cost of the batched F-lowQ Q5 pass against the per-candidate pass, on ONE worker (all thread
counts 1, quiet machine: `preflight.require_quiet()`, JSON saved next to the results).

Per case (25 random cases, seed 0, as `timing.py`) the same 200 harvested candidates (source 0)
are scored by
  * the per-candidate pass `features(R, ...)` (reference), the Q_max-8 local evaluation (the unit),
  * `features_batch` in batches of 1, 50 and 200.
Per-candidate seconds = total / 200, median of REPS repetitions, then the median over cases. All
methods run interleaved inside every repetition. Equality of the batched and the per-candidate
features is checked on the same candidates (exact).
Output (benchmarks/coarse_proxy/): timing_batch.json, timing_batch_preflight.json,
timing_batch_tables.md

  uv run python scripts/coarse_proxy/timing_batch.py
"""

import json
import sys
import time
from pathlib import Path
from typing import Any, Callable, Dict, List

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import base as B  # noqa: E402

C = B.C
import features as F  # noqa: E402

sys.path.insert(0, str(B.ICENINE_PY / "scripts" / "common"))
import doc_tables  # noqa: E402
import preflight  # noqa: E402

N_CASES, N_CAND, REPS = 25, 200, 3
SIZES = (1, 50, 200)


def timed(fn: Callable[[], Any]) -> float:
    t0 = time.perf_counter()
    fn()
    return time.perf_counter() - t0


def main() -> None:
    import torch
    from optimizer_sweep import voxel_context

    torch.set_num_threads(1)
    B.OUT.mkdir(parents=True, exist_ok=True)
    info = preflight.require_quiet()
    (B.OUT / "timing_batch_preflight.json").write_text(json.dumps(info, indent=1, default=str))
    C.init_worker(C.worker_args())
    W = C.get_worker()
    D = B.load_e2_aligned()
    rng = np.random.default_rng(0)
    ts = B.tasks()
    cases = rng.choice(len(ts), N_CASES, replace=False)
    sizes = [len(np.load(C.CACHE_DIR / "e2" / f"v{v}_{var}.npz")["err"]) for v, _, var in ts]
    off = np.concatenate([[0], np.cumsum(sizes)])
    names = ["local_q8", "ref", "ref_geom", *[f"batch{n}" for n in SIZES], "batch200_geom"]
    per_case: Dict[str, List[float]] = {k: [] for k in names}
    n_equal, n_total = 0, 0
    for ci in cases:
        v, vpos, var = ts[ci]
        keys = np.load(C.CACHE_DIR / "images" / f"v{v}_{var}.npz")["keys"]
        C.attach(keys)
        vctx = voxel_context(W.ctx, v)
        vert, ph = vctx.vertices, vctx.voxel.phase
        fe = F.FeatureExtractor(W.local_fn, W.ctx.geo, ph, q_max=5.0)
        fe.set_image(keys)
        idx = np.nonzero(D["source"][off[ci] : off[ci + 1]] == 0)[0]
        R = D["R"][off[ci] + rng.choice(idx, N_CAND, replace=False)]
        Rf = R.astype(np.float32)
        ref = np.stack([fe.features(r, vert, ph) for r in R])  # also the warm-up
        got = fe.features_batch(R, vert, ph)
        n_equal += int((ref == got).all(axis=1).sum())
        n_total += len(R)
        reps: Dict[str, List[float]] = {k: [] for k in names}
        for _ in range(REPS):
            reps["local_q8"].append(timed(lambda: [W.local_fn.evaluate(f, vert, ph) for f in Rf]))
            reps["ref"].append(timed(lambda: [fe.features(r, vert, ph) for r in R]))
            reps["ref_geom"].append(timed(lambda: [fe.features(r, vert, ph, False) for r in R]))
            for n in SIZES:
                reps[f"batch{n}"].append(
                    timed(
                        lambda: [
                            fe.features_batch(R[i : i + n], vert, ph) for i in range(0, N_CAND, n)
                        ]
                    )
                )
            reps["batch200_geom"].append(timed(lambda: fe.features_batch(R, vert, ph, False)))
        for k in names:
            per_case[k].append(float(np.median(reps[k])) / N_CAND)
    med = {k: float(np.median(v)) for k, v in per_case.items()}
    unit = med["local_q8"]
    gbt = json.loads((B.OUT / "timing.json").read_text())["seconds_median"][
        "model_predict_per_candidate"
    ]
    proxy = {k: med[k] + gbt for k in ("ref", *[f"batch{n}" for n in SIZES])}
    out = dict(
        label="single-worker", n_cases=N_CASES, n_candidates_per_case=N_CAND, reps=REPS,
        seconds_median=med, equivalents_of_local_q8={k: v / unit for k, v in med.items()},
        gbt_predict_per_candidate=gbt,
        proxy_equivalents={k: v / unit for k, v in proxy.items()},
        identical_rows=n_equal, total_rows=n_total,
    )  # fmt: skip
    C.save_json(B.OUT / "timing_batch.json", out)
    rows = [
        dict(
            path="per-candidate (reference)",
            batch="1",
            ms=med["ref"] * 1e3,
            speedup=1.0,
            equiv=med["ref"] / unit,
            proxy=proxy["ref"] / unit,
        )
    ]
    for n in SIZES:
        k = f"batch{n}"
        rows.append(
            dict(
                path="batched",
                batch=str(n),
                ms=med[k] * 1e3,
                speedup=med["ref"] / med[k],
                equiv=med[k] / unit,
                proxy=proxy[k] / unit,
            )
        )
    cols = ["path", "batch", "ms", "speedup", "equiv", "proxy"]
    tab = doc_tables.markdown_table(
        rows, cols, formats=dict(ms=".3f", speedup=".1f", equiv=".2f", proxy=".2f")
    )
    doc_tables.write_tables(B.OUT / "timing_batch_tables.md", {"t3_lowq_batch_timing": tab})
    print(f"unit local Q8 evaluation {unit * 1e3:.3f} ms; GBT {gbt * 1e3:.4f} ms/candidate")
    print(f"identical rows {n_equal}/{n_total}")
    print(tab)
    print({k: round(v * 1e3, 4) for k, v in med.items()})


if __name__ == "__main__":
    main()
