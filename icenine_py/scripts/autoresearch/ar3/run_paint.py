"""AR3 step 1: the paint run on the 2,000-voxel pilot region (clean images, Q-max 8), one worker.

  uv run python scripts/autoresearch/ar3/run_paint.py --name paint [--max-voxels N] [--no-side]

Order: the shuffled voxel order of the AR2 BFS runs (``default_rng(0).shuffle(range(n))``; no new
draw). The same ``default_rng(0)`` then serves the local fits and the grid-search fallback.
Writes ``<OUT>/runs/<name>/`` with ``recon.mic`` (the painted map; ``g2_bfs.errors`` reads it),
``records.json`` (source per voxel), ``paint.json`` (entries, stages, evaluations, scores per
grain), ``run.json`` and, with the side arm, ``side/`` (``recon.mic``, ``side.json``).
"""

import argparse
import json
import time
from pathlib import Path
from typing import Any, Dict, List

import ar3_common as C
import numpy as np


def entry_row(e: Any) -> Dict[str, Any]:
    return dict(
        voxel=int(e.voxel),
        kind=e.kind,
        cost=float(e.cost),
        hit=float(e.hit),
        best_hit=float(e.best_hit),
        near_earlier=e.near_earlier,
        near_deg=None if e.near_earlier is None else float(e.near_deg),
        flood_scores=int(e.flood_scores),
        flood_passes=int(e.flood_passes),
        R=np.asarray(e.R).tolist(),
    )


def side_arm(hooks: Any, res: Any, n: int, seed: int = 0) -> Dict[str, Any]:
    """On every contested voxel painted by a grain (two or more grains pass; fit-settled voxels
    are left out): one local fit from each of the two lowest-cost grains, the lower final cost
    kept. Evaluations are counted apart from the paint run's."""
    hooks.rng = np.random.default_rng(seed)
    rows: List[Dict[str, Any]] = []
    R_kept: Dict[int, np.ndarray] = {}
    R_fit = []
    evals = 0
    t0 = time.perf_counter()
    for v in range(n):
        if v in res.fits:
            continue
        r = res.ranked(v)
        if len(r) < 2:
            continue
        fits = []
        for g, cost, hit in r[:2]:
            R, c, h, ne = hooks.local_fit(v, res.entries[g].R)
            fits.append((g, R, c, h, ne))
            evals += ne
        k = int(np.argmin([f[2] for f in fits]))
        R_kept[v] = fits[k][1]
        R_fit.append((v, fits[0][1], fits[1][1]))
        rows.append(
            dict(
                voxel=int(v),
                grains=[int(f[0]) for f in fits],
                paint_costs=[float(r[0][1]), float(r[1][1])],
                margin=float(r[1][1] - r[0][1]),
                n_passing=len(r),
                fit_costs=[float(f[2]) for f in fits],
                fit_hits=[float(f[3]) for f in fits],
                fit_evals=[int(f[4]) for f in fits],
                kept=k,
            )
        )
    return dict(
        rows=rows,
        evals=int(evals),
        wall_s=time.perf_counter() - t0,
        R_kept={v: R for v, R in R_kept.items()},
        R_fits=R_fit,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--max-voxels", type=int, default=None, help="first N voxels of the order")
    ap.add_argument("--no-side", action="store_true")
    ap.add_argument("--config", default=str(C.CONFIG))
    ap.add_argument("--no-timing-gate", action="store_true", help="skip the preflight (tests)")
    a = ap.parse_args()
    out = C.OUT / "runs" / a.name
    out.mkdir(parents=True, exist_ok=True)
    label = "not a timing run" if a.no_timing_gate else C.preflight_label(out)
    cfg = Path(a.config)
    setup, bfs, t_setup = C.build(cfg, out)

    from icenine import grain_paint as GP

    mic = setup.sample.get_mic()
    n = len(mic.voxels)
    rng = np.random.default_rng(0)
    order = list(range(n))
    rng.shuffle(order)  # the BFS runs' order
    if a.max_voxels:
        order = order[: a.max_voxels]
    hooks = GP.ReconstructorHooks(bfs, mic, rng)
    t0 = time.time()
    res = GP.paint_grains(order, n, hooks)
    wall = time.time() - t0
    pairing = float(getattr(bfs.reconstructor, "index_pairing_s", 0.0))

    final: Dict[int, np.ndarray] = {}
    records: Dict[str, Any] = {}
    for v in range(n):
        R, src = res.final_orientation(v)
        if R is None:  # not reached (only with --max-voxels)
            R, src = np.asarray(mic.voxels[v].orientation), "not_visited"
        final[v] = R
        records[str(v)] = dict(source=src, n_evals=0, stage=res.stage.get(v, ""))
    C.write_mic(mic, final, out / "recon.mic")
    (out / "records.json").write_text(json.dumps(records) + "\n")
    paint = dict(
        entries=[entry_row(e) for e in res.entries],
        stage={str(k): v for k, v in res.stage.items()},
        fits={
            str(v): dict(R=f.R.tolist(), cost=f.cost, hit=f.hit, from_grain=f.from_grain)
            for v, f in res.fits.items()
        },
        failed={str(v): dict(cost=f.cost, hit=f.hit) for v, f in res.failed.items()},
        evals=res.evals,
        wall_s=res.wall_s,
        counts=res.counts,
        log=res.log,
        scores=[[[int(v), s[0], s[1], bool(s[2])] for v, s in sc.items()] for sc in res.scores],
    )
    (out / "paint.json").write_text(json.dumps(paint) + "\n")
    info: Dict[str, Any] = dict(
        name=a.name,
        config=cfg.name,
        n_voxels=n,
        n_order=len(order),
        n_entries=len(res.entries),
        evals=res.evals,
        evals_total=int(sum(res.evals.values())),
        counts=res.counts,
        wall_total_s=wall,
        wall_by_stage_s=res.wall_s,
        index_pairing_s=pairing,
        setup_s=t_setup,
        timing_label=label,
        fused_cost=bool(getattr(setup.search_params, "fused_cost", False)),
        torch_threads=C.torch.get_num_threads(),
        cost_fn_eval_count=int(hooks.cost_fn.eval_count),
    )
    if not a.no_side:
        side = side_arm(hooks, res, n)
        sd = out / "side"
        sd.mkdir(exist_ok=True)
        final2 = dict(final)
        final2.update(side["R_kept"])
        mic2 = setup.sample.get_mic()
        C.write_mic(mic2, final2, sd / "recon.mic")
        (sd / "records.json").write_text(json.dumps(records) + "\n")
        (sd / "side.json").write_text(
            json.dumps(dict(rows=side["rows"], evals=side["evals"], wall_s=side["wall_s"])) + "\n"
        )
        np.savez(
            sd / "side_fits.npz",
            voxel=np.array([x[0] for x in side["R_fits"]]),
            R1=np.array([x[1] for x in side["R_fits"]]),
            R2=np.array([x[2] for x in side["R_fits"]]),
        )
        info["side"] = dict(
            n_contested=len(side["rows"]), evals=side["evals"], wall_s=side["wall_s"]
        )
    (out / "run.json").write_text(json.dumps(info, indent=2, default=float) + "\n")
    print(json.dumps(info, indent=2, default=float), flush=True)


if __name__ == "__main__":
    main()
