"""AR3 step 2: one run of the key-on BFS or of the key PaintGrains, on the pilot region or on the
whole layer (clean images, Q-max 8, FusedCost on, one worker, voxel order of ``default_rng(0)``).

  run_layer.py --name A --kind paint --sample pilot     # run A
  run_layer.py --name B --kind bfs   --sample pilot --fit-log   # run B
  run_layer.py --name C --kind paint --sample layer     # run C
  run_layer.py --name D --kind bfs   --sample layer --fit-log   # run D

Writes ``<OUT>/layer/runs/<name>/``: ``recon.mic``, ``records.json``, ``run.json``,
``preflight.json``; for ``paint``: ``paint.json`` (entries, stages, evaluations, scores per grain,
the loop's log) and ``seed_records.json`` (the seed path's result per seed: strict costs, finisher
costs and evaluations, the candidate taken); for ``bfs --fit-log``: ``fit_log.json`` (one row per
local fit of the BFS).

The fit log wraps ``BFSReconstruction._local_fit``. The start's strict cost and hit ratio come from
a second cost function (one extra evaluation per fit, counted apart; it draws nothing from the
run's generator).
"""

import argparse
import json
import time
from pathlib import Path
from typing import Any, Dict, List

import ar3_common as C
import numpy as np

LAYER_OUT = C.OUT / "layer"
LAYER_SAMPLE = C.EX / "SimInput" / "rand_500grains_1mm_neworient_s0_grid.mic"


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


def install_fit_log(bfs: Any, mic: Any) -> List[Dict[str, Any]]:
    """Wrap ``bfs._local_fit`` to log every local fit (see the module docstring)."""
    from icenine.ray_indexing import misorientation_deg, proper_rotations

    sym = proper_rotations(bfs.setup.exp_setup.get_sample_symmetry())
    cost_fn2 = bfs.reconstructor._make_local_cost_fn()
    rows: List[Dict[str, Any]] = []
    orig = bfs._local_fit

    def logged(
        voxel: Any, vertices: Any, rng: Any, rec: Any, phase: str, sigma0_deg: Any = None
    ) -> Any:
        start = np.array(voxel.orientation, dtype=np.float64)
        n0 = rec.n_evals
        out = orig(voxel, vertices, rng, rec, phase, sigma0_deg)
        res, hit, _conf = out
        end = np.asarray(res.orientation, dtype=np.float64)
        info = cost_fn2.evaluate(start, vertices, voxel.phase)
        s_hit, _ = bfs._ratios(info)
        rows.append(
            dict(
                voxel=int(rec.voxel_idx),
                phase=phase,
                angle_deg=float(misorientation_deg(start[None], end[None], sym)[0]),
                exact_same=bool(np.array_equal(start, end)),
                end_cost=float(res.cost),
                end_hit=float(hit),
                evals=int(rec.n_evals - n0),
                start_cost=float(info.cost),
                start_hit=float(s_hit),
                start=start.tolist(),
                end=end.tolist(),
            )
        )
        return out

    bfs._local_fit = logged
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--kind", required=True, choices=["paint", "bfs"])
    ap.add_argument("--sample", required=True, choices=["pilot", "layer"])
    ap.add_argument("--fit-log", action="store_true")
    ap.add_argument("--max-voxels", type=int, default=None)
    ap.add_argument("--config", default=str(C.CONFIG))
    ap.add_argument("--no-timing-gate", action="store_true", help="skip the preflight (tests)")
    a = ap.parse_args()
    out = LAYER_OUT / "runs" / a.name
    out.mkdir(parents=True, exist_ok=True)
    label = "not a timing run" if a.no_timing_gate else C.preflight_label(out)
    sample = LAYER_SAMPLE if a.sample == "layer" else None
    extra = "PaintGrains 1" if a.kind == "paint" else ""
    setup, bfs, t_setup = C.build(Path(a.config), out, sample=sample, extra=extra)
    mic = setup.sample.get_mic()
    n = len(mic.voxels)
    p = setup.search_params
    assert bool(p.paint_grains) == (a.kind == "paint") and p.seed_by_indexing and p.fused_cost
    fit_rows = install_fit_log(bfs, mic) if a.fit_log and a.kind == "bfs" else None

    rng = np.random.default_rng(0)
    t0 = time.time()
    order = bfs.reconstruct_sample(
        output_mic=str(out / "recon.mic"), max_voxels=a.max_voxels, rng=rng
    )
    wall = time.time() - t0
    (out / "records.json").write_text(
        json.dumps({str(k): vars(v) for k, v in bfs.records.items()}) + "\n"
    )
    info: Dict[str, Any] = dict(
        name=a.name,
        kind=a.kind,
        sample=a.sample,
        n_voxels=n,
        n_processed=len(order),
        wall_total_s=wall,
        setup_s=t_setup,
        timing_label=label,
        fused_cost=bool(p.fused_cost),
        seed_by_indexing=bool(p.seed_by_indexing),
        paint_grains=bool(p.paint_grains),
        torch_threads=C.torch.get_num_threads(),
        stats=bfs.stats,
    )
    if a.kind == "paint":
        res, hooks = bfs.paint_result, bfs.paint_hooks
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
        (out / "seed_records.json").write_text(json.dumps(hooks.seed_records) + "\n")
        (out / "fit_records.json").write_text(json.dumps(hooks.fit_records) + "\n")
        info.update(
            n_entries=len(res.entries),
            evals=res.evals,
            evals_total=int(sum(res.evals.values())),
            counts=res.counts,
            index_pairing_s=float(bfs.stats.get("index_pairing_s", 0.0)),
            cost_fn_eval_count=int(hooks.cost_fn.eval_count),
        )
    else:
        info["evals_total"] = int(
            sum(bfs.stats[k] for k in ("n_evals_seed", "n_evals_neighbor", "n_evals_revisit"))
            + bfs.stats.get("n_evals_reseed", 0)
            + bfs.stats.get("n_evals_restart", 0)
        )
        info["index_pairing_s"] = float(bfs.stats.get("index_pairing_s", 0.0))
        if getattr(bfs, "index_log", None):
            (out / "index_log.json").write_text(json.dumps(bfs.index_log) + "\n")
    if fit_rows is not None:
        (out / "fit_log.json").write_text(json.dumps(fit_rows) + "\n")
        info["n_fits_logged"] = len(fit_rows)
    (out / "run.json").write_text(json.dumps(info, indent=2, default=float) + "\n")
    print(json.dumps({k: v for k, v in info.items() if k != "stats"}, indent=2, default=float))


if __name__ == "__main__":
    main()
