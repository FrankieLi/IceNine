"""AR3 step 6, part 2: whole-layer BFS runs on the gradient sample with the picked seed remedy.

  run_seed.py --kind S     # D4's setting (indexed seeds, the step 4 fit wrapper) + picked remedy
  run_seed.py --kind K2    # run K's setting (the library's cma keys, no fit wrapper) + remedy

No library change (``seed_pick.SeedPick`` replaces two instance methods). Writes
``seed/runs/<kind>/``: ``recon.mic.gz``, ``records.json.gz``, ``run.json`` (stats, counters of the
remedy), ``seed_log.json.gz`` (one row per seed step), ``progress.json`` (every 5 minutes),
``preflight.json``, ``run.config``. The boundary pass is made after by ``seed_pass.py``.
"""

import argparse
import json
import time
from pathlib import Path
from typing import Optional

import ar3_common as C
import gradient_common as G
import numpy as np
import run_keys as RK
import run_wrapped as RW
import seed_common as S
from seed_pick import SeedPick

RUNS = S.SEED / "runs"


def run(kind: str, out: Path, max_voxels: Optional[int]) -> None:
    from icenine import ray_indexing as RI

    label = C.preflight_label(out)
    extra = RK.KEYS if kind == "K2" else ""
    setup, bfs, t_setup = C.build(C.CONFIG, out, sample=G.GRID_MIC, stack_name=G.STACK, extra=extra)
    mic = setup.sample.get_mic()
    p = setup.search_params
    assert not p.paint_grains and p.seed_by_indexing and p.fused_cost
    w = None
    if kind in ("S", "S5"):
        assert p.local_optimizer != "cma"
        w = RW.Wrap(bfs)
        orig_lo = w.orig
    else:
        assert p.local_optimizer == "cma" and p.cma_sigma0_deg == 0.05
        orig_lo = bfs.reconstructor.local_optimization
    sym = RI.proper_rotations(setup.exp_setup.get_sample_symmetry())
    pick = SeedPick(bfs, mic, orig_lo, sym, use_tries=kind != "S5")
    pick.install()
    rng = np.random.default_rng(0)
    t0 = time.time()
    RW.progress_thread(out, w, bfs, t0) if w is not None else RK.progress(out, bfs, t0)
    order = bfs.reconstruct_sample(
        output_mic=str(out / "recon.mic"), max_voxels=max_voxels, rng=rng
    )
    wall = time.time() - t0
    (out / "records.json").write_text(
        json.dumps({str(k): vars(v) for k, v in bfs.records.items()}) + "\n"
    )
    RW.gz_file(out / "records.json")
    RW.gz_file(out / "recon.mic")
    G.dump_gz(pick.log, out / "seed_log.json.gz")
    st = bfs.stats
    by = {k: sum(1 for e in pick.log if e["by"] == k) for k in ("tries", "seed_step")}
    info = dict(
        name=kind,
        n_voxels=len(mic.voxels),
        n_processed=len(order),
        wall_total_s=wall,
        setup_s=t_setup,
        timing_label=label,
        stats=st,
        wrapper=None if w is None else w.counts(),
        keys=extra.split("\n") if extra else [],
        pick=dict(
            seed_steps=len(pick.log),
            settled_by_tries=by["tries"],
            went_to_seed_step=by["seed_step"],
            tries_evals=int(sum(e.get("tries_evals", 0) for e in pick.log)),
            evals=int(sum(e["evals"] for e in pick.log)),
        ),
        evals_total=int(
            sum(st[k] for k in ("n_evals_seed", "n_evals_neighbor", "n_evals_revisit"))
            + st.get("n_evals_reseed", 0)
            + st.get("n_evals_restart", 0)
        ),
        index_pairing_s=float(st.get("index_pairing_s", 0.0)),
    )
    (out / "run.json").write_text(json.dumps(info, indent=2, default=float) + "\n")
    print(kind, "done", wall, "s", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kind", required=True, choices=["S", "K2", "S5"])
    ap.add_argument("--smoke", type=int, default=0, help="pilot: first N voxels, own folder")
    a = ap.parse_args()
    out = RUNS / (f"smoke_{a.kind}" if a.smoke else a.kind)
    out.mkdir(parents=True, exist_ok=True)
    run(a.kind, out, a.smoke or None)


if __name__ == "__main__":
    main()
