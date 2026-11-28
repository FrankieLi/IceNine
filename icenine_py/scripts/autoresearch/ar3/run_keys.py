"""AR3 step 5, part 2: whole-layer BFS runs with the library's own keys, no wrapper.

  run_keys.py --kind K     # indexed seeds on the gradient sample (compared with D4 and step 3's D)
  run_keys.py --kind KE    # the same on step 2's sample (compared with E4 and step 2's D)

Keys added to ``key_on_fused.config`` (every other key as in D4 and E4, which use the same file):

    LocalOptimizer cma        # the neighbour, revisit and retry fits, the seed's local
                              # optimization and the grid search's refinement all use CMA
    CMASigma0 0.05            # first step in degrees (the library default is 0.2)
    CMANeighborMaxEvals 60    # budget of a neighbour or revisit fit (start's evaluation included)

``CMARetrySigma0`` and ``CMAMaxEvals`` are left at the library defaults (1.5 degrees, retry on;
1000 evaluations). The retry is a second fit from the same start with a wider first step (same
60-evaluation budget) when the caller rejects the first. D4 and E4 had a different trigger: the
step 4 wrapper ran the stock variance-minimizing fit (box 5 degrees) when the CMA end cost was not
below the start cost. The seed's own fit: with ``SeedByIndexing`` the seed is indexed, the strict
cost is taken on 5 candidates, and the 2 with the lowest cost are refined by the local
optimization (``refine_from_candidates``); the grid search of a failed index seed refines its
candidates there too. In cma mode ``refine_from_candidates`` makes one CMA run (``CMASigma0``,
``CMAMaxEvals`` 1000) per candidate, up to ``MaxDiscreteCandidates`` (30 in this config), stops
when a run that lowers the best cost reaches hit ratio 1, and skips VarianceMinimizing. D4's seed
path used FindOptimal MC plus VarianceMinimizing. So K and KE change the seed path as well as the
neighbour fits; this is what the keys do and is not hidden.

Writes ``boundary/runs/<kind>/``: ``recon.mic.gz``, ``records.json.gz``, ``run.json``,
``progress.json`` (every 5 minutes), ``preflight.json``, ``run.config``.
"""

import argparse
import json
import threading
import time
from pathlib import Path
from typing import Optional

import ar3_common as C
import boundary_common as B
import gradient_common as G
import numpy as np
import run_layer as RL
import run_wrapped as RW

KEYS = "LocalOptimizer cma\nCMASigma0 0.05\nCMANeighborMaxEvals 60"


def progress(out: Path, bfs: object, t0: float, every: float = 300.0) -> None:
    def loop() -> None:
        while True:
            time.sleep(every)
            try:
                stats = bfs.stats  # type: ignore
                st = {k: v for k, v in stats.items() if isinstance(v, (int, float))}
                (out / "progress.json").write_text(
                    json.dumps(dict(elapsed_s=time.time() - t0, stats=st), default=float) + "\n"
                )
            except Exception:  # a progress write must never stop the run
                pass

    threading.Thread(target=loop, daemon=True).start()


def run(kind: str, out: Path, max_voxels: Optional[int]) -> None:
    label = C.preflight_label(out)
    if kind == "K":
        setup, bfs, t_setup = C.build(
            C.CONFIG, out, sample=G.GRID_MIC, stack_name=G.STACK, extra=KEYS
        )
    else:
        setup, bfs, t_setup = C.build(C.CONFIG, out, sample=RL.LAYER_SAMPLE, extra=KEYS)
    mic = setup.sample.get_mic()
    p = setup.search_params
    assert not p.paint_grains and p.seed_by_indexing and p.fused_cost
    assert (
        p.local_optimizer == "cma" and p.cma_sigma0_deg == 0.05 and p.cma_neighbor_max_evals == 60
    )
    rng = np.random.default_rng(0)
    t0 = time.time()
    progress(out, bfs, t0)
    order = bfs.reconstruct_sample(
        output_mic=str(out / "recon.mic"), max_voxels=max_voxels, rng=rng
    )
    wall = time.time() - t0
    (out / "records.json").write_text(
        json.dumps({str(k): vars(v) for k, v in bfs.records.items()}) + "\n"
    )
    RW.gz_file(out / "records.json")
    RW.gz_file(out / "recon.mic")
    st = bfs.stats
    info = dict(
        name=kind,
        n_voxels=len(mic.voxels),
        n_processed=len(order),
        wall_total_s=wall,
        setup_s=t_setup,
        timing_label=label,
        keys=KEYS.split("\n"),
        cma=dict(
            sigma0_deg=p.cma_sigma0_deg,
            neighbor_max_evals=p.cma_neighbor_max_evals,
            retry_sigma0_deg=p.cma_retry_sigma0_deg,
            seed_max_evals=p.cma_max_evals,
        ),
        stats=st,
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
    ap.add_argument("--kind", required=True, choices=["K", "KE"])
    ap.add_argument("--smoke", type=int, default=0, help="pilot: first N voxels, own folder")
    a = ap.parse_args()
    out = B.RUNS / (f"smoke_{a.kind}" if a.smoke else a.kind)
    out.mkdir(parents=True, exist_ok=True)
    run(a.kind, out, a.smoke or None)


if __name__ == "__main__":
    main()
