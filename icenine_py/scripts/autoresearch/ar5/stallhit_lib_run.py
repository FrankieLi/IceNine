"""AR5 step 12, phase B: a library run with the key ``SeedIndexMCStallHit 3000`` on (``LibH``).

  RESEED_SAMPLE=mid uv run python scripts/autoresearch/ar5/stallhit_lib_run.py [--cap-min 240]

Everything is the library's own: the keys are written to the config (``ar3_common.build``
appends them, as ``rescue_lib_run.py`` of step 9 does); no method is replaced and nothing is
counted by patching. ``LibH`` is step 9's ``LibG`` (the keys of run ``S11L`` of step 5 plus
``SeedRescueCMA 0.25 1000``) plus ``SeedIndexMCStallHit 3000``: the library form of the script arm
``Hl`` (``stallhit_hl.py``). ``default_rng(0)``, clean images, Q-max 8, one worker. As in
``rescue_lib_run.py``: a capture of the boundary pass's result (and of the counters before it), a
progress file every 5 minutes and a watchdog at ``--cap-min`` minutes. The map before the pass is
rebuilt from ``bfs.boundary_log``.

Writes ``<cache>/LibH_<sample>/`` (cache = ``STALLHIT_CACHE``, default
``benchmarks/autoresearch/ar5/stallhit/cache``, untracked): ``recon.mic.gz`` (after the pass),
``recon_before_pass.mic.gz``, ``records.json.gz``, ``index_log.json.gz``, ``pass_tries.npz``,
``run.json``, ``preflight.json``, ``run.config``, ``progress.json``; ``run.done`` at the end.
Refuses to start if the folder exists. Wall times are contended, not timing results.
"""

import argparse
import json
import os
import sys
import threading
import time
from pathlib import Path
from typing import Any, Dict

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "ar3"))
sys.path.insert(0, str(HERE))

import reseed_common as RC  # noqa: E402  (first: RESEED_SAMPLE -> AR3_SAMPLE)

import ar3_common as C  # noqa: E402
import gradient_common as G  # noqa: E402
import numpy as np  # noqa: E402
import rescue_lib_run as RLR  # noqa: E402

KEYS = RLR.S11L + ["SeedRescueCMA 0.25 1000", "SeedIndexMCStallHit 3000"]
NAME = "LibH"
CACHE = Path(
    os.environ.get(
        "STALLHIT_CACHE",
        str(C.ICE / "benchmarks" / "autoresearch" / "ar5" / "stallhit" / "cache"),
    )
)


def run(out: Path, cap_min: float) -> None:
    from icenine import boundary_pass as LBP

    label = C.preflight_label(out)
    setup, bfs, t_setup = C.build(
        C.CONFIG, out, sample=G.GRID_MIC, stack_name=G.STACK, extra="\n".join(KEYS)
    )
    mic = setup.sample.get_mic()
    p = setup.search_params
    assert p.seed_by_indexing and p.fused_cost and p.local_optimizer == "mc"
    assert p.bfs_neighbor_cma_first and p.bfs_neighbor_cma_no_fallback and p.boundary_pass
    assert (p.seed_index_n_cost, p.seed_index_n_refine) == (20, 5) and p.seed_index_stop_hit == 0.95
    assert p.cma_sigma0_deg == 0.05 and p.cma_neighbor_max_evals == 60
    assert p.seed_index_cma_first == 200
    assert (p.seed_rescue_sigma_deg, p.seed_rescue_evals) == (0.25, 1000)
    assert p.seed_index_mc_stall_hit == 3000
    captured: Dict[str, Any] = {}
    orig_pass = LBP.run_boundary_pass

    def capture(bfs_: Any, mic_: Any) -> Any:
        captured["stats_before"] = {k: v for k, v in bfs_.stats.items() if k != "reseeds"}
        res = orig_pass(bfs_, mic_)
        captured["result"] = res
        return res

    LBP.run_boundary_pass = capture
    t0 = time.time()

    def progress() -> None:
        while True:
            time.sleep(300.0)
            try:
                (out / "progress.json").write_text(
                    json.dumps(dict(elapsed_s=time.time() - t0, stats=RLR.numeric(bfs.stats)))
                    + "\n"
                )
            except Exception:  # a progress write must never stop the run
                pass

    def watchdog() -> None:
        time.sleep(cap_min * 60.0)
        (out / "capped.json").write_text(
            json.dumps(
                dict(cap_min=cap_min, n_records=len(bfs.records), stats=RLR.numeric(bfs.stats))
            )
            + "\n"
        )
        print(f"CAP of {cap_min} min reached; stopped", flush=True)
        os._exit(3)

    threading.Thread(target=progress, daemon=True).start()
    threading.Thread(target=watchdog, daemon=True).start()
    order = bfs.reconstruct_sample(output_mic=str(out / "recon.mic"), rng=np.random.default_rng(0))
    wall = time.time() - t0
    res = captured["result"]
    np.savez_compressed(
        out / "pass_tries.npz",
        voxel=np.array([t["voxel"] for t in res.tries], np.int32),
        cand=np.array([t["cand"] for t in res.tries], np.int32),
        end=np.array([t["end"] for t in res.tries], np.float64).reshape(-1, 3, 3),
        cost=np.array([t["cost"] for t in res.tries]),
        hit=np.array([t["hit"] for t in res.tries]),
        evals=np.array([t["evals"] for t in res.tries], np.int32),
        boundary=np.array(res.boundary_voxels, np.int32),
        changed=np.array(res.changed, np.int32),
    )
    (out / "records.json").write_text(
        json.dumps({str(k): vars(v) for k, v in bfs.records.items()}) + "\n"
    )
    RLR.gz_file(out / "records.json")
    (out / "index_log.json").write_text(json.dumps(bfs.index_log, default=float) + "\n")
    RLR.gz_file(out / "index_log.json")
    for e in bfs.boundary_log:  # rebuild the map before the pass from the pass's log
        v = mic.voxels[e["voxel"]]
        v.orientation = e["old"].copy()
        v.cost, v.confidence, v.overlap_ratio = (
            e["old_cost"],
            e["old_confidence"],
            e["old_overlap_ratio"],
        )
    mic.write(str(out / "recon_before_pass.mic"))
    RLR.gz_file(out / "recon_before_pass.mic")
    RLR.gz_file(out / "recon.mic")
    st = bfs.stats
    info = dict(
        name=NAME,
        sample=RC.SAMPLE,
        keys=KEYS,
        n_voxels=len(mic.voxels),
        n_processed=len(order),
        wall_total_s=wall,
        setup_s=t_setup,
        timing_label=label,
        wall_label="contended (other jobs on the machine); not a timing result",
        stats={k: v for k, v in st.items() if k != "reseeds"},
        stats_before_pass=RLR.numeric(captured["stats_before"]),
        evals_total_before_pass=RLR.evals_total(captured["stats_before"]),
        evals_total_with_pass=RLR.evals_total(st) + int(st.get("n_evals_boundary_pass", 0)),
        n_boundary_log=len(bfs.boundary_log),
    )
    (out / "run.json").write_text(json.dumps(info, indent=2, default=float) + "\n")
    (out / "run.done").write_text("done\n")
    print(NAME, RC.SAMPLE, "done", wall, "s", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cap-min", type=float, default=240.0)
    a = ap.parse_args()
    assert RC.SAMPLE in ("mid", "grad", "big"), "held and held2 are not run in phase B"
    out = CACHE / f"{NAME}_{RC.SAMPLE}"
    if out.exists():
        raise SystemExit(f"refusing to start: {out} exists")
    out.mkdir(parents=True)
    run(out, a.cap_min)
    os._exit(0)


if __name__ == "__main__":
    main()
