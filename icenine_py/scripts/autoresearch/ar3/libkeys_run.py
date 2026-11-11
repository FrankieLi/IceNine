"""AR3 step 10: a validation run of the library keys on the sample ``big``.

  libkeys_run.py --kind D4L     # BFSNeighborCMAFirst + CMA settings + BoundaryPass
  libkeys_run.py --kind S5L     # + SeedIndexCandidates 20 5
  libkeys_run.py --kind S8L     # + SeedIndexStopHitRatio 0.95
  ... --smoke N                 # pilot: first N voxels (a BFS may still flood the whole layer)

Everything is the library's own: the keys are written to the config (``ar3_common.build`` appends
them to ``key_on_fused.config``); no method is replaced and nothing is counted by patching. The
only things added here are a capture of the pass's result (the function is called through its
module, so a wrapper sees its return value), a progress file every 5 minutes, and a watchdog at
twice the budget (60 minutes). The map before the pass is rebuilt from ``bfs.boundary_log`` and
written too, so it can be compared with the stored run's map.

Writes ``<CACHE>/<kind>/``: ``recon.mic.gz`` (after the pass), ``recon_before_pass.mic.gz``,
``records.json.gz``, ``run.json``, ``pass_tries.npz``, ``preflight.json``, ``run.config``,
``progress.json``; ``run.done`` at the end.
"""

import argparse
import gzip
import json
import os
import shutil
import threading
import time
from pathlib import Path
from typing import Any, Dict, Optional

import libkeys_common as K
import numpy as np

import ar3_common as C  # noqa: E402
import gradient_common as G  # noqa: E402

CAP_MIN = 60  # twice the 30 minute budget


def gz_file(path: Path) -> None:
    with open(path, "rb") as fi, open(str(path) + ".gz", "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as fo:
            shutil.copyfileobj(fi, fo)
    path.unlink()


def progress(out: Path, bfs: Any, t0: float, every: float = 300.0) -> None:
    def loop() -> None:
        while True:
            time.sleep(every)
            try:
                st = {k: v for k, v in bfs.stats.items() if isinstance(v, (int, float))}
                (out / "progress.json").write_text(
                    json.dumps(dict(elapsed_s=time.time() - t0, stats=st), default=float) + "\n"
                )
            except Exception:  # a progress write must never stop the run
                pass

    threading.Thread(target=loop, daemon=True).start()


def watchdog(out: Path, bfs: Any, cap_s: float, t0: float) -> None:
    def loop() -> None:
        time.sleep(cap_s)
        st = {k: v for k, v in bfs.stats.items() if isinstance(v, (int, float))}
        (out / "capped.json").write_text(
            json.dumps(
                dict(cap_s=cap_s, elapsed_s=time.time() - t0, n_records=len(bfs.records), stats=st),
                default=float,
            )
            + "\n"
        )
        print(f"CAP of {cap_s} s reached; stopped", flush=True)
        os._exit(3)

    threading.Thread(target=loop, daemon=True).start()


def run(kind: str, out: Path, max_voxels: Optional[int]) -> None:
    from icenine import boundary_pass as LBP

    label = C.preflight_label(out)
    extra = "\n".join(K.KEYS[kind])
    setup, bfs, t_setup = C.build(C.CONFIG, out, sample=G.GRID_MIC, stack_name=G.STACK, extra=extra)
    mic = setup.sample.get_mic()
    p = setup.search_params
    assert p.seed_by_indexing and p.fused_cost and p.local_optimizer == "mc" and not p.paint_grains
    assert p.bfs_neighbor_cma_first and p.boundary_pass
    assert (p.seed_index_n_cost, p.seed_index_n_refine) == ((5, 2) if kind == "D4L" else (20, 5))
    assert p.seed_index_stop_hit == (0.95 if kind == "S8L" else 0.0)
    assert p.cma_sigma0_deg == 0.05 and p.cma_neighbor_max_evals == 60
    captured: Dict[str, Any] = {}
    orig_pass = LBP.run_boundary_pass

    def capture(bfs_: Any, mic_: Any) -> Any:
        captured["stats_before"] = {k: v for k, v in bfs_.stats.items() if k != "reseeds"}
        res = orig_pass(bfs_, mic_)
        captured["result"] = res
        return res

    LBP.run_boundary_pass = capture
    rng = np.random.default_rng(0)
    t0 = time.time()
    progress(out, bfs, t0)
    watchdog(out, bfs, CAP_MIN * 60.0, t0)
    order = bfs.reconstruct_sample(
        output_mic=str(out / "recon.mic"), max_voxels=max_voxels, rng=rng
    )
    wall = time.time() - t0
    res = captured["result"]
    np.savez_compressed(
        out / "pass_tries.npz",
        voxel=np.array([t["voxel"] for t in res.tries], np.int32),
        cand=np.array([t["cand"] for t in res.tries], np.int32),
        end=np.array([t["end"] for t in res.tries], np.float32).reshape(-1, 3, 3),
        cost=np.array([t["cost"] for t in res.tries]),
        hit=np.array([t["hit"] for t in res.tries]),
        conf=np.array([t["conf"] for t in res.tries]),
        evals=np.array([t["evals"] for t in res.tries], np.int32),
        boundary=np.array(res.boundary_voxels, np.int32),
        changed=np.array(res.changed, np.int32),
        old_score=np.array([res.old_score[v] for v in res.boundary_voxels]).reshape(-1, 3),
    )
    (out / "records.json").write_text(
        json.dumps({str(k): vars(v) for k, v in bfs.records.items()}) + "\n"
    )
    gz_file(out / "records.json")
    # rebuild the map before the pass from the pass's log and write it
    for e in bfs.boundary_log:
        v = mic.voxels[e["voxel"]]
        v.orientation = e["old"].copy()
        v.cost, v.confidence, v.overlap_ratio = (
            e["old_cost"],
            e["old_confidence"],
            e["old_overlap_ratio"],
        )
    mic.write(str(out / "recon_before_pass.mic"))
    gz_file(out / "recon_before_pass.mic")
    gz_file(out / "recon.mic")
    st = bfs.stats
    info = dict(
        name=kind,
        keys=K.KEYS[kind],
        n_voxels=len(mic.voxels),
        n_processed=len(order),
        wall_total_s=wall,
        setup_s=t_setup,
        timing_label=label,
        stats={k: v for k, v in st.items() if k != "reseeds"},
        stats_before_pass={
            k: v for k, v in captured["stats_before"].items() if isinstance(v, (int, float))
        },
        evals_total_before_pass=K.evals_total(st),
        evals_total_with_pass=K.evals_total(st)
        + int(st.get("n_evals_boundary_pass", 0))
        + int(st.get("n_evals_boundary_score", 0)),
        n_boundary_log=len(bfs.boundary_log),
        index_pairing_s=float(st.get("index_pairing_s", 0.0)),
    )
    (out / "run.json").write_text(json.dumps(info, indent=2, default=float) + "\n")
    (out / "run.done").write_text("done\n")
    print(kind, "done", wall, "s", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kind", required=True, choices=sorted(K.KEYS))
    ap.add_argument("--smoke", type=int, default=0, help="pilot: first N voxels, own folder")
    a = ap.parse_args()
    out = K.CACHE / (f"smoke_{a.kind}" if a.smoke else a.kind)
    out.mkdir(parents=True, exist_ok=True)
    if a.smoke:  # a smoke run limits the voxel count through the BFS call
        orig_build = C.build

        def build(*aa: Any, **kk: Any) -> Any:
            setup, bfs, t = orig_build(*aa, **kk)
            rs = bfs.reconstruct_sample

            def limited(*x: Any, **y: Any) -> Any:
                y["max_voxels"] = a.smoke
                return rs(*x, **y)

            bfs.reconstruct_sample = limited
            return setup, bfs, t

        C.build = build
    run(a.kind, out, a.smoke or None)


if __name__ == "__main__":
    main()
