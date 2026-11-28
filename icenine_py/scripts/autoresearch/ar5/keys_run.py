"""AR5 step 5: a library run of the two keys (parity runs; see ``keys_common``).

  AR3_SAMPLE=big uv run python scripts/autoresearch/ar5/keys_run.py --kind S11L   # big
  uv run python scripts/autoresearch/ar5/keys_run.py --kind S11L                  # first gradient
  ... --smoke N                                       # pilot: first N voxels, own folder

Everything is the library's own: the keys are written to the config (``ar3_common.build`` appends
them to ``key_on_fused.config``); no method is replaced and nothing is counted by patching. Added
here: a capture of the boundary pass's result (the function is called through its module), a
progress file every 5 minutes and a watchdog at 90 minutes (a run on `big` takes about
15 minutes; the run on the gradient sample is capped at 90). The map before
the pass is rebuilt from ``bfs.boundary_log`` and written too, as in step 10.

Writes ``<CACHE>/<kind>_<sample>/``: ``recon.mic.gz`` (after the pass),
``recon_before_pass.mic.gz``,
``records.json.gz``, ``run.json``, ``pass_tries.npz``, ``preflight.json``, ``run.config``,
``progress.json``, ``index_log.json.gz``; ``run.done`` at the end. Refuses to start if the folder
exists.
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

import keys_common as K
import numpy as np

import ar3_common as C  # noqa: E402
import gradient_common as G  # noqa: E402

CAP_MIN = 90


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
    assert (p.seed_index_n_cost, p.seed_index_n_refine) == (20, 5) and p.seed_index_stop_hit == 0.95
    assert p.cma_sigma0_deg == 0.05 and p.cma_neighbor_max_evals == 60
    assert p.bfs_neighbor_cma_no_fallback == (kind in ("S10L", "S11L"))
    assert p.seed_index_cma_first == (200 if kind == "S11L" else 0)
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
    (out / "index_log.json").write_text(json.dumps(bfs.index_log, default=float) + "\n")
    gz_file(out / "index_log.json")
    for e in bfs.boundary_log:  # rebuild the map before the pass from the pass's log
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
        sample=K.SAMPLE,
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
        evals_total_before_pass=K.evals_total(captured["stats_before"]),
        evals_total_with_pass=K.evals_total(st) + int(st.get("n_evals_boundary_pass", 0)),
        n_boundary_log=len(bfs.boundary_log),
        index_pairing_s=float(st.get("index_pairing_s", 0.0)),
    )
    (out / "run.json").write_text(json.dumps(info, indent=2, default=float) + "\n")
    (out / "run.done").write_text("done\n")
    print(kind, K.SAMPLE, "done", wall, "s", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kind", required=True, choices=sorted(K.KEYS))
    ap.add_argument("--smoke", type=int, default=0, help="pilot: first N voxels, own folder")
    a = ap.parse_args()
    out = K.CACHE / (f"smoke_{a.kind}_{K.SAMPLE}" if a.smoke else f"{a.kind}_{K.SAMPLE}")
    if out.exists():
        raise SystemExit(f"refusing to start: {out} exists")
    out.mkdir(parents=True)
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
