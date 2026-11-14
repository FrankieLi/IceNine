"""AR5 step 9, phase B: a library run with the key ``SeedRescueCMA 0.25 1000`` on (``LibG``).

  RESEED_SAMPLE=mid uv run python scripts/autoresearch/ar5/rescue_lib_run.py --kind LibG|LibS11

Everything is the library's own: the keys are written to the config (``ar3_common.build``
appends them to ``key_on_fused.config``, as ``keys_run.py`` of AR5 step 5 does); no method is
replaced and nothing is counted by patching. ``LibG`` is run ``S11L`` of step 5 (the keys of
``keys_common.KEYS["S11L"]``) plus ``SeedRescueCMA 0.25 1000``: the library form of the script arm
``Gl`` (``rescue_gl.py``); ``LibS11`` is ``S11L`` with the key off (a control against the stored
script run ``S11`` of step 3). The names differ from ``Gl`` by more than case: macOS file names
are case-insensitive. ``default_rng(0)``, clean images, Q-max 8, one worker. Added here, as
in ``keys_run.py``: a capture of the boundary pass's result, a progress file every 5 minutes and a
watchdog at ``--cap-min`` minutes. The map before the pass is rebuilt from ``bfs.boundary_log``.

Writes ``<cache>/<kind>_<sample>/`` (cache = ``benchmarks/autoresearch/ar5/rescue/cache``,
untracked): ``recon.mic.gz`` (after the pass), ``recon_before_pass.mic.gz``, ``records.json.gz``,
``index_log.json.gz``, ``pass_tries.npz``, ``run.json``, ``preflight.json``, ``run.config``,
``progress.json``; ``run.done`` at the end. Refuses to start if the folder exists. Wall times are
contended, not timing results.
"""

import argparse
import gzip
import json
import os
import shutil
import sys
import threading
import time
from pathlib import Path
from typing import Any, Dict

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "ar3"))

import reseed_common as RC  # noqa: E402  (first: RESEED_SAMPLE -> AR3_SAMPLE)

import ar3_common as C  # noqa: E402
import gradient_common as G  # noqa: E402
import numpy as np  # noqa: E402

S8_KEYS = [
    "BFSNeighborCMAFirst 1",
    "CMASigma0 0.05",
    "CMANeighborMaxEvals 60",
    "BoundaryPass 1",
    "SeedIndexCandidates 20 5",
    "SeedIndexStopHitRatio 0.95",
]
S11L = S8_KEYS + ["BFSNeighborCMANoFallback 1", "SeedIndexCMAFirst 200"]
KEYS = {"LibG": S11L + ["SeedRescueCMA 0.25 1000"], "LibS11": S11L}
CACHE = Path(
    os.environ.get(
        "RESCUE_CACHE", str(C.ICE / "benchmarks" / "autoresearch" / "ar5" / "rescue" / "cache")
    )
)


def gz_file(path: Path) -> None:
    with open(path, "rb") as fi, open(str(path) + ".gz", "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as fo:
            shutil.copyfileobj(fi, fo)
    path.unlink()


def numeric(d: Dict[str, Any]) -> Dict[str, Any]:
    return {k: v for k, v in d.items() if isinstance(v, (int, float)) and not isinstance(v, bool)}


def evals_total(st: Dict[str, Any]) -> int:
    """The run scripts' own total of a BFS (the pass is not in it)."""
    return int(
        sum(st[k] for k in ("n_evals_seed", "n_evals_neighbor", "n_evals_revisit"))
        + st.get("n_evals_reseed", 0)
        + st.get("n_evals_restart", 0)
    )


def run(kind: str, out: Path, cap_min: float) -> None:
    from icenine import boundary_pass as LBP

    label = C.preflight_label(out)
    setup, bfs, t_setup = C.build(
        C.CONFIG, out, sample=G.GRID_MIC, stack_name=G.STACK, extra="\n".join(KEYS[kind])
    )
    mic = setup.sample.get_mic()
    p = setup.search_params
    assert p.seed_by_indexing and p.fused_cost and p.local_optimizer == "mc"
    assert p.bfs_neighbor_cma_first and p.bfs_neighbor_cma_no_fallback and p.boundary_pass
    assert (p.seed_index_n_cost, p.seed_index_n_refine) == (20, 5) and p.seed_index_stop_hit == 0.95
    assert p.cma_sigma0_deg == 0.05 and p.cma_neighbor_max_evals == 60
    assert p.seed_index_cma_first == 200
    want = (0.25, 1000) if kind == "LibG" else (0.0, 0)
    assert (p.seed_rescue_sigma_deg, p.seed_rescue_evals) == want
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
                    json.dumps(dict(elapsed_s=time.time() - t0, stats=numeric(bfs.stats))) + "\n"
                )
            except Exception:  # a progress write must never stop the run
                pass

    def watchdog() -> None:
        time.sleep(cap_min * 60.0)
        (out / "capped.json").write_text(
            json.dumps(dict(cap_min=cap_min, n_records=len(bfs.records), stats=numeric(bfs.stats)))
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
        sample=RC.SAMPLE,
        keys=KEYS[kind],
        n_voxels=len(mic.voxels),
        n_processed=len(order),
        wall_total_s=wall,
        setup_s=t_setup,
        timing_label=label,
        wall_label="contended (other jobs on the machine); not a timing result",
        stats={k: v for k, v in st.items() if k != "reseeds"},
        stats_before_pass=numeric(captured["stats_before"]),
        evals_total_before_pass=evals_total(captured["stats_before"]),
        evals_total_with_pass=evals_total(st) + int(st.get("n_evals_boundary_pass", 0)),
        n_boundary_log=len(bfs.boundary_log),
    )
    (out / "run.json").write_text(json.dumps(info, indent=2, default=float) + "\n")
    (out / "run.done").write_text("done\n")
    print(kind, RC.SAMPLE, "done", wall, "s", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kind", required=True, choices=sorted(KEYS))
    ap.add_argument("--cap-min", type=float, default=240.0)
    a = ap.parse_args()
    assert RC.SAMPLE in ("mid", "grad", "big"), "held is not run in phase B"
    out = CACHE / f"{a.kind}_{RC.SAMPLE}"
    if out.exists():
        raise SystemExit(f"refusing to start: {out} exists")
    out.mkdir(parents=True)
    run(a.kind, out, a.cap_min)
    os._exit(0)


if __name__ == "__main__":
    main()
