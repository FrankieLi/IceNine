"""AR5 step 5, check 5: the library pass on a stored map gives the stored tries and totals.

  AR3_SAMPLE=big uv run python scripts/autoresearch/ar5/keys_pass.py --kind S11L
  uv run python scripts/autoresearch/ar5/keys_pass.py --kind S11L            # first gradient

Reads the map the script's pass read (``before_<name>.npy``: the stored run's ``.mic`` with its
rounded angles) into the BFS's map in memory, with every voxel FITTED, and applies
``icenine.boundary_pass.run_boundary_pass`` to it (the same function as in step 10's
``libkeys_pass.py``). Exact: its tries (voxel, candidate, end as float32, cost, hit ratio,
confidence, evaluations) and its evaluation total equal the stored ``tries.npz``; the total with
the pass is the stored run's BFS total plus the pass's tries. Writes
``benchmarks/autoresearch/ar5/keys/pass_check_<kind>_<sample>.json``. The pass in a library run
(``keys_run.py``) reads the unrounded map and may differ; that is reported by ``keys_check.py``.
"""

import argparse
import json
import time
from typing import Any, Dict

import keys_common as K
import numpy as np

import ar3_common as C  # noqa: E402
import gradient_common as G  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kind", required=True, choices=sorted(K.KEYS))
    a = ap.parse_args()
    from icenine import boundary_pass as LBP
    from icenine.mic_file import ReconstructionState

    key = (a.kind, K.SAMPLE)
    pdir, before_npy, after_npy = K.STORED_PASS[key]
    out = K.CACHE / f"pass_check_{a.kind}_{K.SAMPLE}"
    if out.exists():
        raise SystemExit(f"refusing to start: {out} exists")
    out.mkdir(parents=True)
    label = C.preflight_label(out)
    setup, bfs, _t = C.build(C.CONFIG, out, sample=G.GRID_MIC, stack_name=G.STACK)
    mic = setup.sample.get_mic()
    R = np.load(before_npy).astype(np.float64)
    assert len(R) == len(mic.voxels)
    for v, r in zip(mic.voxels, R):
        v.orientation = r.copy()
        v.reconstruction_id = ReconstructionState.FITTED
    bfs.stats = bfs._new_stats()
    t0 = time.time()
    res = LBP.run_boundary_pass(bfs, mic)
    wall = time.time() - t0
    s = np.load(pdir / "tries.npz")
    n = len(res.tries)
    lib = dict(
        voxel=np.array([t["voxel"] for t in res.tries]),
        cand=np.array([t["cand"] for t in res.tries]),
        end=np.array([t["end"] for t in res.tries], np.float32).reshape(-1, 3, 3),
        cost=np.array([t["cost"] for t in res.tries]),
        hit=np.array([t["hit"] for t in res.tries]),
        conf=np.array([t["conf"] for t in res.tries]),
        evals=np.array([t["evals"] for t in res.tries]),
    )
    stored_run = json.loads((K.STORED[key] / "run.json").read_text())
    stored_bfs = stored_run.get("evals_total", K.evals_total(stored_run["stats"]))
    info: Dict[str, Any] = dict(
        kind=a.kind,
        sample=K.SAMPLE,
        n_tries=(n, int(len(s["voxel"]))),
        evals=(int(res.n_evals_tries), int(s["evals"].sum())),
        total_with_pass=(
            stored_bfs + int(res.n_evals_tries),
            K.EXPECT.get(key, {}).get("with_pass"),
        ),
        wall_s=wall,
        timing_label=label,
    )
    if n == len(s["voxel"]):
        for f in ("voxel", "cand", "evals", "cost", "hit", "conf"):
            bad = np.where(lib[f] != s[f])[0]
            info[f"differ_{f}"] = int(len(bad))
            if len(bad):
                info[f"first_{f}"] = dict(
                    try_index=int(bad[0]), voxel=int(lib["voxel"][bad[0]]),
                    new=float(lib[f][bad[0]]), stored=float(s[f][bad[0]]),
                )  # fmt: skip
        info["differ_end"] = int((~(lib["end"] == s["end"]).reshape(n, -1).all(1)).sum())
    after = np.load(after_npy)
    got = K.orientations(mic)
    changed = np.where((got != R).reshape(len(R), -1).any(1))[0]
    stored_changed = np.where((after != R).reshape(len(R), -1).any(1))[0]
    info["voxels_changed"] = (int(len(changed)), int(len(stored_changed)))
    info["same_changed_voxels"] = bool(set(changed.tolist()) == set(stored_changed.tolist()))
    K.OUT.mkdir(parents=True, exist_ok=True)
    (K.OUT / f"pass_check_{a.kind}_{K.SAMPLE}.json").write_text(json.dumps(info, indent=1) + "\n")
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
