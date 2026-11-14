"""AR3 step 10, bit-identity check (c): the BFS with every new key absent, on the pilot region.

  libkeys_off.py --tag parent      # run with PYTHONPATH = a copy of the parent commit's package
  libkeys_off.py --tag branch      # run with the task branch's package
  libkeys_off.py --compare         # compare the two folders

``key_on_fused.config`` on the 2,000-voxel pilot region (clean images, ``default_rng(0)``), none of
the new keys written. Writes ``<CACHE>/off_<tag>/``: ``recon.mic`` (the written file), ``orient.
npy`` (orientations), ``records.json``, ``stats.json`` and ``logs.json`` (stats without the wall
times; the reseed list and the index log without their seconds). ``--compare`` checks that all
five are equal and writes ``keys/off_compare.json``.
"""

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import libkeys_common as K  # noqa: E402

import ar3_common as C  # noqa: E402


def drop_time(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {
            k: drop_time(v) for k, v in obj.items() if not K.is_time_key(k) and "seconds" not in k
        }
    if isinstance(obj, list):
        return [drop_time(v) for v in obj]
    return obj


def run(tag: str) -> None:
    import icenine

    out = K.CACHE / f"off_{tag}"
    out.mkdir(parents=True, exist_ok=True)
    label = C.preflight_label(out)
    setup, bfs, _t = C.build(C.CONFIG, out)  # the pilot region, clean.npy
    mic = setup.sample.get_mic()
    p = setup.search_params
    assert p.seed_by_indexing and p.fused_cost and p.bfs_reseed
    assert not getattr(p, "boundary_pass", False)  # the parent's package has no such field
    assert not getattr(p, "bfs_neighbor_cma_first", False)
    t0 = time.time()
    bfs.reconstruct_sample(output_mic=str(out / "recon.mic"), rng=np.random.default_rng(0))
    wall = time.time() - t0
    np.save(out / "orient.npy", K.orientations(mic))
    recs = {
        str(k): {f: v for f, v in vars(r).items() if f != "wall_s"} for k, r in bfs.records.items()
    }
    (out / "records.json").write_text(json.dumps(recs, default=float) + "\n")
    stats = drop_time({k: v for k, v in bfs.stats.items() if k != "reseeds"})
    (out / "stats.json").write_text(json.dumps(stats, default=float, sort_keys=True) + "\n")
    logs = dict(
        reseeds=drop_time(bfs.stats.get("reseeds", [])),
        index_log=drop_time(bfs.index_log),
        reseed_log=drop_time(bfs.reseed_log),
    )
    (out / "logs.json").write_text(json.dumps(logs, default=float, sort_keys=True) + "\n")
    (out / "info.json").write_text(
        json.dumps(dict(tag=tag, package_folder=Path(icenine.__file__).parts[-4],
                        has_boundary_pass=hasattr(p, "boundary_pass"),
                        n_voxels=len(mic.voxels), wall_s=wall, timing_label=label), indent=1)
        + "\n"
    )  # fmt: skip
    print(tag, "done", wall, "s", flush=True)


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def compare() -> None:
    a, b = K.CACHE / "off_parent", K.CACHE / "off_branch"
    res: Dict[str, Any] = {}
    for f in ("recon.mic", "orient.npy", "records.json", "stats.json", "logs.json"):
        res[f] = dict(equal=sha(a / f) == sha(b / f), parent_sha=sha(a / f)[:16],
                      branch_sha=sha(b / f)[:16], bytes=(a / f).stat().st_size)  # fmt: skip
    res["all_equal"] = all(v["equal"] for v in res.values())
    res["n_voxels"] = json.loads((a / "info.json").read_text())["n_voxels"]
    (K.OUT / "off_compare.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", choices=["parent", "branch"])
    ap.add_argument("--compare", action="store_true")
    a = ap.parse_args()
    compare() if a.compare else run(a.tag)
