"""AR3 step 10: score the library runs against the truth, after the runs.

The library never sees the truth; only this script does.

  libkeys_score.py [--kinds D4L S5L S8L]

For each run: the voxels whose error is above 1 degree ("wrong", step 7's definition) in the map
before the pass (rebuilt from the pass's log) and after it, and the same for the stored run (its
``recon.mic.gz`` and the stored ``after_<map>.npy``). The library map after the pass is compared
voxel by voxel with the stored one: fixed (wrong in the stored map only), broken (wrong in the
library map only), both wrong. Writes ``keys/score.json`` and ``keys/errors.npz`` (per-voxel
errors, float32; slim).
"""

import argparse
import json
from typing import Any, Dict

import libkeys_common as K  # first: it sets AR3_SAMPLE before the common modules read it
import numpy as np

import gradient_common as G  # noqa: E402
from stats import wilson  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kinds", nargs="+", default=["D4L", "S5L", "S8L"])
    a = ap.parse_args()
    T = G.Truth()
    assert T.n == K.N
    out: Dict[str, Any] = {}
    errs: Dict[str, np.ndarray] = {}
    for kind in a.kinds:
        run = K.CACHE / kind
        _pdir, _before, after_npy = K.STORED_PASS[kind]
        maps = {
            "lib_before": G.read_map(run / "recon_before_pass.mic.gz", K.N)[0],
            "lib_after": G.read_map(run / "recon.mic.gz", K.N)[0],
            "stored_before": G.read_map(K.STORED[kind] / "recon.mic.gz", K.N)[0],
            "stored_after": np.load(after_npy).astype(np.float64),
        }
        e = {k: T.err(R) for k, R in maps.items()}
        for k, v in e.items():
            errs[f"{kind}_{k}"] = v.astype(np.float32)
        w = {k: v > G.WRONG for k, v in e.items()}
        info: Dict[str, Any] = {k: int(v.sum()) for k, v in w.items()}
        info["wrong_after_ci95"] = [round(x, 5) for x in wilson(int(w["lib_after"].sum()), K.N)]
        info["expected"] = dict(
            wrong_before=K.EXPECT[kind]["wrong_before"], wrong_after=K.EXPECT[kind]["wrong_after"]
        )
        for tag in ("before", "after"):
            lib, st = w[f"lib_{tag}"], w[f"stored_{tag}"]
            info[f"paired_{tag}"] = dict(
                fixed_vs_stored=int((st & ~lib).sum()),
                broken_vs_stored=int((~st & lib).sum()),
                both_wrong=int((st & lib).sum()),
                n_voxels_with_different_error=int((e[f"lib_{tag}"] != e[f"stored_{tag}"]).sum()),
            )
        # the pass on the library's own map: fixed / broken against the library map before it
        lb, la = w["lib_before"], w["lib_after"]
        info["pass_on_library_map"] = dict(
            fixed=int((lb & ~la).sum()),
            broken=int((~lb & la).sum()),
            both_wrong=int((lb & la).sum()),
        )
        out[kind] = info
    (K.OUT / "score.json").write_text(json.dumps(out, indent=1) + "\n")
    np.savez_compressed(K.OUT / "errors.npz", **errs)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
