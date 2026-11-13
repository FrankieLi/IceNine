"""AR5 step 9, phase B pilot: the library's ``rescue_seed`` on the replay's cases (``mid``).

  RESEED_SAMPLE=mid uv run python scripts/autoresearch/ar5/rescue_lib_pilot.py [--n N]

For the first ``N`` cases of ``rescue/cases_<sample>.json`` (the failed seed steps of arm ``A``
with their kept ends, as the phase A replay read them) it calls
``AdaptiveVoxelReconstructor.rescue_seed`` (accept at 0.95, 0.25 degree, 1,000 evaluations) and
compares every try (kept end tried, hit ratio, cost, evaluations) and the outcome with the stored
replay (``rescue/replay_<sample>.json``, made by the script rule ``rescue_rules.rescue_choose`` +
``rescue_common.cma_rescue``). Same inputs, so they must be equal. Prints one line per case and
writes ``rescue/cache/lib_pilot_<sample>.json`` (untracked). Wall time: contended, not a timing.
"""

import argparse
import json
import os
import time
from typing import Any, Dict, List

import rescue_common as RCM  # noqa: E402  (first)

import fit_common as FC  # noqa: E402
import numpy as np  # noqa: E402
import reseed_common as RC  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=0, help="first N cases (0: all)")
    a = ap.parse_args()
    s = FC.SAMPLE
    assert s in ("grad", "mid"), "the pilot reads a disclosed sample's replay"
    cases = json.loads((RCM.RESCUE / f"cases_{s}.json").read_text())
    replay = {
        int(r["voxel"]): r
        for r in json.loads((RCM.RESCUE / f"replay_{s}.json").read_text())["rows"]
    }
    if a.n:
        cases = cases[: a.n]
    env = RC.make_env()
    rec = env.hooks.rec
    rows: List[Dict[str, Any]] = []
    t0 = time.time()
    for c in cases:
        v = int(c["voxel"])
        verts = env.hooks.vertices(v)
        ph = int(env.mic.voxels[v].phase)
        kept = c["kept"]
        out = rec.rescue_seed(
            verts,
            ph,
            [np.asarray(k["R"], np.float64).reshape(3, 3) for k in kept],
            [float(k["cost"]) for k in kept],
            0.95,
            0.25,
            1000,
        )
        r = replay[v]
        lib = [(t["k"], t["hit"], t["cost"], t["evals"]) for t in out["tries"]]
        scr = [(t["k"], t["hit"], t["cost"], t["evals"]) for t in r["tries"]]
        same = lib == scr and out["rescued"] == r["rescued"] and out["evals"] == r["evals"]
        rows.append(
            dict(voxel=v, same=same, rescued=out["rescued"], tries=len(lib), evals=out["evals"])
        )
        print(
            f"{s} voxel {v}: same={same} rescued={out['rescued']} tries={len(lib)} "
            f"evals={out['evals']} {time.time() - t0:.0f}s",
            flush=True,
        )
    res = dict(
        sample=s,
        n=len(rows),
        all_same=all(r["same"] for r in rows),
        rows=rows,
        wall_s=time.time() - t0,
        wall_label="contended; not a timing result",
    )
    (RCM.RCACHE / f"lib_pilot_{s}.json").write_text(json.dumps(res, indent=1) + "\n")
    print("all_same", res["all_same"], flush=True)
    os._exit(0)


if __name__ == "__main__":
    main()
