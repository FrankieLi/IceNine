"""AR3 step 8, part 3: slim copies of run S8 (one sample) and of the runs it is compared with.

  RESEED_SAMPLE=big  uv run python scripts/autoresearch/ar3/reseed_s8_slim.py
  RESEED_SAMPLE=grad uv run python scripts/autoresearch/ar3/reseed_s8_slim.py

Needs the full run folders and the pass cache (not in git). Writes, under ``reseed/slim/``:
``s8_<sample>.json`` (scalar run statistics, E1 counters, pass counters),
``baselines_<sample>.json`` (the same for the stored runs S5 and D4, and their passes) and
``errors_<sample>.npz`` (float32
error to the truth of every voxel: S8, S5 and D4, before and after the pass; the voxel's level and
grain).
"""

import json
from typing import Any, Dict

import numpy as np
import reseed_common as RC

import boundary_common as B  # noqa: E402
import gradient_common as G  # noqa: E402

SLIM = RC.OUT / "slim"
N = 24570


def scalars(st: Dict[str, Any]) -> Dict[str, Any]:
    return {k: v for k, v in st.items() if isinstance(v, (int, float))}


def run_info(p: Any) -> Dict[str, Any]:
    r = json.loads((p / "run.json").read_text())
    r["stats"] = scalars(r["stats"])
    for k in ("pick", "wrapper"):
        r.pop(k, None)
    return r


def main() -> None:
    T = G.Truth()
    S8 = RC.CACHE / f"S8_{RC.SAMPLE}"
    PAS = RC.CACHE / f"passes_{RC.SAMPLE}"
    s8 = run_info(S8)
    s8["s8"] = json.loads((S8 / "s8.json").read_text())
    pj = json.loads((PAS / "D4_S8" / "pass.json").read_text())
    log = RC.read_gz(S8 / "e1_log.json.gz")
    s8["pass"] = pj
    s8["e1"] = dict(
        steps=len(log),
        stopped_early=sum(1 for e in log if e["stop"] >= 0),
        evals=int(sum(e.get("evals", 0) for e in log)),
        costing=int(sum(e["costed"] for e in log)),
        refine_evals=int(sum(sum(e["fin_evals"]) for e in log)),
        stop_at=[sum(1 for e in log if e["stop"] == i) for i in range(5)],
    )
    (SLIM / f"s8_{RC.SAMPLE}.json").write_text(json.dumps(s8, indent=1) + "\n")

    base: Dict[str, Any] = {}
    errs: Dict[str, Any] = dict(level=T.level.astype(np.float32), grain=T.gmap)
    R8, _ = G.read_map(S8 / "recon.mic.gz", N)
    errs["S8_before"] = T.err(R8).astype(np.float32)
    errs["S8_after"] = T.err(np.load(PAS / "after_S8.npy")).astype(np.float32)
    if RC.SAMPLE == "big":
        BIG = RC.C.ICE / "benchmarks" / "autoresearch" / "ar3" / "big"
        z = np.load(BIG / "errors.npz")
        for k in ("S5", "D4"):
            errs[f"{k}_before"], errs[f"{k}_after"] = z[f"{k}_before"], z[f"{k}_after"]
            r = json.loads((BIG / "slim" / f"{k}_run.json").read_text())
            r["pass"] = json.loads((BIG / "passes" / f"{k}_pass.json").read_text())
            base[k] = r
    else:
        A = RC._GRAD_WT
        z = np.load(A / "seed" / "runs" / "errors.npz")
        for k in ("S5", "D4"):
            errs[f"{k}_before"], errs[f"{k}_after"] = z[f"{k} before"], z[f"{k} after"]
        base["S5"] = run_info(A / "seed" / "runs" / "S5")
        base["S5"]["pass"] = json.loads((A / "seed/cache/passes/D4_S5/pass.json").read_text())
        base["D4"] = run_info(A / "fine" / "runs" / "D4")
        base["D4"]["pass"] = json.loads((A / "boundary/runs/D4/pass.json").read_text())
    (SLIM / f"baselines_{RC.SAMPLE}.json").write_text(json.dumps(base, indent=1) + "\n")
    np.savez_compressed(SLIM / f"errors_{RC.SAMPLE}.npz", **errs)
    for k in ("S8", "S5", "D4"):
        print(
            k,
            "wrong before/after:",
            int((errs[f"{k}_before"] > G.WRONG).sum()),
            int((errs[f"{k}_after"] > G.WRONG).sum()),
        )
    assert B.G is G


if __name__ == "__main__":
    main()
