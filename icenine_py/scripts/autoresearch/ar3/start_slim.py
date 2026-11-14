"""AR3 step 9, part 3: slim copies of the whole-layer runs made in this task (one sample).

  RESEED_SAMPLE=mid RESEED_CACHE=<start/cache> uv run python scripts/autoresearch/ar3/start_slim.py

Reads the full run folders and the pass cache (not in git); writes, under ``start/``:
``runs_<sample>.json`` (for each run present among D4, S5, S8 and S9: the scalar run statistics, the
pass counters, the seed-step counters of S8 and S9, and for S9 the counters of F1 and F2) and
``errors_<sample>.npz`` (float32 error to the truth of every voxel for each run present, before
and after the pass; the voxel's level and grain). For ``grad`` it also writes
``check_logged_grad.json``: the logged rerun of S5 against the stored run S5 (the baseline check).
"""

import json
from typing import Any, Dict

import numpy as np
import reseed_common as RC

import gradient_common as G  # noqa: E402
import start_part1 as P1  # noqa: E402

N = 24570
KINDS = ("D4", "S5", "S8", "S9")


def scalars(st: Dict[str, Any]) -> Dict[str, Any]:
    return {k: v for k, v in st.items() if isinstance(v, (int, float))}


def run_info(d: Any) -> Dict[str, Any]:
    r = json.loads((d / "run.json").read_text())
    r["stats"] = scalars(r["stats"])
    for k in ("pick", "wrapper"):
        r.pop(k, None)
    return r


def main() -> None:
    T = G.Truth()
    PAS = RC.CACHE / f"passes_{RC.SAMPLE}"
    runs: Dict[str, Any] = {}
    errs: Dict[str, Any] = dict(level=T.level.astype(np.float32), grain=T.gmap)
    for k in KINDS:
        d = RC.CACHE / f"{k}_{RC.SAMPLE}"
        if not (d / "run.json").exists():
            continue
        info = run_info(d)
        for nm in ("s8.json", "s9.json"):
            if (d / nm).exists():
                info[nm.split(".")[0]] = json.loads((d / nm).read_text())
        if (d / "e1_log.json.gz").exists():
            log = RC.read_gz(d / "e1_log.json.gz")
            info["e1"] = dict(
                steps=len(log),
                stopped_early=sum(1 for e in log if e["stop"] >= 0),
                evals=int(sum(e.get("evals", 0) for e in log)),
                costing=int(sum(e["costed"] for e in log)),
                refine_evals=int(sum(sum(e["fin_evals"]) for e in log)),
            )
        pj = PAS / f"D4_{k}" / "pass.json"
        if pj.exists():
            info["pass"] = json.loads(pj.read_text())
        runs[k] = info
        R, _ = G.read_map(d / "recon.mic.gz", N)
        errs[f"{k}_before"] = T.err(R).astype(np.float32)
        if (PAS / f"after_{k}.npy").exists():
            errs[f"{k}_after"] = T.err(np.load(PAS / f"after_{k}.npy")).astype(np.float32)
        print(
            k,
            "evals",
            info["evals_total"],
            "wrong before/after:",
            int((errs[f"{k}_before"] > G.WRONG).sum()),
            int((errs.get(f"{k}_after", np.zeros(1)) > G.WRONG).sum()),
        )
    P1.START.mkdir(parents=True, exist_ok=True)
    (P1.START / f"runs_{RC.SAMPLE}.json").write_text(json.dumps(runs, indent=1) + "\n")
    np.savez_compressed(P1.START / f"errors_{RC.SAMPLE}.npz", **errs)
    if RC.SAMPLE == "grad" and (P1.LOGGED["grad"] / "run.json").exists():
        run = json.loads((P1.LOGGED["grad"] / "run.json").read_text())
        stored = RC.stored_run("S5")
        keys = [
            "n_seeds",
            "n_reseeds",
            "n_evals_seed",
            "n_evals_reseed",
            "n_evals_neighbor",
            "n_evals_revisit",
            "n_neighbor_fits",
            "n_revisit_attempted",
            "index_evals",
        ]
        same = {k: [run["stats"][k], stored["stats"][k]] for k in keys}
        st_n, ref = run["stats"], stored["stats"]
        res_same = [
            (a["voxel"], a["cluster_size"], a["full_search_evaluations"], a["accepted"])
            == (b["voxel"], b["cluster_size"], b["full_search_evaluations"], b["accepted"])
            for a, b in zip(st_n["reseeds"], ref["reseeds"])
        ]
        chk = dict(
            counters_logged_vs_stored=same,
            counters_equal=all(a == b for a, b in same.values()),
            reseeds_logged=len(st_n["reseeds"]),
            reseeds_stored=len(ref["reseeds"]),
            reseeds_equal=bool(all(res_same) and len(st_n["reseeds"]) == len(ref["reseeds"])),
            evals_total_logged=run["evals_total"],
            evals_total_stored=stored["evals_total"],
            wall_logged_s=run["wall_total_s"],
            timing_label=run["timing_label"],
        )
        (P1.START / "check_logged_grad.json").write_text(json.dumps(chk, indent=1) + "\n")
        print(json.dumps(chk)[:600])


if __name__ == "__main__":
    main()
