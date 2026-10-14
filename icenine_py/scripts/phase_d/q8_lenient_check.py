#!/usr/bin/env python3
"""Q8 step 8: why do the lenient wrong counts of the Step 4 re-runs differ from the stored ones?

Step 4 re-ran the lenient rule (arm `off`) on the E0 200-voxel set and got 35 / 27 wrong (clean /
old realistic) at min_sin_eta 0.3 and 41 / 25 at min_sin_eta 0, where the stored runs of the Q3
study and of the ray-overlap study have 40 / 25 and 39 / 18. The C++ rule and ray arms reproduce
their stored counts and evaluations exactly, so only the lenient arm differs. The hypothesis is the
level-0 screen: the stored lenient runs used the vertex-0 screen (the library default before
2026-10-12; MIGRATION_HISTORY 2026-10-12: "Centroid start point ... lenient rule on E0 clean 39 ->
41 wrong"), the re-runs use the centroid screen, which every other arm already used.

Test: arm `offv0` of e0_cpprule.py is the lenient rule with the vertex-0 screen
(scripts/phase_d/run_q8_offv0.sh). This script compares, per (min_sin_eta, images), the stored
count, the `offv0` count and the current `off` count (wrong = error over 1 degree), and pairs
`offv0` against `off` per voxel (fixed = wrong under `off`, right under `offv0`).

  uv run python scripts/phase_d/q8_lenient_check.py

Writes benchmarks/q8_nodrop/lenient_check_summary.json and lenient_check_tables.md.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICE / "scripts" / "findoptimal_robustness"))
sys.path.insert(0, str(ICE / "scripts" / "common"))
import common as C  # noqa: E402
import doc_tables  # noqa: E402
import q8_e0_summary as Q  # noqa: E402
import stats as S  # noqa: E402

OUT = ICE / "benchmarks" / "q8_nodrop"
LABEL = {"clean": "clean", "all": "old realistic (10% drop)"}
# (min_sin_eta, current arm id, vertex-0 arm id, source of the stored count)
SETTINGS = [
    (0.3, "off", "offv0", "Q3 study (summary_on.json, faithful-MC rerun)"),
    (0.0, "off_ms0", "offv0_ms0", "ray-overlap study (e2e_summary.json, arm old)"),
]


def stored(var: str, msin: float) -> Tuple[int, float]:
    """(wrong count, mean evaluations per voxel) of the stored lenient run."""
    if msin == 0.3:
        z = json.loads(Q.EARLIER_Q3.read_text())
        return int(z[var]["wrong_before"]), float(z[var]["mean_evals_before"])
    z = json.loads(Q.EARLIER_RAY.read_text())
    key = f"{'clean' if var == 'clean' else 'realistic'}|{Q.EARLIER_RAY_KEY['off_ms0']}"
    return int(z[key]["wrong"]), float(z[key]["mean_evals"])


HEADER = {
    "msin": "min_sin_eta",
    "images": "images",
    "stored": "stored count (vertex-0 screen)",
    "vertex0": "re-run, vertex-0 screen (arm offv0)",
    "interval": "95% interval",
    "centroid": "re-run, centroid screen (arm off, Step 4)",
    "fixed": "centroid fixes (wrong with vertex-0, right with centroid)",
    "broken": "centroid breaks",
    "both_wrong": "wrong with both",
    "ev_st": "stored evaluations per voxel",
    "ev_v0": "evaluations per voxel, vertex-0",
    "ev_c": "evaluations per voxel, centroid",
}


def main() -> None:
    vox = Q.voxels()
    rows: List[Dict[str, Any]] = []
    summ: Dict[str, Any] = {}
    for msin, cur_id, v0_id, src in SETTINGS:
        for var in ("clean", "all"):
            cur, v0 = Q.load(cur_id, var, vox), Q.load(v0_id, var, vox)
            if cur is None or v0 is None:
                continue
            n = len(vox)
            k_cur, k_v0 = int(cur["wrong"].sum()), int(v0["wrong"].sum())
            sp = Q.split(cur["wrong"], v0["wrong"])  # a = current, ref = vertex-0 start
            lo, hi = S.wilson(k_v0, n)
            st, st_ev = stored(var, msin)
            ev_v0 = float((v0["evals_global"] + v0["evals_local"]).mean())
            ev_c = float((cur["evals_global"] + cur["evals_local"]).mean())
            summ[f"{msin}|{var}"] = dict(
                msin=msin, images=var, stored=st, stored_evals=st_ev, stored_source=src,
                vertex0=k_v0, vertex0_interval=[lo, hi], centroid=k_cur,
                evals_vertex0=ev_v0, evals_centroid=ev_c, paired_centroid_vs_vertex0=sp,
                reproduces_stored=bool(k_v0 == st), evals_equal_stored=bool(abs(ev_v0 - st_ev) < 1e-6),
            )  # fmt: skip
            rows.append(
                dict(
                    msin=msin, images=LABEL[var], stored=st, vertex0=f"{k_v0} of {n}",
                    interval=f"{lo:.3f} to {hi:.3f}", centroid=f"{k_cur} of {n}",
                    fixed=sp["fixed"], broken=sp["broken"], both_wrong=sp["both_wrong"],
                    ev_st=st_ev, ev_v0=ev_v0, ev_c=ev_c,
                )
            )  # fmt: skip
    summ["n_voxels"] = len(vox)
    (OUT / "lenient_check_summary.json").write_text(
        json.dumps(summ, indent=1, default=float) + "\n"
    )
    rows2 = [{HEADER[k]: v for k, v in r.items()} for r in rows]
    fmt = {HEADER[k]: ".0f" for k in ("ev_st", "ev_v0", "ev_c")}
    fmt[HEADER["msin"]] = "g"
    tbl = doc_tables.markdown_table(rows2, list(HEADER.values()), fmt)
    doc_tables.write_tables(OUT / "lenient_check_tables.md", {"q8_lenient_check": tbl})
    print(json.dumps({k: v for k, v in summ.items() if k != "n_voxels"}, indent=1, default=float))


if __name__ == "__main__":
    main()
