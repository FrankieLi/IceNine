#!/usr/bin/env python3
"""Step 3 summary: E0 end to end (200 voxels, per-voxel images, min_sin_eta 0).

Arms old, cpp, R2, A.

  uv run python scripts/phase_d/ray_overlap_e2e_summary.py

old = scripts/phase_d/cache/off_ms0 and cpp = on_ms0 (stored earlier runs, same inputs and seeds);
R2, A = cache/rayR2_ms0, rayA_ms0 (e0_cpprule.py run --arm rayR2|rayA --msin 0). Writes
benchmarks/phase_d_seed_diag/ray_overlap/e2e_summary.json and e2e_tables.md.
Paired splits are counted
per voxel: fixed = wrong under the reference but right under the arm, broken = the reverse.
Also the 6 full-sample voxels (runs/ray_*_full_v*.json from seed_diag.py) when present.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
OUT = ICE / "benchmarks" / "phase_d_seed_diag" / "ray_overlap"
sys.path.insert(0, str(ICE / "scripts" / "findoptimal_robustness"))
sys.path.insert(0, str(ICE / "scripts" / "common"))
import common as C  # noqa: E402
import doc_tables  # noqa: E402
import gzio as GZ  # noqa: E402
import stats as S  # noqa: E402

CACHE = HERE / "cache"
DIRS = {"old": "off_ms0", "cpp": "on_ms0", "R2": "rayR2_ms0", "A": "rayA_ms0"}
RUNS = ICE / "benchmarks" / "phase_d_seed_diag" / "runs"
SIX = [15901, 17034, 18558, 19889, 22867, 5242]
# old = lenientnow_full_*: lenient rule re-run 2026-10-12 in the current build (the earlier
# base_full_* runs of 2026-10-08 predate the MC port and the variance-stage fix: old_stale only)
SIX_FILES = {"old": "lenientnow_full_v{}.json", "cpp": "cpprule_full_v{}.json",
             "R2": "rayR2_full_v{}.json", "A": "rayA_full_v{}.json"}  # fmt: skip


def main() -> None:
    summ: Dict[str, Any] = {}
    rows: List[Dict[str, Any]] = []
    pr: List[Dict[str, Any]] = []
    for var in C.VARIANTS:
        label = "realistic" if var == "all" else var
        names = sorted(p.name for p in (CACHE / DIRS["old"]).glob(f"*_{var}.npz"))
        ok = [n for n in names if all((CACHE / d / n).exists() for d in DIRS.values())]
        z = {a: [np.load(CACHE / d / n) for n in ok] for a, d in DIRS.items()}
        wrong = {}
        for a in DIRS:
            Rt = np.stack([x["R_true"] for x in z[a]])
            err = C.err_deg(np.stack([x["R"] for x in z[a]]), Rt)
            wrong[a] = err > C.WRONG_DEG
            n = len(ok)
            k = int(wrong[a].sum())
            lo, hi = S.wilson(k, n)
            ev = np.array([int(x["evals_global"]) + int(x["evals_local"]) for x in z[a]])
            summ[f"{label}|{a}"] = dict(
                n=n, wrong=k, lo=lo, hi=hi, mean_evals=float(ev.mean()),
                mean_runtime_s=float(np.mean([float(x["runtime"]) for x in z[a]])),
                median_err_right=float(np.median(err[~wrong[a]])),
            )  # fmt: skip
            s = summ[f"{label}|{a}"]
            rows.append(
                dict(variant=label, arm=a, wrong=f"{k}/{n} [{lo:.3f}, {hi:.3f}]",
                     mean_evals=s["mean_evals"], runtime_s=s["mean_runtime_s"],
                     median_err_right=s["median_err_right"])
            )  # fmt: skip
        for a in ("cpp", "R2", "A"):
            for ref in ("cpp", "old"):
                if a == ref or (a == "cpp" and ref == "cpp"):
                    continue
                broken, fixed = S.paired_discordant(wrong[a], wrong[ref])
                both = int((wrong[a] & wrong[ref]).sum())
                summ[f"{label}|{a}_vs_{ref}"] = dict(fixed=fixed, broken=broken, both_wrong=both)
                pr.append(dict(variant=label, arm=a, ref=ref, fixed=fixed, broken=broken,
                               both_wrong=both))  # fmt: skip
    tabs = {
        "ray_e2e_e0": doc_tables.markdown_table(
            rows, list(rows[0].keys()),
            {"mean_evals": ".0f", "runtime_s": ".1f", "median_err_right": ".4f"},
        ),
        "ray_e2e_paired": doc_tables.markdown_table(
            pr, ["variant", "arm", "ref", "fixed", "broken", "both_wrong"], {}
        ),
    }  # fmt: skip
    six: List[Dict[str, Any]] = []
    lab: List[Dict[str, Any]] = []
    for v in SIX:
        row: Dict[str, Any] = dict(voxel=v)
        for a, f in SIX_FILES.items():
            p = RUNS / f.format(v)
            if GZ.exists(p):
                d = GZ.load_json(p)
                row[f"{a}_err"], row[f"{a}_evals"], row[f"{a}_wall"] = (
                    d["err_deg"], d["evals_total"], d["wall_s"],
                )  # fmt: skip
        six.append(row)
        for a, f in {**SIX_FILES, "old_stale": "base_full_v{}.json"}.items():
            p = RUNS / f.format(v)
            if GZ.exists(p):
                d = GZ.load_json(p)
                lab.append(
                    dict(voxel=v, arm=a, quiet=bool(d["quiet_preflight"]),
                         load1=d["preflight"]["loadavg"][0],
                         torch_threads=d["preflight"]["torch_threads"],
                         date=str(d["preflight"]["timestamp"])[:10])
                )  # fmt: skip
    if any(len(r) > 1 for r in six):
        cols = ["voxel"] + [f"{a}_{k}" for a in DIRS for k in ("err", "evals", "wall")]
        for r in six:
            for c in cols:
                r.setdefault(c, float("nan"))
        tabs["ray_e2e_six"] = doc_tables.markdown_table(
            six, cols, {c: (".4f" if c.endswith("err") else ".0f") for c in cols if c != "voxel"}
        )
        summ["six"] = six
        summ["six_labels"] = lab
        for a in list(DIRS) + ["old_stale"]:
            rows_a = [x for x in lab if x["arm"] == a]
            summ[f"six_quiet|{a}"] = dict(n=len(rows_a), n_quiet=sum(x["quiet"] for x in rows_a))
        tabs["ray_e2e_six_labels"] = doc_tables.markdown_table(
            lab, ["voxel", "arm", "quiet", "load1", "torch_threads", "date"], {"load1": ".2f"}
        )
        for a in DIRS:
            ev = [r[f"{a}_evals"] for r in six]
            wl = [r[f"{a}_wall"] for r in six]
            summ[f"six|{a}"] = dict(
                mean_evals=float(np.mean(ev)), mean_wall_s=float(np.mean(wl)),
                n_wrong=int(sum(r[f"{a}_err"] > C.WRONG_DEG for r in six)),
            )  # fmt: skip
    (OUT / "e2e_summary.json").write_text(json.dumps(summ, indent=1, default=float) + "\n")
    doc_tables.write_tables(OUT / "e2e_tables.md", tabs)
    print(json.dumps({k: v for k, v in summ.items() if k != "six"}, indent=1, default=float))


if __name__ == "__main__":
    main()
