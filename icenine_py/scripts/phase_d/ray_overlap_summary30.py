#!/usr/bin/env python3
"""Summary of the 30-voxel coarse-stage rerun (arms old, cpp, f0_R2, f0_A; levels 0-2).

  uv run python scripts/phase_d/ray_overlap_summary30.py

Reads benchmarks/phase_d_seed_diag/ray_overlap/coarse30/*.json, writes summary30.json and
tables30.md (doc_tables). Per variant, arm and level, over the voxels:
  surv3 / surv2 / surv1 = voxels with a candidate kept by the ranking cut that lies within 3 / 2 /
                          1 deg (symmetry reduced) of the truth, after the quick MC (Wilson
                          interval)
  in_disc3 / 2 / 1      = voxels whose discrete-stage list (before the quick MC) has a candidate
                          within the threshold
  rank_disc             = rank of the closest discrete-stage candidate by the discrete score, before
                          the quick MC moves it (tie rule: number of candidates with a strictly
                          smaller score); disc_kept = voxels where that rank < n_keep (the cut size)
  post_kept             = voxels where that same candidate, after the MC, is within the kept set
Paired against the cpp arm and the old arm (3 deg survival, per voxel and level): fixed = survives
under the arm but not under the reference, broken = the reverse, both_fail.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
OUT = ICE / "benchmarks" / "phase_d_seed_diag" / "ray_overlap"
sys.path.insert(0, str(ICE / "scripts" / "common"))
import doc_tables  # noqa: E402
import stats as S  # noqa: E402

VARIANTS = ("clean", "realistic")
ARMS = ["old", "cpp", "f0_R2", "f0_A"]


def load() -> Dict[str, Dict[int, Dict[str, Any]]]:
    d: Dict[str, Dict[int, Dict[str, Any]]] = {v: {} for v in VARIANTS}
    for p in sorted((OUT / "coarse30").glob("*.json")):
        r = json.loads(p.read_text())
        d[r["variant"]].setdefault(r["voxel"], {"arms": {}})["arms"].update(r["arms"])
    return d


def main() -> None:
    data = load()
    summ: Dict[str, Any] = {}
    rows_all: Dict[str, List[Dict[str, Any]]] = {}
    for var in VARIANTS:
        # voxels with all four arms
        vox = sorted(v for v, r in data[var].items() if all(a in r["arms"] for a in ARMS))
        summ[f"{var}|voxels"] = vox
        rows: List[Dict[str, Any]] = []
        fail: Dict[str, Dict[int, np.ndarray]] = {}
        for arm in ARMS:
            for lev in range(3):
                lv = [data[var][v]["arms"][arm]["levels"] for v in vox]
                lv = [x[lev] for x in lv if lev < len(x)]
                n = len(lv)
                if not n:
                    continue
                s: Dict[str, Any] = dict(n=n)
                for th, key in ((3, "basin_survives"), (2, "survives_2"), (1, "survives_1")):
                    k = sum(1 for x in lv if x[key])
                    lo, hi = S.wilson(k, n)
                    s[f"surv{th}"], s[f"surv{th}_lo"], s[f"surv{th}_hi"] = k, lo, hi
                s["in_disc3"] = sum(1 for x in lv if x["n_within_disc"] > 0)
                s["in_disc2"] = sum(1 for x in lv if x["n_within_disc_2"] > 0)
                s["in_disc1"] = sum(1 for x in lv if x["n_within_disc_1"] > 0)
                s["rank_disc_median"] = float(np.median([x["rank_closest_disc"] for x in lv]))
                s["rank_disc_max"] = int(max(x["rank_closest_disc"] for x in lv))
                s["n_keep_median"] = float(np.median([x["n_keep"] for x in lv]))
                s["n_cand_median"] = float(np.median([x["n_cand"] for x in lv]))
                s["disc_kept"] = sum(1 for x in lv if x["rank_closest_disc"] < x["n_keep_disc"])
                s["post_kept"] = sum(1 for x in lv if x["rank_post_of_disc_closest"] < x["n_keep"])
                s["d_min_disc_median"] = float(np.median([x["d_min_disc"] for x in lv]))
                s["evals_median"] = float(np.median([x["evals_end"] for x in lv]))
                summ[f"{var}|{arm}|L{lev}"] = s
                fail.setdefault(arm, {})[lev] = np.array([not x["basin_survives"] for x in lv])
                rows.append(
                    dict(
                        arm=arm, level=lev, n=n, cand=s["n_cand_median"], keep=s["n_keep_median"],
                        surv3=f"{s['surv3']}/{n} [{s['surv3_lo']:.2f}, {s['surv3_hi']:.2f}]",
                        surv2=f"{s['surv2']}/{n} [{s['surv2_lo']:.2f}, {s['surv2_hi']:.2f}]",
                        surv1=f"{s['surv1']}/{n} [{s['surv1_lo']:.2f}, {s['surv1_hi']:.2f}]",
                        in_disc3=s["in_disc3"], in_disc1=s["in_disc1"],
                        rank_disc_med=s["rank_disc_median"], disc_kept=s["disc_kept"],
                        post_kept=s["post_kept"],
                    )
                )  # fmt: skip
        for ref in ("cpp", "old"):
            for arm in ARMS:
                if arm == ref:
                    continue
                for lev in range(3):
                    if lev not in fail.get(arm, {}) or lev not in fail.get(ref, {}):
                        continue
                    broken, fixed = S.paired_discordant(fail[arm][lev], fail[ref][lev])
                    both = int((fail[arm][lev] & fail[ref][lev]).sum())
                    summ[f"{var}|{arm}_vs_{ref}|L{lev}"] = dict(
                        fixed=fixed, broken=broken, both_fail=both
                    )
        rows_all[var] = rows
    tabs: Dict[str, str] = {}
    for var, rows in rows_all.items():
        if rows:
            cols = list(rows[0].keys())
            tabs[f"ray30_coarse_{var}"] = doc_tables.markdown_table(
                rows, cols, {"cand": ".0f", "keep": ".0f", "rank_disc_med": ".0f"}
            )
        prow = []
        for k, v in summ.items():
            if k.startswith(var + "|") and "_vs_" in k:
                _, name, lev = k.split("|")
                prow.append(dict(pair=name, level=lev, **v))
        if prow:
            tabs[f"ray30_paired_{var}"] = doc_tables.markdown_table(
                prow, ["pair", "level", "fixed", "broken", "both_fail"], {}
            )
    (OUT / "summary30.json").write_text(json.dumps(summ, indent=1, default=float) + "\n")
    doc_tables.write_tables(OUT / "tables30.md", tabs)
    print({v: len(summ[f"{v}|voxels"]) for v in VARIANTS})


if __name__ == "__main__":
    main()
