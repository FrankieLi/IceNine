#!/usr/bin/env python3
"""Summaries of the ray-overlap study steps 1-2 -> summary.json and tables.md (doc_tables).

  uv run python scripts/phase_d/ray_overlap_summary.py

Step 1 (peak level). Per variant (clean, realistic), cost function (screen = pixel radius 3,
local = radius 0), arm and orientation class (truth, near 0.5/1/2/5 deg, random):
  kept  = sum of valid overlapping peaks / sum of detector-0 overlaps (pooled peaks, Wilson
          interval; peaks of one voxel are not independent, so the interval is optimistic)
  pass  = share of evaluations with at least one valid overlapping peak (Wilson over evaluations)
The detector-0 overlap count of an arm is that of the f0_R1 arm in the same cost function (the
detector-0 test is arm independent; the old arm's screen projects vertex 0 instead of the
centroid, so its denominator is approximate).
Step 2 (coarse stage). Per arm: candidates per level, evaluations, whether the truth basin (3 deg)
survives each ranking cut, closest candidate misorientation and rank.
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
TOLS = ("R1", "R2", "R4", "R8", "A")
ARMS = ["old", "cpp"] + [f"f{f}_{t}" for t in TOLS for f in (0, 1)]
CLASSES = ("truth", "near0.5", "near1.0", "near2.0", "near5.0", "random")
CASE_NAMES = ("a", "b", "c", "d", "e", "f_dark", "f_bright")
REF = "f0_R1"


def load(sub: str) -> Dict[str, List[Dict[str, Any]]]:
    out: Dict[str, List[Dict[str, Any]]] = {v: [] for v in VARIANTS}
    merged: Dict[Any, Dict[str, Any]] = {}
    for p in sorted((OUT / sub).glob("*_v*.json")):
        d = json.loads(p.read_text())
        key = (d["variant"], d["voxel"])
        if key in merged and "arms" in d:
            merged[key]["arms"].update(d["arms"])
        else:
            merged[key] = d
    for (var, _), d in merged.items():
        out[var].append(d)
    return out


def peak_summary() -> Dict[str, Any]:
    data = load("peak")
    summ: Dict[str, Any] = {}
    for var in VARIANTS:
        recs = data[var]
        if not recs:
            continue
        for fn in ("screen", "local"):
            for arm in ARMS:
                for cls in CLASSES:
                    kept = ref = pas = n_ev = 0
                    cases = np.zeros(8, dtype=np.int64)
                    for r in recs:
                        sel = np.array([lab == cls for lab in r["labels"]])
                        a = np.array(r[f"{fn}|{arm}"])[sel]
                        b = np.array(r[f"{fn}|{REF}"])[sel]
                        kept += int(a[:, 0].sum())
                        ref += int(b[:, 5 + 7].sum())  # det0_bright of the reference arm
                        pas += int((a[:, 0] > 0).sum())
                        n_ev += int(sel.sum())
                        cases += a[:, 5:13].sum(axis=0).astype(np.int64)
                    lo, hi = S.wilson(
                        min(kept, ref), ref
                    )  # old arm screen: denominator approximate
                    plo, phi = S.wilson(pas, n_ev)
                    summ[f"{var}|{fn}|{arm}|{cls}"] = dict(
                        peaks_kept=kept, det0_overlaps=ref, kept_share=kept / ref if ref else None,
                        kept_lo=lo, kept_hi=hi, n_eval=n_ev, n_pass=pas,
                        pass_share=pas / n_ev if n_ev else None, pass_lo=plo, pass_hi=phi,
                        cases=[int(x) for x in cases[:7]], n_voxels=len(recs),
                    )  # fmt: skip
    return summ


def coarse_summary() -> Dict[str, Any]:
    data = load("coarse")
    summ: Dict[str, Any] = {}
    fails = []
    for var in VARIANTS:
        for r in data[var]:
            for arm, a in r["arms"].items():
                for lv in a["levels"]:
                    if not lv["basin_survives"]:
                        fails.append(
                            dict(
                                variant=var, voxel=r["voxel"], arm=arm, level=lv["level"],
                                n_keep=lv["n_keep"], rank_closest=lv["rank_closest"],
                                n_within_disc=lv["n_within_disc"], n_within_qmc=lv["n_within_qmc"],
                                d_min_disc=lv["d_min_disc"], d_min_qmc=lv["d_min_qmc"],
                            )
                        )  # fmt: skip
    summ["basin_failures"] = fails
    for var in VARIANTS:
        recs = data[var]
        if not recs:
            continue
        for arm in ARMS:
            for lev in range(3):
                rows = []
                for r in recs:
                    lv = r["arms"].get(arm, {}).get("levels", [])
                    if lev < len(lv):
                        rows.append(lv[lev])
                if not rows:
                    continue
                n = len(rows)
                surv = sum(1 for x in rows if x["basin_survives"])
                within = sum(1 for x in rows if x["n_within_disc"] > 0)
                lo, hi = S.wilson(surv, n)
                summ[f"{var}|{arm}|L{lev}"] = dict(
                    n_voxels=n,
                    n_cand_median=float(np.median([x["n_cand"] for x in rows])),
                    n_cand_total=int(sum(x["n_cand"] for x in rows)),
                    n_fz_in_median=float(np.median([x["n_fz"] for x in rows])),
                    n_cliques_median=float(np.median([x["n_cliques"] for x in rows])),
                    n_keep_median=float(np.median([x["n_keep"] for x in rows])),
                    evals_level_median=float(
                        np.median(
                            [
                                x["evals_end"] - (r_prev["evals_end"] if r_prev else 0)
                                for x, r_prev in _with_prev(recs, arm, lev)
                            ]
                        )
                    ),
                    basin_in_candidates=within,
                    basin_survives=surv,
                    closest_kept=sum(1 for x in rows if x["rank_closest"] < x["n_keep"]),
                    survive_lo=lo,
                    survive_hi=hi,
                    d_min_qmc_median=float(np.median([x["d_min_qmc"] for x in rows])),
                    d_min_disc_median=float(np.median([x["d_min_disc"] for x in rows])),
                    rank_closest_median=float(np.median([x["rank_closest"] for x in rows])),
                )
            wall = [r["arms"][arm]["wall_s"] for r in recs if arm in r["arms"]]
            if wall:
                summ[f"{var}|{arm}|wall"] = dict(
                    wall_median_s=float(np.median(wall)), n_voxels=len(wall)
                )
    return summ


def _with_prev(recs: List[Dict[str, Any]], arm: str, lev: int) -> List[Any]:
    out = []
    for r in recs:
        lv = r["arms"].get(arm, {}).get("levels", [])
        if lev < len(lv):
            out.append((lv[lev], lv[lev - 1] if lev > 0 else None))
    return out


def tables(peak: Dict[str, Any], coarse: Dict[str, Any]) -> Dict[str, str]:
    t: Dict[str, str] = {}
    for var in VARIANTS:
        for fn in ("screen", "local"):
            rows = []
            for arm in ARMS:
                row: Dict[str, Any] = dict(arm=arm)
                ok = False
                for cls in ("truth", "near1.0", "near5.0", "random"):
                    s = peak.get(f"{var}|{fn}|{arm}|{cls}")
                    if s is None:
                        continue
                    ok = True
                    row[f"kept {cls}"] = s["kept_share"]
                    row[f"pass {cls}"] = s["pass_share"]
                if ok:
                    rows.append(row)
            if rows:
                cols = ["arm"] + [
                    f"{k} {c}"
                    for c in ("truth", "near1.0", "near5.0", "random")
                    for k in ("kept", "pass")
                ]
                t[f"ray_peak_{var}_{fn}"] = doc_tables.markdown_table(
                    rows, cols, {c: ".3f" for c in cols if c != "arm"}
                )
        rows = []
        for arm in ARMS[2:]:
            row = dict(arm=arm)
            for cls in ("truth", "random"):
                s = peak.get(f"{var}|screen|{arm}|{cls}")
                if s is None:
                    continue
                tot = max(1, sum(s["cases"]))
                for name, c in zip(CASE_NAMES, s["cases"]):
                    row[f"{name} {cls}"] = c / tot
                row[f"events {cls}"] = tot
            if len(row) > 1:
                rows.append(row)
        if rows:
            cols = (
                ["arm"]
                + [f"{n} {c}" for c in ("truth", "random") for n in CASE_NAMES]
                + ["events truth", "events random"]
            )
            fm = {c: ".3f" for c in cols}
            fm.update({"arm": "s", "events truth": "d", "events random": "d"})
            t[f"ray_cases_{var}"] = doc_tables.markdown_table(rows, cols, fm)
        rows = []
        for arm in ARMS:
            for lev in range(3):
                s = coarse.get(f"{var}|{arm}|L{lev}")
                if s is None:
                    continue
                rows.append(
                    dict(
                        arm=arm, level=lev, cand_median=s["n_cand_median"],
                        keep_median=s["n_keep_median"], evals_median=s["evals_level_median"],
                        basin_surv=f"{s['basin_survives']}/{s['n_voxels']}",
                        surv_ci=f"{s['survive_lo']:.2f}-{s['survive_hi']:.2f}",
                        closest_kept=f"{s['closest_kept']}/{s['n_voxels']}",
                        d_min_disc=s["d_min_disc_median"], d_min_qmc=s["d_min_qmc_median"],
                        rank_closest=s["rank_closest_median"],
                    )
                )  # fmt: skip
        if rows:
            cols = list(rows[0].keys())
            t[f"ray_coarse_{var}"] = doc_tables.markdown_table(
                rows, cols, {"cand_median": ".0f", "keep_median": ".0f", "evals_median": ".0f",
                             "d_min_disc": ".3f", "d_min_qmc": ".3f", "rank_closest": ".0f"},
            )  # fmt: skip
    return t


def main() -> None:
    peak = peak_summary()
    coarse = coarse_summary()
    (OUT / "summary.json").write_text(json.dumps(dict(peak=peak, coarse=coarse), indent=1) + "\n")
    doc_tables.write_tables(OUT / "tables.md", tables(peak, coarse))
    print(f"summary.json: {len(peak)} peak rows, {len(coarse)} coarse rows")


if __name__ == "__main__":
    main()
