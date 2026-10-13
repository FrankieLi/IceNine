#!/usr/bin/env python3
"""Summary of the centroid-start-point task, step 3 (clean images only) -> summary.json, tables.md.

  uv run python scripts/phase_d/centroid_step3_summary.py

Every comparison is "before" (screen projects triangle vertex 0 for the old rule; stored runs)
against "after" (screen projects the voxel centroid for every rule; this task), same voxels and
seeds.
Only the old (lenient) rule changed its start point; cpp, R2 and A already used the centroid,
so their before/after pairs are a check that nothing else moved.

  a. peak level (30 voxels): kept share and pass share per cost function, arm, orientation class.
  b. coarse stage (30 voxels, levels 0-2): 3 deg basin survival, discrete-stage candidates, evals;
     paired fixed / broken / both fail per level (counted by voxel) after against before.
  c. E0 (200 voxels, min_sin_eta 0): wrong count (> C.WRONG_DEG), paired fixed / broken / both wrong
     after against before, mean evaluations.
  d. six full-sample voxels, old rule: error, evaluations, wall (contended) against lenientnow_*.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
OUT = ICE / "benchmarks" / "phase_d_seed_diag" / "centroid"
RAY = ICE / "benchmarks" / "phase_d_seed_diag" / "ray_overlap"
RUNS = ICE / "benchmarks" / "phase_d_seed_diag" / "runs"
CACHE = HERE / "cache"
sys.path.insert(0, str(ICE / "scripts" / "findoptimal_robustness"))
sys.path.insert(0, str(ICE / "scripts" / "common"))
sys.path.insert(0, str(HERE))
import common as C  # noqa: E402
import doc_tables  # noqa: E402
import ray_overlap_summary as RS  # noqa: E402
import stats as S  # noqa: E402

ARMS = ["old", "cpp", "f0_R2", "f0_A"]
E0_DIRS = {
    "old": ("off_ms0", "off_ms0_centroid"),
    "cpp": ("on_ms0", "on_ms0_centroid"),
    "f0_R2": ("rayR2_ms0", "rayR2_ms0_centroid"),
    "f0_A": ("rayA_ms0", "rayA_ms0_centroid"),
}
SIX = [15901, 17034, 18558, 19889, 22867, 5242]


def peak_part() -> Any:
    RS.ARMS = ["old", "cpp", "f0_R1", "f0_R2", "f0_A"]
    both = ("clean", "realistic")  # the stored peak/ dir also holds realistic files; use clean only
    orig = RS.load
    summ: Dict[str, Any] = {}
    for tag, sub in (("before", "peak"), ("after", "peak_centroid")):

        def load_clean(s: str, sub: str = sub) -> Any:
            RS.VARIANTS = both
            out = orig(sub)
            RS.VARIANTS = ("clean",)
            return out

        RS.load = load_clean  # type: ignore[assignment]
        for k, v in RS.peak_summary().items():
            summ[f"{tag}|{k}"] = v
    RS.load = orig
    rows = []
    for fn in ("screen", "local"):
        for cls in ("truth", "near1.0", "random"):
            for arm in ARMS:
                b = summ[f"before|clean|{fn}|{arm}|{cls}"]
                a = summ[f"after|clean|{fn}|{arm}|{cls}"]
                rows.append(
                    dict(cost=fn, cls=cls, arm=arm,
                         kept_before=f"{b['peaks_kept']}/{b['det0_overlaps']}",
                         kept_after=f"{a['peaks_kept']}/{a['det0_overlaps']}",
                         share_before=b["kept_share"], share_after=a["kept_share"],
                         pass_before=f"{b['n_pass']}/{b['n_eval']}",
                         pass_after=f"{a['n_pass']}/{a['n_eval']}")
                )  # fmt: skip
    tab = doc_tables.markdown_table(
        rows, list(rows[0].keys()), {"share_before": ".4f", "share_after": ".4f"}
    )
    return summ, {"centroid_peak_clean": tab}


def load_coarse(sub: str) -> Dict[int, Dict[str, Any]]:
    d: Dict[int, Dict[str, Any]] = {}
    for p in sorted((RAY / sub).glob("clean_v*.json")):
        r = json.loads(p.read_text())
        d.setdefault(r["voxel"], {"arms": {}})["arms"].update(r["arms"])
    return d


def coarse_part() -> Any:
    bef, aft = load_coarse("coarse30"), load_coarse("coarse30_centroid")
    vox = sorted(
        v
        for v in aft
        if all(a in aft[v]["arms"] and a in bef.get(v, {"arms": {}})["arms"] for a in ARMS)
    )
    if not vox:
        raise ValueError("no voxel has all four arms before and after yet")
    summ: Dict[str, Any] = {"voxels": vox}
    rows: List[Dict[str, Any]] = []
    prow: List[Dict[str, Any]] = []
    fail: Dict[str, Dict[str, Dict[int, np.ndarray]]] = {"before": {}, "after": {}}
    for arm in ARMS:
        for lev in range(3):
            for tag, data in (("before", bef), ("after", aft)):
                lv = [data[v]["arms"][arm]["levels"][lev] for v in vox]
                n = len(lv)
                s3 = sum(1 for x in lv if x["basin_survives"])
                lo, hi = S.wilson(s3, n)
                summ[f"{arm}|L{lev}|{tag}"] = dict(
                    n=n, surv3=s3, lo=lo, hi=hi,
                    surv1=sum(1 for x in lv if x["survives_1"]),
                    surv2=sum(1 for x in lv if x["survives_2"]),
                    n_keep_median=float(np.median([x["n_keep"] for x in lv])),
                    rank_disc_median=float(np.median([x["rank_closest_disc"] for x in lv])),
                    in_disc2=sum(1 for x in lv if x["n_within_disc_2"] > 0),
                    in_disc1=sum(1 for x in lv if x["n_within_disc_1"] > 0),
                    in_disc3=sum(1 for x in lv if x["n_within_disc"] > 0),
                    n_cand_median=float(np.median([x["n_cand"] for x in lv])),
                    evals_end_median=float(np.median([x["evals_end"] for x in lv])),
                    evals_end_mean=float(np.mean([x["evals_end"] for x in lv])),
                )  # fmt: skip
                fail[tag].setdefault(arm, {})[lev] = np.array([not x["basin_survives"] for x in lv])
                rows.append(
                    dict(arm=arm, level=lev, when=tag, surv3=f"{s3}/{n} [{lo:.2f}, {hi:.2f}]",
                         surv2=summ[f"{arm}|L{lev}|{tag}"]["surv2"],
                         surv1=summ[f"{arm}|L{lev}|{tag}"]["surv1"],
                         keep_median=summ[f"{arm}|L{lev}|{tag}"]["n_keep_median"],
                         rank_disc_median=summ[f"{arm}|L{lev}|{tag}"]["rank_disc_median"],
                         in_disc3=summ[f"{arm}|L{lev}|{tag}"]["in_disc3"],
                         cand_median=summ[f"{arm}|L{lev}|{tag}"]["n_cand_median"],
                         evals_median=summ[f"{arm}|L{lev}|{tag}"]["evals_end_median"])
                )  # fmt: skip
            fb, fa = fail["before"][arm][lev], fail["after"][arm][lev]
            broken, fixed = S.paired_discordant(fa, fb)
            summ[f"{arm}|L{lev}|paired"] = dict(
                fixed=fixed, broken=broken, both_fail=int((fa & fb).sum())
            )
            prow.append(dict(arm=arm, level=lev, fixed=fixed, broken=broken,
                             both_fail=int((fa & fb).sum())))  # fmt: skip
    # old after against cpp after (the lenient-versus-C++ gap with equal start points)
    for lev in range(3):
        fo, fc = fail["after"]["old"][lev], fail["after"]["cpp"][lev]
        broken, fixed = S.paired_discordant(fo, fc)
        summ[f"old_vs_cpp|L{lev}|after"] = dict(
            fixed=fixed, broken=broken, both_fail=int((fo & fc).sum())
        )
    tabs = {
        "centroid_coarse_clean": doc_tables.markdown_table(
            rows, list(rows[0].keys()), {"cand_median": ".0f", "evals_median": ".0f"}
        ),
        "centroid_coarse_paired": doc_tables.markdown_table(
            prow, ["arm", "level", "fixed", "broken", "both_fail"], {}
        ),
    }
    return summ, tabs


def e0_part() -> Any:
    summ: Dict[str, Any] = {}
    rows: List[Dict[str, Any]] = []
    prow: List[Dict[str, Any]] = []
    wrong: Dict[str, Dict[str, np.ndarray]] = {"before": {}, "after": {}}
    for arm, (db, da) in E0_DIRS.items():
        names = sorted(p.name for p in (CACHE / da).glob("*_clean.npz"))
        names = [n for n in names if (CACHE / db / n).exists()]
        n = len(names)
        for tag, d in (("before", db), ("after", da)):
            z = [np.load(CACHE / d / nm) for nm in names]
            Rt = np.stack([x["R_true"] for x in z])
            err = C.err_deg(np.stack([x["R"] for x in z]), Rt)
            w = err > C.WRONG_DEG
            wrong[tag][arm] = w
            ev = np.array([int(x["evals_global"]) + int(x["evals_local"]) for x in z])
            lo, hi = S.wilson(int(w.sum()), n)
            summ[f"{arm}|{tag}"] = dict(
                n=n, wrong=int(w.sum()), lo=lo, hi=hi, mean_evals=float(ev.mean()),
                median_err_right=float(np.median(err[~w])),
                mean_runtime_s=float(np.mean([float(x["runtime"]) for x in z])),
            )  # fmt: skip
            s = summ[f"{arm}|{tag}"]
            rows.append(
                dict(arm=arm, when=tag, wrong=f"{s['wrong']}/{n} [{lo:.3f}, {hi:.3f}]",
                     mean_evals=s["mean_evals"], median_err_right=s["median_err_right"])
            )  # fmt: skip
        broken, fixed = S.paired_discordant(wrong["after"][arm], wrong["before"][arm])
        both = int((wrong["after"][arm] & wrong["before"][arm]).sum())
        summ[f"{arm}|paired"] = dict(fixed=fixed, broken=broken, both_wrong=both)
        prow.append(dict(arm=arm, fixed=fixed, broken=broken, both_wrong=both))
    for ref in ("cpp", "f0_R2", "f0_A"):
        broken, fixed = S.paired_discordant(wrong["after"]["old"], wrong["after"][ref])
        both = int((wrong["after"]["old"] & wrong["after"][ref]).sum())
        summ[f"old_vs_{ref}|after"] = dict(fixed=fixed, broken=broken, both_wrong=both)
        # same split, old rule as reference (fixed = wrong under old, right under ref)
        b2, f2 = S.paired_discordant(wrong["after"][ref], wrong["after"]["old"])
        summ[f"{ref}_vs_old|after"] = dict(fixed=f2, broken=b2, both_wrong=both)
        prow.append(dict(arm=f"{ref} vs old (after)", fixed=f2, broken=b2, both_wrong=both))
    allright = ~(
        wrong["after"]["old"]
        | wrong["after"]["cpp"]
        | wrong["after"]["f0_R2"]
        | wrong["after"]["f0_A"]
    )
    summ["right_under_all_after"] = dict(n=int(allright.sum()))
    for arm, (db, da) in E0_DIRS.items():
        names = sorted(p.name for p in (CACHE / da).glob("*_clean.npz"))
        names = [n for n in names if (CACHE / db / n).exists()]
        for tag, d in (("before", db), ("after", da)):
            z = [np.load(CACHE / d / nm) for nm in names]
            err = C.err_deg(np.stack([x["R"] for x in z]), np.stack([x["R_true"] for x in z]))
            summ[f"{arm}|{tag}"]["median_err_right_all_after"] = float(np.median(err[allright]))
    tabs = {
        "centroid_e0_clean": doc_tables.markdown_table(
            rows, list(rows[0].keys()), {"mean_evals": ".0f", "median_err_right": ".4f"}
        ),
        "centroid_e0_paired": doc_tables.markdown_table(
            prow, ["arm", "fixed", "broken", "both_wrong"], {}
        ),
    }
    return summ, tabs


def six_part() -> Any:
    rows: List[Dict[str, Any]] = []
    for v in SIX:
        row: Dict[str, Any] = dict(voxel=v)
        for tag, f in (
            ("before", f"lenientnow_full_v{v}.json"),
            ("after", f"base_centroid_full_v{v}.json"),
        ):
            p = RUNS / f
            if p.exists():
                d = json.loads(p.read_text())
                row[f"err_{tag}"], row[f"evals_{tag}"], row[f"wall_{tag}"] = (
                    d["err_deg"], d["evals_total"], d["wall_s"],
                )  # fmt: skip
                row[f"label_{tag}"] = d["label"]
        rows.append(row)
    cols = [
        "voxel",
        "err_before",
        "err_after",
        "evals_before",
        "evals_after",
        "wall_before",
        "wall_after",
    ]
    for r in rows:
        for c in cols:
            r.setdefault(c, float("nan"))
    summ: Dict[str, Any] = dict(rows=rows)
    ok = [r for r in rows if r["err_after"] == r["err_after"]]
    if ok:
        summ["n"] = len(ok)
        summ["wrong_before"] = int(sum(r["err_before"] > C.WRONG_DEG for r in ok))
        summ["wrong_after"] = int(sum(r["err_after"] > C.WRONG_DEG for r in ok))
        summ["mean_evals_before"] = float(np.mean([r["evals_before"] for r in ok]))
        summ["mean_evals_after"] = float(np.mean([r["evals_after"] for r in ok]))
    for arm, f in (("cpp", "cpprule_full_v{}.json"), ("f0_R2", "rayR2_full_v{}.json"),
                   ("f0_A", "rayA_full_v{}.json")):  # fmt: skip
        ds = [json.loads((RUNS / f.format(v)).read_text()) for v in SIX]
        summ[f"mean_evals_{arm}"] = float(np.mean([d["evals_total"] for d in ds]))
        summ[f"mean_wall_{arm}"] = float(np.mean([d["wall_s"] for d in ds]))
    if ok:
        summ["mean_wall_before"] = float(np.mean([r["wall_before"] for r in ok]))
        summ["mean_wall_after"] = float(np.mean([r["wall_after"] for r in ok]))
        for arm in ("f0_R2", "f0_A"):
            summ[f"evals_ratio_{arm}_vs_cpp"] = summ[f"mean_evals_{arm}"] / summ["mean_evals_cpp"]
        for tag in ("before", "after"):
            summ[f"evals_ratio_old_{tag}_vs_cpp"] = (
                summ[f"mean_evals_{tag}"] / summ["mean_evals_cpp"]
            )
            for arm in ("f0_R2", "f0_A"):
                summ[f"evals_ratio_old_{tag}_vs_{arm}"] = (
                    summ[f"mean_evals_{tag}"] / summ[f"mean_evals_{arm}"]
                )
    tab = doc_tables.markdown_table(
        rows, cols, {c: (".4f" if c.startswith("err") else ".0f") for c in cols if c != "voxel"}
    )
    return summ, {"centroid_six_old": tab}


def main() -> None:
    summ: Dict[str, Any] = {}
    tabs: Dict[str, str] = {}
    for name, fn in (
        ("peak", peak_part),
        ("coarse", coarse_part),
        ("e0", e0_part),
        ("six", six_part),
    ):
        try:
            s, t = fn()
        except (FileNotFoundError, KeyError, ValueError) as e:  # a part that has not finished
            print(f"{name}: skipped ({type(e).__name__}: {e})")
            continue
        summ[name] = s
        tabs.update(t)
    (OUT / "summary.json").write_text(json.dumps(summ, indent=1, default=float) + "\n")
    doc_tables.write_tables(OUT / "tables.md", tabs)
    print("parts:", list(summ))


if __name__ == "__main__":
    main()
