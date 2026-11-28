"""tex1: score the two whole-layer arms against the ``tex24`` truth (light, one thread).

  uv run python scripts/autoresearch/tex1/tex1_score.py [--arms CSP,RCP]

Run it only after both end markers exist (``<cache>/runs/whole/{CSP,RCP}.done``): nothing of
``tex24`` is read before that. Per arm, read the finished map (``tex1_cpp.mic.LBFS`` of ``CSP``,
``recon.mic.gz`` of ``RCP``), match every output voxel to the grid row of the same direction by left
vertex (tolerance 1e-3 mm, exactly one row) and take the cubic-symmetry misorientation to that
voxel's OWN truth. Wrong: over 1 degree; a grid row with no output, or a zero orientation, counts
as wrong and is listed apart ("no output"). Wilson 95% interval (``common/stats.py``).

Per class (sharp, noise, gradient, both) and over all voxels: n, wrong, no output, median and 90th
percentile of the error, the median error of the trivial answer (``R_g`` for every voxel of the
grain) and the share of voxels the arm is closer to the truth than the trivial answer is. In
noise grains the trivial answer is never wrong (a voxel is within 1 degree of the smooth field),
so the 1-degree rule finds gross errors only. Per gradient grain: the recovered spread (largest
and 90th-percentile misorientation to ``R_g`` of the voxels with output) against the true one.
Paired by voxel: both wrong, only one. Evaluations: ``RCP`` from ``run.json``; ``CSP`` the sum of
the ``[cost evaluations] client N: X`` lines of the build with the count line, copied to the
tracked ``cost_evaluation_lines.txt`` at the first call. Wall: ``/usr/bin/time -p`` real, the
launcher's ``<cache>/logs/whole_<arm>.time``. No p-values.

Writes ``benchmarks/autoresearch/tex1/score.json`` and ``score_tables.md``; per-voxel errors go to
``<cache>/score/`` (not tracked).
"""

import argparse
import gzip
import json
import re
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "common"))
sys.path.insert(0, str(HERE.parents[0] / "cpptime"))
import doc_tables as DT  # noqa: E402
import parse_logs as PL  # noqa: E402
import stats as ST  # noqa: E402
import tex1_lib as F  # noqa: E402
import tex1_paths as P  # noqa: E402

ARMS = ("CSP", "RCP")
WORKERS = {"CSP": 12, "RCP": 1}
WRONG_DEG = 1.0
TOL_MM = 1e-3
RUNS = P.RUNS / "whole"
EVAL_RECORD = P.BENCH / "cost_evaluation_lines.txt"
EVAL_LINE = re.compile(r"\[cost evaluations\] client (\d+):\s*(\d+)")
CLASSES = ("all",) + F.CLASS_NAMES
RECOVERED_MEDIAN_DEG = 0.2  # a noise grain is "recovered" if its median error is under this
SPREAD_BAND = (0.8, 1.25)  # recovered / true spread counted as close inside this band


def read_map(arm: str) -> Any:
    from icenine.mic_file import MicFile

    run = RUNS / arm
    if arm == "CSP":
        return MicFile.read(str(run / "tex1_cpp.mic.LBFS"))
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td) / "recon.mic"
        with gzip.open(run / "recon.mic.gz", "rb") as fi, open(tmp, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        return MicFile.read(str(tmp))


def grid_rows() -> Any:
    a = np.loadtxt(str(P.GRID), skiprows=1, ndmin=2)
    return a[:, :2], a[:, 3].astype(int)


def match(
    out_xy: np.ndarray, out_dir: np.ndarray, g_xy: np.ndarray, g_dir: np.ndarray
) -> np.ndarray:
    """Grid row of every output voxel (exactly one candidate of the same direction)."""
    idx = np.full(len(out_xy), -1, dtype=np.int64)
    for code in (1, 2):
        rows = np.nonzero(g_dir == code)[0]
        sel = np.nonzero(out_dir == code)[0]
        if len(sel) == 0:
            continue
        d, j = cKDTree(g_xy[rows]).query(out_xy[sel], k=1)
        assert float(d.max()) < TOL_MM, f"an output voxel is {d.max()} mm from every grid row"
        idx[sel] = rows[j]
    assert (idx >= 0).all() and len(set(idx.tolist())) == len(idx), "match is not one to one"
    return idx


def class_row(
    err: np.ndarray, has: np.ndarray, triv: np.ndarray, sel: np.ndarray
) -> Dict[str, Any]:
    """The numbers of one voxel subset ``sel`` (bool, n): ``err`` misorientation to the voxel's
    truth (inf where no output), ``has`` output exists, ``triv`` error of the trivial answer."""
    n = int(sel.sum())
    out = sel & has
    wrong = sel & (~has | (err > WRONG_DEG))
    lo, hi = ST.wilson(int(wrong.sum()), n) if n else (float("nan"), float("nan"))
    e = err[out]
    t = triv[out]
    closer = int((e < t).sum())
    clo, chi = ST.wilson(closer, int(out.sum())) if out.any() else (float("nan"), float("nan"))
    return dict(
        n=n,
        n_no_output=int((sel & ~has).sum()),
        n_wrong=int(wrong.sum()),
        wrong_fraction=float(wrong.sum() / n) if n else float("nan"),
        wilson95=[lo, hi],
        median_err_deg=float(np.median(e)) if len(e) else float("nan"),
        p90_err_deg=float(np.percentile(e, 90)) if len(e) else float("nan"),
        trivial_median_err_deg=float(np.median(t)) if len(t) else float("nan"),
        trivial_wrong=int((sel & (triv > WRONG_DEG)).sum()),
        n_closer_than_trivial=closer,
        closer_fraction=float(closer / out.sum()) if out.any() else float("nan"),
        closer_wilson95=[clo, chi],
    )


def class_rows(
    err: np.ndarray,
    has: np.ndarray,
    triv: np.ndarray,
    vox_class: np.ndarray,
    subset: np.ndarray,
) -> Dict[str, Any]:
    """``class_row`` over all voxels of ``subset`` and over each class inside it."""
    out = {"all": class_row(err, has, triv, subset)}
    for c, name in enumerate(F.CLASS_NAMES):
        out[name] = class_row(err, has, triv, subset & (vox_class == c))
    return out


def spread_rows(
    R_out: np.ndarray,
    has: np.ndarray,
    ids: np.ndarray,
    R_g: Dict[int, np.ndarray],
    truth_R: np.ndarray,
    grains: List[int],
    sym: Any,
) -> List[Dict[str, Any]]:
    """Per gradient grain: the true spread (largest and 90th-percentile misorientation of the
    truth to ``R_g``) and the recovered one (the same of the voxels with output)."""
    from icenine.ray_indexing import misorientation_deg

    rows = []
    for g in grains:
        sel = np.nonzero(ids == g)[0]
        Rg = np.broadcast_to(R_g[g], (len(sel), 3, 3))
        t = np.asarray(misorientation_deg(truth_R[sel], Rg, sym), dtype=float)
        h = sel[has[sel]]
        r = np.asarray(
            misorientation_deg(R_out[h], np.broadcast_to(R_g[g], (len(h), 3, 3)), sym), dtype=float
        )
        rows.append(
            dict(
                id=int(g),
                n_voxels=int(len(sel)),
                n_with_output=int(len(h)),
                true_max_deg=float(t.max()),
                true_p90_deg=float(np.percentile(t, 90)),
                rec_max_deg=float(r.max()) if len(r) else float("nan"),
                rec_p90_deg=float(np.percentile(r, 90)) if len(r) else float("nan"),
            )
        )
    return rows


def spread_summary(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Median and quartiles over grains of recovered / true 90th-percentile spread."""
    q = np.array([r["rec_p90_deg"] / r["true_p90_deg"] for r in rows if r["true_p90_deg"] > 0])
    q = q[np.isfinite(q)]
    if not len(q):
        return dict(n=0)
    return dict(
        n=int(len(q)),
        n_within_band=int(((q >= SPREAD_BAND[0]) & (q <= SPREAD_BAND[1])).sum()),
        band=list(SPREAD_BAND),
        ratio_p90_median=float(np.median(q)),
        ratio_p90_q25=float(np.percentile(q, 25)),
        ratio_p90_q75=float(np.percentile(q, 75)),
    )


def noise_grain_summary(
    err: np.ndarray,
    ids: np.ndarray,
    cls_of: Dict[int, int],
    grains: List[int],
) -> Dict[str, Any]:
    """Per grain with noise (noise only or both): its median error; a grain is "recovered" if the
    median is under ``RECOVERED_MEDIAN_DEG``. Summary per class: grains, recovered with a Wilson
    interval, the pooled median and the wrong count and share (over 1 degree) of the voxels of the
    recovered grains, and the grain size ranges of the recovered and the other grains."""
    rows = []
    for g in grains:
        e = err[ids == g]
        rows.append(
            dict(
                id=int(g),
                cls=F.CLASS_NAMES[cls_of[g]],
                n_voxels=int(len(e)),
                median_err_deg=float(np.median(e)),
                n_wrong=int((e > WRONG_DEG).sum()),
            )
        )
    out: Dict[str, Any] = dict(
        rule=f"recovered = grain median error under {RECOVERED_MEDIAN_DEG} deg"
    )
    for c in (F.NOISE_ONLY, F.BOTH):
        name = F.CLASS_NAMES[c]
        sub = [r for r in rows if r["cls"] == name]
        rec = [r for r in sub if r["median_err_deg"] < RECOVERED_MEDIAN_DEG]
        oth = [r for r in sub if r["median_err_deg"] >= RECOVERED_MEDIAN_DEG]
        m = np.isin(ids, [r["id"] for r in rec])
        e = err[m]
        lo, hi = ST.wilson(len(rec), len(sub)) if sub else (float("nan"), float("nan"))
        out[name] = dict(
            n_grains=len(sub),
            n_recovered=len(rec),
            recovered_wilson95=[lo, hi],
            recovered_pooled_median_err_deg=float(np.median(e)) if len(e) else float("nan"),
            recovered_n_voxels=int(m.sum()),
            recovered_n_wrong=int((e > WRONG_DEG).sum()),
            recovered_wrong_fraction=float((e > WRONG_DEG).mean()) if len(e) else float("nan"),
            size_range_recovered=(
                [min(r["n_voxels"] for r in rec), max(r["n_voxels"] for r in rec)] if rec else None
            ),
            size_range_other=(
                [min(r["n_voxels"] for r in oth), max(r["n_voxels"] for r in oth)] if oth else None
            ),
            other_median_range_deg=(
                [
                    min(r["median_err_deg"] for r in oth),
                    max(r["median_err_deg"] for r in oth),
                ]
                if oth
                else None
            ),
        )
    out["grains"] = rows
    return out


def pair_counts(wa: np.ndarray, wb: np.ndarray) -> Dict[str, int]:
    """Counts by voxel of two arms' wrong flags (a = CSP, b = RCP): both wrong, only a wrong
    (fixed by b), only b wrong (broken by b), neither."""
    return dict(
        both_wrong=int((wa & wb).sum()),
        only_a_wrong=int((wa & ~wb).sum()),
        only_b_wrong=int((~wa & wb).sum()),
        neither=int((~wa & ~wb).sum()),
    )


def eval_lines_from_log(text: str) -> Dict[str, Any]:
    vals = [(int(m.group(1)), int(m.group(2))) for m in EVAL_LINE.finditer(text)]
    return dict(n_lines=len(vals), total=int(sum(v for _, v in vals)))


def csp_evaluations() -> Dict[str, Any]:
    """The 12 client lines, copied to the tracked record at the first call (the launcher's log
    can be overwritten later); the sum is read from the record."""
    if not EVAL_RECORD.exists():
        log = (P.LOGS / "whole_CSP.log").read_text()
        lines = [ln for ln in log.splitlines() if EVAL_LINE.search(ln)]
        EVAL_RECORD.write_text("# CSP\n" + "\n".join(lines) + "\n")
    return eval_lines_from_log(EVAL_RECORD.read_text())


def time_seconds(arm: str) -> Dict[str, float]:
    t = PL.parse_time_file((P.LOGS / f"whole_{arm}.time").read_text())
    return dict(real=t["real"], user=t["user"], sys=t["sys"])


def marker(arm: str) -> str:
    return (RUNS / f"{arm}.done").read_text().strip()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--arms", default=",".join(ARMS))
    a = ap.parse_args()
    arms = a.arms.split(",")
    for arm in ARMS:  # nothing of tex24 is read before both end markers exist
        if not (RUNS / f"{arm}.done").exists():
            raise SystemExit(f"{RUNS / (arm + '.done')} is missing: score after both runs end")
    from icenine.mic_file import MicFile
    from icenine.ray_indexing import misorientation_deg, proper_rotations
    from icenine.symmetry import create_cubic_symmetry

    sym = proper_rotations(create_cubic_symmetry(3.6))
    truth = MicFile.read(str(P.TRUTH))
    truth_R = np.array([v.orientation for v in truth.voxels], dtype=np.float64)
    info = json.loads(P.SAMPLE_JSON.read_text())
    ids = np.load(str(P.GRAINMAP))["grain"]
    cls_of = {int(g["id"]): int(g["cls"]) for g in info["grains"]}
    R_g = {int(g["id"]): np.array(g["R_centre"]) for g in info["grains"]}
    vox_class = np.array([cls_of[int(g)] for g in ids])
    Rg_v = np.stack([R_g[int(g)] for g in ids])
    triv = np.asarray(misorientation_deg(Rg_v, truth_R, sym), dtype=float)
    P.SCORE.mkdir(parents=True, exist_ok=True)
    np.save(P.SCORE / "trivial_err.npy", triv)
    np.save(P.SCORE / "voxel_class.npy", vox_class)
    g_xy, g_dir = grid_rows()
    n = len(g_xy)
    assert len(truth_R) == n == P.N_FINE
    P.SCORE.mkdir(parents=True, exist_ok=True)
    grad_grains = [g for g, c in cls_of.items() if c in (F.GRADIENT_ONLY, F.BOTH)]
    res: Dict[str, Any] = {}
    maps: Dict[str, Any] = {}
    for arm in arms:
        mic = read_map(arm)
        xy = np.array([v.position[:2] for v in mic.voxels], dtype=np.float64)
        up = np.array([1 if v.points_up else 2 for v in mic.voxels])
        R = np.array([v.orientation for v in mic.voxels], dtype=np.float64)
        idx = match(xy, up, g_xy, g_dir)
        R_all = np.zeros((n, 3, 3))
        has = np.zeros(n, dtype=bool)
        R_all[idx] = R
        has[idx] = np.abs(np.linalg.det(R)) > 1e-6
        err = np.full(n, np.inf)
        err[has] = np.asarray(misorientation_deg(R_all[has], truth_R[has], sym), dtype=float)
        maps[arm] = dict(has=has, wrong=~has | (err > WRONG_DEG), err=err)
        np.save(P.SCORE / f"{arm}_err.npy", err)
        sp = spread_rows(R_all, has, ids, R_g, truth_R, grad_grains, sym)
        noise_grains = [g for g, c in cls_of.items() if c in (F.NOISE_ONLY, F.BOTH)]
        t = time_seconds(arm)
        d: Dict[str, Any] = dict(
            workers=WORKERS[arm],
            marker=marker(arm),
            n_output_voxels=int(len(idx)),
            classes=class_rows(err, has, triv, vox_class, np.ones(n, dtype=bool)),
            spread_rows=sp,
            spread_summary=spread_summary(sp),
            noise_grains=noise_grain_summary(err, ids, cls_of, noise_grains),
            wall_s=t["real"],
            user_s=t["user"],
            sys_s=t["sys"],
        )
        if arm == "CSP":
            ev = csp_evaluations()
            d.update(evals_all=ev["total"], evals_lines=ev["n_lines"])
        else:
            rj = json.loads((RUNS / arm / "run.json").read_text())
            d.update(
                evals_all=rj["evals_all"],
                evals_global=rj["evals_global"],
                evals_by_key=rj["evals_by_key"],
                global_searches=rj["global_searches"],
                wall_inprocess_s=rj["wall_total_s"],
                n_without_orientation=rj["n_without_orientation"],
                index_s=rj["stats"]["index_s"],
                index_pairs=rj["stats"]["index_pairs"],
                index_votes=rj["stats"]["index_votes"],
            )
        res[arm] = d
        print(arm, d["classes"]["all"]["n_wrong"], d["evals_all"], flush=True)
    common = np.ones(n, dtype=bool)
    for k in maps:
        common &= maps[k]["has"]
    for k, d in res.items():
        d["common_n"] = int(common.sum())
        d["common_classes"] = class_rows(
            maps[k]["err"], maps[k]["has"], triv, vox_class, common.copy()
        )
    pairs: Dict[str, Any] = {}
    if "CSP" in maps and "RCP" in maps:
        a_, b_ = maps["CSP"]["wrong"], maps["RCP"]["wrong"]
        pairs["all_voxels"] = pair_counts(a_, b_)
        for c, name in enumerate(F.CLASS_NAMES):
            m = vox_class == c
            pairs[name] = pair_counts(a_[m], b_[m])
        pairs["note"] = "a = CSP, b = RCP; counts by voxel; no output counts as wrong"
    out = dict(
        label="contended and oversubscribed: see the markers; CSP ran 13 MPI ranks on 12 cores "
        "(cpptools holds a core). Wall seconds are /usr/bin/time -p real for both arms.",
        n_voxels=n,
        trivial_answer="R_g of the voxel's grain for every voxel",
        arms=res,
        pairs=pairs,
    )
    (P.BENCH / "score.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    write_tables(res, pairs)
    return 0


def fmt_share(d: Dict[str, Any]) -> str:
    lo, hi = d["wilson95"]
    return f"{100 * d['wrong_fraction']:.3f}% ({100 * lo:.3f}, {100 * hi:.3f})"


def write_tables(res: Dict[str, Any], pairs: Dict[str, Any]) -> None:
    """``score_tables.md``: marker-delimited tables (``doc_tables``) the docs copy from."""
    t: Dict[str, str] = {}
    rows = []
    for arm, d in res.items():
        for c in CLASSES:
            r = d["classes"][c]
            rows.append(
                dict(
                    arm=arm,
                    cls=c,
                    voxels=r["n"],
                    wrong=r["n_wrong"],
                    share=fmt_share(r),
                    no_output=r["n_no_output"],
                    median=r["median_err_deg"],
                    p90=r["p90_err_deg"],
                    trivial_median=r["trivial_median_err_deg"],
                    trivial_wrong=r["trivial_wrong"],
                    closer=(
                        "n/a (pools sharp voxels)"
                        if c == "all"
                        else f"{100 * r['closer_fraction']:.1f}%"
                    ),
                )
            )
    t["tex1_per_class"] = DT.markdown_table(
        rows,
        [
            "arm",
            "cls",
            "voxels",
            "wrong",
            "share",
            "no_output",
            "median",
            "p90",
            "trivial_median",
            "trivial_wrong",
            "closer",
        ],
        formats={
            "voxels": ",d",
            "trivial_wrong": ",d",
            "wrong": ",d",
            "no_output": ",d",
            "median": ".4f",
            "p90": ".4f",
            "trivial_median": ".4f",
        },
    )
    t["tex1_cost_wall"] = DT.markdown_table(
        [
            dict(
                arm=k,
                workers="12 workers, 13 ranks on 12 cores" if d["workers"] > 1 else "1 worker",
                evals=d["evals_all"],
                wall_s=d["wall_s"],
                wall_h=f"{d['wall_s'] / 3600:.2f}",
            )
            for k, d in res.items()
        ],
        ["arm", "workers", "evals", "wall_s", "wall_h"],
        formats={"evals": ",d", "wall_s": ",.0f"},
    )
    if pairs:
        t["tex1_pairs"] = DT.markdown_table(
            [
                dict(subset=k, **{kk: v for kk, v in pairs[k].items()})
                for k in ("all_voxels",) + F.CLASS_NAMES
            ],
            ["subset", "both_wrong", "only_a_wrong", "only_b_wrong", "neither"],
            formats={k: ",d" for k in ("both_wrong", "only_a_wrong", "only_b_wrong", "neither")},
        )
    sp = []
    for k, d in res.items():
        s = d["spread_summary"]
        if s.get("n"):
            sp.append(
                dict(
                    arm=k,
                    grains=s["n"],
                    median=s["ratio_p90_median"],
                    q25=s["ratio_p90_q25"],
                    q75=s["ratio_p90_q75"],
                    within=f"{s['n_within_band']} of {s['n']}",
                )
            )
    if sp:
        t["tex1_spread_ratio"] = DT.markdown_table(
            sp,
            ["arm", "grains", "median", "q25", "q75", "within"],
            formats={"median": ".3f", "q25": ".3f", "q75": ".3f"},
        )
    ng = []
    for k, d in res.items():
        for name in ("noise", "both"):
            r = d["noise_grains"][name]
            lo, hi = r["recovered_wilson95"]
            ng.append(
                dict(
                    arm=k,
                    cls=name,
                    grains=r["n_grains"],
                    recovered=f"{r['n_recovered']} ({100 * lo:.0f}%, {100 * hi:.0f}%)",
                    pooled_median=r["recovered_pooled_median_err_deg"],
                    wrong_in_recovered=(
                        f"{r['recovered_n_wrong']} of {r['recovered_n_voxels']} "
                        f"({100 * r['recovered_wrong_fraction']:.1f}%)"
                        if r["n_recovered"]
                        else "none recovered"
                    ),
                    size_recovered=str(r["size_range_recovered"] or "none"),
                    size_other=str(r["size_range_other"] or "none"),
                )
            )
    t["tex1_noise_grains"] = DT.markdown_table(
        ng,
        [
            "arm",
            "cls",
            "grains",
            "recovered",
            "pooled_median",
            "wrong_in_recovered",
            "size_recovered",
            "size_other",
        ],
        formats={"pooled_median": ".3f"},
    )
    DT.write_tables(P.BENCH / "score_tables.md", t)
    print((P.BENCH / "score_tables.md").read_text())


if __name__ == "__main__":
    sys.exit(main())
