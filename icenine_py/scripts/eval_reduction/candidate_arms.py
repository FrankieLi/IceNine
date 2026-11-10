"""Task C step 3, E0 part: the de-duplication arms (config key CandidateDedup) on the E0 set (200
voxels, per-voxel clean images, seed 0, ``mc`` arm, min_sin_eta 0, default lenient rule, default
batched screen), paired per voxel against the baseline arm (key off).

  uv run python scripts/eval_reduction/candidate_arms.py run --arm base --spec off --workers 10
  uv run python scripts/eval_reduction/candidate_arms.py run --arm c1_t05 --spec pre=0.5
  uv run python scripts/eval_reduction/candidate_arms.py summary

One small json per voxel and arm under ``benchmarks/eval_reduction/C/cache/e0/<arm>/`` (gitignored):
error, evaluations in total and per level (screen + re-cost, quick MC), candidates entering and
leaving each stage, and the smallest symmetry-reduced distance from the truth to the candidates
at each stage (discrete list, quick-MC inputs, post-MC list, kept quarter) for the truth-basin
presence per level. Wrong = error above 1 deg (``common.WRONG_DEG``). The wall time is contended
multi-worker and is not a timing claim.
"""

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICE / "scripts" / "findoptimal_robustness"))
sys.path.insert(0, str(ICE / "scripts" / "common"))
sys.path.insert(0, str(ICE / "scripts" / "phase_d"))

BASE_DIR = ICE / "benchmarks" / "eval_reduction" / "C"
CACHE = BASE_DIR / "cache" / "e0"
OUT = BASE_DIR / "step3"
BASIN = (1.0, 3.0)
# The stored E0 runs of the lenient rule before and after the centroid screen start point
# (scripts/phase_d/cache, gitignored; the centroid study's table centroid_e0_paired) and the
# uninstrumented off arm on today's code: the baseline arm is compared with them run by run.
PHASE_D_CACHE = ICE / "scripts" / "phase_d" / "cache"
REFERENCE_RUNS = {
    "off_ms0": "before the centroid screen start point (2026-10-12) and BatchedScreen",
    "off_ms0_centroid": "centroid start point, per-orientation screen (before BatchedScreen)",
    "off_ms0_c_ref": "e0_cpprule.py off arm on the Task C code, uninstrumented",
}


def level_record(L: Dict[str, Any], Rt: np.ndarray) -> Dict[str, Any]:
    from csl import reduced_misorientation_deg

    def mind(R: np.ndarray) -> float:
        if len(R) == 0:
            return float("inf")
        return float(reduced_misorientation_deg(R, np.broadcast_to(Rt, R.shape).copy()).min())

    Ra = L["R_after_ranked"]
    nk = int(L["n_keep"])
    return dict(
        level=int(L["level"]),
        n_disc=int(L["n_disc"]),
        n_mc=int(L["n_mc"]),
        n_keep=nk,
        screen_recost=int(L["screen_recost_evals"]),
        quick=int(L["quick_evals_total"]),
        d_disc=mind(L["R_disc"]),
        d_mcin=mind(L["R_disc"][L["mc_idx"]]),
        d_after=mind(Ra),
        d_kept=mind(Ra[:nk]),
    )


def task(item: Tuple[Any, ...]) -> str:
    import common as C
    from optimizer_sweep import voxel_context
    from candidate_probe import Probe

    vidx, vpos, spec, path = item
    W = C.get_worker()
    t0 = time.time()
    keys = np.load(C.CACHE_DIR / "images" / f"v{vidx}_clean.npz")["keys"]
    C.attach(keys)
    vctx = voxel_context(W.ctx, vidx)
    rec = W.rec
    Rt = np.asarray(vctx.R_true, dtype=np.float64)
    rec.min_sin_eta = 0.0
    sp = rec.setup.search_params
    sp.candidate_dedup = spec
    probe = Probe(rec, keep_raw=False)
    try:
        res = probe.run(vctx.vertices, vctx.voxel.phase, C.run_seed(vpos, 0))
    finally:
        sp.candidate_dedup = "off"
        rec.min_sin_eta = 0.3
    err = float(C.err_deg(np.asarray(res.orientation, dtype=np.float64), Rt))
    levels = [level_record(L, Rt) for L in probe.levels]
    out = dict(
        voxel=int(vidx), vpos=int(vpos), spec=spec, err_deg=err, evals=probe.total_evals,
        evals_levels=int(sum(l["screen_recost"] + l["quick"] for l in levels)), levels=levels,
        wall_s=time.time() - t0,
    )  # fmt: skip
    Path(path).write_text(json.dumps(out, default=float) + "\n")
    return f"v{vidx} {spec}: evals={probe.total_evals} err={err:.3f} {time.time() - t0:.0f}s"


def run(arm: str, spec: str, workers: int, limit: int) -> None:
    import common as C

    info = dict(np.load(C.OUT_DIR / "voxels.npz"))
    vox = [int(v) for v in info["voxel_indices"]][: limit or C.N_VOXELS]
    root = CACHE / arm
    root.mkdir(parents=True, exist_ok=True)
    items = []
    for p, v in enumerate(vox):
        if not (C.CACHE_DIR / "e0" / f"v{v}_clean.npz").exists():
            continue
        fn = root / f"v{v}.json"
        if not fn.exists():
            items.append((v, p, spec, str(fn)))
    (root / "spec.txt").write_text(spec + "\n")
    print(f"{arm} ({spec}): {len(items)} tasks, {workers} workers", flush=True)
    t0 = time.time()
    C.run_pool(task, items, workers, f"c_e0_{arm}")
    (root / "run_meta.json").write_text(
        json.dumps(
            dict(workers=workers, wall_s=time.time() - t0, tasks=len(items), spec=spec,
                 note="contended multi-worker wall, not a timing claim")
        )  # fmt: skip
    )


def load_arm(arm: str) -> Dict[int, Dict[str, Any]]:
    out = {}
    for fn in sorted((CACHE / arm).glob("v*.json")):
        d = json.loads(fn.read_text())
        out[d["voxel"]] = d
    return out


def reference_rows(
    B: Dict[int, Dict[str, Any]], wrong_deg: float, err_deg: Any
) -> List[Dict[str, Any]]:
    """The baseline arm against the stored E0 runs (REFERENCE_RUNS) run by run: identical errors
    and the paired split of the wrong flags. Empty when the stored runs are not on this machine."""
    out: List[Dict[str, Any]] = []
    vox = sorted(B)
    for run, what in REFERENCE_RUNS.items():
        files = [PHASE_D_CACHE / run / f"v{v}_clean.npz" for v in vox]
        if not all(f.exists() for f in files):
            continue
        z = [np.load(f) for f in files]
        er = err_deg(np.stack([x["R"] for x in z]), np.stack([x["R_true"] for x in z]))
        eb = np.array([B[v]["err_deg"] for v in vox])
        wr, wb = er > wrong_deg, eb > wrong_deg
        out.append(
            dict(
                run=run, what=what, n=len(vox), wrong_run=int(wr.sum()), wrong_base=int(wb.sum()),
                same_error_1e9=int(np.sum(np.abs(er - eb) < 1e-9)),
                max_abs_err_diff=float(np.max(np.abs(er - eb))),
                wrong_only_in_run=int((wr & ~wb).sum()), wrong_only_in_base=int((wb & ~wr).sum()),
                both_wrong=int((wr & wb).sum()),
            )
        )  # fmt: skip
    return out


def summary(arms: List[str], base: str = "base") -> None:
    import common as C
    import doc_tables as DT
    import stats as S

    B = load_arm(base)
    summ: Dict[str, Any] = {}
    rows: List[Dict[str, Any]] = []
    rows_t: List[Dict[str, Any]] = []
    rows_l0: List[Dict[str, Any]] = []
    for arm in arms:
        A = load_arm(arm)
        vox = sorted(set(A) & set(B))
        if not vox:
            continue
        ea = np.array([A[v]["err_deg"] for v in vox])
        eb = np.array([B[v]["err_deg"] for v in vox])
        wa, wb = ea > C.WRONG_DEG, eb > C.WRONG_DEG
        broken, fixed = S.paired_discordant(wa, wb)
        n = len(vox)
        lo, hi = S.wilson(int(wa.sum()), n)
        eva = np.array([A[v]["evals"] for v in vox], dtype=float)
        evb = np.array([B[v]["evals"] for v in vox], dtype=float)
        both = ~wa & ~wb
        spec = (CACHE / arm / "spec.txt").read_text().strip()
        s: Dict[str, Any] = dict(
            spec=spec, n=n, wrong=int(wa.sum()), wrong_base=int(wb.sum()),
            wilson_lo=lo, wilson_hi=hi, broken=broken, fixed=fixed,
            both_wrong=int((wa & wb).sum()), mean_evals=float(eva.mean()),
            mean_evals_base=float(evb.mean()), evals_ratio=float(eva.sum() / evb.sum()),
            evals_total=int(eva.sum()),
            median_err_right_both=float(np.median(ea[both])) if both.any() else float("nan"),
            n_right_both=int(both.sum()),
            err_equal_base=int(np.sum(ea == eb)),
        )  # fmt: skip
        # level 0: runs in which a filter dropped candidates before the quick MC (the later random
        # draws then differ), and the runs that lose the truth basin (3 deg) from the quick-MC
        # inputs against the baseline, with their final outcome (a signal that uses no random
        # numbers: the level-0 discrete stage is deterministic)
        s["L0_filter_drop_runs"] = int(
            sum(A[v]["levels"][0]["n_mc"] < A[v]["levels"][0]["n_disc"] for v in vox)
        )
        b0 = np.array([B[v]["levels"][0]["d_mcin"] < 3.0 for v in vox])
        a0 = np.array([A[v]["levels"][0]["d_mcin"] < 3.0 for v in vox])
        for tag, m in (("lost", b0 & ~a0), ("gained", ~b0 & a0)):
            s[f"L0_mcin3_{tag}"] = dict(
                runs=int(m.sum()),
                broken=int((m & ~wb & wa).sum()), fixed=int((m & wb & ~wa).sum()),
                both_wrong=int((m & wb & wa).sum()), both_right=int((m & ~wb & ~wa).sum()),
            )  # fmt: skip
        # per-level evaluations and truth-basin presence (3 deg)
        nl = max(len(A[v]["levels"]) for v in vox)
        for l in range(nl):
            Ls = [A[v]["levels"][l] for v in vox if l < len(A[v]["levels"])]
            s[f"L{l}_mean_candidates"] = float(np.mean([x["n_disc"] for x in Ls]))
            s[f"L{l}_mean_mc_inputs"] = float(np.mean([x["n_mc"] for x in Ls]))
            s[f"L{l}_mean_kept"] = float(np.mean([x["n_keep"] for x in Ls]))
            s[f"L{l}_mean_screen_recost"] = float(np.mean([x["screen_recost"] for x in Ls]))
            s[f"L{l}_mean_quick"] = float(np.mean([x["quick"] for x in Ls]))
            for b in BASIN:
                for stage in ("disc", "mcin", "after", "kept"):
                    s[f"L{l}_truth{b:g}_{stage}"] = int(sum(x[f"d_{stage}"] < b for x in Ls))
            s[f"L{l}_runs"] = len(Ls)
        # truth basin in the kept quarter at the end of level 2 (the last cut that feeds level 3)
        summ[arm] = s
        rows.append(
            dict(
                arm=arm, spec=spec, n=n, wrong=f"{int(wa.sum())} [{lo:.3f}, {hi:.3f}]",
                broken=broken, fixed=fixed, both_wrong=s["both_wrong"],
                mean_evals=s["mean_evals"], vs_base=f"{100 * s['evals_ratio']:.1f}%",
                n_right_both=s["n_right_both"], med_err_right_both=s["median_err_right_both"],
            )
        )  # fmt: skip
        lost, gained = s["L0_mcin3_lost"], s["L0_mcin3_gained"]
        rows_l0.append(
            dict(
                arm=arm, spec=spec, l0_filter_drop_runs=s["L0_filter_drop_runs"],
                mcin_L0=s["L0_truth3_mcin"], lost_L0=lost["runs"], lost_broken=lost["broken"],
                lost_fixed=lost["fixed"], lost_both_wrong=lost["both_wrong"],
                lost_both_right=lost["both_right"], gained_L0=gained["runs"],
                broken=broken, fixed=fixed,
            )
        )  # fmt: skip
        rows_t.append(
            dict(
                arm=arm,
                **{
                    f"kept_L{l}": f"{s[f'L{l}_truth3_kept']}/{s[f'L{l}_runs']}"
                    for l in range(nl)
                },
                **{f"mcin_L{l}": f"{s[f'L{l}_truth3_mcin']}" for l in range(nl)},
            )
        )  # fmt: skip
    rows_ref = reference_rows(B, C.WRONG_DEG, C.err_deg)
    if rows_ref:
        summ["_reference"] = {r["run"]: r for r in rows_ref}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary_e0.json").write_text(json.dumps(summ, indent=1, default=float) + "\n")
    tabs = {
        "c_e0_arms": DT.markdown_table(
            rows,
            ["arm", "spec", "n", "wrong", "broken", "fixed", "both_wrong", "mean_evals",
             "vs_base", "n_right_both", "med_err_right_both"],
            formats={"mean_evals": ",.0f", "med_err_right_both": ".4f"},
        ),  # fmt: skip
        "c_e0_truth": DT.markdown_table(rows_t, list(rows_t[0].keys())),
        "c_e0_level0": DT.markdown_table(rows_l0, list(rows_l0[0].keys())),
    }
    if rows_ref:
        tabs["c_e0_reference"] = DT.markdown_table(
            rows_ref, list(rows_ref[0].keys()), formats={"max_abs_err_diff": ".2e"}
        )
    DT.write_tables(OUT / "tables_e0.md", tabs)
    print(tabs["c_e0_arms"])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("cmd", choices=["run", "summary"])
    ap.add_argument("--arm", default="base")
    ap.add_argument("--spec", default="off")
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--arms", nargs="*", default=[])
    a = ap.parse_args()
    if a.cmd == "run":
        run(a.arm, a.spec, a.workers, a.limit)
    else:
        arms = a.arms or sorted(p.name for p in CACHE.iterdir() if p.is_dir())
        summary(arms)


if __name__ == "__main__":
    main()
