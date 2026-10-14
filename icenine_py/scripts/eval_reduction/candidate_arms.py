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


def summary(arms: List[str], base: str = "base") -> None:
    import common as C
    import doc_tables as DT
    import stats as S

    B = load_arm(base)
    summ: Dict[str, Any] = {}
    rows, rows_t = [], []
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
                med_err=s["median_err_right_both"],
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
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary_e0.json").write_text(json.dumps(summ, indent=1, default=float) + "\n")
    tabs = {
        "c_e0_arms": DT.markdown_table(
            rows,
            ["arm", "spec", "n", "wrong", "broken", "fixed", "both_wrong", "mean_evals",
             "vs_base", "med_err"],
            formats={"mean_evals": ",.0f", "med_err": ".4f"},
        ),  # fmt: skip
        "c_e0_truth": DT.markdown_table(rows_t, list(rows_t[0].keys())),
    }
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
