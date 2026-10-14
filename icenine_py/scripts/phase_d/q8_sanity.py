#!/usr/bin/env python3
"""Q8 step 2: three image variants side by side (clean, realistic, realistic_nodrop) at the truth.

For 30 full-sample voxels (the ray-overlap study's selection) and each variant, evaluates the hard
cost at the TRUE orientation under four overlap rules (old = lenient default, cpp = the C++ rule,
f0_R2 / f0_A = the ray rule with the 2-pixel and the geometric window) and three settings:

  screen  initial screen (pixel radius 3, Q-max 5)       as ray_overlap_study.py step 1
  local   no search radius (pixel radius 0, Q-max 5)     as ray_overlap_study.py step 1
  full    no search radius, Q-max 8 (the reconstruction's Q-max)

Everything is read at the truth, so each detector-0 overlap is a real peak and every one a rule
does not keep is a rejected real peak: rejected = det0 overlaps (the shared detector-0 test,
counted by the ray-rule arms, case det0_bright) minus peaks kept. The cost at the truth is
1 - quality. The 6-voxel Q-max 8 cost of scripts/phase_d/sanity_checks.py (q_true) is read from
results/sanity_full_<variant>.json.

  uv run python scripts/phase_d/q8_sanity.py run --workers 6      # builds the stacks if absent
  uv run python scripts/phase_d/q8_sanity.py summary

Writes benchmarks/q8_nodrop/sanity_truth.json (per voxel, per variant), sanity_summary.json and
sanity_tables.md (doc_tables blocks). No timing claims.
"""

import argparse
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
EX = HERE.parents[2] / "Examples" / "Example2.ManyGrains"
OUT = ICE / "benchmarks" / "q8_nodrop"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICE / "scripts" / "common"))

VARIANTS = ("clean", "realistic", "realistic_nodrop")
IMG = {
    "clean": "full/clean",
    "realistic": "full/realistic",
    "realistic_nodrop": "full/realistic_nodrop",
}
ARMS = ("old", "cpp", "f0_R2", "f0_A")
SETTINGS = (("screen", 3, 5.0), ("local", 0, 5.0), ("full", 0, 8.0))
COLS = ("peak_overlap", "peak_on_detector", "pixel_overlap", "pixel_on_detector", "quality")
DET0_BRIGHT = 5 + 7  # column of the case count det0_bright in a row (see ray_overlap_study.py)
STACKS = EX / "ScatteringData_PhaseD" / "stacks"


def build_stack(variant: str) -> float:
    """Write the uint8 stack of a variant once (as memmap_check.py); returns the seconds taken."""
    from icenine.experimental_data import write_binary_stack

    stack = STACKS / f"{variant}.npy"
    if stack.exists():
        return 0.0
    t0 = time.time()
    write_binary_stack(
        EX / "ScatteringData_PhaseD" / IMG[variant],
        "500Grains.sim",
        "d",
        5,
        180,
        2,
        2048,
        2048,
        stack,
    )
    return time.time() - t0


def truth_task(job: Tuple[str, int]) -> Dict[str, Any]:
    variant, voxel = job
    import ray_overlap_study as RS
    from icenine import reconstructor as RC
    from icenine.cost_functions import VoxelCostFunction

    RS.init_worker(variant)
    setup = RS.Ctx.setup
    v = setup.sample.get_mic().voxels[voxel]
    verts = RC._get_voxel_vertices(v)
    R = np.asarray(RS.Ctx.truth.voxels[voxel].orientation, dtype=np.float64).astype(np.float32)
    reg = float(setup.exp_setup.config_file.detector_dist_deviation)
    res: Dict[str, Any] = dict(voxel=voxel, variant=variant)
    for fn_name, radius, max_q in SETTINGS:
        for arm in ARMS:
            cpp, ray = RS.arm_spec(arm, reg)
            cf = VoxelCostFunction(
                simulator=setup.simulator,
                detector_list=setup.detector_list,
                range_map=setup.range_map,
                exp_data=setup.exp_data,
                sample=setup.sample,
                structure_list=setup.structure_list,
                mode="hard",
                eta_limit=setup.exp_setup.get_eta_limit(),
                pixel_radius=radius,
                max_q=max_q,
                min_sin_eta=0.0,
                cpp_overlap_rule=cpp,
                center_screen=RS.CENTER_SCREEN,
                ray_rule=ray,
            )
            info = cf.evaluate(R, verts, v.phase)
            res[f"{fn_name}|{arm}"] = [
                info.peak_overlap,
                info.peak_on_detector,
                info.pixel_overlap,
                info.pixel_on_detector,
                info.quality,
            ] + list(info.ray_cases if info.ray_cases is not None else [0] * RS.N_CASES)
    return res


def run(workers: int) -> None:
    import ray_overlap_study as RS

    OUT.mkdir(parents=True, exist_ok=True)
    t_stack = {v: build_stack(v) for v in VARIANTS}
    voxels = RS.select_voxels()
    jobs = [(var, v) for var in VARIANTS for v in voxels]
    t0 = time.time()
    with mp.get_context("spawn").Pool(workers) as pool:
        rows = pool.map(truth_task, jobs, chunksize=1)
    (OUT / "sanity_truth.json").write_text(
        json.dumps(dict(voxels=voxels, stack_write_s=t_stack, wall_s=time.time() - t0, rows=rows))
        + "\n"
    )
    print(f"{len(rows)} (variant, voxel) rows, {time.time() - t0:.0f}s (not a timing claim)")


def summary() -> None:
    import doc_tables
    import stats as S

    d = json.loads((OUT / "sanity_truth.json").read_text())
    rows = d["rows"]
    summ: Dict[str, Any] = {"n_voxels": len(d["voxels"])}
    rej_rows: List[Dict[str, Any]] = []
    cost_rows: List[Dict[str, Any]] = []
    per_var: Dict[str, Dict[int, Dict[str, Any]]] = {v: {} for v in VARIANTS}
    for r in rows:
        per_var[r["variant"]][r["voxel"]] = r
    voxels = d["voxels"]
    for var in VARIANTS:
        rs = [per_var[var][v] for v in voxels]
        for fn, _, _ in SETTINGS[:2]:
            ref = int(sum(np.array(r[f"{fn}|f0_R2"])[DET0_BRIGHT] for r in rs))
            # the shared detector-0 test: the geometric window must count the same overlaps
            assert ref == int(sum(np.array(r[f"{fn}|f0_A"])[DET0_BRIGHT] for r in rs))
            row: Dict[str, Any] = {"variant": var, "setting": fn, "det0_overlaps": ref}
            for arm in ARMS:
                kept = int(sum(r[f"{fn}|{arm}"][0] for r in rs))
                rej = ref - kept
                lo, hi = S.wilson(max(rej, 0), ref)
                summ[f"{var}|{fn}|{arm}"] = dict(kept=kept, det0=ref, rejected=rej, lo=lo, hi=hi)
                row[arm] = f"{rej} of {ref} ({100.0 * rej / ref:.1f}%)" if ref else "n/a"
                row[arm + "_ci"] = f"{100 * lo:.2f} to {100 * hi:.2f}%"
            rej_rows.append(row)
        crow: Dict[str, Any] = {"variant": var}
        for arm in ARMS:
            q = np.array([r[f"full|{arm}"][4] for r in rs], dtype=float)
            c = 1.0 - q
            summ[f"{var}|cost_full|{arm}"] = dict(
                median=float(np.median(c)), mean=float(c.mean()), min=float(c.min()),
                max=float(c.max()), n=len(c),
            )  # fmt: skip
            crow[arm] = float(np.median(c))
        cost_rows.append(crow)
    # paired per voxel: cost at the truth of each variant against clean (same voxels)
    for var in ("realistic", "realistic_nodrop"):
        for arm in ARMS:
            diff = np.array(
                [
                    (1 - per_var[var][v][f"full|{arm}"][4])
                    - (1 - per_var["clean"][v][f"full|{arm}"][4])
                    for v in voxels
                ]
            )
            summ[f"{var}|cost_full_minus_clean|{arm}"] = dict(
                median=float(np.median(diff)), mean=float(diff.mean()), n=len(diff),
                n_higher=int((diff > 1e-12).sum()), n_equal=int((np.abs(diff) <= 1e-12).sum()),
                n_lower=int((diff < -1e-12).sum()),
            )  # fmt: skip
    # the 6-voxel Q-max 8 cost of sanity_checks.py (lenient rule, min_sin_eta default 0.3)
    six: Dict[str, Any] = {}
    for var in VARIANTS:
        p = HERE / "results" / f"sanity_full_{var}.json"
        if p.exists():
            q = [r["q_true"] for r in json.loads(p.read_text())["cost"]["rows"]]
            six[var] = dict(
                voxels=[r["voxel"] for r in json.loads(p.read_text())["cost"]["rows"]],
                q_true=q,
                cost_median=float(np.median(1 - np.array(q))),
                cost_mean=float(np.mean(1 - np.array(q))),
            )
    summ["sanity_checks_six_voxels"] = six
    (OUT / "sanity_summary.json").write_text(json.dumps(summ, indent=1, default=float) + "\n")

    for r in cost_rows:
        r["six_voxel_cost_median"] = six.get(r["variant"], {}).get("cost_median")
    t_cost = doc_tables.markdown_table(
        cost_rows,
        ["variant", "six_voxel_cost_median", "old", "cpp", "f0_R2", "f0_A"],
        formats={k: ".4f" for k in ("six_voxel_cost_median", "old", "cpp", "f0_R2", "f0_A")},
    )
    cols = ["variant", "setting", "det0_overlaps"] + [
        c for a in ARMS if a != "old" for c in (a, a + "_ci")
    ]
    t_rej = doc_tables.markdown_table(
        [r for r in rej_rows if r["setting"] in ("screen", "local")], cols
    )
    doc_tables.write_tables(
        OUT / "sanity_tables.md", {"q8_sanity_cost": t_cost, "q8_sanity_rejected": t_rej}
    )
    print(t_cost)
    print(t_rej)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("cmd", choices=["run", "summary"])
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    if a.cmd == "run":
        run(a.workers)
    else:
        summary()
