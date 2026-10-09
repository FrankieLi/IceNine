"""Collect the before/after/C++ numbers of the VarianceMinimizing restart fix into
benchmarks/phase_d_seed_diag/variance_parity.json and the doc tables.

Inputs: runs/base_full_v*.json (before), runs/varfix_full_v*.json (after), C++ logs with the
patched build's "VT" lines (--cpp-dir, files cpp_<voxel>.log; 15901 from cpp_trace2.log),
variance_parity/*.json from seed_diag_variance_parity.py, run_estimates.json for the seed counts.
Usage (from icenine_py/): uv run python scripts/phase_d/seed_diag_variance_parity_summary.py\
  --cpp-dir DIR
"""

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
OUT = ICE / "benchmarks" / "phase_d_seed_diag"
sys.path.insert(0, str(ICE / "scripts" / "common"))
import doc_tables as DT  # noqa: E402

VOXELS = [15901, 17034, 18558, 19889, 22867, 5242]


def cpp_steps(log: Path):
    steps, evals, us = [], [], []
    for line in log.read_text().splitlines():
        m = re.search(r"\|Step \|\s+(\d+)", line)
        if m:
            steps.append(int(m.group(1)))
        m = re.search(r"adap_evals=(\d+) adap_us_per_eval=([\d.]+)", line)
        if m:
            evals.append(int(m.group(1)))
            us.append(float(m.group(2)))
    return steps, evals, us


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpp-dir", required=True)
    a = ap.parse_args()
    rows, out = [], {"voxels": {}}
    for v in VOXELS:
        b = json.load(open(OUT / "runs" / f"base_full_v{v}.json"))
        f = json.load(open(OUT / "runs" / f"varfix_full_v{v}.json"))
        log = Path(a.cpp_dir) / ("cpp_trace2.log" if v == 15901 else f"cpp_{v}.log")
        steps, ev, us = cpp_steps(log)
        d = dict(
            before=dict(
                evals=b["evals_total"],
                variance_evals=b["variance_evals"],
                err_deg=b["err_deg"],
                q_true=b["q_true"],
                wall_s=b["wall_s"],
                quiet=b["quiet_preflight"],
            ),
            after=dict(
                evals=f["evals_total"],
                variance_evals=f["variance_evals"],
                err_deg=f["err_deg"],
                q_true=f["q_true"],
                cost_final=f["cost_final"],
                wall_s=f["wall_s"],
                us_per_eval=f["mean_us_per_eval"],
                quiet=f["quiet_preflight"],
                label=f["label"],
            ),
            cpp=dict(
                variance_steps=steps,
                evals=ev,
                us_per_eval=float(np.mean(us)),
                n_runs=len(steps),
            ),
        )
        out["voxels"][str(v)] = d
        rows.append(
            dict(
                voxel=v,
                var_before=b["variance_evals"],
                var_after=f["variance_evals"],
                var_cpp=int(np.median(steps)),
                evals_before=b["evals_total"],
                evals_after=f["evals_total"],
                evals_cpp=int(np.median(ev)),
                err_before=b["err_deg"],
                err_after=f["err_deg"],
                wall_after=f["wall_s"],
            )
        )
    m = lambda k: float(np.mean([r[k] for r in rows]))  # noqa: E731
    rows.append(dict(voxel="mean", **{k: m(k) for k in rows[0] if k != "voxel"}))
    fm = {
        k: ",.0f"
        for k in ("var_before", "var_after", "var_cpp", "evals_before", "evals_after", "evals_cpp")
    }
    fm.update(err_before=".3f", err_after=".3f", wall_after=".0f")
    tbl = DT.markdown_table(rows, columns=list(rows[0]), formats=fm)
    # seed-pass estimate: per-seed time = mean evals x quiet us per eval (timing_single_worker.json)
    us = json.load(open(OUT / "timing_single_worker.json"))["us_per_eval"]
    us_mean = float(np.mean(us))
    per_seed_eval = m("evals_after") * us_mean * 1e-6
    per_seed_wall_mean = m("wall_after")
    est = json.load(open(OUT / "run_estimates.json"))["rows"][0]
    seeds = [est["seeds_low"], est["seeds_mid"], est["seeds_high"]]
    out["estimate"] = dict(
        us_per_eval_quiet=us_mean,
        per_seed_s_evals_x_quiet_us=per_seed_eval,
        per_seed_s_measured_mean_wall=per_seed_wall_mean,
        measured_label=sorted({r["after"]["label"] for r in out["voxels"].values()}),
        seeds=seeds,
        hours_from_evals=[s * per_seed_eval / 3600 for s in seeds],
        hours_from_wall=[s * per_seed_wall_mean / 3600 for s in seeds],
    )
    vp = OUT / "variance_parity"
    out["cost_agreement_15901"] = json.load(open(vp / "costs.json"))
    for k in ("before", "after"):
        t = json.load(open(vp / f"trace_{k}.json"))
        out[f"trace_{k}_python_15901"] = t["python"]
        out["trace_cpp_15901"] = t["cpp"]
    (OUT / "variance_parity.json").write_text(json.dumps(out, indent=1) + "\n")
    est_rows = [
        dict(
            basis=name,
            per_seed_s=out["estimate"][key],
            hours_low=out["estimate"][hk][0],
            hours_mid=out["estimate"][hk][1],
            hours_high=out["estimate"][hk][2],
        )
        for name, key, hk in (
            ("evals x quiet us/eval", "per_seed_s_evals_x_quiet_us", "hours_from_evals"),
            ("measured wall (mean of 6)", "per_seed_s_measured_mean_wall", "hours_from_wall"),
        )
    ]
    t2 = DT.markdown_table(
        est_rows,
        columns=list(est_rows[0]),
        formats=dict(per_seed_s=".0f", hours_low=".1f", hours_mid=".1f", hours_high=".1f"),
    )
    DT.write_tables(
        OUT / "variance_parity_tables.md",
        {"variance_parity_seeds": tbl, "variance_parity_estimate": t2},
    )
    print(tbl)
    print(json.dumps(out["estimate"], indent=1))


if __name__ == "__main__":
    main()
