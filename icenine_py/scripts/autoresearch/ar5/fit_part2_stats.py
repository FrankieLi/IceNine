"""AR5 step 2, Part 2: numbers of one finished run against stored S8, on one sample.

  scripts/autoresearch/ar5/fit_bg.sh <sample> stats_<name>_<sample> fit_part2_stats.py --name S10

Reads the run ``<cache>/<name>_<sample>/`` and its pass ``<cache>/passes_<sample>/`` (made by
``fit_run.py`` and ``fit_pass.py``) and the stored S8 run and pass of the same sample; writes
``benchmarks/autoresearch/ar5/fit/part2_<name>_<sample>.json``: evaluations by stage and in total
with the pass, grid fallbacks, the wrapper counts, wrong voxels (more than 1 degree from the truth)
before and after the pass with Wilson 95% intervals, median error, fixed and broken against S8
(paired by voxel, before and after the pass), and the pass rule of the plan
(``fit_rules.s10_pass``; it is evaluated for runs named ``S10``). Truth is used for scoring only.
Also writes ``part2_<name>_<sample>_err.npz`` (float32 errors of the maps before and after the
pass).
A run that reached the cap (``capped.json``) is written as unfinished with no scores.
"""

import argparse
import json
from pathlib import Path
from typing import Any, Dict

import fit_common as FC
import fit_rules as FR
import numpy as np
import reseed_common as RC

import gradient_common as G  # noqa: E402

import stats as ST  # noqa: E402

N = 24570


def stage_evals(run: Dict[str, Any], pass_evals: int) -> Dict[str, int]:
    st = run["stats"]
    d = dict(
        seed_step=int(run["pick"]["evals"]),
        grid_fallback=int(st.get("index_fallback_grid_evals", 0)),
        neighbour=int(st["n_evals_neighbor"]),
        revisit=int(st["n_evals_revisit"]),
        pass_=int(pass_evals),
    )
    assert d["seed_step"] + d["grid_fallback"] + d["neighbour"] + d["revisit"] == int(
        run["evals_total"]
    ), "stages do not add up to evals_total"
    d["before_pass"] = int(run["evals_total"])
    d["total"] = d["before_pass"] + d["pass_"]
    return d


def score(R: np.ndarray, T: Any) -> np.ndarray:
    return np.asarray(T.err(R), dtype=float)


def side(err: np.ndarray) -> Dict[str, Any]:
    k = int((err > G.WRONG).sum())
    lo, hi = ST.wilson(k, len(err))
    return dict(wrong=k, wilson_lo=lo, wilson_hi=hi, median_err_deg=float(np.median(err)))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    a = ap.parse_args()
    s = FC.SAMPLE
    rd = FC.CACHE / f"{a.name}_{s}"
    out: Dict[str, Any] = dict(run=a.name, sample=s)
    if (rd / "capped.json").exists() or not (rd / "run.json").exists():
        out["unfinished"] = True
        out["capped"] = (
            json.loads((rd / "capped.json").read_text()) if (rd / "capped.json").exists() else None
        )
        (FC.FIT / f"part2_{a.name}_{s}.json").write_text(json.dumps(out, indent=1) + "\n")
        print("unfinished")
        return
    pd = FC.CACHE / f"passes_{s}"
    T = G.Truth()
    run = json.loads((rd / "run.json").read_text())
    fit = json.loads((rd / "fit.json").read_text())
    pj = json.loads((pd / f"D4_{a.name}" / "pass.json").read_text())
    old_run = json.loads((FC.STORED_S8[s] / "run.json").read_text())
    old_pj = json.loads((FC.STORED_PASS[s] / "D4_S8" / "pass.json").read_text())
    Rb, _ = G.read_map(rd / "recon.mic.gz", N)
    Ra = np.load(pd / f"after_{a.name}.npy")
    Rb0, _ = G.read_map(FC.STORED_S8[s] / "recon.mic.gz", N)
    Ra0 = np.load(FC.STORED_PASS[s] / "after_S8.npy")
    eb, ea, eb0, ea0 = score(Rb, T), score(Ra, T), score(Rb0, T), score(Ra0, T)
    np.savez_compressed(
        FC.FIT / f"part2_{a.name}_{s}_err.npz",
        before=eb.astype(np.float32),
        after=ea.astype(np.float32),
        s8_before=eb0.astype(np.float32),
        s8_after=ea0.astype(np.float32),
    )
    st = stage_evals(run, int(pj["evals"]))
    st0 = stage_evals(old_run, int(old_pj["evals"]))
    w = run["wrapper"]
    out.update(
        unfinished=False,
        seed_arm=fit["seed_arm"],
        no_mc_fallback=fit["no_mc_fallback"],
        seed_steps=fit["seed_steps"],
        stopped_early=fit["stopped_early"],
        stages=st,
        stages_s8=st0,
        grid_fallbacks=int(run["stats"].get("index_fallback_seed_test", 0)),
        grid_fallbacks_s8=int(old_run["stats"].get("index_fallback_seed_test", 0)),
        wrapper={k: w[k] for k in w if k != "wrapper_calls_by_caller"},
        before=side(eb),
        after=side(ea),
        s8_before=side(eb0),
        s8_after=side(ea0),
        wall_s=fit["wall_s"],
        wall_label=fit["wall_label"],
        wall_pass_s=pj["wall_tries_s"],
    )
    for tag, e, e0 in (("before", eb, eb0), ("after", ea, ea0)):
        wrong, wrong0 = e > G.WRONG, e0 > G.WRONG
        # paired_discordant(a, b) returns (wrong under a only, wrong under b only)
        fixed, broken = ST.paired_discordant(wrong0, wrong)  # fixed: wrong in S8, right in the run
        out[f"paired_{tag}"] = dict(
            fixed=int(fixed), broken=int(broken), both_wrong=int((wrong & wrong0).sum())
        )
    out["rule"] = FR.s10_pass(
        evals=st["total"],
        evals_s8=st0["total"],
        broken=out["paired_after"]["broken"],
        fixed=out["paired_after"]["fixed"],
        med=out["after"]["median_err_deg"],
        med_s8=out["s8_after"]["median_err_deg"],
    )
    (FC.FIT / f"part2_{a.name}_{s}.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1)[:3000])


if __name__ == "__main__":
    main()
