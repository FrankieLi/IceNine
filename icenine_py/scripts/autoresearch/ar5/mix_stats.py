"""AR5 step 3, Part 2: numbers of one finished S11 run on one sample.

  scripts/autoresearch/ar5/mix_bg.sh <sample> stats_<name>_<sample> mix_stats.py --name S11

Step 2's ``fit_part2_stats.py`` (evaluations by stage, grid fallbacks, wrong voxels before and after
the pass with Wilson 95% intervals, fixed and broken against stored S8 paired by voxel, median
error, the pass rule ``fit_rules.s10_pass`` with step 2's numbers) with its output folder set to
``benchmarks/autoresearch/ar5/mix/`` (it writes ``part2_<name>_<sample>.json`` and
``part2_<name>_<sample>_err.npz`` there). Then adds, from the same run: the comparison with run
S10 of step 2 (its stored numbers in ``benchmarks/autoresearch/ar5/fit/``): evaluations, fixed and
broken paired by voxel against S10 (before and after the pass), and the seed-step evaluations of
the run against Part 1's count for the same sample (``part1.json``). A run that reached the cap is
written as unfinished with no scores.
"""

import json

import mix_common as MC  # noqa: E402  (first: sets the cache)

import fit_common as FC  # noqa: E402
import fit_part2_stats as P2  # noqa: E402
import numpy as np  # noqa: E402

import gradient_common as G  # noqa: E402

import stats as ST  # noqa: E402


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    name = ap.parse_known_args()[0].name
    FC.FIT = MC.MIX  # step 2's stats script writes its outputs here, not into step 2's folder
    P2.main()
    s = FC.SAMPLE
    p = MC.MIX / f"part2_{name}_{s}.json"
    out = json.loads(p.read_text())
    if out.get("unfinished"):
        return
    s10_dir = MC.ICE / "benchmarks" / "autoresearch" / "ar5" / "fit"
    d10 = json.loads((s10_dir / f"part2_S10_{s}.json").read_text())
    e10 = np.load(s10_dir / f"part2_S10_{s}_err.npz")
    e11 = np.load(MC.MIX / f"part2_{name}_{s}_err.npz")
    vs: dict = dict(
        evals_total_s10=d10["stages"]["total"],
        evals_before_pass_s10=d10["stages"]["before_pass"],
        stages_s10=d10["stages"],
        grid_fallbacks_s10=d10["grid_fallbacks"],
        after_s10=d10["after"],
        before_s10=d10["before"],
    )
    for tag in ("before", "after"):
        w10 = np.asarray(e10[tag], float) > G.WRONG
        w11 = np.asarray(e11[tag], float) > G.WRONG
        fixed, broken = ST.paired_discordant(w10, w11)
        vs[f"paired_{tag}"] = dict(
            fixed=int(fixed), broken=int(broken), both_wrong=int((w10 & w11).sum())
        )
    vs["fewer_than_s10"] = 1.0 - out["stages"]["total"] / d10["stages"]["total"]
    p1 = json.loads((MC.MIX / "part1.json").read_text())["by_sample"][s]["arms"]["X200"]
    seed = out["stages"]["seed_step"]
    vs["seed_step_evals_part1"] = p1["evals"]
    vs["seed_step_evals_run"] = seed
    vs["seed_step_run_over_part1"] = seed / p1["evals"] - 1.0
    out["vs_s10"] = vs
    mlog = MC.CACHE_DIR / f"{name}_{s}" / "mix_log.json"
    if mlog.exists():
        rows = json.loads(mlog.read_text())
        out["mc_share"] = dict(
            candidates=len(rows),
            mc_ran=sum(1 for r in rows if r["mc_ran"]),
            share=sum(1 for r in rows if r["mc_ran"]) / max(1, len(rows)),
        )
    p.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: out[k] for k in ("rule", "vs_s10", "mc_share") if k in out}, indent=1))


if __name__ == "__main__":
    main()
