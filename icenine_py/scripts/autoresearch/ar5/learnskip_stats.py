"""AR5 step 18: the numbers of a finished arm ``J5L`` run of one sample, against arm ``J5`` (aligned
draws) and against ``S8``, all with the boundary pass, and the pass rules (1) to (5) of the plan.

  RESEED_SAMPLE=<held4|held3|held|mid|grad|big> LEARNSKIP_CACHE=<dir> uv run python \
      scripts/autoresearch/ar5/learnskip_stats.py [--out-dir DIR]

Reads ``<cache>/<arm>_<s>/`` and ``<cache>/passes_<s>/`` for the arms ``S8`` (stored for held, mid,
grad and big; for held3 and held4 run into this cache), ``J5`` (step 16's stored run, linked per
file; for held4 run into this cache) and ``J5L``. Uses step 16's numbers (``high500_stats``: stages,
error vectors, fixed and broken by voxel, the seed steps matched by (voxel, occurrence) with wrong
seeds). Writes ``stats_<s>.json`` (default: ``ar5/learnskip/``). ``held4`` is read only when the
marker ``<cache>/held4_runs_done.json`` exists (``learnskip_guard``). Truth for scoring only; no
p-values.

Rules (1) to (5): ``learnskip_rules.judge`` (pure): total with the pass at least 3% fewer than
``J5``; broken minus fixed (against ``J5``) at most 5; median error at most 0.001 degree above
``J5``'s; at most 1 more wrong seed (a taken seed end wrong where ``J5``'s was right, matched or
unmatched), grid searches within 1, reseeds at most 2% more; the house rule against ``S8`` with the
pass (at least 20% fewer evaluations, broken minus fixed at most 5, median at most 0.005 degree
higher); if ``J5`` fails (5) on the sample, ``J5L`` is judged by (1) to (4).
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict

import learnskip_common as LCM  # noqa: E402  (first)

import fit_common as FC  # noqa: E402
import gradient_common as G  # noqa: E402
import high500_stats as HS  # noqa: E402
import learnskip_guard as LG  # noqa: E402
import learnskip_rules as LR  # noqa: E402
import order_logs as OL  # noqa: E402

N = 24570
S = FC.SAMPLE
CACHE = LCM.LCACHE


def dirs(arm: str) -> Any:
    if arm == "S8":
        if S in ("held3", "held4"):
            return CACHE / f"S8_{S}", CACHE / f"passes_{S}"
        return HS.s8_dirs(S)
    return CACHE / f"{arm}_{S}", CACHE / f"passes_{S}"


def summarize(a: Dict[str, Any], grid: int) -> Dict[str, Any]:
    return dict(
        total=int(a["stages"]["total"]),
        median=float(a["after"]["median_err_deg"]),
        reseeds=int(a["run"]["stats"]["n_reseeds"]),
        grid=grid,
    )


def main() -> None:
    out_dir = Path(sys.argv[sys.argv.index("--out-dir") + 1]) if "--out-dir" in sys.argv else LCM.LS
    if S == "held4" and not LG.held4_open():
        LG.refuse("learnskip guard: held4's results are not read before its runs end")
    T = G.Truth()
    arms: Dict[str, Dict[str, Any]] = {}
    for n in ("S8", "J5", "J5L"):
        rd, pd = dirs(n)
        arms[n] = HS.arm(n, rd, pd, T)
    pj = {n: json.loads((dirs(n)[0] / "polish.json").read_text()) for n in ("J5", "J5L")}
    arms["S8"]["grid_searches"] = arms["S8"]["failed_seed_tests"]
    for n in ("J5", "J5L"):
        arms[n]["grid_searches"] = int(pj[n]["grid_fallbacks_left"])
    logs = {n: LCM.load_gz(dirs(n)[0] / "e1_log.json.gz") for n in ("J5", "J5L")}
    ss = HS.seed_steps(
        OL.index(logs["J5"]), OL.index(logs["J5L"]), T, arms["J5"]["eb"], arms["J5L"]["eb"]
    )
    vs_j5 = HS.pair(arms["J5"], arms["J5L"])
    vs_s8 = HS.pair(arms["S8"], arms["J5L"])
    j5_vs_s8 = HS.pair(arms["S8"], arms["J5"])
    j5l = dict(
        summarize(arms["J5L"], arms["J5L"]["grid_searches"]),
        wrong_seeds=int(ss["wrong_seeds"]["broad"]),
        vs_j5=vs_j5["after"],
        vs_s8=vs_s8["after"],
    )
    j5 = dict(
        summarize(arms["J5"], arms["J5"]["grid_searches"]),
        wrong_seeds=0,
        vs_s8=j5_vs_s8["after"],
    )
    s8 = dict(summarize(arms["S8"], arms["S8"]["grid_searches"]))
    rules = LR.judge(j5l, j5, s8)
    # the learned skip's own counts (the hook's file)
    ls = json.loads((dirs("J5L")[0] / "learnskip.json").read_text())
    res = dict(
        sample=S,
        arms={n: HS.strip(a) for n, a in arms.items()},
        j5l_vs_j5=vs_j5,
        j5l_vs_s8=vs_s8,
        j5_vs_s8=j5_vs_s8,
        seed_steps=dict(ss, unmatched=None),
        learnskip=ls,
        summary=dict(j5l=j5l, j5=j5, s8=s8),
        rules=rules,
        validity=dict(
            finished=all(a["n_processed"] == N for a in arms.values()),
            none_capped=all(not (dirs(n)[0] / "capped.json").exists() for n in arms),
            aligned_draws=bool(
                json.loads((dirs("J5L")[0] / "order.json").read_text())["aligned_draws"]
            ),
        ),
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / f"stats_{S}.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    print(
        json.dumps(
            dict(sample=S, summary=res["summary"], rules=rules, validity=res["validity"]), indent=1
        )
    )


if __name__ == "__main__":
    main()
