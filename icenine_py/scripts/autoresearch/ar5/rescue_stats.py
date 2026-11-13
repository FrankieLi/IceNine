"""AR5 step 9, phase A step 6: the numbers of the finished arm G and arm A runs of one sample.

  RESEED_SAMPLE=<held|mid|grad|big> uv run python scripts/autoresearch/ar5/rescue_stats.py

Runs step 2's ``fit_part2_stats`` for arm ``A`` (the stored step 8 run, linked) and arm ``G``
(each against the stored ``S8``), then adds: the stage split of G with the rescue evaluations
apart from the grid fallbacks; failed seed steps, rescues tried, rescued, grid fallbacks left;
error and hit ratio of every rescued end; wrong rescued seeds; seed steps whose taken end
differs; fixed and broken of G against A (before and after the pass, by voxel); whether the map
and evaluations of G equal A's (``big``). Writes ``rescue/stats_<sample>.json``. Truth is used
for scoring only; no p-values.
"""

import gzip
import json
import sys
from typing import Any, Dict, List

import rescue_common as RCM  # noqa: E402  (first)

import fit_common as FC  # noqa: E402
import fit_part2_stats as P2  # noqa: E402
import numpy as np  # noqa: E402
import order_logs as OL  # noqa: E402
import rescue_judge as RJ  # noqa: E402

import gradient_common as G  # noqa: E402

import stats as ST  # noqa: E402

_HELD_AR3 = RCM.OC.WT / "agent-af8961874e131bda2/icenine_py/benchmarks/autoresearch/ar3/held/cache"
FC.STORED_S8.setdefault("held", _HELD_AR3 / "S8_held")
FC.STORED_PASS.setdefault("held", _HELD_AR3 / "passes_held")
SAMPLE = FC.SAMPLE
N = 24570


def read_gz(p: Any) -> Any:
    with gzip.open(p, "rt") as f:
        return json.load(f)


def main() -> None:
    sd = RCM.RCACHE / "stats"
    sd.mkdir(parents=True, exist_ok=True)
    FC.FIT = sd  # fit_part2_stats writes here (untracked)
    per: Dict[str, Any] = {}
    for arm in ("A", "G"):
        sys.argv = ["fit_part2_stats.py", "--name", arm]
        P2.main()
        d = json.loads((sd / f"part2_{arm}_{SAMPLE}.json").read_text())
        assert not d["unfinished"], arm
        per[arm] = d
    T = G.Truth()
    ga = RCM.RCACHE / f"G_{SAMPLE}"
    resc = json.loads((ga / "rescue.json").read_text())
    resc.pop("rule", None)
    log = read_gz(ga / "rescue_log.json.gz")
    # stages of G with the rescue evaluations apart (P2 counts them inside grid_fallback)
    st = dict(per["G"]["stages"])
    st["rescue"] = int(resc["rescue_evals"])
    st["grid_fallback"] = int(st["grid_fallback"]) - st["rescue"]
    assert st["grid_fallback"] == int(resc["grid_evals"]), "grid evals do not match the log"
    per["G"]["stages_split"] = st
    # rescued ends
    ends: List[Dict[str, Any]] = []
    for r in log:
        if r["rescued"]:
            t = r["tries"][-1]
            R = np.asarray(t["R"], np.float64).reshape(3, 3)
            ends.append(
                dict(
                    voxel=int(r["voxel"]),
                    n_tries=len(r["tries"]),
                    hit=float(t["hit"]),
                    err=float(T.err(R[None], [int(r["voxel"])])[0]),
                )
            )
    wrong_rescued = sum(1 for e in ends if e["err"] > G.WRONG)
    e_a = dict(np.load(sd / f"part2_A_{SAMPLE}_err.npz"))
    e_g = dict(np.load(sd / f"part2_G_{SAMPLE}_err.npz"))
    pair: Dict[str, Any] = {}
    for tag in ("before", "after"):
        wa, wg = e_a[tag] > G.WRONG, e_g[tag] > G.WRONG
        fixed, broken = ST.paired_discordant(wa, wg)  # fixed: wrong in A, right in G
        pair[tag] = dict(fixed=int(fixed), broken=int(broken), both_wrong=int((wa & wg).sum()))
    la = read_gz(RCM.RCACHE / f"A_{SAMPLE}" / "e1_log.json.gz")
    lg = read_gz(ga / "e1_log.json.gz")

    def wrong_fn(v: int, R: List[float]) -> bool:
        Rm = np.asarray(R, np.float64).reshape(3, 3)
        return float(T.err(Rm[None], [v])[0]) > G.WRONG

    seed_cmp = OL.compare_logs(la, lg, wrong_fn)
    Rb_a, _ = G.read_map(RCM.RCACHE / f"A_{SAMPLE}" / "recon.mic.gz", N)
    Rb_g, _ = G.read_map(ga / "recon.mic.gz", N)
    Ra_a = np.load(RCM.RCACHE / f"passes_{SAMPLE}" / "after_A.npy")
    Ra_g = np.load(RCM.RCACHE / f"passes_{SAMPLE}" / "after_G.npy")
    same_map = bool(np.array_equal(Rb_a, Rb_g) and np.array_equal(Ra_a, Ra_g))
    same_evals = per["A"]["stages"]["total"] == per["G"]["stages"]["total"]
    ha = per["A"]
    hg = per["G"]
    house_g = RJ.house(
        hg["stages"]["total"],
        hg["stages_s8"]["total"],
        hg["paired_after"]["broken"],
        hg["paired_after"]["fixed"],
        hg["after"]["median_err_deg"],
        hg["s8_after"]["median_err_deg"],
    )
    house_a = RJ.house(
        ha["stages"]["total"],
        ha["stages_s8"]["total"],
        ha["paired_after"]["broken"],
        ha["paired_after"]["fixed"],
        ha["after"]["median_err_deg"],
        ha["s8_after"]["median_err_deg"],
    )
    rules = RJ.rules(
        SAMPLE,
        hg["stages"]["total"],
        ha["stages"]["total"],
        int(resc["grid_fallbacks_left"]),
        int(ha["grid_fallbacks"]),
        pair["after"]["broken"] - pair["after"]["fixed"],
        hg["after"]["median_err_deg"],
        ha["after"]["median_err_deg"],
        wrong_rescued,
        bool(same_map and same_evals),
        house_g,
        house_a,
    )
    out = dict(
        sample=SAMPLE,
        arms=per,
        rescue=resc,
        rescued_ends=ends,
        rescued_hit_min=min((e["hit"] for e in ends), default=None),
        rescued_err_median=float(np.median([e["err"] for e in ends])) if ends else None,
        rescued_err_max=max((e["err"] for e in ends), default=None),
        wrong_rescued=wrong_rescued,
        g_vs_a=dict(voxels=pair, seed_steps=seed_cmp, same_map=same_map, same_evals=same_evals),
        house_g=house_g,
        house_a=house_a,
        rules=rules,
    )
    (RCM.RESCUE / f"stats_{SAMPLE}.json").write_text(
        json.dumps(out, indent=1, default=float) + "\n"
    )
    print(json.dumps(dict(sample=SAMPLE, rules=rules, pair=pair), indent=1))


if __name__ == "__main__":
    main()
