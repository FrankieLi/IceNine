"""AR5 step 11, phase A step 5: the numbers of the finished arm S run of one sample.

  RESEED_SAMPLE=<held|mid|grad|big> uv run python scripts/autoresearch/ar5/patience_stats.py

Runs step 2's ``fit_part2_stats`` for arm ``G`` (stored step 9 run, linked) and arm ``S`` (each
against the stored ``S8``), then adds: the stage split of S (rescue evaluations apart from the grid
fallbacks); cuts (refinements, cut, evaluations in cut ones, the hit ratio of each cut end = the
running best at the stall); seed steps matched with G's (same voxel, same occurrence): fits both
ran compared exactly (a difference must be a cut MC fit), taken seed ends that changed (hit ratio,
cost, error, angle to G's), wrong seeds (S's taken end wrong where G's was right); failed seed
tests, rescues, grid searches; fixed and broken of S against G and against S8 (before and after the
pass, by voxel); the decision rules (``patience_judge``). Writes ``patience/stats_<sample>.json``.
Truth is used for scoring only; no p-values.
"""

import gzip
import json
import sys
from typing import Any, Dict, List, Set, Tuple

import patience_common as PC  # noqa: E402  (first)

import fit_common as FC  # noqa: E402
import fit_part2_stats as P2  # noqa: E402
import numpy as np  # noqa: E402
import order_logs as OL  # noqa: E402
import patience_judge as PJ  # noqa: E402

import gradient_common as G  # noqa: E402

import stats as ST  # noqa: E402

_HELD_AR3 = (
    PC.RCM.OC.WT / "agent-af8961874e131bda2/icenine_py/benchmarks/autoresearch/ar3/held/cache"
)
FC.STORED_S8.setdefault("held", _HELD_AR3 / "S8_held")
FC.STORED_PASS.setdefault("held", _HELD_AR3 / "passes_held")
SAMPLE = FC.SAMPLE
N = 24570


def read_gz(p: Any) -> Any:
    with gzip.open(p, "rt") as f:
        return json.load(f)


def same_R(a: Any, b: Any) -> bool:
    x = np.round(np.asarray(a, np.float64).ravel(), 12)
    y = np.round(np.asarray(b, np.float64).ravel(), 12)
    return bool(np.max(np.abs(x - y)) <= 1e-9)


def main() -> None:
    sd = PC.PCACHE / "stats"
    sd.mkdir(parents=True, exist_ok=True)
    FC.FIT = sd
    per: Dict[str, Any] = {}
    for arm in ("G", "S"):
        sys.argv = ["fit_part2_stats.py", "--name", arm]
        P2.main()
        d = json.loads((sd / f"part2_{arm}_{SAMPLE}.json").read_text())
        assert not d["unfinished"], arm
        per[arm] = d
    T = G.Truth()
    gd = PC.PCACHE / f"G_{SAMPLE}"
    sdir = PC.PCACHE / f"S_{SAMPLE}"
    rg = json.loads((gd / "rescue.json").read_text())
    rs = json.loads((sdir / "rescue.json").read_text())
    for r in (rg, rs):
        r.pop("rule", None)
    pat = json.loads((sdir / "patience.json").read_text())
    pat.pop("rule", None)
    order_s = json.loads((sdir / "order.json").read_text())
    st = dict(per["S"]["stages"])
    st["rescue"] = int(rs["rescue_evals"])
    st["grid_fallback"] = int(st["grid_fallback"]) - st["rescue"]
    assert st["grid_fallback"] == int(rs["grid_evals"]), "grid evals do not match the log"
    per["S"]["stages_split"] = st
    stg = dict(per["G"]["stages"])
    stg["rescue"] = int(rg["rescue_evals"])
    stg["grid_fallback"] = int(stg["grid_fallback"]) - stg["rescue"]
    per["G"]["stages_split"] = stg

    lg = read_gz(gd / "e1_log.json.gz")
    ls = read_gz(sdir / "e1_log.json.gz")
    ts = read_gz(sdir / "mc_traces.json.gz")
    ig, isx = OL.index(lg), OL.index(ls)
    row_of = {id(e): i for i, e in enumerate(ls)}
    cut_keys: Set[Tuple[int, int]] = {(int(t["step"]) - 1, int(t["j"])) for t in ts if t.get("cut")}
    # cut refinements: the hit ratio of the cut end (the running best at the stall)
    cut_hits: List[float] = []
    for ri, j in cut_keys:
        for f in ls[ri]["fits"]:
            if f[2] == "M" and f[1] == j:
                cut_hits.append(float(f[3]))
    keys = sorted(set(ig) & set(isx))
    fits_cmp = fits_diff = fits_diff_cut = fits_diff_other = 0
    changed: List[Dict[str, Any]] = []
    wrong_g_n = wrong_s_n = extra_wrong = 0
    Rg_all: List[np.ndarray] = []
    Rs_all: List[np.ndarray] = []
    vox: List[int] = []
    meta: List[Tuple[Dict[str, Any], Dict[str, Any]]] = []
    for k in keys:
        eg, es = ig[k], isx[k]
        if not eg["fits"] or not es["fits"]:
            continue
        ri = row_of[id(es)]
        fg = {(f[1], f[2]): f for f in eg["fits"]}
        fs = {(f[1], f[2]): f for f in es["fits"]}
        for key in set(fg) & set(fs):
            fits_cmp += 1
            x, y = fg[key], fs[key]
            if (x[3], x[4], x[5], x[6]) != (y[3], y[4], y[5], y[6]):
                fits_diff += 1
                if key[1] == "M" and (ri, key[0]) in cut_keys:
                    fits_diff_cut += 1
                else:
                    fits_diff_other += 1
        tg, tsx = OL.taken_end(eg), OL.taken_end(es)
        Rg_all.append(np.asarray(tg["R"], np.float64).reshape(3, 3))
        Rs_all.append(np.asarray(tsx["R"], np.float64).reshape(3, 3))
        vox.append(int(k[0]))
        meta.append((tg, tsx))
    Rg_a, Rs_a = np.stack(Rg_all), np.stack(Rs_all)
    err_g = np.array([float(T.err(Rg_a[i][None], [vox[i]])[0]) for i in range(len(vox))])
    err_s = np.array([float(T.err(Rs_a[i][None], [vox[i]])[0]) for i in range(len(vox))])
    mis = np.asarray(T.mis(Rg_a, Rs_a), dtype=float)
    for i, (tg, tsx) in enumerate(meta):
        wg, ws = err_g[i] > G.WRONG, err_s[i] > G.WRONG
        wrong_g_n += int(wg)
        wrong_s_n += int(ws)
        extra_wrong += int(ws and not wg)
        if not same_R(tg["R"], tsx["R"]):
            changed.append(
                dict(
                    voxel=vox[i],
                    stopped_g=bool(tg["hit"] >= 0.95),
                    hit_g=float(tg["hit"]),
                    hit_s=float(tsx["hit"]),
                    cost_g=float(tg["cost"]),
                    cost_s=float(tsx["cost"]),
                    err_g=float(err_g[i]),
                    err_s=float(err_s[i]),
                    mis_to_g_deg=float(mis[i]),
                )
            )
    # voxels, S against G and against S8
    e_g = dict(np.load(sd / f"part2_G_{SAMPLE}_err.npz"))
    e_s = dict(np.load(sd / f"part2_S_{SAMPLE}_err.npz"))
    pair: Dict[str, Any] = {}
    for tag in ("before", "after"):
        wg, ws = e_g[tag] > G.WRONG, e_s[tag] > G.WRONG
        fixed, broken = ST.paired_discordant(wg, ws)  # fixed: wrong in G, right in S
        pair[tag] = dict(fixed=int(fixed), broken=int(broken), both_wrong=int((wg & ws).sum()))
    hs, hg = per["S"], per["G"]
    house_s = PJ.house(
        hs["stages"]["total"],
        hs["stages_s8"]["total"],
        hs["paired_after"]["broken"],
        hs["paired_after"]["fixed"],
        hs["after"]["median_err_deg"],
        hs["s8_after"]["median_err_deg"],
    )
    house_g = PJ.house(
        hg["stages"]["total"],
        hg["stages_s8"]["total"],
        hg["paired_after"]["broken"],
        hg["paired_after"]["fixed"],
        hg["after"]["median_err_deg"],
        hg["s8_after"]["median_err_deg"],
    )
    rules = PJ.rules(
        SAMPLE,
        hs["stages"]["total"],
        hg["stages"]["total"],
        hs["stages"]["seed_step"],
        hg["stages"]["seed_step"],
        pair["after"]["broken"] - pair["after"]["fixed"],
        hs["after"]["median_err_deg"],
        hg["after"]["median_err_deg"],
        extra_wrong,
        int(rs["grid_fallbacks_left"]),
        int(rg["grid_fallbacks_left"]),
        house_s,
        house_g,
    )
    validity = dict(
        s_finished=bool(not (sdir / "capped.json").exists() and (sdir / "run.json").exists()),
        n_processed=int(json.loads((sdir / "run.json").read_text())["n_processed"]),
        aligned_draws=bool(order_s["aligned_draws"]),
        seed_steps_g=len(lg),
        seed_steps_s=len(ls),
        matched_steps=len(keys),
        fits_compared=fits_cmp,
        fits_differ=fits_diff,
        fits_differ_cut_mc=fits_diff_cut,
        fits_differ_other=fits_diff_other,
    )
    cut_n = len(cut_keys)
    out = dict(
        sample=SAMPLE,
        arms=per,
        rescue=dict(G=rg, S=rs),
        cuts=dict(
            **pat,
            refinements_in_traces=len(ts),
            cut_refinements=cut_n,
            cut_end_hit_under_095=int(sum(1 for h in cut_hits if h < 0.95)),
            cut_end_hit_at_least_095=int(sum(1 for h in cut_hits if h >= 0.95)),
        ),
        taken_changed=changed,
        taken_changed_n=len(changed),
        taken_changed_median_mis_deg=(
            float(np.median([c["mis_to_g_deg"] for c in changed])) if changed else None
        ),
        wrong_seeds=dict(g=wrong_g_n, s=wrong_s_n, s_wrong_where_g_right=extra_wrong),
        s_vs_g=dict(voxels=pair),
        validity=validity,
        house_s=house_s,
        house_g=house_g,
        rules=rules,
    )
    (PC.PAT / f"stats_{SAMPLE}.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(json.dumps(dict(sample=SAMPLE, rules=rules, pair=pair, validity=validity), indent=1))


if __name__ == "__main__":
    main()
