"""AR5 step 12, phase A step 6: the numbers of the finished arm H run of one sample.

  RESEED_SAMPLE=<held2|held|mid|grad|big> uv run python scripts/autoresearch/ar5/stallhit_stats.py

Runs step 2's ``fit_part2_stats`` for arm ``G`` (the stored step 9 run, linked; on ``held2`` the
run of this step) and arm ``H`` (each against ``S8``: the stored one, on ``held2`` the run of this
step), then adds: the stage split (rescue evaluations apart from the grid fallbacks); the stall
checks (from the traces: checks, stops, checks that went on, the hit ratio of each); cuts; seed
steps matched with G's (same voxel, same occurrence): fits both ran compared exactly (a difference
must be a cut MC fit, or an MC fit whose only difference is one counted evaluation per stall check
that went on); taken seed ends that changed (hit ratio, cost, error, angle to G's); the seed steps
only one arm has (``unmatched``: the taken end, the other arm's map there); wrong seeds (H's taken
end wrong where G's was right), matched and unmatched; failed seed tests, rescues, grid searches;
fixed and broken of H against G and against S8 (before and after the pass, by voxel); the decision
rules (``stallhit_judge``), rule (4) read broadly (matched plus unmatched). Writes
``stallhit/stats_<sample>.json`` and ``stallhit/unmatched_<sample>.json``. Truth is used for
scoring only; no p-values.
"""

import gzip
import json
import sys
from typing import Any, Dict, List, Set, Tuple

import stallhit_common as SC  # noqa: E402  (first)

import fit_common as FC  # noqa: E402
import fit_part2_stats as P2  # noqa: E402
import numpy as np  # noqa: E402
import order_logs as OL  # noqa: E402
import stallhit_judge as SJ  # noqa: E402

import gradient_common as G  # noqa: E402

import stats as ST  # noqa: E402

_HELD_AR3 = (
    SC.RCM.OC.WT / "agent-af8961874e131bda2/icenine_py/benchmarks/autoresearch/ar3/held/cache"
)
FC.STORED_S8.setdefault("held", _HELD_AR3 / "S8_held")
FC.STORED_PASS.setdefault("held", _HELD_AR3 / "passes_held")
# held2 has no stored runs: its S8 and G are runs of this step, in this step's cache
FC.STORED_S8.setdefault("held2", SC.SCACHE / "S8_held2")
FC.STORED_PASS.setdefault("held2", SC.SCACHE / "passes_held2")
SAMPLE = FC.SAMPLE
N = 24570
SEED_TEST = 0.80


def read_gz(p: Any) -> Any:
    with gzip.open(p, "rt") as f:
        return json.load(f)


def same_R(a: Any, b: Any) -> bool:
    x = np.round(np.asarray(a, np.float64).ravel(), 12)
    y = np.round(np.asarray(b, np.float64).ravel(), 12)
    return bool(np.max(np.abs(x - y)) <= 1e-9)


def unmatched_side(
    keys: List[Any], idx: Dict[Any, Dict[str, Any]], err_other: Dict[str, np.ndarray], T: Any
) -> List[Dict[str, Any]]:
    """One row per unmatched step of one arm: its taken end and the other arm's map there."""
    rows: List[Dict[str, Any]] = []
    for k in keys:
        e = idx[k]
        v = int(k[0])
        row: Dict[str, Any] = dict(voxel=v, occurrence=int(k[1]), reseed=bool(e.get("reseed")))
        if not e["fits"]:
            row.update(has_fits=False)
            rows.append(row)
            continue
        t = OL.taken_end(e)
        R = np.asarray(t["R"], np.float64).reshape(1, 3, 3)
        err = float(T.err(R, [v])[0])
        row.update(
            has_fits=True,
            hit=float(t["hit"]),
            cost=float(t["cost"]),
            passes_seed_test=bool(t["hit"] >= SEED_TEST),
            err_deg=err,
            wrong=bool(err > G.WRONG),
            other_map_err_before_deg=float(err_other["before"][v]),
            other_map_err_after_deg=float(err_other["after"][v]),
            other_map_wrong_before=bool(err_other["before"][v] > G.WRONG),
            other_map_wrong_after=bool(err_other["after"][v] > G.WRONG),
            s8_map_err_after_deg=float(err_other["s8_after"][v]),
        )
        rows.append(row)
    return rows


def main() -> None:
    sd = SC.SCACHE / "stats"
    sd.mkdir(parents=True, exist_ok=True)
    FC.FIT = sd
    per: Dict[str, Any] = {}
    for arm in ("G", "H"):
        sys.argv = ["fit_part2_stats.py", "--name", arm]
        P2.main()
        d = json.loads((sd / f"part2_{arm}_{SAMPLE}.json").read_text())
        assert not d["unfinished"], arm
        per[arm] = d
    T = G.Truth()
    gd = SC.SCACHE / f"G_{SAMPLE}"
    hdir = SC.SCACHE / f"H_{SAMPLE}"
    rg = json.loads((gd / "rescue.json").read_text())
    rh = json.loads((hdir / "rescue.json").read_text())
    for r in (rg, rh):
        r.pop("rule", None)
    stl = json.loads((hdir / "stallhit.json").read_text())
    stl.pop("rule", None)
    order_h = json.loads((hdir / "order.json").read_text())
    st = dict(per["H"]["stages"])
    st["rescue"] = int(rh["rescue_evals"])
    st["grid_fallback"] = int(st["grid_fallback"]) - st["rescue"]
    assert st["grid_fallback"] == int(rh["grid_evals"]), "grid evals do not match the log"
    per["H"]["stages_split"] = st
    stg = dict(per["G"]["stages"])
    stg["rescue"] = int(rg["rescue_evals"])
    stg["grid_fallback"] = int(stg["grid_fallback"]) - stg["rescue"]
    per["G"]["stages_split"] = stg

    lg = read_gz(gd / "e1_log.json.gz")
    lh = read_gz(hdir / "e1_log.json.gz")
    th = read_gz(hdir / "mc_traces.json.gz")
    ig, ihx = OL.index(lg), OL.index(lh)
    row_of = {id(e): i for i, e in enumerate(lh)}
    by_key = {(int(t["step"]) - 1, int(t["j"])): t for t in th}
    cut_keys: Set[Tuple[int, int]] = {k for k, t in by_key.items() if t.get("cut")}
    # stall checks: [n, hit ratio of the running best, stopped]
    checks = [c for t in th for c in t.get("checks", [])]
    stops = [c for c in checks if c[2]]
    went_on = [c for c in checks if not c[2]]
    cut_hits = [float(c[1]) for c in stops]
    keys = sorted(set(ig) & set(ihx))
    fits_cmp = fits_diff = fits_diff_cut = fits_diff_checks = fits_diff_other = 0
    only_one_side = 0
    changed: List[Dict[str, Any]] = []
    wrong_g_n = wrong_h_n = extra_wrong = 0
    Rg_all: List[np.ndarray] = []
    Rh_all: List[np.ndarray] = []
    vox: List[int] = []
    meta: List[Tuple[Dict[str, Any], Dict[str, Any]]] = []
    for k in keys:
        eg, eh = ig[k], ihx[k]
        if not eg["fits"] or not eh["fits"]:
            continue
        ri = row_of[id(eh)]
        fg = {(f[1], f[2]): f for f in eg["fits"]}
        fh = {(f[1], f[2]): f for f in eh["fits"]}
        only_one_side += len(set(fg) ^ set(fh))
        for key in set(fg) & set(fh):
            fits_cmp += 1
            x, y = fg[key], fh[key]
            if (x[3], x[4], x[5], x[6]) != (y[3], y[4], y[5], y[6]):
                fits_diff += 1
                tr = by_key.get((ri, key[0])) if key[1] == "M" else None
                if tr is not None and tr.get("cut"):
                    fits_diff_cut += 1
                elif (
                    tr is not None
                    and (x[3], x[4], x[6]) == (y[3], y[4], y[6])
                    and int(y[5]) == int(x[5]) + int(tr.get("continued", 0))
                ):
                    fits_diff_checks += 1
                else:
                    fits_diff_other += 1
        tg, thx = OL.taken_end(eg), OL.taken_end(eh)
        Rg_all.append(np.asarray(tg["R"], np.float64).reshape(3, 3))
        Rh_all.append(np.asarray(thx["R"], np.float64).reshape(3, 3))
        vox.append(int(k[0]))
        meta.append((tg, thx))
    Rg_a, Rh_a = np.stack(Rg_all), np.stack(Rh_all)
    err_g = np.array([float(T.err(Rg_a[i][None], [vox[i]])[0]) for i in range(len(vox))])
    err_h = np.array([float(T.err(Rh_a[i][None], [vox[i]])[0]) for i in range(len(vox))])
    mis = np.asarray(T.mis(Rg_a, Rh_a), dtype=float)
    for i, (tg, thx) in enumerate(meta):
        wg, wh = err_g[i] > G.WRONG, err_h[i] > G.WRONG
        wrong_g_n += int(wg)
        wrong_h_n += int(wh)
        extra_wrong += int(wh and not wg)
        if not same_R(tg["R"], thx["R"]):
            changed.append(
                dict(
                    voxel=vox[i],
                    stopped_g=bool(tg["hit"] >= 0.95),
                    hit_g=float(tg["hit"]),
                    hit_h=float(thx["hit"]),
                    cost_g=float(tg["cost"]),
                    cost_h=float(thx["cost"]),
                    err_g=float(err_g[i]),
                    err_h=float(err_h[i]),
                    mis_to_g_deg=float(mis[i]),
                )
            )
    # voxels, H against G and against S8
    e_g = dict(np.load(sd / f"part2_G_{SAMPLE}_err.npz"))
    e_h = dict(np.load(sd / f"part2_H_{SAMPLE}_err.npz"))
    pair: Dict[str, Any] = {}
    for tag in ("before", "after"):
        wg, wh = e_g[tag] > G.WRONG, e_h[tag] > G.WRONG
        fixed, broken = ST.paired_discordant(wg, wh)  # fixed: wrong in G, right in H
        pair[tag] = dict(fixed=int(fixed), broken=int(broken), both_wrong=int((wg & wh).sum()))
    # unmatched steps (seed steps only one arm has)
    only_h = sorted(set(ihx) - set(ig))
    only_g = sorted(set(ig) - set(ihx))
    uh = unmatched_side(only_h, ihx, e_g, T)
    ug = unmatched_side(only_g, ig, e_h, T)
    cut_rows = {k[0] for k in cut_keys}
    for k, r in zip(only_h, uh):
        r["cut_in_step"] = bool(row_of[id(ihx[k])] in cut_rows)
    h_wrong_g_right = [r for r in uh if r.get("wrong") and not r["other_map_wrong_before"]]
    g_wrong_h_right = [r for r in ug if r.get("wrong") and not r["other_map_wrong_before"]]
    unmatched = dict(
        sample=SAMPLE,
        matched_steps=len(keys),
        seed_steps_g=len(lg),
        seed_steps_h=len(lh),
        only_h_n=len(only_h),
        only_g_n=len(only_g),
        only_h_wrong=sum(1 for r in uh if r.get("wrong")),
        only_g_wrong=sum(1 for r in ug if r.get("wrong")),
        only_h_wrong_where_g_map_right=len(h_wrong_g_right),
        only_g_wrong_where_h_map_right=len(g_wrong_h_right),
        only_h_with_a_cut=sum(1 for r in uh if r.get("cut_in_step")),
        only_h_wrong_where_g_map_right_voxels=[r["voxel"] for r in h_wrong_g_right],
        note=(
            "unmatched = a seed step (voxel, occurrence) that only one arm has; the other arm's "
            "map is read before the pass for the counts; rule (4) adds them to the matched steps"
        ),
        only_h=uh,
        only_g=ug,
    )
    (SC.SH / f"unmatched_{SAMPLE}.json").write_text(
        json.dumps(unmatched, indent=1, default=float) + "\n"
    )
    broad = extra_wrong + len(h_wrong_g_right)
    hs, hg = per["H"], per["G"]
    house_h = SJ.house(
        hs["stages"]["total"],
        hs["stages_s8"]["total"],
        hs["paired_after"]["broken"],
        hs["paired_after"]["fixed"],
        hs["after"]["median_err_deg"],
        hs["s8_after"]["median_err_deg"],
    )
    house_g = SJ.house(
        hg["stages"]["total"],
        hg["stages_s8"]["total"],
        hg["paired_after"]["broken"],
        hg["paired_after"]["fixed"],
        hg["after"]["median_err_deg"],
        hg["s8_after"]["median_err_deg"],
    )
    rules = SJ.rules(
        SAMPLE,
        hs["stages"]["total"],
        hg["stages"]["total"],
        hs["stages"]["seed_step"],
        hg["stages"]["seed_step"],
        pair["after"]["broken"] - pair["after"]["fixed"],
        hs["after"]["median_err_deg"],
        hg["after"]["median_err_deg"],
        broad,
        int(rh["grid_fallbacks_left"]),
        int(rg["grid_fallbacks_left"]),
        house_h,
        house_g,
    )
    validity = dict(
        h_finished=bool(not (hdir / "capped.json").exists() and (hdir / "run.json").exists()),
        n_processed=int(json.loads((hdir / "run.json").read_text())["n_processed"]),
        aligned_draws=bool(order_h["aligned_draws"]),
        seed_steps_g=len(lg),
        seed_steps_h=len(lh),
        matched_steps=len(keys),
        fits_compared=fits_cmp,
        fits_only_one_side=only_one_side,
        fits_differ=fits_diff,
        fits_differ_cut_mc=fits_diff_cut,
        fits_differ_checks_only=fits_diff_checks,
        fits_differ_other=fits_diff_other,
        every_stop_at_095_or_more=bool(all(h >= 0.95 for h in cut_hits)),
    )
    out = dict(
        sample=SAMPLE,
        arms=per,
        rescue=dict(G=rg, H=rh),
        cuts=dict(
            **stl,
            refinements_in_traces=len(th),
            cut_refinements=len(cut_keys),
            stall_checks=len(checks),
            stall_checks_stopped=len(stops),
            stall_checks_went_on=len(went_on),
            went_on_hit_min=float(min(c[1] for c in went_on)) if went_on else None,
            went_on_hit_median=float(np.median([c[1] for c in went_on])) if went_on else None,
            went_on_hit_max=float(max(c[1] for c in went_on)) if went_on else None,
            stop_hit_min=float(min(cut_hits)) if cut_hits else None,
        ),
        taken_changed=changed,
        taken_changed_n=len(changed),
        taken_changed_median_mis_deg=(
            float(np.median([c["mis_to_g_deg"] for c in changed])) if changed else None
        ),
        wrong_seeds=dict(
            g=wrong_g_n,
            h=wrong_h_n,
            h_wrong_where_g_right=extra_wrong,
            unmatched_h_wrong_where_g_map_right=len(h_wrong_g_right),
            broad=broad,
        ),
        h_vs_g=dict(voxels=pair),
        validity=validity,
        house_h=house_h,
        house_g=house_g,
        rules=rules,
    )
    (SC.SH / f"stats_{SAMPLE}.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(json.dumps(dict(sample=SAMPLE, rules=rules, pair=pair, validity=validity), indent=1))


if __name__ == "__main__":
    main()
