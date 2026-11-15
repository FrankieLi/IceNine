"""AR5 step 14, phase A step 6: the numbers of the finished arm J run of one sample, against arm H.

  RESEED_SAMPLE=<held|mid|grad|big|held2> POLISH_CACHE=<dir> uv run python \
      scripts/autoresearch/ar5/polish_stats.py

Arm ``H`` is step 12's stored run (linked per entry, ``polish_link_cache.sh`` and
``polish_link_stats.sh``); arm ``J`` is this step's run and pass. For ``big``, ``grad``, ``mid``
and ``held``: step 2's ``fit_part2_stats`` for both arms (each against the stored ``S8``); the
stage split (``G1`` rescue and ``G2`` evaluations apart from the grid searches); the seed steps
matched with ``H``'s (same voxel, same occurrence: the fits must be equal, the seed step being
deterministic and independent of the rest of the run); the seed steps only one arm has
(``unmatched``); wrong seeds, matched and unmatched; grid searches left; every ``G2`` end taken
(error, angle to ``H``'s grid end); the grid searches both arms ran (the same voxel: the same
evaluations); fixed and broken of ``J`` against ``H`` and ``S8``; the rules (``polish_judge``).
``held2``: ONLY what the plan allows before step 13's ``held2`` result is on record: totals before
the pass, grid searches, ``G2`` ends, the maps before the pass (error against the truth, paired
against ``H``'s), matched and unmatched seed steps. Nothing under any ``passes_held2`` folder is
opened (the guards refuse it), so no ``held2`` number after the pass, and no ``held2`` rule, is
computed here. Truth is used for scoring only; no p-values. Writes ``polish/stats_<sample>.json``
and ``polish/unmatched_<sample>.json``.
"""

import gzip
import json
import sys
from typing import Any, Dict, List

import polish_common as PC  # noqa: E402  (first)

import fit_common as FC  # noqa: E402
import gradient_common as G  # noqa: E402
import numpy as np  # noqa: E402
import order_logs as OL  # noqa: E402
import polish_judge as PJ  # noqa: E402
import reseed_common as RC  # noqa: E402
import stats as ST  # noqa: E402

N = 24570
SEED_TEST = 0.80
SAMPLE = FC.SAMPLE

# stored S8 of the held sample (as step 12's stats do)
_HELD_AR3 = (
    PC.RCM.OC.WT / "agent-af8961874e131bda2/icenine_py/benchmarks/autoresearch/ar3/held/cache"
)
FC.STORED_S8.setdefault("held", _HELD_AR3 / "S8_held")
FC.STORED_PASS.setdefault("held", _HELD_AR3 / "passes_held")


def same_R(a: Any, b: Any) -> bool:
    x = np.round(np.asarray(a, np.float64).ravel(), 12)
    y = np.round(np.asarray(b, np.float64).ravel(), 12)
    return bool(np.max(np.abs(x - y)) <= 1e-9)


def read_gz(p: Any) -> Any:
    with gzip.open(PC.PG.guarded(p), "rt") as f:
        return json.load(f)


def unmatched_side(
    keys: List[Any], idx: Dict[Any, Dict[str, Any]], err_other: Dict[str, np.ndarray], T: Any
) -> List[Dict[str, Any]]:
    """One row per unmatched step of one arm: its taken end and the other arm's map there (before
    the pass; after the pass only where the sample allows it)."""
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
        err = float(T.err(np.asarray(t["R"], np.float64).reshape(1, 3, 3), [v])[0])
        row.update(
            has_fits=True,
            hit=float(t["hit"]),
            cost=float(t["cost"]),
            passes_seed_test=bool(t["hit"] >= SEED_TEST),
            err_deg=err,
            wrong=bool(err > G.WRONG),
            other_map_err_before_deg=float(err_other["before"][v]),
            other_map_wrong_before=bool(err_other["before"][v] > G.WRONG),
        )
        if "after" in err_other:
            row.update(
                other_map_err_after_deg=float(err_other["after"][v]),
                other_map_wrong_after=bool(err_other["after"][v] > G.WRONG),
            )
        rows.append(row)
    return rows


def seed_steps(
    ig: Dict[Any, Any],
    ij: Dict[Any, Any],
    T: Any,
    err_h: Dict[str, np.ndarray],
    err_j: Dict[str, np.ndarray],
) -> Dict[str, Any]:
    """Matched steps (fits must be equal), unmatched steps, wrong seeds."""
    keys = sorted(set(ig) & set(ij))
    fits_cmp = fits_diff = 0
    taken_changed = extra_wrong = wrong_h = wrong_j = 0
    for k in keys:
        eh, ej = ig[k], ij[k]
        if not eh["fits"] or not ej["fits"]:
            continue
        fh = {(f[1], f[2]): f for f in eh["fits"]}
        fj = {(f[1], f[2]): f for f in ej["fits"]}
        for key in set(fh) & set(fj):
            fits_cmp += 1
            x, y = fh[key], fj[key]
            if (x[3], x[4], x[5], x[6]) != (y[3], y[4], y[5], y[6]):
                fits_diff += 1
        th, tj = OL.taken_end(eh), OL.taken_end(ej)
        v = int(k[0])
        e1 = float(T.err(np.asarray(th["R"], np.float64).reshape(1, 3, 3), [v])[0])
        e2 = float(T.err(np.asarray(tj["R"], np.float64).reshape(1, 3, 3), [v])[0])
        wrong_h += int(e1 > G.WRONG)
        wrong_j += int(e2 > G.WRONG)
        extra_wrong += int(e2 > G.WRONG and not e1 > G.WRONG)
        taken_changed += int(not same_R(th["R"], tj["R"]))
    only_j = sorted(set(ij) - set(ig))
    only_h = sorted(set(ig) - set(ij))
    uj = unmatched_side(only_j, ij, err_h, T)
    uh = unmatched_side(only_h, ig, err_j, T)
    j_wrong_h_right = [r for r in uj if r.get("wrong") and not r["other_map_wrong_before"]]
    h_wrong_j_right = [r for r in uh if r.get("wrong") and not r["other_map_wrong_before"]]
    return dict(
        matched_steps=len(keys),
        fits_compared=fits_cmp,
        fits_differ=fits_diff,
        taken_changed=taken_changed,
        wrong_seeds=dict(
            h=wrong_h,
            j=wrong_j,
            j_wrong_where_h_right=extra_wrong,
            unmatched_j_wrong_where_h_map_right=len(j_wrong_h_right),
            broad=extra_wrong + len(j_wrong_h_right),
        ),
        unmatched=dict(
            only_j_n=len(only_j),
            only_h_n=len(only_h),
            only_j_wrong=sum(1 for r in uj if r.get("wrong")),
            only_h_wrong=sum(1 for r in uh if r.get("wrong")),
            only_j_wrong_where_h_map_right=len(j_wrong_h_right),
            only_h_wrong_where_j_map_right=len(h_wrong_j_right),
            only_j_wrong_where_h_map_right_voxels=[r["voxel"] for r in j_wrong_h_right],
            only_j=uj,
            only_h=uh,
        ),
    )


def polish_rows(T: Any, R_h_before: np.ndarray, jlog: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Every failed rescue of J (G2 accepted or not): evaluations, the end taken, its error and its
    angle to H's map at the voxel before the pass (H's grid end, where H ran a grid search)."""
    rows = []
    for r in jlog:
        if not r["attempted"] or r["rescued"]:
            continue
        v = int(r["voxel"])
        row: Dict[str, Any] = dict(
            voxel=v,
            polished=bool(r["polished"]),
            rescue_evals=int(r["rescue_evals"]),
            polish_evals=int(r["polish_evals"]),
            grid_evals=int(r["grid_evals"]),
            n_polish_tries=len(r["polish_tries"]),
        )
        if r["polished"]:
            R = np.asarray(r["polish_tries"][-1]["R"], np.float64).reshape(3, 3)
            row.update(
                hit=float(r["polish_tries"][-1]["hit"]),
                err_deg=float(T.err(R[None], [v])[0]),
                angle_to_h_map_deg=float(T.mis(R[None], R_h_before[v][None])[0]),
            )
        rows.append(row)
    return rows


def main_full() -> None:
    sd = PC.PCACHE / "stats"
    sd.mkdir(parents=True, exist_ok=True)
    FC.FIT = sd
    per: Dict[str, Any] = {}
    for arm in ("H", "J"):
        sys.argv = ["fit_part2_stats.py", "--name", arm]
        import fit_part2_stats as P2

        P2.main()
        d = json.loads((sd / f"part2_{arm}_{SAMPLE}.json").read_text())
        assert not d["unfinished"], arm
        per[arm] = d
    T = G.Truth()
    hdir, jdir = PC.PCACHE / f"H_{SAMPLE}", PC.PCACHE / f"J_{SAMPLE}"
    rh = json.loads((hdir / "rescue.json").read_text())
    pj = json.loads((jdir / "polish.json").read_text())
    for r in (rh, pj):
        r.pop("rule", None)
    order_j = json.loads((jdir / "order.json").read_text())
    st = dict(per["J"]["stages"])
    st["rescue"] = int(pj["rescue_evals"])
    st["polish"] = int(pj["polish_evals"])
    st["grid_fallback"] = int(st["grid_fallback"]) - st["rescue"] - st["polish"]
    assert st["grid_fallback"] == int(pj["grid_evals"]), "grid evals do not match the log"
    per["J"]["stages_split"] = st
    sh = dict(per["H"]["stages"])
    sh["rescue"] = int(rh["rescue_evals"])
    sh["grid_fallback"] = int(sh["grid_fallback"]) - sh["rescue"]
    assert sh["grid_fallback"] == int(rh["grid_evals"]), "H grid evals do not match its log"
    per["H"]["stages_split"] = sh

    lh = read_gz(hdir / "e1_log.json.gz")
    lj = read_gz(jdir / "e1_log.json.gz")
    e_h = dict(np.load(sd / f"part2_H_{SAMPLE}_err.npz"))
    e_j = dict(np.load(sd / f"part2_J_{SAMPLE}_err.npz"))
    ss = seed_steps(OL.index(lh), OL.index(lj), T, e_h, e_j)
    pair: Dict[str, Any] = {}
    for tag in ("before", "after"):
        wh, wj = e_h[tag] > G.WRONG, e_j[tag] > G.WRONG
        fixed, broken = ST.paired_discordant(wh, wj)  # fixed: wrong in H, right in J
        pair[tag] = dict(fixed=int(fixed), broken=int(broken), both_wrong=int((wh & wj).sum()))
    R_h, _ = G.read_map(PC.PG.guarded(hdir / "recon.mic.gz"), N)
    jlog = read_gz(jdir / "polish_log.json.gz")
    hlog = read_gz(hdir / "rescue_log.json.gz")
    prow = polish_rows(T, R_h, jlog)
    # grid searches both arms ran at the same voxel: the same evaluations (aligned draws)
    hg = {int(r["voxel"]): int(r["grid_evals"]) for r in hlog if r["grid_evals"] > 0}
    jg = {int(r["voxel"]): int(r["grid_evals"]) for r in jlog if r["grid_evals"] > 0}
    both = sorted(set(hg) & set(jg))
    grid_same = sum(1 for v in both if hg[v] == jg[v])
    sj, sh_ = per["J"], per["H"]
    house_j = PJ.house(
        sj["stages"]["total"],
        sj["stages_s8"]["total"],
        sj["paired_after"]["broken"],
        sj["paired_after"]["fixed"],
        sj["after"]["median_err_deg"],
        sj["s8_after"]["median_err_deg"],
    )
    house_h = PJ.house(
        sh_["stages"]["total"],
        sh_["stages_s8"]["total"],
        sh_["paired_after"]["broken"],
        sh_["paired_after"]["fixed"],
        sh_["after"]["median_err_deg"],
        sh_["s8_after"]["median_err_deg"],
    )
    rules = PJ.rules(
        SAMPLE,
        sj["stages"]["total"],
        sh_["stages"]["total"],
        pair["after"]["broken"] - pair["after"]["fixed"],
        sj["after"]["median_err_deg"],
        sh_["after"]["median_err_deg"],
        ss["wrong_seeds"]["broad"],
        house_j,
        house_h,
    )
    identical_map = None
    if SAMPLE == "big":
        import hashlib

        def sha(p: Any) -> str:
            with gzip.open(p) as f:
                return hashlib.sha256(f.read()).hexdigest()

        identical_map = bool(sha(jdir / "recon.mic.gz") == sha(hdir / "recon.mic.gz"))
    validity = dict(
        j_finished=bool(not (jdir / "capped.json").exists() and (jdir / "run.json").exists()),
        n_processed=int(json.loads((jdir / "run.json").read_text())["n_processed"]),
        aligned_draws=bool(order_j["aligned_draws"]),
        seed_steps_h=len(lh),
        seed_steps_j=len(lj),
        matched_steps=ss["matched_steps"],
        fits_compared=ss["fits_compared"],
        fits_differ=ss["fits_differ"],
        grid_searches_both_ran=len(both),
        grid_searches_same_evals=grid_same,
        big_map_identical_to_h=identical_map,
        big_evals_identical_to_h=bool(sj["stages"]["total"] == sh_["stages"]["total"]),
    )
    out = dict(
        sample=SAMPLE,
        arms=per,
        rescue=dict(H=rh, J=pj),
        g2_cases=prow,
        g2_accepted=sum(1 for r in prow if r["polished"]),
        g2_failed_rescues=len(prow),
        seed_steps=ss,
        j_vs_h=dict(voxels=pair),
        validity=validity,
        house_j=house_j,
        house_h=house_h,
        rules=rules,
    )
    (PC.PO / f"stats_{SAMPLE}.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    (PC.PO / f"unmatched_{SAMPLE}.json").write_text(
        json.dumps(ss["unmatched"], indent=1, default=float) + "\n"
    )
    print(json.dumps(dict(sample=SAMPLE, rules=rules, pair=pair, validity=validity), indent=1))


def main_held2() -> None:
    """Before the pass only (see the module docstring)."""
    T = G.Truth()
    hdir, jdir = PC.PCACHE / f"H_{SAMPLE}", PC.PCACHE / f"J_{SAMPLE}"
    rh = json.loads((hdir / "rescue.json").read_text())
    pj = json.loads((jdir / "polish.json").read_text())
    for r in (rh, pj):
        r.pop("rule", None)
    run_h = json.loads((hdir / "run.json").read_text())
    run_j = json.loads((jdir / "run.json").read_text())
    order_j = json.loads((jdir / "order.json").read_text())
    st = dict(
        before_pass=int(run_j["evals_total"]),
        rescue=int(pj["rescue_evals"]),
        polish=int(pj["polish_evals"]),
        grid_fallback_total=int(run_j["stats"].get("index_fallback_grid_evals", 0)),
    )
    st["grid_fallback"] = st["grid_fallback_total"] - st["rescue"] - st["polish"]
    assert st["grid_fallback"] == int(pj["grid_evals"]), "grid evals do not match the log"
    R_h, _ = G.read_map(PC.PG.guarded(hdir / "recon.mic.gz"), N)
    R_j, _ = G.read_map(PC.PG.guarded(jdir / "recon.mic.gz"), N)
    eh, ej = np.asarray(T.err(R_h), float), np.asarray(T.err(R_j), float)
    wh, wj = eh > G.WRONG, ej > G.WRONG
    fixed, broken = ST.paired_discordant(wh, wj)
    lh = read_gz(hdir / "e1_log.json.gz")
    lj = read_gz(jdir / "e1_log.json.gz")
    ss = seed_steps(OL.index(lh), OL.index(lj), T, dict(before=eh), dict(before=ej))
    jlog = read_gz(jdir / "polish_log.json.gz")
    hlog = read_gz(hdir / "rescue_log.json.gz")
    prow = polish_rows(T, R_h, jlog)
    hg = {int(r["voxel"]): int(r["grid_evals"]) for r in hlog if r["grid_evals"] > 0}
    jg = {int(r["voxel"]): int(r["grid_evals"]) for r in jlog if r["grid_evals"] > 0}
    both = sorted(set(hg) & set(jg))
    lo_h, hi_h = ST.wilson(int(wh.sum()), len(wh))
    lo_j, hi_j = ST.wilson(int(wj.sum()), len(wj))
    out = dict(
        sample=SAMPLE,
        scope=(
            "before the boundary pass only; nothing under any passes_held2 folder was opened; "
            "the totals with the pass, the maps after the pass and the rules wait for step 13's "
            "held2 result (plan step 7)"
        ),
        evals_before_pass=dict(j=int(run_j["evals_total"]), h=int(run_h["evals_total"])),
        cut_before_pass_pct=PJ.pct_fewer(int(run_j["evals_total"]), int(run_h["evals_total"])),
        stages_j=st,
        rescue=dict(H=rh, J=pj),
        g2_cases=prow,
        g2_accepted=sum(1 for r in prow if r["polished"]),
        g2_failed_rescues=len(prow),
        before_pass=dict(
            h=dict(
                wrong=int(wh.sum()), wilson_lo=lo_h, wilson_hi=hi_h, median=float(np.median(eh))
            ),
            j=dict(
                wrong=int(wj.sum()), wilson_lo=lo_j, wilson_hi=hi_j, median=float(np.median(ej))
            ),
            fixed=int(fixed),
            broken=int(broken),
            both_wrong=int((wh & wj).sum()),
        ),
        seed_steps=ss,
        validity=dict(
            j_finished=bool(not (jdir / "capped.json").exists() and (jdir / "run.json").exists()),
            n_processed=int(run_j["n_processed"]),
            aligned_draws=bool(order_j["aligned_draws"]),
            seed_steps_h=len(lh),
            seed_steps_j=len(lj),
            matched_steps=ss["matched_steps"],
            fits_compared=ss["fits_compared"],
            fits_differ=ss["fits_differ"],
            grid_searches_both_ran=len(both),
            grid_searches_same_evals=sum(1 for v in both if hg[v] == jg[v]),
        ),
    )
    (PC.PO / f"stats_{SAMPLE}.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    (PC.PO / f"unmatched_{SAMPLE}.json").write_text(
        json.dumps(ss["unmatched"], indent=1, default=float) + "\n"
    )
    print(
        json.dumps({k: out[k] for k in ("evals_before_pass", "before_pass", "validity")}, indent=1)
    )


if __name__ == "__main__":
    RC.SAMPLE  # noqa: B018
    main_held2() if SAMPLE == "held2" else main_full()
