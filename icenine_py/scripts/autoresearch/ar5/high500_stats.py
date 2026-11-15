"""AR5 step 16, phase A step 6: the numbers of the finished arm J5 run of one sample, against arm J
(aligned draws) and against S8, all with the boundary pass.

  RESEED_SAMPLE=<held3|held|mid|grad|big> HIGH500_CACHE=<dir> uv run python \
      scripts/autoresearch/ar5/high500_stats.py

Arms: ``S8`` (stored for held, mid, grad and big; run here for held3), ``J`` (step 14's stored run
for the disclosed samples, linked per file; run here for held3) and ``J5`` (this step). Reads
``<cache>/<arm>_<s>/`` and ``<cache>/passes_<s>/``. Writes ``ar5/high500/stats_<s>.json`` and
``unmatched_<s>.json``: evaluations by stage, stall checks and stops of the MC refinements (from
``mc_traces``), grid searches, wrong voxels before and after the pass (Wilson 95%), median error,
fixed and broken of J5 against J and against S8 and of J against S8 (count by voxel; fixed: wrong
in the older arm, right in the newer), the seed steps matched with J's (same voxel, same
occurrence: the fits must be equal unless a stop changed an MC end), wrong seeds matched and
unmatched, and the rules (``high500_judge``). ``held3`` is read only because
``<cache>/held3_runs_done.json`` exists (``high500_guard``). Truth for scoring only; no p-values.
"""

import json
from typing import Any, Dict, List

import high500_common as BC  # noqa: E402  (first)

import fit_common as FC  # noqa: E402
import fit_part2_stats as P2  # noqa: E402
import gradient_common as G  # noqa: E402
import high500_judge as HJ  # noqa: E402
import high500_rules as HR  # noqa: E402
import numpy as np  # noqa: E402
import order_logs as OL  # noqa: E402
import stats as ST  # noqa: E402

N = 24570
SEED_TEST = 0.80
SAMPLE = FC.SAMPLE
CACHE = BC.HCACHE
_HELD_AR3 = (
    BC.ICE.parents[1] / "agent-af8961874e131bda2/icenine_py/benchmarks/autoresearch/ar3/held/cache"
)


def s8_dirs(s: str) -> Any:
    if s == "held3":
        return CACHE / "S8_held3", CACHE / "passes_held3"
    if s == "held":
        return _HELD_AR3 / "S8_held", _HELD_AR3 / "passes_held"
    return FC.STORED_S8[s], FC.STORED_PASS[s]


def same_R(a: Any, b: Any) -> bool:
    x = np.round(np.asarray(a, np.float64).ravel(), 12)
    y = np.round(np.asarray(b, np.float64).ravel(), 12)
    return bool(np.max(np.abs(x - y)) <= 1e-9)


def arm(name: str, rd: Any, pd: Any, T: Any) -> Dict[str, Any]:
    """One arm's evaluations by stage and its error vectors (before and after the pass)."""
    run = json.loads((rd / "run.json").read_text())
    assert not (rd / "capped.json").exists(), f"{name} reached its cap"
    pj = json.loads((pd / f"D4_{name}" / "pass.json").read_text())
    Rb, _ = G.read_map(rd / "recon.mic.gz", N)
    Ra = np.load(pd / f"after_{name}.npy")
    eb, ea = np.asarray(T.err(Rb), float), np.asarray(T.err(Ra), float)
    return dict(
        run=run,
        stages=P2.stage_evals(run, int(pj["evals"])),
        grid_searches=int(run["stats"].get("index_fallback_seed_test", 0)),
        before=P2.side(eb),
        after=P2.side(ea),
        eb=eb,
        ea=ea,
        n_processed=int(run["n_processed"]),
    )


def pair(old: Dict[str, Any], new: Dict[str, Any]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for tag in ("before", "after"):
        wo, wn = old["e" + tag[0]] > G.WRONG, new["e" + tag[0]] > G.WRONG
        fixed, broken = ST.paired_discordant(wo, wn)  # fixed: wrong in old, right in new
        out[tag] = dict(fixed=int(fixed), broken=int(broken), both_wrong=int((wo & wn).sum()))
    return out


def traces(rd: Any) -> Dict[str, int]:
    tr = BC.load_gz(rd / "mc_traces.json.gz")
    return dict(
        mc_refinements=len(tr),
        mc_evals=int(sum(int(t["evals"]) for t in tr)),
        stops=int(sum(1 for t in tr if t["cut"])),
        checks=int(sum(len(t.get("checks", [])) for t in tr)),
        checks_went_on=int(sum(int(t.get("continued", 0)) for t in tr)),
    )


def unmatched_side(
    keys: List[Any], idx: Dict[Any, Dict[str, Any]], err_other: np.ndarray, T: Any
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for k in keys:
        e = idx[k]
        v = int(k[0])
        row: Dict[str, Any] = dict(voxel=v, occurrence=int(k[1]), reseed=bool(e.get("reseed")))
        if not e["fits"]:
            rows.append(dict(row, has_fits=False))
            continue
        t = OL.taken_end(e)
        err = float(T.err(np.asarray(t["R"], np.float64).reshape(1, 3, 3), [v])[0])
        rows.append(
            dict(
                row,
                has_fits=True,
                hit=float(t["hit"]),
                passes_seed_test=bool(t["hit"] >= SEED_TEST),
                err_deg=err,
                wrong=bool(err > G.WRONG),
                other_map_err_before_deg=float(err_other[v]),
                other_map_wrong_before=bool(err_other[v] > G.WRONG),
            )
        )
    return rows


def seed_steps(
    ij: Dict[Any, Any], i5: Dict[Any, Any], T: Any, e_j: np.ndarray, e_5: np.ndarray
) -> Dict[str, Any]:
    keys = sorted(set(ij) & set(i5))
    fits_cmp = fits_diff = taken_changed = extra_wrong = wrong_j = wrong_5 = 0
    mc_changed = 0
    for k in keys:
        a, b = ij[k], i5[k]
        if not a["fits"] or not b["fits"]:
            continue
        fa = {(f[1], f[2]): f for f in a["fits"]}
        fb = {(f[1], f[2]): f for f in b["fits"]}
        for key in set(fa) & set(fb):
            fits_cmp += 1
            x, y = fa[key], fb[key]
            if (x[3], x[4], x[5], x[6]) != (y[3], y[4], y[5], y[6]):
                fits_diff += 1
                mc_changed += int(key[1] == "M")
        ta, tb = OL.taken_end(a), OL.taken_end(b)
        v = int(k[0])
        e1 = float(T.err(np.asarray(ta["R"], np.float64).reshape(1, 3, 3), [v])[0])
        e2 = float(T.err(np.asarray(tb["R"], np.float64).reshape(1, 3, 3), [v])[0])
        wrong_j += int(e1 > G.WRONG)
        wrong_5 += int(e2 > G.WRONG)
        extra_wrong += int(e2 > G.WRONG and not e1 > G.WRONG)
        taken_changed += int(not same_R(ta["R"], tb["R"]))
    only_5, only_j = sorted(set(i5) - set(ij)), sorted(set(ij) - set(i5))
    u5 = unmatched_side(only_5, i5, e_j, T)
    uj = unmatched_side(only_j, ij, e_5, T)
    w5 = [r for r in u5 if r.get("wrong") and not r["other_map_wrong_before"]]
    wj = [r for r in uj if r.get("wrong") and not r["other_map_wrong_before"]]
    return dict(
        matched_steps=len(keys),
        fits_compared=fits_cmp,
        fits_differ=fits_diff,
        mc_fits_differ=mc_changed,
        taken_changed=taken_changed,
        wrong_seeds=dict(
            j=wrong_j,
            j5=wrong_5,
            j5_wrong_where_j_right=extra_wrong,
            unmatched_j5_wrong_where_j_map_right=len(w5),
            broad=extra_wrong + len(w5),
        ),
        unmatched=dict(
            only_j5_n=len(only_5),
            only_j_n=len(only_j),
            only_j5_wrong=sum(1 for r in u5 if r.get("wrong")),
            only_j_wrong=sum(1 for r in uj if r.get("wrong")),
            only_j5_wrong_where_j_map_right=len(w5),
            only_j_wrong_where_j5_map_right=len(wj),
            only_j5=u5,
            only_j=uj,
        ),
    )


def strip(a: Dict[str, Any]) -> Dict[str, Any]:
    return {k: v for k, v in a.items() if k not in ("eb", "ea", "run")}


def main() -> None:
    s = SAMPLE
    T = G.Truth()
    sd, spd = s8_dirs(s)
    rj, rj5 = CACHE / f"J_{s}", CACHE / f"J5_{s}"
    pd = CACHE / f"passes_{s}"
    arms = dict(
        S8=arm("S8", sd, spd, T),
        J=arm("J", rj, pd, T),
        J5=arm("J5", rj5, pd, T),
    )
    # the stored counts of the plan's stop rule (disclosed samples)
    stored_ok = None
    if s in HR.EXPECT:
        stored_ok = bool(arms["J"]["stages"]["total"] == HR.EXPECT[s]["total"])
    tj, t5 = traces(rj), traces(rj5)
    if s in HR.EXPECT:
        stored_ok = bool(stored_ok and tj["mc_evals"] == HR.EXPECT[s]["mc"])
    logs = {n: BC.load_gz(d / "e1_log.json.gz") for n, d in (("J", rj), ("J5", rj5))}
    ss = seed_steps(OL.index(logs["J"]), OL.index(logs["J5"]), T, arms["J"]["eb"], arms["J5"]["eb"])
    pj = {n: json.loads((d / "polish.json").read_text()) for n, d in (("J", rj), ("J5", rj5))}
    for n in pj:
        pj[n].pop("rule", None)
    stages = {n: dict(arms[n]["stages"]) for n in ("J", "J5")}
    for n in ("J", "J5"):
        st = stages[n]
        st["rescue"] = int(pj[n]["rescue_evals"])
        st["polish"] = int(pj[n]["polish_evals"])
        st["grid_fallback"] = int(st["grid_fallback"]) - st["rescue"] - st["polish"]
        assert st["grid_fallback"] == int(pj[n]["grid_evals"]), f"{n}: grid evals vs its log"
    p_j5_j, p_j5_s8, p_j_s8 = (
        pair(arms["J"], arms["J5"]),
        pair(arms["S8"], arms["J5"]),
        pair(arms["S8"], arms["J"]),
    )

    def h(n: str) -> Dict[str, bool]:
        a = arms[n]
        return HJ.house(
            a["stages"]["total"],
            arms["S8"]["stages"]["total"],
            pair(arms["S8"], a)["after"]["broken"],
            pair(arms["S8"], a)["after"]["fixed"],
            a["after"]["median_err_deg"],
            arms["S8"]["after"]["median_err_deg"],
        )

    house_j5, house_j = h("J5"), h("J")
    rules = HJ.rules(
        s,
        arms["J5"]["stages"]["total"],
        arms["J"]["stages"]["total"],
        p_j5_j["after"]["broken"] - p_j5_j["after"]["fixed"],
        arms["J5"]["after"]["median_err_deg"],
        arms["J"]["after"]["median_err_deg"],
        ss["wrong_seeds"]["broad"],
        arms["J5"]["grid_searches"],
        arms["J"]["grid_searches"],
        house_j5,
        house_j,
    )
    validity = dict(
        finished=all(a["n_processed"] == N for a in arms.values()),
        none_capped=all(not (d / "capped.json").exists() for d in (sd, rj, rj5)),
        stored_j_counts_reproduce=stored_ok,
        aligned_draws=bool(json.loads((rj5 / "order.json").read_text())["aligned_draws"]),
        seed_steps_j=len(logs["J"]),
        seed_steps_j5=len(logs["J5"]),
        matched_steps=ss["matched_steps"],
        fits_compared=ss["fits_compared"],
        fits_differ=ss["fits_differ"],
        mc_fits_differ=ss["mc_fits_differ"],
    )
    out = dict(
        sample=s,
        arms={n: strip(a) for n, a in arms.items()},
        stages_split=stages,
        mc=dict(J=tj, J5=t5),
        polish=pj,
        seed_steps=ss,
        j5_vs_j=p_j5_j,
        j5_vs_s8=p_j5_s8,
        j_vs_s8=p_j_s8,
        house_j5=house_j5,
        house_j=house_j,
        rules=rules,
        validity=validity,
    )
    BC.HD.mkdir(parents=True, exist_ok=True)
    (BC.HD / f"stats_{s}.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    (BC.HD / f"unmatched_{s}.json").write_text(
        json.dumps(ss["unmatched"], indent=1, default=float) + "\n"
    )
    print(
        json.dumps(
            dict(sample=s, rules=rules, validity=validity, j5_vs_j=p_j5_j["after"]), indent=1
        )
    )


if __name__ == "__main__":
    main()
