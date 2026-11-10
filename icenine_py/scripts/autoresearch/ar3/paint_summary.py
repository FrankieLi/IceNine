"""AR3 step 1: summary.json, tables.md and the figure of the paint run on the pilot region.

  uv run python scripts/autoresearch/ar3/paint_summary.py

Reads ``benchmarks/autoresearch/ar3/runs/{paint,paint/side,keyon_profile,keyon_timers}`` and the
stored AR2 runs (key off, key on, key on with FusedCost). Nothing is rerun. Wrong = more than 1
degree from the truth (symmetry-reduced). Counts by voxel; intervals are Wilson; no p-values.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import ar3_common as C
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0] / "ar1"))
sys.path.insert(0, str(HERE.parents[0] / "ar2"))
import g2_bfs  # noqa: E402
import key_summary as KS  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

from icenine import ray_indexing as RI  # noqa: E402

import doc_tables as DT  # noqa: E402
import stats as ST  # noqa: E402

OUT = C.OUT
RUN = OUT / "runs" / "paint"
WRONG = 1.0
STORED_KEY_ON_EVALS = 2_819_906
STORED_KEY_ON_WRONG = 22
GATE_BREAKS = 20
MARGINS = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0]


def wil(k: int, n: int) -> str:
    lo, hi = ST.wilson(k, n)
    return f"{k} of {n} ({100 * lo:.1f} to {100 * hi:.1f}%)"


def wil_dict(k: int, n: int) -> Dict[str, Any]:
    lo, hi = ST.wilson(k, n) if n else (float("nan"), float("nan"))
    return dict(k=int(k), n=int(n), lo=float(lo), hi=float(hi))


def paired(ea: np.ndarray, eb: np.ndarray) -> Dict[str, int]:
    wa, wb = ea > WRONG, eb > WRONG
    return dict(
        wrong_a=int(wa.sum()),
        wrong_b=int(wb.sum()),
        fixed=int((wa & ~wb).sum()),
        broken=int((~wa & wb).sum()),
        both_wrong=int((wa & wb).sum()),
    )


def mis(Ra: np.ndarray, Rb: np.ndarray, sym: np.ndarray) -> float:
    return float(RI.misorientation_deg(Ra[None], Rb[None], sym)[0])


def premise(L: Any) -> Dict[str, Any]:
    """Distinct orientations inside each true grain: the largest misorientation of a grain's
    voxels to its first voxel, from the truth map (the region and the whole layer)."""
    from icenine.mic_file import MicFile

    ex = C.EX / "SimInput"
    truth = MicFile.read(str(ex / "rand_500grains_1mm_neworient_s0.mic"))
    R = np.array([v.orientation for v in truth.voxels], dtype=np.float64)
    out: Dict[str, Any] = {}
    for name, sel in (("region", np.sort(L.region)), ("layer", np.arange(len(R)))):
        g = L.gmap[sel]
        worst = {}
        for gid in np.unique(g):
            vs = sel[g == gid]
            d = RI.misorientation_deg(R[vs], np.broadcast_to(R[vs[0]], R[vs].shape), L.sym)
            worst[int(gid)] = float(d.max())
        w = np.array(list(worst.values()))
        out[name] = dict(
            n_grains=len(w),
            n_voxels=len(sel),
            max_spread_deg=float(w.max()),
            n_grains_over_1e3_deg=int((w > 1e-3).sum()),
            n_grains_over_1_deg=int((w > 1.0).sum()),
        )
    return out


def main() -> None:
    L = KS.Layer()
    sym = L.sym
    run = json.loads((RUN / "run.json").read_text())
    paint = json.loads((RUN / "paint.json").read_text())
    E = g2_bfs.errors(RUN)
    ES = g2_bfs.errors(RUN / "side")
    EK = g2_bfs.errors(C.STORED_KEY_ON)
    EF = g2_bfs.errors(C.STORED_KEY_ON_FUSED)
    EO = g2_bfs.errors(KS.OFF)
    n = len(E["err"])
    full = L.region[E["row"]]
    kind = np.array([L.kind(int(f)) for f in full])
    gid = E["grain"]
    Rtrue = np.stack([L.Rg[L._pos[int(g)]] for g in gid])
    wrong = E["err"] > WRONG
    S: Dict[str, Any] = dict(n_voxels=n)

    # ---- the premise
    S["premise"] = premise(L)

    # ---- entries
    ents = paint["entries"]
    Re = np.array([e["R"] for e in ents])
    seed_v = np.array([e["voxel"] for e in ents])
    e_err = np.array([mis(Re[i], Rtrue[seed_v[i]], sym) for i in range(len(ents))])
    e_grain = [L.grain_of(Re[i]) for i in range(len(ents))]
    per_grain: Dict[int, int] = {}
    for g in e_grain:
        per_grain[-1 if g is None else g] = per_grain.get(-1 if g is None else g, 0) + 1
    near = [e for e in ents if e["near_earlier"] is not None]
    S["entries"] = dict(
        n=len(ents),
        kinds={k: sum(e["kind"] == k for e in ents) for k in ("index", "grid")},
        wrong=int((e_err > WRONG).sum()),
        true_grains_with_entry=len([g for g in per_grain if g != -1]),
        true_grains_in_region=int(len(np.unique(gid))),
        grains_with_two_or_more=int(sum(v >= 2 for g, v in per_grain.items() if g != -1)),
        near_pairs=len(near),
        near_deg=[float(e["near_deg"]) for e in near],
        median_seed_err_deg=float(np.median(e_err)),
        flood_scores_mean=float(np.mean([e["flood_scores"] for e in ents])),
        flood_passes_mean=float(np.mean([e["flood_passes"] for e in ents])),
    )
    # grains of the region without an entry
    have = {g for g in e_grain if g is not None}
    S["entries"]["true_grains_missing_entry"] = sorted(int(g) for g in set(gid) - have)
    S["entries"]["voxels_of_missing_grains"] = int(np.isin(gid, list(set(gid) - have)).sum())

    # ---- scores per voxel
    passing: Dict[int, List[Tuple[int, float, float]]] = {}
    scored: Dict[int, Dict[int, Tuple[float, float, bool]]] = {}
    for g, rows in enumerate(paint["scores"]):
        for v, c, h, ok in rows:
            scored.setdefault(v, {})[g] = (c, h, ok)
            if ok:
                passing.setdefault(v, []).append((g, c, h))
    n_pass = np.array([len(passing.get(v, [])) for v in range(n)])
    stage = {int(k): v for k, v in paint["stage"].items()}
    S["voxels"] = dict(
        by_source={s: int((E["source"] == s).sum()) for s in sorted(set(E["source"]))},
        unpainted=int((n_pass == 0).sum()),
        one_grain=int((n_pass == 1).sum()),
        contested=int((n_pass >= 2).sum()),
        seeds=int(sum(s == "seed" for s in stage.values())),
        fits=int(sum(s == "fit" for s in stage.values())),
        failed=int(sum(s == "failed" for s in stage.values())),
    )

    # ---- evaluations
    ev = dict(run["evals"])
    # The stored paint run predates the split of "seed_rejected" into "index_rejected" and
    # "grid_rejected" (review fix 10). The old bucket mixed the two, so it is only read when 0.
    if "seed_rejected" in ev:
        legacy = int(ev.pop("seed_rejected"))
        assert legacy == 0, f"unsplit seed_rejected = {legacy}: rerun the paint run"
        ev.setdefault("index_rejected", 0)
        ev.setdefault("grid_rejected", 0)
    key_on = json.loads((C.STORED_KEY_ON / "run.json").read_text())["stats"]
    key_f = json.loads((C.STORED_KEY_ON_FUSED / "run.json").read_text())["stats"]
    key_off = KS.baseline()
    ke = KS.evals_of(key_on)
    kf = KS.evals_of(key_f)
    assert ke["total"] == STORED_KEY_ON_EVALS, ke["total"]
    paint_total = int(sum(ev.values()))
    fit_ev = ev["fit_accepted"] + ev["fit_rejected"]

    def on_split(st: Dict[str, Any], e: Dict[str, int]) -> Dict[str, int]:
        return dict(
            total=e["total"],
            seed_paths=int(st["index_evals"]),
            scoring=0,
            local_fits=e["local"],
            grid=int(st.get("index_fallback_grid_evals", 0)),
        )

    S["evals"] = dict(
        key_off=dict(
            total=key_off["evals"]["total"],
            seed_paths=0,
            scoring=0,
            local_fits=key_off["evals"]["local"],
            grid=key_off["evals"]["full"],
            full_searches=key_off["full_searches"],  # grid searches of the key-off run (fix 11)
        ),
        key_on=on_split(key_on, ke),
        key_on_fused=on_split(key_f, kf),
        paint=dict(
            total=paint_total,
            seed_paths=ev["seed_index"] + ev["index_rejected"],
            scoring=ev["score"],
            local_fits=fit_ev,
            grid=ev["seed_grid"] + ev["grid_rejected"],
            index_rejected=ev["index_rejected"],
            grid_rejected=ev["grid_rejected"],
            fit_accepted=ev["fit_accepted"],
            fit_rejected=ev["fit_rejected"],
        ),
        side_arm=dict(extra=int(json.loads((RUN / "side" / "side.json").read_text())["evals"])),
        paint_share_of_key_on=paint_total / ke["total"],
        paint_share_of_key_on_fused=paint_total / kf["total"],
        paint_share_of_key_off=paint_total / key_off["evals"]["total"],
        cost_fn_eval_count=run["cost_fn_eval_count"],
    )

    # ---- wrong voxels, paired
    S["wrong"] = dict(
        paint=wil_dict(int(wrong.sum()), n),
        key_on=wil_dict(int((EK["err"] > WRONG).sum()), n),
        key_on_fused=wil_dict(int((EF["err"] > WRONG).sum()), n),
        key_off=wil_dict(int((EO["err"] > WRONG).sum()), n),
        side=wil_dict(int((ES["err"] > WRONG).sum()), n),
        paint_vs_key_on=paired(EK["err"], E["err"]),
        paint_vs_key_on_fused=paired(EF["err"], E["err"]),
        paint_vs_key_off=paired(EO["err"], E["err"]),
        side_vs_paint=paired(E["err"], ES["err"]),
        side_vs_key_on=paired(EK["err"], ES["err"]),
        key_on_vs_key_on_fused=paired(EK["err"], EF["err"]),
        kinds_wrong={
            "paint": {k: int(((kind == k) & wrong).sum()) for k in ("edge", "corner", "interior")},
            "key_on": {
                k: int(((kind == k) & (EK["err"] > WRONG)).sum())
                for k in ("edge", "corner", "interior")
            },
        },
        kinds_all={k: int((kind == k).sum()) for k in ("edge", "corner", "interior")},
    )

    # ---- trace of every wrong voxel of the paint run (and the breaks against key on)
    def trace(r: int) -> Dict[str, Any]:
        f = int(full[r])
        R = None
        src = str(E["source"][r])
        pg = None
        if src == "grain":
            ranked = sorted(passing.get(r, []), key=lambda t: (t[1], t[0]))
            pg = ranked[0][0]
            R = Re[pg]
        else:
            R = None
        g_taken = L.grain_of(Re[pg]) if pg is not None else None
        # grains whose orientation is the voxel's truth: entries within 1 degree of it
        true_entries = [i for i in range(len(ents)) if mis(Re[i], Rtrue[r], sym) <= WRONG]
        passed_true = [
            i for i in true_entries if i in dict((g, 1) for g, _, _ in passing.get(r, []))
        ]
        if not true_entries:
            why = "no entry of the true grain"
        elif not passed_true:
            why = (
                "true grain's entry scored but did not pass"
                if any(i in scored.get(r, {}) for i in true_entries)
                else "true grain's entry never scored the voxel"
            )
        else:
            why = "true grain's entry passed and lost on cost"
        d: Dict[str, Any] = dict(
            voxel=int(r),
            full_sample_voxel=f,
            kind=str(kind[r]),
            err=float(E["err"][r]),
            source=src,
            painted_entry=pg,
            painted_entry_seed=None if pg is None else int(seed_v[pg]),
            taken_grain=g_taken,
            true_grain=int(gid[r]),
            crossing=bool(g_taken is not None and g_taken != int(gid[r]) and L.touches(f, g_taken)),
            n_passing=int(n_pass[r]),
            why=why,
            cost_painted=None if pg is None else float(scored[r][pg][0]),
            cost_true_entry=(
                None if not passed_true else float(min(scored[r][i][0] for i in passed_true))
            ),
            err_key_on=float(EK["err"][r]),
        )
        return d

    wrong_idx = [int(r) for r in np.nonzero(wrong)[0]]
    tr = [trace(r) for r in wrong_idx]
    brk = [t for t in tr if t["err_key_on"] <= WRONG]
    S["wrong_trace"] = tr
    S["broken_trace"] = brk
    S["wrong_summary"] = dict(
        n=len(tr),
        why={w: sum(t["why"] == w for t in tr) for w in sorted({t["why"] for t in tr})},
        crossing={k: sum(t["crossing"] and t["kind"] == k for t in tr) for k in ("edge", "corner")},
        crossing_total=sum(t["crossing"] for t in tr),
        interior_wrong=sum(t["kind"] == "interior" for t in tr),
        broken_n=len(brk),
        broken_why={w: sum(t["why"] == w for t in brk) for w in sorted({t["why"] for t in brk})},
        broken_crossing={
            k: sum(t["crossing"] and t["kind"] == k for t in brk) for k in ("edge", "corner")
        },
        broken_kind={k: sum(t["kind"] == k for t in brk) for k in ("edge", "corner", "interior")},
    )

    # ---- direct test on the contested voxels
    cont = [int(v) for v in np.nonzero(n_pass >= 2)[0]]
    rows: List[Dict[str, Any]] = []
    for v in cont:
        pg = passing[v]
        by_cost = sorted(pg, key=lambda t: (t[1], t[0]))
        by_hit = sorted(pg, key=lambda t: (-t[2], t[1], t[0]))
        c0, h0 = by_cost[0][0], by_hit[0][0]
        pair_deg = [mis(Re[a[0]], Re[b[0]], sym) for a in pg for b in pg if a[0] < b[0]]
        distinct = any(d > WRONG for d in pair_deg)
        rows.append(
            dict(
                voxel=v,
                kind=str(kind[v]),
                n_passing=len(pg),
                distinct=bool(distinct),
                max_pair_deg=float(max(pair_deg)),
                cost_right=bool(mis(Re[c0], Rtrue[v], sym) <= WRONG),
                hit_right=bool(mis(Re[h0], Rtrue[v], sym) <= WRONG),
                cost_hit_differ=bool(mis(Re[c0], Re[h0], sym) > WRONG),
                true_among=bool(any(mis(Re[g], Rtrue[v], sym) <= WRONG for g, _, _ in pg)),
                margin=float(by_cost[1][1] - by_cost[0][1]),
            )
        )
    direct: Dict[str, Any] = dict(n_contested=len(rows))
    # Contested voxels whose passing entries are all within 1 degree of each other: duplicate
    # entries of one grain. Every rule picks an entry of that grain there, so the direct test cannot
    # fail on them and cannot tell the rules apart (review fix 2).
    dup = [r for r in rows if not r["distinct"]]
    direct["duplicate_entries"] = dict(
        n=len(dup),
        kinds={k: sum(r["kind"] == k for r in dup) for k in ("edge", "corner", "interior")},
        true_among=sum(r["true_among"] for r in dup),
        cost_right=sum(r["cost_right"] for r in dup),
        hit_right=sum(r["hit_right"] for r in dup),
        cost_hit_differ=sum(r["cost_hit_differ"] for r in dup),  # 0 by construction
        max_pair_deg=float(max((r["max_pair_deg"] for r in dup), default=float("nan"))),
        pair_deg_values=sorted({round(r["max_pair_deg"], 5) for r in dup}),
    )
    for sel_name, sel in (
        ("all", lambda r: True),
        ("different_orientations", lambda r: r["distinct"]),
    ):
        sub = [r for r in rows if sel(r)]
        d: Dict[str, Any] = dict(
            n=len(sub),
            true_among=wil_dict(sum(r["true_among"] for r in sub), len(sub)),
            cost_right=wil_dict(sum(r["cost_right"] for r in sub), len(sub)),
            hit_right=wil_dict(sum(r["hit_right"] for r in sub), len(sub)),
            cost_hit_differ=wil_dict(sum(r["cost_hit_differ"] for r in sub), len(sub)),
        )
        for k in ("edge", "corner", "interior"):
            ss = [r for r in sub if r["kind"] == k]
            d[k] = dict(
                n=len(ss),
                cost_right=wil_dict(sum(r["cost_right"] for r in ss), len(ss)),
                hit_right=wil_dict(sum(r["hit_right"] for r in ss), len(ss)),
                true_among=wil_dict(sum(r["true_among"] for r in ss), len(ss)),
                cost_hit_differ=wil_dict(sum(r["cost_hit_differ"] for r in ss), len(ss)),
            )
        direct[sel_name] = d
    direct["rows"] = rows
    S["direct"] = direct
    kw = np.nonzero(EK["err"] > WRONG)[0]
    # Does the key-on run's (wrong) orientation pass in paint? An entry within 1 degree of it that
    # passes on the voxel (review fix 6: where it does not, the fix comes from not fitting).
    from score_bfs import mats_to_quats, pair_misorientation_deg

    q_ent = mats_to_quats(Re)

    def near_key_on(v: int, gs: List[int]) -> bool:
        if not gs:
            return False
        qk = np.broadcast_to(EK["q_rec"][v], (len(gs), 4))
        return bool((pair_misorientation_deg(q_ent[gs], qk, EK["sym"]) <= WRONG).any())

    wrong_passes = [near_key_on(int(v), [g for g, _, _ in passing.get(int(v), [])]) for v in kw]
    wrong_scored = [
        near_key_on(int(v), [g for g, t in scored.get(int(v), {}).items() if not t[2]]) for v in kw
    ]
    S["key_on_wrong_in_paint"] = dict(
        n=len(kw),
        passing_counts={str(k): int((n_pass[kw] == k).sum()) for k in sorted(set(n_pass[kw]))},
        painted_right=int((E["err"][kw] <= WRONG).sum()),
        wrong_orientation_passes=int(sum(wrong_passes)),
        wrong_orientation_does_not_pass=int(len(kw) - sum(wrong_passes)),
        wrong_orientation_scored_not_passed=int(
            sum(s and not p for s, p in zip(wrong_scored, wrong_passes))
        ),
        one_passing_and_wrong_does_not_pass=int(
            sum((n_pass[v] == 1) and not p for v, p in zip(kw, wrong_passes))
        ),
    )

    # ---- how far a wrong grain floods
    trees: Dict[int, Any] = {}

    def tree(g: int) -> Any:
        if g not in trees:
            trees[g] = cKDTree(L.centre[L.gmap == g])
        return trees[g]

    bins = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 1e9)]
    flood_hist = {f"{a:g} to {b:g}" if b < 1e8 else f"{a:g} or more": 0 for a, b in bins}
    flood_out = 0
    flood_in = 0
    flood_no_grain = 0
    per_entry_out: List[int] = []
    for i, rows_i in enumerate(paint["scores"]):
        tg = e_grain[i]
        n_out = 0
        for v, c, h, ok in rows_i:
            if not ok or v == ents[i]["voxel"]:
                continue
            if mis(Re[i], Rtrue[v], sym) <= WRONG:
                flood_in += 1
                continue
            flood_out += 1
            n_out += 1
            if tg is None:
                flood_no_grain += 1
                continue
            dist = float(tree(tg).query(L.centre[int(full[v])])[0] / L.side[int(full[v])])
            for (a, b), key in zip(bins, flood_hist):
                if a <= dist < b:
                    flood_hist[key] += 1
                    break
        per_entry_out.append(n_out)
    S["flood"] = dict(
        passes_inside_true_grain=flood_in,
        passes_outside=flood_out,
        outside_no_true_grain_for_entry=flood_no_grain,
        outside_by_distance_sides=flood_hist,
        entries_with_outside_passes=int(sum(x > 0 for x in per_entry_out)),
        max_outside_per_entry=int(max(per_entry_out)) if per_entry_out else 0,
        outside_voxels_by_kind={
            k: sum(
                1
                for i, rows_i in enumerate(paint["scores"])
                for v, c, h, ok in rows_i
                if ok
                and v != ents[i]["voxel"]
                and kind[v] == k
                and mis(Re[i], Rtrue[v], sym) > WRONG
            )
            for k in ("edge", "corner", "interior")
        },
    )

    # ---- accuracy without a fit
    by_grain = E["source"] == "grain"
    right_both = by_grain & (E["err"] <= WRONG) & (EK["err"] <= WRONG)
    q = [50, 90, 99]
    # A painted voxel carries its entry's orientation and the truth is one orientation per grain,
    # so the paint errors take one value per painting entry (review fix 5).
    rb = [int(v) for v in np.nonzero(right_both)[0]]
    painted_entry = {v: min(passing[v], key=lambda t: (t[1], t[0]))[0] for v in rb}

    def n_orient(e: Any, sel: np.ndarray) -> int:  # distinct stored orientations (to 1e-7)
        q = e["q_rec"][sel]
        q = q * np.where(q[:, :1] < 0, -1.0, 1.0)
        return int(len(np.unique(np.round(q, 7), axis=0)))

    S["accuracy_distinct"] = dict(
        n_entries_painting=len(set(painted_entry.values())),
        n_distinct_paint_err=len({round(float(E["err"][v]), 9) for v in rb}),
        n_distinct_key_on_err=len({round(float(EK["err"][v]), 9) for v in rb}),
        # the key-on map is not one fit result per voxel either: it holds few orientations
        n_orientations_paint_map=n_orient(E, np.arange(n)),
        n_orientations_key_on_map=n_orient(EK, np.arange(n)),
        n_orientations_key_on_fused_map=n_orient(EF, np.arange(n)),
        n_orientations_key_on_right_both=n_orient(EK, right_both),
    )
    S["accuracy"] = dict(
        n_painted_by_grain=int(by_grain.sum()),
        n_right_in_both=int(right_both.sum()),
        paint_pct=dict(zip(map(str, q), map(float, np.percentile(E["err"][right_both], q)))),
        key_on_pct=dict(zip(map(str, q), map(float, np.percentile(EK["err"][right_both], q)))),
        paint_max=float(E["err"][right_both].max()),
        key_on_max=float(EK["err"][right_both].max()),
        over_0p1_paint=int((E["err"][right_both] > 0.1).sum()),
        over_0p1_key_on=int((EK["err"][right_both] > 0.1).sum()),
        over_0p05_paint=int((E["err"][right_both] > 0.05).sum()),
        over_0p05_key_on=int((EK["err"][right_both] > 0.05).sum()),
        entry_seed_err_pct=dict(zip(map(str, q), map(float, np.percentile(e_err, q)))),
    )

    # ---- the side arm
    side = json.loads((RUN / "side" / "side.json").read_text())
    srows = side["rows"]
    curve: List[Dict[str, Any]] = []
    for m in MARGINS + [float("inf")]:
        sel = [r for r in srows if r["margin"] < m]
        fixed = sum(1 for r in sel if wrong[r["voxel"]] and not (ES["err"][r["voxel"]] > WRONG))
        broken = sum(1 for r in sel if (not wrong[r["voxel"]]) and ES["err"][r["voxel"]] > WRONG)
        curve.append(
            dict(
                margin=m if np.isfinite(m) else None,
                voxels=len(sel),
                fits=2 * len(sel),
                evals=int(sum(sum(r["fit_evals"]) for r in sel)),
                fixed=fixed,
                broken=broken,
                wrong_after=int(wrong.sum()) - fixed + broken,
            )
        )
    S["side"] = dict(
        n_contested_fitted=len(srows),
        evals=side["evals"],
        wall_s=side["wall_s"],
        curve=curve,
        kept_lower_cost_first=sum(r["kept"] == 0 for r in srows),
        both_fits_wrong_or_right=dict(
            kept_right_when_paint_wrong=int(
                sum(wrong[r["voxel"]] and ES["err"][r["voxel"]] <= WRONG for r in srows)
            ),
        ),
        total_evals_with_side=paint_total + side["evals"],
        side_share_of_key_on=(paint_total + side["evals"]) / ke["total"],
    )

    # ---- wall time and profile
    prof = json.loads((OUT / "runs" / "keyon_profile" / "run.json").read_text())
    tim_p = OUT / "runs" / "keyon_timers" / "run.json"
    tim = json.loads(tim_p.read_text()) if tim_p.exists() else None
    S["wall"] = dict(
        paint=dict(
            wall_total_s=run["wall_total_s"],
            pairing_s=run["index_pairing_s"],
            by_stage_s=run["wall_by_stage_s"],
            label=run["timing_label"],
            fused_cost=run["fused_cost"],
        ),
        side_wall_s=side["wall_s"],
        key_on_fused_stored=dict(
            wall_total_s=json.loads((C.STORED_KEY_ON_FUSED / "run.json").read_text())[
                "wall_total_s"
            ],
            label=json.loads((C.STORED_KEY_ON_FUSED / "run.json").read_text())["timing_label"],
        ),
        profile=dict(
            wall_total_s=prof["wall_total_s"],
            profiled=prof["profiled"],
            label=prof["timing_label"],
            evals_total=prof["evals_total"],
            parts=prof["parts"],
            pairing_s=prof["pairing_s"],
            top=prof["top_by_tottime"][:12],
        ),
        timers_only=(
            None
            if tim is None
            else dict(
                wall_total_s=tim["wall_total_s"],
                label=tim["timing_label"],
                evals_total=tim["evals_total"],
                parts=tim["parts"],
                pairing_s=tim["pairing_s"],
            )
        ),
    )

    # ---- gate
    brk_n = S["wrong"]["paint_vs_key_on"]["broken"]
    S["gate"] = dict(
        broken_vs_key_on=brk_n,
        broken_limit=GATE_BREAKS,
        evals_share_of_key_on=paint_total / ke["total"],
        evals_limit_share=0.5,
        go=bool(brk_n <= GATE_BREAKS and paint_total <= 0.5 * ke["total"]),
    )
    (OUT / "summary.json").write_text(json.dumps(S, indent=1, default=float) + "\n")
    tables(S)
    figure(S, E, EK, EF, EO, L)
    print(
        json.dumps(
            {k: S[k] for k in ("evals", "wrong", "gate", "wrong_summary")}, indent=1, default=float
        )[:6000]
    )


def tables(S: Dict[str, Any]) -> None:
    T: Dict[str, str] = {}
    ev = S["evals"]
    w = S["wrong"]

    def row(name: str, e: Dict[str, Any], wk: Dict[str, Any]) -> Dict[str, Any]:
        return dict(
            arm=name,
            total=e["total"],
            seed=e["seed_paths"],
            score=e["scoring"],
            fit=e["local_fits"],
            grid=e["grid"],
            wrong=f"{wk['k']} of {wk['n']} ({100 * wk['lo']:.1f} to {100 * wk['hi']:.1f}%)",
        )

    rows = [
        row("key off (stored rerun)", ev["key_off"], w["key_off"]),
        row("key on (stored)", ev["key_on"], w["key_on"]),
        row("key on, FusedCost (stored)", ev["key_on_fused"], w["key_on_fused"]),
        row("paint (FusedCost)", ev["paint"], w["paint"]),
    ]
    sd = dict(
        total=ev["paint"]["total"] + ev["side_arm"]["extra"],
        seed_paths=ev["paint"]["seed_paths"],
        scoring=ev["paint"]["scoring"],
        local_fits=ev["paint"]["local_fits"] + ev["side_arm"]["extra"],
        grid=ev["paint"]["grid"],
    )
    rows.append(row("paint + fit the contested (side arm)", sd, w["side"]))
    T["ar3_evals"] = (
        DT.markdown_table(
            rows,
            ["arm", "total", "seed", "score", "fit", "grid", "wrong"],
            formats={k: "," for k in ("total", "seed", "score", "fit", "grid")},
            labels=None,
        )
        .replace("| seed |", "| seed paths |")
        .replace("| score |", "| scoring |")
    )
    T["ar3_evals"] = (
        T["ar3_evals"].replace("| fit |", "| local fits |").replace("| grid |", "| grid searches |")
    )

    def prow(name: str, p: Dict[str, int]) -> Dict[str, Any]:
        return dict(vs=name, **p)

    pr = [
        prow("key on (stored)", w["paint_vs_key_on"]),
        prow("key on, FusedCost (stored)", w["paint_vs_key_on_fused"]),
        prow("key off (stored)", w["paint_vs_key_off"]),
    ]
    T["ar3_paired"] = (
        DT.markdown_table(
            pr,
            ["vs", "wrong_a", "wrong_b", "fixed", "broken", "both_wrong"],
            labels=None,
        )
        .replace("wrong_a", "wrong there")
        .replace("wrong_b", "wrong in paint")
        .replace("both_wrong", "wrong in both")
    )
    ws = S["wrong_summary"]
    wt = [
        dict(what="wrong voxels of paint", n=ws["n"]),
        dict(what="  of which crossing a boundary, sharing an edge", n=ws["crossing"]["edge"]),
        dict(what="  of which crossing a boundary, corner only", n=ws["crossing"]["corner"]),
        dict(what="  of which interior voxels", n=ws["interior_wrong"]),
    ] + [dict(what=f"  why: {k}", n=v) for k, v in ws["why"].items()]
    T["ar3_wrong"] = DT.markdown_table(wt, ["what", "n"])
    d = S["direct"]

    def drow(label: str, x: Dict[str, Any]) -> Dict[str, Any]:
        f = lambda t: f"{t['k']} of {t['n']} ({100 * t['lo']:.0f} to {100 * t['hi']:.0f}%)"  # noqa
        return dict(
            set=label,
            n=x["n"] if "n" in x else x["true_among"]["n"],
            true_among=f(x["true_among"]),
            cost_right=f(x["cost_right"]),
            hit_right=f(x["hit_right"]),
            cost_hit_differ=f(x["cost_hit_differ"]),
        )

    # Review fixes 1 and 2: the informative set leads (passing grains over 1 degree apart, with
    # intervals); the duplicate-entry voxels, where the test cannot fail, get bare counts only.
    far = "passing grains over 1 degree apart"
    dr = [drow(f"{far}: all kinds", d["different_orientations"])]
    for k in ("edge", "corner", "interior"):
        if d["different_orientations"][k]["n"]:
            dr.append(drow(f"{far}: {k}", d["different_orientations"][k]))
    du = d["duplicate_entries"]
    bare = "{} of {} (any rule; no interval)"
    dr.append(
        dict(
            set=(
                "passing entries all of one grain (duplicates within 1 degree; "
                "rules cannot differ)"
            ),
            n=du["n"],
            true_among=f"{du['true_among']} of {du['n']} (no interval)",
            cost_right=bare.format(du["cost_right"], du["n"]),
            hit_right=bare.format(du["hit_right"], du["n"]),
            cost_hit_differ=f"{du['cost_hit_differ']} of {du['n']} (by construction)",
        )
    )
    a = d["all"]
    dr.append(
        dict(
            set="all contested (the two sets above)",
            n=a["n"],
            true_among=f"{a['true_among']['k']} of {a['n']}",
            cost_right=f"{a['cost_right']['k']} of {a['n']}",
            hit_right=f"{a['hit_right']['k']} of {a['n']}",
            cost_hit_differ=f"{a['cost_hit_differ']['k']} of {a['n']}",
        )
    )
    cols = ["set", "n", "true_among", "cost_right", "hit_right", "cost_hit_differ"]
    head = {
        "true_among": "true grain passes",
        "cost_right": "lowest cost is the true grain",
        "hit_right": "highest hit is the true grain",
        "cost_hit_differ": "lowest cost and highest hit name grains over 1 degree apart",
    }
    md = DT.markdown_table(dr, cols).split("\n")
    md[0] = "| " + " | ".join(head.get(c, c) for c in cols) + " |"  # header row only
    T["ar3_direct"] = "\n".join(md)
    cv = [
        dict(
            m="all" if c["margin"] is None else f"{c['margin']:g}",
            voxels=c["voxels"],
            fits=c["fits"],
            evals=c["evals"],
            fixed=c["fixed"],
            broken=c["broken"],
            wrong=c["wrong_after"],
        )
        for c in S["side"]["curve"]
    ]
    T["ar3_margin"] = DT.markdown_table(
        cv,
        ["m", "voxels", "fits", "evals", "fixed", "broken", "wrong"],
        formats={"evals": ","},
    ).replace("| m |", "| fit if the two lowest costs differ by less than |")
    fl = S["flood"]
    ft = [dict(dist=k, n=v) for k, v in fl["outside_by_distance_sides"].items()]
    T["ar3_flood"] = DT.markdown_table(ft, ["dist", "n"]).replace(
        "| dist |", "| distance from the grain's true region (voxel sides) |"
    )
    ac = S["accuracy"]
    at = [
        dict(
            arm="paint (no fit)",
            p50=ac["paint_pct"]["50"],
            p90=ac["paint_pct"]["90"],
            p99=ac["paint_pct"]["99"],
            mx=ac["paint_max"],
        ),
        dict(
            arm="key on (fitted)",
            p50=ac["key_on_pct"]["50"],
            p90=ac["key_on_pct"]["90"],
            p99=ac["key_on_pct"]["99"],
            mx=ac["key_on_max"],
        ),
    ]
    T["ar3_accuracy"] = (
        DT.markdown_table(
            at,
            ["arm", "p50", "p90", "p99", "mx"],
            formats={k: ".4f" for k in ("p50", "p90", "p99", "mx")},
        )
        .replace("| p50 |", "| median |")
        .replace("| p90 |", "| 90th |")
        .replace("| p99 |", "| 99th |")
        .replace("| mx |", "| largest |")
    )
    wl = S["wall"]
    prof = wl["profile"]
    wt2 = [
        dict(what="paint run (FusedCost, unbatched loop)", s=wl["paint"]["wall_total_s"]),
        dict(what="  of which the one-time pairing", s=wl["paint"]["pairing_s"]),
        dict(
            what="key on, FusedCost, timers only, this session",
            s=None if wl["timers_only"] is None else wl["timers_only"]["wall_total_s"],
        ),
        dict(
            what="key on, FusedCost, stored (AR2 step 2)",
            s=wl["key_on_fused_stored"]["wall_total_s"],
        ),
        dict(what="key on under cProfile (not a timing)", s=prof["wall_total_s"]),
        dict(what="side arm (fits on the contested voxels)", s=wl["side_wall_s"]),
    ]
    T["ar3_wall"] = DT.markdown_table(wt2, ["what", "s"], formats={"s": ".1f"})
    parts = (wl["timers_only"] or prof)["parts"]
    tot = (wl["timers_only"] or prof)["wall_total_s"]
    pt = []
    lf = parts["local_fit"]
    sp = parts["seed_path"]
    pt.append(dict(part="local fits, cost evaluations", s=lf["eval_s"], n=lf["evals"]))
    pt.append(
        dict(part="local fits, optimizer's own work", s=lf["incl_s"] - lf["eval_s"], n=lf["calls"])
    )
    pt.append(
        dict(
            part="seed path: indexing", s=parts["indexing"]["incl_s"], n=parts["indexing"]["calls"]
        )
    )
    pt.append(
        dict(
            part="seed path: finisher, cost evaluations",
            s=parts["finisher"]["eval_s"],
            n=parts["finisher"]["evals"],
        )
    )
    pt.append(
        dict(
            part="seed path: finisher, optimizer's own work",
            s=parts["finisher"]["incl_s"] - parts["finisher"]["eval_s"],
            n=parts["finisher"]["calls"],
        )
    )
    pair_s = (wl["timers_only"] or prof)["pairing_s"]
    pt.append(
        dict(
            part="seed path: strict costs on candidates and the rest of the seed path",
            s=sp["incl_s"] - parts["indexing"]["incl_s"] - parts["finisher"]["incl_s"] - pair_s,
            n=sp["calls"],
        )
    )
    pt.append(
        dict(
            part="overlap evaluations (the BFS's own)",
            s=parts["overlap"]["incl_s"],
            n=parts["overlap"]["calls"],
        )
    )
    pt.append(
        dict(part="grid search", s=parts["grid_search"]["incl_s"], n=parts["grid_search"]["calls"])
    )
    pt.append(
        dict(
            part="one-time pairing (inside the first seed path)",
            s=(wl["timers_only"] or prof)["pairing_s"],
            n=None,
        )
    )
    acc = sum(p["s"] for p in pt)
    pt.append(
        dict(part="everything else (BFS bookkeeping, neighbours, records)", s=tot - acc, n=None)
    )
    for p in pt:
        p["share"] = p["s"] / tot
    T["ar3_profile"] = DT.markdown_table(
        pt, ["part", "s", "n", "share"], formats={"s": ".1f", "share": ".1%", "n": ","}
    ).replace("| n |", "| calls or evaluations |")
    pm = S["premise"]
    T["ar3_premise"] = DT.markdown_table(
        [dict(scope=k, **v) for k, v in pm.items()],
        [
            "scope",
            "n_grains",
            "n_voxels",
            "max_spread_deg",
            "n_grains_over_1e3_deg",
            "n_grains_over_1_deg",
        ],
        formats={"max_spread_deg": ".5f"},
    )
    DT.write_tables(OUT / "tables.md", T)


def figure(S: Dict[str, Any], E: Any, EK: Any, EF: Any, EO: Any, L: Any) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    # left: the effect (wrong voxels by arm, Wilson 95% interval, evaluations under each bar)
    arms = [
        ("key off", "key_off", "C7"),
        ("key on", "key_on", "C0"),
        ("paint", "paint", "C1"),
    ]
    ev = S["evals"]
    for i, (lab, k, col) in enumerate(arms):
        w = S["wrong"][k]
        ax[0].bar(i, w["k"], color=col)
        lo, hi = w["lo"] * w["n"], w["hi"] * w["n"]
        ax[0].errorbar(i, w["k"], yerr=[[w["k"] - lo], [hi - w["k"]]], color="k", capsize=4)
        ax[0].annotate(
            str(w["k"]),
            (i, hi),
            ha="center",
            va="bottom",
            fontsize=9,
            xytext=(0, 2),
            textcoords="offset points",
        )
    ax[0].set_xticks(range(len(arms)))
    ax[0].set_xticklabels([f"{lab}\n{ev[k]['total']:,}\nevaluations" for lab, k, _ in arms])
    ax[0].set_ylim(0, 1.15 * max(S["wrong"][k]["hi"] * S["wrong"][k]["n"] for _, k, _ in arms))
    ax[0].set_ylabel("wrong voxels of 2,000 (Wilson 95% interval)")
    ax[0].set_title("Wrong voxels by arm (one run each)")
    by_grain = E["source"] == "grain"
    ok = by_grain & (E["err"] <= WRONG) & (EK["err"] <= WRONG)
    for e, lab, col in ((E, "paint (no fit)", "C1"), (EK, "key on (every voxel fitted)", "C0")):
        x = np.sort(np.maximum(e["err"][ok], 1e-4))
        ax[1].plot(x, np.arange(1, len(x) + 1) / len(x), label=lab, color=col)
    ax[1].set_xscale("log")
    ax[1].set_xlabel("distance to the truth (degrees)")
    ax[1].set_ylabel("share of voxels")
    nd = S["accuracy_distinct"]["n_entries_painting"]
    ax[1].set_title(
        f"Voxels painted and right in both runs (n={int(ok.sum())});\n"
        f"paint takes {nd} distinct values, one per entry",
        fontsize=10,
    )
    ax[1].legend(fontsize=8)
    fh = S["flood"]["outside_by_distance_sides"]
    ax[2].bar(list(fh), list(fh.values()), color="C3")
    ax[2].set_xlabel("distance from the grain's true region (voxel sides)")
    ax[2].set_ylabel("voxels where a grain passes outside its true grain")
    ax[2].set_title("How far a wrong grain floods")
    fig.tight_layout()
    fig.savefig(OUT / "ar3_paint.png", dpi=130)


if __name__ == "__main__":
    main()
