"""AR3 step 2: summary.json, tables.md and the figure of runs A-D (and the side arm).

  uv run python scripts/autoresearch/ar3/layer_summary.py

Reads ``benchmarks/autoresearch/ar3/layer/runs/{A,B,C,D}`` (+ ``C/side``). Nothing is rerun. Wrong =
more than 1 degree from the truth (symmetry-reduced). Counts by voxel; intervals are Wilson; no
p-values. Every count of C and D is given for the whole layer, the pilot region (2,000 voxels) and
the rest (22,570 voxels). A voxel is the layer's grid index.
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
import key_summary as KS  # noqa: E402
from run_layer import LAYER_OUT, LAYER_SAMPLE  # noqa: E402

from icenine import ray_indexing as RI  # noqa: E402
from icenine.mic_file import MicFile  # noqa: E402

import doc_tables as DT  # noqa: E402
import stats as ST  # noqa: E402

RUNS = LAYER_OUT / "runs"
WRONG = 1.0
KINDS = ("edge", "corner", "interior")
PILOT_MIC = C.REGION / "region_grid.mic"


def wil_dict(k: int, n: int) -> Dict[str, Any]:
    lo, hi = ST.wilson(k, n) if n else (float("nan"), float("nan"))
    return dict(k=int(k), n=int(n), lo=float(lo), hi=float(hi))


def wil_str(d: Dict[str, Any]) -> str:
    if not d["n"]:
        return "0 of 0"
    return f"{d['k']} of {d['n']} ({100 * d['lo']:.1f} to {100 * d['hi']:.1f}%)"


def overlap(a: Dict[str, Any], b: Dict[str, Any]) -> bool:
    return not (a["hi"] < b["lo"] or b["hi"] < a["lo"])


class World:
    """The layer's truth and the maps of the runs, all over the layer's voxel index."""

    def __init__(self) -> None:
        self.L = KS.Layer()
        L = self.L
        self.sym = L.sym
        self.n = len(L.gmap)
        pos = np.array([L._pos[int(g)] for g in L.gmap])
        self.Rtrue = L.Rg[pos]
        truth = MicFile.read(str(C.EX / "SimInput" / "rand_500grains_1mm_neworient_s0.mic"))
        self.xy = np.array([v.position[:2] for v in truth.voxels], dtype=float)
        assert len(self.xy) == self.n
        self.region = np.asarray(L.region)  # region row -> layer voxel (unsorted, as stored)
        self.in_region = np.zeros(self.n, dtype=bool)
        self.in_region[self.region] = True
        self.kind = np.array([L.kind(v) for v in range(self.n)])

    def read_map(self, run: Path) -> Tuple[np.ndarray, np.ndarray]:
        """(orientations (n,3,3), have (n,) bool) of a run's recon.mic over the layer's index."""
        mic = MicFile.read(str(run / "recon.mic"))
        R = np.array([v.orientation for v in mic.voxels], dtype=np.float64)
        xy = np.array([v.position[:2] for v in mic.voxels], dtype=float)
        have = np.zeros(self.n, dtype=bool)
        full = np.zeros((self.n, 3, 3))
        if len(R) == self.n:
            assert np.abs(xy - self.xy).max() < 1e-4
            return R, np.ones(self.n, dtype=bool)
        assert len(R) == len(self.region)
        assert np.abs(xy - self.xy[self.region]).max() < 1e-4
        full[self.region] = R
        have[self.region] = True
        return full, have

    def err(self, R: np.ndarray) -> np.ndarray:
        return np.asarray(RI.misorientation_deg(R, self.Rtrue, self.sym), dtype=float)

    def mis(self, Ra: np.ndarray, Rb: np.ndarray) -> float:
        return float(RI.misorientation_deg(Ra[None], Rb[None], self.sym)[0])

    def crossing(self, v: int, R: np.ndarray) -> Optional[bool]:
        """Does the voxel carry the orientation of a different grain that touches its voxel?"""
        g = self.L.grain_of(R)
        if g is None:
            return None
        return bool(g != int(self.L.gmap[v]) and self.L.touches(v, g))


def paired(ea: np.ndarray, eb: np.ndarray, sel: np.ndarray) -> Dict[str, int]:
    """Baseline ``ea`` (run D), other ``eb`` (run C) over the voxels in ``sel``."""
    wa, wb = (ea > WRONG)[sel], (eb > WRONG)[sel]
    return dict(
        wrong_base=int(wa.sum()),
        wrong_other=int(wb.sum()),
        fixed=int((wa & ~wb).sum()),
        broken=int((~wa & wb).sum()),
        both_wrong=int((wa & wb).sum()),
    )


def pct(x: np.ndarray, q: Tuple[int, ...] = (50, 90, 99)) -> Dict[str, float]:
    return {str(k): float(np.percentile(x, k)) for k in q} if len(x) else {}


def main() -> None:
    W = World()
    L = W.L
    n = W.n
    regions = {
        "layer": np.ones(n, dtype=bool),
        "pilot": W.in_region.copy(),
        "rest": ~W.in_region,
    }
    S: Dict[str, Any] = dict(n_voxels=n, n_pilot=int(W.in_region.sum()))

    # ---- runs and checks
    info = {k: json.loads((RUNS / k / "run.json").read_text()) for k in "ABCD"}
    S["checks"] = {
        k: json.loads((RUNS / k / "check.json").read_text())
        for k in "AB"
        if (RUNS / k / "check.json").exists()
    }
    RC, haveC = W.read_map(RUNS / "C")
    RD, haveD = W.read_map(RUNS / "D")
    assert haveC.all() and haveD.all()
    eC, eD = W.err(RC), W.err(RD)
    paint = json.loads((RUNS / "C" / "paint.json").read_text())
    ents = paint["entries"]
    Re = np.array([e["R"] for e in ents], dtype=np.float64)
    seed_v = np.array([e["voxel"] for e in ents])
    stage = {int(k): v for k, v in paint["stage"].items()}
    scored: Dict[int, Dict[int, Tuple[float, float, bool]]] = {}
    passing: Dict[int, List[Tuple[int, float, float]]] = {}
    for g, rows in enumerate(paint["scores"]):
        for v, c, h, ok in rows:
            scored.setdefault(int(v), {})[g] = (c, h, ok)
            if ok:
                passing.setdefault(int(v), []).append((g, c, h))
    n_pass = np.array([len(passing.get(v, [])) for v in range(n)])

    # ---- evaluations
    evC = {k: int(v) for k, v in paint["evals"].items()}
    kD = KS.evals_of(info["D"]["stats"])
    stD = info["D"]["stats"]
    seed_recs = json.loads((RUNS / "C" / "seed_records.json").read_text())
    fit_recs = json.loads((RUNS / "C" / "fit_records.json").read_text())
    ev_vox_C = np.zeros(n, dtype=np.int64)
    # a seed's own score is its entry (not an evaluation); the flood scores are the evaluations
    for g, rows in enumerate(paint["scores"]):
        for v, c, h, ok in rows:
            if int(v) != ents[g]["voxel"]:
                ev_vox_C[int(v)] += 1
    for r in seed_recs:
        ev_vox_C[r["voxel"]] += r["evals"]
    for v, ne in fit_recs:
        ev_vox_C[v] += ne
    recD = json.loads((RUNS / "D" / "records.json").read_text())
    ev_vox_D = np.array([recD.get(str(v), {}).get("n_evals", 0) for v in range(n)])
    totC = int(sum(evC.values()))
    S["evals"] = dict(
        C=dict(
            total=totC,
            seed_paths=evC["seed_index"] + evC["index_rejected"],
            scoring=evC["score"],
            local_fits=evC["fit_accepted"] + evC["fit_rejected"],
            grid=evC["seed_grid"] + evC["grid_rejected"],
            by_stage=evC,
            by_region={k: int(ev_vox_C[m].sum()) for k, m in regions.items()},
        ),
        D=dict(
            total=kD["total"],
            seed_paths=int(stD["index_evals"]),
            scoring=0,
            local_fits=kD["total"]
            - int(stD["index_evals"])
            - int(stD.get("index_fallback_grid_evals", 0)),
            grid=int(stD.get("index_fallback_grid_evals", 0)),
            by_region={k: int(ev_vox_D[m].sum()) for k, m in regions.items()},
            records_sum=int(ev_vox_D.sum()),
        ),
        share_C_of_D=totC / kD["total"],
        pairing_s=dict(C=info["C"]["index_pairing_s"], D=info["D"]["index_pairing_s"]),
        cost_fn_eval_count_C=info["C"]["cost_fn_eval_count"],
        n_entries=len(ents),
        counts=paint["counts"],
    )

    # ---- wrong voxels
    wC, wD = eC > WRONG, eD > WRONG
    S["wrong"] = {
        k: dict(
            C=wil_dict(int(wC[m].sum()), int(m.sum())),
            D=wil_dict(int(wD[m].sum()), int(m.sum())),
            paired=paired(eD, eC, m),
            kinds_C={q: int((wC & m & (W.kind == q)).sum()) for q in KINDS},
            kinds_D={q: int((wD & m & (W.kind == q)).sum()) for q in KINDS},
            kinds_all={q: int((m & (W.kind == q)).sum()) for q in KINDS},
        )
        for k, m in regions.items()
    }
    for k in regions:
        w = S["wrong"][k]
        w["intervals_overlap"] = overlap(w["C"], w["D"])
    crossC = {int(v): W.crossing(int(v), RC[v]) for v in np.nonzero(wC)[0]}
    crossD = {int(v): W.crossing(int(v), RD[v]) for v in np.nonzero(wD)[0]}
    for k, m in regions.items():
        S["wrong"][k]["crossing_C"] = {
            q: sum(1 for v, c in crossC.items() if c and m[v] and W.kind[v] == q)
            for q in ("edge", "corner")
        }
        S["wrong"][k]["crossing_D"] = {
            q: sum(1 for v, c in crossD.items() if c and m[v] and W.kind[v] == q)
            for q in ("edge", "corner")
        }

    # ---- traces of the breaks (C wrong, D right) and the fixes (D wrong, C right)
    entry_grain = [L.grain_of(Re[i]) for i in range(len(ents))]

    def why(v: int) -> str:
        true_e = [i for i in range(len(ents)) if W.mis(Re[i], W.Rtrue[v]) <= WRONG]
        if not true_e:
            return "no entry of the true grain"
        if any(i in {g for g, _, _ in passing.get(v, [])} for i in true_e):
            return "true grain's entry passed and lost on cost"
        if any(i in scored.get(v, {}) for i in true_e):
            return "true grain's entry scored but did not pass"
        return "true grain's entry never scored the voxel"

    broken = [int(v) for v in np.nonzero(wC & ~wD)[0]]
    fixed = [int(v) for v in np.nonzero(wD & ~wC)[0]]
    brk_rows = []
    for v in broken[:20]:
        pg = min(passing.get(v, []), key=lambda t: (t[1], t[0]), default=None)
        brk_rows.append(
            dict(
                voxel=v,
                in_pilot=bool(W.in_region[v]),
                kind=str(W.kind[v]),
                err_C=float(eC[v]),
                crossing=crossC.get(v),
                stage=stage.get(v, "painted"),
                n_passing=int(n_pass[v]),
                painted_entry=None if pg is None else int(pg[0]),
                why=why(v),
            )
        )
    S["breaks"] = dict(
        n=len(broken),
        listed=brk_rows,
        why_all=(
            {w_: sum(why(v) == w_ for v in broken) for w_ in sorted({why(v) for v in broken})}
            if len(broken) <= 400
            else None
        ),
        kinds={q: sum(W.kind[v] == q for v in broken) for q in KINDS},
        crossing={
            q: sum(1 for v in broken if crossC.get(v) and W.kind[v] == q)
            for q in ("edge", "corner")
        },
    )
    # fixed voxels in two bins: the grain run D gave passes in paint and loses on cost, or not
    bins = {"passes_and_loses_on_cost": 0, "does_not_pass": 0, "no_entry_near_D_orientation": 0}
    for v in fixed:
        near = [i for i in range(len(ents)) if W.mis(Re[i], RD[v]) <= WRONG]
        if not near:
            bins["no_entry_near_D_orientation"] += 1
        elif any(i in {g for g, _, _ in passing.get(v, [])} for i in near):
            bins["passes_and_loses_on_cost"] += 1
        else:
            bins["does_not_pass"] += 1
    S["fixed_bins"] = dict(n=len(fixed), **bins)

    # ---- direct test on the contested voxels
    rows: List[Dict[str, Any]] = []
    for v in (int(x) for x in np.nonzero(n_pass >= 2)[0]):
        pg = passing[v]
        by_cost = sorted(pg, key=lambda t: (t[1], t[0]))
        by_hit = sorted(pg, key=lambda t: (-t[2], t[1], t[0]))
        c0, h0 = by_cost[0][0], by_hit[0][0]
        dmax = max(W.mis(Re[a[0]], Re[b[0]]) for a in pg for b in pg if a[0] < b[0])
        rows.append(
            dict(
                voxel=v,
                kind=str(W.kind[v]),
                pilot=bool(W.in_region[v]),
                distinct=bool(dmax > WRONG),
                cost_right=bool(W.mis(Re[c0], W.Rtrue[v]) <= WRONG),
                hit_right=bool(W.mis(Re[h0], W.Rtrue[v]) <= WRONG),
                cost_hit_differ=bool(W.mis(Re[c0], Re[h0]) > WRONG),
                true_among=bool(any(W.mis(Re[g], W.Rtrue[v]) <= WRONG for g, _, _ in pg)),
                margin=float(by_cost[1][1] - by_cost[0][1]),
                fit_settled=bool(stage.get(v) == "fit"),
            )
        )
    direct: Dict[str, Any] = dict(n_contested=len(rows))
    dup = [r for r in rows if not r["distinct"]]
    direct["same_grain_pairs"] = dict(
        n=len(dup),
        cost_right=sum(r["cost_right"] for r in dup),
        pilot=sum(r["pilot"] for r in dup),
    )
    for rname, rsel in (
        ("layer", lambda r: True),
        ("pilot", lambda r: r["pilot"]),
        ("rest", lambda r: not r["pilot"]),
    ):
        sub = [r for r in rows if r["distinct"] and rsel(r)]
        d: Dict[str, Any] = dict(
            n=len(sub),
            cost_right=wil_dict(sum(r["cost_right"] for r in sub), len(sub)),
            hit_right=wil_dict(sum(r["hit_right"] for r in sub), len(sub)),
            true_among=wil_dict(sum(r["true_among"] for r in sub), len(sub)),
            cost_hit_differ=wil_dict(sum(r["cost_hit_differ"] for r in sub), len(sub)),
        )
        for q in KINDS:
            ss = [r for r in sub if r["kind"] == q]
            d[q] = dict(
                n=len(ss),
                cost_right=wil_dict(sum(r["cost_right"] for r in ss), len(ss)),
                true_among=wil_dict(sum(r["true_among"] for r in ss), len(ss)),
            )
        direct[rname] = d
    direct["wrong_voxels"] = [r["voxel"] for r in rows if r["distinct"] and not r["cost_right"]][
        :40
    ]
    S["direct"] = direct

    # ---- the side arm
    sp = RUNS / "C" / "side" / "side.json"
    if sp.exists():
        side = json.loads(sp.read_text())
        z = np.load(RUNS / "C" / "side" / "side_fits.npz")
        Rs = RC.copy()
        Rs[z["voxel"]] = z["R"]
        eS = W.err(Rs)
        wS = eS > WRONG
        sv = np.array([r["voxel"] for r in side["rows"]], dtype=int)
        margins = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, float("inf")]
        curve = []
        for mg in margins:
            sel = [r for r in side["rows"] if r["margin"] < mg]
            vs = np.array([r["voxel"] for r in sel], dtype=int)
            curve.append(
                dict(
                    margin=mg if np.isfinite(mg) else None,
                    voxels=len(sel),
                    evals=int(sum(sum(r["fit_evals"]) for r in sel)),
                    fixed=int((wC[vs] & ~wS[vs]).sum()) if len(vs) else 0,
                    broken=int((~wC[vs] & wS[vs]).sum()) if len(vs) else 0,
                )
            )
        S["side"] = dict(
            n_voxels=len(sv),
            evals=side["evals"],
            wall_s=side["wall_s"],
            fixed=int((wC & ~wS).sum()),
            broken=int((~wC & wS).sum()),
            wrong_after=int(wS.sum()),
            by_region={
                k: dict(
                    fixed=int((wC & ~wS & m).sum()),
                    broken=int((~wC & wS & m).sum()),
                    voxels=int(m[sv].sum()),
                )
                for k, m in regions.items()
            },
            curve=curve,
            holds=bool((wC & ~wS).sum() <= (~wC & wS).sum()),
        )

    # ---- the grain list
    gids = set(L.gids)
    have_g = {g for g in entry_grain if g is not None}
    per_grain: Dict[int, int] = {}
    for g in entry_grain:
        if g is not None:
            per_grain[g] = per_grain.get(g, 0) + 1
    seed_ev = {r["voxel"]: r["evals"] for r in seed_recs}
    near = [(i, e) for i, e in enumerate(ents) if e["near_earlier"] is not None]
    loop_log = {int(r["voxel"]): r for r in paint["log"]}
    dup_rows = []
    for i, e in near:
        lr = loop_log.get(e["voxel"], {})
        fg = lr.get("fit_from")
        dup_rows.append(
            dict(
                entry=i,
                voxel=e["voxel"],
                in_pilot=bool(W.in_region[e["voxel"]]),
                near_deg=e["near_deg"],
                seed_path_evals=int(seed_ev.get(e["voxel"], 0)),
                flood_scores=e["flood_scores"],
                fit_from_grain=fg,
                fit_hit=lr.get("fit_hit"),
                best_hit_of_that_grain=None if fg is None else ents[fg]["best_hit"],
                fit_ratio=(
                    None
                    if fg is None or not ents[fg]["best_hit"]
                    else lr["fit_hit"] / ents[fg]["best_hit"]
                ),
                seed_hit=e["hit"],
                had_a_fit=fg is not None,
            )
        )
    S["grains"] = dict(
        n_true_grains=len(gids),
        n_entries=len(ents),
        kinds={k: sum(e["kind"] == k for e in ents) for k in ("index", "grid")},
        entries_matching_a_true_grain=len([g for g in entry_grain if g is not None]),
        entries_matching_none=sum(g is None for g in entry_grain),
        true_grains_with_entry=len(have_g),
        true_grains_without_entry=sorted(int(g) for g in gids - have_g),
        voxels_of_grains_without_entry=int(np.isin(L.gmap, list(gids - have_g)).sum()),
        grains_with_two_or_more=int(sum(v >= 2 for v in per_grain.values())),
        duplicates=dict(
            n=len(near),
            n_pilot=sum(r["in_pilot"] for r in dup_rows),
            seed_path_evals=sum(r["seed_path_evals"] for r in dup_rows),
            flood_scores=sum(r["flood_scores"] for r in dup_rows),
            rows=dup_rows,
        ),
    )

    # ---- voxels no grain paints; fallbacks; fit-settled voxels
    unpainted = [v for v in range(n) if n_pass[v] == 0 and v not in stage]
    failed_v = [v for v, s in stage.items() if s == "failed"]
    unp = [v for v in range(n) if n_pass[v] == 0 and stage.get(v) != "fit"]
    grid_rows = [
        r for r in paint["log"] if r.get("seed") in ("grid", "failed") or "grid_evals" in r
    ]
    S["no_grain"] = dict(
        unpainted_not_fit=len(unp),
        unpainted_listed=[
            dict(voxel=v, stage=stage.get(v, ""), kind=str(W.kind[v]), err=float(eC[v]))
            for v in unp[:40]
        ],
        failed=len(failed_v),
        grid_fallbacks_run=len(grid_rows),
        grid_fallbacks=[
            dict(
                voxel=r["voxel"],
                outcome=r.get("seed"),
                index_hit=r.get("index_hit"),
                grid_hit=r.get("grid_hit"),
                grid_evals=r.get("grid_evals"),
                err_C=float(eC[r["voxel"]]),
            )
            for r in grid_rows
        ],
    )
    fits_acc = [int(v) for v, s in stage.items() if s == "fit"]
    fit_rows = []
    for v in fits_acc:
        f = paint["fits"][str(v)]
        Rf = np.array(f["R"])
        lower = [t for t in passing.get(v, []) if t[1] < f["cost"]]
        best = min(passing.get(v, []), key=lambda t: (t[1], t[0]), default=None)
        fit_rows.append(
            dict(
                voxel=v,
                err_fit=W.mis(Rf, W.Rtrue[v]),
                grain_with_lower_cost=bool(lower),
                err_lowest_cost_grain=None if best is None else W.mis(Re[best[0]], W.Rtrue[v]),
            )
        )
    S["fit_settled"] = dict(
        n=len(fits_acc),
        n_rejected=int(paint["counts"]["fits_rejected"]),
        wrong_fit_kept=sum(r["err_fit"] > WRONG for r in fit_rows),
        wrong_lowest_cost_wins=sum(
            (r["err_lowest_cost_grain"] if r["err_lowest_cost_grain"] is not None else r["err_fit"])
            > WRONG
            for r in fit_rows
        ),
        rows=fit_rows,
    )

    # ---- accuracy without a fit
    both = (~wC) & (~wD)
    S["accuracy"] = {
        k: dict(
            n=int((both & m).sum()),
            C=pct(eC[both & m]),
            D=pct(eD[both & m]),
            C_max=float(eC[both & m].max()) if (both & m).any() else None,
            D_max=float(eD[both & m].max()) if (both & m).any() else None,
            D_over_0p1=int((eD[both & m] > 0.1).sum()),
            C_over_0p1=int((eC[both & m] > 0.1).sum()),
        )
        for k, m in regions.items()
    }

    # ---- wall
    S["wall"] = dict(
        C=dict(
            wall_s=info["C"]["wall_total_s"],
            label=info["C"]["timing_label"],
            pairing_s=info["C"]["index_pairing_s"],
            by_stage_s=info["C"]["stats"].get("paint_wall_s"),
        ),
        D=dict(
            wall_s=info["D"]["wall_total_s"],
            label=info["D"]["timing_label"],
            pairing_s=info["D"]["index_pairing_s"],
        ),
        A=dict(wall_s=info["A"]["wall_total_s"], label=info["A"]["timing_label"]),
        B=dict(wall_s=info["B"]["wall_total_s"], label=info["B"]["timing_label"]),
    )

    # ---- the fit log (runs B and D) and the seed record (run C)
    S["fit_log"] = {k: fit_log_summary(W, RUNS / k, e) for k, e in (("B", None), ("D", eD))}
    S["seed_record"] = seed_record_summary(W, seed_recs)

    (LAYER_OUT / "summary.json").write_text(json.dumps(S, indent=1, default=float) + "\n")
    tables(S)
    figure(S, W, eC, eD, both, regions)
    print(
        json.dumps(
            {k: S[k] for k in ("evals", "wrong", "breaks", "fixed_bins")}, indent=1, default=float
        )[:5000]
    )


# ------------------------------------------------------------------------------------------
BINS = [
    ("not at all", 0.0, 0.0),
    ("under 0.01", 0.0, 0.01),
    ("0.01 to 0.1", 0.01, 0.1),
    ("over 0.1", 0.1, 1e9),
]


def fit_log_summary(W: World, run: Path, err: Optional[np.ndarray]) -> Dict[str, Any]:
    rows = json.loads((run / "fit_log.json").read_text())
    if err is None:  # run B: the pilot region, errors from its own map
        R, have = W.read_map(run)
        err = W.err(R)
    ang = np.array([r["angle_deg"] for r in rows])
    same = np.array([r["exact_same"] for r in rows])
    ev = np.array([r["evals"] for r in rows])
    sc = np.array([r["start_cost"] for r in rows])
    ec = np.array([r["end_cost"] for r in rows])
    vox = np.array([r["voxel"] for r in rows])
    # the pilot run's voxels are the run's rows; map to the layer index for B
    if run.name == "B":
        vox = W.region[vox]
    cls = np.where(same, 0, np.where(ang < 0.01, 1, np.where(ang < 0.1, 2, 3)))
    names = [b[0] for b in BINS]
    out: Dict[str, Any] = dict(n_fits=len(rows), evals=int(ev.sum()))
    out["bins"] = {
        names[i]: dict(n=int((cls == i).sum()), evals=int(ev[cls == i].sum())) for i in range(4)
    }
    fall = sc - ec
    qs = np.quantile(sc, [0, 0.25, 0.5, 0.75, 1.0])
    fb = []
    for a, b in zip(qs[:-1], qs[1:]):
        m = (sc >= a) & (sc <= b)
        fb.append(
            dict(
                start_cost_lo=float(a),
                start_cost_hi=float(b),
                n=int(m.sum()),
                mean_fall=float(fall[m].mean()),
                median_fall=float(np.median(fall[m])),
            )
        )
    out["cost_fall_by_start_cost_quartile"] = fb
    out["fall_share_zero"] = float((fall <= 0).mean())
    out["start_hit_median"] = float(np.median([r["start_hit"] for r in rows]))
    moved = cls == 3
    final_wrong = err[vox[moved]] > WRONG
    out["moved_over_0p1"] = dict(
        n=int(moved.sum()), wrong_at_end=wil_dict(int(final_wrong.sum()), int(moved.sum()))
    )
    unmoved = ~moved
    out["not_moved_over_0p1"] = dict(
        n=int(unmoved.sum()),
        wrong_at_end=wil_dict(int((err[vox[unmoved]] > WRONG).sum()), int(unmoved.sum())),
    )
    out["angle_pct"] = pct(ang, (50, 90, 99))
    out["extra_start_evals"] = len(rows)
    return out


def seed_record_summary(W: World, recs: List[Dict[str, Any]]) -> Dict[str, Any]:
    ok = [r for r in recs if r["orientation"] is not None and r.get("candidate_taken") is not None]
    first = [r["refined"] and r["taken"] == r["refined"][0] for r in ok]
    moved = np.array(
        [W.mis(np.array(r["candidate_taken"]), np.array(r["orientation"])) for r in ok]
    )
    i_taken = []
    fall = []
    for r in ok:
        sc = r["strict_costs"]
        t = r["taken"]
        fall.append((sc[t] if t < len(sc) else np.nan) - r["cost"])
        i_taken.append(t < len(sc))
    fall_a = np.array(fall, dtype=float)
    return dict(
        n_seed_paths=len(recs),
        n_with_candidate=len(ok),
        lowest_strict_cost_is_taken=wil_dict(int(sum(first)), len(ok)),
        taken_among_costed=int(sum(i_taken)),
        finisher_moved_deg=pct(moved, (50, 90, 99)),
        finisher_moved_max=float(moved.max()) if len(moved) else None,
        cost_fall_median=float(np.nanmedian(fall_a)) if len(ok) else None,
        cost_fall_pct=pct(fall_a[~np.isnan(fall_a)], (50, 90)),
        evals_per_seed=float(np.mean([r["evals"] for r in recs])) if recs else None,
        finisher_evals_share=float(
            sum(sum(r["finisher_evals"]) for r in recs) / max(1, sum(r["evals"] for r in recs))
        ),
        both_finishers_same_end=int(
            sum(
                1
                for r in ok
                if len(r["finisher_R"]) == 2
                and W.mis(np.array(r["finisher_R"][0]), np.array(r["finisher_R"][1])) <= WRONG
            )
        ),
    )


def tables(S: Dict[str, Any]) -> None:
    T: Dict[str, str] = {}
    ev, w = S["evals"], S["wrong"]
    rows = []
    for name, e in (("C paint (key PaintGrains)", ev["C"]), ("D key-on BFS", ev["D"])):
        k = "C" if name.startswith("C") else "D"
        rows.append(
            dict(
                run=name,
                total=e["total"],
                seed=e["seed_paths"],
                score=e["scoring"],
                fit=e["local_fits"],
                grid=e["grid"],
                pilot=e["by_region"]["pilot"],
                rest=e["by_region"]["rest"],
                wrong=wil_str(w["layer"][k]),
            )
        )
    T["ar3b_evals"] = DT.markdown_table(
        rows,
        ["run", "total", "seed", "score", "fit", "grid", "pilot", "rest", "wrong"],
        formats={k: "," for k in ("total", "seed", "score", "fit", "grid", "pilot", "rest")},
    )
    prow = []
    for reg in ("layer", "pilot", "rest"):
        d = w[reg]
        prow.append(
            dict(
                region=reg,
                voxels=d["C"]["n"],
                C=wil_str(d["C"]),
                D=wil_str(d["D"]),
                fixed=d["paired"]["fixed"],
                broken=d["paired"]["broken"],
                both=d["paired"]["both_wrong"],
                intervals="overlap" if d["intervals_overlap"] else "do not overlap",
            )
        )
    T["ar3b_paired"] = DT.markdown_table(
        prow,
        ["region", "voxels", "C", "D", "fixed", "broken", "both", "intervals"],
        formats={"voxels": ","},
    )
    krow = []
    for reg in ("layer", "pilot", "rest"):
        d = w[reg]
        for q in KINDS:
            krow.append(
                dict(
                    region=reg,
                    kind=q,
                    voxels=d["kinds_all"][q],
                    C=d["kinds_C"][q],
                    D=d["kinds_D"][q],
                    crossC=d["crossing_C"].get(q, 0) if q != "interior" else 0,
                    crossD=d["crossing_D"].get(q, 0) if q != "interior" else 0,
                )
            )
    T["ar3b_kinds"] = DT.markdown_table(
        krow, ["region", "kind", "voxels", "C", "D", "crossC", "crossD"], formats={"voxels": ","}
    )
    dr = []
    for reg in ("layer", "pilot", "rest"):
        d = S["direct"][reg]
        dr.append(
            dict(
                region=reg,
                n=d["n"],
                lowest_cost=wil_str(d["cost_right"]),
                highest_hit=wil_str(d["hit_right"]),
                differ=wil_str(d["cost_hit_differ"]),
            )
        )
        for q in KINDS:
            dq = d[q]
            dr.append(
                dict(
                    region=f"{reg}, {q}",
                    n=dq["n"],
                    lowest_cost=wil_str(dq["cost_right"]),
                    highest_hit="",
                    differ="",
                )
            )
    T["ar3b_direct"] = DT.markdown_table(
        dr, ["region", "n", "lowest_cost", "highest_hit", "differ"]
    )
    if "side" in S:
        sr = [
            dict(
                m="all" if c["margin"] is None else f"{c['margin']:g}",
                voxels=c["voxels"],
                evals=c["evals"],
                fixed=c["fixed"],
                broken=c["broken"],
            )
            for c in S["side"]["curve"]
        ]
        T["ar3b_side"] = DT.markdown_table(
            sr, ["m", "voxels", "evals", "fixed", "broken"], formats={"evals": ","}
        )
    for k in ("B", "D"):
        fl = S["fit_log"][k]
        fr = [
            dict(
                moved=name,
                fits=b["n"],
                share=b["n"] / fl["n_fits"],
                evals=b["evals"],
                evals_share=b["evals"] / max(1, fl["evals"]),
            )
            for name, b in fl["bins"].items()
        ]
        T[f"ar3b_fitlog_{k}"] = DT.markdown_table(
            fr,
            ["moved", "fits", "share", "evals", "evals_share"],
            formats={"fits": ",", "evals": ",", "share": ".3f", "evals_share": ".3f"},
        )
    ac = [
        dict(
            region=reg,
            n=d["n"],
            C_p50=d["C"].get("50"),
            C_p90=d["C"].get("90"),
            C_p99=d["C"].get("99"),
            D_p50=d["D"].get("50"),
            D_p90=d["D"].get("90"),
            D_p99=d["D"].get("99"),
        )
        for reg, d in S["accuracy"].items()
    ]
    T["ar3b_accuracy"] = DT.markdown_table(
        ac,
        ["region", "n", "C_p50", "C_p90", "C_p99", "D_p50", "D_p90", "D_p99"],
        formats={k: ".4f" for k in ("C_p50", "C_p90", "C_p99", "D_p50", "D_p90", "D_p99")},
    )
    DT.write_tables(LAYER_OUT / "tables.md", T)


def figure(
    S: Dict[str, Any],
    W: World,
    eC: np.ndarray,
    eD: np.ndarray,
    both: np.ndarray,
    regions: Dict[str, np.ndarray],
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(15, 4))
    lo = 1e-4
    for e, lab in ((eC, "C paint"), (eD, "D key-on BFS")):
        x = np.sort(np.clip(e[both], lo, None))
        ax[0].plot(x, np.arange(1, len(x) + 1) / len(x), label=lab)
    ax[0].set_xscale("log")
    ax[0].set_xlabel("distance to the truth (degrees), voxels right in both")
    ax[0].set_ylabel("share of voxels")
    ax[0].legend()
    ax[0].set_title("Accuracy without a fit")
    for k in ("B", "D"):
        fl = S["fit_log"][k]
        names = list(fl["bins"])
        ax[1].bar(
            np.arange(4) + (0.0 if k == "B" else 0.4),
            [fl["bins"][nm]["n"] / fl["n_fits"] for nm in names],
            0.4,
            label=f"run {k}",
        )
    ax[1].set_xticks(np.arange(4) + 0.2)
    ax[1].set_xticklabels(["none", "<0.01", "0.01-0.1", ">0.1"])
    ax[1].set_xlabel("how far the BFS local fit moved (degrees)")
    ax[1].set_ylabel("share of fits")
    ax[1].legend()
    sp = S["seed_record"]
    ax[2].bar(
        ["lowest cost = taken", "other"],
        [
            sp["lowest_strict_cost_is_taken"]["k"],
            sp["n_with_candidate"] - sp["lowest_strict_cost_is_taken"]["k"],
        ],
    )
    ax[2].set_title("Seed paths (run C)")
    fig.tight_layout()
    fig.savefig(LAYER_OUT / "ar3_layer.png", dpi=110)


if __name__ == "__main__":
    main()
