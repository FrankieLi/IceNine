"""AR3 step 3: summary of runs P, R, R2, D on the gradient sample (judging points 0-8 of the plan).

  uv run python scripts/autoresearch/ar3/gradient_summary.py

Reads ``gradient/runs/{P,D}``, ``gradient/reading0.json.gz``, ``gradient/grad_stack.json`` (written
by ``grad_stack_note.py``), the render records of the gradient and the stored images
(``scripts/phase_d/results``), step 2's ``layer/summary.json`` and the step-2 maps of runs C and D
(the control). Also measures the local fit's trial steps at the run's settings with the library's
``QuaternionGrid`` (``step_sizes``). Writes ``gradient/summary.json``, ``gradient/tables.md`` and
the figure ``gradient/ar3_gradient.png``. No p-values. A rule's threshold is picked on the pick
set A (grains with an even id) and reported on the test set B (odd id) with the threshold fixed.
Ties in a pick take the more conservative threshold (stated where it applies).
"""

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import ar3_common as C
import gradient_common as G
import numpy as np

import doc_tables as DT  # noqa: E402
import stats as ST  # noqa: E402

RUNS = G.GRAD / "runs"
L2 = C.OUT / "layer" / "runs"
RENDER = C.ICE / "scripts" / "phase_d" / "results"
LEVELS = G.LEVELS
HIT_T = [round(0.90 + 0.01 * i, 2) for i in range(11)]
ZAP = 0.02  # degrees: "within 0.02 degrees of its own end"


def wil(k: int, n: int) -> Dict[str, Any]:
    lo, hi = ST.wilson(int(k), int(n)) if n else (float("nan"), float("nan"))
    return dict(k=int(k), n=int(n), lo=float(lo), hi=float(hi))


def wil_str(d: Dict[str, Any]) -> str:
    if not d["n"]:
        return "0 of 0"
    nd = 2 if d["hi"] < 0.01 else 1
    return f"{d['k']} of {d['n']} ({100 * d['lo']:.{nd}f} to {100 * d['hi']:.{nd}f}%)"


def overlap(a: Dict[str, Any], b: Dict[str, Any]) -> bool:
    return not (a["hi"] < b["lo"] or b["hi"] < a["lo"])


def four(err: np.ndarray) -> Dict[str, Any]:
    """Median, 95th percentile, share above 0.1 degrees, share wrong (error above 1 degree)."""
    n = len(err)
    if not n:
        return dict(n=0, median=float("nan"), p95=float("nan"), gt01=wil(0, 0), wrong=wil(0, 0))
    return dict(
        n=n,
        median=float(np.median(err)),
        p95=float(np.percentile(err, 95)),
        gt01=wil(int((err > 0.1).sum()), n),
        wrong=wil(int((err > G.WRONG).sum()), n),
    )


def four_row(label: str, err: np.ndarray) -> Dict[str, Any]:
    f = four(err)
    return dict(
        group=label,
        n=f["n"],
        median=f["median"],
        p95=f["p95"],
        gt01=wil_str(f["gt01"]),
        wrong=wil_str(f["wrong"]),
    )


FOUR_COLS = ["group", "n", "median", "p95", "gt01", "wrong"]
FOUR_FMT = {"median": ".4f", "p95": ".4f", "n": ","}
FOUR_LABELS = [
    "group",
    "voxels",
    "median error (deg)",
    "95th pct (deg)",
    "above 0.1 deg",
    "wrong (above 1 deg)",
]


def lvl_name(L: float) -> str:
    return f"level {L:g}"


class Data:
    def __init__(self) -> None:
        self.T = G.Truth()
        T = self.T
        self.n = T.n
        self.lev = T.level
        self.isA = (T.gmap % 2) == 0
        P = RUNS / "P"
        mp = np.load(P / "maps.npz")
        self.Pm = mp["P"].astype(np.float64)
        self.Rm = mp["R"].astype(np.float64)
        self.R2m = mp["R2"].astype(np.float64)
        self.paint = G.load_gz(P / "paint.json.gz")
        self.rR = G.load_gz(P / "refine_R.json.gz")
        self.rR2 = G.load_gz(P / "refine_R2.json.gz")
        self.runP = json.loads((P / "run.json").read_text())
        self.trace_check = json.loads((P / "trace_check.json").read_text())
        self.cfgP = (P / "run.config").read_text()
        self.runD = json.loads((RUNS / "D" / "run.json").read_text())
        self.Dlog = G.load_gz(RUNS / "D" / "fit_log.json.gz")
        self.Dm, _ = G.read_map(RUNS / "D" / "recon.mic.gz", self.n)
        self.r0 = np.array(G.load_gz(G.GRAD / "reading0.json.gz")["rows"])
        # strict cost and hit ratio of each voxel's true orientation (Reading 0), by voxel index
        self.cost_true = np.full(self.n, np.nan)
        self.hit_true = np.full(self.n, np.nan)
        self.cost_true[self.r0[:, 0].astype(int)] = self.r0[:, 3]
        self.hit_true[self.r0[:, 0].astype(int)] = self.r0[:, 4]
        assert not np.isnan(self.cost_true).any()
        self.step2 = json.loads((C.OUT / "layer" / "summary.json").read_text())
        self.render_grad = json.loads((RENDER / "render_grad_nodrop.json").read_text())
        self.render_full = json.loads((RENDER / "render_full.json").read_text())
        self.stack_note = json.loads((G.GRAD / "grad_stack.json").read_text())
        self.eP = T.err(self.Pm)
        self.eR = T.err(self.Rm)
        self.eR2 = T.err(self.R2m)
        self.eD = T.err(self.Dm)
        # source of each painted voxel
        pas: Dict[int, List[Tuple[int, float, bool]]] = {}
        for g, sc in enumerate(self.paint["scores"]):
            for v, c, h, ok in sc:
                if ok:
                    pas.setdefault(int(v), []).append((g, c, h))
        self.src = np.empty(self.n, dtype=object)
        stage = {int(k): s for k, s in self.paint["stage"].items()}
        for v in range(self.n):
            if str(v) in self.paint["fits"]:
                self.src[v] = "fit"
            elif v in pas:
                self.src[v] = "seed" if stage.get(v) == "seed" else "grain"
            elif str(v) in self.paint["failed"]:
                self.src[v] = "failed"
            else:
                self.src[v] = "none"
        self.vR = np.array([r["voxel"] for r in self.rR], dtype=int)
        self.painted = np.zeros(self.n, dtype=bool)
        self.painted[self.vR] = True


def sel_level(D: Data, L: float, extra: Optional[np.ndarray] = None) -> np.ndarray:
    m = D.lev == L
    return m if extra is None else (m & extra)


# --------------------------------------------------------------------------------------------
def point0(D: Data, S: Dict[str, Any], T: Dict[str, str]) -> None:
    a = D.r0
    rows = []
    out = {}
    for L in LEVELS:
        s = a[a[:, 2] == L]
        r = s[:, 6] / s[:, 4]
        k = int((r < 0.9).sum())
        out[str(L)] = dict(
            n=len(s),
            hit_true_median=float(np.median(s[:, 4])),
            hit_centre_median=float(np.median(s[:, 6])),
            ratio_median=float(np.median(r)),
            ratio_p05=float(np.percentile(r, 5)),
            below_09=wil(k, len(s)),
            cost_true_median=float(np.median(s[:, 3])),
            cost_centre_median=float(np.median(s[:, 5])),
        )
        rows.append(
            dict(
                level=lvl_name(L),
                n=len(s),
                ht=out[str(L)]["hit_true_median"],
                hc=out[str(L)]["hit_centre_median"],
                rm=out[str(L)]["ratio_median"],
                r5=out[str(L)]["ratio_p05"],
                below=wil_str(out[str(L)]["below_09"]),
                ct=out[str(L)]["cost_true_median"],
                cc=out[str(L)]["cost_centre_median"],
            )
        )
    S["reading0"] = out
    T["ar3c_reading0"] = DT.markdown_table(
        rows,
        ["level", "n", "ht", "hc", "rm", "r5", "below", "ct", "cc"],
        formats={k: ".3f" for k in ("ht", "hc", "rm", "r5", "ct", "cc")} | {"n": ","},
        labels=None,
    )


def point1(D: Data, S: Dict[str, Any], T: Dict[str, str]) -> None:
    rows = [four_row("all voxels", D.eP)]
    for L in LEVELS:
        rows.append(four_row(lvl_name(L), D.eP[D.lev == L]))
    S["start"] = {
        "all": four(D.eP),
        **{str(L): four(D.eP[D.lev == L]) for L in LEVELS},
    }
    T["ar3c_start"] = DT.markdown_table(rows, FOUR_COLS, FOUR_FMT)
    # sources by level
    srows = []
    srcs = ["grain", "fit", "seed", "failed", "none"]
    S["sources"] = {}
    for name, m in [("all voxels", np.ones(D.n, bool))] + [
        (lvl_name(L), D.lev == L) for L in LEVELS
    ]:
        cnt = {s: int(((D.src == s) & m).sum()) for s in srcs}
        S["sources"][name] = cnt
        srows.append(dict(group=name, n=int(m.sum()), **cnt))
    T["ar3c_sources"] = DT.markdown_table(srows, ["group", "n"] + srcs, {"n": ","})
    # entries per true grain
    ent = D.paint["entries"]
    T_ = D.T
    per: Dict[int, int] = {}
    for e in ent:
        per[int(T_.gmap[e["voxel"]])] = per.get(int(T_.gmap[e["voxel"]]), 0) + 1
    erows = []
    S["entries"] = {}
    for L in LEVELS:
        gs = [g for g in T_.gids if T_.level_of[g] == L]
        c = np.array([per.get(g, 0) for g in gs])
        kinds = [
            e for e in ent if T_.level_of[int(T_.gmap[e["voxel"]])] == L and e["kind"] == "grid"
        ]
        S["entries"][str(L)] = dict(
            grains=len(gs),
            entries=int(c.sum()),
            mean=float(c.mean()),
            max=int(c.max()),
            more_than_one=wil(int((c > 1).sum()), len(gs)),
            none=int((c == 0).sum()),
            grid_entries=len(kinds),
        )
        erows.append(
            dict(
                level=lvl_name(L),
                grains=len(gs),
                entries=int(c.sum()),
                mean=float(c.mean()),
                mx=int(c.max()),
                multi=wil_str(S["entries"][str(L)]["more_than_one"]),
                none=int((c == 0).sum()),
                grid=len(kinds),
            )
        )
    c_all = np.array([per.get(g, 0) for g in T_.gids])
    S["entries"]["all"] = dict(
        grains=len(c_all),
        entries=int(c_all.sum()),
        mean=float(c_all.mean()),
        max=int(c_all.max()),
        n_entries_total=len(ent),
    )
    T["ar3c_entries"] = DT.markdown_table(
        erows,
        ["level", "grains", "entries", "mean", "mx", "multi", "none", "grid"],
        {"mean": ".2f"},
    )
    # entries within 1 degree of an earlier entry (paint.json.gz `near_earlier`), the true grain
    # of the two entries' voxels, the flood size, and the fit-settled branch (P/run.json counts)
    near = [e for e in ent if e["near_earlier"] is not None]
    same = sum(
        1 for e in near if T_.gmap[e["voxel"]] == T_.gmap[ent[int(e["near_earlier"])]["voxel"]]
    )
    cnt = D.runP["counts"]
    S["entries_detail"] = dict(
        near_earlier=len(near),
        near_same_true_grain=int(same),
        near_other_true_grain=len(near) - int(same),
        flood_passes_median=float(np.median([e["flood_passes"] for e in ent])),
        fit_branch_tried=int(cnt["fits_accepted"] + cnt["fits_rejected"]),
        fit_branch_accepted=int(cnt["fits_accepted"]),
        seed_share_of_voxels=len(ent) / D.n,
    )
    ev = D.paint["evals"]
    tot = int(sum(ev.values()))
    n_idx = sum(1 for e in ent if e["kind"] == "index")
    n_grid = sum(1 for e in ent if e["kind"] == "grid")
    S["seed_cost"] = dict(
        entries=len(ent),
        index_entries=n_idx,
        grid_entries=n_grid,
        evals_per_index_entry=ev.get("seed_index", 0) / max(n_idx, 1),
        evals_per_grid_entry=ev.get("seed_grid", 0) / max(n_grid, 1),
        seed_share_of_total=(
            ev.get("seed_index", 0)
            + ev.get("seed_grid", 0)
            + ev.get("index_rejected", 0)
            + ev.get("grid_rejected", 0)
        )
        / tot,
        grid_share_of_total=(ev.get("seed_grid", 0) + ev.get("grid_rejected", 0)) / tot,
        per_voxel=tot / D.n,
    )
    sc = S["seed_cost"]
    T["ar3c_seed_cost"] = DT.markdown_table(
        [
            dict(what="entries from the seed path", n=n_idx, per=sc["evals_per_index_entry"]),
            dict(
                what="entries from the grid search (fallback)",
                n=n_grid,
                per=sc["evals_per_grid_entry"],
            ),
        ],
        ["what", "n", "per"],
        {"per": ",.0f", "n": ","},
        labels=None,
    )
    # the 172 grid-fallback entries alone (without the 4 rejected grid searches of failed voxels)
    sc["grid_entries_only_share"] = ev.get("seed_grid", 0) / tot
    st_d = D.runD["stats"]
    sc["D_fallbacks"] = int(st_d["index_fallback_seed_test"])
    sc["D_fallback_grid_evals"] = int(st_d["index_fallback_grid_evals"])
    # P against step 2's whole-layer paint (one orientation per grain), with and without the
    # grid fallback
    s2 = D.step2["evals"]
    no_grid = tot - ev.get("seed_grid", 0) - ev.get("grid_rejected", 0)
    S["vs_step2_paint"] = dict(
        step2_evals=int(s2["C"]["total"]),
        step2_entries=int(s2["n_entries"]),
        step2_key_on_bfs_evals=int(s2["D"]["total"]),
        P_evals=tot,
        P_evals_without_grid=int(no_grid),
        ratio=tot / s2["C"]["total"],
        ratio_without_grid=no_grid / s2["C"]["total"],
        entries_ratio=len(ent) / s2["n_entries"],
    )
    S["evals_P_by_stage"] = dict(ev, total=tot)
    T["ar3c_evals_stage"] = DT.markdown_table(
        [dict(stage=k, evals=v) for k, v in ev.items()] + [dict(stage="total", evals=tot)],
        ["stage", "evals"],
        {"evals": ","},
    )


def paired_vs_start(D: Data, S: Dict[str, Any], T: Dict[str, str]) -> None:
    m = D.painted
    rows = [
        four_row("start, painted voxels", D.eP[m]),
        four_row("after R", D.eR[m]),
        four_row("D (BFS)", D.eD[m]),
    ]
    S["refine_painted"] = {"start": four(D.eP[m]), "R": four(D.eR[m]), "D": four(D.eD[m])}
    T["ar3c_refine_all"] = DT.markdown_table(rows, FOUR_COLS, FOUR_FMT)
    prow = []
    S["refine_by_level"] = {}
    for L in LEVELS + (None,):
        mm = m if L is None else (m & (D.lev == L))
        nm = int(mm.sum())
        closer = int((D.eR[mm] < D.eP[mm] - G.HURT).sum())
        hurt = int((D.eR[mm] > D.eP[mm] + G.HURT).sum())
        wP, wR, wD = (D.eP[mm] > 1).sum(), (D.eR[mm] > 1).sum(), (D.eD[mm] > 1).sum()
        key = "all" if L is None else str(L)
        S["refine_by_level"][key] = dict(
            n=nm,
            closer=wil(closer, nm),
            hurt=wil(hurt, nm),
            start=four(D.eP[mm]),
            R=four(D.eR[mm]),
            D=four(D.eD[mm]),
            wrong_start=int(wP),
            wrong_R=int(wR),
            wrong_D=int(wD),
            R_fixed_vs_D=int(((D.eD[mm] > 1) & (D.eR[mm] <= 1)).sum()),
            R_broken_vs_D=int(((D.eD[mm] <= 1) & (D.eR[mm] > 1)).sum()),
        )
        d = S["refine_by_level"][key]
        prow.append(
            dict(
                level="all" if L is None else lvl_name(L),
                n=nm,
                closer=wil_str(d["closer"]),
                hurt=wil_str(d["hurt"]),
                med_s=d["start"]["median"],
                med_r=d["R"]["median"],
                med_d=d["D"]["median"],
                p95_s=d["start"]["p95"],
                p95_r=d["R"]["p95"],
                p95_d=d["D"]["p95"],
                w_s=int(wP),
                w_r=int(wR),
                w_d=int(wD),
            )
        )
    T["ar3c_refine_level"] = DT.markdown_table(
        prow,
        [
            "level",
            "n",
            "closer",
            "hurt",
            "med_s",
            "med_r",
            "med_d",
            "p95_s",
            "p95_r",
            "p95_d",
            "w_s",
            "w_r",
            "w_d",
        ],
        {k: ".4f" for k in ("med_s", "med_r", "med_d", "p95_s", "p95_r", "p95_d")} | {"n": ","},
    )
    # by the size of the start's error: does the fit help where it is needed?
    edges = [0.0, 0.05, 0.1, 0.2, 0.5, 1e9]
    brow = []
    S["refine_by_start_error"] = {}
    for lo, hi in zip(edges[:-1], edges[1:]):
        mm = m & (D.eP >= lo) & (D.eP < hi)
        nm = int(mm.sum())
        if not nm:
            continue
        closer = int((D.eR[mm] < D.eP[mm] - G.HURT).sum())
        hurt = int((D.eR[mm] > D.eP[mm] + G.HURT).sum())
        moved = int((D.T.mis(D.Pm[mm], D.Rm[mm]) > ZAP).sum())
        name = f"{lo:g} to {hi:g}" if hi < 1e8 else f"{lo:g} and above"
        S["refine_by_start_error"][name] = dict(
            n=nm, moved=wil(moved, nm), closer=wil(closer, nm), hurt=wil(hurt, nm)
        )
        brow.append(
            dict(
                err=name,
                n=nm,
                moved=wil_str(wil(moved, nm)),
                closer=wil_str(wil(closer, nm)),
                hurt=wil_str(wil(hurt, nm)),
            )
        )
    T["ar3c_refine_bins"] = DT.markdown_table(
        brow, ["err", "n", "moved", "closer", "hurt"], {"n": ","}
    )
    evP = D.runP["evals_total"]
    S["evals_total"] = dict(
        P=evP,
        R=D.runP["evals_R"],
        R2=D.runP["evals_R2"],
        P_plus_R=evP + D.runP["evals_R"],
        D=D.runD["evals_total"],
    )
    T["ar3c_evals"] = DT.markdown_table(
        [
            dict(arm="P (paint)", evals=evP),
            dict(arm="R (one fit on each painted voxel)", evals=D.runP["evals_R"]),
            dict(arm="P plus R", evals=evP + D.runP["evals_R"]),
            dict(arm="R2 (second fit, every fourth voxel)", evals=D.runP["evals_R2"]),
            dict(arm="D (BFS, indexed seeds)", evals=D.runD["evals_total"]),
        ],
        ["arm", "evals"],
        {"evals": ","},
    )


ERR_EDGES = [0.0, 0.05, 0.1, 0.2, 0.5, 1e9]


def bin_name(lo: float, hi: float) -> str:
    return f"{lo:g} to {hi:g}" if hi < 1e8 else f"{lo:g} and above"


def fit_cause(D: Data, S: Dict[str, Any], T: Dict[str, str]) -> None:
    """Can the strict cost tell the painted start from the truth, and did the fit get there?

    Per painted voxel: the strict cost of its true orientation (Reading 0) against the cost of its
    start (refine_R `start_cost`), and whether R's fit ended at a lower cost than its start. The
    truth enters this analysis only (no decision)."""
    v = D.vR
    sc = np.array([r["start_cost"] for r in D.rR])
    ec = np.array([r["end_cost"] for r in D.rR])
    sh = np.array([r["start_hit"] for r in D.rR])
    ct, ht = D.cost_true[v], D.hit_true[v]
    lev, es = D.lev[v], D.eP[v]
    groups: List[Tuple[str, np.ndarray]] = [("all painted voxels", np.ones(len(v), bool))]
    groups += [(lvl_name(L), lev == L) for L in LEVELS]
    groups += [
        (f"start error {bin_name(lo, hi)}", (es >= lo) & (es < hi))
        for lo, hi in zip(ERR_EDGES[:-1], ERR_EDGES[1:])
    ]
    rows = []
    S["fit_cause"] = {}
    for name, m in groups:
        nm = int(m.sum())
        lower = int((ct[m] < sc[m]).sum())
        d = dict(
            n=nm,
            truth_lower=wil(lower, nm),
            truth_lower_share=lower / nm,
            equal=int((ct[m] == sc[m]).sum()),
            truth_higher=int((ct[m] > sc[m]).sum()),
            fit_lowered=wil(int((ec[m] < sc[m]).sum()), nm),
            median_gap=float(np.median(sc[m] - ct[m])),
        )
        S["fit_cause"][name] = d
        rows.append(
            dict(
                group=name,
                n=nm,
                truth_lower=wil_str(d["truth_lower"]),
                equal=d["equal"],
                truth_higher=d["truth_higher"],
                fit_lowered=wil_str(d["fit_lowered"]),
                med_gap=d["median_gap"],
            )
        )
    gt = es > 0.1
    S["fit_cause_extra"] = dict(
        truth_hit_higher=int((ht > sh).sum()),
        start_error_median=float(np.median(es)),
        gt01=int(gt.sum()),
        gt01_truth_lower=int((ct[gt] < sc[gt]).sum()),
    )
    T["ar3c_truth_vs_start"] = DT.markdown_table(
        rows,
        ["group", "n", "truth_lower", "equal", "truth_higher", "fit_lowered", "med_gap"],
        {"n": ",", "med_gap": ".4f"},
    )


def cfg_value(text: str, key: str) -> float:
    import re

    m = re.search(rf"(?m)^{key}\s+(\S+)", text)
    assert m, key
    return float(m.group(1))


def step_sizes(D: Data, S: Dict[str, Any], T: Dict[str, str]) -> None:
    """The trial steps of the local fit at the run's settings, measured with the library.

    ``variance_minimizing_optimize`` (C++ ``AdaptiveSamplingZeroTemp``) at box width
    LocalOrientationGridRadius / 2^MinLocalResolution. Block k searches a subregion of radius
    r_k (r_1 = tan(box) / sqrt(48); doubled after a failed block, capped at the box) with
    max(ceil(r_k / 0.5 deg)^2.7, 10) trial steps; a trial is ``get_near_identity_point`` of an
    offset drawn uniformly in +-tan(r_k) / sqrt(12) per axis. Before block k >= 2 the search
    restarts from the start plus ``get_near_identity_point`` of an offset uniform in +-r_k. The
    blocks are listed as in a fit whose blocks all fail (no budget extension), until MaxMCSteps.
    Step and offset angles: 10,000 draws each from default_rng(0) (a measurement of the
    distribution, not part of any run)."""
    import math

    from icenine.sampling import QuaternionGrid

    radius_deg = cfg_value(D.cfgP, "LocalOrientationGridRadius")
    min_res = cfg_value(D.cfgP, "MinLocalResolution")
    max_steps = cfg_value(D.cfgP, "MaxMCSteps")
    box = math.radians(radius_deg) / 2**min_res
    res = math.radians(0.5)
    g = QuaternionGrid()
    rng = np.random.default_rng(0)
    n_draw = 10_000

    def angles(half: float) -> np.ndarray:
        xyz = rng.uniform(-half, half, size=(n_draw, 3))
        w = np.array([abs(g.get_near_identity_point(*p)[0]) for p in xyz])
        return np.degrees(2 * np.arccos(np.minimum(w, 1.0)))

    zero = float(np.degrees(2 * np.arccos(min(abs(g.get_near_identity_point(0, 0, 0)[0]), 1.0))))
    sub = math.tan(box) / math.sqrt(48.0)
    total, ev, k = 0, 1, 0
    rows = []
    blocks = []
    while total < max_steps:
        k += 1
        n = max(int(math.ceil(sub / res) ** 2.7), 10)
        st = angles(math.tan(sub) / math.sqrt(12.0))
        rs = angles(sub) if k > 1 else None
        ev += 1 + n
        b = dict(
            block=k,
            radius_deg=math.degrees(sub),
            steps=n,
            evals_at_end=ev,
            step_median=float(np.median(st)),
            step_p05=float(np.percentile(st, 5)),
            restart_median=None if rs is None else float(np.median(rs)),
        )
        blocks.append(b)
        rows.append(
            dict(
                block=k,
                radius_deg=b["radius_deg"],
                steps=n,
                evals_at_end=ev,
                step_med=b["step_median"],
                step_p05=b["step_p05"],
                restart_med="none" if rs is None else f"{b['restart_median']:.1f}",
            )
        )
        total += n
        sub = min(2.0 * sub, box)
    ev_fit = ev + 1  # the final overlap evaluation
    tr = [r["trace"] for r in D.rR]
    ev_pat = [b["evals_at_end"] for b in blocks]
    four_failed = sum(1 for t in tr if t["ev"] == ev_pat and all(t["fl"]))
    ev_R = np.array([r["evals"] for r in D.rR])
    S["fit_steps"] = dict(
        LocalOrientationGridRadius=radius_deg,
        MinLocalResolution=min_res,
        MaxMCSteps=max_steps,
        box_deg=math.degrees(box),
        identity_offset_deg=zero,
        n_draws_per_row=n_draw,
        blocks=blocks,
        evals_all_blocks_fail=ev_fit,
        R_median_evals=float(np.median(ev_R)),
        R_fits=len(tr),
        R_fits_all_blocks_failed=int(four_failed),
    )
    assert ev_fit == S["fit_steps"]["R_median_evals"]
    T["ar3c_steps"] = DT.markdown_table(
        rows,
        ["block", "radius_deg", "steps", "evals_at_end", "step_med", "step_p05", "restart_med"],
        {"radius_deg": ".3f", "step_med": ".2f", "step_p05": ".2f"},
    )


# --------------------------------------------------------------------------------------------
def pick_skip(D: Data, S: Dict[str, Any], T: Dict[str, str]) -> Dict[str, Any]:
    rows = D.rR
    v = D.vR
    ev = np.array([r["evals"] for r in rows])
    sh = np.array([r["start_hit"] for r in rows])
    eb = np.array([r["entry_best_hit"] for r in rows])
    sc = np.array([r["start_cost"] for r in rows])
    es, ee = D.eP[v], D.eR[v]
    ratio = sh / eb
    hurt_if = es > ee + G.HURT
    isA = D.isA[v]
    res: Dict[str, Any] = {}

    def sweep(kind: str, grid: List[float], skip_fn: Any) -> None:
        sweep_rows = []
        for t in grid:
            sk = skip_fn(t)
            a = sk & isA
            nA = int(isA.sum())
            hA = int((a & hurt_if).sum())
            sweep_rows.append(
                dict(
                    t=float(t),
                    A_skipped_evals=int(ev[a].sum()),
                    A_skipped=int(a.sum()),
                    A_hurt=hA,
                    A_hurt_share=hA / nA,
                )
            )
        feas = [r for r in sweep_rows if r["A_hurt_share"] <= 0.01 and r["A_skipped"] > 0]
        best = None
        if feas:
            mx = max(r["A_skipped_evals"] for r in feas)
            # tie rule: of the thresholds with the most skipped evaluations, the more conservative
            tied = [r for r in feas if r["A_skipped_evals"] == mx]
            best = (max if kind == "hit" else min)(tied, key=lambda r: r["t"])
        res[kind] = dict(sweep=sweep_rows, pick=best)

    sweep("hit", HIT_T, lambda t: ratio >= t)
    cgrid = [float(x) for x in np.percentile(sc[isA], list(range(10, 101, 10)))]
    sweep("cost", cgrid, lambda c: sc <= c)
    # on B
    nB = int((~isA).sum())
    tot_ev_B = int(ev[~isA].sum())
    out_rows = []
    S["skip"] = {}
    for kind, fn in (("hit", lambda t: ratio >= t), ("cost", lambda c: sc <= c)):
        p = res[kind]["pick"]
        if p is None:
            S["skip"][kind] = dict(pick=None)
            continue
        sk = fn(p["t"])
        b = sk & ~isA
        hB = int((b & hurt_if).sum())
        # the resulting map on B's voxels: skipped painted voxels keep the start
        err_map = D.eP.copy()
        err_map[v[~sk]] = ee[~sk]  # fitted painted voxels
        err_map[v[sk]] = es[sk]
        mB = ~D.isA
        # voxels never refined keep P's value (D.eP) already
        S["skip"][kind] = dict(
            pick=p,
            B_painted=nB,
            B_skipped=int(b.sum()),
            B_skipped_evals=int(ev[b].sum()),
            B_total_fit_evals=tot_ev_B,
            B_hurt=wil(hB, nB),
            B_four=four(err_map[mB]),
            B_four_refit_all=four(D.eR[mB]),
            by_level={
                str(L): dict(
                    n=int(((D.lev[v] == L) & ~isA).sum()),
                    skipped=int((b & (D.lev[v] == L)).sum()),
                    hurt=int((b & hurt_if & (D.lev[v] == L)).sum()),
                    skipped_evals=int(ev[b & (D.lev[v] == L)].sum()),
                )
                for L in LEVELS
            },
        )
        d = S["skip"][kind]
        # the most permissive end of the grid: the lowest t; for c the 100th percentile of A's
        # start costs (A's largest start cost)
        edge = res[kind]["sweep"][0]["t"] if kind == "hit" else res[kind]["sweep"][-1]["t"]
        d["pick_at_grid_edge"] = bool(p["t"] == edge)
        out_rows.append(
            dict(
                rule="start hit ratio >= t x entry best" if kind == "hit" else "start cost <= c",
                threshold=p["t"],
                A_share=f"{p['A_skipped_evals'] / ev[isA].sum():.3f}",
                B_skipped=f"{d['B_skipped']} of {nB}",
                B_evals=f"{d['B_skipped_evals']:,} of {tot_ev_B:,}",
                B_hurt=wil_str(d["B_hurt"]),
                B_med=d["B_four"]["median"],
                B_p95=d["B_four"]["p95"],
                B_gt01=wil_str(d["B_four"]["gt01"]),
                B_wrong=wil_str(d["B_four"]["wrong"]),
            )
        )
    S["skip_sweeps"] = {k: res[k]["sweep"] for k in res}
    T["ar3c_skip"] = DT.markdown_table(
        out_rows,
        [
            "rule",
            "threshold",
            "A_share",
            "B_skipped",
            "B_evals",
            "B_hurt",
            "B_med",
            "B_p95",
            "B_gt01",
            "B_wrong",
        ],
        {"threshold": ".4g", "B_med": ".4f", "B_p95": ".4f"},
    )
    lrows = []
    for kind in ("hit", "cost"):
        d = S["skip"].get(kind, {})
        if not d.get("pick"):
            continue
        for L in LEVELS:
            b = d["by_level"][str(L)]
            lrows.append(dict(rule=kind, level=lvl_name(L), **b))
    T["ar3c_skip_level"] = DT.markdown_table(
        lrows, ["rule", "level", "n", "skipped", "hurt", "skipped_evals"], {"skipped_evals": ","}
    )
    return res


def block_stop(D: Data, S: Dict[str, Any], T: Dict[str, str]) -> None:
    rows = D.rR
    v = D.vR
    tr = [r["trace"] for r in rows]
    ev_tot = np.array([r["evals"] for r in rows])
    ee = D.eR[v]
    isA = D.isA[v]
    lev = D.lev[v]
    d0 = D.T.mis(D.Pm[v], D.Rm[v])  # start to end
    # converged point: first block within ZAP of the end; the start counts when it already is
    after = np.zeros(len(rows))
    conv_ev = np.zeros(len(rows))
    for i, t in enumerate(tr):
        if d0[i] <= ZAP:
            ce = 1
        else:
            j = next((j for j, d in enumerate(t["de"]) if d <= ZAP), None)
            ce = t["ev"][j] if j is not None else ev_tot[i]
        conv_ev[i] = ce
        after[i] = max(ev_tot[i] - ce, 0)
    S["converged"] = {
        "all": dict(
            total_evals=int(ev_tot.sum()),
            evals_after=int(after.sum()),
            share=float(after.sum() / ev_tot.sum()),
            median_share=float(np.median(after / ev_tot)),
            start_within=wil(int((d0 <= ZAP).sum()), len(d0)),
            n_fits=len(rows),
        )
    }
    crow = []
    for L in LEVELS:
        m = lev == L
        S["converged"][str(L)] = dict(
            n=int(m.sum()),
            share=float(after[m].sum() / ev_tot[m].sum()),
            median_share=float(np.median(after[m] / ev_tot[m])),
            median_evals=float(np.median(ev_tot[m])),
            start_within=wil(int((d0[m] <= ZAP).sum()), int(m.sum())),
        )
        d = S["converged"][str(L)]
        crow.append(
            dict(
                level=lvl_name(L),
                n=d["n"],
                med_evals=d["median_evals"],
                share=d["share"],
                med_share=d["median_share"],
                start_within=wil_str(d["start_within"]),
            )
        )
    a = S["converged"]["all"]
    crow.append(
        dict(
            level="all",
            n=a["n_fits"],
            med_evals=float(np.median(ev_tot)),
            share=a["share"],
            med_share=a["median_share"],
            start_within=wil_str(a["start_within"]),
        )
    )
    T["ar3c_converged"] = DT.markdown_table(
        crow,
        ["level", "n", "med_evals", "share", "med_share", "start_within"],
        {"med_evals": ".0f", "share": ".3f", "med_share": ".3f", "n": ","},
    )

    # the two stop rules, simulated on the trace
    def stop_first_fail(i: int, _c: float) -> Tuple[float, float]:
        t = tr[i]
        f = next((j for j, x in enumerate(t["fl"]) if x), None)
        if f is None:
            return float(ev_tot[i]), float(ee[i])
        return float(min(t["ev"][f] + 1, ev_tot[i])), float(t["te"][f])

    def stop_cost(i: int, c: float) -> Tuple[float, float]:
        t = tr[i]
        j = next((j for j, x in enumerate(t["bc"]) if x <= c), None)
        if j is None:
            return float(ev_tot[i]), float(ee[i])
        return float(min(t["ev"][j] + 1, ev_tot[i])), float(t["te"][j])

    endc = np.array([r["end_cost"] for r in rows])
    cgrid = [float(x) for x in np.percentile(endc[isA], list(range(10, 101, 10)))]
    S["stop"] = {}
    sweeps: Dict[str, Any] = {}
    prow = []
    for kind, fn, grid in (
        ("first_fail", stop_first_fail, [0.0]),
        ("cost", stop_cost, cgrid),
    ):
        sw = []
        for c in grid:
            ev_used = np.zeros(len(rows))
            err_stop = np.zeros(len(rows))
            for i in range(len(rows)):
                ev_used[i], err_stop[i] = fn(i, c)
            hurt = err_stop > ee + G.HURT
            a_ = isA
            sw.append(
                dict(
                    t=float(c),
                    A_saved=int((ev_tot[a_] - ev_used[a_]).sum()),
                    A_hurt=int(hurt[a_].sum()),
                    A_hurt_share=float(hurt[a_].mean()),
                    ev_used=ev_used,
                    err_stop=err_stop,
                    hurt=hurt,
                )
            )
        feas = [r for r in sw if r["A_hurt_share"] <= 0.01]
        pick = None
        if feas:
            mx = max(r["A_saved"] for r in feas)
            tied = [r for r in feas if r["A_saved"] == mx]
            pick = min(tied, key=lambda r: r["t"])  # tie: the lower (stricter) cost threshold
        sweeps[kind] = [
            {k: x for k, x in r.items() if k not in ("ev_used", "err_stop", "hurt")} for r in sw
        ]
        if pick is None:
            S["stop"][kind] = dict(pick=None)
            continue
        b = ~isA
        saved = ev_tot[b] - pick["ev_used"][b]
        err_map = D.eR.copy()
        err_map[v] = np.where(isA, ee, pick["err_stop"])
        hB = int(pick["hurt"][b].sum())
        mB = ~D.isA
        S["stop"][kind] = dict(
            pick={k: x for k, x in pick.items() if k not in ("ev_used", "err_stop", "hurt")},
            B_fits=int(b.sum()),
            B_saved=int(saved.sum()),
            B_total=int(ev_tot[b].sum()),
            B_hurt=wil(hB, int(b.sum())),
            B_four=four(err_map[mB]),
            by_level={
                str(L): dict(
                    n=int((b & (lev == L)).sum()),
                    saved=int((ev_tot - pick["ev_used"])[b & (lev == L)].sum()),
                    total=int(ev_tot[b & (lev == L)].sum()),
                    hurt=int(pick["hurt"][b & (lev == L)].sum()),
                )
                for L in LEVELS
            },
        )
        d = S["stop"][kind]
        if kind == "cost":
            # which grid point the pick is: percentiles 10 to 100 of A's end costs (100 = A's
            # largest end cost, the grid's edge)
            d["pick_percentile"] = 10 * (grid.index(pick["t"]) + 1)
            d["pick_at_grid_edge"] = bool(pick["t"] == grid[-1])
        prow.append(
            dict(
                rule=(
                    "stop at the first failed block"
                    if kind == "first_fail"
                    else "stop when best cost <= c"
                ),
                threshold="none" if kind == "first_fail" else f"{pick['t']:.4f}",
                A_saved_share=f"{pick['A_saved'] / ev_tot[isA].sum():.3f}",
                B_saved=f"{d['B_saved']:,} of {d['B_total']:,}",
                B_saved_share=d["B_saved"] / d["B_total"],
                B_hurt=wil_str(d["B_hurt"]),
                B_med=d["B_four"]["median"],
                B_p95=d["B_four"]["p95"],
                B_gt01=wil_str(d["B_four"]["gt01"]),
                B_wrong=wil_str(d["B_four"]["wrong"]),
            )
        )
    S["stop_sweeps"] = sweeps
    S["stop_B_refit_all"] = four(D.eR[~D.isA])
    T["ar3c_stop"] = DT.markdown_table(
        prow,
        [
            "rule",
            "threshold",
            "A_saved_share",
            "B_saved",
            "B_saved_share",
            "B_hurt",
            "B_med",
            "B_p95",
            "B_gt01",
            "B_wrong",
        ],
        {"B_saved_share": ".3f", "B_med": ".4f", "B_p95": ".4f"},
    )
    lrows = []
    for kind in ("first_fail", "cost"):
        d = S["stop"].get(kind, {})
        if not d.get("pick"):
            continue
        for L in LEVELS:
            b = d["by_level"][str(L)]
            lrows.append(dict(rule=kind, level=lvl_name(L), **b))
    T["ar3c_stop_level"] = DT.markdown_table(
        lrows, ["rule", "level", "n", "saved", "total", "hurt"], {"saved": ",", "total": ","}
    )


def second_fit(D: Data, S: Dict[str, Any], T: Dict[str, str]) -> None:
    r1 = {r["voxel"]: r for r in D.rR}
    v2 = np.array([r["voxel"] for r in D.rR2], dtype=int)
    if not len(v2):
        return
    e1, e2 = D.eR[v2], D.eR2[v2]
    moved_deg = D.T.mis(D.Rm[v2], D.R2m[v2])
    moved = moved_deg > G.HURT
    helped = e2 < e1 - G.HURT
    hurt = e2 > e1 + G.HURT
    ratio = np.array([r1[int(v)]["end_hit"] / r1[int(v)]["entry_best_hit"] for v in v2])
    lev = D.lev[v2]
    isA = D.isA[v2]
    n = len(v2)
    S["second"] = dict(
        n=n,
        moved=wil(int(moved.sum()), n),
        helped=wil(int(helped.sum()), n),
        hurt=wil(int(hurt.sum()), n),
        evals=int(sum(r["evals"] for r in D.rR2)),
        by_level={
            str(L): dict(
                n=int((lev == L).sum()),
                moved=int((moved & (lev == L)).sum()),
                helped=int((helped & (lev == L)).sum()),
                hurt=int((hurt & (lev == L)).sum()),
            )
            for L in LEVELS
        },
    )
    rows = []
    for L in LEVELS:
        d = S["second"]["by_level"][str(L)]
        rows.append(dict(level=lvl_name(L), **d))
    d = S["second"]
    rows.append(
        dict(
            level="all", n=n, moved=int(moved.sum()), helped=int(helped.sum()), hurt=int(hurt.sum())
        )
    )
    T["ar3c_second"] = DT.markdown_table(rows, ["level", "n", "moved", "helped", "hurt"])
    T["ar3c_second_wil"] = DT.markdown_table(
        [
            dict(what="moved by more than 0.05 deg", v=wil_str(d["moved"])),
            dict(what="brought closer by more than 0.05 deg", v=wil_str(d["helped"])),
            dict(what="pushed away by more than 0.05 deg", v=wil_str(d["hurt"])),
        ],
        ["what", "v"],
    )
    # threshold on the first fit's end hit ratio over the entry's best
    sweep = []
    for t in HIT_T:
        pred = ratio >= t
        a = pred & isA
        sweep.append(dict(t=t, A_pred=int(a.sum()), A_helped=int((a & helped).sum())))
    pick = None
    for r in sweep:  # the lowest t with at most 1% of A's predicted voxels helped
        if r["A_pred"] > 0 and r["A_helped"] / r["A_pred"] <= 0.01:
            pick = r
            break
    S["second"]["sweep"] = sweep
    S["second"]["pick"] = pick
    if pick is not None:
        pred = (ratio >= pick["t"]) & ~isA
        k = int(pred.sum())
        right = int((pred & ~helped).sum())
        S["second"]["B"] = dict(
            predicted=k,
            right=wil(right, k),
            right_share=right / k if k else float("nan"),
            B_total=int((~isA).sum()),
            B_helped=wil(int((helped & ~isA).sum()), int((~isA).sum())),
            predicted_helped=int((pred & helped).sum()),
        )
        d2 = S["second"]["B"]
        T["ar3c_second_rule"] = DT.markdown_table(
            [
                dict(
                    threshold=pick["t"],
                    A_pred=pick["A_pred"],
                    A_helped=pick["A_helped"],
                    B_pred=k,
                    B_right=wil_str(d2["right"]),
                    B_base=wil_str(d2["B_helped"]),
                )
            ],
            ["threshold", "A_pred", "A_helped", "B_pred", "B_right", "B_base"],
            {"threshold": ".2f"},
        )


def point6(D: Data, S: Dict[str, Any], T: Dict[str, str]) -> None:
    lg = D.Dlog
    vox = np.array([r["voxel"] for r in lg], dtype=int)
    ang = np.array([r["angle_deg"] for r in lg])
    ev = np.array([r["evals"] for r in lg])
    same = np.array([r["exact_same"] for r in lg])
    end = np.array([r["end"] for r in lg])
    err_end = D.T.err(end, vox)
    err_start = D.T.err(np.array([r["start"] for r in lg]), vox)
    # `exact_same` is bitwise (np.array_equal); "not lowered" is the fit's own test (end cost not
    # below the start cost), which also counts ends equal to the start up to the quaternion round
    # trip
    sc_d = np.array([r["start_cost"] for r in lg])
    ec_d = np.array([r["end_cost"] for r in lg])
    not_lowered = ec_d >= sc_d
    rt = not_lowered & ~same
    S["D_fit_inert"] = dict(
        fits=len(lg),
        exact_same=int(same.sum()),
        not_lowered=wil(int(not_lowered.sum()), len(lg)),
        not_lowered_share=float(not_lowered.mean()),
        roundtrip_only=int(rt.sum()),
        roundtrip_moved_median_deg=float(np.median(ang[rt])) if rt.any() else float("nan"),
        roundtrip_moved_max_deg=float(np.max(ang[rt])) if rt.any() else float("nan"),
    )
    # the fit's success by the start's error, and whether the truth's strict cost was lower than
    # the start's (Reading 0; truth used in this analysis only)
    ct = D.cost_true[vox]
    phase = np.array([r["phase"] for r in lg])
    edges = [0.0, 0.05, 0.1, 0.2, 0.5, 1.0, 5.0, 1e9]
    drow = []
    S["D_fit_by_start_error"] = {}
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (err_start >= lo) & (err_start < hi)
        nm = int(m.sum())
        if not nm:
            continue
        name = bin_name(lo, hi)
        lowered = int((ec_d[m] < sc_d[m]).sum())
        tl = int((ct[m] < sc_d[m]).sum())
        S["D_fit_by_start_error"][name] = dict(
            fits=nm,
            revisit_fits=int((phase[m] == "revisit").sum()),
            lowered=wil(lowered, nm),
            lowered_share=lowered / nm,
            truth_lower=wil(tl, nm),
            truth_lower_share=tl / nm,
        )
        drow.append(
            dict(
                start_error=name,
                fits=nm,
                lowered=wil_str(wil(lowered, nm)),
                truth_lower=wil_str(wil(tl, nm)),
            )
        )
    T["ar3c_dfit_bins"] = DT.markdown_table(
        drow, ["start_error", "fits", "lowered", "truth_lower"], {"fits": ","}
    )
    rows = []
    S["D_fits"] = {}
    for L in LEVELS + (None,):
        m = np.ones(len(lg), bool) if L is None else (D.lev[vox] == L)
        key = "all" if L is None else str(L)
        nm = int(m.sum())
        S["D_fits"][key] = dict(
            fits=nm,
            evals=int(ev[m].sum()),
            median_evals=float(np.median(ev[m])),
            moved_median=float(np.median(ang[m])),
            moved_p95=float(np.percentile(ang[m], 95)),
            exact_same=wil(int(same[m].sum()), nm),
            not_lowered=wil(int(not_lowered[m].sum()), nm),
            end_four=four(err_end[m]),
            start_four=four(err_start[m]),
        )
        d = S["D_fits"][key]
        rows.append(
            dict(
                level="all" if L is None else lvl_name(L),
                fits=nm,
                evals=d["evals"],
                med_evals=d["median_evals"],
                mv_med=d["moved_median"],
                mv_p95=d["moved_p95"],
                same=wil_str(d["exact_same"]),
                not_lowered=wil_str(d["not_lowered"]),
                s_med=d["start_four"]["median"],
                e_med=d["end_four"]["median"],
                e_p95=d["end_four"]["p95"],
                e_wrong=wil_str(d["end_four"]["wrong"]),
            )
        )
    T["ar3c_dfits"] = DT.markdown_table(
        rows,
        [
            "level",
            "fits",
            "evals",
            "med_evals",
            "mv_med",
            "mv_p95",
            "same",
            "not_lowered",
            "s_med",
            "e_med",
            "e_p95",
            "e_wrong",
        ],
        {
            "fits": ",",
            "evals": ",",
            "med_evals": ".0f",
            "mv_med": ".4f",
            "mv_p95": ".4f",
            "s_med": ".4f",
            "e_med": ".4f",
            "e_p95": ".4f",
        },
    )
    # maps: P plus R against D (all voxels), by level
    prow = []
    S["map_vs_D"] = {}
    for L in LEVELS + (None,):
        m = np.ones(D.n, bool) if L is None else (D.lev == L)
        key = "all" if L is None else str(L)
        wR, wD = D.eR[m] > 1, D.eD[m] > 1
        S["map_vs_D"][key] = dict(
            n=int(m.sum()),
            PR=four(D.eR[m]),
            D=four(D.eD[m]),
            P=four(D.eP[m]),
            fixed=int((wD & ~wR).sum()),
            broken=int((~wD & wR).sum()),
            both=int((wD & wR).sum()),
            intervals_overlap=overlap(four(D.eR[m])["wrong"], four(D.eD[m])["wrong"]),
            gt01_overlap=overlap(four(D.eR[m])["gt01"], four(D.eD[m])["gt01"]),
        )
        d = S["map_vs_D"][key]
        prow.append(
            dict(
                level="all" if L is None else lvl_name(L),
                n=d["n"],
                P_wrong=wil_str(d["P"]["wrong"]),
                PR_wrong=wil_str(d["PR"]["wrong"]),
                D_wrong=wil_str(d["D"]["wrong"]),
                fixed=d["fixed"],
                broken=d["broken"],
                both=d["both"],
                iv="overlap" if d["intervals_overlap"] else "do not overlap",
                iv_gt01="overlap" if d["gt01_overlap"] else "do not overlap",
            )
        )
    T["ar3c_map_vs_D"] = DT.markdown_table(
        prow,
        [
            "level",
            "n",
            "P_wrong",
            "PR_wrong",
            "D_wrong",
            "fixed",
            "broken",
            "both",
            "iv",
            "iv_gt01",
        ],
        {"n": ","},
    )
    # gains under the plan's rule (the two rates' intervals do not overlap), P plus R against D
    S["gains_PR_vs_D"] = [
        dict(
            group=k,
            measure=meas,
            PR=S["map_vs_D"][k]["PR"][meas],
            D=S["map_vs_D"][k]["D"][meas],
            PR_better=S["map_vs_D"][k]["PR"][meas]["k"] < S["map_vs_D"][k]["D"][meas]["k"],
        )
        for k in S["map_vs_D"]
        for meas, flag in (("gt01", "gt01_overlap"), ("wrong", "intervals_overlap"))
        if not S["map_vs_D"][k][flag]
    ]
    frows = [
        four_row("P plus R, all", D.eR),
        four_row("D, all", D.eD),
        four_row("P (start), all", D.eP),
    ]
    for L in LEVELS:
        frows.append(four_row(f"P plus R, {lvl_name(L)}", D.eR[D.lev == L]))
        frows.append(four_row(f"D, {lvl_name(L)}", D.eD[D.lev == L]))
    T["ar3c_map_four"] = DT.markdown_table(frows, FOUR_COLS, FOUR_FMT)


def control(D: Data, S: Dict[str, Any], T: Dict[str, str]) -> None:
    m = D.lev == 0.0
    out = {}
    rows = []
    maps = {}
    for k in ("C", "D"):
        mapk, _ = G.read_map(L2 / k / "recon.mic", D.n)
        maps[k] = D.T.err(mapk)
    sets = [
        ("this run P (paint)", D.eP),
        ("this run P plus R", D.eR),
        ("this run D (BFS)", D.eD),
        ("step 2 C (paint)", maps["C"]),
        ("step 2 D (BFS)", maps["D"]),
    ]
    for name, e in sets:
        out[name] = four(e[m])
        rows.append(four_row(name, e[m]))
    S["control"] = out
    T["ar3c_control"] = DT.markdown_table(rows, FOUR_COLS, FOUR_FMT)


def walls(D: Data, S: Dict[str, Any], T: Dict[str, str]) -> None:
    rp, rd = D.runP, D.runD
    rows = [
        dict(
            run="P (paint, with set-up of pairing)",
            wall=rp["wall_total_s"],
            label=rp["timing_label"],
        ),
        # R and R2 ran in a separate refine-stage process (`--stage refine`, P's stored paint and
        # generator state) with no preflight, about 11:24 to 11:42; run D ran until 11:39 (the
        # mtimes of refine_R.json.gz 11:38, refine_R2.json.gz 11:42, D's fit_log.json.gz 11:39)
        dict(
            run="R (one fit per painted voxel)",
            wall=rp["wall_R_s"],
            label="contended: separate refine-stage process, no preflight; run D running",
        ),
        dict(
            run="R2 (second fit, every fourth voxel)",
            wall=rp["wall_R2_s"],
            label="contended: separate refine-stage process, no preflight; run D running at start",
        ),
        dict(run="D (BFS)", wall=rd["wall_total_s"], label=rd["timing_label"]),
    ]
    S["walls"] = {r["run"]: r for r in rows}
    tot_s = rp["wall_total_s"] + rp["wall_R_s"] + rp["wall_R2_s"] + rd["wall_total_s"]
    S["compute"] = dict(
        P_h=rp["wall_total_s"] / 3600,
        D_h=rd["wall_total_s"] / 3600,
        R_plus_R2_h=(rp["wall_R_s"] + rp["wall_R2_s"]) / 3600,
        total_h=tot_s / 3600,
        budget_h=4.0,
        P_cap_min=90,
        D_cap_min=150,
        caps_coded=False,
    )
    T["ar3c_wall"] = DT.markdown_table(rows, ["run", "wall", "label"], {"wall": ".0f"})


def figure(D: Data, S: Dict[str, Any]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 3, figsize=(16, 9))
    x = np.arange(len(LEVELS))
    labs = [f"{L:g}" for L in LEVELS]
    # (a) the start, the refinement and D by level (95th percentile and median)
    for e, name, c in ((D.eP, "start (P)", "C0"), (D.eR, "after R", "C2"), (D.eD, "D (BFS)", "C3")):
        ax[0, 0].plot(
            x, [np.median(e[D.lev == L]) for L in LEVELS], "o-", c=c, label=name + " median"
        )
        ax[0, 0].plot(
            x, [np.percentile(e[D.lev == L], 95) for L in LEVELS], "s--", c=c, label=name + " 95th"
        )
    ax[0, 0].set_yscale("log")
    ax[0, 0].set_xticks(x, labs)
    ax[0, 0].set_xlabel("level of the grain (degrees)")
    ax[0, 0].set_ylabel("error (degrees)")
    ax[0, 0].set_title("error by level")
    ax[0, 0].legend(fontsize=7)
    # (b) reading 0
    r0 = S["reading0"]
    ax[0, 1].bar(x, [r0[str(L)]["below_09"]["k"] / r0[str(L)]["below_09"]["n"] for L in LEVELS])
    ax[0, 1].set_xticks(x, labs)
    ax[0, 1].set_xlabel("level of the grain (degrees)")
    ax[0, 1].set_ylabel("share of voxels")
    ax[0, 1].set_title("centre orientation hits < 0.9 x the true one")
    # (c) skip rules: A's curves
    for kind, c in (("hit", "C0"), ("cost", "C1")):
        sw = S["skip_sweeps"][kind]
        tot = max(1, sum(r["A_skipped_evals"] for r in sw[-1:]))
        ax[0, 2].plot(
            [100 * r["A_hurt_share"] for r in sw],
            [r["A_skipped_evals"] / 1e6 for r in sw],
            "o-",
            c=c,
            label={"hit": "hit ratio rule", "cost": "start cost rule"}[kind],
        )
    ax[0, 2].axvline(1.0, c="k", ls=":")
    ax[0, 2].set_xlabel("share of pick-set voxels hurt (%)")
    ax[0, 2].set_ylabel("evaluations skipped (millions)")
    ax[0, 2].set_title("skip the fit (pick set A)")
    ax[0, 2].legend(fontsize=7)
    # (d) after-convergence share
    cv = S["converged"]
    ax[1, 0].bar(x, [cv[str(L)]["share"] for L in LEVELS])
    ax[1, 0].set_xticks(x, labs)
    ax[1, 0].set_xlabel("level of the grain (degrees)")
    ax[1, 0].set_ylabel("share of the fit's evaluations")
    ax[1, 0].set_title("evaluations after within 0.02 deg of the end")
    # (e) stop-in-fit rules on A
    for kind, c in (("first_fail", "C0"), ("cost", "C1")):
        sw = S["stop_sweeps"][kind]
        ax[1, 1].plot(
            [100 * r["A_hurt_share"] for r in sw],
            [r["A_saved"] / 1e6 for r in sw],
            "o-",
            c=c,
            label=kind,
        )
    ax[1, 1].axvline(1.0, c="k", ls=":")
    ax[1, 1].set_xlabel("share of pick-set voxels hurt (%)")
    ax[1, 1].set_ylabel("evaluations saved (millions)")
    ax[1, 1].set_title("stop inside the fit (pick set A)")
    ax[1, 1].legend(fontsize=7)
    # (f) D vs P plus R wrong shares by level
    w = 0.25
    for i, (k, nm) in enumerate((("P", "P"), ("PR", "P plus R"), ("D", "D"))):
        ax[1, 2].bar(
            x + (i - 1) * w,
            [
                S["map_vs_D"][str(L)][k]["wrong"]["k"] / S["map_vs_D"][str(L)][k]["wrong"]["n"]
                for L in LEVELS
            ],
            w,
            label=nm,
        )
    ax[1, 2].set_xticks(x, labs)
    ax[1, 2].set_xlabel("level of the grain (degrees)")
    ax[1, 2].set_ylabel("share wrong (above 1 degree)")
    ax[1, 2].set_title("wrong voxels by level")
    ax[1, 2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(G.GRAD / "ar3_gradient.png", dpi=110)


def main() -> None:
    D = Data()
    S: Dict[str, Any] = dict(
        n_voxels=D.n,
        truth_sha256=D.T.sha,
        n_by_level={str(L): int((D.lev == L).sum()) for L in LEVELS},
    )
    T: Dict[str, str] = {}
    sn = D.stack_note
    S["provenance"] = dict(
        render_grad_lit_clean=int(D.render_grad["lit_clean_total"]),
        render_grad_spots=int(D.render_grad["spots_total"]),
        render_stored_lit_clean=int(D.render_full["lit_clean_total"]),
        render_stored_spots=int(D.render_full["spots_total"]),
        stack_nonzero=int(sn["stack_nonzero"]),
        stack_matches_render=bool(sn["stack_matches_render"]),
        render_mic_sha256=sn["mic_sha256"],
        trace_check=D.trace_check,
    )
    assert sn["mic_sha256"] == D.T.sha
    point0(D, S, T)
    point1(D, S, T)
    paired_vs_start(D, S, T)
    fit_cause(D, S, T)
    step_sizes(D, S, T)
    pick_skip(D, S, T)
    block_stop(D, S, T)
    second_fit(D, S, T)
    point6(D, S, T)
    control(D, S, T)
    walls(D, S, T)
    figure(D, S)
    (G.GRAD / "summary.json").write_text(json.dumps(S, indent=1, default=float) + "\n")
    DT.write_tables(G.GRAD / "tables.md", T)
    print("wrote summary.json, tables.md, ar3_gradient.png")


if __name__ == "__main__":
    main()
