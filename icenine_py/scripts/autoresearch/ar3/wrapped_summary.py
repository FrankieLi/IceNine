"""AR3 step 4, part 3: summary of the whole-layer runs D4, P4 (with R4) and E4.

  uv run python scripts/autoresearch/ar3/wrapped_summary.py

Reads ``fine/runs/{D4,P4,E4}`` and the compared runs of steps 2 and 3 (``gradient/runs/{D,P}``,
``layer/runs/D``). Nothing is rerun. Writes the Part 3 entries into ``fine/summary.json`` (key
``part3``) and ``fine/tables.md`` (run ``fine_summary.py`` first; this script adds to its files).
A run that has no ``run.json`` is reported as not finished, with its last ``progress.json``. Wrong =
more than 1 degree from the truth (symmetry-reduced). Counts by voxel, Wilson 95% intervals, no
p-values; a difference is a gain or a loss only when intervals do not overlap. One run per arm,
one voxel order.
"""

import gzip
import json
import shutil
import sys
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import ar3_common as C
import gradient_common as G
import numpy as np
import run_gradient as RG
from fine_summary import FINE, FRUNS, LEVELS, lvl_name, overlap, wil, wil_str

import doc_tables as DT  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0] / "ar1"))
sys.path.insert(0, str(HERE.parents[0] / "ar2"))

L2 = C.OUT / "layer" / "runs"
WRONG = G.WRONG


def four(err: np.ndarray) -> Dict[str, Any]:
    n = len(err)
    return dict(
        n=n,
        median=float(np.median(err)),
        p95=float(np.percentile(err, 95)),
        gt01=wil(int((err > 0.1).sum()), n),
        wrong=wil(int((err > WRONG).sum()), n),
    )


def ungz(path: Path) -> Path:
    """A plain copy of a gzipped file in a temp folder (the caller's ``read_map`` wants a path)."""
    td = Path(tempfile.mkdtemp())
    out = td / "recon.mic"
    with gzip.open(path, "rb") as fi, open(out, "wb") as fo:
        shutil.copyfileobj(fi, fo)
    return td


def mic_path(run: Path) -> Path:
    return run / "recon.mic" if (run / "recon.mic").exists() else run / "recon.mic.gz"


def read_gradient_map(run: Path, n: int) -> np.ndarray:
    R, _ = G.read_map(mic_path(run), n)
    return R


def paired(ea: np.ndarray, eb: np.ndarray, sel: np.ndarray) -> Dict[str, Any]:
    """Baseline ``ea`` (the compared run), other ``eb`` (this step's run) over ``sel``."""
    wa, wb = (ea > WRONG)[sel], (eb > WRONG)[sel]
    return dict(
        n=int(sel.sum()),
        wrong_base=int(wa.sum()),
        wrong_other=int(wb.sum()),
        fixed=int((wa & ~wb).sum()),
        broken=int((~wa & wb).sum()),
        both_wrong=int((wa & wb).sum()),
    )


def finished(kind: str) -> bool:
    return (FRUNS / kind / "run.json").exists()


def status_row(kind: str) -> Dict[str, Any]:
    out = FRUNS / kind
    killed = (FRUNS.parent / "cache" / "logs" / f"{kind}.killed").exists()
    prog = out / "progress.json"
    d: Dict[str, Any] = dict(kind=kind, finished=finished(kind), killed_at_cap=killed)
    if prog.exists():
        d["last_progress"] = json.loads(prog.read_text())
    return d


def err_rows(
    S: Dict[str, Any], key: str, name: str, err: np.ndarray, lev: Optional[np.ndarray]
) -> List[Dict[str, Any]]:
    rows = []
    groups: List[Tuple[str, Optional[float]]] = [("all", None)]
    if lev is not None:
        groups += [(lvl_name(L), L) for L in LEVELS]
    S[key] = {}
    for g, L in groups:
        m = np.ones(len(err), bool) if L is None else (lev == L)
        f = four(err[m])
        S[key][g] = f
        rows.append(
            dict(
                run=name,
                group=g,
                n=f["n"],
                median=f["median"],
                p95=f["p95"],
                gt01=wil_str(f["gt01"]),
                wrong=wil_str(f["wrong"]),
            )
        )
    return rows


SIDE_BY_SIDE = "; D4, P4 and E4 ran side by side (three heavy processes on the machine)"
RATIO_BINS: List[Tuple[str, float, float]] = [
    ("below 0.90", -1.0, 0.90),
    ("0.90 to 0.92", 0.90, 0.92),
    ("0.92 to 0.95", 0.92, 0.95),
    ("0.95 to 0.98", 0.95, 0.98),
    ("0.98 to 1.00", 0.98, 1.0),
]


def part3_label(stored: str) -> str:
    """The stored ``run.json`` label says "single run, one worker"; the three part 3 runs ran side
    by side, so the summary writes the label with that corrected (the stored one is kept too)."""
    return stored.replace("; single run, one worker", "") + SIDE_BY_SIDE


def ratio_rows(
    name: str, ratio: np.ndarray, wrong: np.ndarray, out: Dict[str, Any]
) -> List[Dict[str, Any]]:
    """Voxels and wrong voxels by acceptance ratio bin (lower edge included, upper excluded; a ratio
    of exactly 1 has its own row). A description, not a rule: no threshold is picked."""
    rows = []
    for lab, lo, hi in RATIO_BINS + [("1.00 exactly", 1.0, 1.0)]:
        m = (ratio == 1.0) if lo == hi else ((ratio >= lo) & (ratio < hi))
        if lab == "below 0.90" and not m.any():
            continue
        w = wil(int(wrong[m].sum()), int(m.sum()))
        out[lab] = w
        rows.append(dict(run=name, ratio=lab, n=int(m.sum()), wrong=wil_str(w)))
    return rows


def bfs_source(rec: Dict[str, Any], v: int) -> str:
    return str(rec[str(v)]["source"])


def stage_evals_bfs(stats: Dict[str, Any]) -> Dict[str, int]:
    keys = (
        "n_evals_seed",
        "n_evals_neighbor",
        "n_evals_revisit",
        "n_evals_reseed",
        "n_evals_restart",
    )
    d = {k.replace("n_evals_", ""): int(stats.get(k, 0)) for k in keys}
    d["total"] = int(sum(d.values()))
    return d


def figure3(P3: Dict[str, Any]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pairs = [
        ("D4 vs D (step 3)", "D_error", "D4_error", "D_run", "D4_run"),
        ("P4+R4 vs P+R (step 3)", "PR_error", "P4R4_error", "P_run", "P4_run"),
        ("E4 vs D (step 2)", "E_error", "E4_error", "E_run", "E4_run"),
    ]
    pairs = [p for p in pairs if p[2] in P3 and p[1] in P3]
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    for i, (_name, ka, kb, ra, rb) in enumerate(pairs):
        for j, (k, r, col) in enumerate(((ka, ra, "gray"), (kb, rb, "tab:red"))):
            x = i + (j - 0.5) * 0.35
            ev = P3[r]["evals_total"]
            if "evals_R" in P3[r]:
                ev += P3[r]["evals_R"]
            if "evals_r4" in P3[r]:
                ev += P3[r]["evals_r4"]
            ax[0].bar(x, ev / 1e6, 0.33, color=col)
            for a_, key in ((1, "gt01"), (2, "wrong")):
                d = P3[k]["all"][key]
                m = d["k"] / d["n"]
                ax[a_].errorbar(
                    x,
                    100 * m,
                    yerr=[[100 * (m - d["lo"])], [100 * (d["hi"] - m)]],
                    fmt="o",
                    color=col,
                )
    for a_ in range(3):
        ax[a_].set_xticks(range(len(pairs)))
        ax[a_].set_xticklabels([p[0] for p in pairs], fontsize=8)
    ax[0].set_ylabel("evaluations (millions)")
    ax[0].set_title("Evaluations (gray: compared run, red: wrapped fit)")
    ax[1].set_ylabel("voxels above 0.1 deg (%), Wilson 95%")
    ax[1].set_title("Share above 0.1 degrees")
    ax[2].set_ylabel("wrong voxels, above 1 deg (%), Wilson 95%")
    ax[2].set_title("Share wrong")
    fig.tight_layout()
    fig.savefig(FINE / "ar3_fine_layer.png", dpi=110)


def main() -> None:
    S = json.loads((FINE / "summary.json").read_text())
    T = DT.read_tables(FINE / "tables.md")
    P3: Dict[str, Any] = {}
    P3["status"] = {k: status_row(k) for k in ("D4", "P4", "E4")}
    T_ = G.Truth()
    n = T_.n
    lev = T_.level
    err_rows_all: List[Dict[str, Any]] = []
    ev_rows: List[Dict[str, Any]] = []
    pair_rows: List[Dict[str, Any]] = []
    wr_rows: List[Dict[str, Any]] = []
    ent_rows: List[Dict[str, Any]] = []
    brk_rows: List[Dict[str, Any]] = []  # BFS runs: where the broken voxels were settled
    ratio_rows_all: List[Dict[str, Any]] = []  # acceptance ratio bins (P4, E4)
    brk_paint_rows: List[Dict[str, Any]] = []  # P4 plus R4: where the broken voxels were settled
    brk_paint_cols: List[str] = []

    def wrapper_row(kind: str, w: Dict[str, Any]) -> Dict[str, Any]:
        return dict(
            run=kind,
            calls=w["wrapper_calls"],
            callers=", ".join(f"{k} {v}" for k, v in w["wrapper_calls_by_caller"].items()),
            fell=w["picked_fell"],
            fallback=w["fallback_runs"],
            fallback_won=w["fallback_won"],
            cma_calls=w["cma_calls_total"],
            vmo_in=w["vmo_calls_in_wrapper"],
            vmo_out=w["vmo_calls_outside_wrapper"],
        )

    # ------------------------------------------------------------------ D4 against D
    if finished("D4"):
        run = json.loads((FRUNS / "D4" / "run.json").read_text())
        D = json.loads((G.GRAD / "runs" / "D" / "run.json").read_text())
        eD = T_.err(read_gradient_map(G.GRAD / "runs" / "D", n))
        e4 = T_.err(read_gradient_map(FRUNS / "D4", n))
        err_rows_all += err_rows(P3, "D4_error", "D4", e4, lev)
        err_rows_all += err_rows(P3, "D_error", "D (step 3)", eD, lev)
        P3["D4_run"] = dict(
            evals_total=run["evals_total"],
            wall_total_s=run["wall_total_s"],
            timing_label=part3_label(run["timing_label"]),
            timing_label_run_json=run["timing_label"],
            stages=stage_evals_bfs(run["stats"]),
            n_processed=run["n_processed"],
            wrapper=run["wrapper"],
        )
        P3["D_run"] = dict(
            evals_total=D["evals_total"],
            wall_total_s=D["wall_total_s"],
            stages=stage_evals_bfs(D["stats"]),
        )
        for nm, rr in (("D4", P3["D4_run"]), ("D (step 3)", P3["D_run"])):
            ev_rows.append(dict(run=nm, **rr["stages"]))
        pr = {"all": paired(eD, e4, np.ones(n, bool))}
        for L in LEVELS:
            pr[lvl_name(L)] = paired(eD, e4, lev == L)
        P3["D4_vs_D"] = pr
        for g, d in pr.items():
            pair_rows.append(dict(pair="D4 vs D", group=g, **d))
        wr_rows.append(wrapper_row("D4", run["wrapper"]))
        for nm, st in (("D4", run["stats"]), ("D (step 3)", D["stats"])):
            ent_rows.append(
                dict(
                    run=nm,
                    entries=int(st["n_seeds"]) + int(st.get("n_reseeds", 0)),
                    seeds=int(st["n_seeds"]),
                    reseeds=int(st.get("n_reseeds", 0)),
                    index_used=int(st.get("index_used", 0)),
                    grid_fallbacks=int(st.get("index_fallback_seed_test", 0)),
                    grid_evals=int(st.get("index_fallback_grid_evals", 0)),
                )
            )
        # where D4's saving comes from: the fits (neighbour plus revisit) and the reseeds
        s4, sD = P3["D4_run"]["stages"], P3["D_run"]["stages"]
        P3["D4_saving"] = dict(
            fits_D4=s4["neighbor"] + s4["revisit"],
            fits_D=sD["neighbor"] + sD["revisit"],
            reseed_D4=s4["reseed"],
            reseed_D=sD["reseed"],
        )
        # where the broken voxels (right in D, wrong in D4) were settled; D's records are not
        # stored (step 3 kept only its fit log), so only D4's side can be given
        rec4 = G.load_gz(FRUNS / "D4" / "records.json.gz")
        brk = np.where((eD <= WRONG) & (e4 > WRONG))[0]
        src4 = Counter(bfs_source(rec4, int(v)) for v in brk)
        unres = [int(k) for k, r in rec4.items() if r["source"] == "unresolved"]
        P3["D4_broken_source"] = dict(
            n_broken=int(len(brk)),
            D4=dict(src4),
            broken_with_a_rejected_fit=int(sum(bool(rec4[str(v)]["local_rejected"]) for v in brk)),
            unresolved_all=len(unres),
            unresolved_broken=int(sum(int(v) in set(unres) for v in brk)),
        )
        for s_ in sorted(src4, key=lambda k: -src4[k]):
            brk_rows.append(
                dict(pair="D4 vs D (step 3)", source=s_, this_run=src4[s_], compared="not stored")
            )
        a, b = P3["D4_error"]["all"]["gt01"], P3["D_error"]["all"]["gt01"]
        P3["D4_vs_D_gt01_overlap"] = bool(overlap(a, b))
        a, b = P3["D4_error"]["all"]["wrong"], P3["D_error"]["all"]["wrong"]
        P3["D4_vs_D_wrong_overlap"] = bool(overlap(a, b))

    # ------------------------------------------------------------------ P4 + R4 against P + R
    if finished("P4"):
        run = json.loads((FRUNS / "P4" / "run.json").read_text())
        P = json.loads((G.GRAD / "runs" / "P" / "run.json").read_text())
        base, _ = G.read_map(mic_path(FRUNS / "P4"), n)
        rows4 = G.load_gz(FRUNS / "P4" / "refine_R4.json.gz")
        Rmap = base.copy()
        for r in rows4:
            Rmap[r["voxel"]] = np.array(r["end"])
        Rbase = np.load(G.GRAD / "runs" / "P" / "maps.npz")["R"].astype(np.float64)
        e4 = T_.err(Rmap)
        eR = T_.err(Rbase)
        e4p = T_.err(base)
        err_rows_all += err_rows(P3, "P4R4_error", "P4 plus R4", e4, lev)
        err_rows_all += err_rows(P3, "PR_error", "P plus R (step 3)", eR, lev)
        err_rows(P3, "P4_error", "P4 (paint only)", e4p, lev)
        P3["P4_run"] = {
            k: run[k] for k in ("evals", "evals_total", "evals_r4", "counts", "n_entries")
        }
        P3["P4_run"].update(
            n_refined_r4=run["n_refined_r4"],
            wall_p_s=run["wall_p_s"],
            wall_r4_s=run["wall_r4_s"],
            timing_label=part3_label(run["timing_label"]),
            timing_label_run_json=run["timing_label"],
            wrapper_after_paint=run["wrapper_after_paint"],
            wrapper_final=run["wrapper_final"],
        )
        P3["P_run"] = {k: P[k] for k in ("evals", "evals_total", "evals_R", "counts", "n_entries")}
        P3["P_run"]["wall_total_s"] = P["wall_total_s"]
        for nm, ev, extra in (
            ("P4 plus R4", run["evals"], run["evals_r4"]),
            ("P plus R (step 3)", P["evals"], P["evals_R"]),
        ):
            row = dict(run=nm)
            row.update({k: int(v) for k, v in ev.items()})
            row["refine"] = int(extra)
            row["total"] = int(sum(ev.values()) + extra)
            ev_rows.append(row)
        pr = {"all": paired(eR, e4, np.ones(n, bool))}
        for L in LEVELS:
            pr[lvl_name(L)] = paired(eR, e4, lev == L)
        P3["P4R4_vs_PR"] = pr
        for g, d in pr.items():
            pair_rows.append(dict(pair="P4+R4 vs P+R", group=g, **d))
        wr_rows.append(wrapper_row("P4 (paint)", run["wrapper_after_paint"]))
        for nm, c, e in (("P4", run["counts"], None), ("P (step 3)", P["counts"], None)):
            ent_rows.append(
                dict(
                    run=nm,
                    entries=int(c["seeds"]),
                    seeds=int(c["seeds"]),
                    reseeds=0,
                    index_used=0,
                    grid_fallbacks=int(c.get("failed", 0)),
                    grid_evals=0,
                )
            )
        # where the wrong voxels sit: by the stage that settled the voxel in the paint. A failed
        # voxel that a later grain painted got a start from that grain and was refined (R4, R).
        st_rows = []
        P3["P4R4_wrong_by_stage"] = {}
        kinds: Dict[str, np.ndarray] = {}
        paints: Dict[str, Any] = {}
        stages = (
            "seed",
            "fit-settled",
            "painted from a grain",
            "failed, later painted",
            "failed, not painted",
        )
        for nm, rd, e_ in (
            ("P4 plus R4", FRUNS / "P4", e4),
            ("P plus R (step 3)", G.GRAD / "runs" / "P", eR),
        ):
            res_ = RG.load_paint(rd)
            paints[nm] = res_
            kind = np.full(n, "painted from a grain", dtype=object)
            for v_, st_ in res_.stage.items():
                kind[v_] = "seed" if st_ == "seed" else "fit-settled"
            for v_ in res_.failed:
                later = res_.painted_grain(v_) is not None
                kind[v_] = "failed, later painted" if later else "failed, not painted"
            kinds[nm] = kind
            for k_ in stages:
                m_ = kind == k_
                w_ = wil(int((e_[m_] > WRONG).sum()), int(m_.sum()))
                P3["P4R4_wrong_by_stage"][nm + " / " + k_] = w_
                st_rows.append(dict(run=nm, stage=k_, n=int(m_.sum()), wrong=wil_str(w_)))
        k4 = kinds["P4 plus R4"]
        n_ref = int(((k4 == "painted from a grain") | (k4 == "failed, later painted")).sum())
        assert n_ref == run["n_refined_r4"] == len(rows4), (n_ref, run["n_refined_r4"])
        P3["P4_refined_r4_check"] = dict(
            painted_from_a_grain=int((k4 == "painted from a grain").sum()),
            failed_later_painted=int((k4 == "failed, later painted").sum()),
            n_refined_r4=n_ref,
        )
        # P4's own records agree with the paint: one "fit" record per fit-settled voxel
        rec_p4 = G.load_gz(FRUNS / "P4" / "records.json.gz")
        n_fit_rec = sum(1 for r in rec_p4.values() if r["source"] == "fit")
        assert n_fit_rec == int((k4 == "fit-settled").sum()), n_fit_rec
        P3["P4_records_fit_count"] = n_fit_rec
        S["_stage_rows"] = st_rows
        # the broken voxels (right in P plus R, wrong in P4 plus R4), by P4's stage, with what
        # step 3's P did on the same voxels. Records only, no new run.
        res4, resP = paints["P4 plus R4"], paints["P plus R (step 3)"]
        brk = np.where((eR <= WRONG) & (e4 > WRONG))[0]
        # ``load_paint`` does not rebuild the loop log, so read it from the stored file
        logP = {
            int(r["voxel"]): r for r in G.load_gz(G.GRAD / "runs" / "P" / "paint.json.gz")["log"]
        }
        kP = kinds["P plus R (step 3)"]

        def p_what(v: int) -> str:
            if kP[v] == "seed":
                r = logP.get(v, {})
                return "seed after a rejected fit" if "fit_from" in r else "seed, no fit tried"
            return str(kP[v])

        p_cols = (
            "seed after a rejected fit",
            "seed, no fit tried",
            "fit-settled",
            "painted from a grain",
            "failed, later painted",
            "failed, not painted",
        )
        b_rows = []
        BR: Dict[str, Any] = dict(n_broken=int(len(brk)))
        BR["wrong_after_paint"] = int((e4p[brk] > WRONG).sum())
        for k_ in stages:
            vs = [int(v) for v in brk if k4[v] == k_]
            if not vs:
                continue
            d_: Dict[str, Any] = dict(n=len(vs), wrong_after_paint=int((e4p[vs] > WRONG).sum()))
            other = ""
            if k_ == "fit-settled":
                fr = [res4.fits[v] for v in vs]
                ents = [res4.entries[f.from_grain] for f in fr]
                oth = [int(T_.gmap[e.voxel] != T_.gmap[v]) for e, v in zip(ents, vs)]
                ratio = np.array([f.hit / max(e.best_hit, f.hit) for f, e in zip(fr, ents)])
                to_entry = T_.mis(np.stack([f.R for f in fr]), np.stack([e.R for e in ents]))
                d_.update(
                    entry_seeded_in_other_grain=int(sum(oth)),
                    ratio_min=float(ratio.min()),
                    passed_0_9=int((ratio > 0.9).sum()),
                    end_to_entry_median_deg=float(np.median(to_entry)),
                    end_to_entry_max_deg=float(np.max(to_entry)),
                )
                other = str(int(sum(oth)))
            inP = Counter(p_what(v) for v in vs)
            d_["in_P"] = dict(inP)
            BR[k_] = d_
            row = dict(stage_P4=k_, broken=len(vs), wrong_after_paint=d_["wrong_after_paint"])
            row["entry_other_grain"] = other
            for c_ in p_cols:
                row["P: " + c_] = int(inP.get(c_, 0))
            b_rows.append(row)
        P3["P4R4_broken"] = BR
        brk_paint_rows.extend(b_rows)
        brk_paint_cols.extend(["stage_P4", "broken", "wrong_after_paint", "entry_other_grain"])
        # columns of step 3's P that no broken voxel reached are left out (all zero)
        brk_paint_cols.extend("P: " + c_ for c_ in p_cols if any(r["P: " + c_] for r in b_rows))
        # acceptance ratio of P4's fit-settled voxels: the fit's hit ratio over the larger of the
        # entry's whole-flood best and its own (exactly the ratio the paint tested)
        vf = sorted(res4.fits)
        rf = np.array(
            [
                res4.fits[v].hit
                / max(res4.entries[res4.fits[v].from_grain].best_hit, res4.fits[v].hit)
                for v in vf
            ]
        )
        P3["P4_fit_ratio"] = dict(n=len(vf), min=float(rf.min()), bins={})
        ratio_rows_all += ratio_rows(
            "P4 fit-settled (ratio tested)", rf, e4[vf] > WRONG, P3["P4_fit_ratio"]["bins"]
        )
        a, b = P3["P4R4_error"]["all"]["gt01"], P3["PR_error"]["all"]["gt01"]
        P3["P4R4_vs_PR_gt01_overlap"] = bool(overlap(a, b))
        a, b = P3["P4R4_error"]["all"]["wrong"], P3["PR_error"]["all"]["wrong"]
        P3["P4R4_vs_PR_wrong_overlap"] = bool(overlap(a, b))

    # ------------------------------------------------------------------ E4 against step 2's D
    if finished("E4"):
        import layer_summary as LS

        W = LS.World()
        run = json.loads((FRUNS / "E4" / "run.json").read_text())
        D2 = json.loads((L2 / "D" / "run.json").read_text())
        Re, have = W.read_map(L2 / "D")
        td = ungz(mic_path(FRUNS / "E4")) if mic_path(FRUNS / "E4").suffix == ".gz" else None
        R4, have4 = W.read_map(td if td is not None else FRUNS / "E4")
        eD2 = W.err(Re)
        eE = W.err(R4)
        err_rows_all += err_rows(P3, "E4_error", "E4", eE, None)
        err_rows_all += err_rows(P3, "E_error", "D (step 2)", eD2, None)
        P3["E4_run"] = dict(
            evals_total=run["evals_total"],
            wall_total_s=run["wall_total_s"],
            timing_label=part3_label(run["timing_label"]),
            timing_label_run_json=run["timing_label"],
            stages=stage_evals_bfs(run["stats"]),
            wrapper=run["wrapper"],
        )
        # where the broken voxels (right in step 2's D, wrong in E4) were settled, in each run
        recE = G.load_gz(FRUNS / "E4" / "records.json.gz")
        recD2 = json.loads((L2 / "D" / "records.json").read_text())
        brk = np.where((eD2 <= WRONG) & (eE > WRONG))[0]
        sE = Counter(bfs_source(recE, int(v)) for v in brk)
        sD2 = Counter(bfs_source(recD2, int(v)) for v in brk)
        P3["E4_broken_source"] = dict(n_broken=int(len(brk)), E4=dict(sE), D_step2=dict(sD2))
        for s_ in sorted(set(sE) | set(sD2), key=lambda k: (-sE.get(k, 0), -sD2.get(k, 0))):
            brk_rows.append(
                dict(
                    pair="E4 vs D (step 2)",
                    source=s_,
                    this_run=int(sE.get(s_, 0)),
                    compared=str(int(sD2.get(s_, 0))),
                )
            )
        # acceptance ratio of E4's accepted fits: the BFS records hold only the fit's hit ratio,
        # not the best of the expansion that the test divided by. The stored hit ratio stands in
        # for it; since that best is at most 1, the stored value is at most the tested ratio.
        vacc = [int(k) for k, r in recE.items() if r["source"] in ("neighbor", "revisit")]
        hE = np.array([recE[str(v)]["hit_ratio"] for v in vacc])
        P3["E4_fit_hit"] = dict(n=len(vacc), min=float(hE.min()), bins={})
        ratio_rows_all += ratio_rows(
            "E4 neighbour or revisit (stored hit ratio)",
            hE,
            eE[vacc] > WRONG,
            P3["E4_fit_hit"]["bins"],
        )
        P3["E_run"] = dict(
            evals_total=D2["evals_total"],
            wall_total_s=D2["wall_total_s"],
            stages=stage_evals_bfs(D2["stats"]),
        )
        for nm, rr in (("E4", P3["E4_run"]), ("D (step 2)", P3["E_run"])):
            ev_rows.append(dict(run=nm, **rr["stages"]))
        d = paired(eD2, eE, np.ones(W.n, bool))
        P3["E4_vs_D2"] = {"all": d}
        pair_rows.append(dict(pair="E4 vs D (step 2)", group="all", **d))
        wr_rows.append(wrapper_row("E4", run["wrapper"]))
        for nm, st in (("E4", run["stats"]), ("D (step 2)", D2["stats"])):
            ent_rows.append(
                dict(
                    run=nm,
                    entries=int(st["n_seeds"]) + int(st.get("n_reseeds", 0)),
                    seeds=int(st["n_seeds"]),
                    reseeds=int(st.get("n_reseeds", 0)),
                    index_used=int(st.get("index_used", 0)),
                    grid_fallbacks=int(st.get("index_fallback_seed_test", 0)),
                    grid_evals=int(st.get("index_fallback_grid_evals", 0)),
                )
            )
        a, b = P3["E4_error"]["all"]["gt01"], P3["E_error"]["all"]["gt01"]
        P3["E4_vs_D2_gt01_overlap"] = bool(overlap(a, b))
        a, b = P3["E4_error"]["all"]["wrong"], P3["E_error"]["all"]["wrong"]
        P3["E4_vs_D2_wrong_overlap"] = bool(overlap(a, b))

    S["part3"] = P3
    figure3(P3)
    if err_rows_all:
        T["ar3d_p3_error"] = DT.markdown_table(
            err_rows_all,
            ["run", "group", "n", "median", "p95", "gt01", "wrong"],
            {"n": ",", "median": ".4f", "p95": ".4f"},
        )
    if ev_rows:
        cols: List[str] = []
        for r in ev_rows:
            for k in r:
                if k not in cols:
                    cols.append(k)
        for r in ev_rows:
            for k in cols:
                r.setdefault(k, "")
        T["ar3d_p3_evals"] = DT.markdown_table(ev_rows, cols, {k: "," for k in cols if k != "run"})
    if pair_rows:
        T["ar3d_p3_paired"] = DT.markdown_table(
            pair_rows,
            ["pair", "group", "n", "wrong_base", "wrong_other", "fixed", "broken", "both_wrong"],
            {"n": ","},
        )
    if wr_rows:
        T["ar3d_p3_wrapper"] = DT.markdown_table(
            wr_rows,
            [
                "run",
                "calls",
                "callers",
                "fell",
                "fallback",
                "fallback_won",
                "cma_calls",
                "vmo_in",
                "vmo_out",
            ],
            {k: "," for k in ("calls", "fell", "fallback", "fallback_won", "cma_calls")}
            | {"vmo_in": ",", "vmo_out": ","},
        )
    if ent_rows:
        T["ar3d_p3_entries"] = DT.markdown_table(
            ent_rows,
            ["run", "entries", "seeds", "reseeds", "index_used", "grid_fallbacks", "grid_evals"],
            {k: "," for k in ("entries", "seeds", "reseeds", "index_used", "grid_evals")},
        )
    if "_stage_rows" in S:
        T["ar3d_p3_stage"] = DT.markdown_table(
            S.pop("_stage_rows"), ["run", "stage", "n", "wrong"], {"n": ","}
        )
    if brk_paint_rows:
        T["ar3d_p3_broken_paint"] = DT.markdown_table(brk_paint_rows, brk_paint_cols)
    if brk_rows:
        T["ar3d_p3_broken_bfs"] = DT.markdown_table(
            brk_rows, ["pair", "source", "this_run", "compared"]
        )
    if ratio_rows_all:
        T["ar3d_p3_ratio"] = DT.markdown_table(
            ratio_rows_all, ["run", "ratio", "n", "wrong"], {"n": ","}
        )
    stat_rows = []
    for k, d in P3["status"].items():
        stat_rows.append(
            dict(
                run=k,
                finished="yes" if d["finished"] else "no",
                killed="yes" if d["killed_at_cap"] else "no",
            )
        )
    T["ar3d_p3_status"] = DT.markdown_table(stat_rows, ["run", "finished", "killed"])
    (FINE / "summary.json").write_text(json.dumps(S, indent=1, default=float) + "\n")
    DT.write_tables(FINE / "tables.md", T)
    print("wrote part 3 into summary.json and tables.md")


if __name__ == "__main__":
    main()
