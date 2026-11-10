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
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import ar3_common as C
import gradient_common as G
import numpy as np
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


def err_rows(S: Dict[str, Any], key: str, name: str, err: np.ndarray, lev: Optional[np.ndarray]):
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
            timing_label=run["timing_label"],
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
        Pbase, _ = G.read_map(G.GRAD / "runs" / "P" / "recon.mic", n)
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
            wall_p_s=run["wall_p_s"],
            wall_r4_s=run["wall_r4_s"],
            timing_label=run["timing_label"],
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
        # where the wrong voxels sit: by the stage that settled the voxel in the paint
        st_rows = []
        P3["P4R4_wrong_by_stage"] = {}
        for nm, rd, e_ in (
            ("P4 plus R4", FRUNS / "P4", e4),
            ("P plus R (step 3)", G.GRAD / "runs" / "P", eR),
        ):
            pj = G.load_gz(rd / "paint.json.gz")
            kind = np.full(n, "painted from a grain", dtype=object)
            for v_, st_ in pj["stage"].items():
                kind[int(v_)] = "seed" if st_ == "seed" else "fit-settled"
            for v_ in pj["failed"]:
                kind[int(v_)] = "failed"
            for k_ in ("seed", "fit-settled", "painted from a grain", "failed"):
                m_ = kind == k_
                w_ = wil(int((e_[m_] > WRONG).sum()), int(m_.sum()))
                P3["P4R4_wrong_by_stage"][nm + " / " + k_] = w_
                st_rows.append(dict(run=nm, stage=k_, n=int(m_.sum()), wrong=wil_str(w_)))
        S["_stage_rows"] = st_rows
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
            timing_label=run["timing_label"],
            stages=stage_evals_bfs(run["stats"]),
            wrapper=run["wrapper"],
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
