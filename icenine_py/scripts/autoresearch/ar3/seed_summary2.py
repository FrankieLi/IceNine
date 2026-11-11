"""AR3 step 6, part 2: runs S and K2 (picked seed remedy) with the boundary pass, against D4 and
step 3's D with the same pass.

  uv run python scripts/autoresearch/ar3/seed_summary2.py

Reads (committed) ``seed/runs/<S|K2>/{run.json,seed_log.json.gz}``, ``seed/runs/<kind>/after.npz``
(the maps before and after the pass, float32, written by ``--pack``), the pass records of step 5
(``boundary/runs/{D4,D}/pass.json``) and, for the two comparison maps after the pass, step 5's
untracked ``after_D4.npy`` and ``after_D.npy`` (packed into ``seed/runs/compare.npz`` by
``--pack`` so that this summary can be rerun from committed files). Error is against the truth
(analysis only); wrong is an error above 1 degree. Wilson 95% intervals, no p-values. All wall
times were measured with other jobs running and are not timing results.
"""

import argparse
import json
from typing import Any, Dict, List, Optional

import ar3_common as C
import boundary_common as B
import boundary_lib as BL
import gradient_common as G
import numpy as np
import seed_common as S
from seed_summary import wil, ws

import doc_tables as DT  # noqa: E402

RUNS = S.SEED / "runs"
OTHER = (
    S.WORKTREES / "agent-a4fd6282107a76186/icenine_py/benchmarks/autoresearch/ar3/boundary/cache"
)
N = 24570


def pack() -> None:
    """Pack the error of every voxel against the truth (float32, one array per map) so that the
    summary reads committed files only; the maps themselves are not committed."""
    T = G.Truth()
    out = {}
    for k in ("S", "K2", "S5"):
        if not (RUNS / k / "run.json").exists():
            continue
        before = G.read_map(RUNS / k / "recon.mic.gz", N)[0]
        after = np.load(S.CACHE / "passes" / f"after_{k}.npy").astype(np.float64)
        out[f"{k} before"] = np.asarray(T.err(before), np.float32)
        out[f"{k} after"] = np.asarray(T.err(after), np.float32)
    # step 5's maps after the pass, checked against the committed tries (rule ``hit``)
    chk = {}
    for k in ("D4", "D"):
        R0 = B.load_map(k, N)
        t = BL.Tries(B.RUNS / k)
        redo, _ = BL.apply_rule(t, R0, "hit")
        a = np.load(OTHER / f"after_{k}.npy").astype(np.float64)
        chk[k] = dict(
            max_abs_diff=float(np.abs(a - redo).max()),
            identical_voxels=int((np.abs(a - redo).reshape(N, -1).max(1) < 1e-6).sum()),
        )
        out[f"{k} before"] = np.asarray(T.err(R0), np.float32)
        out[f"{k} after"] = np.asarray(T.err(redo), np.float32)
    np.savez_compressed(RUNS / "errors.npz", **out)
    (RUNS / "compare_check.json").write_text(json.dumps(chk) + "\n")
    print("packed", chk)


def stage_evals(st: Dict[str, Any], pass_evals: int) -> Dict[str, int]:
    grid = int(st.get("index_fallback_grid_evals", 0))
    seed_all = int(st["n_evals_seed"]) + int(st.get("n_evals_reseed", 0))
    return dict(
        seed_step=seed_all - grid,
        grid_fallback=grid,
        neighbour=int(st["n_evals_neighbor"]),
        revisit=int(st["n_evals_revisit"]) + int(st.get("n_evals_restart", 0)),
        pass_=int(pass_evals),
        total=seed_all
        + int(st["n_evals_neighbor"])
        + int(st["n_evals_revisit"])
        + int(st.get("n_evals_restart", 0))
        + int(pass_evals),
    )


def per_level(T: Any, err: np.ndarray) -> Dict[str, Any]:
    e = np.asarray(err, dtype=float)
    out: Dict[str, Any] = dict(
        all=dict(n=len(e), wrong=int((e > G.WRONG).sum()), median=float(np.median(e)))
    )
    for L in G.LEVELS:
        m = T.level == L
        out[f"level {L:g}"] = dict(
            n=int(m.sum()), wrong=int((e[m] > G.WRONG).sum()), median=float(np.median(e[m]))
        )
    return dict(rows=out, err=e)


def paired(ea: np.ndarray, eb: np.ndarray, sel: Optional[np.ndarray] = None) -> Dict[str, int]:
    """a against b (b the comparison): fixed = wrong in b, right in a."""
    wa, wb = ea > G.WRONG, eb > G.WRONG
    if sel is not None:
        wa, wb = wa[sel], wb[sel]
    return dict(
        n=int(len(wa)),
        wrong_a=int(wa.sum()),
        wrong_b=int(wb.sum()),
        fixed=int((wb & ~wa).sum()),
        broken=int((~wb & wa).sum()),
        both_wrong=int((wa & wb).sum()),
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pack", action="store_true")
    a = ap.parse_args()
    if a.pack:
        pack()
    T = G.Truth()
    z = np.load(RUNS / "errors.npz")
    order = [
        "D before",
        "D after",
        "D4 before",
        "D4 after",
        "S before",
        "S after",
        "K2 before",
        "K2 after",
        "S5 before",
        "S5 after",
    ]
    maps: Dict[str, np.ndarray] = {k: z[k].astype(np.float64) for k in order if k in z.files}
    info: Dict[str, Any] = {}
    for k in ("S", "K2", "S5"):
        if f"{k} after" not in maps:
            continue
        info[k] = json.loads((RUNS / k / "run.json").read_text())
    info["D4"] = json.loads((C.OUT / "fine" / "runs" / "D4" / "run.json").read_text())
    info["D"] = json.loads((G.GRAD / "runs" / "D" / "run.json").read_text())
    pass_ev = {
        "D4": json.loads((C.OUT / "boundary" / "runs" / "D4" / "pass.json").read_text())["evals"],
        "D": json.loads((C.OUT / "boundary" / "runs" / "D" / "pass.json").read_text())["evals"],
    }
    for k in ("S", "K2", "S5"):
        if k in info:
            # the pass record is copied from seed/cache/passes/D4_<kind>/ to seed/runs/<kind>/
            pass_ev[k] = json.loads((RUNS / k / "pass.json").read_text())["evals"]
    SUM: Dict[str, Any] = dict(pass_evals=pass_ev)
    T_: Dict[str, str] = {}
    # ---- stages
    names = [k for k in ("D", "D4", "S", "K2", "S5") if k in info]
    rows = []
    for k in names:
        st = info[k]["stats"]
        se = stage_evals(st, pass_ev[k])
        rows.append(
            dict(
                run=k + " + pass",
                seeds=int(st["n_seeds"]),
                reseeds=int(st["n_reseeds"]),
                grid_fallbacks=int(st.get("index_fallback_seed_test", 0)),
                wall_min=info[k].get("wall_total_s", 0.0) / 60.0,
                **se,
            )
        )
    SUM["stages"] = rows
    T_["ar3f_stages"] = DT.markdown_table(
        rows,
        [
            "run",
            "seeds",
            "reseeds",
            "grid_fallbacks",
            "seed_step",
            "grid_fallback",
            "neighbour",
            "revisit",
            "pass_",
            "total",
            "wall_min",
        ],
        formats={
            c: ",d"
            for c in ("seed_step", "grid_fallback", "neighbour", "revisit", "pass_", "total")
        }
        | {"wall_min": ".1f"},
    )
    # ---- errors
    errs: Dict[str, np.ndarray] = {}
    levels = []
    for nm, e0 in maps.items():
        pl = per_level(T, e0)
        errs[nm] = pl["err"]
        for lv, d in pl["rows"].items():
            levels.append(
                dict(
                    map=nm,
                    level=lv,
                    n=d["n"],
                    wrong=wil(d["wrong"], d["n"]),
                    wrong_k=d["wrong"],
                    median=d["median"],
                )
            )
    SUM["levels"] = [{k: v for k, v in r.items()} for r in levels]
    rl = []
    for lv in ["all"] + [f"level {L:g}" for L in G.LEVELS]:
        row: Dict[str, Any] = dict(level=lv)
        for nm in maps:
            r = next(x for x in levels if x["map"] == nm and x["level"] == lv)
            row[nm] = f"{r['wrong_k']} ({r['median']:.3f})"
        rl.append(row)
    T_["ar3f_levels2"] = DT.markdown_table(rl, ["level"] + list(maps))
    # ---- paired, the pass on both sides
    pairs = [
        ("S after", "D4 after"),
        ("S5 after", "D4 after"),
        ("S5 after", "D after"),
        ("S5 after", "S after"),
        ("S5 before", "D4 before"),
        ("S after", "D after"),
        ("K2 after", "S after"),
        ("K2 after", "D4 after"),
        ("K2 after", "D after"),
        ("S before", "D4 before"),
        ("K2 before", "S before"),
    ]
    prow = []
    for a_, b_ in pairs:
        if a_ not in maps or b_ not in maps:
            continue
        p = paired(errs[a_], errs[b_])
        p["a"], p["b"] = a_, b_
        prow.append(p)
        isA = T.gmap % 2 == 0
        for sn, sel in (("A", isA), ("B", ~isA)):
            q = paired(errs[a_], errs[b_], sel)
            p[f"fixed_{sn}"], p[f"broken_{sn}"] = q["fixed"], q["broken"]
    for p in prow:
        wa0 = wil(p["wrong_a"], p["n"])
        p["wrong_a_s"] = f"{wa0['k']} ({1e4 * wa0['lo']:.0f} to {1e4 * wa0['hi']:.0f} per 10,000)"
        wb0 = wil(p["wrong_b"], p["n"])
        p["wrong_b_s"] = f"{wb0['k']} ({1e4 * wb0['lo']:.0f} to {1e4 * wb0['hi']:.0f} per 10,000)"
        wa_, wb_ = wil(p["wrong_a"], p["n"]), wil(p["wrong_b"], p["n"])
        p["intervals_overlap"] = bool(not (wa_["hi"] < wb_["lo"] or wb_["hi"] < wa_["lo"]))
    SUM["paired"] = prow
    # where S's and K2's seed results were made (tries or seed step) and how often they are wrong
    # in the run's own map, against the same voxels in D4's map (before the pass)
    rows_by = []
    for k in ("S", "K2", "S5"):
        lg = (
            G.load_gz(RUNS / k / "seed_log.json.gz")
            if (RUNS / k / "seed_log.json.gz").exists()
            else []
        )
        for by in ("tries", "seed_step"):
            vs = np.array(sorted({int(e["voxel"]) for e in lg if e["by"] == by}), dtype=int)
            if not len(vs):
                continue
            ek = maps[f"{k} before"][vs]
            e4 = maps["D4 before"][vs]
            rows_by.append(
                dict(
                    run=k,
                    settled_by=by,
                    voxels=len(vs),
                    wrong_in_run=int((ek > G.WRONG).sum()),
                    wrong_in_D4=int((e4 > G.WRONG).sum()),
                    run_wrong_D4_right=int(((ek > G.WRONG) & ~(e4 > G.WRONG)).sum()),
                    run_right_D4_wrong=int((~(ek > G.WRONG) & (e4 > G.WRONG)).sum()),
                )
            )
    SUM["by_path"] = rows_by
    if rows_by:
        T_["ar3f_bypath"] = DT.markdown_table(
            rows_by,
            [
                "run",
                "settled_by",
                "voxels",
                "wrong_in_run",
                "wrong_in_D4",
                "run_wrong_D4_right",
                "run_right_D4_wrong",
            ],
        )
    T_["ar3f_paired"] = DT.markdown_table(
        prow,
        [
            "a",
            "b",
            "wrong_a_s",
            "wrong_b_s",
            "fixed",
            "broken",
            "both_wrong",
            "fixed_A",
            "broken_A",
            "fixed_B",
            "broken_B",
        ],
    )
    # ---- the picked remedy's own counters
    prow2 = []
    for k in ("S", "K2", "S5"):
        if k in info:
            pk = info[k]["pick"]
            st = info[k]["stats"]
            prow2.append(
                dict(
                    run=k,
                    seed_steps=pk["seed_steps"],
                    by_tries=pk["settled_by_tries"],
                    to_seed_step=pk["went_to_seed_step"],
                    tries_evals=pk["tries_evals"],
                    grid_fallbacks=int(st.get("index_fallback_seed_test", 0)),
                    grid_evals=int(st.get("index_fallback_grid_evals", 0)),
                    seed_evals=pk["evals"],
                )
            )
    SUM["pick_counters"] = prow2
    if prow2:
        T_["ar3f_pickrun"] = DT.markdown_table(
            prow2,
            [
                "run",
                "seed_steps",
                "by_tries",
                "to_seed_step",
                "tries_evals",
                "grid_fallbacks",
                "grid_evals",
                "seed_evals",
            ],
            formats=dict(tries_evals=",d", grid_evals=",d", seed_evals=",d"),
        )
    DT.write_tables(S.SEED / "tables_part2.md", T_)
    (S.SEED / "summary_part2.json").write_text(json.dumps(SUM, indent=1, default=float) + "\n")
    figure2(rows, maps)
    print("part 2 summary written", [r["run"] for r in rows])


def figure2(rows: List[Dict[str, Any]], maps: Dict[str, np.ndarray]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
    stages = [
        ("seed_step", "seed step"),
        ("grid_fallback", "grid fallbacks"),
        ("neighbour", "neighbour fits"),
        ("revisit", "revisit fits"),
        ("pass_", "boundary pass"),
    ]
    bottom = np.zeros(len(rows))
    for key, lab in stages:
        v = np.array([r[key] for r in rows], dtype=float) / 1e6
        ax[0].bar([r["run"] for r in rows], v, bottom=bottom, label=lab)
        bottom += v
    ax[0].set_ylabel("evaluations (millions)")
    ax[0].set_title("Where the evaluations go (each run plus the boundary pass)")
    ax[0].legend(fontsize=8)
    for i, r in enumerate(rows):
        nm = r["run"].replace(" + pass", "") + " after"
        w = wil(int((maps[nm] > G.WRONG).sum()), len(maps[nm]))
        p = w["k"] / w["n"]
        ax[1].errorbar(i, 1e4 * p, yerr=[[1e4 * (p - w["lo"])], [1e4 * (w["hi"] - p)]], fmt="o")
    ax[1].set_xticks(range(len(rows)))
    ax[1].set_xticklabels([r["run"] for r in rows])
    ax[1].set_ylabel("wrong voxels per 10,000 (Wilson 95%)")
    ax[1].set_title("Wrong voxels after the pass (error above 1 degree)")
    ax[1].set_ylim(bottom=0)
    fig.tight_layout()
    fig.savefig(S.SEED / "ar3_seed_runs.png", dpi=110)


if __name__ == "__main__":
    main()
