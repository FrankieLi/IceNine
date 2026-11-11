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
    """Pack the maps (float32) so that the summary reads committed files only."""
    for k in ("S", "K2"):
        if not (RUNS / k / "run.json").exists():
            continue
        before = G.read_map(RUNS / k / "recon.mic.gz", N)[0]
        after = np.load(S.CACHE / "passes" / f"after_{k}.npy")
        np.savez_compressed(
            RUNS / k / "maps.npz", before=before.astype(np.float32), after=after.astype(np.float32)
        )
    # step 5's maps after the pass, checked against the committed tries (rule ``hit``)
    out = {}
    chk = {}
    for k in ("D4", "D"):
        R0 = B.load_map(k, N)
        t = BL.Tries(B.RUNS / k)
        redo, _ = BL.apply_rule(t, R0, "hit")
        a = np.load(OTHER / f"after_{k}.npy").astype(np.float64)
        d = B.orthonormal(redo)
        chk[k] = dict(
            max_abs_diff=float(np.abs(a - redo).max()),
            identical_voxels=int((np.abs(a - redo).reshape(N, -1).max(1) < 1e-6).sum()),
        )
        out[f"{k}_before"] = R0.astype(np.float32)
        out[f"{k}_after"] = redo.astype(np.float32)
        del d
    np.savez_compressed(RUNS / "compare.npz", **out)
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


def per_level(T: Any, R: np.ndarray) -> Dict[str, Any]:
    e = np.asarray(T.err(R), dtype=float)
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
    cmp = np.load(RUNS / "compare.npz")
    maps: Dict[str, np.ndarray] = {
        "D4 before": cmp["D4_before"].astype(np.float64),
        "D4 after": cmp["D4_after"].astype(np.float64),
        "D before": cmp["D_before"].astype(np.float64),
        "D after": cmp["D_after"].astype(np.float64),
    }
    info: Dict[str, Any] = {}
    for k in ("S", "K2"):
        if not (RUNS / k / "maps.npz").exists():
            continue
        z = np.load(RUNS / k / "maps.npz")
        maps[f"{k} before"], maps[f"{k} after"] = z["before"].astype(np.float64), z["after"].astype(
            np.float64
        )
        info[k] = json.loads((RUNS / k / "run.json").read_text())
    info["D4"] = json.loads((C.OUT / "fine" / "runs" / "D4" / "run.json").read_text())
    info["D"] = json.loads((G.GRAD / "runs" / "D" / "run.json").read_text())
    pass_ev = {
        "D4": json.loads((C.OUT / "boundary" / "runs" / "D4" / "pass.json").read_text())["evals"],
        "D": json.loads((C.OUT / "boundary" / "runs" / "D" / "pass.json").read_text())["evals"],
    }
    for k in ("S", "K2"):
        if k in info:
            p = (
                json.loads((S.CACHE / "passes" / f"D4_{k}" / "pass.json").read_text())
                if (S.CACHE / "passes" / f"D4_{k}" / "pass.json").exists()
                else None
            )
            pass_ev[k] = p["evals"] if p else (info[k].get("pass_evals") or 0)
    SUM: Dict[str, Any] = dict(pass_evals=pass_ev)
    T_: Dict[str, str] = {}
    # ---- stages
    names = [k for k in ("D", "D4", "S", "K2") if k in info]
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
    for nm, R in maps.items():
        pl = per_level(T, R)
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
    SUM["paired"] = prow
    T_["ar3f_paired"] = DT.markdown_table(
        prow,
        [
            "a",
            "b",
            "wrong_a",
            "wrong_b",
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
    for k in ("S", "K2"):
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
    print("part 2 summary written", [r["run"] for r in rows])


if __name__ == "__main__":
    main()
