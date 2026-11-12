"""AR3 step 11: numbers of the S8 baseline on ``held``, next to the three other samples.

  scripts/autoresearch/ar3/held_bg.sh held stats ar3/held_stats.py

Reads ``<cache>/S8_held/`` (``ar5/fit_run.py --name S8``, both changes off) and its pass
``<cache>/passes_held/`` (``ar5/fit_pass.py --name S8``); the other samples' numbers are the
stored S8 numbers in ``benchmarks/autoresearch/ar5/fit/part2_S10_<sample>.json`` and the stored
passes (read only). Writes ``benchmarks/autoresearch/ar3/held/``: ``summary.json``, ``tables.md``
and ``errors.npz`` (float32 errors of the maps before and after the pass). Truth is used for
scoring only. A run at its cap (``capped.json``) is written as unfinished with no scores.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "ar5"))
import fit_common as FC  # noqa: E402
import fit_part2_stats as P2  # noqa: E402

import gradient_common as G  # noqa: E402

import doc_tables as DT  # noqa: E402

N = 24570
HELD = FC.ICE / "benchmarks" / "autoresearch" / "ar3" / "held"
AR5 = FC.ICE / "benchmarks" / "autoresearch" / "ar5" / "fit"
NAMES = {"big": "big", "grad": "first gradient", "mid": "mid", "held": "held"}
GRAINS = {"big": 32, "grad": 497, "mid": 120, "held": 250}


def other(s: str) -> Dict[str, Any]:
    d = json.loads((AR5 / f"part2_S10_{s}.json").read_text())
    pj = json.loads((FC.STORED_PASS[s] / "D4_S8" / "pass.json").read_text())
    return dict(
        stages=d["stages_s8"],
        grid_fallbacks=d["grid_fallbacks_s8"],
        before=d["s8_before"],
        after=d["s8_after"],
        n_boundary=int(pj["n_boundary_voxels"]),
    )


def pct(k: int, n: int) -> str:
    lo, hi = P2.ST.wilson(k, n)
    return f"{k} ({100 * lo:.2f} to {100 * hi:.2f}%)"


def sample_table(sample: Dict[str, Any], check: Dict[str, Any], stack: Dict[str, Any]) -> str:
    sz = check["grain_size_min_median_max"]
    return DT.markdown_table(
        [
            dict(
                grains=sample["n_grains"],
                smin=sz[0],
                smed=sz[1],
                smax=sz[2],
                lit=stack["render_lit_clean_total"],
                spots=stack["render_spots_total"],
                err=check["max_abs_spread_minus_level_deg"],
            )
        ],
        ["grains", "smin", "smed", "smax", "lit", "spots", "err"],
        formats=dict(lit=",d", spots=",d", smed=".1f", err=".5f"),
    )


def first_run_counters() -> Any:
    """Counters of the first, capped run (120 minutes), moved to ``S8_held_cap120`` (D-N17)."""
    p = FC.CACHE / "S8_held_cap120" / "capped.json"
    if not p.exists():
        return None
    c = json.loads(p.read_text())
    s = c["stats"]
    return dict(
        cap_s=c["cap_s"],
        evals_before_pass=int(
            s["index_evals"]
            + s["index_fallback_grid_evals"]
            + s["n_evals_neighbor"]
            + s["n_evals_revisit"]
        ),
        grid_fallbacks=int(s["index_fallback_seed_test"]),
        first_seeds=int(s["n_seeds"]),
        reseeds=int(s["n_reseeds"]),
    )


def unfinished(out: Dict[str, Any], rd: Path) -> None:
    """The run reached its cap (or has no result): counters at the cap, no scores, no pass."""
    cap = json.loads((rd / "capped.json").read_text()) if (rd / "capped.json").exists() else None
    out["unfinished"] = True
    out["capped"] = cap
    sample = json.loads((G.SAMPLE / "sample.json").read_text())
    check = json.loads((G.SAMPLE / "held_check.json").read_text())
    stack = json.loads((G.GRAD / "held_stack.json").read_text())
    tabs: Dict[str, str] = {"ar3_held_sample": sample_table(sample, check, stack)}
    if cap is not None:
        s = cap["stats"]
        st = dict(
            seed_step=int(s["index_evals"]),
            grid_fallback=int(s["index_fallback_grid_evals"]),
            neighbour=int(s["n_evals_neighbor"]),
            revisit=int(s["n_evals_revisit"]),
        )
        st["before_pass"] = sum(st.values())
        assert st["seed_step"] + st["grid_fallback"] == s["n_evals_seed"] + s["n_evals_reseed"]
        out["at_cap"] = dict(
            stages=st,
            grid_fallbacks=int(s["index_fallback_seed_test"]),
            first_seeds=int(s["n_seeds"]),
            reseeds=int(s["n_reseeds"]),
            index_runs=int(s["index_runs"]),
            grid_share_of_evals_pct=round(100.0 * st["grid_fallback"] / st["before_pass"], 1),
            grid_evals_per_fallback=round(st["grid_fallback"] / s["index_fallback_seed_test"]),
            voxels_with_record=int(cap["voxels_with_record"]),
            wall_reseed_s=float(s["wall_reseed_total_s"]),
            wall_neighbour_s=float(s["wall_neighbor_s"]),
            wall_revisit_s=float(s["wall_revisit_s"]),
            wall_grid_s=float(s["index_fallback_grid_wall_s"]),
            cap_s=float(cap["cap_s"]),
            wall_label="contended (other jobs on the machine); not a timing result",
        )
        out["predictions"] = dict(
            p1_evals_20_to_40_million="not met: the count at the cap is already above 40 million",
            p2_wrong_after_30_to_120="not evaluated: no finished map, no pass",
            p3_grid_fallbacks_5_to_40="not met: the count at the cap is already above 40",
            evals_before_pass_at_cap=st["before_pass"],
            grid_fallbacks_at_cap=int(s["index_fallback_seed_test"]),
        )
        rows: List[Dict[str, Any]] = []
        for k in ("big", "grad", "mid"):
            d = other(k)
            x = d["stages"]
            rows.append(
                dict(
                    sample=NAMES[k],
                    state="finished",
                    seed=x["seed_step"],
                    grid=x["grid_fallback"],
                    neighbour=x["neighbour"],
                    revisit=x["revisit"],
                    before=x["before_pass"],
                    fallbacks=d["grid_fallbacks"],
                )
            )
        rows.append(
            dict(
                sample="held",
                state="stopped at 120 minutes",
                seed=st["seed_step"],
                grid=st["grid_fallback"],
                neighbour=st["neighbour"],
                revisit=st["revisit"],
                before=st["before_pass"],
                fallbacks=int(s["index_fallback_seed_test"]),
            )
        )
        tabs["ar3_held_capped"] = DT.markdown_table(
            rows,
            ["sample", "state", "seed", "grid", "neighbour", "revisit", "before", "fallbacks"],
            formats=dict(
                seed=",d", grid=",d", neighbour=",d", revisit=",d", before=",d", state="s"
            ),
        )
    (HELD / "summary.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    DT.write_tables(HELD / "tables.md", tabs)
    print("unfinished", json.dumps(out.get("predictions", {})))


def main() -> None:
    rd = FC.CACHE / "S8_held"
    pd = FC.CACHE / "passes_held"
    out: Dict[str, Any] = dict(run="S8", sample="held")
    if (rd / "capped.json").exists() or not (rd / "run.json").exists():
        unfinished(out, rd)
        return
    T = G.Truth()
    run = json.loads((rd / "run.json").read_text())
    fit = json.loads((rd / "fit.json").read_text())
    pj = json.loads((pd / "D4_S8" / "pass.json").read_text())
    Rb, _ = G.read_map(rd / "recon.mic.gz", N)
    Ra = np.load(pd / "after_S8.npy")
    eb, ea = P2.score(Rb, T), P2.score(Ra, T)
    np.savez_compressed(
        HELD / "errors.npz", before=eb.astype(np.float32), after=ea.astype(np.float32)
    )
    st = P2.stage_evals(run, int(pj["evals"]))
    held = dict(
        stages=st,
        grid_fallbacks=int(run["stats"].get("index_fallback_seed_test", 0)),
        before=P2.side(eb),
        after=P2.side(ea),
        n_boundary=int(pj["n_boundary_voxels"]),
    )
    sample = json.loads((G.SAMPLE / "sample.json").read_text())
    check = json.loads((G.SAMPLE / "held_check.json").read_text())
    stack = json.loads((G.GRAD / "held_stack.json").read_text())
    wb, wa = eb > G.WRONG, ea > G.WRONG
    by_level = {}
    for L in G.LEVELS:
        m = T.level == L
        by_level[str(L)] = dict(
            n_grains=int(sample["n_by_level"][str(L)]),
            n_voxels=int(m.sum()),
            wrong_before=int(wb[m].sum()),
            wrong_after=int(wa[m].sum()),
            fixed_by_pass=int((wb & ~wa & m).sum()),
            broken_by_pass=int((~wb & wa & m).sum()),
        )
    top = str(max(G.LEVELS))
    lower = [v for L, v in by_level.items() if L != top]
    gids, gsize = np.unique(T.gmap, return_counts=True)
    derived = dict(
        note="derived from the stored run files and the maps; quoted in the findings page",
        first_seeds=int(run["stats"]["n_seeds"]),
        reseeds=int(run["stats"]["n_reseeds"]),
        wall_min=round(float(fit["wall_s"]) / 60.0, 1),
        wall_pass_min=round(float(pj["wall_tries_s"]) / 60.0, 1),
        wall_label=(
            "contended: cpptools at 101% CPU at the start of the run and of the pass, and the "
            "three AR5 step 4 runs on held (D-20261111-N18) on the machine from 16:53; "
            "not a timing result"
        ),
        lower_levels=[float(L) for L in G.LEVELS if str(L) != top],
        lower_levels_wrong_before=int(sum(v["wrong_before"] for v in lower)),
        lower_levels_wrong_after=int(sum(v["wrong_after"] for v in lower)),
        top_level=float(top),
        evals_total_minus_grid=int(st["total"] - st["grid_fallback"]),
        grains_with_wrong_before=int(len(np.unique(T.gmap[wb]))),
        grains_with_wrong_after=int(len(np.unique(T.gmap[wa]))),
        small_grain_voxels=40,
        n_grains_under_small=int((gsize < 40).sum()),
        n_grains=int(len(gids)),
    )
    assert derived["first_seeds"] + derived["reseeds"] == fit["seed_steps"]
    out.update(
        unfinished=False,
        derived=derived,
        seed_steps=fit["seed_steps"],
        stopped_early=fit["stopped_early"],
        wrapper={k: v for k, v in run["wrapper"].items() if k != "wrapper_calls_by_caller"},
        held=held,
        by_level=by_level,
        wrong_both=int((wb & wa).sum()),
        wrong_fixed_by_pass=int((wb & ~wa).sum()),
        wrong_broken_by_pass=int((~wb & wa).sum()),
        sample=dict(
            n_grains=sample["n_grains"],
            grain_size_min_median_max=check["grain_size_min_median_max"],
            n_by_level=sample["n_by_level"],
            max_abs_spread_minus_level_deg=check["max_abs_spread_minus_level_deg"],
            lit_pixels=stack["render_lit_clean_total"],
            spots=stack["render_spots_total"],
        ),
        wall_s=fit["wall_s"],
        wall_label=fit["wall_label"],
        wall_pass_s=pj["wall_tries_s"],
        cap_min=fit["cap_min"],
        first_run_capped=first_run_counters(),
    )
    allv = {s: other(s) for s in ("big", "grad", "mid")}
    allv["held"] = held
    ev = float(st["total"])
    out["predictions"] = dict(
        p1_evals_20_to_40_million=bool(20e6 <= ev <= 40e6),
        p2_wrong_after_30_to_120=bool(30 <= held["after"]["wrong"] <= 120),
        p3_grid_fallbacks_5_to_40=bool(5 <= held["grid_fallbacks"] <= 40),
        evals_total=int(ev),
        wrong_after=held["after"]["wrong"],
        grid_fallbacks=held["grid_fallbacks"],
    )
    out["others"] = {s: allv[s] for s in ("big", "grad", "mid")}
    (HELD / "summary.json").write_text(json.dumps(out, indent=1, default=float) + "\n")

    rows_e: List[Dict[str, Any]] = []
    rows_w: List[Dict[str, Any]] = []
    for s in ("big", "grad", "mid", "held"):
        d = allv[s]
        x = d["stages"]
        rows_e.append(
            dict(
                sample=NAMES[s],
                grains=GRAINS[s],
                seed=x["seed_step"],
                grid=x["grid_fallback"],
                neighbour=x["neighbour"],
                revisit=x["revisit"],
                pass_=x["pass_"],
                total=x["total"],
                fallbacks=d["grid_fallbacks"],
            )
        )
        rows_w.append(
            dict(
                sample=NAMES[s],
                before=pct(d["before"]["wrong"], N),
                after=pct(d["after"]["wrong"], N),
                median_b=d["before"]["median_err_deg"],
                median_a=d["after"]["median_err_deg"],
                boundary=f"{d['n_boundary']} ({100 * d['n_boundary'] / N:.1f}%)",
            )
        )
    tabs: Dict[str, str] = {}
    tabs["ar3_held_evals"] = DT.markdown_table(
        rows_e,
        ["sample", "grains", "seed", "grid", "neighbour", "revisit", "pass_", "total", "fallbacks"],
        formats=dict(seed=",d", grid=",d", neighbour=",d", revisit=",d", pass_=",d", total=",d"),
    )
    tabs["ar3_held_wrong"] = DT.markdown_table(
        rows_w,
        ["sample", "before", "after", "median_b", "median_a", "boundary"],
        formats=dict(median_b=".4f", median_a=".4f"),
    )
    lv = [
        dict(
            level=L,
            grains=v["n_grains"],
            voxels=v["n_voxels"],
            wrong_before=v["wrong_before"],
            wrong_after=v["wrong_after"],
            fixed=v["fixed_by_pass"],
            broken=v["broken_by_pass"],
        )
        for L, v in by_level.items()
    ]
    tabs["ar3_held_levels"] = DT.markdown_table(
        lv,
        ["level", "grains", "voxels", "wrong_before", "wrong_after", "fixed", "broken"],
        formats=dict(level="s"),
    )
    tabs["ar3_held_sample"] = sample_table(sample, check, stack)
    DT.write_tables(HELD / "tables.md", tabs)
    print(json.dumps(out["predictions"], indent=1))


if __name__ == "__main__":
    main()
