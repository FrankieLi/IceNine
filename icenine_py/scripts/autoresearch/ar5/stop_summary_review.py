"""AR5 step 1: numbers added in the fixes round (review 20261111-1215, decision D-20261111-N11).

Called by ``stop_summary.py``; writes nothing itself. Every number comes from committed files:
the result files of ``benchmarks/autoresearch/ar5/`` (``starts_*.npz``, ``part1_*.json``,
``part2_cases.json``, ``part3_flags.json``, ``part3_votes.npz``), the sample files of AR3
(``benchmarks/autoresearch/ar3/{big,gradient,mid}/.../sample.json``: grains and gradient levels) and
the S8 totals of AR3 step 9 (``benchmarks/autoresearch/ar3/start/summary.json``). Nothing here is
a new run. The what-ifs are labelled as what-ifs in the returned dict and on the page.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
AR3 = ICE / "benchmarks" / "autoresearch" / "ar3"
sys.path.insert(0, str(HERE.parents[1] / "common"))
sys.path.insert(0, str(HERE))
import doc_tables as DT  # noqa: E402
import stop_part1_lib as P1  # noqa: E402
import stop_rules as SR  # noqa: E402

SAMPLES = ("big", "grad", "mid")
PICK_SAMPLES = ("big", "grad")
SR_WRONG = 1.0  # degrees; as WRONG in stop_common and stop_part1_lib
SR_AR5 = P1.AR5
SAMPLE_FILES = {
    "big": AR3 / "big" / "gradient" / "sample" / "sample.json",
    "grad": AR3 / "gradient" / "sample" / "sample.json",
    "mid": AR3 / "mid" / "gradient" / "sample" / "sample.json",
}


def _wilson(k: int, n: int) -> Dict[str, Any]:
    lo, hi = SR.wilson_row(k, n)[2:]
    return dict(k=int(k), n=int(n), share=k / n, lo=float(lo), hi=float(hi))


def scope() -> Dict[str, Any]:
    """Grains and gradient levels of each sample; the evaluations of the stored S8 runs."""
    out: Dict[str, Any] = {}
    s8 = json.loads((AR3 / "start" / "summary.json").read_text())["part3"]
    for s in SAMPLES:
        d = json.loads(SAMPLE_FILES[s].read_text())
        out[s] = dict(
            n_grains=int(d["n_grains"]),
            n_voxels=int(d["n_voxels"]),
            levels_deg=[float(x) for x in d["levels_deg"]],
            max_level_deg=float(max(d["levels_deg"])),
            s8_total_evals=int(s8[s]["S8"]["total"]),
        )
    return out


def part1(p1: Dict[str, Any], p1p: Dict[str, Any]) -> Dict[str, Any]:
    """The pick's margin; P3 at every d (the one free choice); every plan setting on mid with
    hindsight (a what-if: the pick never read mid)."""
    full_gain = p1p["full_pass"]["net_gain"]
    need = (1.0 - SR.MAX_LOSS) * full_gain
    pick_row = next(r for r in p1p["rows"] if r["name"] == p1p["pick"]["name"])
    loaded = {s: P1.load(s) for s in SAMPLES}
    p3_rows: List[Dict[str, Any]] = []
    for d in SR.D_VALUES:
        for m in SR.M_VALUES:
            per = {s: P1.evaluate(loaded[s], "P3", d, m) for s in SAMPLES}
            pooled_gain = sum(per[s]["net_gain"] for s in PICK_SAMPLES)
            pooled_ev = sum(per[s]["evals"] for s in PICK_SAMPLES)
            mid_full = p1["table"]["full"]["mid"]
            p3_rows.append(
                dict(
                    d=float(d),
                    m=int(m),
                    pick_samples_gain=int(pooled_gain),
                    pick_samples_evals=int(pooled_ev),
                    reaches_needed_gain=bool(pooled_gain >= need - SR.TOL),
                    mid_gain=int(per["mid"]["net_gain"]),
                    mid_lost_pc=100.0 * (1 - per["mid"]["net_gain"] / mid_full["net_gain"]),
                    mid_cut_pc=100.0 * (1 - per["mid"]["evals"] / mid_full["evals"]),
                )
            )
    mid_full = p1["table"]["full"]["mid"]
    hind: Dict[str, Any] = {}
    for name, per in p1["table"].items():
        if name == "full":
            continue
        r = per["mid"]
        cut = 1 - r["evals"] / mid_full["evals"]
        lost = 1 - r["net_gain"] / mid_full["net_gain"]
        hind[name] = dict(
            cut_pc=100.0 * cut,
            lost_pc=100.0 * lost,
            passes=bool(cut >= 0.5 - SR.TOL and lost <= SR.MAX_LOSS + SR.TOL),
        )
    p13 = [v["lost_pc"] for k, v in hind.items() if k.startswith(("P1", "P3"))]
    p13_all_d = [r["mid_lost_pc"] for r in p3_rows] + [
        v["lost_pc"] for k, v in hind.items() if k.startswith("P1")
    ]
    return dict(
        pick_margin=dict(
            full_gain=int(full_gain),
            needed_gain=float(need),
            pick_gain=int(pick_row["net_gain"]),
            margin_voxels=float(pick_row["net_gain"] - need),
        ),
        p3_every_d=p3_rows,
        p3_every_d_max_pick_samples_gain=int(max(r["pick_samples_gain"] for r in p3_rows)),
        p3_every_d_any_reaches=bool(any(r["reaches_needed_gain"] for r in p3_rows)),
        hindsight_mid_whatif=dict(
            settings=hind,
            any_plan_setting_passes=bool(any(v["passes"] for v in hind.values())),
            min_lost_pc_p1_p3_plan=float(min(p13)),
            min_lost_pc_p1_p3_every_d=float(min(p13_all_d)),
        ),
    )


def part2(cases: Dict[str, List[Dict[str, Any]]]) -> Tuple[Dict[str, Any], Dict[str, str]]:
    allc = [c for s in ("grad", "mid") for c in cases[s]]
    n = len(allc)
    n_cand = [int(c["n_cand"]) for c in allc]
    r1_all = sum(1 for c in allc if int(c["R1"]["costed"]) == int(c["n_cand"]))
    repeats = sum(
        1
        for c in allc
        if len(c["R1"]["fin_hit"]) == len(c["stored_fin_hit"])
        and all(abs(a - b) <= SR.TOL for a, b in zip(c["R1"]["fin_hit"], c["stored_fin_hit"]))
    )
    grid_ev = np.array([c["grid_evals"] for c in allc], float)
    lad_ev = float(sum(c["ladder"]["evals"] for c in allc))
    grid_tot = float(grid_ev.sum())
    at_grid = [c for c in allc if c["ladder"]["rung"] == "grid"]
    l0_at_grid = float(sum(c["R3"]["level0_evals"] for c in at_grid))
    r1_tot = float(sum(c["R1"]["evals"] for c in allc))
    reuse = lad_ev - l0_at_grid
    reuse_no_r1 = reuse - r1_tot
    whatif = dict(
        label="what-if, not a result: a ladder whose grid search reuses R3's level 0",
        assumption="the grid search's own level 0 costs what R3's level 0 cost (R3 used "
        "default_rng(0); the stored grid search used another random state)",
        ladder_evals=lad_ev,
        grid_evals=grid_tot,
        two_thirds_limit=2.0 / 3.0 * grid_tot,
        cases_ending_at_grid=len(at_grid),
        r3_level0_evals_in_those=l0_at_grid,
        reuse_level0_evals=reuse,
        reuse_level0_share_pc=100.0 * reuse / grid_tot,
        r1_evals_all_cases=r1_tot,
        reuse_level0_no_r1_evals=reuse_no_r1,
        reuse_level0_no_r1_share_pc=100.0 * reuse_no_r1 / grid_tot,
        under_two_thirds=bool(reuse <= 2.0 / 3.0 * grid_tot),
        under_two_thirds_no_r1=bool(reuse_no_r1 <= 2.0 / 3.0 * grid_tot),
    )
    counts = dict(
        r1_right=_wilson(sum(1 for c in allc if c["R1"]["err"] <= SR_WRONG), n),
        r3_pass=_wilson(sum(1 for c in allc if c["R3"]["hit"] >= SR.SEED_TEST - SR.TOL), n),
        r3_right=_wilson(sum(1 for c in allc if c["R3"]["err"] <= SR_WRONG), n),
        r3_right_when_pass=_wilson(
            sum(
                1
                for c in allc
                if c["R3"]["hit"] >= SR.SEED_TEST - SR.TOL and c["R3"]["err"] <= SR_WRONG
            ),
            sum(1 for c in allc if c["R3"]["hit"] >= SR.SEED_TEST - SR.TOL),
        ),
        ladder_right=_wilson(sum(1 for c in allc if c["ladder"]["right"]), n),
        grid_right=_wilson(sum(1 for c in allc if c["grid"]["err"] <= SR_WRONG), n),
    )
    labels = {
        "r1_right": "R1 end right",
        "r3_pass": "R3 passes the seed test",
        "r3_right": "R3 end right",
        "r3_right_when_pass": "R3 end right, of the cases it passes",
        "ladder_right": "ladder end right",
        "grid_right": "grid search end right",
    }
    rows = [
        dict(count=labels[k], k=v["k"], n=v["n"], lo=100 * v["lo"], hi=100 * v["hi"])
        for k, v in counts.items()
    ]
    tables = {
        "p2_intervals": DT.markdown_table(
            rows, ["count", "k", "n", "lo", "hi"], {"lo": ".0f", "hi": ".0f"}
        )
    }
    out = dict(
        n_cases=n,
        n_cand_min=min(n_cand),
        n_cand_max=max(n_cand),
        r1_costed_every_candidate=r1_all,
        r1_repeats_seed_step_refinements=repeats,
        grid_evals_mean=float(grid_ev.mean()),
        grid_evals_min=int(grid_ev.min()),
        grid_evals_max=int(grid_ev.max()),
        level0_reuse_whatif=whatif,
        intervals=counts,
    )
    return out, tables


def part3(p3f: Dict[str, Any], thr: float) -> Dict[str, Any]:
    key, lower = f"{thr:.2f}", f"{thr - 0.05:.2f}"
    flags: Dict[str, Any] = {}
    for s in SAMPLES:
        t = p3f["tables"][s]
        flags[s] = dict(
            flagged=int(t[key]["all"]["flagged"]),
            flagged_between_lower_and_thr=int(
                t[key]["all"]["flagged"] - t[lower]["all"]["flagged"]
            ),
            lower=float(lower),
            boundary_share_pc=100.0 * t[key]["boundary"]["n"] / t[key]["all"]["n"],
        )
    return flags


def votes() -> Dict[str, Any]:
    z = np.load(SR_AR5 / "part3_votes.npz")

    def split(samples: Tuple[str, ...]) -> Tuple[np.ndarray, np.ndarray]:
        w = np.concatenate([z[f"{s}_votes"][z[f"{s}_flagged"] & z[f"{s}_wrong"]] for s in samples])
        r = np.concatenate([z[f"{s}_votes"][z[f"{s}_flagged"] & ~z[f"{s}_wrong"]] for s in samples])
        return w, r

    cuts = {s: SR.best_cut(*split((s,)))["cut"] for s in SAMPLES}
    wp, rp = split(PICK_SAMPLES)
    pick_cut = SR.best_cut(wp, rp)["cut"]
    wm, rm = split(("mid",))
    wg, rg = int((wm >= pick_cut).sum()), int((rm >= pick_cut).sum())
    wa, ra = split(SAMPLES)
    all_cut = SR.best_cut(wa, ra)
    above = all_cut["wrong_ge"] + all_cut["right_ge"]
    return dict(
        cut_per_sample=cuts,
        cut_min=float(min(cuts.values())),
        cut_max=float(max(cuts.values())),
        pick_sample_cut_whatif=dict(
            label="what-if: the cut chosen on big and grad pooled, read on mid",
            cut=float(pick_cut),
            mid_wrong=_wilson(wg, len(wm)),
            mid_right=_wilson(rg, len(rm)),
            separates=SR.separates(wg, len(wm), rg, len(rm)),
        ),
        all_three=dict(
            cut=float(all_cut["cut"]),
            wrong_among_above=_wilson(all_cut["wrong_ge"], above),
            wrong_among_flagged=_wilson(len(wa), len(wa) + len(ra)),
        ),
    )


def review(
    p1: Dict[str, Any],
    p1p: Dict[str, Any],
    cases: Dict[str, List[Dict[str, Any]]],
    p3f: Dict[str, Any],
    thr: float,
) -> Tuple[Dict[str, Any], Dict[str, str]]:
    p2, tables = part2(cases)
    out = dict(
        scope=scope(),
        part1=part1(p1, p1p),
        part2=p2,
        part3_flags=part3(p3f, thr),
        part3_votes=votes(),
    )
    return out, tables
