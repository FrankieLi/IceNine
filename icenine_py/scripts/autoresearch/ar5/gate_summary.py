"""AR5 step 7, Phase C: the tables, the plot, the verdicts and the pareto rows of the skip gate.

  uv run python scripts/autoresearch/ar5/gate_summary.py

Reads ``benchmarks/autoresearch/ar5/gate/`` (``replay.json``, ``replay_held.json``,
``live_<sample>_<run>_h<h>.json``, ``frozen.json``) and writes there ``tables.md`` (marker blocks
for ``scripts/dev/sync_doc_tables.py``), ``summary.json``, ``ar5_gate_curve.png`` and
``pareto_rows.json`` (rows for the review page; not merged into the shared pareto file here).
Predictions and the decision rule are the frozen ones; verdicts are computed, not typed.
"""

import json
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import gate_common as GC
import gate_rules as GR

from doc_tables import markdown_table, write_tables  # noqa: E402

G = GC.GATE
ENTRIES: List[Tuple[str, str]] = GC.DISCLOSED + [("held", "S11v2")]
LIVE: List[Tuple[str, str]] = [("mid", "S11n"), ("grad", "S11n"), ("big", "S11"), ("held", "S11v2")]
COLORS = {"big": "tab:purple", "grad": "tab:green", "mid": "tab:blue", "held": "tab:red"}
NAME = {"big": "big", "grad": "first gradient", "mid": "mid", "held": "held"}


def load(p: Path) -> Any:
    return json.loads(p.read_text())


def rows_replay() -> Dict[Tuple[str, str], Dict[str, Any]]:
    out = {(r["sample"], r["run"]): r for r in load(G / "replay.json")["rows"]}
    hp = G / "replay_held.json"
    if hp.exists() and load(hp).get("old_score_stored"):
        r = load(hp)["row"]
        out[(r["sample"], r["run"])] = r
    return out


def at(r: Dict[str, Any], h: float = GR.H_S) -> Dict[str, Any]:
    return [x for x in r["by_h"] if abs(x["h"] - h) < 1e-9][0]


def live_rows() -> Dict[Tuple[str, str], Dict[str, Any]]:
    out = {}
    for s, r in LIVE:
        p = G / f"live_{s}_{r}_h0.95.json"
        if p.exists():
            out[(s, r)] = load(p)
    return out


def parity_map_text(v: Dict[str, Any]) -> str:
    if v.get("map_equals_stored_after"):
        return "yes"
    fp = G / "flips.json"
    if fp.exists():
        m = load(fp)["parity_map"]
        return (
            f"no: {m['voxels_differing']} changed voxels differ by at most {m['max_abs_diff']:.1e}"
            f" (matrix entries); {m['differing_among_unchanged']} unchanged voxels differ"
        )
    return "no"


def pct(x: float) -> str:
    return f"{100 * x:.1f}%"


def interval(r: Dict[str, Any]) -> str:
    return f"{100 * r['wilson_lo']:.2f} to {100 * r['wilson_hi']:.2f}%"


def tables() -> Tuple[Dict[str, str], Dict[str, Any]]:
    R = rows_replay()
    L = live_rows()
    T: Dict[str, str] = {}
    S: Dict[str, Any] = {}

    rr = []
    for s, run in ENTRIES:
        if (s, run) not in R:
            continue
        r, x = R[(s, run)], at(R[(s, run)])
        rr.append(
            dict(
                sample=NAME[s],
                run=run,
                disclosed="yes" if s != "held" else "no (test)",
                full=r["pass_evals_full"],
                kept=x["pass_evals"],
                cut=pct(x["cut"]),
                skipped=f"{x['voxels_skipped']} of {r['boundary_voxels']}",
                wf=r["full"]["after"]["wrong"],
                wr=f"{x['wrong']} ({interval(x)})",
                fb=f"{x['fixed_vs_full']} / {x['broken_vs_full']}",
                mf=f"{r['full']['after']['median_err_deg']:.5f}",
                mr=f"{x['median_err_deg']:.5f}",
                mh=f"{r['max_own_hit_wrong_before']:.3f}",
                rule="yes" if x["rule"]["passes"] else "no",
            )
        )
    T["replay"] = markdown_table(
        rr,
        [
            "sample",
            "run",
            "disclosed",
            "full",
            "kept",
            "cut",
            "skipped",
            "wf",
            "wr",
            "fb",
            "mf",
            "mr",
            "mh",
            "rule",
        ],
        formats={"full": ",d", "kept": ",d"},
    ).replace(
        "| sample | run | disclosed | full | kept | cut | skipped | wf | wr | fb | mf | mr | mh "
        "| rule |",
        "| sample | run | disclosed | pass evaluations, full | pass evaluations, skip | cut | "
        "voxels skipped | wrong, full pass | wrong, skip (95% interval) | fixed / broken | "
        "median error, full (deg) | median error, skip (deg) | highest own hit ratio of a voxel "
        "wrong before the pass | rule parts of the pass met |",
    )

    ll = []
    for s, run in LIVE:
        if (s, run) not in L:
            continue
        v = L[(s, run)]
        ll.append(
            dict(
                sample=NAME[s],
                run=run,
                live=v["live_pass_evals"],
                replay=v["replay_pass_evals"],
                tries=f"{v['live_tries']} (replay {v['replay_tries']})",
                skipped=f"{v['live_skipped']} (replay {v['replay_skipped']})",
                wrong=f"{v['live']['wrong']} ({interval(v['live'])})",
                wrep=v["wrong_replay"],
                wfull=v["wrong_full"],
                fbf=f"{v['fixed_vs_full']} / {v['broken_vs_full']}",
                fbs=f"{v['fixed_vs_s8']} / {v['broken_vs_s8']}",
                med=f"{v['live']['median_err_deg']:.5f}",
                total=v["total_with_pass"],
                fewer=pct(1 - v["total_with_pass"] / v["s8_total"]),
            )
        )
    T["live"] = markdown_table(
        ll,
        ["sample", "run", "live", "replay", "tries", "skipped", "wrong", "wrep", "wfull"]
        + ["fbf", "fbs", "med", "total", "fewer"],
        formats={"live": ",d", "replay": ",d", "total": ",d"},
    ).replace(
        "| sample | run | live | replay | tries | skipped | wrong | wrep | wfull | fbf | fbs | med "
        "| total | fewer |",
        "| sample | run | pass evaluations, live | pass evaluations, replay | tries "
        "| voxels skipped | wrong, live (95% interval) | wrong, replay | wrong, full pass "
        "| fixed / broken vs the "
        "full pass | fixed / broken vs S8 | median error, live (deg) | run total, live | fewer "
        "than S8 |",
    )

    rl = []
    for s, run in LIVE:
        if (s, run) not in L:
            continue
        v = L[(s, run)]
        u = v["rule"]
        yn = lambda b: "met" if b else "not met"  # noqa: E731
        dm_full = v["live"]["median_err_deg"] - v["median_full"]
        dm_s8 = v["live"]["median_err_deg"] - v["s8_median"]
        rl.append(
            dict(
                sample=NAME[s],
                run=run,
                cut=f"{yn(u['pass_cut'])} ({pct(v['cut'])})",
                loss=f"{yn(u['pass_loss'])} ({v['broken_vs_full'] - v['fixed_vs_full']} <= "
                f"{u['allowed_loss']})",
                med=f"{yn(u['pass_median'])} ({dm_full:+.5f})",
                s8e=yn(u["pass_s8_evals"]),
                s8w=f"{yn(u['pass_s8_wrong'])} ({v['broken_vs_s8'] - v['fixed_vs_s8']} <= 5)",
                s8m=f"{yn(u['pass_s8_median'])} ({dm_s8:+.5f})",
                passes="yes" if u["passes"] else "no",
            )
        )
    T["rule_live"] = markdown_table(
        rl,
        ["sample", "run", "cut", "loss", "med", "s8e", "s8w", "s8m", "passes"],
    ).replace(
        "| sample | run | cut | loss | med | s8e | s8w | s8m | passes |",
        "| sample | run | pass cut at least 50% | broken minus fixed vs the full pass within the "
        "allowance | median error at most 0.001 deg higher than the full pass | at least 20% "
        "fewer evaluations than S8 | broken minus fixed vs S8 at most 5 | median error at most "
        "0.005 deg above S8 | the skip passes |",
    )

    # replay-level rule (disclosed samples, before any library change)
    rp = []
    for s, run in ENTRIES:
        if (s, run) not in R:
            continue
        u = at(R[(s, run)])["rule"]
        rp.append(
            dict(
                sample=NAME[s],
                run=run,
                cut="met" if u["pass_cut"] else "not met",
                loss="met" if u["pass_loss"] else "not met",
                med="met" if u["pass_median"] else "not met",
                s8=(
                    "met"
                    if (u["pass_s8_evals"] and u["pass_s8_wrong"] and u["pass_s8_median"])
                    else "not met"
                ),
                passes="yes" if u["passes"] else "no",
            )
        )
    T["rule_replay"] = markdown_table(
        rp,
        ["sample", "run", "cut", "loss", "med", "s8", "passes"],
    ).replace(
        "| sample | run | cut | loss | med | s8 | passes |",
        "| sample | run | pass cut at least 50% | broken minus fixed within the allowance | median "
        "within 0.001 deg | house rules vs S8 | the skip passes (replay) |",
    )

    # what-if curve
    hs = [x["h"] for x in R[("mid", "S11n")]["by_h"]]
    cols = ["h"] + [f"{s} {r}" for s, r in ENTRIES if (s, r) in R]
    cv = []
    for i, h in enumerate(hs):
        row: Dict[str, Any] = {"h": f"{h:.2f}"}
        for s, r in ENTRIES:
            if (s, r) in R:
                x = R[(s, r)]["by_h"][i]
                row[f"{s} {r}"] = f"{100 * x['cut']:.1f}% / {x['net_gain_lost']}"
        cv.append(row)
    T["curve"] = markdown_table(cv, cols)

    # own hit ratio of the voxels the full pass fixed
    fh = []
    for s, run in ENTRIES:
        if (s, run) not in R:
            continue
        r = R[(s, run)]
        o = r["own_hit_of_fixed"]
        fh.append(
            dict(
                sample=NAME[s],
                run=run,
                n=len(o),
                lo=f"{min(o):.3f}",
                hi=f"{max(o):.3f}",
                ge=sum(1 for x in o if x >= GR.H_S),
                wrong=f"{r['max_own_hit_wrong_before']:.3f} / {r['max_own_hit_wrong_after']:.3f}",
            )
        )
    T["fixed_hits"] = markdown_table(
        fh,
        ["sample", "run", "n", "lo", "hi", "ge", "wrong"],
    ).replace(
        "| sample | run | n | lo | hi | ge | wrong |",
        "| sample | run | voxels the full pass fixed | lowest own hit ratio | highest own hit "
        "ratio | of them with own hit ratio at least 0.95 | highest own hit ratio of a voxel wrong "
        "before / after the pass |",
    )

    # predictions
    pred = predictions(R, L)
    T["predictions"] = markdown_table(pred, ["no", "sample", "observed", "verdict"]).replace(
        "| no | sample | observed | verdict |",
        "| prediction | sample | observed | verdict |",
    )
    S["predictions"] = pred

    # parity
    pp = G / "live_mid_S11n_h0.json"
    if pp.exists():
        v = load(pp)
        T["parity"] = markdown_table(
            [
                dict(
                    tries=f"{v['live_tries']} of {v['replay_tries']}",
                    evals=f"{v['live_pass_evals']:,} of {v['stored_full_pass_evals']:,}",
                    diffs=" / ".join(
                        str(v.get(f"differ_{k}")) for k in ("voxel", "cand", "evals", "cost")
                    )
                    + " / "
                    + " / ".join(str(v.get(f"differ_{k}")) for k in ("hit", "conf", "end")),
                    mp=parity_map_text(v),
                    ch=f"{v['voxels_changed']} of {v['stored_voxels_changed']}",
                    ss="yes" if v.get("same_changed_voxels") else "no",
                )
            ],
            ["tries", "evals", "diffs", "mp", "ch", "ss"],
        ).replace(
            "| tries | evals | diffs | mp | ch | ss |",
            "| tries, live and stored | pass evaluations, live and stored | tries that differ: "
            "voxel / candidate / evaluations / cost / hit ratio / confidence / end | map equals "
            "the stored map | voxels changed, live and stored | same voxels changed |",
        )
    fp = G / "flips.json"
    if fp.exists():
        fl = load(fp)
        fr_ = []
        for s, run in LIVE:
            f = fl.get(f"{s} {run}")
            if f is None:
                continue
            bo, fo = f["broken_own_hit"], f["fixed_own_hit"]
            fr_.append(
                dict(
                    sample=NAME[s],
                    run=run,
                    kept=f["boundary_kept"],
                    sk=f["skipped"],
                    dsk=f["end_differs_skipped"],
                    dk=f["end_differs_kept"],
                    fb=f"{f['fixed_vs_full']} / {f['broken_vs_full']}",
                    fl=f["flips_skipped"],
                    rg=f"{min(bo + fo):.3f} to {max(bo + fo):.3f}" if bo + fo else "n/a",
                )
            )
        T["flips"] = markdown_table(
            fr_, ["sample", "run", "kept", "sk", "dsk", "dk", "fb", "fl", "rg"]
        ).replace(
            "| sample | run | kept | sk | dsk | dk | fb | fl | rg |",
            "| sample | run | boundary voxels with tries | voxels skipped | skipped voxels whose "
            "end the full pass had changed | kept voxels whose end differs from the full pass "
            "| fixed / broken vs the full pass | of those voxels, skipped ones | own hit ratio "
            "of those voxels |",
        )
        S["flips"] = fl
    S["rows_replay"] = {f"{s} {r}": at(R[(s, r)]) for s, r in ENTRIES if (s, r) in R}
    S["live"] = {f"{s} {r}": L[(s, r)] for s, r in LIVE if (s, r) in L}
    return T, S


def predictions(R: Dict[Tuple[str, str], Dict[str, Any]], L: Dict[Tuple[str, str], Any]) -> List:
    fr = load(G / "frozen.json")["predictions"]
    out: List[Dict[str, Any]] = []
    h = R.get(("held", "S11v2"))
    if h is not None:
        x = at(h)
        ok1 = 55 <= 100 * x["cut"] <= 85 and x["fixed_vs_full"] <= 0 and x["broken_vs_full"] <= 2
        out.append(
            dict(
                no=1,
                text=fr[0],
                sample="held (replay)",
                observed=f"cut {pct(x['cut'])}; fixed / broken {x['fixed_vs_full']} / "
                f"{x['broken_vs_full']}",
                verdict="met" if ok1 else "not met",
            )
        )
    lv = L.get(("held", "S11v2"))
    if lv is not None:
        fewer = 1 - lv["total_with_pass"] / lv["total_with_full_pass"]
        out.append(
            dict(
                no=2,
                text=fr[1],
                sample="held (live)",
                observed=f"{pct(fewer)} fewer ({lv['total_with_pass']:,} against "
                f"{lv['total_with_full_pass']:,}); the pass was {lv['stored_full_pass_evals']:,}",
                verdict="met" if 5 <= 100 * fewer <= 12 else "not met",
            )
        )
    for s, run in LIVE:
        v = L.get((s, run))
        if v is None:
            continue
        out.append(
            dict(
                no=3,
                text=fr[2],
                sample=f"{s} {run}",
                observed=f"live {v['live_pass_evals']:,}, replay {v['replay_pass_evals']:,}; "
                f"{v['live_tries']} tries at "
                f"{v['live_pass_evals'] / max(1, v['live_tries']):.1f} each",
                verdict="met" if v["evals_equal_replay"] else "not met",
            )
        )
    for s, run in LIVE:
        v = L.get((s, run))
        if v is None:
            continue
        d = v["live"]["wrong"] - v["wrong_replay"]
        out.append(
            dict(
                no=4,
                text=fr[3],
                sample=f"{s} {run}",
                observed=f"live {v['live']['wrong']}, replay {v['wrong_replay']} "
                f"(difference {d:+d})",
                verdict="met" if abs(d) <= 2 else "not met",
            )
        )
    for s, run in LIVE:
        v = L.get((s, run))
        if v is None:
            continue
        d = v["live"]["median_err_deg"] - v["median_full"]
        out.append(
            dict(
                no=5,
                text=fr[4],
                sample=f"{s} {run}",
                observed=f"{d:+.5f} degree against the full pass",
                verdict="met" if d <= 0.0005 + 1e-12 else "not met",
            )
        )
    return out


def plot(
    R: Dict[Tuple[str, str], Dict[str, Any]],
    L: Dict[Tuple[str, str], Any],
    L2: Optional[Dict[Tuple[str, str], Any]] = None,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.5, 5))
    for s, run in ENTRIES:
        if (s, run) not in R:
            continue
        by = R[(s, run)]["by_h"]
        x = [100 * b["cut"] for b in by]
        y = [b["net_gain_lost"] for b in by]
        col = COLORS[s]
        (ln,) = ax.plot(
            x,
            y,
            marker=".",
            color=col,
            lw=2.5 if s == "held" else 1.2,
            ls="--" if run == "S11" else "-",
            label=f"{NAME[s]} {run}" + (" (test)" if s == "held" else ""),
        )
        b95 = at(R[(s, run)])
        ax.plot(100 * b95["cut"], b95["net_gain_lost"], marker="o", ms=9, mfc="none", mec="k")
        v = L.get((s, run))
        if v is not None:
            ax.plot(
                100 * v["cut"],
                v["live"]["wrong"] - v["wrong_full"],
                marker="x",
                ms=10,
                color=ln.get_color(),
                ls="none",
            )
        v2 = (L2 or {}).get((s, run))
        if v2 is not None:
            ax.plot(
                100 * v2["cut"],
                v2["live"]["wrong"] - v2["wrong_full"],
                marker="+",
                ms=14,
                mew=2,
                color="k",
                ls="none",
            )
    ax.axvline(50, color="grey", ls=":", lw=1)
    ax.set_xlabel("share of the boundary pass cut (percent of its evaluations)")
    ax.set_ylabel("net gain lost (wrong voxels after the skip minus after the full pass)")
    ax.set_title(
        "What-if curve over h from 0.99 (left) to 0.90 (right); circle = 0.95 (frozen)\n"
        "x = live pass, variant 1 (one draw stream); + = live pass, variant 2 (aligned draws,\n"
        "decided after the first live result)"
    )
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(G / "ar5_gate_curve.png", dpi=130)


def pareto(
    R: Dict[Tuple[str, str], Dict[str, Any]],
    L: Dict[Tuple[str, str], Any],
    L2: Optional[Dict[Tuple[str, str], Any]] = None,
) -> None:
    rows = []
    sid = {"big": "big_24570", "grad": "gradient_24570", "mid": "mid_24570", "held": "held_24570"}
    for s, run in ENTRIES:
        if (s, run) not in R:
            continue
        r = R[(s, run)]
        full_total = r["bfs_evals"] + r["pass_evals_full"]
        base = dict(
            sample=sid[s],
            step="AR5 step 7",
            family="BFS with indexed seeds",
            boundary_pass=True,
            n=GC.N,
            disclosed_sample=s != "held",
            source_json="icenine_py/benchmarks/autoresearch/ar5/gate/summary.json",
            source_page="icenine_py/docs/autoresearch/findings_ar5_gate.md",
        )
        rows.append(
            dict(
                base,
                method=f"run {run} with the full boundary pass (stored)",
                evals_total=full_total,
                wrong=r["full"]["after"]["wrong"],
                kind="stored",
            )
        )
        v = L.get((s, run))
        if v is not None:
            rows.append(
                dict(
                    base,
                    method=f"run {run} with BoundaryPassSkipHitRatio 0.95, variant 1 (live pass)",
                    evals_total=v["total_with_pass"],
                    wrong=v["live"]["wrong"],
                    kind="live",
                )
            )
        v2 = (L2 or {}).get((s, run))
        if v2 is not None:
            rows.append(
                dict(
                    base,
                    method=f"run {run} with BoundaryPassSkipHitRatio 0.95, variant 2 (aligned "
                    "draws; live pass)",
                    evals_total=v2["total_with_pass"],
                    wrong=v2["live"]["wrong"],
                    kind="live",
                    flag="variant 2, decided after the first live result, D-20261112-N5",
                )
            )
        if v is None:
            x = at(r)
            rows.append(
                dict(
                    base,
                    method=f"run {run} with the skip at 0.95 (replay of the stored pass)",
                    evals_total=x["total_with_skip"],
                    wrong=x["wrong"],
                    kind="replay",
                )
            )
    (G / "pareto_rows.json").write_text(json.dumps(rows, indent=1) + "\n")


V2 = "variant 2, aligned draws, decided after the first live result, D-20261112-N5"


def live_rows_v2() -> Dict[Tuple[str, str], Dict[str, Any]]:
    out = {}
    for s, r in LIVE:
        p = G / f"live_{s}_{r}_h0.95_v2.json"
        if p.exists():
            out[(s, r)] = load(p)
    return out


def yn(b: bool) -> str:
    return "met" if b else "not met"


def v2_tables(L1: Dict[Tuple[str, str], Any], L2: Dict[Tuple[str, str], Any]) -> Dict[str, Any]:
    """Tables of variant 2 (aligned draws; decided after the first live result) beside variant 1."""
    T: Dict[str, str] = {}
    rows, rr, vs, nz = [], [], [], []
    fl = load(G / "flips.json") if (G / "flips.json").exists() else {}
    for s, run in LIVE:
        v, v1 = L2.get((s, run)), L1.get((s, run))
        if v is None or v1 is None:
            continue
        u, u1 = v["rule"], v1["rule"]
        fd = v.get("first_difference")
        rows.append(
            dict(
                sample=NAME[s],
                run=run,
                ev=f"{v['live_pass_evals']:,}",
                ident=f"{v['kept_tries_identical']:,} of {v['kept_tries_total']:,}",
                first="none" if fd is None else f"voxel {fd['voxel']}, candidate {fd['cand']}",
                wrong=f"{v['live']['wrong']} ({interval(v['live'])})",
                fbf=f"{v['fixed_vs_full']} / {v['broken_vs_full']}",
                fbs=f"{v['fixed_vs_s8']} / {v['broken_vs_s8']}",
                med=f"{v['live']['median_err_deg']:.5f}",
                total=f"{v['total_with_pass']:,}",
                fewer=pct(1 - v["total_with_pass"] / v["s8_total"]),
            )
        )
        dm_full = v["live"]["median_err_deg"] - v["median_full"]
        dm_s8 = v["live"]["median_err_deg"] - v["s8_median"]
        rr.append(
            dict(
                sample=NAME[s],
                run=run,
                cut=f"{yn(u['pass_cut'])} ({pct(v['cut'])})",
                loss=f"{yn(u['pass_loss'])} ({v['broken_vs_full'] - v['fixed_vs_full']} <= "
                f"{u['allowed_loss']})",
                med=f"{yn(u['pass_median'])} ({dm_full:+.5f})",
                s8e=yn(u["pass_s8_evals"]),
                s8w=f"{yn(u['pass_s8_wrong'])} ({v['broken_vs_s8'] - v['fixed_vs_s8']} <= 5)",
                s8m=f"{yn(u['pass_s8_median'])} ({dm_s8:+.5f})",
                passes="yes" if u["passes"] else "no",
            )
        )
        vs.append(
            dict(
                sample=NAME[s],
                run=run,
                a="yes" if u1["passes"] else "no",
                b=f"{v1['fixed_vs_full']} / {v1['broken_vs_full']}",
                c="yes" if u["passes"] else "no",
                d=f"{v['fixed_vs_full']} / {v['broken_vs_full']}",
            )
        )
        f = fl.get(f"{s} {run}", {})
        nz.append(
            dict(
                sample=NAME[s],
                run=run,
                kept=f.get("boundary_kept", "n/a"),
                diff=f.get("end_differs_kept", "n/a"),
                fb=f"{v1['fixed_vs_full']} / {v1['broken_vs_full']}",
                wrong=f"{v1['live']['wrong']} against {v1['wrong_full']}",
            )
        )
    T["live_v2"] = markdown_table(
        rows,
        ["sample", "run", "ev", "ident", "first", "wrong", "fbf", "fbs", "med", "total"]
        + ["fewer"],
    ).replace(
        "| sample | run | ev | ident | first | wrong | fbf | fbs | med | total | fewer |",
        "| sample | run | pass evaluations, live (equal to the replay) "
        "| kept tries identical to the "
        "stored full pass (voxel, candidate, end, cost, hit ratio, confidence, evaluations) | "
        "first difference | wrong, live (95% interval) | fixed / broken vs the full pass | "
        "fixed / broken vs S8 | median error (deg) | run total | fewer than S8 |",
    )
    T["rule_live_v2"] = markdown_table(
        rr, ["sample", "run", "cut", "loss", "med", "s8e", "s8w", "s8m", "passes"]
    ).replace(
        "| sample | run | cut | loss | med | s8e | s8w | s8m | passes |",
        "| sample | run | pass cut at least 50% | broken minus fixed vs the full pass within the "
        "allowance | median error at most 0.001 deg higher than the full pass | at least 20% "
        "fewer evaluations than S8 | broken minus fixed vs S8 at most 5 | median error at most "
        "0.005 deg above S8 | the skip passes |",
    )
    T["variants"] = markdown_table(vs, ["sample", "run", "a", "b", "c", "d"]).replace(
        "| sample | run | a | b | c | d |",
        "| sample | run | variant 1 (one draw stream): passes the rule | variant 1: fixed / broken "
        "vs the full pass | variant 2 (aligned draws): passes the rule | variant 2: fixed / "
        "broken vs the full pass |",
    )
    T["noise"] = markdown_table(nz, ["sample", "run", "kept", "diff", "fb", "wrong"]).replace(
        "| sample | run | kept | diff | fb | wrong |",
        "| sample | run | boundary voxels kept (with tries) | kept voxels whose end differs from "
        "the full pass | fixed / broken vs the full pass | wrong, live against full pass |",
    )
    pred = predictions(rows_replay(), L2)
    T["predictions_v2"] = markdown_table(pred, ["no", "sample", "observed", "verdict"]).replace(
        "| no | sample | observed | verdict |", "| prediction | sample | observed | verdict |"
    )
    return dict(tables=T, predictions=pred, live={f"{s} {r}": v for (s, r), v in L2.items()})


def main() -> None:
    T, S = tables()
    R, L = rows_replay(), live_rows()
    L2 = live_rows_v2()
    if L2:
        d = v2_tables(L, L2)
        T.update(d["tables"])
        S["variant_2"] = dict(label=V2, predictions=d["predictions"], live=d["live"])
    write_tables(G / "tables.md", T)
    plot(R, L, L2)
    pareto(R, L, L2)
    (G / "summary.json").write_text(json.dumps(S, indent=1, default=float) + "\n")
    print("wrote tables.md, summary.json, ar5_gate_curve.png, pareto_rows.json")
    for p in S["predictions"]:
        print(p["no"], p["sample"], p["verdict"], p["observed"])


if __name__ == "__main__":
    main()
