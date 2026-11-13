"""AR5 step 9, phase A step 6: validity, plot, tables, pareto rows, predictions, verdict.

  uv run python scripts/autoresearch/ar5/rescue_summary.py

Reads ``rescue/stats_<sample>.json`` (``rescue_stats.py``) and ``rescue/gate.json`` and writes
``rescue/summary.json``, ``rescue/tables.md`` (doc-table blocks, a reading line before each table),
``rescue/pareto_rows.json`` (the format of ``ar5/order/pareto_rows.json``) and
``rescue/ar5_rescue_evals_wrong.png`` (the plot is made first). No p-values.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "common"))
RES = HERE.parents[2] / "benchmarks" / "autoresearch" / "ar5" / "rescue"
SAMPLES = ["held", "mid", "grad", "big"]
N = 24570

import doc_tables as DT  # noqa: E402
import rescue_judge as RJ  # noqa: E402
import stats as ST  # noqa: E402

STORED = {  # the plan's stop-rule totals (arm A, with the pass) and fallback evaluations
    "grad": (12_189_777, 3_726_723),
    "mid": (15_645_472, 4_007_241),
    "big": (3_807_747, 0),
}


def f(x: float, d: int = 5) -> str:
    return f"{x:.{d}f}"


def n(x: int) -> str:
    return f"{x:,}"


def wil(k: int, m: int) -> str:
    """``k / m (lo% to hi%)`` with the Wilson 95% interval (``-`` when ``m`` is 0)."""
    if m == 0:
        return "-"
    lo, hi = ST.wilson(k, m)
    return f"{k} / {m} ({100 * lo:.1f}% to {100 * hi:.1f}%)"


def validity(st: Dict[str, Any]) -> Dict[str, Any]:
    v: Dict[str, Any] = {}
    for s in SAMPLES:
        a, g = st[s]["arms"]["A"], st[s]["arms"]["G"]
        ok = dict(
            runs_complete=not a["unfinished"] and not g["unfinished"],
            seed_steps_ok=g["seed_steps"] > 0,
        )
        if s in STORED:
            ok["stored_total_reproduces"] = a["stages"]["total"] == STORED[s][0]
            ok["stored_fallback_evals_reproduce"] = a["stages"]["grid_fallback"] == STORED[s][1]
        if s == "big":
            # no seed test fails on big, so the rescue never runs there: this shows only that the
            # idle wrapper adds no draw or side effect, not that the draws are isolated
            ok["big_idle_path_g_equals_a_evals_and_map"] = bool(
                st[s]["g_vs_a"]["same_map"] and st[s]["g_vs_a"]["same_evals"]
            )
        v[s] = ok
    v["all_ok"] = all(all(x.values()) for k, x in v.items() if k in SAMPLES)
    return v


def isolation(iso: Dict[str, Any]) -> Dict[str, Any]:
    """The check of draw isolation (review fix 6): voxels that differ between A and G."""
    out: Dict[str, Any] = {}
    for s in SAMPLES:
        b, a = iso[s]["differ_before"], iso[s]["differ_after"]
        out[s] = dict(
            rescued=iso[s]["rescued"]["n"],
            differ_before=b["n"],
            differ_before_others=b["others"],
            differ_before_others_in_a_rescued_grain=b["others_in_a_rescued_grain"],
            differ_before_others_max_deg=b["others_max_deg"],
            differ_before_is_the_rescued_set=b["equals_rescued_set"],
            differ_after=a["n"],
            differ_after_others=a["others"],
        )
    out["note"] = (
        "With aligned draws only the rescued voxels and voxels grown from them may change. "
        "mid: exactly the rescued voxels differ. grad, held: the other differing voxels are "
        "mostly in a rescued voxel's grain (consistent with local growth from a changed seed; "
        "the growth link of the rest is not traced). big has no failed seed step and is not a "
        "test of isolation."
    )
    return out


def plot(st: Dict[str, Any]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    col = dict(held="tab:red", mid="tab:blue", grad="tab:green", big="tab:orange")
    for ax, tag, title in zip(axes, ("before", "after"), ("before the pass", "after the pass")):
        for s in SAMPLES:
            a, g = st[s]["arms"]["A"], st[s]["arms"]["G"]
            key = "before_pass" if tag == "before" else "total"
            pts = [
                ("S8", "^", a["stages_s8"][key], a[f"s8_{tag}"]),
                ("A", "o", a["stages"][key], a[tag]),
                ("G", "s", g["stages"][key], g[tag]),
            ]
            for arm, mk, x, side in pts:
                k = side["wrong"]
                lo, hi = N * side["wilson_lo"], N * side["wilson_hi"]
                ax.errorbar(
                    x / 1e6,
                    k,
                    yerr=[[k - lo], [hi - k]],
                    fmt=mk,
                    color=col[s],
                    mfc="none" if arm == "S8" else col[s],
                    ms=6,
                    capsize=2,
                    lw=0.8,
                    label=f"{s} {arm}" if tag == "before" else None,
                )
        ax.set_xscale("log")
        ax.set_xlabel(
            "run evaluations (millions, log%s)" % (", with the pass" if tag == "after" else "")
        )
        ax.set_ylabel("wrong voxels (more than 1 degree), Wilson 95%")
        ax.set_title(title)
        ax.grid(alpha=0.3)
    fig.legend(loc="lower center", fontsize=7, ncol=6, frameon=False)
    fig.suptitle(
        "Stored S8 (triangle), arm A (circle), arm G = A with rescue G1 (square)", fontsize=10
    )
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    fig.savefig(RES / "ar5_rescue_evals_wrong.png", dpi=110)


def pareto(st: Dict[str, Any]) -> None:
    rows: List[Dict[str, Any]] = []
    for s in SAMPLES:
        for arm, order in (("A", "interleaved"), ("G", "interleaved + rescue G1")):
            d = st[s]["arms"][arm]
            for after in (False, True):
                side = d["after" if after else "before"]
                total = d["stages"]["total" if after else "before_pass"]
                rows.append(
                    dict(
                        sample=s,
                        arm=arm,
                        order=order,
                        boundary_pass=after,
                        evals_total=total,
                        evals_per_voxel=round(total / N, 2),
                        wrong=side["wrong"],
                        n=N,
                        interval_pct=[
                            round(100 * side["wilson_lo"], 2),
                            round(100 * side["wilson_hi"], 2),
                        ],
                        source="icenine_py/benchmarks/autoresearch/ar5/rescue/stats_%s.json" % s,
                    )
                )
    (RES / "pareto_rows.json").write_text(json.dumps(rows, indent=1) + "\n")


def main() -> None:
    st = {s: json.loads((RES / f"stats_{s}.json").read_text()) for s in SAMPLES}
    iso = {s: json.loads((RES / f"isolation_{s}.json").read_text()) for s in SAMPLES}
    for s in SAMPLES:  # each try is counted as 1,001 evaluations (checked, not assumed)
        z = st[s]["rescue"]
        assert z["rescue_evals"] == 1001 * z["tries"], (s, z)
    gate = json.loads((RES / "gate.json").read_text())["all"]
    val = validity(st)
    plot(st)
    pareto(st)
    T: Dict[str, str] = {}
    # ---- table 1: evaluations
    r = [
        (
            "| sample | A total | G total | fewer | A seed step | G seed step | G rescue | A grid |"
            " G grid | neighbour+revisit (G) | pass A / G |"
        ),
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for s in SAMPLES:
        a, g = st[s]["arms"]["A"], st[s]["arms"]["G"]
        sp = g["stages_split"]
        cut = RJ.pct_fewer(g["stages"]["total"], a["stages"]["total"])
        r.append(
            f"| {s} | {n(a['stages']['total'])} | {n(g['stages']['total'])} | {cut:.1f}% | "
            f"{n(a['stages']['seed_step'])} | {n(sp['seed_step'])} | {n(sp['rescue'])} | "
            f"{n(a['stages']['grid_fallback'])} | {n(sp['grid_fallback'])} | "
            f"{n(sp['neighbour'] + sp['revisit'])} | {n(a['stages']['pass_'])} / {n(sp['pass_'])} |"
        )
    T["evals"] = (
        "Reading line: evaluations by stage, run plus boundary pass, arm A (stored step 8) and arm "
        "G; `fewer` is G's total against A's. Wall times are not given (contended).\n\n"
        + "\n".join(r)
    )
    # ---- table 2: rescues
    r = [
        (
            "| sample | failed seed steps (A / G) | rescued / tried (Wilson 95%) | tries | grid"
            " fallbacks left | rescued hit ratio min | rescued error median / max | wrong rescued"
            " / rescued (Wilson 95%) |"
        ),
        "|---|---|---|---|---|---|---|---|",
    ]
    for s in SAMPLES:
        a = st[s]["arms"]["A"]
        z = st[s]["rescue"]
        em = st[s]["rescued_err_median"]
        hit = f"{st[s]['rescued_hit_min']:.3f}" if em is not None else "-"
        err = f"{em:.4f} / {st[s]['rescued_err_max']:.4f}" if em is not None else "-"
        r.append(
            f"| {s} | {a['grid_fallbacks']} / {z['rescues_attempted']} |"
            f" {wil(z['rescued'], z['rescues_attempted'])} | {z['tries']} |"
            f" {z['grid_fallbacks_left']} | {hit} | {err} |"
            f" {wil(st[s]['wrong_rescued'], z['rescued'])} |"
        )
    k_all = sum(st[s]["rescue"]["rescued"] for s in SAMPLES)
    w_all = sum(st[s]["wrong_rescued"] for s in SAMPLES)
    r.append(f"| all | - | - | - | - | - | - | {wil(w_all, k_all)} |")
    T["rescues"] = (
        "Reading line: what rule G1 did at the failed seed steps of arm G (errors in degrees"
        " against the truth, scoring only; wrong means over 1 degree). A try is counted as 1,001"
        " evaluations; the one overlap evaluation per try for the hit ratio is not counted (as in"
        " the seed step).\n\n" + "\n".join(r)
    )
    # ---- table 3: accuracy
    r = [
        (
            "| sample | arm | wrong before pass (Wilson 95%) | wrong after pass (Wilson 95%) |"
            " median error after pass |"
        ),
        "|---|---|---|---|---|",
    ]
    for s in SAMPLES:
        a, g = st[s]["arms"]["A"], st[s]["arms"]["G"]
        for name, side_b, side_a in (
            ("S8", a["s8_before"], a["s8_after"]),
            ("A", a["before"], a["after"]),
            ("G", g["before"], g["after"]),
        ):
            r.append(
                f"| {s} | {name} | {side_b['wrong']} ({100*side_b['wilson_lo']:.2f}% to"
                f" {100*side_b['wilson_hi']:.2f}%) | {side_a['wrong']}"
                f" ({100*side_a['wilson_lo']:.2f}% to {100*side_a['wilson_hi']:.2f}%) |"
                f" {f(side_a['median_err_deg'])} |"
            )
    T["accuracy"] = (
        f"Reading line: wrong voxels out of {N:,} (more than 1 degree from the truth), before and "
        "after the boundary pass, and the median error after the pass in degrees.\n\n"
        + "\n".join(r)
    )
    # ---- table 4: paired
    r = [
        (
            "| sample | pair | before pass fixed / broken | after pass fixed / broken | both wrong"
            " after | seed step's own end differs (pairing check) |"
        ),
        "|---|---|---|---|---|---|",
    ]
    for s in SAMPLES:
        a, g = st[s]["arms"]["A"], st[s]["arms"]["G"]
        gv = st[s]["g_vs_a"]
        r.append(
            f"| {s} | G vs A | {gv['voxels']['before']['fixed']} /"
            f" {gv['voxels']['before']['broken']} | {gv['voxels']['after']['fixed']} /"
            f" {gv['voxels']['after']['broken']} | {gv['voxels']['after']['both_wrong']} |"
            f" {gv['seed_steps']['taken_differs']} of {gv['seed_steps']['matched']} |"
        )
        r.append(
            f"| {s} | G vs S8 | {g['paired_before']['fixed']} / {g['paired_before']['broken']} | "
            f"{g['paired_after']['fixed']} / {g['paired_after']['broken']} | "
            f"{g['paired_after']['both_wrong']} | - |"
        )
        r.append(
            f"| {s} | A vs S8 | {a['paired_before']['fixed']} / {a['paired_before']['broken']} | "
            f"{a['paired_after']['fixed']} / {a['paired_after']['broken']} | "
            f"{a['paired_after']['both_wrong']} | - |"
        )
    T["paired"] = (
        "Reading line: the base is the arm named after `vs`; fixed = wrong in the base and right "
        "in the first-named arm, broken = the reverse, counted by voxel. The last column is 0 "
        "by construction: the rescue runs after the seed step, whose own fits are unchanged; it "
        "checks the pairing, not whether a voxel's end changed (see the table of changed "
        "ends).\n\n" + "\n".join(r)
    )
    # ---- table 4b: draw isolation and changed ends (review fixes 5 and 6)
    r = [
        (
            "| sample | rescued | voxels whose end changed, before pass (rescued + others) |"
            " others in a rescued voxel's grain | others, largest change | changed after pass |"
            " G vs A at the rescued voxels: median / largest | largest where A was right |"
            " A wrong (over 1 degree) at a rescued voxel |"
        ),
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for s in SAMPLES:
        b, af, rz = iso[s]["differ_before"], iso[s]["differ_after"], iso[s]["rescued"]
        if rz["n"]:
            med = f"{f(rz['median_deg'], 4)} / {f(rz['max_deg'], 4)}"
            mr = f(rz["max_deg_where_a_right"], 4)
            aw = (
                f"{len(rz['a_wrong'])} (A off "
                + ", ".join(f"{w['err_a']:.0f}" for w in rz["a_wrong"])
                + " degrees; G's end "
                + ("right" if all(w["err_g"] <= 1.0 for w in rz["a_wrong"]) else "NOT right")
                + ")"
                if rz["a_wrong"]
                else "0"
            )
        else:
            med = mr = aw = "-"
        r.append(
            f"| {s} | {rz['n']} | {b['n']} ({b['rescued_in']} + {b['others']}) |"
            f" {b['others_in_a_rescued_grain']} | {f(b['others_max_deg'], 4)} | {af['n']} |"
            f" {med} | {mr} | {aw} |"
        )
    T["isolation"] = (
        "Reading line: per sample, the voxels whose stored orientation differs between arm A and "
        "arm G (angles in degrees). With aligned draws only the rescued voxels and voxels grown "
        "from them may change; this is the check of draw isolation (big has no failed seed "
        "step, so its G = A tests nothing about it).\n\n" + "\n".join(r)
    )
    # ---- rules
    r = [
        (
            "| sample | (1) fallbacks / total | (2) net vs A | (3) median vs A | (4) wrong rescued"
            " | (5) house vs S8 (A passes it?) | G passes |"
        ),
        "|---|---|---|---|---|---|---|",
    ]
    for s in SAMPLES:
        ru, g, a = st[s]["rules"], st[s]["arms"]["G"], st[s]["arms"]["A"]
        gv = st[s]["g_vs_a"]["voxels"]["after"]
        r.append(
            f"| {s} | {'yes' if ru['r1'] else 'NO'}: grid {a['grid_fallbacks']} to"
            f" {st[s]['rescue']['grid_fallbacks_left']}, total {ru['total_cut_pct']:.1f}% fewer |"
            f" {'yes' if ru['r2'] else 'NO'}: {gv['broken'] - gv['fixed']} |"
            f" {'yes' if ru['r3'] else 'NO'}:"
            f" {g['after']['median_err_deg'] - a['after']['median_err_deg']:+.5f} |"
            f" {'yes' if ru['r4'] else 'NO'}: {st[s]['wrong_rescued']} |"
            f" {'yes' if ru['r5'] else 'NO'} (A: {'no' if ru['a_fails_5'] else 'yes'}) |"
            f" {'PASS' if ru['passes'] else 'FAIL'} |"
        )
    T["rules"] = (
        "Reading line: the plan's five frozen rules for arm G against arm A (and the house rule "
        "against S8), per sample; a rule that is not met says NO.\n\n" + "\n".join(r)
    )
    # ---- predictions
    h, m, gr, b = (st[s] for s in ("held", "mid", "grad", "big"))
    fb_cut_h = RJ.pct_fewer(h["rescue"]["grid_fallbacks_left"], h["arms"]["A"]["grid_fallbacks"])
    cut = {s: st[s]["rules"]["total_cut_pct"] for s in SAMPLES}
    nets = {
        s: (
            st[s]["g_vs_a"]["voxels"]["after"]["broken"]
            - st[s]["g_vs_a"]["voxels"]["after"]["fixed"]
        )
        for s in SAMPLES
    }
    meds = {
        s: (
            st[s]["arms"]["G"]["after"]["median_err_deg"]
            - st[s]["arms"]["A"]["after"]["median_err_deg"]
        )
        for s in SAMPLES
    }
    p3 = 40 <= fb_cut_h <= 85 and 20 <= cut["held"] <= 45
    p4 = 12 <= cut["mid"] <= 21 and 15 <= cut["grad"] <= 28 and st["big"]["rules"]["r1"]
    p5 = all(-3 <= nets[s] <= 5 for s in SAMPLES) and all(
        abs(meds[s]) <= 0.001 + 1e-9 for s in SAMPLES
    )
    p2n = gate
    preds = {
        "1": dict(
            met=True,
            text="cost check: grid end cheaper than every kept end in 23 of 23 (needed 21)",
        ),
        "2": dict(
            met=False,
            text=(
                f"replay: rescued {p2n['rescued']} of 23 (12 to 21 met), all right (met), median"
                f" error {p2n['rescued_err_median']:.4f} (met), mean tries {p2n['tries_mean']:.2f}"
                " (at most 2.5 NOT met)"
            ),
        ),
        "3": dict(
            met=bool(p3),
            text=(
                f"held: grid fallbacks {h['arms']['A']['grid_fallbacks']} to"
                f" {h['rescue']['grid_fallbacks_left']} ({fb_cut_h:.1f}% fewer, predicted 40 to"
                f" 85), total {cut['held']:.1f}% fewer (predicted 20 to 45)"
            ),
        ),
        "4": dict(
            met=bool(p4),
            text=(
                f"mid {cut['mid']:.1f}% fewer (12 to 21), grad {cut['grad']:.1f}% (15 to 28), "
                f"big identical: {st['big']['rules']['r1']}"
            ),
        ),
        "5": dict(
            met=bool(p5),
            text="broken minus fixed against A (after pass) "
            + ", ".join(f"{s} {nets[s]}" for s in SAMPLES)
            + " (-3 to 5); median error moves "
            + ", ".join(f"{s} {meds[s]:+.5f}" for s in SAMPLES)
            + " (at most 0.001); held's -4 is outside the range, on the favourable side",
        ),
    }
    r = ["| prediction | met | numbers |", "|---|---|---|"]
    for k, v in preds.items():
        r.append(f"| {k} | {'met' if v['met'] else 'NOT met'} | {v['text']} |")
    T["predictions"] = "Reading line: the five predictions written before any run.\n\n" + "\n".join(
        r
    )
    per = {s: bool(st[s]["rules"]["passes"]) for s in SAMPLES}
    verdict = RJ.step_verdict(per)
    head = (
        "Reading line: the plan's verdict: the step passes if G passes on held and on all "
        "disclosed samples; if it fails on held there is no key.\n\n"
        "| sample | G passes |\n|---|---|\n"
    )
    body = "\n".join(f"| {s} | {'PASS' if per[s] else 'FAIL'} |" for s in SAMPLES)
    tail = f"\n\nStep verdict: {'PASS' if verdict else 'FAIL'}"
    tail += " (held fails: no key)" if not per["held"] else ""
    T["verdict"] = head + body + tail
    DT.write_tables(RES / "tables.md", T)
    (RES / "summary.json").write_text(
        json.dumps(
            dict(
                validity=val,
                draw_isolation=isolation(iso),
                counting=dict(
                    evals_per_try=1001,
                    hit_ratio_overlap_evals_not_counted=sum(
                        st[s]["rescue"]["tries"] for s in SAMPLES
                    ),
                    note=(
                        "each rescue try is 1,001 evaluations in the totals; the one overlap "
                        "evaluation per try for the hit ratio is not counted (as in the seed step)"
                    ),
                ),
                rescued_wilson={
                    s: dict(
                        rescued=st[s]["rescue"]["rescued"],
                        tried=st[s]["rescue"]["rescues_attempted"],
                        wilson_95=list(
                            ST.wilson(
                                st[s]["rescue"]["rescued"], st[s]["rescue"]["rescues_attempted"]
                            )
                        ),
                    )
                    for s in SAMPLES
                    if st[s]["rescue"]["rescues_attempted"]
                },
                predictions=preds,
                per_sample_passes=per,
                step_passes=verdict,
                rules={s: st[s]["rules"] for s in SAMPLES},
            ),
            indent=1,
            default=float,
        )
        + "\n"
    )
    print(
        json.dumps(
            dict(
                validity=val,
                per=per,
                verdict=verdict,
                preds={k: v["met"] for k, v in preds.items()},
            ),
            indent=1,
        )
    )


if __name__ == "__main__":
    main()
