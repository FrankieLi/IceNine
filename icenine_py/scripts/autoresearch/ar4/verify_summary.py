"""AR4 step 1: use 1 and use 2 on ``mid``, the four predictions -> ``summary.json``, ``tables.md``,
the plots.

  uv run python scripts/autoresearch/ar4/verify_summary.py

Reads the pair inputs of the three samples (``cache/inputs_<sample>.npz``), ``models.npz`` and
``thresholds.json`` (committed before ``mid`` was touched), the stored seed-step replays
(``cache/seeds_<sample>.json``), AR5 step 1's ``part2_cases.json`` (the fallback cases) with the
stored S8 log and ``run.json`` (``cache/stored/<sample>/S8_<sample>/``: the seed-step evaluations
of the reseed cases, subtracted to give grid-only evaluations, and the library's grid counter as
the check) and AR5's ``part3_votes.npz`` (only to count which voxels have a stored vote count).
No model is trained and no threshold is picked here.
"""

import gzip
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
AR4 = ICE / "benchmarks" / "autoresearch" / "ar4"
AR5 = ICE / "benchmarks" / "autoresearch" / "ar5"
sys.path.insert(0, str(ICE / "scripts" / "common"))
sys.path.insert(0, str(HERE))
import doc_tables as DT  # noqa: E402
import verify_fit as VF  # noqa: E402
import verify_lib as VL  # noqa: E402

SCORES = VF.SCORES
LABEL = {
    "B0": "B0 (hit ratio)",
    "B1": "B1 (logistic, summary numbers)",
    "N1": "N1 (set network)",
    "lit": "share of peaks with a lit pixel",
}
KINDS = (
    "after",
    "before",
    "try_start",
    "seed_end",
    "seed_cand",
    "turn_01",
    "turn_03",
    "turn_05",
)
SAMPLES = ("big", "grad", "mid")
KIND_SEED_END = 3
TRAIN_SAMPLES = ("big", "grad")


def f3(x: float) -> str:
    return f"{x:.3f}"


def grid_only(cases: List[Dict[str, Any]], s: str) -> Tuple[Dict[int, int], Dict[str, Any]]:
    """Grid-search evaluations per fallback voxel, without the seed step of the reseed cases.

    AR5's ``grid_evals`` of a reseed case is ``full_search_evaluations`` (seed step plus grid
    search); the seed step's evaluations are read from the stored log of S8 (``e1_log``) and
    subtracted. The total is checked against the library counter ``index_fallback_grid_evals``
    of S8's ``run.json``."""
    s8 = VF.VC.CACHE / "stored" / s / f"S8_{s}"
    log = json.loads(gzip.open(s8 / "e1_log.json.gz", "rt").read())
    counter = int(json.loads((s8 / "run.json").read_text())["stats"]["index_fallback_grid_evals"])
    seed_ev: Dict[int, List[int]] = {}
    for e in log:
        if e["stop"] == -1 and e["reseed"]:
            seed_ev.setdefault(int(e["voxel"]), []).append(int(e["evals"]))
    out: Dict[int, int] = {}
    sub = 0
    for c in cases:
        v = int(c["voxel"])
        ev = 0
        if c["reseed"]:
            assert len(seed_ev[v]) == 1, (s, v, seed_ev.get(v))
            ev = seed_ev[v][0]
            sub += ev
        out[v] = VL.grid_only_evals(c["grid_evals"], bool(c["reseed"]), ev)
    total = int(sum(out.values()))
    assert total == counter, (s, total, counter)
    check = dict(
        ar5_field_total=int(sum(int(c["grid_evals"]) for c in cases)),
        seed_step_evals_of_reseed_cases=sub,
        n_reseed_cases=int(sum(bool(c["reseed"]) for c in cases)),
        grid_only_total=total,
        library_counter_index_fallback_grid_evals=counter,
        equal=bool(total == counter),
    )
    return out, check


def main() -> None:
    thr = json.loads((AR4 / "thresholds.json").read_text())
    cuts = {k: v["cut"] for k, v in thr["thresholds"].items()}
    models = VF.load_models(AR4 / "models.npz")
    net = VF.net_from(models)
    inp = {s: VF.load_inputs(s) for s in SAMPLES}
    sc = {s: VF.all_scores(inp[s], models, net) for s in SAMPLES}
    info = {s: json.loads((AR4 / f"inputs_{s}.json").read_text()) for s in SAMPLES}
    out: Dict[str, Any] = dict(
        thresholds=thr["thresholds"],
        training=dict(
            n1_parameters=thr["n1_parameters"],
            n1_epochs=thr["n1_epochs"],
            n1_final_loss=thr["n1_final_loss"],
            n1_wall_s=thr["n1_wall_s"],
            b1_weights=thr["b1_weights"],
            n_train_pairs=thr["n_pairs"],
            n_use1_voxels=thr["n_use1_voxels"],
            n_use1_wrong=thr["n_use1_wrong"],
        ),
        pairs={s: dict(by_kind=info[s]["by_kind"], n_pairs=info[s]["n_pairs"]) for s in SAMPLES},
        sum_back={s: info[s]["sum_back"] for s in SAMPLES},
        seed_replay={
            s: dict(
                n_steps=info[s]["n_seed_steps"],
                match=info[s]["seed_replay_evals_and_hits_match"],
                mismatch=info[s]["seed_replay_mismatch"],
            )
            for s in SAMPLES
        },
    )
    # the pairs whose per-peak sums differ from the fused evaluation, by label (review fix 4)
    sbl: Dict[str, Any] = {}
    for s in SAMPLES:
        d = VL.sumback_differs(inp[s]["peaks"], inp[s]["offsets"], inp[s]["counts"])
        r = VL.right_label(inp[s]["err"])
        sbl[s] = dict(
            n_differ=int(d.sum()),
            n_right=int(r.sum()),
            n_wrong=int((~r).sum()),
            n_differ_right=int(d[r].sum()),
            n_differ_wrong=int(d[~r].sum()),
            share_differ_right=float(d[r].mean()),
            share_differ_wrong=float(d[~r].mean()),
            equal_to_stored_count=bool(
                int(d.sum()) == int(info[s]["sum_back"]["n_pairs_pixel_counts_differ_from_fused"])
            ),
            median_peak_pixels=float(np.median(inp[s]["peaks"][:, 4])),
            min_peak_pixels=float(inp[s]["peaks"][:, 4].min()),
        )
    out["sum_back_by_label"] = sbl

    # ---- use 1: the S8 maps after the pass (pairs of kind "after")
    use1: Dict[str, Any] = {}
    for s in SAMPLES:
        m = inp[s]["kind"] == 0
        wrong = ~VL.right_label(inp[s]["err"][m])
        row: Dict[str, Any] = dict(n_voxels=int(m.sum()), n_wrong=int(wrong.sum()))
        for name in SCORES:
            s_ = sc[s][name][m]
            fl, wf = VL.flag_counts(s_, wrong, cuts[name])
            lo, hi = VL.wilson_interval(wf, int(wrong.sum())) if wrong.sum() else (np.nan, np.nan)
            row[name] = dict(
                flagged=fl,
                wrong_flagged=wf,
                flagged_share=fl / int(m.sum()),
                wrong_flagged_wilson=[lo, hi],
                auc=VL.auc(s_, ~wrong),
            )
        b0 = row["B0"]
        for name in SCORES:
            r = row[name]
            r["pass_rule"] = bool(
                r["wrong_flagged"] >= b0["wrong_flagged"] and r["flagged"] <= b0["flagged"] / 2
            )
        use1[s] = row
    out["use1"] = use1

    # ---- the flag curve on mid (measured points: the cut at each wrong voxel's score)
    m = inp["mid"]["kind"] == 0
    wrong_mid = ~VL.right_label(inp["mid"]["err"][m])
    curve: Dict[str, Any] = {}
    for name in SCORES:
        s_ = sc["mid"][name][m]
        c, fl, wf = VL.curve_at_wrong_scores(s_, wrong_mid)
        curve[name] = dict(
            cut=c.tolist(), flagged=fl.tolist(), wrong_flagged=wf.tolist(), n=int(m.sum())
        )
    out["curve_mid"] = curve

    # ---- boundary voxels of mid (voxels with at least one try): prediction 4
    tries_vox = np.unique(np.load(AR5 / "starts_mid.npz")["voxel"])
    b_mask = np.isin(inp["mid"]["voxel"][m], tries_vox)
    bw = wrong_mid[b_mask]
    bnd: Dict[str, Any] = dict(n_voxels=int(b_mask.sum()), n_wrong=int(bw.sum()))
    for name in SCORES:
        s_ = sc["mid"][name][m][b_mask]
        fl, wf = VL.flag_counts(s_, bw, cuts[name])
        bnd[name] = dict(auc=VL.auc(s_, ~bw), flagged=fl, wrong_flagged=wf)
    # the measured distribution of the lit share at the boundary voxels (why auc is near 0.5)
    ls_b = sc["mid"]["lit"][m][b_mask]
    bnd["lit_detail"] = dict(
        share_right_all_lit=float((ls_b[~bw] == 1.0).mean()),
        share_wrong_all_lit=float((ls_b[bw] == 1.0).mean()),
        lowest_right=float(ls_b[~bw].min()),
        lowest_wrong=float(ls_b[bw].min()),
        n_distinct_right=int(len(np.unique(ls_b[~bw]))),
        n_distinct_wrong=int(len(np.unique(ls_b[bw]))),
        median_peaks_per_pair=float(np.median(inp["mid"]["n_peaks"][m][b_mask])),
    )
    out["boundary_mid"] = bnd
    out["all_wrong_mid_are_boundary"] = bool(int(bw.sum()) == int(wrong_mid.sum()))

    # ---- what-if on mid (cuts picked on mid itself; not a result of the step): how many voxels
    # each score must flag to reach 30 of 31 (95%) and 31 of 31 of the wrong voxels
    hind: Dict[str, Any] = {}
    for name in SCORES:
        s_ = sc["mid"][name][m]
        h: Dict[str, Any] = {}
        for tag, frac in (("95pct", 0.95), ("all", 1.0)):
            cut = VL.pick_threshold(s_, wrong_mid, frac)
            fl, wf = VL.flag_counts(s_, wrong_mid, cut)
            h[tag] = dict(cut=cut, flagged=fl, wrong_flagged=wf, flagged_share=fl / int(m.sum()))
        hind[name] = h
    out["whatif_mid_cut_picked_on_mid"] = hind
    # the wrong voxels of mid that N1 does not flag at its train cut
    vox_m = inp["mid"]["voxel"][m]
    summ_m = inp["mid"]["summary"][m]
    miss = np.where(wrong_mid & (sc["mid"]["N1"][m] > cuts["N1"]))[0]
    out["n1_missed_mid"] = [
        dict(
            voxel=int(vox_m[i]),
            err_deg=float(inp["mid"]["err"][m][i]),
            hit=float(summ_m[i, 0]),
            B0=float(sc["mid"]["B0"][m][i]),
            N1=float(sc["mid"]["N1"][m][i]),
            boundary=bool(vox_m[i] in set(tries_vox.tolist())),
        )
        for i in miss
    ]
    lit: Dict[str, Any] = {}
    for s in SAMPLES:
        mm = inp[s]["kind"] == 0
        ls = inp[s]["summary"][mm, VL.SUMMARY_NAMES.index("lit_share")]
        ww = ~VL.right_label(inp[s]["err"][mm])
        lit[s] = dict(
            share_of_voxels_with_all_peaks_lit=float((ls == 1.0).mean()),
            lowest_lit_share_of_wrong=float(ls[ww].min()),
            median_lit_share_of_wrong=float(np.median(ls[ww])),
            n_wrong_with_all_peaks_lit=int((ls[ww] == 1.0).sum()),
            n_wrong=int(ww.sum()),
        )
    out["lit_share_after_pass"] = lit

    # ---- use 2: the chosen end of each stored seed step
    part2 = json.loads((AR5 / "part2_cases.json").read_text())
    use2: Dict[str, Any] = {}
    for s in ("grad", "mid"):
        seeds = json.loads((VF.VC.CACHE / f"seeds_{s}.json").read_text())["steps"]
        pair_of = {
            (int(st), int(a)): i
            for i, (k, st, a) in enumerate(zip(inp[s]["kind"], inp[s]["step"], inp[s]["aux"]))
            if k == KIND_SEED_END
        }
        grid, grid_check = grid_only(part2[s], s)
        fb_vox = sorted(
            st["voxel"] for st in seeds if not st["passed_seed_test"] and st["taken_pos"] >= 0
        )
        rows_fb: List[Dict[str, Any]] = []
        rows_ps: List[Dict[str, Any]] = []
        for st in seeds:
            if st["taken_pos"] < 0:
                continue
            i = pair_of[(st["step"], st["taken_pos"])]
            r = dict(
                step=st["step"],
                voxel=st["voxel"],
                right=bool(inp[s]["err"][i] <= VL.RIGHT_DEG),
                err=float(inp[s]["err"][i]),
                **{n: float(sc[s][n][i]) for n in SCORES},
            )
            (rows_ps if st["passed_seed_test"] else rows_fb).append(r)
        u: Dict[str, Any] = dict(
            n_fallback=len(rows_fb),
            fallback_voxels_equal_ar5=bool(fb_vox == sorted(c["voxel"] for c in part2[s])),
            n_fb_right=sum(r["right"] for r in rows_fb),
            n_fb_wrong=sum(not r["right"] for r in rows_fb),
            n_passed=len(rows_ps),
            n_passed_right=sum(r["right"] for r in rows_ps),
        )
        grid_all = int(sum(grid[r["voxel"]] for r in rows_fb))
        mean_grid = grid_all / max(len(rows_fb), 1)
        fb_right = np.array([r["right"] for r in rows_fb], bool)
        fb_grid = np.array([grid[r["voxel"]] for r in rows_fb], np.int64)
        ps_right = np.array([r["right"] for r in rows_ps], bool)
        for name in SCORES:
            c = VL.use2_counts(
                np.array([r[name] for r in rows_fb]),
                fb_right,
                fb_grid,
                np.array([r[name] for r in rows_ps]),
                ps_right,
                cuts[name],
            )
            u[name] = dict(
                **c,
                right_accepted_wilson=VL.wilson_interval(
                    c["right_accepted"], max(u["n_fb_right"], 1)
                ),
                wrong_accepted_wilson=VL.wilson_interval(
                    c["wrong_accepted"], max(u["n_fb_wrong"], 1)
                ),
                rules_arithmetic_not_run=VL.use2_rule_costs(c, mean_grid),
            )
        u["rows_fallback"] = rows_fb
        u["grid_evals_all_fallbacks"] = grid_all
        u["grid_evals_mean_per_fallback"] = mean_grid
        u["grid_evals_check"] = grid_check
        use2[s] = u
    out["use2"] = use2

    # ---- predictions
    votes = np.load(AR5 / "part3_votes.npz")
    u1, u2 = use1["mid"], use2["mid"]
    n15 = u2["n_fb_right"]
    out["predictions"] = dict(
        p1_B1_does_not_pass_use1=dict(
            holds=not u1["B1"]["pass_rule"], pass_rule=u1["B1"]["pass_rule"]
        ),
        p2_N1_passes_use1=dict(holds=u1["N1"]["pass_rule"], pass_rule=u1["N1"]["pass_rule"]),
        p3_N1_use2=dict(
            holds=bool(
                n15 > 0
                and u2["N1"]["right_accepted"] >= 2 * n15 / 3
                and u2["N1"]["wrong_accepted"] <= 1
            ),
            right_accepted=u2["N1"]["right_accepted"],
            of_right=n15,
            wrong_accepted=u2["N1"]["wrong_accepted"],
        ),
        p4_lit_better_than_hit_at_boundary=dict(
            holds=bool(bnd["lit"]["auc"] > bnd["B0"]["auc"]),
            auc_lit=bnd["lit"]["auc"],
            auc_hit=bnd["B0"]["auc"],
            auc_N1=bnd["N1"]["auc"],
        ),
    )
    out["derived"] = dict(
        total_pairs=int(sum(out["pairs"][s]["n_pairs"] for s in SAMPLES)),
        half_of_B0_flagged_mid=u1["B0"]["flagged"] / 2,
        two_thirds_of_right_fallback_mid=2 * n15 / 3,
        n1_missed_wrong_mid=int(len(out["n1_missed_mid"])),
        n1_missed_wrong_train=int(
            thr["n_use1_wrong"] - thr["thresholds"]["N1"]["train_wrong_flagged"]
        ),
        n_use1_wrong_train=int(thr["n_use1_wrong"]),
        replace_mean_grid_evals_mid=u2["grid_evals_mean_per_fallback"],
        n1_or_rule_saved_share_of_fallback_grid_mid=u2["N1"]["grid_evals_of_right_accepted"]
        / u2["grid_evals_all_fallbacks"],
        ar5_stored_votes_voxels={s: int(len(votes[f"{s}_voxel"])) for s in SAMPLES},
        n_passed_wrong_mid=u2["n_passed"] - u2["n_passed_right"],
        n1_passed_refused_wrong_mid=u2["N1"]["passed_refused"] - u2["N1"]["passed_refused_right"],
        n_wrong_with_lit_under_1_train=int(
            sum(
                out["lit_share_after_pass"][s]["n_wrong"]
                - out["lit_share_after_pass"][s]["n_wrong_with_all_peaks_lit"]
                for s in TRAIN_SAMPLES
            )
        ),
    )
    (AR4 / "summary.json").write_text(json.dumps(out, indent=1, default=float) + "\n")

    # ---- tables
    DT.write_tables(AR4 / "tables.md", tables(out))
    plots(out, sc, inp, cuts)
    print(json.dumps(out["predictions"], default=float), flush=True)


def tables(o: Dict[str, Any]) -> Dict[str, str]:
    t: Dict[str, str] = {}
    rows = []
    for s in SAMPLES:
        r = dict(sample=s, pairs=o["pairs"][s]["n_pairs"])
        r.update(o["pairs"][s]["by_kind"])
        rows.append(r)
    t["ar4_pairs"] = DT.markdown_table(rows, ["sample", "pairs", *KINDS])
    rows = []
    for s in SAMPLES:
        b = o["sum_back_by_label"][s]
        rows.append(
            dict(
                sample=s,
                right=f"{b['n_differ_right']} of {b['n_right']} "
                f"({100 * b['share_differ_right']:.2f}%)",
                wrong=f"{b['n_differ_wrong']} of {b['n_wrong']} "
                f"({100 * b['share_differ_wrong']:.2f}%)",
            )
        )
    t["ar4_sumback"] = DT.markdown_table(rows, ["sample", "right", "wrong"])
    rows = []
    for name in SCORES:
        d = o["thresholds"][name]
        rows.append(
            dict(
                score=LABEL[name],
                cut=f3(d["cut"]),
                flagged=d["train_voxels_flagged"],
                wrong=d["train_wrong_flagged"],
                share=f"{100 * d['train_flagged_share']:.2f}%",
            )
        )
    t["ar4_cuts"] = DT.markdown_table(
        rows,
        ["score", "cut", "flagged", "wrong", "share"],
        labels=None,
    )
    for s in ("mid", "big", "grad"):
        u = o["use1"][s]
        rows = []
        for name in SCORES:
            r = u[name]
            lo, hi = r["wrong_flagged_wilson"]
            rows.append(
                dict(
                    score=LABEL[name],
                    wrong=f"{r['wrong_flagged']} of {u['n_wrong']} ({lo:.2f} to {hi:.2f})",
                    flagged=f"{r['flagged']} of {u['n_voxels']}",
                    share=f"{100 * r['flagged_share']:.2f}%",
                    auc=f3(r["auc"]),
                    rule="baseline" if name == "B0" else ("yes" if r["pass_rule"] else "no"),
                )
            )
        t[f"ar4_use1_{s}"] = DT.markdown_table(
            rows, ["score", "wrong", "flagged", "share", "auc", "rule"]
        )
    b = o["boundary_mid"]
    rows = [
        dict(
            score=LABEL[n],
            auc=f3(b[n]["auc"]),
            wrong=f"{b[n]['wrong_flagged']} of {b['n_wrong']}",
            flagged=f"{b[n]['flagged']} of {b['n_voxels']}",
        )
        for n in SCORES
    ]
    t["ar4_boundary"] = DT.markdown_table(rows, ["score", "auc", "wrong", "flagged"])
    h = o["whatif_mid_cut_picked_on_mid"]
    rows = [
        dict(
            score=LABEL[n],
            p95=f"{h[n]['95pct']['flagged']} ({100 * h[n]['95pct']['flagged_share']:.2f}%)",
            p100=f"{h[n]['all']['flagged']} ({100 * h[n]['all']['flagged_share']:.2f}%)",
        )
        for n in SCORES
    ]
    t["ar4_whatif"] = DT.markdown_table(rows, ["score", "p95", "p100"], labels=None)
    rows = [
        dict(
            voxel=r["voxel"],
            err=f"{r['err_deg']:.2f}",
            hit=f"{r['hit']:.3f}",
            n1=f"{r['N1']:.3f}",
            bnd="yes" if r["boundary"] else "no",
        )
        for r in o["n1_missed_mid"]
    ]
    t["ar4_missed"] = DT.markdown_table(rows, ["voxel", "err", "hit", "n1", "bnd"])
    rows = [
        dict(
            sample=s,
            all_lit=f"{100 * v['share_of_voxels_with_all_peaks_lit']:.1f}%",
            wrong_all_lit=f"{v['n_wrong_with_all_peaks_lit']} of {v['n_wrong']}",
            lowest=f"{v['lowest_lit_share_of_wrong']:.3f}",
            median=f"{v['median_lit_share_of_wrong']:.3f}",
        )
        for s, v in o["lit_share_after_pass"].items()
    ]
    t["ar4_lit"] = DT.markdown_table(
        rows, ["sample", "all_lit", "wrong_all_lit", "lowest", "median"]
    )
    for s in ("mid", "grad"):
        u = o["use2"][s]
        rows = []
        for name in SCORES:
            r = u[name]
            (rl, rh), (wl, wh) = r["right_accepted_wilson"], r["wrong_accepted_wilson"]
            rows.append(
                dict(
                    score=LABEL[name],
                    right=f"{r['right_accepted']} of {u['n_fb_right']} ({rl:.2f} to {rh:.2f})",
                    wrong=f"{r['wrong_accepted']} of {u['n_fb_wrong']} ({wl:.2f} to {wh:.2f})",
                    saved=f"{r['grid_evals_of_right_accepted']:,}",
                    refused=f"{r['passed_refused']} of {u['n_passed']}"
                    f" ({r['passed_refused_right']} right)",
                )
            )
        t[f"ar4_use2_{s}"] = DT.markdown_table(
            rows, ["score", "right", "wrong", "saved", "refused"]
        )
    # the two readings of use 2 on mid, priced by arithmetic from stored numbers (not a run)
    u = o["use2"]["mid"]
    rows = []
    for name in SCORES:
        r, a = u[name], u[name]["rules_arithmetic_not_run"]
        rows.append(
            dict(
                score=LABEL[name],
                or_saved=f"{a['or_saved']:,.0f}",
                or_wrong=f"{r['wrong_accepted']} of {u['n_fb_wrong']}",
                rep_refused=f"{r['passed_refused_right']} right, {r['passed_refused_wrong']} wrong",
                rep_extra_right=f"{a['replace_extra_right']:,.0f}",
                rep_extra_wrong=f"{a['replace_extra_wrong']:,.0f}",
            )
        )
    t["ar4_use2_rules"] = DT.markdown_table(
        rows,
        ["score", "or_saved", "or_wrong", "rep_refused", "rep_extra_right", "rep_extra_wrong"],
    )
    # grid-only evaluations of the fallback cases against AR5's field
    rows = []
    for s in ("grad", "mid"):
        g = o["use2"][s]["grid_evals_check"]
        rows.append(
            dict(
                sample=s,
                cases=o["use2"][s]["n_fallback"],
                reseed=g["n_reseed_cases"],
                ar5_field=f"{g['ar5_field_total']:,}",
                seed_steps=f"{g['seed_step_evals_of_reseed_cases']:,}",
                grid_only=f"{g['grid_only_total']:,}",
                counter=f"{g['library_counter_index_fallback_grid_evals']:,}",
            )
        )
    t["ar4_grid_evals"] = DT.markdown_table(
        rows, ["sample", "cases", "reseed", "ar5_field", "seed_steps", "grid_only", "counter"]
    )
    return t


def plots(
    o: Dict[str, Any], sc: Dict[str, Any], inp: Dict[str, Any], cuts: Dict[str, float]
) -> None:
    # 1. the flag curve on mid
    fig, ax = plt.subplots(figsize=(7, 4.5))
    n_wrong = o["use1"]["mid"]["n_wrong"]
    colors = dict(B0="tab:gray", B1="tab:blue", N1="tab:red", lit="tab:green")
    for name in SCORES:
        c = o["curve_mid"][name]
        x = np.array(c["wrong_flagged"])
        y = 100 * np.array(c["flagged"]) / c["n"]
        # markers only: the measured points (no line between them)
        ax.plot(x, y, "o", ms=3, color=colors[name], label=LABEL[name])
        r = o["use1"]["mid"][name]
        ax.plot(
            [r["wrong_flagged"]],
            [100 * r["flagged_share"]],
            "*",
            ms=14,
            mfc="none",
            mec=colors[name],
            mew=1.5,
        )
    ax.set_xlabel(f"wrong voxels flagged (of {n_wrong}) on mid")
    ax.set_ylabel("all voxels flagged (%)")
    ax.set_title("Flag curve on mid; open star: the cut picked on the train samples")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(AR4 / "ar4_flag_curve.png", dpi=130)
    plt.close(fig)
    # 2. use 2: scores of the fallback ends of mid
    rows = o["use2"]["mid"]["rows_fallback"]
    fig, axs = plt.subplots(1, 4, figsize=(11, 3.8), sharey=True)
    for ax, name in zip(axs, SCORES):
        for j, r in enumerate(rows):
            ax.plot(
                [0 if r["right"] else 1],
                [r[name]],
                "o",
                color="tab:green" if r["right"] else "tab:red",
                alpha=0.7,
            )
        ax.axhline(cuts[name], color="k", ls="--", lw=1)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["right end", "wrong end"])
        ax.set_title(LABEL[name], fontsize=8)
    axs[0].set_ylabel("score of the refined end")
    fig.suptitle("Use 2 on mid: the 15 seed steps that went to a grid search; dashed: the cut")
    fig.tight_layout()
    fig.savefig(AR4 / "ar4_use2.png", dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    main()
