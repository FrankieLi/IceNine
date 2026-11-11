"""AR4 step 1: use 1 and use 2 on ``mid``, the four predictions -> ``summary.json``, ``tables.md``,
the plots.

  uv run python scripts/autoresearch/ar4/verify_summary.py

Reads the pair inputs of the three samples (``cache/inputs_<sample>.npz``), ``models.npz`` and
``thresholds.json`` (committed before ``mid`` was touched), the stored seed-step replays
(``cache/seeds_<sample>.json``) and AR5 step 1's ``part2_cases.json`` (stored grid-search
evaluations of the fallback cases). No model is trained and no threshold is picked here.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List

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


def f3(x: float) -> str:
    return f"{x:.3f}"


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
    out["boundary_mid"] = bnd

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
        grid = {c["voxel"]: c["grid_evals"] for c in part2[s]}
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
        for name in SCORES:
            acc = [r for r in rows_fb if r[name] >= cuts[name]]
            ref = [r for r in rows_ps if r[name] < cuts[name]]
            u[name] = dict(
                right_accepted=sum(r["right"] for r in acc),
                wrong_accepted=sum(not r["right"] for r in acc),
                grid_evals_of_right_accepted=int(sum(grid[r["voxel"]] for r in acc if r["right"])),
                passed_refused=len(ref),
                passed_refused_right=sum(r["right"] for r in ref),
            )
        u["rows_fallback"] = rows_fb
        u["grid_evals_all_fallbacks"] = int(sum(grid[r["voxel"]] for r in rows_fb))
        use2[s] = u
    out["use2"] = use2

    # ---- predictions
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
                    rule="yes" if r["pass_rule"] else "no",
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
    for s in ("mid", "grad"):
        u = o["use2"][s]
        rows = []
        for name in SCORES:
            r = u[name]
            rows.append(
                dict(
                    score=LABEL[name],
                    right=f"{r['right_accepted']} of {u['n_fb_right']}",
                    wrong=f"{r['wrong_accepted']} of {u['n_fb_wrong']}",
                    saved=f"{r['grid_evals_of_right_accepted']:,}",
                    refused=f"{r['passed_refused']} of {u['n_passed']}"
                    f" ({r['passed_refused_right']} right)",
                )
            )
        t[f"ar4_use2_{s}"] = DT.markdown_table(
            rows, ["score", "right", "wrong", "saved", "refused"]
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
        ax.plot(x, y, "o-", ms=3, color=colors[name], label=LABEL[name])
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
