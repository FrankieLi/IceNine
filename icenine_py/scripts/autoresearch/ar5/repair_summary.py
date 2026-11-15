"""AR5 step 19 wrap-up (the step stopped at its gate, D-20261115-N8): summary, tables, plot and the
manifest of what the repair stored per fit. Reads only files that exist (no new run):
``gate_fail.json`` and ``flags.json`` (tracked) and the repair's ``fits.npz`` per sample (cache).

  uv run python scripts/autoresearch/ar5/repair_summary.py

Writes in ``benchmarks/autoresearch/ar5/repair/``: ``summary.json``, ``tables.md`` (blocks for
``sync_doc_tables.py``), ``ar5_repair_gate.png`` and ``fits_manifest.json``.
"""

import json
from pathlib import Path
from typing import Any, Dict

import shutil
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
RD = ICE / "benchmarks" / "autoresearch" / "ar5" / "repair"
CACHE = RD / "cache"
sys.path.insert(0, str(ICE / "scripts" / "common"))
sys.path.insert(0, str(HERE))
import doc_tables as DT  # noqa: E402
import repair_data as RD_  # noqa: E402
import repair_rules as RR  # noqa: E402
import stats  # noqa: E402

ORDER = ["held3", "held", "mid", "grad", "big", "mid_J5L"]


def main() -> None:
    gate = json.loads((RD / "gate_fail.json").read_text())
    rows = gate["rows"]
    summ: Dict[str, Any] = {}
    manifest_rows: Dict[str, Any] = {}
    for k in ORDER:
        r = rows[k]
        arm, smp = r["arm"], r["sample"]
        fz = np.load(CACHE / f"{arm}R_{smp}" / "fits.npz")
        ev = fz["evals"]
        s = dict(r)
        s.update(
            cap_pct=(r["gate"]["cap_pct"] if "gate" in r else None),
            mean_evals_per_fit=float(ev.mean()),
            min_evals=int(ev.min()),
            max_evals=int(ev.max()),
            q25=float(np.percentile(ev, 25)),
            q75=float(np.percentile(ev, 75)),
            own_end_fits=int((fz["start"] < 0).sum()),
        )
        own = fz["start"] < 0
        s.update(
            wilson_J5_pct=[round(100 * x, 3) for x in stats.wilson(int(r["wrong_J5"]), 24570)],
            wilson_J5R_pct=[round(100 * x, 3) for x in stats.wilson(int(r["wrong_J5R"]), 24570)],
            frac_fits_at_662=float((ev == 662).mean()),
            own_end_frac_at_662=float((ev[own] == 662).mean()),
            cand_frac_at_662=float((ev[~own] == 662).mean()),
            own_end_evals_share_pct=float(100.0 * ev[own].sum() / ev.sum()),
            max_median_err_move_deg=abs(r["median_err_J5"] - r["median_err_J5R"]),
        )
        summ[k] = s
        dst = RD / "fits" / f"{arm}R_{smp}_fits.npz"  # the only per-fit record: tracked copy
        dst.parent.mkdir(exist_ok=True)
        shutil.copyfile(CACHE / f"{arm}R_{smp}" / "fits.npz", dst)
        manifest_rows[k] = dict(
            tracked_copy=f"benchmarks/autoresearch/ar5/repair/fits/{arm}R_{smp}_fits.npz",
            folder=f"benchmarks/autoresearch/ar5/repair/cache/{arm}R_{smp}/",
            n_fits=int(len(ev)),
            n_flagged=r["flagged"],
            n_evals=r["repair_evals"],
        )
    summ["max_median_err_move_deg"] = max(
        summ[k]["max_median_err_move_deg"] for k in ORDER if k != "mid_J5L"
    )
    # J5L repaired against S8 (mid's stored S8 map: the same sample) and against J5R (prediction 6)
    y = RD_.sym()
    dm, dl = RD_.load("mid", "J5"), RD_.load("mid", "J5L")
    Rt = dm["R_t"]
    wl = RR.wrong_mask(np.load(CACHE / "J5LR_mid" / "after.npy"), Rt, y)
    wl0 = RR.wrong_mask(dl["R"], Rt, y)
    w8 = RR.wrong_mask(dm["R_s8"], Rt, y)
    wr = RR.wrong_mask(np.load(CACHE / "J5R_mid" / "after.npy"), Rt, y)
    w5 = RR.wrong_mask(dm["R"], Rt, y)
    f8, b8 = RR.fixed_broken(w8, wl)
    f80, b80 = RR.fixed_broken(w8, wl0)
    fr, br = RR.fixed_broken(wr, wl)
    f5, b5 = RR.fixed_broken(w5, wl)
    summ["J5L_mid"] = dict(
        wrong_J5L=int(wl0.sum()),
        wrong_J5LR=int(wl.sum()),
        wilson_J5L_pct=[round(100 * x, 3) for x in stats.wilson(int(wl0.sum()), 24570)],
        wilson_J5LR_pct=[round(100 * x, 3) for x in stats.wilson(int(wl.sum()), 24570)],
        vs_S8_before=dict(fixed=f80, broken=b80, net=b80 - f80),
        vs_S8_after=dict(fixed=f8, broken=b8, net=b8 - f8),
        vs_J5R=dict(fixed=fr, broken=br),
        vs_J5=dict(fixed=f5, broken=b5),
    )
    # the one voxel J5R broke on big
    db, dbr = RD_.load("big", "J5"), np.load(CACHE / "J5R_big" / "after.npy")
    wb0 = RR.wrong_mask(db["R"], db["R_t"], y)
    wb1 = RR.wrong_mask(dbr, db["R_t"], y)
    fzb = np.load(CACHE / "J5R_big" / "fits.npz")
    broken_v = [int(v) for v in np.where(~wb0 & wb1)[0]]
    big_rows = []
    for v in broken_v:
        err0 = float(RR.misorientation_deg(db["R"][v], db["R_t"][v], y))
        err1 = float(RR.misorientation_deg(dbr[v], db["R_t"][v], y))
        m = np.where(fzb["voxel"] == v)[0]
        best = m[np.argmax(fzb["hit"][m])]
        big_rows.append(
            dict(
                voxel=v,
                err_before_deg=err0,
                err_after_deg=err1,
                hit_before=float(db["hit"][v]),
                hit_after=float(fzb["hit"][best]),
                start_voxel=int(fzb["start"][best]),
                ratio=float(db["ratio"][v]),
                margin=float(db["margin"][v]),
            )
        )
    summ["big_broken_voxels"] = big_rows
    (RD / "summary.json").write_text(json.dumps(summ, indent=1) + "\n")

    t1 = DT.markdown_table(
        [
            dict(
                sample=k if k != "mid_J5L" else "mid, J5L",
                flagged=summ[k]["flagged"],
                fits=summ[k]["fits"],
                changed=summ[k]["changed"],
                evals=summ[k]["repair_evals"],
                per_fit_med=summ[k]["median_evals_per_fit"],
                per_fit_mean=summ[k]["mean_evals_per_fit"],
                share=summ[k]["share_pct"],
                cap=summ[k]["cap_pct"],
                verdict=("cost over the cap" if "gate" in summ[k] else "measurement"),
            )
            for k in ORDER
        ],
        [
            "sample",
            "flagged",
            "fits",
            "changed",
            "evals",
            "per_fit_med",
            "per_fit_mean",
            "share",
            "cap",
            "verdict",
        ],
        formats=dict(evals=",d", per_fit_med=".0f", per_fit_mean=".0f", share=".2f", cap=".1f"),
    )
    t2 = DT.markdown_table(
        [
            dict(
                sample=k if k != "mid_J5L" else "mid, J5L",
                wrong_before=summ[k]["wrong_J5"],
                ci_before="%.3f to %.3f%%" % tuple(summ[k]["wilson_J5_pct"]),
                wrong_after=summ[k]["wrong_J5R"],
                ci_after="%.3f to %.3f%%" % tuple(summ[k]["wilson_J5R_pct"]),
                fixed=summ[k]["fixed_vs_J5"],
                broken=summ[k]["broken_vs_J5"],
                med_before=summ[k]["median_err_J5"],
                med_after=summ[k]["median_err_J5R"],
            )
            for k in ORDER
        ],
        [
            "sample",
            "wrong_before",
            "ci_before",
            "wrong_after",
            "ci_after",
            "fixed",
            "broken",
            "med_before",
            "med_after",
        ],
        formats=dict(med_before=".5f", med_after=".5f"),
    )
    jl = summ["J5L_mid"]
    t4 = DT.markdown_table(
        [
            dict(
                comparison="J5L against S8 (before the repair)",
                fixed=jl["vs_S8_before"]["fixed"],
                broken=jl["vs_S8_before"]["broken"],
                net=jl["vs_S8_before"]["net"],
            ),
            dict(
                comparison="J5L repaired against S8",
                fixed=jl["vs_S8_after"]["fixed"],
                broken=jl["vs_S8_after"]["broken"],
                net=jl["vs_S8_after"]["net"],
            ),
            dict(
                comparison="J5L repaired against J5R (mid)",
                fixed=jl["vs_J5R"]["fixed"],
                broken=jl["vs_J5R"]["broken"],
                net=jl["vs_J5R"]["broken"] - jl["vs_J5R"]["fixed"],
            ),
        ],
        ["comparison", "fixed", "broken", "net"],
    )
    t5 = DT.markdown_table(
        [
            dict(
                sample=k if k != "mid_J5L" else "mid, J5L",
                own_fits=summ[k]["own_end_fits"],
                own_at_662=summ[k]["own_end_frac_at_662"] * 100,
                cand_at_662=summ[k]["cand_frac_at_662"] * 100,
                own_share=summ[k]["own_end_evals_share_pct"],
            )
            for k in ORDER
        ],
        ["sample", "own_fits", "own_at_662", "cand_at_662", "own_share"],
        formats=dict(own_at_662=".0f", cand_at_662=".0f", own_share=".1f"),
    )
    t3 = DT.markdown_table(
        [
            dict(
                sample=k,
                j5_fixed=summ[k]["J5_vs_S8"]["fixed"],
                j5_broken=summ[k]["J5_vs_S8"]["broken"],
                j5_net=summ[k]["J5_vs_S8"]["net"],
                r_fixed=summ[k]["J5R_vs_S8"]["fixed"],
                r_broken=summ[k]["J5R_vs_S8"]["broken"],
                r_net=summ[k]["J5R_vs_S8"]["net"],
            )
            for k in ORDER
            if "J5_vs_S8" in summ[k]
        ],
        ["sample", "j5_fixed", "j5_broken", "j5_net", "r_fixed", "r_broken", "r_net"],
    )
    DT.write_tables(
        RD / "tables.md",
        dict(
            repair_cost=t1,
            repair_accuracy=t2,
            repair_vs_s8=t3,
            repair_j5l=t4,
            repair_own_end=t5,
        ),
    )

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    five = ORDER[:5]
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    x = np.arange(len(five))
    ax[0].bar(x, [summ[k]["share_pct"] for k in five], color="#c0504d", label="repair cost")
    for i, k in enumerate(five):
        ax[0].plot([i - 0.4, i + 0.4], [summ[k]["cap_pct"]] * 2, color="k", lw=2)
    ax[0].plot([], [], color="k", lw=2, label="cap (2.5%; 3.5% on grad)")
    ax[0].set_xticks(x)
    ax[0].set_xticklabels(five)
    ax[0].set_ylabel("repair evaluations, % of the J5 run")
    ax[0].set_title("Cost: over the cap on every sample")
    ax[0].legend(fontsize=8)
    ax[1].bar(x - 0.2, [summ[k]["fixed_vs_J5"] for k in five], 0.4, color="#4f81bd", label="fixed")
    ax[1].bar(
        x + 0.2, [summ[k]["broken_vs_J5"] for k in five], 0.4, color="#c0504d", label="broken"
    )
    ax[1].set_xticks(x)
    ax[1].set_xticklabels(five)
    ax[1].set_ylabel("voxels, against J5 (same run)")
    ax[1].set_title("Accuracy: fixed and broken voxels")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(RD / "ar5_repair_gate.png", dpi=110)

    manifest = dict(
        what="What the repair (arm J5R, step 19) stored per fit. Written after the gate failed; "
        "no new run.",
        per_fit_file="<folder>/fits.npz (untracked cache; one row per MC fit, in the order made)",
        columns=dict(
            voxel="the contested voxel (int32)",
            start="-1 for the fit from the voxel's own end, else the candidate voxel whose "
            "snapshot orientation was the start (int32)",
            end="the end orientation, 3x3 float64",
            cost="strict cost of the end (final overlap evaluation)",
            hit="hit ratio of the end",
            evals="evaluations of the whole fit, final overlap evaluation included "
            "(Reconstructor.last_local_optimization_evals)",
        ),
        also=dict(
            repair_json="<folder>/repair.json: counts, changed voxels, wall (contended), "
            "median evaluations per fit",
            after_npy="<folder>/after.npy: the map after the repair",
        ),
        trace_stored=False,
        trace_note="No per-evaluation trace and no per-gain trace of any MC fit is stored: only "
        "the end and the total evaluations of each fit. A later step that wants to replay a "
        "cheaper stop on these fits would have to re-run them, with one generator "
        "default_rng(0) over the voxels in index order (a fit's draws depend on every earlier "
        "fit, so a replay must keep that order).",
        samples=manifest_rows,
        not_stored="the list of flagged voxels (it is recomputed by repair_rules from the stored "
        "pass tries; flags.json holds only the counts) and the candidates per voxel",
    )
    (RD / "fits_manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
    print("wrote summary.json tables.md ar5_repair_gate.png fits_manifest.json")


if __name__ == "__main__":
    main()
