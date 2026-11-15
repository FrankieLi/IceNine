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

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
RD = ICE / "benchmarks" / "autoresearch" / "ar5" / "repair"
CACHE = RD / "cache"
import sys  # noqa: E402

sys.path.insert(0, str(ICE / "scripts" / "common"))
import doc_tables as DT  # noqa: E402

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
        summ[k] = s
        manifest_rows[k] = dict(
            folder=f"benchmarks/autoresearch/ar5/repair/cache/{arm}R_{smp}/",
            n_fits=int(len(ev)),
            n_flagged=r["flagged"],
            n_evals=r["repair_evals"],
        )
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
                wrong_after=summ[k]["wrong_J5R"],
                fixed=summ[k]["fixed_vs_J5"],
                broken=summ[k]["broken_vs_J5"],
                med_before=summ[k]["median_err_J5"],
                med_after=summ[k]["median_err_J5R"],
            )
            for k in ORDER
        ],
        ["sample", "wrong_before", "wrong_after", "fixed", "broken", "med_before", "med_after"],
        formats=dict(med_before=".5f", med_after=".5f"),
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
    DT.write_tables(RD / "tables.md", dict(repair_cost=t1, repair_accuracy=t2, repair_vs_s8=t3))

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
