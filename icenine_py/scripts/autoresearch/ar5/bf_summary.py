"""AR5 step 24, phase A step 8: the numbers, the plot and the verdict of arm ``J5SDB`` on the four
disclosed samples run live (``mid``, ``held``, ``grad``, ``big``) and the clean test
sample ``held7``, by the frozen decision rules.

  uv run python scripts/autoresearch/ar5/bf_summary.py

Reads (read only) the runs of ``bf_run_all.sh`` in ``benchmarks/autoresearch/ar5/bf/cache``
(``J5DB_<s>``, ``j5sdb/J5S_<s>``; on ``held7`` also ``S8_held7``, ``J5D_held7``, ``j5sd/J5S_held7``
and their passes), step 22's stored arms (links) and tracked ``summary.json``, and the truth maps.
``held7`` is opened only after ``logs/bf_runs.done`` exists (``bf_guard``). Writes
``summary.json``, ``tables.md`` (blocks for ``sync_doc_tables.py``) and the figure
``ar5_bf_evals_wrong.png`` (the plot first: run evaluations against wrong voxels, five panels:
``S8``, ``J5SD``, ``J5SDB``). No p-values: counts, Wilson intervals (``scripts/common/stats.py``)
and paired fixed / broken by voxel (wrong: over 1 degree from the truth, symmetry-reduced). The
decision rules are ``bf_rules.JUDGE`` (frozen); nothing is tuned here.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bf_guard as BG  # noqa: E402
import bf_rules as BR  # noqa: E402
import bf_score as BS  # noqa: E402
import doc_tables as DT  # noqa: E402

RD = BS.ICE / "benchmarks" / "autoresearch" / "ar5" / "bf"
CACHE = BS.CACHE
ORDER = ("held7", "mid", "held", "grad", "big")


def inside(x: float, lo: float, hi: float) -> bool:
    return bool(lo <= x <= hi)


def step_verdict(rows: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """The step passes only if ``J5SDB`` passes on ``held7`` and on the four samples run live."""
    per = {s: bool(rows[s]["vs_j5sd"]["ok"] and rows[s]["vs_s8"]["ok"]) for s in rows}
    need = set(BR.LIVE_SAMPLES) | {BR.NEW_SAMPLE}
    return dict(per_sample=per, step_passes=bool(need <= set(per) and all(per[s] for s in need)))


def predictions(rows: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    bands = dict(held=(4.5, 10), held7=(4.5, 10), mid=(5, 9.5), grad=(2, 5), big=(1.5, 4))
    p3: Dict[str, Any] = {}
    for s, r in rows.items():
        lo, hi = bands[s]
        sav = 100 * r["fewer"]["j5sd"]
        p3[s] = dict(
            saving_pct=sav,
            band=[lo, hi],
            saving_inside=inside(sav, lo, hi),
            net=r["pair"]["j5sd"]["net"],
            net_inside=inside(r["pair"]["j5sd"]["net"], -3, 3),
            grid_within_2=bool(abs(r["grid"]["j5db"] - r["grid"]["j5d"]) <= 2),
            seed_test_within_10=bool(abs(r["seed_test"]["j5db"] - r["seed_test"]["j5d"]) <= 10),
            median_higher_deg=r["pair"]["j5sd"]["median_higher"],
            median_inside=bool(abs(r["pair"]["j5sd"]["median_higher"]) <= 0.0005),
            house_rule=r["vs_s8"]["ok"],
        )
    return p3


def figure(rows: Dict[str, Dict[str, Any]], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(2, 3, figsize=(12, 6.5))
    n = 24570
    for ax, s in zip(axs.ravel(), ORDER):
        r = rows[s]
        for lab, k, col, mk in (
            ("S8", "s8", "tab:gray", "s"),
            ("J5SD", "j5sd", "tab:blue", "o"),
            ("J5SDB", "j5sdb", "tab:red", "D"),
        ):
            x = r["totals"][k] / 1e6
            ax.scatter(x, r["wrong"][k]["n"], c=col, marker=mk, s=60, label=lab)
            lo, hi = r["wrong"][k]["wilson"]
            ax.vlines(x, lo * n, hi * n, color=col, lw=1)
        ax.set_title(f"{s}{' (clean test)' if s == BR.NEW_SAMPLE else ''}")
        ax.set_xlabel("evaluations (millions)")
        ax.set_ylabel("wrong voxels")
        ax.set_ylim(bottom=0)
    axs[0, 0].legend()
    axs[1, 2].axis("off")
    fig.suptitle("Run evaluations against wrong voxels (over 1 degree), Wilson 95% bars")
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def tables(rows: Dict[str, Dict[str, Any]]) -> Dict[str, str]:
    o = ORDER
    ev = [
        dict(
            sample=s,
            s8=rows[s]["totals"]["s8"],
            j5sd=rows[s]["totals"]["j5sd"],
            j5sdb=rows[s]["totals"]["j5sdb"],
            vs_j5sd=100 * rows[s]["fewer"]["j5sd"],
            vs_s8=100 * rows[s]["fewer"]["s8"],
        )
        for s in o
    ]
    acc = []
    for s in o:
        r = rows[s]
        w = r["wrong"]
        lo, hi = w["j5sdb"]["wilson"]
        acc.append(
            dict(
                sample=s,
                w_s8=w["s8"]["n"],
                w_j5sd=w["j5sd"]["n"],
                w_j5sdb=w["j5sdb"]["n"],
                wilson=f"{100 * lo:.2f}-{100 * hi:.2f}%",
                fb_j5sd=f"{r['pair']['j5sd']['fixed']}/{r['pair']['j5sd']['broken']}",
                fb_s8=f"{r['pair']['s8']['fixed']}/{r['pair']['s8']['broken']}",
                med_j5sd=r["median_err_deg"]["j5sd"],
                med_j5sdb=r["median_err_deg"]["j5sdb"],
            )
        )
    dec = []
    for s in o:
        r = rows[s]
        a, b = r["vs_j5sd"], r["vs_s8"]
        dec.append(
            dict(
                sample=s,
                fewer=100 * r["fewer"]["j5sd"],
                r1=a["r1"],
                net=r["pair"]["j5sd"]["net"],
                r2=a["r2"],
                med=r["pair"]["j5sd"]["median_higher"],
                r3=a["r3"],
                house=b["ok"],
                passes=bool(a["ok"] and b["ok"]),
            )
        )
    stg = []
    for s in o:
        d = dict(rows[s]["stages"])
        d["sample"] = s
        stg.append(d)
    bl = []
    for s in o:
        b = rows[s]["ball"]
        bl.append(
            dict(
                sample=s,
                refinements=b["refinements"],
                floor=b["ends_floor"],
                cap=b["ends_cap"],
                proposals=b["proposals"],
                gains=b["gains"],
                halvings=b["halvings"],
                stage1=b["stage1_evals"],
                stage2=b["stage2_evals"],
                grid=f"{rows[s]['grid']['j5d']}/{rows[s]['grid']['j5db']}",
            )
        )
    fb = []
    for s in o:
        b = rows[s]["bf"]
        fb.append(
            dict(
                sample=s,
                steps=b["steps"],
                ball_found=b["ball_found"],
                fallback=b["fallback_steps"],
                fb_found=b["fallback_found"],
                fb_none=b["fallback_none"],
                waste=b["ball_waste_evals"],
                extra=b["fallback_extra_evals"],
                seed_test=f"{rows[s]['seed_test']['j5d']}/{rows[s]['seed_test']['j5db']}",
            )
        )
    keys = [
        "index",
        "rescue",
        "polish",
        "grid",
        "bfs_neighbour",
        "bfs_revisit",
        "screen_pass",
        "repair",
    ]
    return dict(
        ar5_bf_evals=DT.markdown_table(
            ev,
            ["sample", "s8", "j5sd", "j5sdb", "vs_j5sd", "vs_s8"],
            dict(s8=",d", j5sd=",d", j5sdb=",d", vs_j5sd=".2f", vs_s8=".2f"),
        ),
        ar5_bf_accuracy=DT.markdown_table(
            acc,
            [
                "sample",
                "w_s8",
                "w_j5sd",
                "w_j5sdb",
                "wilson",
                "fb_j5sd",
                "fb_s8",
                "med_j5sd",
                "med_j5sdb",
            ],
            dict(med_j5sd=".5f", med_j5sdb=".5f"),
        ),
        ar5_bf_decision=DT.markdown_table(
            dec,
            ["sample", "fewer", "r1", "net", "r2", "med", "r3", "house", "passes"],
            dict(fewer=".2f", med=".5f"),
        ),
        ar5_bf_fallback=DT.markdown_table(
            fb,
            [
                "sample",
                "steps",
                "ball_found",
                "fallback",
                "fb_found",
                "fb_none",
                "waste",
                "extra",
                "seed_test",
            ],
            dict(
                steps=",d",
                ball_found=",d",
                fallback=",d",
                fb_found=",d",
                fb_none=",d",
                waste=",d",
                extra=",d",
            ),
        ),
        ar5_bf_stage2=DT.markdown_table(
            bl,
            [
                "sample",
                "refinements",
                "floor",
                "cap",
                "proposals",
                "gains",
                "halvings",
                "stage1",
                "stage2",
                "grid",
            ],
            dict(
                refinements=",d",
                floor=",d",
                cap=",d",
                proposals=",d",
                gains=",d",
                halvings=",d",
                stage1=",d",
                stage2=",d",
            ),
        ),
        ar5_bf_stages=DT.markdown_table(stg, ["sample"] + keys, {k: ",d" for k in keys}),
    )


def main() -> int:
    marker = CACHE / BG.MARKER
    if not marker.exists():
        print("held7's runs have not ended (no bf_runs.done): nothing is read")
        return 2
    if marker.read_text().splitlines()[0] != "exit=0":
        print("the runs did not end cleanly:", marker.read_text().splitlines()[0])
        return 2
    rows = {s: BS.one(s) for s in ORDER}
    verdict = step_verdict(rows)
    pred = predictions(rows)
    out = dict(rows=rows, verdict=verdict, prediction_3_4=pred, rule=BR.JUDGE)
    (RD / "summary.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    DT.write_tables(RD / "tables.md", tables(rows))
    figure(rows, RD / "ar5_bf_evals_wrong.png")
    for s in ORDER:
        r = rows[s]
        print(
            f"{s:6s} J5SDB {r['totals']['j5sdb']:>11,} vs J5SD {r['totals']['j5sd']:>11,} "
            f"({100 * r['fewer']['j5sd']:.2f}% fewer) wrong {r['wrong']['j5sdb']['n']} vs "
            f"{r['wrong']['j5sd']['n']} net {r['pair']['j5sd']['net']:+d} med "
            f"{r['pair']['j5sd']['median_higher']:+.5f} | vs S8 {100 * r['fewer']['s8']:.1f}% "
            f"net {r['pair']['s8']['net']:+d} | passes {r['vs_j5sd']['ok'] and r['vs_s8']['ok']}"
        )
    print("step passes:", verdict["step_passes"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
