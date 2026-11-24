"""AR5 step 23, information only (the step failed its gate g1): the replay of all eight disclosed
samples under rule ``CB``, the lost and gained 0.95 ends by step kind, and where the lost ends were.

  uv run python scripts/autoresearch/ar5/ball_info.py

Reads ``<cache>/replay_on_<s>/rows.jsonl`` (``ball_replay.py``) and ``<cache>/where_<s>.json``
(``ball_where.py``). Writes ``ar5/ball/info.json``, ``info_tables.md`` (blocks for
``sync_doc_tables.py``) and the figure ``ar5_ball_lost_gained.png`` (lost and gained ends by step
kind; the plot first). Every number is labelled disclosed, for information; the step failed its
gate. Step kinds: first seeds (``reseed`` False) and reseed steps (``reseed`` True). Counts only; no
tests of significance.
"""

import gzip
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ball_guard as BG  # noqa: E402
import ball_replay_collect as BC  # noqa: E402
import ball_rules as BR  # noqa: E402
import doc_tables as DT  # noqa: E402
import stats as ST  # noqa: E402

RD = BC.RD
CACHE = BG.base_cache()
ORDER = BC.ORDER
ZONES = ("stage1", "stage2_first_block", "stage2_later", "candidate_not_refined_by_CB")
BINS = ("0-0.1", "0.1-0.33", "0.33-0.64", "over 0.64")
LABEL = "disclosed, for information; the step failed its gate"


def kind_counts(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for kind, flag in (("first", False), ("reseed", True)):
        sub = [r for r in rows if bool(r["reseed"]) == flag]
        out[kind] = dict(
            steps=len(sub),
            j5d_end=sum(1 for r in sub if r["j5d_stop"] >= 0),
            lost=sum(1 for r in sub if r["lost_end"]),
            gained=sum(1 for r in sub if r["gained_end"]),
        )
    return out


def where_counts(s: str) -> Dict[str, Any]:
    d = json.loads(BG.guarded(CACHE / f"where_{s}.json").read_text())["rows"]
    zone = {z: sum(1 for r in d if r["zone"] == z) for z in ZONES}
    dist = {b: sum(1 for r in d if r.get("dist_bin") == b) for b in BINS}
    cross = {
        z: {b: sum(1 for r in d if r["zone"] == z and r.get("dist_bin") == b) for b in BINS}
        for z in ZONES[:3]
    }
    return dict(n=len(d), zone=zone, dist=dist, cross=cross)


def mid_lost_facts() -> Dict[str, Any]:
    """The facts about the 70 lost ``mid`` steps (from the rows and J5D's stored e1_log)."""
    rows = [
        json.loads(x)
        for x in BG.guarded(CACHE / "replay_on_mid" / "rows.jsonl").read_text().splitlines()
        if x
    ]
    with gzip.open(BG.guarded(CACHE / "J5D_mid" / "e1_log.json.gz"), "rt") as fh:
        ent = json.load(fh)
    lost = [r for r in rows if r["lost_end"]]
    hits = []
    for r in lost:
        e = ent[r["index"]]
        hits.append(
            min(f[3] for f in e["fits"] if f[0] == e["stop"] and f[2] == "M" and f[3] >= 0.95)
        )
    hits.sort()
    n = len(hits)
    med = hits[n // 2] if n % 2 else 0.5 * (hits[n // 2 - 1] + hits[n // 2])
    return dict(
        lost=n,
        median_stored_end_hit_ratio=float(med),
        min_stored_end_hit_ratio=float(hits[0]),
        stored_end_under_0_97=sum(1 for h in hits if h < 0.97),
        stored_end_is_mc_fit=n,
        reseed_steps=sum(1 for r in lost if r["reseed"]),
        j5d_steps_with_a_pause=sum(1 for r in lost if r["j5d_paused"] > 0),
        j5d_evals=sum(r["j5d_evals"] for r in lost),
        cb_evals=sum(r["evals"] for r in lost),
    )


def main() -> int:
    info: Dict[str, Any] = dict(label=LABEL, samples={})
    for s in ORDER:
        rows = [
            json.loads(x)
            for x in BG.guarded(CACHE / f"replay_on_{s}" / "rows.jsonl").read_text().splitlines()
            if x
        ]
        one = BC.one(s)
        t = dict(one)
        t.pop("moved_list", None)
        w = where_counts(s)
        lo, hi = ST.wilson(int(t["lost"]), int(t["steps"]))
        t["lost_share_wilson95"] = [float(lo), float(hi)]
        summ = json.loads(BG.guarded(CACHE / f"replay_on_{s}" / "summary.json").read_text())
        t["pauses"] = int(summ["stats"]["pauses"])
        w["over_0_33"] = int(w["dist"]["0.33-0.64"] + w["dist"]["over 0.64"])
        w["over_0_33_share"] = w["over_0_33"] / max(1, w["n"])
        info["samples"][s] = dict(replay=t, kinds=kind_counts(rows), where=w)
    info["totals"] = dict(
        steps=sum(info["samples"][s]["replay"]["steps"] for s in ORDER),
        lost=sum(info["samples"][s]["where"]["n"] for s in ORDER),
        lost_over_0_33=sum(info["samples"][s]["where"]["over_0_33"] for s in ORDER),
        pauses=sum(info["samples"][s]["replay"]["pauses"] for s in ORDER),
    )
    info["mid_lost_facts"] = mid_lost_facts()
    pp = RD / "replay_partial.json"
    part = json.loads(pp.read_text())
    part["mid_lost_facts"] = info["mid_lost_facts"]
    pp.write_text(json.dumps(part, indent=1) + "\n")
    # predictions of the plan (prediction 2), hit or missed
    pred: Dict[str, Any] = {}
    for s in ORDER:
        r = info["samples"][s]["replay"]
        lo, hi = (
            (0.40, 0.60)
            if s in ("held5", "held4", "held3", "held", "mid", "held6")
            else (0.20, 0.45)
        )
        pred[s] = {
            "index stage fewer": (lo <= r["fewer"] <= hi),
            "250-450 per refinement": (250 <= r["mean_refinement_evals"] <= 450),
            "floor ends >= 90%": (r["floor_share"] >= 0.90),
            "end steps within 3% of J5D": (abs(r["end_steps_rel"]) <= 0.03),
            "g1: fewer >= 30% (held types, mid)": (
                r["fewer"] >= 0.30 if s in BR.GATE["g1_fewer_samples"] else None
            ),
            "g1: lost <= 3% of steps": (r["lost_share"] <= 0.03),
        }
    info["prediction_2"] = pred
    (RD / "info.json").write_text(json.dumps(info, indent=1) + "\n")

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    import numpy as np

    fig, axs = plt.subplots(2, 4, figsize=(14, 6.5))
    for ax, s in zip(axs.ravel(), ORDER):
        k = info["samples"][s]["kinds"]
        x = np.arange(2)
        lost = [k["first"]["lost"], k["reseed"]["lost"]]
        gain = [k["first"]["gained"], k["reseed"]["gained"]]
        ax.bar(x - 0.18, lost, 0.36, color="tab:red", label="lost end")
        ax.bar(x + 0.18, gain, 0.36, color="tab:green", label="gained end")
        ax.set_xticks(x)
        ax.set_xticklabels(
            [f"first seeds\n({k['first']['steps']})", f"reseed\n({k['reseed']['steps']})"]
        )
        ax.set_title(s)
        ax.set_ylabel("steps")
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    axs[0, 0].legend()
    fig.suptitle(
        f"Rule CB replay: 0.95 ends lost and gained, by step kind ({LABEL});"
        " each panel has its own y scale"
    )
    fig.tight_layout()
    fig.savefig(RD / "ar5_ball_lost_gained.png", dpi=110)
    plt.close(fig)

    ev, kd, zn, ds = [], [], [], []
    for s in ORDER:
        S = info["samples"][s]
        r, k, w = S["replay"], S["kinds"], S["where"]
        ev.append(
            dict(
                sample=s,
                steps=r["steps"],
                j5d_evals=r["j5d_evals"],
                cb_evals=r["evals"],
                fewer=100 * r["fewer"],
                per_ref=r["mean_refinement_evals"],
                floor=100 * r["floor_share"],
                j5d_end=r["j5d_end_steps"],
                cb_end=r["end_steps"],
            )
        )
        kd.append(
            dict(
                sample=s,
                first_steps=k["first"]["steps"],
                first_lost=k["first"]["lost"],
                first_gained=k["first"]["gained"],
                reseed_steps=k["reseed"]["steps"],
                reseed_lost=k["reseed"]["lost"],
                reseed_gained=k["reseed"]["gained"],
                lost_pct=100 * r["lost_share"],
            )
        )
        lo, hi = r["lost_share_wilson95"]
        kd[-1]["wilson95"] = f"{100 * lo:.2f}-{100 * hi:.2f}%"
        zn.append(dict(sample=s, lost=w["n"], **{z: w["zone"][z] for z in ZONES}))
        ds.append(
            dict(
                sample=s,
                lost=w["n"],
                **{b: w["dist"][b] for b in BINS},
                over_0_33=w["over_0_33"],
                share=100 * w["over_0_33_share"],
            )
        )
    T = info["totals"]
    ev.append(dict(sample="all eight", steps=T["steps"]))
    kd.append(dict(sample="all eight", lost_pct=100 * T["lost"] / T["steps"]))
    zn.append(dict(sample="all eight", lost=T["lost"]))
    ds.append(
        dict(
            sample="all eight",
            lost=T["lost"],
            over_0_33=T["lost_over_0_33"],
            share=100 * T["lost_over_0_33"] / T["lost"],
        )
    )
    DT.write_tables(
        RD / "info_tables.md",
        dict(
            ar5_ball_info_evals=DT.markdown_table(
                ev,
                [
                    "sample",
                    "steps",
                    "j5d_evals",
                    "cb_evals",
                    "fewer",
                    "per_ref",
                    "floor",
                    "j5d_end",
                    "cb_end",
                ],
                dict(
                    steps=",d",
                    j5d_end=",d",
                    cb_end=",d",
                    j5d_evals=",d",
                    cb_evals=",d",
                    fewer=".1f",
                    per_ref=".0f",
                    floor=".1f",
                ),
            ),
            ar5_ball_info_kinds=DT.markdown_table(
                kd,
                [
                    "sample",
                    "first_steps",
                    "first_lost",
                    "first_gained",
                    "reseed_steps",
                    "reseed_lost",
                    "reseed_gained",
                    "lost_pct",
                    "wilson95",
                ],
                dict(first_steps=",d", reseed_steps=",d", lost_pct=".2f"),
            ),
            ar5_ball_info_zone=DT.markdown_table(zn, ["sample", "lost"] + list(ZONES)),
            ar5_ball_info_dist=DT.markdown_table(
                ds, ["sample", "lost"] + list(BINS) + ["over_0_33", "share"], dict(share=".1f")
            ),
        ),
    )
    for s in ORDER:
        S = info["samples"][s]
        print(s, S["kinds"], S["where"]["zone"], S["where"]["dist"], pred[s])
    return 0


if __name__ == "__main__":
    sys.exit(main())
