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

import json
import sys
from pathlib import Path
from typing import Any, Dict, List

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ball_guard as BG  # noqa: E402
import ball_replay_collect as BC  # noqa: E402
import ball_rules as BR  # noqa: E402
import doc_tables as DT  # noqa: E402

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
        info["samples"][s] = dict(replay=t, kinds=kind_counts(rows), where=where_counts(s))
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
    axs[0, 0].legend()
    fig.suptitle(f"Rule CB replay: 0.95 ends lost and gained, by step kind ({LABEL})")
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
        zn.append(dict(sample=s, lost=w["n"], **{z: w["zone"][z] for z in ZONES}))
        ds.append(dict(sample=s, lost=w["n"], **{b: w["dist"][b] for b in BINS}))
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
                ],
                dict(first_steps=",d", reseed_steps=",d", lost_pct=".2f"),
            ),
            ar5_ball_info_zone=DT.markdown_table(zn, ["sample", "lost"] + list(ZONES)),
            ar5_ball_info_dist=DT.markdown_table(ds, ["sample", "lost"] + list(BINS)),
        ),
    )
    for s in ORDER:
        S = info["samples"][s]
        print(s, S["kinds"], S["where"]["zone"], S["where"]["dist"], pred[s])
    return 0


if __name__ == "__main__":
    sys.exit(main())
