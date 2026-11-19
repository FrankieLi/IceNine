"""AR5 step 23, phase A steps 4 and 5: judge the replay runs (g0 on ``mid``, g1 on all eight) and
write the numbers of prediction 2.

  uv run python scripts/autoresearch/ar5/ball_replay_collect.py g0
  uv run python scripts/autoresearch/ar5/ball_replay_collect.py g1

Reads ``<cache>/replay_<mode>_<s>/rows.jsonl`` (written by ``ball_replay.py``; the replay runs of
``ball_replay_all.sh``). Arithmetic only: no image, no truth. Writes ``ar5/ball/g0.json``, or
``ar5/ball/replay.json`` with ``replay_tables.md`` (blocks for ``sync_doc_tables.py``). Exit 1 when
the gate fails (the plan's stop rule: no tuning, no ``held7``).

g0 (prediction 1, exactly, on ``mid``): with CB off every stored step is reproduced: every fit
equal (hit ratio, cost, end orientation, evaluations), the stop and the step's evaluations equal;
with CB on the first stage of every refinement equals ``J5D``'s stored first stage (the gains of its
trace up to the end of the first stage, evaluation by evaluation), and no stall check or pause
fires. g1: the replay on all eight samples: the index stage (the sum of the seed steps'
evaluations) at least 30% fewer than ``J5D``'s on the held types and ``mid``, and the steps that
lose their 0.95 end at most 3% of the replayed steps (the denominator is every replayed step; the
share of steps that had an end is reported too).
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICE / "scripts" / "common"))
import ball_guard as BG  # noqa: E402
import ball_rules as BR  # noqa: E402
import doc_tables as DT  # noqa: E402

CACHE = BG.base_cache()
RD = ICE / "benchmarks" / "autoresearch" / "ar5" / "ball"
ORDER = ("held5", "held4", "held3", "held", "mid", "grad", "big", "held6")
J5D_REFINEMENT_MEAN = (793, 1341)  # the plan's range of J5D's evaluations per refinement


def read_rows(mode: str, s: str) -> List[Dict[str, Any]]:
    p = CACHE / f"replay_{mode}_{s}" / "rows.jsonl"
    done = CACHE / f"replay_{mode}_{s}" / "replay.done"
    if not done.exists():
        raise SystemExit(f"{s} {mode}: the replay has not finished (no replay.done)")
    return [json.loads(x) for x in BG.guarded(p).read_text().splitlines() if x]


def off_row_ok(r: Dict[str, Any]) -> bool:
    return bool(r["fits_equal"] and r["evals"] == r["j5d_evals"] and r["stop"] == r["j5d_stop"])


def g0() -> int:
    off = read_rows("off", "mid")
    on = read_rows("on", "mid")
    n_cmp = sum(r["n_first_stage_compared"] for r in on)
    n_ref = sum(len(r["refinements"]) for r in on)
    out: Dict[str, Any] = dict(
        sample="mid",
        off_steps=len(off),
        off_fits_equal=sum(1 for r in off if r["fits_equal"]),
        off_evals_equal=sum(1 for r in off if r["evals"] == r["j5d_evals"]),
        off_stop_equal=sum(1 for r in off if r["stop"] == r["j5d_stop"]),
        off_ok=bool(all(off_row_ok(r) for r in off) and off),
        on_steps=len(on),
        on_refinements=n_ref,
        on_first_stage_compared=n_cmp,
        on_first_stage_equal=sum(
            1 for r in on for x in r["refinements"] if x["first_stage_equal"] is True
        ),
        on_first_stage_unequal=sum(
            1 for r in on for x in r["refinements"] if x["first_stage_equal"] is False
        ),
        on_not_compared=n_ref - n_cmp,
        on_checks=sum(r["checks"] for r in on),
        on_ok=False,
    )
    out["on_ok"] = bool(
        n_ref > 0
        and out["on_first_stage_unequal"] == 0
        and out["on_checks"] == 0
        and out["on_first_stage_equal"] == n_cmp
    )
    out["g0_passes"] = bool(out["off_ok"] and out["on_ok"])
    out["note"] = (
        "on_not_compared: refinements the stored J5D run did not make (a step whose end moved, "
        "so another candidate was refined); their first stage has no stored counterpart"
    )
    (RD / "g0.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    return 0 if out["g0_passes"] else 1


def one(s: str) -> Dict[str, Any]:
    on = read_rows("on", s)
    n = len(on)
    refs = [x for r in on for x in r["refinements"]]
    j5d_ev = sum(r["j5d_evals"] for r in on)
    ev = sum(r["evals"] for r in on)
    j5d_end = sum(1 for r in on if r["j5d_stop"] >= 0)
    end = sum(1 for r in on if r["stop"] >= 0)
    lost = sum(1 for r in on if r["lost_end"])
    gained = sum(1 for r in on if r["gained_end"])
    moved = [r for r in on if r["taken_moved_deg"] > 0.05]
    row: Dict[str, Any] = dict(
        sample=s,
        steps=n,
        j5d_evals=j5d_ev,
        evals=ev,
        fewer=BR.fewer(ev, j5d_ev),
        refinements=len(refs),
        mean_refinement_evals=sum(x["evals"] for x in refs) / max(1, len(refs)),
        floor_share=sum(1 for x in refs if x["end"] == "floor") / max(1, len(refs)),
        cap_ends=sum(1 for x in refs if x["end"] == "cap"),
        mean_proposals=sum(x["proposals"] for x in refs) / max(1, len(refs)),
        mean_halvings=sum(x["halvings"] for x in refs) / max(1, len(refs)),
        j5d_end_steps=j5d_end,
        end_steps=end,
        end_steps_rel=(end - j5d_end) / max(1, j5d_end),
        lost=lost,
        lost_share=lost / max(1, n),
        lost_share_of_ended=lost / max(1, j5d_end),
        gained=gained,
        tried_mean=sum(r["tried"] for r in on) / max(1, n),
        j5d_tried_mean=sum(r["j5d_tried"] for r in on) / max(1, n),
        costing=sum(r["costed"] for r in on),
        cma=BR.CMA_EVALS * sum(r["n_cma"] for r in on),
        stage1=sum(x["stage1"] for x in refs),
        stage2=sum(x["stage2"] for x in refs),
        final=sum(x["final"] for x in refs),
        first_stage_unequal=sum(1 for x in refs if x["first_stage_equal"] is False),
        checks=sum(r["checks"] for r in on),
        moved_over_0_05=len(moved),
        moved_list=[
            dict(
                index=r["index"],
                voxel=r["voxel"],
                moved_deg=round(r["taken_moved_deg"], 4),
                j5d_taken=r["j5d_taken"],
                taken=r["taken"],
            )
            for r in moved[:60]
        ],
    )
    row["split_sum_ok"] = bool(
        row["costing"] + row["cma"] + row["stage1"] + row["stage2"] + row["final"] == ev
    )
    return row


def inside(x: float, lo: float, hi: float) -> bool:
    return bool(lo <= x <= hi)


def main() -> int:
    if len(sys.argv) > 1 and sys.argv[1] == "g0":
        return g0()
    rows = {s: one(s) for s in ORDER}
    g = BR.GATE
    g1: Dict[str, Any] = {}
    for s, r in rows.items():
        g1[s] = dict(
            fewer=r["fewer"],
            fewer_ok=(r["fewer"] >= g["g1_fewer_min"]) if s in g["g1_fewer_samples"] else None,
            lost_share=r["lost_share"],
            lost_ok=bool(r["lost_share"] <= g["g1_lost_max"]),
        )
    ok = all(
        (v["fewer_ok"] is not False)
        and v["lost_ok"]
        and rows[s]["checks"] == 0
        and rows[s]["first_stage_unequal"] == 0
        for s, v in g1.items()
    )
    pred: Dict[str, Any] = {}
    for s, r in rows.items():
        lo, hi = (
            (0.40, 0.60)
            if s in ("held5", "held4", "held3", "held", "mid", "held6")
            else (0.20, 0.45)
        )
        pred[s] = dict(
            fewer_inside=inside(r["fewer"], lo, hi),
            band=[lo, hi],
            per_refinement_inside=inside(r["mean_refinement_evals"], 250, 450),
            floor_inside=bool(r["floor_share"] >= 0.90),
            end_steps_inside=bool(abs(r["end_steps_rel"]) <= 0.03),
        )
    out = dict(rows=rows, g1=g1, g1_passes=bool(ok), prediction_2=pred, gate=BR.GATE)
    (RD / "replay.json").write_text(json.dumps(out, indent=1) + "\n")
    cols = [
        "sample",
        "steps",
        "j5d_evals",
        "evals",
        "fewer",
        "mean_ref",
        "floor",
        "j5d_end",
        "end",
        "lost",
        "lost_pct",
        "gained",
        "moved",
    ]
    tab = [
        dict(
            sample=s,
            steps=r["steps"],
            j5d_evals=r["j5d_evals"],
            evals=r["evals"],
            fewer=100 * r["fewer"],
            mean_ref=r["mean_refinement_evals"],
            floor=100 * r["floor_share"],
            j5d_end=r["j5d_end_steps"],
            end=r["end_steps"],
            lost=r["lost"],
            lost_pct=100 * r["lost_share"],
            gained=r["gained"],
            moved=r["moved_over_0_05"],
        )
        for s, r in rows.items()
    ]
    cols2 = ["sample", "costing", "cma", "stage1", "stage2", "final", "tried", "j5d_tried"]
    tab2 = [
        dict(
            sample=s,
            costing=r["costing"],
            cma=r["cma"],
            stage1=r["stage1"],
            stage2=r["stage2"],
            final=r["final"],
            tried=r["tried_mean"],
            j5d_tried=r["j5d_tried_mean"],
        )
        for s, r in rows.items()
    ]
    DT.write_tables(
        RD / "replay_tables.md",
        dict(
            ar5_ball_replay=DT.markdown_table(
                tab,
                cols,
                dict(
                    steps=",d",
                    j5d_evals=",d",
                    evals=",d",
                    fewer=".1f",
                    mean_ref=".0f",
                    floor=".1f",
                    lost_pct=".2f",
                ),
            ),
            ar5_ball_replay_split=DT.markdown_table(
                tab2,
                cols2,
                dict(
                    costing=",d",
                    cma=",d",
                    stage1=",d",
                    stage2=",d",
                    final=",d",
                    tried=".3f",
                    j5d_tried=".3f",
                ),
            ),
        ),
    )
    for s, r in rows.items():
        print(
            f"{s:6s} steps {r['steps']:5d} fewer {100 * r['fewer']:.1f}% mean/refinement "
            f"{r['mean_refinement_evals']:.0f} floor {100 * r['floor_share']:.1f}% lost "
            f"{r['lost']} ({100 * r['lost_share']:.2f}%) gained {r['gained']} moved "
            f"{r['moved_over_0_05']}"
        )
    print("g1 passes:", ok)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
