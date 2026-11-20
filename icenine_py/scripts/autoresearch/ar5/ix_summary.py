"""AR5 step 26: summary of gates e1, w1, e2, w2 (reads the exactness files and the live runs).

  uv run python scripts/autoresearch/ar5/ix_summary.py [--cache DIR]

Reads ``ar5/ix/exact_<sample>.json`` (``ix_exact.py``), the live runs ``<cache>/live_old`` and
``<cache>/live_new`` (``ix_run_all.sh --stage live``) and the stored J5SD rows of step 22 (linked
in ``<cache>/stored_J5SD_<sample>``, ``<cache>/J5D_<sample>``). Writes ``ar5/ix/summary.json``,
``ar5/ix/tables.md`` (``doc_tables.write_tables``) and the plot ``ar5/ix/ix_seconds.png``. A second
counts only in a gate when its run's label is ``single-worker`` (the preflight passed at the start
and the end); otherwise it is listed as contended and the gate is "not judged: machine not quiet".
No truth is read; no p-value is computed.
"""

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
RD = ICE / "benchmarks" / "autoresearch" / "ar5" / "ix"
import sys  # noqa: E402

sys.path.insert(0, str(ICE / "scripts" / "common"))
import doc_tables as DT  # noqa: E402

SAMPLES = ("mid", "held", "held3", "held4", "held5", "held6", "grad", "big")
STORED_EVALS = {"mid": 7383368, "held": 17114589}


def load(p: Path) -> Optional[Dict[str, Any]]:
    return json.loads(p.read_text()) if p.exists() else None


def live_row(cache: Path, arm: str, s: str) -> Optional[Dict[str, Any]]:
    d = cache / f"live_{arm}"
    run = load(d / f"J5D_{s}" / "run.json")
    rep = load(d / "j5sd" / f"J5S_{s}" / "repair.json")
    if run is None or rep is None:
        return None
    st = run["stats"]
    after = d / "j5sd" / f"J5S_{s}" / "after.npy"
    stored_after = cache / f"stored_J5SD_{s}" / "after.npy"
    equal = None
    if after.exists() and stored_after.exists():
        equal = bool(np.array_equal(np.load(after), np.load(stored_after)))
    label_f = cache / "logs" / f"live_{arm}_{s}.label"
    label = label_f.read_text().strip() if label_f.exists() else "unknown"
    return dict(
        sample=s,
        arm=arm,
        evals=int(run["evals_total"]) + int(rep["pass_evals"]) + int(rep["n_evals"]),
        wall_s=float(run["wall_total_s"]) + float(rep["wall_pass_s"]) + float(rep["wall_repair_s"]),
        index_s=float(st["index_s"]),
        n_seed_steps=int(st.get("n_seed_steps", 0)),
        map_equal_stored=equal,
        label=label,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", default=str(RD / "cache"))
    a = ap.parse_args()
    cache = Path(a.cache)
    exact = {
        s: load(RD / f"exact_{s}.json") or load(RD / f"exact_{s}_untimed.json") for s in SAMPLES
    }
    live: Dict[str, Dict[str, Any]] = {}
    for s in ("mid", "held"):
        for arm in ("old", "new"):
            r = live_row(cache, arm, s)
            if r:
                live[f"{arm}_{s}"] = r
    summ: Dict[str, Any] = dict(exact={k: v for k, v in exact.items() if v}, live=live)
    # e1
    e1_rows: List[Dict[str, Any]] = []
    for s in SAMPLES:
        x = exact[s]
        e1_rows.append(
            dict(
                sample=s,
                voxels=x["n_voxels"] if x else None,
                identical=x["tier1_identical"] if x else None,
                differ=x["n_differ"] if x else None,
            )
        )
    # w1 (mid and held)
    w1_rows: List[Dict[str, Any]] = []
    old_t = new_t = 0.0
    quiet_all = True
    for s in ("mid", "held"):
        x = exact[s]
        if not x:
            continue
        old_t += x["old_s"]
        new_t += x["new_s"]
        quiet_all = quiet_all and x["timing_label"] == "single-worker"
        w1_rows.append(
            dict(
                sample=s,
                voxels=x["n_voxels"],
                old_s=x["old_s"],
                new_s=x["new_s"],
                less_pct=x["less_pct"],
                label=x["timing_label"],
            )
        )
    summ["w1"] = dict(
        old_s=old_t,
        new_s=new_t,
        less_pct=100 * (1 - new_t / old_t) if old_t else None,
        quiet=quiet_all,
        gate="not judged: machine not quiet" if not quiet_all else None,
    )
    # e2 / w2
    w2_rows: List[Dict[str, Any]] = []
    for s in ("mid", "held"):
        o, n = live.get(f"old_{s}"), live.get(f"new_{s}")
        if not o or not n:
            continue
        w2_rows.append(
            dict(
                sample=s,
                evals_old=o["evals"],
                evals_new=n["evals"],
                stored=STORED_EVALS[s],
                map_new_equals_stored=n["map_equal_stored"],
                map_old_equals_stored=o["map_equal_stored"],
                wall_old=o["wall_s"],
                wall_new=n["wall_s"],
                less_pct=100 * (1 - n["wall_s"] / o["wall_s"]),
                index_old=o["index_s"],
                index_new=n["index_s"],
                label=f"{o['label']} / {n['label']}",
            )
        )
    pooled = None
    if len(w2_rows) == 2:
        po = sum(r["wall_old"] for r in w2_rows)
        pn = sum(r["wall_new"] for r in w2_rows)
        pooled = 100 * (1 - pn / po)
    summ["w2"] = dict(rows=w2_rows, pooled_less_pct=pooled)
    (RD / "summary.json").write_text(json.dumps(summ, indent=1) + "\n")
    tables = {
        "ix_e1": DT.markdown_table(
            e1_rows,
            ["sample", "voxels", "identical", "differ"],
            labels=["sample", "voxels", "bit-identical", "differ"],
        ),
        "ix_w1": DT.markdown_table(
            w1_rows,
            ["sample", "voxels", "old_s", "new_s", "less_pct", "label"],
            formats={"old_s": ".1f", "new_s": ".1f", "less_pct": ".1f"},
            labels=["sample", "voxels", "old s", "new s", "% less", "timing label"],
        ),
        "ix_w2": DT.markdown_table(
            w2_rows,
            [
                "sample",
                "evals_old",
                "evals_new",
                "map_new_equals_stored",
                "wall_old",
                "wall_new",
                "less_pct",
                "label",
            ],
            formats={
                "evals_old": ",d",
                "evals_new": ",d",
                "wall_old": ".0f",
                "wall_new": ".0f",
                "less_pct": ".1f",
            },
            labels=[
                "sample",
                "evals old",
                "evals new",
                "new map = stored",
                "wall old s",
                "wall new s",
                "% less",
                "labels (old / new)",
            ],
        ),
    }
    DT.write_tables(RD / "tables.md", tables)
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
        names = [r["sample"] for r in w1_rows]
        x = np.arange(len(names))
        if w1_rows:
            ax[0].bar(x - 0.2, [r["old_s"] for r in w1_rows], 0.4, label="old")
            ax[0].bar(x + 0.2, [r["new_s"] for r in w1_rows], 0.4, label="new")
            ax[0].set_xticks(x)
            ax[0].set_xticklabels(names)
        ax[0].set_ylabel("indexing seconds (all voxels)")
        ax[0].set_title("index only (w1)")
        ax[0].legend()
        if w2_rows:
            names2 = [r["sample"] for r in w2_rows]
            x2 = np.arange(len(names2))
            ax[1].bar(x2 - 0.2, [r["wall_old"] for r in w2_rows], 0.4, label="old")
            ax[1].bar(x2 + 0.2, [r["wall_new"] for r in w2_rows], 0.4, label="new")
            ax[1].set_xticks(x2)
            ax[1].set_xticklabels(names2)
        ax[1].set_ylabel("wall seconds (J5SD)")
        ax[1].set_title("whole run (w2)")
        fig.tight_layout()
        fig.savefig(RD / "ix_seconds.png", dpi=110)
    except Exception as e:  # noqa: BLE001
        print("plot skipped:", e)
    print(json.dumps({k: v for k, v in summ.items() if k in ("w1", "w2")}, indent=1))


if __name__ == "__main__":
    main()
