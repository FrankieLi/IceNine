"""Item 14 and 3, phase B: summary.json, tables.md and item14_time.png from the cache.

  uv run python scripts/autoresearch/item14/summary.py

Reads benchmarks/autoresearch/item14/cache (gitignored: E0 arms, runs A1/A2/B1/B2, quiet verdicts,
pair checks, split) and the stored arms (E0 ``e0_old_base`` and ``e0_old_e34t25``; pilot
``keyoff_rerun`` and ``smooth_only``). Writes to benchmarks/autoresearch/item14/. Wrong = error over
1 degree. No p-values: Wilson intervals and paired fixed / broken counts only.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
sys.path.insert(0, str(ICE / "scripts" / "common"))
sys.path.insert(0, str(ICE / "scripts" / "eval_reduction"))

import doc_tables as DT  # noqa: E402
import stats as ST  # noqa: E402

OUT = ICE / "benchmarks" / "autoresearch" / "item14"
CACHE = OUT / "cache"
E_STORED = ICE / "benchmarks" / "eval_reduction" / "E" / "cache" / "rank"
G = ICE / "benchmarks" / "eval_reduction" / "G"
PILOT_STORED = {
    "off": G / "g2" / "bfs" / "runs" / "keyoff_rerun",
    "S": G / "g3" / "bfs" / "runs" / "smooth_only",
}
WRONG = 1.0
ARMS = ["off", "F", "S", "FS"]


# ----------------------------------------------------------------------------------- E0
def e0_arm(path: Path) -> Dict[int, Dict[str, float]]:
    out = {}
    for f in sorted(path.glob("v*.npz")):
        d = np.load(f)
        out[int(d["voxel"])] = dict(
            err=float(d["err_deg"]), evals=int(d["evals"]), smooth=int(d["smooth_evals"])
        )
    return out


def paired(a: List[bool], b: List[bool]) -> Dict[str, int]:
    fixed, broken = ST.paired_discordant(a, b)  # wrong in A only, wrong in B only
    return dict(fixed=fixed, broken=broken)


def e0_summary() -> Dict[str, Any]:
    dirs = {
        "off": E_STORED / "e0_old_base",
        "S": E_STORED / "e0_old_e34t25",
        "F": CACHE / "e0" / "e0_old_off_fused",
        "FS": CACHE / "e0" / "e0_old_e34t25_fused",
    }
    arms = {k: e0_arm(v) for k, v in dirs.items()}
    ids = sorted(set.intersection(*[set(a) for a in arms.values()]))
    res: Dict[str, Any] = dict(n=len(ids), arms={})
    wrong = {}
    for k, a in arms.items():
        w = [a[v]["err"] > WRONG for v in ids]
        wrong[k] = w
        lo, hi = ST.wilson(sum(w), len(ids))
        res["arms"][k] = dict(
            wrong=int(sum(w)), lo=lo, hi=hi,
            evals_per_voxel=float(np.mean([a[v]["evals"] for v in ids])),
            smooth_per_voxel=float(np.mean([a[v]["smooth"] for v in ids])),
        )  # fmt: skip
    res["pairs"] = {
        "F vs off": paired(wrong["off"], wrong["F"]),
        "FS vs S": paired(wrong["S"], wrong["FS"]),
        "FS vs F": paired(wrong["F"], wrong["FS"]),
        "S vs off": paired(wrong["off"], wrong["S"]),
    }
    res["evals_ratio_F_vs_off"] = res["arms"]["F"]["evals_per_voxel"] / (
        res["arms"]["off"]["evals_per_voxel"]
    )
    return res


# --------------------------------------------------------------------------------- pilot
def run_dirs() -> Dict[str, Path]:
    return {
        "off": PILOT_STORED["off"],
        "S": PILOT_STORED["S"],
        "F": CACHE / "runs" / "A1",
        "FS": CACHE / "runs" / "B1",
    }


def pilot_summary() -> Dict[str, Any]:
    import g2_bfs

    errs, runs = {}, {}
    for k, d in run_dirs().items():
        errs[k] = g2_bfs.errors(d)["err"]
        runs[k] = json.loads((d / "run.json").read_text())
    n = len(errs["off"])
    res: Dict[str, Any] = dict(n=n, arms={})
    wrong = {k: list(errs[k] > WRONG) for k in ARMS}
    for k in ARMS:
        st = runs[k]["stats"]
        lo, hi = ST.wilson(sum(wrong[k]), n)
        sm = runs[k].get("smooth_evals_per_full_search", [])
        res["arms"][k] = dict(
            wrong=int(sum(wrong[k])), lo=lo, hi=hi,
            full_searches=int(st["n_seeds"] + st["n_reseeds"]),
            evals_full=int(st["n_evals_seed"] + st["n_evals_reseed"]),
            evals_local=int(st["n_evals_neighbor"] + st["n_evals_revisit"]),
            evals_total=int(
                st["n_evals_seed"] + st["n_evals_reseed"] + st["n_evals_neighbor"]
                + st["n_evals_revisit"]
            ),
            smooth_evals=int(sum(sm)),
            smooth_per_full_search=float(np.mean(sm)) if sm else 0.0,
        )  # fmt: skip
    res["pairs"] = {
        "FS vs S": paired(wrong["S"], wrong["FS"]),
        "FS vs F": paired(wrong["F"], wrong["FS"]),
        "F vs off": paired(wrong["off"], wrong["F"]),
        "S vs off": paired(wrong["off"], wrong["S"]),
    }
    res["r3_wrong_A1_minus_B1"] = res["arms"]["F"]["wrong"] - res["arms"]["FS"]["wrong"]
    return res


# --------------------------------------------------------------------------------- timing
def real_s(name: str) -> float:
    for line in (CACHE / "logs" / f"{name}.time").read_text().splitlines():
        if line.startswith("real"):
            return float(line.split()[1])
    raise ValueError(name)


def timing_summary() -> Dict[str, Any]:
    runs: Dict[str, Any] = {}
    for n in ("A1", "A2", "B1", "B2"):
        r = json.loads((CACHE / "runs" / n / "run.json").read_text())
        st = r["stats"]
        q = json.loads((CACHE / "quiet" / f"{n}.json").read_text())
        pre = json.loads((CACHE / f"preflight_{n}.json").read_text())
        runs[n] = dict(
            real_s=real_s(n), wall_total_s=float(r["wall_total_s"]),
            full_s=float(st["wall_seed_s"] + st["wall_reseed_s"]),
            local_s=float(st["wall_neighbor_s"] + st["wall_revisit_s"]),
            label=q["label"], strict_quiet=bool(q["quiet"]), relaxed_ok=bool(q["relaxed_ok"]),
            max_load=float(q["max_load"]), median_load=float(q["median_load"]),
            n_load=int(q["n_load"]), busy=[b["command"].split("/")[-1] for b in q["busy"]],
            preflight_load=float(pre["loadavg"][0]),
        )  # fmt: skip
    pairs = {a: json.loads((CACHE / f"pair_{a}.json").read_text()) for a in ("A", "B")}
    arm: Dict[str, Any] = {}
    for a, (x, y) in {"A": ("A1", "A2"), "B": ("B1", "B2")}.items():
        for key in ("real_s", "wall_total_s", "full_s", "local_s"):
            v = [runs[x][key], runs[y][key]]
            arm.setdefault(a, {})[key] = dict(median=float(np.median(v)), lo=min(v), hi=max(v))
    ratio = arm["B"]["wall_total_s"]["median"] / arm["A"]["wall_total_s"]["median"]
    sp = json.loads((CACHE / "quiet" / "split.json").read_text())
    split = json.loads((CACHE / "split" / "split.json").read_text())
    return dict(
        runs=runs, pairs=pairs, arm=arm, ratio_B_over_A=ratio,
        smooth_share=(ratio - 1) / ratio,
        split_quiet=dict(strict=bool(sp["quiet"]), max_load=sp["max_load"],
                         median_load=sp["median_load"]),
        split=[
            dict(voxel=r["voxel"], arm=r["arm"], wall_s=r["wall_s"], shares=r["shares"],
                 evals=r["evals"]["total"], smooth_evals=r["smooth_evals"],
                 err_deg=r["err_deg"],
                 all_strict_fused=bool(r["evaluate_fused_n"] == r["evaluate_n"]))
            for r in split["runs"]
        ],
        e0_walls={
            k: json.loads((CACHE / "e0" / d / "run_meta.json").read_text())["wall_s"]
            for k, d in {"F": "e0_old_off_fused", "FS": "e0_old_e34t25_fused"}.items()
        },
    )  # fmt: skip


# ---------------------------------------------------------------------------------- tables
def fmt_ci(a: Dict[str, Any]) -> str:
    return f"{a['wrong']} ({a['lo'] * 100:.1f}% to {a['hi'] * 100:.1f}%)"


def tables(s: Dict[str, Any]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    e, p, t = s["e0"], s["pilot"], s["timing"]
    rows = [
        dict(arm=k, wrong=fmt_ci(e["arms"][k]),
             evals=e["arms"][k]["evals_per_voxel"], smooth=e["arms"][k]["smooth_per_voxel"])
        for k in ARMS
    ]  # fmt: skip
    out["e0_arms"] = (
        DT.markdown_table(
            rows,
            ["arm", "wrong", "evals", "smooth"],
            {"evals": ".0f", "smooth": ".0f"},
        )
        .replace("| evals |", "| strict evaluations per voxel |")
        .replace("| smooth |", "| smooth evaluations per voxel |")
    )
    rows = [dict(pair=k, fixed=v["fixed"], broken=v["broken"]) for k, v in e["pairs"].items()]
    out["e0_pairs"] = DT.markdown_table(rows, ["pair", "fixed", "broken"])
    rows = [
        dict(arm=k, wrong=fmt_ci(p["arms"][k]), full=p["arms"][k]["full_searches"],
             ev=p["arms"][k]["evals_total"], full_ev=p["arms"][k]["evals_full"],
             loc=p["arms"][k]["evals_local"], sm=p["arms"][k]["smooth_evals"])
        for k in ARMS
    ]  # fmt: skip
    out["pilot_arms"] = (
        DT.markdown_table(
            rows,
            ["arm", "wrong", "full", "ev", "full_ev", "loc", "sm"],
            {c: "," for c in ("ev", "full_ev", "loc", "sm")},
        )
        .replace("| full |", "| full searches |")
        .replace("| ev |", "| strict evaluations |")
        .replace("| full_ev |", "| of these, full searches |")
        .replace("| loc |", "| of these, local optimization |")
        .replace("| sm |", "| smooth evaluations |")
    )
    rows = [dict(pair=k, fixed=v["fixed"], broken=v["broken"]) for k, v in p["pairs"].items()]
    out["pilot_pairs"] = DT.markdown_table(rows, ["pair", "fixed", "broken"])
    rows = []
    for n, r in t["runs"].items():
        rows.append(
            dict(run=n, arm="A (F)" if n[0] == "A" else "B (FS)", real=r["real_s"],
                 total=r["wall_total_s"], full=r["full_s"], loc=r["local_s"],
                 label=r["label"] or "strict quiet check failed ()")
        )  # fmt: skip
    out["timing_runs"] = (
        DT.markdown_table(
            rows,
            ["run", "arm", "real", "total", "full", "loc", "label"],
            {c: ".1f" for c in ("real", "total", "full", "loc")},
        )
        .replace("| real |", "| whole process (s) |")
        .replace("| total |", "| BFS wall_total_s |")
        .replace("| full |", "| full searches (s) |")
        .replace("| loc |", "| local optimization (s) |")
    )
    rows = [
        dict(run=n, mx=r["max_load"], med=r["median_load"], n=r["n_load"],
             pre=r["preflight_load"], busy=", ".join(r["busy"]) or "none caught")
        for n, r in t["runs"].items()
    ]  # fmt: skip
    out["quiet_runs"] = (
        DT.markdown_table(
            rows,
            ["run", "mx", "med", "n", "pre", "busy"],
            {"mx": ".2f", "med": ".2f", "pre": ".2f"},
        )
        .replace("| mx |", "| max load |")
        .replace("| med |", "| median load |")
        .replace("| n |", "| load samples |")
        .replace("| pre |", "| preflight load |")
        .replace("| busy |", "| busy processes caught |")
    )
    rows = []
    for a, pr in t["pairs"].items():
        rows.append(
            dict(arm=a, t1=pr["real_s"][0], t2=pr["real_s"][1], d=pr["rel_diff"] * 100,
                 med=t["arm"][a]["wall_total_s"]["median"])
        )  # fmt: skip
    out["pairs_time"] = (
        DT.markdown_table(
            rows,
            ["arm", "t1", "t2", "d", "med"],
            {"t1": ".2f", "t2": ".2f", "d": ".2f", "med": ".1f"},
        )
        .replace("| t1 |", "| first run, whole process (s) |")
        .replace("| t2 |", "| repeat (s) |")
        .replace("| d |", "| difference (%) |")
        .replace("| med |", "| median BFS wall_total_s (s) |")
    )
    rows = [
        dict(voxel=r["voxel"], arm=r["arm"], wall=r["wall_s"], screen=r["shares"]["screen"] * 100,
             smooth=r["shares"]["smooth"] * 100, strict=r["shares"]["strict"] * 100,
             rest=r["shares"]["rest"] * 100, ev=r["evals"], sm=r["smooth_evals"])
        for r in t["split"]
    ]  # fmt: skip
    out["split"] = (
        DT.markdown_table(
            rows,
            ["voxel", "arm", "wall", "screen", "smooth", "strict", "rest", "ev", "sm"],
            {c: ".1f" for c in ("wall", "screen", "smooth", "strict", "rest")}
            | {"ev": ",", "sm": ","},
        )
        .replace("| wall |", "| wall (s) |")
        .replace("| screen |", "| screen (%) |")
        .replace("| smooth |", "| smooth (%) |")
        .replace("| strict |", "| strict (%) |")
        .replace("| rest |", "| rest (%) |")
        .replace("| ev |", "| strict evaluations |")
        .replace("| sm |", "| smooth evaluations |")
    )
    return out


# ----------------------------------------------------------------------------------- figure
def figure(s: Dict[str, Any], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t = s["timing"]
    fig, (ax, bx, cx) = plt.subplots(
        1, 3, figsize=(13, 4.2), gridspec_kw=dict(width_ratios=[1, 1, 1])
    )
    for i, a in enumerate(("A", "B")):
        full, loc = t["arm"][a]["full_s"]["median"], t["arm"][a]["local_s"]["median"]
        ax.bar(i, full, color="#4c72b0", label="full searches" if i == 0 else None)
        ax.bar(i, loc, bottom=full, color="#dd8452", label="local optimization" if i == 0 else None)
        names = ("A1", "A2") if a == "A" else ("B1", "B2")
        for k, n in enumerate(names):
            ax.plot(i + (k - 0.5) * 0.15, t["runs"][n]["wall_total_s"], "ko", ms=5)
        ax.text(i, t["arm"][a]["wall_total_s"]["median"] + 80,
                f"{t['arm'][a]['wall_total_s']['median']:.0f} s", ha="center")  # fmt: skip
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["A: FusedCost", "B: FusedCost + SmoothCostRank"])
    ax.set_ylabel("pilot BFS wall_total_s (s)")
    ax.set_title("Pilot time, 2 repeats (dots)\nsingle-worker, R4b labels (see page)", fontsize=9)
    ax.legend(fontsize=8)
    for a_, key, title in (
        (bx, "e0", "E0: wrong voxels of 200"),
        (cx, "pilot", "Pilot: wrong of 2000"),
    ):
        d = s[key]
        w = [d["arms"][k]["wrong"] for k in ARMS]
        lo = [d["arms"][k]["wrong"] - d["arms"][k]["lo"] * d["n"] for k in ARMS]
        hi = [d["arms"][k]["hi"] * d["n"] - d["arms"][k]["wrong"] for k in ARMS]
        a_.bar(ARMS, w, yerr=[lo, hi], color="#55a868", capsize=3)
        a_.set_title(title + " (Wilson 95%)", fontsize=9)
        for i, v in enumerate(w):
            a_.text(i, v, str(v), ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=130)


def main() -> None:
    s: Dict[str, Any] = dict(e0=e0_summary(), pilot=pilot_summary(), timing=timing_summary())
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(s, indent=1, default=float) + "\n")
    DT.write_tables(OUT / "tables.md", tables(s))
    figure(s, OUT / "item14_time.png")
    print(
        json.dumps({k: s[k]["arms"] if k != "timing" else s[k]["arm"] for k in s}, default=float)[
            :1500
        ]
    )
    for name, p in (("A1", CACHE / "runs" / "A1"), ("B1", CACHE / "runs" / "B1")):
        dst = OUT / "runs" / name
        dst.mkdir(parents=True, exist_ok=True)
        for f in ("recon.mic", "records.json", "run.json"):
            (dst / f).write_bytes((p / f).read_bytes())


if __name__ == "__main__":
    main()
