"""Summary of the C++ pilot time run (cpptime, phase B).

  uv run python scripts/autoresearch/cpptime/summary.py --item14-summary FILE --python-recon FILE

Reads ``benchmarks/autoresearch/cpptime/cache`` (runs C0, C1, C2: parse.json, score.json,
time.txt; quiet/C1.json, C2.json and the .jsonl load samples; pair_C.json) and the Python
numbers of item 14 (arm A: runs A1, A2) from its ``summary.json``. Writes, beside the cache,
``summary.json``, ``tables.md`` and ``cpp_time.png``. Evaluations: not measured (lldb unavailable).
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

ICE = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ICE / "scripts" / "common"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import doc_tables as DT  # noqa: E402
from stats import wilson  # noqa: E402

OUT = ICE / "benchmarks" / "autoresearch" / "cpptime"
CACHE = OUT / "cache"
PY_EVALS = 25_602_383  # item 14 arm A, both runs (summary.json pilot.arms.F)
PY_FULL_SEARCHES = 69


def load_samples(run: str) -> Dict[str, Any]:
    rows = [json.loads(x) for x in (CACHE / "quiet" / f"{run}.jsonl").read_text().splitlines()]
    loads = [r for r in rows if r["kind"] == "load"]
    procs = [r for r in rows if r["kind"] == "procs"]
    return dict(
        n_load=len(loads),
        n_over_2_5=sum(r["load"] > 2.5 for r in loads),
        n_over_3_5=sum(r["load"] > 3.5 for r in loads),
        max_load=max(r["load"] for r in loads),
        n_proc_samples=len(procs),
        load_period_s=round(loads[1]["t"] - loads[0]["t"]),
        proc_period_s=round(procs[1]["t"] - procs[0]["t"]),
    )


def cpp_run(name: str) -> Dict[str, Any]:
    d = CACHE / "runs" / name
    p = json.loads((d / "parse.json").read_text())
    s = json.loads((d / "score.json").read_text())
    start_up = p["time_real"] - p["total_processing_s"]
    out: Dict[str, Any] = dict(
        real_s=p["time_real"],
        user_s=p["time_user"],
        sys_s=p["time_sys"],
        user_over_real=p["time_user"] / p["time_real"],
        user_sys_over_real=(p["time_user"] + p["time_sys"]) / p["time_real"],
        bfs_s=p["total_processing_s"],
        start_up_s=start_up,
        read_s=p["read_s"],
        global_s=p["global_s_sum"],
        local_s=p["total_processing_s"] - p["global_s_sum"],
        global_searches=p["global_searches"],
        local_optimizations=p["local_optimizations"],
        rejected_seeds=p["rejected_seeds"],
        wrong=s["wrong"],
        wrong_ci=s["wrong_wilson"],
        n_output=s["n_output"],
        n_missing=s["n_missing"],
    )
    if name != "C0":
        q = json.loads((CACHE / "quiet" / f"{name}.json").read_text())
        out.update(
            strict_quiet=q["quiet"],
            relaxed_ok=q["relaxed_ok"],
            label=q["label"],
            median_load=q["median_load"],
            busy=[b["command"] for b in q["busy"]],
            samples=load_samples(name),
        )
    return out


def timing_label(r: Dict[str, Any]) -> str:
    s = r["samples"]
    ranks = "2 MPI ranks: 1 fitting client + 1 busy-waiting server"
    if r["strict_quiet"]:
        return f"{ranks}; strictly quiet (limit 3.5)"
    return (
        f"{ranks}; strict quiet check failed (load over 3.5 in {s['n_over_3_5']} of {s['n_load']} "
        f"samples, over 2.5 in {s['n_over_2_5']}; no process over 20% CPU caught by the sampler); "
        "reported under R4b"
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--item14-summary", type=Path, required=True)
    ap.add_argument("--python-recon", type=Path, required=True, help="item 14 A1 recon.mic folder")
    a = ap.parse_args()
    py = json.loads(a.item14_summary.read_text())
    pt = py["timing"]["runs"]
    pil = py["pilot"]["arms"]["F"]
    runs = {n: cpp_run(n) for n in ("C0", "C1", "C2")}
    pair = json.loads((CACHE / "pair_C.json").read_text())
    from score_cpp import score_run

    paired = score_run(CACHE / "runs" / "C1", python_run=a.python_recon)
    pyr: Dict[str, Any] = {}
    for n in ("A1", "A2"):
        r = pt[n]
        pyr[n] = dict(
            wall_s=r["wall_total_s"],
            real_s=r["real_s"],
            full_s=r["full_s"],
            local_s=r["local_s"],
            label=r["label"],
            median_load=r["median_load"],
            n_load=r["n_load"],
        )
    summ = dict(
        cpp=runs,
        pair=pair,
        python=pyr,
        python_arm_a=dict(
            wrong=pil["wrong"],
            full_searches=pil["full_searches"],
            evals_total=pil["evals_total"],
            evals_full=pil["evals_full"],
        ),
        paired_c1_vs_python_a1=paired,
        evaluations_cpp="not measured (lldb cannot launch the binary: developer mode off)",
    )
    (OUT / "summary.json").write_text(json.dumps(summ, indent=1) + "\n")

    # --- tables
    t: Dict[str, str] = {}
    rows: List[Dict[str, Any]] = []
    for n, side in (("A1", "Python"), ("A2", "Python")):
        r = pyr[n]
        rows.append(dict(run=n, side=side, real=r["real_s"], label="single-worker; " + r["label"]))
    for n in ("C1", "C2"):
        r = runs[n]
        rows.append(dict(run=n, side="C++ one client", real=r["real_s"], label=timing_label(r)))
    t["walls"] = DT.markdown_table(rows, ["run", "side", "real", "label"], {"real": ".2f"}).replace(
        "| real |", "| whole process (s) |"
    )
    rows = []
    for n in ("C0", "C1", "C2"):
        r = runs[n]
        rows.append(
            dict(
                run=n + (" (untimed)" if n == "C0" else ""),
                real=r["real_s"],
                su=r["start_up_s"],
                g=r["global_s"],
                loc=r["local_s"],
                u=r["user_over_real"],
                us=r["user_sys_over_real"],
            )
        )
    t["cpp_split"] = DT.markdown_table(
        rows,
        ["run", "real", "su", "g", "loc", "u", "us"],
        {"real": ".2f", "su": ".1f", "g": ".0f", "loc": ".0f", "u": ".2f", "us": ".2f"},
    ).replace(
        "| real | su | g | loc | u | us |",
        "| whole (s) | start-up and exit (s) | global searches (s) | local optimizations (s) "
        "| user / real | (user + sys) / real |",
    )
    c1, a1, a2 = runs["C1"], pyr["A1"], pyr["A2"]
    wrong_lo, wrong_hi = c1["wrong_ci"]
    t["work"] = DT.markdown_table(
        [
            dict(
                side="Python arm A (item 14)",
                g=pil["full_searches"],
                loc="not counted here",
                ev=pil["evals_total"],
                wrong=f"{pil['wrong']} of 2000",
            ),
            dict(
                side="C++ one client (C1; C0 and C2 equal)",
                g=c1["global_searches"],
                loc=c1["local_optimizations"],
                ev="not measured",
                wrong=f"{c1['wrong']} of 2000 ({wrong_lo:.1%} to {wrong_hi:.1%})",
            ),
        ],
        ["side", "g", "loc", "ev", "wrong"],
    ).replace(
        "| side | g | loc | ev | wrong |",
        "| side | global searches | local optimizations | evaluations | wrong voxels |",
    )
    pw = paired["paired_on_cpp_voxels"]
    t["paired"] = DT.markdown_table(
        [dict(b=pw["both"], oc=pw["only_cpp"], op=pw["only_python"], n=pw["neither"])],
        ["b", "oc", "op", "n"],
    ).replace(
        "| b | oc | op | n |",
        "| wrong in both | wrong only in C++ | wrong only in Python | wrong in neither |",
    )
    pc = lambda r, k: r[k] / r["global_searches"]  # noqa: E731
    t["unit"] = DT.markdown_table(
        [
            dict(
                side="Python arm A1",
                gs=a1["full_s"] / pil["full_searches"],
                bfs=a1["wall_s"],
                ev=1e6 * a1["wall_s"] / pil["evals_total"],
            ),
            dict(side="C++ C1", gs=pc(c1, "global_s"), bfs=c1["bfs_s"], ev=None),
        ],
        ["side", "gs", "bfs", "ev"],
        {"gs": ".1f", "bfs": ".0f", "ev": ".1f"},
    ).replace(
        "| side | gs | bfs | ev |",
        "| side | seconds per global search | BFS time (s) | microseconds per evaluation |",
    )
    s1, s2 = runs["C1"]["samples"], runs["C2"]["samples"]
    t["quiet"] = DT.markdown_table(
        [
            dict(
                run="C1",
                n=s1["n_load"],
                o25=s1["n_over_2_5"],
                o35=s1["n_over_3_5"],
                mx=s1["max_load"],
                md=runs["C1"]["median_load"],
                strict=str(runs["C1"]["strict_quiet"]),
            ),
            dict(
                run="C2",
                n=s2["n_load"],
                o25=s2["n_over_2_5"],
                o35=s2["n_over_3_5"],
                mx=s2["max_load"],
                md=runs["C2"]["median_load"],
                strict=str(runs["C2"]["strict_quiet"]),
            ),
        ],
        ["run", "n", "o25", "o35", "mx", "md", "strict"],
        {"mx": ".2f", "md": ".2f"},
    ).replace(
        "| run | n | o25 | o35 | mx | md | strict |",
        "| run | load samples | over 2.5 | over 3.5 (limit) | max load | median load "
        "| strict check passed |",
    )
    DT.write_tables(OUT / "tables.md", t)

    # --- plot
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10, 4.2), gridspec_kw=dict(width_ratios=[3, 2]))
    names = ["Python A1", "Python A2", "C++ C1", "C++ C2"]
    start = [0.1, 0.1, runs["C1"]["start_up_s"], runs["C2"]["start_up_s"]]
    glob = [a1["full_s"], a2["full_s"], runs["C1"]["global_s"], runs["C2"]["global_s"]]
    loc = [a1["local_s"], a2["local_s"], runs["C1"]["local_s"], runs["C2"]["local_s"]]
    ax.bar(names, start, label="start-up", color="#999999")
    ax.bar(names, glob, bottom=start, label="global searches", color="#1f77b4")
    ax.bar(
        names,
        loc,
        bottom=[s + g for s, g in zip(start, glob)],
        label="local optimizations",
        color="#ff7f0e",
    )
    ax.set_ylabel("whole process (s)")
    ax.set_title("Wall time on the 2,000-voxel pilot")
    ax.legend()
    bx.bar(["wrong\nPython", "wrong\nC++"], [pil["wrong"], c1["wrong"]], color="#d62728")
    bx2 = bx.twinx()
    bx2.bar(
        ["global\nPython", "global\nC++"],
        [pil["full_searches"], c1["global_searches"]],
        color="#2ca02c",
        alpha=0.7,
    )
    bx.set_title("Wrong voxels (red) and global searches (green)")
    fig.text(
        0.01,
        0.01,
        "Python: single worker. C++: 2 MPI ranks (1 fitting client). Labels in the page.",
        fontsize=7,
    )
    fig.tight_layout()
    fig.savefig(OUT / "cpp_time.png", dpi=130)
    print("wrote summary.json, tables.md, cpp_time.png")


if __name__ == "__main__":
    main()
