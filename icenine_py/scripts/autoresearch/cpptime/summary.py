"""Summary of the C++ pilot time run (cpptime, phase B).

  uv run python scripts/autoresearch/cpptime/summary.py --item14-commit 07d308f

Reads ``benchmarks/autoresearch/cpptime/cache`` (runs C0, C1, C2: parse.json, score.json,
time.txt; quiet/*.json and the .jsonl load samples; pair_C.json), ``trace_C1.json`` (made by
``trace_cpp_log.py``), ``python_time/A1.time, A2.time`` (copies of item 14's ``/usr/bin/time -p``
files) and the Python numbers of item 14 arm A (runs A1, A2) from
``benchmarks/autoresearch/item14/summary.json`` and ``runs/A1``. Writes, beside the cache,
``summary.json``, ``tables.md`` and ``cpp_time.png``. C++ evaluations: not measured (lldb is
unavailable). Every number the findings page quotes is in ``summary.json``.
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
ITEM14 = ICE / "benchmarks" / "autoresearch" / "item14"
SCRIPTS = Path(__file__).resolve().parent
ITEM14_REL = "benchmarks/autoresearch/item14"


def limit_used() -> float:
    return float((SCRIPTS / "r4_limit.txt").read_text().split()[0])


def load_samples(run: str) -> Dict[str, Any]:
    rows = [json.loads(x) for x in (CACHE / "quiet" / f"{run}.jsonl").read_text().splitlines()]
    loads = [r for r in rows if r["kind"] == "load"]
    procs = [r for r in rows if r["kind"] == "procs"]
    lim = limit_used()
    return dict(
        n_load=len(loads),
        n_over_2_5=sum(r["load"] > 2.5 for r in loads),
        n_over_limit=sum(r["load"] > lim for r in loads),
        max_load=max(r["load"] for r in loads),
        n_proc_samples=len(procs),
        load_period_s=round(loads[1]["t"] - loads[0]["t"]),
        proc_period_s=round(procs[1]["t"] - procs[0]["t"]),
    )


def read_time(path: Path) -> Dict[str, float]:
    t = {k: float(v) for k, v in (ln.split() for ln in path.read_text().splitlines())}
    t["user_sys_over_real"] = (t["user"] + t["sys"]) / t["real"]
    t["user_over_real"] = t["user"] / t["real"]
    return t


def cpp_run(name: str) -> Dict[str, Any]:
    d = CACHE / "runs" / name
    p = json.loads((d / "parse.json").read_text())
    s = json.loads((d / "score.json").read_text())
    out: Dict[str, Any] = dict(
        real_s=p["time_real"],
        user_s=p["time_user"],
        sys_s=p["time_sys"],
        user_over_real=p["time_user"] / p["time_real"],
        user_sys_over_real=(p["time_user"] + p["time_sys"]) / p["time_real"],
        bfs_s=p["total_processing_s"],
        start_up_s=p["time_real"] - p["total_processing_s"],
        read_s=p["read_s"],
        global_s=p["global_s_sum"],
        local_s=p["total_processing_s"] - p["global_s_sum"],
        global_searches=p["global_searches"],
        local_optimizations=p["local_optimizations"],
        wrong=s["wrong"],
        wrong_ci=s["wrong_wilson"],
        n_output=s["n_output"],
        n_missing=s["n_missing"],
        samples=load_samples(name),
    )
    if name != "C0":
        q = json.loads((CACHE / "quiet" / f"{name}.json").read_text())
        out.update(
            strict_quiet=q["quiet"],
            relaxed_ok=q["relaxed_ok"],
            median_load=q["median_load"],
            busy=[b["command"] for b in q["busy"]],
        )
    return out


def timing_label(r: Dict[str, Any]) -> str:
    """Timing label of a timed C++ run: worker count, quiet verdict (D-20261114-N1, N2 wording)."""
    s, lim = r["samples"], limit_used()
    ranks = "2 MPI ranks (1 fitting client + 1 server using about one core)"
    if r["strict_quiet"]:
        return (
            f"{ranks}; strictly quiet at limit {lim:g} only: {s['n_over_2_5']} of {s['n_load']} "
            f"load samples are over 2.5 (max {s['max_load']:.2f}); the limit came from an assumed "
            "idle load of 1.5, not a measured one"
        )
    busy = (
        "no process over 20% CPU caught by the sampler"
        if not r["busy"]
        else ("busy: " + ", ".join(r["busy"]))
    )
    note = "; condition (ii) of N1 held only because nothing was caught" if not r["busy"] else ""
    return (
        f"{ranks}; strict quiet check failed: load samples above 2.5 ({s['n_over_2_5']} of "
        f"{s['n_load']}; above {lim:g}: {s['n_over_limit']}) and {busy} (processes sampled every "
        f"{s['proc_period_s']} s, load every {s['load_period_s']} s){note}; reported under R4b "
        "(D-20261114-N1)"
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--item14-commit", required=True, help="commit of the item 14 inputs")
    a = ap.parse_args()
    py = json.loads((ITEM14 / "summary.json").read_text())
    pt = py["timing"]["runs"]
    pil = py["pilot"]["arms"]["F"]
    runs = {n: cpp_run(n) for n in ("C0", "C1", "C2")}
    pair = json.loads((CACHE / "pair_C.json").read_text())
    trace = json.loads((OUT / "trace_C1.json").read_text())
    from score_cpp import score_run

    paired = score_run(CACHE / "runs" / "C1", python_run=ITEM14 / "runs" / "A1")
    pyr: Dict[str, Any] = {}
    for n in ("A1", "A2"):
        r = pt[n]
        tm = read_time(OUT / "python_time" / f"{n}.time")
        pyr[n] = dict(
            wall_s=r["wall_total_s"],
            real_s=r["real_s"],
            full_s=r["full_s"],
            local_s=r["local_s"],
            label=r["label"],
            median_load=r["median_load"],
            n_load=r["n_load"],
            user_over_real=tm["user_over_real"],
            user_sys_over_real=tm["user_sys_over_real"],
        )
    setup_s = json.loads((ITEM14 / "runs" / "A1" / "run.json").read_text())["setup_s"]
    py_lo, py_hi = wilson(pil["wrong"], 2000)
    c1, c2 = runs["C1"], runs["C2"]
    a1, a2 = pyr["A1"], pyr["A2"]
    ratios = dict(
        whole_process_A1_over_C1=a1["real_s"] / c1["real_s"],
        bfs_A1_over_C1=a1["wall_s"] / c1["bfs_s"],
        per_global_search_A1_over_C1=(a1["full_s"] / pil["full_searches"])
        / (c1["global_s"] / c1["global_searches"]),
        caveat="different jobs: C++ left 416 of 2000 wrong, Python arm A 24; labels in the page",
    )
    gaps = dict(
        cpp_pair_C1_C2=abs(c1["real_s"] - c2["real_s"]) / min(c1["real_s"], c2["real_s"]),
        python_pair_A1_A2=abs(a1["real_s"] - a2["real_s"]) / min(a1["real_s"], a2["real_s"]),
    )
    summ = dict(
        limit_used=limit_used(),
        cpp=runs,
        c0_max_load=runs["C0"]["samples"]["max_load"],
        pair=pair,
        pair_gaps=gaps,
        ratios=ratios,
        python=pyr,
        python_arm_a=dict(
            wrong=pil["wrong"],
            wrong_ci=[py_lo, py_hi],
            full_searches=pil["full_searches"],
            evals_total=pil["evals_total"],
            evals_full=pil["evals_full"],
            setup_s=setup_s,
        ),
        trace_c1=trace,
        paired_c1_vs_python_a1=paired,
        evaluations_cpp="not measured (lldb cannot launch the binary: developer mode off)",
        provenance=dict(
            item14_summary=f"{ITEM14_REL}/summary.json",
            item14_python_recon=f"{ITEM14_REL}/runs/A1/recon.mic",
            item14_time_files="python_time/A1.time, A2.time (copies of item 14 cache/logs)",
            item14_commit=a.item14_commit,
        ),
    )
    (OUT / "summary.json").write_text(json.dumps(summ, indent=1) + "\n")

    # --- tables
    t: Dict[str, str] = {}
    rows: List[Dict[str, Any]] = []
    for n in ("A1", "A2"):
        rows.append(
            dict(
                run=n,
                side="Python",
                real=pyr[n]["real_s"],
                label="single-worker; " + pyr[n]["label"],
            )
        )
    for n in ("C1", "C2"):
        rows.append(
            dict(run=n, side="C++ one client", real=runs[n]["real_s"], label=timing_label(runs[n]))
        )
    t["walls"] = DT.markdown_table(rows, ["run", "side", "real", "label"], {"real": ".2f"}).replace(
        "| real |", "| whole process (s) |"
    )
    t["cpu"] = DT.markdown_table(
        [
            dict(run="Python A1", u=a1["user_over_real"], us=a1["user_sys_over_real"]),
            dict(run="Python A2", u=a2["user_over_real"], us=a2["user_sys_over_real"]),
            dict(run="C++ C1", u=c1["user_over_real"], us=c1["user_sys_over_real"]),
            dict(run="C++ C2", u=c2["user_over_real"], us=c2["user_sys_over_real"]),
        ],
        ["run", "u", "us"],
        {"u": ".2f", "us": ".2f"},
    ).replace("| u | us |", "| user / real | (user + sys) / real |")
    rows = []
    for n in ("C0", "C1", "C2"):
        r = runs[n]
        rows.append(
            dict(
                run=n + (" (not gated: no preflight or quiet check)" if n == "C0" else ""),
                real=r["real_s"],
                su=r["start_up_s"],
                g=r["global_s"],
                loc=r["local_s"],
            )
        )
    t["cpp_split"] = DT.markdown_table(
        rows,
        ["run", "real", "su", "g", "loc"],
        {"real": ".2f", "su": ".1f", "g": ".0f", "loc": ".0f"},
    ).replace(
        "| real | su | g | loc |",
        "| whole (s) | start-up and exit (s) | global searches (s) | local optimizations (s) |",
    )
    wlo, whi = c1["wrong_ci"]
    t["work"] = DT.markdown_table(
        [
            dict(
                side="Python arm A (item 14)",
                g=pil["full_searches"],
                loc="not counted here",
                ev=pil["evals_total"],
                wrong=f"{pil['wrong']} of 2000 ({py_lo:.1%} to {py_hi:.1%})",
            ),
            dict(
                side="C++ one client (C1; C0 and C2 equal)",
                g=c1["global_searches"],
                loc=c1["local_optimizations"],
                ev="not measured",
                wrong=f"{c1['wrong']} of 2000 ({wlo:.1%} to {whi:.1%})",
            ),
        ],
        ["side", "g", "loc", "ev", "wrong"],
    ).replace(
        "| side | g | loc | ev | wrong |",
        "| side | global searches | local optimizations | evaluations | wrong voxels |",
    )
    tr = trace
    t["trace"] = DT.markdown_table(
        [
            dict(what="local-optimization records in the client log", n=tr["records"]),
            dict(what="voxels marked REFIT (hit ratio / best at most 0.9)", n=tr["refit"]),
            dict(what="of those, wrong", n=tr["refit_wrong"]),
            dict(what="voxels accepted (above 0.9) that are wrong", n=tr["accepted_wrong"]),
            dict(
                what="voxels with no record (seeds), of which wrong", n=tr["voxels_not_in_records"]
            ),
            dict(
                what="wrong among the voxels with no record",
                n=tr["wrong_among_voxels_not_in_records"],
            ),
            dict(what="wrong voxels in total", n=tr["wrong_total"]),
            dict(
                what="wrong voxels within 1 degree of another grain's truth",
                n=tr["wrong_near_another_grain_truth"],
            ),
            dict(
                what="records whose Optimized orientation equals Original (printed digits)",
                n=tr["unchanged"],
            ),
            dict(
                what="distinct orientations in the C++ output map",
                n=tr["distinct_output_orientations"],
            ),
        ],
        ["what", "n"],
    ).replace("| what | n |", "| count in run C1 | voxels |")
    pw = paired["paired_on_cpp_voxels"]
    t["paired"] = DT.markdown_table(
        [dict(b=pw["both"], oc=pw["only_cpp"], op=pw["only_python"], n=pw["neither"])],
        ["b", "oc", "op", "n"],
    ).replace(
        "| b | oc | op | n |",
        "| wrong in both | wrong only in C++ | wrong only in Python | wrong in neither |",
    )
    t["unit"] = DT.markdown_table(
        [
            dict(
                side="Python arm A1",
                gs=a1["full_s"] / pil["full_searches"],
                bfs=a1["wall_s"],
                ev=1e6 * a1["wall_s"] / pil["evals_total"],
            ),
            dict(
                side="C++ C1", gs=c1["global_s"] / c1["global_searches"], bfs=c1["bfs_s"], ev=None
            ),
        ],
        ["side", "gs", "bfs", "ev"],
        {"gs": ".1f", "bfs": ".0f", "ev": ".1f"},
    ).replace(
        "| side | gs | bfs | ev |",
        "| side | seconds per global search | BFS time (s) | microseconds per evaluation |",
    )
    qrows = []
    for n in ("C0", "C1", "C2"):
        s = runs[n]["samples"]
        qrows.append(
            dict(
                run=n,
                n=s["n_load"],
                o25=s["n_over_2_5"],
                ol=s["n_over_limit"],
                mx=s["max_load"],
                md=runs[n].get("median_load"),
                strict=str(runs[n].get("strict_quiet", "not checked")),
            )
        )
    t["quiet"] = DT.markdown_table(
        qrows,
        ["run", "n", "o25", "ol", "mx", "md", "strict"],
        {"mx": ".2f", "md": ".2f"},
    ).replace(
        "| run | n | o25 | ol | mx | md | strict |",
        f"| run | load samples | over 2.5 | over {limit_used():g} (limit) | max load "
        "| median load | strict check passed |",
    )
    DT.write_tables(OUT / "tables.md", t)

    # --- plot
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10, 4.2), gridspec_kw=dict(width_ratios=[3, 2]))
    names = ["Python A1", "Python A2", "C++ C1", "C++ C2"]
    start = [setup_s, setup_s, c1["start_up_s"], c2["start_up_s"]]
    glob = [a1["full_s"], a2["full_s"], c1["global_s"], c2["global_s"]]
    loc = [a1["local_s"], a2["local_s"], c1["local_s"], c2["local_s"]]
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
    ax.set_title("Wall time on the 2,000-voxel pilot (different jobs)")
    ax.legend()
    bx.bar(["wrong\nPython", "wrong\nC++"], [pil["wrong"], c1["wrong"]], color="#d62728")
    bx2 = bx.twinx()
    bx2.bar(
        ["global\nPython", "global\nC++"],
        [pil["full_searches"], c1["global_searches"]],
        color="#2ca02c",
        alpha=0.7,
    )
    bx.set_title("Wrong voxels (red), global searches (green)")
    fig.text(
        0.01,
        0.01,
        "Python: single worker. C++: 2 MPI ranks, 1 fitting client. Quiet labels in the page.",
        fontsize=7,
    )
    fig.tight_layout()
    fig.savefig(OUT / "cpp_time.png", dpi=130)
    print("wrote summary.json, tables.md, cpp_time.png")


if __name__ == "__main__":
    main()
