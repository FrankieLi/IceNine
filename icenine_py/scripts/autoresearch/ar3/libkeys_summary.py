"""AR3 step 10: the summary of the validation runs: tables, plot and slim result files.

  libkeys_summary.py

Reads ``keys/check_<kind>.json``, ``keys/score.json``, ``keys/pass_check_S8L.json``,
``keys/off_compare.json`` (all made by the ``libkeys_*.py`` scripts) and the slim copy of each run
(``keys/slim/<kind>_run.json``, made here from ``keys/cache/<kind>/run.json``). The stored runs'
stage numbers and wrong counts are read from ``start/tables.md`` (committed; step 9's tables).
Writes ``keys/tables.md``, ``keys/summary.json`` and ``keys/ar3_keys_evals.png``.
"""

import json
from pathlib import Path
from typing import Any, Dict, List

import libkeys_common as K
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from doc_tables import markdown_table as _markdown_table  # noqa: E402
from doc_tables import write_tables  # noqa: E402


def markdown_table(rows: List[List[Any]], columns: List[str]) -> str:
    """``doc_tables.markdown_table`` for rows given as lists (cells already formatted)."""
    return _markdown_table([dict(zip(columns, r)) for r in rows], columns)


START = K.C.ICE / "benchmarks" / "autoresearch" / "ar3" / "start" / "tables.md"
KINDS = ("D4L", "S5L", "S8L")
STORED_NAME = {"D4L": "D4", "S5L": "S5", "S8L": "S8"}
STAGES = ("seed_step", "grid_fallback", "neighbour", "revisit", "pass_")


def jl(p: Path) -> Any:
    return json.loads(p.read_text())


def table_rows(name: str) -> List[List[str]]:
    """Rows of a committed doc table (cells as strings)."""
    text = START.read_text()
    a = text.index(f"<!-- table:{name} -->")
    b = text.index(f"<!-- /table:{name} -->")
    lines = [ln for ln in text[a:b].splitlines() if ln.startswith("|")][2:]
    return [[c.strip() for c in ln.strip("|").split("|")] for ln in lines]


def stored_stages(sample: str, run: str) -> Dict[str, int]:
    for r in table_rows("p3_stages"):
        if r[0] == sample and r[1] == f"{run} + pass":
            names = ["sample", "run", "seeds", "reseeds", "grid", *STAGES, "total"]
            return {n: int(v) for n, v in zip(names, r) if n not in ("sample", "run")}
    raise KeyError((sample, run))


def lib_stages(st: Dict[str, Any]) -> Dict[str, int]:
    seedlike = int(st["n_evals_seed"] + st["n_evals_reseed"])
    d = dict(
        seeds=int(st["n_seeds"]),
        reseeds=int(st["n_reseeds"]),
        grid=int(st.get("index_fallback_seed_test", 0) + st.get("index_fallback_no_candidate", 0)),
        seed_step=int(st["index_evals"]),
        grid_fallback=seedlike - int(st["index_evals"]),
        neighbour=int(st["n_evals_neighbor"]),
        revisit=int(st["n_evals_revisit"]),
        pass_=int(st["n_evals_boundary_pass"]),
    )
    d["total"] = sum(d[s] for s in STAGES)
    return d


def main() -> None:
    slim = K.OUT / "slim"
    slim.mkdir(parents=True, exist_ok=True)
    runs: Dict[str, Any] = {}
    for kind in KINDS:
        info = jl(K.CACHE / kind / "run.json")
        (slim / f"{kind}_run.json").write_text(json.dumps(info, indent=1, default=float) + "\n")
        runs[kind] = info
    score = jl(K.OUT / "score.json")
    checks = {k: jl(K.OUT / f"check_{k}.json") for k in KINDS}
    S: Dict[str, Any] = dict(runs={}, checks=checks, score=score)
    T: Dict[str, str] = {}

    # ---- stages: script run against library run
    rows = []
    stage_bars: Dict[str, Dict[str, Dict[str, int]]] = {}
    for kind in KINDS:
        name = STORED_NAME[kind]
        lib = lib_stages(runs[kind]["stats"])
        sto = stored_stages("big", name)
        stage_bars[name] = dict(script=sto, library=lib)
        S["runs"][kind] = dict(library=lib, stored=sto)
        wrong = {
            ("script", "before"): score[kind]["stored_before"],
            ("script", "after"): score[kind]["stored_after"],
            ("library", "before"): score[kind]["lib_before"],
            ("library", "after"): score[kind]["lib_after"],
        }
        for src, d in (("script", sto), ("library", lib)):
            rows.append([name, src, d["seeds"], d["reseeds"], d["grid"], d["seed_step"],
                         d["grid_fallback"], d["neighbour"], d["revisit"], d["pass_"], d["total"],
                         wrong[(src, "before")], wrong[(src, "after")]])  # fmt: skip
    T["keys_evals"] = markdown_table(
        rows,
        ["run", "code", "seeds", "reseeds", "grid", "seed_step", "grid_fallback", "neighbour",
         "revisit", "pass", "total", "wrong_before", "wrong_after"],
    )  # fmt: skip

    # ---- exact checks
    rows = []
    for kind in KINDS:
        c = checks[kind]
        pn = c["plan_numbers"]
        for k, (a, b) in pn.items():
            rows.append([kind, k, f"{a:,}", f"{b:,}", "exact" if a == b else "differs"])
        cnt = c["counters"]
        n_diff = len(cnt["differ"])
        rows.append(
            [
                kind,
                "counters (all numeric, no wall times)",
                f"{n_diff} differ",
                f"{cnt['n_stored_counters']} compared",
                "exact" if n_diff == 0 and not cnt["only_in_stored"] else "differs",
            ]
        )

        rec = c["records"]
        rows.append([kind, "per-voxel records", "equal" if rec["equal"] else "differ",
                     "24,570 voxels", "exact" if rec["equal"] else "differs"])  # fmt: skip
        mp = c["map_before_pass_vs_stored_recon_mic"]
        rows.append([kind, "map before the pass vs stored recon.mic", "equal" if mp["equal"]
                     else "differ", f"{mp.get('bytes', 0):,} bytes",
                     "exact" if mp["equal"] else "differs"])  # fmt: skip
        for k, v in c.get("wrapper_counts", {}).items():
            rows.append([kind, f"wrapper: {k}", f"{v[0]:,}", f"{v[1]:,}",
                         "exact" if v[0] == v[1] else "differs"])  # fmt: skip
        for k, v in c.get("early_stops", {}).items():
            rows.append([kind, f"early stop: {k}", f"{v[0]:,}", f"{v[1]:,}",
                         "exact" if v[0] == v[1] else "differs"])  # fmt: skip
    T["keys_exact"] = markdown_table(rows, ["run", "number", "library", "stored", "verdict"])

    # ---- the pass
    rows = []
    for kind in KINDS:
        pf = K.OUT / f"pass_check_{kind}.json"
        if not pf.exists():
            continue
        pc = jl(pf)
        differ = max(
            pc[f"differ_{f}"] for f in ("voxel", "cand", "evals", "cost", "hit", "conf", "end")
        )
        rows.append(
            [
                STORED_NAME[kind],
                f"{pc['n_tries'][0]:,}",
                f"{pc['n_tries'][1]:,}",
                f"{pc['evals'][0]:,}",
                f"{pc['evals'][1]:,}",
                differ,
                pc["voxels_changed"][0],
                pc["changed_voxels_with_a_different_end_after_the_script_rounding"],
            ]
        )
    T["keys_pass1"] = markdown_table(
        rows,
        ["stored map", "tries", "tries_stored", "evals", "evals_stored", "tries_differing",
         "voxels_changed", "other_end"],
    )  # fmt: skip
    rows = []
    for kind in KINDS:
        pr = checks[kind]["pass_in_run"]
        sc = score[kind]
        pa = sc["paired_after"]
        rows.append([kind, f"{pr['tries'][0]:,}", f"{pr['tries'][1]:,}", f"{pr['evals'][0]:,}",
                     f"{pr['evals'][1]:,}", pr.get("same_voxel_and_candidate", "n/a"),
                     pr.get("same_end_exact", "n/a"), sc["lib_after"], sc["stored_after"],
                     pa["fixed_vs_stored"], pa["broken_vs_stored"]])  # fmt: skip
    T["keys_pass2"] = markdown_table(
        rows,
        ["run", "tries", "tries_stored", "evals", "evals_stored", "same_voxel_and_cand",
         "same_end", "wrong_after", "wrong_after_stored", "fixed", "broken"],
    )  # fmt: skip

    # ---- check (c)
    oc = jl(K.OUT / "off_compare.json")
    rows = [
        [f, "equal" if oc[f]["equal"] else "differ", f"{oc[f]['bytes']:,}"]
        for f in ("recon.mic", "orient.npy", "records.json", "stats.json", "logs.json")
    ]
    T["keys_off"] = markdown_table(rows, ["file", "parent vs branch", "bytes"])

    # ---- the measured numbers of steps 7 to 9 (read from start/tables.md)
    wrong = {(r[0], r[1]): r for r in table_rows("p3_wrong")}
    rows = []
    for sample in ("big", "grad", "mid"):
        for name in ("D4", "S5", "S8"):
            s = stored_stages(sample, name)
            w = wrong[(sample, name)]
            rows.append([sample, name, s["total"], w[2], w[4]])
    T["keys_measured"] = markdown_table(
        rows, ["sample", "run", "evals_with_pass", "wrong_before", "wrong_after"]
    )

    write_tables(K.OUT / "tables.md", T)
    S["tables"] = list(T)
    (K.OUT / "summary.json").write_text(json.dumps(S, indent=1, default=float) + "\n")

    # ---- plot: evaluations by stage, script run against library run
    fig, ax = plt.subplots(figsize=(8, 4.5))
    colours = dict(seed_step="#4c72b0", grid_fallback="#c44e52", neighbour="#55a868",
                   revisit="#8172b2", pass_="#ccb974")  # fmt: skip
    labels = dict(seed_step="seed step", grid_fallback="grid fallback", neighbour="neighbour fits",
                  revisit="revisit fits", pass_="boundary pass")  # fmt: skip
    xs, xt = [], []
    x = 0.0
    for name in ("D4", "S5", "S8"):
        for src in ("script", "library"):
            d = stage_bars[name][src]
            bottom = 0.0
            for st in STAGES:
                ax.bar(x, d[st] / 1e6, bottom=bottom, color=colours[st], width=0.8,
                       label=labels[st] if (name, src) == ("D4", "script") else None)  # fmt: skip
                bottom += d[st] / 1e6
            ax.text(x, bottom + 0.15, f"{bottom:.2f}", ha="center", fontsize=8)
            xs.append(x)
            xt.append(f"{name}\n{src}")
            x += 1.0
        x += 0.6
    ax.set_xticks(xs)
    ax.set_xticklabels(xt, fontsize=8)
    ax.set_ylabel("strict-cost evaluations (millions)")
    ax.set_title("Sample big: the stored script runs and the library-key runs, by stage")
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(K.OUT / "ar3_keys_evals.png", dpi=130)
    print("wrote", K.OUT / "tables.md")


if __name__ == "__main__":
    main()
