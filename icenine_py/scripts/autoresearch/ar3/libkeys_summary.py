"""AR3 step 10: the summary of the validation runs: tables, plot and slim result files.

  libkeys_summary.py [--suite-log LOG] [--file-log LOG]

Reads ``keys/check_<kind>.json``, ``keys/score.json``, ``keys/pass_check_<kind>.json``,
``keys/off_compare.json``, ``keys/fresh_check_<kind>.json`` (fixes round, ``libkeys_fresh.py``;
those present) (all made by the ``libkeys_*.py`` scripts) and the slim copy of each run
(``keys/slim/<kind>_run.json``, made here from ``keys/cache/<kind>/run.json``). The stored runs'
stage numbers and wrong counts (with their Wilson intervals) are read from ``start/tables.md``
(committed; step 9's tables). The suite counts come from pytest logs given with ``--suite-log`` /
``--file-log`` (their last lines are kept in ``keys/suite.json``), else from ``keys/suite.json``.
Writes ``keys/tables.md``, ``keys/summary.json`` and ``keys/ar3_keys_evals.png``.
"""

import argparse
import json
import re
from pathlib import Path
from typing import Any, Dict, List, Optional

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


FRESH_ROWS = (
    ("i_inrun_vs_stored_step10_inrun", "(i) rerun's in-run pass vs step 10's in-run pass"),
    ("iii_inproc_rounded_vs_script", "(iii) same process, rounded map vs script tries"),
    ("iv_fresh_unrounded_vs_inrun", "(iv) fresh process, unrounded map vs in-run pass"),
    ("rounding_inproc_rounded_vs_inrun", "(v) same process, rounded map vs in-run pass"),
)
SUITE = K.OUT / "suite.json"
_LAST = re.compile(r"(\d+) (passed|skipped|deselected|failed|errors?|xfailed|xpassed)")


def parse_last_line(log: Path) -> Dict[str, Any]:
    """Counts of pytest's last summary line in ``log`` (the line itself is kept too)."""
    lines = [ln for ln in log.read_text().splitlines() if " in " in ln and "passed" in ln]
    line = lines[-1].strip().strip("=").strip()
    return dict(line=line, **{k: int(n) for n, k in _LAST.findall(line)})


def suite_counts(full: Optional[Path] = None, file_: Optional[Path] = None) -> Dict[str, Any]:
    """The suite counts of this branch: parsed from pytest logs when given (and stored in
    ``keys/suite.json``), else read from ``keys/suite.json``."""
    old = jl(SUITE) if SUITE.exists() else {}
    if full is not None:
        old["full_suite"] = parse_last_line(full)
    if file_ is not None:
        old["test_bfs_ar3_keys"] = parse_last_line(file_)
    if full is not None or file_ is not None:
        SUITE.write_text(json.dumps(old, indent=1) + "\n")
    return old


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite-log", type=Path, help="log of `uv run pytest tests/`")
    ap.add_argument(
        "--file-log", type=Path, help="log of `uv run pytest tests/test_bfs_ar3_keys.py`"
    )
    a = ap.parse_args()
    if a.suite_log or a.file_log:
        suite_counts(a.suite_log, a.file_log)
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
                     pa["fixed_vs_stored"], pa["broken_vs_stored"],
                     pa["both_wrong"]])  # fmt: skip
    T["keys_pass2"] = markdown_table(
        rows,
        ["run", "tries", "tries_stored", "evals", "evals_stored", "same_voxel_and_cand",
         "same_end", "wrong_after", "wrong_after_stored", "fixed", "broken", "both_wrong"],
    )  # fmt: skip

    # ---- fixes round (D-20261111-N9): the reading of the in-run pass, tested
    fresh = {k: jl(K.OUT / f"fresh_check_{k}.json") for k in KINDS
             if (K.OUT / f"fresh_check_{k}.json").exists()}  # fmt: skip
    S["fresh_checks"] = fresh
    rows = []
    for kind, fc in fresh.items():
        for key, what in FRESH_ROWS:
            if key not in fc:
                continue
            c = fc[key]
            rows.append([kind, what, f"{c['n_tries'][0]:,}", f"{c['n_tries'][1]:,}",
                         f"{c['evals'][0]:,}", f"{c['evals'][1]:,}",
                         c["voxels_with_differing_tries"],
                         "equal" if c["equal"] else "differ"])  # fmt: skip
        rows.append(
            [
                kind,
                "(ii) before-pass .mic vs step 10's (bytes)",
                "",
                "",
                "",
                "",
                "",
                "equal" if fc.get("ii_before_mic_equals_stored_step10_bytes") else "differ",
            ]
        )
        rows.append(
            [
                kind,
                "(ii) map read back vs the script's before_*.npy",
                "",
                "",
                "",
                "",
                "",
                "equal" if fc.get("ii_read_back_equals_script_before_npy") else "differ",
            ]
        )
    T["keys_fresh"] = markdown_table(
        rows,
        ["run", "compared", "tries_a", "tries_b", "evals_a", "evals_b", "voxels_differing",
         "verdict"],
    )  # fmt: skip

    # ---- costs left out of the totals, and the wall times (all contended)
    rows = []
    walls = {}
    for kind in KINDS:
        st = runs[kind]["stats"]
        hit = fresh.get(kind, {}).get("early_stop_hit_evals")
        if hit is None:
            hit = 0 if runs[kind]["stats"].get("index_stopped_early") is None else "not stored"
        walls[kind] = dict(
            wall_total_s=float(runs[kind]["wall_total_s"]),
            wall_boundary_pass_s=float(st["wall_boundary_pass_s"]),
            timing_label=runs[kind]["timing_label"],
        )
        rows.append([kind, f"{int(st['n_evals_boundary_score']):,}",
                     f"{hit:,}" if isinstance(hit, int) else hit,
                     f"{walls[kind]['wall_total_s']:,.0f}",
                     f"{walls[kind]['wall_boundary_pass_s']:.1f}"])  # fmt: skip
    S["walls_contended"] = walls
    T["keys_uncounted"] = markdown_table(
        rows, ["run", "boundary_score_evals", "early_stop_hit_evals", "wall_s", "pass_wall_s"]
    )

    # ---- the test suite on this branch (pytest's last lines, kept in keys/suite.json)
    S["suite"] = suite_counts()
    S["notes"] = dict(
        stored_S8_name=(
            "The stored S8 run.json carries the name 'S5' (reseed_s8.py calls run_seed.run('S5')"
            " with the early stop patched in); its s8.json reads 'S8 = S5 + E1' and its totals"
            " are S8's. check_S8L.json's stored_run 'S5' is that name, not a comparison with S5."
        )
    )

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
            rows.append([sample, name, s["total"], w[2], w[3], w[4], w[5]])
    T["keys_measured"] = markdown_table(
        rows,
        ["sample", "run", "evals_with_pass", "wrong_before", "before_ci_pct", "wrong_after",
         "after_ci_pct"],
    )  # fmt: skip

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
