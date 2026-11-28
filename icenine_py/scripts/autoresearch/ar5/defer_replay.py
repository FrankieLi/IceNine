"""AR5 step 22, phase A step 4 (gate g1): the tracked replay of rule ``DEF`` on arm ``J5``'s stored
seed steps of the seven disclosed samples; prediction 1.

  uv run python scripts/autoresearch/ar5/defer_replay.py [--cache DIR]

Arithmetic only: no run, no image, no truth. Reads (read only; linked into the cache by
``defer_link_cache.sh``) ``J5_<s>/e1_log.json.gz`` and ``mc_traces.json.gz`` (every seed step's fits
and every MC refinement's checks) and ``run.json`` (the search's evaluations), and step 21's tracked
``ar5/screen/tables.md`` (``J5S``'s pass and repair evaluations: ``J5S`` total = search + pass +
repair). For each seed step ``defer_rules.walk_step`` gives its type (A, N, C or none), the exact
saving of A steps and C's earlier tails, and the evaluations a type A or N step costs under ``DEF``
(``J5``'s own count minus the tails). ``defer_rules.size_sample`` gives the plan's table row; each
row is compared with ``defer_rules.PLAN_TABLE`` (counts exactly; shares at the table's decimals).

Writes ``ar5/defer/replay.json`` (one row per sample; the per-step walks are not stored: the gate
recomputes them) and ``ar5/defer/replay_tables.md`` (blocks for ``sync_doc_tables.py``). Exit 1 when
a stored count does not reproduce or prediction 1 fails (the plan's stop rule).
"""

import argparse
import gzip
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICE / "scripts" / "common"))
import defer_guard as DG  # noqa: E402
import defer_rules as DR  # noqa: E402
import doc_tables as DT  # noqa: E402

RD = ICE / "benchmarks" / "autoresearch" / "ar5" / "defer"
S21_TABLES = ICE / "benchmarks" / "autoresearch" / "ar5" / "screen" / "tables.md"


def read_gz(p: Path) -> Any:
    with gzip.open(DG.guarded(p), "rt") as f:
        return json.load(f)


def j5s_rows() -> Dict[str, Dict[str, int]]:
    """Step 21's tables: per sample the first row of its J5 / J5S table (12 columns): the search,
    the old pass, the screen pass, the repair and the total of ``J5S``."""
    rows: Dict[str, Dict[str, int]] = {}
    for line in DG.guarded(S21_TABLES).read_text().splitlines():
        if not line.startswith("| "):
            continue
        c = [x.strip() for x in line.strip("|").split("|")]
        name = c[0].split()[0] if c and c[0] else ""
        if name in DR.DISCLOSED and len(c) == 12 and name not in rows:

            def n(s: str) -> int:
                return int(s.replace(",", ""))

            rows[name] = dict(
                search=n(c[2]), old_pass=n(c[3]), pass_=n(c[7]), repair=n(c[8]), total=n(c[9])
            )
    return rows


def walks(entries: List[Dict[str, Any]], traces: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """``walk_step`` of every seed step (``mc_traces`` rows are in the order of the MC fits)."""
    it = iter(traces)
    out: List[Dict[str, Any]] = []
    for i, x in enumerate(entries):
        by: Dict[int, Dict[str, Any]] = {}
        for f in x["fits"]:
            if f[2] == "M":
                t = next(it)
                assert t["voxel"] == x["voxel"] and t["evals"] == f[5], (x["voxel"], t["evals"], f)
                by[int(f[0])] = t
        w = DR.walk_step(x, by)
        w["index"] = i
        out.append(w)
    assert next(it, None) is None, "unused mc_traces rows"
    return out


def one(s: str, cache: Path, tab: Dict[str, Dict[str, int]]) -> Tuple[Dict[str, Any], List[str]]:
    d = cache / f"J5_{s}"
    entries = read_gz(d / "e1_log.json.gz")
    traces = read_gz(d / "mc_traces.json.gz")
    st = json.loads(DG.guarded(d / "run.json").read_text())["stats"]
    search = int(
        st["n_evals_seed"] + st["n_evals_reseed"] + st["n_evals_neighbor"] + st["n_evals_revisit"]
    )
    bad: List[str] = []
    t = tab[s]
    if search != t["search"]:
        bad.append(f"{s}: stored search evaluations {search} vs step 21's table {t['search']}")
    if t["search"] + t["pass_"] + t["repair"] != t["total"]:
        bad.append(f"{s}: step 21's table does not add up")
    w = walks(entries, traces)
    r = DR.size_sample(w, t["total"])
    r.update(
        sample=s,
        steps=len(entries),
        mc_refinements=len(traces),
        search_evals=search,
        screen_pass_evals=t["pass_"],
        repair_evals=t["repair"],
        n_none=sum(1 for x in w if x["type"] == "-"),
        n_c_k={str(k): sum(1 for x in w if x["type"] == "C" and x["stop"] == k) for k in range(4)},
    )
    bad += DR.matches_plan(s, r)
    r["reproduces_plan"] = not bad
    return r, bad


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", default=str(RD / "cache"))
    a = ap.parse_args()
    cache = Path(a.cache)
    tab = j5s_rows()
    missing = [s for s in DR.DISCLOSED if s not in tab]
    if missing:
        print("step 21's tables lack:", missing)
        return 1
    rows: List[Dict[str, Any]] = []
    bad: List[str] = []
    for s in DR.DISCLOSED:
        r, b = one(s, cache, tab)
        rows.append(r)
        bad += b
        print(
            f"{s:6s} A {r['a_steps']:5d} {r['a_pct']:5.2f}% C {r['c_steps']:4d} "
            f"N {r['n_steps']:4d} "
            f"net [{r['lo']:.1f}, {r['hi']:.1f}]% {'OK' if r['reproduces_plan'] else 'DIFFERS'}"
        )
    RD.mkdir(parents=True, exist_ok=True)
    (RD / "replay.json").write_text(
        json.dumps(dict(rows=rows, reproduces_plan=not bad, differences=bad), indent=1) + "\n"
    )
    cols = ["sample", "total", "a_steps", "a_saving", "a_pct", "c_steps", "n_steps", "lo", "hi"]
    fmt = dict(
        total=",d",
        a_steps=",d",
        a_saving=",d",
        a_pct=".2f",
        c_steps=",d",
        n_steps=",d",
        lo=".1f",
        hi=".1f",
    )
    DT.write_tables(
        RD / "replay_tables.md",
        dict(
            ar5_defer_replay=DT.markdown_table(rows, cols, fmt),
        ),
    )
    if bad:
        print("PREDICTION 1 DOES NOT REPRODUCE (stop rule):")
        for b in bad:
            print("  ", b)
        return 1
    print("prediction 1 reproduces exactly on the seven disclosed samples")
    return 0


if __name__ == "__main__":
    sys.exit(main())
