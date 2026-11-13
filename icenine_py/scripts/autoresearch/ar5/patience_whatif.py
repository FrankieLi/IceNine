"""AR5 step 11, phase A step 1: count checks and the what-if tables, from the stored arm-G files.

  uv run python scripts/autoresearch/ar5/patience_whatif.py [--sample big|grad|mid ...]

Reads (read only, nothing rerun) for ``big``, ``grad``, ``mid`` (never ``held``): arm ``G`` of
step 9 as linked in ``ar5/patience/cache`` by ``patience_link_cache.sh``: ``run.json``,
``e1_log.json.gz`` (the fits of each seed step, interleaved order), ``mc_traces.json.gz`` (the
best cost at each improvement of every MC refinement), ``passes_<s>/D4_G/pass.json``. It

1. checks the stored counts of the plan's stop rule (totals with the pass, seed-step evaluations,
   traces) and exits 2 on any difference (``counts.json`` holds the numbers either way);
2. writes the role table (where the MC evaluations go) and the rule table (``P(k)`` for several k,
   and the stall-while-worse rule (iii)) as ``whatif.json`` and ``whatif.md`` (blocks for
   ``sync_doc_tables.py``). Each rule is walked through whole seed steps with ``patience_rules.
   walk_step`` (arithmetic on the traces; kept ends, whether the step still stops, which end is
   taken; an MC end cut before a later gain keeps its stored 0.95 class unless ``--pessimistic``).

The best-of-120 cost-gated rules and the in-sample checkpoint rule of the plan's table (ii) are not
rebuilt (their definitions are not given to the digit); ``whatif.md`` says so.
"""

import argparse
import gzip
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
PAT = ICE / "benchmarks" / "autoresearch" / "ar5" / "patience"
CACHE = Path(os.environ.get("PATIENCE_CACHE", str(PAT / "cache")))
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent.parent / "common"))

import doc_tables as DT  # noqa: E402
import patience_rules as PR  # noqa: E402

SAMPLES = ("mid", "grad", "big")
# stored counts of the plan's stop rule: totals (with the pass), seed steps, traces, steps
EXPECT: Dict[str, Dict[str, int]] = {
    "mid": dict(total=13004457, seed_evals=5417024, traces=1605, steps=1607),
    "grad": dict(total=10196692, seed_evals=1528284, traces=991, steps=1280),
    "big": dict(total=3807747, seed_evals=687902, traces=153, steps=317),
}
KS = (1000, 2000, 3000, 5000)


def _gz(p: Path) -> Any:
    with gzip.open(p, "rt") as f:
        return json.load(f)


def load(s: str) -> Dict[str, Any]:
    d = CACHE / f"G_{s}"
    run = json.loads((d / "run.json").read_text())
    pas = json.loads((CACHE / f"passes_{s}" / "D4_G" / "pass.json").read_text())
    e1 = _gz(d / "e1_log.json.gz")
    traces = _gz(d / "mc_traces.json.gz")
    return dict(run=run, pass_evals=int(pas["evals"]), e1=e1, traces=traces)


def counts(s: str, d: Dict[str, Any]) -> Dict[str, Any]:
    run = d["run"]
    got: Dict[str, Any] = dict(
        total=int(run["evals_total"]) + d["pass_evals"],
        seed_evals=int(run["stats"]["index_evals"]),
        traces=len(d["traces"]),
        steps=len(d["e1"]),
        e1_evals_sum=int(sum(int(x.get("evals", 0)) for x in d["e1"])),
        before_pass=int(run["evals_total"]),
        pass_evals=d["pass_evals"],
    )
    exp = EXPECT[s]
    got["expected"] = exp
    got["ok"] = bool(
        all(got[k] == exp[k] for k in exp) and got["e1_evals_sum"] == exp["seed_evals"]
    )
    return got


def steps_with_mc(d: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Per seed step with an MC fit: the stored fits (dicts) and the traces per position."""
    by_step: Dict[int, Dict[int, Dict[str, Any]]] = {}
    for t in d["traces"]:
        by_step.setdefault(int(t["step"]), {})[int(t["j"])] = t
    out = []
    for si, e in enumerate(d["e1"], start=1):
        if not e["fits"]:
            continue
        fits = [
            dict(pos=f[0], j=f[1], kind=f[2], hit=f[3], cost=f[4], evals=f[5]) for f in e["fits"]
        ]
        if not any(f["kind"] == "M" for f in fits):
            continue
        tr: Dict[int, Dict[str, Any]] = {}
        for f in fits:
            if f["kind"] == "M":
                t = by_step[si][f["j"]]
                assert int(t["voxel"]) == int(e["voxel"]) and int(t["evals"]) == f["evals"]
                tr[f["pos"]] = t
        out.append(dict(step=si, voxel=int(e["voxel"]), entry=e, fits=fits, traces=tr))
    return out


def rule_walk(
    st: Dict[str, Any], k: int, gated: bool, pessimistic: bool
) -> Tuple[Dict[str, Any], Dict[int, Dict[str, Any]]]:
    fits = st["fits"]
    cuts: Dict[int, Dict[str, Any]] = {}
    kept_cost: List[float] = []
    for p in range(1 + max(f["pos"] for f in fits)):
        c = [f for f in fits if f["pos"] == p and f["kind"] == "C"][0]
        m = [f for f in fits if f["pos"] == p and f["kind"] == "M"]
        cc = float(c["cost"])
        if m:
            ref = min([cc] + kept_cost) if gated else None
            cp = PR.cut_point(st["traces"][p]["trace"], int(m[0]["evals"]), k, worse_than=ref)
            cuts[p] = cp
            mc = float(cp["cost"]) if cp["kind"] == "moved" else float(m[0]["cost"])
            kept_cost.append(mc if mc < cc else cc)
        else:
            kept_cost.append(cc)
    n_avail = max(5, 1 + max(f["pos"] for f in fits))
    return PR.walk_step(fits, cuts, pessimistic=pessimistic, n_avail=n_avail), cuts


def evaluate(
    sts: List[Dict[str, Any]], total: int, k: int, gated: bool, pessimistic: bool
) -> Dict[str, Any]:
    saved = changed = stop_changed = lost = cut_after = cut_moved = extra = 0
    for st in sts:
        w, cuts = rule_walk(st, k, gated, pessimistic)
        stored = PR.walk_step(st["fits"], {})
        saved += stored["evals"] - w["evals"]
        extra += w["extra_evals"]
        for cp in cuts.values():
            cut_after += cp["kind"] == "after_last"
            cut_moved += cp["kind"] == "moved"
        t = w["taken"]
        end = w["kept"][t]
        ch = (
            t != stored["taken"]
            or end["which"] != stored["kept"][stored["taken"]]["which"]
            or end["moved"]
        )
        changed += ch
        if ch and stored["stopped"]:
            stop_changed += 1
        lost += w["lost"]
    return dict(
        saved_evals=int(saved),
        saved_share=saved / total,
        extra_evals=int(extra),
        saved_share_without_extra=(saved + extra) / total,
        taken_changed=int(changed),
        stop_type_changed=int(stop_changed),
        stops_lost=int(lost),
        refinements_cut_after_last=int(cut_after),
        refinements_cut_before_a_gain=int(cut_moved),
    )


def roles(sts: List[Dict[str, Any]], total: int) -> Dict[str, Any]:
    names = ("reached_095_taken", "taken_no_095", "not_taken")
    acc: Dict[str, Dict[str, Any]] = {n: dict(n=0, evals=0, tail=0) for n in names}
    for st in sts:
        w = PR.walk_step(st["fits"], {})
        stopped = w["stopped"]
        for p, t in st["traces"].items():
            ev = int(t["evals"])
            tail = ev - int(t["trace"][-1][0])
            if stopped and p == w["taken"]:
                r = names[0]
            elif (not stopped) and p == w["taken"] and w["kept"][p]["which"] == "M":
                r = names[1]
            else:
                r = names[2]
            acc[r]["n"] += 1
            acc[r]["evals"] += ev
            acc[r]["tail"] += tail
    for r in names:
        acc[r]["share"] = acc[r]["evals"] / total
        acc[r]["tail_share"] = acc[r]["tail"] / total
    return acc


def write_md(out: Dict[str, Any], samples: List[str]) -> None:
    def pct(x: float) -> str:
        return f"{100 * x:.1f}%"

    role_rows = []
    for r, label in (
        ("reached_095_taken", "reached 0.95, taken"),
        ("taken_no_095", "taken in a step with no 0.95 end"),
        ("not_taken", "not taken (late gains idle)"),
    ):
        row: Dict[str, Any] = {"role": label}
        for s in samples:
            x = out[s]["roles"][r]
            row[s] = f"{x['n']:,}; {pct(x['share'])} ({pct(x['tail_share'])})"
        role_rows.append(row)
    rule_rows = []
    for name in out[samples[0]]["rules"]:
        row2: Dict[str, Any] = {"rule": name}
        for s in samples:
            x = out[s]["rules"][name]
            row2[s] = f"{pct(x['saved_share'])}; {x['taken_changed']}"
        rule_rows.append(row2)
    tabs = {
        "patience_roles": DT.markdown_table(role_rows, ["role"] + list(samples)),
        "patience_rules": DT.markdown_table(rule_rows, ["rule"] + list(samples)),
    }
    DT.write_tables(PAT / "whatif.md", tabs)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", nargs="*", default=list(SAMPLES))
    ap.add_argument("--pessimistic", action="store_true")
    a = ap.parse_args()
    PAT.mkdir(parents=True, exist_ok=True)
    out: Dict[str, Any] = dict(pessimistic=a.pessimistic)
    cnt: Dict[str, Any] = {}
    bad = []
    for s in a.sample:
        assert s in SAMPLES, "held is not read by this script"
        d = load(s)
        c = counts(s, d)
        cnt[s] = c
        print(s, "counts", {k: v for k, v in c.items() if k != "expected"}, flush=True)
        if not c["ok"]:
            bad.append(s)
        sts = steps_with_mc(d)
        total = c["total"]
        res: Dict[str, Any] = dict(
            steps_with_mc=len(sts), mc_refinements=sum(len(x["traces"]) for x in sts)
        )
        res["roles"] = roles(sts, total)
        res["rules"] = {}
        for k in KS:
            res["rules"][f"P({k})"] = evaluate(sts, total, k, False, a.pessimistic)
        res["rules"]["(iii) 500 while worse"] = evaluate(sts, total, 500, True, a.pessimistic)
        out[s] = res
    (PAT / "counts.json").write_text(json.dumps(cnt, indent=1) + "\n")
    out["counts_ok"] = not bad
    (PAT / "whatif.json").write_text(json.dumps(out, indent=1) + "\n")
    if bad:
        print("STOP: a stored count does not reproduce:", bad, file=sys.stderr)
        sys.exit(2)
    write_md(out, list(a.sample))
    print("whatif written")


if __name__ == "__main__":
    main()
