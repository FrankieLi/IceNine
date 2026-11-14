"""AR5 step 12, phase A step 1: count checks and the what-if tables, from the stored files.

  uv run python scripts/autoresearch/ar5/stallhit_whatif.py [--sample big|grad|mid ...]
  uv run python scripts/autoresearch/ar5/stallhit_whatif.py --md-only   # whatif.md from the json

Reads (read only, nothing rerun; never ``held``) in ``ar5/stallhit/cache`` (linked by
``stallhit_link_cache.sh``): arm ``G`` of step 9 (``run.json``, ``e1_log.json.gz``,
``mc_traces.json.gz``, ``passes_<s>/D4_G/pass.json``) and step 11's arm ``S`` ``e1_log.json.gz``
(rule ``P(3000)``: its cut end is the best at the first 3,000-evaluation stall, so the hit ratio of
the running best at that stall is on record there only). It

1. checks the stored counts of the plan's stop rule (totals with the pass 13,004,457 / 10,196,692 /
   3,807,747; seed steps 5,417,024 / 1,528,284 / 687,902; traces; steps) and exits 2 on a
   difference (``counts.json`` holds the numbers either way);
2. writes ``whatif.json`` and ``whatif.md`` (blocks for ``sync_doc_tables.py``): the rule table
   (``P(3000)`` as step 11 measured it, gated high = this rule, gated low = its mirror), the cases
   behind the ``mid`` row of gated high, and the table of other k.

Hit ratio at a stall, as far as it is on record. The last gain's best is G's own end (the fit's
stored hit). The first stall before a later gain is read from step 11's arm S: the same candidate
(voxel, occurrence of the seed step at that voxel, candidate ``j``) cut at ``n + k + 1`` evaluations
with the cost of that gain. Any other stall is unknown. Gated high resolves an unknown stall of a
refinement whose G end is at 0.95 or more two ways (goes on: the low saving and the known changed
ends; stops: the high saving and the bound on changed ends); an unknown stall of a refinement
whose G end is under 0.95 is read as under 0.95 (the plan's reading: in every known case a stall
best at 0.95 belongs to a G end at 0.95). Gated low reads an unknown stall as the class of G's end.
"""

import argparse
import collections
import gzip
import json
import os
import sys
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
SH = ICE / "benchmarks" / "autoresearch" / "ar5" / "stallhit"
PAT = ICE / "benchmarks" / "autoresearch" / "ar5" / "patience"
CACHE = Path(os.environ.get("STALLHIT_CACHE", str(SH / "cache")))
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent.parent / "common"))

import doc_tables as DT  # noqa: E402
import patience_rules as PR  # noqa: E402  (step 11: cut_point, pure)
import stallhit_rules as SR  # noqa: E402

SAMPLES = ("mid", "grad", "big")
# stored counts of the plan's stop rule: totals (with the pass), seed steps, traces, steps
EXPECT: Dict[str, Dict[str, int]] = {
    "mid": dict(total=13004457, seed_evals=5417024, traces=1605, steps=1607),
    "grad": dict(total=10196692, seed_evals=1528284, traces=991, steps=1280),
    "big": dict(total=3807747, seed_evals=687902, traces=153, steps=317),
}
KS = (2000, 3000, 5000)
K = SR.K
TAU = SR.TAU


def _gz(p: Path) -> Any:
    with gzip.open(p, "rt") as f:
        return json.load(f)


def load(s: str) -> Dict[str, Any]:
    d = CACHE / f"G_{s}"
    run = json.loads((d / "run.json").read_text())
    pas = json.loads((CACHE / f"passes_{s}" / "D4_G" / "pass.json").read_text())
    return dict(
        run=run,
        pass_evals=int(pas["evals"]),
        e1=_gz(d / "e1_log.json.gz"),
        traces=_gz(d / "mc_traces.json.gz"),
        s_e1=_gz(CACHE / f"S_{s}" / "e1_log.json.gz"),
    )


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


def _occ(log: List[Dict[str, Any]]) -> List[Tuple[int, int]]:
    c: collections.Counter = collections.Counter()
    out = []
    for e in log:
        c[int(e["voxel"])] += 1
        out.append((int(e["voxel"]), c[int(e["voxel"])]))
    return out


def s_cut_fits(d: Dict[str, Any]) -> Dict[Tuple[int, int, int], List[Any]]:
    """Step 11's arm-S MC fits by (voxel, occurrence, candidate j): the fit row."""
    out: Dict[Tuple[int, int, int], List[Any]] = {}
    for key, e in zip(_occ(d["s_e1"]), d["s_e1"]):
        for f in e["fits"]:
            if f[2] == "M":
                out[(key[0], key[1], int(f[1]))] = f
    return out


def steps_with_mc(d: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Per seed step with an MC fit: the stored fits (dicts), the traces per position and the
    arm-S fit of each MC candidate (None when the step has no counterpart in S)."""
    by_step: Dict[int, Dict[int, Dict[str, Any]]] = {}
    for t in d["traces"]:
        by_step.setdefault(int(t["step"]), {})[int(t["j"])] = t
    sfit = s_cut_fits(d)
    s_keys = {k[:2] for k in sfit}
    out = []
    for si, (key, e) in enumerate(zip(_occ(d["e1"]), d["e1"]), start=1):
        if not e["fits"]:
            continue
        fits = [
            dict(pos=f[0], j=f[1], kind=f[2], hit=f[3], cost=f[4], evals=f[5]) for f in e["fits"]
        ]
        if not any(f["kind"] == "M" for f in fits):
            continue
        tr: Dict[int, Dict[str, Any]] = {}
        sf: Dict[int, Optional[List[Any]]] = {}
        for f in fits:
            if f["kind"] == "M":
                t = by_step[si][f["j"]]
                assert int(t["voxel"]) == int(e["voxel"]) and int(t["evals"]) == f["evals"]
                tr[f["pos"]] = t
                sf[f["pos"]] = sfit.get(key + (int(f["j"]),)) if key in s_keys else None
        out.append(dict(step=si, voxel=int(e["voxel"]), entry=e, fits=fits, traces=tr, s=sf))
    return out


def hit_oracle(
    trace: List[List[float]], fit: Dict[str, Any], sfit: Optional[List[Any]], k: int
) -> Callable[[int], Optional[float]]:
    """``hit_of(i)`` for ``stallhit_rules.refine_walk``: the hit ratio of the best after the
    ``i``-th gain, as far as it is on record (see the module text); None when unknown."""
    last = len(trace) - 1
    ev = SR.stall_events(trace, int(fit["evals"]), k)
    first = ev[0]["i"] if ev else None

    def hit_of(i: int) -> Optional[float]:
        if i == last:
            return float(fit["hit"])
        if k == K and i == first and sfit is not None:
            # step 11's S cut this refinement at the first stall: its end is the best there
            if (
                int(sfit[5]) == int(trace[i][0]) + k + 1
                and abs(float(sfit[4]) - float(trace[i][1])) <= 1e-9
            ):
                return float(sfit[3])
        return None

    return hit_of


def _resolve(
    hit_of: Callable[[int], Optional[float]], g_hit: float, mode: str
) -> Callable[[int], Optional[float]]:
    """Gated high: an unknown stall of a refinement whose G end is under 0.95 is read as under
    0.95 (hit 0); otherwise left unknown for ``refine_walk``'s ``unknown`` argument."""

    def f(i: int) -> Optional[float]:
        h = hit_of(i)
        if h is None and g_hit < TAU:
            return 0.0
        return h

    return f


def walk_rule(
    st: Dict[str, Any], k: int, unknown: str, gated: str
) -> Tuple[Dict[str, Any], Dict[int, Dict[str, Any]]]:
    """The whole seed step under gated-high (``"high"``, this rule) or gated-low (``"low"``)."""
    cuts: Dict[int, Dict[str, Any]] = {}
    for p, t in st["traces"].items():
        f = [x for x in st["fits"] if x["pos"] == p and x["kind"] == "M"][0]
        ho = hit_oracle(t["trace"], f, st["s"].get(p), k)
        g_hit = float(f["hit"])
        if gated == "high":
            cuts[p] = SR.refine_walk(
                t["trace"], int(f["evals"]), k, _resolve(ho, g_hit, unknown), unknown=unknown
            )
        else:  # low: stop when the best is under 0.95 (mirror); unknown read as G end's class

            def lo(i: int, ho: Callable[[int], Optional[float]] = ho, g: float = g_hit) -> float:
                h = ho(i)
                h = g if h is None else h
                return 1.0 if h < TAU else 0.0

            cuts[p] = SR.refine_walk(t["trace"], int(f["evals"]), k, lo, unknown="go")
            # a stop here is at "under 0.95": the kept end's hit class is under 0.95
            if cuts[p]["kind"] == "moved":
                cuts[p]["hit"] = 0.0
    n_avail = max(5, 1 + max(f["pos"] for f in st["fits"]))
    return SR.walk_step(st["fits"], cuts, n_avail=n_avail), cuts


def evaluate(
    sts: List[Dict[str, Any]], total: int, k: int, unknown: str, gated: str
) -> Dict[str, Any]:
    saved = changed = lost = extra = n_cut = 0
    for st in sts:
        w, cuts = walk_rule(st, k, unknown, gated)
        stored = SR.walk_step(st["fits"], {})
        saved += stored["evals"] - w["evals"]
        extra += w["extra_evals"]
        n_cut += sum(1 for c in cuts.values() if c["kind"] != "none")
        t = w["taken"]
        end = w["kept"][t]
        ch = (
            t != stored["taken"]
            or end["which"] != stored["kept"][stored["taken"]]["which"]
            or end["moved"]
        )
        changed += bool(ch)
        lost += bool(w["lost"])
    return dict(
        saved_evals=int(saved),
        saved_share=saved / total,
        saved_share_without_extra=(saved + extra) / total,
        extra_evals=int(extra),
        taken_changed=int(changed),
        stops_lost=int(lost),
        refinements_cut=int(n_cut),
    )


def cases(sts: List[Dict[str, Any]], total: int) -> Dict[str, Any]:
    """The cases behind gated high at k = 3000 (G end at 0.95 or more): A cut after the last gain
    (the first stall is the last), B a first stall known under 0.95 then a stop at a later stall,
    C a first stall known at 0.95 (stops before a later gain), D a first stall unknown, E a first
    stall known under 0.95 and no stop after it. Each: refinements, evaluations saved, share."""
    acc: Dict[str, Dict[str, Any]] = {
        c: dict(n=0, saved=0, taken=0) for c in ("A", "B", "C", "D", "E")
    }
    for st in sts:
        w0 = SR.walk_step(st["fits"], {})
        for p, t in st["traces"].items():
            f = [x for x in st["fits"] if x["pos"] == p and x["kind"] == "M"][0]
            if float(f["hit"]) < TAU:
                continue
            ev = SR.stall_events(t["trace"], int(f["evals"]), K)
            if not ev:
                continue
            ho = hit_oracle(t["trace"], f, st["s"].get(p), K)
            taken = bool(p == w0["taken"] and w0["kept"][p]["which"] == "M")
            if len(ev) == 1 and ev[0]["last"]:
                c = "A"
                saved = int(f["evals"]) - min(int(f["evals"]), ev[0]["n"] + K + 1)
            else:
                h0 = ho(ev[0]["i"])
                if h0 is None:
                    c = "D"
                elif h0 >= TAU:
                    c = "C"
                else:
                    c = "B"
                if c == "B":
                    r = SR.refine_walk(t["trace"], int(f["evals"]), K, ho)
                    c = "B" if r["kind"] != "none" else "E"
                    saved = int(f["evals"]) - r["evals"] if c == "B" else 0
                elif c == "C":
                    saved = int(f["evals"]) - min(int(f["evals"]), ev[0]["n"] + K + 1)
                else:  # D: saving if it stops at its first stall
                    saved = int(f["evals"]) - min(int(f["evals"]), ev[0]["n"] + K + 1)
            acc[c]["n"] += 1
            acc[c]["saved"] += saved
            acc[c]["taken"] += int(taken)
    for c in acc:
        acc[c]["share"] = acc[c]["saved"] / total
    return acc


def k_table(sts: List[Dict[str, Any]], total: int) -> Dict[str, Any]:
    """Other k, gated high. Exact saving: refinements cut after their last gain (step 11's
    ``cut_point``) whose G end is at 0.95 or more. Bound on changed ends: refinements cut before a
    later gain whose G end is at 0.95 or more."""
    out: Dict[str, Any] = {}
    for k in KS:
        exact = bound = 0
        for st in sts:
            for p, t in st["traces"].items():
                f = [x for x in st["fits"] if x["pos"] == p and x["kind"] == "M"][0]
                if float(f["hit"]) < TAU:
                    continue
                cp = PR.cut_point(t["trace"], int(f["evals"]), k)
                if cp["kind"] == "after_last":
                    exact += int(f["evals"]) - int(cp["evals"])
                elif cp["kind"] == "moved":
                    bound += 1
        out[str(k)] = dict(exact_saved=exact, exact_share=exact / total, bound_changed=bound)
    return out


def write_md(out: Dict[str, Any], samples: List[str]) -> None:
    def pct(x: float, nd: int = 2) -> str:
        return f"{100 * x:.{nd}f}%"

    def rng(lo: float, hi: float, nd: int = 1) -> str:
        a, b = pct(lo, nd), pct(hi, nd)
        return a if a == b else f"{a} to {b}"

    def ichg(lo: int, hi: int) -> str:
        return str(lo) if lo == hi else f"{lo} to {hi}"

    rule_rows = []
    row: Dict[str, Any] = {"rule (k = 3000)": "`P(3000)`, step 11, measured whole run (changed)"}
    for s in samples:
        m = out[s]["step11_measured"]
        row[s] = f"{m['total_cut_pct']:.2f}%; {m['taken_changed_n']}"
    rule_rows.append(row)
    row = {"rule (k = 3000)": "gated high: cut only if the best is at 0.95 or more"}
    for s in samples:
        lo, hi = out[s]["high_go"], out[s]["high_stop"]
        row[s] = (
            f"{rng(lo['saved_share_without_extra'], hi['saved_share_without_extra'])}; "
            f"{ichg(lo['taken_changed'], hi['taken_changed'])}"
        )
    rule_rows.append(row)
    row = {"rule (k = 3000)": "gated low: cut only if the best is under 0.95"}
    for s in samples:
        lw = out[s]["low"]
        row[s] = f"{pct(lw['saved_share_without_extra'], 1)}; {lw['taken_changed']}"
    rule_rows.append(row)

    case_rows = []
    for c, label in (
        ("A", "A: cut after the last gain, G end at 0.95 or more (end unchanged)"),
        ("B", "B: first stall known under 0.95, stops at a later stall"),
        ("C", "C: first stall known at 0.95, stops before a later gain"),
        ("D", "D: first stall unknown, G end at 0.95 or more"),
        ("E", "E: first stall known under 0.95, no stop after it"),
    ):
        r: Dict[str, Any] = {"case": label}
        for s in samples:
            x = out[s]["cases"][c]
            r[s] = f"{x['n']}; {x['saved']:,} ({pct(x['share'])})"
        case_rows.append(r)

    k_rows = []
    for k in KS:
        r2: Dict[str, Any] = {"k": str(k)}
        for s in samples:
            x = out[s]["k_table"][str(k)]
            r2[s] = f"{pct(x['exact_share'])}; {x['bound_changed']}"
        k_rows.append(r2)
    tabs = {
        "stallhit_whatif_rules": DT.markdown_table(rule_rows, ["rule (k = 3000)"] + list(samples)),
        "stallhit_whatif_cases": DT.markdown_table(case_rows, ["case"] + list(samples)),
        "stallhit_whatif_k": DT.markdown_table(k_rows, ["k"] + list(samples)),
    }
    DT.write_tables(SH / "whatif.md", tabs)
    body = (SH / "whatif.md").read_text()
    head = (
        "What-if tables of AR5 step 12 (arithmetic on the stored arm-G traces of step 9 and step\n"
        "11's arm S; a cost estimate, not a measurement). Written by `stallhit_whatif.py`.\n\n"
        "Reading line, `stallhit_whatif_rules`: one row per stop rule; each cell is the share of\n"
        "arm G's run (with the pass) saved WITHOUT the estimated extra candidates of a step that\n"
        "loses its 0.95 stop (the plan's basis; `whatif.json` also holds the share with them),\n"
        "then the seed steps whose taken end changes (a range where some stall hit ratios are not\n"
        "on record: from the plain reading to the reading where every unknown stall at a G end of\n"
        "0.95 or more stops). Step 11's row is the measured whole run\n"
        "(`patience/stats_*.json`).\n\n"
        "Reading line, `stallhit_whatif_cases`: refinements, then the evaluations they save, then\n"
        "that as a share of G's run. Cases are exclusive and cover the refinements whose G end is\n"
        "at 0.95 or more and that have a stall at k = 3000 (see the script text).\n\n"
        "Reading line, `stallhit_whatif_k`: other k, gated high; each cell is the exact saving\n"
        "(refinements cut after their last gain, G end at 0.95 or more), then the refinements cut\n"
        "before a later gain with a G end at 0.95 or more (the count of ends that could change).\n"
        "\n"
    )
    (SH / "whatif.md").write_text(head + body)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", nargs="*", default=list(SAMPLES))
    ap.add_argument("--md-only", action="store_true")
    a = ap.parse_args()
    if a.md_only:
        stored = json.loads((SH / "whatif.json").read_text())
        write_md(stored, [s for s in SAMPLES if s in stored])
        print("whatif.md rewritten from whatif.json")
        return
    SH.mkdir(parents=True, exist_ok=True)
    out: Dict[str, Any] = {}
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
        st11 = json.loads((PAT / f"stats_{s}.json").read_text())
        out[s] = dict(
            steps_with_mc=len(sts),
            mc_refinements=sum(len(x["traces"]) for x in sts),
            step11_measured=dict(
                total_cut_pct=float(st11["rules"]["total_cut_pct"]),
                taken_changed_n=int(st11["taken_changed_n"]),
            ),
            high_go=evaluate(sts, total, K, "go", "high"),
            high_stop=evaluate(sts, total, K, "stop", "high"),
            low=evaluate(sts, total, K, "go", "low"),
            cases=cases(sts, total),
            k_table=k_table(sts, total),
        )
    (SH / "counts.json").write_text(json.dumps(cnt, indent=1) + "\n")
    out["counts_ok"] = not bad
    (SH / "whatif.json").write_text(json.dumps(out, indent=1) + "\n")
    if bad:
        print("STOP: a stored count does not reproduce:", bad, file=sys.stderr)
        sys.exit(2)
    write_md(out, list(a.sample))
    print("whatif written")


if __name__ == "__main__":
    main()
