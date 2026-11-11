"""AR3 step 6, parts 0 and 1: summary of the replay (tables, plots, summary.json).

  uv run python scripts/autoresearch/ar3/seed_summary.py [--slim]

``--slim`` first turns ``seed/cache/part01.jsonl`` into ``seed/part01.json.gz`` (the only per-voxel
file that is committed: floats rounded to 6 digits, the 50 candidate errors dropped). The summary
reads that file, ``seed/voxel_facts.json`` and the stored ``run.json`` of D4. Wilson 95% intervals,
no p-values. Truth is used for analysis only. Sets: A is the fallback voxels whose true grain id is
even, B odd (step 3); the pick is made on A and reported on B.
"""

import argparse
import gzip
import json
import sys
from itertools import permutations
from typing import Any, Callable, Dict, List, Optional, Tuple

import ar3_common as C
import gradient_common as G
import numpy as np
import seed_common as S

sys.path.insert(0, str(C.ICE / "scripts" / "common"))
import doc_tables as DT  # noqa: E402
import stats as ST  # noqa: E402

WRONG = G.WRONG
REMEDIES = ["more20_2", "more20_5", "cma_seed", "wide", "neighbours", "neighbours_only", "outward"]
EXTRAS = ["order2", "order5"]  # added beyond the plan: the indexer's own order, no strict cost
NAME = {
    "order2": "7a order, first 2 (added)",
    "order5": "7b order, first 5 (added)",
    "more20_2": "1a more, best 2",
    "more20_5": "1b more, best 5",
    "cma_seed": "2 cma-seed",
    "wide": "3 wide",
    "neighbours": "4 neighbours",
    "neighbours_only": "5 defer (tries only)",
    "outward": "6 outward",
}
# the remedy numbers the findings page uses for the code names
PAGE = {
    "more20_2": "1a",
    "more20_5": "1b",
    "cma_seed": "2",
    "wide": "3",
    "neighbours": "4",
    "neighbours_only": "5",
    "outward": "6",
    "order2": "7a",
    "order5": "7b",
}
ORDER = ("D4", "D", "P4", "P")
GRID_EVALS = 339_914  # D4: index_fallback_grid_evals / index_fallback_seed_test (25,833,494 / 76)


def slim() -> None:
    src = S.CACHE / "part01.jsonl"
    out = []

    def rnd(x: Any) -> Any:
        if isinstance(x, float):
            return round(x, 6)
        if isinstance(x, list):
            return [rnd(i) for i in x]
        if isinstance(x, dict):
            return {k: rnd(v) for k, v in x.items()}
        return x

    for line in src.read_text().splitlines():
        r = json.loads(line)
        if "cand" in r:
            r["cand"].pop("errs", None)
        out.append(rnd(r))
    with gzip.open(S.SEED / "part01.json.gz", "wt", compresslevel=9) as f:
        json.dump(out, f, separators=(",", ":"))
    print("slim", len(out), "rows")


def wil(k: int, n: int) -> Dict[str, Any]:
    lo, hi = ST.wilson(int(k), int(n)) if n else (float("nan"), float("nan"))
    return dict(k=int(k), n=int(n), lo=float(lo), hi=float(hi))


def ws(d: Dict[str, Any]) -> str:
    if not d["n"]:
        return "0 of 0"
    return f"{d['k']} of {d['n']} ({100 * d['lo']:.0f} to {100 * d['hi']:.0f}%)"


def page_name(chain_name: str) -> str:
    """'neighbours_only then more20_5' -> '5 then 1b' (the remedy numbers of the page)."""
    return " then ".join(PAGE[p] for p in chain_name.split(" then "))


def ref_run(r: Dict[str, Any]) -> str:
    return next(x for x in ORDER if x in r["runs"])


def grid_end(r: Dict[str, Any]) -> Dict[str, Any]:
    return r["stored"][ref_run(r)]


def rem_of(r: Dict[str, Any], name: str) -> Dict[str, Any]:
    x = r["remedies"][name]
    if name in ("neighbours", "neighbours_only"):
        return x[ref_run(r)] if ref_run(r) in x else x[next(iter(x))]
    return x


def grid_cost(r: Dict[str, Any]) -> int:
    """Evaluations the stored run spent at the voxel (seed step plus grid search)."""
    f = r["runs"][ref_run(r)].get("full_evals", 0)
    return int(f) if f else int(r["stock"]["evals"] + GRID_EVALS)


def right(e: Optional[float]) -> bool:
    return e is not None and e <= WRONG


def one(r: Dict[str, Any], name: str) -> Dict[str, Any]:
    """The remedy at a fallback voxel: the end it hands over (its own if it passes the seed test,
    else the grid search's), its evaluations alone and with the grid search when it fails."""
    x = rem_of(r, name)
    g = grid_end(r)
    ev = int(x["evals"])
    ok = bool(x["passed"])
    return dict(
        passed=ok,
        err=x["err"] if ok else g["err"],
        own_err=x["err"],
        evals=ev,
        total=ev + (0 if ok else GRID_EVALS),
        seed_in=bool(x.get("seed_in")),
    )


def chain(r: Dict[str, Any], a: str, b: str) -> Dict[str, Any]:
    xa, xb = one(r, a), one(r, b)
    if xa["passed"]:
        return dict(xa)
    ev = xa["evals"] + xb["evals"] - (r["stock"]["evals"] if xa["seed_in"] and xb["seed_in"] else 0)
    ok = xb["passed"]
    g = grid_end(r)
    return dict(
        passed=ok,
        err=xb["own_err"] if ok else g["err"],
        own_err=xb["own_err"],
        evals=ev,
        total=ev + (0 if ok else GRID_EVALS),
        seed_in=True,
    )


def stat_set(
    rows: List[Dict[str, Any]], get: Callable[[Dict[str, Any]], Dict[str, Any]]
) -> Dict[str, Any]:
    n = len(rows)
    xs = [get(r) for r in rows]
    gr = [right(grid_end(r)["err"]) for r in rows]
    er = [right(x["err"]) for x in xs]
    ev = np.array([x["evals"] for x in xs], dtype=float) if n else np.zeros(0)
    tot = np.array([x["total"] for x in xs], dtype=float) if n else np.zeros(0)
    return dict(
        n=n,
        passed=wil(sum(x["passed"] for x in xs), n),
        need_grid=wil(sum(not x["passed"] for x in xs), n),
        pass_and_right=sum(x["passed"] and right(x["own_err"]) for x in xs),
        right_after=wil(sum(er), n),
        grid_right=wil(sum(gr), n),
        fixed=sum((not g) and e for g, e in zip(gr, er)),
        broken=sum(g and (not e) for g, e in zip(gr, er)),
        both_right=sum(g and e for g, e in zip(gr, er)),
        both_wrong=sum((not g) and (not e) for g, e in zip(gr, er)),
        ev_median=float(np.median(ev)) if n else None,
        ev_mean=float(ev.mean()) if n else None,
        ev_max=float(ev.max()) if n else None,
        total_mean=float(tot.mean()) if n else None,
        total_median=float(np.median(tot)) if n else None,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--slim", action="store_true")
    a = ap.parse_args()
    if a.slim:
        slim()
    with gzip.open(S.SEED / "part01.json.gz", "rt") as f:
        rows = json.load(f)
    facts = json.loads((S.SEED / "voxel_facts.json").read_text())
    fb = [r for r in rows if r["kind"] == "fallback"]
    ctl = [r for r in rows if r["kind"] == "control"]
    small = [r for r in rows if r["kind"] == "small"]
    for r in fb + ctl:
        r.update(facts[str(r["voxel"])])
    A = [r for r in fb if r["isA"]]
    B = [r for r in fb if not r["isA"]]
    SUM: Dict[str, Any] = dict(
        n_fallback_voxels=len(fb), n_A=len(A), n_B=len(B), n_controls=len(ctl)
    )
    T: Dict[str, str] = {}

    # ---------------- part 0 item 1: the lists
    d4 = json.loads((C.OUT / "fine" / "runs" / "D4" / "run.json").read_text())["stats"]
    lists = json.loads((S.SEED / "lists.json").read_text())
    SUM["lists"] = lists
    rows_l = []
    for run in ORDER:
        L = lists[run]
        rows_l.append(
            dict(
                run=run,
                listed=L["listed"],
                counter=L.get("counter"),
                at_seed=L.get("listed_seed"),
                at_reseed=L.get("listed_reseed"),
                at_paint=L.get("listed") if run in ("P", "P4") else None,
            )
        )
    T["ar3f_lists"] = DT.markdown_table(
        rows_l, ["run", "listed", "counter", "at_seed", "at_reseed", "at_paint"]
    )

    # D4 fallbacks by cluster size, and every D4 reseed by cluster size
    def cbin(c: int) -> str:
        return "1" if c == 1 else "2" if c == 2 else "3 to 5" if c <= 5 else "6 or more"

    bins = ["1", "2", "3 to 5", "6 or more"]
    rs = d4["reseeds"]
    rows_c = []
    for b in bins:
        sel = [x for x in rs if cbin(x["cluster_size"]) == b]
        f = [x for x in sel if x["full_search_evaluations"] > S.BIG]
        rows_c.append(
            dict(
                cluster=b,
                reseeds=len(sel),
                fallbacks=len(f),
                evals_all=sum(x["full_search_evaluations"] for x in sel),
                evals_fallbacks=sum(x["full_search_evaluations"] for x in f),
                evals_per_voxel=sum(x["full_search_evaluations"] for x in sel)
                / sum(x["cluster_size"] for x in sel),
            )
        )
    SUM["d4_cluster"] = rows_c
    T["ar3f_cluster"] = DT.markdown_table(
        rows_c,
        ["cluster", "reseeds", "fallbacks", "evals_all", "evals_fallbacks", "evals_per_voxel"],
        formats=dict(evals_all=",d", evals_fallbacks=",d", evals_per_voxel=",.0f"),
    )

    # by level and boundary, unique fallback voxels (all runs)
    rows_lv = []
    for L in G.LEVELS:
        sel = [r for r in fb if r["level"] == L]
        rows_lv.append(
            dict(
                level=f"{L:g}",
                voxels=len(sel),
                at_boundary=sum(r["boundary"] for r in sel),
                inside=sum(not r["boundary"] for r in sel),
            )
        )
    rows_lv.append(
        dict(
            level="all",
            voxels=len(fb),
            at_boundary=sum(r["boundary"] for r in fb),
            inside=sum(not r["boundary"] for r in fb),
        )
    )
    for L in G.LEVELS:
        sel = [r for r in ctl if r["level"] == L]
        rows_lv.append(
            dict(
                level=f"controls, level {L:g}",
                voxels=len(sel),
                at_boundary=sum(r["boundary"] for r in sel),
                inside=sum(not r["boundary"] for r in sel),
            )
        )
    rows_lv.append(
        dict(
            level="controls, all",
            voxels=len(ctl),
            at_boundary=sum(r["boundary"] for r in ctl),
            inside=sum(not r["boundary"] for r in ctl),
        )
    )
    SUM["by_level"] = rows_lv
    T["ar3f_levels"] = DT.markdown_table(rows_lv, ["level", "voxels", "at_boundary", "inside"])

    # ---------------- part 0 item 2: the check
    f_fail = sum(not r["stock"]["passed"] for r in fb)
    c_pass = sum(r["stock"]["passed"] for r in ctl)
    c_evals_same = sum(r["stock"]["evals"] == r["runs"]["D4"]["stored_evals"] for r in ctl)
    chk = dict(
        fallbacks=wil(f_fail, len(fb)),
        controls=wil(c_pass, len(ctl)),
        controls_evals_equal=wil(c_evals_same, len(ctl)),
        disagree_fallbacks=len(fb) - f_fail,
        disagree_controls=len(ctl) - c_pass,
        stop=bool((len(fb) - f_fail) > 0.05 * len(fb) or (len(ctl) - c_pass) > 0.05 * len(ctl)),
        same_first20=wil(sum(bool(r["cand"].get("same_first20")) for r in fb + ctl), len(fb + ctl)),
    )
    SUM["check2"] = chk
    T["ar3f_check"] = DT.markdown_table(
        [
            dict(
                group="fallback voxels (stored: fail)",
                n=len(fb),
                agree=f_fail,
                disagree=len(fb) - f_fail,
            ),
            dict(
                group="control index seeds (stored: pass)",
                n=len(ctl),
                agree=c_pass,
                disagree=len(ctl) - c_pass,
            ),
        ],
        ["group", "n", "agree", "disagree"],
    )

    # ---------------- part 0 item 4: classes
    classes = list("abcdef")
    rows_k = []
    for nm, sel in (("all fallback voxels", fb), ("set A", A), ("set B", B)):
        d = {c: sum(r["class"] == c for r in sel) for c in classes}
        rows_k.append(dict(set=nm, n=len(sel), **d))
    SUM["classes"] = rows_k
    T["ar3f_classes"] = DT.markdown_table(rows_k, ["set", "n"] + classes)
    SUM["gate_a_more_than_half"] = bool(sum(r["class"] == "a" for r in fb) > len(fb) / 2)

    # item 3 and 6: rank of the first candidate within 1 degree; hit ratio of the truth and of the
    # first candidates
    def rank_row(nm: str, sel: List[Dict[str, Any]]) -> Dict[str, Any]:
        rk = [r["cand"]["rank1"] for r in sel]
        found = [x for x in rk if x is not None]
        return dict(
            set=nm,
            n=len(sel),
            none=sum(x is None for x in rk),
            rank0=sum(x == 0 for x in found),
            rank1_4=sum(1 <= x <= 4 for x in found),
            rank5_19=sum(5 <= x <= 19 for x in found),
            rank20_49=sum(x >= 20 for x in found),
            median_cand=float(np.median([r["cand"]["n"] for r in sel])),
            truth_hit_below=sum(r["truth"]["hit"] < S.TEST for r in sel),
            median_truth_hit=float(np.median([r["truth"]["hit"] for r in sel])),
            median_first_hit=float(
                np.median([r["cand"]["first"][0]["hit"] for r in sel if r["cand"]["first"]])
            ),
        )

    rows_r = [rank_row("fallback voxels", fb), rank_row("controls", ctl)]
    SUM["rank"] = rows_r
    T["ar3f_rank"] = DT.markdown_table(
        rows_r,
        [
            "set",
            "n",
            "none",
            "rank0",
            "rank1_4",
            "rank5_19",
            "rank20_49",
            "median_cand",
            "truth_hit_below",
            "median_truth_hit",
            "median_first_hit",
        ],
        formats=dict(median_cand=".0f", median_truth_hit=".3f", median_first_hit=".3f"),
    )
    # hit ratio of the first within-1-degree candidate against the best of the first 5
    gap = []
    for r in fb:
        k = r["cand"]["rank1"]
        if k is not None and k < len(r["cand"]["first"]):
            fh = [c["hit"] for c in r["cand"]["first"][:5]]
            gap.append(
                (
                    r["cand"]["first"][k]["hit"],
                    max(fh),
                    r["cand"]["first"][k]["cost"],
                    min(c["cost"] for c in r["cand"]["first"][:5]),
                )
            )
    SUM["raw_near_candidate"] = dict(
        n=len(gap),
        median_hit_near=float(np.median([g[0] for g in gap])) if gap else None,
        median_hit_best5=float(np.median([g[1] for g in gap])) if gap else None,
        near_is_lowest_cost_of_5=sum(g[2] <= g[3] + 1e-12 for g in gap),
    )
    # where the first candidate within 1 degree of the truth falls in the strict-cost order of the
    # first 5 (0 = lowest cost, so it would be refined): fallback voxels and controls
    for nm, sel in (("fallbacks", fb), ("controls", ctl)):
        pos = []
        for r in sel:
            k = r["cand"]["rank1"]
            if k is not None and k < 5 and len(r["cand"]["first"]) >= 5:
                cs = [c["cost"] for c in r["cand"]["first"][:5]]
                pos.append(int(np.argsort(cs, kind="stable").tolist().index(k)))
        SUM[f"near_cost_pos_{nm}"] = dict(
            n=len(pos),
            refined=sum(p < 2 for p in pos),
            lowest=sum(p == 0 for p in pos),
            # added in the review fixes (D-20261111-N1): Wilson 95% interval of ``refined``
            refined_w=wil(sum(p < 2 for p in pos), len(pos)),
        )
        # error of that candidate (degrees) and its raw hit ratio, by group
        ne = [
            (
                r["cand"]["first"][r["cand"]["rank1"]]["err"],
                r["cand"]["first"][r["cand"]["rank1"]]["hit"],
            )
            for r in sel
            if r["cand"]["rank1"] is not None and r["cand"]["rank1"] < len(r["cand"]["first"])
        ]
        SUM[f"near_err_{nm}"] = dict(
            n=len(ne),
            median_err=float(np.median([x[0] for x in ne])),
            q25_err=float(np.percentile([x[0] for x in ne], 25)),
            q75_err=float(np.percentile([x[0] for x in ne], 75)),
            median_hit=float(np.median([x[1] for x in ne])),
            hit_ge_080=sum(x[1] >= S.TEST for x in ne),
        )
    # the refined ends of the stock replay: how close to the truth, and their hit ratio
    fe = [(e, h) for r in fb for e, h in zip(r["stock"]["fin_err"], r["stock"]["fin_hit"])]
    near_h = [h for e, h in fe if e <= WRONG]  # the refined ends within 1 degree of the truth
    SUM["stock_ends"] = dict(
        n=len(fe),
        within1=sum(e <= WRONG for e, _ in fe),
        median_hit=float(np.median([h for _, h in fe])) if fe else None,  # over all n ends
        # added in the review fixes (D-20261111-N1): the ends within 1 degree alone
        within1_median_hit=float(np.median(near_h)) if near_h else None,
        within1_max_hit=float(max(near_h)) if near_h else None,
        within1_below_test=sum(h < S.TEST for h in near_h),
    )

    # ---------------- part 0 item 5: the grid end and the neighbours
    rows_g = []
    for run in ORDER:
        sel = [r for r in fb if run in r["runs"]]
        st = [r["stored"][run] for r in sel]
        rows_g.append(
            dict(
                run=run,
                n=len(sel),
                right=wil(sum(right(s["err"]) for s in st), len(sel)),
                hit_ge_080=sum(s["hit"] >= S.TEST for s in st),
                median_err=float(np.median([s["err"] for s in st])),
                with_nbr_within1=sum(s["n_nbr_within1"] >= 1 for s in st),
                median_nbr_within1=float(np.median([s["n_nbr_within1"] for s in st])),
                median_nbr=float(np.median([s["n_nbr"] for s in st])),
            )
        )
    SUM["grid_end"] = rows_g
    for d in rows_g:
        d["right_s"] = ws(d["right"])
    T["ar3f_gridend"] = DT.markdown_table(
        rows_g,
        [
            "run",
            "n",
            "right_s",
            "hit_ge_080",
            "median_err",
            "with_nbr_within1",
            "median_nbr_within1",
            "median_nbr",
        ],
        formats=dict(median_err=".4f", median_nbr_within1=".0f", median_nbr=".0f"),
    )

    # ---------------- part 1
    sets = {"A": A, "B": B, "A+B": fb}
    P1: Dict[str, Any] = {}
    for sn, sel in sets.items():
        P1[sn] = {m: stat_set(sel, lambda r, m=m: one(r, m)) for m in REMEDIES + EXTRAS}
        # the stock seed step followed by the grid search (the stored behaviour)
    SUM["remedies"] = P1

    def rem_table(sn: str) -> str:
        out = []
        for m in REMEDIES + EXTRAS:
            s = P1[sn][m]
            out.append(
                dict(
                    remedy=NAME[m],
                    n=s["n"],
                    pass_s=ws(s["passed"]),
                    need_grid=s["need_grid"]["k"],
                    right_s=ws(s["right_after"]),
                    fixed=s["fixed"],
                    broken=s["broken"],
                    both_right=s["both_right"],
                    both_wrong=s["both_wrong"],
                    ev_median=s["ev_median"],
                    ev_mean=s["ev_mean"],
                    ev_max=s["ev_max"],
                )
            )
        return DT.markdown_table(
            out,
            [
                "remedy",
                "n",
                "pass_s",
                "need_grid",
                "right_s",
                "fixed",
                "broken",
                "both_right",
                "both_wrong",
                "ev_median",
                "ev_mean",
                "ev_max",
            ],
            formats=dict(ev_median=",.0f", ev_mean=",.0f", ev_max=",.0f"),
        )

    T["ar3f_rem_A"] = rem_table("A")
    T["ar3f_rem_B"] = rem_table("B")

    # controls
    def ctl_one(r: Dict[str, Any], m: str) -> Dict[str, Any]:
        x = rem_of(r, m)
        ok = bool(x["passed"])
        st = r["stock"]
        end_err = x["err"] if ok else st["err"]
        return dict(
            passed=ok,
            err=end_err,
            ev_add=int(x["evals"]) - int(st["evals"]),
            stock_right=right(st["err"]),
            right=right(end_err),
        )

    rows_ct = []
    for m in REMEDIES + EXTRAS:
        xs = [ctl_one(r, m) for r in ctl]
        ea = np.array([x["ev_add"] for x in xs], dtype=float)
        rows_ct.append(
            dict(
                remedy=NAME[m],
                n=len(xs),
                pass_s=ws(wil(sum(x["passed"] for x in xs), len(xs))),
                stock_right=sum(x["stock_right"] for x in xs),
                broken=sum(x["stock_right"] and not x["right"] for x in xs),
                fixed=sum((not x["stock_right"]) and x["right"] for x in xs),
                add_median=float(np.median(ea)),
                add_mean=float(ea.mean()),
                add_max=float(ea.max()),
                # added in the review fixes (D-20261111-N1): Wilson 95% interval of ``broken``
                broken_w=wil(sum(x["stock_right"] and not x["right"] for x in xs), len(xs)),
            )
        )
    SUM["controls"] = rows_ct
    T["ar3f_controls"] = DT.markdown_table(
        rows_ct,
        [
            "remedy",
            "n",
            "pass_s",
            "stock_right",
            "broken",
            "fixed",
            "add_median",
            "add_mean",
            "add_max",
        ],
        formats=dict(add_median=",.0f", add_mean=",.0f", add_max=",.0f"),
    )

    # small clusters (question 3, remedy 5)
    rows_s = []
    for nm, sel in (
        (
            "cluster 1, fallback",
            [x for x in small if x["cluster"] == 1 and x["stored_evals"] > S.BIG],
        ),
        (
            "cluster 1, index seed used",
            [x for x in small if x["cluster"] == 1 and x["stored_evals"] <= S.BIG],
        ),
        (
            "cluster 2, fallback",
            [x for x in small if x["cluster"] == 2 and x["stored_evals"] > S.BIG],
        ),
        (
            "cluster 2, index seed used",
            [x for x in small if x["cluster"] == 2 and x["stored_evals"] <= S.BIG],
        ),
    ):
        ok = [
            x
            for x in sel
            if x["d_to_stored"] is not None
            and x["d_to_stored"] <= 0.1
            and x["hit"] >= x["stored_hit"] - 1e-9
        ]
        pr = [x for x in sel if x["hit"] >= S.TEST and x["err"] is not None and x["err"] <= WRONG]
        rows_s.append(
            dict(
                group=nm,
                n=len(sel),
                ok_s=ws(wil(len(ok), len(sel))),
                pass_right=ws(wil(len(pr), len(sel))),
                hit_gap_median=(
                    float(np.median([x["stored_hit"] - x["hit"] for x in sel])) if sel else None
                ),
                tries_pass=sum(x["hit"] >= S.TEST for x in sel),
                stored_ev=sum(x["stored_evals"] for x in ok),
                tries_ev=sum(x["evals"] for x in ok),
                saved=sum(x["stored_evals"] - x["evals"] for x in ok),
                stored_ev_all=sum(x["stored_evals"] for x in sel),
            )
        )
    SUM["small"] = rows_s
    T["ar3f_small"] = DT.markdown_table(
        rows_s,
        [
            "group",
            "n",
            "ok_s",
            "pass_right",
            "hit_gap_median",
            "stored_ev",
            "tries_ev",
            "saved",
            "stored_ev_all",
        ],
        formats=dict(
            stored_ev=",d", tries_ev=",d", saved=",d", stored_ev_all=",d", hit_gap_median=".3f"
        ),
    )

    # ---------------- the pick (on A; reported on B)
    gate = SUM["gate_a_more_than_half"]
    mean_ev = {m: P1["A"][m]["ev_mean"] for m in REMEDIES + EXTRAS}

    def make_pick(pool: List[str]) -> List[Dict[str, Any]]:
        allowed = ["neighbours", "neighbours_only"] if gate else list(pool)
        cheap = sorted(allowed, key=lambda m: mean_ev[m])
        cands: List[Tuple[str, Callable[[Dict[str, Any]], Dict[str, Any]]]] = [
            (m, lambda r, m=m: one(r, m)) for m in cheap
        ]
        for x, y in permutations(cheap, 2):
            if mean_ev[x] <= mean_ev[y]:
                cands.append((f"{x} then {y}", lambda r, x=x, y=y: chain(r, x, y)))
        table = []
        for nm, fn in cands:
            table.append(dict(name=nm, A=stat_set(A, fn), B=stat_set(B, fn)))
        ok = [t for t in table if t["A"]["broken"] <= t["A"]["fixed"]]
        ok.sort(key=lambda t: (t["A"]["need_grid"]["k"], t["A"]["total_mean"]))
        return ok

    okc = make_pick(REMEDIES)
    okx = make_pick(REMEDIES + EXTRAS)  # with the two added remedies; not the plan's pick
    table = okc
    SUM["pick_with_extras"] = dict(
        picked=okx[0]["name"] if okx else None,
        A=okx[0]["A"] if okx else None,
        B=okx[0]["B"] if okx else None,
        top=[t["name"] for t in okx[:6]],
    )
    rows_x = [
        dict(
            chain=t["name"],
            remedies=page_name(t["name"]),
            need_A=ws(t["A"]["need_grid"]),
            broken_A=t["A"]["broken"],
            fixed_A=t["A"]["fixed"],
            total_A=t["A"]["total_mean"],
            need_B=ws(t["B"]["need_grid"]),
            broken_B=t["B"]["broken"],
            fixed_B=t["B"]["fixed"],
            total_B=t["B"]["total_mean"],
        )
        for t in okx[:6]
    ]
    T["ar3f_pick_extras"] = DT.markdown_table(
        rows_x,
        [
            "remedies",
            "chain",
            "need_A",
            "broken_A",
            "fixed_A",
            "total_A",
            "need_B",
            "broken_B",
            "fixed_B",
            "total_B",
        ],
        formats=dict(total_A=",.0f", total_B=",.0f"),
    )
    pick = okc[0] if okc else None
    SUM["pick"] = dict(
        gate=gate,
        picked=None if pick is None else pick["name"],
        A=None if pick is None else pick["A"],
        B=None if pick is None else pick["B"],
        runner_up=[t["name"] for t in okc[1:4]],
        n_candidates=len(table),
        # added in the review fixes (D-20261111-N1): the pick in the page's remedy numbers
        picked_page=None if pick is None else page_name(pick["name"]),
    )
    removed = None if pick is None else 1 - pick["A"]["need_grid"]["k"] / max(len(A), 1)
    SUM["pick"]["share_removed_A"] = removed
    SUM["part2_runs"] = bool(removed is not None and removed >= 0.5)
    rows_p = []
    for t in okc[:6]:
        rows_p.append(
            dict(
                chain=t["name"],
                remedies=page_name(t["name"]),
                need_A=ws(t["A"]["need_grid"]),
                broken_A=t["A"]["broken"],
                fixed_A=t["A"]["fixed"],
                total_A=t["A"]["total_mean"],
                need_B=ws(t["B"]["need_grid"]),
                broken_B=t["B"]["broken"],
                fixed_B=t["B"]["fixed"],
                total_B=t["B"]["total_mean"],
            )
        )
    T["ar3f_pick"] = DT.markdown_table(
        rows_p,
        [
            "remedies",
            "chain",
            "need_A",
            "broken_A",
            "fixed_A",
            "total_A",
            "need_B",
            "broken_B",
            "fixed_B",
            "total_B",
        ],
        formats=dict(total_A=",.0f", total_B=",.0f"),
    )
    # stored: mean evaluations per voxel of the seed step followed by the grid search
    SUM["baseline"] = dict(
        mean_A=float(np.mean([grid_cost(r) for r in A])),
        mean_B=float(np.mean([grid_cost(r) for r in B])),
    )
    DT.write_tables(S.SEED / "tables_part01.md", T)
    (S.SEED / "summary_part01.json").write_text(json.dumps(SUM, indent=1, default=float) + "\n")
    figure(SUM, P1)
    print("pick", SUM["pick"]["picked"], "removed", removed, "gate", gate, "check", chk["stop"])


def figure(SUM: Dict[str, Any], P1: Dict[str, Any]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    items = [(NAME[m], P1["A"][m], P1["B"][m]) for m in REMEDIES + EXTRAS]
    pk = SUM["pick"]
    if pk["picked"]:
        # page names: '5 defer (tries only)' then '1b more, best 5' (fix 9 of the review; the old
        # label dropped '_only' and so named the chain of remedy 4, not the pick)
        parts = [NAME[p] for p in pk["picked"].split(" then ")]
        lab = "pick: " + " then\n".join(parts)
        items.append((lab, pk["A"], pk["B"]))
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    for j, (sn, col) in enumerate((("A", "tab:blue"), ("B", "tab:orange"))):
        for i, it in enumerate(items):
            s = it[1 + j]
            n = s["n"]
            p = s["need_grid"]["k"] / n if n else float("nan")
            x = i + (j - 0.5) * 0.3
            ax[0].errorbar(
                x,
                100 * p,
                yerr=[[100 * (p - s["need_grid"]["lo"])], [100 * (s["need_grid"]["hi"] - p)]],
                fmt="o",
                color=col,
                label=f"set {sn} (n={n})" if i == 0 else None,
            )
            ax[1].bar(x, s["total_mean"], width=0.28, color=col)
    ax[0].axhline(100, color="gray", lw=0.6)
    ax[0].set_ylabel("still sent to the grid search (%)")
    ax[0].set_title("Fallback voxels after each remedy (Wilson 95%); stored runs: 100%")
    ax[0].legend(loc="center right")
    ax[1].set_yscale("log")
    ax[1].set_ylim(1e3, 6e5)
    ax[1].set_ylabel("mean evaluations per voxel")
    ax[1].set_title("Cost, grid search included when needed (dashed: stored runs)")
    ax[1].axhline(SUM["baseline"]["mean_A"], color="gray", lw=0.8, ls="--")
    for a_ in ax:
        a_.set_xticks(range(len(items)))
        a_.set_xticklabels([it[0] for it in items], rotation=35, ha="right", fontsize=8)
    fig.tight_layout()
    fig.savefig(S.SEED / "ar3_seed.png", dpi=110)


if __name__ == "__main__":
    main()
