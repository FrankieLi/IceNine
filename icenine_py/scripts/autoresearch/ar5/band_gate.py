"""AR5 step 15, phase A step 5: the gate of both forms of rule ``PB(k)`` on the four disclosed
samples, the predictions 1 and 2, and the choice of the form (by the plan's fixed rule).

  uv run python scripts/autoresearch/ar5/band_gate.py [--no-choose]

Reads (read only) the recording rows (``band_record.py``) and the changed-path rows
(``band_paths.py``) in ``<cache>``, the stored arm ``J`` (``e1_log``, ``polish_log``) and step 14's
tracked ``stats_<s>.json``. Needs ``ar5/band/frozen.json`` committed and agreeing for ``held``
(``band_common``). Writes ``ar5/band/gate.json`` and ``ar5/band/gate_tables.md`` (blocks for
``sync_doc_tables.py``).

Per form and sample: stops (high, low, none), checks that went on, net saving (every extra
evaluation counted: checks, extra candidates, rescues, ``G2b``, grid searches; arithmetic for steps
whose path is unchanged, the live path for the others), changed taken ends (and, against the truth,
broken: ``J``'s end right and the new one wrong; fixed), grid searches, failed seed tests, lost
0.95 ends. The decision rules are ``band_rules.GATE`` (frozen). The net saving is a cost estimate on
``J``'s own steps: growth after a changed voxel is not replayed (the whole layer decides that).
"""

import argparse
import json
import re
import sys
from typing import Any, Dict, List, Tuple

import band_common as BC  # noqa: E402  (first)

import band_paths as BPA  # noqa: E402
import band_record as BREC  # noqa: E402
import band_rules as BR  # noqa: E402

sys.path.insert(0, str(BC.ICE / "scripts" / "common"))
import doc_tables as DT  # noqa: E402

SAMPLES = ("held", "mid", "grad", "big")
TOL = 1e-12


def read_paths(s: str, form: str) -> Dict[int, Dict[str, Any]]:
    rx = re.compile(rf"^paths_{form}_{s}_s\d+of\d+\.jsonl$")
    rows: Dict[int, Dict[str, Any]] = {}
    for p in sorted((BC.BCACHE / "paths").glob("paths_*.jsonl")):
        if rx.match(p.name):
            for line in BC.BG.guarded(p).read_text().splitlines():
                if line.strip():
                    r = json.loads(line)
                    rows[int(r["step"])] = r
    return rows


def prediction1(s: str, recs: List[Dict[str, Any]], d: Dict[str, Any]) -> Dict[str, Any]:
    c = BREC.counts(s, d)
    n_steps = len(BREC.mc_steps(d))
    fits = sum(len(r["refinements"]) for r in recs)
    ok = bool(
        c["ok"]
        and len(recs) == n_steps
        and all(r["p1_ok"] for r in recs)
        and all(x["trace_ok"] and x["checks_ok"] for r in recs for x in r["refinements"])
    )
    return dict(
        sample=s,
        counts=c,
        steps_expected=n_steps,
        steps_recorded=len(recs),
        refinements=fits,
        steps_not_equal=int(sum(1 for r in recs if not r["p1_ok"])),
        worst_abs_diff=float(max([r["worst_d"] for r in recs] + [0.0])),
        ok=ok,
    )


def bracketed(recs: List[Dict[str, Any]], k: int) -> Dict[str, int]:
    """The bracketed refinements at stall size ``k``: the first stall comes before a later gain
    (the best at the stall is not J's end). How many have a stall best under 0.80, in 0.80 to 0.95,
    at 0.95 or more."""
    out = dict(n=0, low=0, mid=0, high=0)
    for r in recs:
        for x in r["refinements"]:
            ev = x["events"].get(str(k), [])
            if not ev or float(ev[0][1]) <= float(x["final_cost"]) + TOL:
                continue
            out["n"] += 1
            h = float(ev[0][2])
            out["low" if h < BR.LOW else ("high" if h >= BR.TAU else "mid")] += 1
    return out


def one(
    s: str,
    form: str,
    recs: List[Dict[str, Any]],
    d: Dict[str, Any],
    plog: List[Dict[str, Any]],
    total: int,
) -> Dict[str, Any]:
    k = BR.FORMS[form]
    paths = read_paths(s, form)
    grid_j = sum(1 for r in plog if int(r["grid_evals"]) > 0)
    out: Dict[str, Any] = dict(
        sample=s,
        form=form,
        j_total_with_pass=total,
        steps=len(recs),
        stops_high=0,
        stops_low=0,
        stops_none=0,
        checks_went_on=0,
        steps_moved=0,
        delta_arith=0,
        delta_live=0,
        changed=0,
        broken=0,
        fixed=0,
        wrong_after=0,
        grid_j=grid_j,
        grid_moved_j=0,
        grid_moved_new=0,
        failed_moved_j=0,
        failed_moved_new=0,
        rescued_new=0,
        g2_accepted_new=0,
        lost_095=0,
        changed_wrong_j_right=[],
    )
    for r in recs:
        cuts = BR.refine_cuts(r["refinements"], k)
        for c in cuts.values():
            out["checks_went_on"] += int(c["w"]["continued"])
            out["stops_" + ("none" if c["w"]["kind"] == "none" else c["w"]["kind"])] += 1
        if not any(c["kind"] == "moved" for c in cuts.values()):
            out["delta_arith"] += BR.step_delta(cuts)
            continue
        out["steps_moved"] += 1
        p = paths.get(int(r["step"]))
        if p is None:
            raise SystemExit(f"{s} {form}: moved step {r['step']} has no path row (run band_paths)")
        out["delta_live"] += int(p["j_total"]) - int(p["new_total"])
        out["grid_moved_j"] += int(p["j_path"] == "grid")
        out["grid_moved_new"] += int(p["path"] == "grid")
        out["failed_moved_j"] += int(p["j_failed"])
        out["failed_moved_new"] += int(p["new_failed"])
        out["rescued_new"] += int(p["path"] == "rescued")
        out["g2_accepted_new"] += int(p["g2_accepted"])
        if int(p["j_stop"]) >= 0 and int(p["new_stop"]) != int(p["j_stop"]) and p["stopped_low"]:
            out["lost_095"] += 1
        if not p["taken_same"]:
            out["changed"] += 1
            jw, nw = float(p["err_j"]) > 1.0, float(p["err_new"]) > 1.0
            out["broken"] += int(nw and not jw)
            out["fixed"] += int(jw and not nw)
            out["wrong_after"] += int(nw)
            if nw and not jw:
                out["changed_wrong_j_right"].append(int(r["step"]))
    out["net_saving"] = int(out["delta_arith"] + out["delta_live"])
    out["net_saving_share"] = out["net_saving"] / total
    out["grid_new"] = out["grid_j"] - out["grid_moved_j"] + out["grid_moved_new"]
    return out


def judge(res: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """The gate of one form: ``res`` maps a sample to ``one``'s dict."""
    g = BR.GATE
    flags: Dict[str, bool] = {}
    flags["g1"] = all(res[s]["broken"] == 0 for s in SAMPLES)
    flags["g2"] = all(res[s]["grid_new"] == res[s]["grid_j"] for s in g["g2_grid_equal"]) and (
        res["held"]["grid_new"] <= res["held"]["grid_j"] + g["g2_grid_extra_held"]
    )
    flags["g3"] = all(res[s]["changed"] <= g["g3_changed_max"][s] for s in SAMPLES)
    flags["g4"] = all(res[s]["net_saving_share"] >= g["g4_saving_min"] for s in g["g4_samples"])
    flags["passes"] = all(flags.values())
    flags["summed_saving_mid_held"] = sum(res[s]["net_saving"] for s in g["g4_samples"])
    return flags


def choose(judged: Dict[str, Dict[str, Any]]) -> str:
    ok = [f for f in BR.FORMS if judged[f]["passes"]]
    if not ok:
        return ""
    if len(ok) == 1:
        return ok[0]
    a, b = (judged[f]["summed_saving_mid_held"] for f in ("PB1000", "PB500"))
    return "PB500" if b > a else "PB1000"  # a tie goes to PB1000


def tables(
    p1: Dict[str, Dict[str, Any]],
    res: Dict[str, Dict[str, Dict[str, Any]]],
    judged: Dict[str, Dict[str, Any]],
    brk: Dict[str, Dict[str, Dict[str, int]]],
) -> Dict[str, str]:
    t: Dict[str, str] = {}
    rows = [
        dict(
            sample=s,
            steps=p1[s]["steps_recorded"],
            refinements=p1[s]["refinements"],
            not_equal=p1[s]["steps_not_equal"],
            worst=p1[s]["worst_abs_diff"],
            ok=p1[s]["ok"],
        )
        for s in SAMPLES
    ]
    t["band_pred1"] = DT.markdown_table(
        rows,
        ["sample", "steps", "refinements", "not_equal", "worst", "ok"],
        {"worst": ".2e"},
    )
    for f in BR.FORMS:
        rr = []
        for s in SAMPLES:
            x = res[f][s]
            rr.append(
                dict(
                    sample=s,
                    stops=f"{x['stops_high']} / {x['stops_low']} / {x['stops_none']}",
                    went_on=x["checks_went_on"],
                    saving=x["net_saving"],
                    share=100.0 * x["net_saving_share"],
                    moved=x["steps_moved"],
                    changed=x["changed"],
                    broken=x["broken"],
                    fixed=x["fixed"],
                    grid=f"{x['grid_j']} -> {x['grid_new']}",
                    failed=f"{x['failed_moved_j']} -> {x['failed_moved_new']}",
                    lost=x["lost_095"],
                )
            )
        t[f"band_gate_{f}"] = DT.markdown_table(
            rr,
            [
                "sample",
                "stops",
                "went_on",
                "saving",
                "share",
                "moved",
                "changed",
                "broken",
                "fixed",
                "grid",
                "failed",
                "lost",
            ],
            {"saving": ",d", "went_on": ",d", "share": ".2f"},
            labels=list(SAMPLES),
        )
        j = judged[f]
        t[f"band_verdict_{f}"] = DT.markdown_table(
            [dict(form=f, g1=j["g1"], g2=j["g2"], g3=j["g3"], g4=j["g4"], passes=j["passes"])],
            ["form", "g1", "g2", "g3", "g4", "passes"],
        )
    rb = []
    for kk in BR.RECORD_KS:
        pooled = dict(n=0, low=0, mid=0, high=0)
        for s in SAMPLES:
            for key in pooled:
                pooled[key] += brk[str(kk)][s][key]
        rb.append(
            dict(
                k=kk,
                n=pooled["n"],
                low=pooled["low"],
                share=100.0 * pooled["low"] / max(1, pooled["n"]),
                mid=pooled["mid"],
                high=pooled["high"],
            )
        )
    t["band_bracketed"] = DT.markdown_table(
        rb, ["k", "n", "low", "share", "mid", "high"], {"share": ".1f"}
    )
    return t


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-choose", action="store_true")
    a = ap.parse_args()
    p1: Dict[str, Dict[str, Any]] = {}
    res: Dict[str, Dict[str, Dict[str, Any]]] = {f: {} for f in BR.FORMS}
    brk: Dict[str, Dict[str, Dict[str, int]]] = {str(k): {} for k in BR.RECORD_KS}
    for s in SAMPLES:
        d = BREC.load(s)
        recs = BPA.read_records(s)
        p1[s] = prediction1(s, recs, d)
        if not p1[s]["ok"]:
            (BC.BD / "gate.json").write_text(json.dumps(dict(prediction1=p1), indent=1) + "\n")
            raise SystemExit(f"STOP: prediction 1 fails on {s}: {p1[s]}")
        plog = BC.load_gz(BREC.stored_dir(s) / "polish_log.json.gz")
        total = int(d["stats"]["arms"]["J"]["stages"]["total"])
        for f in BR.FORMS:
            res[f][s] = one(s, f, recs, d, plog, total)
        for kk in BR.RECORD_KS:
            brk[str(kk)][s] = bracketed(recs, kk)
    judged = {f: judge(res[f]) for f in BR.FORMS}
    chosen = "" if a.no_choose else choose(judged)
    out = dict(
        prediction1=p1,
        forms=res,
        verdict=judged,
        bracketed=brk,
        chosen=chosen,
        rule=BR.RULE,
        gate=BR.GATE,
    )
    (BC.BD / "gate.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    DT.write_tables(BC.BD / "gate_tables.md", tables(p1, res, judged, brk))
    if chosen:
        import hashlib

        sha = hashlib.sha256((BC.BD / "gate.json").read_bytes()).hexdigest()
        (BC.BD / "chosen.json").write_text(
            json.dumps(
                dict(
                    chosen=chosen,
                    k=BR.FORMS[chosen],
                    gate_json_sha256=sha,
                    rule="the plan's fixed choice rule (band_rules.GATE['choice'])",
                ),
                indent=1,
            )
            + "\n"
        )
    for f in BR.FORMS:
        print(f, {k: v for k, v in judged[f].items()})
    print("chosen:", chosen or "NONE (stop)")


if __name__ == "__main__":
    main()
