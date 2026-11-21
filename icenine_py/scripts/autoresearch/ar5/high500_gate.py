"""AR5 step 16, phase A step 4: the gate of rule ``PH(500)`` on the four disclosed samples, and the
predictions 1 and 2.

  uv run python scripts/autoresearch/ar5/high500_gate.py

Reads (read only) step 15's recording rows and arm ``J``'s stored run (linked into the cache by
``high500_link_cache.sh``), the live replay rows (``high500_paths.py``) and step 14's tracked
``stats_<s>.json``. Needs ``ar5/high500/frozen.json`` committed and agreeing (``high500_common``,
sample ``held``). Writes ``ar5/high500/gate.json`` and ``gate_tables.md`` (blocks for
``sync_doc_tables.py``).

Per sample: stored counts; per refinement the walk of the recorded 500-stall events with the high
stop only; high stops, checks that went on, steps with a moved end, steps whose path can change (a
moved end whose own ``J`` end is under 0.95), the net saving of the arithmetic walk of every step
(the plan's table: every check counted, growth after a changed voxel not replayed), and the net
saving with the live replay of every moved step in place of its walk. The decision rules are
``high500_rules.GATE`` (frozen). Exit 1 when the gate fails or a stored count does not reproduce.
"""

import json
import re
import sys
from typing import Any, Dict, List

import high500_common as BC  # noqa: E402  (first)

import high500_rules as HR  # noqa: E402
import high500_stored as HS  # noqa: E402

sys.path.insert(0, str(BC.ICE / "scripts" / "common"))
import doc_tables as DT  # noqa: E402

SAMPLES = HR.GATE["samples"]
TOL = 1e-12


def read_paths(s: str) -> Dict[int, Dict[str, Any]]:
    rx = re.compile(rf"^paths_{s}_s\d+of\d+\.jsonl$")
    rows: Dict[int, Dict[str, Any]] = {}
    for p in sorted((BC.HCACHE / "paths").glob("paths_*.jsonl")):
        if rx.match(p.name):
            for line in BC.HG.guarded(p).read_text().splitlines():
                if line.strip():
                    r = json.loads(line)
                    rows[int(r["step"])] = r
    return rows


def prediction1(s: str, recs: List[Dict[str, Any]], d: Dict[str, Any]) -> Dict[str, Any]:
    c = HS.counts(s, d)
    n_steps = len(HS.mc_steps(d))
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
        refinements=sum(len(r["refinements"]) for r in recs),
        steps_not_equal=int(sum(1 for r in recs if not r["p1_ok"])),
        worst_abs_diff=float(max([r["worst_d"] for r in recs] + [0.0])),
        ok=ok,
    )


def one(
    s: str, recs: List[Dict[str, Any]], d: Dict[str, Any], plog: List[Dict[str, Any]], total: int
) -> Dict[str, Any]:
    paths = read_paths(s)
    out: Dict[str, Any] = dict(
        sample=s,
        j_total_with_pass=total,
        steps=len(recs),
        refinements=0,
        stops_high=0,
        stops_none=0,
        checks_went_on=0,
        moved_steps=[],
        path_change_steps=[],
        saving_walk=0,
        delta_arith=0,
        delta_live=0,
        live_rows=0,
        live_taken_changed=[],
        live_path_changed=[],
        live_broken=[],
        live_fixed=[],
        live_failed_j=0,
        live_failed_new=0,
        live_grid_j=0,
        live_grid_new=0,
        grid_j=sum(1 for r in plog if int(r["grid_evals"]) > 0),
    )
    for r in recs:
        cuts = HR.refine_cuts(r["refinements"])
        for c in cuts.values():
            out["refinements"] += 1
            out["checks_went_on"] += int(c["w"]["continued"])
            out["stops_high" if c["w"]["kind"] == "high" else "stops_none"] += 1
        walk = sum(int(c["base"]) - int(c["evals"]) for c in cuts.values())
        out["saving_walk"] += walk
        moved = any(c["kind"] == "moved" for c in cuts.values())
        if HR.path_change_refs(r["refinements"]):
            out["path_change_steps"].append(int(r["step"]))
        if not moved:
            out["delta_arith"] += walk
            continue
        step = int(r["step"])
        out["moved_steps"].append(step)
        p = paths.get(step)
        if p is None:
            raise SystemExit(f"{s}: moved step {step} has no live row (run high500_paths)")
        out["live_rows"] += 1
        out["delta_live"] += int(p["j_total"]) - int(p["new_total"])
        out["live_failed_j"] += int(p["j_failed"])
        out["live_failed_new"] += int(p["new_failed"])
        out["live_grid_j"] += int(p["j_path"] == "grid")
        out["live_grid_new"] += int(p["path"] == "grid")
        if p["path"] != p["j_path"]:
            out["live_path_changed"].append(step)
        if not p["taken_same"]:
            out["live_taken_changed"].append(step)
            jw, nw = float(p["err_j"]) > 1.0, float(p["err_new"]) > 1.0
            if nw and not jw:
                out["live_broken"].append(step)
            if jw and not nw:
                out["live_fixed"].append(step)
    out["net_saving"] = int(out["delta_arith"] + out["delta_live"])
    out["net_saving_share"] = out["net_saving"] / total
    out["saving_walk_share"] = out["saving_walk"] / total
    out["grid_new"] = out["grid_j"] - out["live_grid_j"] + out["live_grid_new"]
    return out


def judge(res: Dict[str, Dict[str, Any]], p1: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    g = HR.GATE
    flags: Dict[str, Any] = {}
    exp = HR.PRED1
    flags["g0"] = all(
        p1[s]["ok"]
        and len(res[s]["moved_steps"]) == exp[s]["moved"]
        and len(res[s]["path_change_steps"]) == exp[s]["path_change"]
        and res[s]["saving_walk"] == exp[s]["saving"]
        for s in SAMPLES
    )
    flags["g1"] = all(len(res[s]["moved_steps"]) <= g["g1_moved_max"][s] for s in SAMPLES) and all(
        len(res[s]["path_change_steps"]) <= g["g1_path_change_max"] for s in SAMPLES
    )
    flags["g2"] = all(
        not res[s]["live_broken"]
        and res[s]["live_failed_new"] <= res[s]["live_failed_j"]
        and res[s]["live_grid_new"] <= res[s]["live_grid_j"]
        for s in SAMPLES
    )
    flags["g3"] = all(res[s]["net_saving_share"] >= g["g3_saving_min"] for s in g["g3_samples"])
    flags["passes"] = bool(all(flags[k] for k in ("g0", "g1", "g2", "g3")))
    return flags


def tables(
    p1: Dict[str, Dict[str, Any]], res: Dict[str, Dict[str, Any]], judged: Dict[str, Any]
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
    t["high500_pred1"] = DT.markdown_table(
        rows, ["sample", "steps", "refinements", "not_equal", "worst", "ok"], {"worst": ".2e"}
    )
    rr = []
    for s in SAMPLES:
        x = res[s]
        rr.append(
            dict(
                sample=s,
                stops=f"{x['stops_high']} / {x['stops_none']}",
                went_on=x["checks_went_on"],
                moved=len(x["moved_steps"]),
                path=len(x["path_change_steps"]),
                walk=x["saving_walk"],
                walk_share=100.0 * x["saving_walk_share"],
                net=x["net_saving"],
                net_share=100.0 * x["net_saving_share"],
            )
        )
    t["high500_gate"] = DT.markdown_table(
        rr,
        ["sample", "stops", "went_on", "moved", "path", "walk", "walk_share", "net", "net_share"],
        {"went_on": ",d", "walk": ",d", "net": ",d", "walk_share": ".2f", "net_share": ".2f"},
        labels=list(SAMPLES),
    )
    t["high500_verdict"] = DT.markdown_table(
        [
            dict(
                g0=judged["g0"],
                g1=judged["g1"],
                g2=judged["g2"],
                g3=judged["g3"],
                passes=judged["passes"],
            )
        ],
        ["g0", "g1", "g2", "g3", "passes"],
    )
    lr = []
    for s in SAMPLES:
        x = res[s]
        lr.append(
            dict(
                sample=s,
                live=x["live_rows"],
                taken_changed=len(x["live_taken_changed"]),
                path_changed=len(x["live_path_changed"]),
                broken=len(x["live_broken"]),
                fixed=len(x["live_fixed"]),
                failed=f"{x['live_failed_j']} -> {x['live_failed_new']}",
                grid=f"{x['live_grid_j']} -> {x['live_grid_new']}",
            )
        )
    t["high500_live"] = DT.markdown_table(
        lr,
        ["sample", "live", "taken_changed", "path_changed", "broken", "fixed", "failed", "grid"],
        labels=list(SAMPLES),
    )
    return t


def main() -> None:
    p1: Dict[str, Dict[str, Any]] = {}
    res: Dict[str, Dict[str, Any]] = {}
    for s in SAMPLES:
        d = HS.load(s)
        recs = HS.read_records(s)
        p1[s] = prediction1(s, recs, d)
        if not p1[s]["ok"]:
            (BC.HD / "gate.json").write_text(json.dumps(dict(prediction1=p1), indent=1) + "\n")
            raise SystemExit(f"STOP: prediction 1 fails on {s}: {p1[s]}")
        plog = BC.load_gz(HS.stored_dir(s) / "polish_log.json.gz")
        total = int(d["stats"]["arms"]["J"]["stages"]["total"])
        res[s] = one(s, recs, d, plog, total)
    judged = judge(res, p1)
    out = dict(
        prediction1=p1, samples=res, verdict=judged, rule=HR.RULE, gate=HR.GATE, pred1=HR.PRED1
    )
    (BC.HD / "gate.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    DT.write_tables(BC.HD / "gate_tables.md", tables(p1, res, judged))
    print(json.dumps(judged))
    for s in SAMPLES:
        x = res[s]
        print(
            s,
            dict(
                moved=x["moved_steps"],
                path_change=x["path_change_steps"],
                walk=x["saving_walk"],
                net=x["net_saving"],
                share=round(100 * x["net_saving_share"], 3),
            ),
        )
    sys.exit(0 if judged["passes"] else 1)


if __name__ == "__main__":
    main()
