"""C++ vs Python RandomRestartZeroTemp (the MC of the quick stage and of FindOptimal) on the
Phase D full-sample clean images.

C++ side: a scratch build of Src/ with debug prints in COrientationMC::RandomRestartZeroTemp
(patch: benchmarks/phase_d_seed_diag/mc_parity/cpp_debug_patch_mc.py.txt; not committed). It prints
one "MC ..." line per call (steps, blocks, restarts, accepts, evaluations, init/final cost, stop
reason) and, with MBTRACE=1, one "MB ..." line per block and one "MR ..." line per restart draw.
Python side: runs/<tag>_full_v<V>.json of seed_diag.py (calls with blocks, steps, restarts, ...).

  uv run python scripts/phase_d/mc_parity.py summary --cpp-dir DIR --tag mcport

DIR holds cppq_<voxel>_<k>.log (3 plain runs per voxel) and cpp_15901.log (MBTRACE=1).
Writes benchmarks/phase_d_seed_diag/mc_parity/summary.json, tables.md and the trimmed trace.
"""

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
OUT = ICE / "benchmarks" / "phase_d_seed_diag"
PAR = OUT / "mc_parity"
sys.path.insert(0, str(ICE / "scripts" / "common"))
import doc_tables as DT  # noqa: E402

sys.path.insert(0, str(HERE))
import seed_diag_summary as SDS  # noqa: E402

VOXELS = [15901, 17034, 18558, 19889, 22867, 5242]
KV = re.compile(r"(\w+) (\S+)")


def parse_cpp(log: Path) -> Dict[str, Any]:
    """MC calls, MB blocks and MR restarts of the adaptive reconstruction (after the line
    'Normal Reconstruction', i.e. not the basic reconstructor's MC calls)."""
    text = log.read_text()
    calls: List[Dict[str, float]] = []
    blocks: List[List[Dict[str, float]]] = []
    cur: List[Dict[str, float]] = []
    adaptive = False
    ev = sec = None
    for line in log.read_text().splitlines():
        if line.startswith("Normal Reconstruction"):
            adaptive = True
            cur = []
            continue
        m = re.search(r"adap_sec=([\d.]+) adap_evals=(\d+)", line)
        if m:
            sec, ev = float(m.group(1)), int(m.group(2))
        if not adaptive:
            continue
        if line.startswith("MB "):
            t = line.split()
            cur.append(dict(total=int(t[2]), nopt=int(t[4]), step_deg=float(t[6]), fail=int(t[12])))
        elif line.startswith("MR "):
            t = line.split()
            cur[-1]["restart_radius"] = float(t[2])
            cur[-1]["restart_xyz"] = [float(t[4]), float(t[5]), float(t[6])]
        elif line.startswith("MC "):
            d = {k: float(v) for k, v in KV.findall(line[3:])}
            calls.append(d)
            blocks.append(cur)
            cur = []
    cand = [
        int(m.group(1))
        for m in re.finditer(r"Num candiates (\d+)", text.split("Normal Reconstruction")[-1])
    ]
    return dict(calls=calls, blocks=blocks, adap_evals=ev, adap_sec=sec, level_candidates=cand)


def block_rules(blocks: List[Dict[str, float]], nmin: int, step0_deg: float) -> Dict[str, int]:
    """Check one call's blocks against the C++ rules; returns violation counts."""
    v = dict(len=0, halve=0, reset=0)
    for i, b in enumerate(blocks):
        last = i == len(blocks) - 1
        if not last and b["nopt"] != nmin:
            v["len"] += 1
        if i == 0:
            if abs(b["step_deg"] - step0_deg) > 1e-4 * step0_deg:
                v["reset"] += 1
        else:
            p = blocks[i - 1]
            want = 0.5 * p["step_deg"] if not p["fail"] else step0_deg
            if abs(b["step_deg"] - want) > 1e-4 * want:
                v["halve" if not p["fail"] else "reset"] += 1
    return v


def py_block_rules(blk: List[List[Any]], nmin: int, step0_deg: float) -> Dict[str, int]:
    bl = [dict(nopt=b[0], step_deg=b[1], fail=int(b[2] == "f")) for b in blk]
    return block_rules(bl, nmin, step0_deg)


def agg(calls: List[Dict[str, float]], key_steps: int) -> Dict[str, float]:
    sel = [c for c in calls if int(c["maxsteps"]) == key_steps]
    if not sel:
        return {}
    a = lambda k: float(np.mean([c[k] for c in sel]))  # noqa: E731
    return dict(
        n=len(sel), steps=a("steps"), blocks=a("blocks"), restarts=a("restarts"),
        accepts=a("accepts"), evals=a("evals"),
        frac_stop_restarts=float(np.mean([c["stop"] == 1 for c in sel])),
        frac_stop_conv=float(np.mean([c["stop"] == 2 for c in sel])),
        frac_budget=float(np.mean([c["stop"] == 0 for c in sel])),
        frac_any_accept=float(np.mean([c["accepts"] > 0 for c in sel])),
    )  # fmt: skip


CALL_FORMATS = {
    **{k: ".2f" for k in ("q_ev_cpp", "q_ev_py", "f_ev_cpp", "f_ev_py", "f_blk_cpp", "f_blk_py")},
    **{k: ".2f" for k in ("f_stop1_cpp", "f_stop1_py", "f_acc_cpp", "f_acc_py")},
    **{k: ".2f" for k in ("q_blk_cpp", "q_blk_py")},
    **{"n_quick_cpp": ".0f", "n_quick_py": "d", "n_find_cpp": "d", "n_find_py": "d"},
}
QUICK_KEYS = ["evals", "blocks", "steps", "restarts", "accepts", "stop"]


def load_py(tag: str, v: int) -> Dict[str, Any]:
    """runs/<tag>_full_v<V>.json; the quick-MC calls may be stored column-wise (compact)."""
    d = json.load(open(OUT / "runs" / f"{tag}_full_v{v}.json"))
    if "quick_cols" in d:
        cols = d.pop("quick_cols")
        q = [
            dict(zip(QUICK_KEYS, row), phase="quick", quick=True)
            for row in zip(*[cols[k] for k in QUICK_KEYS])
        ]
        d["calls"] = q + d["calls"]
    return d


def compact(tag: str) -> None:
    """Rewrite runs/<tag>_full_v*.json with the quick-MC calls as integer columns (the per-call
    seconds, costs and hit ratios of ~15k calls are not needed here)."""
    for v in VOXELS:
        f = OUT / "runs" / f"{tag}_full_v{v}.json"
        d = json.load(open(f))
        if "quick_cols" in d:
            continue
        quick = [c for c in d["calls"] if c["quick"]]
        d["quick_cols"] = {k: [int(c[k]) for c in quick] for k in QUICK_KEYS}
        d["calls"] = [c for c in d["calls"] if not c["quick"]]
        f.write_text(json.dumps(d, default=float, separators=(",", ":")) + "\n")


def py_agg(calls: List[Dict[str, Any]], quick: bool) -> Dict[str, float]:
    sel = [c for c in calls if c["phase"] != "variance" and bool(c["quick"]) == quick]
    if not sel:
        return {}
    a = lambda k: float(np.mean([c[k] for c in sel]))  # noqa: E731
    return dict(
        n=len(sel), steps=a("steps"), blocks=a("blocks"), restarts=a("restarts"),
        accepts=a("accepts"), evals=a("evals"),
        frac_stop_restarts=float(np.mean([c["stop"] == 1 for c in sel])),
        frac_stop_conv=float(np.mean([c["stop"] == 2 for c in sel])),
        frac_budget=float(np.mean([c["stop"] == 0 for c in sel])),
        frac_any_accept=float(np.mean([c["accepts"] > 0 for c in sel])),
    )  # fmt: skip


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["summary", "compact"])
    ap.add_argument("--cpp-dir", default="")
    ap.add_argument("--tag", default="mcport")
    a = ap.parse_args()
    if a.cmd == "compact":
        compact(a.tag)
        return
    cdir = Path(a.cpp_dir)
    PAR.mkdir(parents=True, exist_ok=True)
    trace = parse_cpp(cdir / "cpp_15901.log")
    res: Dict[str, Any] = {"voxels": {}}
    rows: List[Dict[str, Any]] = []
    rule_rows: List[Dict[str, Any]] = []
    tot = dict(cpp_quick=[], py_quick=[], cpp_find=[], py_find=[])  # per-call evals
    for v in VOXELS:
        logs = sorted(cdir.glob(f"cppq_{v}_*.log"))
        cps = [parse_cpp(p) for p in logs]
        py = load_py(a.tag, v)
        calls_c = [c for p in cps for c in p["calls"]]
        d: Dict[str, Any] = dict(
            cpp_runs=len(cps),
            cpp_adap_evals=[p["adap_evals"] for p in cps],
            cpp_adap_sec=[p["adap_sec"] for p in cps],
            py_evals_total=py["evals_total"],
            py_quick=py_agg(py["calls"], True), py_find=py_agg(py["calls"], False),
            cpp_quick=agg(calls_c, 10), cpp_find=agg(calls_c, 200),
        )  # fmt: skip
        # evaluations by stage (Python) vs the C++ run's totals
        d["py_stage_evals"] = {
            k: py[k] for k in ("quick_evals", "find_evals", "variance_evals") if k in py
        }
        d["cpp_quick_evals_per_run"] = float(
            np.mean([sum(c["evals"] for c in p["calls"] if c["maxsteps"] == 10) for p in cps])
        )
        d["cpp_find_evals_per_run"] = float(
            np.mean([sum(c["evals"] for c in p["calls"] if c["maxsteps"] == 200) for p in cps])
        )
        d["cpp_n_quick_per_run"] = float(
            np.mean([sum(1 for c in p["calls"] if c["maxsteps"] == 10) for p in cps])
        )
        d["py_n_quick"] = py["n_quick_mc"]
        # Python FindOptimal / quick evaluation vs C++ totals
        res["voxels"][v] = d
        rows.append(
            dict(voxel=v, n_quick_cpp=d["cpp_n_quick_per_run"], n_quick_py=d["py_n_quick"],
                 n_find_cpp=d["cpp_find"].get("n"), n_find_py=d["py_find"].get("n"),
                 q_blk_cpp=d["cpp_quick"].get("blocks"), q_blk_py=d["py_quick"].get("blocks"),
                 q_ev_cpp=d["cpp_quick"].get("evals"), q_ev_py=d["py_quick"].get("evals"),
                 f_ev_cpp=d["cpp_find"].get("evals"), f_ev_py=d["py_find"].get("evals"),
                 f_blk_cpp=d["cpp_find"].get("blocks"), f_blk_py=d["py_find"].get("blocks"),
                 f_stop1_cpp=d["cpp_find"].get("frac_stop_restarts"),
                 f_stop1_py=d["py_find"].get("frac_stop_restarts"),
                 f_acc_cpp=d["cpp_find"].get("frac_any_accept"),
                 f_acc_py=d["py_find"].get("frac_any_accept"))
        )  # fmt: skip
        # block-rule violations (Python: every FindOptimal call's blocks)
        nmin = 31
        step0 = 0.131687
        vio = dict(len=0, halve=0, reset=0)
        n_calls = 0
        for c in py["calls"]:
            if c.get("blk"):
                n_calls += 1
                for k, x in py_block_rules(c["blk"], nmin, step0).items():
                    vio[k] += x
        cvio = dict(len=0, halve=0, reset=0)
        cn = 0
        for p in [trace] if v == 15901 else []:  # block lines only in the MBTRACE run
            for call, bl in zip(p["calls"], p["blocks"]):
                if int(call["maxsteps"]) == 200 and bl:
                    cn += 1
                    for k, x in block_rules(bl, nmin, step0).items():
                        cvio[k] += x
        rule_rows.append(dict(voxel=v, py_calls=n_calls, py_viol=sum(vio.values()),
                              cpp_calls=cn, cpp_viol=sum(cvio.values())))  # fmt: skip
        d["rule_violations"] = dict(py=vio, cpp=cvio)
    # stage split: Python by phase; C++ quick and FindOptimal from the MC lines, the rest = discrete
    # stages + variance stage + final evaluation
    st_rows = []
    for v in VOXELS:
        py = load_py(a.tag, v)
        sp = SDS.stage_evals(py)
        cps = [parse_cpp(p) for p in sorted(cdir.glob(f"cppq_{v}_*.log"))]
        cq = float(
            np.mean([sum(c["evals"] for c in p["calls"] if c["maxsteps"] == 10) for p in cps])
        )
        cf = float(
            np.mean([sum(c["evals"] for c in p["calls"] if c["maxsteps"] == 200) for p in cps])
        )
        ct = float(np.mean([p["adap_evals"] for p in cps]))
        lev_c = np.mean(
            [p["level_candidates"][:4] for p in cps if len(p["level_candidates"]) >= 4], axis=0
        )
        lev_p = [lv["n_returned"] for lv in py["levels"]][:4]
        st_rows.append(dict(
            voxel=v, quick_py=sp["quick_L0"] + sp["quick_L1_3"], quick_cpp=cq,
            find_py=sp["find"], find_cpp=cf,
            rest_py=sp["disc_L0"] + sp["disc_L1_3"] + sp["variance"] + 1,
            rest_cpp=ct - cq - cf, total_py=py["evals_total"], total_cpp=ct,
            cand_L0_py=lev_p[0], cand_L0_cpp=float(lev_c[0]),
            cand_all_py=sum(lev_p), cand_all_cpp=float(lev_c.sum()),
        ))  # fmt: skip
    res["stage_split"] = st_rows
    # seed cost: after the MC port (mcport) against the VM-fix-only run (varfix) and C++
    seeds = [350, 395, 504]  # run_estimates.json (oracle fits, p_fail 0 / 0.1 / 0.3)
    sc_rows = []
    for v in VOXELS:
        old = json.load(open(OUT / "runs" / f"varfix_full_v{v}.json"))
        new = load_py(a.tag, v)
        cps = [parse_cpp(p) for p in sorted(cdir.glob(f"cppq_{v}_*.log"))]
        sc_rows.append(dict(
            voxel=v, evals_before=old["evals_total"], evals_after=new["evals_total"],
            evals_cpp=float(np.mean([p["adap_evals"] for p in cps])),
            wall_before=old["wall_s"], wall_after=new["wall_s"],
            wall_cpp=float(np.mean([p["adap_sec"] for p in cps])),
            err_before=old["err_deg"], err_after=new["err_deg"], quiet_after=new["quiet_preflight"],
        ))  # fmt: skip
    mean = lambda k: float(np.mean([r[k] for r in sc_rows]))  # noqa: E731
    sc_rows.append(
        dict(voxel="mean", **{k: mean(k) for k in sc_rows[0] if k not in ("voxel", "quiet_after")})
    )
    res["seed_cost"] = sc_rows
    wall = mean("wall_after")
    res["serial_seed_pass_hours"] = dict(
        per_seed_s=wall, seeds=seeds, hours=[n * wall / 3600.0 for n in seeds],
        label="single-process, quiet (preflight passed in 6/6)",
    )  # fmt: skip
    # whole-run evaluations
    ev_rows = []
    for v in VOXELS:
        d = res["voxels"][v]
        c = float(np.mean(d["cpp_adap_evals"]))
        ev_rows.append(dict(voxel=v, cpp_evals=c, py_evals=d["py_evals_total"],
                            ratio=d["py_evals_total"] / c))  # fmt: skip
    res["evals"] = ev_rows
    # restart draws of the traced run (voxel 15901): base and range are pinned by the unit test;
    # here the printed radius is checked against tan(box)/sqrt(48)
    tr = trace
    rad = [b["restart_radius"] for bl in tr["blocks"] for b in bl if "restart_radius" in b]
    res["trace_15901"] = dict(
        n_calls=len(tr["calls"]), n_restart_draws=len(rad),
        radius_min=float(min(rad)), radius_max=float(max(rad)),
        radius_expected_find=float(np.tan(np.radians(0.329218)) / np.sqrt(48.0)),
        radius_expected_quick=float(np.tan(np.radians(5.0 / 3.0)) / np.sqrt(48.0)),
    )  # fmt: skip
    # trimmed trace: the first 3 FindOptimal calls and the first 5 quick calls
    lines = (cdir / "cpp_15901.log").read_text().splitlines()
    out_lines, adaptive, nq, nf, buf = [], False, 0, 0, []
    for ln in lines:
        if ln.startswith("Normal Reconstruction"):
            adaptive = True
            continue
        if not adaptive or not (ln.startswith(("MB ", "MR ", "MC "))):
            continue
        buf.append(ln)
        if ln.startswith("MC "):
            q = "maxsteps 10 " in ln
            if (q and nq < 5) or (not q and nf < 3):
                out_lines += buf
                nq += q
                nf += not q
            buf = []
    (PAR / "cpp_trace_15901_trimmed.log").write_text("\n".join(out_lines) + "\n")
    (PAR / "summary.json").write_text(json.dumps(res, indent=1, default=float))
    tables = {
        "mc_parity_calls": DT.markdown_table(
            rows,
            ["voxel", "n_quick_cpp", "n_quick_py", "q_blk_cpp", "q_blk_py", "n_find_cpp",
             "n_find_py", "q_ev_cpp", "q_ev_py", "f_ev_cpp", "f_ev_py",
             "f_blk_cpp", "f_blk_py", "f_stop1_cpp", "f_stop1_py", "f_acc_cpp", "f_acc_py"],
            formats=CALL_FORMATS,
        ),
        "mc_parity_evals": DT.markdown_table(
            ev_rows, ["voxel", "cpp_evals", "py_evals", "ratio"],
            formats={"cpp_evals": ",.0f", "py_evals": ",d", "ratio": ".2f"},
        ),
        "mc_parity_stages": DT.markdown_table(
            st_rows,
            ["voxel", "quick_py", "quick_cpp", "find_py", "find_cpp", "rest_py", "rest_cpp",
             "total_py", "total_cpp", "cand_L0_py", "cand_L0_cpp", "cand_all_py", "cand_all_cpp"],
            formats={k: ",.0f" for k in ("quick_py", "quick_cpp", "find_py", "find_cpp", "rest_py",
                     "rest_cpp", "total_py", "total_cpp", "cand_L0_py", "cand_L0_cpp",
                     "cand_all_py", "cand_all_cpp")},
        ),
        "mc_parity_rules": DT.markdown_table(
            rule_rows, ["voxel", "py_calls", "py_viol", "cpp_calls", "cpp_viol"]
        ),
    }  # fmt: skip
    tables["mc_parity_seed_cost"] = DT.markdown_table(
        sc_rows,
        ["voxel", "evals_before", "evals_after", "evals_cpp", "wall_before", "wall_after",
         "wall_cpp", "err_before", "err_after"],
        formats={"evals_before": ",.0f", "evals_after": ",.0f", "evals_cpp": ",.0f",
                 "wall_before": ".0f", "wall_after": ".0f", "wall_cpp": ".1f",
                 "err_before": ".3f", "err_after": ".3f"},
    )  # fmt: skip
    hrs = res["serial_seed_pass_hours"]
    tables["mc_parity_estimate"] = DT.markdown_table(
        [dict(basis="measured wall, mean of 6 (quiet)", per_seed_s=hrs["per_seed_s"],
              hours_low=hrs["hours"][0], hours_mid=hrs["hours"][1], hours_high=hrs["hours"][2])],
        ["basis", "per_seed_s", "hours_low", "hours_mid", "hours_high"],
        formats={"per_seed_s": ".0f", "hours_low": ".1f", "hours_mid": ".1f", "hours_high": ".1f"},
    )  # fmt: skip
    DT.write_tables(PAR / "tables.md", tables)
    print(json.dumps(res["evals"], indent=1, default=float))


if __name__ == "__main__":
    main()
