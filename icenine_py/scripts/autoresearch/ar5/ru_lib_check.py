"""AR5 step 29, phase B: is the key-on library run equal to the script arm ``RUl`` (L1), and what
does the library form save against ``LibPXF`` (L2)?

  RESEED_SAMPLE=mid RU_CACHE=<cache> uv run python scripts/autoresearch/ar5/ru_lib_check.py

L1 (equality, before the boundary pass, ``mid``): ``LibRU`` (``LibPXF``'s keys plus
``SeedReuseKnown 1``) against ``RUl`` (step 28's ``PXFl`` plus the ``RU`` wrapper; ``ru_hl.py``):
the map (decompressed ``.mic`` bytes), the evaluation total, every counter of the script run, all
per-voxel records, the rescues, the MC refinements with their checks, the pauses and resumes, the
PXF steps (the steps a try did not end), and every try row (voxel, tried, starts, costs, each fit's
start, evaluations, hit ratio and cost, success, second). Seconds are not compared.

L2 (library form): ``LibRU`` before the pass at least 10% fewer evaluations than ``LibPXF``
(6,731,510); after the pass, against ``LibPXF``'s map: broken minus fixed at most 3 and the median
error at most 0.001 degree higher (the truth is read only to score). Also reported (no gate): the
try share, second fits, evaluations per try, the lowest accepted hit ratio, and the prefix of rows
that equal phase A's stored ``ru_steps_mid``.

The criteria are those of the plan section "AR5 step 29 phase B" and are never changed after a run.
Writes ``ar5/ru/lib_check.json`` and ``ar5/ru/lib_tables.md``. Does not go through ``ru_common``
(the library changed, so that freeze cannot agree).
"""

import json
import os
import sys
from pathlib import Path
from typing import Any, Callable, Dict, List

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
RD = ICE / "benchmarks" / "autoresearch" / "ar5" / "ru"


def _default_cache() -> Path:
    parts = ICE.parts
    if ".claude" in parts and "worktrees" in parts:
        return Path(*parts[: parts.index(".claude")]) / (
            "icenine_py/benchmarks/autoresearch/ar5/ru/cache"
        )
    return RD / "cache"


SCACHE = Path(os.environ.get("RU_CACHE", str(_default_cache())))
os.environ["HIGH500_CACHE"] = str(SCACHE)
os.environ["RESCUE_CACHE"] = str(SCACHE)
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICE / "scripts" / "common"))

import rescue_lib_check as RLC  # noqa: E402  (first: rescue_common sets the cache)

import gradient_common as G  # noqa: E402
import high500_lib_check as H5C  # noqa: E402
import numpy as np  # noqa: E402
import pxf_lib_check as PLC  # noqa: E402
import reseed_common as RC  # noqa: E402

LIBPXF_EVALS = 6_731_510  # the stored LibPXF run (step 28), before its boundary pass
OLD_CACHE_PXF = SCACHE.parent.parent / "pxf" / "cache"


def _r(x: Any) -> Any:
    return None if x is None else round(float(x), 12)


def row_key_script(r: Dict[str, Any]) -> Any:
    return (
        int(r["voxel"]),
        bool(r["tried"]),
        [int(j) for j in r["starts"]],
        [_r(c) for c in r["costs"]],
        [
            (f["start"], f["evals"], _r(f["hit"]), _r(f["cost"]), bool(f["second"]))
            for f in r["fits"]
        ],
        int(r["evals"]),
        bool(r["success"]),
        bool(r["second"]),
    )


def lib_try_rows(index_log: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    out = []
    for e in index_log:
        r = dict(e["reuse"])
        r["voxel"] = e["voxel"]
        out.append(r)
    return out


def compare_tries(hl: Path, index_log: List[Dict[str, Any]]) -> Dict[str, Any]:
    srows = RLC.read_gz(hl / "ru_rows.json.gz")["rows"]
    lrows = lib_try_rows(index_log)
    s_all = [row_key_script(r) for r in srows]
    l_all = [row_key_script(r) for r in lrows]
    first = next((i for i, (a, b) in enumerate(zip(s_all, l_all)) if a != b), None)
    if first is None and len(s_all) != len(l_all):
        first = min(len(s_all), len(l_all))
    return dict(
        rows=dict(script=len(s_all), library=len(l_all)),
        tried=dict(
            script=sum(1 for r in srows if r["tried"]), library=sum(1 for r in lrows if r["tried"])
        ),
        success=dict(
            script=sum(1 for r in srows if r["success"]),
            library=sum(1 for r in lrows if r["success"]),
        ),
        first_differ=first,
        first_differ_voxel=(
            None if first is None else (int(srows[first]["voxel"]) if first < len(srows) else None)
        ),
        equal=bool(s_all == l_all),
    )


def compare_pxf_steps(hl: Path, index_log: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Step 28's fallback comparison on the seed steps a try did not end; the library's entry
    evaluations include the try's (added after the PXF step), the script's step rows do not."""
    srows = RLC.read_gz(hl / "pxf_steps.json.gz")["rows"]

    def s_key(r: Dict[str, Any]) -> Any:
        f, s = r["first"], r.get("second")
        return (
            bool(r["fallback"]), f["evals"], f["n_pairs"], f["n_votes"], PLC._round(f["hit"]),
            None if s is None else (s["evals"], s["n_pairs"], s["n_votes"]),
        )  # fmt: skip

    def l_key(e: Dict[str, Any]) -> Any:
        ri, ev = e["ray_index"], e["evals"] - e["reuse"]["evals"]
        if not ri["fallback_full"]:
            return (False, ev, ri["n_pairs"], ri["n_votes"], PLC._round(e.get("hit_ratio")), None)
        f = ri["first"]
        return (
            True, f["evals"], f["n_pairs"], f["n_votes"], PLC._round(f["hit"]),
            (ev - f["evals"], ri["n_pairs"], ri["n_votes"]),
        )  # fmt: skip

    lib = [e for e in index_log if not e["reuse"]["success"]]
    s_all = [s_key(r) for r in srows]
    l_all = [l_key(e) for e in lib]
    first = next((i for i, (a, b) in enumerate(zip(s_all, l_all)) if a != b), None)
    return dict(
        steps=dict(script=len(s_all), library=len(l_all)),
        fallbacks=dict(
            script=sum(1 for r in srows if r["fallback"]),
            library=sum(1 for e in lib if e["ray_index"]["fallback_full"]),
        ),
        first_differ=first,
        equal=bool(s_all == l_all),
    )


def accuracy_vs_libpxf(lib_mic: Path, base_mic: Path) -> Dict[str, Any]:
    """Wrong voxels (more than ``G.WRONG`` degrees from the truth; symmetry-reduced), the paired
    split by voxel with ``LibPXF``'s map after its pass as the base (fixed: wrong in the base only;
    broken: wrong in ``LibRU`` only) and the median error over all voxels."""
    import stats as ST

    T = G.Truth()
    n = RLC.N
    Rl, _ = G.read_map(lib_mic, n)
    Rb, _ = G.read_map(base_mic, n)
    el, eb = T.err(Rl), T.err(Rb)
    wl, wb = el > G.WRONG, eb > G.WRONG
    broken, fixed = ST.paired_discordant(wl, wb)
    return dict(
        wrong_threshold_deg=float(G.WRONG),
        wrong=dict(libru=int(wl.sum()), libpxf=int(wb.sum())),
        fixed=int(fixed),
        broken=int(broken),
        both_wrong=int((wl & wb).sum()),
        broken_minus_fixed=int(broken) - int(fixed),
        median_error_deg=dict(libru=float(np.median(el)), libpxf=float(np.median(eb))),
        median_higher_deg=float(np.median(el) - np.median(eb)),
        n_voxels=n,
    )


def try_stats(index_log: List[Dict[str, Any]]) -> Dict[str, Any]:
    rows = lib_try_rows(index_log)
    tried = [r for r in rows if r["tried"]]
    ok = [r for r in tried if r["success"]]
    second = [r for r in tried if len(r["fits"]) == 2]
    hits = [r["fits"][-1]["hit"] for r in ok]
    return dict(
        seed_steps=len(rows),
        tried=len(tried),
        success=len(ok),
        success_share_of_steps=len(ok) / max(1, len(rows)),
        second_fits=len(second),
        second_fits_succeeded=sum(1 for r in second if r["success"]),
        evals_per_try=float(np.mean([r["evals"] for r in tried])) if tried else 0.0,
        evals_in_failed_tries=int(sum(r["evals"] for r in tried if not r["success"])),
        lowest_accepted_hit=float(min(hits)) if hits else None,
    )


def prefix_vs_phase_a(index_log: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Rows from the start that equal phase A's stored ``ru_steps_mid`` (rows only; not a gate)."""
    p = SCACHE / "live_ru" / "ru_steps_mid.json.gz"
    a = RLC.read_gz(p)["rows"]
    lib = lib_try_rows(index_log)

    def prefix(same: Callable[[Dict[str, Any], Dict[str, Any]], bool]) -> int:
        n = 0
        for x, y in zip(a, lib):
            if not same(x, y):
                break
            n += 1
        return n

    return dict(
        phase_a_rows=len(a),
        phase_a_successes=sum(1 for r in a if r["success"]),
        equal_prefix_rows=prefix(lambda x, y: row_key_script(x) == row_key_script(y)),
        same_voxel_order_prefix=prefix(lambda x, y: x["voxel"] == y["voxel"]),
        same_voxel_starts_success_prefix=prefix(
            lambda x, y: (x["voxel"], x["tried"], x["starts"], x["success"])
            == (y["voxel"], y["tried"], y["starts"], y["success"])
        ),
        note="consistent with phase A's arm having G2b and aligned draws, which the library "
        "lacks (starts of the same voxels have other orientations, hence other costs); that "
        "cause is not tested",
    )


def script_hashes() -> Dict[str, Any]:
    """sha256 of ``ru_hl.py`` and of every script the freeze (``ru/frozen.json``) lists, as on disk
    now, and whether each equals its frozen hash (the scripts ``RUl`` imports are among them; the
    library files are not listed: they changed in this phase on purpose)."""
    import hashlib

    frozen = json.loads((RD / "frozen.json").read_text())["files"]
    out: Dict[str, Any] = {}
    for rel, want in sorted(frozen.items()):
        p = ICE / rel
        if rel.startswith("scripts/") and p.exists():
            h = hashlib.sha256(p.read_bytes()).hexdigest()
            out[rel] = dict(sha256=h, equals_frozen=bool(h == want))
    p = HERE / "ru_hl.py"
    out["scripts/autoresearch/ar5/ru_hl.py"] = dict(
        sha256=hashlib.sha256(p.read_bytes()).hexdigest(), equals_frozen=None
    )
    return dict(
        n_scripts=len(out),
        all_frozen_scripts_equal=bool(
            all(v["equals_frozen"] for v in out.values() if v["equals_frozen"] is not None)
        ),
        files=out,
    )


def main() -> None:
    assert RC.SAMPLE == "mid"
    hl, lib = SCACHE / "RUl_mid", SCACHE / "LibRU_mid"
    st_hl = json.loads((hl / "run.json").read_text())["stats"]
    run_lib = json.loads((lib / "run.json").read_text())
    st_lib = RLC.lib_stats(run_lib)
    ok_c, bad_c = RLC.compare_counters(st_hl, st_lib, rescue_split=True)
    ilog = RLC.read_gz(lib / "index_log.json.gz")
    slog = [r for r in RLC.read_gz(hl / "rescue_log.json.gz") if r["attempted"]]
    llog = [e for e in ilog if "rescue" in e]

    def tries(ts: List[Dict[str, Any]]) -> List[Any]:
        return [(t["k"], t["hit"], t["cost"], t["evals"]) for t in ts]

    srows = [(int(r["voxel"]), bool(r["rescued"]), tries(r["tries"])) for r in slog]
    lrows = [
        (int(e["voxel"]), bool(e["rescue"]["rescued"]), tries(e["rescue"]["tries"])) for e in llog
    ]
    rs = H5C.stall_rows_script(RLC.read_gz(hl / "mc_traces.json.gz"))
    rl = H5C.stall_rows_library(ilog)
    first_stall = next((i for i, (a, b) in enumerate(zip(rs, rl)) if a != b), None)
    b = dict(
        map=RLC.compare_maps(hl / "recon.mic.gz", lib / "recon_before_pass.mic.gz"),
        evals_total=dict(script=RLC.evals_total(st_hl), library=run_lib["evals_total_before_pass"]),
        counters_equal=ok_c,
        counters_differ=bad_c,
        records=RLC.compare_records(hl / "records.json.gz", lib / "records.json.gz"),
        rescues=dict(script=len(srows), library=len(lrows), equal=srows == lrows),
        stall=dict(
            refinements=dict(script=len(rs), library=len(rl)),
            first_differ=first_stall,
            equal=bool(rs == rl),
        ),
        defer=PLC.compare_defer(hl, st_lib),
        pxf_steps=compare_pxf_steps(hl, ilog),
        tries=compare_tries(hl, ilog),
    )
    l1 = bool(
        b["map"]["voxels_differ"] == 0
        and b["evals_total"]["script"] == b["evals_total"]["library"]
        and b["counters_equal"]
        and b["records"]["n_differ"] == 0
        and b["rescues"]["equal"]
        and b["stall"]["equal"]
        and b["defer"]["equal"]
        and b["pxf_steps"]["equal"]
        and b["tries"]["equal"]
    )
    pxf = json.loads((OLD_CACHE_PXF / "LibPXF_mid" / "run.json").read_text())
    base_evals = int(pxf["evals_total_before_pass"])
    assert base_evals == LIBPXF_EVALS, "the stored LibPXF total is not the one the plan fixes"
    lib_evals = int(run_lib["evals_total_before_pass"])
    fewer = 1.0 - lib_evals / base_evals
    acc = accuracy_vs_libpxf(lib / "recon.mic.gz", OLD_CACHE_PXF / "LibPXF_mid" / "recon.mic.gz")
    l2 = dict(
        evals_before_pass=dict(libru=lib_evals, libpxf=base_evals, fewer_share=fewer),
        evals_ok=bool(fewer >= 0.10),
        accuracy_after_pass=acc,
        net_ok=bool(acc["broken_minus_fixed"] <= 3),
        median_ok=bool(acc["median_higher_deg"] <= 0.001),
    )
    l2["ok"] = bool(l2["evals_ok"] and l2["net_ok"] and l2["median_ok"])
    res = dict(
        sample="mid",
        L1_equal_before_pass=l1,
        L1_detail=b,
        L2=l2,
        tries=try_stats(ilog),
        prefix_vs_phase_a=prefix_vs_phase_a(ilog),
        scripts_sha256=script_hashes(),
        in_run_pass_evals_library=int(run_lib["evals_total_with_pass"]) - lib_evals,
        walls_s=dict(
            library_run=float(run_lib["wall_total_s"]),
            script_arm_RUl=float(json.loads((hl / "run.json").read_text())["wall_total_s"]),
            libpxf_run=float(pxf["wall_total_s"]),
            label="contended (other jobs on the machine); not timings",
        ),
        predictions=dict(
            L2_fewer="15% to 30%", try_share="75% to 88%", broken_minus_fixed="-3 to +3"
        ),
    )
    (RD / "lib_check.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    import doc_tables as DT

    def row(name: str, script: Any, library: Any, same: bool) -> str:
        return f"| {name} | {script} | {library} | {'equal' if same else 'DIFFERENT'} |"

    t = b["tries"]
    lines = [
        "| L1: compared before the boundary pass | script arm RUl | library, key on | verdict |",
        "|---|---|---|---|",
        row("map (decompressed .mic bytes)", "-", "-", b["map"]["voxels_differ"] == 0),
        row(
            "evaluations",
            f"{b['evals_total']['script']:,}",
            f"{b['evals_total']['library']:,}",
            b["evals_total"]["script"] == b["evals_total"]["library"],
        ),
        row("counters of the script run", b["counters_differ"]["n_script_counters"], "all", ok_c),
        row(
            "per-voxel records",
            b["records"]["n_a"],
            b["records"]["n_b"],
            b["records"]["n_differ"] == 0,
        ),
        row("seed steps (try rows)", t["rows"]["script"], t["rows"]["library"], t["equal"]),
        row("seed steps with a try", t["tried"]["script"], t["tried"]["library"], t["equal"]),
        row("tries that succeeded", t["success"]["script"], t["success"]["library"], t["equal"]),
        row(
            "steps with the PXF step",
            b["pxf_steps"]["steps"]["script"],
            b["pxf_steps"]["steps"]["library"],
            b["pxf_steps"]["equal"],
        ),
        row(
            "MC refinements with checks",
            b["stall"]["refinements"]["script"],
            b["stall"]["refinements"]["library"],
            b["stall"]["equal"],
        ),
    ]
    DT.write_tables(
        RD / "lib_tables.md",
        {
            "ru_lib": "\n".join(lines),
            "ru_lib_l2": "\n".join(
                [
                    "| L2: library form, mid | LibPXF | LibRU | verdict |",
                    "|---|---|---|---|",
                    f"| evaluations before the pass | {base_evals:,} | {lib_evals:,} | "
                    f"{fewer:.1%} fewer ({'met' if l2['evals_ok'] else 'not met'}) |",
                    f"| wrong voxels after the pass | {acc['wrong']['libpxf']} | "
                    f"{acc['wrong']['libru']} | broken {acc['broken']}, fixed {acc['fixed']}, "
                    f"net {acc['broken_minus_fixed']} "
                    f"({'met' if l2['net_ok'] else 'not met'}) |",
                    "| median error after the pass (deg) | "
                    f"{acc['median_error_deg']['libpxf']:.6f} | "
                    f"{acc['median_error_deg']['libru']:.6f} | "
                    f"{acc['median_higher_deg']:+.6f} "
                    f"({'met' if l2['median_ok'] else 'not met'}) |",
                ]
            ),
        },
    )
    print("L1 equal before pass:", l1)
    print("L2:", json.dumps(l2, default=float))


if __name__ == "__main__":
    main()
