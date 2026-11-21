"""AR5 step 28, phase B: is the key-on library run equal to the script arm ``PXFl``? (stop rule)

  RESEED_SAMPLE=mid PXF_CACHE=<cache> uv run python scripts/autoresearch/ar5/pxf_lib_check.py

Compares, on ``mid`` and before the boundary pass, ``LibPXF`` (library: ``SeedRescueCMA 0.25 1000``,
``SeedIndexMCStallHit 500``, ``SeedIndexMCDefer 1``, ``RayIndexMaxPairs 50000``,
``RayIndexFallbackFull 1``) with ``PXFl`` (step 22's ``J5Dl`` with rule ``PXF``; ``pxf_hl.py``),
with the helpers of step 22's ``defer_lib_check.py``: the map (decompressed ``.mic`` bytes), the
evaluation total, every counter of the script run, the per-voxel records, the rescues, the MC
refinements in run order with their checks, the pauses and resumes, and (new) the fallbacks: the
steps that ran twice, both steps' evaluations, hit ratios and pair counts (script:
``pxf_steps.json.gz``; library: ``index_log`` ``ray_index``).

The phase A arm ``J5D`` (rule ``PXF``, ``G2b``, aligned draws) is not the comparison: the library
has neither ``G2b`` nor aligned draws (as in step 22). Writes ``pxf/lib_check.json``. No truth is
read. Does not go through ``pxf_common``.
"""

import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
SD = ICE / "benchmarks" / "autoresearch" / "ar5" / "pxf"
SCACHE = Path(os.environ.get("PXF_CACHE", str(SD / "cache")))
os.environ["HIGH500_CACHE"] = str(SCACHE)
os.environ["RESCUE_CACHE"] = str(SCACHE)
sys.path.insert(0, str(HERE))

import rescue_lib_check as RLC  # noqa: E402  (first: rescue_common sets the cache)

import high500_lib_check as H5C  # noqa: E402
import reseed_common as RC  # noqa: E402


def compare_defer(hl: Path, st_lib: Dict[str, Any]) -> Dict[str, Any]:
    dj = json.loads((hl / "defer.json").read_text())
    counts = dict(
        steps_with_pause=(int(dj["steps_with_pause"]), int(st_lib.get("index_defer_steps", 0))),
        pauses=(int(dj["pauses"]), int(st_lib.get("index_defer_pauses", 0))),
        resumes=(int(dj["resumes"]), int(st_lib.get("index_defer_resumes", 0))),
    )
    return dict(counts_script_library=counts, equal=bool(all(a == b for a, b in counts.values())))


def _round(x: Any) -> Any:
    return None if x is None else round(float(x), 12)


def compare_fallbacks(hl: Path, index_log: List[Dict[str, Any]]) -> Dict[str, Any]:
    """The script's PXF rows against the library's ``ray_index`` entries, step by step in run
    order: whether the step ran twice, both steps' evaluations, hit ratios and pair counts."""
    srows = RLC.read_gz(hl / "pxf_steps.json.gz")["rows"]

    def s_key(r: Dict[str, Any]) -> Any:
        f, s = r["first"], r.get("second")
        return (
            bool(r["fallback"]),
            f["evals"],
            f["n_pairs"],
            f["n_votes"],
            _round(f["hit"]),
            None if s is None else (s["evals"], s["n_pairs"], s["n_votes"]),
        )

    def l_key(e: Dict[str, Any]) -> Any:
        ri = e["ray_index"]
        if not ri["fallback_full"]:
            return (
                False,
                e["evals"],
                ri["n_pairs"],
                ri["n_votes"],
                _round(e.get("hit_ratio")),
                None,
            )
        f = ri["first"]
        return (
            True,
            f["evals"],
            f["n_pairs"],
            f["n_votes"],
            _round(f["hit"]),
            (e["evals"] - f["evals"], ri["n_pairs"], ri["n_votes"]),
        )

    s_all = [s_key(r) for r in srows]
    l_all = [l_key(e) for e in index_log]
    first = next((i for i, (a, b) in enumerate(zip(s_all, l_all)) if a != b), None)
    return dict(
        steps=dict(script=len(s_all), library=len(l_all)),
        fallbacks=dict(
            script=sum(1 for r in srows if r["fallback"]),
            library=sum(1 for e in index_log if e["ray_index"]["fallback_full"]),
        ),
        passed_second=dict(
            script=sum(1 for r in srows if r.get("passed_second")),
            library=sum(1 for e in index_log if e["ray_index"]["fallback_full"] and e.get("used")),
        ),
        first_differ=first,
        equal=bool(s_all == l_all),
    )


def main() -> None:
    assert RC.SAMPLE == "mid"
    hl, lib = SCACHE / "PXFl_mid", SCACHE / "LibPXF_mid"
    st_hl = json.loads((hl / "run.json").read_text())["stats"]
    run_lib = json.loads((lib / "run.json").read_text())
    st_lib = RLC.lib_stats(run_lib)
    ok_c, bad_c = RLC.compare_counters(st_hl, st_lib, rescue_split=True)
    slog = [r for r in RLC.read_gz(hl / "rescue_log.json.gz") if r["attempted"]]
    ilog = RLC.read_gz(lib / "index_log.json.gz")
    llog = [e for e in ilog if "rescue" in e]

    def tries(ts: List[Dict[str, Any]]) -> List[Any]:
        return [(t["k"], t["hit"], t["cost"], t["evals"]) for t in ts]

    srows = [(int(r["voxel"]), bool(r["rescued"]), tries(r["tries"])) for r in slog]
    lrows = [
        (int(e["voxel"]), bool(e["rescue"]["rescued"]), tries(e["rescue"]["tries"])) for e in llog
    ]
    rs = H5C.stall_rows_script(RLC.read_gz(hl / "mc_traces.json.gz"))
    rl = H5C.stall_rows_library(ilog)
    first = next((i for i, (a, b) in enumerate(zip(rs, rl)) if a != b), None)
    b = dict(
        map=RLC.compare_maps(hl / "recon.mic.gz", lib / "recon_before_pass.mic.gz"),
        evals_total=dict(script=RLC.evals_total(st_hl), library=run_lib["evals_total_before_pass"]),
        counters_equal=ok_c,
        counters_differ=bad_c,
        records=RLC.compare_records(hl / "records.json.gz", lib / "records.json.gz"),
        rescues=dict(
            script=len(srows),
            library=len(lrows),
            rescued_script=sum(1 for r in srows if r[1]),
            rescued_library=sum(1 for r in lrows if r[1]),
            equal=srows == lrows,
        ),
        stall=dict(
            refinements=dict(script=len(rs), library=len(rl)),
            first_differ=first,
            equal=bool(rs == rl),
        ),
        defer=compare_defer(hl, st_lib),
        fallbacks=compare_fallbacks(hl, ilog),
    )
    equal = bool(
        b["map"]["voxels_differ"] == 0
        and b["evals_total"]["script"] == b["evals_total"]["library"]
        and b["counters_equal"]
        and b["records"]["n_differ"] == 0
        and b["rescues"]["equal"]
        and b["stall"]["equal"]
        and b["defer"]["equal"]
        and b["fallbacks"]["equal"]
    )
    # for reference only: the phase A arm J5D (rule PXF with G2b, aligned draws; the live run) is
    # not the comparison
    j5d = json.loads((SCACHE / "live_pxf" / "J5D_mid" / "run.json").read_text())["stats"]
    res = dict(
        sample="mid",
        libPXF_vs_PXFl=dict(
            before_pass=b,
            equal_before_pass=equal,
            in_run_pass_evals_library=int(run_lib["evals_total_with_pass"])
            - int(run_lib["evals_total_before_pass"]),
        ),
        reference=dict(
            phase_a_J5D_search_evals=int(RLC.evals_total(j5d)),
            libPXF_evals_before_pass=int(run_lib["evals_total_before_pass"]),
            note="J5D has rule G2b and aligned draws; the library has neither; not compared",
        ),
    )
    (SD / "lib_check.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    sys.path.insert(0, str(ICE / "scripts" / "common"))
    import doc_tables as DT

    def row(name: str, script: Any, library: Any, same: bool) -> str:
        return f"| {name} | {script} | {library} | {'equal' if same else 'DIFFERENT'} |"

    def same(pair: Any) -> bool:
        return bool(pair[0] == pair[1])

    fb, st, df = b["fallbacks"], b["stall"], b["defer"]["counts_script_library"]
    lines = [
        "| compared before the boundary pass | script arm PXFl | library, keys on | verdict |",
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
        row("seed steps", fb["steps"]["script"], fb["steps"]["library"], fb["equal"]),
        row("second tries", fb["fallbacks"]["script"], fb["fallbacks"]["library"], fb["equal"]),
        row(
            "second tries that passed",
            fb["passed_second"]["script"],
            fb["passed_second"]["library"],
            fb["equal"],
        ),
        row(
            "rescues (rescued)",
            f"{b['rescues']['script']} ({b['rescues']['rescued_script']})",
            f"{b['rescues']['library']} ({b['rescues']['rescued_library']})",
            b["rescues"]["equal"],
        ),
        row(
            "MC refinements with checks",
            st["refinements"]["script"],
            st["refinements"]["library"],
            st["equal"],
        ),
        row("steps with a pause", *df["steps_with_pause"], same(df["steps_with_pause"])),
        row("pauses", *df["pauses"], same(df["pauses"])),
        row("resumes", *df["resumes"], same(df["resumes"])),
    ]
    DT.write_tables(SD / "lib_tables.md", {"pxf_lib": "\n".join(lines)})
    print("equal_before_pass:", equal)
    print(json.dumps(b["evals_total"]), json.dumps(b["fallbacks"]))


if __name__ == "__main__":
    main()
