"""AR5 step 22, phase B: is the key-on library run equal to the script arm ``J5Dl``? (stop rule)

  RESEED_SAMPLE=mid DEFER_CACHE=<cache> uv run python scripts/autoresearch/ar5/defer_lib_check.py

Compares, on ``mid`` and before the boundary pass, ``LibJ5D`` (library: ``SeedRescueCMA 0.25 1000``,
``SeedIndexMCStallHit 500``, ``SeedIndexMCDefer 1``) with ``J5Dl`` (step 16's ``J5l`` with rule
``DEF``: ``J5D`` without ``G2b``, with the library's own draws; ``defer_hl.py``), with the helpers
of step 9's ``rescue_lib_check.py`` and step 16's ``high500_lib_check.py``: the map (decompressed
``.mic`` bytes; else the voxels that differ and the largest angle), the evaluation total, every
counter of the script run (the script counts its rescue evaluations inside
``index_fallback_grid_evals``, the library apart), the per-voxel records, the rescues (voxel,
tries), and the MC refinements in run order, each with its checks ``[evaluation, hit ratio,
stopped]``, its checks that went on, whether it stopped and its counted evaluations (script:
``mc_traces.json.gz``,
whose rows include the paused part of a refinement never resumed; library: ``index_log``
``stall_hit``), and the pauses and resumes of the script (``defer.json``) against the library's
``index_defer_*`` counters. The in-run boundary pass is reported apart (it is not part of the
equality: the script arm makes none).

``J5D`` itself is not the comparison: it has rule ``G2b`` and aligned draws, which the library does
not have (as ``J5`` was not in step 16: ``J5l`` was). Writes ``defer/lib_check.json``. No truth is
read. Does not go through ``defer_common`` (see ``defer_lib_run.py``).
"""

import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
SD = ICE / "benchmarks" / "autoresearch" / "ar5" / "defer"
SCACHE = Path(os.environ.get("DEFER_CACHE", str(SD / "cache")))
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


def main() -> None:
    assert RC.SAMPLE == "mid"
    hl, lib = SCACHE / "J5Dl_mid", SCACHE / "LibJ5D_mid"
    st_hl = json.loads((hl / "run.json").read_text())["stats"]
    run_lib = json.loads((lib / "run.json").read_text())
    st_lib = RLC.lib_stats(run_lib)
    ok_c, bad_c = RLC.compare_counters(st_hl, st_lib, rescue_split=True)
    slog = [r for r in RLC.read_gz(hl / "rescue_log.json.gz") if r["attempted"]]
    llog = [e for e in RLC.read_gz(lib / "index_log.json.gz") if "rescue" in e]

    def tries(ts: List[Dict[str, Any]]) -> List[Any]:
        return [(t["k"], t["hit"], t["cost"], t["evals"]) for t in ts]

    srows = [(int(r["voxel"]), bool(r["rescued"]), tries(r["tries"])) for r in slog]
    lrows = [
        (int(e["voxel"]), bool(e["rescue"]["rescued"]), tries(e["rescue"]["tries"])) for e in llog
    ]
    rs = H5C.stall_rows_script(RLC.read_gz(hl / "mc_traces.json.gz"))
    rl = H5C.stall_rows_library(RLC.read_gz(lib / "index_log.json.gz"))
    first = next((i for i, (a, b) in enumerate(zip(rs, rl)) if a != b), None)
    before = dict(
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
    )
    b = before
    res = dict(
        sample="mid",
        libJ5D_vs_J5Dl=dict(
            before_pass=before,
            equal_before_pass=bool(
                b["map"]["voxels_differ"] == 0
                and b["evals_total"]["script"] == b["evals_total"]["library"]
                and b["counters_equal"]
                and b["records"]["n_differ"] == 0
                and b["rescues"]["equal"]
                and b["stall"]["equal"]
                and b["defer"]["equal"]
            ),
            in_run_pass_evals_library=int(run_lib["evals_total_with_pass"])
            - int(run_lib["evals_total_before_pass"]),
        ),
    )
    # what the key gives in library form (against step 16's LibJ5 on mid, before the pass) and why
    # the script arm J5D itself is not the comparison (its search is cheaper by this much)
    lib5 = json.loads((SD.parent / "high500" / "lib_check.json").read_text())
    j5_lib = int(lib5["libJ5_vs_J5l"]["before_pass"]["evals_total"]["library"])
    j5d = SCACHE / "J5D_mid"
    st_d = json.loads((j5d / "run.json").read_text())["stats"]
    pol_d = json.loads((j5d / "polish.json").read_text())
    d_search = RLC.evals_total(st_d)
    d_fb = int(pol_d["rescue_evals"] + pol_d["polish_evals"] + pol_d["grid_evals"])
    l_fb = int(st_hl["index_fallback_grid_evals"])  # the script arm's rescue + grid evaluations
    lib_total = int(run_lib["evals_total_before_pass"])
    res["library_form"] = dict(
        libJ5_evals_before_pass=j5_lib,
        libJ5D_evals_before_pass=lib_total,
        saving_share=1.0 - lib_total / j5_lib,
        j5d_search_evals=int(d_search),
        gap_to_j5d=int(lib_total - d_search),
        gap_in_fallbacks=int(l_fb - d_fb),
        j5dl_fallback_evals=l_fb,
        j5d_fallback_evals=d_fb,
        grid_searches=dict(
            j5dl=int(json.loads((hl / "rescue.json").read_text())["grid_fallbacks_left"]),
            j5d=int(pol_d["grid_fallbacks_left"]),
        ),
        screen_pass_run_in_library="no",
        accuracy_measured="no",
    )
    (SD / "lib_check.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    print(json.dumps(res, indent=1, default=float), flush=True)
    os._exit(0)


if __name__ == "__main__":
    main()
