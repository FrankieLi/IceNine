"""AR5 step 12, phase B: is the key-on library run equal to the script arm ``Hl``? (stop rule)

  RESEED_SAMPLE=mid STALLHIT_CACHE=<cache> RESCUE_CACHE=<cache> \
      uv run python scripts/autoresearch/ar5/stallhit_lib_check.py

Compares, on ``mid``, ``LibH`` (library, ``SeedRescueCMA 0.25 1000`` and
``SeedIndexMCStallHit 3000``) with ``Hl`` (script arm ``H`` with library draws), with the helpers of
step 9's ``rescue_lib_check.py``:

* before the boundary pass: the map (decompressed ``.mic`` bytes; else the voxels that differ and
  the largest angle), the evaluation total, every counter of the script run (time counters left
  out; the script counts its rescue evaluations inside ``index_fallback_grid_evals``, so that count
  is compared with the library's grid plus rescue evaluations), the per-voxel records, the rescues
  (voxel, tries) and, new here, the stall checks: the MC refinements in run order, each with its
  checks ``[evaluation, hit ratio, stopped]``, its checks that went on, whether it stopped and its
  counted evaluations (script: ``mc_traces.json.gz``; library: ``index_log`` ``stall_hit``), and
  the script's ``stallhit.json`` counts against the library's ``index_stallhit_*`` counters;
* after the pass, reported apart (the in-run pass reads the unrounded map in memory, the script
  pass ``order_pass.py --name Hl --aligned 0`` the stored, rounded ``.mic``): evaluations, tries,
  the map, the passes try by try, and the accuracy (truth of ``mid``, a disclosed sample).

Writes ``stallhit/lib_check.json`` and ``stallhit/tables_lib.md``. Only the accuracy after the pass
reads the truth.
"""

import json
import os
from pathlib import Path
from typing import Any, Dict, List

import rescue_lib_check as RLC  # noqa: E402  (first: rescue_common sets the cache)

import numpy as np  # noqa: E402

import reseed_common as RC  # noqa: E402

HERE = Path(__file__).resolve().parent
SH = HERE.parents[2] / "benchmarks" / "autoresearch" / "ar5" / "stallhit"
SCACHE = Path(os.environ.get("STALLHIT_CACHE", str(SH / "cache")))


def stall_rows_script(traces: List[Dict[str, Any]]) -> List[Any]:
    return [
        (
            [list(map(float, c[:2])) + [bool(c[2])] for c in t.get("checks", [])],
            int(t.get("continued", 0)),
            bool(t.get("cut", False)),
            int(t["evals"]),
        )
        for t in traces
    ]


def stall_rows_library(index_log: List[Dict[str, Any]]) -> List[Any]:
    return [
        (
            [list(map(float, c[:2])) + [bool(c[2])] for c in x["checks"]],
            int(x["continued"]),
            bool(x["stopped"]),
            int(x["evals"]),
        )
        for e in index_log
        for x in e.get("stall_hit", [])
    ]


def compare_stall(hl: Path, lib: Path, st_lib: Dict[str, Any]) -> Dict[str, Any]:
    rs = stall_rows_script(RLC.read_gz(hl / "mc_traces.json.gz"))
    rl = stall_rows_library(RLC.read_gz(lib / "index_log.json.gz"))
    first = next((i for i, (a, b) in enumerate(zip(rs, rl)) if a != b), None)
    sj = json.loads((hl / "stallhit.json").read_text())
    counts = dict(
        mc_refinements=(
            int(sj["mc_refinements"]),
            int(st_lib.get("index_stallhit_refinements", 0)),
        ),
        checks=(int(sj["checks_total"]), int(st_lib.get("index_stallhit_checks", 0))),
        stops=(int(sj["cut"]), int(st_lib.get("index_stallhit_stops", 0))),
        went_on=(int(sj["checks_went_on"]), int(st_lib.get("index_stallhit_went_on", 0))),
        evals_in_stopped=(
            int(sj["evals_in_cut"]),
            int(st_lib.get("index_stallhit_evals_in_stopped", 0)),
        ),
    )
    return dict(
        refinements=dict(script=len(rs), library=len(rl)),
        first_differ=first,
        equal=bool(rs == rl),
        counts_script_library=counts,
        counts_equal=bool(all(a == b for a, b in counts.values())),
    )


def lib_vs_hl(C: Path) -> Dict[str, Any]:
    hl, lib = C / "Hl_mid", C / "LibH_mid"
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
            tries_script=sum(len(r[2]) for r in srows),
            tries_library=sum(len(r[2]) for r in lrows),
            equal=srows == lrows,
        ),
        stall=compare_stall(hl, lib, st_lib),
    )
    out: Dict[str, Any] = dict(before_pass=before)
    pdir = C / "passes_mid"
    if (pdir / "D4_Hl" / "pass.json").exists():
        pj = json.loads((pdir / "D4_Hl" / "pass.json").read_text())
        lt = np.load(lib / "pass_tries.npz")
        out["after_pass"] = dict(
            pass_evals=dict(script=int(pj["evals"]), library=int(lt["evals"].sum())),
            pass_tries=dict(script=int(pj["n_tries"]), library=int(len(lt["evals"]))),
            total_with_pass=dict(
                script=RLC.evals_total(st_hl) + int(pj["evals"]),
                library=int(run_lib["evals_total_with_pass"]),
            ),
            map=RLC.compare_maps(lib / "recon.mic.gz", np.load(pdir / "after_Hl.npy")),
            tries=RLC.try_level(lib / "pass_tries.npz", pdir / "D4_Hl" / "tries.npz"),
            accuracy=RLC.accuracy_after(lib / "recon.mic.gz", pdir / "after_Hl.npy"),
        )
        out["compute"] = dict(  # the phase B budget: the two mid runs together
            evals_before_passes=int(RLC.evals_total(st_hl) + run_lib["evals_total_before_pass"]),
            evals_with_passes=int(
                out["after_pass"]["total_with_pass"]["script"]
                + out["after_pass"]["total_with_pass"]["library"]
            ),
        )
    b = before
    out["equal_before_pass"] = bool(
        b["map"]["voxels_differ"] == 0
        and b["evals_total"]["script"] == b["evals_total"]["library"]
        and b["counters_equal"]
        and b["records"]["n_differ"] == 0
        and b["rescues"]["equal"]
        and b["stall"]["equal"]
        and b["stall"]["counts_equal"]
    )
    return out


def tables(res: Dict[str, Any]) -> None:
    import sys

    sys.path.insert(0, str(HERE.parents[1] / "common"))
    import doc_tables as DT

    def yn(x: bool) -> str:
        return "yes" if x else "NO"

    r = res["libH_vs_Hl"]
    b = r["before_pass"]
    ev, z, s = b["evals_total"], b["rescues"], b["stall"]
    cs = s["counts_script_library"]
    es = cs["evals_in_stopped"]
    rows = [
        "| item (before the boundary pass) | library run LibH | script arm Hl | equal |",
        "|---|---|---|---|",
        f"| map: voxels that differ | {b['map']['voxels_differ']} | - |"
        f" {yn(b['map']['voxels_differ'] == 0)} |",
        f"| evaluations | {ev['library']:,} | {ev['script']:,} |"
        f" {yn(ev['library'] == ev['script'])} |",
        f"| counters of the script run ({b['counters_differ']['n_script_counters']}) | - | - |"
        f" {yn(b['counters_equal'])} |",
        f"| per-voxel records ({b['records']['n_a']:,}) | - | - |"
        f" {yn(b['records']['n_differ'] == 0)} |",
        f"| rescues: tried, rescued, tries | {z['library']}, {z['rescued_library']},"
        f" {z['tries_library']} | {z['script']}, {z['rescued_script']}, {z['tries_script']} |"
        f" {yn(z['equal'])} |",
        f"| MC refinements with the stall check (each: checks, stop, evaluations) |"
        f" {s['refinements']['library']:,} | {s['refinements']['script']:,} | {yn(s['equal'])} |",
        f"| stall checks; stopped; went on | {cs['checks'][1]:,}; {cs['stops'][1]:,};"
        f" {cs['went_on'][1]:,} | {cs['checks'][0]:,}; {cs['stops'][0]:,}; {cs['went_on'][0]:,} |"
        f" {yn(s['counts_equal'])} |",
        f"| evaluations in stopped refinements | {es[1]:,} | {es[0]:,} | {yn(es[0] == es[1])} |",
    ]
    T = {
        "stallhit_libcheck": "\n".join(rows)
        + f"\n\nEqual before the pass: {yn(r['equal_before_pass'])}."
    }
    a = r.get("after_pass")
    if a:
        m, acc, t = a["map"], a["accuracy"], a["tries"]
        rows = [
            "| item | library run LibH: pass in the run | script arm Hl: pass after the run |",
            "|---|---|---|",
            f"| pass tries | {a['pass_tries']['library']:,} | {a['pass_tries']['script']:,} |",
            f"| pass evaluations | {a['pass_evals']['library']:,} |"
            f" {a['pass_evals']['script']:,} |",
            f"| total with the pass | {a['total_with_pass']['library']:,} |"
            f" {a['total_with_pass']['script']:,} |",
            "| tries in the same (voxel, candidate) order from the start |"
            f" {t['same_order_tries']:,} | - |",
            f"| map after the pass against the script's: voxels that differ in any bit /"
            f" more than 0.001 degree apart / more than 1 degree apart | {m['voxels_differ']:,} /"
            f" {m['voxels_over_0_001_deg']:,} / {m['voxels_over_1_deg']:,} | - |",
            f"| wrong voxels after the pass (more than {acc['wrong_threshold_deg']:g} degree from"
            f" the truth) | {acc['wrong']['library']:,} | {acc['wrong']['script']:,} |",
            f"| by voxel, script pass as the base: fixed / broken / wrong in both |"
            f" {acc['fixed']} / {acc['broken']} / {acc['both_wrong']} | - |",
            f"| median error after the pass, all {acc['n_voxels']:,} voxels (degrees) |"
            f" {acc['median_error_deg']['library']:.5f} |"
            f" {acc['median_error_deg']['script']:.5f} |",
        ]
        T["stallhit_libpass"] = "\n".join(rows)
    DT.write_tables(SH / "tables_lib.md", T)


def main() -> None:
    assert RC.SAMPLE == "mid"
    C = SCACHE
    res: Dict[str, Any] = dict(sample="mid")
    if (C / "LibH_mid" / "run.done").exists() and (C / "Hl_mid" / "Hl.done").exists():
        res["libH_vs_Hl"] = lib_vs_hl(C)
        tables(res)
    (SH / "lib_check.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    print(json.dumps(res, indent=1, default=float), flush=True)
    os._exit(0)


if __name__ == "__main__":
    main()
