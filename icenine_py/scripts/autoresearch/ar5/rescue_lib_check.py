"""AR5 step 9, phase B: is the key-on library run equal to the script arm ``Gl``? (stop rule)

  RESEED_SAMPLE=mid uv run python scripts/autoresearch/ar5/rescue_lib_check.py

Compares, on ``mid``:

* ``LibG`` (library, ``SeedRescueCMA 0.25 1000``) with ``Gl`` (script arm ``G`` with library
  draws), before the boundary pass: the map (decompressed ``.mic`` bytes; else the voxels that
  differ and the largest angle), the evaluation total and every counter of the run (time counters
  left out; the script counts its rescue evaluations inside ``index_fallback_grid_evals``, the
  library apart in ``index_rescue_evals``, so the script's value is compared with the sum), the
  per-voxel records, and the rescues (voxel, tries: kept end, hit ratio, cost, evaluations).
  After the pass: the library's in-run pass against ``order_pass.py --name Gl --aligned 0``
  (evaluations, tries, map). The in-run pass starts from the unrounded map in memory, the script
  pass from the stored ``.mic`` (rounded angles). Try by try: how long the two passes keep the same
  (voxel, candidate) order, and over those tries how many ends (compared as float32, the script's
  stored type), costs and evaluation counts differ, the first try whose end differs, and the angle
  between the ends (no symmetry: the same candidate's two ends). The accuracy after the pass reads
  the truth of ``mid`` (a disclosed sample): wrong voxels of each pass and the split by voxel.
* the control ``LibS11`` (library, key off) with the stored script run ``S11`` of AR5 step 3 (read
  only, in its kept worktree), before the pass, the same way, and its passes try by try.

Writes ``rescue/lib_check.json`` and ``rescue/tables_lib.md``. Only the accuracy after the pass
reads the truth (of ``mid``).
"""

import gzip
import json
import os
from pathlib import Path
from typing import Any, Dict, List, Tuple

import rescue_common as RCM  # noqa: E402  (first)

import numpy as np  # noqa: E402

import gradient_common as G  # noqa: E402
import reseed_common as RC  # noqa: E402

N = 24570
MIX = RC.WORKTREES / "agent-a5a3601b746680640/icenine_py/benchmarks/autoresearch/ar5/mix/cache"


def read_gz(p: Path) -> Any:
    with gzip.open(p, "rt") as f:
        return json.load(f)


def gz_bytes(p: Path) -> bytes:
    with gzip.open(p, "rb") as f:
        return f.read()


def is_time_key(k: str) -> bool:
    return "wall" in k or k.endswith("_s")


def counters(st: Dict[str, Any]) -> Dict[str, Any]:
    return {
        k: v
        for k, v in st.items()
        if isinstance(v, (int, float)) and not isinstance(v, bool) and not is_time_key(k)
    }


def evals_total(st: Dict[str, Any]) -> int:
    return int(
        sum(st[k] for k in ("n_evals_seed", "n_evals_neighbor", "n_evals_revisit"))
        + st.get("n_evals_reseed", 0)
        + st.get("n_evals_restart", 0)
    )


def mis_deg(Ra: np.ndarray, Rb: np.ndarray) -> np.ndarray:
    """Symmetry-reduced angle (cubic) between matching orientations, degrees."""
    from icenine.ray_indexing import misorientation_deg, proper_rotations
    from icenine.symmetry import create_cubic_symmetry

    if len(Ra) == 0:
        return np.zeros(0)
    sym = proper_rotations(create_cubic_symmetry(3.6))
    return np.asarray(misorientation_deg(Ra, Rb, sym), dtype=float)


def compare_maps(pa: Path, pb_or_R: Any) -> Dict[str, Any]:
    """Voxels whose stored orientation differs (any bit), and among them the symmetry-reduced
    angle: the largest and how many differ by more than 0.001 degree (a .mic stores rounded
    angles, so an unrounded map differs from its stored form by less than that)."""
    Ra, _ = G.read_map(pa, N)
    if isinstance(pb_or_R, Path):
        Rb, _ = G.read_map(pb_or_R, N)
        same_bytes = gz_bytes(pa) == gz_bytes(pb_or_R)
    else:
        Rb, same_bytes = np.asarray(pb_or_R, np.float64).reshape(N, 3, 3), None
    diff = np.where(~np.all(Ra.reshape(N, -1) == Rb.reshape(N, -1), axis=1))[0]
    ang = mis_deg(Ra[diff], Rb[diff])
    return dict(
        same_bytes=same_bytes,
        voxels_differ=int(len(diff)),
        max_angle_deg=float(ang.max()) if len(ang) else 0.0,
        voxels_over_0_001_deg=int((ang > 1e-3).sum()),
        voxels_over_1_deg=int((ang > 1.0).sum()),
    )


# counters the library sets at the end of the run, after the boundary pass (reconstructor.py,
# ``reconstruct_sample``): they are read from the run's final stats, not from the capture made
# when the pass starts
END_KEYS = ("n_unresolved", "n_voxels")


def lib_stats(run_lib: Dict[str, Any]) -> Dict[str, Any]:
    st = dict(run_lib["stats_before_pass"])
    for k in END_KEYS:
        st[k] = run_lib["stats"][k]
    return st


def compare_counters(
    a: Dict[str, Any], b: Dict[str, Any], rescue_split: bool
) -> Tuple[bool, Dict[str, Any]]:
    """``a`` script, ``b`` library; every counter of the script run is compared (as in
    ``keys_check.py``); counters only the library keeps (its own bookkeeping of keys the script
    layer emulated) are listed. With ``rescue_split`` the script's grid evaluations hold the
    rescue evaluations (library: ``index_rescue_evals``)."""
    ca, cb = counters(a), counters(b)
    if rescue_split:
        cb = dict(cb)
        cb["index_fallback_grid_evals"] = cb.get("index_fallback_grid_evals", 0) + cb.get(
            "index_rescue_evals", 0
        )
    differ = {k: (ca[k], cb.get(k)) for k in sorted(ca) if ca[k] != cb.get(k)}
    out = dict(
        n_script_counters=len(ca),
        differ=differ,
        only_in_library={k: cb[k] for k in sorted(set(cb) - set(ca))},
    )
    return not differ, out


def compare_records(pa: Path, pb: Path) -> Dict[str, Any]:
    ra, rb = read_gz(pa), read_gz(pb)

    def strip(r: Dict[str, Any]) -> Dict[str, Any]:
        return {k: v for k, v in r.items() if k != "wall_s"}

    keys = sorted(set(ra) | set(rb), key=int)
    bad = [k for k in keys if strip(ra.get(k, {})) != strip(rb.get(k, {}))]
    return dict(n_a=len(ra), n_b=len(rb), n_differ=len(bad), first=bad[:5])


def try_level(lib_npz: Path, scr_npz: Path) -> Dict[str, Any]:
    """The library's in-run pass tries against the script's, try by try. ``same_order_tries``:
    how many tries from the start have the same (voxel, candidate) in both. Over those tries:
    ends that differ (the library's float64 end cast to the script's stored float32), costs and
    evaluation counts that differ, the first try whose end differs, and the angle between the two
    ends of each try (degrees, no symmetry). Same targets and the same fresh generator would give
    an equal first try if the input maps were equal; a shifted generator would keep the early
    tries equal."""
    lt, st = np.load(lib_npz), np.load(scr_npz)
    nl, ns = len(lt["voxel"]), len(st["voxel"])
    n = min(nl, ns)
    same = (lt["voxel"][:n] == st["voxel"][:n]) & (lt["cand"][:n] == st["cand"][:n])
    m = n if same.all() else int(np.argmin(same))
    le = lt["end"][:m].astype(np.float32).reshape(m, -1)
    se = st["end"][:m].astype(np.float32).reshape(m, -1)
    end_diff = ~np.all(le == se, axis=1)
    Rl = lt["end"][:m].astype(np.float64).reshape(m, 3, 3)
    Rs = st["end"][:m].astype(np.float64).reshape(m, 3, 3)
    tr = np.einsum("nij,nij->n", Rl, Rs)  # trace(Rl Rs^T)
    ang = np.degrees(np.arccos(np.clip((tr - 1.0) / 2.0, -1.0, 1.0)))
    cl = np.bincount(lt["voxel"], minlength=N)
    cs = np.bincount(st["voxel"], minlength=N)
    return dict(
        tries=dict(library=int(nl), script=int(ns)),
        same_order_tries=int(m),
        ends_differ=int(end_diff.sum()),
        costs_differ=int((lt["cost"][:m] != st["cost"][:m]).sum()),
        evals_differ=int((lt["evals"][:m] != st["evals"][:m]).sum()),
        first_end_differ=int(np.argmax(end_diff)) if end_diff.any() else None,
        try0=dict(
            voxel=int(lt["voxel"][0]),
            cand=int(lt["cand"][0]),
            same_cost=bool(lt["cost"][0] == st["cost"][0]),
            angle_deg=float(ang[0]) if m else None,
        ),
        end_angle_deg=dict(
            median=float(np.median(ang)) if m else None, max=float(ang.max()) if m else None
        ),
        voxels_with_other_try_count=int((cl != cs).sum()),
    )


def accuracy_after(lib_mic: Path, scr_after: Path) -> Dict[str, Any]:
    """Accuracy of the map after each pass, against the truth of ``mid``: wrong voxels (more
    than ``G.WRONG`` degrees from the truth), the split by voxel with the script's pass as the
    base (fixed: wrong after the script's pass, right after the in-run pass; broken: the
    reverse), and the median error over all voxels."""
    sys_path_common()
    import stats as ST

    T = G.Truth()
    Rl, _ = G.read_map(lib_mic, N)
    Rs = np.load(scr_after).reshape(N, 3, 3)
    el, es = T.err(Rl), T.err(Rs)
    wl, ws = el > G.WRONG, es > G.WRONG
    broken, fixed = ST.paired_discordant(wl, ws)  # (wrong in-run only, wrong script only)
    return dict(
        wrong_threshold_deg=float(G.WRONG),
        wrong=dict(library=int(wl.sum()), script=int(ws.sum())),
        fixed=int(fixed),
        broken=int(broken),
        both_wrong=int((wl & ws).sum()),
        median_error_deg=dict(library=float(np.median(el)), script=float(np.median(es))),
        n_voxels=N,
        base="script pass (Gl, from the stored map)",
    )


def sys_path_common() -> None:
    """``scripts/common`` (``stats.py``) on the import path."""
    import sys

    common = str(Path(__file__).resolve().parents[2] / "common")
    if common not in sys.path:
        sys.path.insert(0, common)


def lib_vs_gl(C: Path) -> Dict[str, Any]:
    gl, lib = C / "Gl_mid", C / "LibG_mid"
    st_gl = json.loads((gl / "run.json").read_text())["stats"]
    run_lib = json.loads((lib / "run.json").read_text())
    st_lib = lib_stats(run_lib)
    ok_c, bad_c = compare_counters(st_gl, st_lib, rescue_split=True)
    # rescues: script log rows vs library index_log entries with a rescue
    slog = [r for r in read_gz(gl / "rescue_log.json.gz") if r["attempted"]]
    llog = [e for e in read_gz(lib / "index_log.json.gz") if "rescue" in e]

    def tries(ts: List[Dict[str, Any]]) -> List[Any]:
        return [(t["k"], t["hit"], t["cost"], t["evals"]) for t in ts]

    srows = [(int(r["voxel"]), bool(r["rescued"]), tries(r["tries"])) for r in slog]
    lrows = [
        (int(e["voxel"]), bool(e["rescue"]["rescued"]), tries(e["rescue"]["tries"])) for e in llog
    ]
    before = dict(
        map=compare_maps(gl / "recon.mic.gz", lib / "recon_before_pass.mic.gz"),
        evals_total=dict(script=evals_total(st_gl), library=run_lib["evals_total_before_pass"]),
        counters_equal=ok_c,
        counters_differ=bad_c,
        records=compare_records(gl / "records.json.gz", lib / "records.json.gz"),
        rescues=dict(
            script=len(srows),
            library=len(lrows),
            rescued_script=sum(1 for r in srows if r[1]),
            rescued_library=sum(1 for r in lrows if r[1]),
            tries_script=sum(len(r[2]) for r in srows),
            tries_library=sum(len(r[2]) for r in lrows),
            rescue_evals_library=int(st_lib.get("index_rescue_evals", 0)),
            equal=srows == lrows,
        ),
    )
    out: Dict[str, Any] = dict(before_pass=before)
    pdir = C / "passes_mid"
    if (pdir / "D4_Gl" / "pass.json").exists():
        pj = json.loads((pdir / "D4_Gl" / "pass.json").read_text())
        lt = np.load(lib / "pass_tries.npz")
        out["after_pass"] = dict(
            pass_evals=dict(script=int(pj["evals"]), library=int(lt["evals"].sum())),
            pass_tries=dict(script=int(pj["n_tries"]), library=int(len(lt["evals"]))),
            total_with_pass=dict(
                script=evals_total(st_gl) + int(pj["evals"]),
                library=int(run_lib["evals_total_with_pass"]),
            ),
            map=compare_maps(lib / "recon.mic.gz", np.load(pdir / "after_Gl.npy")),
            tries=try_level(lib / "pass_tries.npz", pdir / "D4_Gl" / "tries.npz"),
            accuracy=accuracy_after(lib / "recon.mic.gz", pdir / "after_Gl.npy"),
        )
        lp = RCM.RESCUE / "lib_pass.json"
        if lp.exists():  # the library pass on the stored map (rescue_lib_pass.py)
            out["library_pass_on_stored_map"] = json.loads(lp.read_text())
    b = before
    out["equal_before_pass"] = bool(
        b["map"]["voxels_differ"] == 0
        and b["evals_total"]["script"] == b["evals_total"]["library"]
        and b["counters_equal"]
        and b["records"]["n_differ"] == 0
        and b["rescues"]["equal"]
    )
    return out


def control(C: Path) -> Dict[str, Any]:
    lib, s11 = C / "LibS11_mid", MIX / "S11_mid"
    st_s = json.loads((s11 / "run.json").read_text())["stats"]
    run_lib = json.loads((lib / "run.json").read_text())
    st_lib = lib_stats(run_lib)
    ok_c, bad_c = compare_counters(st_s, st_lib, rescue_split=False)
    b = dict(
        map=compare_maps(s11 / "recon.mic.gz", lib / "recon_before_pass.mic.gz"),
        evals_total=dict(script=evals_total(st_s), library=run_lib["evals_total_before_pass"]),
        counters_equal=ok_c,
        counters_differ=bad_c,
        records=compare_records(s11 / "records.json.gz", lib / "records.json.gz"),
    )
    pj = json.loads((MIX / "passes_mid" / "D4_S11" / "pass.json").read_text())
    lt = np.load(lib / "pass_tries.npz")
    after = dict(  # the in-run pass (unrounded map) against the stored script pass (stored map)
        pass_tries=dict(library=int(len(lt["evals"])), script=int(pj["n_tries"])),
        pass_evals=dict(library=int(lt["evals"].sum()), script=int(pj["evals"])),
        tries=try_level(lib / "pass_tries.npz", MIX / "passes_mid" / "D4_S11" / "tries.npz"),
    )
    return dict(
        before_pass=b,
        after_pass=after,
        equal_before_pass=bool(
            b["map"]["voxels_differ"] == 0
            and b["evals_total"]["script"] == b["evals_total"]["library"]
            and b["counters_equal"]
            and b["records"]["n_differ"] == 0
        ),
    )


def main() -> None:
    assert RC.SAMPLE == "mid"
    C = RCM.RCACHE
    res: Dict[str, Any] = dict(sample="mid")
    if (C / "LibG_mid" / "run.done").exists() and (C / "Gl_mid" / "Gl.done").exists():
        res["libG_vs_Gl"] = lib_vs_gl(C)
    if (C / "LibS11_mid" / "run.done").exists():
        res["control_libS11_vs_stored_S11"] = control(C)
    (RCM.RESCUE / "lib_check.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    if "libG_vs_Gl" in res and "control_libS11_vs_stored_S11" in res:
        tables(res)
    print(json.dumps(res, indent=1, default=float), flush=True)
    os._exit(0)


def tables(res: Dict[str, Any]) -> None:
    """``rescue/tables_lib.md``: the doc table of the check (``doc_tables.write_tables``)."""
    import doc_tables as DT

    def yn(x: bool) -> str:
        return "yes" if x else "NO"

    r = [
        "| pair | map before the pass: voxels that differ | evaluations before the pass |"
        " counters (of the script run) equal | records equal | rescues: tried, rescued, tries |"
        " rescues equal | equal |",
        "|---|---|---|---|---|---|---|---|",
    ]
    g = res["libG_vs_Gl"]["before_pass"]
    c = res["control_libS11_vs_stored_S11"]["before_pass"]
    z = g["rescues"]
    for name, b, resc, eq in (
        ("LibG vs Gl", g, z, res["libG_vs_Gl"]["equal_before_pass"]),
        (
            "LibS11 vs S11 (control)",
            c,
            None,
            res["control_libS11_vs_stored_S11"]["equal_before_pass"],
        ),
    ):
        ev = b["evals_total"]
        r.append(
            f"| {name} | {b['map']['voxels_differ']} | {ev['library']:,} (library),"
            f" {ev['script']:,} (script) | {yn(b['counters_equal'])} |"
            f" {yn(b['records']['n_differ'] == 0)} | "
            + (
                f"{resc['library']}, {resc['rescued_library']}, {resc['tries_library']} (library);"
                f" {resc['script']}, {resc['rescued_script']}, {resc['tries_script']} (script) |"
                f" {yn(resc['equal'])}"
                if resc
                else "- | -"
            )
            + f" | {yn(eq)} |"
        )
    T = {
        "libcheck": "Reading line: one row per pair of runs on `mid`, before the boundary pass; "
        "the library run is compared with the script run it must equal. A counter is compared "
        "when the script run has it; the script counts its rescue evaluations inside the grid "
        "fallback's count, so that count is compared with the library's grid plus rescue "
        "evaluations.\n\n" + "\n".join(r)
    }
    a = res["libG_vs_Gl"].get("after_pass")
    lp = res["libG_vs_Gl"].get("library_pass_on_stored_map")
    if a and lp:
        m = a["map"]
        same = (
            "same tries (every field) and the same changed voxels"
            if lp["exact"] and lp.get("differ_end", 1) == 0
            else "NOT the same (lib_pass.json)"
        )
        rows = [
            "| | library run (LibG): pass in the run | script (Gl): pass after the run |"
            " library pass function on the stored map of Gl |",
            "|---|---|---|---|",
            f"| pass tries | {a['pass_tries']['library']:,} | {a['pass_tries']['script']:,} |"
            f" {lp['n_tries']['library']:,} |",
            f"| pass evaluations | {a['pass_evals']['library']:,} |"
            f" {a['pass_evals']['script']:,} | {lp['evals']['library']:,} |",
            f"| total with the pass | {a['total_with_pass']['library']:,} |"
            f" {a['total_with_pass']['script']:,} | - |",
            "| map after the pass against the script's: voxels more than 0.001 degree apart /"
            " more than 1 degree apart / largest angle |"
            f" {m['voxels_over_0_001_deg']:,} / {m['voxels_over_1_deg']:,} /"
            f" {m['max_angle_deg']:.1f} degrees | - | {same} |",
            "| map after the pass against the script's: voxels that differ in any bit (also"
            f" storage differences below 0.001 degree) | {m['voxels_differ']:,} | - | - |",
        ]
        acc = a.get("accuracy")
        if acc:
            rows += [
                f"| wrong voxels after the pass (more than {acc['wrong_threshold_deg']:g} degree"
                f" from the truth) | {acc['wrong']['library']:,} | {acc['wrong']['script']:,} |"
                " - |",
                "| by voxel, script pass as the base: fixed / broken / wrong in both |"
                f" {acc['fixed']} / {acc['broken']} / {acc['both_wrong']} | - | - |",
                f"| median error after the pass, all {acc['n_voxels']:,} voxels (degrees) |"
                f" {acc['median_error_deg']['library']:.5f} |"
                f" {acc['median_error_deg']['script']:.5f} | - |",
            ]
        ca = res["control_libS11_vs_stored_S11"]["after_pass"]
        rows.append(
            f"| control, key off: LibS11 / stored S11, pass tries |"
            f" {ca['pass_tries']['library']:,} | {ca['pass_tries']['script']:,} | - |"
        )
        rows.append(
            f"| control, key off: LibS11 / stored S11, pass evaluations |"
            f" {ca['pass_evals']['library']:,} | {ca['pass_evals']['script']:,} | - |"
        )
        T["libpass"] = (
            "Reading line: the boundary pass, by where it was made. The library run makes it "
            "inside the run, from the map in memory (unrounded angles). The script makes it "
            "after the run, from the stored map (a .mic file stores rounded angles), with the "
            "library's single generator. The last column is the library's pass function run on "
            "that same stored map.\n\n" + "\n".join(rows)
        )
        t_rows = [
            "| pair | tries in the same (voxel, candidate) order from the start | of these:"
            " ends differ | costs differ | evaluation counts differ | first try whose end"
            " differs | angle between the two ends: median / largest (degrees) | voxels with"
            " another number of tries |",
            "|---|---|---|---|---|---|---|---|",
        ]
        for name, t in (
            ("LibG vs Gl", a.get("tries")),
            ("LibS11 vs S11 (control)", ca.get("tries")),
        ):
            if not t:
                continue
            e = t["end_angle_deg"]
            t_rows.append(
                f"| {name} | {t['same_order_tries']:,} | {t['ends_differ']:,} |"
                f" {t['costs_differ']:,} | {t['evals_differ']:,} | {t['first_end_differ']} |"
                f" {e['median']:.4f} / {e['max']:.2f} | {t['voxels_with_other_try_count']:,} |"
            )
        if len(t_rows) > 2:
            T["libtries"] = (
                "Reading line: the in-run pass of the library against the script's pass, try by "
                "try, one row per pair. Both passes try the same voxels and candidates in the same "
                "order at first. Columns 3 to 7 count, over those tries, where the two differ "
                "(ends compared as float32, the script's stored type; the angle between the two "
                "ends of a try, without symmetry). The last column counts, over the whole pass, "
                "the voxels that got a different number of tries. Same input maps would give an "
                "equal first try (try 0); a shifted generator would instead keep the early tries "
                "equal.\n\n" + "\n".join(t_rows)
            )
    DT.write_tables(RCM.RESCUE / "tables_lib.md", T)


if __name__ == "__main__":
    main()
