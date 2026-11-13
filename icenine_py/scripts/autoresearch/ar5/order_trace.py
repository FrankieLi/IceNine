"""AR5 step 8, review fixes 2 and 3: where the two orders part, and where the flips are made.

  RESEED_SAMPLE=<held|mid|grad|big> uv run python scripts/autoresearch/ar5/order_trace.py
  uv run python scripts/autoresearch/ar5/order_trace.py --tables

Written after the runs, from the finished run files only (no run, nothing tuned; the MC traces are
not read). Per sample it compares the arms ``A`` (present order) and ``Y`` (sweep order):

* the divergence chain: how many seed steps, from the start, visit the same voxels in the same
  order; the first step (at the same voxel, in both logs) whose taken end differs; how many voxels
  differ between the two maps before and after the pass (misorientation above 1e-6 degree);
* the seed steps: matched by (voxel, occurrence), run in one arm only, and the wrong ends (more
  than 1 degree from the truth) of all steps and of the one-arm steps;
* every voxel that is wrong in exactly one arm after the pass ("flip"; fixed = wrong in A, right in
  Y): whether the same flip was already there before the pass, what filled the voxel in each arm
  (the run's record ``source``), the seed steps run at the voxel in each arm, and the region of
  changed voxels around it (connected through the nearest-neighbour voxels).

This is a description of the two paths, not a test of cause: a flip inside a region that holds a
differing seed step is near that step, which does not show the step caused it.

Writes ``benchmarks/autoresearch/ar5/order/trace_<sample>.json``; ``--tables`` reads the four and
writes ``trace_tables.md`` (blocks for ``scripts/dev/sync_doc_tables.py``), with the draw-noise
reference of step 7 (``ar5/gate/flips.json``) and ``Yl`` against ``Y`` on ``big``.
"""

import gzip
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
ORDER = ICE / "benchmarks" / "autoresearch" / "ar5" / "order"
GATE = ICE / "benchmarks" / "autoresearch" / "ar5" / "gate"
SAMPLES = ["held", "mid", "grad", "big"]
GATE_KEY = dict(held="held S11v2", mid="mid S11n", grad="grad S11n", big="big S11")
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICE / "scripts" / "common"))

import order_logs as OL  # noqa: E402

MIS_SAME = 1e-6  # degrees: two maps "differ" at a voxel above this misorientation


def _load_gz(p: Path) -> Any:
    with gzip.open(p, "rt") as f:
        return json.load(f)


def compute() -> Dict[str, Any]:
    import numpy as np
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree

    import order_common as OC  # first of the sample imports (held guard)

    import gradient_common as G

    s = os.environ["RESEED_SAMPLE"]
    assert s in SAMPLES, s
    C = OC.ORDER_CACHE
    T = G.Truth()
    n = T.n
    logs = {a: _load_gz(C / f"{a}_{s}" / "e1_log.json.gz") for a in ("A", "Y")}
    recs = {a: _load_gz(C / f"{a}_{s}" / "records.json.gz") for a in ("A", "Y")}
    errs = {a: dict(np.load(C / "stats" / f"part2_{a}_{s}_err.npz")) for a in ("A", "Y")}
    maps = {
        (a, tag): np.load(C / f"passes_{s}" / f"{tag}_{a}.npy").reshape(n, 3, 3)
        for a in ("A", "Y")
        for tag in ("before", "after")
    }

    def end_err(v: int, R: List[float]) -> float:
        return float(T.err(np.asarray(R, np.float64).reshape(1, 3, 3), [v])[0])

    la, ly = logs["A"], logs["Y"]
    cmp = OL.compare_logs(la, ly, lambda v, R: end_err(v, R) > G.WRONG)
    ia, iy = OL.index(la), OL.index(ly)
    one = {"A": set(ia) - set(iy), "Y": set(iy) - set(ia)}
    wrong_all: Dict[str, int] = {}
    wrong_one: Dict[str, int] = {}
    for a, idx in (("A", ia), ("Y", iy)):
        wa = wo = 0
        for k, e in idx.items():
            if not e["fits"]:
                continue
            w = end_err(k[0], OL.taken_end(e)["R"]) > G.WRONG
            wa += int(w)
            wo += int(w and k in one[a])
        wrong_all[a], wrong_one[a] = wa, wo

    mis_b = T.mis(maps[("A", "before")], maps[("Y", "before")])
    mis_a = T.mis(maps[("A", "after")], maps[("Y", "after")])
    changed = mis_b > MIS_SAME

    # regions of changed voxels (before the pass), linked through nearest-neighbour voxels
    tree = cKDTree(T.xy)
    dist, _ = tree.query(T.xy, k=7)
    d1 = float(np.min(dist[:, 1:][dist[:, 1:] > 1e-9]))
    radius = 1.02 * d1
    pairs = tree.query_pairs(radius, output_type="ndarray")
    keep = changed[pairs[:, 0]] & changed[pairs[:, 1]]
    p = pairs[keep]
    g = coo_matrix((np.ones(len(p)), (p[:, 0], p[:, 1])), shape=(n, n))
    _, lab = connected_components(g, directed=False)
    size: Dict[int, int] = {}
    for v in np.where(changed)[0]:
        size[int(lab[v])] = size.get(int(lab[v]), 0) + 1
    diff_seed_vox = set()
    for k in set(ia) & set(iy):
        ea, ey = ia[k], iy[k]
        if ea["fits"] and ey["fits"]:
            ta, ty = OL.taken_end(ea), OL.taken_end(ey)
            if (ta["j"], ta["kind"]) != (ty["j"], ty["kind"]):
                diff_seed_vox.add(k[0])
    one_vox = {k[0] for a in one for k in one[a]}
    comp_diff = {int(lab[v]) for v in diff_seed_vox if changed[v]}
    comp_one = {int(lab[v]) for v in one_vox if changed[v]}

    def seed_rows(a: str, v: int) -> List[Dict[str, Any]]:
        out = []
        idx = ia if a == "A" else iy
        for k, e in idx.items():
            if k[0] != v or not e["fits"]:
                continue
            t = OL.taken_end(e)
            ends = OL.kept_ends(e)
            best = max(ends, key=lambda x: x["hit"])
            out.append(
                dict(
                    occurrence=k[1],
                    one_arm_only=k in one[a],
                    reseed=bool(e.get("reseed", False)),
                    stop=e["stop"],
                    taken=dict(
                        j=t["j"], kind=t["kind"], hit=round(t["hit"], 3), err=end_err(v, t["R"])
                    ),
                    best_hit=dict(
                        j=best["j"],
                        kind=best["kind"],
                        hit=round(best["hit"], 3),
                        err=end_err(v, best["R"]),
                    ),
                )
            )
        return out

    wA_b, wY_b = errs["A"]["before"] > G.WRONG, errs["Y"]["before"] > G.WRONG
    wA_a, wY_a = errs["A"]["after"] > G.WRONG, errs["Y"]["after"] > G.WRONG
    flips = []
    for v in np.where(wA_a != wY_a)[0]:
        v = int(v)
        if not changed[v]:
            region = "identical before the pass"
        elif int(lab[v]) in comp_diff:
            region = "changed region with a matched seed step whose taken end differs"
        elif int(lab[v]) in comp_one:
            region = "changed region with one-arm seed steps only"
        else:
            region = "changed region without a seed step"
        flips.append(
            dict(
                voxel=v,
                kind="fixed" if wA_a[v] else "broken",
                made=OL.flip_origin(bool(wA_b[v]), bool(wY_b[v]), bool(wA_a[v]), bool(wY_a[v])),
                source_A=recs["A"][str(v)]["source"],
                source_Y=recs["Y"][str(v)]["source"],
                err_before=[round(float(errs[a]["before"][v]), 3) for a in ("A", "Y")],
                err_after=[round(float(errs[a]["after"][v]), 3) for a in ("A", "Y")],
                region=region,
                region_size=size.get(int(lab[v]), 0) if changed[v] else 0,
                seed_steps_A=seed_rows("A", v),
                seed_steps_Y=seed_rows("Y", v),
            )
        )
    first = OL.first_taken_difference(la, ly)
    out = dict(
        sample=s,
        steps_A=len(la),
        steps_Y=len(ly),
        matched=cmp["matched"],
        only_A=cmp["only_a"],
        only_Y=cmp["only_b"],
        taken_differs_matched=cmp["taken_differs"],
        differing_wrong=[cmp["differing_wrong_first"], cmp["differing_wrong_second"]],
        wrong_ends_all=[wrong_all["A"], wrong_all["Y"]],
        wrong_ends_one_arm=[wrong_one["A"], wrong_one["Y"]],
        same_voxel_prefix=OL.same_voxel_prefix(la, ly),
        first_taken_difference_index0=first,
        first_taken_difference_step=first + 1 if first >= 0 else None,
        diff_before_voxels=int(changed.sum()),
        diff_after_voxels=int((mis_a > MIS_SAME).sum()),
        changed_regions=len(size),
        largest_changed_region=max(size.values()) if size else 0,
        neighbour_radius=radius,
        flips_fixed=sum(1 for f in flips if f["kind"] == "fixed"),
        flips_broken=sum(1 for f in flips if f["kind"] == "broken"),
        flips_made_before=sum(1 for f in flips if f["made"] == "before the pass"),
        flips_made_in_pass=sum(1 for f in flips if f["made"] == "in the pass"),
        flips=flips,
    )
    (ORDER / f"trace_{s}.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


def _seed_note(f: Dict[str, Any]) -> str:
    parts = []
    for a in ("A", "Y"):
        for r in f[f"seed_steps_{a}"]:
            t = r["taken"]
            tag = f"{a}-only " if r["one_arm_only"] else ""
            stop = "an end reached 0.95" if r["stop"] >= 0 else "no end reached 0.95"
            parts.append(
                f"{tag}seed step ({stop}): took hit {t['hit']:.3f} at {t['err']:.2f} deg"
                + (
                    f" (highest hit {r['best_hit']['hit']:.3f} at {r['best_hit']['err']:.2f} deg)"
                    if r["best_hit"]["j"] != t["j"] or r["best_hit"]["kind"] != t["kind"]
                    else ""
                )
            )
    return "; ".join(parts) if parts else "no seed step at the voxel"


def tables() -> None:
    import doc_tables as DT

    tr = {s: json.loads((ORDER / f"trace_{s}.json").read_text()) for s in SAMPLES}
    big = json.loads((ORDER / "stats_big.json").read_text())
    gate = json.loads((GATE / "flips.json").read_text())
    stats = {s: json.loads((ORDER / f"stats_{s}.json").read_text()) for s in SAMPLES}
    T: Dict[str, str] = {}
    T["order_trace_chain"] = DT.markdown_table(
        [
            {
                "sample": s,
                "seed steps at the same voxels from the start": t["same_voxel_prefix"],
                "first step whose taken end differs (counting from 1)": t[
                    "first_taken_difference_step"
                ],
                "voxels that differ before the pass": f"{t['diff_before_voxels']:,}",
                "voxels that differ after the pass": f"{t['diff_after_voxels']:,}",
                "matched seed steps": f"{t['matched']:,}",
                "of those, taken end differs": t["taken_differs_matched"],
                "seed steps in one arm only (A / Y)": f"{t['only_A']} / {t['only_Y']}",
                "wrong ends, one-arm steps (A / Y)": f"{t['wrong_ends_one_arm'][0]} / "
                f"{t['wrong_ends_one_arm'][1]}",
                "wrong ends, all seed steps (A / Y)": f"{t['wrong_ends_all'][0]} / "
                f"{t['wrong_ends_all'][1]}",
            }
            for s, t in tr.items()
        ],
        [
            "sample",
            "seed steps at the same voxels from the start",
            "first step whose taken end differs (counting from 1)",
            "voxels that differ before the pass",
            "voxels that differ after the pass",
            "matched seed steps",
            "of those, taken end differs",
            "seed steps in one arm only (A / Y)",
            "wrong ends, one-arm steps (A / Y)",
            "wrong ends, all seed steps (A / Y)",
        ],
    )
    T["order_trace_origin"] = DT.markdown_table(
        [
            {
                "sample": s,
                "flips after the pass, Y against A (fixed / broken)": f"{t['flips_fixed']} / "
                f"{t['flips_broken']}",
                "made before the pass": t["flips_made_before"],
                "made in the pass": t["flips_made_in_pass"],
            }
            for s, t in tr.items()
        ],
        [
            "sample",
            "flips after the pass, Y against A (fixed / broken)",
            "made before the pass",
            "made in the pass",
        ],
    )
    rows = []
    for s in ("held", "mid"):
        for f in tr[s]["flips"]:
            rows.append(
                {
                    "sample": s,
                    "voxel": f["voxel"],
                    "Y against A": f["kind"],
                    "made": f["made"],
                    "filled by, A / Y": f"{f['source_A']} / {f['source_Y']}",
                    "error before the pass, A / Y (deg)": f"{f['err_before'][0]} / "
                    f"{f['err_before'][1]}",
                    "error after the pass, A / Y (deg)": f"{f['err_after'][0]} / "
                    f"{f['err_after'][1]}",
                    "seed steps at the voxel": _seed_note(f),
                    "region before the pass": f["region"]
                    + (
                        f" ({f['region_size']} voxel{'s' if f['region_size'] > 1 else ''})"
                        if f["region_size"]
                        else ""
                    ),
                }
            )
    T["order_trace_flips"] = DT.markdown_table(
        rows,
        [
            "sample",
            "voxel",
            "Y against A",
            "made",
            "filled by, A / Y",
            "error before the pass, A / Y (deg)",
            "error after the pass, A / Y (deg)",
            "seed steps at the voxel",
            "region before the pass",
        ],
    )
    c_ya = "Y against A, after the pass (this step; order changed, draws aligned)"
    c_g7 = "one re-draw of the kept pass tries against the full pass (step 7, variant 1)"
    c_yl = "Yl against Y, after the pass (this step; only the draws changed)"
    nz = []
    for s in SAMPLES:
        va = stats[s]["pairs"]["Y_vs_A"]["voxels"]["after"]
        g = gate[GATE_KEY[s]]
        r = {"sample": s, c_ya: f"{va['fixed']} / {va['broken']}", c_yl: ""}
        r[c_g7] = f"{g['fixed_vs_full']} / {g['broken_vs_full']}"
        if s == "big":
            yl = big["pairs"]["Yl_vs_Y"]["voxels"]["after"]
            r[c_yl] = f"{yl['fixed']} / {yl['broken']}"
        nz.append(r)
    T["order_noise"] = DT.markdown_table(nz, ["sample", c_ya, c_g7, c_yl])
    DT.write_tables(ORDER / "trace_tables.md", T)


def main() -> None:
    if "--tables" in sys.argv[1:]:
        tables()
        return
    out = compute()
    print(json.dumps({k: v for k, v in out.items() if k != "flips"}, indent=1))


if __name__ == "__main__":
    main()
