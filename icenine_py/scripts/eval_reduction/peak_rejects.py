"""Task B step 3: the real peaks the C++ and ray rules reject at the truth on clean images.

The ray-overlap study counted, over its 30 voxels, 524 detector-0 overlaps at the true orientation
(screen: pixel radius 3, Q-max 5, centroid) and found that the C++ rule and both ray windows reject
3
of them (cause not checked). This script finds those peaks with the per-peak trace of cost_edges.py
(lenient arm accepts, C++ arm rejects) and traces each: predicted pixels on both detectors, the lit
pixels around them (distance to the nearest lit pixel per detector), the complete cost's verdict and
the same peak under CppEdgeRules.

  uv run python scripts/eval_reduction/peak_rejects.py --out benchmarks/eval_reduction/B/step3
"""

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "scripts" / "phase_d"))
import cost_edges as CE  # noqa: E402
import parity_common as PC  # noqa: E402

MAXD = 40  # search window (pixels) for the nearest lit pixel


def nearest_lit(img: np.ndarray, col: int, row: int) -> Dict[str, Any]:
    """Distance (Chebyshev, pixels) to the nearest lit pixel and its position, within MAXD."""
    nr, nc = img.shape
    r0, r1 = max(0, row - MAXD), min(nr, row + MAXD + 1)
    c0, c1 = max(0, col - MAXD), min(nc, col + MAXD + 1)
    if r0 >= r1 or c0 >= c1:
        return dict(dist=None, at=None, n_in_5=0)
    sub = img[r0:r1, c0:c1]
    ys, xs = np.nonzero(sub)
    if len(ys) == 0:
        return dict(dist=None, at=None, n_in_5=0)
    d = np.maximum(np.abs(ys + r0 - row), np.abs(xs + c0 - col))
    i = int(d.argmin())
    w = img[max(0, row - 2) : row + 3, max(0, col - 2) : col + 3]
    return dict(dist=int(d[i]), at=[int(xs[i] + c0), int(ys[i] + r0)], n_in_5=int(w.sum()))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    import ray_overlap_study as RS

    ctx = PC.PyContext()
    s = ctx.setup
    base = SimpleNamespace(
        simulator=s.simulator,
        detector_list=s.detector_list,
        range_map=s.range_map,
        exp_data=s.exp_data,
        sample=s.sample,
        structure_list=s.structure_list,
        eta_limit=ctx.eta_limit,
        reg_mm=float(s.exp_setup.config_file.detector_dist_deviation),
    )
    voxels = RS.select_voxels()
    ph0 = int(ctx.mic.voxels[voxels[0]].phase)
    cfs = {
        (st, r): CE.make_cf(base, st, r)
        for st in ("screen", "full")
        for r in ("old", "cpp", "cppE", "r2")
    }
    gs = cfs[("screen", "old")]._phase_recip_vecs[ph0][0].numpy()
    gf = cfs[("full", "old")]._phase_recip_vecs[ph0][0].numpy()
    # the screen's reflection list (Q <= 5) is a subset of the full list: map screen index -> full
    to_full = np.array([int(np.abs(gf - g).sum(1).argmin()) for g in gs])
    assert all(np.allclose(gf[to_full[i]], gs[i]) for i in range(len(gs)))
    edge_stats = [
        dict(
            n_lenient_with_centre_outside=0,
            n_of_those_rejected_by_cpp=0,
            n_lenient_centre_outside_box_reaches_in=0,
            n_of_those_rejected=0,
        )
        for _ in range(2)
    ]
    tot = {k: 0 for k in cfs}
    rej: Dict[Any, int] = {k: 0 for k in cfs}
    rows: List[Dict[str, Any]] = []
    exact_rows: List[Dict[str, Any]] = []
    for v in voxels:
        verts, ph, Rt = ctx.verts(v), ctx.mic.voxels[v].phase, ctx.truth(v).astype(np.float32)
        tr = {k: CE.trace_eval(cf, Rt, verts, ph) for k, cf in cfs.items()}
        for k in cfs:
            tot[k] += int((tr[k]["cnt"][:, 2] > 0).sum())
        for k in cfs:
            if k[1] == "old":
                continue
            lenient = tr[(k[0], "old")]["cnt"][:, 2] > 0
            strict = tr[k]["cnt"][:, 2] > 0
            rej[k] += int((lenient & ~strict).sum())
        lost = (tr[("full", "cpp")]["cnt"][:, 2] > 0) & (tr[("full", "cppE")]["cnt"][:, 2] == 0)
        for key in np.nonzero(lost)[0]:
            ta, tb = tr[("full", "cpp")], tr[("full", "cppE")]
            exact_rows.append(
                dict(
                    voxel=int(v),
                    key=int(key),
                    same_frame=bool(ta["frame"][key] == tb["frame"][key]),
                    same_vertex_pixels=bool((ta["sig"][key] == tb["sig"][key]).all()),
                    detector1_vertex_pixels=tb["sig"][key, 1].tolist(),
                    counts_cpp=ta["cnt"][key].tolist(),
                    counts_cppE=tb["cnt"][key].tolist(),
                )
            )
        lenient = tr[("screen", "old")]["cnt"][:, 2] > 0
        t_c = tr[("screen", "cpp")]
        sg = t_c["sig"][:, :, 0, :]  # (P, n_det, 2) centre (col, row)
        outside = ((sg < 0) | (sg >= 2048)).any(axis=2)  # (P, n_det)
        reach = outside & ((sg >= -3) & (sg < 2048 + 3)).all(axis=2)  # box still touches the image
        for d in range(2):
            edge_stats[d]["n_lenient_centre_outside_box_reaches_in"] += int(
                (lenient & reach[:, d]).sum()
            )
            edge_stats[d]["n_of_those_rejected"] += int(
                (lenient & reach[:, d] & ~(t_c["cnt"][:, 2] > 0)).sum()
            )
            m = lenient & outside[:, d]
            edge_stats[d]["n_lenient_with_centre_outside"] += int(m.sum())
            edge_stats[d]["n_of_those_rejected_by_cpp"] += int((m & ~(t_c["cnt"][:, 2] > 0)).sum())
        for key in np.nonzero(lenient & ~(tr[("screen", "cpp")]["cnt"][:, 2] > 0))[0]:
            t = tr[("screen", "cpp")]
            fr = int(t["frame"][key])
            imgs = [ctx.setup.exp_data.get_image(fr, d).get_binary_numpy() for d in range(2)]
            per_det = []
            for d in range(2):
                c, r = (int(x) for x in t["sig"][key, d, 0])
                per_det.append(dict(center=[c, r], **nearest_lit(imgs[d], c, r)))
            fk = int(2 * to_full[key // 2] + key % 2)
            full = {
                f"full|{r}": dict(
                    pkov=int(tr[("full", r)]["cnt"][fk, 2]),
                    pkon=int(tr[("full", r)]["cnt"][fk, 3]),
                    pixov=int(tr[("full", r)]["cnt"][fk, 0]),
                    pixon=int(tr[("full", r)]["cnt"][fk, 1]),
                )
                for r in ("old", "cpp", "cppE", "r2")
            }
            scr = {
                f"screen|{r}": int(tr[("screen", r)]["cnt"][key, 2])
                for r in ("old", "cpp", "cppE", "r2")
            }
            # full-cost per-detector vertices
            vt = tr[("full", "cpp")]["sig"][fk]  # (n_det, 3, 2) floor of vertex pixel coordinates
            rows.append(
                dict(
                    voxel=int(v),
                    key=int(key),
                    refl=int(key // 2),
                    branch=int(key % 2),
                    omega_deg=float(np.degrees(t["omega"][key])),
                    frame=fr,
                    screen_detectors=per_det,
                    screen_peak_overlap=scr,
                    full=full,
                    full_vertex_pixels=vt.tolist(),
                    cnt_screen=tr[("screen", "cpp")]["cnt"][key].tolist(),
                )
            )
    summary = dict(
        n_voxels=len(voxels),
        lenient_overlaps_at_truth={f"{k[0]}|{k[1]}": v for k, v in tot.items()},
        screen_centre_outside_the_image_by_detector=edge_stats,
        rejected_vs_lenient={f"{k[0]}|{k[1]}": v for k, v in rej.items() if k[1] != "old"},
        complete_cost_peaks_lost_with_cpp_edge_rules=exact_rows,
        peaks=rows,
    )
    (out / "summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps({k: v for k, v in summary.items() if k not in ("peaks",)}, indent=1)[:3000])
    for r in rows:
        print(json.dumps(r)[:700])


if __name__ == "__main__":
    main()
