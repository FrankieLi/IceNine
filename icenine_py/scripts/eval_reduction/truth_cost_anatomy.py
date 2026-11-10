"""Task B (extra): why the cost at the true orientation is not zero.

For every case of cost_edges_run.py (six full-sample voxels, 30 E0 voxels), trace the complete cost
at
the truth (C++ rule) and split the loss of quality, 1 - quality = sum over peaks of (1 - q_i) / N,
by
cause: a peak whose pixels all overlap but that overlaps on only one of the two detectors (q = 1/2),
a peak with a pixel ratio below 1, a peak counted (pixel on detector) but with no qualified overlap.

  uv run python scripts/eval_reduction/truth_cost_anatomy.py --kind full|e0 --out DIR
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cost_edges as CE  # noqa: E402
import cost_edges_run as RUN  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kind", required=True, choices=["full", "e0"])
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    RUN._init(a.kind)
    rows = []
    for rule, v in [(r, v) for r in ("cpp", "cppE") for v in RUN.cases(a.kind)]:
        base, Rt, verts, ph = RUN.case_data(v)
        cf = CE.make_cf(base, "full", rule)
        tr = CE.trace_eval(cf, Rt.astype(np.float32), verts, ph)
        cnt, q = tr["cnt"], tr["q"]
        on = cnt[:, 1] > 0
        n = int(on.sum())
        loss = (1.0 - q[on]) / max(n, 1)
        one_det = on & (cnt[:, 4] == 1) & (cnt[:, 0] == cnt[:, 1])
        ratio = on & (cnt[:, 0] < cnt[:, 1]) & (cnt[:, 4] > 0)
        none = on & (cnt[:, 4] == 0)
        sg = tr["sig"]  # (P, n_det, 3, 2)
        inb = ((sg >= 0) & (sg < 2048)).all(axis=3).any(axis=2)  # (P, n_det): any vertex in bounds
        det_out = ~inb[:, 1] & inb[:, 0]
        rows.append(
            dict(
                rule=rule,
                voxel=int(v),
                n_peaks=n,
                cost=float(1.0 - tr["agg"]["quality"]),
                loss_one_detector=float(((1.0 - q) / max(n, 1))[one_det].sum()),
                loss_pixel_ratio=float(((1.0 - q) / max(n, 1))[ratio].sum()),
                loss_no_overlap=float(((1.0 - q) / max(n, 1))[none].sum()),
                n_one_detector=int(one_det.sum()),
                n_one_detector_det1_off_image=int((one_det & det_out).sum()),
                n_pixel_ratio=int(ratio.sum()),
                n_no_overlap=int(none.sum()),
            )
        )
    out = Path(a.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    (out / f"truth_cost_{a.kind}.json").write_text(json.dumps(rows, indent=1) + "\n")
    summ = {}
    for rule in ("cpp", "cppE"):
        rr = [r for r in rows if r["rule"] == rule]
        c = np.array([r["cost"] for r in rr])
        print(
            a.kind,
            rule,
            "cases",
            len(rr),
            "cost>0:",
            int((c > 1e-9).sum()),
            "median %.4f max %.4f" % (np.median(c), c.max()),
        )
        tot_cost = max(sum(r["cost"] for r in rr), 1e-12)
        summ[rule] = dict(
            cases=len(rr),
            cases_cost_above_zero=int((c > 1e-9).sum()),
            median_cost=float(np.median(c)),
            max_cost=float(c.max()),
            n_one_detector=sum(r["n_one_detector"] for r in rr),
            n_one_detector_det1_off_image=sum(r["n_one_detector_det1_off_image"] for r in rr),
            n_no_overlap=sum(r["n_no_overlap"] for r in rr),
            share_one_detector=sum(r["loss_one_detector"] for r in rr) / tot_cost,
            share_no_overlap=sum(r["loss_no_overlap"] for r in rr) / tot_cost,
        )
        print(
            "   one-detector peaks with every detector-1 vertex off the image:",
            sum(r["n_one_detector_det1_off_image"] for r in rr),
        )
        for k in ("one_detector", "pixel_ratio", "no_overlap"):
            print(
                "  ",
                k,
                "peaks",
                sum(r["n_" + k] for r in rr),
                "share of the total loss %.3f"
                % (sum(r["loss_" + k] for r in rr) / max(sum(r["cost"] for r in rr), 1e-12)),
            )
    (out / f"truth_cost_{a.kind}_summary.json").write_text(json.dumps(summ, indent=1) + "\n")


if __name__ == "__main__":
    main()
