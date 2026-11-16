"""Work item (1): grain grouping, regions and the coverage renderer against the stock render.

For a set of grains of ``grad``:
* grouping of the truth map and of the painted map (flood fill, neighbour misorientation < 2 deg)
  against the true grain map;
* at the TRUE orientations, the share of predicted coverage mass that lands on a lit pixel of the
  stock-rendered image in the predicted frame (pixel = floor of the continuous coordinate);
* the region built at the painted start over the TRUE grain's voxels (windows, touching
  components, ratio filter): spots kept against spots before the ratio filter (the objective
  itself is measured in run_grains.py, over start-map groups);
* timings of region building and of one objective evaluation.

Usage (from icenine_py/): uv run python scripts/grain_opt/check_render.py [--grains 20]
"""

import argparse
import json
import time

import common as CM  # noqa: I001  (first: threads)
import numpy as np
import torch

from icenine import grain_fit as GF


def lit_share(S: CM.Sample, R: np.ndarray, members, bary) -> dict:
    """Share of coverage mass on lit pixels (floor convention, predicted frame), per detector."""
    with torch.no_grad():
        rows = S.model.predict(torch.as_tensor(R[members]), torch.as_tensor(S.V[members]))
    stack = S.images.stack
    H, W = stack.shape[1:]
    res = {}
    for det in range(S.model.n_det):
        pc, pr, m, ok = GF._rows_points(rows, bary, det)
        c = torch.floor(pc).long()
        r = torch.floor(pr).long()
        f = rows.frame[:, None].expand_as(c)
        inb = ok[:, None] & (c >= 0) & (c < W) & (r >= 0) & (r < H)
        cc, rr, ff, mm = c[inb].numpy(), r[inb].numpy(), f[inb].numpy(), m[inb].numpy()
        lit = np.asarray(stack[det * 180 + ff, rr, cc]) > 0
        # also: lit in the frame before / after (frame mismatch diagnostics)
        litm = np.asarray(stack[det * 180 + np.clip(ff - 1, 0, 179), rr, cc]) > 0
        litp = np.asarray(stack[det * 180 + np.clip(ff + 1, 0, 179), rr, cc]) > 0
        res[f"det{det}"] = dict(
            mass_on_lit=float((mm * lit).sum() / mm.sum()),
            mass_on_lit_adjacent_frame_only=float((mm * (~lit & (litm | litp))).sum() / mm.sum()),
            n_points=int(len(mm)),
        )
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--grains", type=int, default=20)
    ap.add_argument("--n-sub", type=int, default=2)
    a = ap.parse_args()
    out = CM.OUT / "check_render"
    S = CM.Sample(out)
    T = S.T
    bary = torch.as_tensor(GF.subtriangle_bary(a.n_sub))
    rec = {"t_setup_s": S.t_setup}

    t0 = time.time()
    lab_truth = GF.group_grains(T.R, S.nb, S.sym, 2.0)
    lab_paint = GF.group_grains(S.R0, S.nb, S.sym, 2.0)
    rec["t_group_s"] = (time.time() - t0) / 2
    gm = np.asarray(T.gmap)

    def purity(lab):
        n_g = int(lab.max()) + 1
        maj = 0
        split = 0
        for g in range(n_g):
            sel = lab == g
            vals, cnt = np.unique(gm[sel], return_counts=True)
            maj += int(cnt.max())
        for tg in np.unique(gm):
            split += len(np.unique(lab[gm == tg])) > 1
        return dict(n_groups=n_g, voxels_in_majority_true_grain=maj, n=len(lab),
                    true_grains_split=int(split), true_grains=int(len(np.unique(gm))))

    rec["group_truth"] = purity(lab_truth)
    rec["group_paint"] = purity(lab_paint)
    print(json.dumps(rec, indent=1), flush=True)

    # pick grains: spread over levels and sizes, from the true grain map
    rng = np.random.default_rng(1)
    picks = []
    levels = sorted(set(T.level_of.values()))
    per = max(a.grains // len(levels), 1)
    for lv in levels:
        gids = [g for g in T.gids if T.level_of[g] == lv]
        sizes = {g: int((gm == g).sum()) for g in gids}
        gids = sorted(gids, key=lambda g: sizes[g])
        qs = np.linspace(0.1, 0.95, per)
        for q in qs:
            picks.append(gids[int(q * (len(gids) - 1))])
    picks = list(dict.fromkeys(picks))[: a.grains]
    rows = []
    for g in picks:
        members = [int(v) for v in np.nonzero(gm == g)[0]]
        shell = GF.shell_of(members, S.nb)
        t1 = time.time()
        reg_t = GF.build_region(S.model, S.images, S.R0, S.V, members, shell, bary)
        n_cand = len(GF.build_region(S.model, S.images, S.R0, S.V, members, [], bary, ratio=None).spots)
        t_reg = time.time() - t1
        ls = lit_share(S, T.R, members, bary)
        err0 = T.err(S.R0[members], members)
        row = dict(
            grain=int(g), level=T.level_of[g], n_vox=len(members), n_shell=len(shell),
            n_spots=len(reg_t.spots), n_spots_before_ratio=n_cand, obs_px_median=float(np.median([len(s.obs) for s in reg_t.spots])),
            t_region_s=t_reg,
            start_err_median=float(np.median(err0)), start_err_max=float(err0.max()),
            lit=ls,
        )
        rows.append(row)
        print(json.dumps(row), flush=True)
    rec["grains"] = rows
    (out / "check_render.json").write_text(json.dumps(rec, indent=1) + "\n")


if __name__ == "__main__":
    main()
