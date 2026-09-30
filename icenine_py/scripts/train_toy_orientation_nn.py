#!/usr/bin/env python3
"""
Train a toy orientation network and evaluate it the Stage 0 way.

Heads
  offset: ToyOffsetNet, rotation-vector offset from nominal (degrees) plus a full
          covariance, trained with Gaussian negative log-likelihood.
  quat:   the v0 ToyOrientationNet, an absolute unit quaternion, trained with
          1 - |q.q_true|. Kept only to re-evaluate v0 under the Stage 0 protocol.

Evaluation (test set = fixed-magnitude offsets in random directions)
  Per magnitude bin, RMS error about the stage axis (z) and perpendicular to it (x, y),
  and the median misorientation angle, for
    - predict-nominal (offset 0),
    - the network,
    - the exact Bayes posterior mean, if --bayes is given (the best possible under
      squared error; its sqrt(trace cov) is the error floor).
  For the offset head, also its predicted per-axis sigma and the mean squared
  Mahalanobis distance of the truth (3.0 if the covariance is calibrated).

Usage:
  cd icenine_py
  uv run python scripts/generate_toy_orientation_dataset.py --n-train 1500 --test-per-bin 30
  uv run python scripts/exact_bayes_baseline.py --test scripts/toy_orientation_stage0_test.pt
  uv run python scripts/train_toy_orientation_nn.py --head offset \
      --train scripts/toy_orientation_stage0_train.pt --test scripts/toy_orientation_stage0_test.pt \
      --bayes scripts/toy_orientation_stage0_test_bayes.npz
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np
import torch


def batches(n, batch_size, rng=None):
    idx = rng.permutation(n) if rng is not None else np.arange(n)
    for a in range(0, n, batch_size):
        yield idx[a : a + batch_size]


def make_prep(meta, no_frame: bool = False):
    """Windows -> network input. Observer-rendered (Stage 1) windows are frame-coded
    uint8 and are decoded to two channels (lit, frame offset), or to the lit channel
    only with no_frame (ablation: hides the frame index); others are cast."""
    from icenine.orientation_eval import decode_windows

    if meta.get("renderer", "simulator") == "observer":
        k = int(meta["frame_half_width"])
        if no_frame:
            return (lambda x: decode_windows(x, k)[..., :1, :, :]), 1
        return (lambda x: decode_windows(x, k)), 2
    return (lambda x: x.float()), 1


def predict(model, head, windows, batch_size, prep=None, forward=None, vid=None):
    """Predicted offsets (N,3 deg) and, for the offset head, Cholesky factors (N,3,3).
    forward(x) or, for multi-voxel data, forward(x, vid_of_batch)."""
    model.eval()
    forward = forward if forward is not None else model
    means, chols = [], []
    with torch.no_grad():
        for b in batches(len(windows), batch_size):
            x = prep(windows[b]) if prep is not None else windows[b].float()
            if head == "offset":
                m, L = forward(x) if vid is None else forward(x, vid[b])
                means.append(m.cpu())
                chols.append(L.cpu())
            else:
                means.append(forward(x).cpu())
    if head == "offset":
        return torch.cat(means).numpy(), torch.cat(chols).numpy()
    return torch.cat(means).numpy(), None


def main():
    from icenine.orientation_eval import (
        error_summary,
        offsets_to_quaternions,
        quaternions_to_offsets_deg,
    )
    from icenine.orientation_nn import gaussian_nll_loss, quaternion_regression_loss
    from icenine.toy_orientation_model import PeakSetNet, ToyOffsetNet, ToyOrientationNet

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--head", choices=["offset", "quat"], default="offset")
    parser.add_argument(
        "--arch",
        choices=["fc", "set"],
        default="fc",
        help="fc: flatten-everything MLP (Stages 0-1); set: shared per-peak encoder + pooling (Stage 3)",
    )
    parser.add_argument(
        "--extra",
        nargs="*",
        default=[],
        metavar="LABEL=NPZ",
        help="extra rows for the table: npz files with pred_deg aligned to the test set",
    )
    parser.add_argument("--train", required=True)
    parser.add_argument("--test", required=True)
    parser.add_argument(
        "--bayes", default=None, help="npz from exact_bayes_baseline.py for the same test set"
    )
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--clip", type=float, default=0.0, help="gradient-norm clip (0 = off)")
    parser.add_argument(
        "--cosine", action="store_true", help="cosine learning-rate decay to 0 over --epochs"
    )
    parser.add_argument(
        "--no-meas", action="store_true", help="set net: drop the explicit measurement features"
    )
    parser.add_argument("--pool", choices=["meanmax", "meansum", "all"], default="meanmax")
    parser.add_argument(
        "--beta-nll", type=float, default=0.0, help="beta-NLL exponent (offset head); 0 = plain NLL"
    )
    parser.add_argument("--val-frac", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--device",
        choices=["cpu", "mps"],
        default="cpu",
        help="training device; mps = Apple GPU (float32 only, so evaluation stays on the CPU)",
    )
    parser.add_argument(
        "--no-frame", action="store_true", help="ablation: hide the frame channel (observer data)"
    )
    parser.add_argument("--results-json", default=None)
    parser.add_argument(
        "--save-predictions",
        default=None,
        help="npz with test predictions (and Cholesky factors) for re-evaluation",
    )
    args = parser.parse_args()

    torch.manual_seed(args.seed)
    rng = np.random.default_rng(args.seed)

    tr = torch.load(args.train)
    te = torch.load(args.test)
    n_peaks, window = tr["n_peaks"], tr["window_size"]
    R_nom = tr["R_nom"].numpy()
    windows, offsets = tr["windows"], tr["offsets_deg"].float()
    prep_cpu, in_channels = make_prep(tr, no_frame=args.no_frame)
    dev = torch.device(args.device)
    prep = lambda x: prep_cpu(x).to(dev)  # noqa: E731
    n = len(windows)
    perm = rng.permutation(n)
    n_val = max(1, int(n * args.val_frac))
    val_idx, train_idx = perm[:n_val], perm[n_val:]
    print(
        f"train {len(train_idx)} / val {n_val} samples, {n_peaks} peaks x {window}x{window}, head={args.head}"
    )

    multi, vid_train = False, None
    if args.head == "quat":
        targets = torch.from_numpy(
            offsets_to_quaternions(offsets.numpy().astype(np.float64), R_nom)
        ).float()
        model = ToyOrientationNet(n_peaks=n_peaks, window_size=window, in_channels=in_channels)
    else:
        targets = offsets
        if args.arch == "set":
            assert "context" in tr, "the set architecture needs a dataset with per-peak context"
            context = tr["context"].float()  # (M, D), or (V, M, D) for multi-voxel data
            multi = bool(tr.get("multi_voxel", False))
            vid_train = tr["voxel_id"].to(torch.long) if multi else None
            model = PeakSetNet(
                window_size=window,
                in_channels=in_channels,
                context_dim=context.shape[-1],
                use_measurements=not args.no_meas,
                frame_half_width=int(tr.get("frame_half_width", 4)),
                pool=args.pool,
            )
        else:
            model = ToyOffsetNet(n_peaks=n_peaks, window_size=window, in_channels=in_channels)
    if args.arch == "set":
        assert args.head == "offset", "the set architecture predicts an offset and covariance"
        if multi:
            forward = lambda x, v: model(x, context[v.to(dev)])  # noqa: E731
        else:
            forward = lambda x: model(x, context)  # noqa: E731
    else:
        forward = model
    n_params = sum(p.numel() for p in model.parameters())
    print(f"architecture {args.arch}: {n_params / 1e6:.2f}M parameters on {dev}", flush=True)
    model.to(dev)
    targets = targets.to(dev)
    if args.arch == "set":
        context = context.to(dev)
    opt = torch.optim.Adam(model.parameters(), lr=args.lr)
    sched = (
        torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=args.epochs) if args.cosine else None
    )

    def loss_fn(x, y, v=None):
        if args.head == "quat":
            return quaternion_regression_loss(forward(x), y)
        mean, chol = forward(x) if v is None else forward(x, v)
        return gaussian_nll_loss(mean, chol, y, beta=args.beta_nll)

    t0 = time.time()
    best_val, best_epoch, best_state = float("inf"), -1, None
    for epoch in range(args.epochs):
        model.train()
        total = 0.0
        for b in batches(len(train_idx), args.batch_size, rng):
            ib = train_idx[b]
            opt.zero_grad()
            loss = loss_fn(prep(windows[ib]), targets[ib], vid_train[ib] if multi else None)
            loss.backward()
            if args.clip > 0:
                torch.nn.utils.clip_grad_norm_(model.parameters(), args.clip)
            opt.step()
            total += loss.item() * len(ib)
        if sched is not None:
            sched.step()
        model.eval()
        with torch.no_grad():
            val = np.mean(
                [
                    loss_fn(
                        prep(windows[val_idx[b]]),
                        targets[val_idx[b]],
                        vid_train[val_idx[b]] if multi else None,
                    ).item()
                    for b in batches(n_val, args.batch_size)
                ]
            )
        if val < best_val:
            best_val, best_epoch = float(val), epoch + 1
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
        print(
            f"epoch {epoch + 1:3d}/{args.epochs}  train {total / len(train_idx):11.6f}  val {val:11.6f}  ({time.time() - t0:.0f}s)",
            flush=True,
        )
    if best_state is None:
        raise RuntimeError("validation loss was never finite; lower --lr or check the data")
    model.load_state_dict(best_state)
    print(f"restored best-validation weights from epoch {best_epoch} (val {best_val:.6f})")

    # ---- evaluation ------------------------------------------------------
    vid_test = te["voxel_id"].to(torch.long) if multi else None
    pred, chol = predict(
        model, args.head, te["windows"], args.batch_size, prep, forward, vid=vid_test
    )
    if args.head == "quat":
        pred = quaternions_to_offsets_deg(pred, R_nom)
    truth = te["offsets_deg"].numpy().astype(np.float64)
    mags = te["magnitudes_deg"].numpy()

    bayes = None
    if args.bayes:
        bz = np.load(args.bayes)
        n_b = len(bz["mean"])
        if n_b != len(truth):
            print(
                f"note: Bayes file covers the first {n_b} of {len(truth)} test cases; using those"
            )
            truth, pred, mags = truth[:n_b], pred[:n_b], mags[:n_b]
            chol = chol[:n_b] if chol is not None else None
        assert np.allclose(
            bz["offsets_deg"], truth[:n_b], atol=1e-5
        ), "Bayes file is for a different test set"
        bayes = bz

    extras = {}
    for item in args.extra:
        label, _, path = item.partition("=")
        ez = np.load(path)
        assert np.allclose(
            ez["truth_deg"], truth[: len(ez["truth_deg"])], atol=1e-5
        ), f"{path} is for a different test set"
        extras[label] = ez["pred_deg"]

    if args.save_predictions:
        np.savez(
            args.save_predictions,
            pred_deg=pred,
            truth_deg=truth,
            magnitudes_deg=mags,
            chol=chol if chol is not None else np.zeros(0),
            best_epoch=best_epoch,
        )
    rows = {}
    print(
        f"\n{'|delta|':>8} {'method':<20} {'n':>3} {'rms_z':>10} {'rms_perp':>10} {'median_ang':>11} {'<0.5deg':>8} {'<0.1deg':>8}   (degrees)"
    )
    if multi:
        held = te["held_out"].numpy()[vid_test.numpy()]
        groups = [("in-dist/", ~held), ("held-out/", held)]
    else:
        groups = [("", np.ones(len(mags), dtype=bool))]
    for tag, gmask in groups:
        if tag:
            print(f"\n=== {tag.rstrip('/')} voxels ({int(gmask.sum())} test cases) ===")
        for mag in sorted(set(mags.tolist())):
            m = (mags == mag) & gmask
            entries = [
                ("predict-nominal", error_summary(np.zeros((m.sum(), 3)), truth[m])),
                (f"net ({args.arch}/{args.head})", error_summary(pred[m], truth[m])),
            ]
            for label, ep in extras.items():
                entries.insert(
                    1, (label, error_summary(ep[m[: len(ep)]], truth[: len(ep)][m[: len(ep)]]))
                )
            if bayes is not None:
                mb = m & np.isfinite(bayes["mean"]).all(
                    axis=1
                )  # cases where the sampler found no members are excluded
                if mb.sum() < m.sum():
                    print(
                        f"note: exact Bayes failed on {m.sum() - mb.sum()} case(s) at |delta| = {mag}; excluded from its row"
                    )
                entries.append(("exact Bayes", error_summary(bayes["mean"][mb], truth[mb])))
            for name, s in entries:
                print(
                    f"{mag:8.2f} {name:<20} {s['n']:3d} {s['rms_z']:10.5f} {s['rms_perp']:10.5f} {s['median_angle']:11.5f} {s['success_0p5']:8.0%} {s['success_0p1']:8.0%}"
                )
            extra = {}
            if bayes is not None:
                floor = float(np.nanmean(np.sqrt(np.trace(bayes["cov"][m], axis1=1, axis2=2))))
                extra["bayes_floor_sqrt_trace"] = floor
                print(f"{'':8} {'Bayes floor':<14} {'':>3} sqrt(tr cov) = {floor:.5f}")
            if chol is not None:
                cov = chol[m] @ chol[m].transpose(0, 2, 1)
                sd = np.sqrt(np.diagonal(cov, axis1=1, axis2=2)).mean(axis=0)
                r = (truth[m] - pred[m])[:, :, None]
                z = np.linalg.solve(chol[m], r)[:, :, 0]
                maha = float((z**2).sum(axis=1).mean())
                extra.update(pred_sigma_xyz=sd.tolist(), mean_mahalanobis_sq=maha)
                print(
                    f"{'':8} {'net sigma xyz':<14} {'':>3} {np.round(sd, 5)}   mean Mahalanobis^2 = {maha:.2f} (3.0 if calibrated)"
                )
            rows[tag + str(mag)] = {name: s for name, s in entries} | extra

    if multi:
        r_perp, n_pk = te["r_perp_um"].numpy(), te["n_peaks_per_voxel"].numpy()
        per_voxel = []
        print(
            f"\n{'voxel':>5} {'r_perp':>7} {'peaks':>5} {'held':>5} {'n':>4}  median angle (deg), all magnitudes"
        )
        for v in range(len(r_perp)):
            mv = vid_test.numpy() == v
            ang = {"net": error_summary(pred[mv], truth[mv])["median_angle"]}
            for label, ep in extras.items():
                ang[label] = error_summary(ep[mv], truth[mv])["median_angle"]
            ang1 = {
                "net": error_summary(pred[mv & (mags == 1.0)], truth[mv & (mags == 1.0)])[
                    "median_angle"
                ]
            }
            per_voxel.append(
                dict(
                    voxel=v,
                    r_perp_um=float(r_perp[v]),
                    n_peaks=int(n_pk[v]),
                    held_out=bool(te["held_out"][v]),
                    median_angle=ang,
                    net_median_angle_1deg=ang1["net"],
                )
            )
            print(
                f"{v:5d} {r_perp[v]:7.0f} {n_pk[v]:5d} {str(bool(te['held_out'][v])):>5} {mv.sum():4d}  "
                + "  ".join(f"{k} {x:.4f}" for k, x in ang.items())
            )
        rows["per_voxel"] = per_voxel

    if args.results_json:
        Path(args.results_json).write_text(json.dumps(rows, indent=2))
        print(f"saved {args.results_json}")


if __name__ == "__main__":
    main()
