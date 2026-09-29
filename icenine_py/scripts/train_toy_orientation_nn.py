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


def predict(model, head, windows, batch_size):
    """Predicted offsets (N,3 deg) and, for the offset head, Cholesky factors (N,3,3)."""
    model.eval()
    means, chols = [], []
    with torch.no_grad():
        for b in batches(len(windows), batch_size):
            x = windows[b].float()
            if head == "offset":
                m, L = model(x)
                means.append(m)
                chols.append(L)
            else:
                means.append(model(x))
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
    from icenine.toy_orientation_model import ToyOffsetNet, ToyOrientationNet

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--head", choices=["offset", "quat"], default="offset")
    parser.add_argument("--train", required=True)
    parser.add_argument("--test", required=True)
    parser.add_argument("--bayes", default=None, help="npz from exact_bayes_baseline.py for the same test set")
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--beta-nll", type=float, default=0.0, help="beta-NLL exponent (offset head); 0 = plain NLL")
    parser.add_argument("--val-frac", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--results-json", default=None)
    parser.add_argument("--save-predictions", default=None, help="npz with test predictions (and Cholesky factors) for re-evaluation")
    args = parser.parse_args()

    torch.manual_seed(args.seed)
    rng = np.random.default_rng(args.seed)

    tr = torch.load(args.train)
    te = torch.load(args.test)
    n_peaks, window = tr["n_peaks"], tr["window_size"]
    R_nom = tr["R_nom"].numpy()
    windows, offsets = tr["windows"], tr["offsets_deg"].float()
    n = len(windows)
    perm = rng.permutation(n)
    n_val = max(1, int(n * args.val_frac))
    val_idx, train_idx = perm[:n_val], perm[n_val:]
    print(f"train {len(train_idx)} / val {n_val} samples, {n_peaks} peaks x {window}x{window}, head={args.head}")

    if args.head == "quat":
        targets = torch.from_numpy(offsets_to_quaternions(offsets.numpy().astype(np.float64), R_nom)).float()
        model = ToyOrientationNet(n_peaks=n_peaks, window_size=window)
    else:
        targets = offsets
        model = ToyOffsetNet(n_peaks=n_peaks, window_size=window)
    opt = torch.optim.Adam(model.parameters(), lr=args.lr)

    def loss_fn(x, y):
        if args.head == "quat":
            return quaternion_regression_loss(model(x), y)
        mean, chol = model(x)
        return gaussian_nll_loss(mean, chol, y, beta=args.beta_nll)

    t0 = time.time()
    best_val, best_epoch, best_state = float("inf"), -1, None
    for epoch in range(args.epochs):
        model.train()
        total = 0.0
        for b in batches(len(train_idx), args.batch_size, rng):
            ib = train_idx[b]
            opt.zero_grad()
            loss = loss_fn(windows[ib].float(), targets[ib])
            loss.backward()
            opt.step()
            total += loss.item() * len(ib)
        model.eval()
        with torch.no_grad():
            val = np.mean([loss_fn(windows[val_idx[b]].float(), targets[val_idx[b]]).item() for b in batches(n_val, args.batch_size)])
        if val < best_val:
            best_val, best_epoch = float(val), epoch + 1
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
        if epoch % max(1, args.epochs // 10) == 0 or epoch == args.epochs - 1:
            print(f"epoch {epoch + 1:3d}/{args.epochs}  train {total / len(train_idx):11.6f}  val {val:11.6f}  ({time.time() - t0:.0f}s)")
    model.load_state_dict(best_state)
    print(f"restored best-validation weights from epoch {best_epoch} (val {best_val:.6f})")

    # ---- evaluation ------------------------------------------------------
    pred, chol = predict(model, args.head, te["windows"], args.batch_size)
    if args.head == "quat":
        pred = quaternions_to_offsets_deg(pred, R_nom)
    truth = te["offsets_deg"].numpy().astype(np.float64)
    mags = te["magnitudes_deg"].numpy()

    bayes = None
    if args.bayes:
        bz = np.load(args.bayes)
        n_b = len(bz["mean"])
        if n_b != len(truth):
            print(f"note: Bayes file covers the first {n_b} of {len(truth)} test cases; using those")
            truth, pred, mags = truth[:n_b], pred[:n_b], mags[:n_b]
            chol = chol[:n_b] if chol is not None else None
        bayes = bz

    if args.save_predictions:
        np.savez(args.save_predictions, pred_deg=pred, truth_deg=truth, magnitudes_deg=mags,
                 chol=chol if chol is not None else np.zeros(0), best_epoch=best_epoch)
    rows = {}
    print(f"\n{'|delta|':>8} {'method':<14} {'n':>3} {'rms_z':>10} {'rms_perp':>10} {'median_ang':>11} {'<0.5deg':>8} {'<0.1deg':>8}   (degrees)")
    for mag in sorted(set(mags.tolist())):
        m = mags == mag
        entries = [
            ("predict-nominal", error_summary(np.zeros((m.sum(), 3)), truth[m])),
            (f"net ({args.head})", error_summary(pred[m], truth[m])),
        ]
        if bayes is not None:
            entries.append(("exact Bayes", error_summary(bayes["mean"][m], truth[m])))
        for name, s in entries:
            print(f"{mag:8.2f} {name:<14} {s['n']:3d} {s['rms_z']:10.5f} {s['rms_perp']:10.5f} {s['median_angle']:11.5f} {s['success_0p5']:8.0%} {s['success_0p1']:8.0%}")
        extra = {}
        if bayes is not None:
            floor = float(np.mean(np.sqrt(np.trace(bayes["cov"][m], axis1=1, axis2=2))))
            extra["bayes_floor_sqrt_trace"] = floor
            print(f"{'':8} {'Bayes floor':<14} {'':>3} sqrt(tr cov) = {floor:.5f}")
        if chol is not None:
            cov = chol[m] @ chol[m].transpose(0, 2, 1)
            sd = np.sqrt(np.diagonal(cov, axis1=1, axis2=2)).mean(axis=0)
            r = (truth[m] - pred[m])[:, :, None]
            z = np.linalg.solve(chol[m], r)[:, :, 0]
            maha = float((z**2).sum(axis=1).mean())
            extra.update(pred_sigma_xyz=sd.tolist(), mean_mahalanobis_sq=maha)
            print(f"{'':8} {'net sigma xyz':<14} {'':>3} {np.round(sd, 5)}   mean Mahalanobis^2 = {maha:.2f} (3.0 if calibrated)")
        rows[str(mag)] = {name: s for name, s in entries} | extra

    if args.results_json:
        Path(args.results_json).write_text(json.dumps(rows, indent=2))
        print(f"saved {args.results_json}")


if __name__ == "__main__":
    main()
