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
from typing import Optional

import numpy as np
import torch


def variant_path(save: str, variant: str) -> str:
    """Prediction file for a test-set variant: the clean variant keeps `save` as given,
    the others get `_<variant>` inserted before the suffix (default .npz)."""
    if variant == "clean":
        return save
    p = Path(save)
    return str(p.with_name(f"{p.stem}_{variant}{p.suffix or '.npz'}"))


def padding_mask(data: dict) -> Optional[torch.Tensor]:
    """(N, M) bool: entry m of sample n is a real peak of its voxel (not zero padding), for
    multi-voxel data; None for single-voxel data (no padding)."""
    npk = data.get("n_peaks_per_voxel")
    if npk is None:
        return None
    vid = data["voxel_id"].long()
    return torch.arange(int(data["n_peaks"]))[None, :] < npk.long()[vid][:, None]


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
        RealismConfig,
        make_realistic_dataset,
        make_realistic_windows,
        error_summary,
        offsets_to_quaternions,
        quaternions_to_offsets_deg,
    )
    from icenine.orientation_nn import (
        decoupled_nll_loss,
        ema_update,
        select_checkpoint_state,
        gaussian_nll_loss,
        mse_deg_loss,
        quaternion_regression_loss,
        split_by_voxel,
    )
    from icenine.toy_orientation_model import (
        FrameProbeNet,
        GNLayerNet,
        PeakSetNet,
        ToyOffsetNet,
        ToyOrientationNet,
    )

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--head", choices=["offset", "quat"], default="offset")
    parser.add_argument(
        "--arch",
        choices=["fc", "set", "probe", "gn"],
        default="fc",
        help="fc: flatten-everything MLP (Stages 0-1); set: shared per-peak encoder + pooling (Stage 3); "
        "probe: frame-only mean-pooled diagnostic net (use with --loss mse); "
        "gn: learned Gauss-Newton layer (needs --aux)",
    )
    parser.add_argument(
        "--aux",
        default=None,
        help="aux table from make_dataset_aux.py (nominal offsets, pair index) for the same voxels",
    )
    parser.add_argument(
        "--subpixel",
        action="store_true",
        help="set net: measurement = centroid/frame minus the exact nominal (needs --aux)",
    )
    parser.add_argument("--gn-iters", type=int, default=1, help="gn arch: unrolled IRLS iterations")
    parser.add_argument(
        "--pairing", action="store_true", help="gn arch: entries see their other-detector partner"
    )
    parser.add_argument(
        "--ema",
        type=float,
        default=0.0,
        help="weight EMA decay per step (0 = off); validation is evaluated on the EMA weights",
    )
    parser.add_argument(
        "--checkpoint",
        choices=["best", "ema", "last"],
        default="best",
        help="final weights: best validation epoch (default), the EMA at the last epoch "
        "(needs --ema), or the last epoch",
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
    parser.add_argument(
        "--loss",
        choices=["nll", "decoupled", "mse", "mse-then-cov", "mse-then-nll"],
        default="nll",
        help="offset head loss. nll: Gaussian NLL (with --beta-nll if > 0); decoupled: MSE on the "
        "mean + NLL of the covariance at stopgrad(mean); mse: mean only; mse-then-cov / "
        "mse-then-nll: MSE for the first half of the epochs, then decoupled / plain NLL",
    )
    parser.add_argument(
        "--mse-scale", type=float, default=0.1, help="degrees; unit of the MSE term (decoupled)"
    )
    parser.add_argument("--val-frac", type=float, default=0.1)
    parser.add_argument(
        "--val-voxels",
        type=int,
        default=0,
        help="multi-voxel data: hold out this many training voxels (one per r_perp stratum, "
        "chosen with --seed) for early stopping instead of a random --val-frac of samples",
    )
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--device",
        choices=["cpu", "mps"],
        default="cpu",
        help="training and inference device; mps = Apple GPU (float32). Predictions are moved to the CPU "
        "and error statistics are computed there in float64",
    )
    parser.add_argument(
        "--no-frame", action="store_true", help="ablation: hide the frame channel (observer data)"
    )
    parser.add_argument(
        "--realistic-train",
        "--corrupt-train",  # historical name
        dest="realistic_train",
        choices=["none", "neighbours", "noise", "all"],
        default="none",
        help="make training windows realistic on the fly (needs dis_windows for neighbours/all; "
        "old name --corrupt-train)",
    )
    parser.add_argument(
        "--mask-padding",
        action="store_true",
        help="multi-voxel data: keep the zero-padded peak entries all-zero when applying the realism layer "
        "(default off: hot pixels/blobs also land in padding, as in the reported runs)",
    )
    parser.add_argument(
        "--eval-variants",
        default="clean",
        help="comma list of test-set variants (clean,neighbours,noise,all); the test windows "
        "are made realistic deterministically (orientation_eval.make_realistic_dataset)",
    )
    parser.add_argument(
        "--extra-variant",
        nargs="*",
        default=[],
        metavar="VARIANT/LABEL=NPZ",
        help="extra rows for a non-clean variant, e.g. neighbours/gn=pred.npz",
    )
    parser.add_argument("--results-json", default=None)
    parser.add_argument(
        "--save-model",
        default=None,
        help="torch file with the final weights (state_dict), the constructor arguments to "
        "rebuild the model and the input settings inference needs (load_model in "
        "scripts/perturbation_sweep.py)",
    )
    parser.add_argument(
        "--save-predictions",
        default=None,
        help="npz with test predictions (and Cholesky factors) for re-evaluation",
    )
    args = parser.parse_args()

    torch.manual_seed(args.seed)
    rng = np.random.default_rng(args.seed)

    assert not (args.subpixel and args.arch == "probe"), "--subpixel is ignored by --arch probe"
    tr = torch.load(args.train)
    az = torch.load(args.aux) if args.aux else None
    te = torch.load(args.test)
    n_peaks, window = tr["n_peaks"], tr["window_size"]
    R_nom = tr["R_nom"].numpy()
    windows, offsets = tr["windows"], tr["offsets_deg"].float()
    prep_cpu, in_channels = make_prep(tr, no_frame=args.no_frame)
    dev = torch.device(args.device)
    prep = lambda x: prep_cpu(x).to(dev)  # noqa: E731
    n = len(windows)
    if args.val_voxels > 0:
        assert "voxel_id" in tr, "--val-voxels needs multi-voxel data"
        train_idx, val_idx, val_voxels = split_by_voxel(
            tr["voxel_id"].numpy(), tr["r_perp_um"].numpy(), args.val_voxels, args.seed
        )
        n_val = len(val_idx)
        print(
            "validation voxels (held out of training, r_perp um): "
            + ", ".join(f"{v} ({tr['r_perp_um'][v]:.0f})" for v in val_voxels)
        )
    else:
        perm = rng.permutation(n)
        n_val = max(1, int(n * args.val_frac))
        val_idx, train_idx = perm[:n_val], perm[n_val:]
    print(
        f"train {len(train_idx)} / val {n_val} samples, {n_peaks} peaks x {window}x{window}, head={args.head}"
    )

    multi, vid_train = False, None
    model_kwargs: dict = {}
    if args.head == "quat":
        targets = torch.from_numpy(
            offsets_to_quaternions(offsets.numpy().astype(np.float64), R_nom)
        ).float()
        model_kwargs = dict(n_peaks=n_peaks, window_size=window, in_channels=in_channels)
        model = ToyOrientationNet(**model_kwargs)
    else:
        targets = offsets
        if args.arch in ("set", "probe", "gn"):
            assert "context" in tr, "the set architecture needs a dataset with per-peak context"
            context = tr["context"].float()  # (M, D), or (V, M, D) for multi-voxel data
            multi = bool(tr.get("multi_voxel", False))
            vid_train = tr["voxel_id"].to(torch.long) if multi else None
            if args.arch == "gn":
                assert args.aux, "--arch gn needs --aux"
                model_kwargs = dict(
                    window_size=window,
                    in_channels=in_channels,
                    context_dim=context.shape[-1],
                    n_iter=args.gn_iters,
                    frame_half_width=int(tr.get("frame_half_width", 4)),
                    frame_width_rad=float(az["frame_width_rad"]),
                    pairing=args.pairing,
                )
                model = GNLayerNet(**model_kwargs)
            elif args.arch == "probe":
                model_kwargs = dict(
                    context_dim=context.shape[-1],
                    frame_half_width=int(tr.get("frame_half_width", 4)),
                )
                model = FrameProbeNet(**model_kwargs)
            else:
                model_kwargs = dict(
                    window_size=window,
                    in_channels=in_channels,
                    context_dim=context.shape[-1],
                    use_measurements=not args.no_meas,
                    frame_half_width=int(tr.get("frame_half_width", 4)),
                    pool=args.pool,
                )
                model = PeakSetNet(**model_kwargs)
        else:
            model_kwargs = dict(n_peaks=n_peaks, window_size=window, in_channels=in_channels)
            model = ToyOffsetNet(**model_kwargs)
    aux_tab = {}
    if az is not None:
        assert (
            bool(az["multi_voxel"]) == multi and az["n_peaks"] == n_peaks
        ), "aux is for other data"
        assert torch.equal(
            az["voxel_indices"].long(),
            (tr["voxel_indices"] if multi else torch.tensor([tr["voxel_index"]])).long(),
        ), "aux is for other voxels"
        aux_tab = {"nom_off": az["nom_off"].float(), "pair_index": az["pair_index"].long()}
        if args.arch == "set" and not args.subpixel:
            aux_tab.pop("nom_off")
    assert args.aux or not args.subpixel, "--subpixel needs --aux"
    if args.arch in ("set", "probe", "gn"):
        assert args.head == "offset", "the set architecture predicts an offset and covariance"

    def make_forward(m):
        if args.arch not in ("set", "probe", "gn"):
            return m
        if multi:
            return lambda x, v: m(  # noqa: E731
                x, context[v.to(dev)], {k: t[v.to(dev)] for k, t in aux_tab.items()}
            )
        return lambda x: m(x, context, dict(aux_tab))  # noqa: E731

    forward = make_forward(model)
    n_params = sum(p.numel() for p in model.parameters())
    print(f"architecture {args.arch}: {n_params / 1e6:.2f}M parameters on {dev}", flush=True)
    model.to(dev)
    targets = targets.to(dev)
    if args.arch in ("set", "probe", "gn"):
        context = context.to(dev)
        aux_tab = {k: t.to(dev) for k, t in aux_tab.items()}
    opt = torch.optim.Adam(model.parameters(), lr=args.lr)
    sched = (
        torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=args.epochs) if args.cosine else None
    )

    two_phase = args.loss in ("mse-then-cov", "mse-then-nll")
    switch_epoch = args.epochs // 2 if two_phase else 0

    def loss_fn(x, y, v=None, phase=1, fwd=None):
        """Returns (loss, mean, chol). phase 0 = MSE-only first half of a two-phase run."""
        fwd = fwd if fwd is not None else forward
        if args.head == "quat":
            return quaternion_regression_loss(fwd(x), y), None, None
        mean, chol = fwd(x) if v is None else fwd(x, v)
        kind = args.loss
        if two_phase:
            kind = "mse" if phase == 0 else ("decoupled" if kind == "mse-then-cov" else "nll")
        if kind == "mse":
            loss = mse_deg_loss(mean, y, args.mse_scale)
        elif kind == "decoupled":
            loss = decoupled_nll_loss(mean, chol, y, args.mse_scale)
        else:
            loss = gaussian_nll_loss(mean, chol, y, beta=args.beta_nll)
        return loss, mean, chol

    import copy

    ema_model = None
    if args.ema > 0:
        ema_model = copy.deepcopy(model)
        for p_ in ema_model.parameters():
            p_.requires_grad_(False)
    eval_fwd = make_forward(ema_model) if ema_model is not None else forward
    eval_model = ema_model if ema_model is not None else model
    if args.checkpoint == "ema":
        assert ema_model is not None, "--checkpoint ema needs --ema"
    realism_cfg = RealismConfig.named(args.realistic_train)
    K_fh = int(tr.get("frame_half_width", 4))
    dis_train = tr.get("dis_windows")
    if realism_cfg is not None and realism_cfg.neighbours:
        assert dis_train is not None, "--realistic-train neighbours/all needs dis_windows"
    realism_gen = torch.Generator().manual_seed(args.seed + 7)

    valid_train = padding_mask(tr) if args.mask_padding else None

    def train_windows(ib):
        w = windows[ib]
        if realism_cfg is None:
            return w
        d = None if dis_train is None else dis_train[ib]
        v = None if valid_train is None else valid_train[ib]
        return make_realistic_windows(w, d, realism_cfg, K_fh, realism_gen, valid=v)

    if realism_cfg is not None:
        val_w = make_realistic_dataset(
            windows[val_idx],
            None if dis_train is None else dis_train[val_idx],
            args.realistic_train,
            K_fh,
            seed=args.seed + 99,
            valid=None if valid_train is None else valid_train[val_idx],
        )
    t0 = time.time()
    best_val, best_epoch, best_state = float("inf"), -1, None
    for epoch in range(args.epochs):
        model.train()
        total = 0.0
        phase = 0 if epoch < switch_epoch else 1
        for b in batches(len(train_idx), args.batch_size, rng):
            ib = train_idx[b]
            opt.zero_grad()
            loss = loss_fn(
                prep(train_windows(ib)), targets[ib], vid_train[ib] if multi else None, phase
            )[0]
            loss.backward()
            if args.clip > 0:
                torch.nn.utils.clip_grad_norm_(model.parameters(), args.clip)
            opt.step()
            if ema_model is not None:
                ema_update(ema_model, model, args.ema)
            total += loss.item() * len(ib)
        if sched is not None:
            sched.step()
        model.eval()
        eval_model.eval()
        vals, errs, sds = [], [], []
        with torch.no_grad():
            for b in batches(n_val, args.batch_size):
                lv, mv, cv = loss_fn(
                    prep(val_w[b] if realism_cfg is not None else windows[val_idx[b]]),
                    targets[val_idx[b]],
                    vid_train[val_idx[b]] if multi else None,
                    phase,
                    eval_fwd,
                )
                vals.append(lv.item())
                if mv is not None:
                    errs.append(((mv - targets[val_idx[b]]) ** 2).cpu())
                    sds.append(
                        torch.sqrt((cv @ cv.transpose(-1, -2)).diagonal(dim1=-2, dim2=-1)).cpu()
                    )
        val = np.mean(vals)
        diag = ""
        if errs:
            rms = torch.cat(errs).mean(0).sqrt()
            sd = torch.cat(sds).mean(0)
            diag = (
                f"  val rms z {rms[2]:.4f} perp {rms[:2].pow(2).mean().sqrt():.4f}"
                f"  sigma z {sd[2]:.4f} perp {sd[:2].mean():.4f}"
            )
        # In a two-phase run the phase-0 loss is not comparable with phase 1's: only
        # phase-1 epochs are eligible for the best-validation checkpoint.
        if phase == 1 and val < best_val:
            best_val, best_epoch = float(val), epoch + 1
            best_state = {k: v.detach().clone() for k, v in eval_model.state_dict().items()}
        print(
            f"epoch {epoch + 1:3d}/{args.epochs}  train {total / len(train_idx):11.6f}  val {val:11.6f}{diag}  ({time.time() - t0:.0f}s)",
            flush=True,
        )
    if best_state is None:
        raise RuntimeError("validation loss was never finite; lower --lr or check the data")
    model.load_state_dict(select_checkpoint_state(args.checkpoint, model, best_state, ema_model))
    if args.checkpoint == "best":
        print(f"restored best-validation weights from epoch {best_epoch} (val {best_val:.6f})")
    elif args.checkpoint == "ema":
        best_epoch = args.epochs
        print(f"using the EMA weights (decay {args.ema}) at the last epoch (val {val:.6f})")
    else:
        best_epoch = args.epochs
        print(f"using the last-epoch weights (val {val:.6f})")

    if args.save_model:
        torch.save(
            {
                "arch": args.arch,
                "head": args.head,
                "model_kwargs": model_kwargs,
                "state_dict": {k: v.detach().cpu() for k, v in model.state_dict().items()},
                "no_frame": args.no_frame,
                "frame_half_width": K_fh,
                "window_size": int(window),
                "n_peaks": int(n_peaks),
                "realistic_train": args.realistic_train,
                "mask_padding": args.mask_padding,
                "seed": args.seed,
                "best_epoch": best_epoch,
                "args": vars(args),
            },
            args.save_model,
        )
        print(f"saved model {args.save_model}")

    # ---- evaluation ------------------------------------------------------
    def evaluate(variant, win_te, extra_items, save_path):
        print(f"\n##### test set variant: {variant} #####")
        vid_test = te["voxel_id"].to(torch.long) if multi else None
        pred, chol = predict(model, args.head, win_te, args.batch_size, prep, forward, vid=vid_test)
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
                    f"note: Bayes file covers the first {n_b} of {len(truth)} test cases; "
                    "using those"
                )
                truth, pred, mags = truth[:n_b], pred[:n_b], mags[:n_b]
                chol = chol[:n_b] if chol is not None else None
            assert np.allclose(
                bz["offsets_deg"], truth[:n_b], atol=1e-5
            ), "Bayes file is for a different test set"
            bayes = bz

        extras = {}
        for item in extra_items:
            label, _, path = item.partition("=")
            ez = np.load(path)
            assert np.allclose(
                ez["truth_deg"], truth[: len(ez["truth_deg"])], atol=1e-5
            ), f"{path} is for a different test set"
            extras[label] = ez["pred_deg"]

        if save_path:
            np.savez(
                save_path,
                pred_deg=pred,
                truth_deg=truth,
                magnitudes_deg=mags,
                chol=chol if chol is not None else np.zeros(0),
                best_epoch=best_epoch,
            )
        rows = {}
        print(
            f"\n{'|delta|':>8} {'method':<20} {'n':>3} {'rms_z':>10} {'rms_perp':>10} "
            f"{'median_ang':>11} {'<0.5deg':>8} {'<0.1deg':>8}   (degrees)"
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
                            f"note: exact Bayes failed on {m.sum() - mb.sum()} case(s) at "
                            f"|delta| = {mag}; excluded from its row"
                        )
                    entries.append(("exact Bayes", error_summary(bayes["mean"][mb], truth[mb])))
                for name, s in entries:
                    print(
                        f"{mag:8.2f} {name:<20} {s['n']:3d} {s['rms_z']:10.5f} "
                        f"{s['rms_perp']:10.5f} {s['median_angle']:11.5f} "
                        f"{s['success_0p5']:8.0%} {s['success_0p1']:8.0%}"
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
                        f"{'':8} {'net sigma xyz':<14} {'':>3} {np.round(sd, 5)}   "
                        f"mean Mahalanobis^2 = {maha:.2f} (3.0 if calibrated)"
                    )
                rows[tag + str(mag)] = {name: s for name, s in entries} | extra

        if multi:
            r_perp, n_pk = te["r_perp_um"].numpy(), te["n_peaks_per_voxel"].numpy()
            per_voxel = []
            print(
                f"\n{'voxel':>5} {'r_perp':>7} {'peaks':>5} {'held':>5} {'n':>4}  "
                "median angle (deg), all magnitudes"
            )
            for v in range(len(r_perp)):
                mv = vid_test.numpy() == v
                ang = {"net": error_summary(pred[mv], truth[mv])["median_angle"]}
                for label, ep in extras.items():
                    assert len(ep) == len(truth), f"--extra {label} does not cover the test set"
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
                    f"{v:5d} {r_perp[v]:7.0f} {n_pk[v]:5d} "
                    f"{str(bool(te['held_out'][v])):>5} {mv.sum():4d}  "
                    + "  ".join(f"{k} {x:.4f}" for k, x in ang.items())
                )
            rows["per_voxel"] = per_voxel
            rv = np.array([d["r_perp_um"] for d in per_voxel])
            held = np.array([d["held_out"] for d in per_voxel])
            summary = {}
            for label in ["net"] + list(extras):
                ev = np.array([d["median_angle"][label] for d in per_voxel])
                summary[label] = dict(
                    corr_err_rperp=float(np.corrcoef(rv, ev)[0, 1]),
                    slope_err_per_100um=float(np.polyfit(rv, ev, 1)[0] * 100.0),
                    median_voxel_err_in_dist=float(np.median(ev[~held])),
                    median_voxel_err_held_out=float(np.median(ev[held])),
                )
                print(
                    f"{label:>14}: corr(voxel median error, r_perp) = "
                    f"{summary[label]['corr_err_rperp']:+.2f}, median per-voxel error in-dist "
                    f"{summary[label]['median_voxel_err_in_dist']:.4f} "
                    f"held-out {summary[label]['median_voxel_err_held_out']:.4f}"
                )
            rows["summary"] = summary

        return rows

    variants = [v for v in args.eval_variants.split(",") if v]
    if any(v in ("neighbours", "all") for v in variants):
        assert (
            te.get("dis_windows") is not None
        ), "--eval-variants neighbours/all needs dis_windows in the test set"
    valid_te = padding_mask(te) if args.mask_padding else None
    all_rows = {}
    for variant in variants:
        win_te = make_realistic_dataset(
            te["windows"],
            te.get("dis_windows"),
            "none" if variant == "clean" else variant,
            int(te.get("frame_half_width", 4)),
            valid=valid_te,
        )
        items = list(args.extra) if variant == "clean" else []
        items += [
            it.split("/", 1)[1] for it in args.extra_variant if it.split("/", 1)[0] == variant
        ]
        save = args.save_predictions
        if save:
            save = variant_path(save, variant)
        all_rows[variant] = evaluate(variant, win_te, items, save)

    if args.results_json:
        # {variant: rows}; summarize_results.py also reads the older single-variant layout
        Path(args.results_json).write_text(json.dumps(all_rows, indent=2))
        print(f"saved {args.results_json}")


if __name__ == "__main__":
    main()
