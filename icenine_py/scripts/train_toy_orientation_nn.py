#!/usr/bin/env python3
"""
Train ToyOrientationNet on a cached windowed-orientation dataset.

Reports success-rate-at-threshold (1/2/5 deg), directly comparable to the
existing HP-sweep baseline (Riemannian Adam 96%/41%/0% @1/2/5deg, MC
92%/40%/6%) from icenine_py/MIGRATION_HISTORY.md.

Usage:
  cd icenine_py
  uv run python scripts/generate_toy_orientation_dataset.py --n-samples 20000
  uv run python scripts/train_toy_orientation_nn.py --dataset scripts/toy_orientation_dataset_threevoxels_v0_n20000.pt
"""

import argparse
from pathlib import Path

import torch
from torch.utils.data import DataLoader, random_split


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", required=True, help="Path to a .pt file from generate_toy_orientation_dataset.py")
    parser.add_argument("--epochs", type=int, default=50)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--val-frac", type=float, default=0.1)
    parser.add_argument("--test-frac", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--checkpoint", default=None, help="Output checkpoint path (default: derived from --dataset)")
    args = parser.parse_args()

    from icenine.orientation_nn import OrientationDataset, quat_misorientation_deg_batch, quaternion_regression_loss
    from icenine.toy_orientation_model import ToyOrientationNet

    dataset = OrientationDataset(args.dataset)
    n = len(dataset)
    n_val = max(1, int(n * args.val_frac))
    n_test = max(1, int(n * args.test_frac))
    n_train = n - n_val - n_test
    if n_train <= 0:
        raise ValueError(f"Dataset too small ({n} samples) for the requested val/test split")

    generator = torch.Generator().manual_seed(args.seed)
    train_set, val_set, test_set = random_split(dataset, [n_train, n_val, n_test], generator=generator)
    print(f"Dataset: {n} samples -> train={n_train} val={n_val} test={n_test}")

    train_loader = DataLoader(train_set, batch_size=args.batch_size, shuffle=True)
    val_loader = DataLoader(val_set, batch_size=args.batch_size, shuffle=False)
    test_loader = DataLoader(test_set, batch_size=args.batch_size, shuffle=False)

    n_peaks = dataset.windows.shape[1]
    window_size = dataset.windows.shape[2]
    model = ToyOrientationNet(n_peaks=n_peaks, window_size=window_size)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)

    def evaluate(loader):
        model.eval()
        total_loss = 0.0
        all_deg = []
        with torch.no_grad():
            for windows, q_true in loader:
                q_pred = model(windows)
                total_loss += quaternion_regression_loss(q_pred, q_true).item() * windows.shape[0]
                all_deg.append(quat_misorientation_deg_batch(q_pred, q_true))
        deg = torch.cat(all_deg)
        avg_loss = total_loss / len(loader.dataset)
        return avg_loss, deg

    for epoch in range(args.epochs):
        model.train()
        train_loss = 0.0
        for windows, q_true in train_loader:
            optimizer.zero_grad()
            q_pred = model(windows)
            loss = quaternion_regression_loss(q_pred, q_true)
            loss.backward()
            optimizer.step()
            train_loss += loss.item() * windows.shape[0]
        train_loss /= n_train

        val_loss, val_deg = evaluate(val_loader)
        print(
            f"Epoch {epoch + 1}/{args.epochs}  train_loss={train_loss:.4f}  "
            f"val_loss={val_loss:.4f}  val_misori_median={val_deg.median():.3f}deg"
        )

    test_loss, test_deg = evaluate(test_loader)
    print(f"\nTest loss: {test_loss:.4f}")
    for threshold in (1.0, 2.0, 5.0):
        success_rate = (test_deg <= threshold).float().mean().item()
        print(f"  Success @ {threshold:g} deg: {success_rate:.1%}")
    print("Baseline (bench_hp_sweep.py, Example2.ThreeVoxels): "
          "Riemannian Adam 96%/41%/0% @1/2/5deg, MC 92%/40%/6%")

    checkpoint_path = args.checkpoint or str(Path(args.dataset).with_suffix("")) + "_model.pt"
    torch.save({"model_state_dict": model.state_dict(), "n_peaks": n_peaks, "window_size": window_size}, checkpoint_path)
    print(f"Saved checkpoint to {checkpoint_path}")


if __name__ == "__main__":
    main()
