"""AR4 step 1: the model definitions, fixed before any training (plan, "Models").

  B0  the hit ratio alone (today's rule): score = hit ratio.
  B1  a logistic regression on the summary numbers (``SUMMARY_USED``), class-balanced.
  N1  a small set network on the per-peak rows plus the summary numbers. PyTorch, CPU, ONE training
      run with generator seed 0 (``SEED``), ``EPOCHS`` epochs. Layers (all sizes below):
        per peak   Linear(5, 64) - ReLU - Linear(64, 64) - ReLU   (the same weights for all peaks)
        pooled     mean and max over the pair's peaks -> 128 numbers
        head       concat(pooled 128, summary 6 standardized) -> Linear(134, 128) - ReLU
                   -> Linear(128, 1) -> sigmoid           (about 22,000 parameters, under 50,000)
      Loss: binary cross-entropy of "right" (error at most 1 degree), each class carrying half of
      the weight (class-balanced). Adam, learning rate 1e-3, weight decay 1e-5, batch 256.
  lit  (prediction 4) the share of predicted peaks with any lit pixel, alone, no training.

The per-peak columns are ``verify_common.PEAK_COLS`` = rank, dist, omega, share, pixels; the model
scales them by the fixed constants in ``scale_peaks`` (no fitted statistics). The summary numbers
are standardized with the mean and standard deviation of the training pairs. ``best_votes`` (the
last summary number) is not computed (see ``verify_inputs.py``) and is not used.
"""

from typing import Tuple

import numpy as np
import torch
from torch import nn

SEED = 0
EPOCHS = 20
BATCH = 256
LR = 1e-3
WEIGHT_DECAY = 1e-5
MAX_PEAKS = 160  # rows per pair kept (a pair has about 65 at Q-max 8)
SUMMARY_USED = (0, 1, 2, 3, 4, 5)  # hit, cost, conf, lit_share, low_quarter, n_explained
HID_PEAK = 64
HID_HEAD = 128
LOGISTIC_EPOCHS = 3000
LOGISTIC_LR = 1.0


def scale_peaks(rows: np.ndarray) -> np.ndarray:
    """Fixed scaling of the per-peak columns (rank, dist, omega, share, pixels)."""
    r = np.asarray(rows, dtype=np.float32)
    out = np.empty_like(r)
    out[..., 0] = r[..., 0] / 10.0
    out[..., 1] = r[..., 1] * 2.0
    out[..., 2] = r[..., 2] / np.pi
    out[..., 3] = r[..., 3]
    out[..., 4] = np.log1p(r[..., 4]) / 6.0
    return out


class SetNet(nn.Module):
    """N1: shared layer per peak, mean and max pooling, two layers with the summary numbers."""

    def __init__(self, n_summary: int = len(SUMMARY_USED)) -> None:
        super().__init__()
        self.peak = nn.Sequential(
            nn.Linear(5, HID_PEAK), nn.ReLU(), nn.Linear(HID_PEAK, HID_PEAK), nn.ReLU()
        )
        self.head = nn.Sequential(
            nn.Linear(2 * HID_PEAK + n_summary, HID_HEAD),
            nn.ReLU(),
            nn.Linear(HID_HEAD, 1),
        )

    def forward(self, x: torch.Tensor, mask: torch.Tensor, summ: torch.Tensor) -> torch.Tensor:
        """x (B, L, 5) scaled rows, mask (B, L) bool, summ (B, n_summary) standardized -> logits."""
        h = self.peak(x)
        m = mask.unsqueeze(-1).to(h.dtype)
        n = m.sum(dim=1).clamp(min=1.0)
        mean = (h * m).sum(dim=1) / n
        mx = (h.masked_fill(~mask.unsqueeze(-1), -1e9)).max(dim=1).values
        mx = torch.where(mask.any(dim=1, keepdim=True), mx, torch.zeros_like(mx))
        return self.head(torch.cat([mean, mx, summ], dim=1)).squeeze(-1)


def n_parameters(model: nn.Module) -> int:
    return int(sum(p.numel() for p in model.parameters()))


def gather_batch(
    flat: np.ndarray, offsets: np.ndarray, idx: np.ndarray, max_len: int = MAX_PEAKS
) -> Tuple[np.ndarray, np.ndarray]:
    """Padded scaled rows (B, max_len, 5) and mask (B, max_len) of the pairs ``idx``."""
    start = offsets[idx]
    cnt = np.minimum(offsets[idx + 1] - start, max_len)
    pos = np.arange(max_len)[None, :]
    mask = pos < cnt[:, None]
    take = np.where(mask, start[:, None] + pos, 0)
    x = scale_peaks(flat[take]) if len(flat) else np.zeros((len(idx), max_len, 5), np.float32)
    x = x * mask[..., None]
    return x.astype(np.float32), mask
