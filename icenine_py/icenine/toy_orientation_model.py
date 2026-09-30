"""
Toy 4-layer fully-connected network for local voxel-orientation refinement.

Given a fixed-size stack of small detector windows around a voxel's known
observable peaks (produced by orientation_nn.render_local_windows), predicts
the voxel's absolute orientation as a unit quaternion -- the same role
RiemannianAdamOptimizer/MCOptimizer play today, but as a single forward pass
instead of an iterative search.

See icenine_py/MIGRATION_HISTORY.md for the approved plan this implements.
"""

from typing import Tuple

import torch
import torch.nn as nn
import torch.nn.functional as F

from .orientation_nn import cholesky_from_raw


class ToyOrientationNet(nn.Module):
    def __init__(
        self,
        n_peaks: int,
        window_size: int = 32,
        hidden: Tuple[int, int, int] = (512, 256, 128),
        in_channels: int = 1,
    ):
        super().__init__()
        in_dim = n_peaks * in_channels * window_size * window_size
        self.fc1 = nn.Linear(in_dim, hidden[0])
        self.fc2 = nn.Linear(hidden[0], hidden[1])
        self.fc3 = nn.Linear(hidden[1], hidden[2])
        self.fc4 = nn.Linear(hidden[2], 4)  # raw quaternion [w, x, y, z]

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        h = F.relu(self.fc1(x.flatten(1)))
        h = F.relu(self.fc2(h))
        h = F.relu(self.fc3(h))
        q = self.fc4(h)
        return q / q.norm(dim=-1, keepdim=True)


class ToyOffsetNet(nn.Module):
    """Same trunk as ToyOrientationNet, but predicts a rotation-vector offset
    from the nominal orientation (degrees) and a full 3x3 covariance.

    The head has 9 outputs: 3 for the mean, 6 for the Cholesky factor of the
    covariance (see orientation_nn.cholesky_from_raw). Train with
    orientation_nn.gaussian_nll_loss.
    """

    def __init__(
        self,
        n_peaks: int,
        window_size: int = 32,
        hidden: Tuple[int, int, int] = (512, 256, 128),
        in_channels: int = 1,
    ):
        super().__init__()
        in_dim = n_peaks * in_channels * window_size * window_size
        self.fc1 = nn.Linear(in_dim, hidden[0])
        self.fc2 = nn.Linear(hidden[0], hidden[1])
        self.fc3 = nn.Linear(hidden[1], hidden[2])
        self.fc4 = nn.Linear(hidden[2], 9)

    def forward(self, x: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        h = F.relu(self.fc1(x.flatten(1)))
        h = F.relu(self.fc2(h))
        h = F.relu(self.fc3(h))
        out = self.fc4(h)
        return out[:, :3], cholesky_from_raw(out[:, 3:])


class PeakSetNet(nn.Module):
    """Shared per-peak encoder with masked pooling: predicts an offset and covariance.

    Every peak's window goes through the same small convolutional encoder (with
    coordinate channels, so absolute position inside the window is available), is
    combined with a description of how that peak responds to an orientation offset
    (its context vector, see BatchedObserver.peak_context), and the per-peak features
    are pooled over the peaks that are present (mean and max). A small MLP head maps
    the pooled features to the offset (3, degrees) and a Cholesky covariance (6).

    The parameter count does not depend on the number of peaks, and the network is
    invariant to peak order, so it can be applied to any peak set.

    forward(x, context): x (B, M, C, W, W) float, context (M, D) or (B, M, D).
    """

    def __init__(
        self,
        window_size: int = 32,
        in_channels: int = 2,
        context_dim: int = 16,
        conv_channels: Tuple[int, int] = (8, 16),
        feat_dim: int = 64,
        hidden: int = 128,
    ):
        super().__init__()
        self.in_channels = in_channels
        c1, c2 = conv_channels
        self.conv1 = nn.Conv2d(in_channels + 2, c1, 3, stride=2, padding=1)
        self.conv2 = nn.Conv2d(c1, c2, 3, stride=2, padding=1)
        side = (window_size + 3) // 4
        self.fc_window = nn.Linear(c2 * side * side, feat_dim)
        self.fc_context1 = nn.Linear(context_dim, feat_dim)
        self.fc_context2 = nn.Linear(feat_dim, feat_dim)
        self.fc_peak1 = nn.Linear(2 * feat_dim, feat_dim)
        self.fc_peak2 = nn.Linear(feat_dim, feat_dim)
        self.head1 = nn.Linear(2 * feat_dim + 1, hidden)
        self.head2 = nn.Linear(hidden, hidden)
        self.head3 = nn.Linear(hidden, 9)
        axis = torch.linspace(-1.0, 1.0, window_size)
        yy, xx = torch.meshgrid(axis, axis, indexing="ij")
        self.register_buffer("coords", torch.stack([xx, yy])[None])  # (1, 2, W, W)

    def forward(self, x: torch.Tensor, context: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        B, M = x.shape[:2]
        present = x[:, :, 0].flatten(2).amax(dim=-1) > 0  # channel 0 = lit
        w = x.reshape(B * M, *x.shape[2:])
        w = torch.cat([w, self.coords.expand(B * M, -1, -1, -1)], dim=1)
        w = F.relu(self.conv1(w))
        w = F.relu(self.conv2(w))
        w = F.relu(self.fc_window(w.flatten(1)))
        ctx = context if context.dim() == 3 else context[None].expand(B, -1, -1)
        c = F.relu(self.fc_context1(ctx.reshape(B * M, -1)))
        c = F.relu(self.fc_context2(c))
        f = F.relu(self.fc_peak1(torch.cat([w, c], dim=-1)))
        f = F.relu(self.fc_peak2(f)).reshape(B, M, -1)
        mask = present[..., None]
        n = present.sum(dim=1, keepdim=True).clamp(min=1).to(f.dtype)
        mean = (f * mask).sum(dim=1) / n
        mx = f.masked_fill(~mask, float("-inf")).amax(dim=1)
        mx = torch.where(torch.isfinite(mx), mx, torch.zeros_like(mx))
        frac = present.float().mean(dim=1, keepdim=True)
        h = F.relu(self.head1(torch.cat([mean, mx, frac], dim=-1)))
        h = F.relu(self.head2(h))
        out = self.head3(h)
        return out[:, :3], cholesky_from_raw(out[:, 3:])
