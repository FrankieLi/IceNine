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
