"""
Toy 4-layer fully-connected network for local voxel-orientation refinement.

Given a fixed-size stack of small detector windows around a voxel's known
observable peaks (produced by orientation_nn.render_local_windows), predicts
the voxel's absolute orientation as a unit quaternion -- the same role
RiemannianAdamOptimizer/MCOptimizer play today, but as a single forward pass
instead of an iterative search.

See icenine_py/MIGRATION_HISTORY.md for the approved plan this implements.
"""

import warnings
from typing import Dict, Optional, Tuple

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


N_MEAS = 5


def measurement_features(
    x: torch.Tensor, frame_half_width: int, nom_off: Optional[torch.Tensor] = None
) -> torch.Tensor:
    """Explicit per-peak measurements from decoded windows x (B, M, C, W, W) -> (B, M, 5).

    Features: present flag, lit-pixel count / 20, mean frame offset over lit pixels (frames),
    and the lit-pixel centroid (col, row) relative to the window centre (px). With a single
    (lit) channel the frame feature is 0. All features are 0 for absent peaks. These are the
    same per-spot quantities as orientation_baselines.extract_measurements.

    nom_off (M, 3) or (B, M, 3), optional (orientation_eval.nominal_offsets): subtracted from
    (centroid col, centroid row, frame) so that they become measurement minus the exact nominal
    prediction instead of measurement minus the window centre / nominal frame index. This
    removes the fixed 0-1 px (and +-0.5 frame) per-peak offset that otherwise has to be learned.
    """
    W = x.shape[-1]
    lit = x[:, :, 0]
    count = lit.sum(dim=(-2, -1))
    present = count > 0
    n = count.clamp(min=1.0)
    axis = torch.arange(W, dtype=x.dtype, device=x.device) + 0.5 - W / 2.0
    cx = (lit.sum(dim=-2) * axis).sum(dim=-1) / n
    cy = (lit.sum(dim=-1) * axis).sum(dim=-1) / n
    if x.shape[2] > 1:
        frame = (x[:, :, 1] * lit).sum(dim=(-2, -1)) / n * frame_half_width
    else:
        frame = torch.zeros_like(cx)
    if nom_off is not None:
        off = nom_off.to(x.dtype)
        if off.dim() == 2:
            off = off[None].expand(x.shape[0], -1, -1)
        cx, cy, frame = cx - off[..., 0], cy - off[..., 1], frame - off[..., 2]
    f = torch.stack([present.to(x.dtype), count / 20.0, frame, cx, cy], dim=-1)
    return f * present[..., None].to(x.dtype)


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
        use_measurements: bool = True,
        frame_half_width: int = 4,
        pool: str = "meanmax",
    ):
        super().__init__()
        assert pool in ("meanmax", "meansum", "all")
        self.in_channels = in_channels
        self.use_measurements = use_measurements
        self.frame_half_width = frame_half_width
        self.pool = pool
        c1, c2 = conv_channels
        self.conv1 = nn.Conv2d(in_channels + 2, c1, 3, stride=2, padding=1)
        self.conv2 = nn.Conv2d(c1, c2, 3, stride=2, padding=1)
        side = (window_size + 3) // 4
        self.fc_window = nn.Linear(c2 * side * side, feat_dim)
        self.fc_context1 = nn.Linear(context_dim, feat_dim)
        self.fc_context2 = nn.Linear(feat_dim, feat_dim)
        if use_measurements:
            self.fc_meas1 = nn.Linear(N_MEAS, feat_dim)
            self.fc_meas2 = nn.Linear(feat_dim, feat_dim)
        n_in = 3 if use_measurements else 2
        self.fc_peak1 = nn.Linear(n_in * feat_dim, feat_dim)
        self.fc_peak2 = nn.Linear(feat_dim, feat_dim)
        self.head1 = nn.Linear((3 if pool == "all" else 2) * feat_dim + 1, hidden)
        self.head2 = nn.Linear(hidden, hidden)
        self.head3 = nn.Linear(hidden, 9)
        axis = torch.linspace(-1.0, 1.0, window_size)
        yy, xx = torch.meshgrid(axis, axis, indexing="ij")
        self.register_buffer("coords", torch.stack([xx, yy])[None])  # (1, 2, W, W)

    def forward(
        self, x: torch.Tensor, context: torch.Tensor, aux: Optional[Dict[str, torch.Tensor]] = None
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        B, M = x.shape[:2]
        nom_off = None if aux is None else aux.get("nom_off")
        present = x[:, :, 0].flatten(2).amax(dim=-1) > 0  # channel 0 = lit
        w = x.reshape(B * M, *x.shape[2:])
        w = torch.cat([w, self.coords.expand(B * M, -1, -1, -1)], dim=1)
        w = F.relu(self.conv1(w))
        w = F.relu(self.conv2(w))
        w = F.relu(self.fc_window(w.flatten(1)))
        ctx = context if context.dim() == 3 else context[None].expand(B, -1, -1)
        c = F.relu(self.fc_context1(ctx.reshape(B * M, -1)))
        c = F.relu(self.fc_context2(c))
        parts = [w, c]
        if self.use_measurements:
            m = measurement_features(x, self.frame_half_width, nom_off).reshape(B * M, -1)
            m = F.relu(self.fc_meas1(m))
            parts.append(F.relu(self.fc_meas2(m)))
        f = F.relu(self.fc_peak1(torch.cat(parts, dim=-1)))
        f = F.relu(self.fc_peak2(f)).reshape(B, M, -1)
        mask = present[..., None]
        n = present.sum(dim=1, keepdim=True).clamp(min=1).to(f.dtype)
        total = (f * mask).sum(dim=1)
        mean = total / n
        if self.pool == "meansum":
            pooled = [mean, total / 30.0]  # sum pooling, scaled to O(1) for ~100 peaks
        else:
            second = f.masked_fill(~mask, float("-inf")).amax(dim=1)
            second = torch.where(torch.isfinite(second), second, torch.zeros_like(second))
            pooled = [mean, second] + ([total / 30.0] if self.pool == "all" else [])
        count = n / 100.0
        h = F.relu(self.head1(torch.cat(pooled + [count], dim=-1)))
        h = F.relu(self.head2(h))
        out = self.head3(h)
        return out[:, :3], cholesky_from_raw(out[:, 3:])


class FrameProbeNet(nn.Module):
    """Diagnostic: can the stage-axis offset survive mean pooling of frame measurements?

    Per peak, only the mean frame offset (measurement_features column 2) and the frame
    gradient d omega*/d delta (context columns 6:9, see BatchedObserver.peak_context) enter
    a small shared MLP; the outputs are mean-pooled over the present peaks and a linear
    head gives the offset. No windows are convolved and there are no pixel features. It
    returns a fixed identity covariance so it fits the same training loop (train it with
    --loss mse).

    forward(x, context): x (B, M, C, W, W) float, context (M, D) or (B, M, D).
    """

    def __init__(self, context_dim: int = 16, frame_half_width: int = 4, feat_dim: int = 64):
        super().__init__()
        assert context_dim >= 9
        self.frame_half_width = frame_half_width
        self.peak1 = nn.Linear(4, feat_dim)  # frame offset + dOmega/d(delta x, y, z)
        self.peak2 = nn.Linear(feat_dim, feat_dim)
        self.head1 = nn.Linear(feat_dim, feat_dim)
        self.head2 = nn.Linear(feat_dim, 3)

    def forward(
        self, x: torch.Tensor, context: torch.Tensor, aux: Optional[Dict[str, torch.Tensor]] = None
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        B, M = x.shape[:2]
        meas = measurement_features(x, self.frame_half_width)
        present = meas[..., 0] > 0
        ctx = context if context.dim() == 3 else context[None].expand(B, -1, -1)
        inp = torch.cat([meas[..., 2:3], ctx[..., 6:9]], dim=-1)
        f = F.relu(self.peak1(inp))
        f = F.relu(self.peak2(f))
        mask = present[..., None].to(f.dtype)
        pooled = (f * mask).sum(dim=1) / mask.sum(dim=1).clamp(min=1.0)
        mean = self.head2(F.relu(self.head1(pooled)))
        chol = torch.eye(3, dtype=mean.dtype, device=mean.device).expand(B, 3, 3)
        return mean, chol


# ---------------------------------------------------------------------------
# Gauss-Newton layer network
# ---------------------------------------------------------------------------

SQRT12 = 12.0**0.5
_SOFTPLUS_ONE = float(torch.log(torch.expm1(torch.tensor(1.0, dtype=torch.float64))))  # 0.5413


def inv3(A: torch.Tensor) -> torch.Tensor:
    """Inverse of batched 3x3 matrices by the adjugate (plain tensor ops: differentiable and
    available on every device, including Apple MPS which lacks linalg.solve)."""
    a, b, c = A[..., 0, 0], A[..., 0, 1], A[..., 0, 2]
    d, e, f = A[..., 1, 0], A[..., 1, 1], A[..., 1, 2]
    g, h, i = A[..., 2, 0], A[..., 2, 1], A[..., 2, 2]
    c00, c01, c02 = e * i - f * h, f * g - d * i, d * h - e * g
    det = a * c00 + b * c01 + c * c02
    adj = torch.stack(
        [
            c00,
            c * h - b * i,
            b * f - c * e,
            c01,
            a * i - c * g,
            c * d - a * f,
            c02,
            b * g - a * h,
            a * e - b * d,
        ],
        dim=-1,
    ).reshape(*A.shape)
    # guard against exactly singular A (no effect when |det| >= 1e-30): keep the sign, clamp |det|
    sign = torch.where(det < 0, -torch.ones_like(det), torch.ones_like(det))
    det = sign * det.abs().clamp(min=1e-30)
    return adj / det[..., None, None]


def chol3(S: torch.Tensor) -> torch.Tensor:
    """Lower Cholesky factor of batched symmetric positive-definite 3x3 matrices (closed form)."""
    l00 = S[..., 0, 0].clamp(min=1e-20).sqrt()
    l10 = S[..., 1, 0] / l00
    l20 = S[..., 2, 0] / l00
    l11 = (S[..., 1, 1] - l10**2).clamp(min=1e-20).sqrt()
    l21 = (S[..., 2, 1] - l20 * l10) / l11
    l22 = (S[..., 2, 2] - l20**2 - l21**2).clamp(min=1e-20).sqrt()
    z = torch.zeros_like(l00)
    return torch.stack([l00, z, z, l10, l11, z, l20, l21, l22], dim=-1).reshape(*S.shape)


def gn_normal_equations(
    J: torch.Tensor, y: torch.Tensor, w: torch.Tensor, ridge: float = 0.0
) -> Tuple[torch.Tensor, torch.Tensor]:
    """Weighted linear least squares min sum_i sum_r w_ir (y_ir - J_ir . d)^2 over d (3).

    J (B, M, R, 3), y (B, M, R), w (B, M, R) >= 0 (zero removes a row). Returns the solution
    d (B, 3) and A^-1 (B, 3, 3) for A = sum J^T W J + ridge I, i.e. the Gauss-Newton step
    from the linearisation point and its covariance (unit-variance rows)."""
    A = torch.einsum("bmrc,bmr,bmrk->bck", J, w, J)
    A = A + ridge * torch.eye(3, dtype=A.dtype, device=A.device)
    b = torch.einsum("bmrc,bmr,bmr->bc", J, w, y)
    Ainv = inv3(A)
    return torch.einsum("bck,bk->bc", Ainv, b), Ainv


class PeakEncoder(nn.Module):
    """Per-peak encoder shared by the set-style networks: window convolutions with coordinate
    channels + context + explicit measurement features -> (B, M, feat_dim)."""

    def __init__(
        self,
        window_size: int = 32,
        in_channels: int = 2,
        context_dim: int = 16,
        conv_channels: Tuple[int, int] = (8, 16),
        feat_dim: int = 64,
    ):
        super().__init__()
        c1, c2 = conv_channels
        self.conv1 = nn.Conv2d(in_channels + 2, c1, 3, stride=2, padding=1)
        self.conv2 = nn.Conv2d(c1, c2, 3, stride=2, padding=1)
        side = (window_size + 3) // 4
        self.fc_window = nn.Linear(c2 * side * side, feat_dim)
        self.fc_context1 = nn.Linear(context_dim, feat_dim)
        self.fc_context2 = nn.Linear(feat_dim, feat_dim)
        self.fc_meas1 = nn.Linear(N_MEAS, feat_dim)
        self.fc_meas2 = nn.Linear(feat_dim, feat_dim)
        self.fc_peak1 = nn.Linear(3 * feat_dim, feat_dim)
        self.fc_peak2 = nn.Linear(feat_dim, feat_dim)
        axis = torch.linspace(-1.0, 1.0, window_size)
        yy, xx = torch.meshgrid(axis, axis, indexing="ij")
        self.register_buffer("coords", torch.stack([xx, yy])[None])

    def forward(self, x: torch.Tensor, ctx: torch.Tensor, meas: torch.Tensor) -> torch.Tensor:
        B, M = x.shape[:2]
        w = x.reshape(B * M, *x.shape[2:])
        w = torch.cat([w, self.coords.expand(B * M, -1, -1, -1)], dim=1)
        w = F.relu(self.conv1(w))
        w = F.relu(self.conv2(w))
        w = F.relu(self.fc_window(w.flatten(1)))
        c = F.relu(self.fc_context1(ctx.reshape(B * M, -1)))
        c = F.relu(self.fc_context2(c))
        m = F.relu(self.fc_meas1(meas.reshape(B * M, -1)))
        m = F.relu(self.fc_meas2(m))
        f = F.relu(self.fc_peak1(torch.cat([w, c, m], dim=-1)))
        return F.relu(self.fc_peak2(f)).reshape(B, M, -1)


class GNLayerNet(nn.Module):
    """Learned Gauss-Newton layer: the physics is in the architecture.

    Per present peak i the shared encoder produces a learned reliability weight W_i (one per
    measured row: centroid column, centroid row, frame) and a correction dy_i to the measurement.
    The measurement y_i (the lit centroid minus the exact nominal centroid, the frame minus the
    exact nominal crossing, in units of their quantisation sigma, see GN in
    orientation_baselines) and the Jacobian J_i = d(y_i)/d(delta) (3x3, from the context: spot
    Jacobian Gamma and d omega*/d delta; this is what carries the r_perp parallax) are combined by
    the pooled normal equations  A = sum J^T W J + ridge,  delta = A^-1 sum J^T W (y + dy).
    The covariance is D A^-1 D with a learned diagonal D (calibration). With weights 1 and
    corrections 0 (the initial state) the layer equals, at ridge 0, one undamped Gauss-Newton
    step from the nominal orientation (CentroidGaussNewton.solve_linear), which on this problem
    equals the converged solution to within quantisation noise for offsets up to 1 degree. The
    default ridge 1e-3 (in sigma units) differs from that step only slightly for well-conditioned A.

    n_iter > 1 unrolls an IRLS-style loop: the weight/correction head additionally sees the
    residual y - J delta and the current delta, so it can down-weight outliers.
    pairing=True: each entry also sees the other-detector entry of the same ray (aux
    "pair_index"): their encodings are mixed before the heads by a residual update
    f += has * pair2(relu(pair1([f, f_partner, has]))) (pair2 zero-initialised and linear, so
    the net starts as the unpaired one and the update trains), and the head sees the partner's
    residual (a pair-consistency check).

    forward(x, context, aux): x (B, M, 2, W, W) decoded windows, context (M, D) or (B, M, D),
    aux dict with nom_off (M|B, M, 3) and, for pairing, pair_index (M|B, M). Returns
    (delta_deg (B, 3), chol (B, 3, 3)).
    """

    def __init__(
        self,
        window_size: int = 32,
        in_channels: int = 2,
        context_dim: int = 16,
        feat_dim: int = 64,
        hidden: int = 64,
        n_iter: int = 1,
        frame_half_width: int = 4,
        frame_width_rad: float = 0.017453292519943295,
        pairing: bool = False,
        ridge: float = 1e-3,
        delta_scale: float = 0.05,
    ):
        super().__init__()
        assert in_channels == 2, "GNLayerNet needs the frame channel"
        self.frame_half_width = frame_half_width
        self.frame_width_rad = frame_width_rad
        self.n_iter = n_iter
        self.pairing = pairing
        self.ridge = ridge
        self.delta_scale = delta_scale  # degrees: unit of delta inside the layer
        self.encoder = PeakEncoder(window_size, in_channels, context_dim, feat_dim=feat_dim)
        if pairing:
            self.pair1 = nn.Linear(2 * feat_dim + 1, feat_dim)
            self.pair2 = nn.Linear(feat_dim, feat_dim)
            nn.init.zeros_(self.pair2.weight)  # start as the unpaired network
            nn.init.zeros_(self.pair2.bias)
        head_in = feat_dim + 6 + (4 if pairing else 0)
        self.head1 = nn.Linear(head_in, hidden)
        self.head2 = nn.Linear(hidden, 6)  # 3 weights + 3 measurement corrections
        nn.init.zeros_(self.head2.weight)
        with torch.no_grad():
            self.head2.bias.copy_(
                torch.tensor([_SOFTPLUS_ONE] * 3 + [0.0] * 3, dtype=self.head2.bias.dtype)
            )
        self.cov1 = nn.Linear(feat_dim, hidden)
        self.cov2 = nn.Linear(hidden, 3)
        nn.init.zeros_(self.cov2.weight)
        nn.init.zeros_(self.cov2.bias)

    def physics(
        self, x: torch.Tensor, context: torch.Tensor, aux: Optional[Dict[str, torch.Tensor]]
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        """(present (B, M), y (B, M, 3), J (B, M, 3, 3), ctx (B, M, D), meas (B, M, 5)):
        the measurement and Jacobian in sigma units, J per delta_scale degrees."""
        B, M = x.shape[:2]
        ctx = context if context.dim() == 3 else context[None].expand(B, -1, -1)
        nom_off = None if aux is None else aux.get("nom_off")
        if nom_off is None:
            warnings.warn(
                "GNLayerNet without aux['nom_off']: measurements fall back to the window centre, "
                "which carries a sub-pixel bias (pass the make_dataset_aux.py table via --aux)",
                stacklevel=3,
            )
        meas = measurement_features(x, self.frame_half_width, nom_off)
        present = meas[..., 0] > 0
        pm = present[..., None].to(x.dtype)
        y = SQRT12 * torch.stack([meas[..., 3], meas[..., 4], meas[..., 2]], dim=-1) * pm
        gamma = ctx[..., 0:6].reshape(B, M, 2, 3) * 20.0  # px / deg
        gomega = ctx[..., 6:9] * (0.017453292519943295 / self.frame_width_rad)  # frames / deg
        J = SQRT12 * torch.cat([gamma, gomega[:, :, None, :]], dim=2) * self.delta_scale
        return present, y, J * pm[..., None], ctx, meas

    def forward(
        self, x: torch.Tensor, context: torch.Tensor, aux: Optional[Dict[str, torch.Tensor]] = None
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        B, M = x.shape[:2]
        present, y, J, ctx, meas = self.physics(x, context, aux)
        pm = present[..., None].to(y.dtype)
        f = self.encoder(x, ctx, meas)
        pidx = has = None
        if self.pairing:
            assert aux is not None and "pair_index" in aux, "pairing needs aux['pair_index']"
            pidx = aux["pair_index"].to(torch.long)
            if pidx.dim() == 1:
                pidx = pidx[None].expand(B, -1)
            safe = pidx.clamp(min=0)
            has = ((pidx >= 0) & torch.gather(present, 1, safe)) & present
            hm = has[..., None].to(f.dtype)
            fp = torch.gather(f, 1, safe[..., None].expand(-1, -1, f.shape[-1])) * hm
            # residual update with a linear output (a ReLU after the zero-initialised layer has
            # ReLU'(0) = 0 and kills the gradient of pair1/pair2 forever); masked so that
            # entries without a partner are unchanged
            f = f + hm * self.pair2(F.relu(self.pair1(torch.cat([f, fp, hm], dim=-1))))
        delta = torch.zeros(B, 3, dtype=y.dtype, device=y.device)
        for _ in range(self.n_iter):
            rho = y - torch.einsum("bmrc,bc->bmr", J, delta)  # residual at current delta
            inp = [f, torch.asinh(rho), delta[:, None, :].expand(-1, M, -1)]
            if self.pairing:
                rp = torch.gather(rho, 1, safe[..., None].expand(-1, -1, 3)) * hm
                inp += [torch.asinh(rp), hm]
            h = F.relu(self.head1(torch.cat(inp, dim=-1)))
            out = self.head2(h)
            w = F.softplus(out[..., :3]) * pm
            delta, Ainv = gn_normal_equations(J, y + out[..., 3:], w, self.ridge)
        n = present.sum(dim=1, keepdim=True).clamp(min=1).to(f.dtype)
        pooled = (f * pm).sum(dim=1) / n
        d = torch.exp(self.cov2(F.relu(self.cov1(pooled))))
        cov = d[:, :, None] * Ainv * d[:, None, :] * self.delta_scale**2
        return delta * self.delta_scale, chol3(cov)
