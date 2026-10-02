"""
Stage 2 baselines that use no learning.

`CentroidGaussNewton` estimates the orientation offset delta from the same
thresholded, frame-coded windows the networks see, by damped Gauss-Newton
(Levenberg-Marquardt) on per-spot measurements:

  - the frame in which the spot is recorded, as the frame's central rotation angle
    (variance dw^2/12: the frame quantisation error), and
  - the spot's pixel centroid, taken as the mean of its lit pixel centres
    (variance 1/12 px^2 per coordinate: the pixel quantisation error).

Predictions come from the batched observer (the exact simulator ray tracing), so
the model is exact and only the measurements are quantised. Jacobians are central
finite differences of the observer (float64); the docs' closed form
(nn_inverse_problem_formulation.md 3.6.3) agrees with them to ~1e-4. The covariance
returned is (J^T W J)^-1 at the solution: the independent-quantisation estimate of
docs 3.6.4, evaluated on the peaks actually recorded.
"""

from dataclasses import dataclass
from typing import Dict, Optional, Sequence, Tuple

import numpy as np
import torch

from .orientation_eval import DEG, BatchedObserver, WindowSpec


@dataclass
class Measurements:
    """Per-spot measurements extracted from one sample's frame-coded windows."""

    used: np.ndarray  # (M,) bool: spot is present in its window
    omega: np.ndarray  # (M,) rad: central rotation angle of the recorded frame
    centroid: np.ndarray  # (M, 2) px: (col, row) mean of lit pixel centres


def frame_center_omega(observer: BatchedObserver, frame_index: np.ndarray) -> np.ndarray:
    """Central rotation angle (rad) of wedge/frame indices; NaN where the frame is unknown."""
    bin_of = {int(w): b for b, w in enumerate(observer.range_index.tolist()) if w >= 0}
    out = np.full(len(frame_index), np.nan)
    for i, j in enumerate(frame_index):
        b = bin_of.get(int(j))
        if b is not None:
            out[i] = observer.range_low + (b + 0.5) * observer.range_width
    return out


def extract_measurements(
    windows: torch.Tensor,
    spec: WindowSpec,
    observer: BatchedObserver,
    detectors: Optional[Sequence[int]] = None,
) -> Measurements:
    """Frame index and lit-pixel centroid of every spot in one sample's windows (M, W, W).

    detectors: keep only spots recorded on these detector indices (default: all), the
    others are marked unused, so the fit sees one detector alone (diagnostics).
    """
    M, K = windows.shape[0], spec.frame_half_width
    w = windows.numpy().astype(np.int64)
    used = np.zeros(M, dtype=bool)
    frame = np.zeros(M, dtype=np.int64)
    centroid = np.zeros((M, 2))
    for m in range(M):
        rows, cols = np.nonzero(w[m])
        if len(rows) == 0:
            continue
        used[m] = True
        frame[m] = spec.frame0[m] + int(w[m, rows[0], cols[0]]) - 1 - K
        centroid[m] = (spec.col0[m] + cols.mean() + 0.5, spec.row0[m] + rows.mean() + 0.5)
    omega = frame_center_omega(observer, frame)
    used &= np.isfinite(omega)
    if detectors is not None:
        used &= np.isin(observer.det_idx.numpy(), list(detectors))
    return Measurements(used=used, omega=omega, centroid=centroid)


def _wrap(a: np.ndarray) -> np.ndarray:
    return (a + np.pi) % (2 * np.pi) - np.pi


class CentroidGaussNewton:
    """Levenberg-Marquardt fit of delta (degrees) to frame and centroid measurements."""

    def __init__(
        self,
        observer: BatchedObserver,
        fd_step_deg: float = 0.002,
        max_iter: int = 30,
        tol_deg: float = 1e-7,
        huber_c: Optional[float] = None,
    ):
        """huber_c: if set, robust fit: Huber loss with this threshold (in units of the
        quantisation sigma, i.e. normalised residuals) fitted by iteratively reweighted
        Levenberg-Marquardt; None = plain least squares (the Stage 2 baseline)."""
        self.huber_c = huber_c
        self.obs = observer
        self.h = fd_step_deg
        self.max_iter = max_iter
        self.tol = tol_deg
        self.sigma_omega = observer.frame_width_rad / np.sqrt(12.0)
        self.sigma_px = 1.0 / np.sqrt(12.0)

    def _residuals(self, deltas: np.ndarray, meas: Measurements, idx: np.ndarray):
        """Weighted residuals (B, n, 3) and validity (B, n) for the used spots idx."""
        o = self.obs.observe(torch.as_tensor(deltas, dtype=self.obs.dtype))
        om = o.omega.numpy()[:, idx]
        cen = o.verts.mean(dim=2).numpy()[:, idx]  # (B, n, 2)
        r_om = _wrap(om - meas.omega[idx][None]) / self.sigma_omega
        r_c = (cen - meas.centroid[idx][None]) / self.sigma_px
        res = np.concatenate([r_om[..., None], r_c], axis=-1)
        return res, o.present.numpy()[:, idx]

    def _weights(self, r: np.ndarray) -> np.ndarray:
        """Per-residual IRLS weights (all ones for plain least squares)."""
        if self.huber_c is None:
            return np.ones_like(r)
        return np.minimum(1.0, self.huber_c / np.maximum(np.abs(r), 1e-12))

    def _cost(self, r: np.ndarray) -> float:
        if self.huber_c is None:
            return 0.5 * float(r @ r)
        c, a = self.huber_c, np.abs(r)
        return float(np.where(a <= c, 0.5 * a * a, c * a - 0.5 * c * c).sum())

    def _linearize(
        self, meas: Measurements, delta: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Linearise at delta by central differences (the used spots present at delta and at
        every probe). Returns the sigma-normalised residual vector r0 (3 per kept spot, no
        Huber weights), the Jacobian J (3 per kept spot x 3, per degree) of r0, and the (n_used,)
        bool mask of kept spots."""
        idx = np.nonzero(meas.used)[0]
        eye = np.eye(3)
        pts = np.stack([delta] + [delta + s * self.h * eye[i] for i in range(3) for s in (1, -1)])
        res, valid = self._residuals(pts, meas, idx)
        ok = valid.all(axis=0)
        J = np.stack(
            [
                (res[1 + 2 * i][ok] - res[2 + 2 * i][ok]).reshape(-1) / (2 * self.h)
                for i in range(3)
            ],
            axis=1,
        )
        return res[0][ok].reshape(-1), J, ok

    def information(self, meas: Measurements, delta: Optional[np.ndarray] = None) -> np.ndarray:
        """J^T J (3x3, per degree^2) of the used spots at delta (default nominal), with
        sigma-normalised residuals: the Fisher information of the independent-quantisation
        model, for conditioning diagnostics. It carries no Huber IRLS weights even when
        huber_c is set, unlike the covariance returned by solve() (J^T W_huber J)."""
        delta = np.zeros(3) if delta is None else np.asarray(delta, dtype=np.float64)
        _r0, J, _ok = self._linearize(meas, delta)
        return J.T @ J

    def solve_linear(self, meas: Measurements) -> np.ndarray:
        """One undamped Gauss-Newton step from the nominal orientation (the linear-model
        least-squares estimate, what a network layer with unit weights computes)."""
        r0, J, _ok = self._linearize(meas, np.zeros(3))
        return -np.linalg.solve(J.T @ J, J.T @ r0)

    def solve(self, meas: Measurements, delta0: Optional[np.ndarray] = None) -> Dict[str, object]:
        idx = np.nonzero(meas.used)[0]
        delta = np.zeros(3) if delta0 is None else np.asarray(delta0, dtype=np.float64).copy()
        eye = np.eye(3)
        lam = 1e-3
        status, n_iter, n_accepted = "max_iter", 0, 0
        for n_iter in range(1, self.max_iter + 1):
            r0, J, ok = self._linearize(meas, delta)  # spots present at delta and every probe
            if ok.sum() < 3:
                status = "too_few_spots"
                break
            wt = self._weights(r0)
            A, g = J.T @ (wt[:, None] * J), J.T @ (wt * r0)
            improved = False
            for _ in range(8):
                step = np.linalg.solve(A + lam * np.diag(np.diag(A)), -g)
                trial, tvalid = self._residuals(delta[None] + step[None], meas, idx)
                common = ok & tvalid[0]
                if common.sum() >= 3:
                    r_new = trial[0][common].reshape(-1)
                    r_old = r0.reshape(-1, 3)[common[ok]].reshape(-1)
                    if self._cost(r_new) < self._cost(r_old):
                        delta = delta + step
                        lam = max(lam / 3.0, 1e-9)
                        improved = True
                        n_accepted += 1
                        break
                lam *= 5.0
            if not improved:
                # No damped step lowers the cost: a stationary point of the (quantised)
                # least-squares objective, unless it fails before any step was accepted.
                status = "no_descent"
                break
            if np.linalg.norm(step) < self.tol:
                status = "step_below_tol"
                break
        converged = status == "step_below_tol" or (status == "no_descent" and n_accepted >= 1)
        # covariance at the solution
        r0, J, ok = self._linearize(meas, delta)
        cov = np.full((3, 3), np.nan)
        chi2 = np.nan
        if ok.sum() >= 3:
            try:
                wt = self._weights(r0)
                cov = np.linalg.inv(J.T @ (wt[:, None] * J))
            except np.linalg.LinAlgError:
                pass
            chi2 = float((r0**2).sum() / max(1, r0.size - 3))
        return dict(
            delta=delta,
            cov=cov,
            n_used=int(ok.sum()),
            n_iter=n_iter,
            converged=bool(converged),
            status=status,
            chi2=chi2,
        )
