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
from typing import Dict, Optional, Tuple

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
    windows: torch.Tensor, spec: WindowSpec, observer: BatchedObserver
) -> Measurements:
    """Frame index and lit-pixel centroid of every spot in one sample's windows (M, W, W)."""
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
    ):
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

    def solve(self, meas: Measurements, delta0: Optional[np.ndarray] = None) -> Dict[str, object]:
        idx = np.nonzero(meas.used)[0]
        delta = np.zeros(3) if delta0 is None else np.asarray(delta0, dtype=np.float64).copy()
        eye = np.eye(3)
        lam = 1e-3
        converged, n_iter = False, 0
        for n_iter in range(1, self.max_iter + 1):
            pts = np.stack(
                [delta] + [delta + s * self.h * eye[i] for i in range(3) for s in (1, -1)]
            )
            res, valid = self._residuals(pts, meas, idx)
            ok = valid.all(axis=0)  # spots present at the point and at every probe
            if ok.sum() < 3:
                converged = False
                break
            r0 = res[0][ok].reshape(-1)
            J = np.stack(
                [
                    (res[1 + 2 * i][ok] - res[2 + 2 * i][ok]).reshape(-1) / (2 * self.h)
                    for i in range(3)
                ],
                axis=1,
            )
            A, g = J.T @ J, J.T @ r0
            cost0 = 0.5 * float(r0 @ r0)
            improved = False
            for _ in range(8):
                step = np.linalg.solve(A + lam * np.diag(np.diag(A)), -g)
                trial, tvalid = self._residuals(delta[None] + step[None], meas, idx)
                common = ok & tvalid[0]
                if common.sum() >= 3:
                    r_new = trial[0][common].reshape(-1)
                    r_old = res[0][common].reshape(-1)
                    if 0.5 * float(r_new @ r_new) < 0.5 * float(r_old @ r_old):
                        delta = delta + step
                        lam = max(lam / 3.0, 1e-9)
                        improved = True
                        break
                lam *= 5.0
            if not improved or np.linalg.norm(step) < self.tol:
                # No damped step lowers the cost, or the step is negligible: a stationary
                # point of the (quantised) least-squares objective.
                converged = True
                break
        # covariance at the solution
        pts = np.stack([delta] + [delta + s * self.h * eye[i] for i in range(3) for s in (1, -1)])
        res, valid = self._residuals(pts, meas, idx)
        ok = valid.all(axis=0)
        cov = np.full((3, 3), np.nan)
        chi2 = np.nan
        if ok.sum() >= 3:
            J = np.stack(
                [
                    (res[1 + 2 * i][ok] - res[2 + 2 * i][ok]).reshape(-1) / (2 * self.h)
                    for i in range(3)
                ],
                axis=1,
            )
            try:
                cov = np.linalg.inv(J.T @ J)
            except np.linalg.LinAlgError:
                pass
            chi2 = float((res[0][ok] ** 2).sum() / max(1, res[0][ok].size - 3))
        return dict(
            delta=delta,
            cov=cov,
            n_used=int(ok.sum()),
            n_iter=n_iter,
            converged=bool(converged),
            chi2=chi2,
        )
