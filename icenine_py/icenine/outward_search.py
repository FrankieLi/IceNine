"""Searching outward from a start orientation (Task F; library only, not wired in).

BFS neighbours, reseeds and the hybrid finisher all start from an orientation that is a few degrees
from the answer. These two searches are built for that case, not for the whole fundamental zone:

* ``RayGaussNewton`` (F4). Gauss-Newton on the observed-ray residual of Task E's arm E4
  (``smooth_cost.SmoothCost``). Every predicted peak is matched to the nearest observed ray attached
  to the voxel (within ``frame_window`` frames, default 1). The residual of a peak is the
  tangent-plane offset of its predicted scattered direction from that ray (two numbers, the sines
  of the angle's components) plus the distance of its predicted omega to the ray's frame (one
  number); the matching and the line search use the angle itself (atan2), which agrees with the
  tangent-plane residual to first order. The Jacobian with respect to a small rotation is exact
  (autograd through ``SmoothCost.predict``, one vmapped backward pass per residual component, 3m
  in all), so one evaluation gives a full correction. The loss is E4's robust
  ``r^2 / (r^2 + tau^2)``, minimised by iteratively reweighted Gauss-Newton with a
  Levenberg-Marquardt damping and a trust region on the step; ``tau`` can shrink over a schedule.
  The strict cost is not used inside the loop (it is flat far from the answer); it is evaluated at
  the start and at the end, and ``accepted`` says whether the end passes the acceptance test.
  The matching is held fixed inside one Jacobian and redone at every trial point.
* ``ShellSearch`` (F5). Expanding shells on SO(3): deterministic points on spheres of radius
  ``k x spacing`` in rotation-vector space about the start (equal-area points per shell, hexagonal
  packing density). Each shell is scored with a coarse screen (callable, batched), the best
  ``n_smooth`` by the screen are re-ranked with a smooth cost (callable), the best ``n_strict`` by
  that are costed with the strict cost; the search stops at the first strict cost that passes the
  acceptance test. The callables keep this module free of the heavy objects; ``screen_fn_from`` and
  ``smooth_fn_from`` build them from the library's cost classes.

Both count their evaluations by kind (strict, screen, smooth), so a result can be put next to
another method's strict evaluations. A "smooth evaluation" is one E4 value or one E4 value with
Jacobian. Measured single-worker (Task F ``timing/e0.json`` and ``timing/full.json``, one voxel
each), a value with its Jacobian costs about 16 single E4 values on E0 images (10.8 ms against
0.66 ms) and about 8 on full-sample images (12.1 ms against 1.50 ms).

Rotations are 3x3 matrices; a step is applied on the left, ``R_new = exp(delta) R`` (the convention
of ``MCOptimizer`` and Task B's walks). Frame and beam assumptions of E4 hold (lab = Rz(omega)
sample, beam along +x, two detectors; see ``smooth_cost``).
"""

import math
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import torch
from scipy.spatial.transform import Rotation

GOLDEN_ANGLE = math.pi * (3.0 - math.sqrt(5.0))
HEX_AREA = math.sqrt(3.0) / 2.0  # area of a hexagonal-packing cell of unit spacing


# --------------------------------------------------------------------------------------------
# rotations
# --------------------------------------------------------------------------------------------
def rotate_left(R: np.ndarray, rotvec_rad: np.ndarray) -> np.ndarray:
    """``exp(rotvec) @ R`` for rotation vectors ``rotvec_rad`` (3,) or (n, 3); R is (3, 3)."""
    return Rotation.from_rotvec(np.asarray(rotvec_rad, dtype=np.float64)).as_matrix() @ R


def angle_deg(R_a: np.ndarray, R_b: np.ndarray) -> float:
    """Rotation angle (deg) of ``R_a R_b^T`` (no crystal-symmetry reduction)."""
    c = (np.trace(R_a @ R_b.T) - 1.0) / 2.0
    return float(np.degrees(np.arccos(np.clip(c, -1.0, 1.0))))


def rodrigues_torch(delta: torch.Tensor) -> torch.Tensor:
    """Rotation matrix of a rotation vector (3,) with a series for the exponential, so that the
    value and the first derivative are exact at ``delta = 0`` (no division by the angle)."""
    zero = torch.zeros((), dtype=delta.dtype)
    K = torch.stack(
        [
            torch.stack([zero, -delta[2], delta[1]]),
            torch.stack([delta[2], zero, -delta[0]]),
            torch.stack([-delta[1], delta[0], zero]),
        ]
    )
    eye = torch.eye(3, dtype=delta.dtype)
    K2 = K @ K
    return eye + K + K2 / 2.0 + K2 @ K / 6.0 + K2 @ K2 / 24.0


# --------------------------------------------------------------------------------------------
# shells (F5)
# --------------------------------------------------------------------------------------------
def shell_size(radius_deg: float, spacing_deg: float) -> int:
    """Points on the shell of ``radius_deg``: the sphere area (in rotation-vector degrees, flat
    metric) over the area of one hexagonal-packing cell of ``spacing_deg``; at least 6."""
    area = 4.0 * math.pi * radius_deg * radius_deg
    return max(6, int(math.ceil(area / (HEX_AREA * spacing_deg * spacing_deg))))


def shell_directions(n: int, shell_index: int = 0) -> np.ndarray:
    """``n`` unit vectors, equal-area (Fibonacci spiral). The spiral is turned by a different
    amount on each shell (``shell_index``) so successive shells do not line up."""
    i = np.arange(n, dtype=np.float64)
    z = 1.0 - 2.0 * (i + 0.5) / n
    phi = GOLDEN_ANGLE * i + shell_index * 1.0
    rho = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    return np.stack([rho * np.cos(phi), rho * np.sin(phi), z], axis=1)


def shell_rotations(
    R0: np.ndarray, radius_deg: float, spacing_deg: float, shell_index: int = 0
) -> np.ndarray:
    """The shell of rotations of rotation angle ``radius_deg`` about ``R0``: (n, 3, 3)."""
    n = shell_size(radius_deg, spacing_deg)
    rv = shell_directions(n, shell_index) * math.radians(radius_deg)
    return rotate_left(R0, rv)


@dataclass
class ShellResult:
    """Outcome of ``ShellSearch.search``. ``radius_deg`` is the radius of the shell that gave the
    accepted candidate (0: the start itself, NaN: none accepted)."""

    R: np.ndarray
    strict_cost: float
    found: bool
    radius_deg: float
    n_strict: int
    n_screen: int
    n_smooth: int
    n_shells: int
    # per shell: (radius deg, best strict cost so far, evaluations so far, best rotation so far)
    trace: List[Tuple[float, float, int, np.ndarray]] = field(default_factory=list)

    @property
    def n_evals(self) -> int:
        return self.n_strict + self.n_screen + self.n_smooth


class ShellSearch:
    """F5: expanding shells about a start, screen then smooth then strict, stop on acceptance.

    ``screen_fn(Rs) -> score`` (higher = better, one number per rotation, batched),
    ``smooth_fn(Rs) -> cost`` (lower = better, may be None: the screen order is then used),
    ``strict_fn(R) -> cost`` (lower = better, one rotation). Ties in the screen score are broken
    by the order of the points on the shell (stable sort), which is deterministic."""

    def __init__(
        self,
        screen_fn: Callable[[np.ndarray], np.ndarray],
        strict_fn: Callable[[np.ndarray], float],
        smooth_fn: Optional[Callable[[np.ndarray], np.ndarray]] = None,
        spacing_deg: float = 0.3,
        r_max_deg: float = 6.0,
        accept_cost: float = 0.8,
        n_smooth: int = 32,
        n_strict: int = 2,
        screen_chunk: int = 256,
    ) -> None:
        if spacing_deg <= 0 or r_max_deg <= 0 or n_strict < 1 or n_smooth < n_strict:
            raise ValueError("ShellSearch: need spacing > 0, r_max > 0, n_smooth >= n_strict >= 1")
        self.screen_fn, self.strict_fn, self.smooth_fn = screen_fn, strict_fn, smooth_fn
        self.spacing_deg = float(spacing_deg)
        self.r_max_deg = float(r_max_deg)
        self.accept_cost = float(accept_cost)
        self.n_smooth = int(n_smooth)
        self.n_strict = int(n_strict)
        self.screen_chunk = int(screen_chunk)

    def n_shells_max(self) -> int:
        return int(math.floor(self.r_max_deg / self.spacing_deg + 1e-9))

    def search(self, R0: np.ndarray) -> ShellResult:
        R0 = np.asarray(R0, dtype=np.float64)
        c0 = float(self.strict_fn(R0))
        res = ShellResult(
            R0, c0, c0 <= self.accept_cost, 0.0 if c0 <= self.accept_cost else math.nan, 1, 0, 0, 0
        )
        if res.found:
            return res
        best_c, best_R = c0, R0
        for k in range(1, self.n_shells_max() + 1):
            r = k * self.spacing_deg
            Rs = shell_rotations(R0, r, self.spacing_deg, k)
            score = np.concatenate(
                [
                    np.asarray(self.screen_fn(Rs[i : i + self.screen_chunk]), dtype=np.float64)
                    for i in range(0, len(Rs), self.screen_chunk)
                ]
            )
            res.n_screen += len(Rs)
            res.n_shells = k
            top = np.argsort(-score, kind="stable")[: self.n_smooth]
            if self.smooth_fn is not None:
                sm = np.asarray(self.smooth_fn(Rs[top]), dtype=np.float64)
                res.n_smooth += len(top)
                top = top[np.argsort(sm, kind="stable")]
            for j in top[: self.n_strict]:
                c = float(self.strict_fn(Rs[j]))
                res.n_strict += 1
                if c < best_c:
                    best_c, best_R = c, Rs[j]
                if c <= self.accept_cost:
                    res.R, res.strict_cost, res.found, res.radius_deg = Rs[j], c, True, r
                    res.trace.append((r, best_c, res.n_evals, best_R))
                    return res
            res.trace.append((r, best_c, res.n_evals, best_R))
        res.R, res.strict_cost = best_R, best_c
        return res


def screen_fn_from(screen_cost: Any, voxel_vertices: Any, phase: int = 0) -> Callable:
    """Batched screen score (number of peaks with an overlap within the pixel radius; higher is
    better) from a pixel-radius ``VoxelCostFunction`` through ``batched_screen.ScreenEvaluator``."""
    from .batched_screen import ScreenEvaluator

    ev = ScreenEvaluator(screen_cost)

    def fn(Rs: np.ndarray) -> np.ndarray:
        po, _ = ev.overlaps_exact(np.asarray(Rs, dtype=np.float64), voxel_vertices, phase)
        screen_cost.eval_count += len(Rs)
        return np.asarray(po, dtype=np.float64)

    return fn


def smooth_fn_from(smooth: Any, voxel_vertices: Any, phase: int = 0, spec: str = "E4,tau=0.25"):
    """Smooth cost (lower is better) of a batch of rotations from a ``SmoothCost``."""
    from .smooth_cost import SmoothSpec

    sp = SmoothSpec.parse(spec)

    def fn(Rs: np.ndarray) -> np.ndarray:
        return smooth.evaluate_many(np.asarray(Rs), voxel_vertices, phase, spec=sp)

    return fn


# --------------------------------------------------------------------------------------------
# Gauss-Newton on the observed-ray residual (F4)
# --------------------------------------------------------------------------------------------
@dataclass
class RayGNResult:
    R: np.ndarray
    strict_cost: float
    start_strict_cost: float
    accepted: bool
    e4_start: float
    e4_end: float
    iterations: int
    n_jac: int  # E4 value-and-Jacobian evaluations
    n_val: int  # E4 value evaluations (trial points)
    n_strict: int
    n_matched: int  # peaks matched to a ray at the last linearisation
    converged: bool
    status: str
    history: List[Tuple[int, float, float]] = field(default_factory=list)  # (iter, tau deg, E4)

    @property
    def n_smooth(self) -> int:
        return self.n_jac + self.n_val

    @property
    def n_evals(self) -> int:
        return self.n_smooth + self.n_strict


class RayGaussNewton:
    """F4: Gauss-Newton on the observed-ray residual, trust region, strict cost for acceptance.

    ``smooth`` is a ``SmoothCost`` of the image set (its observed rays are built on first use),
    ``strict_fn(R) -> cost`` the strict cost (evaluated twice per solve: start and end).
    ``tau_schedule_deg`` is the sequence of robust scales; each stage runs to convergence (a step
    below ``step_tol_deg`` or ``max_iter`` iterations) before the next. ``trust_deg`` is the
    largest step of one iteration, the first step is limited to ``trust0_deg``. A predicted peak is
    matched only to observed rays whose frame is within ``frame_window`` frames of its own (1, as in
    ``SmoothCost`` E4); a start far from the truth moves peaks by more frames than that."""

    def __init__(
        self,
        smooth: Any,
        voxel_vertices: Any,
        strict_fn: Callable[[np.ndarray], float],
        phase: int = 0,
        tau_schedule_deg: Sequence[float] = (0.25,),
        max_iter: int = 30,
        trust_deg: float = 2.0,
        trust0_deg: float = 1.0,
        step_tol_deg: float = 5e-4,
        accept_cost: float = 0.8,
        frame_window: int = 1,
    ) -> None:
        self.sc = smooth
        self.verts = voxel_vertices
        self.strict_fn = strict_fn
        self.phase = int(phase)
        self.taus = [math.radians(t) for t in tau_schedule_deg]
        if not self.taus or min(self.taus) <= 0:
            raise ValueError("RayGaussNewton: tau_schedule_deg must be positive")
        self.max_iter = int(max_iter)
        self.trust = math.radians(trust_deg)
        self.trust0 = math.radians(trust0_deg)
        self.step_tol = math.radians(step_tol_deg)
        self.accept_cost = float(accept_cost)
        self.frame_window = int(frame_window)
        if self.frame_window < 1:
            raise ValueError("RayGaussNewton: frame_window must be at least 1")
        self._ray_tables: Optional[Dict[str, Any]] = None

    # ---------------------------------------------------------------- rays of the voxel
    def _tables(self) -> Dict[str, Any]:
        """Direction, tangent basis and omega interval of every ray attached to the voxel."""
        if self._ray_tables is not None:
            return self._ray_tables
        sc = self.sc
        rays = sc.observed_rays()
        att = sc.attached_rays(self.verts)
        d = rays.dirn[att].astype(np.float64)
        helper = np.where(
            np.abs(d[:, :1]) > 0.9, np.array([[0.0, 1.0, 0.0]]), np.array([[1.0, 0.0, 0.0]])
        )
        t1 = np.cross(helper, d)
        t1 /= np.maximum(np.linalg.norm(t1, axis=1, keepdims=True), 1e-300)
        t2 = np.cross(d, t1)
        lo_bin = np.array(
            [np.nonzero(sc._frame_index[:-1] == w)[0].min() for w in range(sc.n_frames)]
        )
        hi_bin = np.array(
            [np.nonzero(sc._frame_index[:-1] == w)[0].max() + 1 for w in range(sc.n_frames)]
        )
        fr = rays.frame[att]
        self._ray_tables = dict(
            n=len(att),
            d=torch.as_tensor(d),
            t1=torch.as_tensor(t1),
            t2=torch.as_tensor(t2),
            frame=torch.as_tensor(fr),
            lo=torch.as_tensor(sc.low + sc.width * lo_bin[fr], dtype=torch.float64),
            hi=torch.as_tensor(sc.low + sc.width * hi_bin[fr], dtype=torch.float64),
        )
        return self._ray_tables

    # ---------------------------------------------------------------- matching and cost
    def _weights(self, pr: Any) -> torch.Tensor:
        """Per-peak weight: 1 for a peak seen on every detector and with a frame, else 0 (E4)."""
        self._tables()  # builds the observed rays, which sets the detector size used below
        col, row, hit = pr.col[:, :, -1], pr.row[:, :, -1], pr.hit[:, :, -1]
        inb = hit & (col >= 0) & (col < self.sc.cols) & (row >= 0) & (row < self.sc.rows)
        return (inb.all(1) & (pr.frame >= 0)).to(torch.float64)

    def _match(self, pr: Any) -> Tuple[torch.Tensor, torch.Tensor]:
        """(index of the nearest attached ray per peak, squared residual to it); the index is -1
        and the residual +inf when no ray lies within one frame."""
        tb = self._tables()
        n = int(pr.b.shape[0])
        if tb["n"] == 0 or n == 0:
            return torch.full((n,), -1, dtype=torch.long), torch.full(
                (n,), math.inf, dtype=torch.float64
            )
        rs, om = pr.ray.detach(), pr.omega.detach()
        cosang = torch.clamp(rs @ tb["d"].T, -1.0, 1.0)
        cross = torch.linalg.norm(
            torch.cross(
                rs[:, None, :].expand(-1, tb["n"], -1),
                tb["d"][None].expand(n, -1, -1),
                dim=-1,
            ),
            dim=-1,
        )
        ang = torch.atan2(cross, cosang)
        dom = torch.clamp(
            torch.maximum(tb["lo"][None] - om[:, None], om[:, None] - tb["hi"][None]), min=0.0
        )
        r2 = ang * ang + dom * dom
        near = (tb["frame"][None] - pr.frame[:, None]).abs() <= self.frame_window
        r2 = torch.where(near, r2, torch.full_like(r2, math.inf))
        r2min, arg = r2.min(dim=1)
        arg = torch.where(torch.isfinite(r2min), arg, torch.full_like(arg, -1))
        return arg, r2min

    def e4_value(self, R: np.ndarray, tau_rad: float) -> float:
        """E4 cost of ``R`` (the mean over the peaks of ``r^2 / (r^2 + tau^2)``, 1 for a peak
        without a ray), the same definition as ``SmoothCost.e4_cost_torch``."""
        with torch.no_grad():
            pr = self.sc.predict(
                torch.as_tensor(R, dtype=torch.float64)[None], self.verts, self.phase
            )
            W = self._weights(pr)
            _, r2 = self._match(pr)
            rho = torch.where(torch.isfinite(r2), r2 / (r2 + tau_rad**2), torch.ones_like(r2))
            sw = W.sum()
            return float((W * rho).sum() / sw) if float(sw) > 0 else 1.0

    # ---------------------------------------------------------------- linearisation
    def linearize(self, R: np.ndarray, tau_rad: float) -> Dict[str, Any]:
        """Residual vector, exact Jacobian and weights at ``R`` for a left rotation ``delta``:
        ``res`` (3m,), ``J`` (3m, 3), ``w`` (3m,) the reweighting of the robust loss, ``cost`` the
        E4 value, ``m`` the number of matched peaks."""
        tb = self._tables()
        R0 = torch.as_tensor(R, dtype=torch.float64)
        delta = torch.zeros(3, dtype=torch.float64, requires_grad=True)
        pr = self.sc.predict((rodrigues_torch(delta) @ R0)[None], self.verts, self.phase)
        W = self._weights(pr).detach()
        arg, r2 = self._match(pr)
        sw = float(W.sum())
        rho_all = torch.where(torch.isfinite(r2), r2 / (r2 + tau_rad**2), torch.ones_like(r2))
        cost = float((W * rho_all).sum() / sw) if sw > 0 else 1.0
        use = (W > 0) & (arg >= 0)
        m = int(use.sum())
        if m == 0:
            return dict(m=0, cost=cost, res=None, J=None, w=None)
        a = arg[use]
        ray = pr.ray[use]
        om = pr.omega[use]
        ra = (ray * tb["t1"][a]).sum(-1)
        rb = (ray * tb["t2"][a]).sum(-1)
        lo, hi = tb["lo"][a], tb["hi"][a]
        ro = torch.where(om > hi, om - hi, torch.where(om < lo, om - lo, torch.zeros_like(om)))
        res = torch.stack([ra, rb, ro], dim=1).reshape(-1)  # (3m,)
        eye = torch.eye(res.numel(), dtype=torch.float64)
        (J,) = torch.autograd.grad(res, delta, grad_outputs=eye, is_grads_batched=True)
        r2m = (res.detach().reshape(-1, 3) ** 2).sum(1)
        wpeak = tau_rad**2 / (r2m + tau_rad**2) ** 2  # d rho / d r^2 of the robust loss
        w = wpeak[:, None].expand(-1, 3).reshape(-1)
        return dict(m=m, cost=cost, res=res.detach(), J=J.detach(), w=w)

    # ---------------------------------------------------------------- solve
    def solve(
        self, R0: np.ndarray, callback: Optional[Callable[[np.ndarray, int, int], None]] = None
    ) -> RayGNResult:
        """Solve from ``R0``. ``callback(R, n_smooth, n_strict)`` is called at the start and after
        every accepted step (a diagnostic hook: the benchmark reads the error to the truth)."""
        R = np.asarray(R0, dtype=np.float64)
        n_jac = n_val = 0
        c_start = float(self.strict_fn(R))
        e4_0 = self.e4_value(R, self.taus[-1])
        n_val += 1
        if callback is not None:
            callback(R, n_jac + n_val, 1)
        it = 0
        hist: List[Tuple[int, float, float]] = []
        status = "max_iter"
        converged = False
        n_match = 0
        trust = self.trust0
        for tau in self.taus:
            lam = 1e-3
            cost = self.e4_value(R, tau)
            n_val += 1
            stage_done = False
            for _ in range(self.max_iter):
                lin = self.linearize(R, tau)
                n_jac += 1
                it += 1
                n_match = lin["m"]
                if lin["m"] < 3:
                    status = "too_few_peaks"
                    stage_done = True
                    break
                J, r, w = lin["J"], lin["res"], lin["w"]
                H = J.T @ (w[:, None] * J)
                g = J.T @ (w * r)
                step = None
                for _try in range(12):
                    A = (
                        H
                        + lam * torch.diag(torch.diag(H))
                        + 1e-14 * torch.eye(3, dtype=torch.float64)
                    )
                    d = -torch.linalg.solve(A, g).numpy()
                    nd = float(np.linalg.norm(d))
                    if nd > trust:
                        d *= trust / nd
                        nd = trust
                    if nd < self.step_tol:
                        step = ("small", d)
                        break
                    R_try = rotate_left(R, d)
                    c_try = self.e4_value(R_try, tau)
                    n_val += 1
                    if c_try < cost - 1e-12:
                        step = ("ok", d, R_try, c_try)
                        break
                    lam *= 4.0
                    trust = max(0.5 * trust, self.step_tol)
                if step is None or step[0] == "small":
                    stage_done = True
                    converged = True
                    break
                _, d, R, cost = step
                if callback is not None:
                    callback(R, n_jac + n_val, 1)
                lam = max(lam / 3.0, 1e-7)
                trust = min(2.0 * trust, self.trust)
                hist.append((it, math.degrees(tau), cost))
            if not stage_done:
                status = "max_iter"
            elif status == "too_few_peaks":
                break
            else:
                status = "converged"
        c_end = float(self.strict_fn(R))
        e4_1 = self.e4_value(R, self.taus[-1])
        n_val += 1
        if c_end <= c_start:
            R_out, c_out = R, c_end
        else:
            R_out, c_out = np.asarray(R0, dtype=np.float64), c_start
        return RayGNResult(
            R=R_out,
            strict_cost=c_out,
            start_strict_cost=c_start,
            accepted=c_out <= self.accept_cost,
            e4_start=e4_0,
            e4_end=e4_1,
            iterations=it,
            n_jac=n_jac,
            n_val=n_val,
            n_strict=2,
            n_matched=n_match,
            converged=converged and status == "converged",
            status=status,
            history=hist,
        )
