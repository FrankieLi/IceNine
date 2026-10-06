"""Features of one candidate orientation from ONE pass of the forward model against the detector
images (E2). Per eligible (peak, detector) pair the pass records the reflection, the |q| family,
the detector and whether the spot lands on lit pixels (centre pixel; centre or any of the 3
triangle vertices; +-1 and +-3 pixel boxes around the centre). Features are aggregates of those
hits by |q| family, by detector and "CSL-aware" aggregates over the reflections an orientation
shares with its Sigma relatives.

Pure functions of (orientation, images, voxel); deterministic.
"""

import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import csl  # noqa: E402

SIGMAS = (3, 5, 7, 9, 11)
RADII = (1, 3)


class FeatureExtractor:
    """Built once per worker; `set_image(keys)` per case; `features(R, vertices, phase)` per
    candidate."""

    def __init__(self, local_fn: Any, geo: Any, phase: int):
        self.lf = local_fn
        self.geo = geo
        g_hkl, g_mag = local_fn._phase_recip_vecs[phase]
        self.g_hkl, self.g_mag = g_hkl, g_mag
        qm = np.round(g_mag.numpy().astype(np.float64), 3)
        self.q_levels = np.unique(qm)
        self.fam = np.searchsorted(self.q_levels, qm)  # |q| family index per reflection
        self.n_fam = len(self.q_levels)
        # reflections a candidate shares with each of its Sigma relatives (pure crystal-frame
        # property): shared[sigma] (n_rel, n_refl) bool
        hk = g_hkl.numpy().astype(np.float64)
        self.shared: Dict[int, np.ndarray] = {}
        for sg in SIGMAS:
            rel, _ = csl.csl_relatives(np.eye(3), sigmas=[sg])
            self.shared[sg] = csl.invariant_reflection_mask(hk, np.eye(3), rel)
        self._keys: Optional[np.ndarray] = None

    def set_image(self, keys: np.ndarray) -> None:
        self._keys = np.asarray(keys, dtype=np.int64)  # sorted unique

    # -- geometry ---------------------------------------------------------------------------
    def peak_table(self, R: np.ndarray, vertices: torch.Tensor) -> Dict[str, np.ndarray]:
        """Eligible (peak, detector) pairs of orientation R: arrays refl, det, hit0, hit_any,
        hit1, hit3 (bool), n_pairs."""
        from icenine.diffraction_core import get_scattering_omegas_torch

        lf = self.lf
        sim, sample = lf.simulator, lf.sample
        Rt = torch.from_numpy(np.asarray(R, dtype=np.float32))
        g_lab = (Rt @ self.g_hkl.T).T
        br = get_scattering_omegas_torch(
            g_lab, self.g_mag, sim.beam_energy, sim.beam_deflection_chi
        )
        obs = br.observable
        empty = dict(
            refl=np.zeros(0, dtype=int), det=np.zeros(0, dtype=int),
            hit0=np.zeros(0, bool), hit_any=np.zeros(0, bool),
            hit1=np.zeros(0, bool), hit3=np.zeros(0, bool),
        )  # fmt: skip
        if not obs.any():
            return empty
        obs_idx = torch.nonzero(obs).squeeze(1)
        obs_g = g_lab[obs]
        normals = obs_g / torch.norm(obs_g, dim=1, keepdim=True)
        omegas = torch.cat([br.omega1[obs], br.omega2[obs]])
        nrm = torch.cat([normals, normals])
        refl = torch.cat([obs_idx, obs_idx])
        N = len(omegas)
        c, s = torch.cos(omegas), torch.sin(omegas)
        Rz = torch.zeros(N, 3, 3)
        Rz[:, 0, 0], Rz[:, 0, 1], Rz[:, 1, 0], Rz[:, 1, 1], Rz[:, 2, 2] = c, -s, s, c, 1.0
        base = sample.sample_to_lab_matrix
        full = Rz @ base[:3, :3]
        lab_n = torch.bmm(full, nrm.unsqueeze(-1)).squeeze(-1)
        beam = sim.beam_direction
        dot = (beam.unsqueeze(0) * lab_n).sum(dim=1, keepdim=True)
        refl_dir = beam.unsqueeze(0) - 2.0 * dot * lab_n
        rn = torch.norm(refl_dir, dim=1)
        safe = torch.where(rn > 0, rn, torch.ones_like(rn))
        eta = torch.atan2(torch.abs(refl_dir[:, 1]) / safe, torch.abs(refl_dir[:, 2]) / safe)
        valid = (eta < lf.eta_limit) & (rn > 0)
        if lf.min_sin_eta > 0.0:
            valid = valid & (torch.sin(eta) >= lf.min_sin_eta)
        # omega -> wedge (frame) index
        rm = lf.range_map
        b = ((omegas.numpy() - rm.low) / rm.width).astype(int)
        wedge = np.full(N, -1, dtype=np.int64)
        il = rm.index_list
        for i in np.nonzero(valid.numpy())[0]:
            if 0 <= b[i] < len(il) and il[b[i]] is not None:
                wedge[i] = il[b[i]]
        keep = np.nonzero(wedge >= 0)[0]
        if len(keep) == 0:
            return empty
        kt = torch.from_numpy(keep)
        full4 = torch.zeros(len(keep), 4, 4)
        full4[:, :3, :3] = full[kt]
        full4[:, :3, 3] = base[:3, 3]
        full4[:, 3, 3] = 1.0
        v4 = torch.cat([vertices, torch.ones(3, 1)], dim=1)
        lab_v = torch.einsum("mij,vj->mvi", full4, v4)[:, :, :3]  # (M, 3, 3)
        M = len(keep)
        rd = refl_dir[kt]
        out_r, out_d = [], []
        h0, ha, h1, h3 = [], [], [], []
        assert self._keys is not None
        geo = self.geo
        for d_i, det in enumerate(lf.detector_list):
            plane = det._detector_plane
            origins = lab_v.reshape(M * 3, 3)
            dirs = rd.unsqueeze(1).expand(M, 3, 3).reshape(M * 3, 3)
            denom = (dirs * plane.normal).sum(dim=1)
            numer = -((origins * plane.normal).sum(dim=1) + plane.d)
            par = torch.abs(denom) < 1e-8
            t = torch.where(
                par,
                torch.zeros_like(denom),
                numer / torch.where(par, torch.ones_like(denom), denom),
            )
            hit = ((~par) & (t > 0)).reshape(M, 3).all(dim=1)
            pts = origins + t.unsqueeze(1) * dirs
            rel = pts - det._position
            pl = rel - det._lab_frame_coord_origin
            j = (pl * det._lab_frame_basis_j).sum(dim=1)
            k = (pl * det._lab_frame_basis_k).sum(dim=1)
            col = ((j + det.pixel_half_width) / det.pixel_width).reshape(M, 3).numpy()
            row = ((k + det.pixel_half_height) / det.pixel_height).reshape(M, 3).numpy()
            cc, rc = col.mean(axis=1), row.mean(axis=1)  # centroid
            inb = hit.numpy() & (cc >= 0) & (cc < geo.W) & (rc >= 0) & (rc < geo.H)
            sel = np.nonzero(inb)[0]
            if len(sel) == 0:
                continue
            fr = wedge[keep[sel]]
            img = fr * geo.n_det + d_i

            def lit(cx: np.ndarray, ry: np.ndarray) -> np.ndarray:
                ok = (cx >= 0) & (cx < geo.W) & (ry >= 0) & (ry < geo.H)
                if len(self._keys) == 0:
                    return np.zeros(cx.shape, dtype=bool)
                key = (
                    (img[:, None] * geo.H + np.clip(ry, 0, geo.H - 1)) * geo.W
                    + np.clip(cx, 0, geo.W - 1)
                    if cx.ndim == 2
                    else (img * geo.H + np.clip(ry, 0, geo.H - 1)) * geo.W
                    + np.clip(cx, 0, geo.W - 1)
                )
                pos = np.searchsorted(self._keys, key)
                pos = np.clip(pos, 0, len(self._keys) - 1)
                return ok & (self._keys[pos] == key)

            cx0 = np.floor(cc[sel]).astype(np.int64)
            ry0 = np.floor(rc[sel]).astype(np.int64)
            l0 = lit(cx0, ry0)
            vx = np.floor(col[sel]).astype(np.int64)
            vy = np.floor(row[sel]).astype(np.int64)
            lv = lit(vx, vy)  # (n, 3)
            boxes = []
            for r_ in RADII:
                dx = np.arange(-r_, r_ + 1)
                gx = (cx0[:, None, None] + dx[None, :, None]) + 0 * dx[None, None, :]
                gy = (ry0[:, None, None] + dx[None, None, :]) + 0 * dx[None, :, None]
                l = lit(gx.reshape(len(sel), -1), gy.reshape(len(sel), -1))
                boxes.append(l.any(axis=1))
            out_r.append(refl[keep[sel]].numpy())
            out_d.append(np.full(len(sel), d_i))
            h0.append(l0)
            ha.append(l0 | lv.any(axis=1))
            h1.append(boxes[0])
            h3.append(boxes[1])
        if not out_r:
            return empty
        return dict(
            refl=np.concatenate(out_r), det=np.concatenate(out_d), hit0=np.concatenate(h0),
            hit_any=np.concatenate(ha), hit1=np.concatenate(h1), hit3=np.concatenate(h3),
        )  # fmt: skip

    # -- features ---------------------------------------------------------------------------
    names: List[str] = []

    def feature_names(self) -> List[str]:
        n = ["log_n_pairs", "hit0", "hit_any", "hit1", "hit3"]
        n += [f"fam{f}_{k}" for f in range(self.n_fam) for k in ("hit0", "hit3", "frac")]
        n += [f"det{d}_{k}" for d in range(self.geo.n_det) for k in ("hit0", "hit3", "frac")]
        for sg in SIGMAS:
            n += [
                f"S{sg}_{k}"
                for k in (
                    "shared_frac",
                    "hit0_shared",
                    "hit0_nonshared_min",
                    "hit3_nonshared_min",
                    "hit0_nonshared_mean",
                )
            ]
        n += ["cost_local", "cost_global3"]
        return n

    def features(
        self, R: np.ndarray, vertices: torch.Tensor, phase: int, with_cost: bool = True
    ) -> np.ndarray:
        t = self.peak_table(R, vertices)
        n = len(t["refl"])
        f: List[float] = []
        if n == 0:
            f = [0.0] * (len(self.feature_names()) - 2)
        else:
            h0, ha, h1, h3 = (t[k].astype(float) for k in ("hit0", "hit_any", "hit1", "hit3"))
            f += [np.log1p(n), h0.mean(), ha.mean(), h1.mean(), h3.mean()]
            fam = self.fam[t["refl"]]
            for q in range(self.n_fam):
                m = fam == q
                f += [h0[m].mean() if m.any() else 0.0, h3[m].mean() if m.any() else 0.0, m.mean()]
            for d in range(self.geo.n_det):
                m = t["det"] == d
                f += [h0[m].mean() if m.any() else 0.0, h3[m].mean() if m.any() else 0.0, m.mean()]
            for sg in SIGMAS:
                sh = self.shared[sg][:, t["refl"]]  # (n_rel, n)
                any_sh = sh.any(axis=0)
                f.append(float(any_sh.mean()))
                f.append(float(h0[any_sh].mean()) if any_sh.any() else 0.0)
                ns = (~sh).astype(float)  # per relative: peaks NOT shared with it
                cnt = ns.sum(axis=1)
                rate0 = np.where(cnt > 0, (ns @ h0) / np.maximum(cnt, 1), 1.0)
                rate3 = np.where(cnt > 0, (ns @ h3) / np.maximum(cnt, 1), 1.0)
                f += [float(rate0.min()), float(rate3.min()), float(rate0.mean())]
        if with_cost:
            lf = self.lf
            c_loc = lf.evaluate(np.asarray(R, dtype=np.float32), vertices, phase).cost
            old = lf.pixel_radius
            lf.pixel_radius = 3
            c_g = lf.evaluate(np.asarray(R, dtype=np.float32), vertices, phase).cost
            lf.pixel_radius = old
            f += [float(c_loc), float(c_g)]
        else:
            f += [np.nan, np.nan]
        return np.asarray(f, dtype=np.float64)
