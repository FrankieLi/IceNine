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
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import csl  # noqa: E402

SIGMAS = (3, 5, 7, 9, 11)


def _have_c_stage_d() -> bool:
    """True if the C stage-D extension of the cost function is available."""
    from icenine import cost_functions

    return bool(getattr(cost_functions, "_HAS_C_RASTERIZE", False))


RADII = (1, 3)


def feature_names(n_fam: int, n_det: int) -> List[str]:
    """Names of the feature vector columns, in `FeatureExtractor.features` order."""
    n = ["log_n_pairs", "hit0", "hit_any", "hit1", "hit3"]
    n += [f"fam{f}_{k}" for f in range(n_fam) for k in ("hit0", "hit3", "frac")]
    n += [f"det{d}_{k}" for d in range(n_det) for k in ("hit0", "hit3", "frac")]
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


def feature_names_for_width(width: int, n_det: int = 2) -> List[str]:
    """Feature names for a stored feature matrix of `width` columns (the study's geometry has
    2 detectors; the number of |q| families follows from the width)."""
    n_fam = (width - len(feature_names(0, n_det))) // 3
    names = feature_names(n_fam, n_det)
    assert len(names) == width, (len(names), width)
    return names


class FeatureExtractor:
    """Built once per worker; `set_image(keys)` per case; `features(R, vertices, phase)` per
    candidate."""

    def __init__(self, local_fn: Any, geo: Any, phase: int, q_max: Optional[float] = None):
        """q_max (default None = all reflections, the E2 behaviour): restrict the extractor to the
        reflections with |q| <= q_max (Angstrom^-1). The peak table is then the full table
        restricted to those reflections (`refl_full` maps its `refl` indices to the full list), the
        aggregates (families <= q_max, detectors, CSL-aware) are those of the restricted peaks, and
        the two cost columns are the local (pixel radius 0) and pixel-radius-3 costs of cost
        functions built with max_q = q_max, so a low-Q feature vector needs no high-|q| work."""
        self.lf = local_fn
        self.geo = geo
        self.q_max = q_max
        g_hkl, g_mag = local_fn._phase_recip_vecs[phase]
        if q_max is None:
            self.refl_full = np.arange(len(g_mag))
        else:
            self.refl_full = np.nonzero(g_mag.numpy() <= q_max)[0]
            if len(self.refl_full) == 0:
                raise ValueError(f"no reflection with |q| <= {q_max}")
            sel = torch.from_numpy(self.refl_full)
            g_hkl, g_mag = g_hkl[sel], g_mag[sel]
        self.g_hkl, self.g_mag = g_hkl, g_mag
        self._lowq_costs: Optional[Tuple[Any, Any]] = None
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
    def _low_q_cost_fns(self) -> Tuple[Any, Any]:
        """Cost functions restricted to |q| <= q_max at pixel radius 0 and 3, sharing the images of
        the local cost function (re-read at every call, so `attach_images` is followed)."""
        from icenine.cost_functions import VoxelCostFunction

        lf = self.lf
        if self._lowq_costs is None:
            self._lowq_costs = tuple(  # type: ignore[assignment]
                VoxelCostFunction(
                    simulator=lf.simulator, detector_list=lf.detector_list,
                    range_map=lf.range_map, exp_data=lf.exp_data, sample=lf.sample,
                    structure_list=lf.structure_list, mode="hard", eta_limit=lf.eta_limit,
                    pixel_radius=pr, max_q=float(self.q_max), min_sin_eta=lf.min_sin_eta,
                )
                for pr in (0, 3)
            )  # fmt: skip
        assert self._lowq_costs is not None
        for fn in self._lowq_costs:
            fn.exp_data = lf.exp_data
        return self._lowq_costs

    def feature_names(self) -> List[str]:
        return feature_names(self.n_fam, self.geo.n_det)

    # -- batched pass -----------------------------------------------------------------------
    def features_batch(
        self,
        Rs: np.ndarray,
        vertices: torch.Tensor,
        phase: int,
        with_cost: bool = True,
        chunk: int = 256,
    ) -> np.ndarray:
        """`features` of many candidates at once: (B, 3, 3) -> (B, n_features). Equal to stacking
        `features(R, ...)` (the per-candidate path stays the reference; the equality is tested).
        The geometry of all candidates and peaks is one vectorised pass (shared by the two cost
        columns, which the per-candidate path recomputes twice); the aggregates are integer counts
        per candidate (exact in float64), so they equal the per-candidate means. The cost columns
        need the C stage-D extension; without it this falls back to the per-candidate loop."""
        Rs = np.asarray(Rs, dtype=np.float64).reshape(-1, 3, 3)
        if len(Rs) == 0:
            return np.zeros((0, len(self.feature_names())))
        if with_cost and not _have_c_stage_d():
            return np.stack([self.features(R, vertices, phase, True) for R in Rs])
        if len(Rs) > chunk:
            return np.concatenate(
                [
                    self.features_batch(Rs[i : i + chunk], vertices, phase, with_cost, chunk)
                    for i in range(0, len(Rs), chunk)
                ]
            )
        B = len(Rs)
        g = self._batch_geometry(Rs, vertices, with_cost)
        X = np.zeros((B, len(self.feature_names())))
        X[:, -2:] = np.nan
        self._fill_aggregates(X, g["c"], g["refl"], g["det"], g["hits"], B)
        if with_cost:
            lf = self.lf
            if self.q_max is not None:
                radii = (0, 3)
                for fn in self._low_q_cost_fns():
                    fn.eval_count += B
            else:
                radii = (lf.pixel_radius, 3)
                lf.eval_count += 2 * B
            X[:, -2], X[:, -1] = self._batch_costs(g, radii, B)
        return X

    def _wedge_lookup(self) -> np.ndarray:
        rm = self.lf.range_map
        if getattr(self, "_il_src", None) is not rm.index_list:
            self._il_arr = np.array(
                [-1 if w is None else int(w) for w in rm.index_list], dtype=np.int64
            )
            self._il_src = rm.index_list
        return self._il_arr

    def _batch_geometry(
        self, Rs: np.ndarray, vertices: torch.Tensor, with_cost: bool
    ) -> Dict[str, Any]:
        """Eligible (candidate, peak, detector) rows of all candidates; the arithmetic of
        `peak_table`, with the rows of each candidate in the cost function's order (omega1 of its
        observable reflections, then omega2) because the stage-D quality is a running mean."""
        from icenine.diffraction_core import get_scattering_omegas_torch

        lf = self.lf
        sim, sample = lf.simulator, lf.sample
        B, K = len(Rs), len(self.g_mag)
        Rt = torch.from_numpy(Rs.astype(np.float32))
        # one 2-D matmul per candidate: a batched matmul picks another kernel and differs in the
        # last bit, which can flip a pixel (the per-candidate path is the reference)
        gT = self.g_hkl.T
        g_lab = torch.cat([Rt[i] @ gT for i in range(B)], dim=1).T  # (B*K, 3), strided like
        # the per-candidate (R @ g.T).T: a contiguous copy differs in the last bit in the omegas
        br = get_scattering_omegas_torch(
            g_lab, self.g_mag.repeat(B), sim.beam_energy, sim.beam_deflection_chi
        )
        obs = br.observable
        none: Dict[str, Any] = dict(
            c=np.zeros(0, dtype=np.int64), refl=np.zeros(0, dtype=np.int64),
            det=np.zeros(0, dtype=np.int64), hits=np.zeros((4, 0), bool), cost=None,
        )  # fmt: skip
        if not obs.any():
            return none
        flat = torch.nonzero(obs).squeeze(1)
        cand_o = (flat // K).numpy()
        obs_g = g_lab[obs]
        normals = obs_g / torch.norm(obs_g, dim=1, keepdim=True)
        perm = np.argsort(np.concatenate([cand_o * 2, cand_o * 2 + 1]), kind="stable")
        pt = torch.from_numpy(perm)
        omegas = torch.cat([br.omega1[obs], br.omega2[obs]])[pt]
        nrm = torch.cat([normals, normals])[pt]
        refl = torch.cat([flat % K, flat % K])[pt]
        cand = np.concatenate([cand_o, cand_o])[perm]
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
        rm = lf.range_map
        il = self._wedge_lookup()
        b = ((omegas.numpy() - rm.low) / rm.width).astype(int)
        inr = valid.numpy() & (b >= 0) & (b < len(il))
        wedge = np.where(inr, il[np.clip(b, 0, len(il) - 1)], -1)
        keep = np.nonzero(wedge >= 0)[0]
        if len(keep) == 0:
            return none
        kt = torch.from_numpy(keep)
        M = len(keep)
        full4 = torch.zeros(M, 4, 4)
        full4[:, :3, :3] = full[kt]
        full4[:, :3, 3] = base[:3, 3]
        full4[:, 3, 3] = 1.0
        v4 = torch.cat([vertices, torch.ones(3, 1)], dim=1)
        lab_v = torch.einsum("mij,vj->mvi", full4, v4)[:, :, :3]
        rd = refl_dir[kt]
        k_cand, k_refl, k_wedge = cand[keep], refl[kt].numpy(), wedge[keep]
        geo = self.geo
        assert self._keys is not None
        keys = self._keys
        rows_c, rows_r, rows_d, rows_h = [], [], [], []
        pix, allhit = [], []
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
            if with_cost:
                pix.append(np.stack([col, row], axis=2))  # (M, 3, 2) float32
                allhit.append(hit.numpy())
            cc, rc = col.mean(axis=1), row.mean(axis=1)
            inb = hit.numpy() & (cc >= 0) & (cc < geo.W) & (rc >= 0) & (rc < geo.H)
            sel = np.nonzero(inb)[0]
            if len(sel) == 0:
                continue
            img = k_wedge[sel] * geo.n_det + d_i

            def lit(cx: np.ndarray, ry: np.ndarray, img: np.ndarray = img) -> np.ndarray:
                ok = (cx >= 0) & (cx < geo.W) & (ry >= 0) & (ry < geo.H)
                if len(keys) == 0:
                    return np.zeros(cx.shape, dtype=bool)
                im = img[:, None] if cx.ndim == 2 else img
                key = (im * geo.H + np.clip(ry, 0, geo.H - 1)) * geo.W + np.clip(cx, 0, geo.W - 1)
                pos = np.clip(np.searchsorted(keys, key), 0, len(keys) - 1)
                return ok & (keys[pos] == key)

            cx0 = np.floor(cc[sel]).astype(np.int64)
            ry0 = np.floor(rc[sel]).astype(np.int64)
            l0 = lit(cx0, ry0)
            lv = lit(np.floor(col[sel]).astype(np.int64), np.floor(row[sel]).astype(np.int64))
            bx = []
            for r_ in RADII:
                dx = np.arange(-r_, r_ + 1)
                gx = (cx0[:, None, None] + dx[None, :, None]) + 0 * dx[None, None, :]
                gy = (ry0[:, None, None] + dx[None, None, :]) + 0 * dx[None, :, None]
                bx.append(lit(gx.reshape(len(sel), -1), gy.reshape(len(sel), -1)).any(axis=1))
            rows_c.append(k_cand[sel])
            rows_r.append(k_refl[sel])
            rows_d.append(np.full(len(sel), d_i))
            rows_h.append(np.stack([l0, l0 | lv.any(axis=1), bx[0], bx[1]]))
        out = dict(none)
        if rows_c:
            out.update(
                c=np.concatenate(rows_c), refl=np.concatenate(rows_r),
                det=np.concatenate(rows_d), hits=np.concatenate(rows_h, axis=1),
            )  # fmt: skip
        if with_cost:
            out["cost"] = dict(
                cand=k_cand,
                wedge=k_wedge,
                pix=np.stack(pix, axis=1),
                allhit=np.stack(allhit, axis=1),
            )  # pix (M, n_det, 3, 2), allhit (M, n_det)
        return out

    def _batch_costs(
        self, g: Dict[str, Any], radii: Tuple[int, int], B: int
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Costs (1 - quality) of the stage-D overlap at the two pixel radii from the geometry
        already computed; equal to `VoxelCostFunction.evaluate(...).cost` per candidate."""
        from icenine.cost_functions import _c_stage_d_overlap

        cd = g["cost"]
        out = np.ones((2, B))  # OverlapInfo() (no peak): quality 0, cost 1
        if cd is None:
            return out[0], out[1]
        n_det = self.geo.n_det
        exp = self.lf.exp_data
        wedge = cd["wedge"].astype(np.int32)
        uw = np.unique(wedge)
        to_local = np.full(int(uw.max()) + 1, -1, dtype=np.int32)
        to_local[uw] = np.arange(len(uw), dtype=np.int32)
        local = to_local[wedge]
        images = [exp.get_image(int(w), d).get_binary_numpy() for w in uw for d in range(n_det)]
        first = exp.get_image(int(uw[0]), 0)
        det_hit = np.ascontiguousarray(cd["allhit"].astype(np.uint8))
        pix = cd["pix"]
        centers = np.ascontiguousarray(pix[:, :, 0, :].astype(np.int32))  # first vertex
        verts = np.ascontiguousarray(pix.astype(np.float32))
        bounds = np.searchsorted(cd["cand"], np.arange(B + 1))
        for b in range(B):
            s, e = int(bounds[b]), int(bounds[b + 1])
            if e == s:
                continue
            for i, pr in enumerate(radii):
                if pr > 0:
                    res = _c_stage_d_overlap(
                        images, local[s:e], det_hit[s:e], centers[s:e], None,
                        n_det, pr, e - s, first.num_rows, first.num_cols,
                    )  # fmt: skip
                else:
                    res = _c_stage_d_overlap(
                        images, local[s:e], det_hit[s:e], None, verts[s:e],
                        n_det, pr, e - s, first.num_rows, first.num_cols,
                    )  # fmt: skip
                out[i, b] = 1.0 - res[4]
        return out[0], out[1]

    def _fill_aggregates(
        self,
        X: np.ndarray,
        c: np.ndarray,
        refl: np.ndarray,
        det: np.ndarray,
        hits: np.ndarray,
        B: int,
    ) -> None:
        """Columns 0 .. -3 of X from integer counts per candidate (see `features`); candidates
        without any eligible pair keep zeros."""
        n = np.bincount(c, minlength=B)
        has = n > 0
        nf = np.maximum(n, 1).astype(np.float64)

        def cnt(w: np.ndarray, idx: np.ndarray, size: int) -> np.ndarray:
            return np.bincount(idx, weights=w, minlength=size)

        h0, ha, h1, h3 = (hits[i].astype(np.float64) for i in range(4))
        cols: List[np.ndarray] = [np.log1p(n.astype(np.float64))]
        for h in (h0, ha, h1, h3):
            cols.append(cnt(h, c, B) / nf)
        for key, nk in ((self.fam[refl], self.n_fam), (det, self.geo.n_det)):
            idx = c * nk + key
            m = cnt(np.ones(len(c)), idx, B * nk).reshape(B, nk)
            s0 = cnt(h0, idx, B * nk).reshape(B, nk)
            s3 = cnt(h3, idx, B * nk).reshape(B, nk)
            mm = np.maximum(m, 1)
            for q in range(nk):
                cols += [
                    np.where(m[:, q] > 0, s0[:, q] / mm[:, q], 0.0),
                    np.where(m[:, q] > 0, s3[:, q] / mm[:, q], 0.0),
                    m[:, q] / nf,
                ]
        for sg in SIGMAS:
            sh = self.shared[sg][:, refl]  # (n_rel, P)
            any_sh = sh.any(axis=0)
            na = cnt(any_sh.astype(np.float64), c, B)
            cols.append(na / nf)
            cols.append(np.where(na > 0, cnt(h0 * any_sh, c, B) / np.maximum(na, 1), 0.0))
            n_rel = sh.shape[0]
            r0 = np.empty((B, n_rel))
            r3 = np.empty((B, n_rel))
            for r in range(n_rel):
                ns = (~sh[r]).astype(np.float64)
                k_ = cnt(ns, c, B)
                r0[:, r] = np.where(k_ > 0, cnt(ns * h0, c, B) / np.maximum(k_, 1), 1.0)
                r3[:, r] = np.where(k_ > 0, cnt(ns * h3, c, B) / np.maximum(k_, 1), 1.0)
            cols += [r0.min(axis=1), r3.min(axis=1), r0.mean(axis=1)]
        F = np.stack(cols, axis=1)
        X[has, : F.shape[1]] = F[has]

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
        if with_cost and self.q_max is not None:
            c0, c3 = self._low_q_cost_fns()
            Rf = np.asarray(R, dtype=np.float32)
            f += [float(c0.evaluate(Rf, vertices, phase).cost)]
            f += [float(c3.evaluate(Rf, vertices, phase).cost)]
        elif with_cost:
            lf = self.lf
            c_loc = lf.evaluate(np.asarray(R, dtype=np.float32), vertices, phase).cost
            old = lf.pixel_radius
            lf.pixel_radius = 3
            try:
                c_g = lf.evaluate(np.asarray(R, dtype=np.float32), vertices, phase).cost
            finally:
                lf.pixel_radius = old
            f += [float(c_loc), float(c_g)]
        else:
            f += [np.nan, np.nan]
        return np.asarray(f, dtype=np.float64)
