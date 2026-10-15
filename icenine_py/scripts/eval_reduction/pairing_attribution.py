"""Task D step 3: peak-to-voxel attribution and the discriminating test.

(a) Attribution against the truth. For a voxel, the paired observed peaks (the chosen pairs of the
    main preset) are attached by
      A  radius: the pair's origin (sample frame, frame-centre omega) lies within r of the voxel
         centroid (the plan's rule), r over a grid;
      B  footprint: the ray from the voxel centroid along the pair's observed direction hits the
         bounding box of both spots (tolerance 1 pixel plus the arc motion r x d_omega / 2). This
         is the ray-consistency test of the ray rule applied once per image set. It is needed
         because an observed spot is nearly always the union of the spots of all voxels of a
         grain (step 1), so its origin is the grain footprint's centre, not one voxel.
    Truth: the voxel's true peaks are the (frame, reflection) peaks it drew on both detectors; an
    attached pair is true when the voxel is a source of both spots of the pair, with a shared
    reflection (the triple (spot 0, spot 1, voxel) occurs among the true peaks).
(b) Discriminating test. At the true orientation and at wrong orientations (random axes at fixed
    misorientation, and rotations about low-index crystal axes), three scores:
      pixel   1 - cost of the strict C++ overlap rule, radius 0, all reflections (the final cost)
      attrA / attrB   completeness of the attributed peak list: the share of the orientation's
         predicted observable peaks that match an attributed observed peak (same frame, ray
         direction within `tol` degrees)
    Wrong orientations are random-axis perturbations at 0.05 to 20 degrees and rotations of
    10 to 60 degrees about <100>, <110>, <111>.

  uv run python scripts/eval_reduction/pairing_attribution.py attach     # (a), no images needed
  uv run python scripts/eval_reduction/pairing_attribution.py separate   # (b), full-sample stack
"""

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "common"))
import doc_tables as DT  # noqa: E402
import pairing_common as pc  # noqa: E402
import stats as ST  # noqa: E402

from icenine import spot_pairing as sp  # noqa: E402

MAIN = "ring_tight"
OUT = pc.BENCH / "step3"
RADII_UM = [2, 5, 10, 15, 20, 30, 50, 100]
VOXELS_30 = json.loads(
    (
        HERE.parents[1] / "benchmarks" / "phase_d_seed_diag" / "ray_overlap" / "voxels.json"
    ).read_text()
)
SWEEP_DEG = [0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
LATTICE_DEG = [10, 20, 30, 40, 50, 60]
N_DIR = 8
TOL_DEG = 0.05  # one pixel of centroid error on both detectors over the 2 mm separation
TOLS = (0.02, 0.05, 0.1, 0.2)


def wil(k: int, n: int) -> Dict[str, float]:
    lo, hi = ST.wilson(int(k), int(n)) if n else (float("nan"), float("nan"))
    return {
        "k": int(k),
        "n": int(n),
        "rate": float(k / n) if n else float("nan"),
        "lo": lo,
        "hi": hi,
    }


def in_sorted(x: np.ndarray, sorted_vals: np.ndarray) -> np.ndarray:
    """Membership of x in a sorted unique array (searchsorted, no re-sorting of the big side)."""
    i = np.searchsorted(sorted_vals, x)
    i = np.clip(i, 0, len(sorted_vals) - 1)
    return sorted_vals[i] == x


class Attributor:
    """Chosen pairs of the main preset plus the truth lookups."""

    def __init__(self, preset: str = MAIN) -> None:
        self.geo = pc.load_geometry()
        self.gs = pc.detector_geoms(self.geo)
        prov = pc.load_provenance()
        self.ct = pc.component_table(prov, self.geo)
        self.N = len(self.ct["npix"])
        self.peaks = pc.both_detector_peaks(prov)
        r = np.load(pc.CACHE / f"pairs_full_{preset}.npz")
        ch = r["chosen"].astype(bool)
        self.c0, self.c1 = r["c0"][ch], r["c1"][ch]
        self.frame = r["frame"][ch].astype(np.int64)
        self.origin = r["origin_sample"][ch]
        self.dirn = r["dirn"][ch]
        self.omc = 0.5 * (self.geo["omega_lo"] + self.geo["omega_hi"])
        self.cen = self.geo["voxel_centroid"]
        self.grain = self.geo["voxel_grain"]
        self.side = float(self.geo["voxel_side"])
        # truth triples (c0, c1, voxel)
        V = len(self.cen)
        self.V = V
        tri = (self.peaks["c0"].astype(np.int64) * self.N + self.peaks["c1"]) * V + self.peaks[
            "voxel"
        ]
        self.true_triples = np.unique(tri)
        self.tree = cKDTree(self.cen)
        self.om_step = float(np.mean(self.geo["omega_hi"] - self.geo["omega_lo"]))
        self.pair_key = self.c0 * self.N + self.c1

    def true_peaks_of(self, v: int) -> Dict[str, np.ndarray]:
        m = self.peaks["voxel"] == v
        return {k: x[m] for k, x in self.peaks.items()}

    def attach_radius(self, v: int, r_mm: float) -> np.ndarray:
        d = np.linalg.norm(self.origin - self.cen[v], axis=1)
        return np.nonzero(d <= r_mm)[0]

    def attach_footprint(self, v: int, margin_px: float = 1.0) -> np.ndarray:
        """Pairs whose two spots' bounding boxes contain the voxel centroid projected along the
        pair's direction (at the frame-centre omega)."""
        om = self.omc[self.frame]
        o = np.zeros((len(om), 3))
        o[:, :2] = sp.sample_to_lab_xy(np.broadcast_to(self.cen[v], (len(om), 2)), om)
        rad = float(np.linalg.norm(self.cen[v]))
        margin = margin_px + rad * self.om_step / 2.0 / self.gs[0].pw
        ok = np.ones(len(om), dtype=bool)
        for d, cc in ((0, self.c0), (1, self.c1)):
            g = self.gs[d]
            den = self.dirn @ g.n
            t = -(o @ g.n + g.d) / np.where(np.abs(den) > 1e-12, den, 1e-12)
            pt = o + t[:, None] * self.dirn
            col, row = sp.lab_to_pixel(g, pt)
            ct = self.ct
            ok &= (
                (col >= ct["cmin"][cc] - margin)
                & (col <= ct["cmax"][cc] + 1 + margin)
                & (row >= ct["rmin"][cc] - margin)
                & (row <= ct["rmax"][cc] + 1 + margin)
            )
        return np.nonzero(ok)[0]

    def score_attach(self, v: int, idx: np.ndarray) -> Dict[str, int]:
        tp = self.true_peaks_of(v)
        n_true = len(tp["frame"])
        key = self.c0[idx] * self.N + self.c1[idx]
        trip = key * self.V + v
        is_true = in_sorted(trip, self.true_triples)
        found = np.isin(tp["c0"].astype(np.int64) * self.N + tp["c1"], key)
        return {
            "true_peaks": n_true,
            "found": int(found.sum()),
            "attached": int(len(idx)),
            "attached_true": int(is_true.sum()),
            "attached_false": int((~is_true).sum()),
        }


def boundary_flags(att: Attributor, voxels: List[int]) -> np.ndarray:
    """True when a voxel within 1.6 side lengths belongs to another grain."""
    out = np.zeros(len(voxels), dtype=bool)
    for i, v in enumerate(voxels):
        nb = att.tree.query_ball_point(att.cen[v], 1.6 * att.side)
        out[i] = bool((att.grain[nb] != att.grain[v]).any())
    return out


def cmd_attach(args: argparse.Namespace) -> None:
    att = Attributor()
    OUT.mkdir(parents=True, exist_ok=True)
    spread = list(range(0, att.V, max(att.V // 300, 1)))[:300]
    sets = {"voxels30": VOXELS_30, "spread300": spread}
    summary: Dict[str, Any] = {"preset": MAIN, "radii_um": RADII_UM, "sets": {}}
    tables: Dict[str, str] = {}
    for sname, vox in sets.items():
        bflag = boundary_flags(att, vox)
        res: Dict[str, Any] = {"n_voxels": len(vox), "n_boundary": int(bflag.sum())}
        rows = []
        modes = [(f"A r={r} um", ("A", r)) for r in RADII_UM] + [("B footprint", ("B", 0))]
        for label, (mode, r) in modes:
            tot = {
                "true_peaks": 0,
                "found": 0,
                "attached": 0,
                "attached_true": 0,
                "attached_false": 0,
            }
            per_b = {False: dict(tot), True: dict(tot)}
            for v, bf in zip(vox, bflag):
                idx = att.attach_radius(v, r / 1000.0) if mode == "A" else att.attach_footprint(v)
                s = att.score_attach(v, idx)
                for k in tot:
                    tot[k] += s[k]
                    per_b[bool(bf)][k] += s[k]
            res[label] = {
                "all": tot,
                "interior": per_b[False],
                "boundary": per_b[True],
                "recall": wil(tot["found"], tot["true_peaks"]),
                "precision": wil(tot["attached_true"], tot["attached"]),
                "false_per_voxel": tot["attached_false"] / len(vox),
                "true_attached_per_voxel": tot["attached_true"] / len(vox),
                "recall_interior": wil(per_b[False]["found"], per_b[False]["true_peaks"]),
                "recall_boundary": wil(per_b[True]["found"], per_b[True]["true_peaks"]),
                "false_per_voxel_interior": per_b[False]["attached_false"]
                / max(len(vox) - int(bflag.sum()), 1),
                "false_per_voxel_boundary": per_b[True]["attached_false"]
                / max(int(bflag.sum()), 1),
            }
            rr = res[label]
            rows.append(
                {
                    "attachment": label,
                    "true peaks found": f"{rr['recall']['rate']:.3f} [{rr['recall']['lo']:.3f}, "
                    f"{rr['recall']['hi']:.3f}] ({tot['found']}/{tot['true_peaks']})",
                    "attached pairs per voxel": tot["attached"] / len(vox),
                    "false attached per voxel": rr["false_per_voxel"],
                    "share of attached that are false": (
                        1.0 - rr["precision"]["rate"] if tot["attached"] else float("nan")
                    ),
                    "recall, interior": f"{rr['recall_interior']['rate']:.3f}",
                    "recall, boundary": f"{rr['recall_boundary']['rate']:.3f}",
                }
            )
        summary["sets"][sname] = res
        tables[f"step3_attach_{sname}"] = DT.markdown_table(
            rows,
            list(rows[0].keys()),
            {
                "attached pairs per voxel": ".1f",
                "false attached per voxel": ".1f",
                "share of attached that are false": ".3f",
            },
        )
    summary["pairs_available_peaks_ceiling"] = _ceiling(att)
    (OUT / "attach_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    DT.write_tables(OUT / "tables_attach.md", tables)
    for t in tables.values():
        print(t)
        print()


def _ceiling(att: Attributor) -> Dict[str, Any]:
    """Share of all true two-detector peaks whose spot pair is among the chosen pairs."""
    key = att.peaks["c0"].astype(np.int64) * att.N + att.peaks["c1"]
    found = np.isin(key, att.pair_key)
    return wil(int(found.sum()), len(found))


# ------------------------------------------------------------------------------------------
# (b) predicted peaks and the discriminating test
# ------------------------------------------------------------------------------------------
class Predictor:
    """Observable peaks of an orientation (stages 1 to 4 of the forward simulation plus the
    voxel centre inside both detectors), batched over orientations."""

    def __init__(self, ctx: Any, att: Attributor) -> None:
        import torch

        self.torch = torch
        s = ctx.setup
        self.sim = s.simulator
        self.eta_limit = ctx.eta_limit
        self.lookup = s.range_map.to_lookup_tensor()
        self.att = att
        self.gs = att.gs
        self.cf = ctx.cost_fn("full")
        self.ctx = ctx


def predicted_peaks(pr: Predictor, R: np.ndarray, phase: int, centre: np.ndarray):
    """Observable peaks of ONE orientation R (3,3): frame index (n,) and unit ray directions."""
    from icenine.diffraction_core import get_scattering_omegas_torch
    from icenine.peak_filters import batch_eta_filter

    torch = pr.torch
    g_hkl, g_mag = pr.cf._phase_recip_vecs[phase]
    Rt = torch.from_numpy(np.asarray(R, dtype=np.float32))
    g_lab = (Rt @ g_hkl.T).T
    br = get_scattering_omegas_torch(g_lab, g_mag, pr.sim.beam_energy, pr.sim.beam_deflection_chi)
    om = torch.stack([br.omega1, br.omega2], dim=1)  # (R, 2)
    obs = br.observable[:, None].expand_as(om)
    idx_t, low, width, n_bins = pr.lookup
    f = (om.reshape(-1).double() - low) / width
    b = f.long()
    valid = (f >= 0) & (b < n_bins)
    wedge = idx_t[b.clamp(0, n_bins - 1)]
    keep = obs.reshape(-1) & valid & (wedge >= 0)
    sel = torch.nonzero(keep)[:, 0]
    if len(sel) == 0:
        return np.zeros(0, np.int64), np.zeros((0, 3))
    refl = sel // 2
    w = om.reshape(-1)[sel]
    c, s_ = torch.cos(w), torch.sin(w)
    n = g_lab[refl] / (torch.norm(g_lab[refl], dim=1, keepdim=True) + 1e-10)
    ln = torch.stack(
        [c * n[:, 0] - s_ * n[:, 1], s_ * n[:, 0] + c * n[:, 1], n[:, 2]], dim=1
    )  # Rz(w) @ n, base rotation = identity
    beam = pr.sim.beam_direction.float()
    ray = beam[None, :] - 2.0 * (ln @ beam)[:, None] * ln
    acc, _ = batch_eta_filter(ray, pr.eta_limit, torch.ones(len(sel)), torch.ones(len(sel)))
    ray = ray[acc].numpy().astype(np.float64)
    frame = wedge[sel][acc].numpy().astype(np.int64)
    wv = w[acc].numpy().astype(np.float64)
    ray /= np.linalg.norm(ray, axis=1, keepdims=True)
    # voxel centre must project inside both detectors at the exact omega
    o = np.zeros((len(wv), 3))
    o[:, :2] = sp.sample_to_lab_xy(np.broadcast_to(centre, (len(wv), 2)), wv)
    inside = np.ones(len(wv), dtype=bool)
    for g in pr.gs:
        den = ray @ g.n
        t = -(o @ g.n + g.d) / np.where(np.abs(den) > 1e-12, den, 1e-12)
        col, row = sp.lab_to_pixel(g, o + t[:, None] * ray)
        inside &= (t > 0) & (col >= 0) & (col < 2048) & (row >= 0) & (row < 2048)
    return frame[inside], ray[inside]


def completeness(
    frame: np.ndarray,
    ray: np.ndarray,
    obs_frame: np.ndarray,
    obs_dir: np.ndarray,
    tol_deg: float,
    frame_window: int = 0,
) -> Tuple[float, float, int]:
    """(share of predicted peaks matched, share of observed peaks matched, n predicted).

    A predicted peak matches an observed one when the frames differ by at most ``frame_window``
    and the ray directions differ by at most ``tol_deg`` degrees."""
    n = len(frame)
    if n == 0 or len(obs_frame) == 0:
        return 0.0, 0.0, n
    cos_tol = np.cos(np.radians(tol_deg))
    match_p = np.zeros(n, dtype=bool)
    match_o = np.zeros(len(obs_frame), dtype=bool)
    order = np.argsort(obs_frame, kind="stable")
    of, od = obs_frame[order], obs_dir[order]
    lo = np.searchsorted(of, frame - frame_window, side="left")
    hi = np.searchsorted(of, frame + frame_window, side="right")
    for i in np.nonzero(hi > lo)[0]:
        c = od[lo[i] : hi[i]] @ ray[i]
        m = c >= cos_tol
        if m.any():
            match_p[i] = True
            match_o[order[lo[i] : hi[i]][m]] = True
    return float(match_p.sum() / n), float(match_o.sum() / len(obs_frame)), n


def random_axes(rng: np.random.Generator, n: int) -> np.ndarray:
    a = rng.normal(size=(n, 3))
    return a / np.linalg.norm(a, axis=1, keepdims=True)


def rot(axis: np.ndarray, deg: float) -> np.ndarray:
    a = np.asarray(axis, dtype=np.float64)
    a = a / np.linalg.norm(a)
    t = np.radians(deg)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(t) * K + (1 - np.cos(t)) * (K @ K)


def orientation_set(R_true: np.ndarray, voxel: int) -> List[Tuple[str, float, np.ndarray]]:
    """(group, angle in degrees, R) of the wrong orientations of a voxel."""
    rng = np.random.default_rng(1000 + voxel)
    out: List[Tuple[str, float, np.ndarray]] = []
    for deg in SWEEP_DEG:
        for ax in random_axes(rng, N_DIR):
            out.append(("random axis", deg, rot(ax, deg) @ R_true))
    axes = {"100": [1, 0, 0], "110": [1, 1, 0], "111": [1, 1, 1]}
    for name, ax in axes.items():
        for deg in LATTICE_DEG:
            out.append((f"about <{name}>", float(deg), R_true @ rot(np.array(ax, float), deg)))
    return out


def cmd_separate(args: argparse.Namespace) -> None:
    sys.path.insert(0, str(HERE))
    import parity_common as PC  # noqa: E402

    att = Attributor()
    ctx = PC.PyContext()
    pr = Predictor(ctx, att)
    voxels = VOXELS_30[: args.limit] if args.limit else VOXELS_30
    bflag = boundary_flags(att, voxels)
    A_idx_r = args.r_um / 1000.0
    rows_out: List[Dict[str, Any]] = []
    t0 = time.time()
    for vi, v in enumerate(voxels):
        phase = int(ctx.mic.voxels[v].phase)
        verts = ctx.verts(v)
        R_true = ctx.truth(v)
        centre = att.cen[v]
        sets = {
            "A": att.attach_radius(v, A_idx_r),
            "A30": att.attach_radius(v, 0.030),
            "B": att.attach_footprint(v),
        }
        obs = {k: (att.frame[ix], att.dirn[ix]) for k, ix in sets.items()}
        # oracle: the voxel's own true pairs (what a perfect attribution would give)
        tp = att.true_peaks_of(v)
        key = tp["c0"].astype(np.int64) * att.N + tp["c1"]
        orc = np.nonzero(np.isin(att.pair_key, key))[0]
        obs["oracle"] = (att.frame[orc], att.dirn[orc])
        cases = [("truth", 0.0, R_true)] + orientation_set(R_true, v)
        for group, ang, R in cases:
            q = float(pr.cf.evaluate(np.asarray(R, np.float32), verts, phase).quality)
            fr, ry = predicted_peaks(pr, R, phase, centre)
            rec: Dict[str, Any] = {
                "voxel": v,
                "boundary": bool(bflag[vi]),
                "group": group,
                "angle": ang,
                "pixel": q,
                "n_pred": int(len(fr)),
            }
            for name, (of, od) in obs.items():
                for tol in TOLS:
                    for fw in (0, 1):
                        c, p, _ = completeness(fr, ry, of, od, tol, fw)
                        rec[f"{name}_c_{tol}_{fw}"] = c
                        rec[f"{name}_p_{tol}_{fw}"] = p
            rows_out.append(rec)
        print(f"voxel {v} done ({vi + 1}/{len(voxels)}), {time.time() - t0:.0f}s", flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    (pc.CACHE / "separate_records.json").write_text(json.dumps(rows_out) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("attach").set_defaults(fn=cmd_attach)
    p = sub.add_parser("separate")
    p.add_argument("--limit", type=int, default=0)
    p.add_argument("--r-um", type=float, default=15.7, help="attribution radius A (um)")
    p.set_defaults(fn=cmd_separate)
    args = ap.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
