"""AR4 step 1: shared set-up and the evaluation of one pair with its per-peak shares.

The sample is chosen with the environment variable ``RESEED_SAMPLE`` (``big``, ``grad`` or
``mid``), as in the scripts of AR3 steps 8 and 9 and AR5 step 1. The stored runs of S8 are read
through links in ``benchmarks/autoresearch/ar4/cache/stored/<sample>/`` (untracked; they point
into the worktrees of AR3 steps 8 and 9) and are never written.

The library's cost function returns the aggregate counts of one evaluation only. The per-peak
pixel counts are taken here, with no change to the library: ``peak_rows`` wraps the stage-D call
of ``cost_functions.calculate_diffraction_overlap_batched`` and calls the same C function once
per predicted peak (the counts of a peak do not depend on the other peaks). ``check_sum_back``
sums them back to the aggregate counts of the library's own evaluation.
"""

import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]  # icenine_py
AR4 = ICE / "benchmarks" / "autoresearch" / "ar4"
CACHE = AR4 / "cache"
SAMPLE = os.environ.get("RESEED_SAMPLE", "grad")
os.environ.setdefault("RESEED_CACHE", str(CACHE / "stored" / SAMPLE))
sys.path.insert(0, str(ICE / "scripts" / "autoresearch" / "ar3"))
sys.path.insert(0, str(ICE / "scripts" / "common"))
sys.path.insert(0, str(HERE))

import reseed_common as RC  # noqa: E402,F401  (sets AR3_SAMPLE from RESEED_SAMPLE)

import numpy as np  # noqa: E402
import verify_lib as VL  # noqa: E402

STORED = CACHE / "stored" / SAMPLE
S8_DIR = STORED / f"S8_{SAMPLE}"
PASS_DIR = STORED / f"passes_{SAMPLE}"
TRIES_DIR = PASS_DIR / "D4_S8"
N = 24570
TRAIN = ("big", "grad")
TEST = "mid"
# per-peak columns of a pair's rows
PEAK_COLS = ("rank", "dist", "omega", "share", "pixels")


def errors_path(sample: str) -> Path:
    base = ICE / "benchmarks" / "autoresearch" / "ar3"
    if sample == "mid":
        return base / "start" / "errors_mid.npz"
    return base / "reseed" / "slim" / f"errors_{sample}.npz"


class Evaluator:
    """One environment (reconstructor, truth, strict cost function) and its pair evaluations."""

    def __init__(self) -> None:
        self.env = RC.make_env()
        env = self.env
        self.rec = env.bfs.reconstructor
        self.cf = env.hooks.cost_fn
        import icenine.cost_functions as CF
        from icenine import ray_indexing as RI

        self.CF, self.RI = CF, RI
        self.mic = env.mic
        self._ladder: Dict[int, np.ndarray] = {}
        self._att: Dict[int, Any] = {}
        self.n_evals = 0

    # ---- per-peak rows
    def ladder(self, phase: int) -> np.ndarray:
        if phase not in self._ladder:
            _, g_mag = self.cf._phase_recip_vecs[phase]
            self._ladder[phase] = np.unique(np.round(g_mag.numpy().astype(np.float64), 4))
        return self._ladder[phase]

    def peak_rows(self, R: np.ndarray, v: int) -> Dict[str, Any]:
        """One evaluation with the per-peak pixel counts.

        Returns the library's aggregate (``hit``, ``cost``, ``conf``, the counts) from this very
        evaluation, and ``rows`` (n, 5) float32 with columns ``PEAK_COLS`` for every predicted peak
        with at least one pixel on a detector (those peaks make the hit ratio's denominator)."""
        CF, cf = self.CF, self.cf
        phase = self.mic.voxels[v].phase
        vertices = self.env.hooks.vertices(v)
        cap: Dict[str, Any] = {"stage": []}
        orig_bragg = CF.scattering_omegas_fixed_beam
        orig_calc = CF.calculate_diffraction_overlap_batched
        orig_stage = CF._c_stage_d_overlap

        def bragg(g_lab: Any, consts: Any) -> Any:
            res = orig_bragg(g_lab, consts)
            cap["g_lab"], cap["obs"] = g_lab, res.observable
            return res

        def stage(*a: Any) -> Any:
            image_list, wedge, det_hit, centers, verts, n_det, prad, M, rows, cols, *rest = a
            per = []
            for i in range(M):
                r = orig_stage(
                    image_list,
                    wedge[i : i + 1],
                    det_hit[i : i + 1],
                    None if centers is None else centers[i : i + 1],
                    None if verts is None else verts[i : i + 1],
                    n_det,
                    prad,
                    1,
                    rows,
                    cols,
                    *[x[i : i + 1] if hasattr(x, "__len__") else x for x in rest],
                )
                per.append((r[0], r[1], r[2]))
            cap["stage"] = per
            cap["wedge"] = np.asarray(wedge)
            if centers is not None:
                pos = np.asarray(centers)[:, 0, :].astype(np.float64)
            else:
                pos = np.asarray(verts)[:, 0, :, :].astype(np.float64).mean(axis=1)
            cap["pos"], cap["shape"] = pos, (rows, cols)
            return orig_stage(*a)

        def calc(**kw: Any) -> Any:
            cap["omegas"] = kw["peak_omegas"]
            cap["normals"] = kw["peak_normals"]
            cap["kw"] = kw
            return orig_calc(**kw)

        # the pair's score: the library's own (fused) evaluation, as the stored runs and AR5 use it
        info_f = cf.evaluate(np.asarray(R, np.float64), vertices, phase)
        fused = cf.fused_cost
        cf.fused_cost = False  # the per-peak breakdown needs the stage-D path of the same data
        CF.scattering_omegas_fixed_beam = bragg
        CF.calculate_diffraction_overlap_batched = calc
        CF._c_stage_d_overlap = stage
        try:
            info = cf.evaluate(np.asarray(R, np.float64), vertices, phase)
        finally:
            cf.fused_cost = fused
            CF.scattering_omegas_fixed_beam = orig_bragg
            CF.calculate_diffraction_overlap_batched = orig_calc
            CF._c_stage_d_overlap = orig_stage
        self.n_evals += 1
        hit, conf = self.env.bfs._ratios(info_f)
        out: Dict[str, Any] = dict(
            hit=float(hit),
            cost=float(info_f.cost),
            conf=float(conf),
            counts=(info_f.pixel_overlap, info_f.pixel_on_detector, info_f.peak_overlap),
            counts_stage_d=(info.pixel_overlap, info.pixel_on_detector, info.peak_overlap),
            rows=np.zeros((0, len(PEAK_COLS)), np.float32),
        )
        if not cap["stage"]:
            return out
        # which stage-D rows: the rows of ``peak_omegas`` with a frame (as the library selects them)
        kw = cap["kw"]
        fi = kw.get("frame_index")
        if fi is None:
            fi = CF.frame_index_array(kw["range_map"])
        om = kw["peak_omegas"]
        om_np = om.numpy() if hasattr(om, "numpy") else np.asarray(om)
        wedge_all = CF.frame_of_omegas(
            om_np, kw["range_map"].low, kw["range_map"].width, fi, kw.get("cpp_edge_rules", False)
        )
        rows_ok = np.where(wedge_all >= 0)[0]
        assert len(rows_ok) == len(cap["stage"]), (len(rows_ok), len(cap["stage"]))
        # which reflection: match the peak's lab normal to the observable reflections' normals
        g_lab = cap["g_lab"]
        obs_idx = np.where(cap["obs"].numpy())[0]
        unit = g_lab[obs_idx] / np.linalg.norm(g_lab[obs_idx].numpy(), axis=1, keepdims=True)
        nrm = kw["peak_normals"][rows_ok]
        k = obs_idx[np.argmax((nrm @ unit.T).numpy(), axis=1)]
        _, g_mag = cf._phase_recip_vecs[phase]
        rank = VL.family_rank_in(g_mag.numpy()[k], self.ladder(phase))
        st = np.array(cap["stage"], dtype=np.float64)  # (M, 3): po, pod, peak_overlap
        rr, cc = cap["shape"]
        pos = cap["pos"]
        dist = np.hypot((pos[:, 0] - cc / 2) / cc, (pos[:, 1] - rr / 2) / rr)
        has = st[:, 1] > 0
        share = np.where(has, st[:, 0] / np.maximum(st[:, 1], 1), 0.0)
        rows = np.stack([rank, dist, om_np[rows_ok], share, st[:, 1]], axis=1).astype(np.float32)[
            has
        ]
        out["rows"] = rows
        out["counts_rows"] = (int(st[:, 0].sum()), int(st[:, 1].sum()), int(st[:, 2].sum()))
        return out

    # ---- votes (AR5 step 1's number)
    def rays(self, v: int) -> Tuple[Any, Any, Any, Any]:
        if v not in self._att:
            phase = self.mic.voxels[v].phase
            sc, crystal, params, _ = self.rec._index_context(phase)
            att = sc.attached_rays(self.env.hooks.vertices(v))
            g = self.RI.scattering_vectors(sc, att)
            self._att[v] = (g, crystal, params, len(att))
            if len(self._att) > 64:
                self._att.pop(next(iter(self._att)))
        return self._att[v]

    def explained_and_votes_cheap(self, R: np.ndarray, v: int) -> int:
        """Attached rays the orientation explains (no indexing of the rest)."""
        g, crystal, params, _n = self.rays(v)
        ex = self.RI.explained_mask(
            np.asarray(R, np.float64), g, crystal, params.ang_tol_deg, params.len_tol
        )
        return int(ex.sum())

    def explained_and_votes(self, R: np.ndarray, v: int) -> Tuple[int, int]:
        """(attached rays the orientation explains, best vote count among the rest).

        The pilot's function only: it timed the vote count (``best_votes``) on a few voxels of
        ``big``; the pilot's output was not kept. No run of the step calls it; the inputs use
        ``explained_and_votes_cheap``, and ``best_votes`` was not computed."""
        g, crystal, params, _n = self.rays(v)
        ex = self.RI.explained_mask(
            np.asarray(R, np.float64), g, crystal, params.ang_tol_deg, params.len_tol
        )
        gu = g[~ex]
        res = self.RI.index_rays(gu, crystal, params, 20)
        votes = int(np.max(res.votes)) if len(res.votes) else 0
        return int(ex.sum()), votes
