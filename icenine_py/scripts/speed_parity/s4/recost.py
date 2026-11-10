"""S4 Step 4: the re-cost of the screen passers in orientation-major order.

For each level-0 orientation, the local cost function's peaks are taken once (its own float32
operations, read from one ``evaluate`` call: omegas, normals, frames, plane-hit flags), their
float64 affine maps ``u = A x + c`` built, and every voxel that passed the screen is rasterised
exactly through the library's stage-D C call (triangle mode), with its vertices
``A (x_v + delta_i) + c`` (float64, handed to the C call as float32 like ``evaluate`` does).

* ``equal``: on the 20 Phase D reference voxels, every passer's OverlapInfo against the per-voxel
  ``evaluate`` stored in the Step 0 reference (differences listed with the vertex that crossed a
  pixel boundary).
* ``time`` (run alone): per-orientation set-up and per (orientation, voxel) pair time on the
  2,000-voxel patch of Step 3 for a sample of orientations, against per-voxel ``evaluate``.

Usage (from icenine_py/):
  uv run python scripts/speed_parity/s4/recost.py equal
  uv run python scripts/speed_parity/s4/recost.py time [--n-orient 200]
  uv run python scripts/speed_parity/s4/recost.py memory [--n-orient 200]   (bytes, not timed)
"""

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import common as C  # noqa: E402
import orderings as ORD  # noqa: E402
import translation_correctness as TC  # noqa: E402

from icenine import cost_functions as CF  # noqa: E402

FIELDS = (
    "pixel_overlap",
    "pixel_on_detector",
    "peak_overlap",
    "peak_on_detector",
    "quality",
    "n_quality_points",
    "detectors_overlap",
)


class OrientRecost:
    """Per-orientation peaks of the local cost function and the per-voxel C call."""

    def __init__(self, lcf: Any, geoms: List[C.DetGeom], t_s: np.ndarray) -> None:
        self.lcf = lcf
        self.geoms = geoms
        self.t_s = t_s
        self.spy = TC.Spy()
        self.table = lcf.frame_table()

    def setup(self, R: np.ndarray, verts_any: torch.Tensor, phase: int) -> Dict[str, Any]:
        with self.spy:
            c = self.spy.run(self.lcf, R, verts_any, phase)
        if c.get("verts") is None:
            return dict(M=0)
        om, nr = c["omegas"], c["normals"]
        fr = TC.frame_valid_mask(self.lcf, om)
        keep = fr >= 0
        full, refl = TC.f64_geometry(self.lcf, om[keep], nr[keep])
        A, cc, _, _, _ = C.affine_maps(full, refl, self.t_s, self.geoms)
        args = list(c["args"])
        return dict(M=int(keep.sum()), A=A, c=cc, args=args)

    def voxel(self, st: Dict[str, Any], X64: np.ndarray) -> Tuple[float, ...]:
        """OverlapInfo fields of one voxel (vertices X64 (3, 3) float64)."""
        if st["M"] == 0:
            return (0, 0, 0, 0, 0.0, 0, 0)
        U = np.einsum("pdaj,qj->pdqa", st["A"], X64) + st["c"][:, :, None, :]  # (M, nd, 3, 2)
        args = st["args"]
        args[4] = np.ascontiguousarray(U.astype(np.float32))
        return tuple(CF._c_stage_d_overlap(*args))

    def voxels(self, st: Dict[str, Any], X64s: np.ndarray) -> List[Tuple[float, ...]]:
        """The same for many voxels (X64s (n, 3, 3)): one vertex product, then one C call each."""
        if st["M"] == 0:
            return [(0, 0, 0, 0, 0.0, 0, 0)] * len(X64s)
        U = np.einsum("pdaj,nqj->npdqa", st["A"], X64s) + st["c"][None, :, :, None, :]
        U32 = np.ascontiguousarray(U.astype(np.float32))
        args = st["args"]
        out = []
        for i in range(len(X64s)):
            args[4] = U32[i]
            out.append(tuple(CF._c_stage_d_overlap(*args)))
        return out


def cmd_equal(a: argparse.Namespace) -> None:
    ds = C.load_set("phase_d")
    ids = C.reference_voxels(ds)
    cap = C.capture_level0(ds, ds.voxels[ids[0]])
    lcf = cap["local_cost_fn"]
    R_all = C.all_orientations(cap)
    rc = OrientRecost(lcf, C.det_geoms(lcf), C.sample_translation(lcf))
    refs = {vid: np.load(C.CACHE / "ref" / f"phase_d_{vid}.npz") for vid in ids}
    # orientation-major: orientation -> voxels that passed
    by_o: Dict[int, List[Tuple[int, int]]] = {}
    for vid, z in refs.items():
        for j, k in enumerate(z["passers"]):
            by_o.setdefault(int(k), []).append((vid, j))
    got = {vid: np.zeros_like(z["recost"]) for vid, z in refs.items()}
    X64 = {vid: C.voxel_vertices(ds.voxels[vid]).double().numpy() for vid in ids}
    t0 = time.time()
    v_any = C.voxel_vertices(ds.voxels[ids[0]])
    for n, (k, lst) in enumerate(sorted(by_o.items())):
        st = rc.setup(R_all[k], v_any, ds.voxels[ids[0]].phase)
        for vid, j in lst:
            got[vid][j] = rc.voxel(st, X64[vid])
        if n % 5000 == 0:
            print(f"{n}/{len(by_o)} {time.time() - t0:.0f}s", flush=True)
    rows = []
    for vid, z in refs.items():
        ex, gt = z["recost"], got[vid]
        fi = {f: i for i, f in enumerate(z["fields"])}
        diff = np.abs(ex - gt).max(axis=1) > 0

        def recost(m: np.ndarray) -> np.ndarray:
            pon = m[:, fi["peak_on_detector"]]
            return np.where(pon > 0, 1.0 - m[:, fi["peak_overlap"]] / np.maximum(pon, 1), 1.0)

        rows.append(
            dict(
                vid=int(vid),
                n_passers=int(len(ex)),
                any_field_diff=int(diff.sum()),
                recost_diff=int((recost(ex) != recost(gt)).sum()),
                quality_abs_diff_max=float(
                    np.abs(ex[:, fi["quality"]] - gt[:, fi["quality"]]).max()
                ),
                pixel_overlap_abs_diff_max=float(
                    np.abs(ex[:, fi["pixel_overlap"]] - gt[:, fi["pixel_overlap"]]).max()
                ),
            )
        )
        print(json.dumps(rows[-1]), flush=True)
    # explain a sample of the differing pairs: does an integer vertex differ between evaluate's
    # float32 vertices and the float64 affine vertices?
    rng = np.random.default_rng(2)
    diffs = [
        (vid, j)
        for vid, z in refs.items()
        for j in np.nonzero(np.abs(z["recost"] - got[vid]).max(axis=1) > 0)[0]
    ]
    pick = [diffs[i] for i in rng.choice(len(diffs), min(a.n_explain, len(diffs)), replace=False)]
    explained = 0
    dists: List[float] = []
    spy = TC.Spy()
    for vid, j in pick:
        k = int(refs[vid]["passers"][j])
        v = ds.voxels[vid]
        st = rc.setup(R_all[k], C.voxel_vertices(v), v.phase)
        with spy:
            c = spy.run(lcf, R_all[k], C.voxel_vertices(v), v.phase)
        U = np.einsum("pdaj,qj->pdqa", st["A"], X64[vid]) + st["c"][:, :, None, :]
        U32 = U.astype(np.float32).astype(np.float64)
        flip = TC.trunc_c(U32) != TC.trunc_c(c["verts"])
        if flip.any():
            explained += 1
            dists.append(float(np.max(np.abs(U[flip] - np.round(U[flip])))))
    res = dict(
        explain_sampled=len(pick),
        explain_integer_vertex_differs=explained,
        explain_flip_dist_max=max(dists) if dists else None,
        rows=rows,
        pairs=int(sum(r["n_passers"] for r in rows)),
        any_field_diff=int(sum(r["any_field_diff"] for r in rows)),
        recost_diff=int(sum(r["recost_diff"] for r in rows)),
        quality_abs_diff_max=float(max(r["quality_abs_diff_max"] for r in rows)),
        wall_s_untimed=time.time() - t0,
        orientations=len(by_o),
    )
    (C.OUT / "step4_recost_equal.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps({k: v for k, v in res.items() if k != "rows"}, indent=1))


def cmd_time(a: argparse.Namespace) -> None:
    pf = C.preflight_json(C.OUT / "preflight_step4.json")
    label = C.label_of(pf)
    torch.set_num_threads(1)
    ctx = ORD.Ctx("phase_d")
    ds = ctx.ds
    v2 = ctx.patch(2000)
    pm = ORD.run_kind(ctx, "o2p", v2, 4, "u8")  # pass flags of the patch
    lcf = ctx.cap["local_cost_fn"]
    R_all = C.all_orientations(ctx.cap)
    rng = np.random.default_rng(11)
    sample = np.sort(rng.choice(len(R_all), a.n_orient, replace=False))
    rc = OrientRecost(lcf, C.det_geoms(lcf), C.sample_translation(lcf))
    X64 = [C.voxel_vertices(ds.voxels[int(v)]).double().numpy() for v in v2]
    verts32 = [C.voxel_vertices(ds.voxels[int(v)]) for v in v2]
    phase = int(ds.voxels[0].phase)
    t_setup = t_pairs = 0.0
    n_pairs = 0
    t_eval = 0.0
    t_batch = 0.0
    n_batch_diff = 0
    X64a = np.stack(X64)
    n_eval = 0
    for k in sample:
        vs = np.nonzero(pm[k])[0]
        t0 = time.perf_counter()
        st = rc.setup(R_all[k], verts32[0], phase)
        t1 = time.perf_counter()
        one = [rc.voxel(st, X64[j]) for j in vs]
        t2 = time.perf_counter()
        many = rc.voxels(st, X64a[vs])
        t2b = time.perf_counter()
        n_batch_diff += sum(int(x != y) for x, y in zip(one, many))
        t_batch += t2b - t2
        t_setup += t1 - t0
        t_pairs += t2 - t1
        n_pairs += len(vs)
        # per-voxel evaluate of a few of the same pairs (today's re-cost)
        for j in vs[: a.n_eval_per_orient]:
            t3 = time.perf_counter()
            lcf.evaluate(orientation=R_all[k], voxel_vertices=verts32[j], phase_index=phase)
            t_eval += time.perf_counter() - t3
            n_eval += 1
    res = dict(
        label=label,
        preflight_load=pf.get("loadavg"),
        threads=1,
        n_orient=int(a.n_orient),
        patch=int(len(v2)),
        pairs=int(n_pairs),
        setup_us_per_orientation=1e6 * t_setup / a.n_orient,
        pair_us=1e6 * t_pairs / max(n_pairs, 1),
        evaluate_us=1e6 * t_eval / max(n_eval, 1),
        n_evaluate=int(n_eval),
        pair_us_batched_vertices=1e6 * t_batch / max(n_pairs, 1),
        batched_vs_per_voxel_differences=int(n_batch_diff),
        speedup_pair_vs_evaluate=(t_eval / max(n_eval, 1)) / (t_pairs / max(n_pairs, 1)),
        preflight_end=C.preflight_json(C.OUT / "preflight_step4_end.json").get("loadavg"),
    )
    (C.OUT / "step4_recost_timing.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))


def cmd_memory(a: argparse.Namespace) -> None:
    """Bytes a stored layer-wide local-cost table would hold per orientation (not timed): for the
    same seeded sample of orientations as ``time``, the arrays the per-pair call needs (A, c, the
    wedge and detector-hit flags, a 2-byte frame index per peak), scaled to every orientation."""
    ds = C.load_set("phase_d")
    cap = C.capture_level0(ds, ds.voxels[C.reference_voxels(ds)[0]])
    lcf = cap["local_cost_fn"]
    R_all = C.all_orientations(cap)
    rng = np.random.default_rng(11)
    sample = np.sort(rng.choice(len(R_all), a.n_orient, replace=False))
    rc = OrientRecost(lcf, C.det_geoms(lcf), C.sample_translation(lcf))
    v_any = C.voxel_vertices(ds.voxels[0])
    phase = int(ds.voxels[0].phase)
    Ms: List[int] = []
    nbytes: List[int] = []
    for k in sample:
        st = rc.setup(R_all[k], v_any, phase)
        M = int(st["M"])
        Ms.append(M)
        if M == 0:
            nbytes.append(0)
            continue
        wedge, det_hit = np.asarray(st["args"][1]), np.asarray(st["args"][2])
        nbytes.append(int(st["A"].nbytes + st["c"].nbytes + wedge.nbytes + det_hit.nbytes + 2 * M))
    n_all = len(R_all)
    res = dict(
        n_orient_sampled=int(a.n_orient),
        n_orient_layer=int(n_all),
        peaks_per_orientation_mean=float(np.mean(Ms)),
        peaks_per_orientation_min=int(min(Ms)),
        peaks_per_orientation_max=int(max(Ms)),
        bytes_per_orientation_mean=float(np.mean(nbytes)),
        bytes_per_peak=float(sum(nbytes) / max(sum(Ms), 1)),
        table_bytes_layer_estimate=float(np.mean(nbytes) * n_all),
        dtypes=dict(A=str(st["A"].dtype), c=str(st["c"].dtype)),
    )
    (C.OUT / "step4_recost_memory.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("equal")
    p.add_argument("--n-explain", type=int, default=400)
    p = sub.add_parser("time")
    p.add_argument("--n-orient", type=int, default=200)
    p.add_argument("--n-eval-per-orient", type=int, default=20)
    p = sub.add_parser("memory")
    p.add_argument("--n-orient", type=int, default=200)
    a = ap.parse_args()
    dict(equal=cmd_equal, time=cmd_time, memory=cmd_memory)[a.cmd](a)


if __name__ == "__main__":
    main()
