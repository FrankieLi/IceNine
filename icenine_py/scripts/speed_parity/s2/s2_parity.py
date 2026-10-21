"""S2 Step 0d: the fused prototype against the torch ``evaluate``, per intermediate and per
``OverlapInfo`` field, with every differing orientation traced to the peaks that differ.

For each dataset (default: ThreeVoxels voxel 0, S1's ``three0``) and each covered mode, on S1's
orientation set (``s1_parity.orientation_set``: 500 within 3 degrees of the truth + 500 FZ points):

* the library ``VoxelCostFunction.evaluate``;
* ``torch_mirror.mirror`` (the same torch operations with the intermediates kept); its result must
  equal the library's on every orientation, else the run stops;
* ``s2fused.FusedEval.debug`` (one variant per ``--variants``: ``portable`` C, or ``accelerate``,
  g_lab from Accelerate's sgemm as torch on macOS).

Per intermediate (g_lab, omega1/2 of the Bragg solve, the row omega, eta, the frame, the pixel
coordinates of every vertex, the ray parameter t) the difference in float32 ulp over the rows both
sides computed (rows on which stage D runs, for the pixel coordinate and t). Per field the number of
orientations whose value differs. Each differing orientation is traced: the rows whose stage-D input
(present or not, frame, detector hit, the integer pixel coordinates stage D uses) differs, the cause
(observability, eta limit, min_sin_eta, frame edge, pixel integer, t > 0), and the distance of the
torch value from that threshold in ulp. Traced = every such row within ``--trace-ulp`` (8) ulp.

  uv run python scripts/speed_parity/s2/s2_parity.py --out benchmarks/s2_step0/parity.json
  uv run python scripts/speed_parity/s2/s2_parity.py --datasets three0 full15901 --n-near 50 ...
"""

import argparse
import contextlib
import io
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
sys.path.insert(0, str(ICE))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))

import s1_parity as S1  # noqa: E402
import s2fused as SF  # noqa: E402
import torch_mirror as TM  # noqa: E402

# covered modes (S1 names; screen_default is the screen with the default rule)
MODES: Dict[str, Dict[str, Any]] = {
    "local": dict(pixel_radius=0, max_q=8.0),
    "cpp_rule": dict(pixel_radius=0, max_q=8.0, cpp_overlap_rule=True),
    "edge": dict(pixel_radius=0, max_q=8.0, cpp_edge_rules=True),
    "edge_cpp_rule": dict(pixel_radius=0, max_q=8.0, cpp_edge_rules=True, cpp_overlap_rule=True),
    "min_sin_eta": dict(pixel_radius=0, max_q=8.0, min_sin_eta=0.2),
    "screen_default": dict(pixel_radius=3, max_q=5.0, center_screen=True),
    "screen_center": dict(pixel_radius=3, max_q=5.0, cpp_overlap_rule=True, center_screen=True),
    "screen_vertex0": dict(pixel_radius=3, max_q=5.0, cpp_overlap_rule=True, center_screen=False),
    "screen_edge": dict(pixel_radius=3, max_q=5.0, cpp_overlap_rule=True, cpp_edge_rules=True),
}
INTERMEDIATES = ("g_lab", "omega_bragg", "omega", "eta", "pixel", "t")


def ordered(a: np.ndarray) -> np.ndarray:
    """float32 -> int64 whose difference is the distance in ulp (sign-magnitude to two's)."""
    i = np.asarray(a, np.float32).view(np.int32).astype(np.int64)
    return np.where(i < 0, -(i & 0x7FFFFFFF), i)


def ulp_diff(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.abs(ordered(a) - ordered(b))


class Stat:
    """Accumulated ulp differences of one intermediate."""

    def __init__(self) -> None:
        self.n = 0
        self.equal = 0
        self.max_ulp = 0
        self.max_abs = 0.0
        self.hist = np.zeros(5, np.int64)  # 0, 1, 2-8, 9-64, > 64 ulp

    def add(self, a: np.ndarray, b: np.ndarray) -> None:
        a = np.asarray(a, np.float32).reshape(-1)
        b = np.asarray(b, np.float32).reshape(-1)
        if a.size == 0:
            return
        u = ulp_diff(a, b)
        self.n += int(a.size)
        self.equal += int((u == 0).sum())
        self.max_ulp = max(self.max_ulp, int(u.max()))
        self.max_abs = max(self.max_abs, float(np.abs(a.astype(np.float64) - b).max()))
        self.hist += np.array(
            [
                (u == 0).sum(),
                (u == 1).sum(),
                ((u >= 2) & (u <= 8)).sum(),
                ((u >= 9) & (u <= 64)).sum(),
                (u > 64).sum(),
            ]
        )

    def as_dict(self) -> Dict[str, Any]:
        return dict(
            n=self.n,
            equal=self.equal,
            max_ulp=self.max_ulp,
            max_abs=self.max_abs,
            hist_0_1_2to8_9to64_gt64=[int(x) for x in self.hist],
        )


class ErrStat:
    """Absolute error against the float64 reference: count, mean, max."""

    def __init__(self) -> None:
        self.n = 0
        self.total = 0.0
        self.max = 0.0

    def add(self, a: np.ndarray, ref: np.ndarray) -> None:
        e = np.abs(np.asarray(a, np.float64).reshape(-1) - np.asarray(ref).reshape(-1))
        if e.size:
            self.n += int(e.size)
            self.total += float(e.sum())
            self.max = max(self.max, float(e.max()))

    def as_dict(self) -> Dict[str, Any]:
        return dict(n=self.n, mean=self.total / max(self.n, 1), max=self.max)


def stage_d_ints(x: Dict[str, Any], r: int, d: int, radius: int) -> Tuple[int, ...]:
    """The integer pixel coordinates stage D uses for row r on detector d."""
    p = x["pix"][r, d].astype(np.float64)
    if radius > 0:
        return (
            int(np.float32(p[0, 0]).astype(np.int32)),
            int(np.float32(p[0, 1]).astype(np.int32)),
        )
    return tuple(int(-1 if v < 0 else np.floor(v)) for v in p.reshape(-1))


def trace(
    m: Dict[str, Any], d: Dict[str, Any], ref: Dict[str, Any], cf: Any, n_det: int
) -> List[Dict[str, Any]]:
    """Rows whose stage-D input differs between torch (m) and fused (d): the cause, the distance
    of the torch value from the threshold in ulp (``margin_ulp``; for pixel peaks also the fused
    value's, ``fused_margin_ulp``, and both in px, ``*_px_from_edge``), the torch, fused and float64
    reference values and the threshold (``ref_side``: the side of the threshold the float64 value
    is on, "torch" or "fused")."""
    out: List[Dict[str, Any]] = []
    K = m["K"]
    radius = cf.pixel_radius
    eta_lim = np.float32(cf.eta_limit)
    rmap = cf.range_map

    def side(tv: float, fv: float, rv: float, thr: float) -> str:
        return "torch" if (rv - thr >= 0) == (tv - thr >= 0) else "fused"

    for r in range(2 * K):
        k = r % K
        ft, ff = int(m["frame"][r]), int(d["frame"][r])
        if ft < 0 and ff < 0:
            continue
        if (ft < 0) != (ff < 0) or ft != ff:
            if m["obs"][k] != d["obs"][k]:
                e = dict(cause="observable", torch=float(m["bragg_num"][k] - m["bragg_den"][k]))
                e["margin_ulp"] = int(ulp_diff(m["bragg_num"][k], m["bragg_den"][k]))
                e["ref"] = float(ref["obs_margin"][k])
            elif m["valid"][r] != d["valid"][r]:
                if cf.min_sin_eta > 0 and (m["eta"][r] < eta_lim) == (d["eta"][r] < eta_lim):
                    thr = float(np.float32(cf.min_sin_eta))
                    e = dict(cause="min_sin_eta", torch=float(m["sin_eta"][r]), threshold=thr)
                    e["margin_ulp"] = int(ulp_diff(m["sin_eta"][r], np.float32(thr)))
                    e["ref"] = float(np.sin(ref["eta"][r]))
                else:
                    thr = float(eta_lim)
                    e = dict(cause="eta_limit", torch=float(m["eta"][r]), threshold=thr)
                    e["fused"] = float(d["eta"][r])
                    e["margin_ulp"] = int(ulp_diff(m["eta"][r], eta_lim))
                    e["ref"] = float(ref["eta"][r])
            else:
                om = float(m["omega"][r])
                n = round((om - float(rmap.low)) / float(rmap.width))
                thr = float(rmap.low) + n * float(rmap.width)
                e = dict(cause="frame_edge", torch=om, fused=float(d["omega"][r]), threshold=thr)
                e["margin_ulp"] = int(ulp_diff(np.float32(om), np.float32(thr)))
                e["ref"] = float(ref["omega"][r])
            if "threshold" in e and "fused" in e:
                e["ref_side"] = side(e["torch"], e["fused"], e["ref"], e["threshold"])
            e["row"] = r
            out.append(e)
            continue
        for det in range(n_det):
            ht, hf = int(m["hit"][r, det]), int(d["hit"][r, det])
            if ht != hf:
                tt = m["t"][r, det]
                out.append(
                    dict(
                        row=r,
                        det=det,
                        cause="t>0",
                        margin_ulp=int(np.min(ulp_diff(tt, 0.0))),
                        torch=float(np.min(tt)),
                        ref=float(np.min(ref["t"][r, det])),
                    )
                )
                continue
            if not ht:
                continue
            it, iff = stage_d_ints(m, r, det, radius), stage_d_ints(d, r, det, radius)
            if it == iff:
                continue
            tv = m["pix"][r, det].reshape(-1)
            fv = d["pix"][r, det].reshape(-1)
            rv = ref["pix"][r, det].reshape(-1)
            worst: Optional[Dict[str, Any]] = None
            for i in [i for i in range(len(it)) if it[i] != iff[i]]:
                v = np.float32(tv[i])
                if cf.cpp_edge_rules and (it[i] == -1) != (iff[i] == -1) and abs(float(v)) <= 1.0:
                    # -1 rule (d < 0): the threshold is the detector coordinate 0
                    jk = m["jk"][r, det].reshape(-1)[i]
                    e = dict(cause="pixel_negative", margin_ulp=int(ulp_diff(jk, 0.0)))
                    e.update(torch=float(jk))
                else:
                    thr = float(np.round(float(v)))
                    e = dict(cause="pixel_integer", margin_ulp=int(ulp_diff(v, np.float32(thr))))
                    e.update(torch=float(v), fused=float(fv[i]), ref=float(rv[i]), threshold=thr)
                    e["ref_side"] = side(e["torch"], e["fused"], e["ref"], thr)
                    e["torch_px_from_edge"] = abs(float(v) - thr)
                    e["fused_margin_ulp"] = int(ulp_diff(np.float32(fv[i]), np.float32(thr)))
                    e["fused_px_from_edge"] = abs(float(fv[i]) - thr)
                if worst is None or e["margin_ulp"] > worst["margin_ulp"]:
                    worst = e
            assert worst is not None
            worst.update(row=r, det=det)
            out.append(worst)
    return out


def dirty_files() -> Optional[List[str]]:
    """Tracked files that differ from HEAD (``git status --porcelain``, untracked not listed), so a
    run on a dirty tree names what was modified; None if git cannot say."""
    try:
        out = subprocess.check_output(
            ["git", "status", "--porcelain", "--untracked-files=no"], cwd=S1.ROOT, text=True
        )
    except Exception:  # noqa: BLE001
        return None
    return [line[3:] for line in out.splitlines() if line.strip()]


def run_dataset(
    ds: Any, R: np.ndarray, modes: List[str], variants: List[str], trace_ulp: int
) -> Dict[str, Any]:
    res: Dict[str, Any] = {}
    n_det = len(ds.setup.detector_list)
    for mode in modes:
        kw = MODES[mode]
        cf_lib = ds.cost_fn(kw)
        cf_mir = ds.cost_fn(kw)
        fes = {v: SF.FusedEval(ds.cost_fn(kw), v) for v in variants}
        stats = {v: {k: Stat() for k in INTERMEDIATES} for v in variants}
        frame_diff = {v: 0 for v in variants}
        field_diff = {v: {f: 0 for f in SF.FIELDS} for v in variants}
        differing: Dict[str, List[Dict[str, Any]]] = {v: [] for v in variants}
        silent = {v: 0 for v in variants}
        err = {x: dict(pixel=ErrStat(), omega=ErrStat()) for x in ["torch", *variants]}
        t0 = time.perf_counter()
        for i in range(len(R)):
            with contextlib.redirect_stdout(io.StringIO()):
                lib = SF.result_tuple(cf_lib.evaluate(R[i], ds.verts, ds.phase))
                m = TM.mirror(cf_mir, R[i], ds.verts, ds.phase)
            if tuple(m["result"]) != lib:
                raise SystemExit(f"mirror differs from the library: {mode} {i} {m['result']} {lib}")
            ref = TM.reference64(cf_mir, R[i], ds.verts, ds.phase)
            ft = m["frame"] >= 0
            err["torch"]["pixel"].add(m["pix"][ft], ref["pix"][ft])
            err["torch"]["omega"].add(m["omega"][ft], ref["omega"][ft])
            for v in variants:
                d = fes[v].debug(R[i], ds.verts, ds.phase)
                st = stats[v]
                ob = (m["obs"] == 1) & (d["obs"] == 1)
                st["g_lab"].add(m["g_lab"], d["g_lab"])
                st["omega_bragg"].add(m["omega1"][ob], d["omega1"][ob])
                st["omega_bragg"].add(m["omega2"][ob], d["omega2"][ob])
                both = (m["valid"] == 1) & (d["valid"] == 1)
                st["omega"].add(m["omega"][both], d["omega"][both])
                st["eta"].add(m["eta"][both], d["eta"][both])
                fb = (m["frame"] >= 0) & (d["frame"] >= 0)
                frame_diff[v] += int((m["frame"][fb] != d["frame"][fb]).sum())
                st["pixel"].add(m["pix"][fb], d["pix"][fb])
                st["t"].add(m["t"][fb], d["t"][fb])
                ff = d["frame"] >= 0
                err[v]["pixel"].add(d["pix"][ff], ref["pix"][ff])
                err[v]["omega"].add(d["omega"][ff], ref["omega"][ff])
                res_f = tuple(d["result"])
                rows = trace(m, d, ref, cf_mir, n_det)
                if res_f != lib:
                    for j, f in enumerate(SF.FIELDS):
                        field_diff[v][f] += int(res_f[j] != lib[j])
                    differing[v].append(
                        dict(
                            index=i,
                            torch=list(lib),
                            fused=list(res_f),
                            rows=rows,
                            traced=bool(rows) and all(x["margin_ulp"] <= trace_ulp for x in rows),
                        )
                    )
                elif rows:
                    silent[v] += 1
        dt = time.perf_counter() - t0
        res[mode] = {
            v: dict(
                n=len(R),
                intermediates={k: s.as_dict() for k, s in stats[v].items()},
                frame_rows_differing=frame_diff[v],
                field_orientations_differing=field_diff[v],
                n_orientations_differing=len(differing[v]),
                n_traced=sum(x["traced"] for x in differing[v]),
                untraced=[x["index"] for x in differing[v] if not x["traced"]],
                max_margin_ulp=max(
                    [r["margin_ulp"] for x in differing[v] for r in x["rows"]], default=None
                ),
                ref_side={
                    s: sum(r.get("ref_side") == s for x in differing[v] for r in x["rows"])
                    for s in ("torch", "fused")
                },
                error_vs_float64=dict(
                    torch={k: e.as_dict() for k, e in err["torch"].items()},
                    fused={k: e.as_dict() for k, e in err[v].items()},
                ),
                differing=differing[v],
                silent_input_differences=silent[v],
            )
            for v in variants
        }
        line = "  ".join(
            f"{v}: {len(differing[v])} differ ({sum(x['traced'] for x in differing[v])} traced)"
            for v in variants
        )
        print(f"  {mode:15s} n={len(R)} {dt:5.1f}s  {line}", flush=True)
    return res


def main() -> int:
    torch.set_num_threads(1)
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--datasets", nargs="+", default=["three0"])
    ap.add_argument("--modes", nargs="+", default=list(MODES), choices=list(MODES))
    ap.add_argument("--variants", nargs="+", default=["portable", "accelerate"])
    ap.add_argument("--n-near", type=int, default=500)
    ap.add_argument("--n-fz", type=int, default=500)
    ap.add_argument("--trace-ulp", type=int, default=8)
    a = ap.parse_args()
    out_path = Path(a.out).resolve()
    result: Dict[str, Any] = dict(
        meta=dict(
            head=S1.git_head(),
            dirty=S1.git_dirty(),
            dirty_files=dirty_files(),
            torch=torch.__version__,
            numpy=np.__version__,
            threads=torch.get_num_threads(),
            seed=S1.SEED,
            n_near=a.n_near,
            n_fz=a.n_fz,
            trace_ulp=a.trace_ulp,
            variants=a.variants,
            modes={k: MODES[k] for k in a.modes},
        ),
        datasets={},
    )
    for name in a.datasets:
        t0 = time.perf_counter()
        with contextlib.redirect_stdout(io.StringIO()):
            ds = S1.load_dataset(name)
        R = S1.orientation_set(ds, a.n_near, a.n_fz)
        print(f"[{name}] loaded in {time.perf_counter() - t0:.1f} s, {len(R)} orientations")
        result["datasets"][name] = run_dataset(ds, R, a.modes, a.variants, a.trace_ulp)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=1) + "\n")
    bad = sum(
        len(r[v]["untraced"])
        for d in result["datasets"].values()
        for r in d.values()
        for v in a.variants
    )
    print(f"wrote {out_path.name}; untraced differing orientations: {bad}")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
