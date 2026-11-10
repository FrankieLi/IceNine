"""S1 parity harness: every OverlapInfo field of VoxelCostFunction.evaluate, for a fixed set of
orientations, voxels and cost modes, dumped to an .npz and compared exactly (``==``).

S1 (docs/todo_s1_evaluate_loops.md) changes how ``evaluate`` finds the frame of each omega, how it
hands the binary frames to the C stage and how it builds the Bragg constants, and must not change a
single result. Step 0 dumps the base commit; Step 5 dumps the new code and compares.

  uv run python scripts/speed_parity/s1_parity.py dump --out scripts/speed_parity/cache/s1/base.npz
  uv run python scripts/speed_parity/s1_parity.py dump --out scripts/speed_parity/cache/s1/new.npz
  uv run python scripts/speed_parity/s1_parity.py compare base.npz new.npz
  uv run python scripts/speed_parity/s1_parity.py batched \
      --dump scripts/speed_parity/cache/s1/base.npz
  uv run python scripts/speed_parity/s1_parity.py compare-e0 \
      scripts/findoptimal_robustness/cache/e0_s1base scripts/findoptimal_robustness/cache/e0_s1
  uv run python scripts/speed_parity/s1_parity.py dump --fixture \
      --out tests/data/s1_evaluate_ref.npz     # the small fixture of tests/test_s1_evaluate.py

Datasets (``--datasets``): ``three0`` = Example2.ThreeVoxels voxel 0 (forward-simulated images);
``full<N>`` = voxel N of the Phase D clean full-sample images (Q-max 8, 180 frames, 2 detectors;
interior voxels 15901, 17034, 18558). Orientations per dataset: ``--n-near`` within 3 degrees of the
truth (random axis, uniform angle) and ``--n-fz`` drawn without replacement from the fundamental
zone file of the config, seed fixed. The modes marked slow (soft mode, the Python fallback path with
the C extension patched off) use ``--n-slow-near`` + ``--n-slow-fz`` of them (the first of each
group).

Each mode is a separate ``VoxelCostFunction`` (see MODES). ``dump`` also counts, for the local and
the screen cost, the omega rows that reach the edges of the scan (the bin ``[low - width, low)``,
bins above the range, bins of an ``index_list`` gap) so the coverage is shown, not assumed.
Single-threaded (``torch.set_num_threads(1)``) on both sides of a comparison: thread count can
change the rounding of the 3x3 @ 3xK product.
"""

import argparse
import contextlib
import io
import json
import os
import platform
import subprocess
import sys
import time
import zlib
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
ROOT = ICE.parent
sys.path.insert(0, str(ICE))

THREE = ROOT / "Examples" / "Example2.ThreeVoxels"
MANY = ROOT / "Examples" / "Example2.ManyGrains"
STACK = MANY / "ScatteringData_PhaseD" / "stacks" / "clean.npy"
PHASE_D_CONFIG = ICE / "scripts" / "phase_d" / "configs" / "ReconstructPhaseD_mc_clean.config"
FULL_VOXELS = (15901, 17034, 18558)
SEED = 20261014

FIELDS = (
    "pixel_overlap",
    "pixel_on_detector",
    "peak_overlap",
    "peak_on_detector",
    "quality",
    "n_quality_points",
    "detectors_overlap",
)

# name -> (VoxelCostFunction keywords, slow, patch). ``ray`` adds RayRule("R", 2 px); ``patch`` is
# applied around the evaluations (module flags restored afterwards).
MODES: Dict[str, Tuple[Dict[str, Any], bool, Optional[str]]] = {
    "local": (dict(pixel_radius=0, max_q=8.0), False, None),
    "cpp_rule": (dict(pixel_radius=0, max_q=8.0, cpp_overlap_rule=True), False, None),
    "edge": (dict(pixel_radius=0, max_q=8.0, cpp_edge_rules=True), False, None),
    "edge_cpp_rule": (
        dict(pixel_radius=0, max_q=8.0, cpp_edge_rules=True, cpp_overlap_rule=True),
        False,
        None,
    ),
    "screen_center": (
        dict(pixel_radius=3, max_q=5.0, cpp_overlap_rule=True, center_screen=True),
        False,
        None,
    ),
    "screen_vertex0": (
        dict(pixel_radius=3, max_q=5.0, cpp_overlap_rule=True, center_screen=False),
        False,
        None,
    ),
    "screen_edge": (
        dict(pixel_radius=3, max_q=5.0, cpp_overlap_rule=True, cpp_edge_rules=True),
        False,
        None,
    ),
    "ray_c_screen": (dict(pixel_radius=3, max_q=5.0, ray=True), False, None),
    "ray_py_screen": (dict(pixel_radius=3, max_q=5.0, ray=True), False, "ray_python"),
    "ray_c_local": (dict(pixel_radius=0, max_q=8.0, ray=True), False, None),
    "min_sin_eta": (dict(pixel_radius=0, max_q=8.0, min_sin_eta=0.2), False, None),
    "soft": (dict(pixel_radius=0, max_q=8.0, mode="soft"), True, None),
    "fallback_local": (dict(pixel_radius=0, max_q=8.0), True, "no_c_ext"),
    "fallback_cpp_edge": (
        dict(pixel_radius=0, max_q=8.0, cpp_overlap_rule=True, cpp_edge_rules=True),
        True,
        "no_c_ext",
    ),
    "fallback_screen": (dict(pixel_radius=3, max_q=5.0, cpp_overlap_rule=True), True, "no_c_ext"),
}


# ----------------------------------------------------------------------------------------------
# datasets
# ----------------------------------------------------------------------------------------------
class Dataset:
    """Everything a VoxelCostFunction needs plus the voxel, the truth and the FZ orientations."""

    def __init__(
        self,
        name: str,
        setup: Any,
        verts: torch.Tensor,
        phase: int,
        truth: np.ndarray,
        fz: np.ndarray,
    ) -> None:
        self.name, self.setup, self.verts, self.phase = name, setup, verts, phase
        self.truth, self.fz = truth, fz

    def cost_fn(self, kw: Dict[str, Any]) -> Any:
        from icenine.cost_functions import VoxelCostFunction
        from icenine.ray_overlap import RayRule

        kw = dict(kw)
        s = self.setup
        if kw.pop("ray", False):
            reg = float(s.exp_setup.config_file.detector_dist_deviation)
            kw["ray_rule"] = RayRule(arm="R", radius_px=2, registration_mm=reg)
        return VoxelCostFunction(
            simulator=s.simulator,
            detector_list=s.detector_list,
            range_map=s.range_map,
            exp_data=s.exp_data,
            sample=s.sample,
            structure_list=s.structure_list,
            eta_limit=s.exp_setup.get_eta_limit(),
            **{"mode": "hard", **kw},
        )


def load_dataset(name: str, three_mode: str = "sparse") -> Dataset:
    from icenine.config_file import ConfigFile
    from icenine.experimental_data import ExperimentalData
    from icenine.mic_file import MicFile
    from icenine.reconstructor import _get_voxel_vertices, setup_reconstruction
    from icenine.sampling import load_fundamental_zone_file

    if name == "three0":
        os.chdir(THREE)
        cfg = ConfigFile.from_file(str(THREE / "ConfigFiles" / "Example2.Simulation.config"))
        exp = ExperimentalData.from_image_directory(
            directory=str(THREE / "ScatteringData_Python"),
            basename="3Grains.sim",
            ext="d",
            serial_length=5,
            n_omega=180,
            n_detectors=2,
            num_rows=2048,
            num_cols=2048,
            mode=three_mode,
        )
        exp.prepare_for_reconstruction()
        fz = np.asarray(load_fundamental_zone_file(cfg.fundamental_zone_filename))
        setup = setup_reconstruction(cfg, exp_data=exp, fz_orientations=fz)
        vox = MicFile.read(str(THREE / "SimInput" / "three_voxels.mic")).voxels[0]
        return Dataset(
            name,
            setup,
            _get_voxel_vertices(vox),
            int(vox.phase),
            np.asarray(vox.orientation, dtype=np.float64),
            np.asarray(fz),
        )
    if name.startswith("full"):
        v = int(name[4:])
        os.chdir(MANY)
        cfg = ConfigFile.from_file(str(PHASE_D_CONFIG))
        exp = ExperimentalData.from_binary_memmap(STACK, 180, 2)
        setup = setup_reconstruction(cfg, exp_data=exp)
        vox = setup.sample.get_mic().voxels[v]
        truth = MicFile.read("SimInput/rand_500grains_1mm_neworient_s0.mic").voxels[v]
        return Dataset(
            name,
            setup,
            _get_voxel_vertices(vox),
            int(vox.phase),
            np.asarray(truth.orientation, dtype=np.float64),
            np.asarray(setup.fz_orientations),
        )
    raise ValueError(f"unknown dataset {name!r}")


def _rot(axis: np.ndarray, ang: float) -> np.ndarray:
    axis = axis / np.linalg.norm(axis)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * K @ K


def orientation_set(ds: Dataset, n_near: int, n_fz: int, seed: int = SEED) -> np.ndarray:
    """(n_near + n_fz, 3, 3) float64: near the truth (<= 3 degrees), then FZ points."""
    rng = np.random.default_rng([seed, zlib.crc32(ds.name.encode())])
    near = [
        _rot(rng.normal(size=3), np.deg2rad(rng.uniform(0.0, 3.0))) @ ds.truth
        for _ in range(n_near)
    ]
    pick = rng.choice(len(ds.fz), size=n_fz, replace=False)
    fz = [np.asarray(ds.fz[i], dtype=np.float64) for i in pick]
    return np.stack(near + fz)


# ----------------------------------------------------------------------------------------------
# evaluation
# ----------------------------------------------------------------------------------------------
@contextlib.contextmanager
def patched(kind: Optional[str]) -> Any:
    from icenine import cost_functions as CF

    saved = (CF._HAS_C_RASTERIZE, CF.USE_C_RAY_FLAGS)
    try:
        if kind == "no_c_ext":
            CF._HAS_C_RASTERIZE = False
        elif kind == "ray_python":
            CF.USE_C_RAY_FLAGS = False
        yield
    finally:
        CF._HAS_C_RASTERIZE, CF.USE_C_RAY_FLAGS = saved


def run_mode(ds: Dataset, mode: str, R: np.ndarray) -> Dict[str, np.ndarray]:
    kw, _slow, patch = MODES[mode]
    cost_fn = ds.cost_fn(kw)
    n = len(R)
    out = {f: np.zeros(n, dtype=np.float64 if f == "quality" else np.int64) for f in FIELDS}
    rays = np.full((n, 8), -1, dtype=np.int64)
    with patched(patch), contextlib.redirect_stdout(io.StringIO()):
        for i in range(n):
            o = cost_fn.evaluate(R[i], ds.verts, ds.phase)
            for f in FIELDS:
                out[f][i] = getattr(o, f)
            if o.ray_cases is not None:
                rays[i] = np.asarray(o.ray_cases, dtype=np.int64)
    out["ray_cases"] = rays
    return out


def coverage(ds: Dataset, R: np.ndarray) -> Dict[str, Dict[str, int]]:
    """Rows reaching the scan edges, for the local (Q-max 8) and the screen (Q-max 5) cost."""
    from icenine import cost_functions as CF

    res: Dict[str, Dict[str, int]] = {}
    for label, mode in (("local", "local"), ("screen", "screen_center")):
        cost_fn = ds.cost_fn(MODES[mode][0])
        rng = cost_fn.range_map
        n_bins = len(rng.index_list)
        none_bin = np.array([i is None for i in rng.index_list] + [True])
        seen: List[np.ndarray] = []
        orig = CF.calculate_diffraction_overlap_batched

        def spy(**kw: Any) -> Any:
            seen.append(kw["peak_omegas"].numpy().copy())
            return orig(**kw)

        CF.calculate_diffraction_overlap_batched = spy
        try:
            for i in range(len(R)):
                cost_fn.evaluate(R[i], ds.verts, ds.phase)
        finally:
            CF.calculate_diffraction_overlap_batched = orig
        om = np.concatenate(seen)
        f = (om - rng.low) / rng.width
        fl = np.floor(f).astype(int)
        res[label] = dict(
            rows=int(len(om)),
            below_edge_bin=int(((f >= -1) & (f < 0)).sum()),
            above_range=int((fl >= n_bins).sum()),
            far_below=int((f < -1).sum()),
            gap_bins=int(((fl >= 0) & (fl < n_bins) & none_bin[np.clip(fl, 0, n_bins)]).sum()),
            frames_used=int(
                len(
                    np.unique(
                        [
                            rng.index_list[b]
                            for b in fl
                            if 0 <= b < n_bins and rng.index_list[b] is not None
                        ]
                    )
                )
            ),
        )
    return res


def git_head() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, text=True
        ).strip()
    except Exception:  # noqa: BLE001
        return "unknown"


def git_dirty() -> Optional[bool]:
    """True when a tracked file differs from HEAD (untracked files do not count); None if git
    cannot say. A dump of a dirty tree is not a dump of the commit named in ``head``."""
    try:
        out = subprocess.check_output(
            ["git", "status", "--porcelain", "--untracked-files=no"], cwd=ROOT, text=True
        )
        return bool(out.strip())
    except Exception:  # noqa: BLE001
        return None


def cmd_dump(a: argparse.Namespace) -> int:
    from icenine import cost_functions as CF

    out_path = Path(a.out).resolve()  # the dataset loaders change the working directory
    out: Dict[str, Any] = {}
    meta: Dict[str, Any] = dict(
        head=git_head(),
        dirty=git_dirty(),
        torch=torch.__version__,
        platform_system=platform.system(),
        platform_machine=platform.machine(),
        threads=int(torch.get_num_threads()),
        c_ext=bool(CF._HAS_C_RASTERIZE),
        seed=SEED,
        n_near=a.n_near,
        n_fz=a.n_fz,
        n_slow=[a.n_slow_near, a.n_slow_fz],
        datasets=list(a.datasets),
        modes=[],
        coverage={},
    )
    modes = a.modes or list(MODES)
    for name in a.datasets:
        t0 = time.perf_counter()
        ds = load_dataset(name, a.three_mode)
        R = orientation_set(ds, a.n_near, a.n_fz)
        out[f"{name}/R"] = R
        print(f"[{name}] loaded in {time.perf_counter() - t0:.1f} s, {len(R)} orientations")
        if not a.fixture:
            meta["coverage"][name] = coverage(ds, R)
            for label, c in meta["coverage"][name].items():
                print(f"  coverage {label}: {c}")
        for mode in modes:
            slow = MODES[mode][1]
            idx = (
                np.r_[np.arange(a.n_slow_near), a.n_near + np.arange(a.n_slow_fz)]
                if slow
                else np.arange(len(R))
            )
            t0 = time.perf_counter()
            res = run_mode(ds, mode, R[idx])
            for f, v in res.items():
                out[f"{name}/{mode}/{f}"] = v
            out[f"{name}/{mode}/index"] = idx
            print(
                f"  {mode:18s} n={len(idx):5d}  {time.perf_counter() - t0:6.1f} s"
                f"  mean quality {res['quality'].mean():.4f}"
            )
            if mode not in meta["modes"]:
                meta["modes"].append(mode)
    out["meta"] = np.array(json.dumps(meta))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out_path, **out)
    print(f"wrote {out_path.name} ({out_path.stat().st_size / 1024:.0f} KB)")
    return 0


def _row_differs(x: np.ndarray, y: np.ndarray, n: int) -> np.ndarray:
    """Per-row ``x != y`` (exact). A shape or dtype mismatch is a difference in every row."""
    if x.shape != y.shape or x.dtype != y.dtype:
        return np.ones(n, dtype=bool)
    d = x != y
    return d.any(axis=tuple(range(1, d.ndim))) if d.ndim > 1 else d


def compare_dumps(A: Any, B: Any, say: Any = print) -> int:
    """Compare two dumps (open npz files); returns the number of differences found (0 = every
    field of every mode identical). Compared exactly: the key sets, the orientation sets ``R``,
    the orientation indices ``index``, and every field array (shape, dtype and values with
    ``!=``, ``quality`` included)."""
    ka, kb = set(A.files) - {"meta"}, set(B.files) - {"meta"}
    bad = 0
    if ka != kb:
        say(f"KEY MISMATCH: only in A {sorted(ka - kb)[:5]}, only in B {sorted(kb - ka)[:5]}")
        bad += 1
    groups: Dict[Tuple[str, str], List[str]] = {}
    for k in sorted(ka & kb):
        parts = k.split("/")
        if len(parts) == 3 and parts[2] != "index":
            groups.setdefault((parts[0], parts[1]), []).append(parts[2])
        elif len(parts) == 2 and parts[1] == "R" and not np.array_equal(A[k], B[k]):
            say(f"ORIENTATION SET DIFFERS for {parts[0]}")
            bad += 1
        elif len(parts) == 3 and not np.array_equal(A[k], B[k]):
            say(f"ORIENTATION INDEX DIFFERS for {parts[0]}/{parts[1]}")
            bad += 1
    say(f"{'dataset':10s} {'mode':18s} {'n':>5s}  differing orientations per field")
    for (ds, mode), fields in sorted(groups.items()):
        n = len(A[f"{ds}/{mode}/{fields[0]}"])
        diff: Dict[str, int] = {}
        any_row = np.zeros(n, dtype=bool)
        for f in fields:
            x, y = A[f"{ds}/{mode}/{f}"], B[f"{ds}/{mode}/{f}"]
            d = _row_differs(x, y, n)
            if len(d) != n:  # a field array whose length differs from that of the first one
                d = np.ones(n, dtype=bool)
            if d.any():
                diff[f] = int(d.sum())
            any_row |= d
        bad += len(diff)
        msg = "identical" if not diff else f"DIFFER {diff} ({int(any_row.sum())} orientations)"
        say(f"{ds:10s} {mode:18s} {n:5d}  {msg}")
    say("RESULT: " + ("bit-identical (every field of every mode)" if bad == 0 else "DIFFERENCES"))
    return bad


def cmd_compare(a: argparse.Namespace) -> int:
    A, B = np.load(a.a, allow_pickle=False), np.load(a.b, allow_pickle=False)
    return 0 if compare_dumps(A, B) == 0 else 1


BATCHED_MODES = ("local", "cpp_rule", "edge", "edge_cpp_rule", "min_sin_eta")


def cmd_batched(a: argparse.Namespace) -> int:
    """S1 Step 4: ``BatchedCostEvaluator`` (compact and fixed variant, on the CPU) against the
    ``evaluate`` results of a dump, every field exactly, for the modes the batched cost supports
    (it needs pixel_radius 0, hard mode and no ray rule). ``--dump`` is normally the Step 0 dump
    of the base commit, so the reference is the original ``evaluate``."""
    from icenine.batched_cost import BatchedCostEvaluator, unsupported_reason

    ref = np.load(a.dump, allow_pickle=False)
    bad = 0
    print(
        f"{'dataset':10s} {'mode':14s} {'variant':8s} {'n':>5s}  differing orientations per field"
    )
    for name in a.datasets:
        ds = load_dataset(name, a.three_mode)
        R = ref[f"{name}/R"]
        for mode in a.modes:
            cost_fn = ds.cost_fn(MODES[mode][0])
            assert unsupported_reason(cost_fn) is None, (mode, unsupported_reason(cost_fn))
            for variant in ("compact", "fixed"):
                cost_fn = ds.cost_fn(MODES[mode][0])
                ev = BatchedCostEvaluator(cost_fn, "cpu")
                got: List[Any] = []
                with contextlib.redirect_stdout(io.StringIO()):
                    for s0 in range(0, len(R), a.batch):
                        got += ev.evaluate_many(R[s0 : s0 + a.batch], ds.verts, ds.phase, variant)
                diff: Dict[str, int] = {}
                for f in FIELDS:
                    x = np.array([getattr(o, f) for o in got])
                    d = int((x != ref[f"{name}/{mode}/{f}"]).sum())
                    if d:
                        diff[f] = d
                bad += len(diff)
                msg = "identical" if not diff else f"DIFFER {diff}"
                print(f"{name:10s} {mode:14s} {variant:8s} {len(R):5d}  {msg}")
    print("RESULT:", "batched cost equals evaluate (every field)" if bad == 0 else "DIFFERENCES")
    return 0 if bad == 0 else 1


def cmd_compare_e0(a: argparse.Namespace) -> int:
    """Per-voxel records of two E0 runs (cache/e0<tag> directories): every stored array must be
    equal except the wall-clock ``*_runtime``; wrong counts (final error > 1 degree, cubic
    reduced) and evaluation totals are reported side by side."""
    sys.path.insert(0, str(ICE / "scripts" / "findoptimal_robustness"))
    from csl import reduced_misorientation_deg

    da, db = Path(a.dir_a), Path(a.dir_b)
    names = sorted(p.name for p in da.glob("v*.npz"))
    only = set(names) ^ {p.name for p in db.glob("v*.npz")}
    if only:
        print(f"FILES ONLY IN ONE DIRECTORY: {sorted(only)[:5]} ({len(only)})")
    names = [n for n in names if (db / n).exists()]
    n_same = 0
    diff_keys: Dict[str, int] = {}
    tot = {"A": dict(wrong=0, eg=0, el=0, n=0), "B": dict(wrong=0, eg=0, el=0, n=0)}
    for n in names:
        x, y = np.load(da / n, allow_pickle=False), np.load(db / n, allow_pickle=False)
        if "unbuildable" in x.files or "unbuildable" in y.files:
            continue
        bad = [
            k
            for k in x.files
            if "runtime" not in k and (k not in y.files or not np.array_equal(x[k], y[k]))
        ]
        bad += [k for k in y.files if k not in x.files]
        for k in bad:
            diff_keys[k] = diff_keys.get(k, 0) + 1
        n_same += not bad
        for tag, d in (("A", x), ("B", y)):
            err = float(reduced_misorientation_deg(d["s0_R_final"], d["R_true"]))
            tot[tag]["wrong"] += err > 1.0
            tot[tag]["eg"] += int(d["s0_evals_global"])
            tot[tag]["el"] += int(d["s0_evals_local"])
            tot[tag]["n"] += 1
    print(f"voxels compared {len(names)}, records identical (all arrays but runtime): {n_same}")
    if diff_keys:
        print(f"DIFFERING KEYS (records): {dict(sorted(diff_keys.items())[:12])}")
    for tag, d in tot.items():
        print(
            f"  run {tag}: n={d['n']}  wrong (>1 deg) {d['wrong']}  evaluations global"
            f" {d['eg']}  local {d['el']}"
        )
    ok = (not only) and n_same == tot["A"]["n"] == tot["B"]["n"] and not diff_keys
    print("RESULT:", "per-voxel records identical" if ok else "DIFFERENCES")
    return 0 if ok else 1


def main() -> int:
    torch.set_num_threads(1)  # as the fixture and the dumps; tests set it around their own runs
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("dump")
    d.add_argument("--out", required=True)
    d.add_argument("--datasets", nargs="+", default=["three0"] + [f"full{v}" for v in FULL_VOXELS])
    d.add_argument("--modes", nargs="*", default=None, choices=list(MODES))
    d.add_argument("--n-near", type=int, default=500)
    d.add_argument("--n-fz", type=int, default=500)
    d.add_argument("--n-slow-near", type=int, default=25)
    d.add_argument("--n-slow-fz", type=int, default=25)
    d.add_argument("--three-mode", default="sparse", choices=["sparse", "dense"])
    d.add_argument(
        "--fixture",
        action="store_true",
        help="small tests/data fixture: three0 only, 32 + 32 orientations, 8 + 8 for slow modes",
    )
    c = sub.add_parser("compare")
    c.add_argument("a")
    c.add_argument("b")
    b = sub.add_parser(
        "batched", help="BatchedCostEvaluator against the evaluate results of a dump"
    )
    b.add_argument("--dump", required=True)
    b.add_argument("--datasets", nargs="+", default=["three0"] + [f"full{v}" for v in FULL_VOXELS])
    b.add_argument("--modes", nargs="+", default=list(BATCHED_MODES), choices=list(MODES))
    b.add_argument("--batch", type=int, default=64)
    b.add_argument("--three-mode", default="sparse", choices=["sparse", "dense"])
    e = sub.add_parser("compare-e0", help="two E0 cache directories (cache/e0<tag>)")
    e.add_argument("dir_a")
    e.add_argument("dir_b")
    a = ap.parse_args()
    if a.cmd == "compare-e0":
        return cmd_compare_e0(a)
    if a.cmd == "batched":
        return cmd_batched(a)
    if a.cmd == "dump" and a.fixture:
        a.datasets = ["three0"]
        a.n_near = a.n_fz = 32
        a.n_slow_near = a.n_slow_fz = 8
    return cmd_dump(a) if a.cmd == "dump" else cmd_compare(a)


if __name__ == "__main__":
    sys.exit(main())
