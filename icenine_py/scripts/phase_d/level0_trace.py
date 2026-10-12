"""Level-0 candidates, C++ vs Python, search point by search point (Phase D task level0-candidates).

Level 0 of a no-start search has no random draws (fixed FZ grid x fixed local grid), so the two
implementations can be compared per FZ point: which local-grid points pass the screen
(LocalAngularInterpProcess + AcceptAnyPeakOverlap, pixel radius 3), the re-cost 1 - confidence of
the passing ones, and which the spacing filter keeps. Levels 1+ are compared by counts.

  uv run python scripts/phase_d/level0_trace.py py  --voxel 15901 [--full]   # Python trace
  ... py --cpp-rule [--centroid]    # same, with the C++ CountQualifiedPeaks rule (read-only probe)
  uv run python scripts/phase_d/level0_trace.py cpp --voxel 15901 --log LOG  # parse a C++ log
  uv run python scripts/phase_d/level0_trace.py compare --voxel 15901        # S4, S1, S2, S3

C++ side: scratch build with benchmarks/phase_d_seed_diag/level0/cpp_debug_patch_level0.py.txt
(not committed), run with L0TRACE=1 from Examples/Example2.ManyGrains. Output goes to
benchmarks/phase_d_seed_diag/level0/ (cpp_v<V>.json, py_v<V>.json, compare_v<V>.json).
"""

import argparse
import json
import math
import os
import re
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
EX = HERE.parents[2] / "Examples" / "Example2.ManyGrains"
OUT = ICE / "benchmarks" / "phase_d_seed_diag" / "level0"
sys.path.insert(0, str(ICE / "scripts" / "common"))


CPP_COUNT_RULE = """static void count_qualified_peaks_c(
    const int *detector_lit, const int *spot_overlap, int n_detectors,
    int *out_peak_on_det, int *out_peak_overlap, int *out_n_det_overlap
) {
    /* Literal port of C++ CountQualifiedPeaks (OverlapInfo.tmpl.cpp:78-139) */
    *out_peak_on_det = 0;
    *out_peak_overlap = 0;
    *out_n_det_overlap = 0;
    if (n_detectors == 0) return;
    int bValidSimPeak = detector_lit[0];
    int bOnPreviousDetector = detector_lit[0];
    int nDetectorsOverlap = 0;
    for (int i = 1; i < n_detectors && bValidSimPeak; i++) {
        if (detector_lit[i]) {
            bValidSimPeak = 1;
        } else {
            if (bOnPreviousDetector) { bValidSimPeak = 1; bOnPreviousDetector = 0; }
            else bValidSimPeak = 0;
        }
    }
    int bPreviousOverlap = spot_overlap[0];
    int bValidOverlap = spot_overlap[0];
    if (spot_overlap[0]) nDetectorsOverlap++;
    for (int i = 1; i < n_detectors && bValidOverlap; i++) {
        if (spot_overlap[i]) {
            nDetectorsOverlap++;
            bValidOverlap = 1;
        } else {
            if (bPreviousOverlap && !detector_lit[i]) { bValidOverlap = 1; bPreviousOverlap = 0; }
            else bValidOverlap = 0;
        }
    }
    if (bValidSimPeak) *out_peak_on_det = 1;
    if (bValidOverlap) *out_peak_overlap = 1;
    if (bValidSimPeak && bValidOverlap) *out_n_det_overlap = nDetectorsOverlap;
}

"""


def build_cpp_rule_ext() -> Any:
    """Compile a copy of icenine/_rasterize.c whose count_qualified_peaks_c is the literal C++
    CountQualifiedPeaks, as module _rasterize_cpp in a scratch dir (probe only; the library is not
    touched). Returns the imported module."""
    import importlib.util
    import subprocess
    import sysconfig
    import tempfile

    import numpy

    src = (ICE / "icenine" / "_rasterize.c").read_text()
    i = src.index("static void count_qualified_peaks_c(")
    j = src.index("/* ---- Main triangle_overlap function")
    src = src[:i] + CPP_COUNT_RULE + src[j:]
    src = src.replace('"_rasterize",\n    "Fast', '"_rasterize_cpp",\n    "Fast')
    src = src.replace("PyInit__rasterize(void)", "PyInit__rasterize_cpp(void)")
    d = Path(os.environ.get("LEVEL0_SCRATCH", tempfile.gettempdir())) / "rast_cpp_rule"
    d.mkdir(parents=True, exist_ok=True)
    (d / "_rasterize_cpp.c").write_text(src)
    suf = sysconfig.get_config_var("EXT_SUFFIX")
    out = d / f"_rasterize_cpp{suf}"
    cmd = [
        "cc", "-O2", "-shared", "-undefined", "dynamic_lookup", "-fPIC",
        f"-I{sysconfig.get_paths()['include']}", f"-I{numpy.get_include()}",
        str(d / "_rasterize_cpp.c"), "-o", str(out),
    ]  # fmt: skip
    subprocess.run(cmd, check=True, capture_output=True)
    spec = importlib.util.spec_from_file_location("_rasterize_cpp", out)
    mod = importlib.util.module_from_spec(spec)  # type: ignore[arg-type]
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    return mod


class StopTrace(Exception):
    """Raised by the recorder to end a run after level 0 (--level0-only)."""


# --------------------------------------------------------------------------------------------
# Python trace
# --------------------------------------------------------------------------------------------
def run_py(a: argparse.Namespace) -> None:
    os.chdir(EX)
    import torch

    torch.set_num_threads(1)
    from icenine import reconstructor as RC
    from icenine.config_file import ConfigFile
    from icenine.experimental_data import ExperimentalData
    import preflight

    cfg = ConfigFile.from_file(str(HERE / "configs" / "ReconstructPhaseD_mc_clean.config"))
    stack = EX / "ScatteringData_PhaseD" / "stacks" / "clean.npy"
    exp = ExperimentalData.from_binary_memmap(stack, 180, 2)
    setup = RC.setup_reconstruction(cfg, exp_data=exp)
    rec = RC.AdaptiveVoxelReconstructor(setup)
    v = setup.sample.get_mic().voxels[a.voxel]
    verts = RC._get_voxel_vertices(v)

    levels: List[Dict[str, Any]] = []
    fz0: List[Dict[str, Any]] = []

    def hook(ev: str, d: Dict[str, Any]) -> None:
        if ev == "level":
            levels.append(dict(d))
        elif ev == "fz" and len(levels) == 1:
            fz0.append(d)
        elif ev == "level_end":
            levels[-1].update(d)
            if a.level0_only and len(levels) == 1:
                raise StopTrace

    rec.discrete_recorder = hook
    if a.cpp_rule:
        from icenine import cost_functions as CF

        CF._c_stage_d_overlap = build_cpp_rule_ext().stage_d_overlap  # probe only
    if a.centroid:
        # C++ GlobalSearchCostFunctions uses VoxelToCenter: the pixel-radius screen projects the
        # voxel centroid, the library projects vertex 0. Probe: swap in the centroid for r > 0 only.
        from icenine import cost_functions as CF

        orig_eval = CF.VoxelCostFunction.evaluate

        def eval_centroid(self, orientation, voxel_vertices, phase_index=0, **kw):  # type: ignore
            if self.pixel_radius > 0:
                voxel_vertices = voxel_vertices.mean(0, keepdim=True).repeat(3, 1)
            return orig_eval(self, orientation, voxel_vertices, phase_index, **kw)

        CF.VoxelCostFunction.evaluate = eval_centroid  # type: ignore[method-assign]
    pf = preflight.preflight()
    rng = np.random.default_rng(a.seed)
    t0 = time.perf_counter()
    try:
        rec.reconstruct_voxel(voxel_vertices=verts, phase_index=v.phase, rng=rng)
    except StopTrace:
        pass
    wall = time.perf_counter() - t0
    # symmetry the filter used (proper rotations of the configured symmetry)
    out = dict(
        voxel=a.voxel,
        side="python",
        seed=a.seed,
        cpp_rule=a.cpp_rule,
        centroid=a.centroid,
        level0_only=a.level0_only,
        wall_s=wall,
        label="single-process, " + ("quiet" if pf.get("quiet", True) else "contended"),
        levels=levels,
        fz=fz0,
    )
    tag = ("_cpprule" if a.cpp_rule else "") + ("_centroid" if a.centroid else "")
    fn = OUT / f"py{tag}_v{a.voxel}{'_L0' if a.level0_only else ''}.json"
    fn.write_text(json.dumps(out, default=float) + "\n")
    print(f"python v{a.voxel}: {wall:.1f}s levels={len(levels)} -> {fn.name}")


# --------------------------------------------------------------------------------------------
# C++ log parser
# --------------------------------------------------------------------------------------------
def parse_cpp(a: argparse.Namespace) -> None:
    levels: List[Dict[str, Any]] = []
    fz0: Dict[int, Dict[str, Any]] = {}
    ov0: Dict[int, List[int]] = {}
    nobs: Dict[int, int] = {}
    for line in Path(a.log).read_text(errors="replace").splitlines():
        if line.startswith("L0H "):
            t = line.split()
            levels.append(
                dict(
                    n_fz=int(t[4]),
                    n_local=int(t[6]),
                    max_q=float(t[8]),
                    diameter_rad=float(t[10]),
                    n_sym=int(t[14]),
                    n_recip=int(t[16]),
                )
            )
        elif line.startswith("L0S "):
            t = line.split()
            levels[-1].update(n_cliques=int(t[4]), n_candidates=int(t[6]))
        elif line.startswith("L0O "):
            t = line.split()
            ov0[int(t[1])] = [int(x) for x in t[3:]]
            nobs[int(t[1])] = int(t[2])
        elif line.startswith("L0P "):
            head, kept = line.split(" | kept")
            m = re.match(r"L0P (\d+) (\d+) (\d+) j (.*) cost (.*)", head)
            fz = int(m.group(1))
            js = [int(x) for x in m.group(4).split()]
            cs = [float(x) for x in m.group(5).split()]
            fz0[fz] = dict(
                fz=fz, n_obs=int(m.group(2)), j=js, cost=cs, kept=[int(x) for x in kept.split()]
            )
    rows = []
    for fz in sorted(ov0):
        r = fz0.get(fz, dict(fz=fz, j=[], cost=[], kept=[]))
        r = dict(r, peak_overlap=ov0[fz], n_obs=nobs[fz])
        rows.append(r)
    out = dict(voxel=a.voxel, side="cpp", levels=levels, fz=rows)
    fn = OUT / f"cpp_v{a.voxel}.json"
    fn.write_text(json.dumps(out) + "\n")
    print(f"cpp v{a.voxel}: levels={len(levels)} fz={len(rows)} -> {fn.name}")


# --------------------------------------------------------------------------------------------
# Comparison, in order S4 (settings), S1 (screen), S2 (re-cost), S3 (spacing filter)
# --------------------------------------------------------------------------------------------
def compare(a: argparse.Namespace) -> None:
    c = json.loads((OUT / f"cpp_v{a.voxel}.json").read_text())
    pyfn = OUT / f"py{a.py_tag}_v{a.voxel}.json"
    if not pyfn.exists():
        pyfn = OUT / f"py{a.py_tag}_v{a.voxel}_L0.json"
    p = json.loads(pyfn.read_text())
    res: Dict[str, Any] = dict(voxel=a.voxel, py_file=pyfn.name)

    # S4 settings per level
    s4 = []
    for lv, (lc, lp) in enumerate(zip(c["levels"], p["levels"])):
        s4.append(
            dict(
                level=lv,
                n_fz=(lc["n_fz"], lp["n_fz"]),
                n_local=(lc["n_local"], lp["n_local"]),
                n_q_max=(lc["max_q"], None),
                n_recip=(lc["n_recip"], lp["n_recip"]),
                diameter_deg=(math.degrees(lc["diameter_rad"]), math.degrees(lp["diameter_rad"])),
                n_sym=(lc["n_sym"], lp["n_sym"]),
                n_cliques=(lc.get("n_cliques"), lp.get("n_cliques")),
                n_candidates=(lc.get("n_candidates"), lp.get("n_candidates")),
            )
        )
    res["S4_levels"] = s4

    cf = {r["fz"]: r for r in c["fz"]}
    pf = {r["fz"]: r for r in p["fz"]}
    common = sorted(set(cf) & set(pf))
    res["n_fz_common"] = len(common)

    # S1 screen: which j pass
    n_pass_c = np.array([len(cf[i]["j"]) for i in common])
    n_pass_p = np.array([len(pf[i]["j"]) for i in common])
    same_set = np.array([cf[i]["j"] == pf[i]["j"] for i in common])
    only_c = np.array([len(set(cf[i]["j"]) - set(pf[i]["j"])) for i in common])
    only_p = np.array([len(set(pf[i]["j"]) - set(cf[i]["j"])) for i in common])
    ovc = np.array([cf[i]["peak_overlap"] for i in common])
    ovp = np.array([pf[i]["peak_overlap"] for i in common])
    res["S1_screen"] = dict(
        total_pass_cpp=int(n_pass_c.sum()),
        total_pass_py=int(n_pass_p.sum()),
        fz_with_pass_cpp=int((n_pass_c > 0).sum()),
        fz_with_pass_py=int((n_pass_p > 0).sum()),
        fz_identical_pass_set=int(same_set.sum()),
        pairs_pass_cpp_only=int(only_c.sum()),
        pairs_pass_py_only=int(only_p.sum()),
        peak_overlap_equal_frac=float((ovc == ovp).mean()),
        peak_overlap_mean_cpp=float(ovc.mean()),
        peak_overlap_mean_py=float(ovp.mean()),
        mean_n_pass_per_fz=(float(n_pass_c.mean()), float(n_pass_p.mean())),
    )
    # peak_overlap by local-grid index j (does the difference sit on j = 0, the unperturbed point?)
    res["S1_by_j"] = [
        dict(
            j=j,
            pass_cpp=int((ovc[:, j] > 0).sum()),
            pass_py=int((ovp[:, j] > 0).sum()),
            mean_overlap_cpp=float(ovc[:, j].mean()),
            mean_overlap_py=float(ovp[:, j].mean()),
        )
        for j in range(ovc.shape[1])
    ]

    # S2 re-cost: on (fz, j) pairs passing in both
    d = []
    for i in common:
        jc = dict(zip(cf[i]["j"], cf[i]["cost"]))
        jp = dict(zip(pf[i]["j"], pf[i]["cost"]))
        for j in set(jc) & set(jp):
            d.append(jp[j] - jc[j])
    d_arr = np.array(d)
    res["S2_recost"] = dict(
        n_pairs=int(d_arr.size),
        mean_py_minus_cpp=float(d_arr.mean()) if d_arr.size else None,
        median=float(np.median(d_arr)) if d_arr.size else None,
        max_abs=float(np.abs(d_arr).max()) if d_arr.size else None,
        frac_abs_lt_1e_6=float((np.abs(d_arr) < 1e-6).mean()) if d_arr.size else None,
        frac_abs_lt_1e_3=float((np.abs(d_arr) < 1e-3).mean()) if d_arr.size else None,
    )

    # S3 spacing filter: kept count per FZ point, and the filter applied to the OTHER side's input
    kept_c = np.array([len(cf[i]["kept"]) for i in common])
    kept_p = np.array([len(pf[i]["kept"]) for i in common])
    both = [i for i in common if cf[i]["j"] == pf[i]["j"] and len(cf[i]["j"]) > 0]
    same_kept = sum(sorted(cf[i]["kept"]) == sorted(pf[i]["kept"]) for i in both)
    res["S3_filter"] = dict(
        total_kept_cpp=int(kept_c.sum()),
        total_kept_py=int(kept_p.sum()),
        fz_same_kept_set=int(sum(sorted(cf[i]["kept"]) == sorted(pf[i]["kept"]) for i in common)),
        fz_with_identical_screen=len(both),
        fz_same_kept_among_identical_screen=int(same_kept),
        mean_kept_per_fz=(float(kept_c[n_pass_c > 0].mean()), float(kept_p[n_pass_p > 0].mean())),
        frac_kept_of_pass=(float(kept_c.sum() / max(n_pass_c.sum(), 1)),
                           float(kept_p.sum() / max(n_pass_p.sum(), 1))),
    )  # fmt: skip
    # distribution of passing counts per FZ point: C++ vs Python (histogram 0..9)
    res["pass_hist"] = dict(
        cpp=np.bincount(n_pass_c, minlength=10).tolist(),
        py=np.bincount(n_pass_p, minlength=10).tolist(),
    )
    res["n_obs_cpp_mean"] = float(np.mean([cf[i].get("n_obs", 0) for i in common]))
    fn = OUT / f"compare{a.py_tag}_v{a.voxel}.json"
    fn.write_text(json.dumps(res, indent=1, default=float) + "\n")
    print(json.dumps({k: v for k, v in res.items() if k != "S1_by_j"}, indent=1, default=float))


def s3_check(a: argparse.Namespace) -> None:
    """S3 on C++'s own inputs: feed the C++ trace's screened (j, cost) lists, rebuilt as
    local_grid[j] @ fz[i], to the library _spacing_filter and compare with C++'s kept lists."""
    os.chdir(EX)
    from icenine import orientation_search as OS
    from icenine import reconstructor as RC
    from icenine.config_file import ConfigFile
    from icenine.experimental_data import ExperimentalData
    from icenine.sampling import generate_local_grid_multi_level

    cfg = ConfigFile.from_file(str(HERE / "configs" / "ReconstructPhaseD_mc_clean.config"))
    exp = ExperimentalData.from_binary_memmap(
        EX / "ScatteringData_PhaseD" / "stacks" / "clean.npy", 180, 2
    )
    setup = RC.setup_reconstruction(cfg, exp_data=exp)
    sym = OS.get_symmetry_quaternions(setup.exp_setup.get_sample_symmetry())
    c = json.loads((OUT / f"cpp_v{a.voxel}.json").read_text())
    rad = c["levels"][0]["diameter_rad"]
    lg = generate_local_grid_multi_level(rad, 0, 1)
    same_set = same_order = n = 0
    diffs = []
    for r in c["fz"]:
        if len(r["j"]) < 2:
            continue
        n += 1
        cands = [
            OS.SearchCandidate(orientation=lg[j] @ setup.fz_orientations[r["fz"]], cost=cost)
            for j, cost in zip(r["j"], r["cost"])
        ]
        j_of = {id(x): j for x, j in zip(cands, r["j"])}
        kept = [j_of[id(x)] for x in OS._spacing_filter(cands, rad, sym)]
        same_set += sorted(kept) == sorted(r["kept"])
        same_order += kept == r["kept"]
        if sorted(kept) != sorted(r["kept"]) and len(diffs) < 10:
            diffs.append(dict(fz=r["fz"], j=r["j"], cost=r["cost"], cpp=r["kept"], py=kept))
    res = dict(
        voxel=a.voxel,
        n_cliques_with_ge2=n,
        same_kept_set=same_set,
        same_kept_order=same_order,
        examples=diffs,
    )
    (OUT / f"s3_on_cpp_inputs_v{a.voxel}.json").write_text(json.dumps(res, indent=1) + "\n")
    print({k: v for k, v in res.items() if k != "examples"})


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["py", "cpp", "compare", "s3"])
    ap.add_argument("--voxel", type=int, default=15901)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--full", dest="level0_only", action="store_false", help="run all levels")
    ap.add_argument("--log", default=None)
    ap.add_argument("--cpp-rule", action="store_true", help="probe: C++ CountQualifiedPeaks rule")
    ap.add_argument("--centroid", action="store_true", help="probe: screen at the voxel centroid")
    ap.add_argument("--py-tag", default="", help="compare: Python file tag, e.g. _cpprule")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    {"py": run_py, "cpp": parse_cpp, "compare": compare, "s3": s3_check}[a.cmd](a)


if __name__ == "__main__":
    main()
