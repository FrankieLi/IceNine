"""Peak memory (max RSS) of one no-start search with the batched level screen off and on
(Task A fixes, before making ``BatchedScreen`` the default).

One process per run, single worker (torch threads 1), Phase D clean full-sample images, Q-max 8
configuration (the default lenient overlap rule unless ``--opt`` says otherwise). The switch is set
explicitly on ``search_params``, so the run does not depend on the library default.

  uv run python scripts/eval_reduction/screen_memory.py run --voxel 15901 --batched 0
  uv run python scripts/eval_reduction/screen_memory.py run --voxel 15901 --batched 1
  uv run python scripts/eval_reduction/screen_memory.py summary

Peak RSS is ``getrusage(RUSAGE_SELF).ru_maxrss`` (bytes on macOS, KiB on Linux), read after the
setup (images, sample, cost functions) and again after the search; it is a high-water mark of the
whole process, so it includes the Python and torch runtime and the touched pages of the
memory-mapped image stack.
"""

import argparse
import contextlib
import io
import json
import os
import resource
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
EX = ICE.parent / "Examples" / "Example2.ManyGrains"
OUT = ICE / "benchmarks" / "eval_reduction" / "A" / "memory"
CONFIG = ICE / "scripts" / "phase_d" / "configs" / "ReconstructPhaseD_mc_clean.config"
STACK = EX / "ScatteringData_PhaseD" / "stacks" / "clean.npy"
sys.path.insert(0, str(ICE / "scripts" / "common"))
sys.path.insert(0, str(ICE / "scripts" / "findoptimal_robustness"))


def max_rss_mib() -> float:
    """Peak resident set size of this process so far, in MiB."""
    r = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return r / 2**20 if sys.platform == "darwin" else r / 2**10


def run_one(voxel: int, batched: bool, opts: List[str], seed: int, out_dir: Path) -> Dict[str, Any]:
    out_dir = out_dir.resolve()
    os.chdir(EX)
    import torch

    torch.set_num_threads(1)
    from icenine import reconstructor as RC
    from icenine.config_file import ConfigFile
    from icenine.experimental_data import ExperimentalData
    from icenine.mic_file import MicFile
    import preflight
    from csl import reduced_misorientation_deg

    sys.path.insert(0, str(HERE))
    from eval_ledger import apply_opts

    cfg = ConfigFile.from_file(str(CONFIG))
    exp = ExperimentalData.from_binary_memmap(STACK, 180, 2)
    setup = RC.setup_reconstruction(cfg, exp_data=exp)
    apply_opts(setup, opts)
    setup.search_params.batched_screen = batched
    rec = RC.AdaptiveVoxelReconstructor(setup)
    v = setup.sample.get_mic().voxels[voxel]
    verts = RC._get_voxel_vertices(v)
    R_true = np.asarray(
        MicFile.read("SimInput/rand_500grains_1mm_neworient_s0.mic").voxels[voxel].orientation,
        dtype=np.float64,
    )
    rss_setup = max_rss_mib()

    pf = preflight.preflight()
    try:
        preflight.require_quiet(info=pf)
        quiet = True
    except preflight.MachineBusyError as e:
        quiet = False
        pf = dict(pf, not_quiet=str(e)[:300])
    rng = np.random.default_rng(seed)
    t0 = time.perf_counter()
    with contextlib.redirect_stdout(io.StringIO()):
        res = rec.reconstruct_voxel(voxel_vertices=verts, phase_index=v.phase, rng=rng)
    wall = time.perf_counter() - t0
    rss_end = max_rss_mib()
    out = dict(
        voxel=voxel,
        batched_screen=batched,
        screen_chunk=int(setup.search_params.screen_chunk),
        opts=opts,
        seed=seed,
        preflight=pf,
        label="single-worker" + ("" if quiet else ", contended"),
        quiet_preflight=quiet,
        wall_s=wall,
        evals_total=int(rec.last_eval_counts[2]),
        err_deg=float(reduced_misorientation_deg(np.asarray(res.orientation), R_true)),
        max_rss_after_setup_mib=rss_setup,
        max_rss_after_search_mib=rss_end,
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    fn = out_dir / f"batched{int(batched)}_v{voxel}.json"
    fn.write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(
        f"voxel {voxel} batched={batched}: max RSS {rss_setup:.1f} MiB after setup, "
        f"{rss_end:.1f} MiB after search; {wall:.1f} s; {out['label']} -> {fn}"
    )
    return out


def summary(d: Path) -> None:
    sys.path.insert(0, str(ICE / "scripts" / "common"))
    import doc_tables as DT

    runs = [json.loads(p.read_text()) for p in sorted(d.glob("batched*_v*.json"))]
    rows = []
    by: Dict[Any, Dict[str, Any]] = {}
    for r in runs:
        by[(r["voxel"], r["batched_screen"])] = r
        rows.append(
            dict(
                voxel=r["voxel"],
                batched_screen="on" if r["batched_screen"] else "off",
                screen_chunk=r["screen_chunk"] if r["batched_screen"] else 1,
                rss_setup_mib=r["max_rss_after_setup_mib"],
                rss_peak_mib=r["max_rss_after_search_mib"],
                evals=r["evals_total"],
                err_deg=r["err_deg"],
                wall_s=r["wall_s"],
                label=r["label"],
            )
        )
    diffs = {}
    for (vox, b), r in by.items():
        if b and (vox, False) in by:
            off = by[(vox, False)]
            diffs[str(vox)] = dict(
                rss_peak_diff_mib=r["max_rss_after_search_mib"] - off["max_rss_after_search_mib"],
                rss_peak_ratio=r["max_rss_after_search_mib"] / off["max_rss_after_search_mib"],
                evals_equal=r["evals_total"] == off["evals_total"],
                err_equal=r["err_deg"] == off["err_deg"],
            )
    (d / "summary.json").write_text(
        json.dumps(dict(rows=rows, diffs=diffs), indent=1, default=float) + "\n"
    )
    T = {
        "memory_batched_screen": DT.markdown_table(
            rows,
            list(rows[0].keys()),
            formats={
                "voxel": "d",
                "rss_setup_mib": ".1f",
                "rss_peak_mib": ".1f",
                "evals": ",d",
                "err_deg": ".4f",
                "wall_s": ".1f",
            },
        )
    }
    DT.write_tables(d / "tables.md", T)
    print(json.dumps(diffs, indent=1))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["run", "summary"])
    ap.add_argument("--voxel", type=int, default=15901)
    ap.add_argument("--batched", type=int, choices=[0, 1])
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--opt", action="append", default=[])
    ap.add_argument("--dir", default=str(OUT))
    a = ap.parse_args()
    if a.cmd == "run":
        run_one(a.voxel, bool(a.batched), a.opt, a.seed, Path(a.dir))
    else:
        summary(Path(a.dir))


if __name__ == "__main__":
    main()
