"""AR5 step 29: run step 22's run, pass or arm script, or step 21's S8 scripts, on a sample with
this step's set-up (``ru_common``) in front, with rule ``RU``, rule ``PXF`` (step 28's), the plain
indexer or nothing.

  RESEED_SAMPLE=<sample> SCREEN_CACHE=<dir> DEFER_CACHE=<dir> RU_ARM=ru|pxf|old|none \
      PXF_INDEX_LOG=<file.json.gz> PXF_STEP_LOG=<file.json.gz> RU_STEP_LOG=<file.json.gz> \
      uv run python scripts/autoresearch/ar5/ru_launch.py <script.py> [args ...]

The scripts accepted are step 28's (``pxf_launch.ALLOWED``); they run unchanged (``runpy``, as
``__main__``).

``RU_ARM=pxf``: step 28's rule ``PXF`` exactly (``pxf_rules.make_index_voxel`` and
``make_e1_step``). ``RU_ARM=ru``: rule ``PXF`` plus rule ``RU`` in front of every seed step: the
seed step of the pick that the library's ``BFSReconstruction._index_seed`` calls
(``seed_pick.SeedPick._index_seed``, which every run script's pick inherits) is wrapped on the class
(``ru_rules.make_index_seed``). ``RU_ARM=old``: the indexer is unchanged and only logged (the
boundary pass and the screen arm, which must log no call). ``none`` changes nothing (arm ``S8``).
One thread is the caller's business (``OMP_NUM_THREADS=1``, ``MKL_NUM_THREADS=1``).
"""

import atexit
import gzip
import json
import os
import runpy
import sys
from typing import Any, Dict, List, Optional

import ru_common as RC  # noqa: E402  (first: guards, freeze, cache, held8)
import ru_rules as RR  # noqa: E402

import numpy as np  # noqa: E402
import pxf_launch as PL  # noqa: E402
import pxf_rules as PR  # noqa: E402

ALLOWED = PL.ALLOWED


class _Geometry:
    """The voxel centres and a cKDTree of them, built once per run from the pick's sample."""

    def __init__(self, mic: Any) -> None:
        from scipy.spatial import cKDTree

        self.pos = np.asarray(mic.positions.numpy(), np.float64)
        self.tree = cKDTree(self.pos)


def make_ctx() -> Any:
    from icenine import ray_indexing as RI
    from icenine.reconstructor import ReconstructionState

    geo: Dict[str, Any] = {}

    def ctx(self: Any, vertices: Any, phase_index: int) -> Dict[str, Any]:
        if "g" not in geo:
            geo["g"] = _Geometry(self.mic)
        g = geo["g"]
        v = int(self.cur)
        voxels = self.mic.voxels
        near = RR.nearest_accepted(
            g.tree,
            g.pos,
            v,
            lambda j: voxels[j].reconstruction_id == ReconstructionState.FITTED,
        )
        assert all(  # the rule's validity check: a start is a voxel accepted at this moment
            voxels[j].reconstruction_id == ReconstructionState.FITTED and j != v for j in near
        )
        rec = self.rec
        cf = rec._index_context(phase_index)[3]

        def strict(R: np.ndarray) -> float:
            return float(cf.evaluate(np.asarray(R, np.float64), vertices, phase_index).cost)

        def fit(R0: np.ndarray) -> Any:
            r = self.orig_lo(
                voxel_vertices=vertices,
                phase_index=phase_index,
                initial_orientation=np.asarray(R0, np.float64),
                rng=np.random.default_rng(0),
                cma_sigma0_deg=RR.SIGMA0_DEG,
                cma_max_evals=RR.MAX_EVALS,
                optimizer="cma",
            )
            return (
                np.asarray(r.orientation, np.float64),
                int(rec.last_local_optimization_evals),
                float(r.cost),
            )

        def hit_of(R: np.ndarray) -> float:
            info = rec.evaluate_overlap(np.asarray(R, np.float64), vertices, phase_index)
            return float(self.bfs._ratios(info)[0])

        def mis(A: np.ndarray, B: np.ndarray) -> float:
            return float(RI.misorientation_deg(A[None], B[None], self.sym)[0])

        return dict(
            voxel=v,
            reseed=getattr(self, "reseed", False),
            near=near,
            orient_of=lambda j: np.asarray(voxels[j].orientation, np.float64),
            strict_cost=strict,
            fit=fit,
            hit_of=hit_of,
            mis_deg=mis,
        )

    return ctx


def main() -> None:
    if len(sys.argv) < 2 or sys.argv[1] not in ALLOWED:
        raise SystemExit(f"ru_launch.py: the script must be one of {', '.join(ALLOWED)}")
    arm = os.environ.get("RU_ARM", "")
    if arm not in ("ru", "pxf", "old", "none"):
        raise SystemExit("ru_launch.py: set RU_ARM=ru, pxf, old or none")
    if arm != "none":
        from icenine import ray_indexing as RI
        from icenine.ray_indexing import SeedResult

        ctl = PR.Controller()
        steps: List[Dict[str, Any]] = []
        ru_rows: List[Dict[str, Any]] = []
        mode = "pxf" if arm in ("ru", "pxf") else "old"
        RI.index_voxel = PR.make_index_voxel(RI.index_voxel, ctl, mode)  # type: ignore[assignment]
        if arm in ("ru", "pxf"):
            import defer_pick as DP

            DP.DeferPick.e1_step = PR.make_e1_step(  # type: ignore[assignment]
                DP.DeferPick.e1_step, ctl, steps, PL._hit_of, PL._threshold_of
            )
        if arm == "ru":
            import seed_pick as SP

            SP.SeedPick._index_seed = RR.make_index_seed(  # type: ignore[assignment]
                SP.SeedPick._index_seed, ru_rows, make_ctx(), SeedResult
            )
        out = os.environ.get("PXF_INDEX_LOG", "")
        out_steps = os.environ.get("PXF_STEP_LOG", "")
        out_ru: Optional[str] = os.environ.get("RU_STEP_LOG", "")
        script_name = sys.argv[1]

        def dump() -> None:
            if out:
                with gzip.open(out, "wt") as f:
                    json.dump(
                        dict(
                            arm=arm,
                            script=script_name,
                            sample=RC.SAMPLE,
                            summed=PR.summed(ctl.log),
                            rows=ctl.log,
                        ),
                        f,
                    )
            if out_steps and arm in ("ru", "pxf"):
                with gzip.open(out_steps, "wt") as f:
                    json.dump(dict(sample=RC.SAMPLE, script=script_name, rows=steps), f)
            if out_ru and arm == "ru":
                with gzip.open(out_ru, "wt") as f:
                    json.dump(dict(sample=RC.SAMPLE, script=script_name, rows=ru_rows), f)

        atexit.register(dump)
    script = sys.argv[1]
    sys.argv = [script] + sys.argv[2:]
    runpy.run_path(str(RC.PC.HERE / script), run_name="__main__")


if __name__ == "__main__":
    main()
