"""AR5 step 28: run step 22's run, pass or arm script, or step 21's S8 scripts, on a sample with
this step's set-up (``pxf_common``) in front, with rule ``PXF``, the plain indexer or nothing.

  RESEED_SAMPLE=<sample> SCREEN_CACHE=<dir> DEFER_CACHE=<dir> PXF_ARM=pxf|old|none \
      PXF_INDEX_LOG=<file.json.gz> PXF_STEP_LOG=<file.json.gz> uv run python \
      scripts/autoresearch/ar5/pxf_launch.py <script.py> [args ...]

Only ``defer_run.py`` (arm ``J5D``), ``high500_pass.py`` (``--name J5D --aligned 1``),
``screen_run.py`` (step 21's ``--arm J5S``, on a sub-cache of links), ``fit_run.py`` (arm ``S8``)
and ``fit_pass.py`` are accepted. They run unchanged (``runpy``, as ``__main__``).

``PXF_ARM=pxf``: rule ``PXF``. ``icenine.ray_indexing.index_voxel`` (the old indexer) is wrapped
(``pxf_rules.make_index_voxel``: it passes the controller's ``max_pairs``, 50,000, or 600,000 in a
fallback, and logs ``n_pairs``, ``n_votes``, ``n_rays`` and seconds of every call) and the seed
step of the pick (``defer_pick.DeferPick.e1_step``) is wrapped on the class
(``pxf_rules.make_e1_step``): a step that fails the seed test after a subsampled indexing runs
again with the full sample. ``PXF_ARM=old``: the indexer is unchanged and only logged (plain
J5SD; the boundary pass and the screen arm, which must log no call). ``none`` changes nothing
(arm ``S8``). One thread is the caller's business (``OMP_NUM_THREADS=1``, ``MKL_NUM_THREADS=1``).
"""

import atexit
import gzip
import json
import os
import runpy
import sys
from typing import Any, Dict, List

import pxf_common as PC  # noqa: E402  (first)
import pxf_rules as PR  # noqa: E402

ALLOWED = ("defer_run.py", "high500_pass.py", "screen_run.py", "fit_run.py", "fit_pass.py")


def _hit_of(self: Any, out: Any, vertices: Any, phase_index: Any) -> float:
    """The library's seed-test hit ratio of a step's result (``BFSReconstruction._index_seed``)."""
    import numpy as np

    info = self.rec.evaluate_overlap(np.asarray(out.orientation, np.float64), vertices, phase_index)
    return float(self.bfs._ratios(info)[0])


def _threshold_of(self: Any) -> float:
    return float(self.bfs.setup.config.min_acceleration_threshold)


def main() -> None:
    if len(sys.argv) < 2 or sys.argv[1] not in ALLOWED:
        raise SystemExit(f"pxf_launch.py: the script must be one of {', '.join(ALLOWED)}")
    arm = os.environ.get("PXF_ARM", "")
    if arm not in ("pxf", "old", "none"):
        raise SystemExit("pxf_launch.py: set PXF_ARM=pxf, old or none")
    if arm != "none":
        from icenine import ray_indexing as RI

        ctl = PR.Controller()
        steps: List[Dict[str, Any]] = []
        RI.index_voxel = PR.make_index_voxel(RI.index_voxel, ctl, arm)  # type: ignore[assignment]
        if arm == "pxf":
            import defer_pick as DP

            DP.DeferPick.e1_step = PR.make_e1_step(  # type: ignore[assignment]
                DP.DeferPick.e1_step, ctl, steps, _hit_of, _threshold_of
            )
        out = os.environ.get("PXF_INDEX_LOG", "")
        out_steps = os.environ.get("PXF_STEP_LOG", "")
        script_name = sys.argv[1]

        def dump() -> None:
            if out:
                with gzip.open(out, "wt") as f:
                    json.dump(
                        dict(
                            arm=arm,
                            script=script_name,
                            sample=PC.SAMPLE,
                            summed=PR.summed(ctl.log),
                            rows=ctl.log,
                        ),
                        f,
                    )
            if out_steps and arm == "pxf":
                with gzip.open(out_steps, "wt") as f:
                    json.dump(dict(sample=PC.SAMPLE, script=script_name, rows=steps), f)

        atexit.register(dump)
    script = sys.argv[1]
    sys.argv = [script] + sys.argv[2:]
    runpy.run_path(str(PC.HERE / script), run_name="__main__")


if __name__ == "__main__":
    main()
