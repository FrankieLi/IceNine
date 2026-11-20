"""AR5 step 27: run step 22's run, pass or arm script, or step 21's S8 scripts, on a sample with
this step's set-up (``px_common``) in front, with the old or the new indexer setting.

  RESEED_SAMPLE=<sample> SCREEN_CACHE=<dir> DEFER_CACHE=<dir> PX_ARM=old|new|none \
      PX_INDEX_LOG=<file.json.gz> uv run python scripts/autoresearch/ar5/px_launch.py \
      <script.py> [args ...]

Only ``defer_run.py`` (arm ``J5D``), ``high500_pass.py`` (``--name J5D --aligned 1``),
``screen_run.py`` (step 21's ``--arm J5S``, on a sub-cache of links), ``fit_run.py`` (arm ``S8``)
and ``fit_pass.py`` are accepted. They run unchanged (``runpy``, as ``__main__``).

``PX_ARM=old`` and ``new`` wrap ``icenine.ray_indexing.index_voxel`` (the old indexer, in both)
in this process: ``new`` passes ``max_pairs`` 50,000 in the ``IndexParams`` of every call, ``old``
leaves it; both log ``n_pairs``, ``n_votes``, ``n_rays`` and the seconds of every call to
``PX_INDEX_LOG`` at exit. ``none`` changes nothing (arm ``S8``). One thread is the caller's
business (``OMP_NUM_THREADS=1``, ``MKL_NUM_THREADS=1``).
"""

import atexit
import gzip
import json
import os
import runpy
import sys
from typing import Any, Dict, List

import px_common as PC  # noqa: E402  (first)
import px_rules as PR  # noqa: E402

ALLOWED = ("defer_run.py", "high500_pass.py", "screen_run.py", "fit_run.py", "fit_pass.py")


def main() -> None:
    if len(sys.argv) < 2 or sys.argv[1] not in ALLOWED:
        raise SystemExit(f"px_launch.py: the script must be one of {', '.join(ALLOWED)}")
    arm = os.environ.get("PX_ARM", "")
    if arm not in ("old", "new", "none"):
        raise SystemExit("px_launch.py: set PX_ARM=old, new or none")
    if arm != "none":
        from icenine import ray_indexing as RI

        log: List[Dict[str, Any]] = []
        RI.index_voxel = PR.make_index_voxel(RI.index_voxel, arm, log)  # type: ignore[assignment]
        out = os.environ.get("PX_INDEX_LOG", "")
        script_name = sys.argv[1]

        def dump() -> None:
            if not out:
                return
            with gzip.open(out, "wt") as f:
                json.dump(
                    dict(
                        arm=arm,
                        script=script_name,
                        sample=PC.SAMPLE,
                        summed=PR.summed(log),
                        rows=log,
                    ),
                    f,
                )

        atexit.register(dump)
    script = sys.argv[1]
    sys.argv = [script] + sys.argv[2:]
    runpy.run_path(str(PC.HERE / script), run_name="__main__")


if __name__ == "__main__":
    main()
