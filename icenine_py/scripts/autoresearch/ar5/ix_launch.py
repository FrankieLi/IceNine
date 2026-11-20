"""AR5 step 26: run step 22's run, pass or arm script on a sample with this step's set-up
(``ix_common``) in front, the old or the new indexer.

  RESEED_SAMPLE=<sample> SCREEN_CACHE=<dir> DEFER_CACHE=<same> IX_INDEXER=old|new uv run python \
      scripts/autoresearch/ar5/ix_launch.py <script.py> [args ...]

Only ``defer_run.py`` (arm ``J5D``), ``high500_pass.py`` (``--name J5D --aligned 1``) and
``screen_run.py`` (step 21's ``--arm J5S``, on a sub-cache of links) are accepted. They run
unchanged (``runpy``, as ``__main__``). ``IX_INDEXER=new`` replaces ``icenine.ray_indexing.
index_voxel`` by ``index_voxel_fast`` in this process before the script starts (the scripts and
the library call it through the module); ``old`` changes nothing. One thread is the caller's
business (``OMP_NUM_THREADS=1``, ``MKL_NUM_THREADS=1``).
"""

import os
import runpy
import sys

import ix_common as IC  # noqa: E402  (first)

ALLOWED = ("defer_run.py", "high500_pass.py", "screen_run.py")


def main() -> None:
    if len(sys.argv) < 2 or sys.argv[1] not in ALLOWED:
        raise SystemExit(f"ix_launch.py: the script must be one of {', '.join(ALLOWED)}")
    which = os.environ.get("IX_INDEXER", "")
    if which not in ("old", "new"):
        raise SystemExit("ix_launch.py: set IX_INDEXER=old or IX_INDEXER=new")
    from icenine import ray_indexing as RI

    if which == "new":
        RI.index_voxel = RI.index_voxel_fast  # type: ignore[assignment]
    script = sys.argv[1]
    sys.argv = [script] + sys.argv[2:]
    runpy.run_path(str(IC.HERE / script), run_name="__main__")


if __name__ == "__main__":
    main()
