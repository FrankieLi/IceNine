"""AR5 step 15: run an older step's run or pass script on any sample of this step, including
``held3``, with this step's set-up (``band_common``) in front.

  RESEED_SAMPLE=<sample> BAND_CACHE=<dir> uv run python scripts/autoresearch/ar5/band_launch.py \
      <script.py> [args ...]

Only these scripts are accepted: ``fit_run.py`` (arm ``S8``: ``--name S8``), ``fit_pass.py``,
``polish_run.py`` (arm ``J``), ``order_pass.py`` (aligned pass of ``J``), ``polish_pass.py``. They
are run unchanged (``runpy``, as ``__main__``) after ``band_common`` has refused ``held2``, checked
this step's freeze for ``held`` and ``held3``, made ``held3`` known to the AR3 modules and pointed
every cache variable at this step's cache.
"""

import runpy
import sys

import band_common as BC  # noqa: E402  (first)

ALLOWED = ("fit_run.py", "fit_pass.py", "polish_run.py", "order_pass.py", "polish_pass.py")


def main() -> None:
    if len(sys.argv) < 2 or sys.argv[1] not in ALLOWED:
        raise SystemExit(f"band_launch.py: the script must be one of {', '.join(ALLOWED)}")
    script = sys.argv[1]
    sys.argv = [script] + sys.argv[2:]
    runpy.run_path(str(BC.HERE / script), run_name="__main__")


if __name__ == "__main__":
    main()
