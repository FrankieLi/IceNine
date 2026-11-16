"""AR5 step 20: run an older step's run or pass script on any sample of this step, including
``held4``, with this step's set-up (``cmarepair_common``) in front.

  RESEED_SAMPLE=<sample> CMAREPAIR_CACHE=<dir> uv run python \
      scripts/autoresearch/ar5/cmarepair_launch.py <script.py> [args ...]

Only these scripts are accepted: ``fit_run.py`` (arm ``S8``: ``--name S8``), ``fit_pass.py``,
``high500_run.py`` (arm ``J5``) and ``high500_pass.py`` (``--name J5 --aligned 1``). They are run
unchanged (``runpy``, as ``__main__``) after ``cmarepair_common`` has refused ``held2``, checked
this step's freeze for the held samples, made ``held3`` and ``held4`` known to the AR3 modules and
pointed every cache variable at this step's cache.
"""

import runpy
import sys

import cmarepair_common as RCM  # noqa: E402  (first)

ALLOWED = (
    "fit_run.py",
    "fit_pass.py",
    "high500_run.py",
    "high500_pass.py",
)


def main() -> None:
    if len(sys.argv) < 2 or sys.argv[1] not in ALLOWED:
        raise SystemExit(f"cmarepair_launch.py: the script must be one of {', '.join(ALLOWED)}")
    script = sys.argv[1]
    sys.argv = [script] + sys.argv[2:]
    runpy.run_path(str(RCM.HERE / script), run_name="__main__")


if __name__ == "__main__":
    main()
