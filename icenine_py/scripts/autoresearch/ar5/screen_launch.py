"""AR5 step 21: run an older step's run or pass script on any sample of this step, including
``held5``, with this step's set-up (``screen_common``) in front.

  RESEED_SAMPLE=<sample> SCREEN_CACHE=<dir> uv run python \
      scripts/autoresearch/ar5/screen_launch.py <script.py> [args ...]

Only these scripts are accepted: ``fit_run.py`` (arm ``S8``: ``--name S8``), ``fit_pass.py``,
``high500_run.py`` (arm ``J5``) and ``high500_pass.py`` (``--name J5 --aligned 1``). They are run
unchanged (``runpy``, as ``__main__``) after ``screen_common`` has refused ``held2``, checked this
step's freeze for the held samples, made ``held3``, ``held4`` and ``held5`` known to the AR3
modules and pointed every cache variable at this step's cache.
"""

import runpy
import sys

import screen_common as SCM  # noqa: E402  (first)

ALLOWED = (
    "fit_run.py",
    "fit_pass.py",
    "high500_run.py",
    "high500_pass.py",
)


def main() -> None:
    if len(sys.argv) < 2 or sys.argv[1] not in ALLOWED:
        raise SystemExit(f"screen_launch.py: the script must be one of {', '.join(ALLOWED)}")
    script = sys.argv[1]
    sys.argv = [script] + sys.argv[2:]
    runpy.run_path(str(SCM.HERE / script), run_name="__main__")


if __name__ == "__main__":
    main()
