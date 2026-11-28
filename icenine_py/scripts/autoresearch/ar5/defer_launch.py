"""AR5 step 22: run an older step's run, pass or arm script, or this step's ``defer_run.py``, on any
sample of this step, including ``held6``, with this step's set-up (``defer_common``) in front.

  RESEED_SAMPLE=<sample> SCREEN_CACHE=<dir> DEFER_CACHE=<main cache> uv run python \
      scripts/autoresearch/ar5/defer_launch.py <script.py> [args ...]

Only these scripts are accepted: ``fit_run.py`` (arm ``S8``: ``--name S8``), ``fit_pass.py``,
``high500_run.py`` (arm ``J5``), ``high500_pass.py`` (``--name J5`` or ``--name J5D``, ``--aligned
1``), ``screen_run.py`` (step 21's ``--arm J5S``: the screen pass and the wide repair; run on
``J5D``'s map from a sub-cache of links named as ``J5``'s, so that the script is imported
unchanged) and ``defer_run.py`` (arm ``J5D``). They are run unchanged (``runpy``, as ``__main__``)
after ``defer_common`` has refused ``held2``, checked this step's freeze (every script and library
file pinned), made ``held6`` known to the AR3 modules and pointed every cache variable at
``SCREEN_CACHE``.
"""

import runpy
import sys

import defer_common as DC  # noqa: E402  (first)

ALLOWED = (
    "fit_run.py",
    "fit_pass.py",
    "high500_run.py",
    "high500_pass.py",
    "screen_run.py",
    "defer_run.py",
)


def main() -> None:
    if len(sys.argv) < 2 or sys.argv[1] not in ALLOWED:
        raise SystemExit(f"defer_launch.py: the script must be one of {', '.join(ALLOWED)}")
    script = sys.argv[1]
    sys.argv = [script] + sys.argv[2:]
    runpy.run_path(str(DC.HERE / script), run_name="__main__")


if __name__ == "__main__":
    main()
