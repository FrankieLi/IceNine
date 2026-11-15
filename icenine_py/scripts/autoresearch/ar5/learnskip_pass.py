"""AR5 step 18: the boundary pass of a run of this step (``order_pass.py`` after the set-up).

  RESEED_SAMPLE=<sample> LEARNSKIP_CACHE=<dir> uv run python \
      scripts/autoresearch/ar5/learnskip_pass.py --name J5L [--aligned 1]

Step 8's ``order_pass.py`` unchanged (aligned draws, step 5's rule ``hit``) on
``<cache>/<name>_<sample>/recon.mic.gz``; it only sets this step's cache and the guards (``held2``
refused, ``held``, ``held3`` and ``held4`` behind this step's freeze) first. Step 17's
``fitskip_pass.py`` for arm ``J5r``, with this step's set-up.
"""

import learnskip_common as LCM  # noqa: E402,F401  (first: cache, guards)

import order_pass as OP  # noqa: E402

if __name__ == "__main__":
    OP.main()
