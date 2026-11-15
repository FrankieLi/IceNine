"""AR5 step 16: the boundary pass of arm ``J5`` (``order_pass.py`` after ``high500_common``).

  RESEED_SAMPLE=<sample> HIGH500_CACHE=<dir> uv run python \
      scripts/autoresearch/ar5/high500_pass.py \
      --name J5 [--aligned 1]

Step 8's ``order_pass.py`` unchanged (aligned draws, step 5's rule ``hit``) on
``<cache>/J5_<sample>/recon.mic.gz``; it only sets this step's cache and the guards (``held2``
refused, ``held`` and ``held3`` behind this step's freeze) first.
"""

import high500_common as BC  # noqa: E402,F401  (first: cache, guards)

import order_pass as OP  # noqa: E402

if __name__ == "__main__":
    OP.main()
