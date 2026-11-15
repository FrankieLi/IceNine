"""AR5 step 15: the boundary pass of arm ``JB`` (``order_pass.py`` after ``band_common``).

  RESEED_SAMPLE=<sample> BAND_CACHE=<dir> uv run python scripts/autoresearch/ar5/band_pass.py \
      --name JB [--aligned 1]

Step 8's ``order_pass.py`` unchanged (aligned draws, step 5's rule ``hit``) on
``<cache>/JB_<sample>/recon.mic.gz``; it only sets the band cache and the guards (``held2``
refused, ``held`` and ``held3`` behind this step's freeze) first.
"""

import band_common as BC  # noqa: E402,F401  (first: cache, guards)

import order_pass as OP  # noqa: E402

if __name__ == "__main__":
    OP.main()
