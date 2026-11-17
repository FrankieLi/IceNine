"""GO4d (I36): cache the paired observed rays of a sample so the ray-indexing seed path is cheap.

``SmoothCost.observed_rays`` pairs every spot of the two detectors once per process (more than
15 minutes on ``rand``). ``ensure(S, sample)`` loads the pairing from
``<PHASE_D>/<sample>/observed_rays.npz`` into the reconstructor's own cached ``SmoothCost`` (the
object ``Reconstructor._index_context`` gets from ``smooth_cost_cache``), or builds and saves it.

Usage (from icenine_py/), builds the cache and checks a reload:
  uv run python scripts/grain_opt/ray_cache.py --sample corr3
"""

import argparse
import json
import time
from pathlib import Path
from typing import Any, Tuple

import common as CM  # noqa: I001
import numpy as np


def cache_path(sample: str) -> Path:
    return CM.PHASE_D / sample / "observed_rays.npz"


def reconstructor_sc(S: Any) -> Any:
    rec = S.bfs.reconstructor
    assert rec.smooth_cost_cache is not None
    eta = rec.setup.exp_setup.get_eta_limit()
    return rec.smooth_cost_cache.get(rec.setup, eta, rec.min_sin_eta)


def ensure(S: Any, sample: str) -> Tuple[Any, float, bool]:
    """(SmoothCost with rays, seconds spent, loaded from the cache)."""
    from icenine.smooth_cost import ObservedRays

    sc = reconstructor_sc(S)
    p = cache_path(sample)
    t0 = time.perf_counter()
    if sc._rays is not None:
        return sc, 0.0, False
    if p.exists():
        z = np.load(p)
        sc._rays = ObservedRays(
            frame=z["frame"],
            dirn=z["dirn"],
            box0=z["box0"],
            box1=z["box1"],
            n_spots=tuple(int(x) for x in z["n_spots"]),
            n_candidates=int(z["n_candidates"]),
        )
        return sc, time.perf_counter() - t0, True
    r = sc.observed_rays()
    np.savez(
        p,
        frame=r.frame,
        dirn=r.dirn,
        box0=r.box0,
        box1=r.box1,
        n_spots=np.asarray(r.n_spots),
        n_candidates=np.asarray(r.n_candidates),
    )
    return sc, time.perf_counter() - t0, False


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    a = ap.parse_args()
    out = CM.OUT / "go4d" / f"raycache_{a.sample}"
    S = CM.Sample(out, start="centre", intensity=False, sample=a.sample)
    sc, t_build, loaded = ensure(S, a.sample)
    r = sc._rays
    # reload into a fresh holder and compare
    z = np.load(cache_path(a.sample))
    same = all(np.array_equal(getattr(r, k), z[k]) for k in ("frame", "dirn", "box0", "box1"))
    info = dict(sample=a.sample, rays=int(len(r.frame)), n_spots=list(r.n_spots),
                n_candidates=int(r.n_candidates), seconds=t_build, loaded=loaded,
                reload_equal=bool(same), file=str(cache_path(a.sample)))
    (out / "raycache.json").write_text(json.dumps(info, indent=1) + "\n")
    print(json.dumps(info))


if __name__ == "__main__":
    main()
