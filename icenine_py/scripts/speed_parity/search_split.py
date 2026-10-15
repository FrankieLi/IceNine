"""Wall time of a no-start reconstruct_voxel split into evaluate() vs everything else.

Usage (from icenine_py/): uv run python scripts/speed_parity/search_split.py [config name]
"""

import sys
import time
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import prof_eval as P
from icenine import cost_functions as cf
from icenine.reconstructor import setup_reconstruction, AdaptiveVoxelReconstructor
from icenine.sampling import load_fundamental_zone_file

cfgname = sys.argv[1] if len(sys.argv) > 1 else "Example2.Simulation.config"
from icenine.config_file import ConfigFile

cfg = ConfigFile.from_file(str(P.ex / "ConfigFiles" / cfgname))
fz = load_fundamental_zone_file(cfg.fundamental_zone_filename)
setup = setup_reconstruction(cfg, exp_data=P.exp_data, fz_orientations=fz)
rec = AdaptiveVoxelReconstructor(setup)

acc = {"t": 0.0, "n": 0}
orig = cf.VoxelCostFunction.evaluate


def timed(self, *a, **k):
    t0 = time.perf_counter()
    r = orig(self, *a, **k)
    acc["t"] += time.perf_counter() - t0
    acc["n"] += 1
    return r


cf.VoxelCostFunction.evaluate = timed

t0 = time.perf_counter()
res = rec.reconstruct_voxel(P.verts, P.vox.phase, rng=np.random.default_rng(1))
wall = time.perf_counter() - t0
n = acc["n"]
print(f"config {cfgname}: wall {wall:.1f} s, evaluations {n}, cost {res.cost:.4f}")
import hashlib  # noqa: E402

orient = np.ascontiguousarray(res.orientation, dtype=np.float64)
print(
    f"  result: cost {res.cost!r}, orientation sha1 "
    f"{hashlib.sha1(orient.tobytes()).hexdigest()[:12]}"
)
inside, outside = acc["t"], wall - acc["t"]
share = 100 * inside / wall
print(f"  inside evaluate: {inside:.1f} s ({share:.1f}%), {1e6 * inside / n:.0f} us/eval")
print(
    f"  outside evaluate (search machinery): {outside:.1f} s, {1e6 * outside / n:.1f} us per eval"
)
