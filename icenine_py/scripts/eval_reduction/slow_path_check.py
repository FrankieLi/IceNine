"""Task A-exact: Python slow path of the cost against the C path (voxel 5242, CppEdgeRules on/off).

  uv run python scripts/eval_reduction/slow_path_check.py [--stub]

Without ``--stub``: 12 orientations per arm (screen, complete cost) x switch on/off, dense images,
slow path (``_HAS_C_RASTERIZE`` False) against the C path: pixel and peak counts and quality must
be identical. With ``--stub``: the images are dark and the slow path's triangle counter returns
(0 overlap, 5 lit pixels) for every triangle, so with the switch the slow path lights a detector
only through an in-bounds vertex, as the C path does on a dark image; the sums of
peak-on-detector are printed per arm (equal for the switch on; the default differs by
construction, and the pixel counts of the stub differ from the C path, so the "differing" count
is not used).
"""

import argparse
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import torch
from scipy.spatial.transform import Rotation

sys.path.insert(0, str(Path(__file__).resolve().parent))
import parity_common as PC  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stub", action="store_true")
    a = ap.parse_args()
    ctx = PC.PyContext()
    from icenine import cost_functions as CF
    from icenine.image_data import ImageData

    orig = ctx.setup.exp_data.get_image
    cache = {}

    def dense_get(w: int, d: int) -> ImageData:
        if (w, d) not in cache:
            im = ImageData(2048, 2048, mode="dense")
            if a.stub:
                im._pixels_dense = torch.zeros(2048, 2048)
            else:
                im._pixels_dense = torch.from_numpy(orig(w, d).get_binary_numpy().astype("float32"))
            cache[(w, d)] = im
        return cache[(w, d)]

    ctx.setup.exp_data.get_image = dense_get
    if a.stub:
        real = ImageData.get_triangle_overlap_property

        def stub(self, v0, v1, v2, mode="hard", **k):  # type: ignore[no-untyped-def]
            if not CF._HAS_C_RASTERIZE:
                return torch.tensor(0.0), torch.tensor(5.0)
            return real(self, v0, v1, v2, mode=mode, **k)

        ImageData.get_triangle_overlap_property = stub  # type: ignore[assignment]

    v, phase = 5242, 1
    verts, r0 = ctx.verts(v), ctx.truth(v)
    rng = np.random.default_rng(1)
    diffs = tot = 0
    sums = defaultdict(int)
    t = time.time()
    for arm in ("screen", "full"):
        for edge in (True, False):
            cf = ctx.cost_fn(arm, cpp_edge_rules=edge)
            for i in range(12):
                if i >= 4:
                    rot = Rotation.random(random_state=int(rng.integers(1 << 30))).as_matrix()
                else:
                    rot = Rotation.from_rotvec(rng.normal(size=3) * 0.01 * i).as_matrix() @ r0
                res = []
                for c_on in (True, False):
                    CF._HAS_C_RASTERIZE = c_on
                    o = cf.evaluate(np.asarray(rot, dtype=np.float32), verts, phase)
                    CF._HAS_C_RASTERIZE = True
                    res.append(
                        (o.pixel_overlap, o.pixel_on_detector, o.peak_overlap, o.peak_on_detector)
                        + (float(o.quality),)
                    )
                    sums[(arm, edge, "C" if c_on else "slow")] += o.peak_on_detector
                tot += 1
                diffs += res[0] != res[1]
    print(f"compared {tot} evaluations, C vs slow differing: {diffs}, {time.time() - t:.1f} s")
    print({"/".join(map(str, k)): n for k, n in sums.items()})


if __name__ == "__main__":
    main()
