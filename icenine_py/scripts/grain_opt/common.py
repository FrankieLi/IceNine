"""Shared set-up of the grain-fit scripts: the AR3 reconstruction setup of a sample (strict cost),
its truth, a start map, voxel triangles and neighbours, the spot model and the image source.

Import first: ``ar3_common`` sets one thread per process before numpy / torch load.
"""

import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(ICE / "scripts" / "autoresearch" / "ar3"))

import ar3_common as C  # noqa: E402  (sets threads)
import gradient_common as G  # noqa: E402
import numpy as np  # noqa: E402
import torch  # noqa: E402

OUT = ICE / "benchmarks" / "grain_opt"
PHASE_D = C.EX / "ScatteringData_PhaseD"


def wilson(k: int, n: int, z: float = 1.96) -> str:
    if n == 0:
        return "0 of 0"
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    lo, hi = max(c - h, 0.0) * 100, min(c + h, 1.0) * 100
    fmt = "{:.2f}" if hi < 1 else "{:.1f}"
    return f"{k} of {n} ({fmt.format(lo)} to {fmt.format(hi)}%)"


def err_stats(err: np.ndarray) -> Dict[str, Any]:
    n = len(err)
    return dict(
        n=n,
        median=float(np.median(err)) if n else float("nan"),
        p95=float(np.percentile(err, 95)) if n else float("nan"),
        gt01=wilson(int((err > 0.1).sum()), n),
        wrong=wilson(int((err > 1.0).sum()), n),
        n_gt01=int((err > 0.1).sum()),
        n_wrong=int((err > 1.0).sum()),
    )


class Sample:
    """Everything a grain fit needs on one sample (``grad`` by default, AR3_SAMPLE for others)."""

    def __init__(self, out: Path, start: str = "paint", intensity: bool = True) -> None:
        from icenine import grain_fit as GF
        from icenine.reconstructor import _get_voxel_vertices
        from icenine.smooth_cost import _build_smooth_cost

        out.mkdir(parents=True, exist_ok=True)
        cwd = os.getcwd()
        self.setup, self.bfs, self.t_setup = C.build(
            C.CONFIG, out, sample=G.GRID_MIC, stack_name=G.STACK
        )
        os.chdir(cwd)
        self.mic = self.setup.sample.get_mic()
        self.n = len(self.mic.voxels)
        self.T = G.Truth()
        assert self.T.n == self.n
        self.sym = self.T.sym
        self.cost_fn = self.bfs.reconstructor._make_local_cost_fn()
        self.V = np.stack([_get_voxel_vertices(v).numpy().astype(np.float64) for v in self.mic.voxels])
        self.xy = self.V.mean(1)[:, :2]
        side = self.mic.voxels[0].side_length
        self.nb: List[List[int]] = [
            [int(u) for u in self.mic.get_neighbors(v, radius=2.0 * side)] for v in range(self.n)
        ]
        self.sc = _build_smooth_cost(
            self.setup, self.setup.exp_setup.get_eta_limit(), self.bfs.reconstructor.min_sin_eta
        )
        self.phase = int(self.mic.voxels[0].phase)
        assert all(int(v.phase) == self.phase for v in self.mic.voxels)
        self.model = GF.SpotModel(self.sc, phase=self.phase)
        stack = np.load(str(PHASE_D / "stacks" / G.STACK), mmap_mode="r")
        store = None
        tag = G.STACK.replace("_clean.npy", "")
        sp = PHASE_D / tag / "intensity_clean.npz"
        if intensity and sp.exists():
            z = np.load(sp)
            store = {k: z[k] for k in z.files}
        self.images = GF.ImageSource(stack, 180, store)
        self.R0 = self.start_map(start)

    def start_map(self, start: str) -> np.ndarray:
        if start == "paint":
            m = np.load(G.GRAD / "runs" / "P" / "maps.npz")
            return m["P"].astype(np.float64)
        if start == "truth":
            return self.T.R.copy()
        R, _ = G.read_map(Path(start), self.n)
        return R

    def strict(self, R: np.ndarray, voxels: List[int]) -> np.ndarray:
        """Strict (local, hard) cost of each voxel at its orientation."""
        out = np.empty(len(voxels))
        for i, v in enumerate(voxels):
            info = self.cost_fn.evaluate(
                np.asarray(R[i], np.float64), torch.as_tensor(self.V[v], dtype=torch.float32), self.phase
            )
            out[i] = float(info.cost)
        return out
