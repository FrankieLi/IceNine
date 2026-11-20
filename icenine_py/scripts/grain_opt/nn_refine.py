"""GO8b (I25): GNLayerNet as a per-voxel refiner on observed windows.

``NNRefiner(S, weights).propose(v, R_nom)`` returns the network's refined orientation of voxel ``v``
from the nominal ``R_nom`` (or None): peak set, window placement, context and nominal offsets from
``perturbation_sweep.prepare_nominal`` on the ManyGrains geometry (its voxel grid equals the
grain-opt samples' grid); observed windows cut from the sample's binary stack (32 x 32 px at each
peak's detector, frames frame0 +- K; each lit pixel coded 1 + (frame - frame0) + K with the lit frame
nearest frame0), zero-padded to the network's 130 peak slots; R_new = exp([delta]x) R_nom.
"""

import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Optional

import numpy as np
import torch
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))  # scripts/: perturbation_sweep, generate_toy_orientation_dataset


class NNRefiner:
    def __init__(self, S: Any, weights: str) -> None:
        from generate_toy_orientation_dataset import example_dir_for, setup_example
        from perturbation_sweep import load_model

        self.S = S
        self.net, ck = load_model(weights)
        self.M = int(ck.get("n_peaks", 130))
        self.W = int(ck.get("window_size", 32))
        self.K = int(ck.get("frame_half_width", 4))
        ex = example_dir_for("manygrains")
        self.ctx = SimpleNamespace(
            args=SimpleNamespace(max_q=8.0, detectors="all", min_sin_eta=0.3, min_peaks=40,
                                 window_size=self.W, frame_half_width=self.K),
            example_dir=ex, setup=setup_example(ex, max_q=8.0),
        )
        self.stack = S.images.stack
        self.n_omega = int(S.images.n_omega)
        self.n_calls = 0

    def windows(self, prep: Any) -> torch.Tensor:
        W, K, n = self.W, self.K, prep.n
        out = np.zeros((self.M, W, W), dtype=np.uint8)
        dets = prep.obs.det_idx.tolist()
        rows_tot, cols_tot = self.stack.shape[1], self.stack.shape[2]
        for m in range(min(n, self.M)):
            c0, r0, f0 = int(prep.spec.col0[m]), int(prep.spec.row0[m]), int(prep.spec.frame0[m])
            ra, rb, ca, cb = max(r0, 0), min(r0 + W, rows_tot), max(c0, 0), min(c0 + W, cols_tot)
            if ra >= rb or ca >= cb:
                continue
            best = np.full((rb - ra, cb - ca), 99, dtype=np.int64)
            for off in sorted(range(-K, K + 1), key=abs, reverse=True):  # nearest frame wins
                f = f0 + off
                if not 0 <= f < self.n_omega:
                    continue
                lit = np.asarray(self.stack[dets[m] * self.n_omega + f, ra:rb, ca:cb]) > 0
                best[lit] = off
            code = np.where(best != 99, 1 + best + K, 0).astype(np.uint8)
            out[m, ra - r0 : rb - r0, ca - c0 : cb - c0] = code
        return torch.from_numpy(out)[None]

    def propose(self, v: int, R_nom: np.ndarray) -> Optional[np.ndarray]:
        from icenine.orientation_eval import decode_windows
        from perturbation_sweep import prepare_nominal

        prep, why = prepare_nominal(self.ctx, int(v), np.asarray(R_nom, np.float64))
        if prep is None:
            return None
        n = min(prep.n, self.M)
        ctx = torch.zeros(1, self.M, prep.context.shape[-1])
        ctx[0, :n] = prep.context[:n]
        nom = torch.zeros(1, self.M, 3)
        nom[0, :n] = prep.nom_off[:n]
        with torch.no_grad():
            x = decode_windows(self.windows(prep), self.K)
            d, _L = self.net(x, ctx, {"nom_off": nom})
        self.n_calls += 1
        delta = np.radians(d[0].double().numpy())
        return Rotation.from_rotvec(delta).as_matrix() @ np.asarray(R_nom, np.float64)
