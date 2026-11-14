"""AR3 step 3: the gradient sample's paths, its truth and the gzipped record helpers."""

import gzip
import json
from pathlib import Path
from typing import Any, Dict, Tuple

import ar3_common as C
import numpy as np

GRAD = C.OUT / "gradient"
SAMPLE = GRAD / "sample"
TRUTH_MIC = SAMPLE / "gradient_truth.mic"
GRID_MIC = SAMPLE / "gradient_grid.mic"
BIG = C.SAMPLE_KIND in ("big", "mid", "held", "held2")
# AR3 step 7: the 32-grain sample (env AR3_SAMPLE=big); step 11: AR3_SAMPLE=held, the fourth sample
# AR3 step 9: AR3_SAMPLE=mid, the third sample, is built and stored the same way
STACK = f"{C.SAMPLE_KIND}_clean.npy" if BIG else "grad_clean.npy"
WRONG = 1.0  # degrees: "wrong" is an error above this
HURT = 0.05  # degrees: a voxel is hurt / helped by more than this
LEVELS = (0.0, 0.2, 0.5, 1.0, 2.0, 4.0) if BIG else (0.0, 0.1, 0.2, 0.5, 1.0)
GRAINMAP_BIG = SAMPLE / "grainmap.npy"


def dump_gz(obj: Any, path: Path) -> None:
    with gzip.open(path, "wt", compresslevel=6) as f:
        f.write(json.dumps(obj, separators=(",", ":"), default=float))


def load_gz(path: Path) -> Any:
    with gzip.open(path, "rt") as f:
        return json.load(f)


class Truth:
    """The gradient truth over the layer's voxel index, as rendered (read back from the .mic)."""

    def __init__(self) -> None:
        import sys

        sys.path.insert(0, str(C.ICE / "scripts" / "autoresearch" / "ar1"))
        import misses as M
        from icenine.mic_file import MicFile
        from icenine.ray_indexing import proper_rotations
        from icenine.symmetry import create_cubic_symmetry

        self.sym = proper_rotations(create_cubic_symmetry(3.6))
        mic = MicFile.read(str(TRUTH_MIC))
        self.R = np.array([v.orientation for v in mic.voxels], dtype=np.float64)
        self.xy = np.array([v.position[:2] for v in mic.voxels], dtype=float)
        self.gmap = np.load(GRAINMAP_BIG if BIG else M.GRAINMAP)
        info = json.loads((SAMPLE / "sample.json").read_text())
        self.gids = [g["id"] for g in info["grains"]]
        self.level_of = {g["id"]: g["level_deg"] for g in info["grains"]}
        self.Rc = {g["id"]: np.array(g["R_centre"]) for g in info["grains"]}
        self.level = np.array([self.level_of[int(g)] for g in self.gmap])
        self.R_centre = np.stack([self.Rc[int(g)] for g in self.gmap])
        self.n = len(self.R)
        self.sha = info["truth_sha256"]

    def err(self, R: np.ndarray, idx: Any = slice(None)) -> np.ndarray:
        from icenine.ray_indexing import misorientation_deg

        return np.asarray(misorientation_deg(R, self.R[idx], self.sym), dtype=float)

    def mis(self, Ra: np.ndarray, Rb: np.ndarray) -> np.ndarray:
        from icenine.ray_indexing import misorientation_deg

        return np.asarray(misorientation_deg(Ra, Rb, self.sym), dtype=float)


def read_map(path: Path, n: int) -> Tuple[np.ndarray, np.ndarray]:
    """(orientations (n, 3, 3), have (n,)) of a run's recon.mic (plain, or gzipped ``.mic.gz``)."""
    import shutil
    import tempfile

    from icenine.mic_file import MicFile

    if path.suffix == ".gz":
        with tempfile.TemporaryDirectory() as td:
            plain = Path(td) / path.stem
            with gzip.open(path, "rb") as fi, open(plain, "wb") as fo:
                shutil.copyfileobj(fi, fo)
            return read_map(plain, n)
    mic = MicFile.read(str(path))
    R = np.array([v.orientation for v in mic.voxels], dtype=np.float64)
    assert len(R) == n
    return R, np.ones(n, dtype=bool)


def tag_dict(**kw: Any) -> Dict[str, Any]:
    return dict(kw)
