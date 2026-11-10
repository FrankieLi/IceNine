"""Shared set-up of the AR3 scripts: the pilot-region reconstruction setup (as
``scripts/phase_d/pilot_bfs.py`` builds it), the preflight gate and the run-folder writers.

One thread per process is set before numpy / torch are imported (import this module first)."""

import os

for _k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_k] = "1"

import json  # noqa: E402
import re  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any, Dict, Optional, Tuple  # noqa: E402

import numpy as np  # noqa: E402
import torch  # noqa: E402

torch.set_num_threads(1)

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]  # icenine_py
EX = ICE.parent / "Examples" / "Example2.ManyGrains"
REGION = ICE / "benchmarks" / "phase_d_pilot"
OUT = ICE / "benchmarks" / "autoresearch" / "ar3"
CONFIG = ICE / "scripts" / "autoresearch" / "ar2" / "key_on_fused.config"
STORED_KEY_ON = ICE / "benchmarks" / "autoresearch" / "ar2" / "key" / "runs" / "key_on"
STORED_KEY_ON_FUSED = ICE / "benchmarks" / "autoresearch" / "ar2" / "key" / "runs" / "key_on_fused"
sys.path.insert(0, str(ICE / "scripts" / "common"))
sys.path.insert(0, str(ICE / "scripts" / "eval_reduction"))
sys.path.insert(0, str(ICE / "scripts" / "phase_d"))


def preflight_label(out: Path) -> str:
    """Save the preflight JSON beside the results; return the timing label (single-worker when the
    machine is quiet, else the reason it was refused; the run goes on, labelled contended)."""
    import preflight as PF

    info = PF.preflight()
    (out / "preflight.json").write_text(json.dumps(info, indent=2) + "\n")
    try:
        PF.require_quiet(info=info)
        return "single-worker (timing preflight passed)"
    except PF.MachineBusyError as e:
        busy = ", ".join(
            f"{Path(b['command'].split()[0]).name} at {b['cpu']:.0f}% CPU"
            for b in info.get("busy_processes", [])
        )
        return (
            f"contended: timing preflight refused ({busy or str(e)[:80]}); single run, one worker"
        )


def build(
    config: Path,
    out: Path,
    sample: Optional[Path] = None,
    extra: str = "",
    stack_name: str = "clean.npy",
) -> Tuple[Any, Any, float]:
    """(setup, BFSReconstruction, setup seconds) with clean images; the pilot region unless
    ``sample`` is given (AR3 step 2: the whole layer). ``extra`` is appended to the config.
    ``stack_name``: the uint8 stack under ``stacks/`` (AR3 step 3: ``grad_clean.npy``)."""
    cfg_text = config.read_text()
    sample = sample or REGION / "region_grid.mic"
    cfg_text, n = re.subn(r"(?m)^SampleFilename\s+\S+", f"SampleFilename {sample}", cfg_text)
    assert n == 1 and "BFSRevisitRefit 1" in cfg_text
    if extra:
        cfg_text += "\n" + extra + "\n"
    cfg_path = out / "run.config"
    cfg_path.write_text(cfg_text)
    os.chdir(EX)
    from icenine.config_file import ConfigFile
    from icenine.experimental_data import ExperimentalData
    from icenine.reconstructor import BFSReconstruction, setup_reconstruction

    cfg = ConfigFile.from_file(str(cfg_path))
    stack = EX / "ScatteringData_PhaseD" / "stacks" / stack_name
    t0 = time.time()
    setup = setup_reconstruction(cfg, exp_data=ExperimentalData.from_binary_memmap(stack, 180, 2))
    return setup, BFSReconstruction(setup), time.time() - t0


def write_mic(mic: Any, orientations: Dict[int, np.ndarray], path: Path) -> None:
    for v, R in orientations.items():
        mic.voxels[v].orientation = np.asarray(R, dtype=np.float64).copy()
    mic.write(str(path))
