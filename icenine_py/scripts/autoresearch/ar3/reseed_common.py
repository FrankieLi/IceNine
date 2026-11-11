"""AR3 step 8 (task ar3h): shared set-up of the reseed study.

Two samples: ``big`` (32 grains, step 7) and ``grad`` (the first gradient sample, 497 grains, steps
3 to 6). The sample is chosen with the environment variable ``RESEED_SAMPLE`` (``big`` or ``grad``;
default ``grad``), which this module turns into ``AR3_SAMPLE`` before ``ar3_common`` is imported.
Import this module first in every ``reseed_*.py`` script.

Stored runs (untracked; in the worktrees of steps 6 and 7, found from this worktree's parent
folder) are read, never written. Results go to ``benchmarks/autoresearch/ar3/reseed/`` (``cache/``
holds large files that are not committed).

A *seed step* is one call of the library's ``index_seed`` at a voxel (a first seed or a reseed).
Its evaluations are the strict cost calls on the costed candidates plus the cost calls of each
refinement (local optimization); attaching rays and indexing are not cost evaluations.
"""

import os

SAMPLE = os.environ.get("RESEED_SAMPLE", "grad")
assert SAMPLE in ("big", "grad"), SAMPLE
os.environ["AR3_SAMPLE"] = "big" if SAMPLE == "big" else "gradient"

import gzip  # noqa: E402
import json  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any, Dict, List  # noqa: E402

import ar3_common as C  # noqa: E402
import gradient_common as G  # noqa: E402
import numpy as np  # noqa: E402

REPO = C.ICE.parent
WORKTREES = REPO.parent  # sibling worktrees hold the untracked stored runs
OUT = C.ICE / "benchmarks" / "autoresearch" / "ar3" / "reseed"
CACHE = OUT / "cache"
CACHE.mkdir(parents=True, exist_ok=True)

_BIG_WT = WORKTREES / "agent-a38ba6f5a25aa033e/icenine_py/benchmarks/autoresearch/ar3/big"
_GRAD_WT = WORKTREES / "agent-a6db3cabf8ab230e9/icenine_py/benchmarks/autoresearch/ar3"
STORED = {
    # sample -> run -> run folder (run.json, records.json.gz, recon.mic.gz, [seed_log.json.gz])
    "big": dict(S5=_BIG_WT / "seed/runs/S5", D4=_BIG_WT / "fine/runs/D4"),
    "grad": dict(S5=_GRAD_WT / "seed/runs/S5", D4=_GRAD_WT / "fine/runs/D4"),
}
SEED_TEST = 0.80  # the seed test (hit ratio)
RIGHT = G.WRONG  # degrees: a voxel is right when its error is at most this (step 7's page)
NCOST_REFINE = dict(S5=(20, 5), D4=(5, 2))  # (candidates costed, candidates refined)
LEVELS = G.LEVELS


def read_gz(p: Path) -> Any:
    with gzip.open(p, "rt") as f:
        return json.load(f)


def write_gz(obj: Any, p: Path) -> None:
    with gzip.open(p, "wt", compresslevel=6) as f:
        f.write(json.dumps(obj, separators=(",", ":"), default=float))


def stored_run(run: str) -> Dict[str, Any]:
    d = STORED[SAMPLE][run]
    return json.loads((d / "run.json").read_text())


def seed_steps(run: str) -> List[Dict[str, Any]]:
    """The seed steps of a stored run: ``voxel``, ``reseed``, ``stored_evals`` and, for D4,
    whether the step was followed by a grid fallback (``grid``: its stored evaluations then include
    the grid search; the index evaluations are not separable and are not used for the check).

    S5: ``seed_log.json.gz`` (one row per seed step, no grid fallbacks in S5 on either sample is
    checked by the caller). D4: the reseed list of ``run.json`` plus the first seeds of
    ``records.json.gz`` (``source == "seed"``, not reseeded)."""
    d = STORED[SAMPLE][run]
    if run == "S5":
        rows = read_gz(d / "seed_log.json.gz")
        return [
            dict(
                voxel=int(r["voxel"]),
                reseed=bool(r["reseed"]),
                stored_evals=int(r["evals"]),
                grid=not bool(r.get("found", True)),
            )
            for r in rows
        ]
    st = json.loads((d / "run.json").read_text())["stats"]
    rec = read_gz(d / "records.json.gz")
    out = [
        dict(
            voxel=int(x["voxel"]),
            reseed=True,
            stored_evals=int(x["full_search_evaluations"]),
            grid=int(x["full_search_evaluations"]) > 100_000,
        )
        for x in st["reseeds"]
    ]
    for k, x in rec.items():
        if x["source"] == "seed" and not x["reseeded"]:
            out.append(
                dict(
                    voxel=int(k),
                    reseed=False,
                    stored_evals=int(x["n_evals"]),
                    grid=int(x["n_evals"]) > 100_000,
                )
            )
    return out


def make_env() -> Any:
    """The step 5 environment (reconstructor, neighbour relation, truth) of the chosen sample."""
    import boundary_common as B

    return B.Env("grad")
