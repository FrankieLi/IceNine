"""tex1: one run of a Python arm (``RCP``) on the fine layer or on its slice.

  uv run python scripts/autoresearch/tex1/tex1_run_py.py --arm RCP --grid slice --out DIR
  uv run python scripts/autoresearch/tex1/tex1_run_py.py --arm RCP --grid full --out DIR --dry-run

The config is ``tex1_config.python_config`` (the arm's base config with only the grid, the images,
the two side lengths and the output name changed). Clean images from the memory-mapped stack, Q-max
8, one worker, ``default_rng(0)`` for the voxel order, the helpers of AR3 (``ar3_common.build``)
unchanged. The run reads no truth. The freeze is checked first: nothing runs before
``tex1/frozen.json`` is committed.

Writes ``<out>/``: ``arm.config`` and ``run.config``, ``preflight.json``, ``progress.jsonl`` (one
line every 5 minutes: elapsed seconds, voxels with a record and the numeric counters; the launcher's
stall guard watches its size), ``recon.mic``, ``records.json.gz``, ``run.json`` (counts and seconds
only: evaluations in all and in global searches, global searches and their seconds, seconds outside
them, voxels with no output) and last ``run.done``. ``--dry-run`` writes the config under
``<out>`` and prints what would run; it builds nothing and starts nothing.
"""

import argparse
import gzip
import json
import shutil
import sys
import threading
import time
from pathlib import Path
from typing import Any, Dict

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import tex1_config as CF  # noqa: E402
import tex1_freeze as FZ  # noqa: E402
import tex1_paths as P  # noqa: E402

EVAL_KEYS_ALL = (
    "n_evals_seed",
    "n_evals_neighbor",
    "n_evals_revisit",
    "n_evals_restart",
    "n_evals_reseed",
    "n_evals_boundary_pass",
    "n_evals_boundary_score",
    "n_evals_boundary_repair",
)


def numeric(d: Dict[str, Any]) -> Dict[str, Any]:
    return {k: v for k, v in d.items() if isinstance(v, (int, float)) and not isinstance(v, bool)}


def counts_summary(stats: Dict[str, Any], wall_total_s: float) -> Dict[str, Any]:
    """Evaluations (all, and global searches apart), global searches and seconds, from a BFS's
    ``stats`` (counts and seconds only). A global search is a seed or a reseed."""
    n_all = int(sum(int(stats.get(k, 0)) for k in EVAL_KEYS_ALL))
    n_global = int(stats.get("n_evals_seed", 0)) + int(stats.get("n_evals_reseed", 0))
    searches = int(stats.get("n_seeds", 0)) + int(stats.get("n_reseeds", 0))
    t_global = float(stats.get("wall_seed_s", 0.0)) + float(stats.get("wall_reseed_s", 0.0))
    return dict(
        evals_all=n_all,
        evals_by_key={k: int(stats.get(k, 0)) for k in EVAL_KEYS_ALL},
        evals_global=n_global,
        evals_other=n_all - n_global,
        global_searches=searches,
        n_seeds=int(stats.get("n_seeds", 0)),
        n_reseeds=int(stats.get("n_reseeds", 0)),
        wall_global_s=t_global,
        wall_other_s=float(wall_total_s) - t_global,
        wall_per_global_search_s=(t_global / searches) if searches else None,
    )


def gz_file(path: Path) -> None:
    with open(path, "rb") as fi, open(str(path) + ".gz", "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as fo:
            shutil.copyfileobj(fi, fo)
    path.unlink()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--arm", required=True, choices=CF.PYTHON_ARMS)
    ap.add_argument("--grid", required=True, choices=("slice", "full"))
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--max-voxels", type=int, default=0, help="smoke test only")
    ap.add_argument("--no-timing-gate", action="store_true", help="skip the preflight (tests)")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--test-grid", type=Path, default=None, help="code-path test: another grid")
    ap.add_argument("--test-stack", type=Path, default=None, help="code-path test: another stack")
    a = ap.parse_args()
    fine_files = {P.GRID.resolve(), P.SLICE.resolve(), P.STACK.resolve()}
    for t in (a.test_grid, a.test_stack):  # the freeze is skipped only on other data
        if t is not None and t.resolve() in fine_files:
            raise SystemExit("--test-grid and --test-stack are for data other than tex24's")
    grid = a.test_grid or (P.SLICE if a.grid == "slice" else P.GRID)
    stack = a.test_stack or P.STACK
    text = CF.python_config(a.arm, grid, str(P.FRAMES / "500Grains.sim"), f"tex1_{a.arm}")
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "arm.config").write_text(text)
    if a.dry_run:
        print(f"dry run: arm {a.arm} grid {grid} stack {stack} out {a.out}")
        print(
            "config changed keys against the base:",
            CF.changed_keys(CF.BASE_J5S.read_text(), text),
        )
        return 0
    if not (a.test_grid or a.test_stack):  # a code-path test on other data skips the freeze
        FZ.check_or_exit()  # nothing runs before frozen.json is committed
    sys.path.insert(0, str(P.ICE / "scripts" / "autoresearch" / "ar3"))
    import ar3_common as C

    label = "not a timing run" if a.no_timing_gate else C.preflight_label(a.out)
    setup, bfs, t_setup = C.build(a.out / "arm.config", a.out, sample=grid, stack_name=str(stack))
    mic = setup.sample.get_mic()
    p = setup.search_params
    assert p.fused_cost and p.bfs_revisit_refit and p.bfs_reseed and p.local_optimizer == "mc"
    # the library keys of lib_check.json plus steps 22, 28, 29; BoundaryPassScreen's preconditions
    assert p.seed_by_indexing and p.boundary_pass and p.boundary_pass_repair
    assert p.boundary_pass_screen and p.boundary_pass_skip_hit == 0.0
    assert p.bfs_neighbor_cma_first and p.bfs_neighbor_cma_no_fallback
    assert (p.seed_index_n_cost, p.seed_index_n_refine) == (20, 5)
    assert p.seed_index_stop_hit == 0.95 and p.seed_index_cma_first == 200
    assert p.cma_sigma0_deg == 0.05 and p.cma_neighbor_max_evals == 60
    assert (p.seed_rescue_sigma_deg, p.seed_rescue_evals) == (0.25, 1000)
    assert p.seed_index_mc_defer == 1 and p.seed_reuse_known == 1
    assert p.ray_index_max_pairs == 50000 and p.ray_index_fallback_full == 1
    from icenine.reconstructor import check_bfs_keys

    check_bfs_keys(p)  # accepts all 16 keys together
    n = len(mic.voxels)
    if not a.test_grid:
        assert n == (P.SLICE_N if a.grid == "slice" else P.N_FINE), n
    import numpy as np

    t0 = time.time()

    def progress() -> None:
        while True:
            time.sleep(300.0)
            try:
                with open(a.out / "progress.jsonl", "a") as f:
                    f.write(
                        json.dumps(
                            dict(
                                t=time.time(),
                                elapsed_s=time.time() - t0,
                                n_records=len(bfs.records),
                                stats=numeric(bfs.stats),
                            )
                        )
                        + "\n"
                    )
            except Exception:  # a progress write must never stop the run
                pass

    threading.Thread(target=progress, daemon=True).start()
    order = bfs.reconstruct_sample(
        output_mic=str(a.out / "recon.mic"),
        max_voxels=a.max_voxels or None,
        rng=np.random.default_rng(0),
    )
    wall = time.time() - t0
    (a.out / "records.json").write_text(
        json.dumps({str(k): vars(v) for k, v in bfs.records.items()}) + "\n"
    )
    gz_file(a.out / "records.json")
    rec_evals = int(sum(int(getattr(r, "n_evals", 0)) for r in bfs.records.values()))
    zero = int(sum(1 for v in mic.voxels if abs(float(np.linalg.det(v.orientation))) < 1e-6))
    st = {k: v for k, v in bfs.stats.items() if k != "reseeds"}
    info: Dict[str, Any] = dict(
        arm=a.arm,
        grid=a.grid,
        n_voxels=n,
        n_processed=len(order),
        n_without_orientation=zero,
        evals_records_sum=rec_evals,
        wall_total_s=wall,
        setup_s=t_setup,
        timing_label=label,
        torch_threads=C.torch.get_num_threads(),
        stats=st,
        **counts_summary(st, wall),
    )
    (a.out / "run.json").write_text(json.dumps(info, indent=2, default=float) + "\n")
    gz_file(a.out / "recon.mic")
    (a.out / "run.done").write_text("done\n")
    print(a.arm, a.grid, "done", wall, "s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
