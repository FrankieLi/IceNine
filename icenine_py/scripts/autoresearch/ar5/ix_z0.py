"""AR5 step 26, gate z0: is indexing a large enough share of a quiet ``J5SD`` run on ``mid``?

  RESEED_SAMPLE=mid OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 IX_PRE_FREEZE=1 uv run python \
      scripts/autoresearch/ar5/ix_z0.py [--n-evals 2000]

Gate z0 (frozen text): the old indexer on the 200 ``mid`` voxels of ``ix_profile.py``, scaled to the
1,626 stored steps, is at least 30% of a quiet ``J5SD`` estimate on ``mid``. The estimate is the
stored evaluation counts by stage (``J5D``'s ``run.json``, ``polish.json`` and the screen pass and
repair of ``stored_J5SD``) at this session's ms per evaluation: the strict cost of one orientation
at a ``mid`` voxel is timed here (``--n-evals`` evaluations, one thread), and the stage's ratio to
the index stage's ms per evaluation in the stored (contended) run is applied to it (the stages
call the same cost). The 200 voxels are timed here, old indexer only. The result carries the
preflight; with a refused preflight the numbers are "contended" and the gate is "not judged:
machine not quiet". Writes ``ar5/ix/z0.json``. No truth is read.
"""

import argparse
import json
import sys
import time
from typing import Any, Dict

import ix_common as IC
import ix_profile as IP

import numpy as np  # noqa: E402

sys.path.insert(0, str(IC.ICE / "scripts" / "common"))
import preflight as PF  # noqa: E402


def stage_table(cache: Any, s: str) -> Dict[str, Any]:
    r = IC.load_json(cache / f"J5D_{s}" / "run.json")
    pol = IC.load_json(cache / f"J5D_{s}" / "polish.json")
    rep = IC.load_json(cache / f"stored_J5SD_{s}" / "repair.json")
    st = r["stats"]
    idx = st["index_evals"]
    rescue = pol["rescue_evals"] + pol["polish_evals"]
    grid = pol["grid_evals"]
    nbr = st["n_evals_neighbor"] + st["n_evals_revisit"]
    ev = dict(
        index=idx,
        rescue_grid=rescue + grid,
        neighbour=nbr,
        screen=rep["pass_evals"],
        repair=rep["n_evals"],
    )
    ev["other_seed"] = (
        st["n_evals_seed"] + st["n_evals_reseed"] - idx - rescue - grid
    )  # the seed steps' other evaluations (counted with the index stage)
    wall = dict(
        index=st["index_path_wall_s"] - st["index_s"],
        rescue_grid=st.get("index_fallback_grid_wall_s", 0.0),
        neighbour=st["wall_neighbor_s"] + st["wall_revisit_s"],
        screen=rep["wall_pass_s"],
        repair=rep["wall_repair_s"],
    )
    ms = {k: 1e3 * wall[k] / max(1, ev[k]) for k in wall}
    return dict(
        evals=ev,
        stored_ms=ms,
        stored_wall=r["wall_total_s"] + rep["wall_pass_s"] + rep["wall_repair_s"],
        stored_index_s=st["index_s"],
        n_steps=len(IC.stored_steps(s)),
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-evals", type=int, default=2000)
    a = ap.parse_args()
    s = IC.SAMPLE
    if s != "mid":
        raise SystemExit("ix_z0.py: gate z0 is judged on mid")
    start = PF.preflight()
    quiet_start = True
    try:
        PF.require_quiet(info=start)
    except PF.MachineBusyError:
        quiet_start = False
    S = IC.setup()
    env, rec = S["env"], S["rec"]
    from icenine import ray_indexing as RI

    steps = IC.stored_steps(s)
    vox = IP.pick_voxels(steps, 200)
    # (1) the old indexer on the 200 voxels
    t0 = time.perf_counter()
    for v in vox:
        verts = env.hooks.vertices(v)
        ph = int(env.mic.voxels[v].phase)
        sc, crystal, params, _ = rec._index_context(ph)
        RI.index_voxel(sc, verts, crystal, params, RI.SEED_TOP_K)
    t_200 = time.perf_counter() - t0
    # (2) ms per strict evaluation at mid voxels
    rng = np.random.default_rng(0)
    n = 0
    t_eval = 0.0
    voxel_pool = vox[:20]
    per = max(1, a.n_evals // len(voxel_pool))
    for v in voxel_pool:
        verts = env.hooks.vertices(v)
        ph = int(env.mic.voxels[v].phase)
        sc, crystal, params, strict_cf = rec._index_context(ph)
        from scipy.spatial.transform import Rotation

        Rs = Rotation.random(per, random_state=int(rng.integers(1 << 30))).as_matrix()
        strict_cf.evaluate(np.asarray(Rs[0], np.float64), verts, ph)  # warm up
        t1 = time.perf_counter()
        for R in Rs:
            strict_cf.evaluate(np.asarray(R, np.float64), verts, ph)
        t_eval += time.perf_counter() - t1
        n += per
    ms_q = 1e3 * t_eval / n
    tab = stage_table(IC.CACHE, s)
    base = tab["stored_ms"]["index"]
    mult = {k: v / base for k, v in tab["stored_ms"].items()}
    ev = tab["evals"]
    est_eval_s = (
        1e-3
        * ms_q
        * (
            (ev["index"] + ev["other_seed"]) * mult["index"]
            + ev["rescue_grid"] * mult["rescue_grid"]
            + ev["neighbour"] * mult["neighbour"]
            + ev["screen"] * mult["screen"]
            + ev["repair"] * mult["repair"]
        )
    )
    idx_scaled = t_200 * tab["n_steps"] / len(vox)
    end = PF.preflight()
    quiet_end = True
    try:
        PF.require_quiet(info=end)
    except PF.MachineBusyError:
        quiet_end = False
    quiet = quiet_start and quiet_end
    share = idx_scaled / (idx_scaled + est_eval_s)
    out = dict(
        sample=s,
        n_voxels=len(vox),
        old_200_s=t_200,
        n_steps=tab["n_steps"],
        indexing_scaled_s=idx_scaled,
        ms_per_strict_eval_this_session=ms_q,
        stored_stage_ms_contended=tab["stored_ms"],
        stage_multipliers=mult,
        evals_by_stage=ev,
        estimated_evaluation_s=est_eval_s,
        estimate_j5sd_s=idx_scaled + est_eval_s,
        share=share,
        stored_wall_contended_s=tab["stored_wall"],
        stored_index_s_contended=tab["stored_index_s"],
        gate="share >= 0.30",
        contended_share_at_or_above_line=bool(share >= 0.30),
        timing_label="single-worker" if quiet else "contended",
        verdict=(
            ("met" if share >= 0.30 else "not met") if quiet else "not judged: machine not quiet"
        ),
        preflight_start=json.loads(
            json.dumps(start).replace(str(__import__("pathlib").Path.home()), "~")
        ),
        preflight_end=json.loads(
            json.dumps(end).replace(str(__import__("pathlib").Path.home()), "~")
        ),
    )
    (IC.RD / "z0.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: v for k, v in out.items() if not k.startswith("preflight")}, indent=1))


if __name__ == "__main__":
    main()
