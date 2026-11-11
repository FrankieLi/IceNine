"""AR3 step 10 fixes round (decision D-20261111-N9): test the reading of the in-run pass.

  libkeys_fresh.py --kind S8L --stage run     # BFS with the keys, the in-run pass, then the pass
                                              # again in the same process on the .mic-rounded map
  libkeys_fresh.py --kind S8L --stage fresh   # a fresh process: the pass on the unrounded map
  ... --smoke N                               # pilot: first N voxels, own folder

Check (1) of step 10 applied the library pass to the script's rounded map in a fresh process; the
in-run pass reads the unrounded map, after the BFS, in the same process. Two things differ: the
input rounding and the process. The step 10 cache holds no unrounded map (``records.json.gz`` has
no orientations, the ``.mic`` is rounded, ``pass_tries.npz`` holds float32 ends), so the BFS is run
again here and the map the pass reads is saved before the pass (float64, ``.npz``, lossless).

Stage ``run`` (one process): the library BFS with the keys of ``libkeys_common.KEYS``; a wrapper
around ``run_boundary_pass`` saves the map (orientation, reconstruction state, cost, confidence,
hit ratio of each voxel) and runs the in-run pass. After the BFS the map is put back as it was
before the pass, written as a ``.mic``, read back (the script's reading, ``G.read_map``) and the
pass is applied again in the same process to this rounded map. Stage ``fresh`` (a new process): the
same set-up and keys, the saved unrounded map put into the BFS's map, the pass applied once.

Compared try for try (per voxel: candidate, end as float64 where both sides have it, else
float32, cost, hit ratio, confidence, evaluations):
  (i)   the rerun's in-run pass against the stored step 10 in-run pass (``cache/<kind>``);
  (ii)  the before-pass ``.mic`` against the stored one (bytes), the map read back against the
        script's ``before_<map>.npy``;
  (iii) the same-process pass on the rounded map against the script's stored ``tries.npz``;
  (iv)  the fresh-process pass on the unrounded map against the in-run pass;
  and, for the size of the rounding's effect, the same-process pass on the rounded map against the
  in-run pass (unrounded).
Writes ``<CACHE>/fresh_<kind>/`` (large files, ``run.done`` / ``fresh.done``) and
``keys/fresh_check_<kind>.json``. No truth is read.
"""

import argparse
import json
import os
import threading
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

import libkeys_common as K
import numpy as np

import ar3_common as C  # noqa: E402
import gradient_common as G  # noqa: E402

CAP_MIN = 60  # twice the 30 minute budget of a whole-layer run
FIELDS = ("cand", "cost", "hit", "conf", "evals")


def watchdog(out: Path, cap_s: float, t0: float) -> None:
    def loop() -> None:
        time.sleep(cap_s)
        (out / "capped.json").write_text(
            json.dumps(dict(cap_s=cap_s, elapsed_s=time.time() - t0)) + "\n"
        )
        print(f"CAP of {cap_s} s reached; stopped", flush=True)
        os._exit(3)

    threading.Thread(target=loop, daemon=True).start()


def tries_arrays(res: Any) -> Dict[str, np.ndarray]:
    """The tries of a ``BoundaryPassResult`` as arrays (ends kept as float64)."""
    t = res.tries
    return dict(
        voxel=np.array([x["voxel"] for x in t], np.int64),
        cand=np.array([x["cand"] for x in t], np.int64),
        end=np.array([x["end"] for x in t], np.float64).reshape(-1, 3, 3),
        cost=np.array([x["cost"] for x in t], np.float64),
        hit=np.array([x["hit"] for x in t], np.float64),
        conf=np.array([x["conf"] for x in t], np.float64),
        evals=np.array([x["evals"] for x in t], np.int64),
        changed=np.array(res.changed, np.int64),
    )


def load_tries(path: Path) -> Dict[str, np.ndarray]:
    with np.load(path) as z:
        return {k: z[k] for k in z.files}


def compare_tries(a: Dict[str, np.ndarray], b: Dict[str, np.ndarray]) -> Dict[str, Any]:
    """Try-for-try comparison of two passes, by voxel. A voxel differs if its list of tries differs
    in length or in any field of any try. Ends are compared as float64 if both are float64, else
    after a cast to float32."""
    f32 = a["end"].dtype != np.float64 or b["end"].dtype != np.float64
    ea = a["end"].astype(np.float32) if f32 else a["end"]
    eb = b["end"].astype(np.float32) if f32 else b["end"]

    def by_voxel(x: Dict[str, np.ndarray]) -> Dict[int, np.ndarray]:
        out: Dict[int, List[int]] = {}
        for i, v in enumerate(x["voxel"].tolist()):
            out.setdefault(int(v), []).append(i)
        return {v: np.array(ix) for v, ix in out.items()}

    ia, ib = by_voxel(a), by_voxel(b)
    differ: List[int] = []
    differ_tries = 0
    for v in sorted(set(ia) | set(ib)):
        xa, xb = ia.get(v, np.array([], int)), ib.get(v, np.array([], int))
        if len(xa) != len(xb):
            differ.append(v)
            differ_tries += max(len(xa), len(xb))
            continue
        bad = np.zeros(len(xa), bool)
        for f in FIELDS:
            bad |= a[f][xa] != b[f][xb]
        bad |= ~(ea[xa] == eb[xb]).reshape(len(xa), -1).all(1)
        if bad.any():
            differ.append(v)
            differ_tries += int(bad.sum())
    ca = set(a["changed"].tolist()) if "changed" in a else None
    cb = set(b["changed"].tolist()) if "changed" in b else None
    return dict(
        n_tries=[int(len(a["voxel"])), int(len(b["voxel"]))],
        evals=[int(a["evals"].sum()), int(b["evals"].sum())],
        ends_compared_as="float32" if f32 else "float64",
        voxels_with_differing_tries=len(differ),
        differing_tries=differ_tries,
        first_differing_voxels=differ[:10],
        voxels_changed=[len(ca) if ca is not None else None, len(cb) if cb is not None else None],
        same_changed_voxels=(ca == cb) if ca is not None and cb is not None else None,
        equal=not differ,
    )


def save_map(mic: Any, path: Path) -> None:
    vox = mic.voxels
    np.savez(
        path,
        R=np.stack([np.asarray(v.orientation, np.float64) for v in vox]),
        state=np.array([int(v.reconstruction_id) for v in vox], np.int64),
        cost=np.array([v.cost for v in vox], np.float64),
        conf=np.array([v.confidence for v in vox], np.float64),
        hit=np.array([v.overlap_ratio for v in vox], np.float64),
    )


def put_map(mic: Any, z: Dict[str, np.ndarray], R: Optional[np.ndarray] = None) -> None:
    R = z["R"] if R is None else R
    for i, v in enumerate(mic.voxels):
        v.orientation = np.asarray(R[i], np.float64).copy()
        v.reconstruction_id = int(z["state"][i])
        v.cost, v.confidence, v.overlap_ratio = (
            float(z["cost"][i]),
            float(z["conf"][i]),
            float(z["hit"][i]),
        )


def stage_run(kind: str, out: Path, max_voxels: Optional[int]) -> None:
    from icenine import boundary_pass as LBP

    label = C.preflight_label(out)
    setup, bfs, _t = C.build(
        C.CONFIG, out, sample=G.GRID_MIC, stack_name=G.STACK, extra="\n".join(K.KEYS[kind])
    )
    mic = setup.sample.get_mic()
    assert setup.search_params.boundary_pass
    orig = LBP.run_boundary_pass
    got: Dict[str, Any] = {}

    def capture(bfs_: Any, mic_: Any) -> Any:
        save_map(mic_, out / "before_pass_unrounded.npz")
        t = time.time()
        res = orig(bfs_, mic_)
        got["inrun_wall_s"] = time.time() - t
        got["inrun"] = res
        return res

    LBP.run_boundary_pass = capture
    t0 = time.time()
    watchdog(out, CAP_MIN * 60.0, t0)
    bfs.reconstruct_sample(
        output_mic=str(out / "recon.mic"), max_voxels=max_voxels, rng=np.random.default_rng(0)
    )
    wall_bfs = time.time() - t0
    LBP.run_boundary_pass = orig
    np.savez(out / "inrun_tries.npz", **tries_arrays(got["inrun"]))
    # the map before the pass, written as a .mic and read back as the script read it
    z = load_tries(out / "before_pass_unrounded.npz")
    put_map(mic, z)
    mic.write(str(out / "before_pass.mic"))
    Rr = G.read_map(out / "before_pass.mic", len(mic.voxels))[0]
    np.save(out / "before_pass_rounded.npy", Rr)
    put_map(mic, z, R=Rr)
    t = time.time()
    res2 = orig(bfs, mic)
    wall2 = time.time() - t
    np.savez(out / "inproc_rounded_tries.npz", **tries_arrays(res2))
    info = dict(
        kind=kind,
        timing_label=label,
        wall_bfs_with_inrun_pass_s=wall_bfs,
        wall_inrun_pass_s=got["inrun_wall_s"],
        wall_inproc_rounded_pass_s=wall2,
        max_abs_rounding=float(np.abs(Rr - z["R"]).max()),
        voxels_rounded_differently=int((Rr != z["R"]).reshape(len(Rr), -1).any(1).sum()),
        # each refinement end of a seed step with the early stop gets one overlap evaluation (its
        # hit ratio), not counted in the totals; ``refined`` lists the refinements made
        seed_steps=len(bfs.index_log),
        seed_step_refinements=int(sum(len(e["refined"]) for e in bfs.index_log)),
        early_stop_hit_evals=(
            int(sum(len(e["refined"]) for e in bfs.index_log))
            if setup.search_params.seed_index_stop_hit > 0
            else 0
        ),
        n_evals_boundary_score=int(bfs.stats.get("n_evals_boundary_score", 0)),
    )
    (out / "run_info.json").write_text(json.dumps(info, indent=1) + "\n")
    (out / "run.done").write_text("done\n")
    print(json.dumps(info, indent=1), flush=True)


def stage_fresh(kind: str, out: Path, smoke: bool) -> None:
    from icenine import boundary_pass as LBP

    label = C.preflight_label(out)
    setup, bfs, _t = C.build(
        C.CONFIG, out, sample=G.GRID_MIC, stack_name=G.STACK, extra="\n".join(K.KEYS[kind])
    )
    mic = setup.sample.get_mic()
    z = load_tries(out / "before_pass_unrounded.npz")
    assert len(z["R"]) == len(mic.voxels)
    put_map(mic, z)
    t = time.time()
    res = LBP.run_boundary_pass(bfs, mic)
    wall = time.time() - t
    fresh = tries_arrays(res)
    np.savez(out / "fresh_unrounded_tries.npz", **fresh)
    inrun = load_tries(out / "inrun_tries.npz")
    inproc = load_tries(out / "inproc_rounded_tries.npz")
    run_info = json.loads((out / "run_info.json").read_text())
    info: Dict[str, Any] = dict(
        kind=kind,
        smoke=smoke,
        timing_label_run=run_info["timing_label"],
        timing_label_fresh=label,
        wall_bfs_with_inrun_pass_s=run_info["wall_bfs_with_inrun_pass_s"],
        wall_inrun_pass_s=run_info["wall_inrun_pass_s"],
        wall_inproc_rounded_pass_s=run_info["wall_inproc_rounded_pass_s"],
        wall_fresh_unrounded_pass_s=wall,
        max_abs_rounding=run_info["max_abs_rounding"],
        voxels_rounded_differently=run_info["voxels_rounded_differently"],
        seed_steps=run_info.get("seed_steps"),
        seed_step_refinements=run_info.get("seed_step_refinements"),
        early_stop_hit_evals=run_info.get("early_stop_hit_evals"),
        iv_fresh_unrounded_vs_inrun=compare_tries(fresh, inrun),
        # the size of the rounding's effect: the same process, rounded against unrounded input
        rounding_inproc_rounded_vs_inrun=compare_tries(inproc, inrun),
    )
    if not smoke:
        stored_inrun = load_tries(K.CACHE / kind / "pass_tries.npz")
        info["i_inrun_vs_stored_step10_inrun"] = compare_tries(inrun, stored_inrun)
        info["iv_fresh_unrounded_vs_stored_step10_inrun"] = compare_tries(fresh, stored_inrun)
        pdir, before_npy, _after = K.STORED_PASS[kind]
        script = load_tries(pdir / "tries.npz")
        info["iii_inproc_rounded_vs_script"] = compare_tries(inproc, script)
        info["inrun_vs_script"] = compare_tries(inrun, script)
        mic_bytes = (out / "before_pass.mic").read_bytes()
        info["ii_before_mic_equals_stored_step10_bytes"] = mic_bytes == K.gz_bytes(
            K.CACHE / kind / "recon_before_pass.mic.gz"
        )
        Rr = np.load(out / "before_pass_rounded.npy")
        info["ii_read_back_equals_script_before_npy"] = bool(
            np.array_equal(Rr, np.load(before_npy).astype(np.float64))
        )
    dest = (out / "fresh_check.json") if smoke else (K.OUT / f"fresh_check_{kind}.json")
    dest.write_text(json.dumps(info, indent=1) + "\n")
    (out / "fresh.done").write_text("done\n")
    print(json.dumps(info, indent=1), flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kind", required=True, choices=sorted(K.KEYS))
    ap.add_argument("--stage", required=True, choices=("run", "fresh"))
    ap.add_argument("--smoke", type=int, default=0, help="pilot: first N voxels, own folder")
    a = ap.parse_args()
    out = K.CACHE / (f"fresh_smoke_{a.kind}" if a.smoke else f"fresh_{a.kind}")
    out.mkdir(parents=True, exist_ok=True)
    if a.stage == "run":
        stage_run(a.kind, out, a.smoke or None)
    else:
        stage_fresh(a.kind, out, bool(a.smoke))


if __name__ == "__main__":
    main()
