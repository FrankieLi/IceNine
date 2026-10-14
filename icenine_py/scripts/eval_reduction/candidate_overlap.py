"""Task C step 1: how much do the candidates of a level overlap? (measurement only, no library
change.) Runs the no-start search with the recorder hooks (``candidate_probe.Probe``) on

  full   the six interior full-sample voxels of the seed-cost diagnosis (clean Phase D images,
         Q-max 8 configuration, seed 1), lenient and C++ rule
  e0     the first 30 voxels of the E0 set (per-voxel clean images, seed 0, min_sin_eta 0), lenient
         and C++ rule

and stores the per-level candidate arrays; ``analyze`` computes the overlap statistics.

  uv run python scripts/eval_reduction/candidate_overlap.py run-full --voxel 15901 --rule cpp
  uv run python scripts/eval_reduction/candidate_overlap.py run-e0 --rule lenient --workers 10
  uv run python scripts/eval_reduction/candidate_overlap.py analyze --set full
  uv run python scripts/eval_reduction/candidate_overlap.py figure --voxel 15901 --rule cpp

Definitions (all misorientations symmetry-reduced; D = the level's diameter, the local-grid radius
and the spacing-filter radius of that level):
  (a) a quick-MC input is a duplicate when it is not the first of its group in a greedy scan in
      discrete-score order with tolerance t x D, t in {0.25, 0.5, 1}; ``has_better`` = some
      better-ranked candidate lies within t x D (an upper bound: that one may itself be dropped)
  (b) a quick-MC output is a duplicate by the same greedy scan in post-MC rank order, with
      tolerance 0.1 deg or 0.5 x D
  (c) the kept quarter = the first n_keep of the post-MC ranking; distinct basins = leaders of a
      greedy scan with tolerance 0.5 x D; the truth basin = a candidate within 3 deg (and 1 deg)
  (d) next-level screened points (kept x next local grid) that lie within 0.1 deg or
      0.25 x D_next of a point screened at this level (or of an earlier point of their own level)
  (e) cliques = the candidates of one FZ point (C++ grouping); two cliques overlap when a candidate
      of one lies within link x D of a candidate of the other; components = connected components
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICE / "scripts" / "findoptimal_robustness"))
sys.path.insert(0, str(ICE / "scripts" / "common"))
sys.path.insert(0, str(ICE / "scripts" / "phase_d"))

OUT = ICE / "benchmarks" / "eval_reduction" / "C" / "step1"
RAW = ICE / "benchmarks" / "eval_reduction" / "C" / "cache" / "step1_raw"  # gitignored
SIX = [15901, 17034, 18558, 19889, 22867, 5242]
N_E0 = 30
TOLS_IN = (0.25, 0.5, 1.0)
BASIN_DEG = (1.0, 3.0)
LINKS = (0.5, 1.0)


# ---------------------------------------------------------------------------------- running
def flatten(probe: Any, meta: Dict[str, Any]) -> Dict[str, Any]:
    out: Dict[str, Any] = dict(meta=json.dumps(meta, default=float))
    for L in probe.levels:
        for k, v in L.items():
            out[f"L{L['level']}_{k}"] = np.asarray(v)
    return out


def sym_quats(rec: Any) -> np.ndarray:
    from icenine.orientation_search import get_symmetry_quaternions

    return get_symmetry_quaternions(rec.setup.exp_setup.get_sample_symmetry())


def run_full(voxel: int, rule: str, seed: int, out_dir: Path) -> str:
    """One full-sample voxel (the Task A ledger setup)."""
    import eval_ledger as EL

    os.chdir(EL.EX)
    import torch

    torch.set_num_threads(1)
    from icenine import reconstructor as RC
    from icenine.config_file import ConfigFile
    from icenine.experimental_data import ExperimentalData
    from icenine.mic_file import MicFile
    from candidate_probe import Probe

    cfg = ConfigFile.from_file(str(EL.CONFIG))
    exp = ExperimentalData.from_binary_memmap(EL.STACK, 180, 2)
    setup = RC.setup_reconstruction(cfg, exp_data=exp)
    if rule == "cpp":
        setup.search_params.cpp_overlap_rule = True
    rec = RC.AdaptiveVoxelReconstructor(setup)
    v = setup.sample.get_mic().voxels[voxel]
    verts = RC._get_voxel_vertices(v)
    R_true = np.asarray(
        MicFile.read("SimInput/rand_500grains_1mm_neworient_s0.mic").voxels[voxel].orientation,
        dtype=np.float64,
    )
    probe = Probe(rec, keep_raw=True)
    t0 = time.perf_counter()
    res = probe.run(verts, v.phase, np.random.default_rng(seed))
    wall = time.perf_counter() - t0
    from csl import reduced_misorientation_deg

    err = float(reduced_misorientation_deg(np.asarray(res.orientation), R_true))
    meta = dict(
        set="full", voxel=voxel, rule=rule, seed=seed, err_deg=err, evals=probe.total_evals,
        wall_s=wall, label="contended (several processes at once)",
        R_true=R_true.tolist(), sym=sym_quats(rec).tolist(),
    )  # fmt: skip
    out_dir.mkdir(parents=True, exist_ok=True)
    fn = out_dir / f"full_{rule}_v{voxel}.npz"
    np.savez_compressed(fn, **flatten(probe, meta))
    return f"full {rule} v{voxel}: {wall:.0f}s evals={probe.total_evals} err={err:.4f} -> {fn.name}"


def e0_task(item: Any) -> str:
    """One E0 voxel (clean image), as scripts/phase_d/e0_cpprule.py."""
    import common as C
    from optimizer_sweep import voxel_context
    from candidate_probe import Probe

    vidx, vpos, rule, path = item
    W = C.get_worker()
    t0 = time.time()
    keys = np.load(C.CACHE_DIR / "images" / f"v{vidx}_clean.npz")["keys"]
    C.attach(keys)
    vctx = voxel_context(W.ctx, vidx)
    rec = W.rec
    Rt = np.asarray(vctx.R_true, dtype=np.float64)
    rec.min_sin_eta = 0.0
    sp = rec.setup.search_params
    sp.cpp_overlap_rule = rule == "cpp"
    sp.cpp_center_screen = True
    probe = Probe(rec, keep_raw=True)
    try:
        res = probe.run(vctx.vertices, vctx.voxel.phase, C.run_seed(vpos, 0))
    finally:
        sp.cpp_overlap_rule = False
        rec.min_sin_eta = 0.3
    err = float(C.err_deg(np.asarray(res.orientation, dtype=np.float64), Rt))
    meta = dict(
        set="e0", voxel=int(vidx), vpos=int(vpos), rule=rule, seed=0, err_deg=err,
        evals=probe.total_evals, R_true=Rt.tolist(), sym=sym_quats(rec).tolist(),
        label="contended multi-worker",
    )  # fmt: skip
    np.savez_compressed(path, **flatten(probe, meta))
    return f"e0 {rule} v{vidx}: {time.time() - t0:.0f}s evals={probe.total_evals} err={err:.3f}"


def run_e0(rule: str, workers: int, limit: int) -> None:
    import common as C

    info = dict(np.load(C.OUT_DIR / "voxels.npz"))
    vox = [int(v) for v in info["voxel_indices"]][: limit or N_E0]
    out = RAW
    out.mkdir(parents=True, exist_ok=True)
    items = [(v, p, rule, str(out / f"e0_{rule}_v{v}.npz")) for p, v in enumerate(vox)]
    items = [it for it in items if not Path(it[3]).exists()]
    print(f"{len(items)} tasks, {workers} workers", flush=True)
    C.run_pool(e0_task, items, workers, f"c_step1_{rule}")


# ---------------------------------------------------------------------------------- analysis
def load(path: Path) -> Dict[str, Any]:
    z = np.load(path, allow_pickle=False)
    meta = json.loads(str(z["meta"]))
    levels: Dict[int, Dict[str, Any]] = {}
    for k in z.files:
        if k.startswith("L") and "_" in k:
            lv, name = k[1:].split("_", 1)
            levels.setdefault(int(lv), {})[name] = z[k]
    return dict(meta=meta, levels=levels)


def run_metrics(run: Dict[str, Any]) -> List[Dict[str, Any]]:
    from icenine import candidate_dedup as CD
    from csl import reduced_misorientation_deg

    meta = run["meta"]
    sym = np.array(meta["sym"])
    Rt = np.array(meta["R_true"])
    rows: List[Dict[str, Any]] = []
    nlev = len(run["levels"])
    for l in range(nlev):
        L = run["levels"][l]
        D = float(L["diameter"])
        R = L["R_disc"]
        n = len(R)
        qe = L["quick_evals"].astype(np.int64)
        tot_q = int(qe.sum())
        row: Dict[str, Any] = dict(
            level=l, n_disc=n, n_mc=int(L["n_mc"]), n_keep=int(L["n_keep"]),
            diameter_deg=float(np.degrees(D)), quick_evals=tot_q,
            screen_recost_evals=int(L["screen_recost_evals"]), clique_ok=bool(L["clique_ok"]),
            n_cliques=int(L["n_cliques"]), n_fz=int(L["n_fz"]),
        )  # fmt: skip
        # (a) duplicate inputs, in discrete-score order (the list is sorted)
        for t in TOLS_IN:
            nb = CD.neighbour_lists(R, t * D, sym)
            kept, leader = CD.greedy_nms(R, t * D, sym, neighbours=nb)
            dup = np.ones(n, bool)
            dup[kept] = False
            row[f"in_dup_{t}"] = int(dup.sum())
            row[f"in_has_better_{t}"] = int(sum(1 for i, a in enumerate(nb) if len(a) and a[0] < i))
            row[f"in_dup_evals_{t}"] = int(qe[dup].sum())
        # (b) duplicate outputs, in post-MC rank order
        Ra = L["R_after_ranked"]
        order = np.sort(L["mc_idx"]) if False else None  # noqa: F841 (documentation only)
        ev_ranked = qe[L["mc_idx"]]
        for name, tol in (("0.1deg", np.radians(0.1)), ("0.5D", 0.5 * D)):
            nb = CD.neighbour_lists(Ra, tol, sym)
            kept, _ = CD.greedy_nms(Ra, tol, sym, neighbours=nb)
            dup = np.ones(len(Ra), bool)
            dup[kept] = False
            row[f"out_dup_{name}"] = int(dup.sum())
            row[f"out_any_{name}"] = int(sum(1 for a in nb if len(a)))
            row[f"out_dup_evals_{name}"] = int(ev_ranked[dup].sum())
        # (c) the kept quarter and the truth basin
        nk = int(L["n_keep"])
        Rk = Ra[:nk]
        kept, _ = CD.greedy_nms(Rk, 0.5 * D, sym)
        row["kept_distinct_0.5D"] = int(len(kept))
        kept1, _ = CD.greedy_nms(Rk, np.radians(0.1), sym)
        row["kept_distinct_0.1deg"] = int(len(kept1))
        for b in BASIN_DEG:
            d_disc = reduced_misorientation_deg(R, np.broadcast_to(Rt, R.shape).copy())
            d_after = reduced_misorientation_deg(Ra, np.broadcast_to(Rt, Ra.shape).copy())
            row[f"truth_in_disc_{b:g}"] = bool((d_disc < b).any())
            row[f"truth_in_after_{b:g}"] = bool((d_after < b).any())
            row[f"truth_in_kept_{b:g}"] = bool((d_after[:nk] < b).any())
        row["min_dist_kept_deg"] = float(d_after[:nk].min())
        # (e) cliques
        cl = L["clique"].astype(np.int64)
        _, inv, cnt = np.unique(cl, return_inverse=True, return_counts=True)
        row["clique_size_mean"] = float(cnt.mean())
        row["clique_size_max"] = int(cnt.max())
        row["clique_size_1"] = int((cnt == 1).sum())
        # second and later candidates of one clique (discrete rank order within the clique)
        seen: Dict[int, int] = {}
        later = np.zeros(n, bool)
        for i, c in enumerate(cl):
            later[i] = seen.get(int(c), 0) > 0
            seen[int(c)] = seen.get(int(c), 0) + 1
        row["later_in_clique"] = int(later.sum())
        row["later_in_clique_evals"] = int(qe[later].sum())
        # kept-quarter contribution per clique (kept discrete indices = first n_keep of the ranking)
        kept_disc = L["mc_idx"][:nk]
        kc = np.bincount(inv[kept_disc], minlength=len(cnt))
        row["cliques_with_kept"] = int((kc > 0).sum())
        row["kept_per_clique_max"] = int(kc.max()) if len(kc) else 0
        row["kept_in_cliques_with_2plus"] = int(kc[kc >= 2].sum())
        for link in LINKS:
            lab = CD.merged_clique_labels(R, cl, link * D, sym)
            cl_lab = np.zeros(len(cnt), np.int64)
            cl_lab[inv] = lab
            sizes = np.bincount(cl_lab)
            sizes = sizes[sizes > 0]
            row[f"link{link}_components"] = int(len(sizes))
            row[f"link{link}_cliques_overlapping"] = int(sizes[sizes >= 2].sum())
            row[f"link{link}_max_component"] = int(sizes.max())
            row[f"link{link}_median_component"] = float(np.median(sizes))
            cand_per_comp = np.bincount(lab)
            row[f"link{link}_max_component_candidates"] = int(cand_per_comp.max())
        # (d) next-level screened points against this level's
        if l + 1 < nlev and "fz" in run["levels"][l + 1] and "fz" in L:
            Ln = run["levels"][l + 1]
            A = (L["local_grid"][None] @ L["fz"][:, None]).reshape(-1, 3, 3)
            B = (Ln["local_grid"][None] @ Ln["fz"][:, None]).reshape(-1, 3, 3)
            AB = np.concatenate([A, B])
            Dn = float(Ln["diameter"])
            for name, tol in (("0.1deg", np.radians(0.1)), ("0.25Dn", 0.25 * Dn)):
                nb = CD.neighbour_lists(AB, tol, sym)
                nA = len(A)
                inA = sum(1 for a in nb[nA:] if len(a) and a[0] < nA)
                self_dup = sum(1 for i, a in enumerate(nb[nA:], start=nA) if len(a[(a >= nA)]) and
                               a[a >= nA][0] < i)  # fmt: skip
                row[f"next_in_prev_{name}"] = int(inA)
                row[f"next_self_dup_{name}"] = int(self_dup)
            row["next_points"] = int(len(B))
            row["this_points"] = int(len(A))
        rows.append(row)
    return rows


def analyze(which: str) -> None:
    import doc_tables as DT

    files = sorted(RAW.glob(f"{which}_*.npz"))
    print(f"{len(files)} run files", flush=True)
    per_run: List[Dict[str, Any]] = []
    for fn in files:
        run = load(fn)
        for r in run_metrics(run):
            r.update(
                rule=run["meta"]["rule"], voxel=run["meta"]["voxel"], err_deg=run["meta"]["err_deg"]
            )
            per_run.append(r)
        print(f"  {fn.name}", flush=True)
    (OUT / f"metrics_{which}.json").write_text(json.dumps(per_run, indent=1, default=float) + "\n")
    summarize(which, per_run)


def pct(a: float, b: float) -> str:
    return f"{100.0 * a / b:.1f}%" if b else "n/a"


def summarize(which: str, per_run: List[Dict[str, Any]]) -> None:
    import doc_tables as DT

    summ: Dict[str, Any] = {}
    tabs: Dict[str, str] = {}
    rules = sorted({r["rule"] for r in per_run})
    levels = sorted({r["level"] for r in per_run})
    for rule in rules:
        for scope, lv_sel in (("all", levels), ("0-2", [l for l in levels if l <= 2])):
            sel = [r for r in per_run if r["rule"] == rule and r["level"] in lv_sel]
            tot = lambda k: sum(r.get(k, 0) for r in sel)  # noqa: E731
            n = tot("n_disc")
            nmc = tot("n_mc")
            q = tot("quick_evals")
            summ[f"{rule}_levels_{scope}"] = dict(
                n_runs=len({(r["voxel"]) for r in sel}),
                n_disc=n, quick_evals=q,
                **{f"in_dup_{t}": tot(f"in_dup_{t}") for t in TOLS_IN},
                **{f"in_dup_evals_{t}": tot(f"in_dup_evals_{t}") for t in TOLS_IN},
                out_dup_0_1deg=tot("out_dup_0.1deg"), out_dup_0_5D=tot("out_dup_0.5D"),
                out_dup_evals_0_5D=tot("out_dup_evals_0.5D"),
                out_dup_evals_0_1deg=tot("out_dup_evals_0.1deg"),
                later_in_clique=tot("later_in_clique"),
                later_in_clique_evals=tot("later_in_clique_evals"),
            )  # fmt: skip
    for rule in rules:
        rows = []
        for l in levels:
            sel = [r for r in per_run if r["rule"] == rule and r["level"] == l]
            if not sel:
                continue
            tot = lambda k: sum(r.get(k, 0) for r in sel)  # noqa: E731
            n, q = tot("n_disc"), tot("quick_evals")
            rows.append(
                dict(
                    level=l,
                    runs=len(sel),
                    diameter_deg=float(np.mean([r["diameter_deg"] for r in sel])),
                    candidates=n,
                    inputs_dup_025=pct(tot("in_dup_0.25"), n),
                    inputs_dup_05=pct(tot("in_dup_0.5"), n),
                    inputs_dup_1=pct(tot("in_dup_1.0"), n),
                    outputs_dup_01deg=pct(tot("out_dup_0.1deg"), n),
                    outputs_dup_05D=pct(tot("out_dup_0.5D"), n),
                    evals_dup_inputs_05=pct(tot("in_dup_evals_0.5"), q),
                    evals_dup_outputs_05D=pct(tot("out_dup_evals_0.5D"), q),
                    evals_later_in_clique=pct(tot("later_in_clique_evals"), q),
                    quick_evals=q,
                )
            )
        tabs[f"c_step1_{which}_{rule}_dups"] = DT.markdown_table(
            rows,
            ["level", "runs", "diameter_deg", "candidates", "inputs_dup_025", "inputs_dup_05",
             "inputs_dup_1", "outputs_dup_01deg", "outputs_dup_05D", "evals_dup_inputs_05",
             "evals_dup_outputs_05D", "evals_later_in_clique", "quick_evals"],
            formats={"diameter_deg": ".2f", "quick_evals": ",d", "candidates": ",d"},
        )  # fmt: skip
        rows = []
        for l in levels:
            sel = [r for r in per_run if r["rule"] == rule and r["level"] == l]
            if not sel:
                continue
            tot = lambda k: sum(r.get(k, 0) for r in sel)  # noqa: E731
            rows.append(
                dict(
                    level=l,
                    runs=len(sel),
                    kept=tot("n_keep"),
                    kept_distinct_05D=tot("kept_distinct_0.5D"),
                    kept_distinct_01deg=tot("kept_distinct_0.1deg"),
                    truth3_in_candidates=sum(r["truth_in_disc_3"] for r in sel),
                    truth3_after_mc=sum(r["truth_in_after_3"] for r in sel),
                    truth3_in_kept=sum(r["truth_in_kept_3"] for r in sel),
                    truth1_after_mc=sum(r["truth_in_after_1"] for r in sel),
                    truth1_in_kept=sum(r["truth_in_kept_1"] for r in sel),
                )
            )
        tabs[f"c_step1_{which}_{rule}_kept"] = DT.markdown_table(
            rows,
            ["level", "runs", "kept", "kept_distinct_05D", "kept_distinct_01deg",
             "truth3_in_candidates", "truth3_after_mc", "truth3_in_kept", "truth1_after_mc",
             "truth1_in_kept"],
            formats={"kept": ",d"},
        )  # fmt: skip
        rows = []
        for l in levels:
            sel = [r for r in per_run if r["rule"] == rule and r["level"] == l]
            if not sel:
                continue
            tot = lambda k: sum(r.get(k, 0) for r in sel)  # noqa: E731
            nf = sum(r["n_cliques"] for r in sel)
            rows.append(
                dict(
                    level=l,
                    runs=len(sel),
                    cliques=nf,
                    candidates=tot("n_disc"),
                    cand_per_clique=tot("n_disc") / max(nf, 1),
                    cliques_with_kept=tot("cliques_with_kept"),
                    kept=tot("n_keep"),
                    kept_in_cliques_2plus=tot("kept_in_cliques_with_2plus"),
                    overlap_cliques_1D=pct(tot("link1.0_cliques_overlapping"), nf),
                    max_comp_cliques_1D=max(r["link1.0_max_component"] for r in sel),
                    overlap_cliques_05D=pct(tot("link0.5_cliques_overlapping"), nf),
                    max_comp_cliques_05D=max(r["link0.5_max_component"] for r in sel),
                    next_in_prev_01deg=pct(tot("next_in_prev_0.1deg"), tot("next_points")),
                    next_in_prev_025Dn=pct(tot("next_in_prev_0.25Dn"), tot("next_points")),
                    next_self_dup_025Dn=pct(tot("next_self_dup_0.25Dn"), tot("next_points")),
                )
            )
        tabs[f"c_step1_{which}_{rule}_cliques"] = DT.markdown_table(
            rows,
            ["level", "runs", "cliques", "candidates", "cand_per_clique", "cliques_with_kept",
             "kept", "kept_in_cliques_2plus", "overlap_cliques_1D", "max_comp_cliques_1D",
             "overlap_cliques_05D", "max_comp_cliques_05D", "next_in_prev_01deg",
             "next_in_prev_025Dn", "next_self_dup_025Dn"],
            formats={"cand_per_clique": ".2f", "cliques": ",d", "candidates": ",d", "kept": ",d"},
        )  # fmt: skip
    (OUT / f"summary_{which}.json").write_text(json.dumps(summ, indent=1, default=float) + "\n")
    DT.write_tables(OUT / f"tables_{which}.md", tabs)
    print(json.dumps(summ, indent=1, default=float)[:3000])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("cmd", choices=["run-full", "run-e0", "analyze"])
    ap.add_argument("--voxel", type=int, default=0)
    ap.add_argument("--rule", choices=["lenient", "cpp"], default="lenient")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--set", dest="which", default="full", choices=["full", "e0"])
    a = ap.parse_args()
    if a.cmd == "run-full":
        print(run_full(a.voxel, a.rule, a.seed, RAW), flush=True)
    elif a.cmd == "run-e0":
        run_e0(a.rule, a.workers, a.limit)
    else:
        analyze(a.which)


if __name__ == "__main__":
    main()
