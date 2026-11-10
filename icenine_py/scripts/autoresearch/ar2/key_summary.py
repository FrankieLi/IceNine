"""AR2 step 2: the 2,000-voxel BFS pilot with ``SeedByIndexing 1`` against the stored key-off rerun.

  uv run python scripts/autoresearch/ar2/key_summary.py baseline   # reproduce the stored numbers
  uv run python scripts/autoresearch/ar2/key_summary.py all        # summary.json, tables.md, figure

Baseline (from the stored files, before anything else): wrong voxels, evaluations and full searches
of ``benchmarks/eval_reduction/G/g2/bfs/runs/keyoff_rerun``. The stop rule of the plan: if they
differ from 37 wrong, 25,547,935 evaluations and 69 full searches, stop.

Run folders read from ``benchmarks/autoresearch/ar2/key/runs/<name>`` (copied from the gitignored
cache by the run script): ``key_on`` (the accuracy run), ``key_on_fused`` (the wall-time run).
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0] / "ar1"))
sys.path.insert(0, str(HERE))
import ar1_common as AC  # noqa: E402
import g2_bfs  # noqa: E402
import misses as M  # noqa: E402
from holdout_summary import neighbours  # noqa: E402

from icenine import ray_indexing as RI  # noqa: E402

import doc_tables as DT  # noqa: E402
import stats as ST  # noqa: E402

OUT = AC.ICE / "benchmarks" / "autoresearch" / "ar2" / "key"
RUNS = OUT / "runs"
OFF = AC.PILOT_RUN
WRONG_DEG = 1.0
STORED = dict(wrong=37, evals=25_547_935, full_searches=69)
EVAL_KEYS = ("seed", "reseed", "neighbor", "revisit")


def load_run(d: Path) -> Dict[str, Any]:
    run = json.loads((d / "run.json").read_text())
    log_p = d / "index_log.json"
    return dict(
        dir=d,
        run=run,
        stats=run["stats"],
        records=json.loads((d / "records.json").read_text()),
        index_log=json.loads(log_p.read_text()) if log_p.exists() else [],
    )


def evals_of(st: Dict[str, Any]) -> Dict[str, int]:
    e = {k: int(st[f"n_evals_{k}"]) for k in EVAL_KEYS}
    e["full"] = e["seed"] + e["reseed"]
    e["local"] = e["neighbor"] + e["revisit"]
    e["total"] = e["full"] + e["local"]
    return e


def baseline() -> Dict[str, Any]:
    """Wrong voxels, evaluations and full searches of the stored key-off rerun, from the files."""
    E = g2_bfs.errors(OFF)
    run = json.loads((OFF / "run.json").read_text())
    st = run["stats"]
    stored = AC.pilot_stored()
    out = dict(
        wrong=int((E["err"] > WRONG_DEG).sum()),
        evals=evals_of(st),
        seeds=int(st["n_seeds"]),
        reseeds=int(st["n_reseeds"]),
        full_searches=int(st["n_seeds"] + st["n_reseeds"]),
        full_search_voxels_by_evals=len(stored),
        wall_total_s=float(run["wall_total_s"]),
        timing_label=run.get("timing_label"),
        torch_threads=run.get("torch_threads"),
    )
    ok = (
        out["wrong"] == STORED["wrong"]
        and out["evals"]["total"] == STORED["evals"]
        and out["full_searches"] == STORED["full_searches"] == out["full_search_voxels_by_evals"]
    )
    out["reproduced"] = bool(ok)
    return out


# --------------------------------------------------------------------------------------------
# geometry of the layer
# --------------------------------------------------------------------------------------------
class Layer:
    """True grains of the full layer: grain map, one orientation per grain, neighbour sets."""

    def __init__(self) -> None:
        AC.init("full")
        self.ctx = AC.bench().RUN._ENV["ctx"]
        self.gmap = np.load(M.GRAINMAP)
        first: Dict[int, int] = {}
        for v, g in enumerate(self.gmap):
            first.setdefault(int(g), v)
        self.gids = sorted(first)
        self.Rg = np.stack([self.ctx.truth(first[g]) for g in self.gids])
        e0 = AC.env("full", 0, need_sc=False)
        self.sym = RI.proper_rotations(AC.symmetry_of(e0))
        self.vxy, self.side = M.voxel_vertices_xy(M.MIC)
        self.corner, self.edge = neighbours(self.vxy, self.side)
        self.region = np.load(AC.ICE / "benchmarks" / "phase_d_pilot" / "region_index.npy")
        self._far: Dict[Any, bool] = {}
        self._pos = {g: i for i, g in enumerate(self.gids)}

    def far(self, a: int, b: int) -> bool:
        if a == b:
            return False
        k = (min(a, b), max(a, b))
        if k not in self._far:
            ra, rb = self.Rg[self._pos[a]], self.Rg[self._pos[b]]
            self._far[k] = bool(RI.misorientation_deg(ra[None], rb[None], self.sym)[0] > WRONG_DEG)
        return self._far[k]

    def kind(self, v: int) -> str:
        """ "edge" (a voxel of another grain shares an edge), "corner" (another grain only at a
        corner) or "interior"."""
        g = int(self.gmap[v])
        if any(self.far(g, int(self.gmap[u])) for u in self.edge[v]):
            return "edge"
        if any(self.far(g, int(self.gmap[u])) for u in self.corner[v]):
            return "corner"
        return "interior"

    def grain_of(self, R: np.ndarray) -> Optional[int]:
        d = RI.misorientation_deg(R[None], self.Rg, self.sym)
        j = int(np.argmin(d))
        return self.gids[j] if d[j] < WRONG_DEG else None

    def touches(self, v: int, g: int) -> bool:
        return bool(M.touch(v, g, self.gmap, self.vxy, self.side)["touches"])


# --------------------------------------------------------------------------------------------
# paired comparison and the trace of every break
# --------------------------------------------------------------------------------------------
def paired_trace(on: Dict[str, Any], L: Layer) -> Dict[str, Any]:
    A, B = g2_bfs.errors(OFF), g2_bfs.errors(on["dir"])
    ea, eb = A["err"], B["err"]
    wa, wb = ea > WRONG_DEG, eb > WRONG_DEG
    pr = g2_bfs.paired(OFF, on["dir"])
    idx = L.region
    log = on["index_log"]
    used = {int(e["voxel"]) for e in log if e["used"]}
    fallback = {int(e["voxel"]) for e in log if not e["used"]}
    # full-search voxels of the on-run: indexed seeds and reseeds, and the grid fallbacks
    fs = sorted(used | fallback)
    qb = B["q_rec"]
    from score_bfs import pair_misorientation_deg

    def rotmat(i: int) -> np.ndarray:
        return recon_R[i]

    from icenine.mic_file import MicFile

    recon_R = np.array(
        [v.orientation for v in MicFile.read(str(on["dir"] / "recon.mic")).voxels], dtype=np.float64
    )

    def describe(r: int) -> Dict[str, Any]:
        f = int(idx[B["row"][r]])
        R = rotmat(r)
        g_taken = L.grain_of(R)
        d = dict(
            voxel=int(r),
            full_sample_voxel=f,
            grain=int(B["grain"][r]),
            err_on=float(eb[r]),
            err_off=float(ea[r]),
            source_on=str(B["source"][r]),
            source_off=str(A["source"][r]),
            kind=L.kind(f),
            taken_grain=g_taken,
            taken_touches=(None if g_taken is None else L.touches(f, g_taken)),
        )
        # origin of the orientation in the on-run
        if r in used:
            d["origin"] = "own indexed seed" if wb[r] else "own indexed seed (right)"
            d["carriers"] = [int(r)]
        elif r in fallback and str(B["source"][r]) == "seed":
            d["origin"] = "own grid fallback"
            d["carriers"] = [int(r)]
        else:
            dfs = pair_misorientation_deg(np.repeat(qb[r][None], len(fs), axis=0), qb[fs], B["sym"])
            car = [int(fs[i]) for i in np.nonzero(dfs <= WRONG_DEG)[0]]
            d["carriers"] = car
            wrong_car = [c for c in car if wb[c]]
            ind_wrong = [c for c in wrong_car if c in used]
            grid_wrong = [c for c in wrong_car if c not in used]
            if ind_wrong:
                d["origin"] = "carried from a wrong indexed seed"
            elif grid_wrong:
                d["origin"] = "carried from a wrong grid fallback"
            elif car:
                d["origin"] = "carried from a right full search (crossing)"
            else:
                d["origin"] = "origin not found"
            d["carriers_wrong_indexed"] = ind_wrong
        d["against_path"] = d["origin"] in ("own indexed seed", "carried from a wrong indexed seed")
        return d

    broken = [describe(r) for r in pr["broken_voxels"]]
    fixed_rows = []
    for r in pr["fixed_voxels"]:
        f = int(idx[A["row"][r]])
        fixed_rows.append(dict(voxel=int(r), full_sample_voxel=f, kind=L.kind(f)))
    both = [int(r) for r in np.nonzero(wa & wb)[0]]
    both_rows = [describe(r) for r in both]
    # every used indexed seed or reseed: right or wrong
    used_err = {int(r): float(eb[r]) for r in sorted(used)}
    n_wrong_used = int(sum(e > WRONG_DEG for e in used_err.values()))
    # kinds of all wrong voxels of each run
    kinds = {}
    for name, w, E in (("off", wa, A), ("on", wb, B)):
        ks = [L.kind(int(idx[E["row"][r]])) for r in np.nonzero(w)[0]]
        kinds[name] = {k: ks.count(k) for k in ("edge", "corner", "interior")}
    return dict(
        wrong_off=int(wa.sum()),
        wrong_on=int(wb.sum()),
        paired=dict(fixed=pr["fixed"], broken=pr["broken"], both_wrong=pr["both_wrong"]),
        broken=broken,
        fixed=fixed_rows,
        both_wrong=both_rows,
        used_indexed=len(used),
        used_indexed_wrong=n_wrong_used,
        fallback_voxels=sorted(fallback),
        kinds_wrong=kinds,
        median_err_right_on=float(np.median(eb[~wb])),
        median_err_right_off=float(np.median(ea[~wa])),
        unresolved_on=int((B["source"] == "unresolved").sum()),
        unresolved_off=int((A["source"] == "unresolved").sum()),
    )


def trace_candidates(
    rows: List[Dict[str, Any]], on: Dict[str, Any], L: Layer
) -> List[Dict[str, Any]]:
    """For each wrong indexed seed that is the origin of a break (or a break itself): the candidate
    taken, the strict costs before and after the finisher for it and for the truth candidate, and
    the strict cost at the exact truth. The indexing is run again (deterministic) and must give
    the logged strict costs."""
    from icenine.orientation_search import SearchCandidate

    log = {int(e["voxel"]): e for e in on["index_log"] if e["used"]}
    need: Set[int] = set()
    for b in rows:
        if b["against_path"]:
            need.update(b["carriers_wrong_indexed"] if "carriers_wrong_indexed" in b else [])
            if b["origin"] == "own indexed seed":
                need.add(b["voxel"])
    out = []
    p = RI.IndexParams()
    for r in sorted(need):
        f = int(L.region[r])
        e = AC.env("full", f)
        sc, rec, lf, verts, phase, Rt = (
            e["sc"],
            e["rec"],
            e["lf"],
            e["verts"],
            e["phase"],
            e["R_true"],
        )
        crystal = RI.crystal_from_smooth_cost(sc, AC.symmetry_of(e), phase, p.q_max)
        res, att = RI.index_voxel(sc, verts, crystal, p, 20)
        ent = log[r]
        costs = np.array([float(lf.evaluate(res.R[j], verts, phase).cost) for j in range(5)])
        assert np.allclose(costs, ent["strict_costs"], atol=1e-12), (r, costs, ent["strict_costs"])
        d = RI.misorientation_deg(res.R, Rt[None], L.sym)
        hit = np.nonzero(d < 2.0)[0]
        jt = int(hit[0]) if len(hit) else None
        row: Dict[str, Any] = dict(
            voxel=int(r),
            full_sample_voxel=f,
            n_cand=int(len(res.R)),
            truth_rank=None if jt is None else jt + 1,
            truth_in_first5=bool(jt is not None and jt < 5),
            taken=int(ent["taken"]),
            taken_rank=int(ent["taken"]) + 1,
            refined=ent["refined"],
            cost_at_truth=float(lf.evaluate(Rt, verts, phase).cost),
        )
        k = ent["refined"].index(ent["taken"])
        row["taken_cost"] = [
            float(ent["strict_costs"][ent["taken"]]),
            float(ent["finisher_costs"][k]),
        ]
        if jt is not None and jt < 5:
            row["truth_start_deg"] = float(d[jt])
            if jt in ent["refined"]:
                kt = ent["refined"].index(jt)
                row["truth_cost"] = [
                    float(ent["strict_costs"][jt]),
                    float(ent["finisher_costs"][kt]),
                ]
                row["truth_refined_in_run"] = True
            else:
                rr = rec.refine_from_candidates(
                    [SearchCandidate(orientation=res.R[jt].copy(), cost=float(costs[jt]))],
                    verts, phase, rng=np.random.default_rng(0),
                )  # fmt: skip
                row["truth_cost"] = [float(costs[jt]), float(rr.cost)]
                row["truth_refined_in_run"] = False
                row["truth_cost_run_for_the_trace_only"] = True
        out.append(row)
    return out


# --------------------------------------------------------------------------------------------
# tables, figure, summary
# --------------------------------------------------------------------------------------------
def wil(k: int, n: int) -> str:
    lo, hi = ST.wilson(k, n)
    return f"{k} of {n} ({100 * lo:.1f} to {100 * hi:.1f}%)"


def run_block(name: str, d: Path) -> Dict[str, Any]:
    r = load_run(d)
    st, run = r["stats"], r["run"]
    out: Dict[str, Any] = dict(
        name=name,
        evals=evals_of(st),
        seeds=int(st["n_seeds"]),
        reseeds=int(st["n_reseeds"]),
        seed_rejected=int(st["n_seed_rejected"]),
        reseed_rejected=int(st["n_reseed_rejected"]),
        neighbor_fits=int(st["n_neighbor_fits"]),
        revisit_attempted=int(st["n_revisit_attempted"]),
        wall_total_s=float(run["wall_total_s"]),
        timing_label=run.get("timing_label"),
        torch_threads=run.get("torch_threads"),
        fused_cost=run.get("fused_cost"),
        seed_by_indexing=run.get("seed_by_indexing"),
    )
    out["full_searches"] = out["seeds"] + out["reseeds"]
    E = g2_bfs.errors(d)
    out["wrong"] = int((E["err"] > WRONG_DEG).sum())
    if r["index_log"]:
        n = len(r["index_log"])
        ix = [e["index_s"] for e in r["index_log"]]
        path = [e["path_s"] for e in r["index_log"]]
        out["index"] = dict(
            runs=n,
            used=int(st.get("index_used", 0)),
            fallback_no_candidate=int(st.get("index_fallback_no_candidate", 0)),
            fallback_seed_test=int(st.get("index_fallback_seed_test", 0)),
            evals_used=int(st.get("index_evals_used", 0)),
            evals_wasted=int(st.get("index_evals_wasted", 0)),
            evals_fallback_grid=int(st.get("index_fallback_grid_evals", 0)),
            pairing_s=float(st["index_pairing_s"]),
            index_s_total=float(sum(ix)),
            index_s_per_voxel_mean=float(np.mean(ix)),
            index_s_per_voxel_median=float(np.median(ix)),
            index_s_per_voxel_max=float(np.max(ix)),
            path_s_total=float(sum(path)),
            path_s_per_voxel_mean=float(np.mean(path)),
            evals_per_voxel_mean=float(np.mean([e["evals"] for e in r["index_log"]])),
            wall_seed_s=float(st["wall_seed_s"]),
            wall_reseed_s=float(st["wall_reseed_s"]),
            wall_neighbor_s=float(st["wall_neighbor_s"]),
            wall_revisit_s=float(st["wall_revisit_s"]),
        )
    return out


def main() -> None:
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    base = baseline()
    print(json.dumps(base, indent=1), flush=True)
    if not base["reproduced"]:
        raise SystemExit("STOP: the stored baseline does not reproduce")
    if what == "baseline":
        OUT.mkdir(parents=True, exist_ok=True)
        (OUT / "baseline.json").write_text(json.dumps(base, indent=1) + "\n")
        return
    on_dir = RUNS / "key_on"
    on = load_run(on_dir)
    L = Layer()
    pt = paired_trace(on, L)
    cand = trace_candidates(pt["broken"] + pt["both_wrong"], on, L)
    summ: Dict[str, Any] = dict(
        baseline=base,
        off=run_block("key off (stored rerun)", OFF),
        on=run_block("key on", on_dir),
        pair=pt,
        candidate_trace=cand,
    )
    fused = RUNS / "key_on_fused"
    if (fused / "run.json").exists():
        summ["on_fused"] = run_block("key on, FusedCost on", fused)
        pf = g2_bfs.paired(OFF, fused)
        summ["on_fused"]["paired_vs_off"] = {k: pf[k] for k in ("fixed", "broken", "both_wrong")}
        pf2 = g2_bfs.paired(on_dir, fused)
        summ["on_fused"]["paired_vs_on"] = {k: pf2[k] for k in ("fixed", "broken", "both_wrong")}
        summ["on_fused"]["recon_identical_to_on"] = (fused / "recon.mic").read_bytes() == (
            on_dir / "recon.mic"
        ).read_bytes()
    off_fused = RUNS / "key_off_fused"
    if (off_fused / "run.json").exists():
        summ["off_fused"] = run_block("key off, FusedCost on", off_fused)
    (OUT / "summary.json").write_text(json.dumps(summ, indent=1, default=float) + "\n")
    tables(summ)
    figure(summ)


def tables(S: Dict[str, Any]) -> None:
    off, on, pt = S["off"], S["on"], S["pair"]
    n = 2000
    # accuracy
    rows = [
        dict(run="key off (stored rerun)", wrong=wil(off["wrong"], n), fixed="", broken="",
             both_wrong="", unresolved=pt["unresolved_off"]),
        dict(run="key on", wrong=wil(on["wrong"], n), fixed=pt["paired"]["fixed"],
             broken=pt["paired"]["broken"], both_wrong=pt["paired"]["both_wrong"],
             unresolved=pt["unresolved_on"]),
    ]  # fmt: skip
    t_acc = DT.markdown_table(rows, ["run", "wrong", "fixed", "broken", "both_wrong", "unresolved"])
    # evaluations
    ev_rows = []
    for name, r in (("key off (stored rerun)", off), ("key on", on)):
        e = r["evals"]
        n_fits = r["neighbor_fits"] + r["revisit_attempted"]
        ev_rows.append(
            dict(run=name, total=e["total"], full=e["full"], local=e["local"],
                 full_searches=r["full_searches"], local_fits=n_fits)
        )  # fmt: skip
    cols_ev = ["run", "total", "full", "local", "full_searches", "local_fits"]
    t_ev = DT.markdown_table(ev_rows, cols_ev, formats=dict(total=",d", full=",d", local=",d"))
    ix = on["index"]
    fb = ix["fallback_no_candidate"] + ix["fallback_seed_test"]
    t_ix_rows = [
        dict(item="seeds and reseeds", value=on["full_searches"]),
        dict(item="indexed (every seed and reseed)", value=ix["runs"]),
        dict(item="used (passed the seed test)", value=ix["used"]),
        dict(item="fallbacks to the grid search", value=fb),
        dict(item="  no candidate", value=ix["fallback_no_candidate"]),
        dict(item="  failed the seed test", value=ix["fallback_seed_test"]),
        dict(item="evaluations of used indexing paths", value=f"{ix['evals_used']:,d}"),
        dict(item="evaluations of indexing paths that fell back", value=f"{ix['evals_wasted']:,d}"),
        dict(item="evaluations of the grid fallbacks", value=f"{ix['evals_fallback_grid']:,d}"),
    ]  # fmt: skip
    t_ix = DT.markdown_table(t_ix_rows, ["item", "value"])
    # times
    lab = on["timing_label"]
    tm_rows = [
        dict(item="one-time pairing of the rays (s)", value=f"{ix['pairing_s']:.1f}"),
        dict(item="indexing, total (s)", value=f"{ix['index_s_total']:.1f}"),
        dict(item="indexing per voxel, mean (s)", value=f"{ix['index_s_per_voxel_mean']:.3f}"),
        dict(item="indexing per voxel, median (s)", value=f"{ix['index_s_per_voxel_median']:.3f}"),
        dict(item="indexing per voxel, max (s)", value=f"{ix['index_s_per_voxel_max']:.3f}"),
        dict(item="whole path per voxel, mean (s)", value=f"{ix['path_s_per_voxel_mean']:.2f}"),
    ]  # fmt: skip
    t_tm = DT.markdown_table(tm_rows, ["item", "value"])
    wall_rows = [
        dict(run="key off, stored rerun (FusedCost off)", wall_s=f"{off['wall_total_s']:.0f}",
             conditions="contended: 3 BFS runs at once, 1 thread each, 2026-10-29"),
        dict(run="key on (FusedCost off)", wall_s=f"{on['wall_total_s']:.0f}",
             conditions=lab),
    ]  # fmt: skip
    if "on_fused" in S:
        f = S["on_fused"]
        wall_rows.append(
            dict(run="key on (FusedCost on)", wall_s=f"{f['wall_total_s']:.0f}",
                 conditions=f["timing_label"])
        )  # fmt: skip
    if "off_fused" in S:
        f = S["off_fused"]
        wall_rows.append(
            dict(run="key off (FusedCost on)", wall_s=f"{f['wall_total_s']:.0f}",
                 conditions=f["timing_label"])
        )  # fmt: skip
    t_wall = DT.markdown_table(wall_rows, ["run", "wall_s", "conditions"])
    # breaks
    br_rows = []
    for b in pt["broken"]:
        br_rows.append(
            dict(voxel=b["full_sample_voxel"], kind=b["kind"], source_on=b["source_on"],
                 err_on=f"{b['err_on']:.2f}",
                 taken=("no grain" if b["taken_grain"] is None else f"grain {b['taken_grain']}"),
                 crossing=("no" if not b["taken_touches"] else "yes"), origin=b["origin"],
                 against_path=("yes" if b["against_path"] else "no"))
        )  # fmt: skip
    cols_br = [
        "voxel",
        "kind",
        "source_on",
        "err_on",
        "taken",
        "crossing",
        "origin",
        "against_path",
    ]
    t_br = DT.markdown_table(br_rows, cols_br) if br_rows else "no broken voxel"
    ct_rows = []
    for c in S["candidate_trace"]:
        ct_rows.append(
            dict(voxel=c["full_sample_voxel"], truth_rank=c["truth_rank"] or "none",
                 taken_rank=c["taken_rank"],
                 taken_cost=f"{c['taken_cost'][0]:.3f} to {c['taken_cost'][1]:.3f}",
                 truth_cost=("n/a" if "truth_cost" not in c else
                             f"{c['truth_cost'][0]:.3f} to {c['truth_cost'][1]:.3f}"),
                 cost_at_truth=f"{c['cost_at_truth']:.4f}",
                 truth_start_deg=("n/a" if "truth_start_deg" not in c
                                  else f"{c['truth_start_deg']:.3f}"))
        )  # fmt: skip
    cols_ct = ["voxel", "truth_rank", "taken_rank", "taken_cost", "truth_cost", "cost_at_truth",
               "truth_start_deg"]  # fmt: skip
    t_ct = (
        DT.markdown_table(ct_rows, cols_ct) if ct_rows else "no wrong indexed result behind a break"
    )
    # kinds
    k = pt["kinds_wrong"]
    t_kind = DT.markdown_table(
        [dict(run="key off (stored rerun)", **k["off"]), dict(run="key on", **k["on"])],
        ["run", "edge", "corner", "interior"],
    )
    fx = [x["kind"] for x in pt["fixed"]]
    bk = [x["kind"] for x in pt["broken"]]
    t_fxbk = DT.markdown_table(
        [dict(flip="fixed", **{kk: fx.count(kk) for kk in ("edge", "corner", "interior")}),
         dict(flip="broken", **{kk: bk.count(kk) for kk in ("edge", "corner", "interior")})],
        ["flip", "edge", "corner", "interior"],
    )  # fmt: skip
    DT.write_tables(
        OUT / "tables.md",
        {
            "ar2b_accuracy": t_acc,
            "ar2b_evals": t_ev,
            "ar2b_index": t_ix,
            "ar2b_times": t_tm,
            "ar2b_wall": t_wall,
            "ar2b_breaks": t_br,
            "ar2b_candidates": t_ct,
            "ar2b_kinds": t_kind,
            "ar2b_flip_kinds": t_fxbk,
        },
    )


def figure(S: Dict[str, Any]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    off, on, pt = S["off"], S["on"], S["pair"]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ix = on["index"]
    labs = ["key off\n(stored rerun)", "key on"]
    # evaluations by stage, millions
    full_off = off["evals"]["full"] / 1e6
    parts_on = [ix["evals_used"], ix["evals_wasted"], ix["evals_fallback_grid"]]
    local = [off["evals"]["local"] / 1e6, on["evals"]["local"] / 1e6]
    ax[0].bar(0, full_off, color="C0", label="full searches (grid)")
    b = 0.0
    names = ("full searches: indexing, used", "indexing that fell back", "grid fallbacks")
    for v, c, lab in zip(parts_on, ("C1", "C3", "C0"), names):
        ax[0].bar(1, v / 1e6, bottom=b, color=c, label=lab)
        b += v / 1e6
    ax[0].bar([0, 1], local, bottom=[full_off, b], color="C2", label="local fits")
    for x, tot in ((0, off["evals"]["total"]), (1, on["evals"]["total"])):
        ax[0].text(x, tot / 1e6 + 0.3, f"{tot / 1e6:.2f} M", ha="center")
    ax[0].set_xticks([0, 1], labs)
    ax[0].set_ylabel("cost evaluations (millions)")
    ax[0].set_title("Evaluations of the 2,000-voxel pilot")
    ax[0].legend(fontsize=8)
    # paired split
    p = pt["paired"]
    names = ["fixed", "broken", "wrong in both"]
    vals = [p["fixed"], p["broken"], p["both_wrong"]]
    ax[1].bar(names, vals, color=["C2", "C3", "C7"])
    for i, v in enumerate(vals):
        ax[1].text(i, v + 0.5, str(v), ha="center")
    ax[1].set_ylabel("voxels (wrong = more than 1 degree)")
    ax[1].set_title(
        f"Key on against key off, by voxel (wrong: {pt['wrong_off']} off, {pt['wrong_on']} on)"
    )
    fig.tight_layout()
    fig.savefig(OUT / "ar2b_key.png", dpi=130)


if __name__ == "__main__":
    main()
