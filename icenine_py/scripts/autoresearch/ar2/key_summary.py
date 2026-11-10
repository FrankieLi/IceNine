"""AR2 step 2: the 2,000-voxel BFS pilot with ``SeedByIndexing 1`` against the stored key-off rerun.

  uv run python scripts/autoresearch/ar2/key_summary.py baseline   # recompute the stored numbers
  uv run python scripts/autoresearch/ar2/key_summary.py all        # summary.json, tables.md, figure

Baseline (from the stored files, before anything else; nothing is rerun): wrong voxels, evaluations
and full searches of ``benchmarks/eval_reduction/G/g2/bfs/runs/keyoff_rerun``, recomputed from its
``recon.mic``, ``records.json`` and ``run.json``. The stop rule of the plan: if they differ from
37 wrong, 25,547,935 evaluations and 69 full searches, stop.

Run folders read from ``benchmarks/autoresearch/ar2/key/runs/<name>`` (copied from the gitignored
cache by the run script): ``key_on`` (the accuracy run), ``key_on_fused`` and ``key_off_fused``
(the FusedCost-on pair).

The trace (``describe``) works on the stored run files only, the same way for both sides of a pair
(D-20261109-N9): the broken voxels in the key-on run and the fixed voxels in the key-off run. It is
an orientation match. A voxel's "matches" are the voxels of the same run whose stored result is
their own full search (``source == "seed"``: an accepted seed or reseed) and lies within 1 degree
of the voxel's result. The "nearest match" is nearest in orientation. ``records.json`` stores no
parent, so that a match is where a voxel's orientation came from (through local fits) is inferred,
not traced. Seeds the seed test rejected keep no result of their own in ``recon.mic`` and cannot be
matches.
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
# The run script wrote "single run" into every run.json. The key-on (FusedCost off) run did not run
# alone: its log spans 17:14:15 to 17:35:28 and the voxel 17361 run's log 17:18:43 to 17:19:06 (the
# birth and last-change times of both logs in the gitignored cache), and short checks whose times
# were not logged ran beside it. The table carries this corrected label; run.json is left as stored.
WALL_LABEL_FIX = {
    "key_on": "contended: timing preflight refused (an editor process at about 200% CPU); not a "
    "single run: overlapped with the voxel 17361 run (its log spans 23 s) and with short checks "
    "whose times were not logged"
}


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
    """Wrong voxels, evaluations and full searches of the stored key-off rerun, recomputed from its
    stored files (no rerun)."""
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
    out["recomputed_matches_stored"] = bool(ok)
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
        self.centre = self.vxy.mean(axis=1)  # (N, 2) triangle centroids, mm
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

    def dist_sides(self, a: int, b: int) -> float:
        """Distance between the centroids of full-sample voxels a and b, in voxel sides of a."""
        return float(np.linalg.norm(self.centre[a] - self.centre[b]) / self.side[a])


# --------------------------------------------------------------------------------------------
# paired comparison and the trace of every break
# --------------------------------------------------------------------------------------------
class RunView:
    """One stored run, read for the trace: errors, sources, orientations and the index log."""

    def __init__(self, d: Path, index_log: List[Dict[str, Any]]) -> None:
        from icenine.mic_file import MicFile

        self.dir = d
        self.E = g2_bfs.errors(d)
        self.wrong = self.E["err"] > WRONG_DEG
        self.R = np.array(
            [v.orientation for v in MicFile.read(str(d / "recon.mic")).voxels], dtype=np.float64
        )
        self.log = {int(e["voxel"]): e for e in index_log}
        self.indexed = {v for v, e in self.log.items() if e["used"]}
        # voxels whose stored result is their own full search (an accepted seed or reseed); a seed
        # the seed test rejected keeps a later result, so its own full search is not in recon.mic
        self.own = [int(i) for i in np.nonzero(self.E["source"] == "seed")[0]]
        self.own_set = set(self.own)
        assert self.indexed <= self.own_set, "an indexed seed whose result is not its own"
        rec = json.loads((d / "records.json").read_text())
        # seeds and reseeds the seed test rejected (records.json ``seed_rejected``)
        self.rejected = {int(k) for k, v in rec.items() if v.get("seed_rejected")}
        self.full_search = self.own_set | self.rejected


def describe(r: int, V: RunView, L: Layer) -> Dict[str, Any]:
    """Trace of region voxel r in run V (see the module docstring): the grain its result matches,
    whether that grain touches the voxel, and the full searches of the run its result matches."""
    from score_bfs import pair_misorientation_deg

    E = V.E
    q = E["q_rec"]
    f = int(L.region[E["row"][r]])
    g_taken = L.grain_of(V.R[r])
    d: Dict[str, Any] = dict(
        voxel=int(r),
        full_sample_voxel=f,
        grain=int(E["grain"][r]),
        err=float(E["err"][r]),
        source=str(E["source"][r]),
        own_indexed_seed=bool(r in V.indexed),
        full_search=bool(r in V.full_search),
        seed_rejected=bool(r in V.rejected),
        kind=L.kind(f),
        taken_grain=g_taken,
        taken_touches=(None if g_taken is None else L.touches(f, g_taken)),
    )
    if r in V.own_set:
        what = "indexed seed" if r in V.indexed else "grid seed"
        d["origin"] = f"own {what} ({'wrong' if V.wrong[r] else 'right'})"
        car = [int(r)]
    else:
        dfs = pair_misorientation_deg(np.repeat(q[r][None], len(V.own), axis=0), q[V.own], E["sym"])
        car = [int(V.own[i]) for i in np.nonzero(dfs <= WRONG_DEG)[0]]
        wrong_car = [c for c in car if V.wrong[c]]
        if any(c in V.indexed for c in wrong_car):
            d["origin"] = "matches a wrong indexed seed"
        elif wrong_car:
            d["origin"] = "matches a wrong grid seed"
        elif car:
            d["origin"] = "matches a right full search"
        else:
            d["origin"] = "no full search matches"
    d["matches"] = car
    wrong_car = [c for c in car if V.wrong[c]]
    d["matches_wrong_indexed"] = [c for c in wrong_car if c in V.indexed]
    d["matches_wrong_grid"] = [c for c in wrong_car if c not in V.indexed]
    # a wrong full search of the run is behind the voxel's result (its own, or a match): "indexed"
    # (the key's path, N7/N8), "grid" (the grid search), or "no"
    d["wrong_seed_behind"] = (
        "indexed" if d["matches_wrong_indexed"] else "grid" if d["matches_wrong_grid"] else "no"
    )
    d["against_path"] = d["wrong_seed_behind"] == "indexed"
    if car:
        if r in V.own_set:
            c = int(r)
        else:
            dm = pair_misorientation_deg(np.repeat(q[r][None], len(car), axis=0), q[car], E["sym"])
            c = int(car[int(np.argmin(dm))])
        d["nearest_match"] = c
        d["nearest_match_full_sample_voxel"] = int(L.region[E["row"][c]])
        d["nearest_match_err_deg"] = float(E["err"][c])
        d["nearest_match_dist_sides"] = L.dist_sides(f, int(L.region[E["row"][c]]))
    return d


def paired_trace(on: Dict[str, Any], L: Layer, off_dir: Path = OFF) -> Dict[str, Any]:
    """Paired split of the key-on run against a key-off run, with every flip traced the same way on
    its wrong side (D-20261109-N9): the broken voxels in the key-on run, the fixed voxels in the
    key-off run."""
    A, B = RunView(off_dir, []), RunView(on["dir"], on["index_log"])
    ea, eb = A.E["err"], B.E["err"]
    wa, wb = A.wrong, B.wrong
    pr = g2_bfs.paired(off_dir, on["dir"])
    idx = L.region
    used = B.indexed
    fallback = {int(e["voxel"]) for e in on["index_log"] if not e["used"]}

    def flip(
        r: int, wrong_side: RunView, other: RunView, tag: str, other_tag: str
    ) -> Dict[str, Any]:
        d = describe(r, wrong_side, L)
        out = {k: v for k, v in d.items() if k not in ("err", "source")}
        out[f"err_{tag}"] = d["err"]
        out[f"source_{tag}"] = d["source"]
        out[f"err_{other_tag}"] = float(other.E["err"][r])
        out[f"source_{other_tag}"] = str(other.E["source"][r])
        out[f"own_indexed_seed_{other_tag}"] = bool(r in other.indexed)
        for k in ("own_indexed_seed", "full_search", "seed_rejected"):
            out[f"{k}_{tag}"] = d[k]
        out[f"full_search_{other_tag}"] = bool(r in other.full_search)
        out[f"seed_rejected_{other_tag}"] = bool(r in other.rejected)
        return out

    broken = [flip(r, B, A, "on", "off") for r in pr["broken_voxels"]]
    fixed_rows = [flip(r, A, B, "off", "on") for r in pr["fixed_voxels"]]
    both = [int(r) for r in np.nonzero(wa & wb)[0]]
    both_rows = [flip(r, B, A, "on", "off") for r in both]
    # every used indexed seed or reseed: right or wrong
    used_err = {int(r): float(eb[r]) for r in sorted(used)}
    n_wrong_used = int(sum(e > WRONG_DEG for e in used_err.values()))
    # kinds of all wrong voxels of each run
    kinds = {}
    for name, w, V in (("off", wa, A), ("on", wb, B)):
        ks = [L.kind(int(idx[V.E["row"][r]])) for r in np.nonzero(w)[0]]
        kinds[name] = {k: ks.count(k) for k in ("edge", "corner", "interior")}
    # the fixed and broken voxels, counted by what the trace finds on their wrong side
    flips: Dict[str, Any] = {}
    for name, rows in (("broken", broken), ("fixed", fixed_rows)):
        flips[name] = dict(
            n=len(rows),
            crossing=sum(bool(x["taken_touches"]) for x in rows),
            origin={
                o: sum(x["origin"] == o for x in rows) for o in sorted({x["origin"] for x in rows})
            },
            wrong_seed_behind={
                k: sum(x["wrong_seed_behind"] == k for x in rows) for k in ("indexed", "grid", "no")
            },
        )
    # N9: the key's own path at a flipped voxel is the chain on that voxel when it is itself an
    # indexed seed of the key-on run (right for a fix, wrong for a break)
    for name, rows in (("broken", broken), ("fixed", fixed_rows)):
        fz = flips[name]
        fz["indexed_seed_on"] = sum(x["own_indexed_seed_on"] for x in rows)
        fz["grid_seed_off"] = sum(x["full_search_off"] for x in rows)
        fz["grid_seed_off_rejected"] = sum(x["seed_rejected_off"] for x in rows)
        # the key's own path at the voxel (N9): an indexed seed of the key-on run, right for a fix
        # and wrong for a break; listed with what the key-off run did at the same voxel
        fz["key_path_voxels"] = [
            dict(voxel=x["full_sample_voxel"], source_off=x["source_off"],
                 grid_seed_off=x["full_search_off"], rejected_off=x["seed_rejected_off"])
            for x in rows if x["own_indexed_seed_on"]
        ]  # fmt: skip
    for name, rows in (("broken", broken), ("fixed", fixed_rows)):
        for tag in ("off", "on"):
            srcs = [x[f"source_{tag}"] for x in rows]
            flips[name][f"source_{tag}"] = {s: srcs.count(s) for s in sorted(set(srcs))}
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
        flips=flips,
        median_err_right_on=float(np.median(eb[~wb])),
        median_err_right_off=float(np.median(ea[~wa])),
        unresolved_on=int((B.E["source"] == "unresolved").sum()),
        unresolved_off=int((A.E["source"] == "unresolved").sum()),
        own_full_search_voxels_on=len(B.own),
        own_full_search_voxels_off=len(A.own),
    )


def ar1_comparison(on: Dict[str, Any], L: Layer) -> Dict[str, Any]:
    """The indexing paths of this run against AR1's stored ones (``ar1/e2e.py``,
    ``benchmarks/autoresearch/ar1/e2e_pilot_base_part0of1.json``, ``recall_pilot_base.json``) on the
    voxels that are both an indexed seed of this run and one of AR1's 69 full-search voxels.

    AR1 paired the rays with the hull of the full sample; the BFS pilot's sample is the 2,000-voxel
    region, so its hull (and with it the paired rays, ``sample_hull_points`` in ``SmoothCost``) is
    smaller. A voxel is "same rays" when the numbers of attached rays are equal. Only those can be
    expected to give the identical result."""
    rows = json.loads((AC.OUT / "e2e_pilot_base_part0of1.json").read_text())["rows"]
    ar1 = {int(r["voxel"]): r for r in rows if not r.get("no_candidates")}
    rec_rows = json.loads((AC.OUT / "recall_pilot_base.json").read_text())["rows"]
    natt = {int(r["voxel"]): int(r["n_att"]) for r in rec_rows}
    out: Dict[str, Any] = dict(
        n_common=0, same_rays=0, same_rays_identical=0, other_rays=0, other_rays_identical=0
    )
    detail = []
    for e in on["index_log"]:
        f = int(L.region[int(e["voxel"])])
        if f not in ar1:
            continue
        a = ar1[f]
        fins = a["finishers"][:2]
        ident = (
            np.array_equal(np.array(a["strict_costs"][:5]), e["strict_costs"])
            and [f_["cand"] for f_ in fins] == e["refined"]
            and [f_["cost"] for f_ in fins] == e["finisher_costs"]
            and [f_["evals"] for f_ in fins] == e["finisher_evals"]
        )
        same_rays = natt[f] == int(e["n_att"])
        out["n_common"] += 1
        k = "same_rays" if same_rays else "other_rays"
        out[k] += 1
        out[k + "_identical"] += int(ident)
        detail.append(
            dict(voxel=f, n_att_ar1=natt[f], n_att_run=int(e["n_att"]), identical=bool(ident))
        )
    out["detail"] = detail
    return out


def trace_flips(
    rows: List[Dict[str, Any]], run_dir: Path, index_log: List[Dict[str, Any]], tag: str
) -> List[Dict[str, Any]]:
    """Every flipped voxel, in full, on its wrong side (``tag`` "on" for a break in the key-on run,
    "off" for a fix in the key-off run): the strict cost of the voxel's stored result and at its
    exact truth, and its nearest match (nearest in orientation among the run's own full searches
    within 1 degree, see the module docstring): its error, its distance from the voxel in voxel
    sides, and, for an indexed seed, the candidate taken and its strict cost before and after the
    finisher. Read from the stored files; the only computation is the strict cost."""
    from icenine.mic_file import MicFile

    log = {int(e["voxel"]): e for e in index_log}
    recon = MicFile.read(str(run_dir / "recon.mic"))
    out = []
    for b in rows:
        r, f = b["voxel"], b["full_sample_voxel"]
        e = AC.env("full", f, need_sc=False)
        lf, verts, phase, Rt = e["lf"], e["verts"], e["phase"], e["R_true"]
        Rb = np.asarray(recon.voxels[r].orientation, dtype=np.float64)
        row: Dict[str, Any] = dict(
            voxel=f,
            region_voxel=r,
            source=b[f"source_{tag}"],
            result_cost=float(lf.evaluate(Rb, verts, phase).cost),
            cost_at_truth=float(lf.evaluate(Rt, verts, phase).cost),
            taken_grain=b["taken_grain"],
            taken_touches=b["taken_touches"],
            kind=b["kind"],
            err_deg=b[f"err_{tag}"],
        )
        if "nearest_match" in b:
            c = int(b["nearest_match"])
            row.update(
                match=b["nearest_match_full_sample_voxel"],
                match_err_deg=b["nearest_match_err_deg"],
                match_dist_sides=b["nearest_match_dist_sides"],
                match_is_self=bool(c == r),
            )
            if c in log:
                ent = log[c]
                k = ent["refined"].index(ent["taken"])
                row.update(
                    match_phase=ent["phase"],
                    match_taken_rank=int(ent["taken"]) + 1,
                    match_cost=[
                        float(ent["strict_costs"][ent["taken"]]),
                        float(ent["finisher_costs"][k]),
                    ],
                )
        out.append(row)
    return out


def recon_difference(a: Path, b: Path, sym: np.ndarray) -> Dict[str, Any]:
    """How two stored runs of the same arm differ voxel by voxel: orientations (the robust angle of
    ``ray_indexing.misorientation_deg``, from the .mic files, whose Euler angles are written with
    4 decimals for the first angle and 6 for the other two, so smaller differences are not seen),
    the .mic voxel lines (any column), and evaluation counts (``records.json``).

    ``arccos_angle_nonzero`` counts the voxels where the arccos-based angle of ``score_bfs`` is
    above zero. It is above zero also for identical orientations (its floor is about 1.7e-6 deg),
    so it is not a count of differing orientations; it is kept to explain that difference."""
    from icenine.mic_file import MicFile

    Ra = np.array([v.orientation for v in MicFile.read(str(a / "recon.mic")).voxels], dtype=float)
    Rb = np.array([v.orientation for v in MicFile.read(str(b / "recon.mic")).voxels], dtype=float)
    diff = np.nonzero(np.any(Ra != Rb, axis=(1, 2)))[0]
    ang = RI.misorientation_deg(Ra[diff], Rb[diff], sym) if len(diff) else np.zeros(0)
    Ea, Eb = g2_bfs.errors(a), g2_bfs.errors(b)  # puts scripts/phase_d on the path
    from score_bfs import pair_misorientation_deg

    acos = pair_misorientation_deg(Ea["q_rec"], Eb["q_rec"], Ea["sym"])
    same = np.ones(len(Ra), dtype=bool)
    same[diff] = False
    la = (a / "recon.mic").read_text().splitlines()[1:]
    lb = (b / "recon.mic").read_text().splitlines()[1:]
    lines_differ = int(sum(x != y for x, y in zip(la, lb)))
    ra = json.loads((a / "records.json").read_text())
    rb = json.loads((b / "records.json").read_text())
    n = len(Ra)
    ne = [
        int(ra.get(str(i), {}).get("n_evals", 0)) != int(rb.get(str(i), {}).get("n_evals", 0))
        for i in range(n)
    ]
    return dict(
        n_voxels=int(n),
        orientation_differs=int(len(diff)),
        angle_deg_median_of_differing=float(np.median(ang)) if len(ang) else 0.0,
        angle_deg_max=float(ang.max()) if len(ang) else 0.0,
        angle_above_0p001_deg=int((ang > 1e-3).sum()),
        angle_above_0p01_deg=int((ang > 1e-2).sum()),
        mic_lines_differ=lines_differ,
        evals_differ=int(sum(ne)),
        arccos_angle_nonzero=int((acos > 0).sum()),
        arccos_angle_nonzero_identical_orientation=int(((acos > 0) & same).sum()),
        arccos_angle_max_identical_orientation_deg=float(acos[same].max()) if same.any() else 0.0,
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
        need.update(b["matches_wrong_indexed"])
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
        timing_label=WALL_LABEL_FIX.get(d.name, run.get("timing_label")),
        timing_label_stored=run.get("timing_label"),
        torch_threads=run.get("torch_threads"),
        fused_cost=run.get("fused_cost"),
        seed_by_indexing=run.get("seed_by_indexing"),
    )
    out["full_searches"] = out["seeds"] + out["reseeds"]
    out["seeds_accepted"] = out["seeds"] - out["seed_rejected"]
    out["reseeds_accepted"] = out["reseeds"] - out["reseed_rejected"]
    E = g2_bfs.errors(d)
    out["wrong"] = int((E["err"] > WRONG_DEG).sum())
    is_seed = E["source"] == "seed"  # voxels whose final result is their own full search
    out["seed_source"] = int(is_seed.sum())
    out["seed_source_wrong"] = int((E["err"][is_seed] > WRONG_DEG).sum())
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
        # best 2 against best 1, as seen in this run: the second finisher run changed the taken
        # candidate (taken is not the first refined one); the finisher left the cost unchanged
        log = r["index_log"]
        out["index"]["taken_not_first_refined"] = int(
            sum(int(e["taken"]) != int(e["refined"][0]) for e in log)
        )
        out["index"]["finisher_cost_equals_strict_cost"] = int(
            sum(
                float(e["finisher_costs"][e["refined"].index(e["taken"])])
                == float(e["strict_costs"][e["taken"]])
                for e in log
            )
        )
    return out


def main() -> None:
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    base = baseline()
    print(json.dumps(base, indent=1), flush=True)
    if not base["recomputed_matches_stored"]:
        raise SystemExit("STOP: the numbers recomputed from the stored baseline files differ")
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
        trace_method=(
            "orientation match on the stored run files: a voxel's matches are the voxels of the "
            "same run whose stored result is their own full search (source seed) within 1 degree "
            "of its result; nearest match = nearest in orientation; no parent is recorded, so "
            "'reached through local fits from the match' is inferred, not traced; rejected seeds "
            "keep no own result and cannot be matches; breaks traced in the key-on run, fixes in "
            "the key-off run (D-20261109-N9)"
        ),
        off=run_block("key off (stored rerun)", OFF),
        on=run_block("key on", on_dir),
        pair=pt,
        candidate_trace=cand,
        break_trace=trace_flips(pt["broken"], on_dir, on["index_log"], "on"),
        fix_trace=trace_flips(pt["fixed"], OFF, [], "off"),
        ar1_comparison=ar1_comparison(on, L),
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
        # FusedCost is not bit-for-bit (float32 rounding at pixel edges, D-20261027-I2)
        summ["on_fused"]["difference_vs_on"] = recon_difference(on_dir, fused, L.sym)
    off_fused = RUNS / "key_off_fused"
    if (off_fused / "run.json").exists():
        summ["off_fused"] = run_block("key off, FusedCost on", off_fused)
        summ["off_fused"]["difference_vs_off"] = recon_difference(OFF, off_fused, L.sym)
        if "on_fused" in summ:
            onf = load_run(fused)
            ptf = paired_trace(onf, L, off_fused)
            summ["pair_fused"] = ptf
            pn = g2_bfs.paired(OFF, off_fused)
            summ["noise"] = {
                k: pn[k] for k in ("wrong_a", "wrong_b", "fixed", "broken", "both_wrong")
            }
            summ["break_trace_fused"] = trace_flips(ptf["broken"], fused, onf["index_log"], "on")
            summ["fix_trace_fused"] = trace_flips(ptf["fixed"], off_fused, [], "off")
    n_fits_off = summ["off"]["neighbor_fits"] + summ["off"]["revisit_attempted"]
    d: Dict[str, Any] = dict(
        evals_ratio_off_over_on=summ["off"]["evals"]["total"] / summ["on"]["evals"]["total"],
        evals_per_grid_full_search_off=summ["off"]["evals"]["full"] / summ["off"]["full_searches"],
        evals_per_indexed_seed_on=summ["on"]["index"]["evals_used"] / summ["on"]["index"]["runs"],
        evals_per_local_fit_off=summ["off"]["evals"]["local"] / n_fits_off,
        full_search_share_of_evals_off=summ["off"]["evals"]["full"] / summ["off"]["evals"]["total"],
        full_search_share_of_evals_on=summ["on"]["evals"]["full"] / summ["on"]["evals"]["total"],
    )
    # Wilson intervals of the small counts quoted on the page
    ix = summ["on"]["index"]
    fb = ix["fallback_no_candidate"] + ix["fallback_seed_test"]
    d["fallback_wilson"] = list(ST.wilson(fb, ix["runs"]))
    d["own_seed_wrong_wilson_on"] = list(
        ST.wilson(summ["on"]["seed_source_wrong"], summ["on"]["seed_source"])
    )
    d["own_seed_wrong_wilson_off"] = list(
        ST.wilson(summ["off"]["seed_source_wrong"], summ["off"]["seed_source"])
    )
    lo_on, hi_on = d["own_seed_wrong_wilson_on"]
    lo_off, hi_off = d["own_seed_wrong_wilson_off"]
    d["own_seed_wrong_intervals_overlap"] = bool(lo_on <= hi_off and lo_off <= hi_on)
    # seeds and reseeds that failed the seed test, out of all seeds and reseeds
    for name in ("off", "on", "off_fused", "on_fused"):
        if name in summ:
            r = summ[name]
            k = r["seed_rejected"] + r["reseed_rejected"]
            d[f"seed_test_failed_{name}"] = [k, r["full_searches"]]
            d[f"seed_test_failed_wilson_{name}"] = list(ST.wilson(k, r["full_searches"]))
    if "on_fused" in summ and "off_fused" in summ:
        d["evals_ratio_fused_off_over_on"] = (
            summ["off_fused"]["evals"]["total"] / summ["on_fused"]["evals"]["total"]
        )
        d["wall_ratio_fused_off_over_on"] = (
            summ["off_fused"]["wall_total_s"] / summ["on_fused"]["wall_total_s"]
        )
        d["wall_fused_on_seed_reseed_s"] = (
            summ["on_fused"]["index"]["wall_seed_s"] + summ["on_fused"]["index"]["wall_reseed_s"]
        )
        d["wall_fused_on_local_fit_s"] = (
            summ["on_fused"]["index"]["wall_neighbor_s"]
            + summ["on_fused"]["index"]["wall_revisit_s"]
        )
    summ["derived"] = d
    (OUT / "summary.json").write_text(json.dumps(summ, indent=1, default=float) + "\n")
    tables(summ)
    figure(summ)


def flip_tables(rows: List[Dict[str, Any]], trace: List[Dict[str, Any]], tag: str) -> Any:
    """(table of the flipped voxels with the trace's finding, table of each one's trace in full).
    ``tag`` is the wrong side: "on" for breaks (traced in the key-on run), "off" for fixes (traced
    in the key-off run)."""
    other = "off" if tag == "on" else "on"
    empty = "no broken voxel" if tag == "on" else "no fixed voxel"
    fl_rows = []
    for b in rows:
        fl_rows.append(
            {"voxel": b["full_sample_voxel"], "kind": b["kind"],
             f"source_{tag}": b[f"source_{tag}"], f"source_{other}": b[f"source_{other}"],
             f"err_{tag}": f"{b[f'err_{tag}']:.2f}",
             "taken": ("no grain" if b["taken_grain"] is None else f"grain {b['taken_grain']}"),
             "crossing": ("no" if not b["taken_touches"] else "yes"), "match": b["origin"],
             "wrong_seed_behind": b["wrong_seed_behind"],
             "seed_off": ("rejected" if b["seed_rejected_off"] else
                          "accepted" if b["full_search_off"] else "no"),
             "seed_on": ("indexed" if b["own_indexed_seed_on"] else "no")}
        )  # fmt: skip
    cols_fl = ["voxel", "kind", f"source_{tag}", f"source_{other}", f"err_{tag}", "taken",
               "crossing", "match", "wrong_seed_behind", "seed_off", "seed_on"]  # fmt: skip
    t_fl = DT.markdown_table(fl_rows, cols_fl) if fl_rows else empty
    tr_rows = []
    for c in trace:
        tr_rows.append(
            dict(voxel=c["voxel"], source=c["source"], err_deg=f"{c['err_deg']:.1f}",
                 result_cost=f"{c['result_cost']:.3f}", cost_at_truth=f"{c['cost_at_truth']:.3f}",
                 match=("none" if "match" not in c else
                        "itself" if c["match_is_self"] else c["match"]),
                 match_err=("n/a" if "match_err_deg" not in c else f"{c['match_err_deg']:.3f}"),
                 match_dist=("n/a" if "match_dist_sides" not in c
                             else f"{c['match_dist_sides']:.1f}"),
                 match_rank=c.get("match_taken_rank", "n/a"),
                 match_cost=("n/a" if "match_cost" not in c else
                             f"{c['match_cost'][0]:.3f} to {c['match_cost'][1]:.3f}"))
        )  # fmt: skip
    cols_tr = ["voxel", "source", "err_deg", "result_cost", "cost_at_truth", "match", "match_err",
               "match_dist", "match_rank", "match_cost"]  # fmt: skip
    t_tr = DT.markdown_table(tr_rows, cols_tr) if tr_rows else empty
    return t_fl, t_tr


def flips_table(S: Dict[str, Any]) -> str:
    """One row per pair and flip: what the trace finds on the flipped voxels' wrong side."""
    rows = []
    pairs = [("stored pair", S["pair"])]
    if "pair_fused" in S:
        pairs.append(("FusedCost-on pair", S["pair_fused"]))
    for pname, p in pairs:
        for flip, side in (("fixed", "key off"), ("broken", "key on")):
            fz = p["flips"][flip]
            o = fz["origin"]
            n_right = sum(v for k, v in o.items() if k.endswith("(right)")) + o.get(
                "matches a right full search", 0
            )
            rows.append(
                dict(pair=pname, flip=flip, traced_in=side, n=fz["n"], crossing=fz["crossing"],
                     right_seed_behind=n_right,
                     wrong_grid_seed_behind=fz["wrong_seed_behind"]["grid"],
                     wrong_indexed_seed_behind=fz["wrong_seed_behind"]["indexed"],
                     no_match=o.get("no full search matches", 0),
                     grid_seed_off=(f"{fz['grid_seed_off']} ({fz['grid_seed_off_rejected']} "
                                    "rejected)"),
                     indexed_seed_on=fz["indexed_seed_on"])
            )  # fmt: skip
    cols = ["pair", "flip", "traced_in", "n", "crossing", "right_seed_behind",
            "wrong_grid_seed_behind", "wrong_indexed_seed_behind", "no_match", "grid_seed_off",
            "indexed_seed_on"]  # fmt: skip
    return str(DT.markdown_table(rows, cols))


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
    if "pair_fused" in S:
        pf = S["pair_fused"]
        rows += [
            dict(run="key off, FusedCost on", wrong=wil(S["off_fused"]["wrong"], n), fixed="",
                 broken="", both_wrong="", unresolved=pf["unresolved_off"]),
            dict(run="key on, FusedCost on", wrong=wil(S["on_fused"]["wrong"], n),
                 fixed=pf["paired"]["fixed"], broken=pf["paired"]["broken"],
                 both_wrong=pf["paired"]["both_wrong"], unresolved=pf["unresolved_on"]),
        ]  # fmt: skip
    t_acc = DT.markdown_table(rows, ["run", "wrong", "fixed", "broken", "both_wrong", "unresolved"])
    # evaluations
    ev_rows = []
    arms = [("key off (stored rerun)", off), ("key on", on)]
    if "pair_fused" in S:
        arms += [("key off, FusedCost on", S["off_fused"]), ("key on, FusedCost on", S["on_fused"])]
    for name, r in arms:
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
        dict(item="fallbacks to the grid search", value=wil(fb, ix["runs"])),
        dict(item="  no candidate", value=ix["fallback_no_candidate"]),
        dict(item="  failed the seed test", value=ix["fallback_seed_test"]),
        dict(item="evaluations of used indexing paths", value=f"{ix['evals_used']:,d}"),
        dict(item="evaluations of indexing paths that fell back", value=f"{ix['evals_wasted']:,d}"),
        dict(item="evaluations of the grid fallbacks", value=f"{ix['evals_fallback_grid']:,d}"),
        dict(item="key off: grid seeds and reseeds that failed the seed test",
             value=wil(*S["derived"]["seed_test_failed_off"])),
        dict(item="key on: indexed seeds and reseeds that failed the seed test",
             value=wil(*S["derived"]["seed_test_failed_on"])),
        dict(item="key off: voxels whose result is their own grid search, wrong",
             value=wil(off["seed_source_wrong"], off["seed_source"])),
        dict(item="key on: voxels whose result is their own indexed seed, wrong",
             value=wil(on["seed_source_wrong"], on["seed_source"])),
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
    t_br, t_bt = flip_tables(pt["broken"], S["break_trace"], "on")
    t_fx, t_ft = flip_tables(pt["fixed"], S["fix_trace"], "off")
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
    # written only when a wrong indexed seed is behind a break or a voxel wrong in both runs (none
    # in these runs: the block was a placeholder and is no longer written)
    t_ct = DT.markdown_table(ct_rows, cols_ct) if ct_rows else ""
    t_diff = ""
    diffs = []
    if "on_fused" in S:
        diffs.append(("key on: FusedCost on against off", S["on_fused"]["difference_vs_on"]))
    if "off_fused" in S:
        diffs.append(("key off: FusedCost on against off", S["off_fused"]["difference_vs_off"]))
    if diffs:
        t_diff = DT.markdown_table(
            [dict(pair=nm, orientation_differs=f"{x['orientation_differs']} of {x['n_voxels']}",
                  median_deg=f"{x['angle_deg_median_of_differing']:.1e}",
                  above_0_001_deg=x["angle_above_0p001_deg"],
                  above_0_01_deg=x["angle_above_0p01_deg"],
                  max_deg=f"{x['angle_deg_max']:.3f}",
                  mic_lines_differ=f"{x['mic_lines_differ']} of {x['n_voxels']}",
                  evals_differ=f"{x['evals_differ']} of {x['n_voxels']}")
             for nm, x in diffs],
            ["pair", "orientation_differs", "median_deg", "above_0_001_deg", "above_0_01_deg",
             "max_deg", "mic_lines_differ", "evals_differ"],
        )  # fmt: skip
    t_noise = ""
    if "noise" in S:
        nz = S["noise"]
        t_noise = DT.markdown_table(
            [dict(pair="key off with FusedCost on, against key off without", wrong_a=nz["wrong_a"],
                  wrong_b=nz["wrong_b"], fixed=nz["fixed"], broken=nz["broken"],
                  both_wrong=nz["both_wrong"])],
            ["pair", "wrong_a", "wrong_b", "fixed", "broken", "both_wrong"],
        )  # fmt: skip
    t_17361 = ""
    p17 = OUT / "voxel_17361.json"
    if p17.exists():
        v = json.loads(p17.read_text())
        fo, vm = v["stages"]["findoptimal"], v["stages"]["variance_minimizing"]
        rs, gn = v["restart_from_end"], v["ray_gauss_newton_from_end"]
        pert = v["strict_cost_at_truth_perturbed"]
        items = [
            ("truth candidate: rank in the indexer's list", v["truth_rank"]),
            ("truth candidate: distance from the truth (deg)", f"{v['truth_start_deg']:.4f}"),
            ("truth candidate: strict cost", f"{v['truth_candidate_cost']:.4f}"),
            ("strict cost at the exact truth", f"{v['cost_at_exact_truth']:.4f}"),
            ("strict cost 0.02 deg from the truth", f"{pert['0.02deg']:.4f}"),
            ("strict cost 0.045 deg from the truth", f"{pert['0.045deg']:.4f}"),
            ("strict cost 0.1 deg from the truth", f"{pert['0.1deg']:.4f}"),
            ("finisher, total evaluations (truth candidate)", v["finisher_evals"][1]),
            ("stage 1 (FindOptimal, MC blocks): evaluations", fo["evals"]),
            ("stage 1: blocks run, blocks accepted",
             f"{int(fo['n_blocks'])}, {int(fo['n_accept'])}"),
            ("stage 1: first MC step and box (deg)", f"{v['stages']['step_deg']:.3f}, "
                                                    f"{v['stages']['box_deg']:.3f}"),
            ("stage 1: cost before, after", f"{fo['cost_before']:.4f}, {fo['cost_after']:.4f}"),
            ("stage 2 (VarianceMinimizing): evaluations", vm["evals"]),
            ("stage 2: cost before, after", f"{vm['cost_before']:.4f}, {vm['cost_after']:.4f}"),
            ("restart of the finisher from its end point: cost, moved (deg)",
             f"{rs['cost']:.4f}, {rs['moved_deg']:.4f}"),
            ("ray Gauss-Newton from the end point: cost, accepted",
             f"{gn['cost']:.4f}, {gn['accepted']}"),
            ("ray Gauss-Newton: moved (deg), distance from the truth after (deg)",
             f"{gn['moved_deg']:.4f}, {gn['err_deg']:.4f}"),
            ("ray Gauss-Newton: strict, smooth evaluations", f"{gn['n_strict']}, {gn['n_smooth']}"),
        ]  # fmt: skip
        t_17361 = DT.markdown_table([dict(item=a, value=b) for a, b in items], ["item", "value"])
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
    blocks: Dict[str, str] = {
        "ar2b_accuracy": t_acc,
        "ar2b_evals": t_ev,
        "ar2b_index": t_ix,
        "ar2b_times": t_tm,
        "ar2b_wall": t_wall,
        "ar2b_flips": flips_table(S),
        "ar2b_breaks": t_br,
        "ar2b_break_trace": t_bt,
        "ar2b_fixes": t_fx,
        "ar2b_fix_trace": t_ft,
        "ar2b_kinds": t_kind,
        "ar2b_flip_kinds": t_fxbk,
    }
    if t_ct:
        blocks["ar2b_candidates"] = t_ct
    if "pair_fused" in S:
        pf = S["pair_fused"]
        blocks["ar2b_breaks_fused"], blocks["ar2b_break_trace_fused"] = flip_tables(
            pf["broken"], S["break_trace_fused"], "on"
        )
        blocks["ar2b_fixes_fused"], blocks["ar2b_fix_trace_fused"] = flip_tables(
            pf["fixed"], S["fix_trace_fused"], "off"
        )
    if t_diff:
        blocks["ar2b_fused_diff"] = t_diff
    if t_noise:
        blocks["ar2b_noise"] = t_noise
    if t_17361:
        blocks["ar2b_17361"] = t_17361
    DT.write_tables(OUT / "tables.md", blocks)


def figure(S: Dict[str, Any]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    off, on, pt = S["off"], S["on"], S["pair"]
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.6), gridspec_kw=dict(width_ratios=[1, 1.5]))
    ix = on["index"]
    labs = ["key off\n(stored rerun)", "key on"]
    # evaluations by stage, millions
    full_off = off["evals"]["full"] / 1e6
    parts_on = [ix["evals_used"], ix["evals_wasted"], ix["evals_fallback_grid"]]
    local = [off["evals"]["local"] / 1e6, on["evals"]["local"] / 1e6]
    ax[0].bar(0, full_off, color="C0", label="full searches (grid)")
    b = 0.0
    names = ("full searches: indexing, used", "indexing that fell back", "grid fallbacks")
    for v, c, lab in zip(parts_on, ("C1", "C3", "C4"), names):
        if v > 0:  # a part that is zero in this run is neither drawn nor listed
            ax[0].bar(1, v / 1e6, bottom=b, color=c, label=lab)
        b += v / 1e6
    ax[0].text(1.0, b / 2 + 0.6, f"indexing: {ix['evals_used'] / 1e6:.2f} M", ha="center",
               fontsize=8)  # fmt: skip
    ax[0].bar([0, 1], local, bottom=[full_off, b], color="C2", label="local fits")
    for x, tot in ((0, off["evals"]["total"]), (1, on["evals"]["total"])):
        ax[0].text(x, tot / 1e6 + 0.3, f"{tot / 1e6:.2f} M", ha="center")
    ax[0].set_xticks([0, 1], labs)
    ax[0].set_ylabel("cost evaluations (millions)")
    ax[0].set_title("Evaluations of the 2,000-voxel pilot")
    ax[0].legend(fontsize=8)
    # paired splits: the stored pair, the FusedCost-on pair, and the key-off noise pair (two key-off
    # runs that differ only in FusedCost), so the key's split is read beside the noise
    pairs = [
        (f"key on vs key off\n(stored rerun)\nwrong {pt['wrong_off']} -> {pt['wrong_on']}",
         pt["paired"]),
    ]  # fmt: skip
    if "pair_fused" in S:
        pf = S["pair_fused"]
        pairs.append(
            (f"key on vs key off\n(both FusedCost on)\nwrong {pf['wrong_off']} -> {pf['wrong_on']}",
             pf["paired"])
        )  # fmt: skip
    if "noise" in S:
        nz = S["noise"]
        pairs.append(
            ("key off, FusedCost on\nvs key off (stored)\n"
             f"wrong {nz['wrong_a']} -> {nz['wrong_b']}", nz)
        )  # fmt: skip
    names = ["fixed", "broken", "wrong in both"]
    keys = ["fixed", "broken", "both_wrong"]
    w = 0.26
    for j, (nm, kk, col) in enumerate(zip(names, keys, ["C2", "C3", "C7"])):
        xs = np.arange(len(pairs)) + (j - 1) * w
        vals = [int(p[kk]) for _, p in pairs]
        ax[1].bar(xs, vals, width=w, color=col, label=nm)
        for x, v in zip(xs, vals):
            ax[1].text(x, v + 0.4, str(v), ha="center", fontsize=9)
    ax[1].set_xticks(np.arange(len(pairs)), [nm for nm, _ in pairs], fontsize=8)
    ax[1].set_ylabel("voxels (wrong = more than 1 degree)")
    ax[1].set_title("Paired by voxel: the key's two pairs and the key-off noise pair")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "ar2b_key.png", dpi=130)


if __name__ == "__main__":
    main()
