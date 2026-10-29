"""Task G3 summary: the total of the evaluation-reduction plan, summary.json, tables.md and the
figure under benchmarks/eval_reduction/G/g3/, and the maps spec of the pilot region.

Parts (each skipped when its inputs are missing):
  bfs      the 2,000-voxel pilot region (clean images, mc arm, revisit and reseed on, BFS seed
           order 0): the stored pilot, the stored configuration on today's code (``keyoff_rerun``,
           Task G2's baseline), SmoothCostRank E3+E4 tau 0.25 alone (``smooth_only``, this task,
           D-20261029-N3), OutwardSearch alone and both (Task G2); every run paired by voxel against
           the stored pilot and against ``keyoff_rerun``, and ``integrated`` against
           ``smooth_only`` (OutwardSearch added to the smooth key)
  trace    every voxel a run breaks or fixes against ``keyoff_rerun`` (and the stored pilot
           against ``keyoff_rerun``, the flips of a code change alone), classed by the orientation
           it carries in the run that has it wrong:
             another grain   within 1 deg of the orientation of another grain of the sample (the
                             BFS carried a neighbouring grain over the boundary; a neighbour fit
                             and its acceptance test, where SmoothCostRank does not act)
             no grain        more than 1 deg from every grain's orientation
           and by origin: own full search, carried from a wrong full search, boundary crossing
           (a local fit carried a right full search's orientation over a grain boundary),
           untraced (see ``classify``); the run's full-search voxels within 1 deg are listed
  probe    Task E's broken-voxel rule (D-20261014-N13) for the SmoothCostRank-only arm's flips,
           from the probe re-runs of their full searches (``g3_trace_probe.py``)
  total    evaluations by stage before and after: six full-sample voxels and the E0 200 (no-start
           search, Task G2), Task G1's 24 full-sample voxels, and the BFS runs; smooth evaluations
           also in strict-equivalent units (Task E's measured time ratio of an E3+E4 evaluation to
           a strict one, full-sample images, batch of 256; Task G2's weights for the OutwardSearch
           policy)
  levers   one row per lever of Tasks A, C, E, F and G1, every number read from that task's
           tables (``doc_tables_all.md``) or summary, plus the lenient against C++ rule E0 split by
           voxel from Task E's stored baseline runs (E/cache/rank, gitignored; skipped if absent)
  figure   evaluations by stage before and after (G/g3/figures/evals_by_stage.png)
  maps     G/maps/spec.json with the five pilot-region runs (render with mic_maps_report.py)

  uv run python scripts/eval_reduction/g3_summary.py
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICE / "scripts" / "common"))
import doc_tables as DT  # noqa: E402
import g2_bfs  # noqa: E402
import g2_summary as G2S  # noqa: E402
from stats import wilson  # noqa: E402

ER = ICE / "benchmarks" / "eval_reduction"
G = ER / "G"
OUT = G / "g3"
G2 = G / "g2"
STORED = ICE / "benchmarks" / "phase_d_pilot_reseed_v2" / "runs" / "mc_clean"
RUNS2 = G2 / "bfs" / "runs"
RUNS3 = OUT / "bfs" / "runs"
KEYOFF = RUNS2 / "keyoff_rerun"
SMOOTH = RUNS3 / "smooth_only"
OUTWARD = RUNS2 / "outward_only"
INTEGRATED = RUNS2 / "integrated"
EX = ICE.parent / "Examples" / "Example2.ManyGrains" / "SimInput"
WRONG_DEG = 1.0
LABEL = {
    "mc_clean": "stored pilot (2026-10-10)",
    "keyoff_rerun": "stored configuration, today's code",
    "smooth_only": "SmoothCostRank only",
    "outward_only": "OutwardSearch only",
    "integrated": "both (Task G2's integrated configuration)",
}
STAGES = ["screen", "recost", "quick", "find", "variance", "final"]
BFS_STAGES = ["seed", "reseed", "neighbor", "revisit"]


# --------------------------------------------------------------------------------------------
# reading other tasks' tables
# --------------------------------------------------------------------------------------------
def md_rows(md: str) -> List[Dict[str, str]]:
    lines = [ln for ln in md.strip().splitlines() if ln.startswith("|")]
    head = [c.strip() for c in lines[0].strip("|").split("|")]
    out = []
    for ln in lines[2:]:
        cells = [c.strip() for c in ln.strip("|").split("|")]
        out.append(dict(zip(head, cells)))
    return out


def table(task_dir: str, name: str, file: str = "doc_tables_all.md") -> List[Dict[str, str]]:
    return md_rows(DT.read_tables(ER / task_dir / file)[name])


def row(rows: List[Dict[str, str]], **match: str) -> Dict[str, str]:
    hit = [r for r in rows if all(r.get(k) == v for k, v in match.items())]
    assert len(hit) == 1, (match, len(hit))
    return hit[0]


def num(s: str) -> float:
    """First number of a table cell ("41 (15.5-26.6%)", "56/200 [0.222, 0.346]", "48,700")."""
    import re

    m = re.search(r"-?\d[\d,]*(?:\.\d+)?", s)
    assert m, s
    return float(m.group(0).replace(",", ""))


# --------------------------------------------------------------------------------------------
# bfs
# --------------------------------------------------------------------------------------------
def bfs() -> Optional[Dict[str, Any]]:
    if not (SMOOTH / "run.json").exists():
        return None
    others = [KEYOFF, SMOOTH, OUTWARD, INTEGRATED]
    res = G2S.bfs_part(STORED, others, "evalG3_bfs", also_base=KEYOFF)
    assert res is not None
    extra = g2_bfs.paired(SMOOTH, INTEGRATED)  # OutwardSearch added on top of the smooth key
    res["paired"].append(extra)
    res["tables"]["evalG3_bfs_paired"] = G2S.paired_rows(res["paired"])
    run = json.loads((SMOOTH / "run.json").read_text())
    res["smooth_only_run"] = dict(
        timing_label=run["timing_label"],
        smooth_rank=run.get("smooth_rank"),
        outward_search=run.get("outward_search"),
        n_full_searches_counted=len(run.get("smooth_evals_per_full_search", [])),
        setup_s=run["setup_s"],
    )
    return res


# --------------------------------------------------------------------------------------------
# trace
# --------------------------------------------------------------------------------------------
_GRAINS: Dict[str, Any] = {}


def grain_quats() -> Tuple[np.ndarray, np.ndarray]:
    """(grain ids, one quaternion per grain) from the full truth .mic and its grain map; the
    largest spread of a grain's voxels from that quaternion is checked to be below 0.01 deg."""
    if not _GRAINS:
        sys.path.insert(0, str(ICE / "scripts" / "phase_d"))
        from score_bfs import mats_to_quats, pair_misorientation_deg

        from icenine.mic_file import MicFile
        from icenine.orientation_search import get_symmetry_quaternions
        from icenine.symmetry import create_cubic_symmetry

        sym = get_symmetry_quaternions(create_cubic_symmetry(1.0))
        truth = MicFile.read(str(EX / "rand_500grains_1mm_neworient_s0.mic"))
        gmap = np.load(EX / "rand_500grains_1mm_neworient_s0_grainmap.npy")
        q = mats_to_quats(np.array([v.orientation for v in truth.voxels]))
        ids = np.unique(gmap)
        first = np.array([int(np.nonzero(gmap == g)[0][0]) for g in ids])
        qg = q[first]
        spread = pair_misorientation_deg(q, qg[np.searchsorted(ids, gmap)], sym)
        assert spread.max() < 0.01, spread.max()
        xy = np.array([v.position[:2] for v in truth.voxels], dtype=float)
        _GRAINS.update(ids=ids, q=qg, sym=sym, spread_max=float(spread.max()), xy=xy, gmap=gmap)
    return _GRAINS["ids"], _GRAINS["q"]


def full_search_voxels(run_dir: Path) -> List[int]:
    """Voxels of a run whose orientation came from a full search: the seeds (records with source
    "seed") and the reseed voxels of run.json."""
    rec = json.loads((run_dir / "records.json").read_text())
    st = json.loads((run_dir / "run.json").read_text())["stats"]
    seeds = {int(k) for k, v in rec.items() if v.get("source") == "seed"}
    reseeds = {int(e["voxel"]) for e in st.get("reseeds", [])}
    return sorted(seeds | reseeds)


ORIGINS = ["own full search", "carried from a wrong full search", "boundary crossing", "untraced"]


def grain_distance_um(v_xy: np.ndarray, g: int) -> float:
    """Distance (um) from a voxel's vertex to the nearest truth voxel vertex of grain ``g``."""
    from scipy.spatial import cKDTree

    trees = _GRAINS.setdefault("trees", {})
    if g not in trees:
        trees[g] = cKDTree(_GRAINS["xy"][_GRAINS["gmap"] == g])
    return float(trees[g].query(v_xy)[0] * 1000.0)


def classify(run_dir: Path, voxels: List[int]) -> List[Dict[str, Any]]:
    """Class and origin of each voxel's orientation in ``run_dir`` (see the module docstring).

    origin: "own full search" (the voxel's record comes from its own full search, a seed or a
    reseed), "carried from a wrong full search" (a local fit; a full-search voxel of the run within
    1 deg of it is wrong for itself), "boundary crossing" (a local fit; every such full-search
    voxel is right for itself, so the BFS carried a right orientation over a grain boundary),
    "untraced" (no full-search voxel within 1 deg)."""
    from score_bfs import pair_misorientation_deg

    ids, qg = grain_quats()
    sym = _GRAINS["sym"]
    E = g2_bfs.errors(run_dir)
    q, err, grain, src = E["q_rec"], E["err"], E["grain"], E["source"]
    idx = np.load(ICE / "benchmarks" / "phase_d_pilot" / "region_index.npy")
    fs = full_search_voxels(run_dir)
    out = []
    for v in voxels:
        d = pair_misorientation_deg(np.repeat(q[v][None], len(qg), axis=0), qg, sym)
        j = int(np.argmin(d))
        dfs = pair_misorientation_deg(np.repeat(q[v][None], len(fs), axis=0), q[fs], sym)
        carriers = [int(fs[i]) for i in np.nonzero(dfs <= WRONG_DEG)[0]]
        wrong_carriers = [int(c) for c in carriers if err[c] > WRONG_DEG]
        cls = "another grain" if d[j] <= WRONG_DEG and int(ids[j]) != int(grain[v]) else "no grain"
        if str(src[v]) == "seed":
            origin = ORIGINS[0]
        elif wrong_carriers:
            origin = ORIGINS[1]
        elif carriers:
            origin = ORIGINS[2]
        else:
            origin = ORIGINS[3]
        full_idx = int(idx[E["row"][v]])
        out.append(
            dict(
                voxel=int(v),
                full_sample_voxel=full_idx,
                grain=int(grain[v]),
                err_deg=float(err[v]),
                source=str(src[v]),
                cls=cls,
                origin=origin,
                nearest_grain=int(ids[j]),
                nearest_grain_deg=float(d[j]),
                nearest_grain_distance_um=grain_distance_um(_GRAINS["xy"][full_idx], int(ids[j])),
                full_search_carriers=carriers,
                carriers_wrong=wrong_carriers,
            )
        )
    return out


def _split(xs: List[Dict[str, Any]]) -> Dict[str, int]:
    out = {o: sum(1 for x in xs if x["origin"] == o) for o in ORIGINS}
    out["another grain"] = sum(1 for x in xs if x["cls"] == "another grain")
    return out


def trace() -> Optional[Dict[str, Any]]:
    if not (SMOOTH / "run.json").exists():
        return None
    pairs = [
        (KEYOFF, SMOOTH),
        (KEYOFF, OUTWARD),
        (KEYOFF, INTEGRATED),
        (STORED, KEYOFF),
        (SMOOTH, INTEGRATED),
    ]
    res: Dict[str, Any] = {}
    rows = []
    cols = ["n", "another grain", *ORIGINS]
    for a, b in pairs:
        p = g2_bfs.paired(a, b)
        br = classify(b, p["broken_voxels"])  # wrong in b
        fx = classify(a, p["fixed_voxels"])  # wrong in a
        res[f"{a.name}->{b.name}"] = dict(broken=br, fixed=fx)
        for kind, xs in (("broken (wrong in b)", br), ("fixed (wrong in a)", fx)):
            rows.append(dict(a=LABEL[a.name], b=LABEL[b.name], flip=kind, n=len(xs), **_split(xs)))
    res["grain_spread_max_deg"] = _GRAINS["spread_max"]
    res["n_full_searches"] = {
        d.name: len(full_search_voxels(d)) for d in (STORED, KEYOFF, SMOOTH, OUTWARD, INTEGRATED)
    }
    # every wrong voxel of each run (not only the paired flips)
    wrong_rows = []
    for d in (STORED, KEYOFF, SMOOTH, OUTWARD, INTEGRATED):
        err = g2_bfs.errors(d)["err"]
        c = classify(d, [int(v) for v in np.nonzero(err > WRONG_DEG)[0]])
        res[f"wrong_{d.name}"] = c
        dist = [x["nearest_grain_distance_um"] for x in c if x["cls"] == "another grain"]
        wrong_rows.append(
            dict(
                run=LABEL[d.name],
                n=len(c),
                **_split(c),
                max_distance_um=max(dist) if dist else None,
                full_searches=res["n_full_searches"][d.name],
            )
        )
    res["tables"] = {
        "evalG3_trace": DT.markdown_table(rows, ["a", "b", "flip", *cols]),
        "evalG3_wrong_by_origin": DT.markdown_table(
            wrong_rows, ["run", *cols, "max_distance_um", "full_searches"],
            formats=dict(max_distance_um=".1f"),
        ),  # fmt: skip
    }
    return res


def probe(tr: Optional[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    """Task E's trace rule (D-20261014-N13) for the flips of the SmoothCostRank-only arm against
    the key-off rerun, from the probe re-runs of their full searches (``g3_trace_probe.py``).

    Verdict of a broken voxel: "boundary crossing" origin -> draw (the BFS path: a neighbour fit
    took a neighbouring grain from a right full search, where the key does not act); a full-search
    origin -> the probe of that full search: key-on wrong and key-off right -> Task E's ``trace``
    of the key-on run (arm / draw / untraced); otherwise (the probe outcome does not depend on the
    key) -> draw; no probe run -> untraced. Fixed voxels the same way with ``fixed_trace``."""
    sys.path.insert(0, str(HERE))
    import g3_trace_probe as TP

    if tr is None or not TP.VOXELS.exists():
        return None
    from smooth_rank_summary import fixed_trace, load_runs
    from smooth_rank_summary import trace as e_trace

    vox = json.loads(TP.VOXELS.read_text())["voxels"]
    base, arm = load_runs(TP.run_dir("base")), load_runs(TP.run_dir("e34t25"))
    if not all(v in base and v in arm for v in vox):
        return None
    idx = np.load(ICE / "benchmarks" / "phase_d_pilot" / "region_index.npy")
    pv: Dict[int, Dict[str, Any]] = {}
    for v in vox:
        eb, ea = float(base[v]["err_deg"]), float(arm[v]["err_deg"])
        d: Dict[str, Any] = dict(err_key_off=eb, err_key_on=ea, evals_key_off=int(base[v]["evals"]),
                                 evals_key_on=int(arm[v]["evals"]))  # fmt: skip
        if ea > WRONG_DEG and eb <= WRONG_DEG:
            d["broken_trace"] = list(e_trace(arm[v], 3))
        if eb > WRONG_DEG and ea <= WRONG_DEG:
            d["fixed_trace"] = fixed_trace(arm[v], 3)
        pv[v] = d

    def verdict(x: Dict[str, Any], kind: str) -> str:
        if x["origin"] == "boundary crossing":
            return "draw"
        if x["origin"] == "untraced":
            return "untraced"
        src = (
            [int(x["full_sample_voxel"])]
            if x["origin"] == "own full search"
            else [int(idx[c]) for c in x["carriers_wrong"]]
        )
        cls = []
        for s in src:
            if s not in pv:
                cls.append("untraced")
            elif kind == "broken" and "broken_trace" in pv[s]:
                cls.append(pv[s]["broken_trace"][0])
            elif kind == "fixed" and "fixed_trace" in pv[s]:
                cls.append(pv[s]["fixed_trace"])
            else:
                cls.append("draw")
        return "arm" if "arm" in cls else ("untraced" if "untraced" in cls else "draw")

    flips = tr[TP.FLIP]
    out: Dict[str, Any] = dict(voxels=pv)
    rows = []
    for kind in ("broken", "fixed"):
        vs = [dict(voxel=x["voxel"], origin=x["origin"], verdict=verdict(x, kind))
              for x in flips[kind]]  # fmt: skip
        out[kind] = vs
        rows.append(
            dict(
                flip=kind,
                n=len(vs),
                traced_to_arm=sum(1 for x in vs if x["verdict"] == "arm"),
                traced_to_draw=sum(1 for x in vs if x["verdict"] == "draw"),
                untraced=sum(1 for x in vs if x["verdict"] == "untraced"),
            )
        )

    def outcome(d: Dict[str, Any]) -> str:
        off, on = d["err_key_off"] > WRONG_DEG, d["err_key_on"] > WRONG_DEG
        if off == on:
            return "both wrong" if off else "both right"
        if on:
            return f"key breaks it (trace: {d['broken_trace'][0]})"
        return f"key fixes it (trace: {d['fixed_trace']})"

    prow = [
        dict(full_sample_voxel=v, **{k: d[k] for k in ("err_key_off", "err_key_on")},
             outcome=outcome(d))
        for v, d in pv.items()
    ]  # fmt: skip
    out["tables"] = {
        "evalG3_trace_n13": DT.markdown_table(
            rows, ["flip", "n", "traced_to_arm", "traced_to_draw", "untraced"]
        ),
        "evalG3_trace_probe": DT.markdown_table(
            prow, ["full_sample_voxel", "err_key_off", "err_key_on", "outcome"],
            formats=dict(err_key_off=".3f", err_key_on=".3f"),
        ),  # fmt: skip
    }
    return out


# --------------------------------------------------------------------------------------------
# total
# --------------------------------------------------------------------------------------------
def smooth_ratio() -> Dict[str, float]:
    """Task E's time of one E3+E4 evaluation in strict evaluations (batch of 256)."""
    r = row(table("E", "evalE_timing_ms"), cost="E3+E4")
    return dict(
        full=num(r["full-sample: x strict"]),
        e0=num(r["E0: x strict"]),
    )


def total(b: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    s2 = json.loads((G2 / "summary.json").read_text())
    k = smooth_ratio()
    st = s2["ledger"]["stages"]
    nostart = []
    for sname in STAGES:
        x, y = st["task_a_base"][sname], st["g2_integrated"][sname]
        nostart.append(dict(stage=sname, before=x, after=y, change=y / x - 1.0 if x else None))
    tb, ta = s2["ledger"]["totals"]["task_a"], s2["ledger"]["totals"]["integrated"]
    sm = s2["ledger"]["totals"]["smooth"]
    nostart.append(dict(stage="strict, all stages", before=tb, after=ta, change=ta / tb - 1.0))
    nostart.append(dict(stage="smooth key (E3+E4)", before=0, after=sm, change=None))
    eq = ta + sm * k["full"]
    nostart.append(
        dict(
            stage=f"strict-equivalent (smooth x {k['full']:.3f})",
            before=tb,
            after=eq,
            change=eq / tb - 1.0,
        )
    )
    T: Dict[str, str] = {
        "evalG3_total_nostart_stages": DT.markdown_table(
            nostart, ["stage", "before", "after", "change"],
            formats=dict(before=",.0f", after=",.0f", change="+.1%"),
        )  # fmt: skip
    }
    out: Dict[str, Any] = dict(smooth_x_strict=k, nostart_stages=nostart)
    # one row per voxel set: the key off against the configuration that G3 leaves on the table
    # (SmoothCostRank E3+E4 tau 0.25, OutwardSearch off)
    e0 = s2["e0"]
    g1 = row(
        table("G/g1", "evalG1_full24"), rule="lenient", **{"arm (SmoothCostRank)": "E3+E4,tau=0.25"}
    )
    sets = []

    def wl(kk: int, n: int) -> str:
        lo, hi = wilson(kk, n)
        return f"{kk} ({100 * lo:.1f}-{100 * hi:.1f}%)"

    sets.append(
        dict(
            set="six full-sample voxels (sum)",
            n=6,
            wrong_before=wl(s2["ledger"]["totals"]["wrong_task_a"], 6),
            wrong_after=wl(s2["ledger"]["totals"]["wrong_integrated"], 6),
            fixed_broken="0 / 0",
            strict_before=tb,
            strict_after=ta,
            smooth=sm,
            equiv_change=eq / tb - 1.0,
        )
    )
    n24 = int(g1["n"])
    sb24, sa24, sm24 = (num(g1[c]) for c in (
        "strict evaluations per voxel, strict search", "strict evaluations per voxel, arm",
        "smooth evaluations per voxel"))  # fmt: skip
    sets.append(
        dict(
            set="Task E's other 24 full-sample voxels (per voxel)",
            n=n24,
            wrong_before=g1["wrong, strict"],
            wrong_after=g1["wrong, arm"],
            fixed_broken=f"{g1['fixed']} / {g1['broken']}",
            strict_before=sb24,
            strict_after=sa24,
            smooth=sm24,
            equiv_change=(sa24 + sm24 * k["full"]) / sb24 - 1.0,
        )
    )
    eb, ea, es = e0["mean_evals_base"], e0["mean_evals_arm"], e0["mean_smooth_evals_arm"]
    sets.append(
        dict(
            set="E0 200 voxels (per voxel)",
            n=e0["n"],
            wrong_before=wl(e0["wrong_base"], e0["n"]),
            wrong_after=wl(e0["wrong_arm"], e0["n"]),
            fixed_broken=f"{e0['fixed']} / {e0['broken']}",
            strict_before=eb,
            strict_after=ea,
            smooth=es,
            equiv_change=(ea + es * k["e0"]) / eb - 1.0,
        )
    )
    out["sets"] = sets
    if b is not None:
        runs = {r["run"]: r for r in b["runs"]}
        prs = {(p["a"], p["b"]): p for p in b["paired"]}
        w = b["weights_ms"]
        bfs_rows = []
        for name in ("keyoff_rerun", "smooth_only", "outward_only", "integrated"):
            r = runs[name]
            pol = r.get("outward_log", {})
            pol_strict_screen = sum(o["strict"] + o["screen"] for o in pol.values())
            pol_equiv = sum(G2S.strict_equiv(o, w) for o in pol.values())
            equiv = (
                r["evals"]["total"]
                - pol_strict_screen
                + pol_equiv
                + r["smooth_full_searches"] * k["full"]
            )
            r["strict_equivalent"] = equiv
            p = prs.get(("keyoff_rerun", name))
            bfs_rows.append(
                dict(
                    run=LABEL[name],
                    wrong=wl(r["wrong"], r["n"]),
                    fixed_broken="-" if p is None else f"{p['fixed']} / {p['broken']}",
                    grains_lost=len(r["grains_lost"]),
                    **{s: r["evals"][s] for s in BFS_STAGES},
                    total=r["evals"]["total"],
                    smooth=r["smooth_full_searches"],
                    strict_equivalent=equiv,
                    wall_s=r["wall"]["total"],
                    timing=r["timing_label"],
                    # for the figure: strict evaluations outside the policy, by stage, and the
                    # policy's strict-equivalent (OutwardSearch runs)
                    strict_outside_policy={
                        s: r["evals"][s] - (pol[s]["strict"] + pol[s]["screen"] if s in pol else 0)
                        for s in BFS_STAGES
                    },
                    policy_equiv=pol_equiv,
                )
            )
        base_eq = runs["keyoff_rerun"]["strict_equivalent"]
        for x in bfs_rows:
            x["equiv_change"] = x["strict_equivalent"] / base_eq - 1.0
        out["bfs_rows"] = bfs_rows
        T["evalG3_total_bfs"] = DT.markdown_table(
            bfs_rows,
            ["run", "wrong", "fixed_broken", "grains_lost", *BFS_STAGES, "total", "smooth",
             "strict_equivalent", "equiv_change", "wall_s", "timing"],
            formats=dict(strict_equivalent=",.0f", equiv_change="+.1%", wall_s=",.0f",
                         **{s: ",d" for s in BFS_STAGES}, total=",d", smooth=",d"),
        )  # fmt: skip
        # time per full search (seeds and reseeds), and the time per smooth evaluation it implies
        # if a strict evaluation took as long as in the key-off run (an assumption: the runs'
        # contention differs)
        kr0 = runs["keyoff_rerun"]
        us_strict = (kr0["wall"]["seed"] + kr0["wall"]["reseed"]) / (
            kr0["evals"]["seed"] + kr0["evals"]["reseed"]
        )
        fst = []
        for name in ("keyoff_rerun", "smooth_only", "integrated"):
            r = runs[name]
            nf = r["seeds"] + r["reseeds"]
            fw = r["wall"]["seed"] + r["wall"]["reseed"]
            ev = r["evals"]["seed"] + r["evals"]["reseed"]
            sm = r["smooth_full_searches"]
            fst.append(
                dict(
                    run=LABEL[name],
                    full_searches=nf,
                    wall_per_full_search_s=fw / nf,
                    strict_per_full_search=ev / nf,
                    smooth_per_full_search=sm / nf,
                    implied_ms_per_smooth=(fw - ev * us_strict) / sm * 1e3 if sm else None,
                    timing=r["timing_label"],
                )
            )
        out["full_search_time"] = dict(rows=fst, us_per_strict_key_off=us_strict * 1e6)
        T["evalG3_full_search_time"] = DT.markdown_table(
            fst,
            ["run", "full_searches", "wall_per_full_search_s", "strict_per_full_search",
             "smooth_per_full_search", "implied_ms_per_smooth", "timing"],
            formats=dict(wall_per_full_search_s=".0f", strict_per_full_search=",.0f",
                         smooth_per_full_search=",.0f", implied_ms_per_smooth=".2f"),
        )  # fmt: skip
        # full-sample estimate (one process): full searches scale with the grains, local fits with
        # the voxels (an assumption, untested); wall times keep their runs' timing labels
        gm = np.load(EX / "rand_500grains_1mm_neworient_s0_grainmap.npy")
        n_vox, n_grains = int(len(gm)), int(len(np.unique(gm)))
        aff = []
        for name in ("keyoff_rerun", "smooth_only"):
            r = runs[name]
            fw = r["wall"]["seed"] + r["wall"]["reseed"]
            lw = r["wall"]["neighbor"] + r["wall"]["revisit"]
            sg, sv = n_grains / r["grains"], n_vox / r["n"]
            aff.append(
                dict(
                    run=LABEL[name],
                    full_searches=r["seeds"] + r["reseeds"],
                    full_search_wall_s=fw,
                    local_fit_wall_s=lw,
                    grain_scale=sg,
                    voxel_scale=sv,
                    estimate_h=(fw * sg + lw * sv) / 3600.0,
                    timing=r["timing_label"],
                )
            )
        out["full_sample_estimate"] = dict(voxels=n_vox, grains=n_grains, rows=aff)
        T["evalG3_full_sample_estimate"] = DT.markdown_table(
            aff,
            ["run", "full_searches", "full_search_wall_s", "local_fit_wall_s", "grain_scale",
             "voxel_scale", "estimate_h", "timing"],
            formats=dict(full_search_wall_s=",.0f", local_fit_wall_s=",.0f", grain_scale=".2f",
                         voxel_scale=".2f", estimate_h=".1f"),
        )  # fmt: skip
        sm_r = runs["smooth_only"]
        kr = runs["keyoff_rerun"]
        p = prs[("keyoff_rerun", "smooth_only")]
        sets.append(
            dict(
                set="BFS pilot region, 2,000 voxels (whole run)",
                n=sm_r["n"],
                wrong_before=wl(kr["wrong"], kr["n"]),
                wrong_after=wl(sm_r["wrong"], sm_r["n"]),
                fixed_broken=f"{p['fixed']} / {p['broken']}",
                strict_before=kr["evals"]["total"],
                strict_after=sm_r["evals"]["total"],
                smooth=sm_r["smooth_full_searches"],
                equiv_change=sm_r["strict_equivalent"] / kr["strict_equivalent"] - 1.0,
            )
        )
    T["evalG3_total_sets"] = DT.markdown_table(
        sets,
        ["set", "n", "wrong_before", "wrong_after", "fixed_broken", "strict_before",
         "strict_after", "smooth", "equiv_change"],
        formats=dict(strict_before=",.0f", strict_after=",.0f", smooth=",.0f",
                     equiv_change="+.1%"),
    )  # fmt: skip
    out["tables"] = T
    return out


# --------------------------------------------------------------------------------------------
# levers
# --------------------------------------------------------------------------------------------
def rule_split_e0() -> Optional[Dict[str, Any]]:
    """Lenient against C++ rule on the E0 200 (Task E's stored strict-search runs), by voxel."""
    root = ER / "E" / "cache" / "rank"
    a, b = root / "e0_old_base", root / "e0_cpp_base"
    if not (a.exists() and b.exists()):
        return None
    ea = {int(np.load(f)["voxel"]): float(np.load(f)["err_deg"]) for f in a.glob("v*.npz")}
    eb = {int(np.load(f)["voxel"]): float(np.load(f)["err_deg"]) for f in b.glob("v*.npz")}
    vox = sorted(set(ea) & set(eb))
    wa = np.array([ea[v] > WRONG_DEG for v in vox])
    wb = np.array([eb[v] > WRONG_DEG for v in vox])
    return dict(
        n=len(vox),
        wrong_lenient=int(wa.sum()),
        wrong_cpp=int(wb.sum()),
        cpp_fixes=int((wa & ~wb).sum()),
        cpp_breaks=int((~wa & wb).sum()),
        both_wrong=int((wa & wb).sum()),
    )


def levers(b: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    rows = []
    lm = table("A", "evalA_ledger_mean")
    base, cpp = row(lm, arm="base"), row(lm, arm="cpp")
    sv = table("A", "evalA_speed_voxel_quiet")
    rs = rule_split_e0()
    ip = table("A", "evalA_e0_interp_paired")
    e2e = table("E", "evalE_e2e")
    sp = rs or {}
    rows.append(
        dict(
            lever="C++ overlap rule (`CppOverlapRule 1`)",
            task="A",
            evaluations=f"{base['total']} -> {cpp['total']} per voxel (six full-sample voxels)",
            wrong=(
                f"E0: {sp['wrong_lenient']} -> {sp['wrong_cpp']} "
                f"({sp['cpp_fixes']} fixed / {sp['cpp_breaks']} broken)"
                if rs
                else "E0: not paired"
            ),
            time=f"{base['wall_s']} -> {cpp['wall_s']} s per voxel (mixed labels)",
            carried="no (lenient rule stays, D-20261013-I6)",
        )
    )
    q1 = row(sv, arm="cpp_edge", voxel="15901")
    q2 = row(sv, arm="cpp_edge_interp_batched", voxel="15901")
    len_interp = row(ip, new="lenint", base="off")
    rows.append(
        dict(
            lever=(
                "C++ edge rules, omega interpolation, batched screen (`CppEdgeRules`, "
                "`CppOmegaInterp`, `BatchedScreen`)"
            ),
            task="A",
            evaluations=(
                f"{base['total']} -> {row(lm, arm='lenient_interp_batched')['total']} per voxel "
                "(six voxels, lenient rule, all three)"
            ),
            wrong=(
                f"E0 (lenient rule, edge rules and interpolation): "
                f"{len_interp['wrong_base'].split('/')[0]} -> "
                f"{len_interp['wrong_new'].split('/')[0]} ({len_interp['fixed']} fixed / "
                f"{len_interp['broken']} broken)"
            ),
            time=(
                f"{q1['wall_s']} -> {q2['wall_s']} s, voxel 15901, C++ rule with `CppEdgeRules`, "
                "adding the other two (single-worker)"
            ),
            carried="`BatchedScreen` on by default since 2026-10-13; the other two opt-in",
        )
    )
    c_six = row(table("C", "evalC_c_six"), arm="c0a_k1")
    c_e0 = row(table("C", "evalC_c_e0_arms"), arm="c0a_k1")
    c_base = row(table("C", "evalC_c_e0_arms"), arm="base")
    rows.append(
        dict(
            lever="candidate de-duplication, best arm (`CandidateDedup clique=cpp:1`)",
            task="C",
            evaluations=f"{c_six['vs_base']} of the base (six voxels); {c_e0['vs_base']} on E0",
            wrong=(
                f"E0: {int(num(c_base['wrong']))} -> {int(num(c_e0['wrong']))} "
                f"({c_e0['fixed']} fixed / {c_e0['broken']} broken)"
            ),
            time="not measured single-worker",
            carried="no (D-20261014-N4)",
        )
    )
    tv = table("G/g1", "evalG1_cache_timing")
    for spec, carried in (
        ("E3+E4,tau=0.25", "yes (D-20261029-N1)"),
        ("E1+E2+E3+E4", "no (G1 default rule: fewer net fixes on the full sample)"),
    ):
        r0 = row(e2e, **{"voxel set": "E0 200", "rule": "lenient", "arm (SmoothCostRank)": spec})
        r6 = row(
            e2e, **{"voxel set": "six full-sample", "rule": "lenient", "arm (SmoothCostRank)": spec}
        )
        t = row(tv, **{"arm (SmoothCostRank)": spec})
        rows.append(
            dict(
                lever=f"smooth ranking key of levels 0-2 (`SmoothCostRank {spec}`)",
                task="E, G1",
                evaluations=(
                    f"E0 {num(r0['strict evaluations per voxel, baseline']):,.0f} -> "
                    f"{num(r0['strict evaluations per voxel, arm']):,.0f} strict + "
                    f"{num(r0['smooth evaluations per voxel']):,.0f} smooth per voxel; six "
                    f"voxels {num(r6['strict evaluations per voxel, baseline']):,.0f} -> "
                    f"{num(r6['strict evaluations per voxel, arm']):,.0f} + "
                    f"{num(r6['smooth evaluations per voxel']):,.0f}"
                ),
                wrong=(
                    f"E0: {int(num(r0['wrong, baseline']))} -> {int(num(r0['wrong, arm']))} "
                    f"({r0['fixed']} fixed / {r0['broken']} broken)"
                ),
                time=(
                    f"{t['extra time over key off, warm cache, A / B'].split('/')[1].strip()} "
                    "per full-sample voxel with the cache (voxel 17034, single-worker)"
                ),
                carried=carried,
            )
        )
    cb = row(table("G/g1", "evalG1_cache_build"), timing="single-worker")
    saved = row(tv, **{"arm (SmoothCostRank)": "E3+E4,tau=0.25"})
    rows.append(
        dict(
            lever="SmoothCost cache per image set",
            task="G1",
            evaluations="none (the same results)",
            wrong="no change (bit-identical)",
            time=(
                f"saves the {cb['build total (s)']} s build per voxel "
                f"({saved['saved per voxel by the cache (s), A / B']} s measured, single-worker)"
            ),
            carried="yes, under the key",
        )
    )
    if b is not None:
        runs = {r["run"]: r for r in b["runs"]}
        prs = {(p["a"], p["b"]): p for p in b["paired"]}
        o, kr = runs["outward_only"], runs["keyoff_rerun"]
        p = prs[("keyoff_rerun", "outward_only")]
        rows.append(
            dict(
                lever="outward search for BFS neighbours and reseeds (`OutwardSearch 1`)",
                task="F, G2",
                evaluations=(
                    f"BFS region {kr['evals']['total'] / 1e6:.1f}M -> "
                    f"{o['evals']['total'] / 1e6:.1f}M (mostly screen evaluations)"
                ),
                wrong=f"BFS region: {kr['wrong']} -> {o['wrong']} ({p['fixed']} fixed / "
                f"{p['broken']} broken)",
                time=f"{kr['wall']['total']:,.0f} -> {o['wall']['total']:,.0f} s (contended)",
                carried="no (D-20261029-N2)",
            )
        )
    return dict(
        rule_split_e0=rs,
        rows=rows,
        tables={
            "evalG3_levers": DT.markdown_table(
                rows, ["lever", "task", "evaluations", "wrong", "time", "carried"]
            )
        },
    )


# --------------------------------------------------------------------------------------------
# figure
# --------------------------------------------------------------------------------------------
def figure(t: Dict[str, Any]) -> Optional[str]:
    if "bfs_rows" not in t:
        return None
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    k = t["smooth_x_strict"]["full"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.6), gridspec_kw=dict(width_ratios=[1, 2]))
    cols = dict(
        screen="#4c72b0", recost="#8fb0d9", quick="#dd8452", find="#55a868", variance="#8c8c8c",
        final="#000000", seed="#4c72b0", reseed="#8fb0d9", neighbor="#dd8452",
        revisit="#f0b98d", smooth="#c44e52", policy="#937860",
    )  # fmt: skip
    # left: no-start search, six full-sample voxels
    st = {r["stage"]: r for r in t["nostart_stages"]}
    for i, when in enumerate(("before", "after")):
        y = 0.0
        for s in ["screen", "recost", "quick", "find", "variance"]:
            v = st[s][when] / 1e6
            ax1.bar(i, v, bottom=y, color=cols[s], label=s if i == 0 else None, width=0.6)
            y += v
        if when == "after":
            v = st["smooth key (E3+E4)"]["after"] * k / 1e6
            ax1.bar(i, v, bottom=y, color=cols["smooth"], hatch="//", width=0.6,
                    label=f"smooth key (x {k:.2f}, strict-equivalent)")  # fmt: skip
            y += v
        ax1.text(i, y * 1.01, "0 of 6 wrong", ha="center", va="bottom", fontsize=9)
    ax1.set_xticks([0, 1], ["key off\n(Task A)", "SmoothCostRank\nE3+E4 tau 0.25"])
    ax1.set_ylabel("evaluations (millions, strict-equivalent)")
    ax1.set_title("No-start search, six full-sample voxels (sum)", fontsize=10)
    ax1.legend(fontsize=8, loc="lower center", bbox_to_anchor=(0.5, -0.45), ncol=2)
    ax1.set_ylim(0, max(st["strict-equivalent" + f" (smooth x {k:.3f})"]["after"], 1) / 1e6 * 1.15)
    # right: BFS pilot region
    rows = t["bfs_rows"]
    names = {
        "stored configuration, today's code": "key off\n(stored config,\ntoday's code)",
        "SmoothCostRank only": "SmoothCostRank\nonly",
        "OutwardSearch only": "OutwardSearch\nonly",
        "both (Task G2's integrated configuration)": "both\n(G2 integrated)",
    }
    for i, r in enumerate(rows):
        y = 0.0
        for s in BFS_STAGES:
            v = r["strict_outside_policy"][s] / 1e6
            ax2.bar(i, v, bottom=y, color=cols[s], label=s if i == 0 else None, width=0.6)
            y += v
        if r["policy_equiv"]:
            v = r["policy_equiv"] / 1e6
            ax2.bar(i, v, bottom=y, color=cols["policy"], hatch="..", width=0.6,
                    label="outward-search policy (strict-equivalent)")  # fmt: skip
            y += v
        if r["smooth"]:
            v = r["smooth"] * k / 1e6
            ax2.bar(i, v, bottom=y, color=cols["smooth"], hatch="//", width=0.6,
                    label=f"smooth key (x {k:.2f}, strict-equivalent)")  # fmt: skip
            y += v
        ax2.text(i, y * 1.01, f"{r['wrong'].split(' ')[0]} wrong\n{r['fixed_broken']}",
                 ha="center", va="bottom", fontsize=8)  # fmt: skip
    ax2.set_xticks(range(len(rows)), [names.get(r["run"], r["run"]) for r in rows], fontsize=8)
    ax2.set_ylabel("evaluations (millions, strict-equivalent)")
    ax2.set_title(
        "BFS, 2,000-voxel pilot region\n(labels: wrong voxels; fixed / broken against key off)",
        fontsize=10,
    )
    h, lab = ax2.get_legend_handles_labels()
    seen = dict(zip(lab, h))
    ax2.legend(seen.values(), seen.keys(), fontsize=8, loc="lower center",
               bbox_to_anchor=(0.5, -0.42), ncol=3)  # fmt: skip
    ax2.set_ylim(0, max(r["strict_equivalent"] for r in rows) / 1e6 * 1.2)
    fig.tight_layout()
    p = OUT / "figures" / "evals_by_stage.png"
    p.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(p, dpi=130)
    plt.close(fig)
    return str(p.relative_to(ICE))


# --------------------------------------------------------------------------------------------
# maps spec
# --------------------------------------------------------------------------------------------
def maps_spec(b: Dict[str, Any]) -> None:
    runs = {r["run"]: r for r in b["runs"]}
    prs = {(p["a"], p["b"]): p for p in b["paired"]}
    names = ["mc_clean", "keyoff_rerun", "smooth_only", "outward_only", "integrated"]
    wrong = "; ".join(f"{LABEL[n]}: {runs[n]['wrong']}" for n in names)
    pr = "; ".join(
        f"{LABEL[n]}: {prs[('keyoff_rerun', n)]['fixed']} fixed and "
        f"{prs[('keyoff_rerun', n)]['broken']} broken"
        for n in names[2:]
    )
    lost = "; ".join(f"{LABEL[n]}: {len(runs[n]['grains_lost'])}" for n in names)
    spec = {
        "title": "Task G pilot maps",
        "summary": (
            "The 2,000-voxel pilot region reconstructed five times with the same BFS seed order: "
            "the stored pilot (2026-10-10), the same configuration on today's code (the "
            "like-for-like baseline: the centroid screen became the default on 2026-10-12), "
            "SmoothCostRank E3+E4 tau 0.25 alone (Task G3), OutwardSearch alone, and both (Task "
            f"G2's integrated configuration). Wrong voxels (more than 1 degree from the truth, of "
            f"2,000): {wrong}. Grains lost (of {runs['keyoff_rerun']['grains']}): {lost}."
        ),
        "reference": {
            "mic": "../../../../../Examples/Example2.ManyGrains/SimInput/"
            "rand_500grains_1mm_neworient_s0.mic",
            "index": "../../../phase_d_pilot/region_index.npy",
            "label": "Original (truth)",
        },
        "z": None,
        "figures": [
            {
                "name": "pilot_region",
                "title": "2,000-voxel pilot region, clean images, mc arm, reseed on",
                "caption": (
                    "Clean images, mc arm, revisit and reseed on, BFS rng seed 0. Paired by voxel "
                    f"against the stored configuration on today's code: {pr}."
                ),
                "mics": [
                    {"path": "../../../phase_d_pilot_reseed_v2/runs/mc_clean/recon.mic",
                     "label": "Stored pilot"},
                    {"path": "../g2/bfs/runs/keyoff_rerun/recon.mic",
                     "label": "Stored configuration, today's code"},
                    {"path": "../g3/bfs/runs/smooth_only/recon.mic",
                     "label": "SmoothCostRank only"},
                    {"path": "../g2/bfs/runs/outward_only/recon.mic",
                     "label": "OutwardSearch only"},
                    {"path": "../g2/bfs/runs/integrated/recon.mic",
                     "label": "Both (integrated configuration)"},
                ],  # fmt: skip
            }
        ],
        "notes": [
            "Clean simulated images (no realistic noise); one BFS seed order (rng seed 0) per run.",
            "Wrong means more than 1 degree from the truth after cubic symmetry.",
        ],
    }
    (G / "maps" / "spec.json").write_text(json.dumps(spec, indent=2) + "\n")


# --------------------------------------------------------------------------------------------
def main() -> None:
    summary: Dict[str, Any] = {}
    tables: Dict[str, str] = {}
    b = bfs()
    if b is not None:
        tables.update(b.pop("tables"))
        summary["bfs"] = b
        maps_spec(b)
    tr = trace()
    if tr is not None:
        tables.update(tr.pop("tables"))
        summary["trace"] = tr
    pr = probe(tr)
    if pr is not None:
        tables.update(pr.pop("tables"))
        summary["probe"] = pr
    t = total(b)
    tables.update(t.pop("tables"))
    summary["total"] = t
    lv = levers(b)
    tables.update(lv.pop("tables"))
    summary["levers"] = lv
    fig = figure(t)
    if fig:
        summary["figure"] = fig
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, default=float) + "\n")
    DT.write_tables(OUT / "tables.md", tables)
    print(f"parts: {sorted(summary)}; tables: {sorted(tables)}")


if __name__ == "__main__":
    main()
