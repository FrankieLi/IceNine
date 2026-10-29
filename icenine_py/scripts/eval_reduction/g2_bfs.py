"""Task G2: helpers of the BFS pilot runs (scripts/phase_d/pilot_bfs.py does the runs).

subregion  write a sub-region of the pilot region grid (the N voxels nearest the region's
           centre, stable sort, ties by file index) for short pilots and the key-off check
compare    key-off check: two run folders of pilot_bfs.py are bit-identical (recon.mic,
           records.json without the wall times, run.json stats without the wall times)

uv run python scripts/eval_reduction/g2_bfs.py subregion --n 150 --out <grid.mic>
uv run python scripts/eval_reduction/g2_bfs.py compare <run dir A> <run dir B> [--json OUT]
"""

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

ICE = Path(__file__).resolve().parents[2]
REGION = ICE / "benchmarks" / "phase_d_pilot"
CENTRE_XY = (0.17, 0.0)  # region.json centre_xy_mm


def subregion(n: int, out: Path) -> None:
    lines = (REGION / "region_grid.mic").read_text().splitlines()
    head, rows = lines[0], [r for r in lines[1:] if r.strip()]
    xy = np.array([[float(t) for t in r.split()[:2]] for r in rows])
    # triangle centroid: the first two columns are a vertex; the region was selected by centroid,
    # but for a compact patch the vertex distance is enough (the patch is not compared to it)
    d = np.hypot(xy[:, 0] - CENTRE_XY[0], xy[:, 1] - CENTRE_XY[1])
    keep = np.sort(np.argsort(d, kind="stable")[:n])
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join([head] + [rows[i] for i in keep]) + "\n")
    np.save(out.with_suffix(".rows.npy"), keep)  # region row of each sub-region voxel
    print(f"wrote {out}: {len(keep)} of {len(rows)} voxels")


def errors(run_dir: Path) -> Dict[str, np.ndarray]:
    """Per voxel of a pilot_bfs run (full region or a sub-region): the symmetry-reduced
    misorientation to the truth (deg), the region row, the source and the evaluations. The region
    row of a voxel is its index for the full region and comes from the sub-region's ``.rows.npy``
    otherwise (positions alone are ambiguous: an up and a down triangle share a vertex)."""
    import re
    import sys

    sys.path.insert(0, str(ICE / "scripts" / "phase_d"))
    from score_bfs import mats_to_quats, pair_misorientation_deg

    from icenine.mic_file import MicFile
    from icenine.orientation_search import get_symmetry_quaternions
    from icenine.symmetry import create_cubic_symmetry

    sym = get_symmetry_quaternions(create_cubic_symmetry(1.0))
    ex = ICE.parent / "Examples" / "Example2.ManyGrains" / "SimInput"
    idx = np.load(REGION / "region_index.npy")
    truth = MicFile.read(str(ex / "rand_500grains_1mm_neworient_s0.mic"))
    recon = MicFile.read(str(run_dir / "recon.mic"))
    cfg_p = run_dir / "run.config"
    m = re.search(r"(?m)^SampleFilename\s+(\S+)", cfg_p.read_text()) if cfg_p.exists() else None
    sample = Path(m.group(1)) if m else REGION / "region_grid.mic"  # stored runs: the region
    if sample.resolve() == (REGION / "region_grid.mic").resolve():
        row = np.arange(len(idx))
    else:
        row = np.load(sample.with_suffix(".rows.npy"))
    assert len(row) == len(recon.voxels)
    tv = [truth.voxels[idx[r]] for r in row]
    pos = np.array([v.position[:2] for v in tv], dtype=float)
    rpos = np.array([v.position[:2] for v in recon.voxels], dtype=float)
    assert np.abs(pos - rpos).max() < 1e-4, "run voxels do not match the region rows"
    q_true = mats_to_quats(np.array([v.orientation for v in tv]))
    q_rec = mats_to_quats(np.array([v.orientation for v in recon.voxels]))
    err = pair_misorientation_deg(q_rec, q_true, sym)
    records = json.loads((run_dir / "records.json").read_text())
    n = len(recon.voxels)
    src = np.array([records.get(str(i), {}).get("source", "not_visited") for i in range(n)])
    nev = np.array([records.get(str(i), {}).get("n_evals", 0) for i in range(n)])
    grain = np.load(ex / "rand_500grains_1mm_neworient_s0_grainmap.npy")[idx[row]]
    return dict(err=err, row=row, source=src, n_evals=nev, grain=grain)


def _records(p: Path) -> Dict[str, Any]:
    rec = json.loads((p / "records.json").read_text())
    for r in rec.values():
        r.pop("wall_s", None)
    return rec


def _stats(p: Path) -> Dict[str, Any]:
    st = json.loads((p / "run.json").read_text())["stats"]
    out = {k: v for k, v in st.items() if "wall" not in k and k != "reseeds"}
    out["reseeds"] = [
        {k: v for k, v in e.items() if "seconds" not in k} for e in st.get("reseeds", [])
    ]
    return out


def compare(a: Path, b: Path) -> Dict[str, Any]:
    same_mic = (a / "recon.mic").read_bytes() == (b / "recon.mic").read_bytes()
    ra, rb = _records(a), _records(b)
    sa, sb = _stats(a), _stats(b)
    diff_rec: List[str] = [k for k in sorted(set(ra) | set(rb)) if ra.get(k) != rb.get(k)]
    diff_st = [k for k in sorted(set(sa) | set(sb)) if sa.get(k) != sb.get(k)]
    res = dict(
        a=str(a.name),
        b=str(b.name),
        recon_mic_identical=bool(same_mic),
        n_records=len(ra),
        records_differing=diff_rec,
        stats_differing=diff_st,
        n_voxels=sa.get("n_voxels"),
        n_seeds=sa.get("n_seeds"),
        n_reseeds=sa.get("n_reseeds"),
        n_neighbor_fits=sa.get("n_neighbor_fits"),
        n_revisit_attempted=sa.get("n_revisit_attempted"),
        identical=bool(same_mic and not diff_rec and not diff_st),
    )
    return res


def run_summary(run_dir: Path, wrong_deg: float = 1.0) -> Dict[str, Any]:
    """Accuracy, provenance, evaluations by stage and the OutwardSearch log of one run."""
    sys_path = ICE / "scripts" / "common"
    import sys

    sys.path.insert(0, str(sys_path))
    from stats import wilson

    e = errors(run_dir)
    run = json.loads((run_dir / "run.json").read_text())
    st = run["stats"]
    err, src = e["err"], e["source"]
    wrong = err > wrong_deg
    n = len(err)
    lo, hi = wilson(int(wrong.sum()), n)
    out: Dict[str, Any] = dict(
        run=run_dir.name,
        n=n,
        wrong=int(wrong.sum()),
        wrong_wilson95=[lo, hi],
        unresolved=int((src == "unresolved").sum()),
        median_err_right_deg=float(np.median(err[~wrong])) if (~wrong).any() else float("nan"),
        sources={s: int((src == s).sum()) for s in sorted(set(src.tolist()))},
        seeds=st["n_seeds"],
        seed_rejected=st["n_seed_rejected"],
        reseeds=st["n_reseeds"],
        reseed_accepted=st["n_reseed_accepted"],
        neighbor_fits=st["n_neighbor_fits"],
        neighbor_accepted_first=st["n_neighbor_accepted_first"],
        revisit_attempted=st["n_revisit_attempted"],
        revisit_accepted=st["n_revisit_accepted"],
        evals=dict(
            seed=st["n_evals_seed"],
            reseed=st["n_evals_reseed"],
            neighbor=st["n_evals_neighbor"],
            revisit=st["n_evals_revisit"],
        ),
        wall=dict(
            seed=st["wall_seed_s"],
            reseed=st["wall_reseed_s"],
            neighbor=st["wall_neighbor_s"],
            revisit=st["wall_revisit_s"],
            total=run["wall_total_s"],
        ),
        timing_label=run.get("timing_label"),
        outward={k: v for k, v in st.items() if k.startswith("outward_")},
    )
    out["evals"]["total"] = sum(out["evals"].values())
    # grains: a grain is lost when none of its voxels is within wrong_deg
    g = e["grain"]
    lost = [int(k) for k in np.unique(g) if not (err[g == k] <= wrong_deg).any()]
    out["grains"] = int(len(np.unique(g)))
    out["grains_lost"] = lost
    log_p = run_dir / "outward_log.json"
    if log_p.exists():
        log = json.loads(log_p.read_text())
        by: Dict[str, Any] = {}
        for ph in sorted({x["phase"] for x in log}):
            rows = [x for x in log if x["phase"] == ph]
            acc = [x for x in rows if x["accepted"]]
            fb = [x for x in rows if x["fallback"]]
            by[ph] = dict(
                runs=len(rows),
                accepted=len(acc),
                fallback=len(fb),
                fallback_accepted=sum(1 for x in fb if x["accepted"]),
                median_wall_s=float(np.median([x["wall_s"] for x in rows])),
                median_wall_fallback_s=(
                    float(np.median([x["wall_s"] for x in fb])) if fb else float("nan")
                ),
                max_wall_s=float(max(x["wall_s"] for x in rows)),
                total_wall_s=float(sum(x["wall_s"] for x in rows)),
                strict=int(sum(x["n_strict"] for x in rows)),
                screen=int(sum(x["n_screen"] for x in rows)),
                smooth=int(sum(x["n_smooth"] for x in rows)),
                jac=int(sum(x["n_jac"] for x in rows)),
                end_cost_accepted_median=(
                    float(np.median([x["strict_cost"] for x in acc])) if acc else float("nan")
                ),
            )
            if ph == "reseed":
                by[ph]["seed_test"] = [x.get("seed_test") for x in rows]
                by[ph]["used_voxels_err_deg"] = [
                    float(err[x["voxel"]]) for x in rows if x.get("seed_test")
                ]
        out["outward_log"] = by
    return out


def paired(a: Path, b: Path, wrong_deg: float = 1.0) -> Dict[str, Any]:
    """Paired split by voxel (fixed / broken / both wrong) of run b against run a."""
    ea, eb = errors(a)["err"], errors(b)["err"]
    wa, wb = ea > wrong_deg, eb > wrong_deg
    return dict(
        a=a.name,
        b=b.name,
        wrong_a=int(wa.sum()),
        wrong_b=int(wb.sum()),
        fixed=int((wa & ~wb).sum()),
        broken=int((~wa & wb).sum()),
        both_wrong=int((wa & wb).sum()),
        broken_voxels=[int(i) for i in np.nonzero(~wa & wb)[0]],
        fixed_voxels=[int(i) for i in np.nonzero(wa & ~wb)[0]],
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["subregion", "compare", "pilot"])
    ap.add_argument("dirs", nargs="*")
    ap.add_argument("--n", type=int, default=150)
    ap.add_argument("--out", default=None)
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    if a.cmd == "subregion":
        subregion(a.n, Path(a.out))
        return
    if a.cmd == "compare":
        res: Any = compare(Path(a.dirs[0]), Path(a.dirs[1]))
    else:  # pilot: summaries of every run, paired against the first (the key-off run)
        dirs = [Path(d) for d in a.dirs]
        res = dict(
            runs=[run_summary(d) for d in dirs],
            paired=[paired(dirs[0], d) for d in dirs[1:]],
        )
    print(json.dumps(res, indent=1, default=float))
    if a.json:
        Path(a.json).write_text(json.dumps(res, indent=1, default=float) + "\n")


if __name__ == "__main__":
    main()
