"""Task B steps 2 and 4: tables from the ray walks of cost_edges_run.py.

  uv run python scripts/eval_reduction/cost_edges_summary.py --out benchmarks/eval_reduction/B

Reads benchmarks/eval_reduction/B/cache/{full,e0}/*.npz; writes summary.json (compact, NaN as
null) to --out (cost_edges_tables.py turns it into tables.md). Definitions (also in
docs/findings_cost_edges_2026-10.md):
  event            one peak whose evaluation changed between two consecutive points of a ray, by
                   type
  silent event     an event that did not change that peak's contribution to the cost
  gross change     sum over steps and peaks of |change of the peak's term of the cost|, by type; the
                   denominator term (the ratio's N changes) is its own type
  net variation    sum over steps of |change of the cost|
  minimum          a local minimum of the cost along a ray, found with scipy.signal.find_peaks on
                   -cost with the given prominence, truth (r = 0) and the ray's end excluded
Minima walks: "dense" (full-sample: 0 to 5 deg; E0: 0 to 0.5 deg, both at 0.0005 deg),
"dense_to_0.5" (full-sample dense walk cut at 0.5 deg, the like-for-like comparison with E0) and
"profile" (E0, 121 log-spaced points: no width, the grid is not uniform). Ray fractions carry a
95% Wilson interval (scripts/common/stats.py).
"""

import argparse
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
from scipy.signal import find_peaks, peak_widths

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICE / "scripts" / "common"))
import cost_edges as CE  # noqa: E402
import cost_edges_run as RUN  # noqa: E402
import stats  # noqa: E402

CAUSES = RUN.CAUSES
REG = [r[0] for r in RUN.REGIMES]
REG_LABEL = {"fin": "< 0.05 deg", "mc": "0.05-0.5 deg", "grid": "0.5-5 deg"}
PROMS = (0.001, 0.01, 0.05)


def regime_of(rmid: np.ndarray) -> np.ndarray:
    out = np.zeros(len(rmid), np.int8)
    for i, (_, lo, hi) in enumerate(RUN.REGIMES):
        out[(rmid >= lo) & (rmid < hi)] = i
    return out


class Agg:
    def __init__(self) -> None:
        z = lambda *s: np.zeros(s)  # noqa: E731
        self.deg = z(3)
        self.steps = z(3)
        self.steps_changed = z(3)
        self.evn = z(3, len(CE.EVENT_TYPES))
        self.silent = z(3, len(CE.EVENT_TYPES))
        self.gross = z(3, len(CAUSES))
        self.net = z(3)
        self.flipc = z(3, len(CE.EVENT_TYPES) + 1)
        self.den_only = z(3)
        self.dn_steps = z(3)
        self.hist = z(3, len(CAUSES), len(RUN.HBINS))
        self.plateau = z(3)
        self.per_ray: List[np.ndarray] = []  # (3, 2): events/deg all and cost-changing, per ray
        self.rays = 0


def load(kind: str) -> Dict[Tuple[str, str], Agg]:
    aggs = {c: Agg() for c in RUN.COMBOS}
    files = sorted((RUN.CACHE / kind).glob("v*_d*.npz"))
    for f in files:
        z = np.load(f)
        for setting, rule in RUN.COMBOS:
            a = aggs[(setting, rule)]
            a.rays += 1
            ray_acc = np.zeros((3, 3))  # degrees, events, cost-changing events
            for seg in [k.split("|")[1] for k in z.files if k.startswith("radii|")]:
                key = f"{setting}|{rule}|parts|{seg}"
                if key not in z.files:
                    continue
                r = z[f"radii|{seg}"].astype(np.float64)
                rmid = 0.5 * (r[:-1] + r[1:])
                dr = np.diff(r)
                reg = regime_of(rmid)
                cost = z[f"{setting}|{rule}|cost|{seg}"].astype(np.float64)
                dcost = np.diff(cost)
                parts = z[key].astype(np.float64)
                evn = z[f"{setting}|{rule}|evn|{seg}"].astype(np.float64)
                sil = z[f"{setting}|{rule}|silent|{seg}"].astype(np.float64)
                fl = z[f"{setting}|{rule}|flipc|{seg}"].astype(np.float64)
                dn = z[f"{setting}|{rule}|dN|{seg}"]
                do = z[f"{setting}|{rule}|den_only|{seg}"]
                for g in range(3):
                    m = reg == g
                    a.deg[g] += dr[m].sum()
                    a.steps[g] += m.sum()
                    a.steps_changed[g] += (dcost[m] != 0).sum()
                    a.evn[g] += evn[m].sum(0)
                    a.silent[g] += sil[m].sum(0)
                    a.gross[g] += np.abs(parts[m]).sum(0)
                    a.net[g] += np.abs(dcost[m]).sum()
                    a.flipc[g] += fl[m].sum(0)
                    a.plateau[g] += ((cost[:-1] >= 0.9999) & (cost[1:] >= 0.9999))[m].sum()
                    ray_acc[g] += [dr[m].sum(), evn[m].sum(), (evn[m] - sil[m]).sum()]
                    a.den_only[g] += do[m].sum()
                    a.dn_steps[g] += (dn[m] != 0).sum()
                a.hist += z[f"{setting}|{rule}|hist_reg|{seg}"]
            with np.errstate(invalid="ignore", divide="ignore"):
                a.per_ray.append(
                    np.stack([ray_acc[:, 1] / ray_acc[:, 0], ray_acc[:, 2] / ray_acc[:, 0]], 1)
                )
    return aggs


def hist_quantiles(h: np.ndarray, qs: Tuple[float, ...]) -> List[float]:
    """Quantiles of the per-peak change from a histogram over log10 bins (bin lower edges HBINS)."""
    tot = h.sum()
    if tot == 0:
        return [float("nan")] * len(qs)
    cdf = np.cumsum(h) / tot
    out = []
    for q in qs:
        i = int(np.searchsorted(cdf, q))
        out.append(float(10 ** (RUN.HBINS[min(i, len(RUN.HBINS) - 1)] + 0.125)))
    return out


# --------------------------------------------------------------------------------------------
# minima along rays (step 4)
# --------------------------------------------------------------------------------------------
MINIMA_WALKS = {
    "full": [("dense", "dense", None), ("dense_to_0.5", "dense", 0.5)],
    "e0": [("dense", "dense", None), ("profile", "profile", None)],
}


def wilson(k: int, n: int) -> List[float]:
    lo, hi = stats.wilson(k, n)
    return [float(lo), float(hi)]


def minima_rays(kind: str) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    files = sorted((RUN.CACHE / kind).glob("v*_d*.npz"))
    for setting, rule in RUN.COMBOS:
        for name, seg, rmax in MINIMA_WALKS[kind]:
            rec: Dict[Any, Any] = {
                p: dict(
                    n_rays=0,
                    with_min=0,
                    n_min=0,
                    depth=[],
                    width=[],
                    rpos=[],
                    cause_in=[],
                    cause_out=[],
                )
                for p in PROMS
            }
            below_truth = 0
            n_rays = 0
            for f in files:
                z = np.load(f)
                c = z[f"{setting}|{rule}|cost|{seg}"].astype(np.float64)
                r = z[f"radii|{seg}"].astype(np.float64)
                parts = z.get(f"{setting}|{rule}|parts|{seg}")
                if rmax is not None:
                    keep = r <= rmax + 1e-9
                    c, r = c[keep], r[keep]
                    if parts is not None:
                        parts = parts[: int(keep.sum()) - 1]
                n_rays += 1
                below_truth += int(c.min() < c[0] - 1e-9)
                for p in PROMS:
                    rr = rec[p]
                    rr["n_rays"] += 1
                    pk, pr = find_peaks(-c, prominence=p)
                    pk = [k for k in pk if k > 0]
                    rr["with_min"] += int(len(pk) > 0)
                    rr["n_min"] += len(pk)
                    if not len(pk):
                        continue
                    w = peak_widths(-c, np.array(pk), rel_height=0.5)[0]
                    dr = np.median(np.diff(r))
                    for k, wd in zip(pk, w):
                        j = list(find_peaks(-c, prominence=p)[0]).index(k)
                        rr["depth"].append(float(pr["prominences"][j]))
                        if seg != "profile":  # a width in steps has no meaning on a log grid
                            rr["width"].append(float(wd * dr))
                        rr["rpos"].append(float(r[k]))
                        if parts is not None:
                            lo = k
                            while lo > 0 and c[lo - 1] >= c[lo]:
                                lo -= 1
                            hi = k
                            while hi < len(c) - 1 and c[hi + 1] >= c[hi]:
                                hi += 1
                            down = -np.minimum(parts[lo:k], 0).sum(
                                0
                            )  # cost lowered going in, by cause
                            up = np.maximum(parts[k:hi], 0).sum(0)  # cost raised going out
                            rr["cause_in"].append(
                                CAUSES[int(down.argmax())] if down.any() else "none"
                            )
                            rr["cause_out"].append(CAUSES[int(up.argmax())] if up.any() else "none")
            if (setting, rule, name) == ("full", "cpp", "dense"):
                nb: Dict[str, Any] = dict(rays=0, voxels=set(), cause=Counter(), argmin_deg=[])
                for f in files:
                    z = np.load(f)
                    c = z["full|cpp|cost|dense"].astype(float)
                    if c.min() < c[0] - 1e-9:
                        parts = z["full|cpp|parts|dense"].astype(float)
                        k = int(np.argmin(c))
                        nb["rays"] += 1
                        nb["voxels"].add(int(f.name.split("_")[0][1:]))
                        nb["cause"][CAUSES[int(np.minimum(parts[:k], 0).sum(0).argmin())]] += 1
                        nb["argmin_deg"].append(float(z["radii|dense"][k]))
                out["truth_not_lowest"] = dict(
                    rays=nb["rays"], voxels=sorted(nb["voxels"]), cause=dict(nb["cause"]),
                    max_argmin_deg=max(nb["argmin_deg"]) if nb["argmin_deg"] else None,
                )  # fmt: skip
            for v in rec.values():
                v["with_min_wilson95"] = wilson(v["with_min"], v["n_rays"])
            out[f"{setting}|{rule}|{name}"] = dict(
                n_rays=n_rays, rays_min_below_truth=below_truth,
                rays_min_below_truth_wilson95=wilson(below_truth, n_rays),
                by_prominence={str(p): v for p, v in rec.items()},
            )  # fmt: skip
    return out


def trace_check(kind: str) -> Dict[str, Any]:
    """Largest difference, over every evaluation of every walk, between the quality recombined from
    the per-peak trace and the library's own (check_trace, stored per task file as maxdiff)."""
    files = sorted((RUN.CACHE / kind).glob("v*_d*.npz"))
    md = [float(np.load(f)["maxdiff"]) for f in files]
    return dict(n_files=len(files), max_recombination_difference=max(md) if md else None)


def frame_share_per_voxel(kind: str) -> Dict[str, Any]:
    """Frame share of the gross change of the cost per voxel (pooled over its rays), complete cost,
    C++ rule, by regime: the per-voxel spread behind the pooled share."""
    files = sorted((RUN.CACHE / kind).glob("v*_d*.npz"))
    fi = CAUSES.index("frame")
    acc: Dict[int, np.ndarray] = {}
    for f in files:
        v = int(f.name.split("_")[0][1:])
        z = np.load(f)
        r = z["radii|dense"].astype(np.float64)
        reg = regime_of(0.5 * (r[:-1] + r[1:]))
        parts = np.abs(z["full|cpp|parts|dense"].astype(np.float64))
        g = acc.setdefault(v, np.zeros((3, 2)))
        for i in range(3):
            m = reg == i
            g[i] += [parts[m, fi].sum(), parts[m].sum()]
    out: Dict[str, Any] = {}
    for i, reg_name in enumerate(REG):
        sh = {str(v): float(g[i, 0] / g[i, 1]) if g[i, 1] else None for v, g in acc.items()}
        vals = [x for x in sh.values() if x is not None]
        out[reg_name] = dict(
            per_voxel=sh, min=min(vals) if vals else None, max=max(vals) if vals else None
        )
    return out


def nan_to_none(x: Any) -> Any:
    if isinstance(x, dict):
        return {k: nan_to_none(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [nan_to_none(v) for v in x]
    if isinstance(x, (float, np.floating)):
        return None if not np.isfinite(x) else float(x)
    if isinstance(x, np.integer):
        return int(x)
    return x


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out).resolve()
    summary: Dict[str, Any] = {}
    tables: Dict[str, str] = {}
    for kind in ("full", "e0"):
        if not list((RUN.CACHE / kind).glob("v*_d*.npz")):
            continue
        aggs = load(kind)
        summary[kind] = {}
        for (setting, rule), g in aggs.items():
            key = f"{setting}|{rule}"
            s: Dict[str, Any] = dict(n_rays=g.rays)
            for i, reg in enumerate(REG):
                deg = g.deg[i]
                s[reg] = dict(
                    degrees=deg,
                    steps=int(g.steps[i]),
                    steps_cost_changed=int(g.steps_changed[i]),
                    events_per_deg={t: g.evn[i, j] / deg for j, t in enumerate(CE.EVENT_TYPES)},
                    cost_changing_events_per_deg={
                        t: (g.evn[i, j] - g.silent[i, j]) / deg
                        for j, t in enumerate(CE.EVENT_TYPES)
                    },
                    events={t: int(g.evn[i, j]) for j, t in enumerate(CE.EVENT_TYPES)},
                    silent={t: int(g.silent[i, j]) for j, t in enumerate(CE.EVENT_TYPES)},
                    gross={c: float(g.gross[i, j]) for j, c in enumerate(CAUSES)},
                    net_variation=float(g.net[i]),
                    flips={
                        t: int(g.flipc[i, j])
                        for j, t in enumerate(CE.EVENT_TYPES + ["unattributed"])
                    },
                    plateau_steps=int(g.plateau[i]),
                    per_ray_events_per_deg_median_min_max=[
                        [float(np.nanmedian(x)), float(np.nanmin(x)), float(np.nanmax(x))]
                        for x in [np.array([r[i, k] for r in g.per_ray]) for k in (0, 1)]
                    ],
                    steps_denominator_changed=int(g.dn_steps[i]),
                    steps_denominator_only=int(g.den_only[i]),
                    jump_quantiles={
                        c: dict(
                            zip(
                                ("p50", "p90", "max_bin"),
                                hist_quantiles(g.hist[i, j], (0.5, 0.9, 1.0)),
                            ),
                            n=int(g.hist[i, j].sum()),
                        )
                        for j, c in enumerate(CAUSES)
                    },
                )
            summary[kind][key] = s
        summary[kind + "_minima"] = minima_rays(kind)
        summary[kind + "_trace_check"] = trace_check(kind)
        if kind == "full":
            summary["full_frame_share_per_voxel"] = frame_share_per_voxel(kind)
    text = json.dumps(nan_to_none(summary), default=float, separators=(",", ":"), allow_nan=False)
    (out / "summary.json").write_text(text + "\n")
    print("wrote", out / "summary.json")


if __name__ == "__main__":
    main()
