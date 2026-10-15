"""Task E measurement 1, summary: the cost along Task B's rays for the strict cost and each arm.

Per ray and cost (strict = Task B's complete cost under the C++ rule and the lenient rule; the
arms of ``smooth_landscape.ARMS``):

* steps per degree: a step between two neighbouring points (0.0005 deg apart) is an isolated step
  when its change exceeds 1e-6 and three times the change of both neighbouring steps (a jump, not
  a slope). Counted per regime (under 0.05 deg, 0.05 to 0.5 deg, 0.5 to 5 deg; on E0 the last is
  four 0.05-deg windows).
* share of the variation in steps: sum of |change| over the isolated steps / sum over all steps.
* type of the isolated steps: Task B's per-peak events at the same step of the same ray (its cache;
  priority presence > frame > detector > pixel; "none" = no strict event at that step).
* local minima (not the truth) by prominence 0.01 and 0.05 (``scipy.signal.find_peaks``), rays with
  one, with a 95% Wilson interval; rays where the lowest cost is not at the truth.
* basin radius at tolerance eps: the largest radius r such that the cost never falls more than eps
  below its running maximum on [0, r] (the cost is eps-monotone outward); full-sample: the dense
  walk to 5 deg; E0: the log profile to 5 deg (121 points).

  uv run python scripts/eval_reduction/smooth_landscape_summary.py
Writes benchmarks/eval_reduction/E/landscape/{summary.json, tables.md} and the figures.
"""

import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
from scipy.signal import find_peaks

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICE / "scripts" / "common"))
import cost_edges as CE  # noqa: E402
import doc_tables as DT  # noqa: E402
import stats as ST  # noqa: E402
from smooth_landscape import ARMS, CACHE  # noqa: E402

BCACHE = ICE / "benchmarks" / "eval_reduction" / "B" / "cache"
OUT = ICE / "benchmarks" / "eval_reduction" / "E" / "landscape"
FIG = ICE / "benchmarks" / "eval_reduction" / "E" / "figures"
REGIMES = [
    ("under 0.05 deg", 0.0, 0.05),
    ("0.05 to 0.5 deg", 0.05, 0.5),
    ("0.5 to 5 deg", 0.5, 5.01),
]
PROMS = (0.01, 0.05)
EPS = (0.01, 0.05)
COSTS = ["strict_cpp", "strict_old"] + list(ARMS)
LABEL = {
    "strict_cpp": "strict cost, C++ rule",
    "strict_old": "strict cost, lenient rule",
    "S": "S: sampled binary (control)",
    "E1s3": "E1: distance kernel 3 px",
    "E1s12": "E1: distance kernel 12 px",
    "E2": "E2: omega blending, 1 frame",
    "E3": "E3: soft peak weights",
    "E4t05": "E4: ray residual, tau 0.05 deg",
    "E4t25": "E4: ray residual, tau 0.25 deg",
    "E12": "E1 + E2",
    "E123": "E1 + E2 + E3",
    "E23": "E2 + E3",
    "E34": "E3 + E4 (0.25 deg)",
}
TYPES = ["frame", "pixel", "presence", "detector", "none"]


def isolated_steps(c: np.ndarray) -> np.ndarray:
    d = np.abs(np.diff(c.astype(np.float64)))
    left = np.concatenate([[0.0], d[:-1]])
    right = np.concatenate([d[1:], [0.0]])
    return (d > 1e-6) & (d > 3.0 * np.maximum(left, right))


def step_types(evn: np.ndarray) -> np.ndarray:
    """Label of each step from Task B's per-peak event counts (all events, silent ones too)."""
    idx = {t: i for i, t in enumerate(CE.EVENT_TYPES)}
    pres = evn[:, [idx["bragg_limit"], idx["eta_cut"], idx["no_frame"]]].sum(1) > 0
    lab = np.full(len(evn), "none", dtype=object)
    lab[evn[:, idx["pixel"]] > 0] = "pixel"
    lab[evn[:, idx["detector"]] > 0] = "detector"
    lab[evn[:, idx["frame"]] > 0] = "frame"
    lab[pres] = "presence"
    return lab


def basin_radius(r: np.ndarray, c: np.ndarray, eps: float) -> Tuple[float, bool]:
    viol = np.maximum.accumulate(c) - c > eps
    if not viol.any():
        return float(r[-1]), True
    k = int(np.argmax(viol))
    return float(r[max(k - 1, 0)]), False


def ray_costs(kind: str, f: Path) -> Dict[str, Dict[str, np.ndarray]]:
    """cost -> segment -> array, with radii under 'radii'."""
    z = np.load(f)
    b = np.load(BCACHE / kind / f.name)
    segs = [k.split("|")[1] for k in z.files if k.startswith("radii|")]
    out: Dict[str, Dict[str, np.ndarray]] = {"radii": {}, "evn": {}}
    for s in segs:
        out["radii"][s] = z[f"radii|{s}"].astype(np.float64)
        rb = b[f"radii|{s}"].astype(np.float64)
        assert np.allclose(rb, out["radii"][s], atol=1e-6), (f, s)
        if f"full|cpp|evn|{s}" in b.files:
            out["evn"][s] = b[f"full|cpp|evn|{s}"]
    for cname in COSTS:
        out[cname] = {}
        for s in segs:
            if cname == "strict_cpp":
                out[cname][s] = b[f"full|cpp|cost|{s}"].astype(np.float64)
            elif cname == "strict_old":
                out[cname][s] = b[f"full|old|cost|{s}"].astype(np.float64)
            else:
                out[cname][s] = z[f"{cname}|{s}"].astype(np.float64)
    return out


def summarise(kind: str) -> Dict[str, Any]:
    files = sorted((CACHE / kind).glob("v*_d*.npz"))
    acc: Dict[str, Any] = {
        c: dict(
            steps=np.zeros(3),
            deg=np.zeros(3),
            var_step=np.zeros(3),
            var_all=np.zeros(3),
            types=Counter(),
            minima={str(p): dict(rays=0, n=0) for p in PROMS},
            minima05={str(p): dict(rays=0, n=0) for p in PROMS},
            not_lowest=0,
            argmin_deg=[],
            below_truth=[],
            basin={str(e): [] for e in EPS},
            basin_full={str(e): 0 for e in EPS},
            cost_truth=[],
        )
        for c in COSTS
    }
    n_rays = 0
    attached = []
    for f in files:
        R = ray_costs(kind, f)
        n_rays += 1
        attached.append(int(np.load(f)["n_attached"]))
        for cname in COSTS:
            a = acc[cname]
            for s, r in R["radii"].items():
                if s == "profile":
                    continue
                c = R[cname][s]
                iso = isolated_steps(c)
                d = np.abs(np.diff(c))
                mid = 0.5 * (r[:-1] + r[1:])
                for i, (_, lo, hi) in enumerate(REGIMES):
                    m = (mid >= lo) & (mid < hi)
                    a["steps"][i] += int((iso & m).sum())
                    a["deg"][i] += float(np.diff(r)[m].sum())
                    a["var_step"][i] += float(d[iso & m].sum())
                    a["var_all"][i] += float(d[m].sum())
                if s in R["evn"]:
                    lab = step_types(R["evn"][s])
                    a["types"].update(lab[iso].tolist())
            c = R[cname]["dense"]
            r = R["radii"]["dense"]
            a["cost_truth"].append(float(c[0]))
            a["not_lowest"] += int(c.min() < c[0] - 1e-9)
            if c.min() < c[0] - 1e-9:
                a["argmin_deg"].append(float(r[int(np.argmin(c))]))
                a["below_truth"].append(float(c[0] - c.min()))
            for p in PROMS:
                pk, _ = find_peaks(-c, prominence=p)
                pk = pk[pk > 0]
                a["minima"][str(p)]["rays"] += int(len(pk) > 0)
                a["minima"][str(p)]["n"] += int(len(pk))
                k05 = r <= 0.5 + 1e-9
                pk5, _ = find_peaks(-c[k05], prominence=p)
                pk5 = pk5[pk5 > 0]
                a["minima05"][str(p)]["rays"] += int(len(pk5) > 0)
                a["minima05"][str(p)]["n"] += int(len(pk5))
            if kind == "full":
                rb, cb = r, c
            else:
                rp = R["radii"]["profile"]
                rb = np.concatenate([[0.0], rp])
                cb = np.concatenate([[c[0]], R[cname]["profile"]])
            for e in EPS:
                br, full = basin_radius(rb, cb, e)
                a["basin"][str(e)].append(br)
                a["basin_full"][str(e)] += int(full)
    res: Dict[str, Any] = dict(n_rays=n_rays, attached_median=float(np.median(attached)))
    for cname, a in acc.items():
        tot = sum(a["types"].values())
        res[cname] = dict(
            steps_per_deg=[
                float(x / y) if y else float("nan") for x, y in zip(a["steps"], a["deg"])
            ],
            share_variation_in_steps=[
                float(x / y) if y else float("nan") for x, y in zip(a["var_step"], a["var_all"])
            ],
            deg_walked=a["deg"].tolist(),
            step_types={t: int(a["types"].get(t, 0)) for t in TYPES},
            step_types_share={
                t: (a["types"].get(t, 0) / tot if tot else float("nan")) for t in TYPES
            },
            minima={
                p: dict(v, wilson95=list(ST.wilson(v["rays"], n_rays)))
                for p, v in a["minima"].items()
            },
            minima_to_05={
                p: dict(v, wilson95=list(ST.wilson(v["rays"], n_rays)))
                for p, v in a["minima05"].items()
            },
            truth_not_lowest=int(a["not_lowest"]),
            lowest_point_median_deg=(
                float(np.median(a["argmin_deg"])) if a["argmin_deg"] else None
            ),
            lowest_point_max_deg=float(np.max(a["argmin_deg"])) if a["argmin_deg"] else None,
            lowest_below_truth_median=(
                float(np.median(a["below_truth"])) if a["below_truth"] else None
            ),
            basin_median_deg={e: float(np.median(v)) for e, v in a["basin"].items()},
            basin_q25_deg={e: float(np.quantile(v, 0.25)) for e, v in a["basin"].items()},
            basin_to_5deg=a["basin_full"],
            median_cost_at_truth=float(np.median(a["cost_truth"])),
        )
    return res


def pct(x: float) -> str:
    return "n/a" if x != x else f"{100.0 * x:.1f}%"


def wil(v: Dict[str, Any], n: int) -> str:
    lo, hi = v["wilson95"]
    return f"{v['rays']} of {n} ({100 * lo:.1f}-{100 * hi:.1f}%)"


def tables(S: Dict[str, Any]) -> Dict[str, str]:
    T: Dict[str, str] = {}
    for kind in ("full", "e0"):
        R = S[kind]
        n = R["n_rays"]
        g = "0.5 to 5 deg" if kind == "full" else "0.5 to 5 deg (4 windows of 0.05 deg)"
        rows = []
        for c in COSTS:
            x = R[c]
            rows.append(
                {
                    "cost": LABEL[c],
                    "cost at the truth (median)": x["median_cost_at_truth"],
                    "steps per deg, under 0.05 deg": x["steps_per_deg"][0],
                    "steps per deg, 0.05 to 0.5 deg": x["steps_per_deg"][1],
                    f"steps per deg, {g}": x["steps_per_deg"][2],
                    "variation in steps, under 0.05 deg": pct(x["share_variation_in_steps"][0]),
                    "variation in steps, 0.05 to 0.5 deg": pct(x["share_variation_in_steps"][1]),
                    f"variation in steps, {g}": pct(x["share_variation_in_steps"][2]),
                }
            )
        T[f"evalE_land_steps_{kind}"] = DT.markdown_table(
            rows,
            list(rows[0].keys()),
            {k: ".1f" for k in rows[0] if k.startswith("steps")}
            | {"cost at the truth (median)": ".4f"},
        )
        rows = []
        for c in COSTS:
            x = R[c]
            tot = sum(x["step_types"].values())
            row: Dict[str, Any] = {"cost": LABEL[c], "isolated steps": tot}
            for t in TYPES:
                row[f"at a {t} event" if t != "none" else "at no strict event"] = pct(
                    x["step_types_share"][t]
                )
            rows.append(row)
        T[f"evalE_land_types_{kind}"] = DT.markdown_table(rows, list(rows[0].keys()), {})
        rows = []
        walk = "dense, 0 to 5 deg" if kind == "full" else "dense, 0 to 0.5 deg"
        for c in COSTS:
            x = R[c]
            row = {"cost": LABEL[c]}
            for p in ("0.01", "0.05"):
                row[f"rays with a minimum of depth {p}, {walk}"] = wil(x["minima"][p], n)
            if kind == "full":
                for p in ("0.01", "0.05"):
                    row[f"rays with a minimum of depth {p}, 0 to 0.5 deg"] = wil(
                        x["minima_to_05"][p], n
                    )
            row["rays where the truth is not the lowest point"] = f"{x['truth_not_lowest']} of {n}"
            lp = x["lowest_point_median_deg"]
            row["lowest point when not the truth: median deg (max)"] = (
                "" if lp is None else f"{lp:.4f} ({x['lowest_point_max_deg']:.4f})"
            )
            lb = x["lowest_below_truth_median"]
            row["... and its depth below the truth (median)"] = "" if lb is None else f"{lb:.4f}"
            rows.append(row)
        T[f"evalE_land_minima_{kind}"] = DT.markdown_table(rows, list(rows[0].keys()), {})
        rows = []
        for c in COSTS:
            x = R[c]
            rows.append(
                {
                    "cost": LABEL[c],
                    "basin radius, eps 0.01 (median deg)": x["basin_median_deg"]["0.01"],
                    "basin radius, eps 0.01 (lower quartile deg)": x["basin_q25_deg"]["0.01"],
                    "rays monotone to 5 deg, eps 0.01": f"{x['basin_to_5deg']['0.01']} of {n}",
                    "basin radius, eps 0.05 (median deg)": x["basin_median_deg"]["0.05"],
                    "rays monotone to 5 deg, eps 0.05": f"{x['basin_to_5deg']['0.05']} of {n}",
                }
            )
        T[f"evalE_land_basin_{kind}"] = DT.markdown_table(
            rows,
            list(rows[0].keys()),
            {k: ".3f" for k in rows[0] if "median" in k or "quartile" in k},
        )
    return T


def figures() -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIG.mkdir(parents=True, exist_ok=True)
    f = CACHE / "full" / "v15901_d0.npz"
    R = ray_costs("full", f)
    r = R["radii"]["dense"]
    show = [
        "strict_cpp",
        "S",
        "E1s3",
        "E1s12",
        "E2",
        "E3",
        "E4t05",
        "E4t25",
        "E12",
        "E123",
        "E23",
        "E34",
    ]
    cmap = plt.get_cmap("tab20")
    fig, axs = plt.subplots(1, 3, figsize=(16, 5.2))
    for ax, hi in zip(axs, (0.1, 0.5, 5.0)):
        m = r <= hi + 1e-9
        for i, c in enumerate(show):
            st = dict(color="k", lw=1.6) if c == "strict_cpp" else dict(color=cmap(i), lw=0.9)
            ax.plot(r[m], R[c]["dense"][m], label=LABEL[c], **st)
        ax.set_xlabel("rotation away from the true orientation (deg)")
        ax.set_title(f"voxel 15901, one ray (Task B's ray plot): 0 to {hi} deg", fontsize=9)
        ax.set_ylim(-0.02, 1.02)
    axs[0].set_ylabel("cost (lower = better)")
    axs[2].legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG / "ray_strict_vs_arms.png", dpi=140)
    plt.close(fig)
    # small multiples: each arm against the strict cost, 0 to 0.5 deg
    arms = show[1:]
    fig, axs = plt.subplots(3, 4, figsize=(15, 9), sharex=True, sharey=True)
    m = r <= 0.5 + 1e-9
    for ax, c in zip(axs.ravel(), arms):
        ax.plot(r[m], R["strict_cpp"]["dense"][m], color="k", lw=1.2, label="strict, C++ rule")
        ax.plot(r[m], R[c]["dense"][m], color="tab:red", lw=1.0, label=LABEL[c])
        ax.set_title(LABEL[c], fontsize=9)
    for ax in axs.ravel()[len(arms) :]:
        ax.axis("off")
    for ax in axs[-1]:
        ax.set_xlabel("deg from the truth")
    axs[0, 0].legend(fontsize=7)
    fig.suptitle("voxel 15901, Task B's ray 0: each arm (red) against the strict cost (black)")
    fig.tight_layout()
    fig.savefig(FIG / "ray_each_arm.png", dpi=130)
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    S = {k: summarise(k) for k in ("full", "e0")}
    (OUT / "summary.json").write_text(json.dumps(S, indent=1, default=float) + "\n")
    T = tables(S)
    DT.write_tables(OUT / "tables.md", T)
    figures()
    for k, t in T.items():
        print(k)
        print(t)
        print()


if __name__ == "__main__":
    main()
