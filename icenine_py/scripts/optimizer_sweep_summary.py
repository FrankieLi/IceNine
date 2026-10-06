"""Tables, JSON and plot for scripts/optimizer_sweep.py: the existing optimizers next to the
network on the perturbation sweep's cases (see that script's docstring for the protocol)."""

import json
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

import perturbation_sweep as ps

VARIANTS = ps.VARIANTS
# (label, method index, pass index) of the optimizer columns; GN passes: 0 = one-shot, 2 = x3
OPT_COLUMNS: List[Tuple[str, str, int]] = [
    ("MC", "mc", 0),
    ("Adam", "adam", 0),
    ("GN", "gn", 0),
    ("GN x3", "gn", 2),
    ("Huber", "huber", 0),
    ("Huber x3", "huber", 2),
]
PLOT_STYLE = {  # Okabe-Ito colours (colour-blind safe); labels are also in the legend
    "MC": ("#0072B2", "o", "-"),
    "Adam": ("#009E73", "s", "-"),
    "GN x3": ("#CC79A7", "^", "-"),
    "Huber x3": ("#E69F00", "v", "-"),
}
NET_COLOR = "#D55E00"


def _stats(err: np.ndarray, r: float) -> Dict[str, float]:
    e = err.astype(np.float64)
    if len(e) == 0:
        return {k: float("nan") for k in ("n", "median", "rms", "frac_lt_0p1", "frac_improved")}
    return dict(
        n=int(len(e)),
        median=float(np.median(e)),
        rms=float(np.sqrt(np.mean(e**2))),
        frac_lt_0p1=float(np.mean(e < 0.1)),
        frac_improved=float(np.mean(e < r)),
    )


def compute(
    opt: Dict[str, np.ndarray], sweep: Dict[str, np.ndarray], net_summary: Dict[str, Any]
) -> Dict[str, Any]:
    """Summary dict: per variant, per radius, per column statistics (optimizers, on the cases where
    the net could run at pass 1 and the optimizer produced an answer), the net's numbers from the
    sweep summary (seed means), and paired comparisons on identical cases."""
    radii = [float(r) for r in opt["radii"]]
    models = [str(m) for m in sweep["models"]]
    groups = ps.seed_groups(models)
    methods = [str(m) for m in opt["methods"]]
    net_ok = sweep["fail_pass1"] == 0  # (V, R, D, 2)
    out: Dict[str, Any] = {"radii": radii, "variants": {}}
    for vi, v in enumerate(VARIANTS):
        vo: Dict[str, Any] = {}
        for ri, r in enumerate(radii):
            key = f"{r:g}"
            row: Dict[str, Any] = {}
            for label, m, p in OPT_COLUMNS:
                mi = methods.index(m)
                err = opt["err_angle"][:, ri, :, vi, mi, p]
                ran = net_ok[:, ri, :, vi] & np.isfinite(err)
                s = _stats(err[ran], r)
                rt = opt["runtime"][:, ri, :, vi, mi, p][ran]
                s["median_runtime_s"] = float(np.median(rt)) if len(rt) else float("nan")
                row[label] = s
            for tname, idx in groups.items():
                for p, tag in ((0, "one-shot"), (2, "x3")):
                    ms = [net_summary["metrics"]["all"][v][models[i]][key][p] for i in idx]
                    row[f"net {tname} {tag}"] = dict(
                        n=int(np.mean([m["n_ok"] for m in ms])),
                        median=float(np.mean([m["median_angle"] for m in ms])),
                        rms=float(np.mean([m["rms_angle"] for m in ms])),
                        frac_lt_0p1=float(np.mean([m["frac_lt_0p1"] for m in ms])),
                        frac_improved=float(np.mean([m["frac_improved"] for m in ms])),
                    )
            # paired: fraction of identical cases where the net's error is smaller
            paired: Dict[str, Any] = {}
            for tname, idx in groups.items():
                for p, tag in ((0, "one-shot"), (2, "x3")):
                    for label, m, op in OPT_COLUMNS:
                        mi = methods.index(m)
                        oe = opt["err_angle"][:, ri, :, vi, mi, op]
                        wins, ns = [], []
                        for i in idx:
                            ne = sweep["err_angle"][:, ri, :, vi, i, p]
                            ok = net_ok[:, ri, :, vi] & np.isfinite(oe) & np.isfinite(ne)
                            ns.append(int(ok.sum()))
                            wins.append(float(np.mean(ne[ok] < oe[ok])) if ok.any() else np.nan)
                        paired[f"net {tname} {tag} < {label}"] = dict(
                            n=ns[0], win_per_seed=wins, win_mean=float(np.mean(wins))
                        )
            row["paired"] = paired
            vo[key] = row
        out["variants"][v] = vo
    return out


def _fmt(x: float, f: str = ".4f") -> str:
    return "   nan" if x != x else format(x, f)


def format_tables(summ: Dict[str, Any], opt: Dict[str, np.ndarray]) -> str:
    radii = [f"{r:g}" for r in summ["radii"]]
    cfg = json.loads(str(opt["config_json"]))
    n_vox, n_dir = opt["err_angle"].shape[0], opt["err_angle"].shape[2]
    lines = [
        f"Existing optimizers vs the network on the perturbation sweep: {n_vox} voxels x "
        f"{len(radii)} radii x {n_dir} directions; MC / Adam on the first {cfg['opt_dirs']} "
        "directions per radius, GN / Huber GN on all.",
        "Error = angle (deg) between the result and the true orientation, over the cases where the "
        "net could run at pass 1 (the same cases for every column). Net columns: seed means "
        "(2 seeds) from perturbation_sweep_summary.json. 'x3' = iterated with re-centring.",
        "MC is told r (search box 1.5 r); the net is not.",
        "",
    ]
    cols = [c[0] for c in OPT_COLUMNS]
    nets = [k for k in next(iter(summ["variants"]["clean"].values())) if k.startswith("net ")]
    allc = cols + nets
    blocks = [
        ("median error (deg)", "median", ".4f"),
        ("RMS error (deg)", "rms", ".4f"),
        ("fraction of cases < 0.1 deg", "frac_lt_0p1", ".3f"),
        ("fraction improved (error < r)", "frac_improved", ".3f"),
    ]
    for v in VARIANTS:
        lines.append(f"=== {v} data ===")
        vo = summ["variants"][v]
        for title, key, f in blocks:
            lines.append(f"-- {title}")
            lines.append(f"{'r':>5} {'n':>5} " + " ".join(f"{c:>20}" for c in allc))
            for r in radii:
                row = vo[r]
                lines.append(
                    f"{r:>5} {row['MC']['n']:>5} "
                    + " ".join(f"{_fmt(row[c][key], f):>20}" for c in allc)
                )
            lines.append("")
        lines.append(
            "-- median runtime per case (s), data preparation excluded (GN x3: cumulative)"
        )
        lines.append(f"{'r':>5} " + " ".join(f"{c:>9}" for c in cols))
        for r in radii:
            lines.append(
                f"{r:>5} "
                + " ".join(f"{_fmt(vo[r][c]['median_runtime_s'], '.3f'):>9}" for c in cols)
            )
        lines.append("")
        for tname in ("realistic", "clean"):
            for tag in ("x3", "one-shot"):
                lines.append(
                    f"-- paired: fraction of identical cases where the {tname}-trained net "
                    f"({tag}) has the smaller error (seed mean; per seed in the JSON)"
                )
                lines.append(f"{'r':>5} " + " ".join(f"{c:>9}" for c in cols))
                for r in radii:
                    pr = vo[r]["paired"]
                    lines.append(
                        f"{r:>5} "
                        + " ".join(
                            f"{_fmt(pr[f'net {tname} {tag} < {c}']['win_mean'], '.3f'):>9}"
                            for c in cols
                        )
                    )
                lines.append("")
    return "\n".join(lines)


def plot(summ: Dict[str, Any], path: Path, n_vox: int, opt_dirs: int) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    radii = np.array(summ["radii"])
    keys = [f"{r:g}" for r in radii]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.9), sharey=True)
    grid, ink = "#d9d9d9", "#333333"
    for ax, v in zip(axes, VARIANTS):
        vo = summ["variants"][v]
        ax.plot(radii, radii, ":", color="0.45", lw=1.2, label="y = x (no correction)")
        ax.axvline(1.0, color="0.7", lw=1.2)
        for label, (c, mk, ls) in PLOT_STYLE.items():
            y = [vo[k][label]["median"] for k in keys]
            ax.plot(radii, y, ls, color=c, marker=mk, ms=4.5, lw=1.8, label=label)
        for tname, a_, lw in (("clean", 0.3, 1.4), ("realistic", 1.0, 2.2)):
            for tag, ls in (("x3", "-"), ("one-shot", "--")):
                if tname == "clean" and tag == "one-shot":
                    continue
                y = [vo[k][f"net {tname} {tag}"]["median"] for k in keys]
                lab = f"net, {tname}-trained, " + ("iterated x3" if tag == "x3" else "one-shot")
                ax.plot(
                    radii, y, ls, color=NET_COLOR, alpha=a_, lw=lw, marker="o", ms=3.5, label=lab
                )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("perturbation radius r (deg)")
        ax.set_title("clean data" if v == "clean" else "realistic data (all)", color=ink)
        ax.grid(True, which="both", color=grid, lw=0.6)
        ax.set_axisbelow(True)
    axes[0].set_ylabel("median angular error (deg)")
    axes[0].legend(fontsize=7.5, loc="upper left", frameon=False)
    fig.suptitle(
        f"Existing optimizers vs the network: {n_vox} voxels, 20 directions per radius "
        f"(MC/Adam: {opt_dirs}); grey line = end of the net's training prior (1 deg); "
        "MC is told r",
        fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def do_summarize(out_dir: Path) -> None:
    opt = dict(np.load(out_dir / "optimizer_sweep_raw.npz", allow_pickle=False))
    sweep = dict(np.load(out_dir / "perturbation_sweep_raw.npz", allow_pickle=False))
    net_summary = json.loads((out_dir / "perturbation_sweep_summary.json").read_text())
    summ = compute(opt, sweep, net_summary)
    cfg = json.loads(str(opt["config_json"]))
    # image-based cost at the truth and at the end (diagnostics)
    diag: Dict[str, Any] = {}
    for vi, v in enumerate(VARIANTS):
        for k, name in enumerate(("mc", "adam")):
            qt = opt["q_true"][..., vi, k]
            qt = qt[np.isfinite(qt)]
            diag[f"{v}_{name}_quality_at_truth_median"] = float(np.median(qt)) if len(qt) else None
    summ["diagnostics"] = diag
    (out_dir / "optimizer_sweep_summary.json").write_text(json.dumps(summ, indent=1))
    text = format_tables(summ, opt)
    text += "\n\nimage cost at the true orientation (median over cases): " + json.dumps(diag)
    (out_dir / "optimizer_sweep_summary.txt").write_text(text + "\n")
    plot(
        summ,
        out_dir / "perturbation_sweep_vs_optimizers.png",
        opt["err_angle"].shape[0],
        cfg["opt_dirs"],
    )
    print(text)
