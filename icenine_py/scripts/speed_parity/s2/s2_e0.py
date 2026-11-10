"""S2 Step 6: the E0 A/B of the key FusedCost (E0 clean, 200 voxels, seed 0), from the per-voxel
records of ``scripts/findoptimal_robustness/e0_run.py``.

Arms (cache directories ``cache/e0<tag>``): ``off`` (_s2off, default keys), ``on`` (_s2on,
``--fused-cost``), ``bqmc`` (_s2bq, ``--batched-quick-mc``), ``bqmc_on`` (_s2bqon, both); and the
key-off arm of the base commit 832db36 (``--base``, a ``git archive`` copy run with its own
build), compared record by record with ``off`` (``s1_parity.compare-e0`` rule: every stored array
equal except the wall-clock ``*_runtime``). Wrong = reduced misorientation of the final orientation
to the truth > 1 degree. Paired by voxel: fixed = wrong under a, right under b; broken = the
reverse; ``same_final`` = the same final orientation (bit for bit).

  uv run python scripts/speed_parity/s2/s2_e0.py --base BASE/cache/e0_s2base \
      --out benchmarks/s2_results/e0.json
"""

import argparse
import json
import statistics
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
sys.path.insert(0, str(ICE / "scripts" / "findoptimal_robustness"))
sys.path.insert(0, str(ICE / "scripts" / "common"))
sys.path.insert(0, str(HERE.parent))

import common as C  # noqa: E402
import s1_parity as S1  # noqa: E402
import stats  # noqa: E402

ARMS = {
    "off": ("_s2off", "key off"),
    "on": ("_s2on", "FusedCost 1"),
    "bqmc": ("_s2bq", "BatchedQuickMC 1"),
    "bqmc_on": ("_s2bqon", "BatchedQuickMC 1, FusedCost 1"),
}
PAIRS = (
    ("off", "on", "key off vs FusedCost 1"),
    ("bqmc", "bqmc_on", "BatchedQuickMC 1 vs BatchedQuickMC 1 + FusedCost 1"),
    ("off", "bqmc", "key off vs BatchedQuickMC 1 (another draw, for scale)"),
)


def records(d: Path, vox: List[int], variant: str) -> Dict[int, Dict[str, Any]]:
    out: Dict[int, Dict[str, Any]] = {}
    for v in vox:
        p = d / f"v{v}_{variant}.npz"
        if not p.exists():
            continue
        x = np.load(p)
        if "unbuildable" in x.files:
            continue
        out[v] = dict(
            wrong=bool(float(C.err_deg(x["s0_R_final"], x["R_true"])) > C.WRONG_DEG),
            err=float(C.err_deg(x["s0_R_final"], x["R_true"])),
            R=x["s0_R_final"],
            g=int(x["s0_evals_global"]),
            loc=int(x["s0_evals_local"]),
            t=float(x["s0_runtime"]),
            fused=bool(x["fused_cost"]) if "fused_cost" in x.files else False,
            bqmc=bool(x["batched_quick_mc"]),
        )
    return out


def identical(da: Path, db: Path, names: List[str]) -> Dict[str, Any]:
    """compare-e0's rule over the given files."""
    n, same, keys_differ = 0, 0, {}
    for nm in names:
        if not (da / nm).exists() or not (db / nm).exists():
            continue
        x, y = np.load(da / nm), np.load(db / nm)
        n += 1
        bad = sorted(
            {
                k
                for k in x.files
                if "runtime" not in k and (k not in y.files or not np.array_equal(x[k], y[k]))
            }
            | {k for k in y.files if k not in x.files}
        )  # fmt: skip
        same += not bad
        for k in bad:
            keys_differ[k] = keys_differ.get(k, 0) + 1
    return dict(n=n, identical=same, keys_differ=keys_differ)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cache", default=str(C.CACHE_DIR))
    ap.add_argument("--base", required=True, help="the base commit's key-off cache directory")
    ap.add_argument("--out", required=True)
    ap.add_argument("--variant", default="clean")
    ap.add_argument("--n-voxels", type=int, default=200)
    ap.add_argument("--workers", type=int, required=True, help="the e0_run.py --workers used")
    a = ap.parse_args()
    cache = Path(a.cache).resolve()
    info = dict(np.load(C.OUT_DIR / "voxels.npz"))
    vox = [int(v) for v in info["voxel_indices"]][: a.n_voxels]
    recs = {arm: records(cache / f"e0{tag}", vox, a.variant) for arm, (tag, _) in ARMS.items()}
    common = sorted(set.intersection(*[set(r) for r in recs.values()]))
    res: Dict[str, Any] = dict(
        meta=dict(
            head=S1.git_head(),
            variant=a.variant,
            seed=0,
            n_voxels=a.n_voxels,
            n_common=len(common),
            workers=a.workers,
            arms={k: v[0] for k, v in ARMS.items()},
        ),
        arms={},
        pairs=[],
    )
    for arm, (tag, label) in ARMS.items():
        r = [recs[arm][v] for v in common]
        res["arms"][arm] = dict(
            label=label,
            n=len(r),
            wrong=sum(x["wrong"] for x in r),
            evals_local_mean=statistics.mean(x["loc"] for x in r),
            evals_global_mean=statistics.mean(x["g"] for x in r),
            evals_local_total=sum(x["loc"] for x in r),
            evals_global_total=sum(x["g"] for x in r),
            runtime_median=statistics.median(x["t"] for x in r),
            flags_ok=all(
                x["fused"] == arm.endswith("on") and x["bqmc"] == arm.startswith("bqmc") for x in r
            ),
            wrong_voxels=[v for v in common if recs[arm][v]["wrong"]],
        )  # fmt: skip
        lo, hi = stats.wilson(res["arms"][arm]["wrong"], len(r))
        res["arms"][arm]["wilson"] = [lo, hi]
    for x, y, label in PAIRS:
        wa = [recs[x][v]["wrong"] for v in common]
        wb = [recs[y][v]["wrong"] for v in common]
        fixed, broken = stats.paired_discordant(wa, wb)
        res["pairs"].append(
            dict(
                label=label,
                a=x,
                b=y,
                n=len(common),
                wrong_a=sum(wa),
                wrong_b=sum(wb),
                fixed=fixed,
                broken=broken,
                both_wrong=sum(p and q for p, q in zip(wa, wb)),
                both_right=sum(not p and not q for p, q in zip(wa, wb)),
                same_final=sum(
                    bool(np.array_equal(recs[x][v]["R"], recs[y][v]["R"])) for v in common
                ),
                same_evals=sum(
                    recs[x][v]["loc"] == recs[y][v]["loc"] and recs[x][v]["g"] == recs[y][v]["g"]
                    for v in common
                ),
                fixed_voxels=[v for v, p, q in zip(common, wa, wb) if p and not q],
                broken_voxels=[v for v, p, q in zip(common, wa, wb) if q and not p],
            )
        )
    base = Path(a.base).resolve()
    names = [f"v{v}_{a.variant}.npz" for v in vox]
    res["keyoff_vs_base"] = identical(cache / f"e0{ARMS['off'][0]}", base, names)
    bw = records(base, vox, a.variant)
    res["keyoff_vs_base"]["base_wrong"] = sum(x["wrong"] for x in bw.values())
    res["keyoff_vs_base"]["base_n"] = len(bw)
    out = Path(a.out).resolve()
    out.write_text(json.dumps(res, indent=1, default=float) + "\n")
    print(json.dumps({k: (v["wrong"], v["n"]) for k, v in res["arms"].items()}))
    for p in res["pairs"]:
        print(p["label"], {k: p[k] for k in ("fixed", "broken", "both_wrong", "same_final")})
    print("key off vs base:", {k: v for k, v in res["keyoff_vs_base"].items()})
    return 0


if __name__ == "__main__":
    sys.exit(main())
