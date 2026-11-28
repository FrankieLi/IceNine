"""AR5 step 17, phase A step 4: prediction 1 for the recording run ``J5r`` (and the stop rule's
stored counts).

  RESEED_SAMPLE=<big|grad|mid|held|held3> FITSKIP_CACHE=<dir> uv run python \
      scripts/autoresearch/ar5/fitskip_p1.py [--name J5r] [--stored J5]

Prediction 1: ``J5r`` gives ``J5``'s map before the pass at every voxel, its fit counts (the plan's
table) and its evaluations. Compares ``<cache>/<name>_<s>/`` (the recording run) with the stored
``<cache>/<stored>_<s>/`` (step 16's run, linked read only by ``fitskip_link_cache.sh``):

* every integer entry of ``run.json``'s ``stats`` (counts and evaluations; walls are floats and are
  not compared), the integer entries of ``polish.json`` and ``high500.json``;
* the map before the pass (``recon.mic.gz``) at every voxel (rotation matrices, to 1e-9) and the
  map after the pass (``passes_<s>/after_<name>.npy``), the pass's ``n_tries`` and ``evals``;
* the plan's fit counts (``fitskip_rules.counts_ok``: revisits, rejected revisits, neighbour fits,
  first-rejected neighbour fits).

Writes ``<cache>/p1_<name>_<s>.json`` (``ok``) and prints one line; exit 1 when it fails (the plan's
stop rule). No row of the recording is read here.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import fitskip_common as FCM  # noqa: E402  (first)

import fit_common as FC  # noqa: E402
import fitskip_rules as FR  # noqa: E402
import numpy as np  # noqa: E402
import reseed_common as RC  # noqa: E402

N_VOXELS = 24570
TOL = 1e-9


def ints(d: Dict[str, Any]) -> Dict[str, int]:
    return {k: int(v) for k, v in d.items() if isinstance(v, int) and not isinstance(v, bool)}


def diff_ints(a: Dict[str, Any], b: Dict[str, Any]) -> List[str]:
    x, y = ints(a), ints(b)
    out = [k for k in sorted(set(x) | set(y)) if x.get(k) != y.get(k)]
    return out


def main() -> None:
    args = sys.argv[1:]
    name = args[args.index("--name") + 1] if "--name" in args else "J5r"
    stored = args[args.index("--stored") + 1] if "--stored" in args else "J5"
    suffix = args[args.index("--suffix") + 1] if "--suffix" in args else ""
    s = RC.SAMPLE
    assert s in FR.EXPECT, "prediction 1 is checked on the five disclosed samples"
    import gradient_common as G

    c = FC.CACHE
    a_dir, b_dir = c / f"{name}_{s}{suffix}", c / f"{stored}_{s}{suffix}"
    ra = FCM.load_json(a_dir / "run.json")
    rb = FCM.load_json(b_dir / "run.json")
    res: Dict[str, Any] = dict(sample=s, name=name, stored=stored)
    res["stats_differ"] = diff_ints(ra["stats"], rb["stats"])
    for f in ("polish.json", "high500.json"):
        key = f.split(".")[0] + "_differ"
        pa, pb = FCM.load_json(a_dir / f), FCM.load_json(b_dir / f)
        res[key] = diff_ints(pa, pb)
    ma, _ = G.read_map(FCM.FG.guarded(a_dir / "recon.mic.gz"), N_VOXELS)
    mb, _ = G.read_map(FCM.FG.guarded(b_dir / "recon.mic.gz"), N_VOXELS)
    res["map_before_voxels_differ"] = int(
        (np.abs(ma.reshape(N_VOXELS, 9) - mb.reshape(N_VOXELS, 9)).max(axis=1) > TOL).sum()
    )
    pa_dir, pb_dir = c / f"passes_{s}", c / f"passes_{s}"
    if suffix:  # a smoke compare: both runs are the same short run, no pass, no table
        res["ok"] = bool(
            not res["stats_differ"]
            and not res["polish_differ"]
            and not res["high500_differ"]
            and res["map_before_voxels_differ"] == 0
        )
        (c / f"p1_{name}_{s}{suffix}.json").write_text(json.dumps(res, indent=1) + "\n")
        print(json.dumps(res))
        if not res["ok"]:
            raise SystemExit(1)
        return
    aa = np.load(FCM.FG.guarded(pa_dir / f"after_{name}.npy"))
    ab = np.load(FCM.FG.guarded(pb_dir / f"after_{stored}.npy"))
    res["map_after_voxels_differ"] = int(
        (np.abs(aa.reshape(N_VOXELS, 9) - ab.reshape(N_VOXELS, 9)).max(axis=1) > TOL).sum()
    )
    qa = FCM.load_json(pa_dir / f"D4_{name}" / "pass.json")
    qb = FCM.load_json(pb_dir / f"D4_{stored}" / "pass.json")
    res["pass_differ"] = [
        k for k in ("n_tries", "evals", "score_evals", "n_boundary_voxels") if qa[k] != qb[k]
    ]
    total = int(
        FCM.load_json(FR_STATS(s))["arms"]["J5"]["stages"]["total"]
    )  # the stored total with the pass (step 16's tracked stats)
    res["counts"] = FR.counts_ok(s, ra["stats"], total)
    res["ok"] = bool(
        not res["stats_differ"]
        and not res["polish_differ"]
        and not res["high500_differ"]
        and res["map_before_voxels_differ"] == 0
        and res["map_after_voxels_differ"] == 0
        and not res["pass_differ"]
        and res["counts"]["ok"]
    )
    (c / f"p1_{name}_{s}.json").write_text(json.dumps(res, indent=1) + "\n")
    print(
        json.dumps({k: v for k, v in res.items() if k != "counts"}),
        "counts_ok=",
        res["counts"]["ok"],
    )
    if not res["ok"]:
        raise SystemExit(1)


def FR_STATS(s: str) -> Path:
    """Step 16's tracked stats file (stored J5 totals with the pass)."""
    return FCM.FS.parent / "high500" / f"stats_{s}.json"


if __name__ == "__main__":
    main()
