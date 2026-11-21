"""AR5 step 8, Phase A step 4: are two runs identical, and does ``Al`` reproduce step 3's ``S11``?

  uv run python scripts/autoresearch/ar5/order_compare.py --a chk1 --b chk2 [--stored-al chkAl]

Compares two run folders of ``order_run.py`` on ``big`` (default ``chk1``, ``chk2``: arm ``A`` twice
with aligned draws): the decompressed ``recon.mic`` (the map), the seed-step log with every fit
(``e1_log.json.gz``), the MC traces, and the evaluation counts of ``run.json``. With
``--stored-al`` it also compares the run ``Al`` (interleaved order, the library's own draws) with
the stored ``S11`` of step 3 on ``big``: the map, the seed-step rows (voxel, evaluations, hit
ratios) and the total before the pass (2,879,705 in step 3). Writes ``order/check_runs.json``.
"""

import argparse
import gzip
import json
from pathlib import Path
from typing import Any, Dict

import order_common as OC  # noqa: E402  (first)


def read(p: Path) -> bytes:
    with gzip.open(p, "rb") as f:
        return f.read()


def jgz(p: Path) -> Any:
    return json.loads(read(p))


def totals(run: Dict[str, Any]) -> Dict[str, Any]:
    st = run["stats"]
    keys = ("n_evals_seed", "n_evals_neighbor", "n_evals_revisit", "n_evals_reseed")
    return dict(
        evals_total=int(run["evals_total"]), **{k: int(st.get(k, 0)) for k in keys if k in st}
    )


def same_runs(a: Path, b: Path, with_fits: bool = True) -> Dict[str, Any]:
    ra, rb = json.loads((a / "run.json").read_text()), json.loads((b / "run.json").read_text())
    out: Dict[str, Any] = dict(
        map_identical=read(a / "recon.mic.gz") == read(b / "recon.mic.gz"),
        totals_a=totals(ra),
        totals_b=totals(rb),
        totals_identical=totals(ra) == totals(rb),
    )
    ea, eb = jgz(a / "e1_log.json.gz"), jgz(b / "e1_log.json.gz")
    if not with_fits:  # stored S11 rows have no per-fit records and no 'order' field
        key = ("voxel", "reseed", "n_cand", "stop", "fin_evals", "fin_hit", "costed", "evals")
        ea = [{k: r.get(k) for k in key} for r in ea]
        eb = [{k: r.get(k) for k in key} for r in eb]
    out["seed_log_identical"] = ea == eb
    out["n_seed_steps"] = [len(ea), len(eb)]
    if with_fits:
        out["traces_identical"] = jgz(a / "mc_traces.json.gz") == jgz(b / "mc_traces.json.gz")
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", default="chk1")
    ap.add_argument("--b", default="chk2")
    ap.add_argument("--stored-al", default="")
    a = ap.parse_args()
    res: Dict[str, Any] = {}
    res["twice"] = same_runs(OC.ORDER_CACHE / f"{a.a}_big", OC.ORDER_CACHE / f"{a.b}_big")
    pa, pb = (OC.ORDER_CACHE / "passes_big" / f"after_{x}.npy" for x in (a.a, a.b))
    if pa.exists() and pb.exists():
        import numpy as np

        ja = json.loads((pa.parent / f"D4_{a.a}" / "pass.json").read_text())
        jb = json.loads((pb.parent / f"D4_{a.b}" / "pass.json").read_text())
        res["twice"]["pass_map_identical"] = bool(np.array_equal(np.load(pa), np.load(pb)))
        res["twice"]["pass_evals"] = [ja["evals"], jb["evals"]]
        res["twice"]["pass_tries"] = [ja["n_tries"], jb["n_tries"]]
    if a.stored_al:
        res["al_vs_stored_S11"] = same_runs(
            OC.ORDER_CACHE / f"{a.stored_al}_big", OC.STORED_S11["big"], with_fits=False
        )
    (OC.ORDER / "check_runs.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
