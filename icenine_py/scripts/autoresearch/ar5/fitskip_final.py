"""AR5 step 17, phase A: which recorded fit end is the voxel's final answer, before and after the
pass (the last field of a recording row).

  RESEED_SAMPLE=<big|grad|mid|held|held3> FITSKIP_CACHE=<dir> uv run python \
      scripts/autoresearch/ar5/fitskip_final.py [--name J5r]

The final answer of a voxel before the pass is the orientation in the run's map
(``<name>_<s>/recon.mic.gz``); after the pass, the orientation in ``passes_<s>/after_<name>.npy``.
A fit's end is the voxel's final answer when it is the voxel's LAST fit that the push rule kept
(``kept``) and its end rotation matrix equals the map's to ``TOL`` (max absolute difference of the
nine entries). A voxel whose map orientation matches no kept fit end (a full search, a rescue, a
grid search, a pass try) has no fit that is final. Writes ``<cache>/record/final_<name>_<s>.json``:
the row numbers ``i`` that are final before and after the pass, the counts, and how many voxels with
a kept fit matched (the rest were replaced later). ``fitskip_rules.replay`` reads the row field
``final_before`` that ``merge_final`` sets from this file. No other row field is read here.
"""

import gzip
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Set

import fitskip_common as FCM  # noqa: E402  (first)

import fit_common as FC  # noqa: E402
import numpy as np  # noqa: E402
import reseed_common as RC  # noqa: E402

N_VOXELS = 24570
TOL = 1e-6


def read_rows(path: Path) -> List[Dict[str, Any]]:
    op = gzip.open if str(path).endswith(".gz") else open
    with op(FCM.FG.guarded(path), "rt") as f:  # type: ignore[operator]
        return [json.loads(x) for x in f if x.strip()]


def last_kept(rows: List[Dict[str, Any]]) -> Dict[int, Dict[str, Any]]:
    """The last row with ``kept`` per voxel (run order)."""
    out: Dict[int, Dict[str, Any]] = {}
    for r in rows:
        if r["kept"]:
            out[int(r["voxel"])] = r
    return out


def final_rows(rows: List[Dict[str, Any]], R_map: np.ndarray, tol: float = TOL) -> Set[int]:
    """The row numbers ``i`` whose end equals the map's orientation at their voxel, among the
    last kept row of each voxel."""
    out: Set[int] = set()
    for v, r in last_kept(rows).items():
        d = float(np.abs(np.asarray(r["end_R"], np.float64) - R_map[v].reshape(9)).max())
        if d <= tol:
            out.add(int(r["i"]))
    return out


def merge_final(rows: List[Dict[str, Any]], fin: Dict[str, Any]) -> None:
    """Set ``final_before`` and ``final_after`` on every row (in place)."""
    b, a = set(fin["before"]), set(fin["after"])
    for r in rows:
        r["final_before"] = int(r["i"]) in b
        r["final_after"] = int(r["i"]) in a


def main() -> None:
    args = sys.argv[1:]
    name = args[args.index("--name") + 1] if "--name" in args else "J5r"
    s = RC.SAMPLE
    import gradient_common as G

    c = FC.CACHE
    rows = read_rows(c / "record" / f"{name}_{s}.jsonl.gz")
    mb, _ = G.read_map(FCM.FG.guarded(c / f"{name}_{s}" / "recon.mic.gz"), N_VOXELS)
    ma = np.load(FCM.FG.guarded(c / f"passes_{s}" / f"after_{name}.npy")).reshape(N_VOXELS, 3, 3)
    before = final_rows(rows, mb)
    after = final_rows(rows, ma)
    lk = last_kept(rows)
    out = dict(
        sample=s,
        name=name,
        tol=TOL,
        n_rows=len(rows),
        n_voxels_with_kept_fit=len(lk),
        final_before=len(before),
        final_after=len(after),
        before=sorted(before),
        after=sorted(after),
    )
    (c / "record" / f"final_{name}_{s}.json").write_text(json.dumps(out) + "\n")
    print({k: v for k, v in out.items() if k not in ("before", "after")})


if __name__ == "__main__":
    main()
