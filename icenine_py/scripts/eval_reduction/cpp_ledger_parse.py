"""Parse the `LG` lines of the scratch C++ ledger build (benchmarks/eval_reduction/A/step2/
cpp_ledger_patch.py.txt) into cpp_ledger.json.

  uv run python scripts/eval_reduction/cpp_ledger_parse.py --logs DIR \
      --out benchmarks/eval_reduction/A/step2

DIR holds cpp_<voxel>.log, each from `LEDGER=1 <scratch>/IceNine r <config>` run in
Examples/Example2.ManyGrains (one-voxel mic, clean full-sample frames, Q-max 8 config).
"""

import argparse
import json
import re
from pathlib import Path
from typing import Any, Dict

SIX = [15901, 17034, 18558, 19889, 22867, 5242]


def parse(path: Path) -> Dict[str, Any]:
    levels = []
    find: Dict[str, Any] = {}
    for line in path.read_text().splitlines():
        if not line.startswith("LG "):
            continue
        t = line.split()
        if t[1] == "L":
            kv = {t[i]: int(t[i + 1]) for i in range(3, len(t), 2)}
            levels.append(
                dict(
                    level=int(t[2]),
                    n_fz=kv["n_fz"],
                    n_local=kv["n_local"],
                    disc=kv["disc"],
                    rerecost=kv["rerecost"],
                    quick=kv["quick"],
                    cand=kv["cand"],
                    kept=kv["kept"],
                )
            )
        elif t[1] == "FIND":
            find = {t[i].lower(): int(t[i + 1]) for i in range(1, len(t), 2)}
    for lv in levels:
        lv["screen"] = lv["n_fz"] * lv["n_local"]
        lv["recost"] = lv["disc"] - lv["screen"]  # first re-cost: every locally passing pair
    tot = dict(
        screen=sum(l["screen"] for l in levels),
        recost=sum(l["recost"] for l in levels),
        rerecost=sum(l["rerecost"] for l in levels),
        quick=sum(l["quick"] for l in levels),
        find=find["find"],
        variance=find["var"],
        final=find["final"],
        candidates=sum(l["cand"] for l in levels),
    )
    tot["sum"] = (
        tot["screen"]
        + tot["recost"]
        + tot["rerecost"]
        + tot["quick"]
        + tot["find"]
        + tot["variance"]
        + tot["final"]
    )
    return dict(levels=levels, find=find, totals=tot, total_counter=find["total"])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--logs", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    res = {str(v): parse(Path(a.logs) / f"cpp_{v}.log") for v in SIX}
    for v, r in res.items():
        assert r["totals"]["sum"] == r["total_counter"], (v, r["totals"]["sum"], r["total_counter"])
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "cpp_ledger.json").write_text(json.dumps(res, indent=1) + "\n")
    print({v: r["totals"]["sum"] for v, r in res.items()})


if __name__ == "__main__":
    main()
