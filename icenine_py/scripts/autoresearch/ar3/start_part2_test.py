"""AR3 step 9, part 2: the test set. Needs the pick (``start/f_pick.json``) committed first.

  uv run python scripts/autoresearch/ar3/start_part2_test.py

Reads the fits of the test set (true grains with an odd index) of the two samples, applies the
setting that ``f_pick.json`` chose for each remedy and writes ``start/part2_verdict.json``: the
pass rule of the plan. A remedy passes if on the test set (pooled over the two samples) it cuts
the evaluations of its target by at least a third and breaks no more cases than it fixes. The
tables of the other settings are written too, for the reader, and are not used for the verdict.
The same table for each sample alone is given; the verdict is on the pooled set.
"""

import json
import subprocess
from typing import Any, Dict

import start_part2_pick as PK

ROOT = PK.HERE.parents[3]
REL = "icenine_py/benchmarks/autoresearch/ar3/start/f_pick.json"
CUT = 1.0 / 3.0


def pick_is_committed() -> bool:
    tracked = subprocess.run(
        ["git", "-C", str(ROOT), "ls-files", "--error-unmatch", REL], capture_output=True
    )
    dirty = subprocess.run(
        ["git", "-C", str(ROOT), "status", "--porcelain", "--", REL], capture_output=True, text=True
    ).stdout.strip()
    return tracked.returncode == 0 and not dirty


def main() -> None:
    assert pick_is_committed(), "commit start/f_pick.json before the test set is read"
    pick = json.loads((PK.START / "f_pick.json").read_text())
    by_sample = {}
    for s in PK.SAMPLES:
        f = PK.load_feat(s)
        by_sample[s] = PK.f1_table(f, (f["grain"] % 2) == 1)
    f1 = PK.pool_f1(list(by_sample.values()))
    out: Dict[str, Any] = dict(sets="test = odd true-grain index", F1=dict(by_sample=by_sample))
    out["F1"]["pooled"] = f1
    t = pick["F1"]["chosen_threshold"]
    if t is None:
        out["F1"]["verdict"] = dict(passes=False, reason="no threshold qualified on the pick set")
    else:
        r = f1[str(t)]
        out["F1"]["verdict"] = dict(
            threshold=t,
            cut=r["cut"],
            fixed=r["fixed"],
            broken=r["broken"],
            passes=bool(r["cut"] >= CUT and r["broken"] <= r["fixed"]),
        )
    f2_files = []
    for s in PK.SAMPLES:
        p = PK.START / f"part2_test_{s}_f2.json"
        if p.exists():
            f2_files.append(json.loads(p.read_text()))
    if len(f2_files) == len(PK.SAMPLES):
        f2 = PK.pool_f2(f2_files)
        out["F2"] = dict(pooled=f2, by_sample={x["sample"]: x for x in f2_files})
        sg = pick["F2"].get("chosen_sigma_deg")
        if sg is None:
            out["F2"]["verdict"] = dict(passes=False, reason="no setting qualified on the pick set")
        else:
            r = f2[str(sg)]
            out["F2"]["verdict"] = dict(
                sigma_deg=sg,
                cut=r["cut"],
                fixed=r["fixed"],
                broken=r["broken"],
                passes=bool(r["cut"] >= CUT and r["broken"] <= r["fixed"]),
            )
    else:
        out["F2"] = dict(note="F2 test files missing")
    (PK.START / "part2_verdict.json").write_text(json.dumps(out, indent=1) + "\n")
    print(
        json.dumps({k: v.get("verdict") for k, v in out.items() if isinstance(v, dict)}, indent=1)
    )


if __name__ == "__main__":
    main()
