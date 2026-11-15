"""AR5 step 17, phase A step 5: the selection, exactly as frozen (``fitskip_rules``).

  uv run python scripts/autoresearch/ar5/fitskip_select.py [--cache DIR]

Replays the 27 rules on the recording of the five disclosed samples (``<cache>/record/J5r_<s>.
jsonl.gz`` with the final-answer flags of ``fitskip_final``), applies the frozen admissibility
tests and the frozen ranking, and writes (tracked, small) ``ar5/fitskip/select.json`` (every rule,
every sample), ``pick.json`` (the pick and the numbers of the five samples) and ``tables.md``
(``doc_tables``). No rule is tuned and there is no second pick. Needs ``ar5/fitskip/frozen.json``
committed and agreeing (``fitskip_freeze``). Symmetry for rule ``D``: the proper rotations of the
fcc cubic group (copper, as the samples).
"""

import gzip
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import fitskip_freeze as FF  # noqa: E402

import fitskip_rules as FR  # noqa: E402

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
FS = ICE / "benchmarks" / "autoresearch" / "ar5" / "fitskip"
sys.path.insert(0, str(ICE / "scripts" / "common"))


def load_rows(cache: Path, s: str) -> List[Dict[str, Any]]:
    with gzip.open(cache / "record" / f"J5r_{s}.jsonl.gz", "rt") as f:
        rows = [json.loads(x) for x in f if x.strip()]
    fin = json.loads((cache / "record" / f"final_J5r_{s}.json").read_text())
    b = set(fin["before"])
    for r in rows:
        r["final_before"] = int(r["i"]) in b
    return rows


def main() -> None:
    cache = FS / "cache"
    if "--cache" in sys.argv:
        cache = Path(sys.argv[sys.argv.index("--cache") + 1])
    FF.check_or_exit()
    from icenine import ray_indexing as RI
    from icenine.symmetry import create_fcc_symmetry

    sym = RI.proper_rotations(create_fcc_symmetry(3.6, "Cu"))
    results: Dict[str, Dict[str, Dict[str, Any]]] = {r.name: {} for r in FR.family()}
    for s in FR.SAMPLES:
        rows = load_rows(cache, s)
        assert len(rows) == FR.EXPECT[s]["revisits"] + FR.EXPECT[s]["nbr"], s
        for r in FR.family():
            results[r.name][s] = FR.replay(r, rows, sym)
        print(s, len(rows), flush=True)
    sel = FR.select(results)
    out = dict(
        rules=results,
        n_admissible=sel["n_admissible"],
        ranking=sel["admissible"],
        pick=sel["pick"],
    )
    (FS / "select.json").write_text(json.dumps(out, indent=1) + "\n")
    pick = sel["pick"]
    pj: Dict[str, Any] = dict(pick=pick, n_admissible=sel["n_admissible"])
    if pick:
        per = results[pick]
        pj["per_sample"] = per
        pj["net_share"] = {s: per[s]["net"] / FR.EXPECT[s]["total"] for s in FR.SAMPLES}
        pj["g2"] = FR.gate_g2(per)
    (FS / "pick.json").write_text(json.dumps(pj, indent=1) + "\n")
    import doc_tables

    ok = {x["rule"]: i for i, x in enumerate(sel["admissible"], start=1)}
    lines = [
        "Reading line: each row is one rule replayed on the recording of the five samples; a wrong "
        "skip is a skipped fit that J5 accepted; net share is the saved evaluations of skipped "
        "fits that J5 rejected, minus every check, as a share of J5's total with the pass; "
        "admissible means at most 0.1% wrong skips and at most 2 changed final answers on every "
        "sample; rank is by the sum of net shares among the admissible rules.",
        "",
        "| rule | admissible | rank | sum of net shares (5 samples) "
        "| wrong skips (held3 / held / mid / grad / big) "
        "| net share % (held3 / held / mid / grad / big) |",
        "|---|---|---|---|---|---|",
    ]
    for r in FR.family():
        per = results[r.name]
        adm = r.name in ok
        wrong = " / ".join(str(per[s]["wrong"]) for s in FR.SAMPLES)
        net = " / ".join(f"{100 * per[s]['net'] / FR.EXPECT[s]['total']:.2f}" for s in FR.SAMPLES)
        tot = float(FR.net_share_sum(per))
        lines.append(
            f"| {r.name} | {'yes' if adm else 'no'} | {ok.get(r.name, '')} | {tot:.4f} | {wrong} |"
            f" {net} |"
        )
    doc_tables.write_tables(FS / "tables.md", {"rules": "\n".join(lines)})
    print("pick:", pick, "admissible:", sel["n_admissible"])


if __name__ == "__main__":
    main()
