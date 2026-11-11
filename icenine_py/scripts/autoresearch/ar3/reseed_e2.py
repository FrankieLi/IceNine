"""AR3 step 8, part 2, remedy E2: keep a rejected fit that passes a second test.

  uv run python scripts/autoresearch/ar3/reseed_e2.py

Stored cases only (the logged rerun of S5 on ``big``; the first gradient sample has no per-fit log,
so E2 and the classes are not evaluated there). Reads the slim files ``rejected_big``,
``members_big``, ``clusters_big`` and ``big/errors.npz`` (the stored S5 map's error per voxel before
the pass).

Cases. For every voxel with at least one rejected neighbour fit (phase "neighbor"), its first
rejected fit. E2 keeps that fit when the second test passes; its end is then the voxel's
orientation.
Judged with the truth: the kept end is right when its error is at most 1 degree. Paired with the
stored run on the same voxel (its final map before the pass): fixed = stored wrong, E2 right;
broken = stored right, E2 wrong. "Settled wrong" = kept and wrong.

Evaluations (an estimate from stored cases, not a rerun). A reseed cluster would not exist if every
voxel in it had been kept by E2 (judged on each member's last rejected fit before that reseed);
its reseed (the seed step and the flood's fits) would not be spent. The cut is the evaluations
of such clusters over the reseed total of the set. A cluster that is only smaller still costs a
reseed, so the estimate leaves its saving out. The sets are by true grain parity (pick: even; test:
odd); a case belongs to the set of its voxel's grain, a cluster to the set of its reseed voxel's
grain (46 of 314 clusters reach into two grains, 3 into three).

Second tests tried on the pick set (all use the fit's own hit ratio ``h``):
  abs-T     h >= T, T in 0.60, 0.65, 0.70, 0.75, 0.80, 0.85 (0.80 is the seed test)
  nb-f      h >= f * (hit ratio of the neighbour the fit started from), f in 0.7, 0.8, 0.9
Pick rule (fixed before the test set is summed): among the tests with broken <= fixed on the pick
set, the one with the largest cut; ties go to the stricter test. If none qualifies, E2 fails.
"""

import gzip
import json
from pathlib import Path
from typing import Any, Callable, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE.parents[2] / "benchmarks" / "autoresearch" / "ar3" / "reseed"
BIG = HERE.parents[2] / "benchmarks" / "autoresearch" / "ar3" / "big"
SLIM = OUT / "slim"
RIGHT = 1.0
ABS = (0.60, 0.65, 0.70, 0.75, 0.80, 0.85)
REL = (0.7, 0.8, 0.9)


def rj(name: str) -> Any:
    with gzip.open(SLIM / name, "rt") as f:
        return json.load(f)


def tests() -> Dict[str, Callable[[float, Any], bool]]:
    d: Dict[str, Callable[[float, Any], bool]] = {}
    for t in ABS:
        d[f"abs-{t:.2f}"] = lambda h, e, t=t: h >= t
    for f in REL:
        d[f"nb-{f:.1f}"] = lambda h, e, f=f: e is not None and h >= f * e
    return d


def main() -> None:
    rej = rj("rejected_big.json.gz")
    mem = rj("members_big.json.gz")
    clu = rj("clusters_big.json.gz")
    z = np.load(BIG / "errors.npz")
    stored_err = z["S5_before"].astype(float)
    # ---- cases: first rejected neighbour fit of each voxel
    first: Dict[int, Dict[str, Any]] = {}
    for r in sorted(rej, key=lambda r: r["row"]):
        if r["phase"] == 0:
            first.setdefault(r["voxel"], r)
    cases = list(first.values())
    by_cluster: Dict[int, List[Dict[str, Any]]] = {}
    for m in mem:
        by_cluster.setdefault(m["k"], []).append(m)
    T = tests()

    def summarize(name: str, parity: int) -> Dict[str, Any]:
        f = T[name]
        cs = [c for c in cases if c["grain"] % 2 == parity]
        kept = [c for c in cs if f(c["hit"], c["exp_hit"])]
        right = np.array([c["err_end"] <= RIGHT for c in kept], bool)
        srt = np.array([stored_err[c["voxel"]] <= RIGHT for c in kept], bool)
        cl = [c for c in clu if c["grain"] % 2 == parity]
        saved = 0
        n_el = 0
        for c in cl:
            ms = by_cluster[c["k"]]
            if all(m["cls"] != "N" and f(m["hit"], m["exp_hit"]) for m in ms):
                saved += c["seed_evals"] + c["flood_evals"]
                n_el += 1
        tot = sum(c["seed_evals"] + c["flood_evals"] for c in cl)
        return dict(
            cases=len(cs),
            kept=len(kept),
            settled_wrong=int((~right).sum()),
            fixed=int((~srt & right).sum()),
            broken=int((srt & ~right).sum()),
            stored_wrong_among_kept=int((~srt).sum()),
            clusters=len(cl),
            clusters_gone=n_el,
            reseed_evals=int(tot),
            saved_evals=int(saved),
            cut=float(saved / tot) if tot else 0.0,
        )

    pick = {n: summarize(n, 0) for n in T}
    ok = [n for n in T if pick[n]["broken"] <= pick[n]["fixed"]]
    chosen = max(ok, key=lambda n: (pick[n]["cut"], list(T).index(n))) if ok else None
    # tie rule: the stricter test = the one with fewer kept
    if ok:
        best = max(pick[n]["cut"] for n in ok)
        tied = [n for n in ok if pick[n]["cut"] == best]
        chosen = min(tied, key=lambda n: pick[n]["kept"])
    pk = json.loads((OUT / "remedies_pick.json").read_text())
    pk["E2"] = dict(
        second_tests="abs-T: h >= T for T in %s; nb-f: h >= f * neighbour's hit for f in %s"
        % (list(ABS), list(REL)),
        pick_set=pick,
        chosen=chosen,
        chosen_by="largest estimated cut among tests with broken <= fixed on the pick set; ties "
        "to the test that keeps fewer; none: E2 fails",
    )
    (OUT / "remedies_pick.json").write_text(json.dumps(pk, indent=1) + "\n")
    out: Dict[str, Any] = dict(chosen=chosen, pick=pick)
    if chosen is not None:
        out["test_set"] = {n: summarize(n, 1) for n in T}
        out["all"] = {
            n: {k: (pick[n][k] + out["test_set"][n][k]) if k != "cut" else None for k in pick[n]}
            for n in T
        }
    (OUT / "remedies_e2.json").write_text(json.dumps(out, indent=1) + "\n")
    print("cases", len(cases), "chosen", chosen)
    for n in T:
        p = pick[n]
        t = out.get("test_set", {}).get(n)
        print(
            n,
            "pick",
            {k: p[k] for k in ("kept", "settled_wrong", "fixed", "broken", "clusters_gone", "cut")},
        )
        if t:
            print(
                "   test",
                {
                    k: t[k]
                    for k in ("kept", "settled_wrong", "fixed", "broken", "clusters_gone", "cut")
                },
            )


if __name__ == "__main__":
    main()
