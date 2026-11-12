"""AR5 step 3, Part 1: the mixed seed refinement from the replay cache (no new evaluation).

  uv run python scripts/autoresearch/ar5/mix_part1.py --pick   # big, grad: writes pick.json
  uv run python scripts/autoresearch/ar5/mix_part1.py --all    # all three: writes part1.json

Reads the stored replay of AR5 step 2 (rows per seed step with the CMA arms ``C60``/``C200`` and the
MC arm ``M`` for the same candidates; linked into this step's cache) and builds ``X60`` and ``X200``
with ``mix_rules.mixed_step``. Per arm and sample: seed-test passes, ends within 1 degree,
seed-step evaluations, the paired counts against ``M`` (passes and right ends: with ``M`` only,
with the arm only), the number of candidates refined and those where MC still ran, and the steps
that needed an MC end the replay did not compute (``missing``).

``--pick`` reads ``big`` and ``grad`` only and applies step 2's rule (``fit_rules.pick_arm``, with
``M`` as the reference): the arm with the fewer seed-step evaluations among those whose passes (a)
and right ends (b) are at most 1% of the steps below ``M``'s. The pick file is committed before any
whole-layer run; ``--all`` (which reads ``mid``) refuses to run before that commit.
"""

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[2]
MIX = ICE / "benchmarks" / "autoresearch" / "ar5" / "mix"
os.environ.setdefault("RESEED_CACHE", str(MIX / "cache"))
os.environ.setdefault("RESEED_SAMPLE", "big")
sys.path.insert(0, str(HERE))

import fit_common as FC  # noqa: E402
import fit_part1_lib as L  # noqa: E402
import fit_rules as FR  # noqa: E402
import mix_rules as MR  # noqa: E402

PICK_FILE = MIX / "pick.json"
PICK_SAMPLES = ("big", "grad")
ALL_SAMPLES = ("big", "grad", "mid")
MIXED = {"X60": "C60", "X200": "C200"}


def pick_commit() -> str:
    import subprocess

    out = subprocess.run(
        ["git", "log", "-1", "--format=%H %cI", "--", str(PICK_FILE)],
        cwd=str(ICE),
        capture_output=True,
        text=True,
    )
    return out.stdout.strip()


def steps_of(rows: List[Dict[str, Any]], arm: str) -> List[Dict[str, Any]]:
    """One ``mixed_step`` result per replay row (arm ``M``/``C*``: the stored values)."""
    out = []
    for r in rows:
        a = r["arms"]
        if arm in MIXED:
            s = MR.mixed_step(r["costed"], a.get(MIXED[arm]), a.get("M"))
        else:
            d = a.get(arm)
            s = dict(
                evals=int(d["evals"]) if d else 0,
                mc_ran=0,
                n_refined=len(d["fin_hit"]) if d else 0,
                missing=False,
                passed=bool(d) and FR.seed_test(d["hit"]),
                right=bool(d) and d["err"] <= MR.RIGHT_DEG,
            )
            if arm == "M" and d:
                s["mc_ran"] = len(d["fin_hit"])
        out.append(s)
    return out


def summarize(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    arms = ("M", "X60", "X200")
    st = {a: steps_of(rows, a) for a in arms}
    out: Dict[str, Any] = dict(n_steps=len(rows), arms={})
    for a in arms:
        s = st[a]
        n_ref = sum(x["n_refined"] for x in s)
        mc = sum(x["mc_ran"] for x in s)
        out["arms"][a] = dict(
            n_pass=int(sum(x["passed"] for x in s)),
            n_right=int(sum(x["right"] for x in s)),
            evals=int(sum(x["evals"] for x in s)),
            candidates_refined=int(n_ref),
            candidates_mc_ran=int(mc),
            share_mc_ran=float(mc / n_ref) if n_ref else 0.0,
            steps_mc_ran=int(sum(x["mc_ran"] > 0 for x in s)),
            steps_missing=int(sum(x["missing"] for x in s)),
        )
    for a in ("X60", "X200"):
        for key, name in (("passed", "paired_passed"), ("right", "paired_right")):
            d = FR.paired_split([x[key] for x in st["M"]], [x[key] for x in st[a]])
            out["arms"][a][name] = dict(
                both=d["both"], m_only=d["only_a"], arm_only=d["only_b"], neither=d["neither"]
            )
        out["arms"][a]["evals_fall_vs_M"] = FR.fall(
            out["arms"][a]["evals"], out["arms"]["M"]["evals"]
        )
    return out


PLAN_WHATIF_X200 = dict(pick=2_501_991, mid=5_643_374)  # PLAN.md, AR5 step 3, what-if table


def whatif_variant(rows: List[Dict[str, Any]], arm: str = "X200") -> Dict[str, Any]:
    """The plan's what-if variant (stop when either end of a candidate reaches 0.95,
    ``mix_rules.mixed_step(stop_on_either=True)``) against the rule as defined (stop test on the
    kept lower-cost end), with every step where they differ (review fix 6)."""
    diffs: List[Dict[str, Any]] = []
    tot_rule = tot_var = pass_var = 0
    for i, r in enumerate(rows):
        c, m = r["arms"].get(MIXED[arm]), r["arms"].get("M")
        a = MR.mixed_step(r["costed"], c, m)
        b = MR.mixed_step(r["costed"], c, m, stop_on_either=True)
        tot_rule += a["evals"]
        tot_var += b["evals"]
        pass_var += int(b["passed"])
        if a["evals"] == b["evals"] and a["taken"] == b["taken"]:
            continue
        assert c is not None and m is not None
        p = int(b["taken"])
        nxt = p + 1 if p + 1 < len(c["fin_hit"]) else None
        diffs.append(
            dict(
                replay_index=i,
                step=int(r["step"]),
                candidate=p,
                cma_hit=float(c["fin_hit"][p]),
                cma_cost=float(c["fin_cost"][p]),
                cma_evals=int(c["fin_evals"][p]),
                mc_hit=float(m["fin_hit"][p]),
                mc_cost=float(m["fin_cost"][p]),
                mc_evals=int(m["fin_evals"][p]),
                next_cma_evals=(int(c["fin_evals"][nxt]) if nxt is not None else None),
                next_cma_hit=(float(c["fin_hit"][nxt]) if nxt is not None else None),
                evals_rule=int(a["evals"]),
                evals_variant=int(b["evals"]),
                passed_rule=bool(a["passed"]),
                passed_variant=bool(b["passed"]),
            )
        )
    return dict(
        evals_rule=int(tot_rule),
        evals_variant=int(tot_var),
        n_pass_variant=int(pass_var),
        steps_differ=diffs,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--pick", action="store_true")
    g.add_argument("--all", action="store_true")
    a = ap.parse_args()
    samples = PICK_SAMPLES if a.pick else ALL_SAMPLES
    if a.all and not pick_commit():
        raise SystemExit("order rule: mid is not read before pick.json is committed")
    per = {s: L.load(s) for s in samples}
    for s in samples:
        assert per[s], f"no replay rows for {s}"
    by_sample = {s: summarize(per[s]) for s in samples}
    pick_rows = [r for s in PICK_SAMPLES for r in per[s]]
    pooled = summarize(pick_rows)
    n = pooled["n_steps"]
    arms = {k: {f: v[f] for f in ("n_pass", "n_right", "evals")} for k, v in pooled["arms"].items()}
    pick_arm = FR.pick_arm(arms, n)
    rule_a = {k: FR.mid_pass(arms[k], arms["M"], n)["a"] for k in ("X60", "X200")}
    rule_b = {k: FR.mid_pass(arms[k], arms["M"], n)["b"] for k in ("X60", "X200")}
    out = dict(
        rule="the arm with the fewer seed-step evaluations among those whose seed-test passes (a) "
        "and ends within 1 degree (b) are each at most 1% of all steps below M's (net counts, "
        "pooled over the pick samples big and grad); M if neither qualifies",
        pick_samples=list(PICK_SAMPLES),
        n_steps=n,
        one_percent_of_steps=n / 100.0,
        pooled=pooled,
        rule_a=rule_a,
        rule_b=rule_b,
        pick=pick_arm,
        stop_needed=bool(pick_arm == "M"),
    )
    if a.all:
        out["by_sample"] = by_sample
        wi = {s: whatif_variant(per[s]) for s in samples}
        for s in samples:
            assert wi[s]["evals_rule"] == by_sample[s]["arms"]["X200"]["evals"]
        pick_var = sum(wi[s]["evals_variant"] for s in PICK_SAMPLES)
        assert pick_var == PLAN_WHATIF_X200["pick"], pick_var
        assert wi["mid"]["evals_variant"] == PLAN_WHATIF_X200["mid"], wi["mid"]["evals_variant"]
        out["plan_whatif"] = dict(
            note="the plan's what-if count for X200 comes from a variant that stops when either "
            "end of a candidate reaches 0.95; the rule as defined (and run) tests the stop on the "
            "kept lower-cost end. Both counts and every step where they differ.",
            plan_counts=PLAN_WHATIF_X200,
            pick_evals_variant=int(pick_var),
            by_sample=wi,
        )
        dest = MIX / "part1.json"
        old = json.loads(PICK_FILE.read_text())
        assert old["pick"] == pick_arm and old["pooled"] == pooled, "pick file does not match"
    else:
        dest = PICK_FILE
    dest.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: out[k] for k in ("n_steps", "rule_a", "rule_b", "pick")}, indent=1))
    for k, v in pooled["arms"].items():
        print(k, {f: v[f] for f in ("n_pass", "n_right", "evals", "share_mc_ran", "steps_missing")})


if __name__ == "__main__":
    main()
