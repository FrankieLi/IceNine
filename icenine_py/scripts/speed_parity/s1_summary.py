"""S1 Step 5 summary: timings (run_s1_timing.sh), the parity dumps (s1_parity.py dump), the E0 A/B
(cache/e0_s1base against cache/e0_s1) and the batched-cost check, to benchmarks/s1_loops/
summary.json and tables.md (doc tables, ``doc_tables.write_tables``).

  uv run python scripts/speed_parity/s1_summary.py --out benchmarks/s1_loops

Inputs (all under scripts/speed_parity/cache/s1 unless noted): timing/{base,s1}_{prof_eval,
stage_split,search_split}_r<round>.txt and timing/preflight_*.json, base.npz and new.npz,
batched.txt (the output of ``s1_parity.py batched``), and the two E0 caches of
scripts/findoptimal_robustness/cache. No p-values: counts and ratios of medians only; the timings
are labelled by the preflight state (contended unless preflight_exit.txt holds 0). The E0 runtimes
(10 workers, two runs at different times) are kept in summary.json with a note and are not a
timing comparison.
"""

import argparse
import json
import re
import shutil
import statistics
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(ICE / "scripts" / "common"))
sys.path.insert(0, str(ICE / "scripts" / "findoptimal_robustness"))
import doc_tables as DT  # noqa: E402
from csl import reduced_misorientation_deg  # noqa: E402

CACHE = HERE / "cache" / "s1"
TIMING = CACHE / "timing"
E0 = ICE / "scripts" / "findoptimal_robustness" / "cache"
ARMS = ("screen", "recost", "final")
ARM_LABEL = {
    "screen": "screen (pixel radius 3, Q-max 5, C++ rule)",
    "recost": "re-cost arm (pixel radius 3, Q-max 8, C++ rule)",
    "final": "final (triangle, Q-max 8, C++ rule)",
}
FIELDS = (
    "pixel_overlap",
    "pixel_on_detector",
    "peak_overlap",
    "peak_on_detector",
    "quality",
    "n_quality_points",
    "detectors_overlap",
    "ray_cases",
)


def med(xs: List[float]) -> float:
    return float(statistics.median(xs))


def files(arm: str, tool: str) -> List[Path]:
    return sorted(TIMING.glob(f"{arm}_{tool}_r*.txt"), key=lambda p: int(p.stem.rsplit("_r", 1)[1]))


def prof_eval_ops(arm: str) -> Dict[str, int]:
    """aten ops per call of each arm (the same in every round)."""
    keys = {"screen (r=3": "screen", "recost/quickMC": "recost", "final (triangle": "final"}
    out: Dict[str, int] = {}
    for line in files(arm, "prof_eval")[0].read_text().splitlines():
        for k, a in keys.items():
            if line.startswith(k):
                out[a] = int(re.search(r"aten ops/eval=(\d+)", line).group(1))  # type: ignore
    return out


def prof_eval(arm: str) -> Dict[str, List[float]]:
    out: Dict[str, List[float]] = {a: [] for a in ARMS}
    keys = {"screen (r=3": "screen", "recost/quickMC": "recost", "final (triangle": "final"}
    for f in files(arm, "prof_eval"):
        for line in f.read_text().splitlines():
            for k, a in keys.items():
                if line.startswith(k):
                    m = re.search(r"([\d.]+) us/eval", line)
                    out[a].append(float(m.group(1)))  # type: ignore[union-attr]
    return out


def num(pattern: str, text: str) -> Optional[float]:
    m = re.search(pattern, text)
    return float(m.group(1)) if m else None


def stage_split(arm: str) -> Dict[str, List[float]]:
    pats = {
        "evaluate_total": r"evaluate\(\) total\s+([\d.]+) us",
        "bragg_torch": r"Bragg solve \(torch\)\s+([\d.]+) us",
        "rest_of_body": r"rest of evaluate body \(eta \.\.\.\)\s+([\d.]+) us",
        "overlap_batched": r"calculate_overlap_batched\s+([\d.]+) us",
        "stage_a_loop": r"stage A Python loop\s+([\d.]+) us",
        "image_list_loop": r"image-list Python loop\s+([\d.]+) us",
        "stage_d_c_call": r"stage D C call\s+([\d.]+) us\s*$",
        "screen_evaluate": r"screen arm: evaluate ([\d.]+) us",
        "bragg_original": r"Bragg solve, original function\s+([\d.]+) us",
        "bragg_hoisted": r"Bragg solve, hoisted constants\s+([\d.]+) us",
        "stage_a_lookup": r"stage A, frame_of_omegas\s+([\d.]+) us",
        "frame_table_check": r"frame table: identity \+ filled\s+([\d.]+) us",
        "c_call_table": r"table of \d+ entries ([\d.]+) us",
        "c_call_compact": r"compact list of \d+ entries ([\d.]+) us",
    }
    out: Dict[str, List[float]] = {k: [] for k in pats}
    for f in files(arm, "stage_split"):
        text = f.read_text()
        for k, p in pats.items():
            v = None
            for line in text.splitlines():
                m = re.search(p, line)
                if m:
                    v = float(m.group(1))
                    break
            if v is not None:
                out[k].append(v)
    return out


def search_split(arm: str) -> Dict[str, Any]:
    wall, evals, cost, sha, share, per = [], [], [], [], [], []
    decimals = []
    for f in files(arm, "search_split"):
        t = f.read_text()
        m = re.search(r"wall ([\d.]+) s, evaluations (\d+)", t)
        wall.append(float(m.group(1)))  # type: ignore[union-attr]
        decimals.append(len(m.group(1).partition(".")[2]))  # type: ignore[union-attr]
        evals.append(int(m.group(2)))  # type: ignore[union-attr]
        r = re.search(r"result: cost (\S+), orientation sha1 (\w+)", t)
        cost.append(r.group(1))  # type: ignore[union-attr]
        sha.append(r.group(2))  # type: ignore[union-attr]
        s = re.search(r"inside evaluate: [\d.]+ s \(([\d.]+)%\), (\d+) us/eval", t)
        share.append(float(s.group(1)))  # type: ignore[union-attr]
        per.append(float(s.group(2)))  # type: ignore[union-attr]
    return dict(
        wall_s=wall,
        wall_decimals=min(decimals),
        evaluations=evals,
        cost=cost,
        sha1=sha,
        evaluate_share_pct=share,
        us_per_eval=per,
    )


def parity() -> Dict[str, Any]:
    a = np.load(CACHE / "base.npz", allow_pickle=False)
    b = np.load(CACHE / "new.npz", allow_pickle=False)
    meta = json.loads(str(a["meta"]))
    meta_new = json.loads(str(b["meta"]))
    keys = set(a.files) - {"meta"}
    if keys != set(b.files) - {"meta"}:
        raise ValueError("the two dumps do not have the same arrays; run s1_parity.py compare")
    # the same orientations (R) and the same orientation subset of each mode (index)
    same_sets = all(
        np.array_equal(a[k], b[k]) for k in keys if k.endswith("/R") or k.endswith("/index")
    )
    n_rows = n_fields = n_diff = 0
    per_mode: Dict[str, Dict[str, int]] = {}
    for k in sorted(keys):
        parts = k.split("/")
        if len(parts) != 3 or parts[2] == "index":
            continue
        ds, mode, f = parts
        x, y = a[k], b[k]
        # a shape or dtype mismatch counts every value as different
        n_bad = int(x.size) if x.shape != y.shape or x.dtype != y.dtype else int((x != y).sum())
        n = len(x)
        m = per_mode.setdefault(mode, dict(orientations=0, fields_compared=0, differing_values=0))
        if f == "quality":
            m["orientations"] += n
        m["fields_compared"] += 1
        m["differing_values"] += n_bad
        n_fields += 1
        n_diff += n_bad
        n_rows += n if f == "quality" else 0
    return dict(
        datasets=meta["datasets"],
        modes=meta["modes"],
        n_orientation_evaluations=n_rows,
        n_field_arrays=n_fields,
        n_differing_values=n_diff,
        per_mode=per_mode,
        coverage=meta["coverage"],
        base_head=meta["head"],
        base_dirty=meta.get("dirty"),
        new_head=meta_new["head"],
        new_dirty=meta_new.get("dirty"),
        new_torch=meta_new["torch"],
        same_orientation_sets_and_indices=same_sets,
    )


def e0() -> Dict[str, Any]:
    sets = {t: E0 / t for t in ("e0_s1base", "e0_s1")}
    names = sorted(p.name for p in sets["e0_s1base"].glob("v*.npz"))
    first200 = sorted(p.name for p in (E0 / "e0_s3off").glob("v*.npz"))
    same = 0
    diff_keys: Dict[str, int] = {}
    tot: Dict[str, Dict[str, Any]] = {
        t: dict(n=0, wrong=0, evals_global=0, evals_local=0, n200=0, wrong200=0, runtime_s=0.0)
        for t in sets
    }
    for n in names:
        d = {t: np.load(p / n, allow_pickle=False) for t, p in sets.items()}
        x, y = d["e0_s1base"], d["e0_s1"]
        bad = [
            k
            for k in x.files
            if "runtime" not in k and (k not in y.files or not np.array_equal(x[k], y[k]))
        ] + [k for k in y.files if k not in x.files]
        for k in bad:
            diff_keys[k] = diff_keys.get(k, 0) + 1
        same += not bad
        for t, z in d.items():
            wrong = bool(reduced_misorientation_deg(z["s0_R_final"], z["R_true"]) > 1.0)
            tt = tot[t]
            tt["n"] += 1
            tt["wrong"] += wrong
            tt["evals_global"] += int(z["s0_evals_global"])
            tt["evals_local"] += int(z["s0_evals_local"])
            tt["runtime_s"] += float(z["s0_runtime"])
            if n in first200:
                tt["n200"] += 1
                tt["wrong200"] += wrong
    return dict(
        voxels=len(names),
        records_identical=same,
        differing_keys=diff_keys,
        runs=tot,
        runtime_note=(
            "runtime_s is the sum of the per-voxel wall-clock s0_runtime of a 10-worker run; the "
            "two runs were made at different times, not interleaved: contended, not a timing "
            "comparison (the timings are in the interleaved single-worker rounds)"
        ),
    )


def batched_check() -> Dict[str, Any]:
    f = CACHE / "batched.txt"
    if not f.exists():
        return {}
    lines = [ln for ln in f.read_text().splitlines() if ln.split()[:1] and "RESULT" not in ln]
    rows = [ln.split() for ln in lines if ln.endswith("identical") or "DIFFER" in ln]
    return dict(
        rows=len(rows),
        identical=sum(1 for r in rows if r[-1] == "identical"),
        modes=sorted({r[1] for r in rows}),
        variants=sorted({r[2] for r in rows}),
        result=[ln for ln in f.read_text().splitlines() if ln.startswith("RESULT")],
    )


def derived(s: Dict[str, Any]) -> Dict[str, Any]:
    """Numbers quoted in prose that are computed from the tables (so the audit finds them)."""
    ss = s["stage_split_us"]

    def m(arm: str, k: str) -> float:
        return med(ss[arm][k])

    stage_a = m("base", "stage_a_loop") - m("s1", "stage_a_lookup")
    frames = m("base", "image_list_loop") - m("s1", "frame_table_check")
    bragg = m("s1", "bragg_original") - m("s1", "bragg_hoisted")
    total = m("base", "evaluate_total") - m("s1", "evaluate_total")
    se = s["search"]
    tab, cmp_ = ss["s1"]["c_call_table"], ss["s1"]["c_call_compact"]
    dec = min(se[arm]["wall_decimals"] for arm in ("base", "s1"))
    half = 0.5 * 10 ** (-dec)  # the walls are printed to ``dec`` decimals
    w_base, w_s1 = med(se["base"]["wall_s"]), med(se["s1"]["wall_s"])
    pe = s["prof_eval_us"]
    recost_saved = med(pe["base"]["recost"]) - med(pe["s1"]["recost"])
    search_saved = med(se["base"]["us_per_eval"]) - med(se["s1"]["us_per_eval"])
    a = np.load(CACHE / "base.npz", allow_pickle=False)
    meta = json.loads(str(a["meta"]))
    edge_changed = {
        ds: int((a[f"{ds}/local/quality"] != a[f"{ds}/edge/quality"]).sum())
        for ds in meta["datasets"]
    }
    fx = np.load(ICE / "tests" / "data" / "s1_evaluate_ref.npz", allow_pickle=False)
    fixture_edge = int(
        (fx["three0/local/peak_on_detector"] != fx["three0/edge/peak_on_detector"]).sum()
    )
    n_unique = None
    for line in (TIMING / "s1_stage_split_r1.txt").read_text().splitlines():
        mm = re.search(r"unique wedges=(\d+)", line)
        if mm:
            n_unique = int(mm.group(1))
    return dict(
        parts_removed_saved_us=stage_a + frames + bragg,
        evaluate_total_saved_us=total,
        unattributed_us=total - (stage_a + frames + bragg),
        search_wall_ratio=w_base / w_s1,
        search_wall_decimals=dec,
        search_wall_ratio_rounding_low=(w_base - half) / (w_s1 + half),
        search_wall_ratio_rounding_high=(w_base + half) / (w_s1 - half),
        search_saved_us_per_evaluate=search_saved,
        recost_arm_saved_us_per_call=recost_saved,
        search_saving_over_recost_saving=search_saved / recost_saved,
        c_call_table_rounds=len(tab),
        c_call_table_rounds_higher=sum(t > c for t, c in zip(tab, cmp_)),
        c_call_table_rounds_equal=sum(t == c for t, c in zip(tab, cmp_)),
        c_call_table_rounds_lower=sum(t < c for t, c in zip(tab, cmp_)),
        c_call_table_minus_compact_us=med(tab) - med(cmp_),
        edge_rules_change_quality_orientations_of_1000=edge_changed,
        fixture_edge_changes_peak_on_detector_of_64=fixture_edge,
        fixture_size_bytes=(ICE / "tests" / "data" / "s1_evaluate_ref.npz").stat().st_size,
        frames_in_use=n_unique,
        dump_design=dict(
            near_truth=meta["n_near"],
            fz=meta["n_fz"],
            slow_near_fz=meta["n_slow"],
            seed=meta["seed"],
        ),
        orientations_per_voxel=meta["n_near"] + meta["n_fz"],
        slow_orientations_per_voxel=sum(meta["n_slow"]),
        torch_version=meta["torch"],
        plan_estimate_us="55 to 60 (plan)",
        plan_search_ratio=1.15,
        preflight_max_load=1.5,
        prof_eval_orientations=400,
    )


def quiet_repeat() -> Dict[str, Any]:
    """The preflight checks made for the quiet-machine repeat that the owner asked for after the
    review (D-20261014-I16): run the timing again only if one of them passes."""
    files_ = sorted(CACHE.glob("review_preflight_*.json"))
    checks = []
    for f in files_:
        info = json.loads(f.read_text())
        checks.append(
            dict(
                timestamp=info["timestamp"],
                load_1min=info["loadavg"][0],
                busy_processes=info["busy_processes"],
                quiet=info["loadavg"][0] <= 1.5 and not info["busy_processes"],
            )
        )
    passed = bool(files_) and any(c["quiet"] for c in checks)
    return dict(
        quiet_check_passed=passed,
        note=(
            "quiet-machine repeat run (see label)"
            if passed
            else (
                "quiet-machine repeat not done (preflight not quiet), "
                "owner decision D-20261014-I16"
            )
        ),
        max_load=1.5,
        checks=checks,
    )


def build() -> Dict[str, Any]:
    pre = json.loads((TIMING / "preflight_start.json").read_text())
    exit_code = int((TIMING / "preflight_exit.txt").read_text().strip())
    label = "single-worker, quiet" if exit_code == 0 else "single-worker, contended"
    s: Dict[str, Any] = dict(label=label, preflight_start=pre, preflight_exit=exit_code)
    end = TIMING / "preflight_end.json"
    if end.exists():
        s["preflight_end"] = json.loads(end.read_text())
    s["prof_eval_us"] = {arm: prof_eval(arm) for arm in ("base", "s1")}
    s["prof_eval_aten_ops"] = {arm: prof_eval_ops(arm) for arm in ("base", "s1")}
    s["stage_split_us"] = {arm: stage_split(arm) for arm in ("base", "s1")}
    s["search"] = {arm: search_split(arm) for arm in ("base", "s1")}
    s["quiet_repeat_preflights"] = quiet_repeat()
    s["parity"] = parity()
    s["e0"] = e0()
    s["batched_check"] = batched_check()
    s["derived"] = derived(s)
    return s


def tables(s: Dict[str, Any]) -> Dict[str, str]:
    t: Dict[str, str] = {}
    pe = s["prof_eval_us"]
    rows = []
    for a in ARMS:
        b, n = pe["base"][a], pe["s1"][a]
        rows.append(
            {
                "arm": ARM_LABEL[a],
                "base commit (us)": med(b),
                "S1 (us)": med(n),
                "saved (us)": med(b) - med(n),
                "base / S1": med(b) / med(n),
                "base min": min(b),
                "S1 min": min(n),
                "rounds": len(b),
                "aten ops per call, base / S1": (
                    f"{s['prof_eval_aten_ops']['base'][a]} / {s['prof_eval_aten_ops']['s1'][a]}"
                ),
            }
        )
    cols = list(rows[0])
    t["s1_per_call"] = DT.markdown_table(
        rows,
        cols,
        {c: ".1f" for c in cols if c not in ("arm", "rounds", "base / S1")} | {"base / S1": ".2f"},
    )
    ss = s["stage_split_us"]

    def m(arm: str, k: str) -> float:
        return med(ss[arm][k])

    parts = [
        (
            "evaluate() total (Q-max 8, 110 valid peaks)",
            m("base", "evaluate_total"),
            m("s1", "evaluate_total"),
        ),
        (
            "stage A: loop over the peaks / frame_of_omegas",
            m("base", "stage_a_loop"),
            m("s1", "stage_a_lookup"),
        ),
        (
            "frame list: 174 get_image calls / table check",
            m("base", "image_list_loop"),
            m("s1", "frame_table_check"),
        ),
        (
            "Bragg solve on g: original function / hoisted constants",
            m("s1", "bragg_original"),
            m("s1", "bragg_hoisted"),
        ),
        (
            "stage D C call: compact list of 174 / table of 360 entries",
            m("s1", "c_call_compact"),
            m("s1", "c_call_table"),
        ),
    ]
    t["s1_stage_split"] = DT.markdown_table(
        [{"part": p, "before (us)": b, "after (us)": a, "saved (us)": b - a} for p, b, a in parts],
        ["part", "before (us)", "after (us)", "saved (us)"],
        {"before (us)": ".1f", "after (us)": ".1f", "saved (us)": ".1f"},
    )
    se = s["search"]
    rows = []
    for arm, lab in (("base", "base commit"), ("s1", "S1")):
        x = se[arm]
        rows.append(
            {
                "arm": lab,
                "wall (s)": med(x["wall_s"]),
                "evaluate calls": x["evaluations"][0],
                "us per evaluate": med(x["us_per_eval"]),
                "evaluate share (%)": med(x["evaluate_share_pct"]),
                "cost": x["cost"][0],
                "orientation sha1": x["sha1"][0],
                "same in all rounds": len(set(x["evaluations"])) == 1 and len(set(x["sha1"])) == 1,
            }
        )
    t["s1_search"] = DT.markdown_table(
        rows,
        list(rows[0]),
        {
            "wall (s)": f".{min(se[a]['wall_decimals'] for a in ('base', 's1'))}f",
            "us per evaluate": ".0f",
            "evaluate share (%)": ".1f",
        },
    )
    p = s["parity"]
    rows = [
        {
            "mode": k,
            "orientations (4 voxels)": v["orientations"],
            "field arrays": v["fields_compared"],
            "differing values": v["differing_values"],
        }
        for k, v in p["per_mode"].items()
    ]
    t["s1_parity"] = DT.markdown_table(rows, list(rows[0]))
    cov = []
    for ds, c in p["coverage"].items():
        for lab in ("local", "screen"):
            v = c[lab]
            cov.append(
                {
                    "voxel": ds,
                    "cost": lab,
                    "omega rows": v["rows"],
                    "in [low - width, low)": v["below_edge_bin"],
                    "above the range": v["above_range"],
                    "further below": v["far_below"],
                    "in gap bins": v["gap_bins"],
                    "frames used": v["frames_used"],
                }
            )
    t["s1_coverage"] = DT.markdown_table(cov, list(cov[0]))
    e = s["e0"]
    rows = []
    for tag, lab in (("e0_s1base", "base commit"), ("e0_s1", "S1")):
        r = e["runs"][tag]
        rows.append(
            {
                "run": lab,
                "voxels": r["n"],
                "wrong (> 1 deg)": r["wrong"],
                "wrong, the 200 voxels of e0_s3off": r["wrong200"],
                "evaluations, global": r["evals_global"],
                "evaluations, local": r["evals_local"],
            }
        )
    t["s1_e0"] = DT.markdown_table(rows, list(rows[0]))
    bc = s.get("batched_check")
    if bc:
        t["s1_batched"] = DT.markdown_table(
            [
                {
                    "check": (
                        "BatchedCostEvaluator (CPU compact and fixed) against the base evaluate"
                    ),
                    "modes": ", ".join(bc["modes"]),
                    "(voxel, mode, variant) rows": bc["rows"],
                    "identical": bc["identical"],
                }
            ],
            ["check", "modes", "(voxel, mode, variant) rows", "identical"],
        )
    return t


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = ICE / a.out if not Path(a.out).is_absolute() else Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    s = build()
    (out / "summary.json").write_text(json.dumps(s, indent=1, default=str) + "\n")
    DT.write_tables(out / "tables.md", tables(s))
    shutil.copy(TIMING / "preflight_start.json", out / "preflight.json")
    print(f"wrote {out}/summary.json, tables.md ({s['label']})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
