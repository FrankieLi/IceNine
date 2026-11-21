"""AR5 step 18, phase A steps 4 and 6: the replay of the learned skip on the recording ``J5r``
(leave one sample out), and the export of the gate model and the final model.

  uv run python scripts/autoresearch/ar5/learnskip_replay.py loso     # g0, g1: the table
  uv run python scripts/autoresearch/ar5/learnskip_replay.py final    # final model (after the gate)

Applies the frozen modules (``learnskip_rules``, ``_features``, ``_model``); nothing is chosen here.
Reads the five disclosed samples' rows of ``J5r`` (``<cache>/record/``, linked read only from step
17's worktree by ``learnskip_link_cache.sh``) through ``learnskip_guard.guarded``; ``held2`` and
``held4`` are never opened. ``loso``:

* g0: per sample the row counts, the kinds and the rejected counts against the stored table of arm
  ``J5``, and step 17's fit-count check (``p1_J5r_<s>.json``: maps before and after the pass,
  counts, totals, all equal to the stored run);
* for each held-out sample: the threshold is chosen from the four others (out-of-fold scores:
  each of the four scored by a model trained on the other three), a model is trained on the four
  and scores the fifth; the counts, the net share, wrong skips, changed final answers are written
  with ``S(0.3, both)`` and the ceiling, and compared with the plan's table (prediction 1);
* the gate model for ``mid`` (the one trained on the other four) is exported to
  ``<cache>/model_gate_mid.json``.

``final``: out-of-fold scores of all five (each by a model on the other four), the threshold by the
same procedure, the model on all five, exported to ``<out>/model_final.json`` (tracked, under 3 MB)
with the training hash.
"""

import gzip
import hashlib
import json
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import learnskip_features as LFE  # noqa: E402
import learnskip_guard as LG  # noqa: E402
import learnskip_model as LM  # noqa: E402
import learnskip_rules as LR  # noqa: E402

ICE = HERE.parents[2]
OUT = ICE / "benchmarks" / "autoresearch" / "ar5" / "learnskip"
CACHE = LG.cache_dir()
SAMPLES = LR.SAMPLES


def read_rows(s: str) -> List[Dict[str, Any]]:
    assert s in LG.DISCLOSED
    rows = []
    with gzip.open(LG.guarded(CACHE / "record" / f"J5r_{s}.jsonl.gz"), "rt") as f:
        for line in f:
            rows.append(json.loads(line))
    return rows


def g0(s: str, rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    """The recording's counts against the stored arm J5 (the plan's table) and step 17's check."""
    e = LR.EXPECT[s]
    rev = [r for r in rows if r["kind"] == "revisit"]
    nbr = [r for r in rows if r["kind"] == "neighbour"]
    got = dict(
        revisits=len(rev),
        rej_revisits=sum(not r["accepted"] for r in rev),
        nbr=len(nbr),
        nbr_rej=sum(not r["accepted"] for r in nbr),
    )
    p1 = json.loads(LG.guarded(CACHE / f"p1_J5r_{s}.json").read_text())
    return dict(
        sample=s,
        n_fits=len(rows),
        n_fits_ok=len(rows) == LR.FITS[s],
        counts=got,
        counts_ok=all(got[k] == e[k] for k in got),
        p1_ok=bool(p1["ok"]),
        p1_total=p1["counts"]["got"]["total"],
        total_ok=p1["counts"]["got"]["total"] == e["total"],
    )


def load(s: str) -> LR.SampleData:
    rows = read_rows(s)
    arr = LFE.rows_to_arrays(rows)
    fin = json.loads(LG.guarded(CACHE / "record" / f"final_J5r_{s}.json").read_text())
    before = set(int(i) for i in fin["before"])
    fb = np.array([int(r["i"]) in before for r in rows])
    d = LR.SampleData(s, arr["X"], arr["accepted"], arr["evals"], fb)
    d.g0 = g0(s, rows)
    d.rows_hash = hashlib.sha256(
        np.ascontiguousarray(arr["X"]).tobytes()
        + arr["accepted"].tobytes()
        + arr["evals"].tobytes()
    ).hexdigest()
    return d


def fit_score(train: List[str], data: Dict[str, LR.SampleData], test: str) -> np.ndarray:
    X = np.vstack([data[t].X for t in train])
    y = np.concatenate([data[t].accepted for t in train])
    m = LM.train(X, y)
    return m.predict_proba(data[test].X)[:, 1]


def oof_scores(train: List[str], data: Dict[str, LR.SampleData]) -> Dict[str, np.ndarray]:
    return {u: fit_score([t for t in train if t != u], data, u) for u in train}


def pct(x: Fraction) -> float:
    return round(float(x) * 100, 2)


def loso() -> Dict[str, Any]:
    t0 = time.time()
    data = {s: load(s) for s in SAMPLES}
    print("loaded", {s: len(d.accepted) for s, d in data.items()}, flush=True)
    out: Dict[str, Any] = dict(g0={s: data[s].g0 for s in SAMPLES}, samples={}, wall_s=None)
    for h in SAMPLES:
        tr = [s for s in SAMPLES if s != h]
        oof = oof_scores(tr, data)
        pick = LR.choose_threshold(oof, {s: data[s] for s in tr})
        X = np.vstack([data[t].X for t in tr])
        y = np.concatenate([data[t].accepted for t in tr])
        m = LM.train(X, y)
        scores = m.predict_proba(data[h].X)[:, 1]
        a = LR.account(LR.skip_vector(scores, pick["threshold"]), data[h])
        sa = LR.account(LR.s_rule_skip(data[h]), data[h])
        ceil = LR.ceiling(data[h])
        cmp = LR.matches_plan_table(h, a, sa, ceil)
        out["samples"][h] = dict(
            threshold=pick["threshold"],
            train_share_sum=pick["share_sum"],
            rounds=int(m.n_iter_),
            n=a["n"],
            skipped=a["skipped"],
            wrong=a["wrong"],
            allowed=int(LR.WRONG_MAX_SHARE * a["n"]),
            changed=a["changed"],
            saved=a["saved"],
            net=a["net"],
            share_pct=pct(a["share"]),
            s_net=sa["net"],
            s_share_pct=pct(sa["share"]),
            s_wrong=sa["wrong"],
            s_changed=sa["changed"],
            ceiling_pct=pct(ceil),
            matches_plan_table=cmp,
            within_limits=LR.within_limits(a),
            g1=bool(LR.gate_g1({h: a})[h]),
        )
        print(h, json.dumps(out["samples"][h]), flush=True)
        if h == "mid":
            model = LM.export(
                m,
                pick["threshold"],
                extra=dict(
                    features=list(LFE.FEATURES),
                    trained_on=tr,
                    role="gate model for mid (trained on the other four disclosed samples)",
                    training_hash=hashlib.sha256(
                        "".join(data[t].rows_hash for t in tr).encode()
                    ).hexdigest(),
                ),
            )
            CACHE.mkdir(parents=True, exist_ok=True)
            LM.dump(model, CACHE / "model_gate_mid.json")
    out["wall_s"] = round(time.time() - t0, 1)
    out["prediction_1"] = bool(
        all(all(v["matches_plan_table"].values()) for v in out["samples"].values())
        and all(
            g["n_fits_ok"] and g["counts_ok"] and g["p1_ok"] and g["total_ok"]
            for g in out["g0"].values()
        )
    )
    out["g0_ok"] = out["prediction_1"]
    out["g1_ok"] = bool(all(v["g1"] for v in out["samples"].values()))
    return out


def final() -> Dict[str, Any]:
    data = {s: load(s) for s in SAMPLES}
    oof = oof_scores(list(SAMPLES), data)
    pick = LR.choose_threshold(oof, data)
    X = np.vstack([data[t].X for t in SAMPLES])
    y = np.concatenate([data[t].accepted for t in SAMPLES])
    m = LM.train(X, y)
    model = LM.export(
        m,
        pick["threshold"],
        extra=dict(
            features=list(LFE.FEATURES),
            trained_on=list(SAMPLES),
            role="final model for held4 and the whole layers (trained on the five samples)",
            training_hash=hashlib.sha256(
                "".join(data[t].rows_hash for t in SAMPLES).encode()
            ).hexdigest(),
        ),
    )
    OUT.mkdir(parents=True, exist_ok=True)
    LM.dump(model, OUT / "model_final.json")
    res = dict(
        threshold=pick["threshold"],
        train_share_sum=pick["share_sum"],
        rounds=int(m.n_iter_),
        plan_sizing=dict(
            threshold=LR.SIZING_FINAL_THRESHOLD, rounds=LR.SIZING_FINAL_ROUNDS, json_mb=0.63
        ),
        json_bytes=(OUT / "model_final.json").stat().st_size,
        per_sample_out_of_fold={
            s: {
                k: v
                for k, v in LR.account(LR.skip_vector(oof[s], pick["threshold"]), data[s]).items()
                if k != "share"
            }
            | dict(
                share_pct=pct(
                    LR.account(LR.skip_vector(oof[s], pick["threshold"]), data[s])["share"]
                )
            )
            for s in SAMPLES
        },
    )
    return res


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else ""
    if what == "loso":
        r = loso()
        (OUT / "loso.json").write_text(json.dumps(r, indent=1) + "\n")
        print("prediction 1:", r["prediction_1"], "g1:", r["g1_ok"], "wall_s", r["wall_s"])
    elif what == "final":
        r = final()
        (OUT / "final.json").write_text(json.dumps(r, indent=1) + "\n")
        print(json.dumps({k: v for k, v in r.items() if k != "per_sample_out_of_fold"}))
    else:
        raise SystemExit("usage: learnskip_replay.py loso|final")
