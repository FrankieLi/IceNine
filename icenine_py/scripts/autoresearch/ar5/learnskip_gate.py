"""AR5 step 18, phase A step 5: the gate of the live run on ``mid`` (g2).

  uv run python scripts/autoresearch/ar5/learnskip_gate.py predict   # before the live run
  uv run python scripts/autoresearch/ar5/learnskip_gate.py compare   # after it

The replay of the gate model (``<cache>/model_gate_mid.json``, trained on the four other disclosed
samples) on the recording ``J5r`` of ``mid`` decides every fit with the features of the recording
(the open-loop replay: each feature is computed from what ``J5`` did). A live run decides with
features computed from what the live run did. Up to the first wrong skip in run order the run of
fits is the same as ``J5``'s (an accepted fit is never skipped; a skipped rejected fit leaves the
expanding voxel's orientation alone), but a skipped fit has no fit end, so a feature that reads a
voxel's kept fit end (own hit, parent hit, parent cost, start over parent) can differ from the
recording's.

``predict``: walks the recording with the live state (``learnskip_features.FeatureState`` and the
live update of a skipped fit) and reports the first fit at which a feature differs from the
open-loop replay by more than 1e-9, and which feature; it runs nothing live and decides nothing.

``compare`` (g2): the live run's rows (``<cache>/J5L_mid/fits.jsonl.gz``) against the replay, fit
by fit, up to the first wrong skip (inclusive) and over the whole run if there is none: voxel, kind,
the 16 features (1e-9), the score (1e-9), the skip, and for a fit that ran its end hit ratio, end
cost, accepted flag and evaluations. Writes ``g2.json``.
"""

import gzip
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import learnskip_features as LFE  # noqa: E402
import learnskip_guard as LG  # noqa: E402
import learnskip_model as LM  # noqa: E402
import learnskip_replay as RP  # noqa: E402
import learnskip_rules as LR  # noqa: E402

OUT = RP.OUT
CACHE = RP.CACHE
TOL = 1e-9


def replay_mid() -> Dict[str, Any]:
    rows = RP.read_rows("mid")
    arr = LFE.rows_to_arrays(rows)
    model = LM.PlainModel.load(CACHE / "model_gate_mid.json")
    scores = model.score_many(arr["X"])
    skip = np.array([model.skips(s) for s in scores])
    w = LR.first_wrong_skip(skip, arr["accepted"])
    return dict(rows=rows, X=arr["X"], accepted=arr["accepted"], scores=scores, skip=skip, w=w)


def predict() -> Dict[str, Any]:
    r = replay_mid()
    rows, X, skip = r["rows"], r["X"], r["skip"]
    w = r["w"]
    last = len(rows) - 1 if w is None else w
    st = LFE.FeatureState(LFE.cubic_rotations())
    first_diff: Optional[Dict[str, Any]] = None
    n_diff = 0
    flips = 0
    model = LM.PlainModel.load(CACHE / "model_gate_mid.json")
    feats_hit: Dict[str, int] = {}
    n_skips = 0
    for k in range(last + 1):
        row = rows[k]
        v = int(row["voxel"])
        x = st.features(
            v,
            str(row["kind"]),
            int(row["n_revisits"]),
            float(row["best"]),
            row["start_R"],
            float(row["start_hit"]),
            float(row["start_cost"]),
            float(row["start_conf"]),
            int(row["parent"]),
            int(row["pos"]),
        )
        d = np.abs(np.array(x) - X[k])
        bad = np.where(d > TOL)[0]
        if model.skips(model.score(x)) != bool(skip[k]):
            flips += 1  # the live decision would differ from the replay's
        if len(bad):
            n_diff += 1
            for j in bad:
                feats_hit[LFE.FEATURES[j]] = feats_hit.get(LFE.FEATURES[j], 0) + 1
            if first_diff is None:
                first_diff = dict(
                    fit=k, features=[LFE.FEATURES[j] for j in bad], max_abs=float(d.max())
                )
        if skip[k]:
            n_skips += 1
            st.update(v, row["start_R"], bool(row["accepted"]), False, None, None)
        else:
            st.update(
                v,
                row["start_R"],
                bool(row["accepted"]),
                bool(row["kept"]),
                float(row["end_hit"]),
                float(row["end_cost"]),
            )
    return dict(
        fits_walked=last + 1,
        first_wrong_skip=w,
        skips_before=n_skips,
        fits_with_a_feature_difference=n_diff,
        decisions_that_would_differ=flips,
        first_difference=first_diff,
        features_hit=feats_hit,
    )


def read_live(path: Path) -> List[Dict[str, Any]]:
    out = []
    with gzip.open(LG.guarded(path), "rt") as f:
        for line in f:
            out.append(json.loads(line))
    return out


def compare() -> Dict[str, Any]:
    r = replay_mid()
    rows, X, skip, scores, w = r["rows"], r["X"], r["skip"], r["scores"], r["w"]
    live = read_live(CACHE / "J5L_mid" / "fits.jsonl.gz")
    last = len(rows) - 1 if w is None else w
    n = min(len(live), last + 1)
    first_bad: Optional[Dict[str, Any]] = None
    worst = dict(feature=0.0, score=0.0, end_hit=0.0, end_cost=0.0)
    n_ok = 0
    for k in range(n):
        a, b = rows[k], live[k]
        why = []
        if int(a["voxel"]) != int(b["voxel"]) or str(a["kind"]) != str(b["kind"]):
            why.append("voxel or kind")
        else:
            df = float(np.abs(np.array(b["x"]) - X[k]).max())
            worst["feature"] = max(worst["feature"], df)
            if df > TOL:
                why.append(f"features (max abs {df:.3g})")
            ds = abs(float(b["score"]) - float(scores[k]))
            worst["score"] = max(worst["score"], ds)
            if ds > TOL:
                why.append(f"score ({ds:.3g})")
            if bool(b["skipped"]) != bool(skip[k]):
                why.append("skip")
            elif not skip[k]:
                dh = abs(float(b["end_hit"]) - float(a["end_hit"]))
                dc = abs(float(b["end_cost"]) - float(a["end_cost"]))
                worst["end_hit"] = max(worst["end_hit"], dh)
                worst["end_cost"] = max(worst["end_cost"], dc)
                if dh > TOL or dc > TOL:
                    why.append(f"end ({dh:.3g}, {dc:.3g})")
                if bool(b["accepted"]) != bool(a["accepted"]):
                    why.append("accepted flag")
                if int(b["evals"]) != int(a["evals"]) + 0:
                    why.append(f"evaluations ({b['evals']} against {a['evals']})")
        if why:
            first_bad = first_bad or dict(fit=k, why=why)
        else:
            n_ok += 1
    res = dict(
        replay_fits=len(rows),
        live_fits=len(live),
        first_wrong_skip=w,
        compared=n,
        equal=n_ok,
        first_difference=first_bad,
        worst=worst,
        tol=TOL,
        g2=bool(first_bad is None and n == last + 1 and n > 0),
    )
    return res


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else ""
    if what == "predict":
        print(json.dumps(predict(), indent=1))
    elif what == "compare":
        r = compare()
        (OUT / "g2.json").write_text(json.dumps(r, indent=1) + "\n")
        print(json.dumps(r, indent=1))
    else:
        raise SystemExit("usage: learnskip_gate.py predict|compare")
