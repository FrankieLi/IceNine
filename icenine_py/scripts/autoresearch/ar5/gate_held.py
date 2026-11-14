"""AR5 step 7, Phase A step 4: the replay of the skip on ``held`` (the test sample).

  uv run python scripts/autoresearch/ar5/gate_held.py

Run only after ``ar5/gate/frozen.json`` is committed. It records the time it first opens a file of
``held`` and the SHA-256 of the files it reads (``held_files.json``), then replays the frozen
``h_s`` (and, as a what-if, the 0.90 to 0.99 curve) on the stored pass of ``S11v2`` and scores it
against the truth (``replay_held.json``). If the stored tries have no ``old_score``, it says so and
writes nothing else (``held`` is then tested live only). No value is tuned here.
"""

import datetime
import json
from typing import Any, Dict

import gate_common as GC
import gate_replay as RP
import gate_rules as GR
import numpy as np

SAMPLE, RUN = "held", "S11v2"


def main() -> None:
    frozen = json.loads((GC.GATE / "frozen.json").read_text())
    assert frozen["h_s"] == GR.H_S and not frozen["held_read"]
    t0 = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    files: Dict[str, Any] = {k: GC.sha256(p) for k, p in GC.pass_files(SAMPLE, RUN).items()}
    s8d = GC.HELD_S8 / "passes_held"
    s8_files = dict(
        run_json=GC.HELD_S8 / "S8_held" / "run.json",
        pass_json=s8d / "D4_S8" / "pass.json",
        after=s8d / "after_S8.npy",
    )
    files.update({f"s8/{k}": GC.sha256(p) for k, p in s8_files.items()})
    (GC.GATE / "held_files.json").write_text(
        json.dumps(dict(first_opened_at=t0, h_s_frozen=GR.H_S, sha256=files), indent=1) + "\n"
    )
    keys = list(np.load(GC.pass_files(SAMPLE, RUN)["tries"]).files)
    if "old_score" not in keys:
        (GC.GATE / "replay_held.json").write_text(
            json.dumps(dict(old_score_stored=False, keys=keys), indent=1) + "\n"
        )
        print("no old_score stored: held is tested live only")
        return
    P = GC.load_pass(SAMPLE, RUN)
    T = GC.Truth(SAMPLE)
    s8 = GC.s8_from_maps(T, s8_files["run_json"], s8_files["pass_json"], s8_files["after"])
    r = RP.run_rows(SAMPLE, RUN, P, T.err, s8, RP.H_CURVE)
    r["pass_json_matches_tries"] = r["pass_evals_full"] == r["pass_json_evals"]
    (GC.GATE / "replay_held.json").write_text(
        json.dumps(dict(h_s=GR.H_S, old_score_stored=True, first_opened_at=t0, row=r), indent=1)
        + "\n"
    )
    x = [b for b in r["by_h"] if abs(b["h"] - GR.H_S) < 1e-9][0]
    print(
        f"held S11v2 pass {r['pass_evals_full']:,} cut@0.95 {100 * x['cut']:.1f}% wrong "
        f"{r['full']['after']['wrong']} -> {x['wrong']} fixed/broken {x['fixed_vs_full']}/"
        f"{x['broken_vs_full']} median {r['full']['after']['median_err_deg']:.5f} -> "
        f"{x['median_err_deg']:.5f} total {x['total_with_skip']:,} vs full {x['total_full']:,}"
    )
    print("rule at 0.95:", x["rule"])
    print("map check:", r["map_check"])


if __name__ == "__main__":
    main()
