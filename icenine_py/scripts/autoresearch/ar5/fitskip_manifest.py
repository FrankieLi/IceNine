"""AR5 step 17: write the tracked manifest of the recording (sha256 and row count per file).

  uv run python scripts/autoresearch/ar5/fitskip_manifest.py

Reads the files of ``<cache>/record/`` (the recording stays in the untracked cache) and writes
``ar5/fitskip/manifest.json`` with the fields, units, which fields are recording-only, the row
counts per sample and the checksums. No run; nothing of ``held2`` or ``held4`` is touched.
"""

import gzip
import hashlib
import json
from pathlib import Path

ICE = Path(__file__).resolve().parents[3]
FS = ICE / "benchmarks" / "autoresearch" / "ar5" / "fitskip"
SAMPLES = ("held3", "held", "mid", "grad", "big")


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    rec = FS / "cache" / "record"
    files = {}
    for s in SAMPLES:
        rp, fp = rec / f"J5r_{s}.jsonl.gz", rec / f"final_J5r_{s}.json"
        with gzip.open(rp, "rt") as f:
            n = sum(1 for x in f if x.strip())
        fin = json.loads(fp.read_text())
        files[s] = dict(
            rows_file=rp.name,
            rows=n,
            rows_sha256=sha(rp),
            rows_bytes=rp.stat().st_size,
            final_file=fp.name,
            final_sha256=sha(fp),
            final_before=fin["final_before"],
            final_after=fin["final_after"],
            voxels_with_kept_fit=fin["n_voxels_with_kept_fit"],
        )
    out = dict(
        what="recording J5r: one row per BFS local optimization of arm J5 "
        "(neighbour fit or revisit)",
        location="untracked cache of the worktree that ran it: icenine_py/benchmarks/autoresearch/"
        "ar5/fitskip/cache/record/ (not in git; archive it before the worktree is removed)",
        counted_fields=[
            "evals",
            "runbest",
            "start_cost",
            "end_R",
            "end_hit",
            "end_cost",
            "end_conf",
        ],
        recording_only_fields=["start_hit", "start_conf", "hit20", "rec_evals"],
        recording_only_note="two overlap evaluations per fit (start_hit, start_conf, hit20) are "
        "made "
        "outside the counted run and are in no stats; a live rule that uses start_hit pays one "
        "counted evaluation; hit20 has no live form (needs a library change)",
        fields=dict(
            i="run-order row number from 0",
            exp="expansion number",
            pos="fits already made in this expansion",
            voxel="voxel index",
            kind="neighbour or revisit",
            parent="voxel that passed its orientation on, -1 if unknown",
            best="expansion's running best hit ratio before the fit, ratio in [0, 1]",
            start_R="start rotation matrix, row-major 3x3, 12 decimals",
            start_hit="hit ratio of the start, ratio in [0, 1] (recording-only)",
            start_conf="confidence of the start (recording-only)",
            start_cost="cost of the start (first evaluation of the fit)",
            runbest="running minimum cost after each evaluation, 9 decimals, 61 values",
            n_calls="cost calls seen in the fit (61)",
            hit20="hit ratio of the running best after 20 evaluations (recording-only)",
            evals="counted evaluations of the fit (61)",
            end_R="end rotation matrix, row-major 3x3, 12 decimals",
            end_hit="hit ratio of the end, ratio in [0, 1]",
            end_cost="cost of the end",
            end_conf="confidence of the end",
            accepted="hit ratio over 0.9 times max(best, hit)",
            kept="the hook's model of the push rule (a first fit is kept; a revisit is kept if "
            "accepted or its confidence is not under the old one)",
            n_revisits="revisit count of the voxel at this fit",
            rec_evals="recording-only evaluations of this row (2)",
        ),
        final_flags="final_J5r_<s>.json: row numbers i that are the last kept fit of their voxel "
        "and whose end matches the run's map (before the pass) or the map after the pass, max "
        "absolute difference of the nine entries at most 1e-6",
        files=files,
    )
    (FS / "manifest.json").write_text(json.dumps(out, indent=1) + "\n")
    print("ok", {s: files[s]["rows"] for s in SAMPLES})


if __name__ == "__main__":
    main()
