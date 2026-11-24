"""AR5 step 28, phase B: one whole-layer run of the script arm ``PXFl`` (step 22's ``J5Dl`` with
rule ``PXF``: arm ``J5D`` plus the small pair sample first, without ``G2b``, with the library's own
draws).

  RESEED_SAMPLE=mid PXF_CACHE=<cache> uv run python scripts/autoresearch/ar5/pxf_hl.py

The script arm that the library run with ``RayIndexMaxPairs 50000`` and ``RayIndexFallbackFull 1``
must equal before the boundary pass. The phase A arm ``J5D`` (rule ``PXF``, ``G2b``, aligned
draws) is not the comparison: the library has neither ``G2b`` nor aligned draws (as in step 22).
Writes ``<cache>/PXFl_<sample>/`` (the files of ``defer_hl.py``, plus ``pxf_index.json.gz`` and
``pxf_steps.json.gz``) and ``PXFl.done``. ``mid`` only. Goes through ``pxf_common`` (the freeze must
agree: run it BEFORE the library keys exist).
"""

import gzip
import json
import sys

import pxf_common as PC  # noqa: E402  (first: guards, freeze, cache)
import pxf_launch as PL  # noqa: E402
import pxf_rules as PR  # noqa: E402

import defer_pick as DP  # noqa: E402
import fit_common as FC  # noqa: E402
import high500_pick as HP  # noqa: E402
import high500_rules as HR  # noqa: E402
import order_run as ORUN  # noqa: E402
import reseed_common as RC  # noqa: E402
import rescue_run as RRUN  # noqa: E402


def main() -> None:
    args = sys.argv[1:]
    assert "--smoke" not in args, "no smoke run: a truncated map does not end quickly on mid"
    assert FC.SAMPLE == "mid", "the library check is made on mid"
    from icenine import ray_indexing as RI

    ctl = PR.Controller()
    steps: list = []
    RI.index_voxel = PR.make_index_voxel(RI.index_voxel, ctl, "pxf")  # type: ignore[assignment]
    DP.DeferPick.e1_step = PR.make_e1_step(  # type: ignore[assignment]
        DP.DeferPick.e1_step, ctl, steps, PL._hit_of, PL._threshold_of
    )
    DP.DEFER["on"] = True
    HP.H5.update(k=HR.K, aux=())
    sys.argv = ["order_run.py", "--arm", "Al", "--name", "PXFl"] + args
    ORUN.OrderPick = DP.DeferPick  # type: ignore[misc]
    ORUN.OrderWrap = RRUN.RescueWrap  # type: ignore[misc]
    ORUN.main()
    if "--guards-only" in args:
        return
    out = FC.CACHE / f"PXFl_{FC.SAMPLE}"
    RC.write_gz(RRUN.LOG, out / "rescue_log.json.gz")
    (out / "rescue.json").write_text(json.dumps(RRUN.summary(), indent=1) + "\n")
    (out / "stallhit.json").write_text(json.dumps(HP.summary(), indent=1) + "\n")
    (out / "defer.json").write_text(json.dumps(DP.summary(), indent=1) + "\n")
    RC.write_gz(DP.LOG, out / "defer_steps.json.gz")
    with gzip.open(out / "pxf_index.json.gz", "wt") as f:
        json.dump(dict(summed=PR.summed(ctl.log), rows=ctl.log), f)
    with gzip.open(out / "pxf_steps.json.gz", "wt") as f:
        json.dump(dict(sample=PC.SAMPLE, rows=steps), f)
    (out / "PXFl.done").write_text("done\n")


if __name__ == "__main__":
    main()
