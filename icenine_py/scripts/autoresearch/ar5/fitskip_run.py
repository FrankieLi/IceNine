"""AR5 step 17, phase A: one whole-layer run of arm ``J5`` with the hook on the BFS local
optimizations (script layer; library unchanged).

  RESEED_SAMPLE=<big|grad|mid|held|held3> FITSKIP_CACHE=<dir> uv run python \
      scripts/autoresearch/ar5/fitskip_run.py [--name J5r] [--rule S(0.3,both)] \
      [--cap-min 300] [--smoke N]

Without ``--rule``: the recording run ``J5r``. Arm ``J5`` (step 16: arm ``J`` with ``PH(500)``;
``high500_run.py``) is run unchanged, with the hook of ``fitskip_hook`` in RECORDING mode: after
every BFS local optimization (neighbour fit or revisit) one row goes to
``<cache>/record/<name>_<sample>.jsonl`` (then gzipped). The hook draws no random number and adds
no counted evaluation; its two overlap evaluations per fit are listed apart (``rec_evals``). The
counted run must equal ``J5`` (``fitskip_p1.py`` checks it against the stored run). Aligned draws,
as ``J5``.

With ``--rule`` (``S(t,scope)`` or ``D(d)``): arm ``J5C``, the same run with the rule acting (no
recording). A skipped fit still takes and discards its aligned draw. ``C`` has no live form.

Writes, besides the files of ``order_run.py`` and ``high500_run.py``: ``fitskip.json`` (the hook's
counts). The boundary pass is made after by ``high500_pass.py --name <name>`` (``fitskip_run_all``
and ``fitskip_record_all.sh`` do it). ``held`` and ``held3`` are refused unless ``ar5/fitskip/
frozen.json`` is committed and agrees (``fitskip_common``).
"""

import gzip
import json
import shutil
import sys
from typing import Any, Dict, List, Optional

import fitskip_common as FCM  # noqa: E402  (first: sets the cache, guards held2, held, held3)

import fit_common as FC  # noqa: E402
import fitskip_hook as FH  # noqa: E402
import fitskip_rules as FR  # noqa: E402
import high500_pick as BP  # noqa: E402
import high500_rules as HR  # noqa: E402
import high500_run as H5R  # noqa: E402
import order_draws as OD  # noqa: E402
import order_run as ORUN  # noqa: E402
import reseed_common as RC  # noqa: E402

CFG: Dict[str, Any] = dict(rule=None, name="J5r", record=True)
HOOKS: List[FH.FitHook] = []


class FitskipWrap(H5R.High500Wrap):
    """Arm ``J5``'s wrapper with the fit hook installed on the BFS instance."""

    def __init__(self, bfs: Any) -> None:
        super().__init__(bfs)
        rule: Optional[FR.Rule] = CFG["rule"]
        sym = None
        if rule is not None and rule.form == "D":
            from icenine import ray_indexing as RI

            sym = RI.proper_rotations(bfs.setup.exp_setup.get_sample_symmetry())
        out_path = None
        if CFG["record"]:
            d = FCM.FCACHE / "record"
            d.mkdir(parents=True, exist_ok=True)
            out_path = d / f"{CFG['name']}_{FC.SAMPLE}{CFG['suffix']}.jsonl"
            if out_path.exists():
                raise SystemExit(f"refusing to start: {out_path} already exists")
        HOOKS.append(
            FH.FitHook(
                bfs,
                wrap=self,
                rule=rule,
                sym=sym,
                out_path=out_path,
                book=self.book if ORUN.ORD["aligned"] else None,
                voxel_key=OD.voxel_key,
            )
        )


def main() -> None:
    args = sys.argv[1:]
    name = "J5r"
    rule: Optional[FR.Rule] = None
    rest: List[str] = []
    i = 0
    while i < len(args):
        if args[i] == "--name":
            name = args[i + 1]
            i += 2
        elif args[i] == "--rule":
            rule = FR.rule_by_name(args[i + 1])
            i += 2
        else:
            rest.append(args[i])
            i += 1
    smoke = "--smoke" in rest and int(rest[rest.index("--smoke") + 1]) > 0
    CFG.update(rule=rule, name=name, record=rule is None, suffix="_smoke" if smoke else "")
    BP.H5.update(k=HR.K, aux=())
    sys.argv = ["order_run.py", "--arm", "A", "--name", name] + rest
    ORUN.OrderPick = BP.High500Pick  # type: ignore[misc]
    ORUN.OrderWrap = FitskipWrap  # type: ignore[misc]
    ORUN.main()
    if "--guards-only" in rest:
        return
    out = FC.CACHE / f"{name}_{FC.SAMPLE}{'_smoke' if smoke else ''}"
    RC.write_gz(H5R.LOG, out / "polish_log.json.gz")
    (out / "polish.json").write_text(json.dumps(H5R.summary(), indent=1) + "\n")
    (out / "high500.json").write_text(json.dumps(BP.summary(), indent=1) + "\n")
    hook = HOOKS[0]
    hook.close()
    (out / "fitskip.json").write_text(json.dumps(hook.summary(), indent=1) + "\n")
    rows = FCM.FCACHE / "record" / f"{name}_{FC.SAMPLE}{CFG['suffix']}.jsonl"
    if CFG["record"] and rows.exists():
        with open(rows, "rb") as fi, open(str(rows) + ".gz", "wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as fo:
                shutil.copyfileobj(fi, fo)
        rows.unlink()


if __name__ == "__main__":
    main()
