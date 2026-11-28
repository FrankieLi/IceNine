"""AR5 step 18, phase A: one whole-layer run of arm ``J5L`` (arm ``J5`` with the learned skip of
the BFS local optimizations; script layer, library unchanged).

  RESEED_SAMPLE=<big|grad|mid|held|held3|held4> LEARNSKIP_CACHE=<dir> uv run python \
      scripts/autoresearch/ar5/learnskip_run.py --model <model.json> [--name J5L] \
      [--cap-min 300] [--smoke N]

Arm ``J5`` (step 16: arm ``J`` with ``PH(500)``; ``high500_run.py``) is run unchanged with the hook
of ``learnskip_hook`` acting: every neighbour fit and revisit is checked (one counted evaluation of
the start), scored by the exported trees and skipped if the score is under the model's threshold
(see ``learnskip_hook``). A skipped fit takes and discards its aligned draw. Aligned draws, as
``J5``. Writes, besides the files of ``order_run.py`` and ``high500_run.py``: ``learnskip.json``
(the hook's counts, the model and feature wall times) and ``fits.jsonl.gz`` (one row per fit). The
boundary pass is made after by ``learnskip_pass.py --name <name>``. ``held``, ``held3`` and
``held4`` are refused unless ``ar5/learnskip/frozen.json`` is committed and agrees
(``learnskip_common``).
"""

import gzip
import json
import shutil
import sys
from pathlib import Path
from typing import Any, Dict, List

import learnskip_common as LCM  # noqa: E402  (first: sets the cache, guards held2, held*)

import fit_common as FC  # noqa: E402
import high500_pick as BP  # noqa: E402
import high500_rules as HR  # noqa: E402
import high500_run as H5R  # noqa: E402
import learnskip_hook as LH  # noqa: E402
import learnskip_model as LM  # noqa: E402
import order_draws as OD  # noqa: E402
import order_run as ORUN  # noqa: E402
import reseed_common as RC  # noqa: E402

CFG: Dict[str, Any] = dict(model=None, name="J5L", suffix="")
HOOKS: List[LH.LearnHook] = []


class LearnWrap(H5R.High500Wrap):
    """Arm ``J5``'s wrapper with the learned-skip hook installed on the BFS instance."""

    def __init__(self, bfs: Any) -> None:
        super().__init__(bfs)
        from icenine import ray_indexing as RI

        sym = RI.proper_rotations(bfs.setup.exp_setup.get_sample_symmetry())
        d = LCM.LCACHE / "learn"
        d.mkdir(parents=True, exist_ok=True)
        out_path = d / f"{CFG['name']}_{FC.SAMPLE}{CFG['suffix']}.jsonl"
        if out_path.exists():
            raise SystemExit(f"refusing to start: {out_path} already exists")
        HOOKS.append(
            LH.LearnHook(
                bfs,
                CFG["model"],
                sym,
                wrap=self,
                out_path=out_path,
                book=self.book if ORUN.ORD["aligned"] else None,
                voxel_key=OD.voxel_key,
            )
        )


def main() -> None:
    args = sys.argv[1:]
    name = "J5L"
    model_path = ""
    rest: List[str] = []
    i = 0
    while i < len(args):
        if args[i] == "--name":
            name = args[i + 1]
            i += 2
        elif args[i] == "--model":
            model_path = args[i + 1]
            i += 2
        else:
            rest.append(args[i])
            i += 1
    if not model_path:
        raise SystemExit("learnskip_run.py: --model <model.json> is required")
    smoke = "--smoke" in rest and int(rest[rest.index("--smoke") + 1]) > 0
    CFG.update(
        model=LM.PlainModel.load(Path(model_path)), name=name, suffix="_smoke" if smoke else ""
    )
    BP.H5.update(k=HR.K, aux=())
    sys.argv = ["order_run.py", "--arm", "A", "--name", name] + rest
    ORUN.OrderPick = BP.High500Pick  # type: ignore[misc]
    ORUN.OrderWrap = LearnWrap  # type: ignore[misc]
    ORUN.main()
    if "--guards-only" in rest:
        return
    out = FC.CACHE / f"{name}_{FC.SAMPLE}{'_smoke' if smoke else ''}"
    RC.write_gz(H5R.LOG, out / "polish_log.json.gz")
    (out / "polish.json").write_text(json.dumps(H5R.summary(), indent=1) + "\n")
    (out / "high500.json").write_text(json.dumps(BP.summary(), indent=1) + "\n")
    hook = HOOKS[0]
    hook.close()
    (out / "learnskip.json").write_text(json.dumps(hook.summary(), indent=1) + "\n")
    rows = LCM.LCACHE / "learn" / f"{name}_{FC.SAMPLE}{CFG['suffix']}.jsonl"
    if rows.exists():
        with open(rows, "rb") as fi, open(out / "fits.jsonl.gz", "wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as fo:
                shutil.copyfileobj(fi, fo)
        rows.unlink()


if __name__ == "__main__":
    main()
