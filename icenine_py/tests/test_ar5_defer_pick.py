"""AR5 step 22 (task ar5v): the live pick's pause and resume of one MC refinement, with a scripted
refinement and a fake tracer (no sample, no image, no reconstructor). Run in a subprocess: importing
``defer_pick`` sets the study's caches and replaces the older freeze checks in its process, which
must not leak into other tests. The freeze check is switched off in that process only (this test
checks the pick's logic, not the freeze)."""

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
AR5 = ROOT / "scripts" / "autoresearch" / "ar5"

CODE = r'''
import sys
sys.path.insert(0, {ar5!r})
import defer_freeze as DF
DF.problems = lambda: ""
import numpy as np
import defer_pick as DP
import mix_common as MC
import stallhit_rules as SR

class FT:
    """A fake tracer with the fields the pick reads."""
    def __init__(self):
        self.on = False; self.hit_fn = None; self.stop_kind = ""; self.events = {{}}
        self.n = 0; self.best = float("inf"); self.trace = []; self.best_R = None
        self.checks = []; self.continued = 0
    def begin(self, k, low=None, aux=()):
        self.on = True; self.n = 0; self.best = float("inf"); self.trace = []; self.best_R = None
        self.checks = []; self.continued = 0; self.stop_kind = ""; self.events = {{}}

class Rec:
    def evaluate_overlap(self, R, v, ph):
        return float(np.asarray(R).ravel()[0])
    def _make_local_cost_fn(self):
        class F:
            def evaluate(self_, orientation, voxel_vertices, phase_index):
                class O: cost = 1.25
                return O()
        return F()
class Bfs:
    def _ratios(self, info):
        return (info,)

HIT = {{"at_check": 0.6}}
SEEN = []
def refine(arm, rec, lo, R, c, v, ph):
    """A refinement: a gain at 100 (cost 2.0), a stall check at 600 (hit HIT['at_check']); if it
    goes on, a gain at 1000 (cost 1.5) and the end at 1500 evaluations."""
    tr = DP._BASE.ptracer
    tr.n = 100; tr.best = 2.0; tr.trace = [[100, 2.0]]
    tr.best_R = np.array([HIT["at_check"] + 0.0] * 9)
    tr.n = 600
    SEEN.append(tr.n)
    h = tr.hit_fn(tr.best_R)
    if h >= 0.95:
        tr.checks.append([600, h, True]); tr.stop_kind = "high"
        raise SR.StallHitStop()
    tr.checks.append([600, h, False]); tr.continued += 1
    tr.n = 1000; tr.best = 1.5; tr.trace.append([1000, 1.5]); tr.best_R = np.array([0.99] * 9)
    return tr.best_R, 1.5, 1500

MC._orig_refine_arm = refine
DP._BASE.ptracer = FT()
DP.ORUN.OrderPick.tracer = DP._BASE.ptracer
p = DP.DeferPick.__new__(DP.DeferPick)
p.rec = Rec(); p.bfs = Bfs(); p.orig_lo = None; p.cur = 3; p.reseed = False
R0 = np.zeros(9)

# 1. a first run that may pause: paused at the check that goes on, 600 evaluations plus the check
f = p._fit_m(p.rec, 2, R0, 3.0, None, 0, True)
assert f["paused"] and f["evals"] == 601 and abs(f["hit"] - 0.6) < 1e-12 and f["cost"] == 2.0, f
assert f["state"]["n"] == 600 and f["state"]["trace"] == [[100, 2.0]]
assert DP.STATS["pauses"] == 1 and DP.STATS["rerun_evals"] == 0
# 2. not allowed to pause (the last candidate): the refinement runs through, J5's count
g = p._fit_m(p.rec, 4, R0, 3.0, None, 0, False)
assert not g.get("paused") and g["evals"] == 1501 and g["cost"] == 1.5, g
# 3. the resume: the first part is re-run and checked, not counted; the fit equals the unpaused one
r = p._fit_m(p.rec, 2, R0, 3.0, None, 0, False, f["state"])
assert r["evals"] == g["evals"] and r["cost"] == g["cost"] and abs(r["hit"] - g["hit"]) < 1e-12
assert np.array_equal(r["R"], g["R"])
assert DP.STATS["rerun_evals"] == 601 and DP.STATS["rerun_checked"] == 1
# 4. a refinement that would stop at its check is never paused (hit at 0.95 or more)
HIT["at_check"] = 0.97
s = p._fit_m(p.rec, 0, R0, 3.0, None, 0, True)
assert not s.get("paused") and s["evals"] == 600 + 1 + 0, s
# 5. a re-run that differs from the first run is a stop (mismatch)
HIT["at_check"] = 0.6
f2 = p._fit_m(p.rec, 1, R0, 3.0, None, 0, True)
bad = dict(f2["state"]); bad["hit"] = 0.61
try:
    p._fit_m(p.rec, 1, R0, 3.0, None, 0, False, bad)
    raise SystemExit("a differing re-run was accepted")
except DP.DeferRerunMismatch:
    pass
print("PICK-OK")
'''


def test_pause_and_resume_of_one_refinement_with_a_scripted_refinement() -> None:
    env = {"RESEED_SAMPLE": "mid", "PATH": "/usr/bin:/bin:/usr/local/bin:/opt/homebrew/bin"}
    import os

    env = dict(os.environ, RESEED_SAMPLE="mid")
    env.pop("VIRTUAL_ENV", None)
    r = subprocess.run(
        [sys.executable, "-c", CODE.format(ar5=str(AR5))],
        capture_output=True,
        text=True,
        env=env,
        cwd=str(ROOT),
    )
    assert "PICK-OK" in r.stdout, r.stdout[-500:] + r.stderr[-1500:]
