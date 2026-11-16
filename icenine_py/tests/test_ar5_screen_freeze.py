"""AR5 step 21 (task ar5u): the freeze lists of the screen step and the appliers' refusal before the
freeze is committed (plan "AR5 step 21", phase A step 2)."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
AR5 = ROOT / "scripts" / "autoresearch" / "ar5"
sys.path.insert(0, str(AR5))

import screen_freeze as SF  # noqa: E402


def test_freeze_lists_exist_and_constants_match_rules() -> None:
    for name in SF.OWN:
        assert (AR5 / name).exists(), name
    cur = SF.current()
    assert "scripts/autoresearch/ar5/screen_rules.py" in cur["files"]
    assert "scripts/autoresearch/ar5/cmarepair_rules.py" in cur["files"]
    assert "icenine/reconstructor.py" in cur["library"]
    c = SF.CONSTANTS
    assert c["screen_max_evals"] == 1 and c["evals_per_screen"] == 2
    assert c["repair_cma_max_evals"] == 120 and c["wide_deg"] == 1.0 and c["ratio_max"] == 0.95
    # no data file of held5 or held2 is hashed (that would read it): scripts and library only
    assert all(
        k.startswith(("scripts/", "icenine/")) or k in ("pyproject.toml", "uv.lock")
        for part in cur.values()
        for k in part
    )
    assert not any("held5" in k for k in cur["files"])


def test_the_appliers_check_the_freeze_first() -> None:
    for src in ("screen_pass.py", "screen_run.py"):
        text = (AR5 / src).read_text()
        assert "check_or_exit()" in text, src
        assert text.index("import screen_common") < text.index("check_or_exit()"), src


def test_the_freeze_does_not_depend_on_step_20s_check() -> None:
    text = (AR5 / "screen_common.py").read_text()
    assert "cmarepair_freeze" not in text.replace("``cmarepair_freeze``", "")
    assert "import cmarepair" not in text


def test_held5_is_built_as_held_with_generator_seed_5(monkeypatch) -> None:  # type: ignore
    import screen_samples as SS

    assert SS.SEEDS["held5"] == 5 and SS.SEEDS["held3"] == 3 and SS.SEEDS["held4"] == 4
    src = SS.patched_source("make", "held5")
    assert "SEED = 5" in src and '{"held": 1, "held2": 2}' not in src
    assert 'K = "held5"' in SS.patched_source("stack", "held5")
    monkeypatch.setenv("AR3_SAMPLE", "held5")
    monkeypatch.setenv("RESEED_SAMPLE", "held5")
    assert SS.wanted() and SS.which_sample() == "held5"
