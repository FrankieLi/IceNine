"""AR5 step 12, review fix 1 of the second review: the run and pass scripts refuse to write when
the real path of an output lies outside their own cache (e.g. a whole-folder link into another
worktree's cache)."""

import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
A5 = ROOT / "scripts" / "autoresearch" / "ar5"
sys.path.insert(0, str(A5))

import stallhit_cache_guard as SCG  # noqa: E402

PASS_FILES = [
    "passes_mid",
    "passes_mid/D4_Hl",
    "passes_mid/before_Hl.npy",
    "passes_mid/after_Hl.npy",
]


def _caches(tmp_path: Path) -> "tuple[Path, Path]":
    cache = tmp_path / "cache"
    other = tmp_path / "other_worktree_cache"
    (cache / "logs").mkdir(parents=True)
    (other / "passes_mid" / "D4_A").mkdir(parents=True)
    (other / "passes_mid" / "after_A.npy").write_bytes(b"")
    return cache, other


def test_a_whole_folder_link_out_of_the_cache_is_refused(tmp_path: Path) -> None:
    cache, other = _caches(tmp_path)
    (cache / "passes_mid").symlink_to(other / "passes_mid")
    bad = SCG.outside(cache, PASS_FILES)
    assert len(bad) == 4 and all("outside the cache" in b for b in bad)
    with pytest.raises(SystemExit) as e:
        SCG.require_inside(cache, ["Hl_mid"] + PASS_FILES)
    assert e.value.code == 2
    assert SCG.main(["--cache", str(cache), "logs"] + PASS_FILES) == 2


def test_a_real_folder_with_per_entry_links_accepts_new_names(tmp_path: Path) -> None:
    cache, other = _caches(tmp_path)
    (cache / "passes_mid").mkdir()
    (cache / "passes_mid" / "D4_A").symlink_to(other / "passes_mid" / "D4_A")
    (cache / "passes_mid" / "after_A.npy").symlink_to(other / "passes_mid" / "after_A.npy")
    assert SCG.outside(cache, ["Hl_mid", "logs"] + PASS_FILES) == []
    SCG.require_inside(cache, ["Hl_mid"] + PASS_FILES)  # no exit
    # writing into a linked stored entry is refused
    assert len(SCG.outside(cache, ["passes_mid/D4_A", "passes_mid/D4_A/pass.json"])) == 2
    assert len(SCG.outside(cache, ["passes_mid/after_A.npy"])) == 1


def test_absolute_paths_dangling_links_and_a_linked_cache(tmp_path: Path) -> None:
    cache, other = _caches(tmp_path)
    assert SCG.outside(cache, [cache / "LibH_mid"]) == []
    assert len(SCG.outside(cache, [other / "LibH_mid"])) == 1
    assert len(SCG.outside(cache, ["../other_worktree_cache/x"])) == 1
    (cache / "LibH_mid").symlink_to(other / "not_there_yet")  # dangling link out
    assert len(SCG.outside(cache, ["LibH_mid", "LibH_mid/run.json"])) == 2
    # the cache itself may be a link: its real path is the root
    alias = tmp_path / "cache_alias"
    alias.symlink_to(cache)
    assert SCG.outside(alias, ["Hl_mid", "logs", cache / "x"]) == []


def test_run_all_dry_run_refuses_a_linked_pass_folder_first(tmp_path: Path) -> None:
    cache, other = _caches(tmp_path)
    (other / "passes_big").mkdir()
    (cache / "passes_big").symlink_to(other / "passes_big")
    env = {k: v for k, v in os.environ.items() if k not in ("STALLHIT_CACHE", "VIRTUAL_ENV")}
    cmd = [str(A5 / "stallhit_run_all.sh"), "--dry-run", "--smoke-cap", "1", "--cache", str(cache)]
    r = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=300)
    assert r.returncode == 2
    assert "outside the cache" in r.stderr and "passes_big" in r.stderr
    assert sorted(p.name for p in cache.iterdir()) == ["logs", "passes_big"]
    assert list((cache / "logs").iterdir()) == [] and list((other / "passes_big").iterdir()) == []
    # the same dry run without the link is not refused by this guard (it may refuse later)
    (cache / "passes_big").unlink()
    r = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=300)
    assert "outside the cache" not in r.stderr


def test_the_runners_and_unfrozen_run_scripts_call_the_guard() -> None:
    pb = (A5 / "stallhit_phase_b.sh").read_text()
    assert pb.index("stallhit_cache_guard.py") < pb.index('mkdir -p "$logs"')
    assert "passes_mid/D4_Hl" in pb and "LibH_mid" in pb
    ra = (A5 / "stallhit_run_all.sh").read_text()
    assert ra.index("stallhit_cache_guard.py") < ra.index("# ---- guard 1:")
    for name in ("stallhit_hl.py", "stallhit_lib_run.py"):
        assert "SCG.require_inside(" in (A5 / name).read_text()
