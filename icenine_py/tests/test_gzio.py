"""scripts/common/gzio.py: transparent reads of plain or gzip-compressed result files."""

import gzip
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "common"))

import gzio  # noqa: E402


def _gz(path: Path, obj: object) -> None:
    with gzip.open(str(path) + ".gz", "wt", encoding="utf-8") as fh:
        json.dump(obj, fh)


def test_load_plain_and_gz(tmp_path: Path) -> None:
    (tmp_path / "a.json").write_text(json.dumps({"x": 1}))
    _gz(tmp_path / "b.json", {"y": 2})
    assert gzio.load_json(tmp_path / "a.json") == {"x": 1}
    assert gzio.load_json(tmp_path / "b.json") == {"y": 2}
    assert gzio.load_json(tmp_path / "b.json.gz") == {"y": 2}
    assert gzio.resolve(tmp_path / "b.json").name == "b.json.gz"
    assert gzio.exists(tmp_path / "b.json") and not gzio.exists(tmp_path / "c.json")
    with pytest.raises(FileNotFoundError):
        gzio.load_json(tmp_path / "c.json")


def test_plain_wins_over_gz(tmp_path: Path) -> None:
    (tmp_path / "a.json").write_text(json.dumps({"new": True}))
    _gz(tmp_path / "a.json", {"new": False})
    assert gzio.load_json(tmp_path / "a.json") == {"new": True}


def test_glob_paths_merges_and_dedupes(tmp_path: Path) -> None:
    (tmp_path / "r_v1.json").write_text("{}")
    _gz(tmp_path / "r_v2.json", {})
    _gz(tmp_path / "r_v1.json", {})
    (tmp_path / "other.txt").write_text("")
    got = gzio.glob_paths(tmp_path / "r_v*.json")
    assert [p.name for p in got] == ["r_v1.json", "r_v2.json"]


def test_write_json_like_keeps_compression(tmp_path: Path) -> None:
    _gz(tmp_path / "a.json", {"v": 0})
    out = gzio.write_json_like(tmp_path / "a.json", '{"v": 1}\n')
    assert out.name == "a.json.gz" and not (tmp_path / "a.json").exists()
    assert gzio.load_json(tmp_path / "a.json") == {"v": 1}
    first = out.read_bytes()
    gzio.write_json_like(tmp_path / "a.json", '{"v": 1}\n')
    assert out.read_bytes() == first  # mtime 0, no name: reproducible bytes
    out2 = gzio.write_json_like(tmp_path / "b.json", "{}\n")
    assert out2.name == "b.json" and out2.read_text() == "{}\n"
