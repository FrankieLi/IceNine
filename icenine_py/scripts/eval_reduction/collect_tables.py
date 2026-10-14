"""Merge the generated tables of a task into one file with a name prefix, for
scripts/dev/sync_doc_tables.py.

  uv run python scripts/eval_reduction/collect_tables.py --root benchmarks/eval_reduction/A \
      --prefix evalA_ --out benchmarks/eval_reduction/A/doc_tables_all.md
"""

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "common"))
import doc_tables as DT  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    root = Path(a.root)
    merged = {}
    for f in sorted(root.glob("*/tables*.md")) + sorted(root.glob("*/*/tables*.md")):
        if "doc_tables_all" in f.name:
            continue
        for name, text in DT.read_tables(f).items():
            key = a.prefix + name
            if key in merged:
                raise SystemExit(f"duplicate table name {key} ({f})")
            merged[key] = text
    DT.write_tables(Path(a.out), merged)
    print(f"{len(merged)} tables -> {a.out}")


if __name__ == "__main__":
    main()
