"""Mechanical half of the claims audit: do the numbers in a doc section occur in the sources?

Usage: audit_numbers.py --doc FILE.md --section "Heading text" --sources FILE_OR_GLOB... [--strict]

Numbers (decimals, percentages, ranges, "a / b" pairs, integers >= 10) are extracted from the
section; a number matches when some source value, or that value x100 or /100, equals it at the
printed precision. Unmatched numbers are listed with line numbers. Exit 0 unless --strict.
"""

import argparse
import glob
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

TEXT_SUFFIXES = {".txt", ".md", ".json", ".csv", ".log", ".tsv", ".yaml", ".yml", ".out"}
MAX_ARRAY = 200
RANGE_PCT_RE = re.compile(r"\s*[–-]\s*\d+(?:\.\d+)?%")
NUM_RE = re.compile(r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")
DATE_RE = re.compile(r"\b\d{4}-\d{2}-\d{2}\b")
URL_RE = re.compile(r"https?://\S+")
CODE_SPAN_RE = re.compile(r"`([^`]*)`")
PATH_LIKE_RE = re.compile(r"[/\\]|\.(?:py|json|txt|md|npz|npy|png|sh|cpp|h|config|csv|log)\b")
PAIR_RE = re.compile(r"(?<![\w.])(\d+)\s*/\s*(\d+)(?!\w|\.\d)")
TOKEN_RE = re.compile(r"(?<![\w.])(\d{1,3}(?:,\d{3})+|\d+)(?:\.(\d+))?(%|k)?(?![\w])")
LIST_MARKER_RE = re.compile(r"^\s*(?:[-*+]\s+)?\d+[.)]\s")
TABLE_SEP_RE = re.compile(r"^\s*\|?[\s:|-]+\|?\s*$")
CONTEXT_SKIP_RE = re.compile(
    r"(?:Sigma|Σ)\s*(?:<=|≤|<|=)?\s*$|Q[_-]?max\s*=?\s*$|(?:Step|Stage|Task|Section|Phase|PR|#)\s*$",
    re.IGNORECASE,
)


def find_section(lines: Sequence[str], heading: str) -> Tuple[int, int]:
    """(start, end) line indexes (0-based, end exclusive) of the section whose heading line
    contains ``heading``; it runs to the next heading of the same or a higher level."""
    for i, line in enumerate(lines):
        m = re.match(r"^(#+)\s+(.*)$", line)
        if m and heading in m.group(2):
            level = len(m.group(1))
            for j in range(i + 1, len(lines)):
                m2 = re.match(r"^(#+)\s", lines[j])
                if m2 and len(m2.group(1)) <= level:
                    return i, j
            return i, len(lines)
    raise ValueError(f"section not found: {heading!r}")


def _clean(line: str) -> str:
    def span(m: "re.Match[str]") -> str:
        return " " if PATH_LIKE_RE.search(m.group(1)) else m.group(0)

    line = CODE_SPAN_RE.sub(span, line)
    line = URL_RE.sub(" ", line)
    return DATE_RE.sub(" ", line)


def extract_numbers(lines: Sequence[str], offset: int = 0) -> List[Dict[str, Any]]:
    """Doc numbers as dicts: line (1-based), text, value, decimals, kind ('num' or 'pair')."""
    found: List[Dict[str, Any]] = []
    for idx, raw in enumerate(lines):
        lineno = offset + idx + 1
        if TABLE_SEP_RE.match(raw) or re.match(r"^#+\s", raw) and DATE_RE.search(raw):
            continue
        line = _clean(raw)
        if LIST_MARKER_RE.match(line):
            line = re.sub(r"^(\s*(?:[-*+]\s+)?)\d+[.)]", r"\1", line)
        for m in PAIR_RE.finditer(line):
            a, b = int(m.group(1)), int(m.group(2))
            if b >= 10 and not CONTEXT_SKIP_RE.search(line[: m.start()]):
                found.append({"line": lineno, "text": m.group(0), "kind": "pair", "a": a, "b": b})
        line = PAIR_RE.sub(" ", line)
        for m in TOKEN_RE.finditer(line):
            if CONTEXT_SKIP_RE.search(line[: m.start()]):
                continue
            whole = m.group(1).replace(",", "")
            frac = m.group(2)
            suffix = m.group(3)
            if frac is None and suffix != "%" and suffix != "k":
                if int(whole) < 10 or (1900 <= int(whole) <= 2100):
                    continue
            decimals = len(frac) if frac else 0
            range_pct = RANGE_PCT_RE.match(line[m.end() :]) is not None
            value = float(whole + ("." + frac if frac else ""))
            if suffix == "k":
                value *= 1000
                decimals = -3
            found.append(
                {
                    "line": lineno,
                    "text": m.group(0),
                    "kind": "num",
                    "value": value,
                    "decimals": decimals,
                    "pct": suffix == "%" or range_pct,
                }
            )
    return found


def _walk(obj: Any, out: List[float]) -> None:
    if isinstance(obj, bool):
        return
    if isinstance(obj, (int, float)):
        out.append(float(obj))
    elif isinstance(obj, dict):
        _walk_group(list(obj.values()), out)
        for v in obj.values():
            _walk(v, out)
    elif isinstance(obj, list):
        if len(obj) > MAX_ARRAY and all(isinstance(v, (int, float)) for v in obj):
            return  # raw data array, not a reported number
        _walk_group(obj, out)
        for v in obj:
            _walk(v, out)
    elif isinstance(obj, str):
        out.extend(float(x) for x in NUM_RE.findall(obj))


def _walk_group(items: Sequence[Any], out: List[float]) -> None:
    """Add 100*k/n for sibling numbers (a count out of a total stored next to it)."""
    nums = [float(v) for v in items if isinstance(v, (int, float)) and not isinstance(v, bool)]
    if 2 <= len(nums) <= 12:
        for k in nums:
            for n in nums:
                if n > 0 and 0 <= k <= n and k == int(k):
                    out.append(k / n)


def source_values(paths: Iterable[Path]) -> Tuple[np.ndarray, np.ndarray]:
    """(raw values, raw values with their x100 and /100 variants), each sorted and unique."""
    vals: List[float] = []
    for p in paths:
        if p.suffix.lower() not in TEXT_SUFFIXES:
            continue
        try:
            text = p.read_text(errors="ignore")
        except OSError:
            continue
        if p.suffix.lower() == ".json":
            try:
                _walk(json.loads(text), vals)
                continue
            except ValueError:
                pass
        vals.extend(float(x) for x in NUM_RE.findall(text))
        for m in PAIR_RE.finditer(text):
            if int(m.group(2)) > 0:
                vals.append(int(m.group(1)) / int(m.group(2)))
    arr = np.asarray(vals, dtype=np.float64)
    arr = arr[np.isfinite(arr)]
    return np.unique(arr), np.unique(np.concatenate([arr, arr * 100.0, arr / 100.0]))


def _hit(sorted_vals: np.ndarray, value: float, decimals: int) -> bool:
    tol = 0.5 * 10.0 ** (-decimals) + 1e-9
    lo = np.searchsorted(sorted_vals, value - tol, side="left")
    hi = np.searchsorted(sorted_vals, value + tol, side="right")
    return bool(hi > lo)


Values = Tuple[np.ndarray, np.ndarray]


def _ok(d: Dict[str, Any], vals: Values) -> bool:
    raw, scaled = vals
    if d["kind"] == "pair":
        if _hit(raw, d["a"], 0) and _hit(raw, d["b"], 0):
            return True
        return _hit(scaled, d["a"] / d["b"], 3) or _hit(scaled, 100.0 * d["a"] / d["b"], 1)
    plain_int = d["decimals"] == 0 and not d["pct"]
    return _hit(raw if plain_int else scaled, d["value"], d["decimals"])


def audit(doc_numbers: Sequence[Dict[str, Any]], vals: Values) -> List[Dict[str, Any]]:
    """The doc numbers with no matching source value. Plain integers need an exact source
    value; decimals and percentages may also match a source value x100 or /100."""
    return [d for d in doc_numbers if not _ok(d, vals)]


def chance_rate(doc_numbers: Sequence[Dict[str, Any]], vals: Values) -> float:
    """Share of doc numbers that still match after being shifted by 3 in the last printed digit:
    how often the check passes by chance (low-precision decimals match almost anything)."""
    shifted = []
    for d in doc_numbers:
        if d["kind"] == "num":
            e = dict(d)
            e["value"] = d["value"] + 3 * 10.0 ** (-d["decimals"])
            shifted.append(e)
    if not shifted:
        return float("nan")
    return 1.0 - len(audit(shifted, vals)) / len(shifted)


def expand_sources(patterns: Sequence[str]) -> List[Path]:
    paths: List[Path] = []
    for pat in patterns:
        hits = [Path(h) for h in sorted(glob.glob(pat))]
        paths.extend(h for h in hits if h.is_file())
    return paths


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--doc", required=True, type=Path)
    ap.add_argument("--section", required=True)
    ap.add_argument("--sources", required=True, nargs="+")
    ap.add_argument("--strict", action="store_true")
    args = ap.parse_args(argv)

    lines = args.doc.read_text().splitlines()
    start, end = find_section(lines, args.section)
    nums = extract_numbers(lines[start:end], offset=start)
    paths = expand_sources(args.sources)
    if not paths:
        print("error: no source files found", file=sys.stderr)
        return 2
    vals = source_values(paths)
    unmatched = audit(nums, vals)
    total = len(nums)
    print(f"{total - len(unmatched)}/{total} numbers matched ({len(paths)} source files)")
    print(
        f"chance-match rate (numbers shifted by 3 in the last digit): {chance_rate(nums, vals):.0%}"
    )
    for d in unmatched:
        print(f"{args.doc}:{d['line']}: unmatched {d['text']!r}")
    return 1 if (unmatched and args.strict) else 0


if __name__ == "__main__":
    sys.exit(main())
