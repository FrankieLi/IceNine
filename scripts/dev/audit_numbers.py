"""Mechanical half of the claims audit: do the numbers in a doc section occur in the sources?

Usage: audit_numbers.py --doc FILE.md --section "Heading text" --sources FILE_OR_GLOB... [--strict] [-v] [--per-file]

Numbers (decimals, percentages, ranges, "a / b" pairs, integers >= 10) are extracted from the
section; a number matches when a source value equals it at the printed precision. Only for a
percentage or a fraction <= 1 are the x100 and /100 variants and sibling k/n ratios also tried.
Unmatched numbers are listed with line numbers; -v prints the nearest source value and its
file:key; --per-file requires each paragraph or table row to match within one source file. Exit 0
unless --strict. A match is not verification: it only shows that the number occurs somewhere.
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
                    "pct": suffix == "%"
                    or range_pct
                    or (decimals >= 1 and value <= 100 and "%" in raw),
                }
            )
    return found


class Store:
    """Source numbers: absolute values sorted, with the file index and key of each."""

    def __init__(self) -> None:
        self.files: List[str] = []
        self._v: List[float] = []
        self._f: List[int] = []
        self._k: List[str] = []
        self._rv: List[float] = []  # derived k/n ratios of sibling numbers
        self._rf: List[int] = []
        self._rk: List[str] = []
        self.sealed = False

    def add(self, value: float, fid: int, key: str) -> None:
        if np.isfinite(value):
            self._v.append(abs(value))
            self._f.append(fid)
            self._k.append(key)

    def add_ratio(self, value: float, fid: int, key: str) -> None:
        self._rv.append(value)
        self._rf.append(fid)
        self._rk.append(key)

    def seal(self) -> "Store":
        order = np.argsort(self._v, kind="stable")
        self.v = np.asarray(self._v, dtype=np.float64)[order]
        self.f = np.asarray(self._f, dtype=np.int64)[order]
        self.k = [self._k[i] for i in order]
        order = np.argsort(self._rv, kind="stable")
        self.rv = np.asarray(self._rv, dtype=np.float64)[order]
        self.rf = np.asarray(self._rf, dtype=np.int64)[order]
        self.rk = [self._rk[i] for i in order]
        self.sealed = True
        return self


def _walk(obj: Any, st: Store, fid: int, key: str) -> None:
    if isinstance(obj, bool):
        return
    if isinstance(obj, (int, float)):
        st.add(float(obj), fid, key)
    elif isinstance(obj, dict):
        _walk_group(list(obj.values()), st, fid, key)
        for k, v in obj.items():
            _walk(v, st, fid, f"{key}.{k}" if key else str(k))
    elif isinstance(obj, list):
        if len(obj) > MAX_ARRAY and all(isinstance(v, (int, float)) for v in obj):
            return  # raw data array, not a reported number
        _walk_group(obj, st, fid, key)
        for i, v in enumerate(obj):
            _walk(v, st, fid, f"{key}[{i}]")
    elif isinstance(obj, str):
        for x in NUM_RE.findall(obj):
            st.add(float(x), fid, key)


def _walk_group(items: Sequence[Any], st: Store, fid: int, key: str) -> None:
    """Add k/n for sibling numbers (a count out of a total stored next to it)."""
    nums = [float(v) for v in items if isinstance(v, (int, float)) and not isinstance(v, bool)]
    if 2 <= len(nums) <= 12:
        for k in nums:
            for n in nums:
                if n > 0 and 0 <= k <= n and k == int(k):
                    st.add_ratio(k / n, fid, f"{key} (k/n of siblings)")


def source_values(paths: Iterable[Path]) -> Store:
    st = Store()
    for p in paths:
        if p.suffix.lower() not in TEXT_SUFFIXES:
            continue
        try:
            text = p.read_text(errors="ignore")
        except OSError:
            continue
        fid = len(st.files)
        st.files.append(str(p))
        if p.suffix.lower() == ".json":
            try:
                _walk(json.loads(text), st, fid, "")
                continue
            except ValueError:
                pass
        for ln, line in enumerate(text.splitlines(), 1):
            for x in NUM_RE.findall(line):
                st.add(float(x), fid, f"line {ln}")
            for m in PAIR_RE.finditer(line):
                if int(m.group(2)) > 0:
                    st.add_ratio(int(m.group(1)) / int(m.group(2)), fid, f"line {ln} (a/b)")
    return st.seal()


Hit = Tuple[float, int, str]  # source value, file index, key


def _search(vals: np.ndarray, value: float, tol: float) -> List[int]:
    lo = int(np.searchsorted(vals, value - tol, side="left"))
    hi = int(np.searchsorted(vals, value + tol, side="right"))
    return list(range(lo, hi))


def _candidates(value: float, decimals: int, pct: bool) -> List[Tuple[float, float, bool]]:
    """(target, tolerance, ratio-array?) pairs to look for in the sources.

    Scaled variants (x100, /100) and the sibling-ratio expansion apply only when the doc number
    is a percentage or a fraction <= 1; other numbers must match a source value exactly.
    """
    tol = 0.5 * 10.0 ** (-decimals) + 1e-9
    out = [(value, tol, False)]
    if pct:
        out += [(value / 100.0, tol / 100.0, False), (value / 100.0, tol / 100.0, True)]
        out += [(value, tol, True)]
    elif value <= 1.0:
        out += [(value * 100.0, tol * 100.0, False), (value, tol, True)]
    return out


def find_matches(d: Dict[str, Any], st: Store) -> List[Hit]:
    """All source hits for a doc number (empty if unmatched)."""
    if d["kind"] == "pair":
        a, b = d["a"], d["b"]
        ha = [(st.v[i], int(st.f[i]), st.k[i]) for i in _search(st.v, a, 0.5)]
        hb = [(st.v[i], int(st.f[i]), st.k[i]) for i in _search(st.v, b, 0.5)]
        files_b = {h[1] for h in hb}
        hits = [h for h in ha if h[1] in files_b]
        if hits:
            return hits
        cands = _candidates(a / b, 3, True) if a < b else []
    else:
        cands = _candidates(d["value"], d["decimals"], d["pct"])
    hits = []
    for target, tol, ratio in cands:
        vals, fs, ks = (st.rv, st.rf, st.rk) if ratio else (st.v, st.f, st.k)
        hits += [(vals[i], int(fs[i]), ks[i]) for i in _search(vals, target, tol)]
    return hits


def audit(doc_numbers: Sequence[Dict[str, Any]], st: Store) -> List[Dict[str, Any]]:
    """The doc numbers with no matching source value."""
    return [d for d in doc_numbers if not find_matches(d, st)]


def _bucket(d: Dict[str, Any]) -> str:
    return (
        "0 dp"
        if d["decimals"] <= 0
        else f"{min(d['decimals'], 3)}{'+' if d['decimals'] >= 3 else ''} dp"
    )


def chance_rates(doc_numbers: Sequence[Dict[str, Any]], st: Store) -> Dict[str, Tuple[float, int]]:
    """Per precision bucket: the share of doc numbers that still match after being shifted by 3
    in the last printed digit (how often the check passes by chance), and the bucket size."""
    groups: Dict[str, List[Dict[str, Any]]] = {}
    for d in doc_numbers:
        if d["kind"] == "num":
            e = dict(d)
            e["value"] = d["value"] + 3 * 10.0 ** (-d["decimals"])
            groups.setdefault(_bucket(d), []).append(e)
    return {k: (1.0 - len(audit(v, st)) / len(v), len(v)) for k, v in sorted(groups.items())}


def chance_rate(doc_numbers: Sequence[Dict[str, Any]], st: Store) -> float:
    """Overall chance-match rate (all precisions)."""
    shifted = []
    for d in doc_numbers:
        if d["kind"] == "num":
            e = dict(d)
            e["value"] = d["value"] + 3 * 10.0 ** (-d["decimals"])
            shifted.append(e)
    if not shifted:
        return float("nan")
    return 1.0 - len(audit(shifted, st)) / len(shifted)


def group_ids(lines: Sequence[str]) -> List[int]:
    """Paragraph / table-row group id of each line (blank lines end a paragraph)."""
    gids, gid, in_para = [], 0, False
    for line in lines:
        if not line.strip():
            in_para = False
            gids.append(-1)
        elif line.lstrip().startswith("|"):
            gid += 1
            in_para = False
            gids.append(gid)
        else:
            if not in_para:
                gid += 1
                in_para = True
            gids.append(gid)
    return gids


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
    ap.add_argument("-v", "--verbose", action="store_true", help="print the nearest source hit")
    ap.add_argument(
        "--per-file",
        action="store_true",
        help="numbers in one paragraph or table row must all match within a single source file",
    )
    args = ap.parse_args(argv)

    lines = args.doc.read_text().splitlines()
    start, end = find_section(lines, args.section)
    nums = extract_numbers(lines[start:end], offset=start)
    paths = expand_sources(args.sources)
    if not paths:
        print("error: no source files found", file=sys.stderr)
        return 2
    st = source_values(paths)
    hits = {id(d): find_matches(d, st) for d in nums}
    unmatched = [d for d in nums if not hits[id(d)]]
    total = len(nums)
    print(f"{total - len(unmatched)}/{total} numbers matched ({len(paths)} source files)")
    rates = chance_rates(nums, st)
    print(
        f"chance-match rate: {chance_rate(nums, st):.0%} overall; "
        + "; ".join(f"{k}: {r:.0%} (n={n})" for k, (r, n) in rates.items())
    )
    bad = len(unmatched)
    for d in unmatched:
        print(f"{args.doc}:{d['line']}: unmatched {d['text']!r}")
    if args.verbose:
        for d in nums:
            h = hits[id(d)]
            if h:
                best = min(h, key=lambda x: abs(x[0] - d.get("value", d.get("a", 0.0))))
                print(
                    f"{args.doc}:{d['line']}: {d['text']!r} <- {best[0]:g} at "
                    f"{Path(st.files[best[1]]).name}:{best[2]}"
                )
    if args.per_file:
        gids = group_ids(lines[start:end])
        by_group: Dict[int, List[Dict[str, Any]]] = {}
        for d in nums:
            if hits[id(d)]:
                by_group.setdefault(gids[d["line"] - start - 1], []).append(d)
        for g, ds in by_group.items():
            common = set.intersection(*({h[1] for h in hits[id(d)]} for d in ds))
            if not common:
                texts = ", ".join(repr(d["text"]) for d in ds)
                print(f"{args.doc}:{ds[0]['line']}: no single source file matches all of {texts}")
                bad += 1
    return 1 if (bad and args.strict) else 0


if __name__ == "__main__":
    sys.exit(main())
