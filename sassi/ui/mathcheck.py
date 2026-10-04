"""Check the LaTeX formulas of the lessons and Help pages with the KaTeX that the GUI uses.

The GUI typesets the formulas of a Markdown text (``$...$``, ``$$...$$``, fenced ``math`` blocks; see
:mod:`sassi.ui.markdown`) in the browser with KaTeX (``sassi/ui/static/katex/katex.min.js``).  This module
runs the same ``katex.min.js`` in Node.js (``katex.renderToString`` with ``throwOnError``) so a formula that
would show as a red error in the GUI fails the tests instead (``tests/unit/test_lessons.py``) -- when
Node.js is not installed the check is skipped.

Command line (for lesson authors)::

    .venv/bin/python -m sassi.ui.mathcheck sassi/ui/lessons/*.md docs/theory/THEORY_MANUAL.md
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence

from .markdown import math_spans

KATEX_JS = Path(__file__).resolve().parent / "static" / "katex" / "katex.min.js"

#: the KaTeX options of the GUI (static/learn.js renderMath): keep the two in step
_NODE_SCRIPT = r"""
const katex = require(process.argv[1]);
let data = "";
process.stdin.on("data", (c) => { data += c; });
process.stdin.on("end", () => {
  const out = JSON.parse(data).map((f) => {
    try {
      katex.renderToString(f.tex, {displayMode: !!f.display, throwOnError: true, strict: "ignore", trust: false});
      return null;
    } catch (e) { return String(e.message || e); }
  });
  process.stdout.write(JSON.stringify(out));
});
"""


def node_path() -> Optional[str]:
    """The Node.js executable, or None."""
    return shutil.which("node")


def check_formulas(formulas: Sequence[Dict]) -> List[Optional[str]]:
    """KaTeX error message (or None) of every ``{"tex", "display"}``; raises RuntimeError without Node.js."""
    node = node_path()
    if node is None:
        raise RuntimeError("Node.js is not installed")
    if not formulas:
        return []
    r = subprocess.run([node, "-e", _NODE_SCRIPT, str(KATEX_JS)], input=json.dumps(list(formulas)),
                       capture_output=True, text=True, timeout=120)
    if r.returncode != 0:
        raise RuntimeError(f"node failed: {r.stderr.strip()[:500]}")
    return json.loads(r.stdout)


def check_text(text: str, where: str = "") -> List[str]:
    """``"<where>:<line>: <error> in <tex>"`` for every formula of Markdown ``text`` that KaTeX rejects."""
    spans = math_spans(text)
    errs = check_formulas(spans)
    return [f"{where}:{s['line']}: {e} -- in: {s['tex'][:120]}" for s, e in zip(spans, errs) if e]


def main(argv: Optional[Sequence[str]] = None) -> int:
    paths = list(argv if argv is not None else sys.argv[1:])
    if not paths:
        print(__doc__)
        return 2
    bad = n = 0
    for p in paths:
        text = Path(p).read_text(encoding="utf-8")
        n += len(math_spans(text))
        for e in check_text(text, p):
            print(e)
            bad += 1
    print(f"{n} formulas, {bad} errors")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
