"""Guided-course lessons: parser unit tests and validation of every lesson in sassi/ui/lessons/
(static checks, and a headless run of every step: no command may produce an error message, and every
file an action names must exist after its step).  Format: docs/internal/lesson_format.md."""
from __future__ import annotations

import pytest

from sassi.ui import lessons as L

SAMPLE = """---
id: 99-sample
title: Sample lesson
part: Fundamentals
order: 99
minutes: 3
summary: a tiny lesson used by the parser tests
objectives: [one, two]
prerequisites: []
---
Intro text.

```sassi-setup
TIT,sample
```

## Define two nodes
Narrative of step 1.

```sassi
* two nodes
N,1,0,0,0
N,2,1,0,0
```

### What this does
`N,<nd>,<x>,<y>,<z>` defines a node.

```sassi-show
N,<nd>,<x>,<y>,<z>
```

### Check yourself
How many nodes?
Answer: two.

```action
plot-nodes
```

## Group
```sassi
GROUP,1,2
```

### Something else
Plain section with a ## inside a fence:

```text
## not a heading
```
"""


def test_parse_sample():
    les = L.parse_lesson(SAMPLE, "sample.md")
    assert les.id == "99-sample" and les.objectives == ["one", "two"] and les.prerequisites == []
    assert les.setup == ["TIT,sample"]
    assert [s.title for s in les.steps] == ["Define two nodes", "Group"]
    s1, s2 = les.steps
    assert s1.commands == ["* two nodes", "N,1,0,0,0", "N,2,1,0,0"]
    assert [x.key for x in s1.sections] == ["what", "check"]
    assert s1.actions == [("plot-nodes", "")]
    assert L.split_answer(s1.sections[1].markdown) == ("How many nodes?", "two.")
    assert s2.commands == ["GROUP,1,2"] and s2.sections[0].key == "other"
    assert "## not a heading" in s2.sections[0].markdown
    d = les.to_dict()
    assert d["nsteps"] == 2 and d["steps"][0]["actions"] == [{"verb": "plot-nodes", "args": ""}]


def test_front_matter_errors():
    with pytest.raises(L.LessonError):
        L.parse_lesson("no front matter\n## step\n")
    with pytest.raises(L.LessonError):
        L.parse_lesson("---\ntitle: x\n---\n## s\n")


def test_headless_run_of_sample(tmp_path):
    les = L.parse_lesson(SAMPLE, "sample.md")
    rep = L.run_lesson_headless(les, tmp_path)
    assert rep.ok, rep
    les.steps[1].blocks[0].text = "NOSUCHCOMMAND,1"
    rep = L.run_lesson_headless(les, tmp_path)
    assert not rep.ok and rep.errors


# ---------------------------------------------------------------------- the real course
LESSONS = L.load_lessons()


def test_course_ids_unique_and_prerequisites_exist():
    ids = [x.id for x in LESSONS]
    assert len(ids) == len(set(ids))
    for les in LESSONS:
        for p in les.prerequisites:
            assert p in ids, f"{les.id}: prerequisite {p} not found"


@pytest.mark.parametrize("lesson", LESSONS, ids=[x.id for x in LESSONS])
def test_lesson_static_checks(lesson):
    assert L.validate_lesson(lesson) == []


@pytest.mark.slow
@pytest.mark.parametrize("lesson", LESSONS, ids=[x.id for x in LESSONS])
def test_lesson_runs_without_errors(lesson, tmp_path):
    rep = L.run_lesson_headless(lesson, tmp_path)
    msg = "\n".join(f"{w}: {c!r} -> {e}" for w, c, e in rep.errors[:20])
    msg += "\n" + "\n".join(f"{w}: missing {p}" for w, p in rep.missing[:20])
    assert rep.ok, msg


# ---------------------------------------------------------------------- LaTeX formulas
def test_math_rendering_rules():
    from sassi.ui.markdown import math_spans, render
    t = "Text $\\omega^2 M_s$, a \\$5 fee, `$code$` and US$5.\n\n```math\na = b\n```\n\n$$ x_i $$\n"
    spans = math_spans(t)
    assert [(s["tex"], s["display"]) for s in spans] == [("\\omega^2 M_s", False), ("a = b", True), ("x_i", True)]
    html = render(t).html
    assert "$5 fee" in html and "<code>$code$</code>" in html and "US$5" in html
    assert '<div class="math-display" data-tex="a = b">' in html


@pytest.mark.parametrize("lesson", LESSONS, ids=[x.id for x in LESSONS])
def test_lesson_formulas_typeset_with_katex(lesson):
    from sassi.ui import mathcheck
    if mathcheck.node_path() is None:
        pytest.skip("Node.js is not installed (the GUI's KaTeX check runs in Node)")
    from pathlib import Path
    errs = mathcheck.check_text(Path(lesson.path).read_text(encoding="utf-8"), lesson.id)
    assert errs == [], "\n".join(errs)
