"""Guided course: lesson files (parser, workspace preparation, headless runner, validator).

The lesson format is specified in ``docs/internal/lesson_format.md``.  Lessons are Markdown files in
``sassi/ui/lessons/`` with a small front-matter header; every ``##`` heading is a step, ``###``
headings are the step's sections ("What this does", "Why it matters", "Technical basis", ...),
```` ```sassi ```` blocks hold the commands a step runs, ```` ```action ```` blocks the GUI actions
offered after it.

The GUI (``sassi.ui.api``) serves the parsed lessons and runs their commands through the session
interpreter; ``tests/unit/test_lessons.py`` validates every lesson by running it headlessly with
:func:`run_lesson_headless`.
"""
from __future__ import annotations

import re
import shutil
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple, Union

PathLike = Union[str, Path]

LESSON_DIR = Path(__file__).resolve().parent / "lessons"
PROJECT_ROOT = Path(__file__).resolve().parents[2]
EXAMPLES_DIR = PROJECT_ROOT / "examples"

#: recognised step-section titles (lower case) -> key used by the GUI
SECTION_KEYS = {
    "what this does": "what",
    "why it matters": "why",
    "technical basis": "basis",
    "in ansys terms": "ansys",
    "try this": "try",
    "check yourself": "check",
}
#: GUI action verbs -> True when the argument list names files that must exist after the step
ACTION_VERBS = {
    "plot-model": False, "plot-nodes": False, "plot-layers": False, "plot-soilprops": False,
    "plot-spectrum": True, "plot-history": True, "open-file": True, "open-listing": False,
    "open-dialog": False, "open-doc": False, "explain": False,
}
RUN_BLOCKS = ("sassi",)
SHOW_BLOCKS = ("sassi-show",)
SETUP_BLOCKS = ("sassi-setup",)
FENCE_RE = re.compile(r"^(`{3,})\s*([\w-]*)\s*(.*)$")
ID_RE = re.compile(r"^[0-9a-z][0-9a-z-]*$")


@dataclass
class Block:
    """A fenced block (``kind`` = 'sassi', 'sassi-show', 'sassi-setup', 'action' or 'code')."""
    kind: str
    lang: str
    text: str

    @property
    def lines(self) -> List[str]:
        return [ln.rstrip() for ln in self.text.split("\n") if ln.strip()]


@dataclass
class Section:
    key: str          # 'what', 'why', 'basis', 'ansys', 'try', 'check' or 'other'
    title: str
    markdown: str


@dataclass
class Step:
    index: int                      # 1-based
    title: str
    narrative: str                  # markdown before the first section
    sections: List[Section] = field(default_factory=list)
    blocks: List[Block] = field(default_factory=list)   # every fenced block of the step, in order

    @property
    def commands(self) -> List[str]:
        """Command lines the step runs (the ``sassi`` blocks, in order)."""
        out: List[str] = []
        for b in self.blocks:
            if b.kind == "sassi":
                out.extend(b.lines)
        return out

    @property
    def actions(self) -> List[Tuple[str, str]]:
        """``(verb, arguments)`` of the step's ``action`` blocks."""
        out: List[Tuple[str, str]] = []
        for b in self.blocks:
            if b.kind == "action":
                for ln in b.lines:
                    verb, _, args = ln.partition(":")
                    out.append((verb.strip().lower(), args.strip()))
        return out


@dataclass
class Lesson:
    id: str
    title: str
    part: str = ""
    order: int = 0
    minutes: int = 0
    example: str = ""
    summary: str = ""
    objectives: List[str] = field(default_factory=list)
    prerequisites: List[str] = field(default_factory=list)
    intro: str = ""
    setup: List[str] = field(default_factory=list)
    steps: List[Step] = field(default_factory=list)
    path: str = ""

    def to_dict(self, full: bool = True) -> Dict[str, Any]:
        d = {k: getattr(self, k) for k in ("id", "title", "part", "order", "minutes", "example", "summary",
                                            "objectives", "prerequisites")}
        d["nsteps"] = len(self.steps)
        d["step_titles"] = [s.title for s in self.steps]
        if full:
            d["intro"] = self.intro
            d["setup"] = list(self.setup)
            d["steps"] = [{"index": s.index, "title": s.title, "narrative": s.narrative,
                           "sections": [asdict(x) for x in s.sections],
                           "blocks": [{"kind": b.kind, "lang": b.lang, "text": b.text} for b in s.blocks],
                           "commands": s.commands,
                           "actions": [{"verb": v, "args": a} for v, a in s.actions]} for s in self.steps]
        return d


class LessonError(ValueError):
    """A lesson file that cannot be parsed."""


# --------------------------------------------------------------------------------------
# parsing
# --------------------------------------------------------------------------------------
def _parse_front_matter(lines: List[str], where: str) -> Tuple[Dict[str, Any], int]:
    if not lines or lines[0].strip() != "---":
        raise LessonError(f"{where}: missing front matter (first line must be ---)")
    meta: Dict[str, Any] = {}
    for i in range(1, len(lines)):
        ln = lines[i]
        if ln.strip() == "---":
            return meta, i + 1
        if not ln.strip() or ln.lstrip().startswith("#"):
            continue
        key, sep, val = ln.partition(":")
        if not sep:
            raise LessonError(f"{where}: front matter line {i + 1} is not 'key: value'")
        val = val.strip()
        if val.startswith("[") and val.endswith("]"):
            meta[key.strip()] = [v.strip() for v in val[1:-1].split(",") if v.strip()]
        else:
            meta[key.strip()] = val
    raise LessonError(f"{where}: front matter not closed with ---")


def _split_blocks(text: str) -> Tuple[str, List[Block]]:
    """Return (text with fences kept, blocks in order) for one chunk of markdown."""
    blocks: List[Block] = []
    lines = text.split("\n")
    i = 0
    while i < len(lines):
        m = FENCE_RE.match(lines[i].strip())
        if m and lines[i].lstrip().startswith("```"):
            fence, lang = m.group(1), m.group(2).lower()
            body: List[str] = []
            i += 1
            while i < len(lines) and not lines[i].strip().startswith(fence):
                body.append(lines[i])
                i += 1
            kind = lang if lang in RUN_BLOCKS + SHOW_BLOCKS + SETUP_BLOCKS + ("action",) else "code"
            blocks.append(Block(kind, lang, "\n".join(body)))
        i += 1
    return text, blocks


def _split_on_heading(lines: List[str], level: int) -> List[Tuple[Optional[str], List[str]]]:
    """Split lines into (heading title or None for the preamble, body lines), ignoring fenced code."""
    out: List[Tuple[Optional[str], List[str]]] = [(None, [])]
    prefix = "#" * level + " "
    in_fence = None
    for ln in lines:
        s = ln.strip()
        if s.startswith("```"):
            tick = re.match(r"^`{3,}", s).group(0)
            if in_fence is None:
                in_fence = tick
            elif s.startswith(in_fence) and s.strip("`").strip() == "":
                in_fence = None
            out[-1][1].append(ln)
            continue
        if in_fence is None and ln.startswith(prefix) and not ln.startswith(prefix + "#"):
            out.append((ln[len(prefix):].strip(), []))
        else:
            out[-1][1].append(ln)
    return out


def parse_lesson(text: str, path: PathLike = "") -> Lesson:
    """Parse one lesson file (see docs/internal/lesson_format.md)."""
    where = str(path) or "<lesson>"
    lines = text.replace("\r\n", "\n").split("\n")
    meta, start = _parse_front_matter(lines, where)
    for req in ("id", "title"):
        if not meta.get(req):
            raise LessonError(f"{where}: front matter needs '{req}'")
    try:
        order = int(meta.get("order", 0) or 0)
        minutes = int(meta.get("minutes", 0) or 0)
    except ValueError as exc:
        raise LessonError(f"{where}: order/minutes must be integers ({exc})")
    lesson = Lesson(id=str(meta["id"]), title=str(meta["title"]), part=str(meta.get("part", "")), order=order,
                    minutes=minutes, example=str(meta.get("example", "") or ""),
                    summary=str(meta.get("summary", "")),
                    objectives=list(meta.get("objectives", []) or []),
                    prerequisites=list(meta.get("prerequisites", []) or []), path=str(path))
    chunks = _split_on_heading(lines[start:], 2)
    intro_text = "\n".join(chunks[0][1]).strip("\n")
    _, intro_blocks = _split_blocks(intro_text)
    lesson.intro = intro_text
    lesson.setup = [ln for b in intro_blocks if b.kind == "sassi-setup" for ln in b.lines]
    for k, (title, body) in enumerate(chunks[1:], start=1):
        sec_chunks = _split_on_heading(body, 3)
        narrative = "\n".join(sec_chunks[0][1]).strip("\n")
        step = Step(index=k, title=title or f"Step {k}", narrative=narrative)
        _, nb = _split_blocks(narrative)
        step.blocks.extend(nb)
        for stitle, sbody in sec_chunks[1:]:
            md = "\n".join(sbody).strip("\n")
            key = SECTION_KEYS.get((stitle or "").strip().lower(), "other")
            step.sections.append(Section(key, stitle or "", md))
            _, sb = _split_blocks(md)
            step.blocks.extend(sb)
        lesson.steps.append(step)
    return lesson


def load_lessons(directory: PathLike = LESSON_DIR) -> List[Lesson]:
    """All lessons of ``directory`` sorted by (part order, order, id)."""
    d = Path(directory)
    if not d.is_dir():
        return []
    lessons = [parse_lesson(p.read_text(encoding="utf-8"), p) for p in sorted(d.glob("*.md"))]
    parts = {"fundamentals": 0, "design applications": 1, "advanced": 2}
    return sorted(lessons, key=lambda l: (parts.get(l.part.lower(), 9), l.order, l.id))


def lesson_by_id(lesson_id: str, directory: PathLike = LESSON_DIR) -> Lesson:
    for lesson in load_lessons(directory):
        if lesson.id == lesson_id:
            return lesson
    raise KeyError(f"no lesson '{lesson_id}'")


def strip_blocks(markdown: str, kinds: Iterable[str] = ("action", "sassi-setup")) -> str:
    """``markdown`` without the fenced blocks of the given kinds (rendered as buttons, not text)."""
    kinds = set(kinds)
    out: List[str] = []
    lines = markdown.split("\n")
    i = 0
    while i < len(lines):
        m = FENCE_RE.match(lines[i].strip())
        if m and lines[i].lstrip().startswith("```") and m.group(2).lower() in kinds:
            fence = m.group(1)
            i += 1
            while i < len(lines) and not lines[i].strip().startswith(fence):
                i += 1
            i += 1
            continue
        out.append(lines[i])
        i += 1
    return "\n".join(out).strip("\n")


def split_answer(markdown: str) -> Tuple[str, str]:
    """``(question, answer)`` of a 'Check yourself' section (answer starts at a line 'Answer:');
    action and setup blocks are removed first."""
    lines = strip_blocks(markdown).split("\n")
    for i, ln in enumerate(lines):
        if ln.strip().lower().startswith("answer:"):
            first = ln.strip()[len("answer:"):].strip()
            return "\n".join(lines[:i]).strip(), "\n".join([first] + lines[i + 1:]).strip()
    return markdown.strip(), ""


# --------------------------------------------------------------------------------------
# workspace and headless run
# --------------------------------------------------------------------------------------
#: repository files an example reads through ``../sassi/...`` (e.g. ``INP,../sassi/data/dynp_library.pre``
#: of examples 4 and 6): the path is relative to the model directory, one level below the examples
#: directory, so in a copied workspace the file must exist at ``<workspace>/sassi/...``
SUPPORT_RE = re.compile(r"(?<![\w.])\.\./(sassi/data/[\w./-]+)")


def example_support_files(pre_text: str) -> List[str]:
    """Repository-relative paths (``sassi/data/...``) of the support files a ``.pre`` file reads."""
    out: List[str] = []
    for m in SUPPORT_RE.finditer(pre_text):
        rel = m.group(1).rstrip(".")
        parts = rel.split("/")
        if ".." in parts or not rel.startswith("sassi/data/"):
            continue
        if rel not in out and (PROJECT_ROOT / rel).is_file():
            out.append(rel)
    return out


def copy_example(name: str, ws: PathLike, examples_dir: PathLike = EXAMPLES_DIR) -> Path:
    """Copy ``examples/<name>.pre``, the ``examples/data/`` folder and the repository support files the
    example reads (:func:`example_support_files`) into the directory ``ws``; returns the copied ``.pre``."""
    ex, ws = Path(examples_dir), Path(ws)
    src = ex / f"{name}.pre"
    if not src.is_file():
        raise FileNotFoundError(f"example {src} not found")
    ws.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, ws / src.name)
    if (ex / "data").is_dir():
        shutil.copytree(ex / "data", ws / "data", dirs_exist_ok=True)
    for rel in example_support_files(src.read_text(encoding="utf-8", errors="replace")):
        dst = ws / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(PROJECT_ROOT / rel, dst)
    return ws / src.name


def prepare_workspace(lesson: Lesson, root: PathLike, examples_dir: PathLike = EXAMPLES_DIR,
                      fresh: bool = True) -> Path:
    """Create ``<root>/<lesson id>/`` (emptied when ``fresh``) and copy the lesson's example files
    (:func:`copy_example`: the ``.pre``, ``data/`` and the support files it reads)."""
    ws = Path(root) / lesson.id
    if fresh and ws.exists():
        shutil.rmtree(ws)
    ws.mkdir(parents=True, exist_ok=True)
    if lesson.example:
        try:
            copy_example(lesson.example, ws, examples_dir)
        except FileNotFoundError:
            raise FileNotFoundError(f"lesson {lesson.id}: example {Path(examples_dir) / (lesson.example + '.pre')} "
                                    f"not found") from None
    return ws


@dataclass
class RunReport:
    lesson: str
    workspace: str
    errors: List[Tuple[str, str, str]] = field(default_factory=list)     # (where, command, message)
    missing: List[Tuple[str, str]] = field(default_factory=list)         # (where, path)

    @property
    def ok(self) -> bool:
        return not self.errors and not self.missing


def _action_files(args: str) -> List[str]:
    files = args.split("|")[0]
    return [f.strip() for f in files.split(",") if f.strip()]


def run_lesson_headless(lesson: Lesson, root: PathLike, examples_dir: PathLike = EXAMPLES_DIR) -> RunReport:
    """Run the lesson like a learner pressing every 'Run step' in order (setup first); collect errors
    and action files that do not exist after their step."""
    from ..prep import Interpreter, Kind
    ws = prepare_workspace(lesson, root, examples_dir)
    rep = RunReport(lesson.id, str(ws))
    ui = Interpreter(cwd=ws)

    def run(where: str, lines: Iterable[str]) -> None:
        for ln in lines:
            before = len(ui.sink.texts(Kind.ERROR))
            ui.execute(ln)
            for e in ui.sink.texts(Kind.ERROR)[before:]:
                rep.errors.append((where, ln, e))

    run("setup", lesson.setup)
    for step in lesson.steps:
        where = f"step {step.index} ({step.title})"
        run(where, step.commands)
        for verb, args in step.actions:
            if verb == "open-listing":
                m = ui.model
                p = Path(m.path or ws) / f"{m.name}_{args.strip().upper()}.out"
                if not p.is_file():
                    rep.missing.append((where, str(p)))
            elif ACTION_VERBS.get(verb):
                for f in _action_files(args):
                    if not (ws / f).is_file():
                        rep.missing.append((where, f))
    return rep


def validate_lesson(lesson: Lesson) -> List[str]:
    """Static checks: front matter, section titles, action verbs, command names, doc links."""
    from ..prep.registry import lookup
    problems: List[str] = []
    if not ID_RE.match(lesson.id):
        problems.append(f"id '{lesson.id}' must match {ID_RE.pattern}")
    if not lesson.steps:
        problems.append("lesson has no steps (## headings)")
    if lesson.example and not (EXAMPLES_DIR / f"{lesson.example}.pre").is_file():
        problems.append(f"example '{lesson.example}' not found in examples/")
    for step in lesson.steps:
        where = f"step {step.index} ({step.title})"
        for verb, _ in step.actions:
            if verb not in ACTION_VERBS:
                problems.append(f"{where}: unknown action verb '{verb}'")
        for ln in step.commands:
            s = ln.strip()
            if s.startswith("*"):
                continue
            name = re.split(r"[,\s]", s, maxsplit=1)[0]
            if lookup(name) is None:
                problems.append(f"{where}: unknown command '{name}'")
        md = step.narrative + "\n" + "\n".join(x.markdown for x in step.sections)
        for target in re.findall(r"\]\(([^)\s]+)\)", md):
            if target.startswith(("http://", "https://", "#", "mailto:")):
                continue
            path = target.split("#", 1)[0]
            if path and not (PROJECT_ROOT / path).exists():
                problems.append(f"{where}: link target '{path}' does not exist (paths are relative to the repository root)")
    return problems
