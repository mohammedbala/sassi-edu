"""Learn: the server side of the guided course (Learn tab, lesson player, examples gallery).

* **Course outline** -- the lessons of ``sassi/ui/lessons/`` (:mod:`sassi.ui.lessons`, format
  ``docs/internal/lesson_format.md``) grouped by part, for the start page and the Learn menu.
* **Lesson rendering** -- every Markdown piece of a lesson (introduction, step narrative, the labelled
  sections) is rendered with the safe renderer of the Help pages (:mod:`sassi.ui.markdown`: raw HTML never
  passes).  Fenced ``sassi`` / ``sassi-show`` / ``sassi-setup`` blocks become placeholders
  (``<div class="lesson-block" data-block="k">``) that the browser fills with the interactive command
  block (explain affordance, run status); ``action`` blocks become buttons; the answer of "Check
  yourself" is split off (:func:`sassi.ui.lessons.split_answer`).  Links to ``docs/...`` (repository-root
  paths) open in the Help tab.
* **Workspaces** -- a lesson runs in ``<course root>/<lesson id>/``, an example in
  ``<course root>/examples/<name>/`` (course root: Learn > Course Workspace Folder, default
  ``<GUI start directory>/sassi-course``).  A workspace is recreated from scratch only when it carries
  the marker file this module writes (:data:`MARKER`), so a folder the course did not create is never
  deleted.
* **Commands** -- opening a lesson or loading an example returns the command text the browser submits
  through the session interpreter (rule L17): ``CD`` to the workspace, a fresh model (``ACTM,0`` and
  ``DMODEL`` of the models in memory), the lesson's ``sassi-setup`` commands, or ``INP`` of the example.
* **Actions** -- the verbs of the lesson ``action`` blocks resolved to command text (``plot-spectrum`` ->
  ``READSPEC`` + ``SPECPLOT`` ...) or to a GUI request (open a file, a dialog, a help page, the explainer).
* **Examples gallery** -- one entry per ``examples/*.pre``: the title and "What you learn" of its comment
  header and the topic / modules / run time row of ``examples/README.md``.
"""
from __future__ import annotations

import os
import posixpath
import re
import shutil
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from ..prep.lexer import join_command
from . import lessons as L
from .files import file_info
from .markdown import render

PARTS = ("Fundamentals", "Design applications", "Advanced")
DEFAULT_DIRNAME = "sassi-course"
EXAMPLES_SUBDIR = "examples"
#: marker file of a workspace created by the course (only such a folder is ever emptied)
MARKER = ".sassi-course-workspace"
EXAMPLE_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_\-]*$")
MODULE_RE = re.compile(r"^[A-Za-z][A-Za-z0-9]*$")
DIALOGS = ("ANALYSIS", "MODEL", "WRITE", "CHECK", "LOADGEN", "LOADGENDYN")
SECTION_TITLES = {"what": "What this does", "why": "Why it matters", "basis": "Technical basis",
                  "ansys": "In ANSYS terms", "try": "Try this", "check": "Check yourself"}
_PLACE = "XSASSIBLOCK{}X"
_PLACE_RE = re.compile(r"<p>XSASSIBLOCK(\d+)X</p>")


class LearnError(ValueError):
    """A request of the Learn features that cannot be served (``status``: HTTP status)."""

    def __init__(self, message: str, status: int = 400):
        super().__init__(message)
        self.status = status


# ======================================================================================
# workspaces
# ======================================================================================
def course_root(setting: str, start_dir: Path) -> Path:
    """The folder of the course workspaces: the setting (absolute, ``~`` expanded; a relative path is
    relative to the GUI start directory) or ``<start_dir>/sassi-course``."""
    s = (setting or "").strip()
    if not s:
        return (Path(start_dir) / DEFAULT_DIRNAME).resolve()
    p = Path(s).expanduser()
    if not p.is_absolute():
        p = Path(start_dir) / p
    return p.resolve()


def check_root(path: str, start_dir: Path) -> Path:
    """Validate a new course root (Learn > Course Workspace Folder): its parent folder must exist and it
    must not be an existing file; the filesystem root and the home folder itself are refused."""
    s = (path or "").strip()
    if not s:
        return course_root("", start_dir)
    if "\x00" in s:
        raise LearnError("invalid folder name")
    p = course_root(s, start_dir)
    if p.exists() and not p.is_dir():
        raise LearnError(f"{p} is a file, not a folder")
    if not p.parent.is_dir():
        raise LearnError(f"the folder {p.parent} does not exist")
    if p == Path(p.anchor) or p == Path.home().resolve():
        raise LearnError(f"choose a dedicated folder, not {p}")
    return p


def _fresh_dir(ws: Path) -> None:
    """Empty ``ws`` for a fresh start -- only when the course created it (marker) or it is empty."""
    if ws.is_symlink():
        raise LearnError(f"{ws} is a symbolic link; the course does not use it", 409)
    if ws.exists():
        if not ws.is_dir():
            raise LearnError(f"{ws} exists and is not a folder", 409)
        if not (ws / MARKER).is_file() and any(ws.iterdir()):
            raise LearnError(f"{ws} exists and was not created by the course; it is left untouched. Choose "
                             f"another course workspace folder (Learn > Course Workspace Folder).", 409)
        shutil.rmtree(ws)


def _mark(ws: Path, what: str) -> None:
    (ws / MARKER).write_text(f"SASSI-EDU guided course workspace ({what}); the GUI may delete and recreate "
                             f"this folder.\n", encoding="utf-8")


def prepare_lesson(lesson: L.Lesson, root: Path, examples_dir: Path) -> Path:
    """Fresh ``<root>/<lesson id>/`` with the lesson's example files (:func:`sassi.ui.lessons.prepare_workspace`)."""
    if not L.ID_RE.match(lesson.id):
        raise LearnError(f"invalid lesson id {lesson.id!r}")
    ws = Path(root) / lesson.id
    _fresh_dir(ws)
    try:
        ws = L.prepare_workspace(lesson, root, examples_dir, fresh=True)
    except FileNotFoundError as exc:
        raise LearnError(str(exc), 404) from None
    except OSError as exc:
        raise LearnError(f"cannot prepare {ws}: {exc}", 500) from None
    _mark(ws, f"lesson {lesson.id}")
    return ws


def prepare_example(name: str, root: Path, examples_dir: Path) -> Tuple[Path, Path]:
    """Fresh ``<root>/examples/<name>/`` with the example, ``data/`` and its support files; returns
    ``(workspace, copied .pre)``."""
    if not EXAMPLE_RE.match(name or "") or not (Path(examples_dir) / f"{name}.pre").is_file():
        raise LearnError(f"no example {name!r}", 404)
    ws = Path(root) / EXAMPLES_SUBDIR / name
    _fresh_dir(ws)
    try:
        pre = L.copy_example(name, ws, examples_dir)
    except OSError as exc:
        raise LearnError(f"cannot prepare {ws}: {exc}", 500) from None
    _mark(ws, f"example {name}")
    return ws, pre


def _tree_bytes(path: Path) -> int:
    total = 0
    for dirpath, _dirs, files in os.walk(path):
        for f in files:
            try:
                total += os.lstat(os.path.join(dirpath, f)).st_size
            except OSError:
                pass
    return total


def _contains(folder: Path, path: Optional[Path]) -> bool:
    if path is None:
        return False
    try:
        Path(path).resolve().relative_to(folder.resolve())
        return True
    except (ValueError, OSError):
        return False


def course_workspaces(root: Path, in_use: Sequence[Optional[Path]] = ()) -> List[Dict[str, Any]]:
    """The workspaces the course created under ``root`` (marker file): ``<root>/<lesson id>/`` and
    ``<root>/examples/<name>/``, with their size on disk; ``in_use`` when one of the ``in_use`` paths (the
    interpreter's working directory, the open lesson) lies inside."""
    root = Path(root)
    cands = [p for p in sorted(root.iterdir()) if p.is_dir()] if root.is_dir() else []
    ex = root / EXAMPLES_SUBDIR
    if ex.is_dir() and not (ex / MARKER).is_file():
        cands = [p for p in cands if p != ex] + [p for p in sorted(ex.iterdir()) if p.is_dir()]
    out = []
    for p in cands:
        if p.is_symlink() or not (p / MARKER).is_file():
            continue
        out.append({"name": p.relative_to(root).as_posix(), "path": str(p), "bytes": _tree_bytes(p),
                    "in_use": any(_contains(p, u) for u in in_use)})
    return out


def delete_workspaces(root: Path, names: Optional[Sequence[str]] = None,
                      in_use: Sequence[Optional[Path]] = ()) -> Dict[str, Any]:
    """Delete course workspaces under ``root`` (all of them, or those ``names`` of
    :func:`course_workspaces`) to free disk space.  Only folders carrying the course marker are deleted, and
    never one that is in use; returns ``{"deleted": [...], "skipped": [...], "freed": bytes}``."""
    ws = {w["name"]: w for w in course_workspaces(root, in_use)}
    want = list(ws) if names is None else list(names)
    deleted, skipped, freed = [], [], 0
    for n in want:
        w = ws.get(n)
        if w is None:
            skipped.append({"name": n, "reason": "not a course workspace"})
        elif w["in_use"]:
            skipped.append({"name": n, "reason": "in use (the working directory or the open lesson)"})
        else:
            try:
                shutil.rmtree(w["path"])
            except OSError as exc:
                skipped.append({"name": n, "reason": str(exc)})
                continue
            deleted.append(n)
            freed += w["bytes"]
    return {"deleted": deleted, "skipped": skipped, "freed": freed}


def fresh_model_lines(models: Sequence[int], active: int) -> List[str]:
    """Commands that leave one empty active model 0 (the learner's other models are removed from
    memory; their files are untouched): ``ACTM,0``, ``DMODEL,<n>`` ..., ``DMODEL,0``."""
    lines = [] if active == 0 else ["ACTM,0"]
    for n in sorted(set(models)):
        if n != 0:
            lines.append(f"DMODEL,{n}")
    lines.append("DMODEL,0")
    return lines


def cd_line(path: Path) -> str:
    return join_command("CD", [str(path)])


# ======================================================================================
# examples gallery
# ======================================================================================
def _header_lines(text: str) -> List[str]:
    out = []
    for ln in text.splitlines():
        if not ln.startswith("*"):
            break
        body = ln[1:]
        out.append(body[1:] if body.startswith(" ") else body)
    return out


def parse_example_header(text: str) -> Dict[str, Any]:
    """Title, number, "What you learn" bullets and the physics / model paragraph of an example's header."""
    lines = _header_lines(text)
    sep = [i for i, ln in enumerate(lines) if re.fullmatch(r"\s*=+\s*", ln)]
    title, number = "", None
    if len(sep) >= 2:
        title = " ".join(ln.strip() for ln in lines[sep[0] + 1:sep[1]]).strip()
        body = lines[sep[1] + 1:]
    else:
        body = lines
    m = re.match(r"^SASSI-EDU tutorial example (\d+)\s*:\s*(.*)$", title)
    if m:
        number, title = int(m.group(1)), m.group(2).strip()
    # sections: a non-indented line starts one ("What you learn", "The model: a rigid ...")
    sections: List[Tuple[str, List[str]]] = []
    for ln in body:
        if re.fullmatch(r"\s*=+\s*", ln):
            break
        if ln and not ln[0].isspace():
            head, colon, rest = ln.partition(":")
            if colon and len(head) <= 40:
                sections.append((head.strip(), [rest.strip()] if rest.strip() else []))
            else:
                sections.append((ln.strip(), []))
        elif sections:
            sections[-1][1].append(ln)
    learn: List[str] = []
    physics = ""
    for name, sl in sections:
        low = name.lower()
        if low.startswith("what you learn"):
            for ln in sl:
                s = ln.strip()
                if not s:
                    continue
                mm = re.match(r"^(?:-|\d+\.)\s+(.*)$", s)
                if mm or not learn:
                    learn.append(mm.group(1) if mm else s)
                else:
                    learn[-1] += " " + s
        elif not physics and (low.startswith("physics") or low.startswith("the model") or low.startswith("the structure")):
            physics = " ".join(ln.strip() for ln in sl if ln.strip())
    learn = [re.sub(r"\s+", " ", x).strip().rstrip(";") for x in learn]
    mdl = re.search(r"^MDL\s*,\s*([^,\s]+)\s*,\s*([^\s,]+)", text, re.M)
    return {"number": number, "title": title, "learn": learn, "physics": re.sub(r"\s+", " ", physics).strip(),
            "model": mdl.group(1) if mdl else "", "results": mdl.group(2) if mdl else ""}


def readme_rows(text: str) -> Dict[str, Dict[str, str]]:
    """``examples/README.md`` table: file stem -> {topic, modules, runtime}."""
    out: Dict[str, Dict[str, str]] = {}
    for ln in text.splitlines():
        m = re.match(r"^\|\s*`([^`]+)\.pre`\s*\|(.*)\|\s*$", ln.strip())
        if not m:
            continue
        cells = [c.strip() for c in m.group(2).split("|")]
        if len(cells) >= 3:
            out[m.group(1)] = {"topic": cells[0], "modules": cells[1], "runtime": cells[2]}
    return out


def examples(examples_dir: Path, lessons: Sequence[L.Lesson]) -> List[Dict[str, Any]]:
    """The gallery: one entry per ``examples/*.pre`` (sorted by name)."""
    d = Path(examples_dir)
    try:
        rows = readme_rows((d / "README.md").read_text(encoding="utf-8", errors="replace"))
    except OSError:
        rows = {}
    out = []
    for p in sorted(d.glob("*.pre")):
        if not EXAMPLE_RE.match(p.stem):
            continue
        try:
            text = p.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        h = parse_example_header(text)
        r = rows.get(p.stem, {})
        out.append({"name": p.stem, "file": p.name, "number": h["number"], "title": h["title"] or p.stem,
                    "learn": h["learn"], "physics": h["physics"], "model": h["model"], "results": h["results"],
                    "topic": r.get("topic", ""), "modules": r.get("modules", ""), "runtime": r.get("runtime", ""),
                    "lines": sum(1 for ln in text.splitlines() if ln.strip() and not ln.lstrip().startswith("*")),
                    "lessons": [{"id": les.id, "title": les.title, "part": les.part}
                                for les in lessons if les.example == p.stem]})
    return out


# ======================================================================================
# course outline and lesson rendering
# ======================================================================================
def outline(lessons: Sequence[L.Lesson]) -> Dict[str, Any]:
    parts: List[Dict[str, Any]] = []
    for part in list(PARTS) + sorted({x.part for x in lessons if x.part not in PARTS}):
        items = [x.to_dict(full=False) for x in lessons if x.part == part]
        if items:
            parts.append({"part": part, "lessons": items})
    return {"parts": parts, "count": len(lessons)}


def load_all(directory: Path) -> Tuple[List[L.Lesson], List[Dict[str, str]]]:
    """Every lesson of ``directory`` that parses; the files that do not are reported (not raised), so one
    lesson being written cannot hide the others."""
    d = Path(directory)
    good: List[L.Lesson] = []
    bad: List[Dict[str, str]] = []
    if not d.is_dir():
        return good, bad
    for p in sorted(d.glob("*.md")):
        try:
            good.append(L.parse_lesson(p.read_text(encoding="utf-8"), p))
        except (L.LessonError, OSError, UnicodeDecodeError) as exc:
            bad.append({"file": p.name, "error": str(exc)})
    order = {k.lower(): i for i, k in enumerate(PARTS)}
    good.sort(key=lambda x: (order.get(x.part.lower(), 9), x.order, x.id))
    return good, bad


def _fence_pass(md: str, start: int, keep_setup: bool = False) -> Tuple[str, List[Dict[str, Any]]]:
    """Replace the fenced blocks of ``md`` (numbered from ``start`` in document order, as
    :func:`sassi.ui.lessons._split_blocks` numbers them): command blocks -> placeholder paragraphs, action
    blocks -> removed, ``sassi-setup`` -> placeholder (``keep_setup``) or removed, other fences kept."""
    out: List[str] = []
    blocks: List[Dict[str, Any]] = []
    lines = md.split("\n")
    i, k = 0, start
    while i < len(lines):
        m = L.FENCE_RE.match(lines[i].strip())
        if m and lines[i].lstrip().startswith("```"):
            fence, lang = m.group(1), m.group(2).lower()
            j = i + 1
            body: List[str] = []
            while j < len(lines) and not lines[j].strip().startswith(fence):
                body.append(lines[j])
                j += 1
            kind = lang if lang in L.RUN_BLOCKS + L.SHOW_BLOCKS + L.SETUP_BLOCKS + ("action",) else "code"
            blk = {"k": k, "kind": kind, "lang": lang, "lines": [ln.rstrip() for ln in body if ln.strip()]}
            blocks.append(blk)
            if kind in ("sassi", "sassi-show") or (kind == "sassi-setup" and keep_setup):
                out.extend(["", _PLACE.format(k), ""])
            elif kind == "code":
                out.extend(lines[i:j + 1])
            k += 1
            i = j + 1
            continue
        out.append(lines[i])
        i += 1
    return "\n".join(out), blocks


class _Renderer:
    """Markdown of a lesson -> HTML (links to docs/... open in Help; images under docs/)."""

    def __init__(self, helpdocs):
        self.link, self.image = helpdocs._resolver("LESSON.md")      # repository-root relative paths

    def __call__(self, md: str) -> str:
        if not md.strip():
            return ""
        html = render(md, link=self.link, image=self.image).html
        return _PLACE_RE.sub(lambda m: f'<div class="lesson-block" data-block="{m.group(1)}"></div>', html)


def _line_names(files: Sequence[str]) -> List[str]:
    """Names of the plotted lines: the file names, or the workspace-relative paths when two files have the
    same name (``ex01/00085TR_X.TFI`` against ``ex01_fixed/00085TR_X.TFI``)."""
    base = [posixpath.basename(f.replace("\\", "/")) for f in files]
    if len(set(base)) == len(base):
        return base
    return [posixpath.normpath(f.replace("\\", "/")) for f in files]


def action_label(verb: str, args: str) -> str:
    """The text of an action button."""
    a = (args or "").strip()
    files = [f.strip() for f in a.split("|")[0].split(",") if f.strip()]
    names = ", ".join(_line_names(files))
    log = "|" in a and "log" in a.split("|", 1)[1].lower()
    return {
        "plot-model": "Plot the model (elements)",
        "plot-nodes": "Plot the nodes",
        "plot-layers": "Plot the soil layers",
        "plot-soilprops": f"Plot soil curves {a}".strip(),
        "plot-spectrum": f"Plot {names}" + (" (log f)" if log else ""),
        "plot-history": f"Plot history {names}",
        "open-file": f"Open {posixpath.basename(a) or a}",
        "open-listing": f"Open the {a.upper()} listing",
        "open-dialog": "Open " + ("Options > Analysis > " + a.split("/", 1)[1].upper() if "/" in a
                                  else {"ANALYSIS": "Options > Analysis", "MODEL": "Options > Model",
                                        "WRITE": "Options > Write", "CHECK": "Options > Check"}.get(a.upper(), a)),
        "open-doc": "Read: " + (a.split("#", 1)[0].rsplit("/", 1)[-1] + (" §" + a.split("#", 1)[1] if "#" in a else "")),
        "explain": f"Explain {a}",
    }.get(verb, f"{verb} {a}".strip())


def _doc_label(args: str, helpdocs) -> str:
    """``Read: <document title> -- <section heading>`` of an ``open-doc`` action (the raw path and anchor
    when the document cannot be read)."""
    path, _, anchor = (args or "").strip().partition("#")
    try:
        doc = helpdocs.document(path)
    except Exception:                                    # noqa: BLE001 - fall back to the plain label
        return action_label("open-doc", args)
    head = next((h["text"] for h in doc.get("toc", []) if h.get("slug") == anchor), "") if anchor else ""
    return f"Read: {doc['title']}" + (f" — {head}" if head else "")


def render_lesson(lesson: L.Lesson, helpdocs) -> Dict[str, Any]:
    """The lesson as the lesson player shows it (see the module docstring)."""
    R = _Renderer(helpdocs)
    d = lesson.to_dict(full=False)
    intro_md, intro_blocks = _fence_pass(lesson.intro, 0, keep_setup=True)
    d["intro_html"] = R(intro_md)
    d["intro_blocks"] = [b for b in intro_blocks if b["kind"] in ("sassi", "sassi-show", "sassi-setup")]
    d["setup"] = list(lesson.setup)
    steps = []
    for st in lesson.steps:
        k = 0
        narr_md, blocks = _fence_pass(st.narrative, k)
        k += len(blocks)
        sections = []
        for sec in st.sections:
            smd, sb = _fence_pass(sec.markdown, k)
            k += len(sb)
            blocks.extend(sb)
            item = {"key": sec.key, "title": sec.title or SECTION_TITLES.get(sec.key, ""), "html": ""}
            if sec.key == "check":
                q, a = L.split_answer(smd)
                item["html"], item["answer_html"] = R(q), R(a)
            else:
                item["html"] = R(smd)
            sections.append(item)
        steps.append({"index": st.index, "title": st.title, "narrative_html": R(narr_md), "sections": sections,
                      "blocks": [b for b in blocks if b["kind"] in ("sassi", "sassi-show")],
                      "commands": st.commands,
                      "actions": [{"verb": v, "args": a, "label": _doc_label(a, helpdocs) if v == "open-doc" else action_label(v, a),
                                   "known": v in L.ACTION_VERBS} for v, a in st.actions]})
    d["steps"] = steps
    return d


# ======================================================================================
# actions
# ======================================================================================
def _ws_file(ws: Path, rel: str) -> Path:
    """``rel`` (relative to the lesson workspace) as an absolute path inside it."""
    s = (rel or "").strip().replace("\\", "/")
    if not s or "\x00" in s or s.startswith("/") or re.match(r"^[A-Za-z]:", s):
        raise LearnError(f"action path {rel!r} must be relative to the lesson workspace")
    norm = posixpath.normpath(s)
    if norm == ".." or norm.startswith("../"):
        raise LearnError(f"action path {rel!r} leaves the lesson workspace")
    return Path(ws) / norm


def _cmd_path(interp, p: Path) -> str:
    """The shortest spelling of file ``p`` that the interpreter resolves to it (L15): relative to the
    working directory when that resolves to the same file, else absolute."""
    try:
        rel = os.path.relpath(str(p), str(interp.cwd))
        if not rel.startswith("..") and Path(interp.resolve_path(rel)).resolve() == p.resolve():
            return rel
    except Exception:            # noqa: BLE001 -- any doubt: the absolute path
        pass
    return str(p)


def _listing(interp, module: str) -> Path:
    m = interp.model
    if not m.name or not m.path:
        raise LearnError("the active model has no name and folder yet (MDL): run the step first", 404)
    d = Path(m.path)
    for cand in (f"{m.name}_{module.upper()}.out", f"{m.name}_{module.lower()}.out"):
        if (d / cand).is_file():
            return d / cand
    if d.is_dir():
        want = f"{m.name}_{module}.out".lower()
        for p in d.iterdir():
            if p.name.lower() == want:
                return p
    raise LearnError(f"{m.name}_{module.upper()}.out not found in {d}: run the step (RUN{module.upper()}) first", 404)


def resolve_action(verb: str, args: str, ws: Path, interp, lines_in_memory: Sequence[int], helpdocs) -> Dict[str, Any]:
    """An action of a lesson step -> ``{"kind": "commands", "lines": [...]}`` (command text, L17) or a GUI
    request: ``{"kind": "file", "path"}``, ``{"kind": "dialog", "name", "tab"}``, ``{"kind": "doc", "doc",
    "anchor"}``, ``{"kind": "explain", "line"}``."""
    verb = (verb or "").strip().lower()
    a = (args or "").strip()
    if verb not in L.ACTION_VERBS:
        raise LearnError(f"unknown action {verb!r}")
    if verb == "plot-model":
        return {"kind": "commands", "lines": ["MODELPLOT"]}
    if verb == "plot-nodes":
        return {"kind": "commands", "lines": ["NODEPLOT"]}
    if verb == "plot-layers":
        return {"kind": "commands", "lines": ["LAYERPLOT"]}
    if verb == "plot-soilprops":
        if not a:
            raise LearnError("plot-soilprops needs the DYNP label")
        return {"kind": "commands", "lines": [join_command("SOILPROPPLOT", [a])]}
    if verb in ("plot-spectrum", "plot-history"):
        files_part, _, opt = a.partition("|")
        files = [f.strip() for f in files_part.split(",") if f.strip()]
        if not files:
            raise LearnError(f"{verb} needs at least one file")
        nxt = max([0] + [int(n) for n in lines_in_memory]) + 1
        lines: List[str] = []
        nums: List[int] = []
        names = _line_names(files)
        rename = names != [posixpath.basename(f.replace("\\", "/")) for f in files]
        for f, label in zip(files, names):
            p = _ws_file(ws, f)
            if not p.is_file():
                raise LearnError(f"{f} not found in the lesson workspace {ws}: run the step first", 404)
            spell = _cmd_path(interp, p)
            if verb == "plot-spectrum":
                lines.append(join_command("READSPEC", [spell, "1", str(nxt)]))
            else:
                pair = int(file_info(p).get("pair", 0) or 0)
                lines.append(join_command("READTH", [spell, str(pair), str(nxt)]))
            if rename:                     # files of the same name: the legend shows where each comes from
                lines.append(join_command("LINENAME", [str(nxt), label]))
            nums.append(nxt)
            nxt += 1
        lines.append(("SPECPLOT," if verb == "plot-spectrum" else "THPLOT,") + ",".join(str(n) for n in nums))
        opts = {o.strip().lower() for o in opt.split(",") if o.strip()}
        if "log" in opts or "logx" in opts:
            lines.append("AXES,,,,,1")
        if "loglog" in opts:
            lines.append("AXES,,,,,1,1")
        return {"kind": "commands", "lines": lines}
    if verb == "open-file":
        p = _ws_file(ws, a)
        if not p.is_file():
            raise LearnError(f"{a} not found in the lesson workspace {ws}: run the step first", 404)
        return {"kind": "file", "path": str(p)}
    if verb == "open-listing":
        if not MODULE_RE.match(a):
            raise LearnError(f"open-listing needs a module name (got {a!r})")
        return {"kind": "file", "path": str(_listing(interp, a))}
    if verb == "open-dialog":
        name, _, tab = a.upper().partition("/")
        name = name.strip()
        from .dialogs import TABS
        if name in TABS and not tab:
            name, tab = "ANALYSIS", name
        if name not in DIALOGS or (tab and tab not in TABS):
            raise LearnError(f"open-dialog: unknown dialog {a!r}")
        return {"kind": "dialog", "name": name, "tab": tab.strip()}
    if verb == "open-doc":
        doc, _, anchor = a.partition("#")
        rel = helpdocs.normalize(doc)
        if rel is None or not helpdocs.is_readable(rel):
            raise LearnError(f"open-doc: {doc} is not a help document", 404)
        return {"kind": "doc", "doc": rel, "anchor": anchor}
    if verb == "explain":
        if not a:
            raise LearnError("explain needs a command line")
        return {"kind": "explain", "line": a}
    raise LearnError(f"unknown action {verb!r}")       # pragma: no cover
