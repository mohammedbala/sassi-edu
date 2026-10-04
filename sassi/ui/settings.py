"""Persistent GUI settings: ``SASSIini.xml`` and ``SASSIdb.xml`` (requirements 5.10, UI-07, D-UI-07).

Both files live in the per-user settings directory (:func:`sassi.plotting.state.settings_dir`:
``~/Library/Application Support/SASSI-EDU`` on macOS; ``SASSI_EDU_SETTINGS_DIR`` overrides it).

``SASSIini.xml`` -- written by OK of the Modules > Location and Modules > Extension dialogs, by the
View > Command Display toggles and at Model > Exit:

* module locations (Modules > Location, spec 04 section 15.1): ``built-in`` (default, the Python
  module of this package) or the path of an external executable that reads the three-line batch
  protocol (for cross-checks against the commercial modules);
* file extensions (Modules > Extension, spec 04 section 15.2);
* Command Display filters (View > Command Display, persisted, D-UI-07) and the six message colours
  (requirements 5.6);
* the folder of the guided-course workspaces (Learn > Course Workspace Folder; empty = the default
  ``<GUI start directory>/sassi-course``);
* the shader options (Options > Shader Options = ``SHADEROPTIONS``: node size, element outline,
  element shrink, vector scale; UI-07, D-UI-07, 5.10) -- written when they change (OK of the Shader
  dialog submits SHADEROPTIONS) and at Exit, restored into the plot state at start-up;

Check Options and toolbar visibility are **not** persisted (manual fidelity, UI-07).

``SASSIdb.xml`` -- the two-level Load Model tree (Model > Open, spec 04 section 3.3): groups and,
per model, only its name, location and title (D-UI-05).  No model data is stored there.
"""
from __future__ import annotations

import os
import xml.etree.ElementTree as ET
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, List, Optional

from ..plotting.state import ShaderOptions, settings_dir

#: Modules > Location rows (spec 04 section 15.1); disabled rows are the modules not in V3
LOCATION_ROWS = ("EQUAKE", "SOIL", "LIQUEF", "SITE", "POINT", "HOUSE", "PINT", "FORCE", "ANALYS", "COMBIN",
                 "MOTION", "STRESS", "RELDISP", "LOADGEN", "NONLINEAR", "SASSIANSYS")
LOCATION_DISABLED = ("LIQUEF", "PINT")
BUILT_IN = "built-in"

#: Modules > Extension defaults (spec 04 section 15.2): module -> (input ext, output ext)
DEFAULT_EXTENSIONS: Dict[str, List[str]] = {
    "EQUAKE": [".equ", "_equake.out"], "SOIL": [".soi", "_soil.out"], "SITE": [".sit", "_site.out"],
    "POINT": [".poi", "_point.out"], "HOUSE": [".hou", "_house.out"], "FORCE": [".frc", "_force.out"],
    "ANALYS": [".anl", "_analys.out"], "COMBIN": ["", "_combin.out"], "MOTION": [".mot", "_motion.out"],
    "STRESS": [".str", "_stress.out"], "RELDISP": [".rdi", "_reldisp.out"], "NONLINEAR": [".eql", "_nonlinear.out"],
    "LOADGEN": [".lgn", "_loadgen.out"],
}

#: View > Command Display groups (sassi.prep.messages.MessageSink.FILTER_GROUPS); info is always shown
DISPLAY_GROUPS = ("echo", "confirm", "comments", "warnerr")
#: default Command History colours (requirements 5.6)
DEFAULT_COLOURS = {"ECHO": "#000000", "CONFIRM": "#0000c0", "COMMENT": "#008000", "INFO": "#800080",
                   "WARNING": "#b07800", "ERROR": "#d00000"}
#: shader options (sassi.plotting.state.ShaderOptions: SHADEROPTIONS,<points>,<linew>,<shrink>,<scale>)
SHADER_FIELDS = ("points", "linew", "shrink", "scale")


def shader_value_ok(name: str, value: float) -> bool:
    """The SHADEROPTIONS ranges: outline and shrink are fractions in [0, 0.5); size and scale > 0."""
    if value != value:                                   # NaN
        return False
    if name in ("linew", "shrink"):
        return 0.0 <= value < 0.5
    return value > 0.0


def _atomic_write(path: Path, root: ET.Element) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    ET.ElementTree(root).write(str(tmp), encoding="utf-8", xml_declaration=True)
    os.replace(tmp, path)


class IniSettings:
    """``SASSIini.xml``: module locations, extensions, Command Display filters, colours, shader options."""

    def __init__(self, directory: Optional[Path] = None):
        self.directory = Path(directory) if directory is not None else settings_dir()
        self.locations: Dict[str, str] = {m: BUILT_IN for m in LOCATION_ROWS}
        self.extensions: Dict[str, List[str]] = {k: list(v) for k, v in DEFAULT_EXTENSIONS.items()}
        self.display: Dict[str, bool] = {g: True for g in DISPLAY_GROUPS}
        self.colours: Dict[str, str] = dict(DEFAULT_COLOURS)
        self.shader: Dict[str, float] = {k: float(v) for k, v in asdict(ShaderOptions()).items()}
        #: Learn (guided course): the folder of the lesson / example workspaces ("" = <GUI start dir>/sassi-course)
        self.course: Dict[str, str] = {"root": ""}
        self.load()

    @property
    def path(self) -> Path:
        return self.directory / "SASSIini.xml"

    # ------------------------------------------------------------------ I/O
    def load(self) -> None:
        if not self.path.exists():
            return
        try:
            root = ET.parse(str(self.path)).getroot()
        except (ET.ParseError, OSError):
            return                      # a corrupt file never stops the GUI; defaults are used
        for el in root.findall("ModuleLocations/Module"):
            name = el.get("name", "").upper()
            if name in self.locations:
                self.locations[name] = el.get("path", BUILT_IN) or BUILT_IN
        for el in root.findall("Extensions/Module"):
            name = el.get("name", "").upper()
            if name in self.extensions:
                self.extensions[name] = [el.get("input", self.extensions[name][0]),
                                         el.get("output", self.extensions[name][1])]
        cd = root.find("CommandDisplay")
        if cd is not None:
            for g in DISPLAY_GROUPS:
                if cd.get(g) is not None:
                    self.display[g] = cd.get(g) not in ("0", "false", "False")
        for el in root.findall("Colors/Color"):
            k = el.get("name", "").upper()
            if k in self.colours and el.get("value"):
                self.colours[k] = el.get("value")
        sh = root.find("Shader")
        if sh is not None:
            for k in SHADER_FIELDS:
                try:
                    v = float(sh.get(k, ""))
                except ValueError:
                    continue                    # a damaged value keeps the default
                if shader_value_ok(k, v):
                    self.shader[k] = v
        co = root.find("Course")
        if co is not None:
            self.course["root"] = co.get("root", "") or ""

    def save(self) -> Path:
        root = ET.Element("SASSIini", {"version": "1"})
        ml = ET.SubElement(root, "ModuleLocations")
        for m in LOCATION_ROWS:
            ET.SubElement(ml, "Module", {"name": m, "path": self.locations.get(m, BUILT_IN)})
        ex = ET.SubElement(root, "Extensions")
        for m, (i, o) in self.extensions.items():
            ET.SubElement(ex, "Module", {"name": m, "input": i, "output": o})
        ET.SubElement(root, "CommandDisplay", {g: "1" if self.display[g] else "0" for g in DISPLAY_GROUPS})
        co = ET.SubElement(root, "Colors")
        for k, v in self.colours.items():
            ET.SubElement(co, "Color", {"name": k, "value": v})
        ET.SubElement(root, "Shader", {k: repr(float(self.shader[k])) for k in SHADER_FIELDS})
        ET.SubElement(root, "Course", {"root": self.course.get("root", "")})
        _atomic_write(self.path, root)
        return self.path

    # ------------------------------------------------------------------ edits
    def external_module(self, module: str) -> Optional[str]:
        """Path of an external executable configured for ``module`` (None = built-in)."""
        p = self.locations.get(module.upper(), BUILT_IN)
        return None if not p or p == BUILT_IN else p

    def set_locations(self, rows: Dict[str, str]) -> List[str]:
        """Modules > Location OK: returns messages for rejected entries (missing executables)."""
        errs = []
        for m, p in rows.items():
            m = m.upper()
            if m not in self.locations or m in LOCATION_DISABLED:
                continue
            p = (p or "").strip() or BUILT_IN
            if p != BUILT_IN and not Path(p).expanduser().is_file():
                errs.append(f"{m}: {p} is not an existing file (kept {self.locations[m]})")
                continue
            self.locations[m] = p if p == BUILT_IN else str(Path(p).expanduser())
        return errs

    def set_extensions(self, rows: Dict[str, List[str]]) -> None:
        """Modules > Extension OK: all edits of the dialog are saved together (spec 04 section 15.2)."""
        for m, v in rows.items():
            m = m.upper()
            if m in self.extensions and isinstance(v, (list, tuple)) and len(v) == 2:
                self.extensions[m] = [str(v[0]).strip(), str(v[1]).strip()]

    def to_dict(self) -> Dict[str, Any]:
        return {"path": str(self.path), "locations": dict(self.locations), "location_rows": list(LOCATION_ROWS),
                "location_disabled": list(LOCATION_DISABLED), "extensions": {k: list(v) for k, v in self.extensions.items()},
                "display": dict(self.display), "colours": dict(self.colours), "shader": dict(self.shader),
                "course": dict(self.course)}


class ModelDatabase:
    """``SASSIdb.xml``: groups of models (name, location, title) for the Load Model dialog."""

    def __init__(self, directory: Optional[Path] = None):
        self.directory = Path(directory) if directory is not None else settings_dir()
        self.groups: List[Dict[str, Any]] = []
        self.load()

    @property
    def path(self) -> Path:
        return self.directory / "SASSIdb.xml"

    def load(self) -> None:
        self.groups = []
        if not self.path.exists():
            return
        try:
            root = ET.parse(str(self.path)).getroot()
        except (ET.ParseError, OSError):
            return
        for g in root.findall("Group"):
            models = [{"name": m.get("name", ""), "path": m.get("path", ""), "title": m.get("title", "")}
                      for m in g.findall("Model")]
            self.groups.append({"name": g.get("name", ""), "models": models})

    def save(self) -> Path:
        root = ET.Element("SASSIdb", {"version": "1"})
        for g in self.groups:
            ge = ET.SubElement(root, "Group", {"name": g["name"]})
            for m in g["models"]:
                ET.SubElement(ge, "Model", {"name": m["name"], "path": m["path"], "title": m.get("title", "")})
        _atomic_write(self.path, root)
        return self.path

    def group(self, name: str) -> Optional[Dict[str, Any]]:
        for g in self.groups:
            if g["name"] == name:
                return g
        return None

    def add_group(self, name: str) -> None:
        name = name.strip()
        if not name:
            raise ValueError("group name required")
        if self.group(name) is not None:
            raise ValueError(f"group {name} already exists")
        self.groups.append({"name": name, "models": []})        # appended at the bottom (spec 04 3.3)
        self.save()

    def remove_group(self, name: str) -> Dict[str, Any]:
        g = self.group(name)
        if g is None:
            raise ValueError(f"group {name} not found")
        self.groups.remove(g)
        self.save()
        return g

    def add_model(self, group: str, name: str, path: str, title: str = "") -> Dict[str, Any]:
        """Add Model (D-UI-05): group, model name, model directory (created if missing), title."""
        g = self.group(group)
        if g is None:
            raise ValueError("select a group first (a model must be in a group)")
        name = name.strip()
        if not name or any(c in name for c in "/\\,"):
            raise ValueError("model name required (no '/', '\\' or ',')")
        p = Path(path).expanduser()
        if not p.is_absolute():
            raise ValueError("the model directory must be an absolute path")
        p.mkdir(parents=True, exist_ok=True)
        if any(m["name"] == name and Path(m["path"]) == p for m in g["models"]):
            raise ValueError(f"model {name} in {p} is already in group {group}")
        entry = {"name": name, "path": str(p), "title": title}
        g["models"].append(entry)
        self.save()
        return entry

    def remove_model(self, group: str, name: str, path: str) -> Dict[str, Any]:
        g = self.group(group)
        if g is None:
            raise ValueError(f"group {group} not found")
        for m in g["models"]:
            if m["name"] == name and m["path"] == path:
                g["models"].remove(m)
                self.save()
                return m
        raise ValueError(f"model {name} not found in group {group}")


def model_files(name: str, path: str) -> List[Path]:
    """Files of model ``name`` in its directory that the second Remove prompt deletes.

    Only files named after the model (``<name>.*``, ``<name>_*``, ``<name>-Sim.pre``) are deleted;
    other files and sub-directories of the folder are left on disk (SASSI-EDU safety choice: model
    folders are often shared with other data).
    """
    d = Path(path)
    if not name or not d.is_dir():
        return []
    out = []
    for p in sorted(d.iterdir()):
        if p.is_file() and (p.name.startswith(name + ".") or p.name.startswith(name + "_")
                            or p.name == f"{name}-Sim.pre"):
            out.append(p)
    return out
