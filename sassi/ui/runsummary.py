"""Run summary: the key inputs, the key outputs and their graphs of a model after a run (the GUI's summary window).

``GET /api/run_summary?model=<name or number>`` (default: the active model) returns what the window shows.  It is
read from the model's folder, so it describes what the modules actually used and wrote:

* **modules** -- the listings ``<model>_<MODULE>.out``: status, time, warnings and errors of each module run;
* **inputs** -- from the decks AFWRITE wrote (the resolved values every module read, :mod:`sassi.io.decks`):
  the model (nodes, elements by type, interaction nodes, embedment, masses from the HOUSE listing), the soil
  profile, the frequencies, the analysis, the input motion (file, scaling, peak, duration, spectrum damping), the
  loads and the output requests;
* **outputs** -- from the result files: the transfer functions (``.TFI``: the largest amplitude and its
  frequency, i.e. the SSI resonance), the peak accelerations (``.ACC``), the in-structure response spectra
  (``.RS``: peak and zero-period acceleration), the relative displacements (``.THD``) and the element stress
  histories (``.THS``), the SOIL module's spectra and histories;
* **charts** -- the transfer-function amplitudes, the response spectra (with the spectrum of the input motion,
  computed here), the acceleration histories (input and the largest response) and the soil profile; each lists
  the result files it draws, so the window can open them as plots (``POST /api/run_summary/plot``: the READSPEC /
  READTH / SPECPLOT / THPLOT command text of a lesson's plot actions, rule L17).

Nothing here changes the model or the files.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..conventions import ELEMENT_TYPE_NAMES
from ..io import decks as D
from ..io.textfiles import read_history, read_tf, read_xy

#: the order in which the modules run (the summary lists them so)
MODULE_ORDER = ("EQUAKE", "SOIL", "SITE", "POINT", "HOUSE", "FORCE", "ANALYS", "COMBIN", "MOTION", "STRESS",
                "RELDISP", "NONLINEAR")
STATUS = re.compile(r"^\s*([A-Z]+) finished with status (\w+) in ([\d.]+) s; (\d+) warning\(s\), (\d+) error\(s\)")
#: nodal result files of MOTION / RELDISP: <node>TR_<dir>[<damping index>].<ext> (RO_ for rotations)
NODAL = re.compile(r"^(\d{5,6})(TR|RO)_([XYZ]{1,2})(\d\d)?\.(TFI|TFU|ACC|RS|THD|TFD)$", re.IGNORECASE)
#: element stress histories of STRESS: <TYPE>_<group>_<element>_<component>.THS
ELEMENT = re.compile(r"^([A-Z]+)_(\d+)_(\d+)_(\w+?)\.THS$", re.IGNORECASE)
MASS = re.compile(r"Total (structural|excavated soil) mass(?: X Y Z)?\s*:\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)"
                  r"(?:\s+\(weight ([-\d.eE+]+)\))?")
WAVES = {1: "Rayleigh", 2: "SV", 3: "P", 4: "SH", 5: "Love"}
DIRS = "XYZ"
MAX_LINES = 6                    # lines per chart
MAX_POINTS = 1500                # points per history line (min / max decimation keeps the peaks)


class SummaryError(ValueError):
    """The model or the request cannot be served (status 404 / 400)."""

    def __init__(self, message: str, status: int = 404):
        super().__init__(message)
        self.status = status


def _g(v: float, digits: int = 4) -> str:
    """Short number (4 significant digits)."""
    v = float(v)
    if v != 0 and (abs(v) < 1e-3 or abs(v) >= 1e6):
        return f"{v:.{digits - 1}e}"
    return f"{v:.{digits}g}"


def _units(gravity: float) -> str:
    """The length unit the gravity implies (the model's consistent system)."""
    for g, u in ((9.81, "m"), (32.2, "ft"), (386.4, "in"), (981.0, "cm"), (9810.0, "mm")):
        if abs(gravity - g) <= 0.01 * g:
            return u
    return "length units"


def find_model(interp, model: Optional[str] = None) -> Tuple[int, Any]:
    """The model by name or number (default: the active model)."""
    if model in (None, "", "-1"):
        n = interp.active_model
        return n, interp.models[n]
    s = str(model).strip()
    if s.lstrip("-").isdigit() and int(s) in interp.models:
        return int(s), interp.models[int(s)]
    for n, m in sorted(interp.models.items()):
        if m.name == s:
            return n, m
    raise SummaryError(f"model {s} is not in memory")


def _folder(m) -> Optional[Path]:
    p = Path(m.path) if m.path else None
    return p if p is not None and p.is_dir() else None


def module_runs(folder: Path, name: str) -> List[Dict[str, Any]]:
    """The module listings of the model's folder: status, seconds, warnings, errors, time of the run."""
    out = []
    for mod in MODULE_ORDER:
        p = folder / f"{name}_{mod}.out"
        if not p.is_file():
            continue
        status = None
        try:
            for ln in p.read_text(encoding="utf-8", errors="replace").splitlines()[-40:]:
                m = STATUS.match(ln)
                if m:
                    status = m
        except OSError:
            continue
        out.append({"module": mod, "listing": p.name, "mtime": p.stat().st_mtime,
                    "ok": bool(status and status.group(2).upper() == "OK"),
                    "seconds": float(status.group(3)) if status else None,
                    "warnings": int(status.group(4)) if status else None,
                    "errors": int(status.group(5)) if status else None})
    return out


def _deck(folder: Path, name: str, module: str) -> Optional[D.Deck]:
    p = D.deck_path(folder, name, module)
    if not p.is_file():
        return None
    try:
        return D.read(p, module)
    except Exception:                       # noqa: BLE001 -- an unreadable deck is left out of the summary
        return None


def _resolve(folder: Path, name: str) -> Optional[Path]:
    """A file named in a deck (relative to the model folder, as the modules read it)."""
    if not name:
        return None
    p = Path(name)
    p = p if p.is_absolute() else folder / p
    try:
        p = p.resolve()
    except OSError:
        return None
    return p if p.is_file() else None


def input_motion(folder: Path, deck: Optional[D.Deck]) -> Optional[Dict[str, Any]]:
    """The control motion of a MOTION / STRESS / RELDISP deck, scaled as the module scales it (g)."""
    if deck is None:
        return None
    pr = deck.params
    p = _resolve(folder, str(pr.get("thfile") or ""))
    if p is None:
        return None
    from ..io.thfile import read_history as read_th
    try:
        acc, dt = read_th(p, int(pr.get("fopt", 0) or 0), int(pr.get("rec1", 1) or 1), int(pr.get("rec2", 0) or 0))
    except (OSError, ValueError):
        return None
    if dt is None or not len(acc):
        return None
    peak = float(np.max(np.abs(acc)))
    mx, mult = float(pr.get("max", 0) or 0), float(pr.get("mult", 0) or 0)
    scale = (mx / peak if peak > 0 else 1.0) if mx > 0 else (mult if mult else 1.0)
    acc = acc * scale
    return {"file": p.name, "path": str(p), "title": str(pr.get("thtit") or ""), "dt": float(dt), "acc": acc,
            "pga": float(np.max(np.abs(acc))), "duration": float(dt) * len(acc), "scale": scale}


# ======================================================================================
# inputs
# ======================================================================================
def _seismic(decks: Dict[str, Optional[D.Deck]]) -> bool:
    anl = decks.get("ANALYS")
    return anl is None or int(anl.params.get("type", 0) or 0) == 0


def _inputs(m, folder: Path, decks: Dict[str, Optional[D.Deck]], motion: Optional[Dict[str, Any]],
            units: str) -> List[Dict[str, Any]]:
    sections: List[Dict[str, Any]] = []
    hou, sit, anl, mot = decks.get("HOUSE"), decks.get("SITE"), decks.get("ANALYS"), decks.get("MOTION")
    seismic = _seismic(decks)

    # ---- the model
    rows: List[List[str]] = []
    if m.title:
        rows.append(["Title", m.title])
    c = m.counts()
    types: Dict[str, int] = {}
    for g, e in m.iter_elements():
        t = ELEMENT_TYPE_NAMES.get(g.type, str(g.type))
        types[t] = types.get(t, 0) + 1
    if c.get("elements", 0):
        rows.append(["Nodes / elements", f"{c.get('nodes', 0)} / {c.get('elements', 0)}"
                     + (f" ({', '.join(f'{n} {t}' for t, n in sorted(types.items()))})" if types else "")])
    else:
        rows.append(["Structure", "none: a free-field (site response) analysis"])
    inter = [n for n, nd in m.nodes.items() if 0 in nd.flags]
    if inter:
        _, xyz = m.global_coordinates(inter)
        xyz = np.asarray(xyz, float).reshape(-1, 3)
        gelev = float(m.ground_elevation)
        embed = max(0.0, gelev - float(xyz[:, 2].min()))
        rows.append(["Interaction nodes", f"{len(inter)} ({3 * len(inter)} DOFs)"
                     + (f"; embedment {_g(embed)} {units}" if embed > 1e-9 else "; surface foundation")])
    if hou is not None:
        meth = {0: "FV (flexible volume)", 1: "FFV", 2: "FI (subtraction)"}.get(int(hou.params.get("imp", 0) or 0), "")
        if meth:
            rows.append(["Method (HOUSE <imp>)", meth])
        rows.append(["Gravity, ground elevation", f"{_g(hou.params.get('gravity', 0))}, {_g(hou.params.get('gelev', 0))}"])
    lst = folder / f"{m.name}_HOUSE.out"
    if lst.is_file():
        try:
            for ln in lst.read_text(encoding="utf-8", errors="replace").splitlines():
                mm = MASS.search(ln)
                if mm:
                    what = "Structure mass X Y Z" if mm.group(1) == "structural" else "Excavated soil mass"
                    vals = f"{_g(mm.group(2))} {_g(mm.group(3))} {_g(mm.group(4))}"
                    rows.append([what, vals + (f" (weight {_g(mm.group(5))})" if mm.group(5) else "")])
        except OSError:
            pass
    sections.append({"title": "Model", "rows": rows})

    # ---- the soil profile and the frequencies (SITE deck)
    if sit is not None:
        lay = sit.table("layers").rows if sit.table("layers") is not None else []
        hs = sit.table("halfspace").rows if sit.table("halfspace") is not None else []
        rows = []
        if lay:
            depth = sum(float(r[1]) for r in lay)
            vs = [float(r[4]) for r in lay]
            ds = [float(r[6]) for r in lay]
            rows.append(["Layers", f"{len(lay)} to depth {_g(depth)} {units}"])
            rows.append(["Vs of the layers", f"{_g(min(vs))}" + (f" to {_g(max(vs))}" if max(vs) > min(vs) else "")
                         + f" {units}/s; damping {_g(100 * min(ds))}" + (f" to {_g(100 * max(ds))}" if max(ds) > min(ds) else "")
                         + " %"])
        if hs:
            r = hs[0]
            rows.append(["Half-space", f"Vs {_g(r[4])} {units}/s, damping {_g(100 * float(r[6]))} %"
                         + ("" if int(sit.params.get("nl", 20) or 0) else " (rigid base)")])
        waves = [WAVES.get(int(w[0]), str(w[0])) for w in (sit.table("waves").rows if sit.table("waves") is not None else [])
                 if int(w[1])]
        rows.append(["Waves; control point", f"{', '.join(waves) or 'none'}; top of layer {int(sit.params.get('cl', 1))}, "
                     f"direction {DIRS[int(sit.params.get('cm', 0) or 0)]}"])
        sections.append({"title": "Soil profile (SITE)", "rows": rows})
        fr = sit.table("freqs").rows if sit.table("freqs") is not None else []
        df = float(sit.params.get("df", 0) or 0) or 1.0 / (float(sit.params.get("delt", 0.005)) * int(sit.params.get("nft", 4096)))
        rows = []
        if fr:
            nums = [int(r[0]) for r in fr]
            rows.append(["SSI frequencies", f"{len(nums)} from {_g(min(nums) * df)} to {_g(max(nums) * df)} Hz"])
        rows.append(["NFFT, Δt, Δf", f"{int(sit.params.get('nft', 0))}, {_g(sit.params.get('delt', 0))} s, {_g(df)} Hz"])
        sections.append({"title": "Frequencies", "rows": rows})

    # ---- the analysis, the loads, the input motion and the output requests
    rows = []
    if anl is not None:
        seismic = int(anl.params.get("type", 0) or 0) == 0
        rows.append(["Analysis", "seismic (free-field motion)" if seismic else "foundation vibration (loads)"])
        sim = int(anl.params.get("simul", 0) or 0)
        if sim:
            rows.append(["Simultaneous cases", str(sim)])
    frc = decks.get("FORCE")
    if frc is not None and frc.table("loads") is not None and frc.table("loads").rows:
        loads = frc.table("loads").rows
        rows.append(["Loads (FORCE)", ", ".join(f"node {int(r[0])} {('FX', 'FY', 'FZ', 'MX', 'MY', 'MZ')[int(r[1]) - 1]} × {_g(r[2])}"
                                                for r in loads[:4]) + (" ..." if len(loads) > 4 else "")])
    if motion is not None:
        what, unit = ("Input motion", " g") if seismic else ("Load history", " (load factor)")
        rows.append([what, motion["file"] + (f" ({motion['title']})" if motion["title"] else "")])
        rows.append(["Peak, duration", f"{_g(motion['pga'])}{unit}, {_g(motion['duration'])} s (Δt {_g(motion['dt'])} s)"
                     + (f"; scaled × {_g(motion['scale'])}" if abs(motion["scale"] - 1) > 1e-9 else "")])
    if mot is not None:
        damp = [float(r[0]) for r in (mot.table("damp").rows if mot.table("damp") is not None else [])]
        if damp:
            rows.append(["Spectrum damping", ", ".join(f"{_g(100 * d)} %" for d in damp)])
        nout = mot.table("nout").rows if mot.table("nout") is not None else []
        if nout:
            rows.append(["Output requests (NOUT)", f"{len({int(r[0]) for r in nout})} nodes, {len(nout)} node-directions"])
    if rows:
        sections.append({"title": "Analysis and input motion", "rows": rows})

    # ---- the free-field modules
    soi = decks.get("SOIL")
    if soi is not None:
        prof = soi.table("profile").rows if soi.table("profile") is not None else []
        rows = [["Sublayers (SOIL)", f"{max(0, len(prof) - 1)} over a half-space"],
                ["Input history", f"{Path(str(soi.params.get('thfile') or '')).name or '-'}, peak {_g(soi.params.get('max', 0))} g"
                 if float(soi.params.get("max", 0) or 0) else Path(str(soi.params.get("thfile") or "")).name or "-"],
                ["Iterations", f"{int(soi.params.get('iter', 0))} (strain ratio {_g(soi.params.get('ratio', 0))})"]]
        sections.append({"title": "Site response (SOIL)", "rows": rows})
    equ = decks.get("EQUAKE")
    if equ is not None:
        sp = equ.table("spectra").rows if equ.table("spectra") is not None else []
        rows = [["Target spectra", ", ".join(Path(str(r[1])).name for r in sp) or "-"],
                ["Duration, damping", f"{_g(equ.params.get('dur', 0))} s, {_g(100 * float(equ.params.get('damp', 0)))} %"]]
        sections.append({"title": "Input motion generation (EQUAKE)", "rows": rows})
    return sections


# ======================================================================================
# outputs and charts
# ======================================================================================
def _decimate(t: np.ndarray, y: np.ndarray, n: int = MAX_POINTS) -> Tuple[np.ndarray, np.ndarray]:
    """At most about n points: the minimum and the maximum of every bucket (the peaks stay)."""
    if len(y) <= n:
        return t, y
    k = int(np.ceil(len(y) / (n / 2)))
    idx: List[int] = []
    for s in range(0, len(y), k):
        seg = y[s:s + k]
        a, b = s + int(np.argmin(seg)), s + int(np.argmax(seg))
        idx.extend(sorted({a, b}))
    idx = np.asarray(idx)
    return t[idx], y[idx]


def _label(node: str, kind: str, comp: str) -> str:
    return f"node {int(node)} {'' if kind.upper() == 'TR' else 'R'}{comp.upper()}"


def _outputs(folder: Path, decks: Dict[str, Optional[D.Deck]], motion: Optional[Dict[str, Any]],
             units: str) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    files = sorted(p for p in folder.iterdir() if p.is_file())
    nodal: Dict[str, List[Tuple[Path, re.Match]]] = {}
    for p in files:
        m = NODAL.match(p.name)
        if m:
            nodal.setdefault(m.group(5).upper(), []).append((p, m))
    sections: List[Dict[str, Any]] = []
    charts: List[Dict[str, Any]] = []
    anl = decks.get("ANALYS")
    seismic = anl is None or int(anl.params.get("type", 0) or 0) == 0

    # ---- transfer functions: the SSI resonance
    tfs = []
    for p, m in nodal.get("TFI", []) or nodal.get("TFU", []):
        try:
            f, H, _ = read_tf(p)
        except (OSError, ValueError):
            continue
        if not len(f):
            continue
        a = np.abs(H)
        sel = f > 0
        if not np.any(sel):
            continue
        k = int(np.argmax(np.where(sel, a, -1)))
        tfs.append({"file": p.name, "label": _label(m.group(1), m.group(2), m.group(3)), "f": f, "a": a,
                    "peak": float(a[k]), "fpeak": float(f[k])})
    if tfs:
        tfs.sort(key=lambda r: -r["peak"])
        top = tfs[0]
        rows = [["Largest amplitude", f"{_g(top['peak'])} at {_g(top['fpeak'])} Hz ({top['label']})"]]
        rows += [[r["label"], f"{_g(r['peak'])} at {_g(r['fpeak'])} Hz"] for r in tfs[1:MAX_LINES]]
        sections.append({"title": "Transfer functions" + (" (per unit control motion)" if seismic else " (per unit load)"),
                         "rows": rows})
        show = tfs[:MAX_LINES]
        fmax = max(float(r["f"].max()) for r in show)
        charts.append({"id": "tf", "title": "Transfer-function amplitude |H|", "x": "Frequency (Hz)",
                       "y": "|H|", "logx": False,
                       "range": [0, fmax], "kind": "spec", "files": [r["file"] for r in show],
                       "lines": [{"name": r["label"], "x": r["f"].tolist(), "y": r["a"].tolist()} for r in show]})

    # ---- peak accelerations
    accs = []
    for p, m in nodal.get("ACC", []):
        try:
            a, dt = read_history(p)
        except (OSError, ValueError):
            continue
        if not len(a) or not dt:
            continue
        k = int(np.argmax(np.abs(a)))
        accs.append({"file": p.name, "label": _label(m.group(1), m.group(2), m.group(3)), "a": a, "dt": float(dt),
                     "peak": float(abs(a[k])), "tpeak": k * float(dt)})
    if accs:
        accs.sort(key=lambda r: -r["peak"])
        rows = []
        if motion is not None and seismic:
            rows.append(["Input motion (peak)", f"{_g(motion['pga'])} g"])
        rows += [[r["label"], f"{_g(r['peak'])} g at {_g(r['tpeak'])} s"
                  + (f" ({_g(r['peak'] / motion['pga'], 3)} × input)" if motion and motion["pga"] > 0 and seismic else "")]
                 for r in accs[:MAX_LINES + 2]]
        sections.append({"title": "Peak accelerations", "rows": rows})
        # the largest response, and the input motion drawn over it
        r = accs[0]
        tt, yy = _decimate(np.arange(len(r["a"])) * r["dt"], r["a"])
        lines = [{"name": f"{r['label']} (largest)", "x": tt.tolist(), "y": yy.tolist()}]
        if motion is not None and seismic:                # a load history is not an acceleration
            tt, yy = _decimate(np.arange(len(motion["acc"])) * motion["dt"], motion["acc"])
            lines.append({"name": f"input motion ({motion['file']})", "x": tt.tolist(), "y": yy.tolist(),
                          "color": "#3a4350"})
        charts.append({"id": "acc", "title": "Acceleration: largest response and input", "x": "Time (s)",
                       "y": "Acceleration (g)", "logx": False, "kind": "th", "files": [r["file"]], "lines": lines})

    # ---- in-structure response spectra (the first damping of each node)
    mot = decks.get("MOTION")
    damps = [float(r[0]) for r in (mot.table("damp").rows if mot is not None and mot.table("damp") is not None else [])]
    rss = []
    for p, m in nodal.get("RS", []):
        di = int(m.group(4) or 1)
        try:
            xy = read_xy(p)
        except (OSError, ValueError):
            continue
        if xy.shape[0] < 2 or xy.shape[1] < 2:
            continue
        f, sa = xy[:, 0], xy[:, 1]
        k = int(np.argmax(sa))
        rss.append({"file": p.name, "label": _label(m.group(1), m.group(2), m.group(3)), "di": di, "f": f, "sa": sa,
                    "peak": float(sa[k]), "fpeak": float(f[k]), "zpa": float(sa[int(np.argmax(f))])})
    if rss:
        # the damping of design spectra (5 %) when it was computed, else the first one
        idxs = sorted({r["di"] for r in rss})
        d0 = min(idxs, key=lambda i: abs(damps[i - 1] - 0.05) if 0 < i <= len(damps) else 1.0 + i)
        first = sorted((r for r in rss if r["di"] == d0), key=lambda r: -r["peak"])
        zeta = damps[d0 - 1] if 0 < d0 <= len(damps) else None
        zt = f"{_g(100 * zeta)} % damping" if zeta is not None else f"damping no. {d0}"
        rows = [[r["label"], f"peak {_g(r['peak'])} g at {_g(r['fpeak'])} Hz; ZPA {_g(r['zpa'])} g"]
                for r in first[:MAX_LINES + 2]]
        sections.append({"title": f"In-structure response spectra ({zt})", "rows": rows})
        lines = [{"name": r["label"], "x": r["f"].tolist(), "y": r["sa"].tolist()} for r in first[:MAX_LINES]]
        if motion is not None and zeta is not None and seismic:
            try:
                from ..core.spectra import response_spectrum
                fr = first[0]["f"]
                sa_in = response_spectrum(motion["acc"], motion["dt"], fr, [zeta])["SA"][0]
                lines.append({"name": "input motion", "x": fr.tolist(), "y": sa_in.tolist(), "dash": True,
                              "color": "#3a4350"})
            except Exception:              # noqa: BLE001 -- the input spectrum is an extra
                pass
        charts.append({"id": "rs", "title": f"Response spectra ({zt})", "x": "Frequency (Hz)",
                       "y": "Spectral acceleration (g)", "logx": True, "kind": "spec",
                       "files": [r["file"] for r in first[:MAX_LINES]], "lines": lines})

    # ---- relative displacements and element stresses
    rel = []
    for p, m in nodal.get("THD", []):
        try:
            d, dt = read_history(p)
        except (OSError, ValueError):
            continue
        if len(d) and np.any(d != 0):              # the reference node itself moves 0
            k = int(np.argmax(np.abs(d)))
            rel.append((float(abs(d[k])), k * float(dt or 0), _label(m.group(1), m.group(2), m.group(3))))
    if rel:
        rel.sort(reverse=True)
        sections.append({"title": "Relative displacements (RELDISP)",
                         "rows": [[lab, f"{_g(v)} {units} at {_g(t)} s"] for v, t, lab in rel[:MAX_LINES]]})
    ths = []
    for p in files:
        m = ELEMENT.match(p.name)
        if not m:
            continue
        try:
            s, dt = read_history(p)
        except (OSError, ValueError):
            continue
        if len(s):
            k = int(np.argmax(np.abs(s)))
            ths.append((float(abs(s[k])), f"{m.group(1).upper()} group {int(m.group(2))} element {int(m.group(3))} "
                        f"{m.group(4).upper()}", k * float(dt or 0)))
    if ths:
        ths.sort(reverse=True)
        sections.append({"title": "Element forces and stresses (STRESS, largest)",
                         "rows": [[lab, f"{_g(v)} at {_g(t)} s"] for v, lab, t in ths[:MAX_LINES]]})

    # ---- the free-field modules (SOIL, EQUAKE): spectra and histories under their own names
    if decks.get("SOIL") is not None or decks.get("EQUAKE") is not None:
        rs_ff, th_ff = [], []
        for p in files:
            if NODAL.match(p.name):
                continue
            ext = p.suffix.lower()
            if ext in (".rs", ".rso"):
                try:
                    xy = read_xy(p)
                except (OSError, ValueError):
                    continue
                if xy.shape[0] >= 2 and xy.shape[1] >= 2:
                    rs_ff.append((p.name, xy[:, 0], xy[:, 1]))
            elif (ext == ".th" and p.name.upper().startswith("ACC")) or ext == ".acc":
                try:
                    a, dt = read_history(p)
                except (OSError, ValueError):
                    continue
                if len(a) and dt:
                    th_ff.append((p.name, a, float(dt)))
        def ff_label(n: str) -> str:
            mm = re.match(r"^(RS|ACC)(\d{3})", n, re.IGNORECASE)
            sf = re.match(r"^SAF(\d{3})\w*?_(\d{3})", n, re.IGNORECASE)
            if sf:
                return f"{n} (amplification, layer {int(sf.group(1))} / {int(sf.group(2))})"
            return f"{n} (layer {int(mm.group(2))})" if mm else n

        rows = [[ff_label(n), f"peak {_g(float(sa.max()))}{'' if n.upper().startswith('SAF') else ' g'} at "
                 f"{_g(float(f[int(np.argmax(sa))]))} Hz"] for n, f, sa in rs_ff[:MAX_LINES]]
        rows += [[ff_label(n), f"peak {_g(float(np.max(np.abs(a))))} g"] for n, a, dt in th_ff[:MAX_LINES]]
        rs_ff = [r for r in rs_ff if not r[0].upper().startswith("SAF")] or rs_ff    # g and ratios: not one axis
        if rows:
            sections.append({"title": "Free-field results (SOIL / EQUAKE)", "rows": rows})
        if rs_ff:
            charts.append({"id": "ffrs", "title": "Free-field response spectra", "x": "Frequency (Hz)",
                           "y": "Spectral acceleration (g)", "logx": True, "kind": "spec",
                           "files": [n for n, _, _ in rs_ff[:MAX_LINES]],
                           "lines": [{"name": n, "x": f.tolist(), "y": sa.tolist()} for n, f, sa in rs_ff[:MAX_LINES]]})
        if th_ff:
            lines = []
            for n, a, dt in th_ff[:3]:
                tt, yy = _decimate(np.arange(len(a)) * dt, a)
                lines.append({"name": n, "x": tt.tolist(), "y": yy.tolist()})
            charts.append({"id": "ffacc", "title": "Free-field acceleration histories", "x": "Time (s)",
                           "y": "Acceleration (g)", "logx": False, "kind": "th",
                           "files": [n for n, _, _ in th_ff[:3]], "lines": lines})
    return sections, charts


def _profile_chart(sit: Optional[D.Deck], units: str) -> Optional[Dict[str, Any]]:
    """The soil profile: Vs against depth (a step line), the half-space below."""
    if sit is None or sit.table("layers") is None or not sit.table("layers").rows:
        return None
    x: List[float] = []
    y: List[float] = []
    z = 0.0
    for r in sit.table("layers").rows:
        vs, h = float(r[4]), float(r[1])
        x += [vs, vs]
        y += [z, z + h]
        z += h
    hs = sit.table("halfspace").rows if sit.table("halfspace") is not None else []
    lines = [{"name": "layers", "x": x, "y": y}]
    if hs:
        vs = float(hs[0][4])
        lines.append({"name": "half-space", "x": [x[-1], vs, vs], "y": [z, z, z * 1.25 if z > 0 else 1.0],
                      "dash": True, "color": "#8a6d3b"})
    vmax = max(max(ln["x"]) for ln in lines)
    return {"id": "soil", "title": "Soil profile: shear-wave velocity", "x": f"Vs ({units}/s)",
            "y": f"Depth ({units})", "logx": False, "reversey": True, "range": [0, 1.1 * vmax], "kind": None,
            "files": [], "lines": lines}


# ======================================================================================
# entry points
# ======================================================================================
def summary(interp, model: Optional[str] = None) -> Dict[str, Any]:
    """Everything the summary window shows for one model (module docstring)."""
    number, m = find_model(interp, model)
    folder = _folder(m)
    others = [{"number": n, "name": mm.name, "title": mm.title} for n, mm in sorted(interp.models.items())
              if _folder(mm) is not None and module_runs(_folder(mm), mm.name)]
    head = {"number": number, "name": m.name, "title": m.title, "path": str(folder) if folder else ""}
    if folder is None or not m.name:
        return {"model": head, "modules": [], "inputs": [], "outputs": [], "charts": [], "others": others,
                "note": "this model has no folder yet (MDL): nothing has run"}
    runs = module_runs(folder, m.name)
    decks = {mod: _deck(folder, m.name, mod) for mod in D.SCHEMAS}
    gravity = 0.0
    for mod in ("HOUSE", "SITE", "MOTION", "ANALYS"):
        if decks.get(mod) is not None and decks[mod].params.get("gravity"):
            gravity = float(decks[mod].params["gravity"])
            break
    units = _units(gravity) if gravity else "length units"
    motion = input_motion(folder, decks.get("MOTION") or decks.get("STRESS") or decks.get("RELDISP"))
    inputs = _inputs(m, folder, decks, motion, units)
    outputs, charts = _outputs(folder, decks, motion, units)
    prof = _profile_chart(decks.get("SITE"), units)
    if prof is not None:
        charts.append(prof)
    if motion is not None:
        motion = {k: v for k, v in motion.items() if k != "acc"}
    return {"model": head, "modules": runs, "inputs": inputs, "outputs": outputs, "charts": charts,
            "others": others, "units": units, "motion": motion,
            "note": "" if runs else "no module has run in this model's folder yet"}


def plot_lines(interp, model: Optional[str], chart: str, helpdocs=None, lines_in_memory: Sequence[int] = ()
               ) -> List[str]:
    """The command text that opens a chart of the summary as a plot (READSPEC / READTH and SPECPLOT / THPLOT)."""
    s = summary(interp, model)
    ch = next((c for c in s["charts"] if c["id"] == chart), None)
    if ch is None or not ch.get("files") or ch.get("kind") not in ("spec", "th"):
        raise SummaryError(f"chart {chart!r} has no result files to plot", 400)
    from . import learn
    verb = "plot-spectrum" if ch["kind"] == "spec" else "plot-history"
    opt = " | log" if ch.get("logx") else ""
    r = learn.resolve_action(verb, ", ".join(ch["files"]) + opt, Path(s["model"]["path"]), interp,
                             list(lines_in_memory), helpdocs)
    return list(r["lines"])
