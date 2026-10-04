"""NONLINEAR module (Option NON): equivalent-linear iterations of nonlinear wall panels and springs, and the
auxiliary program COMB_XYZ_THD (requirements 4.15, 2.4; spec 05d section 3; decisions D-NON-01 ... D-NON-12).

What NONLINEAR does (for the structural engineer)
-------------------------------------------------
The SSI solution of SASSI is linear in the frequency domain.  Option NON treats nonlinear RC wall
panels (groups of coplanar SHELL elements) and nonlinear translational springs by **iterative
equivalent linearisation** (manual Fig. 1.2)::

    elastic SSI:  SITE, POINT, HOUSE, ANALYS (save restart files), MOTION, RELDISP (X, Y, Z input)
                  COMB_XYZ_THD -> NONLINEAR  (no PANEL.NON: elastic run, writes PANEL.NON = 1)
    iteration k:  <model>_new.hou -> <model>.hou, HOUSE, ANALYS ("New Structure" restart), MOTION, RELDISP
                  COMB_XYZ_THD -> NONLINEAR  (PANEL.NON = 1: iteration) ... until converged

For every nonlinear element NONLINEAR reads the relative-displacement histories ``.THD`` of its nodes
(panel corners: ``nnnnnTR_X/Y/Z.THD``; spring ends: the spring DOF), forms the deformation history x(t)
(panel shear strain, :func:`sassi.core.panels.panel_deformations`; spring elongation u_J - u_I), runs
the hysteresis model of the element (:mod:`sassi.core.hysteresis`: CMS, TAK or GMR) to get the force
history F(t), and linearises it at ``x_eq = EDF max|x|``:

* ``E_new = E_el K_sec(x_eq)/K_el`` (panels; Poisson's ratio constant, so shear, bending and axial
  stiffness degrade together, manual section 1.5.4) or ``k_new = k_el K_sec/K_el`` (springs), with
  ``K_el`` the BBC first slope Y1/X1 (D-NON-11);
* ``xi = min(cutoff, scale xi_h + [xi_el])`` with the loop-area damping ``xi_h(x_eq)`` (D-NON-04);
* ductility ``mu = max|x|/x_cr`` (x_cr = BBC point 1) and force-reduction factor
  ``F_mu = K_el max|x| / |F(t*)|`` at the time t* of max|x|.

Inputs (model directory): ``<model>.eql`` (this module's deck, D-NON-01: EQL header, BBC points, panel
and spring records resolved from the model by RUNNONLINEAR), ``<model>.hou`` (the HOUSE deck of the
analysis just made), the ``.THD`` files (RELDISP, combined by COMB_XYZ_THD for X/Y/Z input) and the state
files below.

State files (manual section 6.5.4): ``PANEL.NON`` / ``SPRING.NON`` absent -> elastic run (the SSI
results are those of the elastic model), the file is then created with ``1``; present with ``1`` ->
nonlinear iteration: the properties of the analysis just made are read from ``Panel_EQL_Matl_Prop.txt``
/ ``SPRING_EQL_Matl_Prop.txt`` ("used" for the convergence measure and checked against ``<model>.hou``).
Deleting both files (NONLINRESET) restarts from the elastic analysis.

Outputs: ``PanelNNNN.thd`` (strain history), ``PanelNNNN_AXIAL.thd`` (vertical axial strain),
``PanelNNNN.ths`` (nonlinear force history), ``PanelNNNN.crv`` (equivalent-linear E/E_el and damping
curves versus amplitude, written by the elastic run), ``Panel_EQL_Matl_Prop.txt`` (new properties),
``Panel.fmu`` (ductility and F_mu), the same with ``SPRING`` for springs, ``<model>_new.hou`` (the
HOUSE deck for the next iteration), ``NONLINEAR_CONVERGENCE.TXT`` (one row per run: max |dE/E| and max
|d xi|; converged when < 2 % and < 0.5 %, D-NON-06) and the listing ``<model>_NONLINEAR.out``.

COMB_XYZ_THD (D-NON-08): ``COMB_XYZ_THD.inp`` holds ``n`` on line 1, then n lines ``outfile fileX fileY
fileZ`` (``-`` for a missing direction); every output is the algebraic sum of the three histories,
time step by time step (:func:`comb_xyz_thd`).

    python -m sassi.modules.nonlinear < NONLINEAR.inp          # batch use (model, deck, listing)
    python -m sassi.modules.nonlinear COMB_XYZ_THD [file.inp]   # COMB_XYZ_THD in the current directory
"""
from __future__ import annotations

import math
import sys
import time
import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, TextIO, Tuple, Union

import numpy as np

from ..conventions import nodal_result_name
from ..core import hysteresis as HY
from ..core.panels import CORNERS, panel_deformations
from ..io import decks
from ..io.deckfmt import Table, read_raw, write_raw
from .base import Listing, ModuleContext, ModuleError

NAME = "NONLINEAR"
DECK_EXT = ".eql"
NON_FILES = {"panel": "PANEL.NON", "spring": "SPRING.NON"}
PROP_FILES = {"panel": "Panel_EQL_Matl_Prop.txt", "spring": "SPRING_EQL_Matl_Prop.txt"}
FMU_FILES = {"panel": "Panel.fmu", "spring": "Spring.fmu"}
PREFIX = {"panel": "Panel", "spring": "SPRING"}
CONV_FILE = "NONLINEAR_CONVERGENCE.TXT"
COMB_INP = "COMB_XYZ_THD.inp"
#: EQL <NonLinOpts> bits (D-NON-02)
OPT_PANELS, OPT_SPRINGS, OPT_BEAMS = 1, 2, 4

PathLike = Union[str, Path]


# ======================================================================================
# The .eql deck (D-NON-01): own schema on the keyword deck format of sassi.io.deckfmt
# ======================================================================================
EQL_PARAMS: Tuple[Tuple[str, type, Any, str], ...] = (
    ("title", str, "", "model title"),
    ("model", str, "", "model name"),
    ("edf", float, 0.8, "EQL <disp>: equivalent-linear displacement factor EDF"),
    ("nonlinopts", int, 0, "EQL <NonLinOpts> resolved: 1 panels, 2 springs (4 beams: not available)"),
    ("dampcutoff", float, 0.0, "EQL <dampCutoff>: damping cut-off in percent (0 = none)"),
    ("dampscale", float, 0.0, "EQL <dampScale>: damping scale factor (0 = 1)"),
    ("elasticd", int, 0, "EQL <ElasicD>: 1 = add the elastic damping to the hysteretic damping"),
    ("gravity", float, 32.2, "HOUSE gravity"),
    ("nonext", int, 0, "EDUOPT,NONEXT: 1 unlocks TAK and bending panels (experimental, D-NON-05)"),
    ("maxit", int, HY.MAX_ITERATIONS, "maximum number of iterations (D-NON-06)"),
    ("tol_e", float, HY.TOL_E, "convergence: max relative change of E"),
    ("tol_xi", float, HY.TOL_XI, "convergence: max absolute change of the damping ratio"),
    ("maxnode", int, 0, "largest node number of the model (width of the .THD file names)"),
    ("crvpoints", int, 60, "number of amplitudes of the .crv curves"),
)
EQL_TABLES: Dict[str, Tuple[Tuple[str, type], ...]] = {
    "bbc": tuple((c, t) for c, t in (("bbc", int), ("type", int), ("yield", int), ("point", int), ("x", float),
                                       ("y", float))),
    "panels": tuple([("num", int), ("group", int), ("bbc", int), ("disp", int), ("force", int), ("mat", int),
                     ("mtype", int), ("val1", float), ("val2", float), ("weight", float), ("pdamp", float),
                     ("sdamp", float), ("e_el", float), ("nu_el", float), ("g_el", float), ("thick", float),
                     ("length", float), ("height", float), ("area", float), ("inertia", float), ("bl", int),
                     ("br", int), ("tr", int), ("tl", int), ("ex", float), ("ey", float), ("nelem", int)]),
    "springs": tuple([("num", int), ("group", int), ("elem", int), ("bbc", int), ("disp", int), ("force", int),
                      ("prop", int), ("node_i", int), ("node_j", int), ("scx", float), ("scy", float), ("scz", float),
                      ("scxx", float), ("scyy", float), ("sczz", float), ("damp", float), ("k_el", float)]),
}


@dataclass
class EqlDeck:
    """Typed content of ``<model>.eql``."""
    params: Dict[str, Any] = field(default_factory=lambda: {p[0]: p[2] for p in EQL_PARAMS})
    tables: Dict[str, List[Dict[str, Any]]] = field(default_factory=lambda: {k: [] for k in EQL_TABLES})

    def __getitem__(self, k: str) -> Any:
        return self.params[k]

    def rows(self, name: str) -> List[Dict[str, Any]]:
        return self.tables[name]

    def add(self, name: str, **row) -> None:
        cols = [c for c, _ in EQL_TABLES[name]]
        missing = [c for c in cols if c not in row]
        extra = [c for c in row if c not in cols]
        if missing or extra:
            raise KeyError(f".eql table {name}: missing {missing}, unknown {extra}")
        self.tables[name].append({c: row[c] for c in cols})

    def backbones(self) -> Dict[int, Tuple[int, int, np.ndarray, np.ndarray]]:
        """``{bbc: (type, yield, x, y)}`` from the bbc table."""
        out: Dict[int, List] = {}
        for r in sorted(self.tables["bbc"], key=lambda r: (r["bbc"], r["point"])):
            b = out.setdefault(int(r["bbc"]), [int(r["type"]), int(r["yield"]), [], []])
            b[2].append(float(r["x"]))
            b[3].append(float(r["y"]))
        return {k: (v[0], v[1], np.asarray(v[2]), np.asarray(v[3])) for k, v in out.items()}


def _coerce(v: Any, typ: type) -> Any:
    if typ is str:
        return "" if v is None else str(v)
    t = str(v).strip().replace("D", "E").replace("d", "e") if isinstance(v, str) else v
    if typ is int:
        f = float(t) if t != "" else 0.0
        if f != int(round(f)):
            raise ValueError(f"expected an integer, got {v!r}")
        return int(round(f))
    return float(t) if t != "" else 0.0


def write_eql(path: PathLike, deck: EqlDeck, comments: Sequence[str] = ()) -> Path:
    """Write ``<model>.eql`` (keyword deck, D-NON-01)."""
    params = {name: _coerce(deck.params.get(name, dflt), typ) for name, typ, dflt, _ in EQL_PARAMS}
    tables: Dict[str, Table] = {}
    for name, cols in EQL_TABLES.items():
        tab = Table(columns=[c for c, _ in cols])
        for r in deck.tables.get(name, []):
            tab.rows.append([_coerce(r[c], t) for c, t in cols])
        tables[name] = tab
    return write_raw(path, NAME, params, tables, comments=comments)


def read_eql(path: PathLike) -> EqlDeck:
    """Read ``<model>.eql``; missing parameters take their defaults."""
    raw = read_raw(path)
    if raw.module != NAME:
        raise ModuleError(f"{Path(path).name}: expected a NONLINEAR deck, found {raw.module}")
    d = EqlDeck()
    for name, typ, dflt, _ in EQL_PARAMS:
        if name in raw.params:
            v = raw.params[name]
            if typ is str and v.startswith('"'):
                import json
                v = json.loads(v)
            d.params[name] = _coerce(v, typ)
    for name, cols in EQL_TABLES.items():
        rt = raw.tables.get(name)
        if rt is None:
            continue
        want = [c for c, _ in cols]
        if rt.columns != want:
            raise ModuleError(f"{Path(path).name}: table {name} has columns {rt.columns}, expected {want}")
        for row in rt.rows:
            d.tables[name].append({c: _coerce(v, t) for (c, t), v in zip(cols, row)})
    return d


# ======================================================================================
# Small text files
# ======================================================================================
def write_history(path: Path, values: np.ndarray, dt: float, header: str = "") -> Path:
    """``.thd`` / ``.ths``: optional '#' header, first line dt, then one value per line (D-FIL-02)."""
    lines = [f"# {header}"] if header else []
    lines.append(f"{float(dt):.16g}")
    lines += [f"{v:.16e}" for v in np.asarray(values, dtype=float)]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def read_history(path: Path) -> Tuple[np.ndarray, float]:
    """``(values, dt)`` of a history file (first number dt; '#' and '*' lines skipped)."""
    vals: List[float] = []
    for ln in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        s = ln.strip()
        if not s or s[0] in "#*":
            continue
        for t in s.replace(",", " ").split():
            try:
                vals.append(float(t.replace("D", "E")))
            except ValueError:
                raise ModuleError(f"{Path(path).name}: '{t}' is not a number") from None
    if not vals:
        raise ModuleError(f"{Path(path).name} is empty")
    return np.asarray(vals[1:], dtype=float), float(vals[0])


def find_file(workdir: Path, name: str) -> Optional[Path]:
    """``workdir/name`` or a case variant of it (the manual's names differ in case between chapters)."""
    p = workdir / name
    if p.exists():
        return p
    low = name.lower()
    for q in workdir.iterdir():
        if q.name.lower() == low:
            return q
    return None


def read_state(workdir: Path, kind: str) -> Optional[int]:
    """Content of PANEL.NON / SPRING.NON (None when absent)."""
    p = find_file(workdir, NON_FILES[kind])
    if p is None:
        return None
    txt = p.read_text(encoding="utf-8", errors="replace").split()
    try:
        return int(float(txt[0])) if txt else 0
    except ValueError:
        return 0


@dataclass
class PropRow:
    num: int
    group: int
    ident: int            # material (panels) or SC property (springs)
    k_el: float           # elastic E (panels) or spring constant (springs)
    k_new: float
    ratio: float
    xi_el: float
    xi_h: float
    xi_new: float
    x_max: float
    x_eq: float


PROP_COLUMNS = ("element", "group", "mat_or_prop", "elastic", "new", "ratio", "xi_elastic", "xi_hyst", "xi_new",
                "x_max", "x_eq")


def write_props(path: Path, kind: str, model: str, iteration: int, rows: Sequence[PropRow], note: str = "") -> Path:
    what = "wall-panel material (E)" if kind == "panel" else "spring constant (k)"
    lines = [f"# SASSI-EDU NONLINEAR equivalent-linear {what} properties (Option NON), model {model}",
             f"# properties of SSI iteration {iteration} (written after SSI analysis {iteration - 1}"
             f"{', the elastic analysis' if iteration == 1 else ''}); damping as ratios",
             f"ITERATION {iteration}"]
    if note:
        lines.append(f"# {note}")
    lines.append("# " + " ".join(f"{c:>23s}" for c in PROP_COLUMNS))
    for r in rows:
        lines.append("  " + " ".join([f"{r.num:>23d}", f"{r.group:>23d}", f"{r.ident:>23d}"]
                                    + [f"{v:>23.16e}" for v in (r.k_el, r.k_new, r.ratio, r.xi_el, r.xi_h, r.xi_new,
                                                               r.x_max, r.x_eq)]))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def read_props(path: Path) -> Tuple[int, Dict[int, PropRow]]:
    """``(iteration, {element: row})`` of a ``*_EQL_Matl_Prop.txt`` file."""
    it = 0
    rows: Dict[int, PropRow] = {}
    for ln in path.read_text(encoding="utf-8", errors="replace").splitlines():
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        if s.upper().startswith("ITERATION"):
            it = int(s.split()[1])
            continue
        t = s.split()
        if len(t) < 11:
            continue
        r = PropRow(int(t[0]), int(t[1]), int(t[2]), *[float(v) for v in t[3:11]])
        rows[r.num] = r
    return it, rows


def read_convergence(workdir: Path) -> List[Dict[str, float]]:
    """Rows of NONLINEAR_CONVERGENCE.TXT: iteration, elements, max dE/E %, max dxi %, converged."""
    p = find_file(Path(workdir), CONV_FILE)
    out: List[Dict[str, float]] = []
    if p is None:
        return out
    for ln in p.read_text(encoding="utf-8", errors="replace").splitlines():
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        t = s.split()
        if len(t) < 5:
            continue
        try:
            out.append(dict(iteration=int(t[0]), elements=int(t[1]), max_de_pct=float(t[2]),
                            max_dxi_pct=float(t[3]), converged=int(t[4])))
        except ValueError:
            continue
    return out


def _append_convergence(workdir: Path, iteration: int, n: int, conv: HY.Convergence) -> None:
    p = Path(workdir) / CONV_FILE
    new = not p.exists()
    with open(p, "a", encoding="utf-8") as fh:
        if new:
            fh.write("# SASSI-EDU NONLINEAR convergence history (D-NON-06: max |dE/E| < 2 %, max |d xi| < 0.5 %)\n")
            fh.write("# analysis = SSI analysis processed (0 = elastic); the properties of analysis+1 are compared "
                     "with those of analysis\n")
            fh.write(f"# {'analysis':>8s} {'elements':>9s} {'max_dE/E_%':>12s} {'max_dxi_%':>12s} converged\n")
        fh.write(f"  {iteration:>8d} {n:>9d} {100 * conv.max_de:>12.5f} {100 * conv.max_dxi:>12.5f} "
                 f"{1 if conv.converged else 0:>9d}\n")


# ======================================================================================
# The nonlinear elements of one run
# ======================================================================================
@dataclass
class NLElement:
    kind: str               # 'panel' or 'spring'
    num: int
    group: int
    ident: int              # material (panel) or SC property (spring)
    code: int               # hysteresis model (Force Opt)
    disp: int               # panel: 1 shear, 2 bending; spring: DOF 1..3
    bbc: int
    backbone: HY.Backbone
    k_el: float             # elastic E (panel) or spring constant (spring)
    xi_el: float            # elastic damping
    row: Dict[str, Any]     # the deck row
    x: Optional[np.ndarray] = None
    axial: Optional[np.ndarray] = None
    result: Optional[HY.EquivalentLinear] = None
    k_cur: float = 0.0      # property of the analysis just made
    xi_cur: float = 0.0

    @property
    def label(self) -> str:
        return f"{PREFIX[self.kind]}{self.num:04d}"


def _node_history(workdir: Path, node: int, dof: int, maxnode: int, cache: Dict[Tuple[int, int], Tuple[np.ndarray, float]],
                  ) -> Tuple[np.ndarray, float]:
    key = (node, dof)
    if key in cache:
        return cache[key]
    names = dict.fromkeys([nodal_result_name(node, dof, "THD", maxnode), nodal_result_name(node, dof, "THD", 0),
                           nodal_result_name(node, dof, "THD", 100000)])
    for nm in names:
        p = find_file(workdir, nm)
        if p is not None:
            cache[key] = read_history(p)
            return cache[key]
    raise ModuleError(f"{nodal_result_name(node, dof, 'THD', maxnode)} not found: node {node} must be in the MOTION "
                      "and RELDISP output requests (NONLINMOTDISP before AFWRITE) and RELDISP (and COMB_XYZ_THD for "
                      "X/Y/Z input) must run before NONLINEAR")


def _elastic_e(row: Dict[str, Any]) -> float:
    return float(row["e_el"])


def scale_material(mtype: int, val1: float, val2: float, ratio: float) -> Tuple[float, float]:
    """M-table values of a material whose moduli are multiplied by ``ratio`` (Poisson's ratio constant)."""
    if mtype == 1:
        return val1 * ratio, val2
    if mtype == 2:
        return val1 * ratio, val2 * ratio
    if mtype == 3:
        r = math.sqrt(ratio)
        return val1 * r, val2 * r
    raise ModuleError(f"material type {mtype} is not 1, 2 or 3")


# ======================================================================================
# Module run
# ======================================================================================
def _refuse_renumbered_model(ctx: ModuleContext) -> None:
    """Option NON works in the model numbering: the panel corners and spring ends of the ``.eql`` deck and
    the ``<model>_new.hou`` it writes are model node numbers, and the ``.THD`` files are read by those
    numbers.  When the HOUSE node optimizer changed the numbering (a non-identity ``<model>.map``; HOUSE
    renames a stale map ``.prev``), FILE8 and therefore the MOTION/RELDISP files use the new numbers, so a
    ``.THD`` file named after a model number holds another node's history: stop instead of computing panel
    strains from the wrong nodes (final audit, docs/verification/AUDIT_REPORT.md F-04)."""
    p = ctx.path(f"{ctx.model}.map")
    if not p.exists():
        return
    try:
        from ..core.renumber import read_map
        moved = sum(1 for o, n in read_map(p).items() if o != n)
    except (OSError, ValueError) as exc:
        raise ModuleError(f"{p.name} (HOUSE node optimizer map) unreadable: {exc}") from None
    if moved:
        raise ModuleError(f"the HOUSE node optimizer renumbered {moved} nodes ({p.name}): FILE8 and the MOTION/RELDISP "
                          "files use the new numbers, while Option NON reads the .THD histories of its panel and "
                          "spring nodes by their model numbers -- run HOUSE without the optimizer (HOUSEX "
                          "<optimize> = 0) for Option NON")


def run(ctx: ModuleContext) -> int:
    """NONLINEAR module entry point (requirements 4.15)."""
    lst = ctx.listing
    wd = ctx.workdir
    deck = read_eql(ctx.deck_path)
    if deck["title"]:
        lst.write(f" {deck['title']}")
    hou_path = ctx.require(f"{ctx.model}.hou", "AFWRITE (HOUSE deck)")
    with _quiet_warnings():
        hou = decks.read(hou_path, "HOUSE")
    _refuse_renumbered_model(ctx)
    opts = int(deck["nonlinopts"])
    kinds = [k for k, bit in (("panel", OPT_PANELS), ("spring", OPT_SPRINGS)) if opts & bit]
    edf = float(deck["edf"])
    cutoff, scale, elasticd = float(deck["dampcutoff"]), float(deck["dampscale"]), bool(int(deck["elasticd"]))
    maxit = int(deck["maxit"])
    tol_e, tol_xi = float(deck["tol_e"]), float(deck["tol_xi"])

    lst.section("NONLINEAR options (EQL)")
    lst.write(f"   equivalent-linear displacement factor EDF : {edf:g}")
    lst.write(f"   nonlinear elements (NonLinOpts {opts})       : "
              + (", ".join("wall panels" if k == "panel" else "springs" for k in kinds) or "none"))
    lst.write(f"   damping cut-off                          : {'none' if cutoff <= 0 else f'{cutoff:g} %'}")
    lst.write(f"   damping scale factor                     : {scale if scale != 0 else 1.0:g}"
              + ("  (0 given: 1.0, D-NON-04)" if scale == 0 else ""))
    lst.write(f"   include elastic damping                  : {'yes' if elasticd else 'no'}")
    lst.write(f"   convergence (D-NON-06)                   : |dE/E| < {100 * tol_e:g} %, |d xi| < {100 * tol_xi:g} %"
              f", at most {maxit} iterations")
    if opts & OPT_BEAMS:
        lst.warning("nonlinear beams (NonLinOpts 4) are not available in this version: ignored")
    if not kinds:
        raise ModuleError("no nonlinear element type selected (EQL <NonLinOpts>: 1 panels, 2 springs)")

    # ---------------------------------------------------------------- backbone curves
    bbcs = deck.backbones()
    lst.section("Backbone curves (BBC; origin implied, point 1 = cracking point)")
    lst.write(f"{'BBC':>6s} {'type':>6s} {'points':>7s} {'yield':>6s} {'X1':>13s} {'Y1':>13s} {'K_el=Y1/X1':>13s} "
              f"{'X_yield':>13s} {'Y_yield':>13s} {'X_last':>13s}")
    backbones: Dict[int, HY.Backbone] = {}
    for b, (typ, yi, x, y) in sorted(bbcs.items()):
        try:
            bb = HY.Backbone(x, y, yi)
        except HY.HysteresisError as exc:
            raise ModuleError(f"BBC {b}: {exc}") from None
        backbones[b] = bb
        lst.write(f"{b:>6d} {HY.MODEL_CODES.get(typ, str(typ)):>6s} {bb.n:>7d} {yi:>6d} {bb.x_cr:>13.5e} {bb.f_cr:>13.5e} "
                  f"{bb.k_el:>13.5e} {bb.x_y:>13.5e} {bb.f_y:>13.5e} {bb.x[-1]:>13.5e}")

    # ---------------------------------------------------------------- elements
    elements: List[NLElement] = []
    if "panel" in kinds:
        for r in deck.rows("panels"):
            bb = backbones.get(int(r["bbc"]))
            if bb is None:
                raise ModuleError(f"panel {r['num']}: BBC {r['bbc']} is not defined")
            elements.append(NLElement("panel", int(r["num"]), int(r["group"]), int(r["mat"]), int(r["force"]),
                                      int(r["disp"]), int(r["bbc"]), bb, _elastic_e(r), float(r["sdamp"]), r))
    if "spring" in kinds:
        for r in deck.rows("springs"):
            bb = backbones.get(int(r["bbc"]))
            if bb is None:
                raise ModuleError(f"spring {r['num']}: BBC {r['bbc']} is not defined")
            elements.append(NLElement("spring", int(r["num"]), int(r["group"]), int(r["prop"]), int(r["force"]),
                                      int(r["disp"]), int(r["bbc"]), bb, float(r["k_el"]), float(r["damp"]), r))
    if not elements:
        raise ModuleError("no nonlinear panels (P) or springs (S) in the deck")
    for el in elements:
        if el.code == 2:
            raise ModuleError(f"{el.kind} {el.num}: the Cheng-Mertz Bending model (Force Opt 2) is not included in this "
                              "version")
        if el.kind == "panel" and (el.code != 1 or el.disp != 1) and not int(deck["nonext"]):
            raise ModuleError(f"panel {el.num}: Force Opt {el.code} / Disp. Opt {el.disp} is not supported (Errors "
                              "126/128: set 1; EDUOPT,NONEXT,1 unlocks TAK and bending panels as experimental)")
        if el.kind == "spring" and el.code != 4 and not int(deck["nonext"]):
            raise ModuleError(f"spring {el.num}: Force Opt {el.code} is not supported (Error 127: set 4, GMR)")
    _list_elements(lst, elements)

    # ---------------------------------------------------------------- state machine
    states = {k: read_state(wd, k) for k in kinds}
    present = [k for k in kinds if states[k] is not None]
    if present and len(present) != len(kinds):
        raise ModuleError("inconsistent state files: " + ", ".join(f"{NON_FILES[k]} {'present' if states[k] is not None else 'absent'}"
                                                                    for k in kinds)
                          + " -- delete them (NONLINRESET) to restart from the elastic analysis")
    iterating = bool(present) and all(states[k] == 1 for k in kinds)
    if present and not iterating:
        raise ModuleError("a .NON state file does not contain 1: delete it (NONLINRESET) to restart from the elastic "
                          "analysis")
    lst.section("Run state (PANEL.NON / SPRING.NON)")
    analysis = 0
    current: Dict[Tuple[str, int], PropRow] = {}
    if iterating:
        its = []
        for k in kinds:
            p = find_file(wd, PROP_FILES[k])
            if p is None:
                lst.warning(f"{NON_FILES[k]} = 1 but {PROP_FILES[k]} is missing: the properties of the HOUSE deck "
                            f"{hou_path.name} are taken as those of the analysis just made")
                continue
            it, rows = read_props(p)
            its.append(it)
            for num, r in rows.items():
                current[(k, num)] = r
            lst.write(f"   {NON_FILES[k]} = 1: nonlinear iteration; {p.name} holds the properties of SSI iteration {it} "
                      "(used as the properties of the analysis just made)")
        analysis = max(its) if its else 1
    else:
        lst.write("   " + " and ".join(NON_FILES[k] for k in kinds) + " absent: ELASTIC run -- the SSI results are those of "
                  "the elastic model; the .NON file(s) will be created with 1")
    lst.write(f"   SSI analysis processed: {analysis} ({'elastic' if analysis == 0 else 'iteration ' + str(analysis)})")
    _current_properties(elements, current, hou, lst, iterating)

    # ---------------------------------------------------------------- deformation histories
    cache: Dict[Tuple[int, int], Tuple[np.ndarray, float]] = {}
    dt_ref: Optional[float] = None
    nt_ref: Optional[int] = None
    maxnode = int(deck["maxnode"])

    def hist(node: int, dof: int) -> np.ndarray:
        nonlocal dt_ref, nt_ref
        v, dt = _node_history(wd, node, dof, maxnode, cache)
        if dt_ref is None:
            dt_ref, nt_ref = dt, len(v)
        elif abs(dt - dt_ref) > 1e-9 * max(dt_ref, 1e-30) or len(v) != nt_ref:
            raise ModuleError(f"node {node}: its .THD history (dt {dt:g}, {len(v)} values) differs from the others "
                              f"(dt {dt_ref:g}, {nt_ref} values): run RELDISP for all nodes with the same options")
        return v

    for el in elements:
        if ctx.cancelled():
            raise ModuleError("run cancelled")
        r = el.row
        if el.kind == "panel":
            disp = {}
            for lab, key in zip(CORNERS, ("bl", "br", "tr", "tl")):
                n = int(r[key])
                disp[lab] = (hist(n, 1), hist(n, 2), hist(n, 3))
            g, k, e = panel_deformations(disp, (float(r["ex"]), float(r["ey"]), 0.0), float(r["length"]),
                                         float(r["height"]))
            el.x = g if el.disp == 1 else k
            el.axial = e
        else:
            dof = el.disp
            el.x = hist(int(r["node_j"]), dof) - hist(int(r["node_i"]), dof)
    dt = float(dt_ref or 0.0)

    # ---------------------------------------------------------------- equivalent linearisation
    lst.section("Equivalent-linear properties for the next SSI iteration")
    lst.write(f"{'element':>12s} {'model':>5s} {'max|x|':>12s} {'x_eq':>12s} {'x_cr':>12s} {'K_sec/K_el':>11s} "
              f"{'xi_h':>8s} {'xi_new':>8s} {'mu':>8s} {'F_mu':>8s}")
    for i, el in enumerate(elements):
        res = HY.equivalent_linear(el.code, el.backbone, el.x, edf, el.xi_el, cutoff, scale, elasticd)
        el.result = res
        lst.write(f"{el.label:>12s} {HY.MODEL_CODES[el.code]:>5s} {res.x_max:>12.5e} {res.x_eq:>12.5e} "
                  f"{el.backbone.x_cr:>12.5e} {res.ratio:>11.5f} {res.xi_h:>8.4f} {res.xi:>8.4f} {res.ductility:>8.3f} "
                  f"{res.f_mu:>8.3f}")
        for n in res.notes:
            lst.write(f"      note ({el.label}): {n}")
        if res.beyond_bbc:
            lst.warning(f"{el.label}: the deformation passed the last BBC point (x = {el.backbone.x[-1]:g}): failure range "
                        "of the backbone, the force is held at its last value")
        if res.xi < el.xi_el and int(el.code) != 3:
            lst.write(f"      note ({el.label}): damping {res.xi:.4f} below the elastic damping {el.xi_el:.4f} (manual: "
                      "hysteretic-only damping should not be smaller; consider Include Elastic Damping)")
        ctx.progress((i + 1) / len(elements), f"NONLINEAR {i + 1}/{len(elements)}")

    # ---------------------------------------------------------------- convergence
    e_old = [el.k_cur for el in elements]
    e_new = [el.k_el * el.result.ratio for el in elements]
    x_old = [el.xi_cur for el in elements]
    x_new = [el.result.xi for el in elements]
    conv = HY.convergence(e_old, e_new, x_old, x_new, tol_e, tol_xi)
    lst.section("Convergence (D-NON-06)")
    lst.write(f"{'element':>12s} {'E or k used':>14s} {'E or k new':>14s} {'dE/E %':>9s} {'xi used':>9s} {'xi new':>9s} "
              f"{'d xi %':>8s}")
    for el, a, b, c, d in zip(elements, e_old, e_new, x_old, x_new):
        lst.write(f"{el.label:>12s} {a:>14.6e} {b:>14.6e} {100 * (b - a) / a:>9.3f} {c:>9.5f} {d:>9.5f} {100 * (d - c):>8.3f}")
    lst.write(f"   SSI analysis {analysis}: " + conv.text())
    if conv.converged:
        lst.write(f"   CONVERGED: the properties used in SSI analysis {analysis} are consistent with its response")
    elif analysis >= maxit:
        lst.warning(f"not converged after {analysis} iterations (limit {maxit}, D-NON-06)")

    # ---------------------------------------------------------------- outputs
    written = _write_outputs(ctx, deck, hou, elements, dt, analysis, kinds, iterating, conv)
    _append_convergence(wd, analysis, len(elements), conv)
    lst.section("Files written")
    for name in written:
        lst.write(f"   {name}")
    return 0


class _quiet_warnings:
    def __enter__(self):
        import warnings
        self._cm = warnings.catch_warnings()
        self._cm.__enter__()
        warnings.simplefilter("ignore")
        return self

    def __exit__(self, *exc):
        return self._cm.__exit__(*exc)


def _list_elements(lst: Listing, elements: Sequence[NLElement]) -> None:
    panels = [e for e in elements if e.kind == "panel"]
    springs = [e for e in elements if e.kind == "spring"]
    if panels:
        lst.section("Wall panels (P; corners BL BR TR TL, local axis e_h, width L, height H)")
        lst.write(f"{'panel':>6s} {'group':>6s} {'mat':>5s} {'BBC':>4s} {'model':>5s} {'disp':>5s} {'corners':>27s} "
                  f"{'L':>10s} {'H':>10s} {'t':>8s} {'E_el':>12s} {'G A':>12s} {'K_el(BBC)':>12s}")
        for e in panels:
            r = e.row
            corners = " ".join(f"{int(r[k])}" for k in ("bl", "br", "tr", "tl"))
            ref = float(r["g_el"]) * float(r["area"]) if e.disp == 1 else float(r["e_el"]) * float(r["inertia"])
            lst.write(f"{e.num:>6d} {e.group:>6d} {e.ident:>5d} {e.bbc:>4d} {HY.MODEL_CODES[e.code]:>5s} {e.disp:>5d} "
                      f"{corners:>27s} {float(r['length']):>10.4g} {float(r['height']):>10.4g} {float(r['thick']):>8.4g} "
                      f"{e.k_el:>12.5e} {ref:>12.5e} {e.backbone.k_el:>12.5e}")
            if ref > 0 and abs(e.backbone.k_el / ref - 1.0) > 0.01:
                lst.warning(f"panel {e.num}: BBC initial slope {e.backbone.k_el:.5g} differs by "
                            f"{100 * (e.backbone.k_el / ref - 1):.2f} % from {'G A_shear' if e.disp == 1 else 'E I'} = "
                            f"{ref:.5g} (EDU-09, D-NON-11)")
    if springs:
        lst.section("Nonlinear springs (S)")
        lst.write(f"{'spring':>6s} {'group':>6s} {'elem':>5s} {'SC':>5s} {'BBC':>4s} {'model':>5s} {'dof':>4s} {'node I':>7s} "
                  f"{'node J':>7s} {'k_el':>12s} {'K_el(BBC)':>12s}")
        for e in springs:
            r = e.row
            lst.write(f"{e.num:>6d} {e.group:>6d} {int(r['elem']):>5d} {e.ident:>5d} {e.bbc:>4d} {HY.MODEL_CODES[e.code]:>5s} "
                      f"{'XYZ'[e.disp - 1]:>4s} {int(r['node_i']):>7d} {int(r['node_j']):>7d} {e.k_el:>12.5e} "
                      f"{e.backbone.k_el:>12.5e}")
            if abs(e.backbone.k_el / e.k_el - 1.0) > 0.01:
                lst.warning(f"spring {e.num}: BBC initial slope {e.backbone.k_el:.5g} differs by "
                            f"{100 * (e.backbone.k_el / e.k_el - 1):.2f} % from the spring constant {e.k_el:.5g} "
                            "(EDU-09, D-NON-11)")


def _current_properties(elements: Sequence[NLElement], current: Dict[Tuple[str, int], PropRow], hou, lst: Listing,
                        iterating: bool) -> None:
    """Properties of the SSI analysis just made: elastic (elastic run) or *_EQL_Matl_Prop.txt (iteration);
    cross-check with the HOUSE deck."""
    mats = {int(r["id"]): r for r in hou.rows("materials")}
    scs = {int(r["id"]): r for r in hou.rows("springprops")}
    mismatch = []
    for el in elements:
        if iterating and (el.kind, el.num) in current:
            r = current[(el.kind, el.num)]
            el.k_cur, el.xi_cur = r.k_new, r.xi_new
        else:
            el.k_cur, el.xi_cur = el.k_el, el.xi_el
        # the HOUSE deck should hold these properties
        if el.kind == "panel":
            m = mats.get(el.ident)
            if m is None:
                raise ModuleError(f"panel {el.num}: material {el.ident} is not in {hou.params.get('model', '')}.hou")
            row = el.row
            v1, v2 = scale_material(int(row["mtype"]), float(row["val1"]), float(row["val2"]), el.k_cur / el.k_el)
            if (abs(float(m["val1"]) - v1) > 1e-6 * abs(v1) or abs(float(m["val2"]) - v2) > 1e-6 * max(abs(v2), 1e-30)
                    or abs(float(m["sdamp"]) - el.xi_cur) > 1e-9 * max(abs(el.xi_cur), 1.0)):
                mismatch.append(f"panel {el.num} (material {el.ident})")
        else:
            s = scs.get(el.ident)
            if s is None:
                raise ModuleError(f"spring {el.num}: SC property {el.ident} is not in the HOUSE deck")
            col = ("scx", "scy", "scz")[el.disp - 1]
            if (abs(float(s[col]) - el.k_cur) > 1e-6 * abs(el.k_cur)
                    or abs(float(s["damp"]) - el.xi_cur) > 1e-9 * max(abs(el.xi_cur), 1.0)):
                mismatch.append(f"spring {el.num} (SC {el.ident})")
    if mismatch:
        lst.warning("the HOUSE deck does not hold the properties of the analysis just made for "
                    + ", ".join(mismatch[:8]) + (" ..." if len(mismatch) > 8 else "")
                    + (": copy <model>_new.hou to <model>.hou before HOUSE (FCOPY) in every iteration" if iterating
                       else ": the elastic HOUSE deck should hold the M/SC values of the model (run AFWRITE)"))


def _write_outputs(ctx: ModuleContext, deck: EqlDeck, hou, elements: Sequence[NLElement], dt: float, analysis: int,
                   kinds: Sequence[str], iterating: bool, conv: HY.Convergence) -> List[str]:
    wd, model = ctx.workdir, ctx.model
    written: List[str] = []
    # histories and curves
    for el in elements:
        res = el.result
        what = "shear strain" if (el.kind == "panel" and el.disp == 1) else (
            "curvature" if el.kind == "panel" else f"relative displacement u_J - u_I ({'XYZ'[el.disp - 1]})")
        write_history(wd / f"{el.label}.thd", el.x, dt, f"{el.label} equivalent-linear deformation history: {what}, "
                                                        f"SSI analysis {analysis}")
        write_history(wd / f"{el.label}.ths", res.force, dt, f"{el.label} nonlinear force history "
                                                             f"({HY.MODEL_NAMES[el.code]}), SSI analysis {analysis}")
        written += [f"{el.label}.thd", f"{el.label}.ths"]
        if el.kind == "panel" and el.axial is not None:
            write_history(wd / f"{el.label}_AXIAL.thd", el.axial, dt, f"{el.label} uniform vertical axial strain history")
            written.append(f"{el.label}_AXIAL.thd")
        crv = wd / f"{el.label}.crv"
        if not iterating or not crv.exists():
            amps = HY.crv_amplitudes(el.backbone, int(deck["crvpoints"]))
            tab = HY.eql_curve(el.code, el.backbone, amps)
            lines = [f"# {el.label}: equivalent-linear curves of BBC {el.bbc} with the {HY.MODEL_NAMES[el.code]} model "
                     f"(stabilised cycles, {HY.LOOP_CYCLES} cycles per amplitude)",
                     f"# K_el = {el.backbone.k_el:.8e} (BBC first slope); E/E_el = K_sec/K_el on the model backbone; "
                     "xi_h = E_D/(4 pi E_S)",
                     f"# {'amplitude':>16s} {'E/E_el':>16s} {'xi_h':>16s} {'K_sec':>16s} {'F':>16s}"]
            lines += ["  " + " ".join(f"{v:>16.8e}" for v in row) for row in tab]
            crv.write_text("\n".join(lines) + "\n", encoding="utf-8")
            written.append(crv.name)
    # properties, ductility and the new HOUSE deck
    for kind in kinds:
        els = [e for e in elements if e.kind == kind]
        if not els:
            continue
        rows = [PropRow(e.num, e.group, e.ident, e.k_el, e.k_el * e.result.ratio, e.result.ratio, e.xi_el, e.result.xi_h,
                        e.result.xi, e.result.x_max, e.result.x_eq) for e in els]
        write_props(wd / PROP_FILES[kind], kind, model, analysis + 1, rows,
                    note=("CONVERGED: " if conv.converged else "") + conv.text())
        written.append(PROP_FILES[kind])
        lines = [f"# SASSI-EDU NONLINEAR ductility and force-reduction factors ({'wall panels' if kind == 'panel' else 'springs'}), "
                 f"SSI analysis {analysis}",
                 "# mu = max|x| / x_cr (x_cr = BBC point 1); F_mu = initial elastic force K_el max|x| / nonlinear force "
                 "|F| at the time of max|x| (ASCE 43-05 inelastic absorption factor)",
                 f"# {'element':>8s} {'max|x|':>15s} {'x_cr':>15s} {'mu':>15s} {'F_elastic':>15s} {'F_nonlinear':>15s} {'F_mu':>15s}"]
        for e in els:
            r = e.result
            fe = e.backbone.k_el * r.x_max
            lines.append(f"  {e.num:>8d} {r.x_max:>15.8e} {e.backbone.x_cr:>15.8e} {r.ductility:>15.8e} {fe:>15.8e} "
                         f"{(fe / r.f_mu if r.f_mu else 0.0):>15.8e} {r.f_mu:>15.8e}")
        (wd / FMU_FILES[kind]).write_text("\n".join(lines) + "\n", encoding="utf-8")
        written.append(FMU_FILES[kind])
        if not iterating:
            (wd / NON_FILES[kind]).write_text("1\n", encoding="utf-8")
            written.append(NON_FILES[kind])
    new = _new_house_deck(hou, elements)
    path = wd / f"{model}_new.hou"
    decks.write(path, new, comments=[f"written by NONLINEAR after SSI analysis {analysis}: equivalent-linear properties "
                                     f"of SSI iteration {analysis + 1} (copy to {model}.hou before HOUSE)"])
    written.append(path.name)
    return written


def _new_house_deck(hou, elements: Sequence[NLElement]):
    """Copy of the HOUSE deck with the panel materials and spring properties of the next iteration."""
    new = decks.new("HOUSE")
    new.params.update(hou.params)
    for name, tab in hou.tables.items():
        new.tables[name] = Table(columns=list(tab.columns), rows=[list(r) for r in tab.rows])
    mat_cols = new.tables["materials"].columns
    sc_cols = new.tables["springprops"].columns
    by_mat = {int(r[mat_cols.index("id")]): r for r in new.tables["materials"].rows}
    by_sc = {int(r[sc_cols.index("id")]): r for r in new.tables["springprops"].rows}
    for el in elements:
        ratio, xi = el.result.ratio, el.result.xi
        if el.kind == "panel":
            row = by_mat[el.ident]
            v1, v2 = scale_material(int(el.row["mtype"]), float(el.row["val1"]), float(el.row["val2"]), ratio)
            row[mat_cols.index("type")] = int(el.row["mtype"])
            row[mat_cols.index("val1")] = v1
            row[mat_cols.index("val2")] = v2
            row[mat_cols.index("weight")] = float(el.row["weight"])
            row[mat_cols.index("pdamp")] = xi
            row[mat_cols.index("sdamp")] = xi
        else:
            row = by_sc[el.ident]
            for c in ("scx", "scy", "scz", "scxx", "scyy", "sczz"):
                row[sc_cols.index(c)] = float(el.row[c])
            col = ("scx", "scy", "scz")[el.disp - 1]
            row[sc_cols.index(col)] = float(el.row["k_el"]) * ratio
            row[sc_cols.index("damp")] = xi
    return new


# ======================================================================================
# COMB_XYZ_THD (D-NON-08)
# ======================================================================================
def read_comb_inp(path: Path) -> List[Tuple[str, List[Optional[str]]]]:
    """``[(outfile, [fileX, fileY, fileZ])]`` of a COMB_XYZ_THD.inp ('-' / NONE = missing direction)."""
    lines = []
    for ln in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        s = ln.strip()
        if s and s[0] not in "#*!":
            lines.append(s)
    if not lines:
        raise ModuleError(f"{Path(path).name} is empty")
    try:
        n = int(float(lines[0].split()[0]))
    except ValueError:
        raise ModuleError(f"{Path(path).name}: line 1 must be the number of combinations") from None
    if len(lines) - 1 < n:
        raise ModuleError(f"{Path(path).name}: {n} combinations announced, {len(lines) - 1} given")
    out = []
    for k, ln in enumerate(lines[1:n + 1], start=2):
        t = ln.replace(",", " ").split()
        if len(t) < 2:
            raise ModuleError(f"{Path(path).name}: combination line {k - 1} needs 'outfile fileX fileY fileZ'")
        srcs = [None if (v in ("-", "0") or v.upper() == "NONE") else v for v in (t[1:4] + ["-"] * (3 - len(t[1:4])))]
        if not any(srcs):
            raise ModuleError(f"{Path(path).name}: combination {t[0]} has no input file")
        out.append((t[0], srcs))
    if len(lines) - 1 > n:
        out_extra = len(lines) - 1 - n
        out.append(("", [f"{out_extra} extra line(s) ignored", None, None]))
    return out


def comb_xyz_thd(workdir: PathLike, inp: str = COMB_INP) -> Tuple[List[str], List[str]]:
    """COMB_XYZ_THD: algebraic sum of the X, Y and Z histories of every combination (D-NON-08).

    Returns ``(files written, warnings)``.  All inputs of one combination must have the same time step
    and length."""
    wd = Path(workdir)
    p = wd / inp if not Path(inp).is_absolute() else Path(inp)
    if not p.exists():
        raise ModuleError(f"{p.name} not found in {wd} (line 1: n; then n lines 'outfile fileX fileY fileZ')")
    combos = read_comb_inp(p)
    written: List[str] = []
    warns: List[str] = []
    for out, srcs in combos:
        if not out:
            warns.append(srcs[0])
            continue
        total = None
        dt0 = None
        for s in srcs:
            if s is None:
                continue
            q = wd / s if not Path(s).is_absolute() else Path(s)
            if not q.exists():
                raise ModuleError(f"COMB_XYZ_THD: {s} (input of {out}) not found")
            v, dt = read_history(q)
            if total is None:
                total, dt0 = v.copy(), dt
            else:
                if abs(dt - dt0) > 1e-9 * max(dt0, 1e-30):
                    raise ModuleError(f"COMB_XYZ_THD: {s} has dt {dt:g}, the first input of {out} {dt0:g}")
                if len(v) != len(total):
                    raise ModuleError(f"COMB_XYZ_THD: {s} has {len(v)} values, the first input of {out} {len(total)}")
                total += v
        q = wd / out if not Path(out).is_absolute() else Path(out)
        write_history(q, total, float(dt0), header="COMB_XYZ_THD: " + " + ".join(s for s in srcs if s))
        written.append(q.name)
    return written, warns


# ======================================================================================
# Runner (sassi.modules.base.run_module knows the 11 P0 modules only) and batch entry
# ======================================================================================
def run_nonlinear(model: str, workdir: PathLike, deck_path: Optional[PathLike] = None,
                  listing_path: Optional[PathLike] = None, echo: Optional[Callable[[str], None]] = None,
                  progress: Optional[Callable[[float, str], None]] = None,
                  cancelled: Optional[Callable[[], bool]] = None) -> int:
    """Run NONLINEAR for ``model`` in ``workdir`` exactly as :func:`sassi.modules.base.run_module` runs the
    other modules (listing ``<model>_NONLINEAR.out``, status line); returns 0 on success."""
    workdir = Path(workdir)
    deck_path = Path(deck_path) if deck_path else workdir / f"{model}{DECK_EXT}"
    listing = Listing(Path(listing_path) if listing_path else workdir / f"{model}_{NAME}.out", echo=echo)
    ctx = ModuleContext(module=NAME, model=model, workdir=workdir, deck_path=deck_path, listing=listing,
                        progress=progress or (lambda f, m: None), cancelled=cancelled or (lambda: False))
    t0 = time.time()
    rc = 2
    try:
        listing.header(NAME)
        if not deck_path.exists():
            raise ModuleError(f"input deck {deck_path.name} not found -- run RUNNONLINEAR (it writes the deck)")
        rc = int(run(ctx) or 0)
    except ModuleError as exc:
        listing.error(str(exc))
        rc = 1
    except Exception as exc:   # unexpected: keep the traceback in the listing
        listing.error(f"unexpected failure: {exc!r}")
        listing.write(traceback.format_exc())
        rc = 2
    finally:
        listing.write("")
        listing.write(f" {NAME} finished with status {'OK' if rc == 0 else 'FAILED'} in {time.time() - t0:.2f} s; "
                      f"{len(listing.warnings)} warning(s), {len(listing.errors)} error(s)")
        listing.close()
    return rc


def main(argv: Optional[Sequence[str]] = None, stdin: Optional[TextIO] = None) -> int:
    """``python -m sassi.modules.nonlinear < NONLINEAR.inp`` or ``... COMB_XYZ_THD [file.inp]``."""
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0].upper() in ("COMB_XYZ_THD", "COMBXYZTHD"):
        try:
            written, warns = comb_xyz_thd(Path.cwd(), argv[1] if len(argv) > 1 else COMB_INP)
        except ModuleError as exc:
            print(f"COMB_XYZ_THD: {exc}", file=sys.stderr)
            return 1
        for w in warns:
            print(f"COMB_XYZ_THD: warning: {w}")
        print(f"COMB_XYZ_THD: {len(written)} files written")
        return 0
    stdin = stdin or sys.stdin
    lines = [ln.strip() for ln in stdin.read().splitlines() if ln.strip()]
    if len(lines) < 3:
        print(f"{NAME}: expected three input lines (model, deck, listing)", file=sys.stderr)
        return 2
    model, deck, out = lines[:3]
    wd = Path.cwd()
    return run_nonlinear(model, wd, deck_path=wd / deck, listing_path=wd / out, echo=print)


if __name__ == "__main__":
    raise SystemExit(main())
