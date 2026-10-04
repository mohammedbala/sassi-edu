"""Verification problems of LOADGEN (Option A, the ACS SASSI-ANSYS two-step approach; manual 6.4.15,
spec 04 section 15.6, :mod:`sassi.modules.loadgen`).

* VP-LA1  equivalent static loads.  A three-mass BEAMS stick on a rigid 4 m x 4 m surface mat of a stiff
          site (Vs = 1000 m/s, practically a fixed base) runs SITE -> POINT -> HOUSE -> ANALYS -> MOTION ->
          RELDISP -> STRESS under an earthquake-like record; LOADGEN ("ANSYS Eq. Static Load", Disp. and
          Accel., lumped masses generated from the HOUSE mass matrix, critical time = maximum base shear)
          writes the APDL.  By the dynamic equilibrium of the stick above the base the inertia forces
          ``F = -m a`` sum to minus the STRESS base shear (force on the base element at I, D-STR-05), and
          their moment about the base to minus the STRESS base moment: both within 1 % at the critical time
          (the residual is the different interpolation of the nodal TFs in MOTION and of the element STF in
          STRESS, D-STR-01).  The interface displacements equal the RELDISP histories (started at rest) to
          1e-10, the generated masses equal the MT masses to 1e-12, and the APDL parses.
* VP-LA2  dynamic loads.  The same stick on a 3 x 3 mat on frequency-independent soil springs (the
          flexible-volume equation with ``X = k I``, solved at every Fourier frequency, so MOTION does not
          interpolate): LOADGEN "ANSYS Dynamic Load" writes TABLE arrays.  The acceleration tables (ACEL of
          the ACC method = the mat-centre ``.ACC``; check tables of the stick nodes) reproduce the MOTION
          ``.ACC`` histories to 1e-10, the REL method's ACEL table the control motion and its D tables the
          RELDISP ``.THD`` (started at rest) to 1e-10, the FILE8 source the RESULTS source to 1e-9, and the
          APDL parses with :func:`parse_apdl` (``*DIM`` / assignments / ``*VFILL`` / ``*IF``, D, F, ACEL,
          solution commands).  Finally the *second step itself* is replayed: the structure with the exported
          boundary conditions (interface D tables + ACEL) is solved at every Fourier frequency with the HOUSE
          matrices (COOSK/COOSM) and reproduces the SSI absolute accelerations of the stick (MOTION ``.ACC``)
          to 1e-6 -- the exactness of the decomposition ``u = u_r + iota u_g`` behind Option A (same damping).

:func:`parse_apdl` is a deliberately small APDL reader: it executes only what LOADGEN writes (scalar and
array parameters, TABLE fills, ``*VFILL ... RAMP``, ``*IF`` blocks) and records every other command; it
reports lines longer than 640 characters, names longer than 32 characters, more than 10 values per array
assignment, unfilled table cells and references to undefined tables.  The unit tests use it too.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

from ...conventions import nodal_result_name
from ...io import decks, textfiles
from ...io.container import read_container
from ...modules import loadgen as LG
from ...modules.base import run_module
from .. import VPResult, problem, worse

G = 9.81
KNOWN_COMMANDS = {"/SOLU", "/PREP7", "/POST1", "FINISH", "ANTYPE", "TRNOPT", "TIMINT", "AUTOTS", "KBC", "DELTIM",
                  "NSUBST", "TIME", "ALPHAD", "BETAD", "OUTRES", "ACEL", "D", "F", "SOLVE", "LSWRITE", "LSSOLVE"}
D_LABELS = {"UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ", "ALL"}
F_LABELS = {"FX", "FY", "FZ", "MX", "MY", "MZ"}
_NAME_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")
_ASSIGN_RE = re.compile(r"^([A-Za-z][A-Za-z0-9_]*)\s*(?:\(\s*([^)]*)\))?\s*=\s*(.*)$")
_TABLE_REF = re.compile(r"^%([A-Za-z][A-Za-z0-9_]*)%$")


# =======================================================================================
# A small APDL reader for the LOADGEN output
# =======================================================================================
@dataclass
class ApdlArray:
    """``*DIM`` parameter: TABLE arrays keep row 0 and column 0 (index values) in ``data``."""
    name: str
    kind: str
    imax: int
    jmax: int
    vars: List[str]
    data: np.ndarray

    def column(self, j: int = 1) -> np.ndarray:
        return self.data[1:, j] if self.kind == "TABLE" else self.data[:, j - 1]

    @property
    def index(self) -> np.ndarray:
        """Row index values (column 0) of a TABLE: the times."""
        return self.data[1:, 0]


@dataclass
class ApdlCommand:
    line: int
    name: str
    fields: List[str]


@dataclass
class ApdlProgram:
    scalars: Dict[str, float] = field(default_factory=dict)
    arrays: Dict[str, ApdlArray] = field(default_factory=dict)
    commands: List[ApdlCommand] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)
    nlines: int = 0
    max_values: int = 0

    def find(self, name: str) -> List[ApdlCommand]:
        return [c for c in self.commands if c.name == name.upper()]

    def loads(self, name: str) -> List[Tuple[int, str, object]]:
        """``D`` / ``F`` commands as (node, label, value) with value a float or ``'%TABLE%'`` name."""
        out = []
        for c in self.find(name):
            node = int(float(c.fields[0]))
            lab = c.fields[1].strip().upper()
            v = c.fields[2].strip()
            m = _TABLE_REF.match(v)
            out.append((node, lab, m.group(1).upper() if m else self.value(v, c.line)))
        return out

    def value(self, tok: str, line: int = 0) -> float:
        t = tok.strip()
        if not t:
            return 0.0
        try:
            return float(t)
        except ValueError:
            pass
        if t.upper() in self.scalars:
            return self.scalars[t.upper()]
        self.errors.append(f"line {line}: cannot evaluate '{t}'")
        return float("nan")

    def table(self, name: str) -> ApdlArray:
        return self.arrays[name.upper()]


def _strip_comment(line: str) -> str:
    k = line.find("!")
    return line if k < 0 else line[:k]


def _split_index(text: str) -> List[str]:
    """Split on the commas that are not inside parentheses (``*VFILL,T(1,0),RAMP,0,DT``)."""
    out, depth, cur = [], 0, []
    for ch in text:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch == "," and depth == 0:
            out.append("".join(cur).strip())
            cur = []
        else:
            cur.append(ch)
    out.append("".join(cur).strip())
    return out


def parse_apdl(text: str, max_line: int = 640) -> ApdlProgram:
    """Parse/execute the subset of APDL that LOADGEN writes (see the module docstring)."""
    prog = ApdlProgram()
    stack: List[bool] = []           # *IF states (True = executing)

    def active() -> bool:
        return all(stack)

    def fill(name: str, idx: List[str], values: List[str], ln: int) -> None:
        arr = prog.arrays.get(name.upper())
        if arr is None:
            prog.errors.append(f"line {ln}: assignment to undefined array {name}")
            return
        if len(values) > 10:
            prog.errors.append(f"line {ln}: {len(values)} values in one assignment (APDL allows 10)")
        prog.max_values = max(prog.max_values, len(values))
        i = int(prog.value(idx[0], ln)) if idx and idx[0] else 1
        j = int(prog.value(idx[1], ln)) if len(idx) > 1 and idx[1] else (1 if arr.kind != "TABLE" else 1)
        for off, v in enumerate(values):
            r = i + off
            if arr.kind == "TABLE":
                if not (0 <= r <= arr.imax and 0 <= j <= arr.jmax):
                    prog.errors.append(f"line {ln}: {name}({r},{j}) outside the table")
                    return
                arr.data[r, j] = prog.value(v, ln)
            else:
                if not (1 <= r <= arr.imax and 1 <= j <= arr.jmax):
                    prog.errors.append(f"line {ln}: {name}({r},{j}) outside the array")
                    return
                arr.data[r - 1, j - 1] = prog.value(v, ln)

    for ln, raw in enumerate(text.splitlines(), start=1):
        prog.nlines = ln
        if len(raw) > max_line:
            prog.errors.append(f"line {ln}: {len(raw)} characters (ANSYS reads at most {max_line})")
        body = _strip_comment(raw).strip()
        if not body:
            continue
        head = body.split(",")[0].strip().upper()
        if head == "*IF":
            f = _split_index(body)
            if len(f) < 5 or f[4].upper() != "THEN":
                prog.errors.append(f"line {ln}: unsupported *IF form")
                stack.append(False)
                continue
            a, op, b = prog.value(f[1], ln), f[2].upper(), prog.value(f[3], ln)
            ok = {"EQ": a == b, "NE": a != b, "LT": a < b, "GT": a > b, "LE": a <= b, "GE": a >= b}.get(op)
            if ok is None:
                prog.errors.append(f"line {ln}: unknown *IF operator {op}")
                ok = False
            stack.append(bool(ok))
            continue
        if head == "*ELSE":
            if stack:
                stack[-1] = not stack[-1]
            continue
        if head == "*ENDIF":
            if not stack:
                prog.errors.append(f"line {ln}: *ENDIF without *IF")
            else:
                stack.pop()
            continue
        if not active():
            continue
        m = _ASSIGN_RE.match(body)
        if m and not body.startswith("*") and not body.startswith("/"):
            name, idx, rhs = m.group(1), m.group(2), m.group(3)
            if len(name) > 32:
                prog.errors.append(f"line {ln}: parameter name {name} longer than 32 characters")
            vals = _split_index(rhs)
            if idx is None:
                if len(vals) != 1:
                    prog.errors.append(f"line {ln}: scalar {name} assigned {len(vals)} values")
                prog.scalars[name.upper()] = prog.value(vals[0], ln)
            else:
                fill(name, _split_index(idx), vals, ln)
            continue
        f = _split_index(body)
        cmd = f[0].upper()
        if cmd == "*DIM":
            name = f[1]
            if not _NAME_RE.match(name) or len(name) > 32:
                prog.errors.append(f"line {ln}: invalid parameter name {name!r}")
            kind = f[2].upper()
            imax, jmax = int(prog.value(f[3], ln)), int(prog.value(f[4], ln) if len(f) > 4 and f[4] else 1)
            kmax = int(prog.value(f[5], ln)) if len(f) > 5 and f[5] else 1
            if kmax != 1:
                prog.errors.append(f"line {ln}: 3-D arrays are not used by LOADGEN")
            if kind == "TABLE":
                data = np.full((imax + 1, jmax + 1), np.nan)
                data[0, :] = 0.0
                data[0, 0] = 0.0
            elif kind == "ARRAY":
                data = np.zeros((imax, jmax))
            else:
                prog.errors.append(f"line {ln}: *DIM type {kind} not supported")
                continue
            prog.arrays[name.upper()] = ApdlArray(name.upper(), kind, imax, jmax, [v.upper() for v in f[6:] if v], data)
            continue
        if cmd == "*SET":
            fill(f[1], [], f[2:], ln) if f[1].upper() in prog.arrays else prog.scalars.__setitem__(
                f[1].upper(), prog.value(f[2], ln))
            continue
        if cmd == "*VFILL":
            mm = re.match(r"^([A-Za-z][A-Za-z0-9_]*)\s*\(\s*([^)]*)\)$", f[1].strip())
            if not mm or f[2].upper() != "RAMP":
                prog.errors.append(f"line {ln}: unsupported *VFILL form")
                continue
            arr = prog.arrays.get(mm.group(1).upper())
            if arr is None:
                prog.errors.append(f"line {ln}: *VFILL of undefined array {mm.group(1)}")
                continue
            idx = _split_index(mm.group(2))
            i0, j = int(prog.value(idx[0], ln)), int(prog.value(idx[1], ln))
            c1, c2 = prog.value(f[3], ln), prog.value(f[4], ln)
            if arr.kind == "TABLE":
                rows = np.arange(i0, arr.imax + 1)
                arr.data[rows, j] = c1 + (rows - i0) * c2
            else:
                rows = np.arange(i0, arr.imax + 1)
                arr.data[rows - 1, j - 1] = c1 + (rows - i0) * c2
            continue
        if cmd.startswith("*"):
            prog.errors.append(f"line {ln}: unsupported command {cmd}")
            continue
        if cmd not in KNOWN_COMMANDS:
            prog.errors.append(f"line {ln}: unexpected command {cmd}")
        args = f[1:]
        prog.commands.append(ApdlCommand(ln, cmd, args))
        if cmd in ("D", "F"):
            lab = args[1].strip().upper() if len(args) > 1 else ""
            if lab not in (D_LABELS if cmd == "D" else F_LABELS):
                prog.errors.append(f"line {ln}: invalid {cmd} label {lab!r}")
            refs = [args[2]] if len(args) > 2 else []
        elif cmd == "ACEL":
            refs = args[:3]
        else:
            refs = []
        for r_ in refs:
            mt = _TABLE_REF.match(r_.strip())
            if mt and (mt.group(1).upper() not in prog.arrays or prog.arrays[mt.group(1).upper()].kind != "TABLE"):
                prog.errors.append(f"line {ln}: %{mt.group(1)}% is not a defined TABLE")
    if stack:
        prog.errors.append("unterminated *IF block")
    for arr in prog.arrays.values():
        if arr.kind == "TABLE" and np.any(np.isnan(arr.data[1:, :])):
            prog.errors.append(f"table {arr.name}: {int(np.isnan(arr.data[1:, :]).sum())} cells not filled")
    return prog


def parse_file(path: Path) -> ApdlProgram:
    return parse_apdl(Path(path).read_text(encoding="utf-8"))


# =======================================================================================
# Model and chain helpers
# =======================================================================================
NFFT, DT = 1024, 0.02
STICK = dict(heights=(4.0, 8.0, 12.0), masses=(2.0, 2.0, 1.5), E=3.0e7, A=0.5, I=0.02, nu=0.2, beta=0.05, J=0.04)


def _stick(B):
    return B.Stick(heights=list(STICK["heights"]), masses=list(STICK["masses"]), E=STICK["E"], A=STICK["A"],
                   I=STICK["I"], nu=STICK["nu"], beta=STICK["beta"], J=STICK["J"])


def control_motion(wd: Path, n: int = 600, peak: float = 0.3, seed: int = 5, dt: float = DT) -> np.ndarray:
    """Earthquake-like control motion ``eq.acc`` (g, peak ``peak``) in the model directory."""
    from .vp_motion import synthetic_motion, write_motion_file
    acc = peak * synthetic_motion(n, dt, seed=seed, fmax=15.0)
    write_motion_file(wd / "eq.acc", acc, dt)
    return acc


def run_motion_reldisp(wd: Path, B, fs, nout_nodes: Sequence[Tuple[int, Sequence[int]]], rd_nodes: Sequence[int],
                       model: str = "m", rd_dofs: Sequence[int] = (1, 2, 3), savetf: int = 1) -> Tuple[int, int]:
    """MOTION (.ACC and .TFI of the listed nodes/DOFs; ``savetf`` = every translational TF) and RELDISP (free
    field) of ``rd_dofs`` at ``rd_nodes``."""
    rows = [(n, k, 1, 1, 0, 0, 0, 1) for n, dofs in nout_nodes for k in dofs]
    md = B.motion_deck(fs, "eq.acc", rows, gravity=G, model=model, cplx=1, savetf=savetf)
    B.write_deck(wd, model, md)
    rc1 = run_module("MOTION", model, wd)
    rd = decks.new("RELDISP")
    rd["model"], rd["thfile"], rd["nft"], rd["delt"] = model, "eq.acc", fs.nft, fs.delt
    rd["df"], rd["gravity"], rd["mult"], rd["max"] = 1.0 / (fs.nft * fs.delt), G, 1.0, 0.0
    for n in rd_nodes:
        rd.table("rdnd").append([int(n)] + [1 if k in rd_dofs else 0 for k in range(1, 7)])
    rd["numfiles"] = len(rd_nodes)
    B.write_deck(wd, model, rd)
    rc2 = run_module("RELDISP", model, wd)
    return rc1, rc2


def spring_model(wd: Path, B, nfft: int = NFFT, dt: float = DT, k_spring: float = 5.0e3, optimize: bool = False,
                 model: str = "m"):
    """The VP stick on a 3 x 3 rigid mat whose interaction translations sit on frequency-independent soil springs
    ``k_spring``: HOUSE (optionally with the node optimizer) and FILE8 of the flexible-volume equation with
    ``X = k I`` at every Fourier frequency below the Nyquist frequency (vp_stress.stiff_site_file8), so MOTION
    reproduces the exact TFs without interpolation.  Returns (FrequencySet, builders Model)."""
    from .vp_stress import stiff_site_file8
    site = B.uniform_site(depth=20.0, nsub=4, vs=300.0, rho=2.0, nu=1.0 / 3.0, beta=0.05, gravity=G)
    fnum = np.arange(1, nfft // 2)
    fs = B.FrequencySet.fourier(dt, nfft, fnum)
    mdl = B.stick_on_mat(site, _stick(B), half_width=2.0, ndiv=2)
    B.write_deck(wd, model, B.site_deck(site, fs, model=model))
    hd = mdl.house.deck(model)
    if optimize:
        hd["optimize"] = 1
    B.write_deck(wd, model, hd)
    B.run("HOUSE", wd, model)
    stiff_site_file8(wd, fnum, 1.0 / (nfft * dt), k_spring, nfft, dt, cm=0, model=model)
    return fs, mdl


def run_loadgen(wd: Path, model: str = "m", **params) -> Tuple[int, str]:
    """Write ``<model>.lgn`` with ``params`` (deck parameters; ``tables`` = {name: nodes}) and run LOADGEN."""
    d = LG.new_deck()
    d["model"] = model
    tables = params.pop("tables", {})
    for k, v in params.items():
        d[k] = v
    for k, v in tables.items():
        d.tables[k] = list(v)
    LG.write_deck(LG.deck_path(wd, model), d)
    rc = LG.run_loadgen(model, wd)
    return rc, (wd / f"{model}_LOADGEN.out").read_text(encoding="utf-8")


def history(wd: Path, node: int, dof: int, ext: str) -> np.ndarray:
    v, _ = textfiles.read_history(wd / nodal_result_name(node, dof, ext))
    return v


def _relmax(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    den = float(np.max(np.abs(b))) or 1.0
    return float(np.max(np.abs(a - b))) / den


# =======================================================================================
# VP-LA1: equivalent static loads vs STRESS base shear
# =======================================================================================
FNUM_LA1 = np.unique(np.concatenate([np.arange(1, 40), np.arange(40, 120, 2), np.arange(120, 300, 5)]))


@problem("VP-LA1", "LOADGEN equivalent static loads: inertia forces of the stick sum to the STRESS base shear",
         tier="P2", modules=["LOADGEN", "STRESS", "MOTION", "RELDISP", "ANALYS", "HOUSE", "SITE", "POINT"],
         source="manual 6.4.15, spec 04 15.6; wave-3 package option_a_loadgen (VP-LA1)")
def vp_la1(workdir) -> VPResult:
    from .. import builders as B
    from .vp_stress import add_eout, read_ths, run_stress, stress_deck
    r = VPResult()
    wd = Path(workdir)
    wd.mkdir(parents=True, exist_ok=True)
    site = B.uniform_site(depth=20.0, nsub=4, vs=1000.0, rho=2.0, nu=1.0 / 3.0, beta=0.02, gravity=G)
    fs = B.FrequencySet.fourier(DT, NFFT, FNUM_LA1)
    mdl = B.stick_on_mat(site, _stick(B), half_width=2.0, ndiv=2)
    try:
        B.run_soil(wd, "m", site, fs, layer=mdl.layer, rad=mdl.rad)
        f4 = B.run_house(wd, "m", mdl)
        B.run_analys(wd, "m", fs, gravity=G)
    except B.ChainError as exc:
        r.require("SITE, POINT, HOUSE and ANALYS runs succeed", False, str(exc)[-600:])
        return r
    control_motion(wd)
    base, stick, mat = mdl["centre"], list(mdl["stick"]), list(mdl["mat"])
    rc1, rc2 = run_motion_reldisp(wd, B, fs, [(n, (1, 2, 3)) for n in [base] + stick], mat)
    r.require("MOTION and RELDISP runs succeed", rc1 == 0 and rc2 == 0)
    if rc1 or rc2:
        return r
    hit = [k for k in range(len(f4["elem_id"])) if int(f4["elem_type"][k]) == 2
           and list(f4["elem_nodes"][k][:2]) == [base, stick[0]]]
    r.require("base element of the stick found in FILE4", len(hit) == 1)
    if len(hit) != 1:
        return r
    g, e = int(f4["elem_group"][hit[0]]), int(f4["elem_id"][hit[0]])
    dk = stress_deck(NFFT, DT, thfile="eq.acc", gravity=G)
    add_eout(dk, g, [e], [2] * 12)
    rc, out = run_stress(wd, dk)
    r.require("STRESS run succeeds", rc == 0, out[-600:] if rc else "")
    if rc:
        return r
    fyi, dts = read_ths(wd, "BEAMS", g, e, "FYI")
    mzi, _ = read_ths(wd, "BEAMS", g, e, "MZI")
    rc, lst = run_loadgen(wd, analysis="STATIC", data=3, genmass=1, crit="V", ncrit=1)
    r.require("LOADGEN static run succeeds", rc == 0, lst[-800:] if rc else "")
    if rc:
        return r
    # ---- generated lumped masses = the MT masses of the stick (massless beams and links)
    masl = LG.read_lumped(wd / "m.masl")
    for n, m in zip(stick, STICK["masses"]):
        r.check(f"generated lumped mass of stick node {n} (X) = MT mass", masl.get(n, np.zeros(6))[0], m, rtol=1e-12)
    r.check("total generated mass (X)", sum(v[0] for v in masl.values()), sum(STICK["masses"]), rtol=1e-12)
    # ---- the APDL file
    prog = parse_file(wd / "m_LGS.inp")
    r.require("APDL of the static load parses without errors", not prog.errors, "; ".join(prog.errors[:5]))
    res = np.loadtxt(wd / "m_LGS_res.txt")
    t, VX, MY = res[:, 0], res[:, 1], res[:, 5]
    k = int(np.argmax(np.abs(VX)))
    ks = int(np.argmax(np.abs(fyi)))
    r.check("critical time (max |base shear| of the inertia forces) = time of max |FYI| of STRESS", k * DT, ks * dts,
            atol=1.01 * DT)
    forces = prog.loads("F")
    fx = sum(v for (_, lab, v) in forces if lab == "FX")
    r.check("APDL: sum of the FX inertia forces = - STRESS base shear FYI at the critical time", fx, -fyi[k], rtol=0.01)
    h = {n: z for n, z in zip(stick, STICK["heights"])}
    mom = sum(h[nd] * v for (nd, lab, v) in forces if lab == "FX")
    r.check("APDL: moment of the inertia forces about the base = - STRESS base moment MZI at the critical time", mom,
            -mzi[k], rtol=0.01)
    r.check("resultant history: max |VX| = max |FYI| of STRESS", np.max(np.abs(VX)), np.max(np.abs(fyi)), rtol=0.01)
    r.check("resultant history: max |MY| about the base = max |MZI| of STRESS", np.max(np.abs(MY)), np.max(np.abs(mzi)),
            rtol=0.01)
    r.inform("resultant history: max_t |VX + FYI| / max |FYI| (interpolation of nodal TFs vs element STFs)",
             _relmax(-VX[:len(fyi)], fyi), 0.0)
    # ---- interface displacements = RELDISP .THD started at rest
    dmax = 0.0
    nd = 0
    for (node, lab, v) in prog.loads("D"):
        dof = {"UX": 1, "UY": 2, "UZ": 3}[lab]
        thd = history(wd, node, dof, "THD")
        ref = thd[k] - thd[0]
        dmax = worse(dmax, abs(v - ref) / max(np.max(np.abs(thd - thd[0])), 1e-30))
        nd += 1
    r.require("APDL: D at the 9 interaction nodes x 3 translations", nd == 27, f"{nd} D commands")
    r.check("APDL: D values = RELDISP .THD(t*) - .THD(0), max error / peak", dmax, 0.0, atol=1e-10)
    # ---- FILE8 source = RESULTS source
    rc, lst8 = run_loadgen(wd, analysis="STATIC", data=3, genmass=1, crit="V", ncrit=1, source="FILE8",
                           apdlfile="m_LGS8.inp")
    r.require("LOADGEN static run with source FILE8 succeeds", rc == 0, lst8[-800:] if rc else "")
    if rc == 0:
        p8 = parse_file(wd / "m_LGS8.inp")
        a = np.array([v for (_, _, v) in prog.loads("F")])
        b = np.array([v for (_, _, v) in p8.loads("F")])
        r.check("source FILE8 = source RESULTS: inertia forces, max difference / max", _relmax(b, a) if len(a) == len(b)
                else 1.0, 0.0, atol=1e-9)
    r.notes.append(f"stick: heights {STICK['heights']} m, masses {STICK['masses']} t, EI {STICK['E'] * STICK['I']:g} kN m2 "
                   f"on a 4 m x 4 m rigid mat, site Vs 1000 m/s; record {600 * DT:g} s, PGA 0.3 g; {len(FNUM_LA1)} SSI "
                   f"frequencies, NFFT {NFFT}, dt {DT} s")
    r.notes.append(f"critical time {k * DT:g} s: sum FX {fx:.6g}, STRESS FYI {fyi[k]:.6g}; moment {mom:.6g}, MZI "
                   f"{mzi[k]:.6g} (kN, kN m)")
    return r


# =======================================================================================
# VP-LA2: dynamic loads (tables) and the replay of the second step
# =======================================================================================
K_SPRING = 5.0e3          # frequency-independent soil spring per interaction translation, kN/m


def replay_second_step(wd: Path, prog: ApdlProgram, dnodes: Sequence[int], model: str = "m") -> Dict[Tuple[int, int], np.ndarray]:
    """Solve the ANSYS second step of an REL dynamic export in the frequency domain with the HOUSE matrices:
    interface translations prescribed by the D tables, ``ACEL`` = reference-frame acceleration, every other
    equation free; returns the absolute accelerations (g) of every free equation (node, dof)."""
    f4 = read_container(wd / f"{model}.N4", "FILE4")
    Ks = read_container(wd / "COOSK", "COOSK").sparse("Ks").toarray()
    Ms = read_container(wd / "COOSM", "COOSM").sparse("Ms").toarray()
    eqn, eqd = np.asarray(f4["eq_node"]), np.asarray(f4["eq_dof"])
    dtabs = {(n, lab): v for (n, lab, v) in prog.loads("D")}
    b = [i for i, (n, k) in enumerate(zip(eqn, eqd)) if int(n) in set(dnodes) and k <= 3]
    s = [i for i in range(len(eqn)) if i not in set(b)]
    lab = {1: "UX", 2: "UY", 3: "UZ"}
    nt = int(prog.scalars["LG_NT"])
    dt = prog.scalars["LG_DT"]
    if nt != NFFT:
        raise ValueError("the replay needs the full Fourier period (dur = 0)")
    ub = np.zeros((NFFT // 2 + 1, len(b)), dtype=complex)
    for c, i in enumerate(b):
        ub[:, c] = np.fft.rfft(prog.table(dtabs[(int(eqn[i]), lab[int(eqd[i])])]).column(1))
    acel = [c.fields for c in prog.find("ACEL")][0]
    ag = np.zeros((NFFT // 2 + 1, 3), dtype=complex)
    for k in range(3):
        m = _TABLE_REF.match(acel[k].strip())
        if m:
            ag[:, k] = np.fft.rfft(prog.table(m.group(1)).column(1))
    iota = np.zeros((len(eqn), 3))
    for i, k in enumerate(eqd):
        if k <= 3:
            iota[i, k - 1] = 1.0
    w = 2.0 * np.pi * np.fft.rfftfreq(NFFT, dt)
    a_abs = np.zeros((NFFT // 2 + 1, len(s)), dtype=complex)
    Kss, Ksb = Ks[np.ix_(s, s)], Ks[np.ix_(s, b)]
    Mss, Msb = Ms[np.ix_(s, s)], Ms[np.ix_(s, b)]
    with np.errstate(divide="ignore", over="ignore", invalid="ignore", under="ignore"):   # Accelerate BLAS flags
        Mi = Ms[s, :] @ iota                               # (ns, 3)
        for q in range(1, NFFT // 2):                      # DC and Nyquist: rigid body / zero (MOTION's rules)
            w2 = w[q] ** 2
            rhs = -(Ksb - w2 * Msb) @ ub[q] - Mi @ ag[q]   # ANSYS: -M iota ACEL (frame acceleration)
            ur = np.linalg.solve(Kss - w2 * Mss, rhs)
            a_abs[q] = -w2 * ur + iota[s] @ ag[q]
        a_abs[0] = iota[s] @ ag[0]
    if not np.all(np.isfinite(a_abs)):
        raise ValueError("non-finite replay result")
    acc = np.fft.irfft(a_abs, n=NFFT, axis=0) / G
    return {(int(eqn[i]), int(eqd[i])): acc[:, c] for c, i in enumerate(s)}


@problem("VP-LA2", "LOADGEN dynamic loads: APDL tables reproduce the SSI histories and the second step replays the SSI",
         tier="P2", modules=["LOADGEN", "MOTION", "RELDISP", "HOUSE"],
         source="manual 6.4.15, spec 04 15.6; wave-3 package option_a_loadgen (VP-LA2)")
def vp_la2(workdir) -> VPResult:
    from .. import builders as B
    r = VPResult()
    wd = Path(workdir)
    wd.mkdir(parents=True, exist_ok=True)
    try:
        fs, mdl = spring_model(wd, B, k_spring=K_SPRING)
    except B.ChainError as exc:
        r.require("HOUSE run succeeds", False, str(exc)[-600:])
        return r
    acc = control_motion(wd)
    base, stick, mat = mdl["centre"], list(mdl["stick"]), list(mdl["mat"])
    rc1, rc2 = run_motion_reldisp(wd, B, fs, [(n, (1, 2, 3)) for n in [base] + stick], mat)
    r.require("MOTION and RELDISP runs succeed", rc1 == 0 and rc2 == 0)
    if rc1 or rc2:
        return r
    # ---- (a) REL method, RESULTS source, check tables of the stick nodes
    rc, lst = run_loadgen(wd, analysis="DYNAMIC", method="REL", alpha=0.4, beta=0.002, tables={"checknodes": stick})
    r.require("LOADGEN dynamic run (REL, RESULTS) succeeds", rc == 0, lst[-800:] if rc else "")
    if rc:
        return r
    prog = parse_file(wd / "m_LGD.inp")
    r.require("APDL of the dynamic load parses without errors", not prog.errors, "; ".join(prog.errors[:5]))
    r.require("APDL: at most 10 values per array assignment", 0 < prog.max_values <= 10, f"{prog.max_values}")
    tabs = [a for a in prog.arrays.values() if a.kind == "TABLE"]
    tgrid = np.arange(NFFT) * DT
    r.check("APDL: every TABLE time column = (i - 1) dt, max |error| (s)", max(float(np.max(np.abs(a.index - tgrid)))
                                                                           for a in tabs), 0.0, atol=1e-12)
    r.require("APDL: transient solution commands (ANTYPE,TRANS; TRNOPT,FULL; DELTIM = dt; TIME = (n-1) dt; SOLVE)",
              [c.fields[0].upper() for c in prog.find("ANTYPE")] == ["TRANS"] and
              [c.fields[0].upper() for c in prog.find("TRNOPT")] == ["FULL"] and
              abs(prog.value(prog.find("DELTIM")[0].fields[0]) - DT) < 1e-15 and
              abs(prog.value(prog.find("TIME")[0].fields[0]) - (NFFT - 1) * DT) < 1e-9 and len(prog.find("SOLVE")) == 1)
    r.check("APDL: ALPHAD", prog.value(prog.find("ALPHAD")[0].fields[0]), 0.4, rtol=1e-15)
    r.check("APDL: BETAD", prog.value(prog.find("BETAD")[0].fields[0]), 0.002, rtol=1e-15)
    ax = prog.table("LG_ACX").column(1)
    ag = np.zeros(NFFT)
    ag[:len(acc)] = acc
    r.check("REL: ACEL table / g = control motion, max error / peak", _relmax(ax / G, ag), 0.0, atol=1e-10)
    errs = []
    for n in stick:
        errs.append(_relmax(prog.table(f"LGA_{n}_UX").column(1) / G, history(wd, n, 1, "ACC")))
    r.check("check tables / g = MOTION .ACC of the stick nodes (X), max error / peak", worse(*errs), 0.0, atol=1e-10)
    errs = []
    for (node, lab, name) in prog.loads("D"):
        dof = {"UX": 1, "UY": 2, "UZ": 3}[lab]
        thd = history(wd, node, dof, "THD")
        errs.append(_relmax(prog.table(name).column(1), thd - thd[0]) if np.max(np.abs(thd)) > 0 else 0.0)
    r.require("REL: D tables at the 9 interaction nodes x 3 translations", len(errs) == 27, f"{len(errs)}")
    r.check("REL: D tables = RELDISP .THD - .THD(0), max error / peak", worse(*errs), 0.0, atol=1e-10)
    # ---- (b) FILE8 source = RESULTS source
    rc, lst8 = run_loadgen(wd, analysis="DYNAMIC", method="REL", alpha=0.4, beta=0.002, source="FILE8",
                           apdlfile="m_LGD8.inp", tables={"checknodes": stick})
    r.require("LOADGEN dynamic run (REL, FILE8) succeeds", rc == 0, lst8[-800:] if rc else "")
    if rc == 0:
        p8 = parse_file(wd / "m_LGD8.inp")
        diffs = [_relmax(p8.table(nm).column(1), a.column(1)) for nm, a in prog.arrays.items() if a.kind == "TABLE"
                 and nm in p8.arrays]
        r.require("FILE8 source writes the same tables", set(p8.arrays) == set(prog.arrays))
        r.check("source FILE8 = source RESULTS: every table, max difference / peak", max(diffs), 0.0, atol=1e-9)
    # ---- (c) ACC method: ACEL = absolute acceleration of the mat centre
    rc, lsta = run_loadgen(wd, analysis="DYNAMIC", method="ACC", refnode=base, alpha=0.4, beta=0.002,
                           apdlfile="m_LGDA.inp")
    r.require("LOADGEN dynamic run (ACC, reference node = mat centre) succeeds", rc == 0, lsta[-800:] if rc else "")
    if rc == 0:
        pa = parse_file(wd / "m_LGDA.inp")
        r.require("APDL (ACC) parses without errors", not pa.errors, "; ".join(pa.errors[:5]))
        errs = []
        for lab, dof in (("X", 1), ("Y", 2), ("Z", 3)):
            a_ref = history(wd, base, dof, "ACC")
            if f"LG_AC{lab}" in pa.arrays:
                errs.append(_relmax(pa.table(f"LG_AC{lab}").column(1) / G, a_ref))
            else:
                errs.append(float(np.max(np.abs(a_ref)) > 0) * 1.0)
        r.check("ACC: ACEL tables / g = MOTION .ACC of the reference node (X, Y, Z), max error / peak", worse(*errs), 0.0,
                atol=1e-10)
        dz = pa.loads("D")
        r.require("ACC: the interface nodes are fixed (D = 0)", len(dz) > 0 and all(v == 0.0 for (_, _, v) in dz))
    # ---- (d) replay of the second step (REL export, HOUSE matrices, same hysteretic damping)
    rep = replay_second_step(wd, prog, mat)
    errs = [_relmax(rep[(n, 1)], history(wd, n, 1, "ACC")) for n in stick]
    r.check("replay: 2nd step with the exported D + ACEL reproduces the SSI absolute accelerations of the stick "
            "(MOTION .ACC), max error / peak", worse(*errs), 0.0, atol=1e-6)
    peak_rel = max(float(np.max(np.abs(prog.table(nm).column(1)))) for (_, _, nm) in prog.loads("D"))
    r.notes.append(f"soil springs k = {K_SPRING:g} kN/m per interaction translation (rocking and sliding of the mat); "
                   f"{NFFT // 2 - 1} SSI frequencies = every Fourier frequency below Nyquist (no interpolation); largest "
                   f"interface relative displacement {peak_rel:.4g} m")
    return r
