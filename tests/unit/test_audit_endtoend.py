"""Final audit: independent end-to-end checks of cross-module conventions (docs/verification/AUDIT_REPORT.md).

1. **Units.**  The same model is run through the command language in SI units (m, kN, t) and in British
   units (ft, kip, kip s^2/ft; g = 9.80665/0.3048 ft/s^2) through SOIL, SITE, POINT, HOUSE, ANALYS,
   MOTION, STRESS, RELDISP and LOADGEN (Option A APDL: forces scale with force, D values and the
   ACEL / displacement tables with length).  Every dimensionless result must agree to round-off and every
   dimensional one must scale exactly with the unit factors: transfer functions of translations,
   accelerations and response spectra in g, SOIL strains in % are identical; rotations scale with
   1/length, relative displacements with length, beam forces with force and moments with
   force x length, SOIL stresses with force/length^2.  A hard-coded 9.81 (or a module using another
   gravity than the HOUSE gravity or SOIL ``grav``) anywhere in the chain breaks this.
2. **Phase convention through to the time domain.**  With the time factor exp(+i w t) a wave that
   reaches a point later carries exp(-i w tau).  (a) Wave passage of HOUSE (FILE77 factors
   exp(-i w x/V_app)) applied as free-field motion (ANALYSX <ffm> = 1) and (b) an inclined SV wave of
   SITE (k = w sin(theta)/V_hs) must both *delay* the MOTION acceleration history of a node at
   distance x by exactly x/V_app (x sin(theta)/V_hs) -- checked sample by sample against the history
   of the reference node with integer-sample delays, so no interpolation is involved.
3. **In-structure response spectra** of MOTION against an independent frequency-domain oscillator
   (exact SDOF transfer function on the zero-padded MOTION history), SA(100 Hz) = ZPA = PGA.
4. **Translation invariance**: an embedded FV model under an inclined wave moved with its ground elevation
   and control point gives the same FILE8 (gelev and the control point are used consistently).
5. **EDU-06** (finding F-03): CHECK reports the free drilling rotations of shells on excavated SOLIDs, which
   FIXROT does not treat (manual 9.7) and which made ANALYS stop on a singular system after a clean CHECK.
"""
from __future__ import annotations

import math
import re
import shutil
from pathlib import Path

import numpy as np
import pytest

from sassi.io import textfiles
from sassi.io.container import read_container
from sassi.prep import Interpreter
from sassi.verify import builders as B

EXAMPLES = Path(__file__).resolve().parents[2] / "examples"
DYNP = Path(__file__).resolve().parents[2] / "sassi" / "data" / "dynp_library.pre"
FT = 0.3048                      # m per ft (exact)
KIP = 4.4482216152605            # kN per kip (exact)
G_SI = 9.80665


# ======================================================================================
# 1. SI vs British units through the whole chain
# ======================================================================================
def _units_model(units: str) -> str:
    """Two-layer site, 3 x 3 surface SHELL mat, two-mass BEAMS stick; SOIL (3 iterations) on the same
    profile.  Every dimensional input is converted with the unit factors of ``units``."""
    if units == "SI":
        L, F = 1.0, 1.0
    else:
        L, F = 1.0 / FT, 1.0 / KIP
    g = G_SI * L
    w, s = F / L ** 3, F / L ** 2          # specific weight, stress / modulus

    def n(x):
        return f"{x:.15g}"
    lines = [
        "MDL,u,u", "TIT,units audit",
        f"L,1,{n(1.0 * L)},{n(19.0 * w)},{n(600 * L)},{n(300 * L)},0.05,0.05",
        f"L,2,{n(1.0 * L)},{n(20.0 * w)},{n(1000 * L)},{n(500 * L)},0.04,0.04",
        f"L,3,{n(1.0 * L)},{n(21.0 * w)},{n(2000 * L)},{n(1000 * L)},0.02,0.02",
        "TOPL,1,1,1,1", "TOPL,2,2,2,2",
        "FREQ,1,4,20,41,82,123,164,205,246,328,410",
        "SITE,0,1,0,20,3,1,0,1,4096,1,0,0.005,8192,1", "WAVE,2,1,1,1,0",
        f"N,1,{n(-5 * L)},{n(-5 * L)},0", f"N,3,{n(5 * L)},{n(-5 * L)},0", "FILL,1,3",
        f"NGEN,2,3,1,3,1,0,{n(5 * L)},0",
        f"N,10,0,0,{n(5 * L)}", f"N,11,0,0,{n(10 * L)}",
        f"N,12,{n(10 * L)},0,{n(2.5 * L)}", f"N,13,{n(10 * L)},0,{n(7.5 * L)}",
        f"M,1,{n(3.0e8 * s)},0.2,{n(24.0 * w)},0.05,0.05,1",
        f"M,2,{n(3.0e7 * s)},0.2,0.0,0.05,0.05,1",
        f"R,1,{n(20.0 * L ** 2)},{n(10.0 * L ** 2)},{n(10.0 * L ** 2)},{n(400.0 * L ** 4)},{n(200.0 * L ** 4)},"
        f"{n(200.0 * L ** 4)}",
        "GROUP,1,SHELL", "MACT,1", "E,1,1,2,5,4", "EGEN,1,1,1", "EGEN,1,3,1,2", f"THICK,1,4,1,{n(1.5 * L)}",
        "GROUP,2,BEAMS", "MACT,2", "RACT,1", "E,1,5,10,12", "E,2,10,11,13",
        "FIXROT",
        f"MT,10,{n(9810 * F)},{n(9810 * F)},{n(9810 * F)}",
        f"MT,11,{n(9810 * F)},{n(9810 * F)},{n(9810 * F)}",
        "INT,1,9,1,1",
        f"POINT,0,0,{n(4.5 * L)}",
        f"HOUSE,{n(g)},0,0,2,0,0,0,0,0",
        "ANALYS,0,0,0,0,1,0,0,0,0,0,0",
        "MOTION,0,0,0,20,0,0.1,100,101,1,0,1,0,0,0,0,1,0,0,1",
        "DAMP,0.05",
        "THFILE,../motion.acc",
        "NOUT,1,1,1,0,0,1,1,5,10,11",
        "NOUT,3,1,1,0,0,1,1,4,6",
        "STRESS,0,0,1,1,1",
        "EOUT,1,2,1,1,1,2,1,1,1,1,1,1,2,1-2",
        "RELD,1,0,0", "RELFILE,00005TR_X.TFI", "RDND,10,1,0,0,0,0,0", "RDND,11,1,0,0,0,0,0",
        # SOIL on the same profile (sublayers 1-8, 9 = half-space), rock outcrop input, 3 iterations
        f"INP,{DYNP}",
        "SPRO,1,1,Sand", "SPRO,2,1,Sand", "SPRO,3,1,Sand", "SPRO,4,1,Sand",
        "SPRO,5,2,Clay", "SPRO,6,2,Clay", "SPRO,7,2,Clay", "SPRO,8,2,Clay", "SPRO,9,3",
        f"SOIL,4000,{n(g)},1,1,0,3,0.65,1,0", "SOILX,0,1,0,9",
        "SACC,1,2,0", "SSTR,2,1,1,1,1",
        "AOPT,0,1,0,1,1,1,0,0,1,0,1,1,1,0",
        "CHECK", "AFWRITE",
        "RUNSOIL", "RUNSITE", "RUNPOINT", "RUNHOUSE", "RUNANALYS", "RUNMOTION", "RUNSTRESS", "RUNRELDISP",
        # Option A: equivalent static loads at the 2 largest base-shear peaks and the dynamic tables
        "LOADGEN,3,0,0,0,1,1,FILE8", "LGTIME,V,2", "RUNLOADGEN,STATIC",
        "LOADGENDYN,0,0,REL,0,FILE8", "RUNLOADGEN,DYNAMIC",
    ]
    return "\n".join(lines) + "\n"


def _run_units(root: Path, units: str) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(EXAMPLES / "data" / "rg160h_030g.acc", root / "motion.acc")
    ui = Interpreter(cwd=root)
    pre = root / f"{units}.pre"
    pre.write_text(_units_model(units), encoding="utf-8")
    summary = ui.run_file(str(pre), resolve=False)
    assert summary.errors == 0, (root / "u" / "u.err").read_text() if (root / "u" / "u.err").exists() else summary
    return root / "u"


def _history(p: Path) -> np.ndarray:
    v, _ = textfiles.read_history(p)
    return np.asarray(v, float)


@pytest.fixture(scope="module")
def unit_runs(tmp_path_factory):
    root = tmp_path_factory.mktemp("audit_units")
    return _run_units(root / "SI", "SI"), _run_units(root / "BS", "BS")


def _rel(a, b) -> float:
    a, b = np.asarray(a), np.asarray(b)
    return float(np.abs(a - b).max() / np.abs(b).max())


def test_units_transfer_functions_are_dimensionless(unit_runs):
    si, bs = unit_runs
    a, b = read_container(si / "FILE8", "FILE8"), read_container(bs / "FILE8", "FILE8")
    assert np.array_equal(np.asarray(a["eq_node"]), np.asarray(b["eq_node"]))
    assert _rel(b["freq"], a["freq"]) < 1e-12
    dof = np.asarray(a["eq_dof"])
    Ha, Hb = np.asarray(a["H"]), np.asarray(b["H"])
    for d in range(1, 6):
        m = dof == d
        if m.any():
            scale = 1.0 if d <= 3 else FT                     # rotation per unit control displacement
            assert _rel(Hb[:, m], scale * Ha[:, m]) < 1e-9, d


def test_units_histories_spectra_and_forces_scale_exactly(unit_runs):
    si, bs = unit_runs
    names = sorted(p.name for p in si.iterdir())
    checked = {".ACC": 0, ".RS": 0, ".THD": 0, ".THS": 0, ".TFD": 0, ".TH": 0}
    for name in names:
        ext = Path(name).suffix.upper()
        p1, p2 = si / name, bs / name
        if ext not in checked:
            continue
        assert p2.exists(), name
        if ext == ".RS":
            a, b = textfiles.read_xy(p1), textfiles.read_xy(p2)
            assert _rel(b[:, 1:], a[:, 1:]) < 1e-9, name             # SA in g
        elif ext == ".TFD":
            f1, h1, _ = textfiles.read_tf(p1)
            f2, h2, _ = textfiles.read_tf(p2)
            assert _rel(h2, np.asarray(h1) / FT) < 1e-9, name          # length per g
        else:
            a, b = _history(p1), _history(p2)
            if ext in (".ACC",) or name.startswith(("ACC", "SN")):
                scale = 1.0                                           # g, strain in %
            elif ext == ".THD":
                scale = 1.0 / FT                                      # relative displacement
            elif name.startswith("SS"):
                scale = (1.0 / KIP) * FT ** 2                         # SOIL stress, force / length^2
            else:
                comp = name.rsplit("_", 1)[-1].split(".")[0]          # beam force FYI / moment MZI
                scale = (1.0 / KIP) if comp.startswith("F") else (1.0 / KIP / FT)
            assert np.abs(a).max() > 0, name
            assert _rel(b, scale * a) < 1e-7, name
        checked[ext] += 1
    assert all(v > 0 for v in checked.values()), checked


def test_isrs_equals_an_independent_frequency_domain_oscillator(unit_runs):
    """In-structure response spectra of MOTION (Nigam-Jennings with sub-stepping, SA in g) vs an
    independent computation from the MOTION acceleration history: the absolute acceleration of each
    5 % oscillator from its exact transfer function (w0^2 + 2 i z w0 w)/(w0^2 - w^2 + 2 i z w0 w) applied
    to the zero-padded record.  Observed: <= 0.47 % over 0.1-100 Hz; SA(100 Hz) = ZPA = PGA."""
    si, _ = unit_runs
    for name in ("00011TR_X", "00005TR_X", "00004TR_Z"):
        acc, dt = textfiles.read_history(si / f"{name}.ACC")
        acc = np.asarray(acc, float)
        rs = textfiles.read_xy(si / f"{name}01.RS")
        nfft = 1 << int(np.ceil(np.log2(16 * len(acc))))
        A = np.fft.rfft(acc, nfft)
        w = 2 * np.pi * np.fft.rfftfreq(nfft, dt)
        ref = []
        for f0 in rs[:, 0]:
            w0 = 2 * np.pi * f0
            H = (w0 ** 2 + 0.1j * w0 * w) / (w0 ** 2 - w ** 2 + 0.1j * w0 * w)
            ref.append(np.abs(np.fft.irfft(H * A, nfft)).max())
        assert np.max(np.abs(rs[:, 1] - ref) / ref) < 0.01, name
        assert abs(rs[-1, 1] - np.abs(acc).max()) / np.abs(acc).max() < 0.01, name     # ZPA at 100 Hz


def _apdl_values(path: Path):
    """F and D values {(node, label): value} and table values {name: array} of a LOADGEN APDL file."""
    loads, tables = {}, {}
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.split("!")[0].strip()
        m = re.match(r"^(F|D),(\d+),(\w+),([-+0-9.Ee]+)$", line)
        if m:
            loads[(m.group(1), int(m.group(2)), m.group(3))] = float(m.group(4))
            continue
        m = re.match(r"^(LG\w*?)\((\d+),(\d+)\)=(.*)$", line)
        if m:
            tables.setdefault(m.group(1), []).extend(float(x) for x in m.group(4).split(","))
    return loads, {k: np.asarray(v) for k, v in tables.items()}


def test_units_loadgen_apdl_scales_exactly(unit_runs):
    si, bs = unit_runs
    for name in ("u_LGS.inp", "u_LGD.inp"):
        (la, ta), (lb, tb) = _apdl_values(si / name), _apdl_values(bs / name)
        assert set(la) == set(lb) and set(ta) == set(tb)
        assert la or ta, name
        top = {k: max(abs(v) for key, v in la.items() if key[0] == k) for k in {key[0] for key in la}}
        for key, v in la.items():
            scale = 1.0 / KIP if key[0] == "F" else 1.0 / FT             # forces / interface displacements
            # values that are zero in exact arithmetic (e.g. UY under X input) carry FFT round-off of
            # ~1e-12 of the field maximum, which depends on the BLAS summation order: absolute floor
            # 1e-10 of the maximum (the APDL file prints 12 significant digits)
            assert abs(lb[key] - scale * v) <= 1e-9 * abs(scale * v) + 1e-10 * scale * top[key[0]], key
        peak = max((float(np.abs(v).max()) for v in ta.values()), default=0.0) / FT
        for key, v in ta.items():                                       # ACEL (length/s^2), D tables (length)
            assert float(np.abs(tb[key] - v / FT).max()) <= 1e-9 * max(float(np.abs(v).max()) / FT, 1e-3 * peak), key
    loads, _ = _apdl_values(si / "u_LGS.inp")
    assert sum(1 for k in loads if k[0] == "F") > 0 and sum(1 for k in loads if k[0] == "D") > 0


# ======================================================================================
# 2. Delays through to the time domain
# ======================================================================================
DT, NFT = 0.005, 2048


def _isolated_surface_nodes(site: B.Site, xs):
    """Surface interaction nodes with no structure: each node hangs on a 1e-6 kN/m spring to a fixed
    point, so the SSI solution is the free field (C = X_ff) to ~1e-12."""
    hb = B.HouseBuilder(site)
    nodes = [hb.add_node(x, 0.0, 0.0, fix=(0, 0, 0, 1, 1, 1)) for x in xs]
    hb.set_interaction(nodes)
    for nd, x in zip(nodes, xs):
        g = hb.add_node(x, 0.0, 1.0, fix=(1,) * 6)
        hb.spring(nd, g, (1e-6, 1e-6, 1e-6, 0.0, 0.0, 0.0))
    return B.Model(hb, {"nodes": nodes}, rad=2.0, layer=0), nodes


def _pulse_motion(wd: Path, fs: B.FrequencySet, nodes):
    t = np.arange(0.0, 6.0, DT)
    acc = np.exp(-((t - 1.0) / 0.05) ** 2) * np.cos(2 * np.pi * 5 * (t - 1.0))
    B.write_history(wd / "pulse.acc", acc, DT)
    md = B.motion_deck(fs, thfile="pulse.acc", nout=[(n, 1, 1, 1, 0, 0, 0, 1) for n in nodes], damp=[0.05],
                       dur=8.0, interp=1)
    B.write_deck(wd, "m", md)
    B.run("MOTION", wd, "m")
    return [_history(wd / f"{n:05d}TR_X.ACC") for n in nodes]


def _assert_delayed(hist, delays_s):
    base = hist[0]
    ref = np.abs(base).max()
    for a, tau in zip(hist[1:], delays_s[1:]):
        lag = int(round(tau / DT))
        assert lag > 0 and abs(lag * DT - tau) < 1e-12
        n = len(a)
        delayed = np.abs(a[lag:] - base[:n - lag]).max() / ref
        advanced = np.abs(a[:n - lag] - base[lag:]).max() / ref
        assert delayed < 1e-9 and advanced > 0.5, (tau, delayed, advanced)


def test_wave_passage_delays_the_motion_downstream(tmp_path):
    """HOUSE WPASS (coherent, V_app = 400 m/s along x) -> FILE77 -> ANALYS FFM -> MOTION."""
    site = B.layered_site([(4.0, 300.0, 600.0, 2.0, 0.05), (4.0, 400.0, 800.0, 2.0, 0.05)],
                          (500.0, 1000.0, 2.1, 0.02))
    fs = B.FrequencySet.fourier(DT, NFT, range(1, int(15.0 * DT * NFT) + 1))   # every Fourier line to 15 Hz
    xs, V = [0.0, 40.0, 100.0], 400.0
    mdl, nodes = _isolated_surface_nodes(site, xs)
    B.run_soil(tmp_path, "m", site, fs, layer=0, rad=mdl.rad)
    d = mdl.deck("m")
    d["wpass"], d["appv"], d["wang"], d["coh"] = 1, V, 0.0, 0
    B.write_deck(tmp_path, "m", fs.fill(d))
    B.run("HOUSE", tmp_path, "m")
    B.write_deck(tmp_path, "m", B.analys_deck(fs, wpass=1, ffm=1))
    B.run("ANALYS", tmp_path, "m")
    f8 = B.read_file8(tmp_path)
    for nd, x in zip(nodes, xs):                               # frequency domain: exp(-i w x / V)
        ref = np.exp(-2j * np.pi * np.asarray(f8["freq"]) * x / V)
        assert np.abs(B.tf(f8, nd, 1) - ref).max() < 1e-9
    _assert_delayed(_pulse_motion(tmp_path, fs, nodes), [x / V for x in xs])


def test_inclined_sv_arrives_later_downstream(tmp_path):
    """SITE inclined SV (30 deg in an undamped uniform site, V_hs = 1000 m/s) -> ANALYS -> MOTION."""
    vs, theta = 1000.0, 30.0
    site = B.layered_site([(5.0, vs, 2 * vs, 2.0, 0.0), (5.0, vs, 2 * vs, 2.0, 0.0)], (vs, 2 * vs, 2.0, 0.0))
    fs = B.FrequencySet.fourier(DT, NFT, range(1, int(15.0 * DT * NFT) + 1))
    xs = [0.0, 40.0, 100.0]
    mdl, nodes = _isolated_surface_nodes(site, xs)
    B.run_soil(tmp_path, "m", site, fs, layer=0, rad=mdl.rad, waves=[(2, 1, 1.0, 1.0, theta)], cm=0)
    B.run_house(tmp_path, "m", mdl)
    B.write_deck(tmp_path, "m", B.analys_deck(fs))
    B.run("ANALYS", tmp_path, "m")
    sin = math.sin(math.radians(theta))
    f8 = B.read_file8(tmp_path)
    for nd, x in zip(nodes, xs):
        ref = np.exp(-2j * np.pi * np.asarray(f8["freq"]) * x * sin / vs)
        assert np.abs(B.tf(f8, nd, 1) - ref).max() < 1e-9
    _assert_delayed(_pulse_motion(tmp_path, fs, nodes), [x * sin / vs for x in xs])


# ======================================================================================
# 3. Translation invariance (gelev and control point handled consistently by every module)
# ======================================================================================
def _shifted_model(dx: float, dy: float, dz: float, drilling: bool = True) -> str:
    """Embedded 2 x 2 x 1 excavation (FV) with a SHELL mat on its top face and a two-mass BEAMS stick,
    under an inclined SV wave (30 deg, k != 0); the whole model is translated by (dx, dy, dz), the ground
    elevation by dz and the ANALYS control point by (dx, dy)."""
    def n(v):
        return f"{v:.12g}"
    nodes = []
    k = 0
    for iz, z in enumerate((-1.0, 0.0)):                 # bottom level first (bottom-up numbering)
        for iy in range(3):
            for ix in range(3):
                k += 1
                nodes.append(f"N,{k},{n(-2 + 2 * ix + dx)},{n(-2 + 2 * iy + dy)},{n(z + dz)}")
    lines = [
        "MDL,t,t",
        "L,1,0.5,19.0,600,300,0.05,0.05", "L,2,1.0,21.0,2000,1000,0.02,0.02",
        "TOPL,1,1,1,1",
        "FREQ,1,20,41,82,123,164,205",
        "SITE,0,1,0,20,2,1,0,1,4096,1,0,0.005,8192,1", "WAVE,2,1,1,1,30",
        *nodes,
        f"N,20,{n(dx)},{n(dy)},{n(4 + dz)}", f"N,21,{n(dx)},{n(dy)},{n(8 + dz)}",
        f"N,22,{n(5 + dx)},{n(dy)},{n(2 + dz)}", f"N,23,{n(5 + dx)},{n(dy)},{n(6 + dz)}",
        "M,1,3.0E8,0.2,24.0,0.05,0.05,1", "M,2,3.0E7,0.2,0.0,0.05,0.05,1",
        "R,1,20.0,10.0,10.0,400.0,200.0,200.0",
        "GROUP,1,SOLID", "MACT,1", "E,1,1,2,5,4,10,11,14,13", "EGEN,1,1,1", "EGEN,1,3,1,2", "ETYPE,1,4,1,2",
        "GROUP,2,SHELL", "MACT,1", "E,1,10,11,14,13", "EGEN,1,1,1", "EGEN,1,3,1,2", "THICK,1,4,1,1.0",
        "GROUP,3,BEAMS", "MACT,2", "RACT,1", "E,1,14,20,22", "E,2,20,21,23",
        "FIXROT",
        # the mat nodes are shared with the excavated SOLIDs: FIXROT does not treat them (shell-only nodes),
        # so their drilling rotation is fixed with D (node 14 carries the stick)
        *(["D,10,13,1,1,ROTZ", "D,15,18,1,1,ROTZ"] if drilling else []),
        "MT,20,9810,9810,9810", "MT,21,9810,9810,9810",
        "INT,1,18,1,1",
        "POINT,0,2,1.8",
        f"HOUSE,9.81,{n(dz)},0,2,0,0,0,0,0",
        f"ANALYS,0,0,0,0,1,0,0,{n(dx)},{n(dy)},0,0",
        "AOPT,0,0,0,1,1,1,0,0,1,0,0,0,0,0",
        "CHECK", "AFWRITE", "RUNSITE", "RUNPOINT", "RUNHOUSE", "RUNANALYS",
    ]
    return "\n".join(lines) + "\n"


def test_translation_of_model_ground_and_control_point_leaves_the_solution_unchanged(tmp_path):
    out = {}
    for tag, shift in (("ref", (0.0, 0.0, 0.0)), ("moved", (100.0, -50.0, 7.0))):
        root = tmp_path / tag
        root.mkdir()
        pre = root / "t.pre"
        pre.write_text(_shifted_model(*shift), encoding="utf-8")
        ui = Interpreter(cwd=root)
        summary = ui.run_file(str(pre), resolve=False)
        assert summary.errors == 0, (root / "t" / "t.err").read_text()
        out[tag] = read_container(root / "t" / "FILE8", "FILE8")
    a, b = out["ref"], out["moved"]
    assert np.array_equal(np.asarray(a["eq_node"]), np.asarray(b["eq_node"]))
    Ha, Hb = np.asarray(a["H"]), np.asarray(b["H"])
    assert np.abs(Ha).max() > 0.1
    assert _rel(Hb, Ha) < 1e-9


def test_check_reports_free_drilling_rotations_of_shells_on_solids(tmp_path):
    """EDU-06 in CHECK covers shell nodes shared with SOLID elements (final audit): before, only HOUSE
    reported them and ANALYS stopped on a singular system after a clean CHECK."""
    from sassi.prep.check import run_check
    ui = Interpreter(cwd=tmp_path)
    text = _shifted_model(0.0, 0.0, 0.0, drilling=False)
    ui.run_text(text[:text.index("CHECK")])
    rep = run_check(ui.model)
    flagged = sorted(int(m.text.split()[-1]) for m in rep.messages if m.number == "EDU-06")
    assert flagged == [10, 11, 12, 13, 15, 16, 17, 18], rep.format()
    assert all("with D" in m.detail for m in rep.messages if m.number == "EDU-06")
    ui.run_text("D,10,13,1,1,ROTZ\nD,15,18,1,1,ROTZ\n")
    assert not any(m.number == "EDU-06" for m in run_check(ui.model).messages)
