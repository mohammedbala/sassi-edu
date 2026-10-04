"""Verification problems of the STRESS module (requirements section 4.10 and 6.3).

* VP-35  stress recovery (spec 05d section 8): a rigid-body nodal TF gives a zero STF (SOLID, BEAMS,
         SHELL, SPRING in 3D, PLANE in 2D; translations, and rotations for the elements that carry
         them); pure shear gives ``tau_oct = (sqrt 6/3) tau = 0.8165 tau`` and uniaxial stress
         ``(sqrt 2/3) sigma``; a cantilever BEAMS element under a quasi-static (low-frequency) tip
         load shape has end forces in equilibrium (``F_I + F_J = 0``, ``M_I + M_J + L e1 x F_J = 0``)
         and ``M3_I = k_L u = -P2 L``.  Tolerance 1e-10.
* VP-S1  a fixed-base cantilever stick (BEAMS, tip mass) excited through a stiff site by a harmonic
         control motion well below its natural frequency: the STRESS base moment equals the
         quasi-static ``m a h`` within 1 %, and the ``.THS`` history equals ``irfft(STF U_g)`` of the
         single harmonic (closed form) within 1e-10.
* VP-S2  a SOLID column in a uniform layer under vertically propagating SV: the STRESS shear stress at
         each depth equals ``G* gamma_ff`` of the exact continuous solution within 2 % (peak values and
         STFs up to 15 Hz).

VP-S1 and VP-S2 run two variants with the same checks:

* a self-contained one that isolates STRESS -- VP-S1 solves the flexible-volume equation of R1
  section 4.4, ``[(K*_s - w^2 M_s) - (K*_e - w^2 M_e) + X] U = X U'_f``, from the HOUSE matrices
  (COOSK/COOSM) with a frequency-independent stiff impedance ``X = k_site I`` at the interaction node
  (rigid-site limit); VP-S2 takes FILE8 = the SITE free field (FILE1) at the column nodes (a column
  with the free-field properties does not change the motion: the zero-SSI identity of VP-16);
* the whole chain SITE -> POINT -> HOUSE -> ANALYS -> STRESS built with :mod:`sassi.verify.builders`
  (VP-S1: a massless stick with a tip mass on a rigid surface mat of a stiff site; VP-S2: an FV
  excavation whose structure is the excavated soil), when ANALYS and the builders are available.

VP-35 synthesises FILE8 from closed-form nodal fields.  Every problem writes HOUSE and STRESS decks
with :mod:`sassi.io.decks` and runs the modules through :func:`sassi.modules.base.run_module`, so the
chain deck -> HOUSE (FILE4 recovery operators) -> STRESS (STF, interpolation, convolution, files) is
exercised.  The builders are public so that the unit tests reuse them.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, Sequence, Tuple

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from ...conventions import element_result_name
from ...core import signal as S
from ...core import stress_lib as SL
from ...core.freefield import free_field_at_nodes
from ...elements import MaterialProps, material_from_M
from ...elements.beam import transformation
from ...io import decks, textfiles
from ...io.files import read_container
from ...modules.base import run_module
from .. import VPResult, problem, worse
from .vp_house import G_SI, add_element, add_node, house_deck, run_house
from .vp_motion import synthetic_motion, write_file8, write_motion_file

# =======================================================================================
# Builders (also used by tests/unit/test_stress*.py)
# =======================================================================================
MAT_STEEL = (1, 1, 2.0e8, 0.3, 78.0, 0.04, 0.04)         # M row: id, type 1 (E, nu), weight, damping
MAT_CONC = (2, 1, 3.0e7, 0.2, 24.0, 0.05, 0.05)
MAT_PD = (3, 3, 900.0, 400.0, 20.0, 0.06, 0.03)           # type 3 (Vp, Vs), beta_p != beta_s
SECTION = (1, 0.25, 0.2, 0.18, 0.012, 0.004, 0.006)       # R row: A As2 As3 J I2 I3
SPRING = (1, 1.0e5, 2.0e5, 3.0e5, 10.0, 20.0, 30.0, 0.03)


def material(row) -> MaterialProps:
    """MaterialProps of an M-table row as HOUSE builds it (gravity G_SI)."""
    _, t, v1, v2, w, pd, sd = row
    return material_from_M(t, v1, v2, w, pd, sd, G_SI)


def stress_deck(nft: int, delt: float, thfile: str = "eq.acc", gravity: float = G_SI, **params) -> decks.Deck:
    """STRESS deck with the control-motion data and ``params`` set (``eout`` rows via :func:`add_eout`)."""
    d = decks.new("STRESS")
    d["model"] = "m"
    d["thfile"] = thfile
    d["nft"] = int(nft)
    d["delt"] = float(delt)
    d["df"] = 1.0 / (nft * delt)
    d["gravity"] = gravity
    d["mult"], d["max"] = 1.0, 0.0
    for k, v in params.items():
        if k not in d.params:
            raise KeyError(f"unknown STRESS parameter {k}")
        d[k] = v
    return d


def add_eout(d: decks.Deck, group: int, elems: Iterable[int], codes: Sequence[int]) -> None:
    codes = (list(codes) + [0] * 12)[:12]
    for e in elems:
        d.table("eout").append([int(group), int(e)] + [int(c) for c in codes])


def run_stress(workdir: Path, d: decks.Deck, model: str = "m") -> Tuple[int, str]:
    decks.write(Path(workdir) / f"{model}.str", d)
    rc = run_module("STRESS", model, workdir)
    return rc, (Path(workdir) / f"{model}_STRESS.out").read_text(encoding="utf-8")


def read_tfu(workdir: Path, etype: str, group: int, elem: int, comp: str, ext: str = "TFU"):
    f, H, _ = textfiles.read_tf(Path(workdir) / element_result_name(etype, group, elem, comp, ext))
    return f, H


def read_ths(workdir: Path, etype: str, group: int, elem: int, comp: str):
    return textfiles.read_history(Path(workdir) / element_result_name(etype, group, elem, comp, "THS"))


def file4_maps(workdir: Path, model: str = "m"):
    """(FILE4 container, {(node, dof): eq}, {node: xyz})."""
    f4 = read_container(Path(workdir) / f"{model}.N4", "FILE4")
    eq = {(int(n), int(k)): i for i, (n, k) in enumerate(zip(f4["eq_node"], f4["eq_dof"]))}
    xyz = {int(n): np.asarray(p, float) for n, p in zip(f4["node_id"], f4["node_xyz"])}
    return f4, eq, xyz


def field_file8(workdir: Path, fnum, df: float, field, nfft: int, delt: float, model: str = "m",
                **kw) -> np.ndarray:
    """FILE8 whose TF at every equation is ``field(node, xyz, f) -> (6,)`` complex (UX..ROTZ)."""
    f4, eq, xyz = file4_maps(workdir, model)
    fnum = np.asarray(list(fnum), dtype=np.int64)
    H = np.zeros((len(fnum), len(eq)), dtype=complex)
    for q, n in enumerate(fnum):
        for (node, k), i in eq.items():
            H[q, i] = field(node, xyz[node], n * df)[k - 1]
    kw.setdefault("model_hash", str(f4.meta.get("model_hash", "")))
    write_file8(Path(workdir) / "FILE8", fnum, df, f4["eq_node"], f4["eq_dof"], H, nfft=nfft, delt=delt, **kw)
    return H


def mixed_model_deck() -> decks.Deck:
    """3D model for VP-35: a distorted SOLID, an inclined BEAMS element with a K node, an inclined flat
    SHELL quad (not parallel to any global plane) and a SPRING, connected at shared nodes."""
    d = house_deck()
    pts = [(0.0, 0.0, 0.0), (1.1, 0.1, 0.05), (1.2, 1.0, -0.1), (-0.1, 0.9, 0.0),
           (0.05, -0.1, 1.0), (1.0, 0.0, 1.1), (1.15, 1.1, 0.95), (0.0, 1.05, 1.05)]
    for i, p in enumerate(pts):
        add_node(d, i + 1, *p)
    add_node(d, 9, 1.8, 2.5, 3.2)                       # beam J
    add_node(d, 10, 6.0, -1.0, 0.5, fix=(1,) * 6)       # K node (orientation only)
    # shell: a flat parallelogram through node 6 in an inclined plane
    o, a, b = np.array(pts[5]), np.array([1.2, 0.3, 0.4]), np.array([-0.2, 1.1, 0.5])
    for n, p in ((11, o + a), (12, o + a + b), (13, o + b)):
        add_node(d, n, *p)
    for r in (MAT_STEEL, MAT_CONC, MAT_PD):
        d.table("materials").append(list(r))
    d.table("beamprops").append(list(SECTION))
    d.table("springprops").append(list(SPRING))
    for gid, typ, title in ((1, 1, "block"), (2, 2, "column"), (3, 3, "slab"), (4, 7, "spring")):
        d.table("groups").append([gid, typ, title])
    add_element(d, 1, 1, list(range(1, 9)), etype=1, mat=3)
    add_element(d, 2, 1, [7, 9, 10], mat=1)
    add_element(d, 3, 1, [6, 11, 12, 13], mat=2, thick=0.25)
    add_element(d, 4, 1, [9, 12], prop=1)
    return d


def plane_model_deck() -> decks.Deck:
    """2D (X-Z) model for VP-35: two distorted PLANE quadrilaterals."""
    d = house_deck(dim=1)
    pts = {1: (0.0, 0.0), 2: (1.0, 0.1), 3: (2.1, 0.0), 4: (0.1, 1.0), 5: (1.1, 1.2), 6: (2.0, 0.9)}
    for n, (x, z) in pts.items():
        add_node(d, n, x, 0.0, z)
    d.table("materials").append(list(MAT_PD))
    d.table("groups").append([1, 4, "plane"])
    add_element(d, 1, 1, [1, 2, 5, 4], etype=1, mat=3)
    add_element(d, 1, 2, [2, 3, 6, 5], etype=1, mat=3)
    return d


def cantilever_deck(L: float = 4.0) -> Tuple[decks.Deck, np.ndarray]:
    """One BEAMS element, I fixed (all six DOFs), inclined, K node off the axis; returns (deck, xyz I/J/K)."""
    d = house_deck()
    xi = np.array([0.5, -0.2, 0.3])
    e1 = np.array([2.0, 1.0, 2.0]) / 3.0
    xj = xi + L * e1
    xk = np.array([3.0, -2.0, 0.0])
    add_node(d, 1, *xi, fix=(1,) * 6)
    add_node(d, 2, *xj)
    add_node(d, 3, *xk, fix=(1,) * 6)
    d.table("materials").append(list(MAT_STEEL))
    d.table("beamprops").append(list(SECTION))
    d.table("groups").append([1, 2, "cantilever"])
    add_element(d, 1, 1, [1, 2, 3], mat=1)
    return d, np.vstack([xi, xj, xk])


def harmonic_file(path: Path, nfft: int, dt: float, k0: int, a0: float) -> np.ndarray:
    """Control motion ``a0 sin(2 pi k0 n / nfft)`` over the whole Fourier period (a single rfft bin)."""
    a = a0 * np.sin(2.0 * np.pi * k0 * np.arange(nfft) / nfft)
    write_motion_file(path, a, dt)
    return a


def _rel_zero(stf: np.ndarray, scale: np.ndarray) -> float:
    m = scale > 0
    return float(np.max(np.abs(stf[m]) / scale[m])) if np.any(m) else 0.0


def _stf_scales(f4, H: np.ndarray) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """Per type: (STF from the FILE4 operators, sum_d |S_cd U_d|) at every frequency -- the size of the
    terms that must cancel for a rigid-body motion."""
    out = {}
    for t in f4.meta["components"]:
        Sx, eqx = f4[f"rec_{t}_S"], f4[f"rec_{t}_eq"]
        U = SL.element_dof_tf(H, eqx)
        out[t] = (np.einsum("ecd,fed->fec", Sx, U), np.einsum("ecd,fed->fec", np.abs(Sx), np.abs(U)))
    return out


# =======================================================================================
# VP-35
# =======================================================================================
NFFT35, DT35 = 512, 0.01
DF35 = 1.0 / (NFFT35 * DT35)
FNUM35 = np.array([1, 2, 3, 5, 8, 12, 17, 23, 30, 40])


def _vp35_rigid(r: VPResult, wd: Path, deck, label: str, motions) -> None:
    rc, out = run_house(wd, deck)
    r.require(f"{label}: HOUSE run succeeds", rc == 0, out[-500:] if rc else "")
    if rc:
        return
    write_motion_file(wd / "eq.acc", 0.2 * synthetic_motion(300, DT35, seed=5, fmax=20.0), DT35)
    f4, _, _ = file4_maps(wd)
    for name, field, types in motions:
        H = field_file8(wd, FNUM35, DF35, field, NFFT35, DT35)
        d = stress_deck(NFFT35, DT35, itran=1)
        groups = sorted(set(int(g) for g in f4["elem_group"]))
        for g in groups:
            ids = sorted(int(e) for gg, e in zip(f4["elem_group"], f4["elem_id"]) if gg == g)
            add_eout(d, g, ids, [2] * 12)
        rc, out = run_stress(wd, d)
        r.require(f"{label}, {name}: STRESS run succeeds", rc == 0, out[-800:] if rc else "")
        if rc:
            continue
        sc = _stf_scales(f4, H)
        ug = np.max(np.abs(np.fft.irfft(S.displacement_spectrum(
            np.fft.rfft(S.pad_record(textfiles.read_history(wd / "eq.acc")[0], NFFT35)),
            S.fourier_grid(NFFT35, DT35), G_SI), NFFT35)))
        for t in types:
            if t not in sc:
                continue
            stf, scale = sc[t]
            comps = f4.meta["components"][t]
            err_tf, err_th = 0.0, 0.0
            for e, gidx in enumerate(f4[f"rec_{t}_idx"]):
                gi, ei = int(f4["elem_group"][gidx]), int(f4["elem_id"][gidx])
                for c, comp in enumerate(comps):
                    v, _ = read_ths(wd, t, gi, ei, comp)
                    err_th = worse(err_th, float(np.max(np.abs(v))) / max(float(np.max(scale)) * ug, 1e-300))
                    if comp == "SOCT":                     # derived in the time domain: no TF file
                        continue
                    _, h = read_tfu(wd, t, gi, ei, comp)   # STF written by STRESS
                    err_tf = worse(err_tf, _rel_zero(np.abs(h), scale[:, e, c]))
            r.check(f"{label}, {name}: {t} max |STF| / sum_d |S_cd U_d| (.TFU)", err_tf, 0.0, atol=1e-10)
            r.check(f"{label}, {name}: {t} max |history| / (max sum_d |S_cd U_d| x max |u_g(t)|) (.THS)",
                    err_th, 0.0, atol=1e-10)


def _rigid_translation(v):
    v = np.asarray(v, dtype=complex)
    return lambda node, x, f: np.concatenate([v, np.zeros(3)])


def _rigid_rotation(theta, x0):
    th = np.asarray(theta, dtype=float)
    return lambda node, x, f: np.concatenate([np.cross(th, x - x0), th]).astype(complex)


def _vp35_states(r: VPResult, wd: Path) -> None:
    """Pure shear and uniaxial stress fields imposed on a distorted SOLID (linear displacement fields are
    reproduced exactly by the trilinear element: patch test)."""
    d = mixed_model_deck()
    rc, out = run_house(wd, d)
    r.require("states: HOUSE run succeeds", rc == 0, out[-500:] if rc else "")
    if rc:
        return
    mat = material(MAT_PD)
    G, E, nu = mat.G, mat.E, mat.nu
    write_motion_file(wd / "eq.acc", 0.25 * synthetic_motion(300, DT35, seed=7, fmax=20.0), DT35)
    gam, eps = 1.0e-4, 2.0e-4
    states = {
        "pure shear tau_xy": (lambda n, x, f: np.array([0.5 * gam * x[1], 0.5 * gam * x[0], 0, 0, 0, 0], complex),
                              {"SXY": G * gam}, SL.OCT_STRESS_PURE_SHEAR, "SXY", 0.8165),
        "uniaxial sigma_xx": (lambda n, x, f: np.array([eps * x[0], -nu * eps * x[1], -nu * eps * x[2], 0, 0, 0],
                                                       complex),
                              {"SXX": E * eps}, SL.OCT_STRESS_UNIAXIAL, "SXX", np.sqrt(2.0) / 3.0),
    }
    comps = ["SXX", "SYY", "SZZ", "SXY", "SXZ", "SYZ"]
    for name, (field, nonzero, ratio, main, quoted) in states.items():
        field_file8(wd, FNUM35, DF35, field, NFFT35, DT35)
        dk = stress_deck(NFFT35, DT35, itran=1)
        add_eout(dk, 1, [1], [2] * 7)
        rc, out = run_stress(wd, dk)
        r.require(f"{name}: STRESS run succeeds", rc == 0, out[-800:] if rc else "")
        if rc:
            continue
        ref_main = abs(nonzero[main])
        for c in comps:
            _, h = read_tfu(wd, "SOLID", 1, 1, c)
            if c in nonzero:
                r.check(f"{name}: STF {c} = {'G* gamma' if c == 'SXY' else 'E* eps'} (complex, all f)",
                        float(np.max(np.abs(h - nonzero[c]))) / ref_main, 0.0, atol=1e-10)
            else:
                r.check(f"{name}: STF {c} = 0", float(np.max(np.abs(h))) / ref_main, 0.0, atol=1e-10)
        s_main, _ = read_ths(wd, "SOLID", 1, 1, main)
        soct, _ = read_ths(wd, "SOLID", 1, 1, "SOCT")
        r.check(f"{name}: max SOCT / max |{main}| (reference {quoted:.4f})", np.max(soct) / np.max(np.abs(s_main)),
                ratio, rtol=1e-10)
        r.check(f"{name}: SOCT(t) = {ratio:.6f} |{main}(t)| at every step",
                float(np.max(np.abs(soct - ratio * np.abs(s_main))) / np.max(np.abs(s_main))), 0.0, atol=1e-10)


def _vp35_cantilever(r: VPResult, wd: Path) -> None:
    """Quasi-static tip-load shape of a cantilever (I fixed): recovered end forces."""
    L = 4.0
    d, X = cantilever_deck(L)
    rc, out = run_house(wd, d)
    r.require("cantilever: HOUSE run succeeds", rc == 0, out[-500:] if rc else "")
    if rc:
        return
    mat = material(MAT_STEEL)
    E, G = mat.beam_moduli()
    _, A, As2, As3, J, I2, I3 = SECTION
    P1, P2, P3, T = 30.0, 12.0, -7.0, 2.5
    uL = np.zeros(6, dtype=complex)
    uL[0] = P1 * L / (E * A)
    uL[1] = P2 * L ** 3 / (3 * E * I3) + P2 * L / (G * As2)
    uL[2] = P3 * L ** 3 / (3 * E * I2) + P3 * L / (G * As3)
    uL[3] = T * L / (G * J)
    uL[4] = -P3 * L ** 2 / (2 * E * I2)                 # t2 = -du3/dx
    uL[5] = P2 * L ** 2 / (2 * E * I3)                  # t3 = +du2/dx
    _, Lam, _ = transformation(X)
    uJ = np.concatenate([Lam.T @ uL[:3], Lam.T @ uL[3:]])
    nfft, dt = 1024, 0.02
    df = 1.0 / (nfft * dt)
    field_file8(wd, np.arange(1, 13), df, lambda n, x, f: uJ if n == 2 else np.zeros(6, complex), nfft, dt)
    write_motion_file(wd / "eq.acc", 0.3 * synthetic_motion(400, dt, seed=3, fmax=0.6, f_hp=0.05), dt)
    dk = stress_deck(nfft, dt, itran=1)
    add_eout(dk, 1, [1], [2] * 12)
    rc, out = run_stress(wd, dk)
    r.require("cantilever: STRESS run succeeds", rc == 0, out[-800:] if rc else "")
    if rc:
        return
    names = ["FXI", "FYI", "FZI", "MXI", "MYI", "MZI", "FXJ", "FYJ", "FZJ", "MXJ", "MYJ", "MZJ"]
    stf = {c: read_tfu(wd, "BEAMS", 1, 1, c)[1] for c in names}
    th = {c: read_ths(wd, "BEAMS", 1, 1, c)[0] for c in names}
    Fs = P2 * L
    exact = {"FXJ": P1, "FYJ": P2, "FZJ": P3, "MXJ": T, "MYJ": 0.0, "MZJ": 0.0}
    for c, v in exact.items():
        r.check(f"cantilever STF {c} = {v:g} (tip load, all f)", float(np.max(np.abs(stf[c] - v))) / Fs, 0.0,
                atol=1e-10)
    for a, b in (("FXI", "FXJ"), ("FYI", "FYJ"), ("FZI", "FZJ"), ("MXI", "MXJ")):
        r.check(f"cantilever equilibrium {a} + {b} = 0 (STF)", float(np.max(np.abs(stf[a] + stf[b]))) / Fs, 0.0,
                atol=1e-10)
    r.check("cantilever M3_I = k_L u = -P2 L (STF, all f)", float(np.max(np.abs(stf["MZI"] + P2 * L))) / Fs, 0.0,
            atol=1e-10)
    r.check("cantilever M2_I = +P3 L (moment equilibrium about axis 2)",
            float(np.max(np.abs(stf["MYI"] - P3 * L))) / Fs, 0.0, atol=1e-10)
    hmax = max(float(np.max(np.abs(v))) for v in th.values())
    r.check("cantilever histories: max |FYI + FYJ| / max |r|", float(np.max(np.abs(th["FYI"] + th["FYJ"]))) / hmax,
            0.0, atol=1e-10)
    r.check("cantilever histories: max |MZI + L FYJ + MZJ| / max |r| (moment equilibrium in time)",
            float(np.max(np.abs(th["MZI"] + L * th["FYJ"] + th["MZJ"]))) / hmax, 0.0, atol=1e-10)


@problem("VP-35", "STRESS recovery: rigid body, octahedral stress, cantilever end forces", tier="P0",
         modules=["HOUSE", "STRESS"], source="requirements 6.3; spec 05d section 8")
def vp35(workdir) -> VPResult:
    r = VPResult()
    wd = Path(workdir)
    x0 = np.array([0.3, 0.4, 0.5])
    _vp35_rigid(r, wd / "rigid3d", mixed_model_deck(), "3D mixed model", [
        ("rigid translation (1, -0.6, 0.35)", _rigid_translation((1.0, -0.6, 0.35)), ("SOLID", "BEAMS", "SHELL", "SPRING")),
        ("rigid rotation (0.02, -0.03, 0.05) rad", _rigid_rotation((0.02, -0.03, 0.05), x0), ("SOLID", "BEAMS", "SHELL")),
    ])
    _vp35_rigid(r, wd / "rigid2d", plane_model_deck(), "2D PLANE model", [
        ("rigid translation (1, 0, -0.7)", _rigid_translation((1.0, 0.0, -0.7)), ("PLANE",)),
        ("rigid rotation 0.03 rad about Y", _rigid_rotation((0.0, 0.03, 0.0), x0), ("PLANE",)),
    ])
    _vp35_states(r, wd / "states")
    _vp35_cantilever(r, wd / "cantilever")
    r.notes.append("Rigid-body errors are |S.U| relative to the sum of the magnitudes of the cancelling terms "
                   "sum_d |S_cd U_d| (FILE4 operators, FILE8 TFs); SPRING forces are not invariant under rigid "
                   "rotation (F = k (u_J - u_I) per global component, D-STR-05) and are checked for translation only.")
    r.notes.append("FILE8 is synthetic (closed-form nodal fields); the STFs are read back from the STRESS .TFU files "
                   "and the octahedral relations from the .THS histories.")
    return r


# =======================================================================================
# VP-S1: cantilever stick through a stiff site, harmonic control motion
# =======================================================================================
STICK1 = dict(h=10.0, ne=10, E=3.0e7, nu=0.2, beta=0.05, A=0.5, J=0.08, I2=0.04, I3=0.05, fn=5.0, beam_mass=0.01)


def stick_site_deck(p=STICK1) -> Tuple[decks.Deck, float, float]:
    """Vertical cantilever stick (z up) of BEAMS elements on one interaction node at the surface
    (translations free, rotations fixed: a rigid foundation), tip mass m; K node on +X so that local
    axis 2 = X (deflection along X bends about local 3 = Y with I3).  Returns (deck, m, rho)."""
    h, ne = p["h"], p["ne"]
    k = 3.0 * p["E"] * p["I3"] / h ** 3
    m = k / (2.0 * np.pi * p["fn"]) ** 2
    rho = p["beam_mass"] * m / (p["A"] * h)
    d = house_deck(gravity=G_SI)
    for i in range(ne + 1):
        add_node(d, i + 1, 0.0, 0.0, h * i / ne, fix=(0, 0, 0, 1, 1, 1) if i == 0 else (0,) * 6)
    add_node(d, 99, 1.0, 0.0, 0.0, fix=(1,) * 6)
    d.table("interaction").append([1])
    d.table("materials").append([1, 1, p["E"], p["nu"], rho * G_SI, p["beta"], p["beta"]])
    d.table("beamprops").append([1, p["A"], 0.0, 0.0, p["J"], p["I2"], p["I3"]])
    d.table("groups").append([1, 2, "stick"])
    for i in range(ne):
        add_element(d, 1, i + 1, [i + 1, i + 2, 99])
    d.table("masses").append([ne + 1, m, m, m, 0.0, 0.0, 0.0, 0])
    return d, m, rho


def stiff_site_file8(workdir: Path, fnum, df: float, k_site: float, nfft: int, delt: float, cm: int = 0,
                     model: str = "m") -> np.ndarray:
    """FILE8 of the flexible-volume equation (R1 section 4.4) with a rigid-site impedance.

    ``[(K*_s - w^2 M_s) - (K*_e - w^2 M_e) + X] U = X U'_f`` with ``X = k_site`` on the interaction
    translations and the free field ``U'_f`` = unit motion in the control direction ``cm`` at every
    interaction node (vertical incidence on a rigid site).  Complex-symmetric system: sparse LU."""
    wd = Path(workdir)
    f4 = read_container(wd / f"{model}.N4", "FILE4")
    ck, cmm = read_container(wd / "COOSK", "COOSK"), read_container(wd / "COOSM", "COOSM")
    Ks, Ke, Ms, Me = ck.sparse("Ks"), ck.sparse("Ke"), cmm.sparse("Ms"), cmm.sparse("Me")
    neq = Ks.shape[0]
    ie = np.asarray(f4["int_eq"]).reshape(-1, 3)
    xd = np.zeros(neq)
    xd[ie[ie >= 0]] = k_site
    X = sp.diags(xd)
    uf = np.zeros(neq, dtype=complex)
    uf[ie[:, cm][ie[:, cm] >= 0]] = 1.0
    fnum = np.asarray(list(fnum), dtype=np.int64)
    H = np.zeros((len(fnum), neq), dtype=complex)
    for q, n in enumerate(fnum):
        w2 = (2.0 * np.pi * n * df) ** 2
        Amat = (Ks - w2 * Ms - (Ke - w2 * Me) + X).tocsc()
        H[q] = spla.spsolve(Amat, X @ uf)
    write_file8(wd / "FILE8", fnum, df, f4["eq_node"], f4["eq_dof"], H, nfft=nfft, delt=delt, cm=cm,
                model_hash=str(f4.meta.get("model_hash", "")))
    return H


NFFT_S1, DT_S1, K0_S1, A0_S1 = 1024, 0.04, 10, 0.2
FNUM_S1 = np.array([1, 2, 4, 6, 8, 10, 12, 16, 20, 30, 40, 60, 80, 100, 140, 180, 205, 220, 260, 300, 400, 512])


def _chain_modules():
    """The SITE -> POINT -> HOUSE -> ANALYS chain helpers (sassi.verify.builders) or None when this
    build has no ANALYS module / builders (the self-contained variants then run alone)."""
    try:
        from .. import builders as B
        import sassi.modules.analys  # noqa: F401
    except ImportError:
        return None
    return B


def _s1_checks(r: VPResult, wd: Path, label: str, H: np.ndarray, base: int, group: int, elem: int,
               m_qs: float, fn: float) -> None:
    """Checks of VP-S1 on one FILE8 (``H``) of a stick whose base element is (``group``, ``elem``)."""
    df = 1.0 / (NFFT_S1 * DT_S1)
    f4, eq, _ = file4_maps(wd)
    q0 = int(np.where(FNUM_S1 == K0_S1)[0][0])
    r.check(f"{label}: stiff site, |H_base,x - 1| at the excitation frequency", abs(H[q0, eq[(base, 1)]] - 1.0), 0.0,
            atol=1e-3)
    a = harmonic_file(wd / "eq.acc", NFFT_S1, DT_S1, K0_S1, A0_S1)
    dk = stress_deck(NFFT_S1, DT_S1, itran=1)
    add_eout(dk, group, [elem], [2] * 12)
    rc, out = run_stress(wd, dk)
    r.require(f"{label}: STRESS run succeeds", rc == 0, out[-800:] if rc else "")
    if rc:
        return
    mzi, dt = read_ths(wd, "BEAMS", group, elem, "MZI")
    r.check(f"{label}: THS time step", dt, DT_S1, rtol=1e-12)
    # (1) quasi-static base moment m a h
    f0 = K0_S1 * df
    r.check(f"{label}: max |base moment MZI| = m a h quasi-static (f0/fn = {f0 / fn:.3f})", np.max(np.abs(mzi)),
            m_qs, rtol=0.01)
    # (2) the history is irfft(STF U_g) of the single harmonic: closed form
    names = list(f4.meta["components"]["BEAMS"])
    row = [k for k, g in enumerate(f4["rec_BEAMS_idx"])
           if int(f4["elem_group"][g]) == group and int(f4["elem_id"][g]) == elem][0]
    S1, eq1 = f4["rec_BEAMS_S"][row], f4["rec_BEAMS_eq"][row]
    U = np.where(eq1 >= 0, H[q0, np.clip(eq1, 0, None)], 0.0)
    # S annihilates the rigid translation r (unit X at every element DOF): S.U = S.(U - r); the second
    # form avoids cancelling the O(1) rigid-body terms of the total-motion TF at this low frequency
    rig = SL.rigid_body_dofs((1, 2, 3, 4, 5, 6), 2, cm=0)
    c = names.index("MZI")
    stf0 = (S1 @ (U - rig))[c]
    plain = (S1 @ U)[c]
    cancel = float(np.sum(np.abs(S1[c] * U)) / abs(stf0))
    w0 = 2.0 * np.pi * f0
    t = np.arange(NFFT_S1) * DT_S1
    ref = (G_SI * A0_S1 / w0 ** 2) * np.real(1j * stf0 * np.exp(1j * w0 * t))
    r.check(f"{label}: .THS MZI = irfft(STF U_g) (closed form of the harmonic), max error / max", float(
        np.max(np.abs(mzi - ref)) / np.max(np.abs(ref))), 0.0, atol=1e-10)
    Ug = S.displacement_spectrum(np.fft.rfft(a), S.fourier_grid(NFFT_S1, DT_S1), G_SI)
    R = np.zeros(NFFT_S1 // 2 + 1, dtype=complex)
    R[K0_S1] = stf0 * Ug[K0_S1]
    r.check(f"{label}: .THS MZI = irfft(STF U_g) at the excitation bin, max error / max",
            float(np.max(np.abs(mzi - np.fft.irfft(R, NFFT_S1))) / np.max(np.abs(ref))), 0.0, atol=1e-10)
    _, stf_file = read_tfu(wd, "BEAMS", group, elem, "MZI")
    r.check(f"{label}: .TFU MZI at f0 = S . U_e(f0) from FILE4 and FILE8", abs(stf_file[q0] - stf0) / abs(stf0), 0.0,
            atol=1e-12)
    lst = (wd / "m_STRESS.out").read_text(encoding="utf-8")
    r.require(f"{label}: listing reports the maximum of MZI with its time", "MZI" in lst and "time (s)" in lst)
    r.notes.append(f"{label}: quasi-static M = {m_qs:.6g}, computed {np.max(np.abs(mzi)):.6g} (ratio "
                   f"{np.max(np.abs(mzi)) / m_qs:.5f}: dynamic amplification ~ 1/(1 - (f0/fn)^2) and the site "
                   f"flexibility); S.U_e at f0 sums terms {cancel:.3g} times larger than the result, the plain product "
                   f"differs from S.(U_e - r) by {abs(plain - stf0) / abs(stf0):.2e} (STRESS subtracts the rigid-body "
                   "motion before the recovery, sassi.core.stress_lib.element_stf)")


def _s1_chain(r: VPResult, wd: Path, B) -> None:
    """VP-S1 through SITE -> POINT -> HOUSE -> ANALYS (sassi.verify.builders): massless BEAMS stick with a
    tip mass on a rigid 4 m x 4 m surface mat of a stiff uniform site (Vs = 1000 m/s)."""
    p = STICK1
    h, ne = p["h"], p["ne"]
    k = 3.0 * p["E"] * p["I3"] / h ** 3
    m = k / (2.0 * np.pi * p["fn"]) ** 2
    site = B.uniform_site(depth=20.0, nsub=4, vs=1000.0, rho=2.0, nu=1.0 / 3.0, beta=0.02, gravity=G_SI)
    fs = B.FrequencySet.fourier(DT_S1, NFFT_S1, FNUM_S1)
    stick = B.Stick(heights=[h * (i + 1) / ne for i in range(ne)], masses=[0.0] * (ne - 1) + [m], E=p["E"],
                    A=p["A"], I=p["I3"], nu=p["nu"], beta=p["beta"], J=p["J"])
    mdl = B.stick_on_mat(site, stick, half_width=2.0, ndiv=2)
    wd.mkdir(parents=True, exist_ok=True)
    try:
        B.run_soil(wd, "m", site, fs, layer=mdl.layer, rad=mdl.rad)
        f4 = B.run_house(wd, "m", mdl)
        B.run_analys(wd, "m", fs, gravity=G_SI)
    except B.ChainError as exc:
        r.require("chain: SITE, POINT, HOUSE and ANALYS runs succeed", False, str(exc)[-600:])
        return
    base, first = mdl["centre"], mdl["stick"][0]
    hit = [k for k in range(len(f4["elem_id"])) if int(f4["elem_type"][k]) == 2
           and list(f4["elem_nodes"][k][:2]) == [base, first]]
    r.require("chain: base element of the stick found in FILE4", len(hit) == 1)
    if len(hit) != 1:
        return
    g, e = int(f4["elem_group"][hit[0]]), int(f4["elem_id"][hit[0]])
    H = np.asarray(B.read_file8(wd)["H"])
    _s1_checks(r, wd, "chain", H, base, g, e, A0_S1 * G_SI * m * h, p["fn"])


@problem("VP-S1", "STRESS base moment of a cantilever stick on a stiff site under a harmonic motion", tier="P0",
         modules=["HOUSE", "STRESS", "SITE", "POINT", "ANALYS"], source="requirements 4.10, D-CNV-06, D-STR-02, D-STR-05")
def vp_s1(workdir) -> VPResult:
    r = VPResult()
    wd = Path(workdir)
    p = STICK1
    # ---- (a) flexible-volume equation with a stiff impedance at the interaction node ---------------------
    wa = wd / "impedance"
    d, m, rho = stick_site_deck(p)
    rc, out = run_house(wa, d)
    r.require("HOUSE run of the stick succeeds", rc == 0, out[-600:] if rc else "")
    if rc:
        return r
    k_struct = 3.0 * p["E"] * p["I3"] / p["h"] ** 3
    df = 1.0 / (NFFT_S1 * DT_S1)
    H = stiff_site_file8(wa, FNUM_S1, df, 1.0e5 * k_struct, NFFT_S1, DT_S1)
    m_qs = A0_S1 * G_SI * (m * p["h"] + rho * p["A"] * p["h"] ** 2 / 2.0)    # tip mass + light beam (q h^2/2)
    _s1_checks(r, wa, "stiff impedance", H, 1, 1, 1, m_qs, p["fn"])
    r.notes.append(f"stick: h = {p['h']} m, k = 3EI/h^3 = {k_struct:.6g}, tip mass m = {m:.6g}, fn ~ {p['fn']} Hz, "
                   f"harmonic f0 = {K0_S1 * df:.4f} Hz, a0 = {A0_S1} g")
    r.notes.append("variant 'stiff impedance': the flexible-volume equation of R1 4.4 solved from COOSK/COOSM with a "
                   "frequency-independent stiff impedance (k_site = 1e5 k) at the interaction node (rigid-site limit)")
    # ---- (b) the whole chain SITE -> POINT -> HOUSE -> ANALYS -> STRESS ---------------------------------
    B = _chain_modules()
    if B is None:
        r.notes.append("variant 'chain' skipped: sassi.verify.builders / ANALYS are not available in this build")
    else:
        _s1_chain(r, wd / "chain", B)
    return r


# =======================================================================================
# VP-S2: SOLID column in a uniform layer under vertical SV (SITE free field)
# =======================================================================================
SITE_S2 = dict(nlay=30, h=1.0, rho=2.0, vs=200.0, vp=400.0, beta=0.05, hs=(0.0, 2.2, 1000.0, 2000.0, 0.01, 0.01))
NCOL_S2 = 10
NFFT_S2, DT_S2 = 2048, 0.01
FNUM_S2 = np.unique(np.concatenate([[2, 4, 6, 8, 10, 13], np.arange(16, 420, 7)])).astype(np.int64)
F_CHECK_S2 = 15.0          # Hz: element height h <= lambda/13 below this frequency


def column_deck(p=SITE_S2, ncol: int = NCOL_S2) -> decks.Deck:
    """1 m x 1 m SOLID column from the surface down ``ncol`` user layers; material = the soil layer
    (M type 3: Vp, Vs, weight, damping) so the column does not disturb the free field."""
    d = house_deck(gravity=G_SI)
    n = 0
    for lev in range(ncol + 1):
        z = -lev * p["h"]
        for (x, y) in ((0, 0), (1, 0), (1, 1), (0, 1)):
            n += 1
            add_node(d, n, x, y, z)
    d.table("materials").append([1, 3, p["vp"], p["vs"], p["rho"] * G_SI, p["beta"], p["beta"]])
    d.table("groups").append([1, 1, "soil column"])
    for e in range(ncol):
        bot = 4 * (e + 1)             # level e+1 is below level e
        top = 4 * e
        add_element(d, 1, e + 1, [bot + 1, bot + 2, bot + 3, bot + 4, top + 1, top + 2, top + 3, top + 4], etype=1,
                    mat=1)
    return d


def _exact_shear_tf(f: np.ndarray, depth: float, p=SITE_S2) -> np.ndarray:
    """``tau_xz / u_surface = G* k* sin(k* d)``: exact continuous SV field in the uniform top layer with the
    control point at the surface (``u(d) = cos(k* d)``, ``gamma_xz = du/dz = k* sin(k* d)`` with z up)."""
    from ...conventions import cfactor
    Gs = p["rho"] * p["vs"] ** 2 * cfactor(p["beta"])
    vss = p["vs"] * np.sqrt(cfactor(p["beta"]))
    k = 2.0 * np.pi * np.asarray(f, dtype=float) / vss
    return Gs * k * np.sin(k * depth)


def _s2_checks(r: VPResult, wd: Path, label: str, group: int, acc: np.ndarray) -> Tuple[float, float]:
    """STRESS run on the column elements of ``group`` and the comparison with ``G* gamma_ff``; returns the
    worst STF and peak errors."""
    p = SITE_S2
    df = 1.0 / (NFFT_S2 * DT_S2)
    f4, _, xyz = file4_maps(wd)
    elems = [(int(f4["elem_id"][k]), float(np.mean([xyz[int(n)][2] for n in f4["elem_nodes"][k] if n > 0])))
             for k in range(len(f4["elem_id"])) if int(f4["elem_group"][k]) == group]
    elems.sort()
    dk = stress_deck(NFFT_S2, DT_S2, itran=1)
    add_eout(dk, group, [e for e, _ in elems], [1, 0, 0, 0, 2, 0, 1])
    rc, out = run_stress(wd, dk)
    r.require(f"{label}: STRESS run succeeds", rc == 0, out[-800:] if rc else "")
    if rc:
        return 1.0, 1.0
    f_grid = S.fourier_grid(NFFT_S2, DT_S2)
    fN = FNUM_S2[-1] * df
    Ug = S.displacement_spectrum(np.fft.rfft(S.pad_record(acc, NFFT_S2)), f_grid, G_SI)
    f_ssi = FNUM_S2 * df
    sel = f_ssi <= F_CHECK_S2
    depths = np.array([-zc for _, zc in elems])                                  # gelev = 0
    ref_cols = np.array([_exact_shear_tf(f_ssi, dc) for dc in depths])           # (ncol, nF)
    col_max = np.max(np.abs(ref_cols), axis=0)
    worst_f, worst_t = 0.0, 0.0
    for (e, _), dc, ref in zip(elems, depths, ref_cols):
        _, stf = read_tfu(wd, "SOLID", group, e, "SXZ")
        err_f = float(np.max(np.abs(stf[sel] - ref[sel]) / col_max[sel]))
        worst_f = worse(worst_f, err_f)
        r.check(f"{label}: STF SXZ = G* gamma_ff at depth {dc:.1f} m, f <= {F_CHECK_S2:g} Hz (error / column max)",
                err_f, 0.0, atol=0.02)
        tau, _ = read_ths(wd, "SOLID", group, e, "SXZ")
        R = np.where(f_grid <= fN * (1 + 1e-12), _exact_shear_tf(f_grid, dc), 0.0) * Ug
        tref = np.fft.irfft(R, NFFT_S2)
        r.check(f"{label}: peak |tau_xz| at depth {dc:.1f} m vs G* gamma_ff (exact continuous)", np.max(np.abs(tau)),
                np.max(np.abs(tref)), rtol=0.02)
        worst_t = worse(worst_t, abs(np.max(np.abs(tau)) / np.max(np.abs(tref)) - 1.0))
    lst = (wd / "m_STRESS.out").read_text(encoding="utf-8")
    r.require(f"{label}: listing has the maxima of SXX, SXZ and SOCT of the column",
              all(c in lst for c in ("SXX", "SXZ", "SOCT")))
    return worst_f, worst_t


def _s2_chain(r: VPResult, wd: Path, B, acc: np.ndarray) -> None:
    """VP-S2 through SITE -> POINT -> HOUSE -> ANALYS: a 1 m x 1 m FV excavation of the top NCOL_S2 layers whose
    structure is the excavated soil itself (sassi.verify.builders.embedded_box, structure 'soil')."""
    p = SITE_S2
    site = B.layered_site([(p["h"], p["vs"], p["vp"], p["rho"], p["beta"], p["beta"])] * p["nlay"],
                          (p["hs"][2], p["hs"][3], p["hs"][1], p["hs"][4], p["hs"][5]), gravity=G_SI)
    fs = B.FrequencySet.fourier(DT_S2, NFFT_S2, FNUM_S2)
    mdl = B.embedded_box(site, 0.5, NCOL_S2, ndiv=1, method="FV", structure="soil")
    wd.mkdir(parents=True, exist_ok=True)
    try:
        B.run_soil(wd, "m", site, fs, layer=mdl.layer, rad=mdl.rad)
        f4 = B.run_house(wd, "m", mdl)
        B.run_analys(wd, "m", fs, gravity=G_SI)
    except B.ChainError as exc:
        r.require("chain: SITE, POINT, HOUSE and ANALYS runs succeed", False, str(exc)[-600:])
        return
    groups = sorted(set(int(g) for g, x in zip(f4["elem_group"], f4["elem_excav"]) if int(x) == 0))
    r.require("chain: one structural SOLID group (the soil column)", len(groups) == 1)
    if len(groups) != 1:
        return
    write_motion_file(wd / "eq.acc", acc, DT_S2)
    wf, wt = _s2_checks(r, wd, "chain", groups[0], acc)
    r.notes.append(f"variant 'chain' (SITE, POINT, HOUSE, ANALYS FV with structure = excavated soil): worst STF error "
                   f"{100 * wf:.3f} %, worst peak error {100 * wt:.3f} %")


@problem("VP-S2", "STRESS shear stress in a SOLID soil column under vertical SV = G* x free-field strain",
         tier="P0", modules=["SITE", "POINT", "HOUSE", "ANALYS", "STRESS"],
         source="requirements 4.10, R1 section 2 / 4.4 (VP-16 identity)")
def vp_s2(workdir) -> VPResult:
    from .vp_site_point import run_site
    r = VPResult()
    wd = Path(workdir)
    p = SITE_S2
    df = 1.0 / (NFFT_S2 * DT_S2)
    acc = 0.2 * synthetic_motion(1500, DT_S2, seed=21, fmax=12.5)
    # ---- (a) SITE free field at the column nodes (the zero-SSI identity, no ANALYS) ---------------------
    wa = wd / "freefield"
    layers = [(p["h"], p["rho"], p["vs"], p["vp"], p["beta"], p["beta"])] * p["nlay"]
    files = run_site(wa / "site", "col", layers, p["hs"], FNUM_S2, df, nl=20, cl=1, cm=0)
    f1 = files["FILE1"]
    rc, out = run_house(wa, column_deck())
    r.require("HOUSE run of the column succeeds", rc == 0, out[-600:] if rc else "")
    if rc:
        return r
    f4, eq, xyz = file4_maps(wa)
    nodes = sorted(xyz)
    iface = np.array([int(round(-xyz[n][2] / p["h"])) + 1 for n in nodes])
    P = np.array([xyz[n] for n in nodes])
    H = np.zeros((len(FNUM_S2), len(eq)), dtype=complex)
    for q in range(len(FNUM_S2)):
        U = free_field_at_nodes(f1, q, P, iface, 0.0, 0.0, 0.0)
        for i, n in enumerate(nodes):
            for k in (1, 2, 3):
                if (n, k) in eq:
                    H[q, eq[(n, k)]] = U[i, k - 1]
    write_file8(wa / "FILE8", FNUM_S2, df, f4["eq_node"], f4["eq_dof"], H, nfft=NFFT_S2, delt=DT_S2,
                model_hash=str(f4.meta.get("model_hash", "")))
    write_motion_file(wa / "eq.acc", acc, DT_S2)
    wf, wt = _s2_checks(r, wa, "free field", 1, acc)
    r.notes.append(f"column h = {p['h']} m in a uniform layer Vs = {p['vs']} m/s, beta = {p['beta']}; variant 'free field' "
                   f"(FILE8 = SITE FILE1 at the column nodes: a column with the free-field properties does not change "
                   f"the motion, VP-16): worst STF error {100 * wf:.3f} % (f <= {F_CHECK_S2:g} Hz), worst peak error "
                   f"{100 * wt:.3f} %; the errors are the element's centroid finite difference (kh)^2/24 and the SITE "
                   "thin-layer discretisation")
    # ---- (b) the whole chain -----------------------------------------------------------------------------
    B = _chain_modules()
    if B is None:
        r.notes.append("variant 'chain' skipped: sassi.verify.builders / ANALYS are not available in this build")
    else:
        _s2_chain(r, wd / "chain", B, acc)
    return r
