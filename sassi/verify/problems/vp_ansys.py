"""Verification of the ANSYS interfaces: VP-A1 (requirements 3.4.K, 6.6; D-ANS-03, D-ANS-06).

VP-A1  A cantilever modelled in ANSYS with BEAM188 elements and a RECT section (B = 0.2 along the
       element y axis, H = 0.4 along z, so the two bending planes have different stiffness) is
       read from a ``.cdb`` file by CONVERT,ANSYS, once with an orientation node K and once with
       ANSYS's default orientation (no K node).  Checks:

       1. the converted models have the same first 8 fixed-base frequencies as the native
          SASSI-EDU model of the same cantilever (rtol 1e-8);
       2. the native model agrees with Euler-Bernoulli theory, ``f_n = (b_n L)^2/(2 pi L^2)
          sqrt(EI/(rho A))``, in both bending planes (1 %; the residual is Timoshenko shear);
       3. the beam-axis mapping itself (D-ANS-03: SASSI axis 2 = ANSYS z, I2 = Izz, I3 = Iyy).
          Checks 1 and 4 compare *sorted* frequency lists, which cannot see a swap of the two
          bending planes of a straight prismatic cantilever (the swapped model has the same
          spectrum).  So, for both converted files and both re-imported models: (a) the bending
          frequencies classified by mode-shape direction (deflection along global Y or Z) equal
          the native ones (1e-8), and (b) the static tip deflections under unit loads along global
          Y and Z equal ``P L^3/(3 E I) + P L/(G As)`` with ``I = Izz = H B^3/12`` (Y) and
          ``Iyy = B H^3/12`` (Z), the ANSYS section convention (1e-8);
       4. the native model exported with the ANSYS command (legacy BEAM4 and modern BEAM188/ASEC,
          i.e. the export mapping IZZ = I2, IYY = I3) and read back by CONVERT,ANSYS keeps its
          frequencies (1e-8): the export is the inverse of the import.

The modal analysis uses the HOUSE element library directly (:func:`fixed_base_frequencies`
builds the element records of a model the way HOUSE does: global coordinates, D fixities,
M/R/SC/MX tables, MT/MR masses in mass units) with the undamped stiffness.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from sassi.elements import (ElemRecord, assemble, build_dofmap, material_from_M, natural_frequencies)
from sassi.verify import VPResult, problem

DATA = Path(__file__).resolve().parents[2] / "data" / "ansys"      # sassi/data/ansys (package data)

#: cantilever of VP-A1 (SI units): length, elements, RECT B x H, steel
L_CANT, NE_CANT, B_RECT, H_RECT = 20.0, 20, 0.2, 0.4
E_STEEL, NU_STEEL, RHO_STEEL, G_SI = 2.0e11, 0.3, 7850.0, 9.81
#: Euler-Bernoulli cantilever eigenvalues b_n L
BETA_L = (1.875104068711961, 4.694091132974175, 7.854757438237613)


# ======================================================================================
# helpers (also used by the unit tests)
# ======================================================================================
def element_records(model, undamped: bool = True):
    """(node_xyz, ids, fix (n, 6), records, masses) of the *structure* of a model, as HOUSE builds
    them; excavated-soil SOLID/PLANE elements (ETYPE resolved by D-HOU-01) are left out."""
    from sassi.prep.check import ModelView
    excavated = {(r.group, r.elem.id) for r in ModelView(model).elems if r.type in (1, 4) and r.etype == 2}
    ids, xyz = model.global_coordinates()
    node_xyz = {int(i): xyz[k] for k, i in enumerate(ids)}
    fix = np.array([model.nodes[int(i)].fix for i in ids], dtype=int)
    g = model.gravity
    incompatible = int(model.mopt.get("incomp")) == 0
    mats = {}
    for k, mt in model.materials.items():
        mats[k] = material_from_M(mt.mtype, mt.val1, mt.val2, mt.weight, mt.pdamp, mt.sdamp, g)
    recs: List[ElemRecord] = []
    for grp, e in model.iter_elements():
        if (grp.id, e.id) in excavated:
            continue
        code = grp.type
        nodes = tuple(e.nodes)
        props: Dict = {}
        mat = mats.get(e.mat)
        if code == 1:
            props = dict(eint=e.eint, incompatible=incompatible)
        elif code == 4:
            props = dict(incompatible=incompatible)
        elif code == 3:
            props = dict(thick=e.thick)
        elif code == 2:
            s = model.sections[e.prop]
            props = dict(section=dict(A=s.axial, As2=s.shear2, As3=s.shear3, J=s.tors, I2=s.flex2, I3=s.flex3),
                         ki="".join(str(v) for v in e.ki), kj="".join(str(v) for v in e.kj))
        elif code == 7:
            sp = model.springs[e.prop]
            props = dict(k=sp.k, damp=sp.damp)
            mat = None
        elif code == 9:
            p = model.matrices[e.prop]
            props = dict(KR=p.full("R"), KI=p.full("I"), MM=p.full("M"), munits=int(model.mopt.get("matrix")),
                         gravity=g)
            mat = None
        recs.append(ElemRecord(grp.id, e.id, code, nodes, False, mat, props))
    masses = {n: v for n, v in model.nodal_masses().items()}
    return node_xyz, [int(i) for i in ids], fix, recs, masses


def assemble_model(model, undamped: bool = True):
    """(dofmap, AssembledModel) of the structure of ``model`` (undamped stiffness by default)."""
    node_xyz, ids, fix, recs, masses = element_records(model)
    dm = build_dofmap(ids, fix, recs, extra_dofs={n: range(1, 7) for n in masses})
    return dm, assemble(node_xyz, recs, dm, masses=masses, undamped=undamped)


def fixed_base_frequencies(model, nmodes: Optional[int] = None, return_modes: bool = False):
    """Undamped fixed-base natural frequencies (Hz) of a model (and the DOF map / mode shapes)."""
    dm, am = assemble_model(model)
    out = natural_frequencies(am.Ks, am.Ms, nmodes, return_modes=return_modes)
    return (out, dm) if return_modes else out


def native_cantilever_text(with_k: bool = True) -> str:
    """``.pre`` text of the native SASSI-EDU cantilever of VP-A1 (K node above: axis 2 = +Z)."""
    from sassi.model.materials import section_rectangle
    from sassi.model.values import fmt_num as f
    s = section_rectangle(B_RECT, H_RECT)              # b along axis 3 (= -Y), h along axis 2 (= +Z)
    lines = ["TIT,VP-A1 native cantilever", f"GRAVITY,{f(G_SI)}"]
    for i in range(NE_CANT + 1):
        lines.append(f"N,{i + 1},{f(L_CANT * i / NE_CANT)},0,0")
    lines.append(f"N,{NE_CANT + 2},0,0,1")
    lines.append(f"M,1,{f(E_STEEL)},{f(NU_STEEL)},{f(RHO_STEEL * G_SI)},0.02,0.02,1")
    lines.append("R,1," + ",".join(f(s[k]) for k in ("axial", "shear2", "shear3", "tors", "flex2", "flex3")))
    lines.append("GROUP,1,BEAMS")
    lines += [f"E,{i + 1},{i + 1},{i + 2},{NE_CANT + 2}" for i in range(NE_CANT)]
    lines.append("D,1,1,1,1,ALL")
    return "\n".join(lines) + "\n"


def native_cantilever():
    from sassi.prep import Interpreter, Kind
    ui = Interpreter()
    ui.run_text(native_cantilever_text())
    errs = ui.sink.texts(Kind.ERROR)
    if errs:
        raise RuntimeError("native cantilever: " + "; ".join(errs))
    return ui.model


def bending_modes(model, n_each: int = 3) -> Dict[str, List[float]]:
    """Frequencies of the first bending modes deflecting along global Y and Z (mode-shape test)."""
    (f, phi), dm = fixed_base_frequencies(model, return_modes=True)
    out: Dict[str, List[float]] = {"Y": [], "Z": []}
    eq_dof = np.asarray(dm.eq_dof)
    for k in range(len(f)):
        v = np.abs(phi[:, k])
        e = {d: float(np.sum(v[eq_dof == d] ** 2)) for d in (1, 2, 3)}
        tot = sum(e.values()) or 1.0
        for d, lab in ((2, "Y"), (3, "Z")):
            if e[d] / tot > 0.9 and len(out[lab]) < n_each:
                out[lab].append(float(f[k]))
    return out


def tip_deflection(model, tip: int, direction: int, P: float = 1.0) -> float:
    """Static tip deflection of a fixed-base model under a unit load (global direction 1..3)."""
    import scipy.sparse.linalg as spla
    from sassi.elements.base import blas_quiet
    dm, am = assemble_model(model)
    F = np.zeros(dm.neq)
    q = dm.eq(tip, direction)
    F[q] = P
    with blas_quiet():
        K = am.Ks.tocsc()
        K = type(K)((np.ascontiguousarray(K.data.real), K.indices, K.indptr), shape=K.shape)
        u = spla.spsolve(K, F)
    return float(u[q])


def converted(path: Path, gravity: float = G_SI):
    from sassi.io.ansys_cdb import convert_cdb, read_cdb
    return convert_cdb(read_cdb(path), gravity).model


def reimported(model, modern: bool = False):
    """``model`` -> APDL (ANSYS command) -> CONVERT,ANSYS reader -> model."""
    from sassi.io.ansys_cdb import convert_cdb, read_cdb_text
    from sassi.io.apdl import ApdlOptions, export_apdl
    text = export_apdl(model, ApdlOptions(modern=modern)).text
    return convert_cdb(read_cdb_text(text), model.gravity).model


# ======================================================================================
# VP-A1
# ======================================================================================
@problem("VP-A1", "ANSYS BEAM188 cantilever (.cdb) vs native model and Euler-Bernoulli", tier="P1",
         modules=["HOUSE", "CONVERT", "ANSYS"], source="requirements 6.6; D-ANS-03; R2 I.1 (beam theory)")
def vpa1(workdir):
    r = VPResult()
    missing = [f for f in ("beam188_cantilever.cdb", "beam188_default_orientation.cdb") if not (DATA / f).exists()]
    if missing:
        r.require("sample .cdb files present", False, note=f"missing {missing} in {DATA}")
        r.notes.append(f"VP-A1 needs the sample files {missing} of sassi/data/ansys")
        return r
    native = native_cantilever()
    f_nat = fixed_base_frequencies(native)
    nm = 8
    models = {"BEAM188 (K node)": converted(DATA / "beam188_cantilever.cdb"),
              "BEAM188 (no K node)": converted(DATA / "beam188_default_orientation.cdb"),
              "native -> APDL (legacy BEAM4)": reimported(native, False),
              "native -> APDL (modern BEAM188 ASEC)": reimported(native, True)}
    # ---- 1. converted .cdb models vs the native model ---------------------------------------
    for label in ("BEAM188 (K node)", "BEAM188 (no K node)"):
        f = fixed_base_frequencies(models[label])
        for k in range(nm):
            r.check(f"{label} -> SASSI: f{k + 1} vs native [Hz]", f[k], f_nat[k], rtol=1e-8)
    # ---- 2. native model vs Euler-Bernoulli in both bending planes ----------------------------
    A = B_RECT * H_RECT
    I_z_bend = B_RECT * H_RECT ** 3 / 12.0      # deflection along Z: ANSYS Iyy = SASSI I3
    I_y_bend = H_RECT * B_RECT ** 3 / 12.0      # deflection along Y: ANSYS Izz = SASSI I2
    modes_nat = bending_modes(native)
    for lab, I in (("Z", I_z_bend), ("Y", I_y_bend)):
        for k, bl in enumerate(BETA_L):
            ref = bl ** 2 / (2.0 * math.pi * L_CANT ** 2) * math.sqrt(E_STEEL * I / (RHO_STEEL * A))
            comp = modes_nat[lab][k] if k < len(modes_nat[lab]) else float("nan")
            r.check(f"native: bending along {lab}, mode {k + 1} vs Euler-Bernoulli [Hz]", comp, ref, rtol=0.01)
    # ---- 3. axis mapping: bending planes by mode shape, static tip deflections --------------
    G = E_STEEL / (2.0 * (1.0 + NU_STEEL))
    As = A / 1.2
    tip = NE_CANT + 1
    for label, m in models.items():
        modes = bending_modes(m)
        for lab in ("Y", "Z"):
            for k in range(len(BETA_L)):
                comp = modes[lab][k] if k < len(modes[lab]) else float("nan")
                ref = modes_nat[lab][k] if k < len(modes_nat[lab]) else float("nan")
                r.check(f"{label}: bending along {lab}, mode {k + 1} vs native [Hz]", comp, ref, rtol=1e-8)
        for d, lab, I in ((2, "Y", I_y_bend), (3, "Z", I_z_bend)):
            ref = L_CANT ** 3 / (3.0 * E_STEEL * I) + L_CANT / (G * As)
            r.check(f"{label}: tip deflection under P_{lab} = 1 [m]", tip_deflection(m, tip, d), ref, rtol=1e-8)
    # ---- 4. export (ANSYS command) -> import round trip ----------------------------------------
    for label in ("native -> APDL (legacy BEAM4)", "native -> APDL (modern BEAM188 ASEC)"):
        f = fixed_base_frequencies(models[label])
        for k in range(nm):
            r.check(f"{label} -> CONVERT: f{k + 1} [Hz]", f[k], f_nat[k], rtol=1e-8)
    r.notes.append(f"cantilever L = {L_CANT} m, {NE_CANT} elements, RECT B = {B_RECT} (y) x H = {H_RECT} (z), "
                   f"E = {E_STEEL:g}, nu = {NU_STEEL}, rho = {RHO_STEEL}; native f1..f4 = "
                   + ", ".join(f"{x:.6f}" for x in f_nat[:4]) + " Hz")
    r.notes.append("SASSI BEAMS are Timoshenko beams with As = 5A/6: the Euler-Bernoulli differences (largest for "
                   "mode 3, a few 0.1 %) are the shear deformation of the converted RECT section")
    r.notes.append("sorted frequency lists (checks 1 and 4) are blind to a swap of the bending planes of this "
                   "prismatic cantilever; the direction-classified frequencies and the static tip deflections "
                   "(check 3) verify the axis mapping of every converted and re-imported model")
    r.notes.append("not checked here (needs ANSYS): ANSYS BEAM188 itself; requirements 6.6 lists the ANSYS-side "
                   "cross-checks for the user")
    return r
