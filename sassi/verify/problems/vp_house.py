"""Verification problems of the HOUSE module: VP-H1 and VP-H2.

* VP-H1  HOUSE matrices of a small mixed model (SOLID + BEAMS + SHELL + SPRING + GENERAL + nodal
         masses) written as a HOUSE deck equal the element-library assembly of the same element records
         (identity 1e-12: it verifies the deck -> element-record translation of HOUSE -- materials, sections,
         releases, ETYPE, masses, DOF map -- not the element formulations, which VP-37/38/39 verify), the
         total masses equal hand computations, and the fixed-base natural frequencies of a cantilever stick
         written through a HOUSE deck match Euler-Bernoulli theory (1 %, R2 I.1).
* VP-H2  excavated soil: the total excavated mass equals sum(rho V) of the layers (1e-12), the
         excavated stiffness/mass of a soil block equal the SOLID element matrices with the layer
         properties (material_from_layer, D-ELM-12) scattered by hand (identity 1e-12), the complex moduli
         G* = rho Vs^2 c(beta_s) and M* = rho Vp^2 c(beta_p) are recovered from HOUSE's Ke with affine fields
         (closed form, 1e-10); the interaction nodes are placed on the SITE interfaces.

The decks are built with :func:`sassi.io.decks.new` exactly as AFWRITE writes them and the module
runs through :func:`sassi.modules.base.run_module`, so these problems exercise the deck -> FILE4 /
COOSK / COOSM path (requirements 4.4).  The builders are public so that the unit tests reuse them.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Sequence, Tuple

import numpy as np
import scipy.sparse as sp

from sassi.elements import (ElemRecord, assemble, build_dofmap, material_from_layer, material_from_M,
                            natural_frequencies)
from sassi.elements import solid
from sassi.elements.base import blas_quiet
from sassi.io import decks
from sassi.io.files import read_container
from sassi.modules.base import run_module
from sassi.verify import VPResult, problem

G_SI = 9.81
#: site profile used by the examples: (L no, thick, weight, vp, vs, dp, ds); last row = half-space
SITE_LAYERS = [(1, 2.0, 18.0, 400.0, 200.0, 0.05, 0.05),
               (2, 3.0, 19.0, 600.0, 300.0, 0.04, 0.04),
               (3, 0.0, 20.0, 1000.0, 500.0, 0.02, 0.02)]


# ======================================================================================
# deck builders (also used by tests/unit/test_house*.py)
# ======================================================================================
def write_site_deck(workdir: Path, model: str = "m", layers=SITE_LAYERS, gravity: float = G_SI, nl: int = 20,
                    cmodform: int = 0) -> Path:
    """Minimal SITE deck ``<model>.sit`` (HOUSE always requires it, D-HOU-02)."""
    d = decks.new("SITE")
    d["model"] = model
    d["gravity"] = gravity
    d["nl"] = nl
    d["cmodform"] = cmodform
    d["hs"] = int(layers[-1][0])
    for row in layers[:-1]:
        d.table("layers").append(list(row))
    d.table("halfspace").append(list(layers[-1]))
    return decks.write(Path(workdir) / f"{model}.sit", d)


def house_deck(gravity: float = G_SI, gelev: float = 0.0, layers=SITE_LAYERS, **params) -> decks.Deck:
    """Empty HOUSE deck with the site table (TOPL rows + half-space row of thickness 0) filled."""
    d = decks.new("HOUSE")
    d["model"] = "m"
    d["gravity"] = gravity
    d["gelev"] = gelev
    for k, v in params.items():
        d[k] = v
    for r in layers[:-1]:
        d.table("sitelayers").append(list(r))
        d.table("layers").append(list(r))
    d.table("sitelayers").append(list(layers[-1][:1]) + [0.0] + list(layers[-1][2:]))
    d.table("layers").append(list(layers[-1]))
    return d


def add_node(d: decks.Deck, nid: int, x: float, y: float, z: float, fix: Sequence[int] = (0,) * 6) -> None:
    d.table("nodes").append([nid, float(x), float(y), float(z)] + [int(v) for v in fix])


def add_element(d: decks.Deck, group: int, eid: int, nodes: Sequence[int], etype: int = 1, mat: int = 1,
                prop: int = 1, eint: int = 0, thick: float = 0.0, ki: str = "000000", kj: str = "000000") -> None:
    n = (list(nodes) + [0] * 8)[:8]
    d.table("elements").append([group, eid, etype, mat, prop, eint, thick] + n + [ki, kj])


def add_matrix_property(d: decks.Deck, prop: int, kind: str, A: np.ndarray) -> None:
    """Rows of the upper triangle of a symmetric 12x12 matrix (MXR/MXI/MXM layout)."""
    for r in range(12):
        terms = list(A[r, r:]) + [0.0] * r
        d.table("matrices").append([prop, kind, r + 1] + [float(t) for t in terms])


def run_house(workdir: Path, d: decks.Deck, model: str = "m", site: bool = True) -> Tuple[int, str]:
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    if site:
        write_site_deck(workdir, model, gravity=float(d["gravity"]), cmodform=int(d["cmodform"]))
    decks.write(workdir / f"{model}.hou", d)
    rc = run_module("HOUSE", model, workdir)
    return rc, (workdir / f"{model}_HOUSE.out").read_text()


# ---- the VP-H1 mixed model -------------------------------------------------------------------
MIX_MATERIALS = {   # id: (type, val1, val2, weight, pdamp, sdamp)
    1: (1, 3.0e7, 0.2, 24.0, 0.05, 0.05),
    2: (1, 2.0e7, 0.25, 25.0, 0.04, 0.04),
    3: (2, 4.0e7, 1.5e7, 23.0, 0.03, 0.02),
}
MIX_SECTION = dict(A=0.25, As2=0.2, As3=0.18, J=0.01, I2=0.004, I3=0.006)
MIX_SPRING = ((1.0e5, 2.0e5, 3.0e5, 10.0, 20.0, 30.0), 0.03)


def _gm_matrices(seed: int = 3) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((12, 6))
    KR = 1.0e4 * (A @ A.T)
    KI = 0.04 * KR
    MM = np.diag(rng.uniform(0.5, 2.0, 12)) * G_SI          # weight units (MOPT <matrix> = 1)
    return KR, KI, MM


def _nid(ix: int, iy: int, iz: int) -> int:
    """Excavation grid node: x in {0,1,2}, y in {0,1}, z in {-5,-2,0}; numbered bottom-up."""
    return 1 + ix + 3 * iy + 6 * iz


EXC_Z = (-5.0, -2.0, 0.0)


def mixed_house_deck(incomp: int = 0) -> decks.Deck:
    """VP-H1 model: a 4-element excavated block (FV interaction set, two soil layers, ETYPE 0 and 2),
    a structural SOLID, a SHELL slab, two BEAMS (one with an end release), a SPRING, two GENERAL
    elements (2-node global and 3-node local, weight units) and nodal masses in both units."""
    d = house_deck(incomp=incomp, gmunits=1, df=0.5)
    for n in (1, 2, 3):
        d.table("freqs").append([n])
    for iz, z in enumerate(EXC_Z):
        for iy in range(2):
            for ix in range(3):
                add_node(d, _nid(ix, iy, iz), ix, iy, z)
    extra = {19: (0, 0, 0.5), 20: (1, 0, 0.5), 21: (1, 1, 0.5), 22: (0, 1, 0.5), 23: (2, 0, 0.5),
             24: (2, 1, 0.5), 25: (2, 1, 3.5), 27: (2, 1, 5.0), 28: (2.5, 1.2, 3.5)}
    for n, p in extra.items():
        add_node(d, n, *p)
    add_node(d, 26, 3.0, 1.0, 0.5, fix=(1,) * 6)                 # K node only
    for mid, (t, v1, v2, w, pd, sd) in MIX_MATERIALS.items():
        d.table("materials").append([mid, t, v1, v2, w, pd, sd])
    s = MIX_SECTION
    d.table("beamprops").append([1, s["A"], s["As2"], s["As3"], s["J"], s["I2"], s["I3"]])
    k, damp = MIX_SPRING
    d.table("springprops").append([1] + list(k) + [damp])
    KR, KI, MM = _gm_matrices()
    add_matrix_property(d, 1, "R", KR)
    add_matrix_property(d, 1, "I", KI)
    add_matrix_property(d, 1, "M", MM)
    for gid, typ, title in ((1, 1, "basemat block"), (2, 3, "slab"), (3, 2, "column"), (4, 7, "contact spring"),
                            (5, 9, "equipment"), (6, 1, "excavated soil")):
        d.table("groups").append([gid, typ, title])
    add_element(d, 1, 1, [13, 14, 17, 16, 19, 20, 21, 22], etype=1, mat=3, eint=1)
    add_element(d, 2, 1, [20, 23, 24, 21], mat=2, thick=0.2)
    add_element(d, 3, 1, [24, 25, 26], mat=1)
    add_element(d, 3, 2, [25, 27, 26], mat=1, ki="000011")
    add_element(d, 4, 1, [15, 23], prop=1)
    add_element(d, 5, 1, [25, 28], prop=1)
    add_element(d, 5, 2, [28, 27, 26], prop=1)
    e = 0
    for iz in range(2):
        for ix in range(2):
            e += 1
            ns = [_nid(ix, 0, iz), _nid(ix + 1, 0, iz), _nid(ix + 1, 1, iz), _nid(ix, 1, iz),
                  _nid(ix, 0, iz + 1), _nid(ix + 1, 0, iz + 1), _nid(ix + 1, 1, iz + 1), _nid(ix, 1, iz + 1)]
            # bottom layer (z -5..-2) is layer 2 with explicit ETYPE 2; top layer ETYPE 0 (resolved)
            add_element(d, 6, e, ns, etype=2 if iz == 0 else 0, mat=2 if iz == 0 else 1)
    d.table("masses").append([27, 2 * G_SI, 2 * G_SI, 2 * G_SI, 0.3 * G_SI, 0.3 * G_SI, 0.1 * G_SI, 1])
    d.table("masses").append([19, 0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 0])
    d.table("masses").append([26, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0])          # fixed node: ignored (W5)
    for n in range(1, 19):
        d.table("interaction").append([n])
    return d


def mixed_reference(d: decks.Deck, incomp: int = 0):
    """Independent element-library assembly of :func:`mixed_house_deck` from the raw numbers.

    Returns (dofmap, AssembledModel, records)."""
    g = G_SI
    ids = [int(r["id"]) for r in d.rows("nodes")]
    xyz = {int(r["id"]): np.array([r["x"], r["y"], r["z"]], float) for r in d.rows("nodes")}
    fix = np.array([[r[c] for c in ("fx", "fy", "fz", "frx", "fry", "frz")] for r in d.rows("nodes")], int)
    mats = {k: material_from_M(t, v1, v2, w, pd, sd, g) for k, (t, v1, v2, w, pd, sd) in MIX_MATERIALS.items()}
    lay = {r[0]: material_from_layer(r[3], r[4], r[2], r[5], r[6], g) for r in SITE_LAYERS}
    KR, KI, MM = _gm_matrices()
    k, damp = MIX_SPRING
    recs = [ElemRecord(1, 1, 1, (13, 14, 17, 16, 19, 20, 21, 22), False, mats[3],
                       dict(eint=1, incompatible=(incomp == 0))),
            ElemRecord(2, 1, 3, (20, 23, 24, 21), False, mats[2], dict(thick=0.2)),
            ElemRecord(3, 1, 2, (24, 25, 26), False, mats[1], dict(section=dict(MIX_SECTION))),
            ElemRecord(3, 2, 2, (25, 27, 26), False, mats[1], dict(section=dict(MIX_SECTION), ki=(0, 0, 0, 0, 1, 1))),
            ElemRecord(4, 1, 7, (15, 23), False, None, dict(k=k, damp=damp)),
            ElemRecord(5, 1, 9, (25, 28), False, None, dict(KR=KR, KI=KI, MM=MM / g)),
            ElemRecord(5, 2, 9, (28, 27, 26), False, None, dict(KR=KR, KI=KI, MM=MM / g))]
    e = 0
    for iz in range(2):
        for ix in range(2):
            e += 1
            ns = (_nid(ix, 0, iz), _nid(ix + 1, 0, iz), _nid(ix + 1, 1, iz), _nid(ix, 1, iz),
                  _nid(ix, 0, iz + 1), _nid(ix + 1, 0, iz + 1), _nid(ix + 1, 1, iz + 1), _nid(ix, 1, iz + 1))
            recs.append(ElemRecord(6, e, 1, ns, True, lay[2 if iz == 0 else 1], dict(eint=0)))
    masses = {27: np.array([2.0, 2.0, 2.0, 0.3, 0.3, 0.1]), 19: np.array([0.5, 0.5, 0.5, 0, 0, 0]),
              26: np.array([1.0, 1.0, 1.0, 0, 0, 0])}
    dm = build_dofmap(ids, fix, recs, extra_dofs={n: (1, 2, 3) for n in range(1, 19)})
    return dm, assemble(xyz, recs, dm, masses=masses), recs


def _rel(A, B) -> float:
    A = A.toarray() if sp.issparse(A) else np.asarray(A)
    B = B.toarray() if sp.issparse(B) else np.asarray(B)
    nb = max(np.abs(B).max(), 1e-300)
    return float(np.abs(A - B).max() / nb)


# ---- the VP-H1 cantilever stick ------------------------------------------------------------------
STICK = dict(L=10.0, ne=20, E=2.0e8, nu=0.3, rho=7.8, A=0.02, J=3.0e-5, I2=2.5e-5, I3=1.0e-4)


def stick_deck() -> decks.Deck:
    """Vertical cantilever (z up, fixed at z = 0) of 20 BEAMS elements, K node on +X: local axis 2 =
    X (bending with deflection along X uses I3), axis 3 = Y (I2)."""
    s = STICK
    d = house_deck(df=1.0)
    d.table("freqs").append([1])
    ne, L = s["ne"], s["L"]
    for i in range(ne + 1):
        add_node(d, i + 1, 0.0, 0.0, L * i / ne, fix=(1,) * 6 if i == 0 else (0,) * 6)
    add_node(d, 99, 1.0, 0.0, 0.0, fix=(1,) * 6)
    d.table("materials").append([1, 1, s["E"], s["nu"], s["rho"] * G_SI, 0.0, 0.0])
    d.table("beamprops").append([1, s["A"], 0.0, 0.0, s["J"], s["I2"], s["I3"]])
    d.table("groups").append([1, 2, "stick"])
    for i in range(ne):
        add_element(d, 1, i + 1, [i + 1, i + 2, 99])
    return d


EB_BETA_L = (1.875104, 4.694091, 7.854757, 10.995541)


@problem("VP-H1", "HOUSE matrices of a mixed model and a cantilever stick through the HOUSE deck", tier="P0",
         modules=["HOUSE"], source="requirements 4.4; R2 I.1 (Euler-Bernoulli)")
def vp_h1(workdir):
    r = VPResult()
    wd = Path(workdir)
    # ---- (a) mixed model: deck -> COOSK/COOSM/FILE4 equals the element-library assembly ----------
    d = mixed_house_deck()
    rc, out = run_house(wd / "mixed", d)
    r.require("HOUSE run of the mixed model succeeds", rc == 0, out[-600:] if rc else "")
    if rc != 0:
        return r
    dm, am, _ = mixed_reference(d)
    f4 = read_container(wd / "mixed" / "m.N4", "FILE4")
    ck = read_container(wd / "mixed" / "COOSK", "COOSK")
    cm = read_container(wd / "mixed" / "COOSM", "COOSM")
    r.require("equation map (eq_node, eq_dof) equals build_dofmap",
              np.array_equal(f4["eq_node"], dm.eq_node) and np.array_equal(f4["eq_dof"], dm.eq_dof))
    for name, A, B in (("Ks", ck.sparse("Ks"), am.Ks), ("Ke", ck.sparse("Ke"), am.Ke),
                       ("Ms", cm.sparse("Ms"), am.Ms), ("Me", cm.sparse("Me"), am.Me)):
        r.check(f"max |{name}(HOUSE) - {name}(element library)| / max |{name}|", _rel(A, B), 0.0, atol=1e-12)
    for t, rd in am.recovery.items():
        r.check(f"recovery operator {t}: max relative difference", _rel(f4[f"rec_{t}_S"].reshape(-1),
                                                                         rd["S"].reshape(-1)), 0.0, atol=1e-12)
        r.require(f"recovery equations {t} identical", np.array_equal(f4[f"rec_{t}_eq"], rd["eq"]))
    # total structural mass per direction (rigid-body mass incl. the supports), from the raw numbers
    g = G_SI
    beams = MIX_MATERIALS[1][3] / g * MIX_SECTION["A"] * (3.0 + 1.5)          # rho A L of both beams
    struct = MIX_MATERIALS[3][3] / g * 0.5 + MIX_MATERIALS[2][3] / g * 0.2 * 1.0 + beams   # SOLID rho V, SHELL rho t A
    MMg = np.diag(_gm_matrices()[2]) / g                                      # GENERAL masses (weight / g)
    pI, pJ, pK = np.array([2.5, 1.2, 3.5]), np.array([2.0, 1.0, 5.0]), np.array([3.0, 1.0, 0.5])
    e1 = (pJ - pI) / np.linalg.norm(pJ - pI)
    v = (pK - pI) - np.dot(pK - pI, e1) * e1
    e2 = v / np.linalg.norm(v)
    Lam = np.vstack([e1, e2, np.cross(e1, e2)])                               # local axes of GENERAL 5/2
    m91 = read_container(wd / "mixed" / "FILE91", "FILE91").meta
    for c, lab in ((0, "X"), (1, "Y"), (2, "Z")):
        gm = MMg[c] + MMg[6 + c] + sum(MMg[off + a] * Lam[a, c] ** 2 for off in (0, 6) for a in range(3))
        ref = struct + gm + 2.0 + 0.5                       # nodal masses: node 27 (weight / g), node 19 (mass)
        r.check(f"total structural mass {lab} (elements + GENERAL + MT; mass on the K-only node ignored)",
                m91["mass_structure"][c], ref, rtol=1e-12)
    # ---- (b) cantilever stick: fixed-base frequencies vs Euler-Bernoulli ---------------------------
    rc, out = run_house(wd / "stick", stick_deck())
    r.require("HOUSE run of the cantilever stick succeeds", rc == 0, out[-600:] if rc else "")
    if rc != 0:
        return r
    ck = read_container(wd / "stick" / "COOSK", "COOSK")
    cm = read_container(wd / "stick" / "COOSM", "COOSM")
    f4 = read_container(wd / "stick" / "m.N4", "FILE4")
    K = ck.sparse("Ks")
    r.check("stick: max |Im K| (undamped material)", float(np.abs(K.imag).max()) if K.nnz else 0.0, 0.0, atol=0.0)
    f, phi = natural_frequencies(K.real, cm.sparse("Ms"), return_modes=True)
    dof = f4["eq_dof"]
    s = STICK
    mbar = s["rho"] * s["A"]
    found = {"X": [], "Y": []}
    for k in range(f.size):
        e = np.array([np.sum(phi[dof == c, k] ** 2) for c in (1, 2, 3, 6)])
        if e.sum() <= 0:
            continue
        frac = e / e.sum()
        if frac[0] > 0.9:
            found["X"].append(f[k])
        elif frac[1] > 0.9:
            found["Y"].append(f[k])
    for direction, I, label in (("X", s["I3"], "I3 (deflection along local axis 2 = X)"),
                                ("Y", s["I2"], "I2 (deflection along local axis 3 = Y)")):
        fx = found[direction]
        for n, ref in enumerate(EB_BETA_L):
            if n >= len(fx):
                r.require(f"bending mode {n + 1} along {direction} found", False)
                continue
            bl = s["L"] * ((2 * math.pi * fx[n]) ** 2 * mbar / (s["E"] * I)) ** 0.25
            r.check(f"cantilever beta_{n + 1} L, bending along {direction}, {label}", bl, ref, rtol=0.01)
    r.notes.append("The stick is written as a HOUSE deck (20 BEAMS, consistent mass, As = 0) and solved from "
                   "COOSK/COOSM with the undamped real stiffness; modes are classified by their UX/UY content.")
    return r


# ---- VP-H2: excavated soil block -------------------------------------------------------------------
def soil_block_deck(distort: Tuple[float, float] = (0.2, -0.1), incomp: int = 0) -> Tuple[decks.Deck, dict]:
    """2 x 2 x 2 excavated block of SOLIDs, x, y in [0, 4], z in [-5, 0] (layer 1 over layer 2), with
    the centre column of nodes moved in plan by ``distort`` (prismatic distortion: element volumes
    are plan area x height), every node an interaction node (FV).  Nodes are numbered bottom-up
    (R4) and the excavation groups from the surface down, one group per embedment layer (R6)."""
    d = house_deck(incomp=incomp)
    xs, ys = (0.0, 2.0, 4.0), (0.0, 2.0, 4.0)
    nid = {}
    n = 0
    for iz, z in enumerate(EXC_Z[0:1] + (-2.0, 0.0)):
        for iy in range(3):
            for ix in range(3):
                n += 1
                x, y = xs[ix], ys[iy]
                if ix == 1 and iy == 1:
                    x, y = x + distort[0], y + distort[1]
                nid[(ix, iy, iz)] = n
                add_node(d, n, x, y, z)
                d.table("interaction").append([n])
    d.table("groups").append([1, 1, "excavated soil, layer 1"])
    d.table("groups").append([2, 1, "excavated soil, layer 2"])
    e = {1: 0, 2: 0}
    for iz in range(2):
        g = 2 if iz == 0 else 1
        for iy in range(2):
            for ix in range(2):
                e[g] += 1
                ns = [nid[(ix, iy, iz)], nid[(ix + 1, iy, iz)], nid[(ix + 1, iy + 1, iz)], nid[(ix, iy + 1, iz)],
                      nid[(ix, iy, iz + 1)], nid[(ix + 1, iy, iz + 1)], nid[(ix + 1, iy + 1, iz + 1)],
                      nid[(ix, iy + 1, iz + 1)]]
                add_element(d, g, e[g], ns, etype=2, mat=2 if iz == 0 else 1)
    return d, nid


@problem("VP-H2", "HOUSE excavated soil: mass rho V and stiffness identity with the layer properties", tier="P0",
         modules=["HOUSE"], source="requirements 4.4 item 2, D-ELM-12, D-CNV-05")
def vp_h2(workdir):
    r = VPResult()
    wd = Path(workdir)
    d, nid = soil_block_deck()
    rc, out = run_house(wd, d)
    r.require("HOUSE run of the excavated block succeeds", rc == 0, out[-600:] if rc else "")
    if rc != 0:
        return r
    f4 = read_container(wd / "m.N4", "FILE4")
    ck = read_container(wd / "COOSK", "COOSK")
    cm = read_container(wd / "COOSM", "COOSM")
    Ke, Me = ck.sparse("Ke").toarray(), cm.sparse("Me").toarray()
    # ---- total mass = sum rho V per layer (plan area 16 is independent of the interior distortion) ------
    lay = {row[0]: row for row in SITE_LAYERS}
    rho1, rho2 = lay[1][2] / G_SI, lay[2][2] / G_SI
    m_ref = 16.0 * (2.0 * rho1 + 3.0 * rho2)
    eqd = f4["eq_dof"]
    for c, lab in ((1, "X"), (2, "Y"), (3, "Z")):
        rv = (eqd == c).astype(float)
        with blas_quiet():                    # spurious Accelerate matmul flags (sassi.elements.base)
            mtot = float(rv @ (Me @ rv))
        r.check(f"excavated mass r_{lab}^T Me r_{lab} = sum rho V", mtot, m_ref, rtol=1e-12)
    # ---- element-by-element identity with the layer material (hand scatter, no assembler) -------------
    eq = {(int(n), int(k)): i for i, (n, k) in enumerate(zip(f4["eq_node"], f4["eq_dof"]))}
    xyz = {int(n): p for n, p in zip(f4["node_id"], f4["node_xyz"])}
    Kref = np.zeros_like(Ke)
    Mref = np.zeros_like(Me)
    mats = {l: material_from_layer(row[3], row[4], row[2], row[5], row[6], G_SI) for l, row in lay.items()}
    vol = 0.0
    for row in d.rows("elements"):
        ns = [int(row[f"n{k}"]) for k in range(1, 9)]
        X = np.array([xyz[n] for n in ns])
        Ke_e, Me_e = solid.matrices(X, mats[int(row["mat"])], incompatible=False, eint=0)
        idx = np.array([eq[(n, c)] for n in ns for c in (1, 2, 3)])
        Kref[np.ix_(idx, idx)] += Ke_e
        Mref[np.ix_(idx, idx)] += Me_e
        vol += Me_e[0::3, 0::3].sum() / mats[int(row["mat"])].rho
    r.check("Ke(HOUSE) = sum of SOLID K*(layer properties): max relative difference", _rel(Ke, Kref), 0.0,
            atol=1e-12)
    r.check("Me(HOUSE) = sum of SOLID M(layer properties): max relative difference", _rel(Me, Mref), 0.0,
            atol=1e-12)
    r.check("excavated volume from the element masses", vol, 16.0 * 5.0, rtol=1e-12)
    Ks = ck.sparse("Ks")
    r.check("no structure: max |Ks|", float(abs(Ks).max()) if Ks.nnz else 0.0, 0.0, atol=0.0)
    # ---- complex moduli of HOUSE's Ke from affine fields, independent of the element and material code
    # (final audit: the former check compared |G*| of the VP's own material object with rho Vs^2, which
    # |c(beta)| = 1 makes blind to the damping).  For an affine field the SOLID reproduces the constant
    # strain exactly, so u^T Ke u = sum_l G*_l gamma^2 V_l (pure shear u_x = gamma z) and
    # sum_l M*_l eps^2 V_l (uniaxial strain u_z = eps z), with G* = rho Vs^2 c(beta_s), M* = rho Vp^2 c(beta_p)
    # and c(b) = 1 - 2 b^2 + 2 i b sqrt(1 - b^2) written out here.
    def c_sassi(b):
        return 1.0 - 2.0 * b * b + 2j * b * math.sqrt(1.0 - b * b)
    zeq = np.array([xyz[int(n)][2] for n in f4["eq_node"]])
    vols = {1: 16.0 * 2.0, 2: 16.0 * 3.0}
    for dof, lab, vel, damp in ((1, "G* (pure shear u_x = z)", 4, 6), (3, "M* (uniaxial strain u_z = z)", 3, 5)):
        u = np.where(eqd == dof, zeq, 0.0)
        with blas_quiet():
            energy = complex(u @ (Ke @ u))
        ref = sum(lay[l][2] / G_SI * lay[l][vel] ** 2 * c_sassi(lay[l][damp]) * vols[l] for l in (1, 2))
        r.check(f"HOUSE Ke: sum of {lab} V over the layers = sum rho V^2 c(beta) V (complex, relative)",
                abs(energy - ref) / abs(ref), 0.0, atol=1e-10)
    # ---- interfaces of the interaction nodes -----------------------------------------------------------
    z = np.array([xyz[int(n)][2] for n in f4["int_node"]])
    expect = np.select([z == 0.0, z == -2.0, z == -5.0], [1, 2, 3], -1)
    r.require("interaction-node interfaces (z = 0, -2, -5 -> 1, 2, 3)", np.array_equal(f4["int_iface"], expect))
    r.require("incompatible modes suppressed on excavated soil although MOPT <incomp> = 0 (D-ELM-02)",
              _rel(Ke, Kref) <= 1e-12)
    r.notes.append("Reference: sassi.elements.solid.matrices with material_from_layer(Vp, Vs, weight, dp, ds, g) "
                   "scattered by hand on the FILE4 equation numbers.")
    return r
