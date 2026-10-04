"""Final audit: complex-modulus and damping conventions are the same in SOIL, SITE and HOUSE
(requirements 0.3, 4.0.2, D-CNV-04: beta_s on G, beta_p on M = lambda + 2G, c(beta) = 1 - 2 beta^2 +
2 i beta sqrt(1 - beta^2), CMODFORM 1: 1 + 2 i beta; mass density = weight / g).

The verification problems use layers with beta_p = beta_s, so an exchange of the two damping ratios, or a
module that ignores CMODFORM, would go unnoticed (VP audit, docs/verification/AUDIT_REPORT.md).  Here the
soil has beta_p = 0.03 and beta_s = 0.08 and every module is compared with one closed form written out in
this file (no sassi helper on the reference side):

* SITE, vertical SV and P (thin layers, 20 sublayers in 2 m): within-layer motion
  u(z)/u(0) = cos(k* z), k* = w / (V sqrt(c(beta))) -- the free-surface solution of a uniform layer (observed
  7.5e-10: the mixed mass cancels the O((k h)^2) thin-layer error of vertical waves),
  independent of the half-space below;
* SOIL (SHAKE, continuous), horizontal (SH, beta_s) and vertical (P, beta_p, SOILX <indir> = 1) input:
  surface / base-within amplification 1/cos(k* H);
* HOUSE: G* and M* recovered from the assembled Ke (excavated soil, L table) and Ks (structure SOLIDs,
  M table type 1 E-nu and type 3 Vp-Vs) with affine fields: u^T K u = G* gamma^2 V (u_x = gamma z) and
  M* eps^2 V (u_z = eps z); the masses r^T M r = rho V.
"""
from __future__ import annotations

import cmath
import math
import shutil
from pathlib import Path

import numpy as np
import pytest

from sassi.io import textfiles
from sassi.io.container import read_container
from sassi.prep import Interpreter

EXAMPLES = Path(__file__).resolve().parents[2] / "examples"
DYNP = Path(__file__).resolve().parents[2] / "sassi" / "data" / "dynp_library.pre"
G = 9.80665
W, VP, VS, BP, BS = 18.0, 400.0, 200.0, 0.03, 0.08            # soil layer (weight, velocities, damping)
H_SUB, N_SUB = 0.1, 20                                        # 20 sublayers of 0.1 m: a uniform 2 m layer
E1, NU1, W1, BP1, BS1 = 3.0e7, 0.2, 24.0, 0.02, 0.06          # structure M type 1 (E, nu)
VP3, VS3, W3, BP3, BS3 = 900.0, 450.0, 21.0, 0.01, 0.04       # structure M type 3 (Vp, Vs)
FNUM = [41, 123, 205, 287, 369, 410]                           # 2 .. 20 Hz (df = 1/(8192 x 0.005))


def c_of(beta: float, form: int) -> complex:
    return 1 + 2j * beta if form == 1 else 1 - 2 * beta * beta + 2j * beta * math.sqrt(1 - beta * beta)


def _model(form: int, wave: str) -> str:
    topl = ",".join(["1"] * N_SUB)
    wave_cmd = "WAVE,2,1,1,1,0" if wave == "SV" else "WAVE,3,1,1,1,0"
    cm = 0 if wave == "SV" else 2
    soil = ([f"INP,{DYNP}"] + [f"SPRO,{k},1,Sand" for k in range(1, N_SUB + 1)]   # curves unused: 0 iterations
            + [f"SPRO,{N_SUB + 1},2"])
    lines = [
        "MDL,mat,mat", f"CMODFORM,{form}",
        f"L,1,{H_SUB},{W},{VP},{VS},{BP},{BS}",
        "L,2,1.0,20.0,1600,800,0.02,0.02",
        f"TOPL,{topl}",
        "FREQ,1," + ",".join(str(n) for n in FNUM),
        f"SITE,0,1,0,20,2,1,{1 if wave == 'SH' else 0},1,4096,1,{cm},0.005,8192,1", wave_cmd,
        f"HOUSE,{G},0,0,2,0,0,0,0,0", "POINT,0,1,0.9",
        # excavated hex (L layer 1, below grade) and two isolated structure hexes above grade
        "N,1,0,0,-0.1", "N,2,2,0,-0.1", "N,3,2,2,-0.1", "N,4,0,2,-0.1",
        "N,5,0,0,0", "N,6,2,0,0", "N,7,2,2,0", "N,8,0,2,0",
        "N,11,10,0,1", "N,12,12,0,1", "N,13,12,2,1", "N,14,10,2,1",
        "N,15,10,0,2", "N,16,12,0,2", "N,17,12,2,2", "N,18,10,2,2",
        "N,21,20,0,1", "N,22,22,0,1", "N,23,22,2,1", "N,24,20,2,1",
        "N,25,20,0,3", "N,26,22,0,3", "N,27,22,2,3", "N,28,20,2,3",
        f"M,1,{E1},{NU1},{W1},{BP1},{BS1},1",
        f"M,2,{VP3},{VS3},{W3},{BP3},{BS3},3",
        "GROUP,1,SOLID", "MACT,1", "E,1,1,2,3,4,5,6,7,8", "ETYPE,1,1,1,2",
        "GROUP,2,SOLID", "MACT,1", "E,1,11,12,13,14,15,16,17,18", "ETYPE,1,1,1,1",
        "GROUP,3,SOLID", "MACT,2", "E,1,21,22,23,24,25,26,27,28", "ETYPE,1,1,1,1",
        "MOPT,1",
        "INT,1,8,1,1",
        # SOIL on the same profile: input within at the top of the half-space (sublayer 21), linear
        *soil,
        f"SOIL,4000,{G},1,0,0,0,0.65,1,0", f"SOILX,{1 if wave == 'P' else 0},1,0,{N_SUB + 1}",
        "THFILE,../motion.acc",
        f"SSAF,1,1,0,0,{N_SUB + 1},0.48828125,surface / base within",
        "DAMP,0.05",
        "AOPT,0,1,0,1,0,1,0,0,0,0,0,0,0,0",
        "CHECK", "AFWRITE", "RUNSOIL", "RUNSITE", "RUNHOUSE",
    ]
    return "\n".join(lines) + "\n"


def _run(root: Path, form: int, wave: str) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(EXAMPLES / "data" / "rg160h_030g.acc", root / "motion.acc")
    pre = root / "m.pre"
    pre.write_text(_model(form, wave), encoding="utf-8")
    ui = Interpreter(cwd=root)
    summary = ui.run_file(str(pre), resolve=False)
    err = (root / "mat" / "mat.err").read_text() if (root / "mat" / "mat.err").exists() else ""
    assert summary.errors == 0, err
    return root / "mat"


@pytest.fixture(scope="module", params=[0, 1], ids=["sassi-form", "cmodform-1"])
def runs(request, tmp_path_factory):
    form = request.param
    root = tmp_path_factory.mktemp(f"audit_mat{form}")
    return form, _run(root / "SV", form, "SV"), _run(root / "P", form, "P")


def test_site_free_field_uses_beta_s_for_sv_and_beta_p_for_p(runs):
    form, sv, p = runs
    for md, comp, v, beta in ((sv, 0, VS, BS), (p, 2, VP, BP)):
        f1 = read_container(md / "FILE1", "FILE1")
        U = np.asarray(f1["U"])[0]                                    # (nF, nI, 3)
        freq = np.asarray(f1["freq"])
        z = H_SUB * np.arange(N_SUB + 1)
        worst, swapped = 0.0, 0.0
        for q, f in enumerate(freq):
            ks = 2 * math.pi * f / (v * cmath.sqrt(c_of(beta, form)))
            ref = np.array([cmath.cos(ks * zz) for zz in z])
            got = U[q, : N_SUB + 1, comp] / U[q, 0, comp]
            worst = max(worst, float(np.abs(got - ref).max()))
            kw = 2 * math.pi * f / (v * cmath.sqrt(c_of(BP + BS - beta, form)))      # the other damping ratio
            swapped = max(swapped, float(np.abs(got - np.array([cmath.cos(kw * zz) for zz in z])).max()))
        assert swapped > 10 * worst, (md.parent.name, worst, swapped)    # the damping assignment is resolved
        # mixed mass (D-CNV-05) cancels the O((k h)^2) error of a vertically propagating wave: observed 7.5e-10
        assert worst < 1e-7, (md.parent.name, worst)


def test_soil_amplification_uses_beta_s_horizontally_and_beta_p_vertically(runs):
    form, sv, p = runs
    for md, v, beta in ((sv, VS, BS), (p, VP, BP)):
        f, H, cplx = textfiles.read_tf(md / f"SAF001W_{N_SUB + 1:03d}W.TFU")
        assert cplx
        f, H = np.asarray(f), np.asarray(H)
        keep = (f > 0) & (f < 30)
        ref = np.array([1.0 / cmath.cos(2 * math.pi * fi / (v * cmath.sqrt(c_of(beta, form))) * H_SUB * N_SUB)
                        for fi in f[keep]])
        # the file prints the frequency with 6 decimals (0.488281 Hz) and H with 9 digits: ~1e-6 near resonance
        assert float(np.max(np.abs(H[keep] - ref) / np.abs(ref))) < 2e-6, md.parent.name


def _energy(K, eq_node, eq_dof, xyz, nodes, dof, along):
    u = np.zeros(K.shape[0])
    for i, (n, d) in enumerate(zip(eq_node, eq_dof)):
        if int(n) in nodes and int(d) == dof:
            u[i] = xyz[int(n)][along]
    return complex(u @ (K @ u))


def test_house_moduli_masses_and_damping_assignment(runs):
    form, sv, _ = runs
    f4 = read_container(sv / "mat.N4")
    ck, cm = read_container(sv / "COOSK", "COOSK"), read_container(sv / "COOSM", "COOSM")
    Ks, Ke = ck.sparse("Ks").tocsr(), ck.sparse("Ke").tocsr()
    Ms, Me = cm.sparse("Ms").tocsr(), cm.sparse("Me").tocsr()
    eq_node, eq_dof = np.asarray(f4["eq_node"]), np.asarray(f4["eq_dof"])
    xyz = {int(n): p for n, p in zip(f4["node_id"], f4["node_xyz"])}
    G1, M1 = E1 / (2 * (1 + NU1)), E1 * (1 - NU1) / ((1 + NU1) * (1 - 2 * NU1))
    cases = [
        ("excavated (L)", Ke, Me, set(range(1, 9)), 2 * 2 * H_SUB, W / G, W / G * VS ** 2 * c_of(BS, form),
         W / G * VP ** 2 * c_of(BP, form)),
        ("structure M type 1", Ks, Ms, set(range(11, 19)), 4.0, W1 / G, G1 * c_of(BS1, form), M1 * c_of(BP1, form)),
        ("structure M type 3", Ks, Ms, set(range(21, 29)), 8.0, W3 / G, W3 / G * VS3 ** 2 * c_of(BS3, form),
         W3 / G * VP3 ** 2 * c_of(BP3, form)),
    ]
    for name, K, M, nodes, V, rho, Gs, Ms_ in cases:
        g = _energy(K, eq_node, eq_dof, xyz, nodes, 1, 2) / V            # u_x = z: gamma_xz = 1
        m = _energy(K, eq_node, eq_dof, xyz, nodes, 3, 2) / V            # u_z = z: eps_zz = 1
        assert abs(g - Gs) / abs(Gs) < 1e-10, (name, g, Gs)
        assert abs(m - Ms_) / abs(Ms_) < 1e-10, (name, m, Ms_)
        r = np.array([1.0 if (int(n) in nodes and int(d) == 1) else 0.0 for n, d in zip(eq_node, eq_dof)])
        assert abs(float(r @ (M @ r)) - rho * V) < 1e-10 * rho * V, name
