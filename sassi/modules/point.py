"""POINT module: point-load solutions of the layered soil -> FILE3 (requirements §4.3, R1 §3).

POINT3 (3D, HOUSE <dim> = 2): for every analysis frequency and every load interface
n = 1 .. <layer>+1, a unit point load is applied on the axis of an axisymmetric central zone of
radius R0 = <rad> (one radial element per SITE sublayer, D-PNT-03) surrounded by the consistent
cylindrical transmitting boundary built from the SITE Mode 1 modes (FILE2).  FILE3 stores the
exterior mode amplitudes (mu = 0 vertical load, mu = 1 horizontal load), the core-axis
displacements and the mode rows at the observation interfaces, so that ANALYS evaluates the
flexibility exactly at any distance (D-PNT-02, :func:`sassi.core.flexibility.flexibility_matrix`).

POINT2 (2D, <dim> = 1, tier P1): plane-strain strip |x| <= R0 with Waas-Lysmer boundaries
(R1 §3.5), see :mod:`sassi.core.strip2d`.

    python -m sassi.modules.point < POINT.inp      # batch use
"""
from __future__ import annotations

import numpy as np

from ..core import axisym, strip2d, tlm
from ..io import decks
from ..io.files import read_container, validate, write_container
from .base import ModuleContext, ModuleError, batch_main

NAME = "POINT"


def _frequency_rows(d: decks.Deck, f2, lst) -> np.ndarray:
    """Rows of FILE2 for the deck frequencies (all of FILE2 when the deck lists none)."""
    fn = [int(r["number"]) for r in d.rows("freqs")]
    if not fn:
        lst.write(" No frequency numbers in the deck: all FILE2 frequencies are used")
        return np.arange(len(f2["fnum"]))
    if len(set(fn)) != len(fn):
        raise ModuleError("frequency set contains duplicate frequency numbers")
    rows = []
    for n in sorted(fn):
        idx = np.flatnonzero(f2["fnum"] == n)
        if idx.size == 0:
            raise ModuleError(f"frequency number {n} not in FILE2 -- re-run SITE Mode 1 with this frequency")
        rows.append(int(idx[0]))
    return np.asarray(rows, int)


def run(ctx: ModuleContext) -> int:
    d = decks.read(ctx.deck_path, NAME)
    lst = ctx.listing
    L = int(d["layer"])
    R0 = float(d["rad"])
    dim = int(d["dim"])
    if L < 0:
        raise ModuleError("Error 56: last layer number in near field zone < 0")
    if R0 <= 0:
        raise ModuleError("Error 57: radius of central zone must be positive")
    if dim not in (1, 2):
        raise ModuleError(f"<dim> = {dim}: POINT needs 1 (2D, POINT2) or 2 (3D, POINT3)")
    f2 = read_container(ctx.require("FILE2", "SITE (Mode 1)"), "FILE2")
    probs = validate(f2)
    if probs:
        raise ModuleError("FILE2: " + "; ".join(probs))
    df = float(f2.meta["df"])
    if float(d["df"]) > 0 and abs(float(d["df"]) - df) > 1e-9 * df:
        raise ModuleError(f"frequency step {d['df']} differs from FILE2 ({df})")
    rows = _frequency_rows(d, f2, lst)
    nI = len(f2["depth_user"])
    iface_user = np.asarray(f2["iface_user"], int)
    load_iface = np.arange(1, L + 2)
    if L + 1 > nI:
        raise ModuleError(f"POINT <layer> = {L} needs interface {L + 1}, but the site has {nI} interfaces")
    load_idx = iface_user[load_iface - 1]
    if np.any(load_idx < 0):
        raise ModuleError(f"interface {L + 1} is the rigid base (SITE nl = 0): POINT <layer> must be < {nI - 1}")
    mass = str(f2.meta.get("mass", tlm.MASS_MIXED))
    program = "POINT3 (axisymmetric central zone)" if dim == 2 else "POINT2 (plane-strain strip)"
    lst.section(f"{program} input")
    lst.write(f" Title                          : {d['title']}")
    lst.write(f" Radius of central zone R0      : {R0:.6g}")
    lst.write(f" Last layer in near field zone  : {L}  -> loads at interfaces 1..{L + 1}")
    lst.write(f" Frequencies                    : {rows.size} of FILE2 (df = {df:.9g} Hz)")
    lst.write(f" Column                         : {nI - 1} user layers + {int(f2.meta['nl'])} generated sublayers, "
              f"{f2.meta['base']} base, {mass} mass")
    depth = np.asarray(f2["depth_user"])
    lst.table(["load iface", "depth", "free index"], [[str(i), depth[i - 1], str(load_idx[k])]
                                                       for k, i in enumerate(load_iface)], fmt="{:>14.6g}")
    if int(d["opmode"]) == 1:
        lst.write("")
        lst.write(" Data check only (<opmode> = 1): FILE3 not written")
        return 0
    if dim == 1:
        arrays, meta = _point2(f2, rows, R0, load_iface, load_idx, mass, ctx)
    else:
        arrays, meta = _point3(f2, rows, R0, load_iface, load_idx, mass, ctx)
    write_container(ctx.path("FILE3"), "FILE3", arrays, meta, module=NAME)
    _listing_results(lst, arrays, R0, dim)
    ctx.progress(1.0, "POINT done")
    return 0


def _point3(f2, rows, R0, load_iface, load_idx, mass, ctx):
    nF = rows.size
    nL = load_iface.size
    out = {"alpha0": [], "alpha1": [], "axis0": [], "axis1": [], "phix_obs": [], "phiz_obs": [], "phiy_obs": [],
           "kR": [], "kL": [], "x_rim0": [], "x_rim1": [], "x_loss": []}
    for i, q in enumerate(rows):
        if ctx.cancelled():
            raise ModuleError("run cancelled")
        f = float(f2["freq"][q])
        col = tlm.column_from_file2(f2, int(q))
        md = tlm.modes_from_file2(f2, int(q))
        try:
            sol = axisym.solve_point3(col, md, 2 * np.pi * f, R0, load_idx, mass)
        except (FloatingPointError, np.linalg.LinAlgError) as exc:
            raise ModuleError(f"POINT3 failed at {f:.6g} Hz: {exc}") from None
        out["alpha0"].append(sol.alpha[0])
        out["alpha1"].append(sol.alpha[1])
        out["axis0"].append(sol.axis[0])
        out["axis1"].append(sol.axis[1])
        out["x_rim0"].append(sol.rim[0])
        out["x_rim1"].append(sol.rim[1])
        out["x_loss"].append(tlm.column_loss(col))     # propagating sector of the table grid (D-SIT-09)
        out["phix_obs"].append(md.phix[load_idx])
        out["phiz_obs"].append(md.phiz[load_idx])
        out["phiy_obs"].append(md.phiy[load_idx])
        out["kR"].append(md.kR)
        out["kL"].append(md.kL)
        ctx.progress((i + 1) / nF, f"POINT3: frequency {i + 1}/{nF} ({f:.4g} Hz)")
    arrays = {k: np.stack(v) for k, v in out.items()}
    arrays.update({"fnum": np.asarray(f2["fnum"])[rows], "freq": np.asarray(f2["freq"])[rows],
                   "depth_user": np.asarray(f2["depth_user"]), "load_iface": load_iface})
    meta = {"R0": R0, "dim": 2, "df": float(f2.meta["df"]), "nl": int(f2.meta["nl"]), "base": f2.meta["base"],
            "mass": mass, "program": "POINT3",
            "alpha_basis": "column j of Psi_mu(rho) uses H_mu(k_j rho) exp(i k_j R0) (scaled to R0)",
            "axis_index": "[frequency, load interface, observation interface, (u_rho, u_theta, u_z)]",
            "load": "unit point load: mu = 0 vertical (upward, +z), mu = 1 horizontal (+x)"}
    return arrays, meta


def _point2(f2, rows, R0, load_iface, load_idx, mass, ctx):
    """POINT2 (P1): plane-strain line-load solutions stored in the FILE3 layout (dim = 1).

    ``alpha1[q, n] = [P-SV amplitudes of the x line load (2Nf); SH amplitudes of the y line load
    (Nf)]``, ``alpha0[q, n] = [P-SV amplitudes of the z line load; 0]``;
    ``axis1[..., :] = (u_x | x load, u_y | y load, u_z | x load)`` and
    ``axis0[..., :] = (u_x | z load, 0, u_z | z load)`` at the centre node (index [q, load, obs]).
    """
    nF = rows.size
    out = {k: [] for k in ("alpha0", "alpha1", "axis0", "axis1", "phix_obs", "phiz_obs", "phiy_obs", "kR", "kL",
                           "x_loss")}
    for i, q in enumerate(rows):
        if ctx.cancelled():
            raise ModuleError("run cancelled")
        f = float(f2["freq"][q])
        col = tlm.column_from_file2(f2, int(q))
        md = tlm.modes_from_file2(f2, int(q))
        try:
            sol = strip2d.solve_point2(col, md, 2 * np.pi * f, R0, load_idx, mass)
        except (FloatingPointError, np.linalg.LinAlgError) as exc:
            raise ModuleError(f"POINT2 failed at {f:.6g} Hz: {exc}") from None
        nL = load_idx.size
        out["alpha1"].append(np.hstack([sol.alpha_x, sol.alpha_y]))
        out["alpha0"].append(np.hstack([sol.alpha_z, np.zeros_like(sol.alpha_y)]))
        c = sol.centre
        out["axis1"].append(np.stack([c["xx"], c["yy"], c["zx"]], -1))
        out["axis0"].append(np.stack([c["xz"], np.zeros((nL, nL), complex), c["zz"]], -1))
        out["phix_obs"].append(md.phix[load_idx])
        out["phiz_obs"].append(md.phiz[load_idx])
        out["phiy_obs"].append(md.phiy[load_idx])
        out["kR"].append(md.kR)
        out["kL"].append(md.kL)
        out["x_loss"].append(tlm.column_loss(col))
        ctx.progress((i + 1) / nF, f"POINT2: frequency {i + 1}/{nF} ({f:.4g} Hz)")
    arrays = {k: np.stack(v) for k, v in out.items()}
    arrays.update({"fnum": np.asarray(f2["fnum"])[rows], "freq": np.asarray(f2["freq"])[rows],
                   "depth_user": np.asarray(f2["depth_user"]), "load_iface": load_iface})
    meta = {"R0": R0, "dim": 1, "df": float(f2.meta["df"]), "nl": int(f2.meta["nl"]), "base": f2.meta["base"],
            "mass": mass, "program": "POINT2",
            "alpha_basis": "right-boundary modes exp(-i k_j (|x| - R0)); P-SV then SH",
            "axis_index": "[frequency, load interface, observation interface, component]",
            "load": "unit line loads per unit length at x = 0"}
    return arrays, meta


def _listing_results(lst, arrays, R0, dim) -> None:
    lst.section("Self-flexibility on the load axis (core-axis displacement under a unit load)")
    nL = arrays["load_iface"].size
    rows = []
    for q in range(arrays["fnum"].size):
        for k in range(nL):
            uh = arrays["axis1"][q, k, k, 0]
            uv = arrays["axis0"][q, k, k, 2]
            rows.append([float(arrays["freq"][q]), str(int(arrays["load_iface"][k])), uh.real, uh.imag, uv.real, uv.imag])
    lst.table(["f (Hz)", "iface", "Re u_h", "Im u_h", "Re u_v", "Im u_v"], rows[:200], fmt="{:>14.6g}")
    if len(rows) > 200:
        lst.write(f" ... {len(rows) - 200} more rows in FILE3")
    if dim == 2:
        lst.write(f" u_h: horizontal displacement per unit horizontal load; u_v: vertical per unit vertical load"
                  f" (central zone R0 = {R0:.6g})")


if __name__ == "__main__":
    raise SystemExit(batch_main(NAME))
