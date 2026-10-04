"""Normative contents of the binary inter-module files (array names, dtypes and shapes).

Every file is a :mod:`sassi.io.container` archive.  Symbols used in shapes:

* ``nF``   number of SSI frequencies in the file (sorted ascending frequency numbers)
* ``nI``   number of *user* interfaces = number of TOPL layers + 1 (the last one is the
           top of the half-space or the rigid base); interface 1 is the ground surface
* ``Nf``   number of free interfaces of the discretised column (user + generated
           sublayer interfaces, including the dashpot base; a rigid base is removed)
* ``nW``   number of active wave fields; ``nL`` number of POINT load interfaces
           (POINT <layer> + 1); ``nEq`` number of active equations (DOFs)
* depths are positive down from the top of TOPL layer 1 (= ground surface); a model
  elevation z maps to depth ``gelev - z`` (HOUSE <gelev>).

Producers must write at least the listed arrays and meta keys; consumers may rely on
them.  Additional arrays are allowed (prefix them with ``x_``).  ``validate`` checks
presence and leading dimensions.
"""
from __future__ import annotations

from typing import Dict, List, Tuple

from .container import Container, read_container, write_container  # noqa: F401  (re-export)

# kind -> (required arrays {name: (dtype kind, shape description)}, required meta keys)
SCHEMAS: Dict[str, Tuple[Dict[str, Tuple[str, str]], List[str]]] = {
    # ---------------- SITE Mode 1 -> SITE Mode 2, POINT ------------------------------------
    "FILE2": ({
        "fnum": ("i", "(nF,) frequency numbers"),
        "freq": ("f", "(nF,) Hz"),
        "depth_user": ("f", "(nI,) depth of user interfaces"),
        "layer_thick": ("f", "(nI-1,) user layer thicknesses"),
        "layer_rho": ("f", "(nI,) mass density of user layers + half-space (last)"),
        "layer_G": ("c", "(nI,) complex shear modulus G* of user layers + half-space"),
        "layer_M": ("c", "(nI,) complex constrained modulus M* = lam*+2G*"),
        "h_gen": ("f", "(nF, nl) generated half-space sublayer thicknesses per frequency (nl may be 0)"),
        "iface_user": ("i", "(nI,) index of each user interface among the Nf free interfaces (-1 if removed rigid base)"),
        "kR": ("c", "(nF, 2Nf) Rayleigh wavenumbers, Im k < 0 (or Re k > 0 if real)"),
        "phix": ("c", "(nF, Nf, 2Nf) Rayleigh mode shapes, x component (Kausel normalisation)"),
        "phiz": ("c", "(nF, Nf, 2Nf) Rayleigh mode shapes, i-scaled z component (u_z = -i phi_z)"),
        "kL": ("c", "(nF, Nf) Love wavenumbers"),
        "phiy": ("c", "(nF, Nf, Nf) Love mode shapes"),
    }, ["df", "nl", "base", "mass", "cmodform"]),
    # ---------------- SITE Mode 2 -> ANALYS ----------------------------------------------
    "FILE1": ({
        "fnum": ("i", "(nF,)"),
        "freq": ("f", "(nF,) Hz"),
        "depth_user": ("f", "(nI,)"),
        "U": ("c", "(nW, nF, nI, 3) unit-normalised free field per wave, components x', y', z' (physical, z' up)"),
        "k": ("c", "(nW, nF) horizontal wavenumber of each wave along +x' (0 for vertical incidence)"),
        "ratio": ("f", "(nW, nF) participation ratio of each wave"),
        "wave_type": ("i", "(nW,) 1 R, 2 SV, 3 P, 4 SH, 5 L"),
    }, ["df", "cl", "cm", "wopt"]),
    # ---------------- POINT -> ANALYS ----------------------------------------------------
    "FILE3": ({
        "fnum": ("i", "(nF,)"),
        "freq": ("f", "(nF,) Hz"),
        "depth_user": ("f", "(nI,)"),
        "load_iface": ("i", "(nL,) user interface numbers (1-based) carrying loads: 1..layer+1"),
        "kR": ("c", "(nF, 2Nf)"),
        "kL": ("c", "(nF, Nf)"),
        "phix_obs": ("c", "(nF, nL, 2Nf) Rayleigh x mode rows at the observation (= load) interfaces"),
        "phiz_obs": ("c", "(nF, nL, 2Nf) Rayleigh i-scaled z mode rows"),
        "phiy_obs": ("c", "(nF, nL, Nf) Love mode rows"),
        "alpha0": ("c", "(nF, nL, 3Nf) exterior mode amplitudes, vertical unit load (mu = 0) at each load interface"),
        "alpha1": ("c", "(nF, nL, 3Nf) exterior mode amplitudes, horizontal unit load (mu = 1)"),
        "axis0": ("c", "(nF, nL, nL, 3) core-axis displacements (u_rho~, u_theta~, u_z~) at observation interfaces, mu = 0"),
        "axis1": ("c", "(nF, nL, nL, 3) core-axis displacements, mu = 1"),
    }, ["R0", "dim", "df"]),
    # ---------------- HOUSE -> ANALYS, STRESS (FILE4 = <model>.N4) ------------------------
    "FILE4": ({
        "node_id": ("i", "(nN,) node ids of the analysis model"),
        "node_xyz": ("f", "(nN, 3) global coordinates"),
        "eq_node": ("i", "(nEq,) node id of each equation"),
        "eq_dof": ("i", "(nEq,) dof 1..6 of each equation"),
        "int_node": ("i", "(nInt,) interaction node ids (interaction order)"),
        "int_iface": ("i", "(nInt,) 1-based user interface number of each interaction node"),
        "int_eq": ("i", "(nInt, 3) equation numbers of UX, UY, UZ (-1 if fixed/absent)"),
        "elem_group": ("i", "(nE,) group of each output-capable element"),
        "elem_id": ("i", "(nE,) element number within its group"),
        "elem_type": ("i", "(nE,) group type code (1 SOLID, 2 BEAMS, ...)"),
        "elem_nodes": ("i", "(nE, 8) node ids (0 = unused)"),
        "elem_excav": ("i", "(nE,) 1 excavated soil element, 0 structure"),
    }, ["gravity", "gelev", "dim", "nEq", "cmodform", "model_hash"]),
    # element recovery operators (STRESS) are stored in FILE4 per element type T as
    #   rec_<T>_S  (nE_T, nComp_T, nDofE_T) complex  -- component TF = S @ u_e
    #   rec_<T>_eq (nE_T, nDofE_T) int               -- equation numbers (-1 = fixed)
    #   rec_<T>_idx (nE_T,) int                      -- row in elem_* arrays
    #   meta["components"][T] = list of component codes (sassi.conventions.ELEMENT_COMPONENTS)
    "COOSK": ({
        # sparse (CSR) complex stiffness matrices on the FILE4 equation numbering
    }, ["nEq"]),          # sparse entries: Ks (structure incl. near field), Ke (excavated soil)
    "COOSM": ({}, ["nEq"]),  # sparse entries: Ms, Me (real)
    # ---------------- HOUSE -> ANALYS restart bookkeeping (requirements §2.1, §2.5) ----------
    "DOFSMAP": ({
        "eq_node": ("i", "(nEq,) node id of each equation (same as FILE4)"),
        "eq_dof": ("i", "(nEq,) dof 1..6"),
        "int_node": ("i", "(nInt,) interaction nodes in interaction order"),
    }, ["nEq"]),
    # FILE90: hashes used to validate restarts; meta keys int_hash (interaction-node ids, coordinates,
    # interfaces), layer_hash (site layer table + gravity + damping form), file4_hash (FILE4 + COOSK/COOSM)
    "FILE90": ({}, ["int_hash", "layer_hash", "file4_hash"]),
    # FILE91: free-form run metadata of the HOUSE run that produced FILE4 (model, title, counts, time)
    "FILE91": ({}, ["model"]),
    # COOXqqq: impedance X_ff (nInt*3, nInt*3) complex for one frequency; COOTKqqq: factorised
    # condensed system (layout private to ANALYS; meta keys fnum, freq and the FILE90 hashes)
    "COOX": ({"X": ("c", "(3nInt, 3nInt) impedance at interaction DOFs, node-major x,y,z")},
             ["fnum", "freq", "int_hash", "layer_hash"]),
    # ---------------- FORCE -> ANALYS -----------------------------------------------------
    "FILE9": ({
        "fnum": ("i", "(nF,)"),
        "freq": ("f", "(nF,)"),
        "load_node": ("i", "(nLd,)"),
        "load_dof": ("i", "(nLd,) 1..6"),
        "P": ("c", "(nF, nLd) complex load amplitude a*exp(-i w t_arrival)"),
    }, ["df"]),
    # ---------------- ANALYS / COMBIN -> MOTION, STRESS, COMBIN ----------------------------
    "FILE8": ({
        "fnum": ("i", "(nF,)"),
        "freq": ("f", "(nF,) Hz"),
        "eq_node": ("i", "(nEq,)"),
        "eq_dof": ("i", "(nEq,)"),
        "H": ("c", "(nF, nEq) transfer functions: seismic U/U_cp (dimensionless), vibration U per unit load factor"),
    }, ["df", "type", "case", "ang", "cm", "nfft", "delt", "model_hash"]),
    # ---------------- HOUSE -> ANALYS (incoherency / multiple excitation, P1) --------------
    "FILE77": ({
        "fnum": ("i", "(nF,)"),
        "freq": ("f", "(nF,)"),
        "int_node": ("i", "(nInt,)"),
        "s": ("c", "(nF, 3, nInt) motion factors per frequency, direction X/Y/Z and interaction node"),
    }, ["df", "method"]),
}


def validate(c: Container) -> List[str]:
    """Return a list of problems (empty if the container satisfies its schema)."""
    problems: List[str] = []
    spec = SCHEMAS.get(c.kind)
    if spec is None:
        return problems
    arrays, meta_keys = spec
    for name in arrays:
        if name not in c.arrays:
            problems.append(f"{c.kind}: missing array '{name}'")
    for k in meta_keys:
        if k not in c.meta:
            problems.append(f"{c.kind}: missing meta key '{k}'")
    if c.kind in ("FILE1", "FILE2", "FILE3", "FILE8", "FILE9", "FILE77") and "fnum" in c.arrays:
        nF = len(c.arrays["fnum"])
        for name in ("freq",):
            if name in c.arrays and len(c.arrays[name]) != nF:
                problems.append(f"{c.kind}: '{name}' length != nF")
    return problems
