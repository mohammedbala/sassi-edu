"""ANSYS APDL export of a SASSI-EDU model: the ``ANSYS`` command (requirements 3.4.K, D-ANS-06;
spec 04 section 12.3; spec 09 section 6.1; R1 section 1.2 for the damping mapping).

The file ``<model>.inp`` is plain ``/PREP7`` input that ANSYS reads with ``/INPUT``: ET/KEYOPT, MP,
R/RMORE (or SECTYPE/SECDATA), N, TYPE/MAT/REAL/SECNUM + EN, D, CMBLOCK, with comments that say
where every item comes from.  It is meant for the cross-checks of requirements 6.6 (modal
frequencies of the fixed-base structure, static checks) and as a starting point for ANSYS work.

Element mapping (D-ANS-06; legacy ANSYS V11-15 elements by default, ``EDUOPT,ANSYSMODERN,1``
writes current elements):

==========================  ===============================  ======================================
SASSI-EDU                   legacy (default)                 modern
==========================  ===============================  ======================================
SOLID (structure)           SOLID45, KEYOPT(1) from MOPT     SOLID185 (KEYOPT(2) = 3 with
                            <incomp> (0 extra shapes)        incompatible modes, else 0)
PLANE (structure)           PLANE42 plane strain             PLANE182 plane strain
SHELL                       SHELL63, R = THICK               SHELL181 + SECTYPE,SHELL
TSHELL                      SHELL181 + SECTYPE,SHELL         SHELL181
BEAMS                       BEAM4 (no releases) / BEAM44     BEAM188 + SECTYPE ASEC (releases are
                            (KEYOPT(7)/(8) = KI/KJ)          not written: ENDRELEASE warning)
BEAMS with I2 = I3 = 0      LINK8 (axial)                    LINK180 (axial)
SPRING                      COMBIN14, KEYOPT(2) = 1..6, one element per non-zero SC component
GENERAL                     MATRIX27 KEYOPT(3) = 4 (K_R) and KEYOPT(3) = 2 (mass), global axes
MT / MR                     MASS21 KEYOPT(3) = 0 (mass units)
M (E, nu, weight, damping)  MP,EX / PRXY / DENS = weight/g / DMPR (see below)
D                           D (only the DOFs the node has in ANSYS; ALL when they are all fixed)
interaction nodes, groups   CMBLOCK node component SSI_INT, element components SSI_G<group>
excavated soil, L layers    not exported (structure only)
==========================  ===============================  ======================================

Beam axes: ANSYS puts the orientation node K in the element x-z plane (z towards K) and SASSI
puts it in the 1-2 plane (axis 2 towards K), so with the *same* K node ANSYS z = SASSI 2 and
ANSYS y = -SASSI 3: ``IZZ = I2``, ``IYY = I3``, ``SHEARZ = A/As2``, ``SHEARY = A/As3`` and the
releases map as UX UY UZ ROTX ROTY ROTZ = P1 P3 P2 M1 M3 M2 (inverse of :mod:`sassi.io.ansys_cdb`).

Damping (the ``<dmap>`` argument of the ANSYS command):

* 0 (default, D-ANS-06): ``MP,DMPR,mat,2*beta_s`` -- the structural damping coefficient
  g = 2 beta, i.e. the ``1 + 2 i beta`` form; E unchanged.
* 1 (R1 section 1.2, exact): ANSYS gives ``K (1 + i g)``; to reproduce the SASSI factor
  ``1 - 2 b^2 + 2 i b sqrt(1 - b^2)`` write ``E_ANSYS = E (1 - 2 b^2)`` and
  ``g = 2 b sqrt(1 - b^2) / (1 - 2 b^2)``, so that ``E_ANSYS (1 + i g) = E c(b)`` exactly.  A model
  with ``CMODFORM,1`` uses ``c = 1 + 2 i b`` (D-CNV-03), whose exact mapping is ``E_ANSYS = E``,
  ``g = 2 b``.
* 2: no damping data.

Mass is always real; SASSI uses 1/2 lumped + 1/2 consistent SOLID mass and lumped SHELL mass,
ANSYS consistent mass by default (``LUMPM``): small differences in modal frequencies are expected
from the mass schemes, not from the converter.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from .. import PRODUCT, __version__
from ..model import SSIModel
from ..model.materials import elastic_constants
from ..model.values import fmt_num

DOF_LABELS = ("UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ")
DAMPING_MODES = {0: "ratio", 1: "exact", 2: "none"}

# SASSI group type codes
SOLID, BEAMS, SHELL, PLANE, TSHELL, SPRING, GENERAL = 1, 2, 3, 4, 5, 7, 9


@dataclass
class ApdlOptions:
    """Export options: ``modern`` elements (EDUOPT,ANSYSMODERN,1), ``damping`` mapping
    ('ratio' = DMPR 2 beta, 'exact' = R1 1.2, 'none'), element components per group."""
    modern: bool = False
    damping: str = "ratio"
    group_components: bool = True


@dataclass
class ApdlExport:
    """Result of :func:`export_apdl`: the file text, messages and the element map."""
    text: str
    warnings: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)
    element_map: List[Tuple[int, int, int]] = field(default_factory=list)   # (ANSYS elem, group, elem)
    counts: Dict[str, int] = field(default_factory=dict)


def ansys_damping(E: float, beta: float, mode: str = "ratio", cmodform: int = 0) -> Tuple[float, float]:
    """``(E_ANSYS, DMPR)`` of a material with Young's modulus E and damping ratio beta.

    'ratio' (D-ANS-06): (E, 2 beta).  'exact' (R1 1.2): ANSYS structural damping multiplies the
    stiffness by ``1 + i g``, so ``E_A (1 + i g)`` is made equal to the model's complex modulus
    ``E c(beta)``: with the SASSI form (``cmodform`` 0, D-CNV-03) ``c = 1 - 2b^2 + 2ib sqrt(1-b^2)``,
    hence ``E_A = E (1 - 2 b^2)`` and ``g = 2 b sqrt(1 - b^2) / (1 - 2 b^2)``; with CMODFORM,1
    ``c = 1 + 2ib``, hence ``E_A = E`` and ``g = 2 b`` (the same numbers as 'ratio').
    'none': (E, 0).
    """
    if mode == "none" or beta == 0.0:
        return E, 0.0
    if mode == "exact" and cmodform == 0:
        c = 1.0 - 2.0 * beta * beta
        if c <= 0.0:
            raise ValueError(f"damping ratio {beta} too large for the exact mapping (1 - 2 b^2 <= 0)")
        return E * c, 2.0 * beta * math.sqrt(1.0 - beta * beta) / c
    return E, 2.0 * beta


def _n(x) -> str:
    return fmt_num(float(x))


def _release_keyopt(code: Sequence[int]) -> int:
    """SASSI KI/KJ (P1 P2 P3 M1 M2 M3) -> BEAM44 KEYOPT(7)/(8) digits UX UY UZ ROTX ROTY ROTZ."""
    p = [1 if v else 0 for v in code]
    digits = [p[0], p[2], p[1], p[3], p[5], p[4]]
    return int("".join(str(d) for d in digits))


def cmblock_lines(name: str, kind: str, ids: Sequence[int]) -> List[str]:
    """CMBLOCK (blocked component definition, format (8i10)); ranges as ``first, -last``."""
    ids = sorted(set(int(i) for i in ids))
    items: List[int] = []
    i = 0
    while i < len(ids):
        j = i
        while j + 1 < len(ids) and ids[j + 1] == ids[j] + 1:
            j += 1
        items.append(ids[i])
        if j > i:
            items.append(-ids[j])
        i = j + 1
    out = [f"CMBLOCK,{name},{kind},{len(items):8d}", "(8i10)"]
    for k in range(0, len(items), 8):
        out.append("".join(f"{v:10d}" for v in items[k:k + 8]))
    return out


def _r_lines(setno: int, values: Sequence[float], comment: str = "") -> List[str]:
    """R,set,v1..v6 then RMORE lines of 6 values each (R7..R12, R13..R18, ...)."""
    vals = [float(v) for v in values]
    first = vals[:6]
    out = ["R," + ",".join([str(setno)] + [_n(v) for v in first]) + (f"   ! {comment}" if comment else "")]
    rest = vals[6:]
    while rest:
        chunk = rest[:6]
        rest = rest[6:]
        out.append("RMORE," + ",".join(_n(v) for v in chunk))
    return out


class _Exporter:
    def __init__(self, model: SSIModel, opts: ApdlOptions, etype_resolver: Optional[Callable] = None):
        self.m = model
        self.o = opts
        self.res = ApdlExport(text="")
        self.etype_of = etype_resolver or _default_resolver(model)
        self.g = model.gravity
        rec = model.options.record("CMODFORM")           # complex-modulus form of the model (D-CNV-03)
        self.cmodform = (rec.integer(1, 0) or 0) if rec is not None else 0
        self.et_lines: List[str] = []
        self.r_lines: List[str] = []
        self.el_lines: List[str] = []
        self.n_et = 0
        self.n_real = 0
        self.n_sec = 0
        self.n_el = 0
        self.used_mats: Dict[int, str] = {}
        self.node_dofs: Dict[int, set] = {}
        self.group_elems: Dict[int, List[int]] = {}
        self.state = {"TYPE": None, "MAT": None, "REAL": None, "SECNUM": None}
        self.warned: set = set()
        self.real_cache: Dict[Tuple, int] = {}
        self.sec_cache: Dict[Tuple, int] = {}
        ids, xyz = model.global_coordinates()
        self.ids = [int(i) for i in ids]
        self.P = {int(i): xyz[k] for k, i in enumerate(ids)}
        types = {g.type for g in model.groups.values() if g.elements}
        self.two_d = PLANE in types
        # x_A = R x_S: 2-D models go from the SASSI X-Z plane to the ANSYS X-Y plane (inverse of
        # the CONVERT,ANSYS rotation); 3-D models are unchanged
        self.R = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]]) if self.two_d else np.eye(3)
        self.axis = [int(np.argmax(np.abs(self.R[:, s]))) for s in range(3)]   # SASSI axis -> ANSYS axis
        if self.two_d and types - {PLANE, SPRING, GENERAL}:
            self.warn("2-D model with PLANE elements and 3-D element types: every node is rotated from the "
                      "SASSI X-Z plane to the ANSYS X-Y plane")

    # ------------------------------------------------------------------ messages
    def warn(self, text: str, key: Optional[str] = None) -> None:
        k = key or text
        if k not in self.warned:
            self.warned.add(k)
            self.res.warnings.append(text)

    def note(self, text: str) -> None:
        self.res.notes.append(text)

    # ------------------------------------------------------------------ helpers
    def ansys_dof(self, s: int) -> int:
        """ANSYS DOF index of SASSI DOF index ``s`` (0..5) under the 2-D rotation."""
        return self.axis[s % 3] + (3 if s >= 3 else 0)

    def A(self, n: int) -> np.ndarray:
        return self.R @ self.P[n]

    def new_et(self, ename: str, keyopts: Sequence[Tuple[int, int, str]], comment: str) -> int:
        self.n_et += 1
        t = self.n_et
        self.et_lines.append(f"ET,{t},{ename}" + " " * max(1, 22 - len(f"ET,{t},{ename}")) + f"! {comment}")
        for k, v, why in keyopts:
            line = f"KEYOPT,{t},{k},{v}"
            self.et_lines.append(line + (" " * max(1, 22 - len(line)) + f"! {why}" if why else ""))
        return t

    def real_set(self, key: Tuple, values: Sequence[float], comment: str) -> int:
        rid = self.real_cache.get(key)
        if rid is None:
            self.n_real += 1
            rid = self.n_real
            self.real_cache[key] = rid
            self.r_lines += _r_lines(rid, values, comment)
        return rid

    def attr(self, **kw) -> None:
        for k, v in kw.items():
            if v is not None and self.state.get(k) != v:
                self.el_lines.append(f"{k},{v}")
                self.state[k] = v

    def element(self, gid: int, eid: int, nodes: Sequence[int], dofs: Sequence[int],
                n_dof_nodes: Optional[int] = None) -> int:
        """``EN`` line; the first ``n_dof_nodes`` nodes (all by default) get ``dofs`` (a beam K node
        is orientation only and gets none, so no D is written for it)."""
        self.n_el += 1
        self.el_lines.append("EN," + ",".join(str(int(x)) for x in [self.n_el] + list(nodes)))
        self.res.element_map.append((self.n_el, gid, eid))
        self.group_elems.setdefault(gid, []).append(self.n_el)
        k = len(nodes) if n_dof_nodes is None else n_dof_nodes
        for n in list(nodes)[:k]:
            self.node_dofs.setdefault(int(n), set()).update(dofs)
        return self.n_el

    def mat(self, mid: int, by: str) -> int:
        if mid not in self.m.materials:
            self.warn(f"material {mid} (used by {by} elements) is not defined; MAT,{mid} written without MP data",
                      key=f"nomat{mid}")
        self.used_mats.setdefault(mid, by)
        return mid

    # ------------------------------------------------------------------ groups
    def solid(self, g, elems) -> None:
        inc = int(self.m.mopt.get("incomp")) == 0
        if self.o.modern:
            t = self.new_et("SOLID185", [(2, 3 if inc else 0, "simplified enhanced strain (incompatible modes)"
                                          if inc else "full integration (B-bar)")],
                            f"group {g.id} SOLID ({g.title})" if g.title else f"group {g.id} SOLID")
        else:
            t = self.new_et("SOLID45", [(1, 0 if inc else 1, "extra displacement shapes " +
                                         ("included (MOPT <incomp> = 0)" if inc else "suppressed (MOPT <incomp> = 1)"))],
                            f"group {g.id} SOLID" + (f" ({g.title})" if g.title else ""))
        dofs = (0, 1, 2) if not self.two_d else (0, 1)
        if any(e.eint for e in elems):
            self.note(f"group {g.id}: SASSI EINT integration orders are not exported (ANSYS uses 2x2x2)")
        for e in elems:
            self.attr(TYPE=t, MAT=self.mat(e.mat, "SOLID"))
            self.element(g.id, e.id, (list(e.nodes) + [0] * 8)[:8], dofs)
        self.res.counts["SOLID"] = self.res.counts.get("SOLID", 0) + len(elems)

    def plane(self, g, elems) -> None:
        inc = int(self.m.mopt.get("incomp")) == 0
        if self.o.modern:
            t = self.new_et("PLANE182", [(1, 3 if inc else 0, "simplified enhanced strain" if inc else "full integration"),
                                         (3, 2, "plane strain")], f"group {g.id} PLANE")
        else:
            t = self.new_et("PLANE42", [(2, 0 if inc else 1, "extra displacement shapes " +
                                         ("included" if inc else "suppressed")), (3, 2, "plane strain")],
                            f"group {g.id} PLANE")
        for e in elems:
            n = [x for x in e.nodes if x]
            if len(n) == 3:
                n = n + [n[2]]
            self.attr(TYPE=t, MAT=self.mat(e.mat, "PLANE"))
            self.element(g.id, e.id, n[:4], (0, 1))
        self.res.counts["PLANE"] = self.res.counts.get("PLANE", 0) + len(elems)

    def shell(self, g, elems) -> None:
        use181 = self.o.modern or g.type == TSHELL
        if g.type == TSHELL:
            self.warn("TSHELL groups exported as SHELL181 (thick shell)")
        if use181:
            t = self.new_et("SHELL181", [(1, 0, "bending and membrane"), (3, 2, "full integration with "
                                                                          "incompatible modes")],
                            f"group {g.id} {'TSHELL' if g.type == TSHELL else 'SHELL'} (Mindlin in ANSYS; "
                            "SASSI SHELL is Kirchhoff)")
        else:
            t = self.new_et("SHELL63", [], f"group {g.id} SHELL (Kirchhoff, R1 = THICK)")
        for e in elems:
            n = [x for x in e.nodes if x]
            if len(n) == 3:
                n = n + [n[2]]
            if e.thick <= 0.0:
                self.warn(f"group {g.id}: shell elements without THICK exported with thickness 0", key=f"thk{g.id}")
            mid = self.mat(e.mat, "SHELL")
            if use181:
                key = ("SHELL", float(e.thick), mid)
                sid = self.sec_cache.get(key)
                if sid is None:
                    self.n_sec += 1
                    sid = self.n_sec
                    self.sec_cache[key] = sid
                    self.r_lines += [f"SECTYPE,{sid},SHELL" + " " * 8 + f"! thickness {_n(e.thick)}, material {mid}",
                                     f"SECDATA,{_n(e.thick)},{mid},0,3", "SECOFFSET,MID"]
                self.attr(TYPE=t, MAT=mid, SECNUM=sid)
            else:
                rid = self.real_set(("SHELL63", float(e.thick)), [e.thick], f"SHELL63 thickness {_n(e.thick)}")
                self.attr(TYPE=t, MAT=mid, REAL=rid)
            self.element(g.id, e.id, n[:4], range(6))
        self.res.counts["SHELL"] = self.res.counts.get("SHELL", 0) + len(elems)

    def beams(self, g, elems) -> None:
        # one ANSYS type per (kind, release pattern): KEYOPT(7)/(8) are type properties (ANSYSREFORMAT)
        patterns: Dict[Tuple, List] = {}
        for e in elems:
            s = self.m.sections.get(e.prop)
            axial = s is not None and s.flex2 == 0.0 and s.flex3 == 0.0
            rel = (tuple(int(bool(v)) for v in e.ki), tuple(int(bool(v)) for v in e.kj))
            patterns.setdefault((axial, rel), []).append(e)
        rel_patterns = {k[1] for k in patterns if not k[0]}
        if len(rel_patterns) > 1 and not self.o.modern:
            self.warn(f"group {g.id}: {len(rel_patterns)} beam release patterns written as {len(rel_patterns)} ANSYS "
                      "element types (run ANSYSREFORMAT to regroup the beams by release pattern)", key=f"rel{g.id}")
        for (axial, rel), els in sorted(patterns.items(), key=lambda kv: (kv[0][0], kv[0][1])):
            released = any(rel[0]) or any(rel[1])
            if axial:
                ename = "LINK180" if self.o.modern else "LINK8"
                t = self.new_et(ename, [], f"group {g.id} BEAMS with I2 = I3 = 0 (axial only)")
            elif self.o.modern:
                t = self.new_et("BEAM188", [(3, 3, "cubic shape functions (one element per SASSI beam)")],
                                f"group {g.id} BEAMS")
                if released:
                    self.warn(f"group {g.id}: beam end releases are not written for BEAM188 (ENDRELEASE adds couplings); "
                              "use the legacy export (EDUOPT,ANSYSMODERN,0) for BEAM44 releases", key=f"er{g.id}")
            elif released:
                k7, k8 = _release_keyopt(rel[0]), _release_keyopt(rel[1])
                t = self.new_et("BEAM44", [(7, k7, f"releases at I (UX UY UZ ROTX ROTY ROTZ) from KI "
                                            f"{''.join(map(str, rel[0]))}"),
                                           (8, k8, f"releases at J from KJ {''.join(map(str, rel[1]))}")],
                                f"group {g.id} BEAMS with end releases")
            else:
                t = self.new_et("BEAM4", [], f"group {g.id} BEAMS")
            for e in els:
                s = self.m.sections.get(e.prop)
                if s is None:
                    self.warn(f"group {g.id}: beam section R,{e.prop} is not defined; elements not exported",
                              key=f"nor{e.prop}")
                    continue
                mid = self.mat(e.mat, "BEAMS")
                n = (list(e.nodes) + [0, 0, 0])[:3]
                if axial:
                    if self.o.modern:
                        key = ("LINK", e.prop)
                        sid = self.sec_cache.get(key)
                        if sid is None:
                            self.n_sec += 1
                            sid = self.sec_cache[key] = self.n_sec
                            self.r_lines += [f"SECTYPE,{sid},LINK" + " " * 8 + f"! R,{e.prop} (axial only)",
                                             f"SECDATA,{_n(s.axial)}"]
                        self.attr(TYPE=t, MAT=mid, SECNUM=sid)
                    else:
                        rid = self.real_set(("LINK8", e.prop), [s.axial], f"LINK8 AREA from R,{e.prop}")
                        self.attr(TYPE=t, MAT=mid, REAL=rid)
                    self.element(g.id, e.id, n[:2], (0, 1, 2))
                    continue
                tkz = math.sqrt(12.0 * s.flex3 / s.axial) if s.axial > 0 else 0.0     # rectangle with IYY = I3
                tky = math.sqrt(12.0 * s.flex2 / s.axial) if s.axial > 0 else 0.0
                shz = s.axial / s.shear2 if s.shear2 > 0 else 0.0
                shy = s.axial / s.shear3 if s.shear3 > 0 else 0.0
                if self.o.modern:
                    G = self._shear_modulus(mid)
                    key = ("ASEC", e.prop, mid if (s.shear2 > 0 or s.shear3 > 0) else 0)
                    sid = self.sec_cache.get(key)
                    if sid is None:
                        self.n_sec += 1
                        sid = self.sec_cache[key] = self.n_sec
                        self.r_lines += [f"SECTYPE,{sid},BEAM,ASEC" + " " * 4 + f"! R,{e.prop}: Iyy = I3, Izz = I2",
                                         "SECDATA," + ",".join(_n(v) for v in (s.axial, s.flex3, 0.0, s.flex2, 0.0,
                                                                               s.tors, 0.0, 0.0, 0.0, 0.0, tkz, tky))]
                        if (s.shear2 > 0 or s.shear3 > 0) and G > 0:
                            self.r_lines.append(f"SECCONTROL,{_n(G * s.shear2)},{_n(G * s.shear3)}"
                                                "   ! transverse shear stiffness k G A (xz, xy) = G As2, G As3")
                    self.attr(TYPE=t, MAT=mid, SECNUM=sid)
                elif released:
                    vals = [s.axial, s.flex2, s.flex3, tkz, tky, s.tors] * 2 + [0.0] * 6 + [shz, shy]
                    rid = self.real_set(("BEAM44", e.prop), vals, f"BEAM44 from R,{e.prop}: IZ = I2, IY = I3, IX = J")
                    self.attr(TYPE=t, MAT=mid, REAL=rid)
                else:
                    vals = [s.axial, s.flex2, s.flex3, tkz, tky, 0.0, 0.0, s.tors, shz, shy]
                    rid = self.real_set(("BEAM4", e.prop), vals,
                                        f"BEAM4 from R,{e.prop}: IZZ = I2, IYY = I3, IXX = J, SHEAR = A/As")
                    self.attr(TYPE=t, MAT=mid, REAL=rid)
                if not n[2]:
                    self.warn(f"group {g.id} element {e.id}: no K node", key=f"nok{g.id}")
                self.element(g.id, e.id, n if n[2] else n[:2], range(6), n_dof_nodes=2)
            self.res.counts["BEAMS"] = self.res.counts.get("BEAMS", 0) + len(els)

    def _shear_modulus(self, mid: int) -> float:
        mt = self.m.materials.get(mid)
        if mt is None:
            return 0.0
        try:
            c = elastic_constants(mt.mtype, mt.val1, mt.val2, mt.weight, self.g)
        except ValueError:
            return 0.0
        E, _ = ansys_damping(c.E, mt.sdamp, self.o.damping, self.cmodform)
        return E / (2.0 * (1.0 + c.nu))

    def springs(self, g, elems) -> None:
        types: Dict[int, int] = {}
        damped = False
        for e in elems:
            sp = self.m.springs.get(e.prop)
            if sp is None:
                self.warn(f"group {g.id}: spring property SC,{e.prop} is not defined; elements not exported",
                          key=f"nosc{e.prop}")
                continue
            damped = damped or sp.damp != 0.0
            n = [x for x in e.nodes if x][:2]
            for s, k in enumerate(sp.k):
                if k == 0.0:
                    continue
                a = self.ansys_dof(s)
                t = types.get(a)
                if t is None:
                    t = types[a] = self.new_et("COMBIN14", [(2, a + 1, f"1-D spring in {DOF_LABELS[a]}")],
                                               f"group {g.id} SPRING component {DOF_LABELS[s]} (SASSI)")
                rid = self.real_set(("C14", float(k)), [k], f"COMBIN14 K = SC {DOF_LABELS[s]}")
                self.attr(TYPE=t, REAL=rid)
                self.element(g.id, e.id, n, (a,))
        if damped:
            self.warn("SC spring damping ratios have no frequency-independent COMBIN14 counterpart; not exported "
                      "(use material or global structural damping in ANSYS)", key="scdamp")
        self.res.counts["SPRING"] = self.res.counts.get("SPRING", 0) + len(elems)

    def general(self, g, elems) -> None:
        from ..elements.base import blas_quiet, frame_from_three_points
        tK = tM = None
        TR = np.kron(np.eye(4), self.R)
        for e in elems:
            if e.prop not in self.m.matrices:
                self.warn(f"group {g.id}: matrix property {e.prop} is not defined; elements not exported",
                          key=f"nomx{e.prop}")
                continue
            K, M = self.m.general_matrices(e.prop)
            if np.any(K.imag != 0.0):
                self.warn("GENERAL imaginary stiffness (MXI) has no MATRIX27 counterpart; only K_R is exported",
                          key="mxi")
            Kr, Mm = K.real, M
            n = [x for x in e.nodes if x]
            local = len(n) >= 3
            with blas_quiet():
                if local:
                    lam, _ = frame_from_three_points(self.P[n[0]], self.P[n[1]], self.P[n[2]])
                    T = np.kron(np.eye(4), lam)
                    Kr, Mm = T.T @ Kr @ T, T.T @ Mm @ T
                if self.two_d:
                    Kr, Mm = TR @ Kr @ TR.T, TR @ Mm @ TR.T
            for kind, A in (("stiffness", Kr), ("mass", Mm)):
                if not np.any(A != 0.0):
                    continue
                if kind == "stiffness" and tK is None:
                    tK = self.new_et("MATRIX27", [(3, 4, "stiffness matrix (K_R)")], f"group {g.id} GENERAL stiffness")
                if kind == "mass" and tM is None:
                    tM = self.new_et("MATRIX27", [(3, 2, "mass matrix (mass units)")], f"group {g.id} GENERAL mass")
                vals = []
                for i in range(12):
                    vals.extend(float(v) for v in A[i, i:])
                key = ("M27", kind, e.prop, tuple(n) if (local or self.two_d) else ())
                rid = self.real_set(key, vals, f"MATRIX27 {kind} of matrix property {e.prop}"
                                    + (" (local axes rotated to global)" if local else ""))
                self.attr(TYPE=tK if kind == "stiffness" else tM, REAL=rid)
                self.element(g.id, e.id, n[:2], range(6), n_dof_nodes=2)
        self.res.counts["GENERAL"] = self.res.counts.get("GENERAL", 0) + len(elems)

    def masses(self) -> None:
        nm = self.m.nodal_masses()
        if not nm:
            return
        if self.two_d:
            t = self.new_et("MASS21", [(3, 4, "2-D mass without rotary inertia")], "nodal masses MT (mass units)")
        else:
            t = self.new_et("MASS21", [(3, 0, "3-D mass with rotary inertia")], "nodal masses MT/MR (mass units)")
        for n in sorted(nm):
            v = nm[n]
            if not np.any(v != 0.0):
                continue
            if self.two_d:
                vals = [float(v[0])]
                if v[0] != v[2] or v[1] != 0.0:
                    self.warn("2-D MASS21 carries one mass for X and Y: the SASSI X mass is used", key="m2d")
                dofs = (0, 1)
            else:
                vals = [float(x) for x in v]
                dofs = range(6)
            rid = self.real_set(("MASS21", tuple(vals)), vals, f"MASS21 at node {n}")
            self.attr(TYPE=t, REAL=rid)
            self.n_el += 1
            self.el_lines.append(f"EN,{self.n_el},{n}")
            self.node_dofs.setdefault(n, set()).update(dofs)
        self.res.counts["MASS21"] = sum(1 for v in nm.values() if np.any(v != 0.0))

    # ------------------------------------------------------------------ materials, constraints
    def materials(self) -> List[str]:
        out: List[str] = []
        mode = self.o.damping
        for mid in sorted(self.used_mats):
            mt = self.m.materials.get(mid)
            if mt is None:
                continue
            try:
                c = elastic_constants(mt.mtype, mt.val1, mt.val2, mt.weight, self.g)
            except ValueError as exc:
                self.warn(f"material {mid}: {exc}; not exported")
                continue
            beta = mt.sdamp
            if mt.pdamp != mt.sdamp and mode != "none":
                self.warn(f"material {mid}: beta_p = {_n(mt.pdamp)} != beta_s = {_n(mt.sdamp)}; ANSYS has one "
                          "damping value per material, beta_s is used", key=f"bp{mid}")
            E, dmpr = ansys_damping(c.E, beta, mode, self.cmodform)
            exact0 = mode == "exact" and self.cmodform == 0
            out.append(f"! material {mid}: type {mt.mtype}, weight {_n(mt.weight)}, beta_p {_n(mt.pdamp)}, "
                       f"beta_s {_n(mt.sdamp)}")
            out.append(f"MP,EX,{mid},{_n(E)}" + (f"   ! E (1 - 2 beta^2), E = {_n(c.E)}" if exact0 and beta else ""))
            out.append(f"MP,PRXY,{mid},{_n(c.nu)}")
            out.append(f"MP,DENS,{mid},{_n(c.rho)}   ! weight / g")
            if dmpr:
                why = ("g = 2 beta sqrt(1 - beta^2)/(1 - 2 beta^2) (R1 1.2)" if exact0
                       else "g = 2 beta: exact for CMODFORM,1 (E (1 + 2 i beta))" if mode == "exact"
                       else "2 beta (D-ANS-06)")
                out.append(f"MP,DMPR,{mid},{_n(dmpr)}   ! {why}")
        return out

    def constraints(self) -> List[str]:
        out: List[str] = []
        for n in self.ids:
            nd = self.m.nodes[n]
            if not any(nd.fix):
                continue
            have = self.node_dofs.get(n, set())
            if not have:
                continue
            fixed = {self.ansys_dof(s) for s in range(6) if nd.fix[s]} & have
            if not fixed:
                continue
            if all(nd.fix):
                out.append(f"D,{n},ALL,0")         # ALL = every DOF the node has in ANSYS
            else:
                out.extend(f"D,{n},{DOF_LABELS[a]},0" for a in sorted(fixed))
        return out

    # ------------------------------------------------------------------ driver
    def run(self) -> ApdlExport:
        m = self.m
        n_exc = 0
        for gid in sorted(m.groups):
            g = m.groups[gid]
            elems = g.sorted_elements()
            if not elems:
                continue
            if g.type in (SOLID, PLANE):
                keep = [e for e in elems if self.etype_of(g, e) != 2]
                n_exc += len(elems) - len(keep)
                elems = keep
                if not elems:
                    continue
            fn = {SOLID: self.solid, PLANE: self.plane, SHELL: self.shell, TSHELL: self.shell, BEAMS: self.beams,
                  SPRING: self.springs, GENERAL: self.general}.get(g.type)
            if fn is None:
                self.warn(f"group {gid}: element type {g.type} not exported")
                continue
            fn(g, elems)
        if n_exc:
            self.note(f"{n_exc} excavated-soil SOLID/PLANE elements and the soil layers are not exported "
                      "(structure only, D-ANS-06)")
        self.masses()
        mp = self.materials()
        dl = self.constraints()
        return self._assemble(mp, dl)

    def _assemble(self, mp: List[str], dl: List[str]) -> ApdlExport:
        m = self.m
        modern = "modern (SOLID185, SHELL181, BEAM188, LINK180)" if self.o.modern else \
            "legacy ANSYS V11-15 (SOLID45, SHELL63, BEAM4/BEAM44, LINK8, PLANE42)"
        damp = {"ratio": "MP,DMPR = 2 beta_s (D-ANS-06: structural damping g = 2 beta)",
                "exact": ("EX = E (1 - 2 beta^2), DMPR = 2 beta sqrt(1 - beta^2)/(1 - 2 beta^2) (R1 1.2, exact "
                          "SASSI complex modulus)" if self.cmodform == 0 else
                          "EX = E, DMPR = 2 beta (exact for the model's CMODFORM,1 modulus E (1 + 2 i beta))"),
                "none": "no damping written"}[self.o.damping]
        L: List[str] = [
            "! " + "=" * 92,
            f"! {PRODUCT} {__version__}: ANSYS APDL export of model '{m.name or 'unnamed'}' (ANSYS command, D-ANS-06)",
            f"! Title    : {m.title}" if m.title else "! Title    : (none)",
            f"! Elements : {modern}",
            f"! Units    : the model's consistent units; DENS = weight / g with g = {_n(self.g)}",
            f"! Damping  : {damp}",
            "!            ANSYS releases differ in the meaning of DMPR (structural damping coefficient vs damping",
            "!            ratio) and newer releases offer MP,DMPS: check your release (docs/user/ANSYS.md)",
            "! Beams    : the SASSI K node is the ANSYS K node (ANSYS z = SASSI axis 2, ANSYS y = -SASSI axis 3)",
        ]
        if self.two_d:
            L.append("! 2-D      : SASSI X-Z plane -> ANSYS X-Y plane (X = x, Y = z), PLANE elements in plane strain")
        L += ["! " + "=" * 92, "/PREP7"]
        if m.title:
            L.append(f"/TITLE,{m.title}")
        L += ["!", "! ---- Element types"] + self.et_lines
        if mp:
            L += ["!", "! ---- Materials (M table)"] + mp
        if self.r_lines:
            L += ["!", "! ---- Real constants and sections"] + self.r_lines
        L += ["!", f"! ---- Nodes ({len(self.ids)}, global coordinates)"]
        for n in self.ids:
            x = self.A(n) + 0.0
            L.append(f"N,{n},{_n(x[0])},{_n(x[1])},{_n(x[2])}")
        if self.el_lines:
            L += ["!", f"! ---- Elements ({self.n_el}); ANSYS element -> SASSI group/element:"]
            for gid in sorted(self.group_elems):
                ids = self.group_elems[gid]
                L.append(f"!      group {gid}: ANSYS elements {_ranges(ids)}")
            L += self.el_lines
        if dl:
            L += ["!", "! ---- Constraints (D)"] + dl
        comp: List[str] = []
        inter = [n for n in self.ids if 0 in m.nodes[n].flags]
        if inter:
            comp += ["! interaction nodes (INT code 0)"] + cmblock_lines("SSI_INT", "NODE", inter)
        if self.o.group_components:
            for gid in sorted(self.group_elems):
                comp += [f"! elements of SASSI group {gid}"] + cmblock_lines(f"SSI_G{gid}", "ELEM",
                                                                               self.group_elems[gid])
        if comp:
            L += ["!", "! ---- Components"] + comp
        L += ["!", "FINISH", ""]
        self.res.text = "\n".join(L)
        self.res.counts["nodes"] = len(self.ids)
        self.res.counts["elements"] = self.n_el
        return self.res


def _ranges(ids: Sequence[int]) -> str:
    ids = sorted(ids)
    parts = []
    i = 0
    while i < len(ids):
        j = i
        while j + 1 < len(ids) and ids[j + 1] == ids[j] + 1:
            j += 1
        parts.append(str(ids[i]) if i == j else f"{ids[i]}-{ids[j]}")
        i = j + 1
    return ", ".join(parts)


def _default_resolver(model: SSIModel) -> Callable:
    """ETYPE of SOLID/PLANE elements resolved with the CHECK rule D-HOU-01 (ETYPE 0 above grade is
    structure, D-ANS-06).  Imported lazily: :mod:`sassi.prep.check` holds the authoritative rule."""
    from ..prep.check import ModelView
    view = ModelView(model)
    table = {(r.group, r.elem.id): r.etype for r in view.elems}

    def resolve(g, e) -> int:
        return table.get((g.id, e.id), 1)
    return resolve


def export_apdl(model: SSIModel, options: Optional[ApdlOptions] = None,
                etype_resolver: Optional[Callable] = None) -> ApdlExport:
    """APDL text of ``model`` (see the module docstring for the mapping)."""
    opts = options or ApdlOptions()
    if opts.damping not in ("ratio", "exact", "none"):
        raise ValueError(f"damping mapping {opts.damping!r} is not 'ratio', 'exact' or 'none'")
    return _Exporter(model, opts, etype_resolver).run()


def write_apdl(path: Union[str, Path], model: SSIModel, options: Optional[ApdlOptions] = None) -> ApdlExport:
    """Write :func:`export_apdl` to ``path``."""
    res = export_apdl(model, options)
    Path(path).write_text(res.text, encoding="utf-8")
    return res
