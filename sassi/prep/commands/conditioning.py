"""Model conditioning: the FIXROT family (requirements section 3.4.I P0 items; spec 09 sections 2.7-2.10).

In ACS SASSI every node carries six DOFs, so DOFs that no element stiffens must be restrained:

* FIXSLDROT -- nodes connected only to SOLID elements get ROTX, ROTY, ROTZ fixed (``D,n,n,1,1,ROT``);
* FIXSHLROT,[stiff] -- the Kirchhoff SHELL has no in-plane ("drilling") rotational stiffness, so at
  every node connected only to coplanar shells a soft rotational spring about the shell normal
  ``n`` is added: a new coincident node with every DOF fixed, a SPRING element in a new group and an
  SC property with ``scxx = stiff nx^2``, ``scyy = stiff ny^2``, ``sczz = stiff nz^2`` (SC springs are
  global and uncoupled; exact ``stiff n n^T`` for axis-parallel shells).  Recommended stiffness:
  no more than 10 % of the shell bending stiffness (default 10);
* FIXSPRROT -- nodes connected only to springs (or springs and solids) get the DOFs without spring
  stiffness fixed (rotations; translations only where there is no mass either);
* FIXROT,[stiff] -- the combination: D on the rotations of solid-only nodes, D on the drilling
  rotation of shell-only nodes whose normal is a global axis, soft springs for oblique shells (they
  cannot be fixed with D), and FIXSPRROT for spring nodes.

TSHELL elements have a small rotational stiffness inside HOUSE and need no fix.  Nodes where
non-coplanar shells meet (wall-floor corners) are left alone: each shell's drilling rotation is
restrained by the other's bending.  The generation commands of section 3.4.I that are tier P1
keep the registry placeholders (``<CMD> is not available in this build (tier P1)``).
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np

from ...model.entities import Element, Group, SpringProp
from ..check import SHELL, SOLID, SPRING, ModelView
from ..registry import command

EPS = 1e-4          # coplanarity / axis-parallel tolerance on unit normals (spec 09 section 2.8)


def _solid_only(v: ModelView) -> List[int]:
    return sorted(n for n, lst in v.node_elems.items() if lst and all(r.type == SOLID for r in lst))


def _shell_only(v: ModelView) -> List[int]:
    return sorted(n for n, lst in v.node_elems.items() if lst and all(r.type == SHELL for r in lst))


def fix_solid_rotations(m, v: ModelView) -> int:
    """D on ROTX/ROTY/ROTZ of nodes connected only to SOLID elements; returns the node count."""
    k = 0
    for n in _solid_only(v):
        fx = m.nodes[n].fix
        if not all(fx[3:6]):
            fx[3:6] = [1, 1, 1]
            k += 1
    return k


class _SpringAdder:
    """Creates the drilling springs of FIXSHLROT / FIXROT (one group, one SC property per normal)."""

    def __init__(self, m, stiff: float, title: str):
        self.m = m
        self.stiff = float(stiff)
        self.title = title
        self.group: Optional[Group] = None
        self.props: Dict[Tuple[float, float, float], int] = {}
        self.count = 0
        self.next_node = (max(m.nodes) if m.nodes else 0) + 1

    def _prop(self, nrm: np.ndarray) -> int:
        k = self.stiff * nrm ** 2
        key = tuple(round(float(x), 12) for x in k)
        if key not in self.props:
            pid = (max(self.m.springs) if self.m.springs else 0) + 1
            self.m.springs[pid] = SpringProp(pid, 0.0, 0.0, 0.0, key[0], key[1], key[2], 0.0)
            self.props[key] = pid
        return self.props[key]

    def add(self, node: int, xyz: np.ndarray, nrm: np.ndarray) -> None:
        m = self.m
        if self.group is None:
            gid = (max(m.groups) if m.groups else 0) + 1
            self.group = Group(gid, SPRING, self.title)
            m.groups[gid] = self.group
        new = self.next_node
        self.next_node += 1
        m.define_node(new, [float(x) for x in xyz], 0)
        m.nodes[new].fix = [1] * 6
        eid = len(self.group.elements) + 1
        self.group.elements[eid] = Element(eid, [node, new], mat=1, prop=self._prop(nrm))
        self.count += 1


def shell_drilling(m, v: ModelView, stiff: float, use_d: bool, title: str) -> Tuple[int, int]:
    """Drilling restraints at coplanar shell-only nodes: D when ``use_d`` and the normal is a global
    axis, otherwise a soft spring.  Returns (nodes fixed with D, springs added)."""
    adder = _SpringAdder(m, stiff, title)
    nd = 0
    for n in _shell_only(v):
        nrm = v.coplanar_normal(n, EPS)
        if nrm is None:
            continue                       # non-coplanar shells restrain each other
        fx = m.nodes[n].fix
        axis = [k for k in range(3) if abs(nrm[k]) >= 1.0 - EPS]
        if use_d and axis:
            if not fx[3 + axis[0]]:
                fx[3 + axis[0]] = 1
                nd += 1
            continue
        if all(fx[3:6]) or (axis and fx[3 + axis[0]]):
            continue                       # already restrained
        adder.add(n, v.P[n], nrm)
    return nd, adder.count


def fix_spring_rotations(m, v: ModelView) -> int:
    """FIXSPRROT: nodes connected only to springs (or springs and solids): fix the DOFs that get no
    spring stiffness (rotations; spring-only translations only when they carry no mass)."""
    k = 0
    for n, lst in v.node_elems.items():
        types = {r.type for r in lst}
        if SPRING not in types or not types <= {SPRING, SOLID}:
            continue
        ksum = np.zeros(6)
        for r in lst:
            if r.type == SPRING:
                sc = m.springs.get(r.elem.prop)
                if sc is not None:
                    ksum += np.abs(np.asarray(sc.k, float))
        fx = m.nodes[n].fix
        dofs = [3, 4, 5]
        if types == {SPRING} and 0 not in m.nodes[n].flags:
            mass = m.tmass.get(n, [0.0, 0.0, 0.0])
            dofs = [d for d in range(3) if mass[d] == 0] + dofs
        changed = False
        for d in dofs:
            if ksum[d] == 0 and not fx[d]:
                fx[d] = 1
                changed = True
        k += int(changed)
    return k


def _stiff(c) -> float:
    s = c.float(1, default=10.0)
    if s <= 0:
        c.fail("<stiff> must be > 0")
    return s


@command("FIXSLDROT", max_args=0)
def cmd_fixsldrot(c):
    """FIXSLDROT: fix the rotations of nodes connected only to SOLID elements."""
    m = c.model
    n = fix_solid_rotations(m, ModelView(m))
    c.confirm(f"FIXSLDROT: rotations fixed at {n} solid-only nodes")


@command("FIXSHLROT", max_args=1)
def cmd_fixshlrot(c):
    """FIXSHLROT,[stiff]: soft drilling springs (default 10) at coplanar-shell nodes."""
    m = c.model
    stiff = _stiff(c)
    _, ns = shell_drilling(m, ModelView(m), stiff, use_d=False, title="FIXSHLROT drilling springs")
    c.confirm(f"FIXSHLROT: {ns} drilling springs of stiffness {stiff:g} added")
    if ns:
        c.info("FIXSHLROT: keep the spring stiffness below 10 % of the shell bending stiffness")


@command("FIXSPRROT", max_args=0)
def cmd_fixsprrot(c):
    """FIXSPRROT (alias FIXSPROT): fix the unstiffened DOFs of spring-only and spring+solid nodes."""
    m = c.model
    n = fix_spring_rotations(m, ModelView(m))
    c.confirm(f"FIXSPRROT: DOFs fixed at {n} spring nodes")


@command("FIXROT", max_args=1)
def cmd_fixrot(c):
    """FIXROT,[Stiff]: D on solid-only rotations and axis-parallel shell drilling rotations; soft
    springs (default 10) for oblique shells; FIXSPRROT for spring nodes."""
    m = c.model
    stiff = _stiff(c)
    v = ModelView(m)
    ns = fix_solid_rotations(m, v)
    nd, nsp = shell_drilling(m, v, stiff, use_d=True, title="FIXROT drilling springs")
    nspr = fix_spring_rotations(m, ModelView(m) if nsp else v)
    c.confirm(f"FIXROT: {ns} solid-only nodes, {nd} shell nodes fixed with D, {nsp} drilling springs "
              f"(stiffness {stiff:g}), {nspr} spring nodes")
