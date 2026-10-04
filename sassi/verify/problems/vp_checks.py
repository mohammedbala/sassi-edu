"""VP-52: model checks with planted faults (requirements section 6.4; spec 01 section 13 item 6;
spec 09 section 3; D-CHK-06, D-CHK-09).

Each check is run on a small model with exactly one planted fault and on the same model without
it (no false positive):

* shared interior excavation node -> EXCSTRCHK (G-15);
* fixed interaction node -> FIXEDINT, and CHECK Error 124 when all translations are fixed;
* beam meeting solids at one node -> HINGED;
* beam K node that is an interaction node -> KINT;
* free spring node -> FREESPRING;
* oblique Kirchhoff shell without FIXROT -> Warning EDU-06; FIXROT removes it with drilling
  springs, and on a flat XY plate FIXROT fixes ROTZ only (spec 09 section 8 item 8);
* TOPL naming an undefined soil layer -> Error 19 under SITE and HOUSE, which both copy the layer
  profile (D-CHK-01), blocking both decks but not POINT.

The reference values are exact counts and node numbers of the planted faults.
"""
from __future__ import annotations

import math
from typing import List

from .. import VPResult, problem


def _block(etype: int, z0: float = -10.0, h: float = 5.0) -> List[str]:
    """2 x 2 x 2 hex block, nodes 1..27 numbered bottom-up (id = 1 + i + 3 j + 9 k)."""
    out = []
    for k in range(3):
        for j in range(3):
            for i in range(3):
                out.append(f"N,{1 + i + 3 * j + 9 * k},{5.0 * i},{5.0 * j},{z0 + h * k}")

    def n(i, j, k):
        return 1 + i + 3 * j + 9 * k
    out.append("GROUP,1,SOLID")
    e = 0
    for k in range(2):
        for j in range(2):
            for i in range(2):
                e += 1
                ids = [n(i, j, k), n(i + 1, j, k), n(i + 1, j + 1, k), n(i, j + 1, k),
                       n(i, j, k + 1), n(i + 1, j, k + 1), n(i + 1, j + 1, k + 1), n(i, j + 1, k + 1)]
                out.append("E," + ",".join(str(x) for x in [e] + ids))
    out.append(f"ETYPE,1,8,1,{etype}")
    return out


SOIL = ["L,1,5,0.12,1500,800,0.05,0.05", "L,2,10,0.13,2400,1200,0.02,0.02", "TOPL,1,1",
        "SITE,0,1,0,20,2,1,0,1,2048,1,0,0.005,4096,1", "FREQ,1,2,4", "POINT,0,2,4.5",
        "M,1,4.32e5,0.25,0.15,0.05,0.05", "R,1,1,0.8,0.8,0.2,0.1,0.1", "SC,1,100,100,100,0,0,0,0.02",
        "AOPT,0,0,0,1,1,1,0,0,0,0,0,0,0,0"]


def _ui(lines: List[str]):
    from ...prep import Interpreter
    ui = Interpreter()
    ui.run_text("\n".join(lines))
    return ui


def _excavated(beam_from: int) -> List[str]:
    """Excavated block (FV interaction nodes) with a structural beam from node ``beam_from`` to node
    30 above grade (K node 31)."""
    return SOIL + _block(2) + ["INT,1,27,1,1", "N,30,5,5,5", "N,31,10,5,5", "GROUP,2,BEAMS",
                               f"E,1,{beam_from},30,31"]


@problem("VP-52", "Model checks with planted faults (EXCSTRCHK, FIXEDINT, HINGED, KINT, FREESPRING, EDU-06)",
         tier="P0", modules=["UI"], source="requirements 6.4; spec 01 section 13 item 6")
def vp52(workdir):
    from ...prep.check import (ModelView, drilling_unrestrained, excstrchk, fixed_interaction, free_springs,
                               hinges, kint, run_check)
    r = VPResult()

    # ---- 1. EXCSTRCHK: beam attached to the interior excavation node 14
    ui = _ui(_excavated(14))
    v = ModelView(ui.model)
    exc = v.excavation()
    r.check("excavation interior nodes of a 2x2x2 block", len(exc["interior"]), 1, atol=0)
    res = excstrchk(v)
    r.check("EXCSTRCHK shared interior nodes (planted)", len(res), 1, atol=0)
    r.require("EXCSTRCHK reports node 14", list(res) == [14])
    rep = run_check(ui.model, modules=["HOUSE"])
    r.require("CHECK reports the EXCSTRCHK error", rep.has("Error", "EXCSTRCHK", "MODEL"))
    ui.execute("EXCSTRCHK")
    r.require("EXCSTRCHK command lists node 14", any("node 14:" in t for t in ui.sink.texts()))
    ui = _ui(_excavated(23))                      # beam on the top-face centre node: boundary, allowed
    r.check("EXCSTRCHK shared interior nodes (clean model)", len(excstrchk(ModelView(ui.model))), 0, atol=0)

    # ---- 2. FIXEDINT and Error 124
    ui = _ui(_excavated(23) + ["D,1,1,1,1,DISP", "D,2,2,1,1,UZ"])
    v = ModelView(ui.model)
    fi = dict(fixed_interaction(v))
    r.check("FIXEDINT fixed interaction nodes (planted 2)", len(fi), 2, atol=0)
    r.require("FIXEDINT: node 1 UX UY UZ, node 2 UZ", fi.get(1) == ["UX", "UY", "UZ"] and fi.get(2) == ["UZ"])
    rep = run_check(ui.model, modules=["HOUSE"])
    e124 = [m for m in rep.errors("MODEL") if m.number == 124]
    r.check("Error 124 count (all translations fixed, D-CHK-06)", len(e124), 1, atol=0)
    r.require("Error 124 text", e124 and e124[0].text == "Node 1 is a fixed interaction node")
    r.check("FIXEDINT on the clean model", len(fixed_interaction(ModelView(_ui(_excavated(23)).model))), 0, atol=0)

    # ---- 3. HINGED: a beam meeting a structural solid block at one node
    base = SOIL + _block(1, z0=0.0) + ["N,30,5,5,15", "N,31,10,5,15", "GROUP,2,BEAMS"]
    ui = _ui(base + ["E,1,23,30,31"])
    h = hinges(ModelView(ui.model))
    r.check("HINGED single-point beam-solid connections (planted)", len(h), 1, atol=0)
    r.require("HINGED reports node 23", h and h[0][0] == 23)
    ui = _ui(base + ["E,1,23,30,31", "E,2,23,22,31", "E,3,23,24,31"])    # massless beams embedded in the solid
    r.check("HINGED with beams penetrating the solid along its edges", len(hinges(ModelView(ui.model))), 0, atol=0)

    # ---- 4. KINT: K node 31 flagged as an interaction node
    ui = _ui(_excavated(23) + ["INT,31,31,1,1"])
    k = kint(ModelView(ui.model))
    r.check("KINT K nodes that are interaction nodes (planted)", len(k), 1, atol=0)
    r.require("KINT reports node 31 (K-only)", k == [(31, True)])
    r.check("KINT on the clean model", len(kint(ModelView(_ui(_excavated(23)).model))), 0, atol=0)

    # ---- 5. FREESPRING: spring from node 30 to an unconstrained node 32
    spring = _excavated(23) + ["N,32,5,5,6", "GROUP,3,SPRING", "RACT,1", "E,1,30,32"]
    fs = free_springs(ModelView(_ui(spring).model))
    r.check("FREESPRING free spring nodes (planted)", len(fs), 1, atol=0)
    r.require("FREESPRING reports node 32 without mass", fs == [(32, False)])
    r.check("FREESPRING with node 32 fixed", len(free_springs(ModelView(_ui(spring + ["D,32,32,1,1,ALL"]).model))),
            0, atol=0)

    # ---- 6. oblique Kirchhoff shell plate (30 deg about X) without FIXROT -> EDU-06
    c, s = math.cos(math.radians(30)), math.sin(math.radians(30))
    plate = []
    for j in range(3):
        for i in range(3):
            y = 5.0 * j
            plate.append(f"N,{1 + i + 3 * j},{5.0 * i},{y * c},{10 + y * s}")
    plate += ["GROUP,1,SHELL", "E,1,1,2,5,4", "E,2,2,3,6,5", "E,3,4,5,8,7", "E,4,5,6,9,8", "THICK,1,4,1,0.3",
              "M,1,4.32e5,0.25,0.15,0.05,0.05", "AOPT,0,0,0,0,0,1,0,0,0,0,0,0,0,0"]
    ui = _ui(plate)
    r.check("EDU-06 unrestrained drilling nodes (oblique plate)", len(drilling_unrestrained(ModelView(ui.model))),
            9, atol=0)
    rep = run_check(ui.model)
    r.check("CHECK Warning EDU-06 count", sum(1 for m in rep.warnings("MODEL") if m.number == "EDU-06"), 9, atol=0)
    ui.execute("FIXROT")
    g = ui.model.groups.get(2)
    r.check("FIXROT drilling springs added (oblique plate)", len(g.elements) if g else 0, 9, atol=0)
    sc = ui.model.springs[max(ui.model.springs)]
    # drilling spring k n_i^2 along the shell normal n = (0, -sin 30, cos 30): scyy = 10 sin^2(30) = 2.5
    r.check("FIXROT spring scyy = 10 sin^2(30) = 2.5 (normal (0, -sin 30, cos 30))", sc.scyy, 10 * s * s,
            rtol=1e-12)
    r.check("FIXROT spring sczz = 10 cos^2(30) = 7.5", sc.sczz, 10 * c * c, rtol=1e-12)
    r.check("EDU-06 after FIXROT", len(drilling_unrestrained(ModelView(ui.model))), 0, atol=0)
    flat = [ln for ln in plate if not ln.startswith("N,")] + [f"N,{1 + i + 3 * j},{5.0 * i},{5.0 * j},10"
                                                            for j in range(3) for i in range(3)]
    ui = _ui(flat)
    ui.execute("FIXROT")
    fixes = {tuple(nd.fix) for nd in ui.model.nodes.values()}
    r.require("flat XY plate: FIXROT fixes ROTZ only, no springs",
              fixes == {(0, 0, 0, 0, 0, 1)} and len(ui.model.groups) == 1)

    # ---- 7. shared free-field data (D-CHK-01): TOPL names an undefined layer 3
    ui = _ui(_excavated(23) + ["TOPL,0", "TOPL,1,3"])
    rep = run_check(ui.model)
    e19 = {m.module for m in rep.errors() if m.number == 19 and m.text == "Soil Layer 3 Is not Defined"}
    r.require("Error 19 for TOPL layer 3 under SITE and HOUSE (D-CHK-01)", e19 == {"SITE", "HOUSE"})
    r.require("SITE and HOUSE blocked, POINT not", rep.blocked("SITE") and rep.blocked("HOUSE")
              and not rep.blocked("POINT"))
    rep = run_check(_ui(_excavated(23)).model)
    r.check("Error 19 on the clean model", sum(1 for m in rep.errors() if m.number == 19), 0, atol=0)
    r.notes.append("planted-fault counts are exact; each check also runs on the model without the fault")
    return r
