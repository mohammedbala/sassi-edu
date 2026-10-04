"""VP-49: section-cut integration (requirements section 6.3, UT-21; spec 10 section 6 tests T-S1, T-S2;
requirements 4.14, D-SEC-01 ... D-SEC-05).

The two tests of spec 10 section 6 are run through the complete command workflow of manual 5.8.2:
``READSTR`` of an ``.ess`` frame on the original model, a cut (``CUTADD`` / ``SLICE``), ``CSECT`` into a
cross-section model, ``ACTM`` and ``CALCPAR``; the same resultants are then recomputed from the original
model by ``CALCSECTHIST`` over a list of the frames (7-column CSV with the ``MAX`` row).

T-S1 (solid section): two unit hexahedra, element 1 at x in [0, 1] and element 2 at x in [1, 2],
y, z in [0, 1]; plane through (0, 0, 0.5), n = (0, 0, 1), r = (1, 0, 0) (so ex = X, ey = Y, ez = Z).
Exact references: A = 2, C = (1, 0.5, 0.5), Ixx = 2 x 1^3/12 = 1/6, Iyy = 1 x 2^3/12 = 2/3, Ixy = 0;

* Szz = +10 (el. 1), -10 (el. 2): Fz = 0, My = +10, Mx = Mz = 0
  (each half carries a force of 10 at the lever arm 0.5: M = sum (c_i - C) x f_i);
* the same plus Sxz = 3 in both: Fx = 3 x 2 = 6;
* Syz = +1 (el. 1), -1 (el. 2): Mz = -1 (torsion of two opposite in-plane shears).

T-S2 (shell wall): two vertical shells in the plane y = 0, x in [0, 2] and [2, 4], z in [0, 2],
thickness 0.5, local x' = X and y' = Z; plane z = 1, n = (0, 0, 1), r = (1, 0, 0).  References:
A = 2 x 2 x 0.5 = 2, Iyy = 0.5 x 4^3/12 = 8/3;

* Sy'y' = +100 (el. 1), -100 (el. 2): Fz = 0, My = +200;
* Sy'y' = +100 in both: Fz = 200.

All references are exact (the stresses are uniform over each piece), so the tolerances are round-off.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Sequence, Tuple

from .. import VPResult, problem

#: round-off tolerances of the exact references
RTOL = 1e-9
ATOL = 1e-9


# ======================================================================================
# Model builders (also used by the unit tests)
# ======================================================================================
def two_hex_lines() -> List[str]:
    """T-S1 model: two unit hexahedra side by side along X (group 1, elements 1 and 2)."""
    pts: Dict[Tuple[int, int, int], int] = {}
    out = []
    nid = 0
    for z in (0, 1):
        for y in (0, 1):
            for x in (0, 1, 2):
                nid += 1
                pts[(x, y, z)] = nid
                out.append(f"N,{nid},{x},{y},{z}")
    out.append("GROUP,1,SOLID")
    for e, x0 in ((1, 0), (2, 1)):
        ids = [pts[(x0, 0, 0)], pts[(x0 + 1, 0, 0)], pts[(x0 + 1, 1, 0)], pts[(x0, 1, 0)],
               pts[(x0, 0, 1)], pts[(x0 + 1, 0, 1)], pts[(x0 + 1, 1, 1)], pts[(x0, 1, 1)]]
        out.append(f"E,{e}," + ",".join(str(k) for k in ids))
    return out


def shell_wall_lines(thick: float = 0.5) -> List[str]:
    """T-S2 model: two 2 x 2 shells in the plane y = 0 (local x' = X, y' = Z), thickness 0.5."""
    return ["N,1,0,0,0", "N,2,2,0,0", "N,3,4,0,0", "N,4,0,0,2", "N,5,2,0,2", "N,6,4,0,2",
            "GROUP,1,SHELL", "E,1,1,2,5,4", "E,2,2,3,6,5", f"THICK,1,2,1,{thick}"]


def _run(ui, lines: Sequence[str], res: VPResult, what: str) -> bool:
    from sassi.prep import Kind
    n0 = len(ui.sink.texts(Kind.ERROR))
    for ln in lines:
        ui.execute(ln)
    errs = ui.sink.texts(Kind.ERROR)[n0:]
    if errs:
        res.notes.append(f"{what}: errors {errs}")
    return res.require(f"{what}: commands run without error", not errs)


def _close(res: VPResult, label: str, computed: float, reference: float) -> None:
    if reference == 0.0:
        res.check(label, computed, reference, atol=ATOL)
    else:
        res.check(label, computed, reference, rtol=RTOL)


def _calcpar(ui, res: VPResult, ess: Path, cut_cmd: str, plane: str, calcpar: str, what: str) -> Dict[str, float]:
    """READSTR -> cut -> CSECT -> ACTM -> CALCPAR (manual 5.8.2) and the CALCPAR values."""
    ok = _run(ui, ["ACTM,0", "CUTCLR,1", f"READSTR,{ess.name}", cut_cmd, f"CSECT,9,1,{plane}", "ACTM,9", calcpar,
                   "ACTM,0"], res, what)
    return dict(ui.session.get("calcpar", {})) if ok else {}


@problem("VP-49", "Section-cut integration T-S1 (solid) and T-S2 (shell wall)", tier="P1",
         modules=["READSTR", "CSECT", "CALCPAR", "CALCSECTHIST"], source="spec 10 section 6 (T-S1, T-S2)")
def vp49(workdir: Path) -> VPResult:
    from sassi.prep import Interpreter
    from sassi.prep.cuts_lib import write_ess
    res = VPResult()
    wd = Path(workdir)
    # ---------------------------------------------------------------- T-S1
    ui = Interpreter(cwd=wd)
    ui.run_text("\n".join(two_hex_lines()), name="ts1.pre")
    states = {
        "Szz +10/-10": {(1, 1): [0, 0, 10, 0, 0, 0], (1, 2): [0, 0, -10, 0, 0, 0]},
        "Szz +10/-10, Sxz 3": {(1, 1): [0, 0, 10, 0, 3, 0], (1, 2): [0, 0, -10, 0, 3, 0]},
        "Syz +1/-1": {(1, 1): [0, 0, 0, 0, 0, 1], (1, 2): [0, 0, 0, 0, 0, -1]},
    }
    expect = {
        "Szz +10/-10": {"Fz": 0.0, "My": 10.0, "Mx": 0.0, "Mz": 0.0},
        "Szz +10/-10, Sxz 3": {"Fx": 6.0, "Fz": 0.0, "My": 10.0},
        "Syz +1/-1": {"Mz": -1.0, "Fz": 0.0},
    }
    files = []
    for k, (name, vals) in enumerate(states.items(), start=1):
        f = write_ess(wd / f"TS1_{k}.ess", {1: "SOLID"}, vals)
        files.append(f)
        par = _calcpar(ui, res, f, "CUTADD,1,1,1-2", "0,0,0.5,0,0,1", "CALCPAR,0,0,1,1,0,0,11", f"T-S1 {name}")
        if not par:
            continue
        if k == 1:
            for q, ref in (("Area", 2.0), ("Xc", 1.0), ("Yc", 0.5), ("Zc", 0.5), ("Ixx", 1.0 / 6.0),
                           ("Iyy", 2.0 / 3.0), ("Ixy", 0.0)):
                _close(res, f"T-S1 {q}", par[q], ref)
        for q, ref in expect[name].items():
            _close(res, f"T-S1 [{name}] {q}", par[q], ref)
    # the same three states as a history from the original model (CALCSECTHIST, ts = 0.01)
    (wd / "TS1.lst").write_text("T-S1 frames\n" + "\n".join(f.name for f in files) + "\n", encoding="utf-8")
    if _run(ui, ["ACTM,0", "CUTCLR,1", "SLICE,1,0,0,0.5,0,0,1",
                 "CALCSECTHIST,TS1.lst,1,0,0,0.5,0,0,1,1,0,0,12,0.01,TS1.csv"], res, "T-S1 CALCSECTHIST"):
        rows = [ln.split(",") for ln in (wd / "TS1.csv").read_text().splitlines()]
        res.require("T-S1 CALCSECTHIST header Time,Fx,Fy,Fz,Mx,My,Mz",
                    rows[0] == ["Time", "Fx", "Fy", "Fz", "Mx", "My", "Mz"])
        res.require("T-S1 CALCSECTHIST 3 frames + MAX row", len(rows) == 5 and rows[-1][0] == "MAX")
        col = {n: i for i, n in enumerate(rows[0])}
        for k, name in enumerate(states, start=1):
            for q, ref in expect[name].items():
                _close(res, f"T-S1 CALCSECTHIST frame {k} {q}", float(rows[k][col[q]]), ref)
        _close(res, "T-S1 CALCSECTHIST time of frame 3 = 2 ts", float(rows[3][0]), 0.02)
        _close(res, "T-S1 CALCSECTHIST MAX Mz (signed)", float(rows[-1][col["Mz"]]), -1.0)
    # ---------------------------------------------------------------- T-S2
    ui2 = Interpreter(cwd=wd)
    ui2.run_text("\n".join(shell_wall_lines()), name="ts2.pre")
    s2 = {"Sy'y' +100/-100": ({(1, 1): [0, 100, 0, 0, 0, 0], (1, 2): [0, -100, 0, 0, 0, 0]}, {"Fz": 0.0, "My": 200.0}),
          "Sy'y' +100 both": ({(1, 1): [0, 100, 0, 0, 0, 0], (1, 2): [0, 100, 0, 0, 0, 0]}, {"Fz": 200.0, "My": 0.0})}
    for k, (name, (vals, exp)) in enumerate(s2.items(), start=1):
        f = write_ess(wd / f"TS2_{k}.ess", {1: "SHELL"}, vals)
        par = _calcpar(ui2, res, f, "SLICE,1,0,0,1,0,0,1", "0,0,1,0,0,1", "CALCPAR,0,0,1,1,0,0,11", f"T-S2 {name}")
        if not par:
            continue
        if k == 1:
            _close(res, "T-S2 Area", par["Area"], 2.0)
            _close(res, "T-S2 Iyy", par["Iyy"], 8.0 / 3.0)
        for q, ref in exp.items():
            _close(res, f"T-S2 [{name}] {q}", par[q], ref)
    res.notes.append("T-S1/T-S2 run through READSTR -> CUTADD/SLICE -> CSECT -> CALCPAR (manual 5.8.2); T-S1 also "
                     "through CALCSECTHIST on the original model. References are exact: tolerances are round-off.")
    return res
