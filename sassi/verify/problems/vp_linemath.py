"""VP-47, VP-48, VP-53: line mathematics, spectrum broadening, line / frame files and CRITFREQ.

The problems run the *commands* through the interpreter (READSPEC ... BROADEN, CRITFREQ,
FRAMECOMBIN) so that the parsing, the line store and the maths are verified together; the
references are the closed-form numbers of requirements section 6.4 (spec 10 section 4.1 tests T-L1,
T-B1..T-B3; spec 06 section 12; spec 09 section 8 item 10).

* VP-47 (T-L1): A x = [0, 1, 2], y = [0, 1, 2]; B x = [0.5, 1.5], y = [10, 20].  Union grid
  [0, 0.5, 1, 1.5, 2] with B held constant outside its range.  Plus the spec 06 case
  L1 = {(1,1), (3,3)}, L2 = {(2,10), (4,20)}: ADDITION = {11, 12, 18, 23}.
* VP-48 (T-B1, T-B2, T-B3): BROADEN of a spike (plateau 8.5..11.5 at 5, ramps (7.65, 1)-(8.5, 5)
  and (11.5, 5)-(12.65, 1)), bridging (flat 3.8 from 5.2 to 8; no change at Smooth1 = 10 with
  the valley-depth criterion), envelope of two lines.  In addition the result of a random
  spectrum is compared with a brute-force evaluation of the D-LIN-01 definition
  ``B(f) = max{E(f') : f/(1+b) <= f' <= f/(1-b)}`` on a dense grid.
* VP-53: WRITESPEC -> READSPEC and WRITETH -> READTH round trips (exact), READSPEC column
  counting (1 + 3 columns, numLines = 3 -> 3 lines), FRAMECOMBIN SRSS / sum / average, CRITFREQ
  flagging a 40 % TFI overshoot at tol 20 % but not at 50 %.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from .. import VPResult, problem


def _interp(workdir: Path):
    from ...prep import Interpreter
    from ...plotting.state import plot_state
    ui = Interpreter(cwd=workdir)
    plot_state(ui).auto_render = False           # no images needed for the numbers
    return ui, plot_state(ui)


def _errors(ui):
    from ...prep.messages import Kind
    return ui.sink.texts(Kind.ERROR)


def _write_xy(path: Path, x, y) -> None:
    path.write_text("\n".join(f"{float(a)!r} {float(b)!r}" for a, b in zip(x, y)) + "\n", encoding="utf-8")


def _value(line, x: float) -> float:
    return float(np.interp(x, line.x, line.y))


@problem("VP-47", "Line mathematics on the union of abscissas (T-L1)", tier="P1", modules=["UI"],
         source="requirements 6.4; spec 10 4.1.7 T-L1; spec 06 12.1")
def vp47(workdir: Path) -> VPResult:
    r = VPResult()
    ui, st = _interp(workdir)
    _write_xy(workdir / "A.txt", [0, 1, 2], [0, 1, 2])
    _write_xy(workdir / "B.txt", [0.5, 1.5], [10, 20])
    ui.run_text("READSPEC,A.txt,1,1\nREADSPEC,B.txt,1,2\nADDITION,3,1,2\nAVERAGE,4,1,2\nSUBTRACTION,5,1,2\n"
                "LINECOMBIN,6,1,2,2,0.5\nSRSS,7,1,2\n")
    r.require("no command errors", not _errors(ui), "; ".join(_errors(ui)))
    grid = [0.0, 0.5, 1.0, 1.5, 2.0]
    refs = {3: ("ADDITION", [10, 10.5, 16, 21.5, 22], 1e-12),
            4: ("AVERAGE", [5, 5.25, 8, 10.75, 11], 1e-12),
            5: ("SUBTRACTION A-B", [-10, -9.5, -14, -18.5, -18], 1e-12),
            6: ("LINECOMBIN 2A+0.5B", [5, 6, 9.5, 13, 14], 1e-12),
            7: ("SRSS", [10, 10.01249, 15.03330, 20.05617, 20.09975], None)}
    for num, (name, ref, tol) in refs.items():
        L = st.lines.get(num)
        if L is None:
            r.require(f"{name} line created", False)
            continue
        r.check(f"{name}: number of grid points", L.n, 5, atol=0)
        for xg, (xv, yv, yr) in zip(grid, zip(L.x, L.y, ref)):
            r.check(f"{name}: abscissa {xg:g}", xv, xg, atol=1e-15)
            if tol is None:
                r.check(f"{name}: y({xg:g})", yv, yr, atol=5e-6, note="reference printed to 5 decimals")
            else:
                r.check(f"{name}: y({xg:g})", yv, yr, atol=tol)
    r.require("default names ('Linear Combin.', 'Average Line', 'SRSS Line')",
              st.lines[3].name == "Linear Combin." and st.lines[4].name == "Average Line"
              and st.lines[7].name == "SRSS Line")
    # spec 06 section 12 item 1
    _write_xy(workdir / "L1.txt", [1, 3], [1, 3])
    _write_xy(workdir / "L2.txt", [2, 4], [10, 20])
    ui.run_text("READSPEC,L1.txt,1,11\nREADSPEC,L2.txt,1,12\nADDITION,13,11,12\n")
    L = st.lines[13]
    for xg, yr in zip([1, 2, 3, 4], [11, 12, 18, 23]):
        r.check(f"spec 06 case: ADDITION y({xg})", _value(L, xg), yr, atol=1e-12)
    r.check("spec 06 case: grid size", L.n, 4, atol=0)
    # spec 06 section 12 items 2 and 3
    _write_xy(workdir / "one.txt", [0, 1, 2], [1, 1, 1])
    _write_xy(workdir / "three.txt", [0, 1, 2], [3, 3, 3])
    ui.run_text("READSPEC,one.txt,1,21\nREADSPEC,one.txt,1,22\nREADSPEC,one.txt,1,23\nSRSS,24,21,22,23\n"
                "READSPEC,three.txt,1,25\nAVERAGE,26,21,25\nLINECOMBIN,27,25,1.0,21,0.4,21,0.4\n"
                "SUBTRACTION,28,25,21,21\n")
    r.check("SRSS of three lines y = 1", float(st.lines[24].y.max()), math.sqrt(3.0), rtol=1e-14)
    r.check("AVERAGE of y = 1 and y = 3", float(st.lines[26].y.mean()), 2.0, rtol=1e-14)
    r.check("LINECOMBIN 1.0/0.4/0.4 (100-40-40 rule)", float(st.lines[27].y[0]), 3.8, rtol=1e-14)
    r.check("SUBTRACTION y1 - y2 - y3", float(st.lines[28].y[0]), 1.0, rtol=1e-14)
    r.require("no command errors (spec 06 cases)", not _errors(ui), "; ".join(_errors(ui)))
    return r


def _brute_broaden(x, y, b, f):
    """Independent dense evaluation of max{E(f') : f/(1+b) <= f' <= f/(1-b)} (E = the line)."""
    out = np.empty(len(f))
    for i, fi in enumerate(f):
        lo, hi = fi / (1 + b), fi / (1 - b)
        fp = np.concatenate([np.linspace(lo, hi, 4001), x[(x > lo) & (x < hi)]])
        out[i] = np.interp(fp, x, y, left=y[0], right=y[-1]).max()
    return out


@problem("VP-48", "BROADEN: envelope, peak broadening and peak bridging (T-B1..T-B3)", tier="P1", modules=["UI"],
         source="requirements 4.13, 6.4; D-LIN-01; spec 10 4.1.3")
def vp48(workdir: Path) -> VPResult:
    r = VPResult()
    ui, st = _interp(workdir)
    # ---- T-B1: broadening of a spike by +-15 %
    _write_xy(workdir / "tb1.rs", [1, 5, 9, 10, 11, 15, 20], [1, 1, 1, 5, 1, 1, 1])
    ui.run_text("READSPEC,tb1.rs,1,1\nBROADEN,2,0,15,1\n")
    B = st.lines[2]
    for xg in (8.5, 9.0, 10.0, 11.0, 11.5):
        r.check(f"T-B1: B({xg:g}) on the plateau", _value(B, xg), 5.0, atol=1e-12)
    for xg, yr in ((7.65, 1.0), (8.075, 3.0), (11.5, 5.0), (12.075, 3.0), (12.65, 1.0), (1.0, 1.0), (5.0, 1.0),
                   (7.0, 1.0), (14.0, 1.0), (20.0, 1.0)):
        r.check(f"T-B1: B({xg:g})", _value(B, xg), yr, atol=1e-12)
    for xg in (7.65, 8.5, 9.35, 10.35, 11.5, 12.65):
        r.require(f"T-B1: augmented grid contains {xg:g}", bool(np.any(np.abs(B.x - xg) < 1e-12)))
    r.check("T-B1: B = 1 outside [7.65, 12.65] (max deviation)",
            float(np.max(np.abs(B.y[(B.x < 7.65 - 1e-9) | (B.x > 12.65 + 1e-9)] - 1.0))), 0.0, atol=1e-12)
    # ---- T-B2: bridging with Smooth1 = 15 % (valley-depth criterion)
    xb = [1, 4, 5, 6, 6.5, 7, 8, 9, 12, 20]
    yb = [0.5, 2, 4, 3, 3.5, 3.0, 3.8, 2, 1, 0.5]
    _write_xy(workdir / "tb2.rs", xb, yb)
    ui.run_text("READSPEC,tb2.rs,1,3\nBROADEN,4,15,0,3\nBROADEN,5,10,0,3\n")
    B2 = st.lines[4]
    for xg, yr in ((1, 0.5), (4, 2.0), (5, 4.0), (5.1, 3.9), (9, 2.0), (12, 1.0), (20, 0.5)):
        r.check(f"T-B2: unchanged y({xg:g})", _value(B2, xg), yr, atol=1e-12)
    for xg in (5.2, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0):
        r.check(f"T-B2: flat 3.8 at {xg:g}", _value(B2, xg), 3.8, atol=1e-12)
    r.require("T-B2: crossing point 5.2 in the grid", bool(np.any(np.abs(B2.x - 5.2) < 1e-12)))
    B3 = st.lines[5]
    r.check("T-B2 with Smooth1 = 10 %: unchanged (valley-depth criterion), max |B - y|",
            float(np.max(np.abs(np.interp(xb, B3.x, B3.y) - np.asarray(yb)))), 0.0, atol=1e-12)
    # ---- T-B3: envelope of two lines (b = 0, p = 0) = pointwise maximum
    _write_xy(workdir / "e1.rs", [1, 2, 3, 4, 5], [1, 4, 1, 1, 1])
    _write_xy(workdir / "e2.rs", [1.5, 2.5, 3.5, 4.5], [1, 1, 3, 1])
    ui.run_text("READSPEC,e1.rs,1,6\nREADSPEC,e2.rs,1,7\nBROADEN,8,0,0,6,7\n")
    E = st.lines[8]
    grid = np.union1d(st.lines[6].x, st.lines[7].x)
    ref = np.maximum(st.lines[6].at(grid), st.lines[7].at(grid))
    r.check("T-B3: max |E - max(y1, y2)| on the union grid", float(np.max(np.abs(E.at(grid) - ref))), 0.0,
            atol=1e-14)
    dense = np.linspace(1, 5, 2001)
    r.check("T-B3: envelope >= each line everywhere (min of E - max y_i, dense)",
            float(np.min(E.at(dense) - np.maximum(st.lines[6].at(dense), st.lines[7].at(dense)))), 0.0,
            atol=1e-12, note="crossing points of the lines are inserted (exact upper envelope)")
    ui.run_text("BROADEN,9,0,0,6\n")
    r.check("single line, Smooth1 = Smooth2 = 0: unchanged", float(np.max(np.abs(st.lines[9].y
                                                                             - st.lines[6].y))), 0.0, atol=0)
    # ---- exactness against the definition (random spectrum, 15 % broadening)
    rng = np.random.default_rng(47)
    xs = np.sort(rng.uniform(0.5, 30.0, 60))
    ys = np.abs(rng.normal(1.0, 0.6, 60)) + 0.2
    _write_xy(workdir / "rand.rs", xs, ys)
    ui.run_text("READSPEC,rand.rs,1,10\nBROADEN,11,0,15,10\n")
    Bx = st.lines[11]
    f = np.linspace(xs[0], xs[-1], 1500)
    brute = _brute_broaden(xs, ys, 0.15, f)
    err = float(np.max(Bx.at(f) - brute))
    r.check("random spectrum: max(B_code - B_definition) (code never above the definition)", max(err, 0.0), 0.0,
            atol=1e-9)
    r.check("random spectrum: max |B_code - B_definition| (brute force over the window ends and every vertex)",
            float(np.max(np.abs(Bx.at(f) - brute))), 0.0, atol=1e-9 * float(ys.max()),
            note="the maximum of a piecewise-linear line over a window lies at a vertex or a window end, both "
                 "evaluated by the brute force, so it is exact (final audit: the tolerance was 2e-3 max(y))")
    r.require("no command errors", not _errors(ui), "; ".join(_errors(ui)))
    r.notes.append(f"T-B1 result: {B.n} points; random-spectrum max deviation from the dense brute force "
                   f"{float(np.max(np.abs(Bx.at(f) - brute))):.2e}")
    return r


@problem("VP-53", "Line / frame files, READSPEC column counting, FRAMECOMBIN and CRITFREQ", tier="P1",
         modules=["UI", "MOTION"], source="requirements 6.4; spec 10 4.2; spec 09 2.4, 2.11, 8.10; D-FIL-02/03")
def vp53(workdir: Path) -> VPResult:
    from ...plotting.state import read_frame, write_frame
    r = VPResult()
    ui, st = _interp(workdir)
    rng = np.random.default_rng(53)
    # ---- WRITESPEC -> READSPEC round trip on the union grid
    x1 = np.sort(rng.uniform(0.1, 50, 37))
    x2 = np.sort(rng.uniform(0.2, 40, 23))
    _write_xy(workdir / "s1.rs", x1, rng.normal(size=37) * np.pi)
    _write_xy(workdir / "s2.rs", x2, rng.normal(size=23) / 3.0)
    ui.run_text("READSPEC,s1.rs,1,1\nREADSPEC,s2.rs,1,2\nWRITESPEC,both.rs,1,2\nREADSPEC,both.rs,2,11,12\n")
    X = np.union1d(x1, x2)
    for a, b in ((1, 11), (2, 12)):
        La, Lb = st.lines[a], st.lines[b]
        r.check(f"WRITESPEC/READSPEC line {a}: grid = union of abscissas (max |dx|)",
                float(np.max(np.abs(Lb.x - X))) if Lb.n == len(X) else 1.0, 0.0, atol=0)
        r.check(f"WRITESPEC/READSPEC line {a}: values (max |dy|, exact)", float(np.max(np.abs(Lb.y - La.at(X)))),
                0.0, atol=0)
    # ---- READSPEC column counting: 1 + 3 columns, numLines = 3
    f = np.linspace(0.5, 10, 20)
    cols = [np.sin(f), np.cos(f), f ** 2]
    (workdir / "cols.tfu").write_text("# f c1 c2 c3\n" + "\n".join(
        " ".join(repr(float(v)) for v in (f[i], cols[0][i], cols[1][i], cols[2][i])) for i in range(len(f))) + "\n")
    ui.run_text("READSPEC,cols.tfu,3,21,22,23\n")
    made = [n for n in (21, 22, 23) if n in st.lines]
    r.check("READSPEC 1 + 3 columns, numLines 3: lines created", len(made), 3, atol=0)
    if len(made) == 3:
        for k, n in enumerate((21, 22, 23)):
            r.check(f"column {k + 1} -> line {n} (frequency column not counted)",
                    float(np.max(np.abs(st.lines[n].y - cols[k]))), 0.0, atol=0)
        r.require("line names cols.tfu, cols.tfu[2], cols.tfu[3] (D-UI-16)",
                  [st.lines[n].name for n in (21, 22, 23)] == ["cols.tfu", "cols.tfu[2]", "cols.tfu[3]"])
    ok4 = ui.execute("READSPEC,cols.tfu,4,31")
    r.require("READSPEC numLines 4 on a 1 + 3 column file is an error", not ok4)
    # ---- WRITETH -> READTH round trip
    dt = 0.005
    acc = rng.normal(size=400) * 0.1
    (workdir / "in.acc").write_text(repr(dt) + "\n" + "\n".join(repr(float(v)) for v in acc) + "\n")
    ui.run_text("READTH,in.acc,0,41\nWRITETH,out.acc,41\nREADTH,out.acc,0,42\n")
    L1, L2 = st.lines[41], st.lines[42]
    r.check("READTH Pair 0: t_k = (k-1) dt (last time)", float(L1.x[-1]), 399 * dt, rtol=1e-14)
    r.check("WRITETH/READTH: times (max |dt|)", float(np.max(np.abs(L1.x - L2.x))), 0.0, atol=0)
    r.check("WRITETH/READTH: values (max |da|, exact)", float(np.max(np.abs(L1.y - L2.y))), 0.0, atol=0)
    (workdir / "pairs.th").write_text("\n".join(f"{k * dt!r} {float(v)!r}" for k, v in enumerate(acc[:50])) + "\n")
    ui.run_text("READTH,pairs.th,1,43\n")
    r.check("READTH Pair 1 (time/value pairs) = Pair 0 data", float(np.max(np.abs(st.lines[43].y - acc[:50]))),
            0.0, atol=0)
    # ---- FRAMECOMBIN SRSS / sum / average
    nodes = np.array([3, 7, 11, 15])
    fa = np.array([[3.0, -1.0, 0.5], [0.0, 2.0, 1.0], [4.0, 0.0, -2.0], [1.0, 1.0, 1.0]])
    fb = np.array([[4.0, 1.0, 0.5], [1.0, -2.0, 3.0], [3.0, 0.0, 2.0], [1.0, 1.0, 1.0]])
    write_frame(workdir / "ACC_00.000_00001", nodes, fa)
    write_frame(workdir / "ACC_00.005_00002", nodes, fb)
    for op, name, ref in ((0, "SRSS", np.sqrt(fa ** 2 + fb ** 2)), (1, "sum", fa + fb), (2, "average", (fa + fb) / 2)):
        ok = ui.execute(f"FRAMECOMBIN,{op},2,ACC_00.000_00001,ACC_00.005_00002,comb{op}.txt")
        r.require(f"FRAMECOMBIN op {op} ({name}) succeeded", ok)
        if ok:
            fr = read_frame(workdir / f"comb{op}.txt")
            r.require(f"FRAMECOMBIN {name}: node ids copied", bool(np.array_equal(fr.nodes, nodes)))
            r.check(f"FRAMECOMBIN {name}: max |value - reference|", float(np.max(np.abs(fr.values - ref))), 0.0,
                    atol=1e-9, note="frames are written with 11 significant digits")
    r.check("FRAMECOMBIN SRSS (3, 4) -> 5", float(read_frame(workdir / "comb0.txt").values[0, 0]), 5.0, rtol=1e-12)
    # ---- CRITFREQ: 40 % overshoot of the interpolated TF between two SSI frequencies
    F = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    AU = np.array([1.0, 1.2, 2.0, 2.5, 1.5, 1.0])
    fi = np.round(np.arange(1.0, 6.0 + 1e-9, 0.05), 10)
    ai = np.interp(fi, F, AU) + 1.25 * np.clip(1.0 - np.abs(fi - 3.5) / 0.5, 0.0, None)
    (workdir / "00010TR_X.TFU").write_text("# computed\n" + "\n".join(f"{float(a)!r} {float(b)!r} 0.0" for a, b in zip(F, AU)))
    (workdir / "00010TR_X.TFI").write_text("# interpolated\n" + "\n".join(f"{float(a)!r} {float(b)!r} 0.0" for a, b in zip(fi, ai)))
    r.check("CRITFREQ test data: TFI peak / max computed neighbour", float(ai.max() / 2.5), 1.4, rtol=1e-12)
    ui.run_text("CRITFREQ,20,50,00010TR_X,FR\nCRITFREQ,50,50,00010TR_X,FR50\n")
    fr20 = ui.variables.get("fr")
    fr50 = ui.variables.get("fr50")
    r.require("CRITFREQ tol 20 %: the 3.5 Hz peak is flagged (frequency number 70 at df 0.05 Hz)",
              fr20 is not None and fr20.items == ["70"], f"items {fr20.items if fr20 else None}")
    r.require("CRITFREQ tol 50 %: nothing flagged", fr50 is not None and fr50.items == [],
              f"items {fr50.items if fr50 else None}")
    r.require("no command errors", not [e for e in _errors(ui) if "READSPEC" not in e], "; ".join(_errors(ui)))
    return r
