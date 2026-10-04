"""Programming commands: macros (UT-09), variables and loops (UT-10), REDUCESET and random lists (UT-11),
substitution grammar (D-PAR-15/16) and the other commands of section 3.4.M."""
from __future__ import annotations

import numpy as np
import pytest

from sassi.prep import Interpreter, Kind


class Recorder(Interpreter):
    """Interpreter that records every line reaching execute() (after macro substitution)."""

    def __init__(self, *a, **k):
        super().__init__(*a, **k)
        self.lines = []

    def execute(self, line):
        self.lines.append(line)
        return super().execute(line)


@pytest.fixture
def ui(tmp_path):
    return Interpreter(cwd=tmp_path)


def run(ui, text):
    ui.run_text(text)
    return ui.model


def errors(ui):
    return ui.sink.texts(Kind.ERROR)


def warnings(ui):
    return ui.sink.texts(Kind.WARNING)


def var(ui, name):
    return ui.variables[name.lower()]


# ------------------------------------------------------------------ UT-09 macros (manual 5.6.1-5.6.3)
def test_ut09_basic_macro(ui, tmp_path):
    (tmp_path / "Node-Macro.pre").write_text("N,1,$1$,$2$,$3$\n")
    m = run(ui, "LoadMacro,Node,.\\Node-Macro.pre\nMacro,Node,13.52,15,100.25\n")
    assert m.nodes[1].xyz == (13.52, 15.0, 100.25)
    assert "NODE" in ui.macros                      # names are upper-cased


SRSS_MACRO = """READSPEC,$1$,1,1
READSPEC,$2$,1,2
READSPEC,$3$,1,3
SRSS,4,1,2,3
WRITESPEC,$4$,4
SPECPLOT,1,2,3,4
CAPTUREPLOT,$4$.png
CLOSEPLOT
"""


def test_ut09_graphing_and_nested_drivers_identical(tmp_path):
    (tmp_path / "SRSS-macro.pre").write_text(SRSS_MACRO)
    (tmp_path / "Nested-macro.pre").write_text("MACRO,SRSS,Node$1$x.rs,Node$1$y.rs,Node$1$z.rs,Node$1$srss.rs\n")
    a = Recorder(cwd=tmp_path)
    a.run_text("LOADMACRO,SRSS,SRSS-macro.pre\n" + "".join(
        f"MACRO,SRSS,Node{n}x.rs,Node{n}y.rs,Node{n}z.rs,Node{n}srss.rs\n" for n in (1, 2, 17)))
    b = Recorder(cwd=tmp_path)
    b.run_text("LOADMACRO,SRSS,SRSS-macro.pre\nLOADMACRO,NESTED,Nested-macro.pre\n"
               "MACRO,NESTED,1\nMACRO,NESTED,2\nMACRO,NESTED,17\n")
    body = lambda r: [ln for ln in r.lines if not ln.upper().startswith(("MACRO", "LOADMACRO"))]
    assert body(a) == body(b)
    assert "WRITESPEC,Node17srss.rs,4" in body(a) and "CAPTUREPLOT,Node17srss.rs.png" in body(a)
    assert len(body(a)) == 3 * 8


def test_macro_missing_argument_namespace_and_recursion(ui, tmp_path):
    (tmp_path / "m.pre").write_text("N,$1$,$2$,0,$3$\n")
    (tmp_path / "self.pre").write_text("MACRO,SELF\n")
    m = run(ui, "LOADMACRO,N,m.pre\nMACRO,N,5,1\nN,6,1,1,1\nLOADMACRO,SELF,self.pre\nMACRO,SELF\nMACRO,NOPE\n")
    assert m.nodes[5].xyz == (1.0, 0.0, 0.0)        # missing $3$ -> empty -> default 0
    assert any("argument 3 missing" in w for w in warnings(ui))
    assert m.nodes[6].xyz == (1.0, 1.0, 1.0)        # a macro named N does not override the command
    assert any("nesting deeper than 64" in e for e in errors(ui))
    assert any("NOPE is not loaded" in e for e in errors(ui))
    run(ui, "LOADMACRO," + "A" * 51 + ",m.pre\nMACROLIST")
    assert any("at most 50 characters" in e for e in errors(ui))
    assert any(t.startswith("N ") and t.endswith("m.pre") for t in ui.sink.texts(Kind.INFO))


def test_macro_lines_substitute_variables_after_placeholders(ui, tmp_path):
    (tmp_path / "v.pre").write_text("N,@C++,$1$,@C,0\n")
    m = run(ui, "VAR,C\nLOADMACRO,V,v.pre\nMACRO,V,7\nMACRO,V,8\n")
    assert m.nodes[1].xyz == (7.0, 1.0, 0.0) and m.nodes[2].xyz == (8.0, 2.0, 0.0)


# ------------------------------------------------------------------ UT-10 variables and loops (manual 5.9)
LOOP_PRE = """Var,NNUM
Var,X,1,2,3,4,5
Var,Y,1,2,3,4,5
Var,Z,1,2,3,4,5
ForEach,Z,ForEach,Y,ForEach,X,N,@NNUM++,@X[#],@Y[#],@Z[#]
"""


def test_ut10_loop_pre(ui, tmp_path):
    (tmp_path / "Loop.pre").write_text(LOOP_PRE)
    ui.execute("INP,Loop.pre")
    m = ui.model
    assert len(m.nodes) == 125 and not errors(ui)
    for iz in range(1, 6):
        for iy in range(1, 6):
            for ix in range(1, 6):
                n = 25 * (iz - 1) + 5 * (iy - 1) + ix
                assert m.nodes[n].xyz == (float(ix), float(iy), float(iz))
    assert var(ui, "NNUM").counter == 125
    # the EOF summary counts the file's command lines, not the loop iterations
    assert any(t.startswith("INP Loop.pre: 5 lines, 5 commands, 0 warnings, 0 errors")
               for t in ui.sink.texts(Kind.INFO))


def test_ut10_counter_operators(ui):
    run(ui, "VAR,X,1.23,2.83,3\nSETVAR,@X+5\n")
    assert var(ui, "x").counter == 5
    run(ui, "SETVAR,@X++")
    assert var(ui, "X").counter == 6
    run(ui, "SETVAR,@X=-15")
    assert var(ui, "X").counter == -15
    run(ui, "SETVAR,@X--,@X-4,@X")
    assert ui.sink.texts(Kind.INFO)[-3] == "SETVAR,-16,-20,-20"
    run(ui, "VAR,X,9")
    assert var(ui, "X").counter == 0 and var(ui, "X").items == ["9"]
    run(ui, "var,x,a,b")                             # names are case-insensitive
    assert len(ui.variables) == 1 and var(ui, "X").items == ["a", "b"]


def test_ut10_nested_foreach_same_variable_is_an_error(ui):
    m = run(ui, "VAR,X,1,2\nVAR,C\nFOREACH,X,FOREACH,X,N,@C++,0,0,0\n")
    assert not m.nodes
    assert sum("cannot reuse the variable X" in e for e in errors(ui)) == 2   # once per outer iteration


def test_items_index_and_hash_grammar(ui):
    m = run(ui, "VAR,X,10,20,30\nVAR,C\nN,1@X[2],@X[3],0,0\nFOREACH,X,N,@X[#],#,#+1,@X[#-0]\n")
    assert m.nodes[120].xyz == (30.0, 0.0, 0.0)
    assert m.nodes[10].xyz == (1.0, 2.0, 10.0) and m.nodes[30].xyz == (3.0, 4.0, 30.0)
    run(ui, "N,@X[4],0,0,0\n")                       # D-PAR-16: out of range -> error, command skipped
    assert any("index out of range 1..3" in e for e in errors(ui))
    run(ui, "TIT,Run #1 of @NOTDEFINED\n")           # '#' outside FOREACH is data; unknown name kept
    assert m.title == "Run #1 of @NOTDEFINED"
    assert any("'NOTDEFINED' is not defined" in w for w in warnings(ui))
    run(ui, "VAR,I,2\nSETVAR,@I=3\nN,99,@X[@I],0,0\n")
    assert m.nodes[99].xyz == (30.0, 0.0, 0.0)


def test_variable_inside_longer_token(ui):
    run(ui, "VAR,X,17,18\nVAR,F\nFOREACH,X,VAR,F,Node@X[#].rs\n")
    assert var(ui, "F").items == ["Node18.rs"]
    run(ui, "VAR,XY,a\nVAR,G,@XY@X\n")               # longest defined name first
    assert var(ui, "G").items == ["00"]


def test_foreach_over_macro(ui, tmp_path):
    (tmp_path / "n.pre").write_text("N,$1$,$1$,0,0\n")
    m = run(ui, "VAR,L,3,5,8\nLOADMACRO,NODE,n.pre\nFOREACH,L,MACRO,NODE,@L[#]\n")
    assert sorted(m.nodes) == [3, 5, 8]


# ------------------------------------------------------------------ UT-11 REDUCESET and random lists
def test_ut11_reduceset(ui):
    run(ui, "VAR,S,10,9,1,9\nREDUCESET,S\n")
    assert var(ui, "S").items == ["1", "10", "9"]
    run(ui, "REDUCESET,S,INT")
    assert var(ui, "S").items == ["1", "9", "10"]
    run(ui, "VAR,T,1,01,1.0,2.5\nREDUCESET,T,FLOAT\n")
    assert var(ui, "T").items == ["1.0", "2.5"]
    run(ui, "VAR,U,1,a\nREDUCESET,U,INT\nREDUCESET,U,WORD\n")
    assert var(ui, "U").items == ["1", "a"]          # unchanged on error
    assert any("variable unchanged" in e for e in errors(ui))
    assert any("STRING, INT or FLOAT" in e for e in errors(ui))


def test_ut11_rnd_reproducible(tmp_path):
    vals = []
    for _ in range(2):
        ui = Interpreter(cwd=tmp_path)
        ui.run_text("RNDSEED,7\nRND,R,1000,UNI,2,4\n")
        vals.append([float(v) for v in var(ui, "R").items])
    a = np.array(vals[0])
    assert vals[0] == vals[1] and len(a) == 1000
    assert a.min() >= 2.0 and a.max() <= 4.0
    assert abs(a.mean() - 3.0) < 0.1
    ui.run_text("RNDSEED,8\nRND,R,1000,UNI,2,4\n")
    assert [float(v) for v in var(ui, "R").items] != vals[0]


def test_rnd_distributions_and_addrnd(ui):
    run(ui, "RNDSEED,11\nRND,I,500,UNINT,1,3\nRND,N,4000,NORM,10,2\nRND,G,4000,LOGNORM,5,1\nRND,P,2000,POSSION,4\n")
    i = [int(v) for v in var(ui, "I").items]
    assert set(i) == {1, 2, 3}
    n = np.array(var(ui, "N").items, float)
    assert abs(n.mean() - 10) < 0.15 and abs(n.std() - 2) < 0.15
    g = np.array(var(ui, "G").items, float)
    assert abs(g.mean() - 5) < 0.1 and abs(g.std() - 1) < 0.1 and g.min() > 0
    p = np.array(var(ui, "P").items, int)
    assert abs(p.mean() - 4) < 0.2
    run(ui, "SETVAR,@I=7\nADDRND,I,503,UNIINT,7,7\nADDRND,I,10,UNI,0,1\nADDRND,NEW,2,UNI,0,1\n")
    assert len(var(ui, "I").items) == 503 and var(ui, "I").items[-3:] == ["7", "7", "7"]
    assert var(ui, "I").counter == 7                 # ADDRND keeps the counter
    assert any("nothing appended" in w for w in warnings(ui))
    assert len(var(ui, "NEW").items) == 2
    run(ui, "RNDSEED,0\nRND,Q,3,CAUCHY,1,1\n")
    assert any("positive integer" in e for e in errors(ui))
    assert any("distribution must be" in e for e in errors(ui))


def test_loadvar_showvar_varlist(ui, tmp_path):
    (tmp_path / "nodes.txt").write_text("  101 \n\n102\n103\n")
    run(ui, "LOADVAR,nodes.txt\nLOADVAR,nodes.txt,other\nVAR,ONE,x\nSHOWVAR,NODES\nVARLIST\nSHOWVAR,ZZ\n")
    assert var(ui, "NODES").items == ["101", "102", "103"] and var(ui, "OTHER").items == ["101", "102", "103"]
    info = ui.sink.texts(Kind.INFO)
    assert "  [2] 102" in info
    assert any(t.startswith("ONE = x") for t in info)
    assert any(t.startswith("NODES : 3 elements") for t in info)
    assert any("ZZ not defined" in e for e in errors(ui))


def test_loadmacro_and_loadvar_ignore_utf8_bom(ui, tmp_path):
    """Spec 04 4.1: macro and variable files saved with a UTF-8 byte-order mark load unchanged."""
    (tmp_path / "m.pre").write_bytes("N,$1$,1,0,0\r\n".encode("utf-8-sig"))
    (tmp_path / "xs.txt").write_bytes("10\r\n20\r\n".encode("utf-8-sig"))
    m = run(ui, "LOADMACRO,MK,m.pre\nMACRO,MK,7\nLOADVAR,xs.txt\nVAR,K\nFOREACH,XS,N,@K++,@XS[#]")
    assert not ui.sink.texts(Kind.ERROR)
    assert ui.variables["xs"].items == ["10", "20"]
    assert m.nodes[7].x == 1.0 and m.nodes[1].x == 10.0 and m.nodes[2].x == 20.0
