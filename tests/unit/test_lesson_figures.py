"""Concept figures of the guided course: the ```figure blocks of the lessons (sassi/ui/lessons.py, learn.py),
their registry in sassi/ui/static/figures.js, the web build that ships them, and -- in Node.js, when it is
installed -- the physics the figures draw, checked against the Python modules of SASSI-EDU (SHAKE recursion,
response spectra, transfer-function interpolation, Masing loops) and against the numbers the lessons quote.
Format: docs/internal/lesson_format.md ("Figures")."""
from __future__ import annotations

import importlib.util
import json
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from sassi.ui import lessons as L

ROOT = Path(__file__).resolve().parents[2]
STATIC = ROOT / "sassi" / "ui" / "static"
FIGJS = STATIC / "figures.js"

SAMPLE = """---
id: 99-figures
title: Figure sample
part: Fundamentals
order: 99
---
Intro.

```figure
substructuring case=surface
Caption of the intro figure with $X_{ff}$.
```

## A step
Narrative.

```figure
soil-column f=7.3 vhs=1000
The column; resonance near $V_s/(4H)$.
```

### Technical basis
Text.

```figure
nosuch x=1
bad
```

```figure
soil-column f=abc wrong=1 vhs
```
"""


# ---------------------------------------------------------------------- parsing and validation
def test_parse_figure_block():
    fig = L.parse_figure("\nsoil-column f=7.3 vhs=1000\nLine one of the caption,\nline two $x$.\n")
    assert fig.name == "soil-column" and fig.params == {"f": "7.3", "vhs": "1000"}
    assert fig.caption == "Line one of the caption,\nline two $x$." and fig.problems == []
    assert L.parse_figure("interaction-sets method=evbn").problems == []
    assert L.parse_figure("interaction-sets").caption == ""


@pytest.mark.parametrize("text, needle", [
    ("", "without a figure name"),
    ("nosuch", "unknown figure 'nosuch'"),
    ("soil-column wrong=1", "unknown parameter 'wrong'"),
    ("soil-column f=abc", "f=abc is not a number"),
    ("soil-column f", "'f' is not key=value"),
    ("interaction-sets method=all", "method must be one of fv, fsin, evbn, ffv"),
])
def test_parse_figure_problems(text, needle):
    assert any(needle in p for p in L.parse_figure(text).problems), L.parse_figure(text).problems


def test_figure_blocks_in_lessons_and_validation():
    les = L.parse_lesson(SAMPLE, "sample.md")
    kinds = [b.kind for b in les.steps[0].blocks]
    assert kinds == ["figure", "figure", "figure"]
    assert les.steps[0].commands == [] and les.steps[0].actions == []           # figures are neither
    probs = L.validate_lesson(les)
    assert any("unknown figure 'nosuch'" in p for p in probs)
    assert any("f=abc is not a number" in p for p in probs) and any("unknown parameter 'wrong'" in p for p in probs)
    assert not [p for p in probs if "introduction" in p]                         # the intro figure is valid
    good = SAMPLE.split("### Technical basis")[0]
    assert L.validate_lesson(L.parse_lesson(good, "good.md")) == []


def test_render_lesson_places_figures():
    from sassi.ui import learn
    from sassi.ui.helpdocs import HelpDocs
    les = L.parse_lesson(SAMPLE.split("### Technical basis")[0], "sample.md")
    d = learn.render_lesson(les, HelpDocs(ROOT))
    figs = [b for b in d["intro_blocks"] if b["kind"] == "figure"]
    assert len(figs) == 1 and figs[0]["figure"] == "substructuring" and figs[0]["params"] == {"case": "surface"}
    assert f'<div class="lesson-figure" data-figure-block="{figs[0]["k"]}"></div>' in d["intro_html"]
    assert 'data-tex="X_{ff}"' in figs[0]["caption_html"] and "```" not in d["intro_html"]
    st = d["steps"][0]
    sf = [b for b in st["blocks"] if b["kind"] == "figure"]
    assert len(sf) == 1 and sf[0]["figure"] == "soil-column" and sf[0]["params"] == {"f": "7.3", "vhs": "1000"}
    assert f'data-figure-block="{sf[0]["k"]}"' in st["narrative_html"] and "V_s/(4H)" in sf[0]["caption_html"]
    assert "soil-column f=7.3" not in st["narrative_html"]                       # not shown as code


# ---------------------------------------------------------------------- the course
LESSONS = L.load_lessons()


def _figure_blocks(les):
    blocks = list(L._split_blocks(les.intro)[1]) + [b for s in les.steps for b in s.blocks]
    return [L.parse_figure(b.text) for b in blocks if b.kind == "figure"]


USED = sorted({(les.id, f.name) for les in LESSONS for f in _figure_blocks(les)})


def test_course_uses_every_figure_and_only_known_ones():
    names = {n for _, n in USED}
    assert names == set(L.FIGURES), (names ^ set(L.FIGURES))
    for les in LESSONS:
        for f in _figure_blocks(les):
            assert f.problems == [], (les.id, f.problems)
            assert f.caption.strip(), f"{les.id}: figure {f.name} needs a caption"
            assert 1 <= len(re.findall(r"[.!?](?:\s|$)", f.caption)) <= 4, (les.id, f.name, "captions are 1-3 sentences")


def test_figure_captions_typeset_with_katex():
    from sassi.ui import mathcheck
    if mathcheck.node_path() is None:
        pytest.skip("Node.js is not installed")
    errs = []
    for les in LESSONS:
        for f in _figure_blocks(les):
            errs += mathcheck.check_text(f.caption, f"{les.id}:{f.name}")
    assert errs == [], "\n".join(errs)


def test_front_end_wiring():
    js, learn, index = FIGJS.read_text(encoding="utf-8"), (STATIC / "learn.js").read_text(encoding="utf-8"), (STATIC / "index.html").read_text(encoding="utf-8")
    css = (STATIC / "styles.css").read_text(encoding="utf-8")
    order = [index.index(f'"static/{f}"') for f in ("plots.js", "figures.js", "learn.js", "main.js")]
    assert order == sorted(order)                                   # figures.js before learn.js (relative URL)
    assert "S.Figures.fill(root, blocks" in learn
    assert "fill, mount, names" in js and not re.search(r"\.(innerHTML|outerHTML)\s*=|insertAdjacentHTML", js)
    assert "IntersectionObserver" in js and "prefers-reduced-motion" in js
    for name in L.FIGURES:
        assert f'FIG["{name}"] = {{' in js or f"FIG.{name} = {{" in js, name
    for var in ("--fig-ink", "--fig-s1", "--fig-soil1", ".lfig-cap", ".lfig-ctrls"):
        assert var in css, var


def test_web_build_ships_figures(tmp_path):
    pytest.importorskip("plotly")
    spec = importlib.util.spec_from_file_location("sassi_web_build_figs", ROOT / "web" / "build.py")
    b = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(b)
    out = tmp_path / "site"
    summary = b.build(out)
    assert "static/figures.js" in summary["files"]
    assert (out / "static" / "figures.js").read_bytes() == FIGJS.read_bytes()
    v = b.file_version(FIGJS.read_bytes())                     # the version of the file itself (D-W6-11)
    assert summary["versions"]["static/figures.js"] == v
    assert f'src="static/figures.js?v={v}"' in (out / "index.html").read_text(encoding="utf-8")
    shutil.rmtree(out, ignore_errors=True)


# ---------------------------------------------------------------------- the figures in Node.js
NODE = shutil.which("node")
needs_node = pytest.mark.skipif(NODE is None, reason="Node.js is not installed")

_RUNNER = r"""
globalThis.SASSI = {};
require(process.argv[1]);
const F = SASSI.Figures, P = F.physics, C = P.C;
let data = "";
process.stdin.on("data", (c) => { data += c; });
process.stdin.on("end", () => {
  const job = JSON.parse(data);
  const out = new Function("F", "P", "C", job.code)(F, P, C);
  process.stdout.write(JSON.stringify(out));
});
"""


def node(code: str):
    """Run JavaScript `code` (a function body returning a JSON-able value) with figures.js loaded."""
    r = subprocess.run([NODE, "-e", _RUNNER, str(FIGJS)], input=json.dumps({"code": code}), capture_output=True,
                       text=True, timeout=120)
    assert r.returncode == 0, r.stderr[-2000:]
    return json.loads(r.stdout)


@needs_node
def test_registry_matches_python():
    reg = node("return Object.fromEntries(Object.entries(F.registry).map(([k, d]) => [k, {params: d.params, "
               "tex: (d.tex || []).concat(d.texParts || []), animated: !!d.animated, title: d.title}]));")
    assert set(reg) == set(L.FIGURES)
    for name, spec in L.FIGURES.items():
        js = reg[name]["params"]
        assert set(js) == set(spec), name
        for k, kind in spec.items():
            if kind == "num":
                assert js[k]["min"] <= js[k]["def"] <= js[k]["max"], (name, k)
            else:
                assert sorted(js[k]["values"]) == sorted(kind) and js[k]["def"] in kind, (name, k)
        assert reg[name]["title"]
    # every formula a figure shows typesets with the GUI's KaTeX
    tex = [t for d in reg.values() for t in d["tex"]]
    assert len(tex) >= 20
    from sassi.ui import mathcheck
    errs = mathcheck.check_formulas([{"tex": t, "display": False} for t in tex])
    assert [(t, e) for t, e in zip(tex, errs) if e] == []


@needs_node
def test_soil_column_against_shake_and_lesson2():
    from sassi.core.shake import complex_velocity, wave_amplitudes
    res = node("""
      const fs = [0.5, 3, 7.3, 12, 20], col = P.COLUMN;
      const waves = fs.map((f) => { const w = P.columnWaves(col, f); return {E: w.E, F: w.F}; });
      let io = 0, iw = 0, fo = 0, fw = 0;
      for (let f = 0.05; f < 11; f += 0.01) {
        const r = P.columnRatios(col, f), o = C.abs(r.outcrop), w = C.abs(r.within);
        if (o > io) { io = o; fo = f; } if (w > iw) { iw = w; fw = f; }
      }
      const w20 = P.columnWaves(col, 20), u = (z) => C.abs(P.columnU(w20, col, z));
      return {fs, waves, io, fo, iw, fw, u20: [u(0), u(1.6), u(12.8), u(55)]};""")
    g = 32.2                                                    # ft, kip, s
    th, rho = np.array([16.0, 39.0, 1.0]), np.array([0.120, 0.125, 0.130]) / g
    vs, beta = np.array([1000.0, 1650.0, 3300.0]), np.array([0.05, 0.04, 0.02])
    E, F = wave_amplitudes(np.array(res["fs"]), th, rho, complex_velocity(rho * vs ** 2, rho, beta))
    Ej = np.array([[complex(*e) for e in w["E"]] for w in res["waves"]])
    Fj = np.array([[complex(*e) for e in w["F"]] for w in res["waves"]])
    assert np.abs(E - Ej).max() < 1e-10 and np.abs(F - Fj).max() < 1e-10
    # lesson 2: outcrop amplification 2.16 at 7.4 Hz, within ratio 17.1 at 7.3 Hz
    assert res["io"] == pytest.approx(2.163, abs=0.01) and res["fo"] == pytest.approx(7.45, abs=0.05)
    assert res["iw"] == pytest.approx(17.07, abs=0.1) and res["fw"] == pytest.approx(7.31, abs=0.05)
    # lesson 2, SITE listing at 20 Hz: 1 at the surface, 0.980 at 1.6 ft, 0.088 at 12.8 ft, 0.370 at 55 ft
    assert res["u20"] == pytest.approx([1.0, 0.98, 0.088, 0.37], abs=0.006)


@needs_node
def test_response_spectrum_against_motion_module():
    from sassi.core.spectra import response_spectrum
    res = node("""
      const mo = P.synthMotion(), fs = [0.5, 1, 2, 3.5, 5, 8, 10];
      return {a: Array.from(mo.a), dt: mo.dt, n: mo.n, fs, sa: P.spectrum(mo.a, mo.dt, fs, 0.05),
              pga: Math.max(...Array.from(mo.a).map(Math.abs))};""")
    assert res["pga"] == pytest.approx(0.3) and res["n"] == 1000
    py = response_spectrum(np.array(res["a"]), res["dt"], res["fs"], [0.05])["SA"][0]
    assert np.allclose(res["sa"], py, rtol=1e-9)


@needs_node
def test_interpolation_and_critfreq_against_motion_module():
    from sassi.core.interp import interpolate_tf
    res = node("""
      const DF = 1 / (8192 * 0.005), out = {};
      const sets = {lesson: [4,20,41,61,82,102,123,143,164,184,205,225,246,287,328,369,410,492,573,655,737,819],
                    coarse: [4,41,82,123,164,205,287,369,492,655,819]};
      const modes = {"2": [{f: 3.47, beta: 0.052, P: 1.34}, {f: 13.0, beta: 0.10, P: -0.34}],
                     "3": [{f: 3.47, beta: 0.052, P: 1.0}, {f: 4.6, beta: 0.04, P: 0.34}, {f: 13.0, beta: 0.10, P: -0.34}]};
      for (const m of ["2", "3"]) for (const s of ["lesson", "coarse"]) {
        const fu = sets[s].map((n) => n * DF), Hu = fu.map((f) => P.modalTF(f, modes[m]));
        const fi = []; for (let q = 0; q <= 819; q++) fi.push(q * DF);
        const Hi = P.interpTF(fu, Hu, fi, [1, 0]);
        out[m + s] = {fu, Hu, Hi, crit: P.critfreq(fu, Hu.map(C.abs), fi, Hi.map(C.abs), 5, 50, DF)};
      }
      return out;""")
    for key, v in res.items():
        fu = np.array(v["fu"])
        Hu = np.array([complex(*h) for h in v["Hu"]])
        fi = np.arange(820) / (8192 * 0.005)
        Hp = interpolate_tf(fu, Hu, fi, option=1, h0=1.0)
        Hj = np.array([complex(*h) for h in v["Hi"]])
        assert np.abs(Hp - Hj).max() < 1e-8 * np.abs(Hp).max(), key
    # lesson 10, experiment 1: the coarse set misses 3.47 Hz, the 2-DOF interpolant recovers the peak, CRITFREQ flags it
    c = res["2coarse"]["crit"]
    assert [q["n"] for q in c if q["flag"]] == [142] and c[0]["d"] > 150
    assert not [q for q in res["2lesson"]["crit"] if q["flag"]]


@needs_node
def test_masing_loops_against_hysteresis_module():
    from sassi.core.hysteresis import Backbone, hysteretic_damping, secant_ratio
    from sassi.core.panels import bbcgen_curve
    amps = [5e-5, 2e-4, 1e-3, 4e-3, 1e-2]
    res = node(f"""
      const Pn = P.PANEL, bb = P.bbcgen(Pn.vu, Pn.vcr, Pn.vcr / Pn.gcr), pl = P.polyline(bb.xs, bb.ys);
      const amps = {json.dumps(amps)};
      const hy = P.hyperbolic(1, 0.1);
      return {{bb, xi: amps.map((a) => P.masing(pl, a).xi), ratio: amps.map((a) => P.masing(pl, a).ratio),
              hyp: [0.01, 0.1, 1].map((g) => P.masing(hy, g).xi), gr: P.sandRefStrain()}};""")
    g, v, yi = bbcgen_curve(2821, 846, 846 / 8.82e-5)
    assert np.allclose(res["bb"]["xs"], g, rtol=1e-12) and np.allclose(res["bb"]["ys"], v, rtol=1e-12)
    bk = Backbone(g, v, yi)
    for a, xi, r in zip(amps, res["xi"], res["ratio"]):
        assert xi == pytest.approx(hysteretic_damping(4, bk, a), abs=1e-9)
        assert r == pytest.approx(secant_ratio(4, bk, a), rel=1e-9)
    # hyperbolic backbone: Theory Eq. 18.2 with s = 1
    for gam, xi in zip([0.01, 0.1, 1], res["hyp"]):
        z = gam / 0.1
        assert xi == pytest.approx(2 / np.pi * (2 * (1 + z) * (z - np.log1p(z)) / z ** 2 - 1), rel=1e-12)
    assert res["gr"] == pytest.approx(0.056, abs=0.001)          # G/Gmax = 0.5 of the library Sand curve


@needs_node
def test_ssi_oscillator_and_impedance_numbers():
    res = node("""
      const fs = []; for (let f = 0.1; f <= 10; f += 0.001) fs.push(f);
      const pk = (m, key) => P.peakInfo(fs, fs.map((f) => C.abs(P.ssiResponse(m, f)[key])));
      const m = P.ssiModel(1), mm = P.ssiModel(1, {m0: 0});
      const k = m.k, vm = P.SSI.ffb / Math.sqrt(1 + k / m.Kx + k * P.SSI.h ** 2 / m.Kt);
      return {fb: pk(m, "Hfb"), ssi: pk(m, "Ht"), massless: pk(mm, "Ht"), vm, stiff: pk(P.ssiModel(1000), "Ht")};""")
    b = 0.05
    # fixed base: Theory Eq. 2.5, |H|max = 1/(2 beta sqrt(1 - beta^2)) at f1 sqrt(1 - 2 beta^2)
    assert res["fb"]["amp"] == pytest.approx(1 / (2 * b * np.sqrt(1 - b * b)), rel=1e-4)
    assert res["fb"]["f"] == pytest.approx(4.987 * np.sqrt(1 - 2 * b * b), abs=0.002)
    # on the springs: the Veletsos-Meek estimate of lesson 4 (about 3.9 Hz) for a massless mat
    assert res["vm"] == pytest.approx(3.90, abs=0.01)
    assert res["massless"]["f"] == pytest.approx(res["vm"], rel=0.01) and res["ssi"]["f"] == pytest.approx(3.9, abs=0.05)
    assert 0.04 < res["ssi"]["zeta"] < 0.07                       # about 5 %, as lesson 1 finds
    assert res["stiff"]["f"] == pytest.approx(res["fb"]["f"], abs=0.01)   # a rigid site is a fixed base


@needs_node
def test_interaction_sets_and_broadening():
    res = node("""
      const counts = Object.fromEntries(["fv", "fsin", "evbn", "ffv"].map((m) => [m, P.setCount(m)]));
      const X = P.logspace(1, 20, 200), E = X.map((f) => 1 / (1 + ((f - 5) / 0.4) ** 2));
      const B = P.broaden(X, E, 0.15), at = (f) => P.interpLin(B.x, B.y, f);
      return {counts, plateau: [at(5 * 0.86), at(5), at(5 * 1.14)], outside: at(5 * 0.75), peak: Math.max(...E)};""")
    assert res["counts"] == {"fv": 150, "fsin": 105, "evbn": 114, "ffv": 132}
    text = (L.LESSON_DIR / "05_embedded.md").read_text(encoding="utf-8")
    assert "FV 150, FI-FSIN 105, FI-EVBN 114, FFV\n132" in text or "FV 150, FI-FSIN 105, FI-EVBN 114, FFV 132" in text
    assert res["plateau"] == pytest.approx([res["peak"]] * 3, rel=1e-3) and res["outside"] < 0.9 * res["peak"]


@needs_node
def test_figures_js_syntax():
    r = subprocess.run([NODE, "--check", str(FIGJS)], capture_output=True, text=True, timeout=60)
    assert r.returncode == 0, r.stderr
