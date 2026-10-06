"""AFWRITE: CHECK, then write the input deck of every AOPT-enabled module without errors
(requirements sections 1.9, 2.6, 3.3; D-AFW-01..06; ARCHITECTURE section 3; spec 07 section 9.2.3).

Each deck ``<model>.<ext>`` is created with :func:`sassi.io.decks.new`, every schema parameter and
table is filled with **resolved** data, and it is written with :func:`sassi.io.decks.write`:

* global node coordinates (local systems resolved), ETYPE 0 resolved to 1/2 (D-HOU-01);
* SITE: the TOPL layers with their L properties and the half-space row; HOUSE: the same layers as
  ``sitelayers`` (with a half-space row of thickness 0);
* the frequency numbers of the selected set (SITE ``<freq>``) and the resolved ``df``; NFFT is the
  nearest power of 2 when the model value is not one (Warning 9, D-CNV-10);
* shared variables copied from their one storage location into every deck that needs them
  (spec 07 section 3): HOUSE gravity, SITE delt/nft/fstep/freq/cl/cm, ANALYS type/ang, MOTION
  history data, THFILE/THTIT, DAMP;
* write-time fix-ups (D-AFW-06): gap nodes (Warning 1) and unused or K-only nodes (Warning 4,
  D-CHK-08) are written with every DOF fixed; empty groups (Warning 7) are skipped; only the first
  100 top layers with ``EDUOPT,LIMITS,PREP`` (Warning 8); duplicate MOTION output requests are
  merged with OR of the flags (D-MOT-10, EDU-16).

SOIL-NON (P2): the SOIL deck schema has no nonlinear time-domain fields, so the NLSOIL record and the
NLSLAYER sets are written to the side file ``<model>.nls`` next to ``<model>.soi`` (read by SOIL; a stale
``.nls`` of a model without NLSOIL is renamed ``.nls.bak``), see
:func:`sassi.prep.commands.soilnon_cmds.write_nls_file`.

A module with CHECK errors (or gated by model errors, D-CHK-01) gets no deck; an existing stale
deck is renamed ``<model>.<ext>.bak`` with a warning (D-AFW-01).  A deck whose data cannot be
resolved although CHECK passed blocks only its own module the same way (the reason is listed in the
window and in the ``.err`` file); the other modules still get their decks.  COMBIN has no deck
(D-CMB-02); NONLINEAR (Option NON) is P2.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

from ..io import deckfmt, decks
from . import defaults as DEF
from .check import PLANE, SOLID, CheckOptions, CheckReport, Checker, ModelView, Resolved, write_err
from .options import eduopt, eduopt_float, get_entries

#: modules with an input deck, in run order (COMBIN has none, D-CMB-02)
DECK_MODULES = ("EQUAKE", "SOIL", "SITE", "POINT", "HOUSE", "FORCE", "ANALYS", "MOTION", "STRESS", "RELDISP")


class AfwriteError(RuntimeError):
    """AFWRITE cannot run (no model name / path, D-MDL-03)."""


class DeckBuildError(ValueError):
    """The data a deck needs is not resolvable (normally caught by CHECK first); the module is
    then treated as blocked: no deck, stale deck renamed (D-AFW-01)."""


@dataclass
class AfwriteResult:
    report: CheckReport
    written: Dict[str, Path] = field(default_factory=dict)
    blocked: Dict[str, str] = field(default_factory=dict)
    renamed: Dict[str, Path] = field(default_factory=dict)
    notes: List[str] = field(default_factory=list)      # write-time fix-ups and information
    warnings: List[str] = field(default_factory=list)
    err_file: Optional[Path] = None
    sim_file: Optional[Path] = None
    model_hash: str = ""


# ======================================================================================
# Deck builders
# ======================================================================================
class DeckBuilder:
    """Builds the decks of one model from its resolved analysis data."""

    def __init__(self, model, dirs: Sequence[Union[str, Path]] = (), checker: Optional[Checker] = None):
        self.m = model
        self.chk = checker or Checker(model, dirs=dirs)
        self.r: Resolved = self.chk.r
        self.notes: List[str] = []

    @property
    def view(self) -> ModelView:
        return self.chk.view

    # ---------------------------------------------------------------- common pieces
    def _common(self, d: decks.Deck, opmode: Optional[int] = None) -> None:
        d["title"] = self.m.title
        d["model"] = self.m.name
        if opmode is not None:
            d["opmode"] = int(opmode)

    def _fft(self, d: decks.Deck) -> None:
        r = self.r
        d["delt"] = r.delt
        d["nft"] = r.nft
        d["df"] = r.df if r.df is not None else 0.0
        if r.nft != r.nft_raw:
            self._note(f"NFFT {r.nft_raw} written as {r.nft} (Warning 9, nearest power of 2, D-CNV-10)")

    def _freqs(self, d: decks.Deck) -> None:
        for n in sorted(set(self.r.fnums or [])):
            d.table("freqs").append([n])

    def _hist(self, d: decks.Deck, module: str = "") -> None:
        """Control-motion data shared by MOTION, STRESS and RELDISP (spec 07 section 3).  A blank THFILE is
        written as its built-in default when the module reads the history (D-W5-03 ... D-W5-05, D-W5-13)."""
        r, mo = self.r, self.r.motion
        hc = r.history_choice
        use = hc.default is not None and (not module or r.history_needed(module))
        d["thfile"] = r.deck_file_name(r.thfile if use or hc.default is None else r.thfile_given)
        if use:
            self._note(f"{hc.default.text} (EDU-29); {hc.name} is written to the decks that read the "
                       f"control motion")
        d["thtit"] = r.thtit
        d["mult"] = float(mo.mult)
        d["max"] = float(mo.get("max"))
        d["rec1"] = int(mo.rec1)
        d["rec2"] = int(mo.rec2)
        d["fopt"] = int(mo.fopt)
        d["dur"] = float(mo.dur)
        d["type"] = r.type
        d["gravity"] = r.gravity
        d["file8"] = eduopt(self.m, "TFFILE") or "FILE8"
        d["maxnode"] = max(self.m.nodes) if self.m.nodes else 0

    def _note(self, text: str) -> None:
        if text not in self.notes:
            self.notes.append(text)

    #: simultaneous-case solution files of ANALYS <simul> = 1 and their control direction (D-ANL-06)
    XYZ_FILE8 = {"FILE8X": 0, "FILE8Y": 1, "FILE8Z": 2}

    def environment(self) -> Tuple[int, float]:
        """Control direction ``<cm>`` and angle ``<ang>`` of the solution file read by MOTION, STRESS
        and RELDISP (deck parameters ``cm`` / ``ang``; requirements 2.4, D-ANL-06).

        With ANALYS ``<simul>`` = 1 a simultaneous-case file FILE8X / FILE8Y / FILE8Z
        (``EDUOPT,TFFILE``) holds the response to the x' / y' / z' control motion of FILE1X / FILE1Y /
        FILE1Z -- the SV x', SH y' and P z' SITE runs at angle 0 (EDU-26 forbids an angle with
        simultaneous cases) -- whatever the options of the *last* SITE run, so the direction is taken
        from the file name (a Windows path ``C:\\SSI\\FILE8Y`` included).  The zero-frequency anchor of
        the interpolation (MOTION, STRESS) and the free-field reference of RELDISP (D-RDP-04) depend on
        it.  Any other case -- another file name, or a single-case FILE8 that the user copied to a name
        like FILE8X while ``<simul>`` = 0 -- uses SITE ``<cm>`` and ANALYS ``<ang>``.  MOTION, STRESS
        and RELDISP finally use the direction stored in the solution file itself when they find it.
        """
        raw = (eduopt(self.m, "TFFILE") or "FILE8").strip().replace("\\", "/")
        name = raw.rsplit("/", 1)[-1].upper()
        site_cm, ang = int(self.r.site.cm), float(self.r.analys.ang)
        if name in self.XYZ_FILE8 and int(self.r.analys.simul) == 1:
            cm = self.XYZ_FILE8[name]
            if cm != site_cm or ang != 0.0:
                self._note(f"{name}: control direction {'xyz'[cm]}' at angle 0 written to the MOTION/STRESS/RELDISP "
                           f"decks (simultaneous X/Y/Z case, D-ANL-06; the last SITE run has <cm> = {site_cm}, "
                           f"ANALYS <ang> = {ang:g})")
            return cm, 0.0
        return site_cm, ang

    def _layer_row(self, l: int, thick: Optional[float] = None) -> List[Any]:
        L = self.m.layers.get(l)
        if L is None:
            raise DeckBuildError(f"soil layer {l} is not defined (Error 19)")
        return [l, L.thick if thick is None else thick, L.weight, L.vp, L.vs, L.pdamp, L.sdamp]

    def site_layers(self) -> List[int]:
        topl = list(self.m.topl)
        if len(topl) > 100 and self.r.limits == "PREP":
            self._note(f"only the first 100 of {len(topl)} top layers written (Warning 8, LIMITS,PREP)")
            topl = topl[:100]
        return topl

    def halfspace_layer(self) -> int:
        """L layer written as the half-space row (SITE ``halfspace``, last HOUSE ``sitelayers`` row).

        * No top layers: no layer profile, no half-space row (0; SITE has Error 46).
        * SITE ``<hs>`` defined: that layer.
        * Rigid base without a half-space layer (``<hs>`` = 0 and ``<nl>`` = 0): the deck still
          needs a row; the last top layer is written (its properties are not used by a rigid base).
        * Otherwise (``<hs>`` undefined, or ``<hs>`` = 0 with generated half-space sublayers) there
          is no valid half-space: CHECK reports Error 19 (SITE and HOUSE) and the deck is not written.
        """
        s = self.r.site
        hs = int(s.hs)
        if not self.m.topl:
            return 0
        if hs in self.m.layers:
            return hs
        if hs == 0 and int(s.nl) == 0:
            self._note(f"rigid base (SITE <nl> = 0) without a half-space layer: the last top layer "
                       f"{self.m.topl[-1]} is written as the half-space row (properties not used)")
            return self.m.topl[-1]
        raise DeckBuildError(f"SITE half-space layer <hs> = {hs} is not defined (Error 19)")

    # ================================================================== EQUAKE
    def equake(self) -> decks.Deck:
        m, r = self.m, self.r
        eq = r.rec["EQUAKE"]
        d = decks.new("EQUAKE")
        self._common(d)
        d["opmode"] = 0
        files, uses = r.equake_files()                 # blank inputs with their defaults (D-W5-07, D-W5-08)
        for u in uses:
            self._note(f"EQUAKE: {u.text} (EDU-29)")
        nrfreq = int(eq.nrfreq)
        if not eq.given(2):
            from .check import count_xy_rows
            rs = [k for k, f in sorted(files["RSIN"].items()) if f]
            p = r.file(files["RSIN"][rs[0]]) if rs else None
            nrfreq = (count_xy_rows(p) or 0) if p is not None else 0
            if rs:
                self._note(f"EQUAKE <nrfreq> = {nrfreq} (records of RSIN {rs[0]}, dialog default)")
        d["accopt"] = int(eq.accopt)
        d["nrfreq"] = nrfreq
        d["rand"] = int(eq.rand)
        d["damp"] = float(eq.damp)
        d["dur"] = float(eq.dur)
        d["corr"] = int(eq.corr)
        d["seeds"] = int(eq.seeds)
        d["tpsd"] = int(eq.tpsd)
        d["delt"] = r.delt
        d["gravity"] = r.gravity
        d["eqtit"] = m.options.string("EQTIT")
        comps = sorted({k for dd in files.values() for k in dd if isinstance(k, int)})
        for i in comps:
            def f(name: str, inp: bool) -> str:
                v = files[name].get(i, "")
                return r.deck_file_name(v) if inp else v
            d.table("spectra").append([i, f("RSIN", True), f("RSOUT", False), f("ACCIN", True), f("ACCOUT", False),
                                       f("TPSD", True)])
        for k, rec in get_entries(m, "CORR"):
            d.table("corr").append([int(rec.no), float(rec.time), float(rec.val)])
        return d

    # ================================================================== SOIL
    def soil(self) -> decks.Deck:
        m, r = self.m, self.r
        so, sx = r.rec["SOIL"], r.rec["SOILX"]
        d = decks.new("SOIL")
        self._common(d)
        d["opmode"] = 0
        for name in ("nrval", "header", "outcrop", "save", "iter"):
            d[name] = int(so.get(name))
        for name in ("grav", "ratio", "gravmult"):
            d[name] = float(so.get(name))
        d["cof"] = 1.0 / (2.0 * r.delt) if r.delt > 0 else 0.0          # cut-off forced to Nyquist (manual)
        d["soilcutoff"] = eduopt_float(m, "SOILCUTOFF", 0.0)          # Fortran D exponents (L6, D-PAR-07)
        d["vppolicy"] = (eduopt(m, "VPPOLICY") or "NU").upper()         # D-SOL-06
        d["delt"] = r.delt
        d["nft"] = r.nft
        d["cl"] = int(sx.cl) if int(sx.cl) > 0 else int(r.site.cl)
        hc = r.soil_history()                          # SOILX file, THFILE or the built-in record (D-W5-06)
        d["thfile"] = r.deck_file_name(hc.name)
        for u in r.defaults_used("SOIL"):
            self._note(f"SOIL: {u.text} (EDU-29)")
        d["thtit"] = r.thtit
        d["mult"] = float(sx.mult)
        d["max"] = float(sx.get("max"))
        d["indir"] = int(sx.indir)
        d["cmodform"] = r.cmodform
        prof = self.chk.soil_profile()                 # SPRO, or the default profile (D-W5-10)
        for k, rec in prof:
            row = self._layer_row(int(rec.prop))
            d.table("profile").append([int(rec.layer)] + row[1:] + [rec.dynprop])
        used = list(dict.fromkeys(rec.dynprop for _, rec in prof if rec.dynprop))
        # model curves, and the built-in curve of a used label the model does not define (D-W5-09)
        rows = DEF.dynp_table(m, used)
        for lab in dict.fromkeys(row[0] for row in rows):
            pts = [row for row in rows if row[0] == lab]
            for row in pts:
                if row[6]:
                    d.table("dynp").append([lab, "G", row[2], row[3]])
            for row in pts:
                if row[7]:
                    d.table("dynp").append([lab, "D", row[4], row[5]])
        for k, rec in get_entries(m, "SACC"):
            d.table("sacc").append([int(rec.layer), int(rec.opt), int(rec.outcrop)])
        for k, rec in get_entries(m, "SRS"):
            d.table("srs").append([int(rec.layer), int(rec.save), int(rec.outcrop)])
        for k, rec in get_entries(m, "SSTR"):
            d.table("sstr").append([int(rec.layer), int(rec.opt1), int(rec.opt2), int(rec.opt3), int(rec.opt4)])
        for k, rec in get_entries(m, "SSAF"):
            d.table("ssaf").append([int(rec.layer), int(rec.save), int(rec.outcrop1), int(rec.outcrop2),
                                    int(rec.layer2), float(rec.freqstep), rec.title])
        if get_entries(m, "SFOU"):
            self._note("SFOU requests are not usable in this version and are not written (compute Fourier "
                       "spectra with EQUAKE)")
        for v in m.damp:
            d.table("damp").append([float(v)])
        return d

    # ================================================================== SITE
    def site(self) -> decks.Deck:
        m, r = self.m, self.r
        s = r.site
        d = decks.new("SITE")
        self._common(d, s.opmode)
        for name in ("mode1", "nl", "mode2", "wopt", "freq1", "cl", "cm", "freq"):
            d[name] = int(s.get(name))
        d["fstep"] = r.fstep
        d["freq2"] = int(s.freq2) if s.given(9) else max(r.nft // 2, 1)
        hs = self.halfspace_layer()
        d["hs"] = int(s.hs)
        d["soilmode"] = int(r.rec["SITEX"].soilmode)
        d["gravity"] = r.gravity
        d["cmodform"] = r.cmodform
        d["hslaw"] = eduopt(m, "HSLAW").lower()
        self._fft(d)
        for l in self.site_layers():
            d.table("layers").append(self._layer_row(l))
        if hs:
            d.table("halfspace").append(self._layer_row(hs))
        if not get_entries(m, "WAVE"):
            self._note("no WAVE command: the new-model wave field is written (vertical "
                       f"{'SV' if int(s.wopt) == 0 else 'SH'}, ratios 1 and 1)")
        for w in self.chk.effective_waves():
            d.table("waves").append([int(w.type), int(w.opt), float(w.ratio1), float(w.ratio2), float(w.angle)])
        self._freqs(d)
        return d

    # ================================================================== POINT
    def point(self) -> decks.Deck:
        r = self.r
        p = r.rec["POINT"]
        d = decks.new("POINT")
        self._common(d, p.opmode)
        d["layer"] = int(p.layer)
        d["rad"] = float(p.rad)
        d["dim"] = int(r.house.dim)
        d["freq"] = r.freq_set
        d["df"] = r.df if r.df is not None else 0.0
        self._freqs(d)
        return d

    # ================================================================== HOUSE
    def house(self) -> decks.Deck:
        m, r = self.m, self.r
        v = self.view
        h = r.house
        d = decks.new("HOUSE")
        self._common(d, h.opmode)
        d["gravity"] = r.gravity
        d["gelev"] = r.gelev
        for name in ("dim", "imp", "coh", "wpass", "me", "cmplxspec"):
            d[name] = int(h.get(name))
        hx = r.rec["HOUSEX"]
        for name in ("optimize", "supmode", "nsim", "nlssi", "ansys"):
            d[name] = int(hx.get(name))
        mo = m.mopt
        d["incomp"] = int(mo.incomp)
        d["gmunits"] = int(mo.matrix)
        d["cmodform"] = r.cmodform
        inc, wp, an = r.rec["INCOH"], r.rec["WPASS"], r.analys
        for name in ("gammax", "gammay", "gammaz", "alpha", "randphz"):
            d[name] = float(inc.get(name))
        for name in ("ngp", "ipr", "nmodes", "met", "hseed", "vseed"):
            d[name] = int(inc.get(name))
        d["appv"] = float(wp.appv)
        d["wang"] = float(wp.ang)
        d["cohf"] = int(wp.cohf)
        d["xc"], d["yc"], d["zc"] = float(an.xc), float(an.yc), float(an.zc)
        self._fft(d)
        # ---- nodes (global coordinates); gap / unused / K-only nodes fixed (W1, W4, D-CHK-08)
        inter = set(v.interaction_nodes())
        used = v.used_dof_nodes()
        ids = sorted(m.nodes)
        gaps = [n for n in range(1, ids[-1]) if n not in m.nodes] if ids else []
        rows: Dict[int, List[Any]] = {}
        nfix = 0
        for n in ids:
            p = v.P[n]
            fx = list(m.nodes[n].fix)
            if n not in used and n not in inter:
                fx = [1] * 6
                nfix += 1
            rows[n] = [n, float(p[0]), float(p[1]), float(p[2])] + [int(x) for x in fx]
        for n in gaps:
            rows[n] = [n, 0.0, 0.0, 0.0] + [1] * 6
        for n in sorted(rows):
            d.table("nodes").append(rows[n])
        if gaps:
            self._note(f"{len(gaps)} gap nodes written with all DOFs fixed (Warning 1)")
        if nfix:
            self._note(f"{nfix} unused or K-only nodes written with all DOFs fixed (Warning 4, D-CHK-08)")
        for n in sorted(inter):
            d.table("interaction").append([n])
        # ---- groups and elements (ETYPE resolved; empty groups skipped)
        skipped = []
        by_group: Dict[int, List[Any]] = {}
        for ref in v.elems:
            by_group.setdefault(ref.group, []).append(ref)
        for gid in sorted(m.groups):
            g = m.groups[gid]
            if not g.elements:
                skipped.append(gid)
                continue
            d.table("groups").append([gid, g.type, g.title])
            for ref in by_group.get(gid, []):
                e = ref.elem
                nodes = (list(e.nodes) + [0] * 8)[:8]
                d.table("elements").append([gid, e.id, int(ref.etype), int(e.mat), int(e.prop), int(e.eint),
                                            float(e.thick)] + [int(x) for x in nodes]
                                           + ["".join(str(int(bool(x))) for x in e.ki),
                                              "".join(str(int(bool(x))) for x in e.kj)])
        if skipped:
            self._note(f"empty groups {skipped} skipped (Warning 7)")
        n_res = sum(1 for ref in v.elems if ref.type in (SOLID, PLANE) and ref.elem.etype == 0)
        if n_res:
            n_exc = sum(1 for ref in v.elems if ref.type in (SOLID, PLANE) and ref.elem.etype == 0 and ref.etype == 2)
            self._note(f"ETYPE 0 resolved for {n_res} SOLID/PLANE elements: {n_exc} excavated, "
                       f"{n_res - n_exc} structure (D-HOU-01)")
        # ---- property tables
        for k, mt in sorted(m.materials.items()):
            d.table("materials").append([k, int(mt.mtype), mt.val1, mt.val2, mt.weight, mt.pdamp, mt.sdamp])
        for k in sorted(m.layers):
            d.table("layers").append(self._layer_row(k))
        for l in self.site_layers():
            d.table("sitelayers").append(self._layer_row(l))
        hs = self.halfspace_layer()
        if hs:
            d.table("sitelayers").append(self._layer_row(hs, thick=0.0))
        for k, s in sorted(m.sections.items()):
            d.table("beamprops").append([k, s.axial, s.shear2, s.shear3, s.tors, s.flex2, s.flex3])
        for k, s in sorted(m.springs.items()):
            d.table("springprops").append([k, s.scx, s.scy, s.scz, s.scxx, s.scyy, s.sczz, s.damp])
        for k, p in sorted(m.matrices.items()):
            for kind in ("R", "I", "M"):
                for row, terms in sorted(p.rows.get(kind, {}).items()):
                    t = (list(terms) + [0.0] * 12)[:12]
                    d.table("matrices").append([k, kind, int(row)] + [float(x) for x in t])
        # ---- masses (raw values with their units flag; HOUSE converts with the gravity)
        for n in sorted(set(m.tmass) | set(m.rmass)):
            t = m.tmass.get(n, [0.0, 0.0, 0.0])
            q = m.rmass.get(n, [0.0, 0.0, 0.0])
            d.table("masses").append([n] + [float(x) for x in t] + [float(x) for x in q] + [m.mass_unit(n)])
        # ---- multiple excitation, amplification ratios, symmetry
        for k, z in get_entries(m, "ME"):
            d.table("me").append([int(z.no), int(z.nfirst), int(z.nlast)])
        cplx = int(h.cmplxspec) == 1
        for no, vals in sorted(m.amp.items()):
            if cplx:
                for i in range(0, len(vals), 2):
                    d.table("amp").append([no, i // 2 + 1, float(vals[i]),
                                           float(vals[i + 1]) if i + 1 < len(vals) else 0.0])
            else:
                for i, a in enumerate(vals):
                    d.table("amp").append([no, i + 1, float(a), 0.0])
        for k, s in get_entries(m, "SYMM"):
            d.table("symm").append([int(s.no), int(s.type), int(s.node1), int(s.node2), int(s.node3)])
        self._freqs(d)
        return d

    # ================================================================== FORCE
    def force(self) -> decks.Deck:
        m, r = self.m, self.r
        f = r.rec["FORCE"]
        d = decks.new("FORCE")
        self._common(d, f.opmode)
        d["gravity"] = r.gravity
        d["mforce"] = int(m.mopt.force)
        d["freq"] = r.freq_set
        d["fstep"] = r.fstep
        self._fft(d)
        dropped = 0
        for table, off in ((m.forces, 0), (m.moments, 3)):
            for n, ld in sorted(table.items()):
                nd = m.nodes.get(n)
                for k in range(3):
                    if ld.factor[k] == 0:
                        continue
                    if nd is not None and nd.fix[off + k]:
                        dropped += 1
                        continue
                    d.table("loads").append([n, off + k + 1, float(ld.factor[k]), float(ld.arrival[k])])
        if dropped:
            self._note(f"{dropped} load components on fixed DOFs ignored (Warnings 10/11)")
        self._freqs(d)
        return d

    # ================================================================== ANALYS
    def analys(self) -> decks.Deck:
        r = self.r
        a = r.analys
        ax = r.rec["ANALYSX"]
        d = decks.new("ANALYS")
        self._common(d, a.opmode)
        for name in ("type", "mode", "save", "prnt", "fopt", "impe", "simul"):
            d[name] = int(a.get(name))
        for name in ("ang", "xc", "yc", "zc"):
            d[name] = float(a.get(name))
        d["ffm"] = int(ax.ffm)
        d["delrst"] = int(ax.delrst)
        d["coh"], d["wpass"], d["me"] = r.coh, r.wpass, r.me
        d["freq"] = r.freq_set
        d["gravity"] = r.gravity
        self._fft(d)
        self._freqs(d)
        return d

    # ================================================================== MOTION
    def motion(self) -> decks.Deck:
        m, r = self.m, self.r
        mo, mx = r.motion, r.rec["MOTIONX"]
        d = decks.new("MOTION")
        self._common(d, mo.opmode)
        for name in ("out", "step", "res", "fstep", "bl", "cplx", "cnvrt", "pzadj", "interp"):
            d[name] = int(mo.get(name))
        for name in ("freq1", "freq2", "smo"):
            d[name] = float(mo.get(name))
        for name in ("f1213", "resp", "srss", "savetf", "saveacc", "savers", "saverot", "rsttf", "rstacc", "rstrs"):
            d[name] = int(mx.get(name))
        d["cm"], d["ang"] = self.environment()
        self._fft(d)
        self._hist(d, "MOTION")
        for v in m.damp:
            d.table("damp").append([float(v)])
        merged: Dict[Tuple[int, int], List[int]] = {}
        dups = 0
        for req in m.nout:
            codes = [1 if int(c) else 0 for c in (list(req.codes) + [0] * 6)[:6]]
            for n in req.nodes:
                key = (n, int(req.dir))
                if key in merged:
                    dups += 1
                    merged[key] = [a | b for a, b in zip(merged[key], codes)]
                else:
                    merged[key] = list(codes)
        for (n, dr), codes in merged.items():
            d.table("nout").append([n, dr] + codes)
        if dups:
            self._note(f"{dups} duplicate MOTION output requests merged, flags OR'ed (D-MOT-10, EDU-16)")
        return d

    # ================================================================== STRESS
    def stress(self) -> decks.Deck:
        m, r = self.m, self.r
        st, sx = r.rec["STRESS"], r.rec["STRESSX"]
        d = decks.new("STRESS")
        self._common(d, st.opmode)
        for name in ("iter", "save", "itran", "interopt"):
            d[name] = int(st.get(name))
        for name in ("pzadj", "skip", "savemax", "saveth", "rstns", "rstsp"):
            d[name] = int(sx.get(name))
        d["smo"] = float(sx.smo)
        sec = m.options.record("SECDATAOPT")
        d["secdataopt"] = sec.integer(1, 0) if sec is not None else 0
        tsr = m.options.record("THSHLSTR")                         # D-TSH-01 (wave-3 lead fix)
        d["thshlstr"] = 1 if (tsr is not None and tsr.integer(1, 0) == 1) else 0
        d["cm"], d["ang"] = self.environment()
        self._fft(d)
        self._hist(d, "STRESS")
        for req in m.eout:
            codes = [int(c) for c in (list(req.codes) + [0] * 12)[:12]]
            for e in req.elements:
                d.table("eout").append([int(req.group), int(e)] + codes)
        return d

    # ================================================================== RELDISP
    def reldisp(self) -> decks.Deck:
        m, r = self.m, self.r
        rd, rx = r.rec["RELD"], r.rec["RELDX"]
        d = decks.new("RELDISP")
        self._common(d)
        d["opmode"] = 0
        d["relfile"] = m.options.string("RELFILE")
        d["reldisoutput"] = int(rd.reldisoutput)
        d["reldispsall"] = int(rd.reldispsall)
        d["numfiles"] = len(m.rdnd)                     # D-RDP-03: written = RDND count
        d["saverot"] = int(rx.saverot)
        d["rstframes"] = int(rx.rstframes)
        d["cm"], d["ang"] = self.environment()
        self._fft(d)
        self._hist(d, "RELDISP")
        for q in m.rdnd:
            d.table("rdnd").append([int(q.node)] + [1 if int(f) >= 1 else 0 for f in (list(q.flags) + [0] * 6)[:6]])
        return d

    def build(self, module: str) -> decks.Deck:
        return getattr(self, module.lower())()


# ======================================================================================
# Deck writing
# ======================================================================================
#: columns whose schema type is int although the quantity is real (reported to the lead):
#: SOIL ssaf.freqstep is a frequency step in Hz (manual example 0.1)
_FLOAT_OVERRIDES = {("SOIL", "ssaf", "freqstep")}


def _coerce(v: Any, typ: type) -> Any:
    if typ is str:
        return "" if v is None else str(v)
    f = float(v)
    if typ is int:
        if f != round(f):
            raise ValueError(f"expected an integer, got {v!r}")
        return int(round(f))
    return f


def write_deck(path: Union[str, Path], deck: decks.Deck, comments: Sequence[str] = ()) -> Path:
    """Write a deck with :func:`sassi.io.decks.write`; decks with a schema type override
    (:data:`_FLOAT_OVERRIDES`) are written with the same layout through :mod:`sassi.io.deckfmt`."""
    mod = deck.module
    if not any(o[0] == mod for o in _FLOAT_OVERRIDES):
        return decks.write(path, deck, comments=comments)
    s = decks.SCHEMAS[mod]
    unknown = set(deck.params) - {p.name for p in s.params}
    if unknown:
        raise KeyError(f"{mod} deck: unknown parameters {sorted(unknown)}")
    params = {}
    for p in s.params:
        val = deck.params.get(p.name, p.default)
        params[p.name] = [(_coerce(x, p.item)) for x in val] if p.type is list else _coerce(val, p.type)
    tables = {}
    for t in s.tables:
        cols = [(c, float if (mod, t.name, c) in _FLOAT_OVERRIDES else typ) for c, typ in t.columns]
        out = deckfmt.Table(columns=[c for c, _ in cols])
        tab = deck.tables.get(t.name)
        if tab is not None:
            for row in tab.rows:
                out.rows.append([_coerce(v, typ) for v, (_, typ) in zip(row, cols)])
        tables[t.name] = out
    return deckfmt.write_raw(path, mod, params, tables, version=decks.DECK_VERSION, comments=comments)


def _stale(path: Path) -> Optional[Path]:
    """Rename a stale deck to ``<deck>.bak`` (D-AFW-01); returns the new path."""
    if not path.exists():
        return None
    bak = path.with_name(path.name + ".bak")
    os.replace(path, bak)
    return bak


def afwrite(model, dirs: Sequence[Union[str, Path]] = (), check_options: Optional[CheckOptions] = None,
            write_options: Optional[Dict[str, Any]] = None, modules: Optional[Sequence[str]] = None) -> AfwriteResult:
    """CHECK, then write the decks of the enabled modules without errors (D-AFW-01).

    ``dirs`` are the directories for relative input files (model directory first).  ``write_options``
    is the interpreter's Options > Write state: ``sim`` (Simulation Commands, D-AFW-04) writes the
    RUN commands of the written decks in run order to ``<model>-Sim.pre`` (``sim_location`` 'sim',
    default) or to the end of ``<model>.pre`` (``sim_location`` 'pre').  ``modules`` overrides AOPT.
    """
    if not model.name or not model.path:
        raise AfwriteError("model name and path are not defined -- use MDL,<Model>,<Path> (D-MDL-03)")
    mdir = Path(model.path)
    if not mdir.is_dir():
        raise AfwriteError(f"model directory {mdir} does not exist")
    dirs = [mdir] + [Path(d) for d in dirs if d and Path(d) != mdir]
    chk = Checker(model, modules, dirs)
    report = chk.run()
    res = AfwriteResult(report=report, model_hash=model.model_hash())
    builder = DeckBuilder(model, dirs, checker=chk)
    comments = [f"written by AFWRITE from model '{model.name}'", f"model_hash = {res.model_hash}"]

    def block(mod: str, path: Path, why: str) -> None:
        res.blocked[mod] = why
        bak = _stale(path)
        if bak is not None:
            res.renamed[mod] = bak
            res.warnings.append(f"{mod}: stale deck {path.name} renamed {bak.name} (D-AFW-01)")

    try:
        for mod in report.modules:
            if mod == "COMBIN":
                # no deck (D-CMB-02); the duplicate-frequency policy goes to COMBIN.opt (sassi.modules.combin)
                if report.blocked(mod):
                    res.blocked[mod] = report.blocking_reason(mod)
                    continue
                policy = eduopt(model, "COMBINDUP")
                (mdir / "COMBIN.opt").write_text(f"EDUOPT,COMBINDUP,{policy}\n", encoding="utf-8")
                res.notes.append(f"COMBIN has no input deck (D-CMB-02): it reads FILE81 and FILE82; COMBIN.opt "
                                 f"written (COMBINDUP {policy}, D-CMB-01)")
                continue
            if mod == "NONLINEAR":
                # Option NON (P2): keyword .eql deck (sassi.modules.nonlinear, D-NON-01)
                from .commands.nonlinear_cmds import build_eql
                from ..modules.nonlinear import write_eql
                path = mdir / f"{model.name}.eql"
                eql, errs, warns = build_eql(model)
                res.warnings += [f"NONLINEAR: {w}" for w in warns]
                if report.blocked(mod) or errs:
                    block(mod, path, report.blocking_reason(mod) if report.blocked(mod) else "; ".join(errs))
                else:
                    res.written[mod] = write_eql(path, eql, comments=comments)
                continue
            path = decks.deck_path(mdir, model.name, mod)
            if report.blocked(mod):
                block(mod, path, report.blocking_reason(mod))
                continue
            # A deck that cannot be built blocks only its own module: the other decks are still
            # written and the reason goes to the window and the .err file (D-AFW-01).
            try:
                deck = builder.build(mod)
                res.written[mod] = write_deck(path, deck, comments=comments)
                if mod == "HOUSE" and int(deck["nlssi"]) == 1:
                    write_pin_file(model, mdir, res, view=chk.view)
                if mod == "HOUSE" and int(deck["coh"]) == 1:
                    write_house_options(model, mdir, res)
                if mod == "SOIL":
                    # SOIL-NON (P2): NLSOIL / NLSLAYER go to <model>.nls (no SOIL deck fields for them)
                    from .commands.soilnon_cmds import write_nls_file
                    write_nls_file(model, mdir, res)
            except Exception as exc:          # noqa: BLE001 -- isolate one module's failure
                why = (str(exc) if isinstance(exc, DeckBuildError)
                       else f"internal error while building the deck: {exc!r}")
                block(mod, path, why)
                report.notes.append(f"AFWRITE: {mod} input deck not written -- {why}")
    finally:
        res.err_file = write_err(report, mdir / f"{model.name}.err", check_options)
    res.notes.extend(builder.notes)
    wo = write_options or {}
    if wo.get("sim"):
        res.sim_file = write_simulation_commands(model, report, res, str(wo.get("sim_location", "sim")))
    return res


def write_house_options(model, mdir: Path, res: AfwriteResult) -> Path:
    """Incoherent analysis (HOUSE <coh> = 1, requirements 4.4 item 6): the EDUOPT switches of HOUSE that
    have no deck parameter -- the mode sign convention ``INCOHSIGN`` (D-INC-03) and the merging of
    coincident plan positions ``INCOHMERGE`` (D-INC-09) -- go to ``HOUSE.opt`` in the model directory
    (EDUOPT lines, as COMBIN.opt), which HOUSE reads (:func:`sassi.core.incoherency.read_options`)."""
    from ..core.incoherency import OPTION_FILE, options_text
    path = mdir / OPTION_FILE
    path.write_text(options_text(eduopt(model, "INCOHSIGN"), eduopt(model, "INCOHMERGE")), encoding="utf-8")
    res.notes.append(f"{OPTION_FILE} written (HOUSE incoherency switches EDUOPT INCOHSIGN {eduopt(model, 'INCOHSIGN')}, "
                     f"INCOHMERGE {eduopt(model, 'INCOHMERGE')}; D-INC-03, D-INC-09)")
    return path


def write_pin_file(model, mdir: Path, res: AfwriteResult, view=None) -> Optional[Path]:
    """Non-linear soil SSI (HOUSEX <nlssi> = 1, requirements 4.4 item 7): write ``<model>.pin`` from the
    PIN / PINGRP / PINMAT records (:func:`sassi.prep.commands.nlsoil_cmds.pin_from_model`).  Without
    PINGRP records a hand-edited ``.pin`` in the model directory is left as it is (HOUSE reads it).
    An existing ``.pin`` with another content is kept as ``.pin.bak``.  The iteration state
    (``<model>.liq``) is never changed here: NLSSIRESET starts a new analysis.  ``view`` is the CHECK
    model view (ETYPE 0 resolved as in CHECK, D-HOU-01)."""
    from .commands.nlsoil_cmds import has_pin_data, pin_text
    path = mdir / f"{model.name}.pin"
    if not has_pin_data(model):
        if not path.exists():
            res.warnings.append(f"HOUSE: Non-Linear SSI is on but {path.name} does not exist and no PINGRP is "
                                "defined (edit the .pin with the HOUSE dialog Input Data, or use PIN/PINGRP/PINMAT)")
        return None
    text, probs, notes = pin_text(model, view=view)
    if text is None:
        res.warnings.append(f"HOUSE: {path.name} not written -- " + "; ".join(t for _, t in probs))
        return None
    if path.exists():
        old = path.read_text(encoding="utf-8", errors="replace")
        if old == text:
            return path
        bak = path.with_name(path.name + ".bak")
        os.replace(path, bak)
        res.notes.append(f"{path.name}: previous content kept as {bak.name}")
    path.write_text(text, encoding="utf-8")
    res.notes.append(f"{path.name} written (non-linear soil SSI input, {text.count(chr(10))} lines)")
    res.notes.extend(f"{path.name}: {t}" for t in notes)
    return path


SIM_BEGIN = "* ---- Simulation commands written by AFWRITE (D-AFW-04): module runs in run order"
SIM_END = "* ---- end of simulation commands"


def write_simulation_commands(model, report: CheckReport, res: AfwriteResult, location: str = "sim") -> Path:
    """Options > Write "Simulation Commands" (D-AFW-04): RUN commands of the AOPT-enabled modules whose
    decks were written, in run order, to ``<model>-Sim.pre`` (location ``'sim'``, default) or as a
    block at the end of ``<model>.pre`` (location ``'pre'``; an earlier block is replaced)."""
    mdir = Path(model.path)
    runs = [f"RUN{mod}" for mod in report.modules
            if mod in res.written or (mod == "COMBIN" and not report.blocked(mod))]
    block = [SIM_BEGIN] + runs + [SIM_END]
    if location.lower() == "pre":
        p = mdir / f"{model.name}.pre"
        old = p.read_text(encoding="utf-8").splitlines() if p.exists() else []
        if SIM_BEGIN in old:
            i = old.index(SIM_BEGIN)
            j = old.index(SIM_END, i) + 1 if SIM_END in old[i:] else len(old)
            old = old[:i] + old[j:]
        p.write_text("\n".join(old + block) + "\n", encoding="utf-8")
        return p
    p = mdir / f"{model.name}-Sim.pre"
    p.write_text("\n".join(block) + "\n", encoding="utf-8")
    return p


__all__ = ["afwrite", "AfwriteResult", "AfwriteError", "DeckBuildError", "DeckBuilder", "write_deck",
           "DECK_MODULES"]
