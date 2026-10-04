"""RELDISP: relative displacements from complex interpolated transfer functions.

Physics (requirements section 4.9, D-RDP-01 ... D-RDP-04, spec 05d section 2)
---------------------------------------------------------------------------
The relative displacement between an output node and a reference location, in one DOF d,
is computed analytically in the frequency domain from the complex interpolated TFs
(``.TFI``) that MOTION wrote with ``<cplx> = 1``::

    H_rel(f) = H_node,d(f) - H_ref,d(f)                       (on the Fourier grid)

* seismic analysis: FILE8 TFs are total-motion ratios ``U/U_cp``, so the relative
  displacement spectrum is ``D(f) = H_rel(f) U_g(f)`` with the control *displacement*
  spectrum ``U_g = -g A/w^2`` (A = FFT of the control acceleration in g, the f = 0 term zero;
  D-CNV-06);
* vibration analysis: the TFs are displacements per unit load factor and
  ``D(f) = H_rel(f) F(f)`` with F the FFT of the reference load history;

and ``d(t) = irfft(D)``.  This avoids the drift of double-integrating accelerations (MOTION
baseline correction) and is the method the manual recommends for relative displacements.

Reference: the RELFILE ``.TFI`` (its DOF tag, e.g. ``00415TR_X.TFI``, fixes the DOF of the run:
one DOF per run), or the *free field* (RELFILE empty or ``FREEFIELD``, D-RDP-04): a synthetic
reference of amplitude 1 and phase 0 in the input direction, i.e. the rigid-body projection
of the control direction on each DOF (0 for the other directions and rotations).  The
free-field reference holds at *every* Fourier frequency (it is the control motion itself),
whereas MOTION's node TFIs are zero above the SSI cut-off f_N (D-MOT-03); above f_N the
relative displacement therefore contains the (small) ground displacement that the absolute
response omits, exactly ``irfft((H - 1) U_g)`` of VP-34.

Outputs: ``.TFD`` = complex relative-displacement TF per unit control acceleration (length per
g, D-RDP-02; vibration: H_rel), ``.THD`` = relative-displacement history (first line dt), the
listing with the maxima, and ``THD``/``THDR`` frames with RELD <RelDispSAll> or RELDX
<rstframes>.  RELD <RelDisOutput> bit 0 = complex ``.TFD`` (else amplitude only), bit 1 = use
``.TFU`` instead of ``.TFI`` (discouraged; the TFU is interpolated here with option 1).
Node files are found with 5- or 6-digit node numbers (D-FIL-05; see :func:`resolve_node_file`)
and the outputs keep the width of the MOTION files they come from.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from .. import conventions as C
from ..core import interp as I
from ..core import signal as S
from ..io import decks, textfiles
from ..io.container import read_container
from .base import ModuleContext, ModuleError, batch_main
from .motion import (MAX_NFFT, read_control_motion, write_frame, write_history_file, write_tf_file,
                     frame_time, optimizer_numbering_note)

NAME = "RELDISP"

_TAG_RE = re.compile(r"(TR_X|TR_Y|TR_Z|R_XX|R_YY|R_ZZ)(?:\d\d)?\.TF[IU]$", re.IGNORECASE)
_NODE_RE = {e: re.compile(r"^(\d{5,6})(TR_X|TR_Y|TR_Z|R_XX|R_YY|R_ZZ)\.%s$" % e, re.IGNORECASE)
            for e in ("TFI", "TFU")}
TAG_DOF = {v: k for k, v in C.DOF_TAGS.items()}
_DOF_COLS = ("x", "y", "z", "xx", "yy", "zz")


def dof_from_name(name: str) -> Optional[int]:
    """DOF of a TF file from its tag (``00415TR_X.TFI`` -> 1); None if no tag (spec 05d section 2.1)."""
    m = _TAG_RE.search(Path(name).name)
    return TAG_DOF[m.group(1).upper()] if m else None


def resolve_node_file(workdir: Path, node: int, dof: int, ext: str, maxnode: int = 0) -> Optional[Path]:
    """Find the MOTION output file of (node, dof) whatever node-number width MOTION used (D-FIL-05).

    File names carry 5 digits, or 6 when the largest node number exceeds 99,999.  MOTION takes the
    width from max(deck <maxnode>, largest node in FILE8), so a model with a node >= 100000 gets
    6-digit names even for small node numbers while the RELDISP deck may hold the default
    <maxnode> = 0.  The name of the deck's width is tried first, then the 5- and 6-digit names.
    """
    names = [C.nodal_result_name(node, dof, ext, maxnode), C.nodal_result_name(node, dof, ext, 0),
             C.nodal_result_name(node, dof, ext, 100000)]
    for nm in dict.fromkeys(names):
        if (Path(workdir) / nm).exists():
            return Path(workdir) / nm
    return None


def _width_maxnode(path: Path, maxnode: int) -> int:
    """A <maxnode> value that reproduces the node-number width of an existing file name."""
    digits = len(path.name) - len(path.name.lstrip("0123456789"))
    return max(int(maxnode), 100000 if digits >= 6 else 0)


def read_tf_on_grid(path: Path, nK: int, df: float, need_complex: bool = True,
                    tfu: bool = False, option: int = 1, h0: complex = 0.0, mode: str = "seismic") -> np.ndarray:
    """Read a ``.TFI`` (or ``.TFU``) and return it on the rfft grid k = 0 .. nK-1 (zero above f_N).

    A ``.TFI`` must lie on the Fourier grid f_k = k df (checked); a ``.TFU`` (computed
    frequencies) is interpolated with ``option``.
    """
    f, H, cplx = textfiles.read_tf(path)
    if need_complex and not cplx:
        raise ModuleError(f"{path.name} holds amplitudes only: run MOTION with Save Complex TF (<cplx> = 1)")
    H = np.asarray(H, dtype=complex)
    out = np.zeros(nK, dtype=complex)
    if len(f) == 0:
        return out
    if tfu:
        grid = np.arange(nK) * df
        return I.interpolate_tf(f, H, grid, option, h0=h0, mode=mode)
    k = np.rint(f / df).astype(np.int64)
    tol = 1e-6 + 1e-8 * np.abs(f)
    if np.any(np.abs(f - k * df) > tol) or np.any(k < 0) or np.any(np.diff(k) <= 0):
        raise ModuleError(f"{path.name}: frequencies are not on the Fourier grid k*df (df = {df:.9g} Hz)")
    if k[-1] >= nK:
        raise ModuleError(f"{path.name}: frequencies beyond the Nyquist frequency of NFFT")
    out[k] = H
    return out


def solution_direction(workdir: Path, d, lst=None) -> Tuple[int, float]:
    """Control direction ``<cm>`` and angle ``<ang>`` of the solution file the TFs come from.

    The free-field reference and the zero-frequency anchor are the rigid-body projection of the
    control motion on each DOF (D-RDP-04), so they need the direction and the angle of the FILE8 that
    MOTION interpolated (deck ``<file8>``: FILE8, or FILE8X / FILE8Y / FILE8Z of simultaneous cases,
    D-ANL-06).  As in MOTION and STRESS, the values stored by ANALYS in that file win over the deck
    (a deck written after another SITE run, or for a copied FILE8, may carry other values); the deck
    values are used when the file is not in the model directory or cannot be read.
    """
    cm, ang = int(d["cm"]), float(d["ang"])
    name = str(d["file8"] or "FILE8").strip().replace("\\", "/")
    p = Path(name)
    if not p.is_absolute():
        p = Path(workdir) / p
    if not p.is_file():
        return cm, ang
    try:
        meta = read_container(p).meta
    except (OSError, ValueError, KeyError):
        return cm, ang
    fcm, fang = int(meta.get("cm", cm)), float(meta.get("ang", ang))
    if (fcm != cm or abs(fang - ang) > 1e-9) and lst is not None:
        lst.warning(f"control direction/angle of {p.name} (cm {fcm}, ang {fang:g}) differ from the deck "
                    f"(cm {cm}, ang {ang:g}); the {p.name} values are used")
    return fcm, fang


def run(ctx: ModuleContext) -> int:
    """RELDISP module entry point (requirements section 4.9)."""
    lst = ctx.listing
    d = decks.read(ctx.deck_path, NAME)
    if d["title"]:
        lst.write(f" {d['title']}")
    seismic = int(d["type"]) == 0
    nfft, dt = int(d["nft"]), float(d["delt"])
    if not C.is_power_of_two(nfft) or nfft > MAX_NFFT:
        raise ModuleError(f"NFFT = {nfft} must be a power of 2 not above {MAX_NFFT}")
    if dt <= 0:
        raise ModuleError("time step <delt> must be > 0")
    df = 1.0 / (nfft * dt)
    if float(d["df"]) > 0 and abs(float(d["df"]) - df) > 1e-6 * df:
        raise ModuleError(f"deck frequency step {d['df']:.9g} Hz differs from 1/(NFFT dt) = {df:.9g} Hz")
    gravity = float(d["gravity"])
    out_code = int(d["reldisoutput"])
    cplx_tfd = bool(out_code & 1)
    use_tfu = bool(out_code & 2)
    allnodes = bool(int(d["reldispsall"]))
    frames = allnodes or bool(int(d["rstframes"]))
    cm, ang = solution_direction(ctx.workdir, d, lst)
    ext = "TFU" if use_tfu else "TFI"
    relfile = (d["relfile"] or "").strip()
    freefield = relfile == "" or relfile.upper() == "FREEFIELD"
    maxnode = int(d["maxnode"])

    lst.section("RELDISP options")
    lst.write(f"   type of analysis : {'seismic' if seismic else 'foundation vibration'}")
    lst.write(f"   reference        : {'free field (unit amplitude, zero phase, input direction)' if freefield else relfile}")
    if freefield and not seismic:
        lst.write("   note: a vibration analysis has no free-field motion; the free-field reference is 0, so the "
                  "results are absolute displacements")
    lst.write(f"   NFFT, dt, df     : {nfft}, {dt:g} s, {df:.9g} Hz")
    lst.write(f"   TF files         : .{ext}" + ("  (WARNING: TFU use is discouraged)" if use_tfu else ""))
    lst.write(f"   .TFD             : {'complex (amplitude + phase)' if cplx_tfd else 'amplitude only'}")
    if use_tfu:
        lst.warning("RELD <RelDisOutput> bit 1: relative displacements from .TFU (computed frequencies) "
                    "interpolated with option 1; MOTION's .TFI is the recommended input")
    if int(d["numfiles"]) != len(d.rows("rdnd")):
        lst.write(f"   note: <RelDispNumFiles> = {d['numfiles']} ignored (RDND list has {len(d.rows('rdnd'))} entries)")

    if not allnodes:                                    # requirements 2.6: HOUSE optimizer numbering
        ref_node = None if freefield else re.match(r"^(\d{5,6})", Path(relfile).name)
        optimizer_numbering_note(ctx, [int(r["node"]) for r in d.rows("rdnd")]
                                 + ([int(ref_node.group(1))] if ref_node else []), lst,
                                 "RDND nodes and the RELFILE reference")

    # ---------------------------------------------------------------- control motion
    thpath = Path(d["thfile"]) if d["thfile"] else None
    if thpath is None:
        raise ModuleError("Error 73: no time-history file (THFILE) given")
    if not thpath.is_absolute():
        thpath = ctx.workdir / thpath
    cmot = read_control_motion(thpath, int(d["fopt"]), int(d["rec1"]), int(d["rec2"]), dt, nfft,
                               float(d["mult"]), float(d["max"]))
    nK = nfft // 2 + 1
    f = S.fourier_grid(nfft, dt)
    lst.write(f"   control motion   : {cmot.path.name}, {len(cmot.acc)} records, peak "
              f"{np.max(np.abs(cmot.acc)):.6g} {'g' if seismic else '(load factor)'}")

    # ---------------------------------------------------------------- reference
    ref_dof = None
    H_ref_file = None
    if not freefield:
        rp = Path(relfile) if Path(relfile).is_absolute() else ctx.workdir / relfile
        if not rp.exists():
            raise ModuleError(f"reference TF file {relfile} not found: run MOTION with Save Complex TF for it")
        ref_dof = dof_from_name(rp.name)
        if ref_dof is None:
            lst.warning(f"reference file {rp.name} carries no DOF tag (TR_X ... R_ZZ): it is used for every "
                        "requested DOF")
        H_ref_file = read_tf_on_grid(rp, nK, df, tfu=use_tfu, mode="seismic" if seismic else "vibration",
                                     h0=I.rigid_body_anchor(ref_dof or 1, cm, ang) if seismic else 0.0)

    # ---------------------------------------------------------------- requests
    req: List[Tuple[int, int]] = []
    if allnodes:
        rot = bool(int(d["saverot"]))
        for p in sorted(ctx.workdir.iterdir()):
            m = _NODE_RE[ext].match(p.name)
            if m:
                k = TAG_DOF[m.group(2).upper()]
                if k <= 3 or rot:
                    req.append((int(m.group(1)), k))
        lst.write(f"   all nodes        : {len(req)} node/DOF .{ext} files found"
                  + ("" if rot else " (translations; RELDX <saverot> adds rotations)"))
    else:
        for r in d.rows("rdnd"):
            for k, col in enumerate(_DOF_COLS, start=1):
                if int(r[col]) >= 1:
                    req.append((int(r["node"]), k))
    req = sorted(set(req))
    if ref_dof is not None:
        skip = [q for q in req if q[1] != ref_dof]
        if skip:
            lst.warning(f"one DOF per run (reference DOF {C.DOF_TAGS[ref_dof]}): requests ignored for "
                        + ", ".join(f"{n} {C.DOF_TAGS[k]}" for n, k in skip[:20]) + (" ..." if len(skip) > 20 else ""))
        req = [q for q in req if q[1] == ref_dof]
    if not req:
        raise ModuleError("no relative-displacement request (RDND) for the reference DOF")
    files = {q: resolve_node_file(ctx.workdir, q[0], q[1], ext, maxnode) for q in req}
    missing = [q for q in req if files[q] is None]
    if missing:
        raise ModuleError("Run MOTION with Save Complex TF for nodes " +
                          ", ".join(f"{n} {C.DOF_TAGS[k]}" for n, k in missing[:30]) +
                          (" ..." if len(missing) > 30 else "") + f" (.{ext} files missing)")

    # ---------------------------------------------------------------- convolution kernel
    if seismic:
        kernel_tfd = S.displacement_spectrum(np.ones(nK, dtype=complex), f, gravity)    # -g/w^2 per unit A
        drive = S.displacement_spectrum(cmot.A, f, gravity)                            # U_g = -g A / w^2
    else:
        kernel_tfd = np.ones(nK, dtype=complex)
        drive = cmot.A                                                                  # F_ref
    nout = S.output_length(nfft, dt, float(d["dur"]))
    lst.section("Relative displacements")
    unit = "length" if seismic else "displacement units"
    lst.write(f"   history length {nout} samples ({nout * dt:g} s)")
    lst.write(f"{'node':>8s} {'dof':>6s} {'max |d|':>16s} {'time (s)':>10s}")
    results: Dict[Tuple[int, int], np.ndarray] = {}
    wd = ctx.workdir
    for i, (node, dof) in enumerate(req):
        if ctx.cancelled():
            raise ModuleError("run cancelled")
        path = files[(node, dof)]
        wmax = _width_maxnode(path, maxnode)          # outputs keep the width of the MOTION files
        H_node = read_tf_on_grid(path, nK, df, tfu=use_tfu, mode="seismic" if seismic else "vibration",
                                 h0=I.rigid_body_anchor(dof, cm, ang) if seismic else 0.0)
        if freefield:
            H_ref = np.full(nK, I.rigid_body_anchor(dof, cm, ang) if seismic else 0.0, dtype=complex)
        else:
            H_ref = H_ref_file
        H_rel = H_node - H_ref
        D = H_rel * drive
        dth = np.fft.irfft(S.hermitian_bins(D, nfft), n=nfft)
        tag = C.DOF_TAGS[dof]
        hdr = f"node {node} dof {tag} relative to {'free field' if freefield else Path(relfile).name}"
        write_tf_file(wd / C.nodal_result_name(node, dof, "TFD", wmax), f, H_rel * kernel_tfd, cplx_tfd,
                      header=hdr + (" (length per g)" if seismic else " (per unit load)"))
        write_history_file(wd / C.nodal_result_name(node, dof, "THD", wmax), dth[:nout], dt)
        k = int(np.argmax(np.abs(dth[:nout])))
        lst.write(f"{node:8d} {tag:>6s} {abs(dth[k]):16.6e} {k * dt:10.4f}")
        if frames:
            results[(node, dof)] = dth[:nout]
        ctx.progress((i + 1) / len(req), f"RELDISP {i + 1}/{len(req)}")
    lst.write(f"   units: {unit} (gravity {gravity:g})" if seismic else f"   units: {unit}")

    if frames and results:
        _write_frames(wd, results, dt, nout, bool(int(d["saverot"])), lst)
    return 0


def _write_frames(wd: Path, results: Dict[Tuple[int, int], np.ndarray], dt: float, nout: int, rot: bool,
                  lst) -> None:
    """``THD``/``THDR`` frames: rows ``node dx dy dz`` per time step (D-FIL-03, Table 3.2)."""
    for dofs, sub in (((1, 2, 3), "THD"), ((4, 5, 6), "THDR")):
        if sub == "THDR" and not rot:
            continue
        nodes = sorted({n for (n, k) in results if k in dofs})
        if not nodes:
            continue
        pos = {n: r for r, n in enumerate(nodes)}
        vals = np.zeros((nout, len(nodes), 3))
        for (n, k), v in results.items():
            if k in dofs:
                vals[:, pos[n], dofs.index(k)] = v
        out = wd / sub
        out.mkdir(exist_ok=True)
        for s in range(nout):
            write_frame(out / f"{sub}_{frame_time(s * dt)}_{s + 1:05d}", np.array(nodes), vals[s])
        lst.write(f"   {sub} frames: {nout} files in {sub}/ ({len(nodes)} nodes)")


if __name__ == "__main__":
    raise SystemExit(batch_main(NAME))
