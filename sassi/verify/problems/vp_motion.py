"""Verification problems of the TF interpolation core, MOTION, RELDISP and COMBIN.

VP-25 COMBIN merge, VP-28 interpolation exactness, VP-29 identity convolution, VP-31 baseline
correction, VP-32 SRSS TF (P1), VP-33 phase adjustment (P1), VP-34 relative displacement
(requirements section 6.3; R1 V5; spec 05c Part C; spec 05d section 8).

ANALYS is not needed: the problems synthesise FILE8 containers (schema of
``sassi.io.files``) from closed-form transfer functions of hysteretic shear chains on a moving
base, write the module decks with ``sassi.io.decks`` and run the modules exactly as RUNMOTION /
RUNRELDISP / RUNCOMBIN do (``sassi.modules.base.run_module``).

Reference system: a fixed-base hysteretic shear chain (masses m_i, storey stiffnesses k_i,
complex modulus factor c(beta) of D-CNV-03) driven by base motion u_g.  With the relative
displacements u = U - u_g::

    (K c(beta) - w^2 M) u = w^2 M 1 u_g        ->     H = U/u_g = 1 + u

which is the total-motion transfer function MOTION receives from ANALYS for a model on a rigid
base (H(0) = 1, the rigid-body limit of D-MOT-03).  The 2-DOF chain is exactly the "two-SDOF
five-parameter" form behind interpolation options 0-5 (R1 section 5.1).
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional, Sequence

import numpy as np

from ...conventions import cfactor, nodal_rs_name, nodal_result_name
from ...core import interp as I
from ...core import signal as S
from ...io import decks, textfiles
from ...io.container import read_container, write_container
from ...modules.base import run_module
from .. import VPResult, problem, worse

# ---------------------------------------------------------------------------------------
# Synthesis helpers (also used by the unit tests)
# ---------------------------------------------------------------------------------------
TWO_DOF = dict(masses=(1.5, 1.0), stiffs=(2.0 * (2 * np.pi * 5.0) ** 2, (2 * np.pi * 3.0) ** 2), beta=0.05)
THREE_DOF = dict(masses=(1.2, 1.0, 0.8), stiffs=(3.0e3, 2.2e3, 1.4e3), beta=0.04)


def shear_chain_tf(f, masses: Sequence[float], stiffs: Sequence[float], beta: float) -> np.ndarray:
    """Total-motion TFs ``H = U/u_g`` (nF, n) of a hysteretic shear chain on a moving base.

    ``masses``/``stiffs`` are ordered bottom to top; storey i joins mass i-1 (or the base) to
    mass i.  Damping is hysteretic with the SASSI complex modulus ``K* = K c(beta)``.
    """
    w = 2.0 * np.pi * np.atleast_1d(np.asarray(f, dtype=float))
    m = np.asarray(masses, dtype=float)
    k = np.asarray(stiffs, dtype=float)
    n = len(m)
    K = np.zeros((n, n))
    for i in range(n):
        K[i, i] += k[i]
        if i + 1 < n:
            K[i, i] += k[i + 1]
            K[i, i + 1] -= k[i + 1]
            K[i + 1, i] -= k[i + 1]
    Ks = K * cfactor(beta)
    A = Ks[None, :, :] - (w ** 2)[:, None, None] * np.diag(m)[None, :, :]
    rhs = (w ** 2)[:, None] * m[None, :]
    u = np.linalg.solve(A, rhs[..., None])[..., 0]
    return 1.0 + u


def write_file8(path: Path, fnum: Iterable[int], df: float, eq_node, eq_dof, H: np.ndarray, type: int = 0,
                cm: int = 0, ang: float = 0.0, nfft: int = 4096, delt: float = 0.005, case: int = 0,
                model_hash: str = "VP-synthetic") -> Path:
    """Write a FILE8 container exactly as the FILE8 schema of ``sassi.io.files`` requires."""
    fnum = np.asarray(list(fnum), dtype=np.int64)
    return write_container(path, "FILE8",
                           {"fnum": fnum, "freq": fnum * df, "eq_node": np.asarray(eq_node, dtype=np.int64),
                            "eq_dof": np.asarray(eq_dof, dtype=np.int64), "H": np.asarray(H, dtype=complex)},
                           meta={"df": df, "type": type, "case": case, "ang": ang, "cm": cm, "nfft": nfft,
                                 "delt": delt, "model_hash": model_hash}, module="ANALYS")


def synthetic_motion(n: int, dt: float, seed: int = 11, fmax: float = 15.0, f_hp: float = 0.25) -> np.ndarray:
    """Deterministic earthquake-like acceleration record (peak 1).

    White noise shaped by a Kanai-Tajimi spectrum (f_g = 2.5 Hz, zeta_g = 0.6), high-passed
    below ``f_hp`` and tapered to zero between ``fmax`` and 1.2 ``fmax``, times a build-up /
    strong-motion / decay envelope.
    """
    rng = np.random.default_rng(seed)
    F = np.fft.rfft(rng.standard_normal(n))
    f = np.fft.rfftfreq(n, dt)
    r = f / 2.5
    kt = np.sqrt((1 + 4 * 0.36 * r ** 2) / ((1 - r ** 2) ** 2 + 4 * 0.36 * r ** 2))
    hp = (f / f_hp) ** 4 / (1 + (f / f_hp) ** 4)
    taper = np.clip((1.2 * fmax - f) / (0.2 * fmax), 0.0, 1.0)
    a = np.fft.irfft(F * kt * hp * taper, n)
    t = np.arange(n) * dt
    T = n * dt
    env = np.where(t < 0.15 * T, (t / (0.15 * T)) ** 2, np.where(t < 0.55 * T, 1.0, np.exp(-(t - 0.55 * T) / (0.12 * T))))
    a = a * env
    return a / np.max(np.abs(a))


def write_motion_file(path: Path, acc: np.ndarray, dt: float) -> Path:
    """Control-motion file in MOTION format fopt 0 (dt on the first line), 17 digits."""
    lines = [f"{dt:.16g}"] + [f"{v:.16e}" for v in acc]
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return Path(path)


def motion_deck(nout: Sequence[Sequence[int]] = (), damp: Sequence[float] = (), **params) -> decks.Deck:
    """MOTION deck with ``params`` set; ``nout`` rows (node, dir, c1..c6); damping list."""
    d = decks.new("MOTION")
    d["mult"], d["max"] = 1.0, 0.0
    for k, v in params.items():
        if k not in d.params:
            raise KeyError(f"unknown MOTION parameter {k}")
        d[k] = v
    if "df" not in params and d["nft"] * d["delt"] > 0:
        d["df"] = 1.0 / (d["nft"] * d["delt"])
    for r in nout:
        d.table("nout").append(list(r))
    for z in damp:
        d.table("damp").append([z])
    return d


def run_deck(workdir: Path, model: str, deck: decks.Deck) -> int:
    """Write ``deck`` as ``<model>.<ext>`` in ``workdir`` and run its module."""
    path = decks.deck_path(workdir, model, deck.module)
    decks.write(path, deck)
    return run_module(deck.module, model, workdir)


def listing_text(workdir: Path, model: str, module: str) -> str:
    return (Path(workdir) / f"{model}_{module}.out").read_text(encoding="utf-8")


def _relmax(a, b) -> float:
    a, b = np.asarray(a), np.asarray(b)
    den = np.max(np.abs(b))
    return float(np.max(np.abs(a - b)) / (den if den > 0 else 1.0))


def _pointwise_rel(a, b) -> float:
    a, b = np.asarray(a), np.asarray(b)
    return float(np.max(np.abs(a - b) / np.abs(b)))


# Common grid: dt = 0.005 s, NFFT = 4096 -> df = 0.048828125 Hz, Nyquist 100 Hz
NFFT, DT = 4096, 0.005
DF = 1.0 / (NFFT * DT)
#: typical SSI frequency set: 60-70 frequency numbers, log spaced, cut-off 25 Hz (fnum 512)
SSI_FNUM = np.unique(np.round(np.geomspace(4, 512, 64)).astype(np.int64))


# =======================================================================================
# VP-28 interpolation exactness
# =======================================================================================
@problem("VP-28", "TF interpolation exactness (options 0-6) and ISRS from interpolated TFs", tier="P0",
         modules=["MOTION", "STRESS"], source="R1 V5; spec 05c C8/C9/C14; requirements 6.3")
def vp28(workdir: Path) -> VPResult:
    r = VPResult()
    workdir = Path(workdir)
    fk = S.fourier_grid(NFFT, DT)
    fs = SSI_FNUM * DF
    ins = (fk >= fs[0]) & (fk <= fs[-1])
    h0 = np.ones(2)
    # (a) exact 2-DOF hysteretic TF -> options 0-5 reproduce it between the computed points
    H2 = shear_chain_tf(fs, **TWO_DOF)
    ex2 = shear_chain_tf(fk[ins], **TWO_DOF)
    for opt in range(6):
        Hi = I.interpolate_tf(fs, H2, fk, opt, h0=h0)
        r.check(f"(a) 2-DOF hysteretic TF, option {opt}: max pointwise rel. error on [f1, fN]",
                _pointwise_rel(Hi[ins], ex2), 0.0, atol=1e-10)
    # (a') hysteretic SDOF: every window is rank deficient (rank 4) but still exact (R1 V5)
    w0 = 2 * np.pi * 4.0
    ks = w0 ** 2 * cfactor(0.05)
    sdof = lambda f: ks / (ks - (2 * np.pi * np.asarray(f)) ** 2)
    fit = I.fit_windows(fs, sdof(fs)[:, None])
    r.require("(a') hysteretic SDOF: all 5x5 window systems rank deficient (rank 4)",
              bool(np.all(fit["rank"] == 4)), note=f"ranks found {sorted(set(fit['rank'].ravel().tolist()))}")
    for opt in range(6):
        Hi = I.interpolate_tf(fs, sdof(fs), fk, opt, h0=1.0)
        r.check(f"(a') hysteretic SDOF, option {opt}: max pointwise rel. error",
                _pointwise_rel(Hi[ins], sdof(fk[ins])), 0.0, atol=1e-10)
    # (a'') lightly damped low-frequency SDOFs on uniform coarse SSI grids: the genuine resonance pole
    # fails the literal D-MOT-01 ratio test min|D| < 1e-3 max|D| (wide window, small x0), but it is a
    # damped (physical) pole, so the refined guard keeps the exact rational form (review defect fix)
    n_trig = 0
    for f0, beta, step in ((0.5, 0.02, 8), (0.2, 0.05, 4), (0.3, 0.02, 4), (0.6, 0.01, 8), (1.0, 0.01, 8)):
        ksu = (2 * np.pi * f0) ** 2 * cfactor(beta)
        sdu = lambda f, ksu=ksu: ksu / (ksu - (2 * np.pi * np.asarray(f)) ** 2)
        fu = np.arange(4, 513, step) * DF
        iu = (fk >= fu[0]) & (fk <= fu[-1])
        fitu = I.fit_windows(fu, sdu(fu)[:, None])
        n_trig += int(fitu["trigger"].sum())
        worst_u = worse(*(_pointwise_rel(I.interpolate_tf(fu, sdu(fu), fk, opt, h0=1.0)[iu], sdu(fk[iu]))
                          for opt in range(6)))      # NaN-propagating (final audit; built-in max drops NaN)
        r.check(f"(a'') hysteretic SDOF f0 = {f0} Hz, beta = {beta}, uniform spacing {step} df: max pointwise "
                "rel. error, options 0-5", worst_u, 0.0, atol=1e-10)
    r.notes.append(f"(a''): {n_trig} windows trip the literal D-MOT-01 ratio test on genuine damped poles; with the "
                   "literal rule (cubic there) these five cases had max pointwise errors 1.0 to 10 (100-1000 %)")
    # (b) option 6 exact on cubic data (with the f = 0 anchor p(0))
    cub = lambda f: (0.3 + 0.1j) + (0.5 - 0.2j) * f + (-0.04 + 0.01j) * f ** 2 + (0.002 + 0.0005j) * f ** 3
    Hi = I.interpolate_tf(fs, cub(fs), fk, 6, h0=cub(0.0))
    r.check("(b) option 6 on complex cubic data: max rel. error", _relmax(Hi[ins], cub(fk[ins])), 0.0, atol=1e-10)
    # (b') option 6 is O(df^4) otherwise: uniform SSI spacings h = 0.05, 0.025, 0.0125 Hz (well below the
    # 0.26 Hz half-power band of the first mode, i.e. in the asymptotic regime), 2-DOF TF up to 12 Hz,
    # errors measured on a 1 mHz grid
    fo = np.arange(0.2, 12.0, 0.001)
    errs = []
    hs = (0.05, 0.025, 0.0125)
    for h in hs:
        fu = np.arange(h, 12.0 + 1e-9, h)
        iu = (fo >= fu[0]) & (fo <= fu[-1])
        Hu = I.interpolate_tf(fu, shear_chain_tf(fu, **TWO_DOF), fo, 6, h0=h0)
        errs.append(_relmax(Hu[iu], shear_chain_tf(fo[iu], **TWO_DOF)))
    order = np.log2(errs[-2] / errs[-1])
    r.check("(b') option 6 convergence order (h = 0.025 -> 0.0125 Hz), O(df^4) expected", order, 4.0, atol=0.5,
            note="max rel. errors " + ", ".join(f"{e:.3e}" for e in errs) + " for h = " + ", ".join(map(str, hs)) + " Hz")
    # (c) TFI = TFU at the computed frequencies, every option and smoothing S
    kf = SSI_FNUM
    worst = 0.0
    for opt in range(7):
        for Sm in (0.0, 10.0, 1000.0):
            Hi = I.interpolate_tf(fs, H2, fk, opt, smooth=Sm, h0=h0)
            worst = worse(worst, float(np.max(np.abs(Hi[kf] - H2))))
    r.check("(c) |TFI - TFU| at computed frequencies, options 0-6, S = 0/10/1000 (in memory)", worst, 0.0, atol=0.0)
    # (c') the same through the MOTION module and its text files (option 0, S = 10)
    eq_node, eq_dof = [1, 2, 3], [1, 1, 1]
    H8 = np.column_stack([np.ones(len(fs)), H2])
    write_file8(workdir / "FILE8", SSI_FNUM, DF, eq_node, eq_dof, H8, nfft=NFFT, delt=DT)
    n = 2400
    acc = 0.3 * synthetic_motion(n, DT)
    write_motion_file(workdir / "eq.acc", acc, DT)
    deck = motion_deck(nout=[(3, 1, 1, 1, 0, 0, 1, 1)], damp=[0.05], thfile="eq.acc", nft=NFFT, delt=DT, interp=0,
                       smo=10.0, cplx=1)
    rc = run_deck(workdir, "vp28c", deck)
    r.require("(c') MOTION run (option 0, S = 10) succeeded", rc == 0)
    if rc == 0:
        fu, Hu, _ = textfiles.read_tf(workdir / nodal_result_name(3, 1, "TFU"))
        fi, Hti, _ = textfiles.read_tf(workdir / nodal_result_name(3, 1, "TFI"))
        r.check("(c') .TFI rows at the computed frequencies equal the .TFU rows (files)",
                float(np.max(np.abs(Hti[kf] - Hu))), 0.0, atol=0.0)
    # (d) ISRS from interpolated vs dense TF (3-DOF chain: not of the 2-DOF form)
    fd = np.arange(1, 513)                                     # every Fourier bin up to 25 Hz
    Hd = np.column_stack([np.ones(len(fd)), shear_chain_tf(fd * DF, **THREE_DOF)])
    eqn, eqd = [1, 2, 3, 4], [1, 1, 1, 1]
    rs_dense = _isrs_run(workdir / "dense", fd, Hd, eqn, eqd, acc, option=1)
    dense_grid = np.unique(np.round(np.geomspace(4, 512, 220)).astype(np.int64))
    for opt in range(7):
        fn = SSI_FNUM if opt < 6 else dense_grid
        Hc = np.column_stack([np.ones(len(fn)), shear_chain_tf(fn * DF, **THREE_DOF)])
        rs = _isrs_run(workdir / f"opt{opt}", fn, Hc, eqn, eqd, acc, option=opt)
        if rs is None or rs_dense is None:
            r.require(f"(d) option {opt}: MOTION runs succeeded", False)
            continue
        err = float(np.max(np.abs(rs - rs_dense) / rs_dense))
        r.check(f"(d) option {opt}: ISRS (top node, 5 %) max rel. difference interpolated ({len(fn)} SSI freq.) vs dense",
                err, 0.0, atol=0.02)
    r.notes.append("SSI frequency set: %d log-spaced frequency numbers 4..512 (0.195-25 Hz), df = %.9g Hz; option 6 "
                   "ISRS uses %d frequencies (manual: the spline needs a denser grid)" % (len(SSI_FNUM), DF, len(dense_grid)))
    return r


def _isrs_run(wd: Path, fnum, H, eq_node, eq_dof, acc, option: int) -> Optional[np.ndarray]:
    """Run MOTION on a synthetic FILE8 and return the 5 % RS of the top node (node 4)."""
    wd.mkdir(parents=True, exist_ok=True)
    write_file8(wd / "FILE8", fnum, DF, eq_node, eq_dof, H, nfft=NFFT, delt=DT)
    write_motion_file(wd / "eq.acc", acc, DT)
    deck = motion_deck(nout=[(4, 1, 0, 0, 0, 0, 1, 1)], damp=[0.05], thfile="eq.acc", nft=NFFT, delt=DT,
                       interp=option, freq1=0.1, freq2=100.0, fstep=301)
    if run_deck(wd, "isrs", deck) != 0:
        return None
    return textfiles.read_xy(wd / nodal_rs_name(4, 1, 1))[:, 1]


# =======================================================================================
# VP-29 identity convolution
# =======================================================================================
@problem("VP-29", "Identity convolution (H = 1 up to Nyquist)", tier="P0", modules=["MOTION"],
         source="spec 05c C7; requirements 6.3")
def vp29(workdir: Path) -> VPResult:
    r = VPResult()
    workdir = Path(workdir)
    nfft, dt = 2048, 0.005
    df = 1.0 / (nfft * dt)
    acc = 0.25 * synthetic_motion(1500, dt, seed=3, fmax=40.0)
    a_pad = S.pad_record(acc, nfft)
    # core: irfft(1 * rfft(a)) = a
    r.check("core convolve(H = 1): max |r - a| / max|a|", _relmax(S.convolve(np.ones(nfft // 2 + 1), np.fft.rfft(a_pad)),
                                                                 a_pad), 0.0, atol=1e-12)
    # module: FILE8 with H = 1 at every frequency number 1..NFFT/2 (Nyquist included), H(0) = 1 anchor
    fnum = np.arange(1, nfft // 2 + 1)
    write_file8(workdir / "FILE8", fnum, df, [1], [1], np.ones((len(fnum), 1)), nfft=nfft, delt=dt)
    write_motion_file(workdir / "eq.acc", acc, dt)
    deck = motion_deck(nout=[(1, 1, 0, 1, 0, 0, 0, 1)], thfile="eq.acc", nft=nfft, delt=dt, interp=1, dur=0.0)
    rc = run_deck(workdir, "vp29", deck)
    r.require("MOTION run succeeded", rc == 0)
    if rc == 0:
        out, dto = textfiles.read_history(workdir / nodal_result_name(1, 1, "ACC"))
        r.require("output length = NFFT (dur = 0)", len(out) == nfft)
        r.check("MOTION .ACC vs zero-padded input: max |r - a| / max|a|", _relmax(out, a_pad), 0.0, atol=1e-12)
    return r


# =======================================================================================
# VP-31 baseline correction
# =======================================================================================
@problem("VP-31", "Baseline correction (Hudson-Housner) removes a constant acceleration offset", tier="P0",
         modules=["MOTION"], source="spec 05c C11; requirements 6.3, D-MOT-08")
def vp31(workdir: Path) -> VPResult:
    r = VPResult()
    workdir = Path(workdir)
    nfft, dt, g = 2048, 0.01, 9.81
    df = 1.0 / (nfft * dt)
    a0 = 0.01                                             # offset, g
    T = (nfft - 1) * dt
    # (a) through MOTION: identity FILE8, a record that is the pure offset over the whole Fourier period
    fnum = np.arange(1, nfft // 2 + 1)
    write_file8(workdir / "FILE8", fnum, df, [1], [1], np.ones((len(fnum), 1)), nfft=nfft, delt=dt)
    write_motion_file(workdir / "offset.acc", np.full(nfft, a0), dt)
    deck = motion_deck(nout=[(1, 1, 0, 1, 0, 0, 0, 1)], thfile="offset.acc", nft=nfft, delt=dt, interp=1, bl=1,
                       f1213=1, gravity=g)
    rc = run_deck(workdir, "vp31", deck)
    r.require("MOTION run with baseline correction and FILE13 succeeded", rc == 0)
    if rc == 0:
        F13 = textfiles.read_xy(workdir / "FILE13")
        t, a_c, v_c, d_c = F13.T
        v_unc_end = a0 * g * T
        d_unc_end = 0.5 * a0 * g * T ** 2
        B = np.stack([t, 0.5 * t * t], axis=1)
        c, *_ = np.linalg.lstsq(B, v_c, rcond=None)
        r.check("(a) corrected acceleration max |a_c| / a0", np.max(np.abs(a_c)) / a0, 0.0, atol=1e-6)
        r.check("(a) velocity trend of the corrected record |c1 T + c2 T^2/2| / (a0 g T)",
                abs(c[0] * T + 0.5 * c[1] * T * T) / v_unc_end, 0.0, atol=1e-6)
        r.check("(a) final displacement |d_c(T)| / (a0 g T^2/2)", abs(d_c[-1]) / d_unc_end, 0.0, atol=1e-6)
    # (b) core: an earthquake record plus the offset -- the offset-induced drift is removed exactly
    aq = 0.2 * synthetic_motion(1200, dt, seed=5)
    aq = S.pad_record(aq, nfft)
    b1 = S.baseline_correction(aq + a0, dt, scale=g)
    b0 = S.baseline_correction(aq, dt, scale=g)
    v_raw, d_raw = S.integrate_trapz(aq + a0, dt)
    r.check("(b) offset removed: |c1[a+a0] - c1[a] - a0| / a0", abs(b1.c1 - b0.c1 - a0) / a0, 0.0, atol=1e-6)
    r.check("(b) offset-induced final displacement |d_c[a+a0](T) - d_c[a](T)| / (a0 g T^2/2)",
            abs(b1.final_disp - b0.final_disp) / (0.5 * a0 * g * T ** 2), 0.0, atol=1e-6)
    tt = np.arange(nfft) * dt
    cc, *_ = np.linalg.lstsq(np.stack([tt, 0.5 * tt * tt], axis=1), b1.vel, rcond=None)
    r.check("(b) velocity trend of the corrected record |c1 T + c2 T^2/2| / (a0 g T)",
            abs(cc[0] * T + 0.5 * cc[1] * T * T) / (a0 * g * T), 0.0, atol=1e-6)
    r.notes.append(f"uncorrected final displacement of (b): {d_raw[-1] * g:.4g} (offset alone {0.5 * a0 * g * T * T:.4g}); "
                   f"corrected: {b1.final_disp:.4g}")
    return r


# =======================================================================================
# VP-34 relative displacement (RELDISP)
# =======================================================================================
def _vp34_setup(workdir: Path, nfft: int, dt: float, g: float):
    """2-DOF chain on a rigid base: nodes 1 (base / control point), 2, 3; MOTION with complex TFI."""
    df = 1.0 / (nfft * dt)
    fnum = np.unique(np.concatenate([np.arange(1, 8), np.round(np.geomspace(8, 512, 64))]).astype(np.int64))
    chain = dict(masses=(1.5, 1.0), stiffs=(2.0 * (2 * np.pi * 4.0) ** 2, (2 * np.pi * 2.2) ** 2), beta=0.05)
    H = np.column_stack([np.ones(len(fnum)), shear_chain_tf(fnum * df, **chain)])
    write_file8(workdir / "FILE8", fnum, df, [1, 2, 3], [1, 1, 1], H, nfft=nfft, delt=dt)
    acc = 0.3 * synthetic_motion(int(0.55 * nfft), dt, seed=21, fmax=12.0)
    write_motion_file(workdir / "eq.acc", acc, dt)
    deck = motion_deck(nout=[(n, 1, 1, 1, 0, 0, 0, 1) for n in (1, 2, 3)], thfile="eq.acc", nft=nfft, delt=dt,
                       interp=0, cplx=1, gravity=g)
    return run_deck(workdir, "vp34", deck), acc


def reldisp_deck(relfile: str, rdnd: Sequence[Sequence[int]], **params) -> decks.Deck:
    d = decks.new("RELDISP")
    d["mult"], d["max"] = 1.0, 0.0
    d["relfile"] = relfile
    for k, v in params.items():
        if k not in d.params:
            raise KeyError(f"unknown RELDISP parameter {k}")
        d[k] = v
    if "df" not in params and d["nft"] * d["delt"] > 0:
        d["df"] = 1.0 / (d["nft"] * d["delt"])
    for row in rdnd:
        d.table("rdnd").append(list(row))
    d["numfiles"] = len(rdnd)
    return d


@problem("VP-34", "Relative displacement (RELDISP)", tier="P0", modules=["RELDISP", "MOTION"],
         source="spec 05d section 8 item 5; requirements 6.3, D-RDP-01, D-CNV-06")
def vp34(workdir: Path) -> VPResult:
    r = VPResult()
    workdir = Path(workdir)
    nfft, dt, g = 8192, 0.0025, 9.81
    rc, acc = _vp34_setup(workdir, nfft, dt, g)
    r.require("MOTION run (complex TFI for nodes 1-3) succeeded", rc == 0)
    if rc != 0:
        return r
    common = dict(thfile="eq.acc", nft=nfft, delt=dt, gravity=g)
    # free-field unit reference
    rc = run_deck(workdir, "ff", reldisp_deck("FREEFIELD", [(3, 1, 0, 0, 0, 0, 0), (2, 1, 0, 0, 0, 0, 0)], **common))
    r.require("RELDISP run with free-field reference succeeded", rc == 0)
    if rc != 0:
        return r
    d_ff, _ = textfiles.read_history(workdir / nodal_result_name(3, 1, "THD"))
    ft, Htd, _ = textfiles.read_tf(workdir / nodal_result_name(3, 1, "TFD"))
    # (1) reference = the node itself -> zero
    rc = run_deck(workdir, "self", reldisp_deck(nodal_result_name(3, 1, "TFI"), [(3, 1, 0, 0, 0, 0, 0)], **common))
    r.require("RELDISP run with reference = node 3 itself succeeded", rc == 0)
    if rc == 0:
        d_self, _ = textfiles.read_history(workdir / nodal_result_name(3, 1, "THD"))
        r.check("(1) reference = node itself: max|d| / max|d_free-field|", np.max(np.abs(d_self)) / np.max(np.abs(d_ff)),
                0.0, atol=1e-10)
    # (2) free-field reference equals irfft((H - 1) U_g), U_g = -g A / w^2 (H from MOTION's .TFI)
    f = S.fourier_grid(nfft, dt)
    nK = len(f)
    fi, Hi, cplx = textfiles.read_tf(workdir / nodal_result_name(3, 1, "TFI"))
    Hg = np.zeros(nK, dtype=complex)
    Hg[np.rint(fi / (f[1] - f[0])).astype(int)] = Hi
    A = np.fft.rfft(S.pad_record(acc, nfft))
    Ug = np.zeros(nK, dtype=complex)
    Ug[1:] = -g * A[1:] / (2 * np.pi * f[1:]) ** 2
    d_ref = np.fft.irfft((Hg - 1.0) * Ug, n=nfft)
    r.require("(2) node .TFI is complex", bool(cplx))
    r.check("(2) free-field reference: max|d - irfft((H-1) U_g)| / max|d|", _relmax(d_ff, d_ref), 0.0, atol=1e-10)
    # (3) double integration of the acceleration difference (MOTION .ACC minus the control motion)
    a_node, _ = textfiles.read_history(workdir / nodal_result_name(3, 1, "ACC"))
    a_rel = g * (a_node - S.pad_record(acc, nfft))
    v, _ = S.integrate_trapz(a_rel, dt)
    v = v - v.mean()                     # periodic solution: zero-mean velocity (f = 0 term of D is 0)
    dd = np.concatenate([[0.0], np.cumsum(0.5 * (v[1:] + v[:-1]) * dt)])
    dd = dd - dd.mean()                  # zero-mean displacement (f = 0 term set to 0)
    r.check("(3) free-field reference vs double integration of a_node - a_g: max|diff| / max|d|",
            _relmax(dd, d_ff), 0.0, atol=1e-3)
    # (4) .TFD of the free-field run = H_rel * (-g/w^2) (length per g)
    k = np.rint(ft / (f[1] - f[0])).astype(int)
    tfd_ref = np.zeros(len(k), dtype=complex)
    nz = k > 0
    tfd_ref[nz] = (Hg[k[nz]] - 1.0) * (-g / (2 * np.pi * f[k[nz]]) ** 2)
    r.check("(4) .TFD = (H - 1)(-g/w^2): max rel. difference", _relmax(Htd, tfd_ref), 0.0, atol=1e-10)
    r.notes.append(f"NFFT {nfft}, dt {dt} s, g {g}; max relative displacement of node 3: {np.max(np.abs(d_ff)):.5g}")
    return r


# =======================================================================================
# VP-25 COMBIN
# =======================================================================================
@problem("VP-25", "COMBIN merge of odd and even frequency sets", tier="P0", modules=["COMBIN"],
         source="spec 02 Part 6; requirements 6.3, D-CMB-01")
def vp25(workdir: Path) -> VPResult:
    r = VPResult()
    workdir = Path(workdir)
    fnum = np.arange(1, 121)
    eq_node = np.repeat([1, 2, 3], 3)
    eq_dof = np.tile([1, 2, 3], 3)
    H = np.empty((len(fnum), 9), dtype=complex)
    Hc = shear_chain_tf(fnum * DF, **TWO_DOF)
    H[:, 0], H[:, 3], H[:, 6] = 1.0, Hc[:, 0], Hc[:, 1]
    H[:, [1, 2, 4, 5, 7, 8]] = 1e-3 * (1 + 0.5j) * Hc[:, [0, 1, 0, 1, 0, 1]] * np.arange(1, 7)
    full = workdir / "full"
    full.mkdir(exist_ok=True)
    write_file8(full / "FILE8", fnum, DF, eq_node, eq_dof, H)
    odd, even = fnum % 2 == 1, fnum % 2 == 0
    write_file8(workdir / "FILE81", fnum[odd], DF, eq_node, eq_dof, H[odd])
    write_file8(workdir / "FILE82", fnum[even], DF, eq_node, eq_dof, H[even])
    rc = run_module("COMBIN", "vp25", workdir)
    r.require("COMBIN run succeeded", rc == 0)
    if rc == 0:
        a = read_container(workdir / "FILE8", "FILE8")
        b = read_container(full / "FILE8", "FILE8")
        r.require("merged frequency numbers equal the full run", np.array_equal(a["fnum"], b["fnum"]))
        r.require("merged DOF map equals the full run",
                  np.array_equal(a["eq_node"], b["eq_node"]) and np.array_equal(a["eq_dof"], b["eq_dof"]))
        r.check("merged frequencies: max |f - f_full|", float(np.max(np.abs(a["freq"] - b["freq"]))), 0.0, atol=0.0)
        r.check("merged H: max |H - H_full|", float(np.max(np.abs(a["H"] - b["H"]))), 0.0, atol=0.0)
        r.require("merged meta df/type equal", a.meta["df"] == b.meta["df"] and a.meta["type"] == b.meta["type"])
    # duplicate frequency -> error
    dup = workdir / "dup"
    dup.mkdir(exist_ok=True)
    sel2 = even | (fnum == 3)
    H2 = H.copy()
    H2[2] *= 2.0                                     # distinguishable FILE82 value at fnum 3
    write_file8(dup / "FILE81", fnum[odd], DF, eq_node, eq_dof, H[odd])
    write_file8(dup / "FILE82", fnum[sel2], DF, eq_node, eq_dof, H2[sel2])
    rc = run_module("COMBIN", "vp25d", dup)
    txt = listing_text(dup, "vp25d", "COMBIN")
    r.require("duplicate frequency number -> COMBIN error", rc != 0 and "present in both" in txt)
    r.require("no FILE8 written after the error", not (dup / "FILE8").exists())
    # EDUOPT,COMBINDUP,PREFER82 keeps the FILE82 value
    (dup / "COMBIN.opt").write_text("EDUOPT,COMBINDUP,PREFER82\n", encoding="utf-8")
    rc = run_module("COMBIN", "vp25p", dup)
    r.require("PREFER82: COMBIN run succeeded", rc == 0)
    if rc == 0:
        a = read_container(dup / "FILE8", "FILE8")
        r.require("PREFER82: frequency numbers equal the full run", np.array_equal(a["fnum"], fnum))
        r.check("PREFER82: duplicate row taken from FILE82", float(np.max(np.abs(a["H"][2] - H2[2]))), 0.0, atol=0.0)
    return r


# =======================================================================================
# VP-32 SRSS TF (P1)
# =======================================================================================
@problem("VP-32", "SRSS TF with one modal FILE8 equal to the coherent FILE8", tier="P1", modules=["MOTION"],
         source="spec 05c C12; requirements 6.3, D-MOT-06")
def vp32(workdir: Path) -> VPResult:
    r = VPResult()
    workdir = Path(workdir)
    fs = SSI_FNUM * DF
    H = np.column_stack([np.ones(len(fs)), shear_chain_tf(fs, **TWO_DOF)])
    acc = 0.3 * synthetic_motion(2400, DT)
    write_motion_file(workdir / "eq.acc", acc, DT)
    for nm in ("FILE8", "FILE8_coh", "FILE8_01"):
        write_file8(workdir / nm, SSI_FNUM, DF, [1, 2, 3], [1, 1, 1], H, nfft=NFFT, delt=DT)
    base = dict(nout=[(3, 1, 1, 1, 0, 0, 0, 1)], thfile="eq.acc", nft=NFFT, delt=DT, interp=1, cplx=1)
    rc = run_deck(workdir, "coh", motion_deck(**base))
    r.require("coherent MOTION run succeeded", rc == 0)
    _, Hcoh_i, _ = textfiles.read_tf(workdir / nodal_result_name(3, 1, "TFI"))
    _, Hcoh_u, _ = textfiles.read_tf(workdir / nodal_result_name(3, 1, "TFU"))
    a_coh, _ = textfiles.read_history(workdir / nodal_result_name(3, 1, "ACC"))
    for p, lines in ((1, "1 1\nFILE8_coh\nFILE8_01\n"), (0, "1 0\nFILE8_01\n")):
        (workdir / "SRSSTF.txt").write_text(lines, encoding="utf-8")
        rc = run_deck(workdir, f"srss{p}", motion_deck(srss=1, **base))
        r.require(f"SRSS MOTION run, phase option {p}, succeeded", rc == 0)
        if rc != 0:
            continue
        _, Hs_i, _ = textfiles.read_tf(workdir / nodal_result_name(3, 1, "TFI"))
        _, Hs_u, _ = textfiles.read_tf(workdir / nodal_result_name(3, 1, "TFU"))
        if p == 1:
            r.check("p = 1: SRSS .TFI equals the coherent .TFI (max rel.)", _relmax(Hs_i, Hcoh_i), 0.0, atol=1e-14)
            r.check("p = 1: SRSS .TFU equals the coherent .TFU (max rel.)", _relmax(Hs_u, Hcoh_u), 0.0, atol=1e-14)
            a_s, _ = textfiles.read_history(workdir / nodal_result_name(3, 1, "ACC"))
            r.check("p = 1: SRSS acceleration history equals the coherent one (max rel.)", _relmax(a_s, a_coh), 0.0,
                    atol=1e-12)
        else:
            r.check("p = 0: SRSS .TFI equals |coherent .TFI| (max rel.)", _relmax(Hs_i, np.abs(Hcoh_i)), 0.0, atol=1e-14)
            r.check("p = 0: SRSS .TFI phase (max |phase|, rad)", float(np.max(np.abs(np.angle(Hs_i)))), 0.0, atol=0.0)
            r.check("p = 0: SRSS .TFU equals |coherent .TFU| (max rel.)", _relmax(Hs_u, np.abs(Hcoh_u)), 0.0, atol=1e-14)
    r.notes.append("'exact' is checked to 1e-14 relative: sqrt(|H|^2) and |H| e^{i arg H} differ from H by a few ulp")
    return r


# =======================================================================================
# VP-33 phase adjustment (P1)
# =======================================================================================
@problem("VP-33", "Phase adjustment (option 6 zero phase; option 2 phases scaled by 1/(1+S))", tier="P1",
         modules=["MOTION"], source="spec 05c C13; requirements 6.3, D-MOT-05")
def vp33(workdir: Path) -> VPResult:
    r = VPResult()
    workdir = Path(workdir)
    fs = SSI_FNUM * DF
    fk = S.fourier_grid(NFFT, DT)
    H = np.column_stack([np.ones(len(fs)), shear_chain_tf(fs, **TWO_DOF)])
    write_file8(workdir / "FILE8", SSI_FNUM, DF, [1, 2, 3], [1, 1, 1], H, nfft=NFFT, delt=DT)
    write_motion_file(workdir / "eq.acc", 0.3 * synthetic_motion(2400, DT), DT)
    base = dict(nout=[(3, 1, 1, 0, 0, 0, 0, 0)], thfile="eq.acc", nft=NFFT, delt=DT, cplx=1, pzadj=1)
    # option 6: |H| exactly (zero phase)
    rc = run_deck(workdir, "o6", motion_deck(interp=6, **base))
    r.require("option 6 + phase adjustment: MOTION run succeeded", rc == 0)
    H6 = I.interpolate_tf(fs, H[:, 2], fk, 6, h0=1.0)
    _, Hm, _ = textfiles.read_tf(workdir / nodal_result_name(3, 1, "TFI"))
    K = len(Hm)
    r.check("option 6: max |phase| of .TFI (rad)", float(np.max(np.abs(np.angle(Hm)))), 0.0, atol=0.0)
    r.check("option 6: .TFI amplitude equals |H| of the spline (max rel.)", _relmax(np.abs(Hm), np.abs(H6[:K])), 0.0,
            atol=1e-14)
    # option 2, S = 1000: unwrapped phase scaled by 1/(1+S)
    Sm = 1000.0
    rc = run_deck(workdir, "o2", motion_deck(interp=2, smo=Sm, **base))
    r.require("option 2, S = 1000 + phase adjustment: MOTION run succeeded", rc == 0)
    H2 = I.interpolate_tf(fs, H[:, 2], fk, 2, smooth=Sm, h0=1.0)
    ref = np.abs(H2) * np.exp(1j * np.unwrap(np.angle(H2)) / (1.0 + Sm))
    _, Hm, _ = textfiles.read_tf(workdir / nodal_result_name(3, 1, "TFI"))
    r.check("option 2, S = 1000: .TFI = |H| exp(i phi/(1+S)) (max rel.)", _relmax(Hm, ref[:K]), 0.0, atol=1e-12)
    # options 0-5 with S = 0: phase adjustment leaves H unchanged (rho = 1)
    worst = 0.0
    for opt in range(6):
        Ha = I.interpolate_tf(fs, H[:, 2], fk, opt, pzadj=1, h0=1.0)
        Hb = I.interpolate_tf(fs, H[:, 2], fk, opt, pzadj=0, h0=1.0)
        worst = worse(worst, _relmax(Ha, Hb))
    r.check("options 0-5, S = 0: phase adjustment leaves H unchanged (max rel.)", worst, 0.0, atol=1e-12)
    # reversed control direction (ang = 180 deg: H -> -H, anchor -1): the phase is referred to f = 0, so the
    # adjusted TF is reversed too and keeps the rigid-body value -1 (D-MOT-05 'reference f = 0', D-MOT-03)
    for opt, Sm in ((6, 0.0), (2, 10.0), (2, 1000.0)):
        pos = I.interpolate_tf(fs, H[:, 2], fk, opt, smooth=Sm, pzadj=1, h0=1.0)
        neg = I.interpolate_tf(fs, -H[:, 2], fk, opt, smooth=Sm, pzadj=1, h0=-1.0)
        r.check(f"option {opt}, S = {Sm:g}, reversed input: max |PA(-H) + PA(H)| / max|H|", _relmax(neg, -pos), 0.0,
                atol=1e-12)
        r.check(f"option {opt}, S = {Sm:g}, reversed input: |H'(0) - (-1)|", abs(neg[0] + 1.0), 0.0, atol=1e-12)
    return r
