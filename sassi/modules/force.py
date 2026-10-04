"""FORCE module: frequency-domain load vectors of external (vibration) loads -> FILE9
(requirements 2.1, 4.5; D-FRC-01..04; spec 05b section 9).

The load model (for the structural engineer)
--------------------------------------------
All nodal loads share one *reference* load history f_ref(t) (given later as the MOTION/STRESS
THFILE).  A load component k (a force F or a moment MM at node n, in global axes) is the reference
history scaled by a **factor** a_k and delayed by an **arrival time** t_k::

    p_k(t) = a_k f_ref(t - t_k)

With the harmonic factor exp(+i w t) a delay is a phase lag (D-CNV-01), so for every SSI frequency
w_j = 2 pi f_j FORCE stores the complex amplitude (D-FRC-01)::

    P_k(w_j) = a_k exp(-i w_j t_k)

ANALYS solves the SSI system with P in place of the seismic load and writes the transfer functions
per unit reference load; MOTION then forms ``IFFT[H(w) F_ref(w)]``.  Moving loads are modelled by
giving the loads along the path increasing arrival times (VP-40).

Repeated definitions of the same (node, DOF) follow MOPT <force> (deck ``mforce``): 1 overwrite
(the last definition wins), 0 add (D-FRC-02 / D-MDL-08): the factors are added and the arrival time
of the *latest* (most recent, in row order) definition with a non-zero factor is kept, with a warning
when it differs from the earlier one.  This is exactly the rule of the F/MM commands of the
interpreter (sassi/prep/commands/loads.py; spec 08 section 9.5 "the arrival time of the new command
replaces the old one"), so a load sequence gives the same FILE9 whether it reaches FORCE through the
interpreter or as repeated deck rows.  A row with a zero factor adds nothing and never changes the
phase of an existing load.  (A sum of two differently delayed loads a1 f(t-t1) + a2 f(t-t2) cannot be
stored as one (factor, arrival) pair; the add mode is therefore an approximation the user is warned
about.)

FILE9 holds ``P[nF, nLd]`` for the loaded (node, dof) pairs; the user copies it to FILE9001 ...
FILE9500 for simultaneous load cases (D-FRC-03).  FREAD/MREAD are not available (D-FRC-04).

    python -m sassi.modules.force < FORCE.inp      # batch use
"""
from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np

from ..conventions import frequency_step, is_power_of_two
from ..core.house_lib import DOF_NAMES
from ..io import decks
from ..io.files import write_container
from .base import ModuleContext, ModuleError, batch_main

NAME = "FORCE"

#: maximum NFFT of FORCE (D-CNV-11)
MAX_NFFT = 32768


def combine_loads(rows: List[dict], mforce: int) -> Tuple[Dict[Tuple[int, int], Tuple[float, float]], List[str]]:
    """One (factor, arrival) pair per loaded (node, dof) from the deck rows, in row order.

    ``mforce`` = MOPT <force>: 1 overwrite (last definition wins), 0 add (D-FRC-02 / D-MDL-08):
    factors are added; the arrival time of the most recent definition with a non-zero factor
    replaces the earlier one (warning when an earlier non-zero load had a different time); a
    zero-factor row never changes an existing load.  Same rule as the F/MM interpreter commands.
    Returns ``({(node, dof): (factor, arrival)}, warnings)``."""
    loads: Dict[Tuple[int, int], Tuple[float, float]] = {}
    warns: List[str] = []
    for r in rows:
        key = (int(r["node"]), int(r["dof"]))
        a, t = float(r["factor"]), float(r["arrival"])
        if key not in loads:
            loads[key] = (a, t)
            continue
        a0, t0 = loads[key]
        lab = f"node {key[0]} {DOF_NAMES[key[1]]}"
        if int(mforce) == 0:
            # add mode: only a definition that actually carries load (a != 0) sets the arrival time
            tn = t if a != 0.0 else t0
            loads[key] = (a0 + a, tn)
            msg = f"{lab} defined more than once: factors added ({a0:.6g} + {a:.6g} = {a0 + a:.6g})"
            if a != 0.0 and a0 != 0.0 and t != t0:
                msg += (f", arrival time {t0:.6g} s replaced by {t:.6g} s of the latest definition "
                        "(a sum of differently delayed loads is approximated by one arrival time)")
            elif a == 0.0 and t != t0:
                msg += f"; the zero-factor definition does not change the arrival time {t0:.6g} s"
            warns.append(msg + " (MOPT <force> = 0, D-FRC-02)")
        else:
            loads[key] = (a, t)
            warns.append(f"{lab} defined more than once: the last definition (factor {a:.6g}, arrival {t:.6g} s) "
                         "replaces the earlier one (MOPT <force> = 1)")
    return loads, warns


def load_amplitudes(freq: np.ndarray, factor: np.ndarray, arrival: np.ndarray) -> np.ndarray:
    """``P[j, k] = a_k exp(-i 2 pi f_j t_k)`` (requirements 4.5, D-FRC-01), shape (nF, nLd)."""
    f = np.asarray(freq, dtype=float)[:, None]
    return np.asarray(factor, float)[None, :] * np.exp(-2j * np.pi * f * np.asarray(arrival, float)[None, :])


def run(ctx: ModuleContext) -> int:
    lst = ctx.listing
    d = decks.read(ctx.deck_path, NAME)
    # ---- frequencies (same rules as SITE: f_j = n_j df) -------------------------------------
    delt, nft, fstep = float(d["delt"]), int(d["nft"]), float(d["fstep"])
    df = float(d["df"])
    if df <= 0:
        if fstep < 0:
            raise ModuleError(f"Error 48: frequency step {fstep} < 0")
        if fstep == 0 and (delt <= 0 or nft <= 0):
            raise ModuleError(f"Errors 49/50: time step {delt} and NFFT {nft} must be > 0 when the frequency "
                              "step is 0")
        df = frequency_step(delt, nft, fstep)
    warns: List[str] = []
    if fstep <= 0:                       # time-history analysis: df = 1/(delt NFFT)
        if nft > MAX_NFFT:
            raise ModuleError(f"NFFT {nft} exceeds the FORCE maximum {MAX_NFFT} (D-CNV-11)")
        if nft > 0 and not is_power_of_two(nft):
            warns.append(f"NFFT {nft} is not a power of 2 (Warning 9; AFWRITE writes the nearest power of 2)")
    fnum = [int(r["number"]) for r in d.rows("freqs")]
    if not fnum:
        raise ModuleError(f"Error 120: frequency set {int(d['freq'])} is empty")
    if len(set(fnum)) != len(fnum):
        raise ModuleError("frequency set contains duplicate frequency numbers")
    if min(fnum) < 1:
        raise ModuleError("frequency numbers must be positive integers")
    fnum_a = np.array(sorted(fnum), dtype=np.int64)
    freq = fnum_a * df
    if nft > 0 and fnum_a[-1] > nft // 2:
        warns.append(f"frequency number {int(fnum_a[-1])} > NFFT/2 = {nft // 2} (G-02, D-CNV-12)")
    # ---- loads ------------------------------------------------------------------------------
    rows = d.rows("loads")
    if not rows:
        raise ModuleError("Error 61: no forces or moments defined (F / MM)")
    for r in rows:
        if not 1 <= int(r["dof"]) <= 6:
            raise ModuleError(f"load on node {int(r['node'])}: dof {int(r['dof'])} must be 1..6 "
                              "(1-3 forces F, 4-6 moments MM)")
        if int(r["node"]) <= 0:
            raise ModuleError(f"load on undefined node {int(r['node'])} (Error 62)")
        if float(r["arrival"]) < 0:
            warns.append(f"node {int(r['node'])} {DOF_NAMES[int(r['dof'])]}: negative arrival time "
                         f"{float(r['arrival']):.6g} s (the load starts before the reference load)")
    mforce = int(d["mforce"])
    if mforce not in (0, 1):
        raise ModuleError(f"MOPT <force> = {mforce} must be 0 (add) or 1 (overwrite)")
    loads, w2 = combine_loads(rows, mforce)
    warns.extend(w2)
    keys = sorted(loads)
    zero = [k for k in keys if loads[k][0] == 0.0]
    keys = [k for k in keys if loads[k][0] != 0.0]
    if zero:
        warns.append(f"{len(zero)} load components with a zero total factor are not written: "
                     + ", ".join(f"node {n} {DOF_NAMES[c]}" for n, c in zero[:10]))
    if not keys:
        raise ModuleError("Error 61: every load factor is zero")
    node = np.array([k[0] for k in keys], dtype=np.int64)
    dof = np.array([k[1] for k in keys], dtype=np.int64)
    fac = np.array([loads[k][0] for k in keys], dtype=float)
    tarr = np.array([loads[k][1] for k in keys], dtype=float)
    P = load_amplitudes(freq, fac, tarr)
    # ---- listing ------------------------------------------------------------------------------
    lst.section("FORCE input")
    lst.write(f" Title                         : {d['title']}")
    lst.write(f" Frequency step                : {df:.9g} Hz"
              + (" (1/(delt NFFT))" if fstep <= 0 and float(d["df"]) <= 0 else ""))
    lst.write(f" Time step / NFFT              : {delt:.6g} s / {nft}")
    lst.write(f" Frequency set                 : {int(d['freq'])}: {fnum_a.size} frequencies, "
              f"{freq[0]:.6g} .. {freq[-1]:.6g} Hz")
    rule = "overwrite (last wins)" if mforce == 1 else "add factors, arrival of the latest non-zero definition"
    lst.write(f" Repeated definitions          : {rule}  (MOPT <force> = {mforce})")
    lst.section("Loads: P(f) = factor * exp(-i 2 pi f t_arrival)")
    lst.write(f"{'node':>8s}{'load':>6s}{'factor':>14s}{'arrival (s)':>14s}{'|P|':>12s}"
              f"{'phase f1 (deg)':>16s}{'phase fN (deg)':>16s}")
    for k in range(len(keys)):
        lst.write(f"{node[k]:>8d}{DOF_NAMES[int(dof[k])]:>6s}{fac[k]:>14.6g}{tarr[k]:>14.6g}{abs(fac[k]):>12.6g}"
                  f"{np.degrees(np.angle(P[0, k])):>16.4f}{np.degrees(np.angle(P[-1, k])):>16.4f}")
    lst.write(f" {len(keys)} loaded DOFs on {np.unique(node).size} nodes")
    for w in warns:
        lst.warning(w)
    if int(d["opmode"]) == 1:
        lst.write("")
        lst.write(" Data check only (<opmode> = 1): FILE9 not written")
        return 0
    arrays = {"fnum": fnum_a, "freq": freq, "load_node": node, "load_dof": dof, "P": P,
              "x_factor": fac, "x_arrival": tarr}
    meta = {"df": df, "mforce": mforce, "nft": nft, "delt": delt, "model": ctx.model, "title": str(d["title"]),
            "convention": "P = factor * exp(-i w t_arrival), harmonic factor exp(+i w t) (D-FRC-01)"}
    write_container(ctx.path("FILE9"), "FILE9", arrays, meta, module=NAME)
    lst.section("Files written")
    lst.write(" FILE9  (copy to FILE9001 ... FILE9500 for simultaneous load cases, D-FRC-03)")
    ctx.progress(1.0, "FORCE done")
    return 0


if __name__ == "__main__":
    raise SystemExit(batch_main(NAME))
