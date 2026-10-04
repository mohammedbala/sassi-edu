"""AFWRITE and the module run commands RUN<MODULE> (requirements sections 2.6, 3.4.A, 3.4.G;
D-AFW-01..06, D-RUN-01..05; ARCHITECTURE section 8).

``AFWRITE`` runs CHECK and writes the deck of every AOPT-enabled module without errors
(:func:`sassi.prep.afwrite.afwrite`).  ``RUN<MODULE>,[model]`` (model = number of a model in
memory, default -1 = the active one) runs :func:`sassi.modules.base.run_module` synchronously in
the model directory and streams the listing lines as information messages (D-RUN-05).  Before the
run the prerequisites of requirements section 2.6 are checked (D-RUN-03): model name and path, the
deck written by AFWRITE, no CHECK error for the module in the last AFWRITE, the input files of
the upstream modules (section 2.1 table) and, for ANALYS, the frequency survey (every frequency to
solve present in FILE1/FILE9 and FILE3).  RUNSITE is implemented with the same semantics (D-RUN-01).
"""
from __future__ import annotations

import importlib.util
import warnings
from pathlib import Path
from typing import List, Optional, Sequence

from ...conventions import DOF_TAGS, nodal_result_name
from ...io import decks
from ...io.container import read_container
from ..afwrite import AfwriteError, afwrite
from ..registry import CommandReported, command
from .checks import check_options, search_dirs

RUN_MODULES = ("EQUAKE", "SOIL", "SITE", "POINT", "HOUSE", "FORCE", "ANALYS", "COMBIN", "MOTION", "RELDISP",
               "STRESS")


def _afw_state(interp, number: int) -> dict:
    return interp.session.setdefault("afwrite", {}).setdefault(number, {})


@command("AFWRITE", max_args=0)
def cmd_afwrite(c):
    """AFWRITE: CHECK, then write <model>.<ext> for every AOPT-enabled module without errors (D-AFW-01)."""
    m = c.model
    opt = check_options(c.interp)
    try:
        res = afwrite(m, dirs=search_dirs(c), check_options=opt, write_options=c.interp.write_options)
    except AfwriteError as exc:
        c.fail(str(exc))
    rep = res.report
    c.interp.session["check_report"] = rep
    st = _afw_state(c.interp, c.interp.active_model)
    st.update(hash=res.model_hash, written={k: str(v) for k, v in res.written.items()},
              blocked=dict(res.blocked), report=rep)
    if rep.messages and not opt.suppress_window:
        for ln in rep.format(opt).rstrip("\n").splitlines():
            c.info(ln)
    for note in res.notes:
        c.info(f"AFWRITE: {note}")
    for mod, path in res.written.items():
        c.confirm(f"{mod}: {Path(path).name} written")
    for mod, why in res.blocked.items():
        c.warn(f"{mod}: input deck not written ({why}; see {res.err_file.name if res.err_file else '.err'})")
    for w in res.warnings:
        c.warn(w)
    if res.sim_file:
        c.info(f"AFWRITE: simulation commands written to {res.sim_file.name} (D-AFW-04)")
    c.info(f"AFWRITE: {len(res.written)} decks written, {len(res.blocked)} modules with errors; "
           f"{rep.summary()}")


# ======================================================================================
# RUN<MODULE>
# ======================================================================================
def _module_available(module: str) -> bool:
    try:
        return importlib.util.find_spec(f"sassi.modules.{module.lower()}") is not None
    except (ImportError, ValueError):
        return False


def _read_deck(path: Path, module: str):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if module == "SOIL":
            try:
                from ...modules.soil import read_deck     # schema workaround of the SOIL package
                return read_deck(path)
            except ImportError:
                pass
        return decks.read(path, module)


def _fnums(path: Path) -> Optional[List[int]]:
    """Frequency numbers ``fnum`` of a binary inter-module file (ARCHITECTURE section 4), or None
    when the file cannot be read as a container (e.g. a legacy or foreign file)."""
    try:
        c = read_container(path)
    except (OSError, ValueError):
        return None
    f = c.get("fnum")
    return None if f is None else [int(v) for v in f]


def _numbers(nums: Sequence[int]) -> str:
    nums = list(nums)
    return ", ".join(str(n) for n in nums[:10]) + (" ..." if len(nums) > 10 else "")


def frequency_survey(mdir: Path, deck, load_files: Sequence[str], producer: str) -> List[str]:
    """ANALYS frequency survey (requirements section 2.6, D-RUN-03).

    Every SSI frequency ANALYS will solve must be present in the free-field / load file it reads
    (FILE1 or FILE1X/Y/Z from SITE; FILE9 or FILE9nnn from FORCE) **and** in FILE3 (POINT), else
    the run stops: "Frequency n not in FILE1/FILE3 -- re-run SITE/POINT".  With ``<fopt>`` = 0 the
    frequencies are those of the deck (the selected frequency set); with ``<fopt>`` = 1 they are all
    the frequencies of the FILE1/FILE9 files, which must then all be in FILE3.  A file that is not a
    readable container is skipped here (ANALYS reports it).
    """
    out: List[str] = []
    why = "ANALYS frequency survey: every frequency to solve must be in {} and FILE3"
    requested = [int(r["number"]) for r in deck.rows("freqs")] if int(deck["fopt"]) == 0 else None
    solved: List[int] = list(requested or [])
    for name in load_files:
        have = _fnums(mdir / name)
        if have is None:
            continue
        if requested is None:                                   # <fopt> = 1: all frequencies of the file
            solved.extend(v for v in have if v not in solved)
            continue
        missing = sorted(set(requested) - set(have))
        if missing:
            out.append(f"Frequency {_numbers(missing)} not in {name} -- re-run {producer} ({why.format(name)})")
    have3 = _fnums(mdir / "FILE3")
    if have3 is not None and solved:
        missing = sorted(set(solved) - set(have3))
        if missing:
            src = "/".join(load_files) if load_files else "FILE1"
            out.append(f"Frequency {_numbers(missing)} not in FILE3 -- re-run POINT ({why.format(src)})")
    return out


def prerequisites(module: str, model_name: str, mdir: Path, deck) -> List[str]:
    """Missing inputs of a module run (requirements sections 2.1 and 2.6 tables), as messages."""
    out: List[str] = []

    def need(name: str, msg: str) -> None:
        if not (mdir / name).exists():
            out.append(msg)

    if module in ("POINT",):
        need("FILE2", "FILE2 missing -- run SITE Mode 1")
    elif module == "SITE":
        if deck is not None and not int(deck["mode1"]) and int(deck["mode2"]):
            need("FILE2", "FILE2 missing -- run SITE Mode 1 (SITE Mode 2 alone needs FILE2)")
        if deck is not None and int(deck["soilmode"]) == 1:
            need("FILE88", "FILE88 missing -- run SOIL first (SITEX Non-Linear Soil)")
    elif module == "HOUSE":
        need(f"{model_name}.sit", f"HOUSE needs {model_name}.sit (run AFWRITE with SITE enabled)")
    elif module == "ANALYS" and deck is not None:
        need("FILE3", "FILE3 missing -- run POINT")
        for f in (f"{model_name}.N4", "COOSK", "COOSM", "DOFSMAP"):     # FILE4 family (section 2.1)
            need(f, f"{f} missing -- run HOUSE")
        simul, coh, me = int(deck["simul"]), int(deck["coh"]), int(deck["me"])
        load_files: List[str] = []
        if int(deck["type"]) == 0:
            producer = "SITE"
            xyz = all((mdir / f).exists() for f in ("FILE1X", "FILE1Y", "FILE1Z"))
            if coh == 1:
                if not xyz and not (mdir / "FILE1").exists():
                    out.append("FILE1 (or FILE1X/FILE1Y/FILE1Z) missing -- run SITE and copy with FCOPY")
                load_files = ["FILE1X", "FILE1Y", "FILE1Z"] if xyz else ["FILE1"]
            elif simul >= 1:
                for f in ("FILE1X", "FILE1Y", "FILE1Z"):
                    need(f, f"{f} missing -- run SITE (SV x' / SH y' / P z') and FCOPY,FILE1,{f}")
                load_files = ["FILE1X", "FILE1Y", "FILE1Z"]
            else:
                need("FILE1", "FILE1 missing -- run SITE Mode 2")
                load_files = ["FILE1"]
        else:
            producer = "FORCE"
            if simul > 1:
                for k in range(1, simul + 1):
                    need(f"FILE9{k:03d}", f"FILE9{k:03d} missing -- run FORCE and FCOPY,FILE9,FILE9{k:03d}")
                load_files = [f"FILE9{k:03d}" for k in range(1, simul + 1)]
            else:
                need("FILE9", "FILE9 missing -- run FORCE")
                load_files = ["FILE9"]
        if not out:
            out.extend(frequency_survey(mdir, deck, load_files, producer))
        if coh == 1 or me == 1:
            if simul > 1 and coh == 1:
                for k in range(1, simul + 1):
                    need(f"FILE77{k:03d}", f"FILE77{k:03d} missing -- run HOUSE (incoherent simulations)")
            else:
                need("FILE77", "FILE77 missing -- run HOUSE (incoherency / multiple excitation)")
    elif module == "COMBIN":
        need("FILE81", "FILE81 missing -- rename the first FILE8 (FMOVE,FILE8,FILE81)")
        need("FILE82", "FILE82 missing -- rename the second FILE8 (FMOVE,FILE8,FILE82)")
    elif module in ("MOTION", "STRESS") and deck is not None:
        f8 = deck["file8"] or "FILE8"
        only_cnvrt = module == "MOTION" and int(deck["cnvrt"]) and not deck.rows("nout")
        if not only_cnvrt:
            need(f8, f"{f8} missing -- run ANALYS (or COMBIN)")
        if module == "STRESS":
            need(f"{model_name}.N4", f"{model_name}.N4 (FILE4) missing -- run HOUSE")
    elif module == "RELDISP" and deck is not None:
        rel = (deck["relfile"] or "").strip()
        if rel and rel.upper() != "FREEFIELD":
            p = Path(rel)
            if not (p if p.is_absolute() else mdir / p).exists():
                out.append(f"reference TF {rel} missing -- run MOTION with Save Complex TF")
        maxnode = int(deck["maxnode"])
        missing = []
        for r in deck.rows("rdnd"):
            for dof, col in enumerate(("x", "y", "z", "xx", "yy", "zz"), start=1):
                if int(r[col]) >= 1:
                    name = nodal_result_name(int(r["node"]), dof, "TFI", maxnode)
                    if not (mdir / name).exists():
                        missing.append(f"{int(r['node'])} {DOF_TAGS[dof]}")
        if missing:
            out.append("Run MOTION with Save Complex TF for nodes " + ", ".join(missing[:10])
                       + (" ..." if len(missing) > 10 else ""))
    return out


def run_command(c, module: str) -> None:
    """RUN<MODULE>,[model]: run one module synchronously in the model directory (D-RUN-05)."""
    interp = c.interp
    num = c.int(1, default=-1)
    if num is None or num < 0:
        num = interp.active_model
    if num not in interp.models:
        c.fail(f"model {num} is not in memory")
    m = interp.models[num]
    if not m.name or not m.path:
        c.fail("Model name/path not defined -- use MDL")
    mdir = Path(m.path)
    if not _module_available(module):
        c.fail(f"module {module} is not available in this build yet (sassi/modules/{module.lower()}.py "
               f"not found)")
    deck = None
    if module != "COMBIN":
        dp = decks.deck_path(mdir, m.name, module)
        if not dp.exists():
            c.fail(f"{dp.name} not found -- Run AFWRITE first (with {module} enabled in AOPT)")
        st = interp.session.get("afwrite", {}).get(num, {})
        if module in st.get("blocked", {}):
            c.fail(f"{module} has CHECK errors in the last AFWRITE ({st['blocked'][module]}); correct them and "
                   f"run AFWRITE again")
        if st.get("hash") and st["hash"] != m.model_hash():
            c.warn("the model changed since the last AFWRITE: the deck may be outdated (run AFWRITE)")
        try:
            deck = _read_deck(dp, module)
        except (OSError, ValueError, KeyError) as exc:
            c.fail(f"cannot read {dp.name}: {exc}")
    missing = prerequisites(module, m.name, mdir, deck)
    if missing:
        for msg in missing:
            c.error(msg)
        raise CommandReported(f"{module}: prerequisites missing")
    from ...modules.base import run_module
    c.info(f"RUN{module}: model {m.name} in {mdir}")
    rc = run_module(module, m.name, mdir, echo=lambda line: c.info(line))
    listing = mdir / f"{m.name}_{module}.out"
    if rc != 0:
        c.fail(f"{module} finished with status FAILED ({rc}); see {listing.name}")
    c.confirm(f"{module} finished (listing {listing.name})")


def _make(module: str):
    def handler(c):
        run_command(c, module)
    handler.__name__ = handler.__qualname__ = f"cmd_run{module.lower()}"
    handler.__doc__ = f"RUN{module},[model]: run module {module} in the model directory (requirements 2.6)."
    return handler


for _mod in RUN_MODULES:
    command(f"RUN{_mod}", max_args=1)(_make(_mod))
