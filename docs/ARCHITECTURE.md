# SASSI-EDU Architecture and Implementation Contracts

SASSI-EDU is an educational Python re-implementation of the ACS SASSI Version 3
methodology (flexible-volume SSI after Lysmer et al., 1981). This document is the
**binding contract** between the parts of the code base. The requirements are in
`docs/spec/00_requirements.md`, which wins on *what* to build; this document wins on
*how the parts connect*. The theory is in `docs/spec/R1_sassi_theory.md`, with
prototype numerics in `docs/spec/R1_checks/`. The verification references are in
`docs/spec/R2_benchmarks.md`.

## 1. Principles

1. **Same nomenclature as ACS SASSI.** Module names (EQUAKE, SOIL, SITE, POINT, HOUSE,
   FORCE, ANALYS, COMBIN, MOTION, RELDISP, STRESS), command names and argument orders,
   option codes, file names (FILE1, FILE2, FILE3, FILE4 = `<model>.N4`, COOSK, COOSM,
   FILE8, FILE9, FILE77, `.TFU`, `.TFI`, `.ACC`, `.RS`, ...), dialog labels and CHECK
   messages are as in the manual (see the glossary in requirements §1).
2. **Same architecture as the original.** A command interpreter ("UI") holds the model.
   AFWRITE writes one input **deck** per module. Each module is an independent program
   that reads its deck plus binary files written by upstream modules, and writes its
   own files into the model directory. Modules never import the interpreter or the
   model database.
3. **One interpreter.** The command console, `.pre` files, macros and the web GUI all
   submit command text to the same `sassi.prep.Interpreter` (requirements UI-01, L17).
4. **Verifiable.** Every physics module ships with verification problems (VPs) that
   compare computed values with published or exact references and print computed,
   reference, error and tolerance (requirements §6). These are registered in
   `sassi.verify` and run by pytest and by the `VERIFY` command.
5. **Pure Python + NumPy/SciPy.** No compiled extensions. Must run on Python 3.9.

## 2. Package layout and ownership

```
sassi/
  __init__.py, conventions.py          (shared; owner: lead)   damping, units, naming formulas
  cli.py                               console entry `sassi` (interpreter REPL / `sassi run file.pre`)
  io/
    container.py, files.py             (lead) binary containers + FILEnn schemas      §4
    deckfmt.py, decks.py               (lead) module deck format + schemas             §3
    textfiles.py, thfile.py            (lead) text result files, control-motion files §5
    ansys_cdb.py, apdl.py              ANSYS .cdb import (CONVERT,ANSYS), APDL export (ANSYS command)
  model/                               interpreter's model database (SSIModel, tables, options)
  prep/                                lexer, interpreter, registry (command catalogue), WRITE, CHECK, AFWRITE
    commands/*.py                      command handlers by category (session, nodes, elements, tables, loads,
                                       options, checks, conditioning, generation, cuts, conversion, plotting,
                                       linemath, frames, program, modules_cmd (RUNxxx), extensions,
                                       incoherency_cmds, nlsoil_cmds, soilnon_cmds, nonlinear_cmds (Option NON),
                                       loadgen_cmds (Option A), water, thickshell, afwrbat, verifyreport)
    generation_lib.py cuts_lib.py water_lib.py   algorithms of INTGEN/EXCAV/SOILMESH..., cuts, FILLPOOL...
  core/
    spectra.py                         (lead) Nigam-Jennings RS, Arias, PSD, integration
    tlm.py                             thin-layer layer matrices, column assembly, VDM half-space, eigenproblems
    freefield.py                       SITE Mode 2 free field + ANALYS evaluation at interaction nodes
    axisym.py                          POINT3 central zone (core element, exterior modes, boundary stiffness)
    strip2d.py                         POINT2 (plane-strain core + Waas-Lysmer boundaries)
    flexibility.py                     F_ff assembly from FILE3 (used by ANALYS)
    greens_tlm.py                      exact Kausel thin-layer Green functions (reference/tests)
    ssi_solver.py                      ANALYS numerics: impedance, assembly, Schur solve, restart factors
    ssi2d.py, symmetry.py              2D plane-strain SSI helpers; SYMM image superposition
    coherency.py, incoherency.py       coherency models (data in sassi/data/coherency), FILE77 synthesis
    house_lib.py, renumber.py          HOUSE helpers; node-numbering optimizer (RCM)
    interp.py                          TF interpolation options 0-6, smoothing, phase adjustment
    signal.py                          FFT helpers, baseline correction (Hudson-Housner), convolution
    stress_lib.py                      STRESS recovery, derived quantities, frames
    shake.py                           SHAKE equivalent-linear recursion (SOIL)
    soilnon.py                         SOIL-NON hyperbolic/Masing time-domain site response
    nlsoil.py                          near-field nonlinear-soil SSI iterations (.pin, FILE74/78)
    hysteresis.py, panels.py           Option NON hysteresis models (GMR/CMS/TAK), wall panels, SHEAR/BBCGEN
    equake_lib.py                      EQUAKE spectral matching (LW + wavelets), SRP 3.7.1 checks
  elements/
    base.py                            MaterialProps, ElementSpec registry
    solid.py plane.py beam.py shell.py tshell.py spring.py general.py
    assemble.py                        DOF map + sparse assembly (structure / excavated parts)
  modules/
    base.py                            (lead) ModuleContext, Listing, run_module, batch_main
    equake.py soil.py site.py point.py house.py force.py analys.py combin.py motion.py reldisp.py stress.py
    nonlinear.py (Option NON)  loadgen.py (Option A)
  plotting/                            plot/line state shared by console and GUI; matplotlib rendering
  verify/
    __init__.py                        (lead) VPResult, @problem registry, run_problem, worse()
    builders.py                        canonical SSI test models written as decks
    report.py                          verification manual and command reference generators
    problems/vp_*.py                   one file per VP group
  ui/                                  local web GUI (server.py, api.py, dialogs.py, helpdocs.py + static/)
  data/                                coherency coefficients (with citations), DYNP curve library
tests/unit, tests/verification, tests/integration   pytest (verification tests call sassi.verify problems)
```

Each source file has one owner per implementation wave; never edit a file owned by
another work package. If you need something from another package that does not exist
yet, write a small private helper in your own file and note it in your report.

## 3. Module input decks (AFWRITE -> module)

The deck format and the schema of every deck are code: `sassi/io/deckfmt.py`
(format) and `sassi/io/decks.py` (`SCHEMAS`). AFWRITE must create decks with
`decks.new(module)`, fill **every** parameter and table, and call `decks.write`.
Modules read them with `decks.read(path, module)`. Rules:

* The deck carries *resolved* data. Examples: SITE gets the TOPL layers with their L
  properties and the half-space row; HOUSE gets global node coordinates (local
  coordinate systems resolved), the element table with ETYPE 0 resolved to 1/2 (D-AFW-06),
  the material, layer, section, spring and matrix tables, nodal masses and the
  interaction-node list; every deck that needs frequencies gets the frequency numbers of
  the selected set and the resolved `df`.
* Shared variables (requirements §3 table "Shared analysis variables") are copied into
  every deck that needs them (`delt`, `nft`, `df`, `gravity`, `type`, history data...).
* Units: the model's consistent units. `weight` is a specific weight; mass density is
  `weight / gravity` (HOUSE/SITE/POINT/FORCE/ANALYS use the HOUSE gravity; SOIL uses its own
  `grav`). Damping ratios are fractions.
* Schema changes are made only by the lead. Missing parameters read back as defaults
  with a warning; unknown parameters are an error on write.

## 4. Binary inter-module files

Containers: `sassi/io/container.py` (`write_container`, `read_container`, sparse helpers).
Schemas: `sassi/io/files.py` (`SCHEMAS`, `validate`). Producers must write every listed
array and meta key. Conventions:

* `fnum` (frequency numbers) are sorted ascending; `freq = fnum * df`.
* Depths are positive down from the top of TOPL layer 1. User interface `i` (1-based) is
  the top of TOPL layer `i`; interface `nI = nTOPL + 1` is the top of the half-space (or the
  rigid base). Interaction nodes must lie on user interfaces `1 .. POINT<layer>+1`
  (EDU-01). A node at elevation z has depth `gelev - z`.
* FILE1 stores, per wave, the unit-normalised free field at user interfaces in the SITE
  axes x'y'z' (physical components, z' up) and the horizontal wavenumber. The motion at a
  point is `sum_w ratio_w * U_w(depth) * exp(-i k_w (x' - x'_cp))` where
  `x' = (x - xc) cos(ang) + (y - yc) sin(ang)` and the result is rotated to global axes by
  `ang` (ANALYS). Vertical incidence has k = 0.
* FILE3 stores, per frequency and load interface, the POINT3 exterior mode amplitudes
  (`alpha0` for a vertical unit load, mu = 0; `alpha1` for a horizontal unit load, mu = 1),
  the core-axis displacements, and the Rayleigh/Love mode rows at the observation
  interfaces, so ANALYS can evaluate `u(rho) = Psi_mu(rho) alpha` at any distance
  (R1 §3.3, §3.6) without FILE2.
* FILE4 (`<model>.N4`) holds the DOF map (`eq_node`, `eq_dof`), interaction nodes and
  their interfaces, element topology and the STRESS recovery operators; COOSK holds the
  sparse complex stiffness `Ks` (structure incl. near-field soil) and `Ke` (excavated
  soil), COOSM the real masses `Ms`, `Me`; all on the FILE4 equation numbering.
* FILE8 holds `H[nF, nEq]` on the FILE4 equation numbering (`eq_node`, `eq_dof` repeated
  so MOTION/STRESS can work from FILE8 and FILE4 alone). Seismic: `H = U/U_cp`
  (dimensionless, equal for displacement and acceleration). Vibration: displacement per
  unit load factor.
* FILE9 holds load amplitudes per (node, dof) and frequency.

## 5. Text files

`sassi/io/textfiles.py` implements D-FIL-02 formats (.TFU/.TFI/.TFD, .ACC/.THD/.THS,
.RS/.RSO, .psd, .fft). `sassi/io/thfile.py` reads control-motion files (MOTION `<fopt>`,
records, scaling) and SOIL input files. File names follow `sassi.conventions`
(`nodal_result_name`, `nodal_rs_name`, `element_result_name`, `layer_th_name` ...).

## 6. Core numerical APIs (signatures are binding)

### 6.1 Damping, materials (`sassi.conventions`, `sassi.elements.base`)

```python
cfactor(beta, form=0) -> complex            # 1-2b^2+2ib*sqrt(1-b^2)  (form 1: 1+2ib)
complex_lame(G, M, beta_s, beta_p, form)    # -> (G*, M*, lam*)

@dataclass
class MaterialProps:          # sassi.elements.base
    rho: float                # mass density
    G: complex                # G*
    M: complex                # M* = lam* + 2G*
    beta_s: float; beta_p: float
    G0: float; M0: float      # real (undamped) moduli
    @property lam -> complex; E -> complex; nu -> complex (from G*, M*); E0, nu0 real
material_from_M(mtype, val1, val2, weight, pdamp, sdamp, gravity, form=0) -> MaterialProps
material_from_layer(vp, vs, weight, dp, ds, gravity, form=0) -> MaterialProps
```

### 6.2 Elements (`sassi.elements`)

```python
@dataclass
class ElementSpec:
    code: int                 # GROUP type: 1 SOLID, 2 BEAMS, 3 SHELL, 4 PLANE, 5 TSHELL, 7 SPRING, 9 GENERAL
    name: str                 # 'SOLID'
    nnodes: int               # nodes carrying DOFs (BEAMS: 2, the K node is orientation only)
    dofs: tuple               # dof labels per node, e.g. (1,2,3) or (1,2,3,4,5,6); PLANE (1,3)
    components: list          # sassi.conventions.ELEMENT_COMPONENTS[name]
    matrices: Callable        # (xyz, mat, **props) -> (K complex (nd,nd), M real (nd,nd))
    recovery: Callable        # (xyz, mat, **props) -> S complex (ncomp, nd); component TF = S @ u_e
ELEMENTS: dict[int, ElementSpec]
```

Element DOF vector order: node-major, `[node1 dofs..., node2 dofs..., ...]` in `ElementSpec.dofs`
order. Element-specific keyword props:

* SOLID: `incompatible: bool` (structure only), `eint: int` (0/1/2 -> 2/3/4 Gauss points),
  repeated nodes for prisms/pyramids. Mass = 1/2 lumped + 1/2 consistent. Recovery: stresses at
  the centroid in global axes (SXX SYY SZZ SXY SXZ SYZ) + SOCT computed by STRESS in time domain
  (rows of S for SOCT may be zero).
* PLANE: X-Z plane, `incompatible`, unit thickness, plane strain.
* BEAMS: `xyz` has 3 rows (I, J, K orientation node); `section=dict(A, As2, As3, J, I2, I3)`,
  `ki`, `kj` release tuples of 6 bools (P1 P2 P3 M1 M2 M3). Recovery: local end forces
  FXI..MZJ = forces exerted on the element (D-STR-05).
* SHELL: 3 or 4 nodes, `thick`; flat facet, membrane + DKQ/DKT bending, zero drilling stiffness;
  lumped mass, no rotary inertia. Recovery: FXX FYY FXY (membrane **stresses**, D-STR-10/D-W1-14) and
  MXX MYY MXY (moments per unit length) in local x'y' at the centroid.
* SPRING: `k=(kx,ky,kz,kxx,kyy,kzz)`, `damp`, `form`; K* = k * cfactor(damp); no mass. Recovery
  F = k*(u_J - u_I) per global component.
* GENERAL: `KR`, `KI`, `MM` 12x12 (full symmetric, built from upper-triangle rows), optional
  `xyz` with 3 rows for local input (transformation as BEAMS).

Assembly (`sassi.elements.assemble`):

```python
@dataclass
class ElemRecord:
    group: int; id: int; code: int; nodes: tuple; excavated: bool
    mat: Optional[MaterialProps]; props: dict
build_dofmap(node_ids, node_fix (nN,6), records, extra_dofs=None) -> DofMap
    # DofMap.eq_node, .eq_dof (arrays), .eq(node, dof) -> int or -1, .neq
    # active dofs of a node = union of dofs of its elements (and masses), minus fixed
assemble(node_xyz: dict[int, (3,)], records, dofmap, masses=None) ->
    AssembledModel(Ks, Ms, Ke, Me, recovery)   # scipy.sparse CSR on dofmap numbering
    # recovery: dict[type name] -> dict(S=(nE,nc,nd), eq=(nE,nd), group=(nE,), id=(nE,))
```

### 6.3 Layered soil (`sassi.core.tlm`, `freefield`, `axisym`, `flexibility`)

Owned by the SITE/POINT package. The functions ANALYS needs are binding:

```python
# sassi.core.freefield
free_field_at_nodes(file1: Container, q: int, xyz: (n,3), depth_iface: (n,) int (1-based),
                    ang_deg: float, xc: float, yc: float) -> (n,3) complex   # global x,y,z components
# sassi.core.flexibility
flexibility_matrix(file3: Container, q: int, xy: (n,2), iface: (n,) int (1-based)) -> (3n,3n) complex
    # node-major [ux1,uy1,uz1, ux2,...]; symmetrised (F + F.T)/2  (D-ANL-02)
```

### 6.4 Post-processing (`sassi.core.interp`, `sassi.core.signal`, `sassi.core.spectra`)

```python
# sassi.core.interp
interpolate_tf(f_ssi: (nF,), H: (nF, m) complex, f_out: (nK,), option: int, smooth: float = 0.0,
               pzadj: int = 0, h0: Optional[(m,) complex] = None, mode: str = 'seismic') -> (nK, m) complex
    # exact at f_ssi; zero above f_ssi[-1]; below f_ssi[0]: linear from h0 (seismic) or H_1 (vibration)
# sassi.core.signal
fourier_grid(nfft, delt) -> f (nfft//2+1,)
convolve(H_grid, A) -> time history via irfft
hudson_baseline(acc, dt) -> corrected acc   (D-MOT-08)
# sassi.core.spectra  (implemented)
response_spectrum(acc, dt, freqs, dampings) -> {'SA','SV','SD','PSA','PSV'} (n_damp, n_freq)
log_frequencies(f1, f2, n); integrate(acc, dt); arias_intensity; strong_motion_window; band_averaged_psd
```

## 7. Modules

Every module lives in `sassi/modules/<name>.py`, defines `NAME` and `run(ctx)` and ends with

```python
if __name__ == "__main__":
    raise SystemExit(batch_main(NAME))
```

`run(ctx)` reads `decks.read(ctx.deck_path, NAME)`, uses `ctx.require(file, producer)` for
inputs, writes outputs to `ctx.workdir`, writes a readable listing (`ctx.listing`), reports
progress with `ctx.progress(fraction, message)` and raises `ModuleError` on fatal input errors.
Algorithms are specified in requirements §4 (normative) and R1. The run sequences and restart
semantics are in requirements §2.

## 8. Interpreter and model (`sassi.prep`, `sassi.model`)

* `sassi.model.SSIModel` is the in-memory database of one numbered model: nodes (with the
  coordinate system they were defined in), coordinate systems, groups/elements with their
  per-element attributes (ETYPE, EINT, MSET, RSET, THICK, KI/KJ), tables (M, L, R, SC, MX),
  loads and masses, fixities, interaction flags, frequency sets, lists (DAMP, TOPL, AMP),
  output requests (NOUT, EOUT, RDND), module option records and extension options. The
  interpreter holds several models (ACTM).
* `sassi.prep.Interpreter.execute(line) / run_file(path)` implements the lexical rules
  L1-L18, abbreviations, variables/macros/FOREACH and dispatch to handlers in
  `sassi/prep/commands/`. Every message goes through a message sink with the six classes
  (echo, confirmation, comment, info, warning, error) so the console and GUI can filter them.
* WRITE produces a `.pre` that INP reads back to an identical state (UT-03).
* CHECK implements the Chapter 10 catalogue and writes `<model>.err`; AFWRITE runs CHECK and
  writes the decks of AOPT-enabled modules without errors (D-AFW-01).
* `RUN<MODULE>` calls `sassi.modules.base.run_module` in the model directory.

## 9. GUI (`sassi.ui`)

A local web application (decision D-GEN-05 revised: the system Tk is 8.5.9, which is
deprecated on current macOS, so the GUI is browser-based):

* `sassi-gui` (`python -m sassi.ui.server`) starts a stdlib `http.server` on
  `127.0.0.1:<port>` and opens the browser. One `Interpreter` per server session.
* JSON API: `POST /api/command {"line": ...}` (returns messages), `GET /api/model` (nodes,
  elements, groups, interaction flags for plotting), `GET/POST /api/options/<MODULE>` (dialog
  fields <-> option records; POST emits the equivalent command text), `POST /api/run/<MODULE>`
  (background thread, streamed listing via `GET /api/jobs/<id>`), `GET /api/files` and
  `GET /api/file?name=` (results), `GET /api/plot/...` (data for 2D plots).
* Front end: one HTML page with the ACS SASSI menu bar (Model, File, Plot, Modules, Options,
  View, Help), toolbar, tabbed document area (Command History first, plot tabs), the
  Command Entry line, and the Options > Analysis dialog with tabs EQUAKE ... AFWRITE. Plotly.js
  (served from the installed `plotly` package, so it works offline) draws 3D models,
  transfer functions, spectra, histories, soil layer and soil property plots.
* Every GUI action submits command text (L17), so a GUI session can be replayed as a `.pre`.

## 10. Verification and tests

* Unit tests: `tests/unit/test_<package>_*.py`.
* Verification problems: functions in `sassi/verify/problems/vp_NN_<topic>.py` decorated with
  `@problem("VP-NN", ...)` returning `VPResult` with `check(quantity, computed, reference, rtol|atol)`.
  `tests/verification/test_vp_<topic>.py` calls `sassi.verify.run_problem(id, tmp_path)` and
  asserts `.passed`. IDs and reference values are those of requirements §6.3.
* A VP must never be weakened to pass: if a tolerance cannot be met, report the observed error
  and the reason (`VPResult.notes`) and leave it failing (or mark the tolerance *provisional* as
  requirements §6.1 allows, with the observed value recorded).
* Run everything with `.venv/bin/python -m pytest -q`.

## 11. Coding standards

* Python 3.9 syntax (`from __future__ import annotations`, no `match`, no `X | Y` at runtime).
* Docstrings state the equation/decision implemented (e.g. "R1 §3.3 (iv)", "D-MOT-04").
* No global mutable state in modules; deterministic results (D-GEN-09).
* Complex-symmetric systems: never use Cholesky/Hermitian solvers.
* Keep functions vectorised with NumPy; the P0 performance target is 2,000 interaction
  nodes and 20,000 DOFs per frequency on a laptop (D-GEN-10).
