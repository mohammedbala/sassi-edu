# SASSI-EDU Command Reference

> Generated on 2026-10-06 by `python -m sassi.verify.report --commands-only` from the command catalogue of the interpreter (`sassi/prep/registry.py`) and the docstrings of the command handlers (`sassi/prep/commands/`). Do not edit by hand.

Every command of the ACS SASSI V3 manual and of the SASSI-EDU dialect is listed with its syntax, its documented abbreviation or alias, its priority tier and whether this build implements it. The rules of the command language (comma-separated fields, blank fields take defaults, `*` comment lines, abbreviations, variables and loops) are explained in the [User Guide](../user/USER_GUIDE.md#4-the-command-language); the normative detail of every command is in [requirements section 3](../spec/00_requirements.md) and the manual-derived specifications `docs/spec/07` to `11`.

**Syntax notation.** `<arg>` is a required argument, `[arg]` an optional one; a blank field takes the documented default. Names are case-insensitive. Only the full name or the listed abbreviation is accepted (no prefix matching).

**Tiers.** P0 = core linear SSI chain; P1 = important engineering features; P2 = advanced or cosmetic (requirements section 0.1).

**Implemented.** *yes*: the command works in this build. *stored only*: the command is parsed and stored (WRITE writes it back, so `.pre` files of the original program load without loss) but the feature behind it is not available. *no*: the command prints `<CMD> is not available in this build (tier Pn)` or the manual's own message.

## Summary

| Category | Commands | Implemented | Stored only | Not available |
|---|---:|---:|---:|---:|
| [3.4.A Session, files, models and global options](#34a-session-files-models-and-global-options) | 23 | 21 | 0 | 2 |
| [3.4.B Frequency sets](#34b-frequency-sets) | 2 | 2 | 0 | 0 |
| [3.4.C Nodes and coordinate systems](#34c-nodes-and-coordinate-systems) | 18 | 18 | 0 | 0 |
| [3.4.D Groups, elements and properties](#34d-groups-elements-and-properties) | 36 | 36 | 0 | 0 |
| [3.4.E Loads and masses](#34e-loads-and-masses) | 20 | 18 | 0 | 2 |
| [3.4.F Module analysis options](#34f-module-analysis-options) | 38 | 38 | 0 | 0 |
| [3.4.G Module run commands](#34g-module-run-commands) | 12 | 12 | 0 | 0 |
| [3.4.H Model checking](#34h-model-checking) | 7 | 7 | 0 | 0 |
| [3.4.I Model conditioning and generation](#34i-model-conditioning-and-generation) | 24 | 24 | 0 | 0 |
| [3.4.J Cuts, submodels and section calculations](#34j-cuts-submodels-and-section-calculations) | 20 | 19 | 0 | 1 |
| [3.4.K File conversion](#34k-file-conversion) | 5 | 5 | 0 | 0 |
| [3.4.L Plotting and line mathematics](#34l-plotting-and-line-mathematics) | 52 | 52 | 0 | 0 |
| [3.4.M Programming (variables, loops, macros)](#34m-programming-variables-loops-macros) | 13 | 13 | 0 | 0 |
| [3.4.N Water modelling](#34n-water-modelling) | 5 | 5 | 0 | 0 |
| [3.4.O Option NON and nonlinear soil](#34o-option-non-and-nonlinear-soil) | 32 | 29 | 0 | 3 |
| [3.4.P Binary databases](#34p-binary-databases) | 15 | 0 | 1 | 14 |
| [3.4.Q Thick shell](#34q-thick-shell) | 2 | 2 | 0 | 0 |
| [3.4.R Extension commands of SASSI-EDU](#34r-extension-commands-of-sassi-edu) | 41 | 41 | 0 | 0 |
| **Total** | 365 | 342 | 1 | 22 |

| Tier | Commands | Implemented |
|---|---:|---:|
| P0 | 162 | 161 |
| P1 | 125 | 125 |
| P2 | 78 | 56 |

## 3.4.A Session, files, models and global options

| Command | Syntax | Abbreviation / alias | Tier | Implemented | Meaning |
|---|---|---|---|---|---|
| **ACTM** | `ACTM,<Model>` | - | P0 | yes | Activate a model (an empty one is created when absent). |
| **AFWRBAT** | `AFWRBAT,<splits>` | - | P1 | yes | Split the SSI frequency set into &lt;splits> folders with run and combine scripts. |
| **AFWRITE** | `AFWRITE` | AFWR | P0 | yes | CHECK, then write &lt;model>.&lt;ext> for every AOPT-enabled module without errors (D-AFW-01). |
| **AOPT** | `AOPT,<EQUAKE>,<SOIL>,<DEP>,<SITE>,<POINT>,<HOUSE>,<DEP>,<FORCE>,<ANALYS>,<COMBIN>,<MOTION>,<STRESS>,<RELDISP>,<PANEL>` | - | P0 | yes | Modules processed by CHECK and AFWRITE (a non-zero flag includes the module). |
| **CD** | `CD,<dir>` | - | P0 | yes | Change the working directory (it must exist), spec 10 section 5.2. |
| **CHECK** | `CHECK` | CHEC | P0 | yes | Check the model and the AOPT-enabled modules; write &lt;model>.err (manual Chapter 10). |
| **CPMODEL** | `CPMODEL,<Mdl>` | - | P0 | yes | Copy the active model (name and path included) to model &lt;Mdl>. |
| **DMODEL** | `DMODEL,<Mdl>` | - | P0 | yes | Remove a model from memory (files untouched). |
| **GETENV** | `GETENV` | - | P2 | no | Show solver settings -- prints 'GETENV is not available in this build (tier P2)' |
| **GRAVITY** | `GRAVITY,<grav>` | - | P0 | yes | Set HOUSE &lt;gravity> only (other HOUSE options untouched). |
| **GROUNDELEV** | `GROUNDELEV,<elev>` | - | P0 | yes | Set HOUSE &lt;gelev> only. |
| **INP** | `INP,<filename>` | - | P0 | yes | Execute the commands of a .pre file (spec 07 section 9.2.18). |
| **MDL** | `MDL,<Model>,<Path>` | - | P0 | yes | Set the model name and directory; also the working directory. |
| **MDLNAME** | `MDLNAME,<name>` | - | P0 | yes | Change the model name only (path and title unchanged). |
| **MKDIR** | `MKDIR,<dir>` | - | P0 | yes | Create a directory (intermediate ones too); warning if it exists. |
| **MODELLIST** | `MODELLIST` | - | P0 | yes | List the models in memory. |
| **MOPT** | `MOPT,<incomp>,<matrix>,<mass>,<force>` | - | P0 | yes | Model options (record setter; blank = default). |
| **RESUME** | `RESUME` | RESU | P0 | yes | Reload the last SAVE of the active model (D-MDL-02; unknown major versions refused). |
| **SAVE** | `SAVE` | - | P0 | yes | Save the active model to &lt;path>/&lt;model>.sdb (D-MDL-02). |
| **SETENV** | `SETENV,<mem>` | - | P2 | no | Solver memory limit (MB) -- prints 'SETENV is not available in this build (tier P2)' |
| **STATUS** | `STATUS` | STAT | P0 | yes | Global model information (spec 07 section 9.2.37). |
| **TIT** | `TIT,<title>` | - | P0 | yes | Model title (the rest of the line, commas included). |
| **WRITE** | `WRITE,[<file>],[<path>]` | WRIT | P0 | yes | Write the model as commands (default &lt;model>.pre in the model path). |

## 3.4.B Frequency sets

| Command | Syntax | Abbreviation / alias | Tier | Implemented | Meaning |
|---|---|---|---|---|---|
| **FREQ** | `FREQ,<ndx>,<f1>,...,<f10>` | - | P0 | yes | Append frequency numbers to set ndx; f1 = 0 deletes the set. |
| **LFREQ** | `LFREQ,[start],[end],[step]` | LFRE | P0 | yes | List frequency sets with numbers and Hz values (current df). |

## 3.4.C Nodes and coordinate systems

| Command | Syntax | Abbreviation / alias | Tier | Implemented | Meaning |
|---|---|---|---|---|---|
| **CSYS** | `CSYS,<ns>` | - | P0 | yes | Activate a coordinate system (0 = global). |
| **D** | `D,<n1>,<n2>,[<inc>],[<val>],<label1>,...,<label6>` | - | P0 | yes | Set fixity codes (0 free, default; 1 fixed). |
| **FILL** | `FILL,[<n1>],[<n2>],[<nr>]` | - | P0 | yes | Nr equally spaced nodes on the line n1-n2 (numbers n1 + k\*d). |
| **GLOBAL** | `GLOBAL,<n1>,<n2>,<inc>` | GLOB | P0 | yes | Convert stored coordinates to global (the active system is unchanged). |
| **INT** | `INT,<n1>,<n2>,[<inc>],<set>,[<code>]` | - | P0 | yes | Set (1) / reset (0) node flag code 0 interaction (default), 1 intermediate, 2 interface, 3 internal (D-MDL-07). |
| **INTLIST** | `INTLIST,[<n1>],[<n2>],[<step>],[<c1>],[<c2>],[<c3>],[<c4>]` | INTL | P0 | yes | List classified nodes. |
| **LMOVE** | `LMOVE,[<dx>],[<dy>],[<dz>],<nd>,<l1>,...,<l15>` | LMOV | P0 | yes | Node nd+i-1 = node l_i + (dx, dy, dz). |
| **LOC** | `LOC,<ns>,<type>,<x0>,<y0>,<z0>,<txy>,<tyz>,<txz>` | - | P0 | yes | Cartesian system, R = Rz(txy) Rx(tyz) Ry(txz). |
| **LOCAL** | `LOCAL,<ns>,<type>,<n1>,<n2>,<n3>` | LOCA | P0 | yes | System with origin n1, x toward n2, n3 in the +y half plane. |
| **N** | `N,<nd>,[<x>],[<y>],[<z>]` | - | P0 | yes | Define a node in the active coordinate system (legacy fields 5-6 ignored). |
| **NDEL** | `NDEL,<n1>,[<n2>],[<inc>]` | - | P0 | yes | Delete nodes with their loads, masses, fixities and flags (D-MDL-11). |
| **NGEN** | `NGEN,[itim],[step],[n1],[n2],[inc],[dx],[dy],[dz]` | - | P0 | yes | Copy a node pattern itim times. |
| **NLIST** | `NLIST,[<n1>],[<n2>],[<inc>]` | NLIS | P0 | yes | List nodes (stored coordinates, system, fixities, INT flags I/M/F/N). |
| **NMED** | `NMED,<nd>,<n1>,[<n2>,...,<n8>]` | - | P0 | yes | Node at the mean of the global coordinates of 1 to 8 nodes. |
| **NMOVE** | `NMOVE,[<dx>],[<dy>],[<dz>],<nd>,<l1>,...,<l15>` | NMOV | P0 | yes | Node nd+i-1 = (x dx, y dy, z dz) of node l_i (defaults 1). |
| **NSCALE** | `NSCALE,<n1>,<n2>,[<inc>],[<sfx>],[<sfy>],[<sfz>]` | NSCA | P0 | yes | Scale coordinates in place (factor 0 -> 1). |
| **SDEL** | `SDEL,<s1>,[<s2>],[<inc>]` | - | P0 | yes | Delete systems; their nodes are converted to global first. |
| **SLIST** | `SLIST,[<s1>],[<s2>],[<inc>]` | SLIS | P0 | yes | List coordinate systems. |

## 3.4.D Groups, elements and properties

| Command | Syntax | Abbreviation / alias | Tier | Implemented | Meaning |
|---|---|---|---|---|---|
| **DELL** | `DELL,<m1>,[<m2>],[<inc>]` | - | P0 | yes | Delete soil layers. |
| **DELM** | `DELM,<m1>,[<m2>],[<inc>]` | - | P0 | yes | Delete materials. |
| **DELR** | `DELR,<r1>,[<r2>],[<inc>]` | - | P0 | yes | Delete beam sections. |
| **DELSC** | `DELSC,<r1>,[<r2>],[<inc>]` | DELS | P0 | yes | Delete spring properties. |
| **E** | `E,<ne>,<n1>,...,<n8>` | - | P0 | yes | Element of the active group (node count by group type). |
| **ECOMPR** | `ECOMPR` | ECOM | P0 | yes | Renumber the elements of the active group 1..N (order kept; EOUT remapped). |
| **EDEL** | `EDEL,<e1>,[<e2>],[<inc>]` | - | P0 | yes | Delete elements of the active group (compress with ECOMPR). |
| **EGEN** | `EGEN,[<itim>],<ninc1>,<e1>,[<e2>],[<inc>],[<ee>]` | - | P0 | yes | Copy an element pattern itim times. |
| **EINT** | `EINT,<e1>,<e2>,[<inc>],<order>` | - | P0 | yes | SOLID 0/1/2 (2, 3, 4 Gauss points); TSHELL 0 reduced / 1 selective. |
| **ELIST** | `ELIST,[<e1>],[<e2>],[<inc>]` | ELIS | P0 | yes | List the elements of the active group. |
| **ETYPE** | `ETYPE,<e1>,<e2>,[<inc>],<type>` | ETYP | P0 | yes | 0 implicit (by ground elevation), 1 structure, 2 excavated soil (SOLID/PLANE) or embedded/buried shell (SHELL/TSHELL). |
| **GDEL** | `GDEL,<g1>,[<g2>],[<inc>]` | - | P0 | yes | Delete groups and their elements (their EOUT requests are dropped). |
| **GLIST** | `GLIST,[<g1>],[<g2>],[<inc>]` | GLIS | P0 | yes | List groups. |
| **GROUP** | `GROUP,<ng>,<type>` | GROU | P0 | yes | Create a group of the given type, or activate an existing group. |
| **GTIT** | `GTIT,[gr],<title>` | - | P0 | yes | Group title (rest of the line); without a group number the active group. |
| **KI** | `KI,<e1>,[<e2>],[<inc>],<k1>,...,<k6>` | - | P0 | yes | Beam end releases at node I (1 = released). |
| **KJ** | `KJ,<e1>,[<e2>],[<inc>],<k1>,...,<k6>` | - | P0 | yes | Beam end releases at node J (1 = released). |
| **L** | `L,<nm>,<thick>,<weight>,<pveloc>,<sveloc>,<pdamp>,<sdamp>` | - | P0 | yes | Soil layer (SITE and excavated soil). |
| **LLIST** | `LLIST,<m1>,[<m2>],[<step>]` | LLIS | P0 | yes | List soil layers with derived G, M, rho. |
| **M** | `M,<nm>,<val1>,<val2>,<weight>,<pdamp>,<sdamp>,<type>` | - | P0 | yes | Material; type 1 (E, nu), 2 (M, G), 3 (Vp, Vs). |
| **MACT** | `MACT,<index>` | - | P0 | yes | Active material / soil-layer index for new elements (all groups). |
| **MLIST** | `MLIST,<m1>,[<m2>],[<step>]` | MLIS | P0 | yes | List materials with derived E, nu, G, M, Vp, Vs, rho. |
| **MSET** | `MSET,<e1>,[<e2>],[<inc>],<index>` | - | P0 | yes | Material (M) index, or soil-layer (L) index for excavated SOLID/PLANE elements (interpreted at CHECK/AFWRITE, D-MDL-10). |
| **MTYPE** | `MTYPE,[<gr>],<type>` | MTYP | P0 | yes | Change the type of a group (default the active group). |
| **MXDEL** | `MXDEL,<p1>,[<p2>],[<step>]` | MXDE | P0 | yes | Delete matrix properties. |
| **MXI** | `MXI,<p>,<row>,<t1>,...,<t12>` | - | P0 | yes | Imaginary stiffness row of GENERAL matrix property p. |
| **MXLIST** | `MXLIST,<p>` | MXLI | P0 | yes | List the real stiffness, imaginary stiffness and mass/weight matrices of property p. |
| **MXM** | `MXM,<p>,<row>,<t1>,...,<t12>` | - | P0 | yes | Mass (or weight, MOPT &lt;matrix> = 1) row of GENERAL matrix property p. |
| **MXR** | `MXR,<p>,<row>,<t1>,...,<t12>` | - | P0 | yes | Real stiffness row of GENERAL matrix property p (upper triangle). |
| **R** | `R,<nm>,<axial>,<shear2>,<shear3>,<tors>,<flex2>,<flex3>` | - | P0 | yes | Beam section (A, As2, As3, J, I2, I3). |
| **RACT** | `RACT,<index>` | - | P0 | yes | Active real / spring / matrix property index for new elements (all groups). |
| **RLIST** | `RLIST,<r1>,[<r2>],[<step>]` | RLIS | P0 | yes | List beam sections. |
| **RSET** | `RSET,<e1>,[<e2>],[<inc>],<index>` | - | P0 | yes | Property index (R for BEAMS, SC for SPRING, MX for GENERAL). |
| **SC** | `SC,<nm>,<scx>,<scy>,<scz>,<scxx>,<scyy>,<sczz>,<damp>` | - | P0 | yes | Spring constants (global axes) and damping. |
| **SCLIST** | `SCLIST,<r1>,[<r2>],[<step>]` | SCLI | P0 | yes | List spring properties. |
| **THICK** | `THICK,<e1>,<e2>,[<inc>],<thick>` | THIC | P0 | yes | Shell thickness (SHELL and TSHELL groups). |

## 3.4.E Loads and masses

| Command | Syntax | Abbreviation / alias | Tier | Implemented | Meaning |
|---|---|---|---|---|---|
| **F** | `F,<n>,<fx>,<fy>,<fz>,<tx>,<ty>,<tz>` | - | P0 | yes | Force factors and arrival times (global X, Y, Z). |
| **FDEL** | `FDEL,<n1>,[<n2>],[<inc>]` | - | P0 | yes | Delete forces. |
| **FLIST** | `FLIST,[<n1>],[<n2>],[<inc>]` | FLIS | P0 | yes | List forces. |
| **FREAD** | `FREAD` | - | P2 | no | FREAD is not available (syntax not documented, D-FRC-04) -- FREAD is not available (syntax not documented, D-FRC-04) |
| **FSCALE** | `FSCALE,[<n1>],[<n2>],[<inc>],[<sx>],[<sy>],[<sz>]` | FSCA | P0 | yes | Scale force factors (0 -> 1; arrival times kept). |
| **MM** | `MM,<n>,<fxx>,<fyy>,<fzz>,<txx>,<tyy>,<tzz>` | - | P0 | yes | Moment factors and arrival times (about X, Y, Z). |
| **MMDEL** | `MMDEL,<n1>,[<n2>],[<inc>]` | MMDE | P0 | yes | Delete moments. |
| **MMLIST** | `MMLIST,[<n1>],[<n2>],[<inc>]` | MMLI | P0 | yes | List moments. |
| **MR** | `MR,<n>,<mxx>,<myy>,<mzz>` | - | P0 | yes | Rotational masses about global X, Y, Z (units by MUNITS). |
| **MRDEL** | `MRDEL,<n1>,[<n2>],[<inc>]` | MRDE | P0 | yes | Delete rotational masses. |
| **MREAD** | `MREAD` | - | P2 | no | MREAD is not available (syntax not documented, D-FRC-04) -- MREAD is not available (syntax not documented, D-FRC-04) |
| **MRGEN** | `MRGEN,[<itim>],[<ninc>],<n1>,<n2>,[<inc>],[<mxx>],[<myy>],[<mzz>]` | MRGE | P0 | yes | Generate rotational masses. |
| **MRSCALE** | `MRSCALE,[<n1>],[<n2>],[<inc>],[<sx>],[<sy>],[<sz>]` | MRSC | P0 | yes | Scale rotational masses (0 -> 1). |
| **MSCALE** | `MSCALE,[<n1>],[<n2>],[<inc>],[<sx>],[<sy>],[<sz>]` | MSCA | P0 | yes | Scale moment factors (0 -> 1). |
| **MT** | `MT,<n>,<mx>,<my>,<mz>` | - | P0 | yes | Translational masses in global X, Y, Z (units by MUNITS). |
| **MTDEL** | `MTDEL,<n1>,[<n2>],[<inc>]` | MTDE | P0 | yes | Delete translational masses. |
| **MTGEN** | `MTGEN,[<itim>],[<ninc>],<n1>,<n2>,[<inc>],[<mx>],[<my>],[<mz>]` | MTGE | P0 | yes | Generate translational masses. |
| **MTLIST** | `MTLIST,[<n1>],[<n2>],[<inc>]` | MTLI | P0 | yes | List translational and rotational masses with units and mass values. |
| **MTSCALE** | `MTSCALE,[<n1>],[<n2>],[<inc>],[<sx>],[<sy>],[<sz>]` | MTSC | P0 | yes | Scale translational masses (0 -> 1). |
| **MUNITS** | `MUNITS,<n1>,[<n2>],[<step>],<units>` | MUNI | P0 | yes | Nodal mass units, 0 mass / 1 weight (default 1, D-MDL-08). |

## 3.4.F Module analysis options

| Command | Syntax | Abbreviation / alias | Tier | Implemented | Meaning |
|---|---|---|---|---|---|
| **ACCIN** | `ACCIN,<no>,<file>` | ACCI | P0 | yes | Seed / external acceleration input file of spectrum &lt;no>. |
| **ACCOUT** | `ACCOUT,<no>,<file>` | ACCO | P0 | yes | Generated acceleration output file of spectrum &lt;no>. |
| **AMP** | `AMP,<no>,<a1>,...,<a100>` | - | P1 | yes | Append spectral amplification ratios of motion &lt;no>; a1 = 0 clears. |
| **ANALYS** | `ANALYS,<opmode>,<type>,<mode>,<save>,<prnt>,<fopt>,<ang>,<xc>,<yc>,<zc>,<impe>,[simul]` | ANAL | P0 | yes | ANALYS options |
| **CORR** | `CORR,<no>,<time>,<val>` | - | P1 | yes | X-Y correlation pair &lt;no> (P1). |
| **DAMP** | `DAMP,<d1>,...,<d10>` | - | P0 | yes | Append RS damping ratios (fractions); d1 = 0 clears the list. |
| **DYNP** | `DYNP,<no>,<sg>,<g>,<sd>,<d>,<label>` | - | P0 | yes | Curve point (strain %, G/Gmax; strain %, damping %). |
| **EOUT** | `EOUT,<code1>,...,<code12>,<group>,<element list>` | - | P0 | yes | Element output request; EOUT,0 clears. |
| **EQTIT** | `EQTIT,<title>` | EQTI | P0 | yes | Spectra title (rest of the line). |
| **EQUAKE** | `EQUAKE,<accopt>,<nrfreq>,<rand>,<damp>,<dur>,<corr>,<seeds>,[tpsd]` | EQUA | P0 | yes | EQUAKE options (D-EQK-01). |
| **FORCE** | `FORCE,<opmode>` | FORC | P0 | yes | FORCE options (gravity, df, delt, NFFT and the frequency set are shared). |
| **HOUSE** | `HOUSE,<gravity>,<gelev>,<opmode>,<dim>,<imp>,<coh>,<wpass>,<me>,<cmplxspec>` | HOUS | P0 | yes | HOUSE options (&lt;coh> 1 incoherent motion, &lt;wpass> 1 wave passage, &lt;me> 1 multiple excitation). |
| **INCOH** | `INCOH,<gammax>,<gammay>,<gammaz>,<alpha>,<ngp>,<ipr>,<nmodes>,<met>,<HSeed>,<VSeed>,<RandPhz>` | INCO | P1 | yes | Incoherency options of HOUSE (requirements 4.4 item 6). |
| **ME** | `ME,<no>,<nfirst>,<nlast>,<xc>,<yc>,<zc>` | - | P1 | yes | Multiple-excitation zone &lt;no> = the interaction nodes &lt;nfirst>..&lt;nlast> (control point stored, not used). |
| **MOTION** | `MOTION,<opmode>,<out>,<step>,<dur>,<res>,<freq1>,<freq2>,<fstep>,<mult>,<max>,<rec1>,<rec2>,<fopt>,<bl>,<smo>,<cplx>,<cnvrt>,<pzadj>,<interp>` | MOTI | P0 | yes | MOTION options (19 arguments). |
| **NOUT** | `NOUT,<dir>,<code1>,...,<code6>,<node list>` | - | P0 | yes | Nodal output request; NOUT,0 clears the list. |
| **POINT** | `POINT,<opmode>,<layer>,<rad>` | POIN | P0 | yes | Embedment layers (0 surface) and central-zone radius. |
| **RDND** | `RDND,<NodeNum>,<X>,<Y>,<Z>,<XX>,<YY>,<ZZ>` | - | P0 | yes | RELDISP node request (flag >= 1 = on); RDND,0 clears. |
| **RELD** | `RELD,<RelDisOutput>,<RelDispSAll>,<RelDispNumFiles>` | - | P0 | yes | RELDISP options (D-RDP-03). |
| **RELFILE** | `RELFILE,<FileName>` | - | P0 | yes | Reference-node complex .TFI file of RELDISP. |
| **RSIN** | `RSIN,<no>,<file>` | - | P0 | yes | Target response-spectrum file of spectrum &lt;no> (blank file clears it). |
| **RSOUT** | `RSOUT,<no>,<file>` | RSOU | P0 | yes | Response-spectrum output file of spectrum &lt;no>. |
| **SACC** | `SACC,<layer>,<opt>,<outcrop>` | - | P0 | yes | Acceleration output at the top of a SOIL sublayer. |
| **SFOU** | `SFOU,<layer>,<out>,<save>,<outcrop>,<smooth>,<nrval>` | - | P0 | yes | Parsed and stored only. |
| **SITE** | `SITE,<opmode>,<mode1>,<fstep>,<nl>,<hs>,<mode2>,<wopt>,<freq1>,<freq2>,<cl>,<cm>,<delt>,<nft>,<freq>` | - | P0 | yes | SITE options; owns Δt, NFFT, Δf, frequency set, control layer |
| **SOIL** | `SOIL,<nrval>,<grav>,<header>,<outcrop>,<save>,<iter>,<ratio>,<gravmult>,<cof>` | - | P0 | yes | SOIL options. |
| **SPRO** | `SPRO,<layer>,<prop>,<dynprop>` | - | P0 | yes | SOIL sublayer -> L property and DYNP label (last = half-space). |
| **SRS** | `SRS,<layer>,<save>,<outcrop>` | - | P0 | yes | Response-spectrum output at the top of a SOIL sublayer. |
| **SSAF** | `SSAF,<layer>,<save>,<outcrop1>,<outcrop2>,<layer2>,<freqstep>,<title>` | - | P0 | yes | Spectral amplification. |
| **SSTR** | `SSTR,<layer>,<opt1>,<opt2>,<opt3>,<opt4>` | - | P0 | yes | Stress/strain output at the centre of a sublayer. |
| **STRESS** | `STRESS,<opmode>,<iter>,<save>,<itran>,<interopt>` | STRE | P0 | yes | STRESS options. |
| **SYMM** | `SYMM,<no>,[<type>],[<node1>],[<node2>],[<node3>]` | - | P1 | yes | Symmetry plane; &lt;node1> = 0 resets plane &lt;no>. |
| **THFILE** | `THFILE,<file>` | THFI | P0 | yes | Control-motion acceleration file (shared by SOIL, MOTION, STRESS, RELDISP). |
| **THTIT** | `THTIT,<title>` | THTI | P0 | yes | Control-motion title. |
| **TOPL** | `TOPL,<l1>,...` | - | P0 | yes | Append L-layer numbers (top first) to the free-field list; l1 = 0 clears it. |
| **TPSD** | `TPSD,<num>,<file>` | - | P0 | yes | Target PSD file of spectrum &lt;num>. |
| **WAVE** | `WAVE,<type>,<opt>,<ratio1>,<ratio2>,<angle>` | - | P0 | yes | Wave field (1 R, 2 SV, 3 P, 4 SH, 5 L). |
| **WPASS** | `WPASS,<appv>,<ang>,<cohf>` | WPAS | P1 | yes | Apparent velocity and angle of Line D; &lt;cohf> = unlagged coherency model 1-7 (D-INC-01). |

## 3.4.G Module run commands

| Command | Syntax | Abbreviation / alias | Tier | Implemented | Meaning |
|---|---|---|---|---|---|
| **RUNANALYS** | `RUNANALYS,[model]` | - | P0 | yes | Run module ANALYS in the model directory (requirements 2.6). |
| **RUNCOMBIN** | `RUNCOMBIN,[model]` | - | P0 | yes | Run module COMBIN in the model directory (requirements 2.6). |
| **RUNEQUAKE** | `RUNEQUAKE,[model]` | - | P0 | yes | Run module EQUAKE in the model directory (requirements 2.6). |
| **RUNFORCE** | `RUNFORCE,[model]` | - | P0 | yes | Run module FORCE in the model directory (requirements 2.6). |
| **RUNHOUSE** | `RUNHOUSE,[model]` | - | P0 | yes | Run module HOUSE in the model directory (requirements 2.6). |
| **RUNMOTION** | `RUNMOTION,[model]` | - | P0 | yes | Run module MOTION in the model directory (requirements 2.6). |
| **RUNNONLINEAR** | `RUNNONLINEAR,[model]` | - | P2 | yes | Write &lt;model>.eql from the EQL / P / S / BBC data (Errors 121-128 stop the run) and run the NONLINEAR module in the model directory (it reads &lt;model>.hou and the .THD files of RELDISP). |
| **RUNPOINT** | `RUNPOINT,[model]` | - | P0 | yes | Run module POINT in the model directory (requirements 2.6). |
| **RUNRELDISP** | `RUNRELDISP,[model]` | - | P0 | yes | Run module RELDISP in the model directory (requirements 2.6). |
| **RUNSITE** | `RUNSITE,[model]` | - | P0 | yes | Run module SITE in the model directory (requirements 2.6). |
| **RUNSOIL** | `RUNSOIL,[model]` | - | P0 | yes | Run module SOIL in the model directory (requirements 2.6). |
| **RUNSTRESS** | `RUNSTRESS,[model]` | - | P0 | yes | Run module STRESS in the model directory (requirements 2.6). |

## 3.4.H Model checking

| Command | Syntax | Abbreviation / alias | Tier | Implemented | Meaning |
|---|---|---|---|---|---|
| **EXCSTRCHK** | `EXCSTRCHK` | - | P0 | yes | Excavation interior nodes shared with structure, beam, spring or GM elements (G-15). |
| **FIXEDINT** | `FIXEDINT` | - | P0 | yes | Interaction nodes with any fixed translational DOF. |
| **FREESPRING** | `FREESPRING` | - | P0 | yes | Unconstrained nodes connected only to springs. |
| **HINGED** | `HINGED` | - | P0 | yes | Possible unintended hinges (6-DOF elements meeting SOLID/PLANE at one node; drilling joints). |
| **INTCOUNT** | `INTCOUNT` | - | P0 | yes | Number of interaction nodes and the ANALYS memory estimate (D-ANL-10, UT-18). |
| **KINT** | `KINT` | - | P0 | yes | Beam K nodes that are interaction nodes. |
| **USED** | `USED` | - | P0 | yes | Fix all DOFs of the nodes not used by any element (D,n,n,1,1,ALL); unused interaction nodes are reported, not fixed (fixing them would trip FIXEDINT). |

## 3.4.I Model conditioning and generation

| Command | Syntax | Abbreviation / alias | Tier | Implemented | Meaning |
|---|---|---|---|---|---|
| **CRITFREQ** | `CRITFREQ,<tol>,<minfilter>,<TF>,<Var>` | - | P1 | yes | Frequencies where interpolated TF peaks deviate from TFU. |
| **ETYPEGEN** | `ETYPEGEN,<type>,[RESET]` | - | P1 | yes | ETYPE of every SOLID/PLANE/SHELL/TSHELL element: 0 by location (explicit 1/2), 1 structure, 2 excavated soil / buried shell; ``0,RESET`` restores the implicit ETYPE 0. |
| **EXCAV** | `EXCAV,<model>,[delta]` | - | P1 | yes | Excavation volume of the active model's basement, stored in &lt;model>. |
| **FIXROT** | `FIXROT,[Stiff]` | - | P0 | yes | D on solid-only rotations and axis-parallel shell drilling rotations; soft springs (default 10) for oblique shells; FIXSPRROT for spring nodes. |
| **FIXSHLROT** | `FIXSHLROT,[stiff]` | - | P0 | yes | Soft drilling springs (default 10) at coplanar-shell nodes. |
| **FIXSLDROT** | `FIXSLDROT` | - | P0 | yes | Fix the rotations of nodes connected only to SOLID elements. |
| **FIXSPRROT** | `FIXSPRROT` | FIXSPROT | P0 | yes | Fix the unstiffened DOFs of spring-only and spring+solid nodes. |
| **FRAMECOMBIN** | `FRAMECOMBIN,<op>,<num>,<InFile1>,...,<InFileX>,<Outfile>` | - | P1 | yes | Combine frames (0 SRSS, 1 sum, 2 average). |
| **FRAMESEL** | `FRAMESEL,<tol>,<Acc>,<Var>` | - | P1 | yes | Critical frame numbers (local extrema >= tol % of max\|a\|) -> variable. |
| **GCOM** | `GCOM` | - | P1 | yes | Renumber groups 1..G without gaps (order kept); element numbers unchanged. |
| **GLB2LOC** | `GLB2LOC,<Start>,<End>,<Stride>,<Sysno>` | - | P1 | yes | Store nodes in Cartesian system Sysno (0 = global). |
| **GROUPMAT** | `GROUPMAT` | - | P2 | yes | One new material per group (copy of the material of its first element); the original material table is replaced (Option NON preparation). |
| **INTGEN** | `INTGEN,<type>,[level skip]` | - | P1 | yes | Generate interaction nodes: 0 clear, 1 FV, 2 FI-EVBN (MSM), 3 FI-FSIN (SM), 4 surface, 5 FFV (internal levels every level skip + 1, default skip 1); 1-5 add. |
| **MERGE** | `MERGE,<Mdl1>,<Mdl2>,<X>,<Y>,<Z>` | - | P1 | yes | Active model = Mdl1 + Mdl2 (numbers offset, translated by X,Y,Z). |
| **MERGEGROUP** | `MERGEGROUP,<dest>,[G1],...,[G10]` | - | P1 | yes | Merge same-type groups into dest (then GCOM). |
| **MERGESOIL** | `MERGESOIL,<Struct>,<Soil>,[Mode],[Stiff],[Stiff2],[SepLevel],[Mapping]` | - | P1 | yes | Join a structure model and an excavation model into the active model (Mode 0 unbonded, 1 merged nodes, 2 stiff springs, 3 stiff below / soft above SepLevel). |
| **MODFRAMES** | `MODFRAMES,<cols>,<framelist>` | - | P1 | yes | Set the column count (second header number) of every listed frame. |
| **NCOM** | `NCOM,[MapFile]` | - | P1 | yes | Renumber nodes 1..N without gaps (order kept) and remap every node reference. |
| **RADIUS** | `RADIUS,[Scale],[FileName]` | - | P1 | yes | POINT central-zone radius r_e = Scale sqrt(A_plan) of every excavation element, with min / average / max (D-PNT-05; Scale default 0.9). |
| **RMVUNUSED** | `RMVUNUSED` | - | P1 | yes | Delete the nodes that are used by no element and are not interaction nodes. |
| **ROTATE** | `ROTATE,<x>,<y>,<z>,<rxy>,<ryz>,<rzx>` | - | P1 | yes | Rotate the model about (x,y,z): rxy about Z, then ryz about X, then rzx about Y (degrees, right-hand rule, D-MDL-14). |
| **SOILMESH** | `SOILMESH,<dest>,<sX>,<sY>,<hori>,<vert>,<xAdj>,<yAdj>,<Zdepth>,<contact>,<rNum>` | - | P1 | yes | Near-field soil mesh around the basement (soil-pressure models), stored in &lt;dest> numbered after the active model. |
| **TRANSLATE** | `TRANSLATE,<x>,<y>,<z>` | - | P1 | yes | Move every node (defaults 0). |
| **WELD** | `WELD,[FORCE]` | - | P1 | yes | Merge coincident nodes in the connectivity (lowest number kept); the two ends of a zero-length spring are never welded together unless FORCE (other nodes at that point still are). |

## 3.4.J Cuts, submodels and section calculations

| Command | Syntax | Abbreviation / alias | Tier | Implemented | Meaning |
|---|---|---|---|---|---|
| **CALCC** | `CALCC` | - | P1 | yes | Volume centroid of the active model (SOLID, SHELL, PLANE, BEAMS elements; on a cross-section model the BEAMS, which are not section pieces, are left out so that it equals the area centroid). |
| **CALCM** | `CALCM` | - | P1 | yes | Element mass (material weight x volume / g) and lumped MT/MR masses of the active model. |
| **CALCMOI** | `CALCMOI,<nx>,<ny>,<nz>,<rx>,<ry>,<rz>,<sysno>` | - | P1 | yes | Section moments of inertia of the active cross-section model (on any other model: the base section of the whole structure, Q-11) about its area centroid, in the local cut axes. |
| **CALCPAR** | `CALCPAR,<nx>,<ny>,<nz>,<rx>,<ry>,<rz>,<sysno>,[verbose]` | - | P1 | yes | Area, centroid, moments of inertia and the six section resultants (D-SEC-02 order) of the active cross-section model (on any other model: the base section of the whole structure, Q-11). |
| **CALCSECTHIST** | `CALCSECTHIST,<infile>,<cutnum>,<px>,<py>,<pz>,<nx>,<ny>,<nz>,<rx>,<ry>,<rz>,<sysno>,[ts],<outfile>` | - | P1 | yes | Section resultants of a cut of the active (original) model for every .ess frame of a list file (7-column CSV with a final MAX row). |
| **CALCSECTHISTDB** | `CALCSECTHISTDB,<cutnum>,<px>,<py>,<pz>,<nx>,<ny>,<nz>,<rx>,<ry>,<rz>,<sys>,<ts>,[start],[end],[Stride],<outfile>` | - | P2 | no | Same from the binary stress DB -- prints 'CALCSECTHISTDB is not available in this build (tier P2)' |
| **CSECT** | `CSECT,<dest>,<cutnum>,<px>,<py>,<pz>,<nx>,<ny>,<nz>` | - | P1 | yes | Cross-section model (unit thickness along n, same group and element numbers) of the cut's elements in model &lt;dest>. |
| **CUT2SUB** | `CUT2SUB,<cutnum>,<dest>,[solid]` | - | P1 | yes | Submodel of the cut's elements in model &lt;dest> (same group and element numbers); solid >= 1 turns shells into thick-shell solids for plotting. |
| **CUTADD** | `CUTADD,<cutnum>,<group>,<e1>,...,<eN> \| CUTADD,<cutnum>,<group>,RANGE,<start>,[end],[stride]` | - | P1 | yes | Add elements of one group to a cut (created if needed; duplicates ignored). |
| **CUTCLR** | `CUTCLR,<first>,[last],[step]` | - | P1 | yes | Delete cuts (they can be rebuilt by adding elements again). |
| **CUTRMV** | `CUTRMV,<cutnum>,<group>,<e1>,...,<eN> \| CUTRMV,<cutnum>,<group>,RANGE,<start>,[end],[stride]` | - | P1 | yes | Remove elements from a cut. |
| **CUTVOL** | `CUTVOL,<cutnum>,[Xmin],[Xmax],[Ymin],[Ymax],[Zmin],[Zmax]` | - | P1 | yes | Add the elements whose nodes all lie in the closed box (blank bounds = the active model's extents; D-MDL-13). |
| **EXTRACTEXCAV** | `EXTRACTEXCAV,<Model>` | - | P1 | yes | Submodel of the excavation elements (explicit ETYPE 2 SOLID/PLANE) with their nodes, interaction flags and soil layers; the active model is unchanged. |
| **READSTR** | `READSTR,<Filename>,[Dir]` | - | P1 | yes | Attach the element stresses of an .ess file to the active model. |
| **SECDATAOPT** | `SECDATAOPT,<flag>` | - | P1 | yes | STRESS writes NSTRESS/ESTRESS\_&lt;n>.ess frames for the whole history (1) or not (0). |
| **SHEAR** | `SHEAR,<panel>,[fc],[fy],[P],[Nu],[Abe],[fybe]` | - | P2 | yes | Shear capacities of wall panel &lt;panel> (0 = all) by ACI 318-08, Wood 1990 (upper/lower bound), Barda 1977 and Gulec-Whittaker 2009; fc, fy in ksi (kN/m2), P = web reinforcement ratio, Nu axial force (kips, kN), Abe/fybe ... |
| **SLICE** | `SLICE,<cutnum>,<px>,<py>,<pz>,<nx>,<ny>,<nz>` | - | P1 | yes | Add the elements that cross (or touch) the plane. |
| **SPLITGROUP** | `SPLITGROUP,<group>,<split>,[dir]` | - | P1 | yes | Split a group by the plane of the first shell of SHELL group &lt;split>, or (dir = X, Y, Z) by the plane &lt;axis> = &lt;split>; elements on the + side move to a new group. |
| **TRANELEM** | `TRANELEM,<dest>,<group>,<begin>,<end>,<stride>` | - | P1 | yes | Copy elements of a group into model &lt;dest> (overwrites elements, group type and nodes there). |
| **TRANVOL** | `TRANVOL,<dest>,[Xmin],[Xmax],[Ymin],[Ymax],[Zmin],[Zmax]` | - | P1 | yes | Copy the elements inside the box (CUTVOL test) into model &lt;dest>. |

## 3.4.K File conversion

| Command | Syntax | Abbreviation / alias | Tier | Implemented | Meaning |
|---|---|---|---|---|---|
| **ANSYS** | `ANSYS,[FileName],[Dir],[<dmap>]` | - | P1 | yes | Write the active model as ANSYS APDL input (&lt;model>.inp). |
| **ANSYSMODELTYPE** | `ANSYSMODELTYPE,<type>` | - | P2 | yes | Option AA operation mode 1 embedded, 2 surface (P2: stored only). |
| **ANSYSREFORMAT** | `ANSYSREFORMAT,<Org>,<Map>` | - | P1 | yes | Copy model &lt;Org> into the empty active model with one BEAMS group per release pattern (ANSYS KEYOPT(7)/(8) are per element type) and write the group map. |
| **CONVERT** | `CONVERT,<ConSel>,<model>,<filename>,...` | - | P1 | yes | Run a file converter (SSI, ANSYS; STRUDL not applicable). |
| **GENMATRIXDAMP** | `GENMATRIXDAMP,<begin>,[end],[stride],<damp>` | - | P2 | yes | Not usable in this version (manual). |

## 3.4.L Plotting and line mathematics

| Command | Syntax | Abbreviation / alias | Tier | Implemented | Meaning |
|---|---|---|---|---|---|
| **ADDITION** | `ADDITION,<dest>,<source1>,...,<source100>` | - | P1 | yes | Sum of lines ('Linear Combin.'). |
| **AVERAGE** | `AVERAGE,<dest>,<source1>,...,<source100>` | - | P1 | yes | Mean of lines ('Average Line'). |
| **AXES** | `AXES,<MaxTickX>,<MaxTickY>,<MinTickX>,<MinTickY>,<LogX>,<LogY>` | - | P1 | yes | Grids and log axes (0/1; blank unchanged). |
| **BROADEN** | `BROADEN,<Dest>,<Smooth1>,<Smooth2>,<source1>,...` | - | P1 | yes | Envelope, +-Smooth2 % peak broadening, Smooth1 % bridging. |
| **BUBBLEPLOT** | `BUBBLEPLOT,<BufferDir>,<MnF>,<MxF>,<ST>,<MnR>,<MxR>,<Col>` | - | P1 | yes | Bubble animation of nodal values. |
| **CAPTUREPLOT** | `CAPTUREPLOT,<FileName>` | - | P1 | yes | Save the active plot (PNG if '.png' is in the name, else BMP). |
| **CLOSEPLOT** | `CLOSEPLOT` | - | P1 | yes | Close the active plot (the most recent remaining plot becomes active). |
| **CNGCENTER** | `CNGCENTER,<X>,<Y>,<Z>` | - | P1 | yes | Centre of rotation of the 3D plot (global coordinates). |
| **CNGVIEW** | `CNGVIEW,<rX>,<rY>,<rZ>,<px>,<py>,<zoom>` | - | P1 | yes | 3D view (degrees, pan in model units, zoom 1 = fit). |
| **COLOR** | `COLOR,<Palette>,<Num>,<R>,<G>,<B>` | - | P1 | yes | Set colour Num (1-based) of a palette (0-255 values). |
| **CONTOURPLOT** | `CONTOURPLOT,<BufferDir>,<MnF>,<MxF>,<ST>,<MnR>,<MxR>,<Col>` | - | P1 | yes | Contour animation on element faces. |
| **CUTPLOT** | `CUTPLOT,[Cut],[Model]` | - | P1 | yes | Wireframe of a model with the elements of a cut filled (any model). |
| **DEBUG** | `DEBUG,[switch]` | - | P1 | yes | View values / animation information on 3D plots (0 off, 1 on, 2 toggle). |
| **DEFORMPLOT** | `DEFORMPLOT,<BufferDir>,<MnF>,<MxF>,<ST>,<Scale>` | - | P1 | yes | Animated deformed shape x' = x + Scale u. |
| **ELECOLOR** | `ELECOLOR,<val>` | - | P1 | yes | Element colours by 1 group, 2 material, 3 property (ElemPalette, 128 colours). |
| **ELENUM** | `ELENUM,[opt]` | - | P1 | yes | Element number labels (-1 toggle, 0 off, 1 on); turns GROUPNUM off. |
| **GROUPNUM** | `GROUPNUM,[opt]` | - | P1 | yes | Group number labels (-1 toggle, 0 off, 1 on); turns ELENUM off. |
| **LAYERPLOT** | `LAYERPLOT` | - | P1 | yes | Soil layer column and property table of the active model (TOPL / L, SITE half-space). |
| **LBINCORS** | `LBINCORS,<out>,<in>` | - | P2 | yes | Lower-bound incoherent RS -- no algorithm in the manual (D-LIN-03). |
| **LINECOMBIN** | `LINECOMBIN,<Dest>,<Line1>,<Coeff1>,...,[Line100],[Coeff100]` | - | P1 | yes | Sum of c_i y_i ('Linear Combin.'). |
| **LINENAME** | `LINENAME,<Num>,<Name>` | - | P1 | yes | Rename a line object (global: every graph shows the new name). |
| **MARKERS** | `MARKERS,<Mark>,<Ln1>,...,<Ln50>` | - | P1 | yes | Data-point markers off (0) / on (1) for lines (global). |
| **MODELPLOT** | `MODELPLOT` | - | P1 | yes | Element plot of the active model (colours by ELECOLOR). |
| **NODENUM** | `NODENUM,[opt]` | - | P1 | yes | Node number labels (-1 toggle, 0 off, 1 on). |
| **NODEPLOT** | `NODEPLOT` | - | P1 | yes | Node plot of the active model (element-connected nodes; interaction nodes red). |
| **NODESEL** | `NODESEL,<N1>,...,<N20>` | - | P1 | yes | Toggle the selection of nodes (model level: all current and future plots). |
| **PAUSE** | `PAUSE,[pz]` | - | P1 | yes | Animation -1 toggle start/stop (default), 0 start, 1 stop. |
| **PLOTRANGE** | `PLOTRANGE,<Xmin>,<Xmax>,<Ymin>,<Ymax>` | - | P1 | yes | Extent of the active 2D plot (blank = data extent). |
| **PLOTTITLE** | `PLOTTITLE,<Title>` | - | P1 | yes | Title of the active plot (2D and 3D; ignored by the soil layer plot). |
| **PROCFRAME** | `PROCFRAME,[AniFile],[BufferDIR],[Data],[anitype]` | - | P1 | yes | Frame list -> frame store + SASSIani.xml. |
| **READSPEC** | `READSPEC,<SpecFile>,<numLines>,<Line1>,...,<LineN>` | - | P1 | yes | Load spectrum columns (frequency column not counted). |
| **READTH** | `READTH,<THFile>,<Pair>,<Num>` | - | P1 | yes | Load a time history (Pair 0 one column with dt first, 1 time/value pairs). |
| **RSTCENTER** | `RSTCENTER` | - | P1 | yes | Rotation centre = centre of the bounding box of the element-connected nodes. |
| **RSTVIEW** | `RSTVIEW` | - | P1 | yes | Reset the view (default isometric view, zoom and pan fitting the model). |
| **SHADEROPTIONS** | `SHADEROPTIONS,[points],[linew],[shrink],[scale]` | - | P1 | yes | Node size, outline, shrink fraction, scale factor. |
| **SHOWDOF** | `SHOWDOF,[label1],...,[label6]` | - | P1 | yes | Mark nodes with fixed DOFs (X Y Z XX YY ZZ DISP ROT ALL; NONE off). |
| **SHOWMASS** | `SHOWMASS,[opt]` | - | P1 | yes | Lumped-mass markers (-1 toggle, 0 off, 1 on); red X, green Y, blue Z. |
| **SHRINK** | `SHRINK,[switch]` | - | P1 | yes | Shrink the elements of the element plot (1 on, 0 off, -1 toggle). |
| **SOILPROPPLOT** | `SOILPROPPLOT,<PropName>` | - | P1 | yes | G/Gmax and damping vs shear strain % of a DYNP property (case-sensitive). |
| **SPECPLOT** | `SPECPLOT,<Line1>,...,<Line50>` | - | P1 | yes | Spectrum plot of line objects (-1 ends the list). |
| **SRSS** | `SRSS,<dest>,<source1>,...,<source100>` | - | P1 | yes | Square root of the sum of squares ('SRSS Line'). |
| **STIPPLE** | `STIPPLE,<switch>` | - | P1 | yes | Dash patterns on the lines of the active plot (-1 toggle, 0 off, 1 on). |
| **SUBTRACTION** | `SUBTRACTION,<dest>,<source1>,...,<source100>` | - | P1 | yes | First line minus the others ('Linear Combin.'). |
| **THPLOT** | `THPLOT,<Line1>,...,<Line50>` | - | P1 | yes | Time-history plot of line objects (-1 ends the list). |
| **VECTORPLOT** | `VECTORPLOT,<BufferDir>,<MnF>,<MxF>,<ST>,<Scale>` | - | P1 | yes | Animated complex-TF vectors (X red, Y green, Z blue). |
| **WINDOWSETTINGS** | `WINDOWSETTINGS,[<field>,<values>]` | - | P1 | yes | Settings window of the active plot. |
| **WIREFRAME** | `WIREFRAME,<Switch>` | - | P1 | yes | Wireframe mode of the element plot (0 off, 1 on, -1 toggle). |
| **WRITESPEC** | `WRITESPEC,<SpecFile>,<Num1>,...,[Num50]` | - | P1 | yes | Write lines on the union of their abscissas. |
| **WRITETH** | `WRITETH,<THFile>,<Num>` | - | P1 | yes | Write a line as a one-column history (dt = x2 - x1). |
| **XTITLE** | `XTITLE,<label>` | - | P1 | yes | X-axis title of the active 2D plot. |
| **YTITLE** | `YTITLE,<label>` | - | P1 | yes | Left Y-axis title of the active 2D plot. |
| **YTITLE2** | `YTITLE2,<label>` | - | P1 | yes | Right Y-axis title (soil property plot: shear modulus axis). |

## 3.4.M Programming (variables, loops, macros)

| Command | Syntax | Abbreviation / alias | Tier | Implemented | Meaning |
|---|---|---|---|---|---|
| **ADDRND** | `ADDRND,<var>,<numsamples>,<dist>,<param1>,...` | - | P1 | yes | Append random values until the list has numsamples items. |
| **FOREACH** | `FOREACH,<Var>,<Command>` | - | P1 | yes | Run Command once per item of Var ('#' = 1-based loop index). |
| **LOADMACRO** | `LOADMACRO,<Name>,<MACROFILE>` | - | P1 | yes | Load and cache a macro file under an upper-case name. |
| **LOADVAR** | `LOADVAR,<filename>,[name]` | - | P1 | yes | Variable from a text file, one item per line (name = file stem). |
| **MACRO** | `MACRO,<Name>,<1>,...,<N>` | - | P1 | yes | Run a loaded macro, replacing `$k$` with argument k. |
| **MACROLIST** | `MACROLIST` | - | P1 | yes | List the loaded macros and their files. |
| **REDUCESET** | `REDUCESET,<var>,[sorttype]` | - | P1 | yes | Sort ascending and remove duplicates (STRING default, INT, FLOAT). |
| **RND** | `RND,<var>,<numsamples>,<dist>,<param1>,...` | - | P1 | yes | Fill a variable with random values (overwrites it). |
| **RNDSEED** | `RNDSEED,<seed>` | - | P1 | yes | Seed the random generator of RND/ADDRND (positive integer). |
| **SETVAR** | `SETVAR,...` | - | P1 | yes | No effect on the model; its arguments are substituted, so counters change. |
| **SHOWVAR** | `SHOWVAR,<varname>` | - | P1 | yes | Show the items (1-based) and the counter of a variable. |
| **VAR** | `VAR,<Name>,<X1>,...,<Xn>` | - | P1 | yes | Define or redefine a variable (the counter is reset to 0). |
| **VARLIST** | `VARLIST` | - | P1 | yes | List the variables (value of single-item variables, item count otherwise). |

## 3.4.N Water modelling

| Command | Syntax | Abbreviation / alias | Tier | Implemented | Meaning |
|---|---|---|---|---|---|
| **FILLPOOL** | `FILLPOOL,<Stiff>,<Sensitivity>,<EmptyLevels>,<ShellArea>,<offset>,<stiff2>` | - | P2 | yes | Fill the pool sub-model with water SOLIDs, wall springs and optional area shells. |
| **LISTPOOLINTER** | `LISTPOOLINTER,<Pool>` | - | P2 | yes | List the interface nodes of FILLPOOL pool model &lt;Pool> that the active (original) model has with the same number and position. |
| **MERGEPOOL** | `MERGEPOOL,<Pool>` | - | P2 | yes | (SASSI-EDU) import the water, interface springs and area shells of FILLPOOL pool model &lt;Pool> into the active model. |
| **POOLDATA** | `POOLDATA,<water>,<springs>,<shells>,<stiff>,<sens>,<empty>,<shellarea>,<offset>,<stiff2>,<watermat>,<zfloor>,<zsurface>` | - | P2 | yes | (SASSI-EDU) FILLPOOL pool data. |
| **REFINEMODEL** | `REFINEMODEL` | - | P2 | yes | Split every quadrilateral SHELL, TSHELL and PLANE element into 4 and every hexahedral SOLID into 8. |

## 3.4.O Option NON and nonlinear soil

| Command | Syntax | Abbreviation / alias | Tier | Implemented | Meaning |
|---|---|---|---|---|---|
| **B** | `B,<num>,<group>,<spgroup>,<bbc>,<force>,<end1>,<end2>` | - | P2 | yes | Nonlinear beam -- not usable in this version (stored only, as the manual). |
| **BBC** | `BBC,<num>,<type>,<points>,<yield>,<file>` | - | P2 | yes | Backbone curve &lt;num> from a file of X Y pairs (one pair per line, the origin not given; exactly &lt;points> pairs are read); &lt;type> 1 CMS, 2 CMB, 3 TAK, 4 GMR (titles only); &lt;yield> = number of the yield point. |
| **BBCGEN** | `BBCGEN,<Panel>,<ShearModel>,[fc],[fy],[Pn],[Nu],[bre],[bys],[CrackingForceLevel]` | - | P2 | yes | 22-point backbone curve of panel &lt;Panel> (0 = all) from its ultimate shear by ShearModel 1 ACI 318-08, 2 Wood 1990, 3 Barda 1977, 4 Gulec-Whittaker 2009 (D-NON-09): cracking point V_cr = 3 sqrt(f'c) A_W ... |
| **BBCI** | `BBCI,<num>,<yield>,<type>` | - | P2 | yes | Set the yield point number and the type (1 CMS, 2 CMB, 3 TAK, 4 GMR) of backbone curve &lt;num>. |
| **BBCP** | `BBCP,<num>,<point>,<X>,<Y>` | - | P2 | yes | Set point &lt;point> of backbone curve &lt;num> to (X, Y), as one row of the BBC grid of the NONLINEAR dialog; &lt;point> = number of points + 1 appends a point (a curve can be built point by point after BBCI). |
| **BBCX** | `BBCX,<num>,<points>,<yield>,<X1>,...,<Xn>` | - | P2 | yes | X values (strain or displacement, the origin not given) of backbone curve &lt;num>. |
| **BBCY** | `BBCY,<num>,<points>,<yield>,<Y1>,...,<Yn>` | - | P2 | yes | Y values (force or moment) of backbone curve &lt;num>. |
| **BEAMPILE** | `BEAMPILE` | - | P2 | no | Not applicable / not usable (message) -- BEAMPILE is not applicable in this version |
| **DCOUPLEBEAM** | `DCOUPLEBEAM` | - | P2 | no | Not applicable / not usable (message) -- DCOUPLEBEAM is not usable in this version |
| **DELBBC** | `DELBBC,<start>,[<end>],[<stride>]` | - | P2 | yes | Delete backbone curves (all their BBCI/BBCX/BBCY data). |
| **DELBM** | `DELBM,<start>,[<end>],[<stride>]` | - | P2 | yes | Delete nonlinear beam records (the beams stay). |
| **DELNLS** | `DELNLS,<start>,[end],[stride]` | - | P2 | yes | Delete NLSLAYER sets start, start+stride, ... <= end (end = start by default; a blank stride or a stride < 1 uses 1 with a warning); linear soil data are not changed. |
| **DELSPR** | `DELSPR,<start>,[<end>],[<stride>]` | - | P2 | yes | Delete nonlinear spring records (the springs stay as linear elements). |
| **DGRDFLR** | `DGRDFLR,<scale>` | - | P2 | yes | Multiply Young's modulus of the materials of the floor panels (horizontal shell groups) by &lt;scale> (after WALLFLR and GROUPMAT); materials also used by other groups are skipped. |
| **EDGE** | `EDGE,<panel>,[X],[Y],[Z]` | - | P2 | yes | Split wall group &lt;panel> along the lines of its outer and opening edges (flags: 0 use the edges parallel to that global axis, 1 ignore them; oblique edges are always used); the cells (top first, then left to right) ... |
| **EDGEMODEL** | `EDGEMODEL,[x],[y],[z]` | - | P2 | yes | EDGE on every wall (vertical shell) group of the model with the same flags. |
| **EQL** | `EQL,<disp>,<NonLinOpts>,<dampCutoff>,<dampScale>,<ElasicD>` | - | P2 | yes | NONLINEAR header -- EDF (0.7-0.9, default 0.8), element types (1 panels + 2 springs; 4 beams not available; blank = the types with P/S records), damping cut-off % (0 none), damping scale (0 = 1), include elastic damping ... |
| **MERGEPANEL** | `MERGEPANEL,<Panel>` | - | P2 | yes | Append the groups and materials of panel model &lt;Panel> to the active (original) model whose shell groups were deleted; nodes are matched by number (coordinates checked), missing nodes are added; element materials are ... |
| **NLSLAYER** | `NLSLAYER,<Num>,[curvefit],[B],[S],[refStrain],[Vis]` | - | P2 | yes | SOIL-NON data of soil sublayer &lt;Num> -- curvefit 1 fits the MKZ parameters to the sublayer's DYNP G/Gmax curve (B, S, refStrain, Vis ignored), 0 uses Beta, S exponent, Reference Strain (in %) and Viscosity of tau = G0 g ... |
| **NLSOIL** | `NLSOIL,<Opt>,<NSTimeSunInc>,<DispConv>,<ForceConv>,<EqualIt>,<BedInt>,<NLDampType>,<MMmult>,<SMmult>` | - | P2 | yes | Global SOIL-NON options (Opt 1 = nonlinear time-domain SOIL; sub-increments, convergence errors, equilibrium iterations, bedrock 0 rigid / 1 viscoelastic, damping type 1 frequency independent / 2 visco-elastic / 3 ... |
| **NONLINBAT** | `NONLINBAT,<Sel>` | - | P2 | yes | Write &lt;model>\_NONLINBAT.pre, a generic batch script of the whole Option NON analysis (elastic SSI, then iterations with NONLINITER): Sel 0 one input direction, 1 X/Y/Z input with COMB_XYZ_THD (a COMB_XYZ_THD.inp listing ... |
| **NONLINMOTDISP** | `NONLINMOTDISP` | NONLINMODISP | P2 | yes | Add the four corner nodes of every panel (node-connection counting) and the end nodes of every S spring to the MOTION (NOUT, transfer functions) and RELDISP (RDND) output requests, without duplicates. |
| **P** | `P,<num>,<group>,<bbc>,<disp>,<force>` | - | P2 | yes | Wall panel &lt;num> = SHELL group &lt;group> (vertical, coplanar) with backbone curve &lt;bbc>, Disp. |
| **PANELIZE** | `PANELIZE` | - | P2 | yes | Split every shell group of a panel model (after WALLFLR) along its intersections with the other shell groups (edges whose two nodes also belong to another group); each connected part becomes a group (the first keeps the ... |
| **PDEL** | `PDEL,<start>,[<end>],[<stride>]` | - | P2 | yes | Delete wall panel records (the elements stay). |
| **PLIST** | `PLIST,[start],[end],[stride]` | - | P2 | yes | List / check the wall panels: group, elements, corners, width and height, coplanarity, material, BBC and the supported options. |
| **PNLGEN** | `PNLGEN (PANELGEN)` | PANELGEN | P2 | yes | One panel record per vertical shell group, numbered 1, 2 ... in group order, with BBC = panel number, Disp. |
| **S** | `S,<num>,<group>,<elem>,<bbc>,<disp>,<force>` | - | P2 | yes | Nonlinear spring &lt;num> = SPRING element &lt;elem> of group &lt;group> with backbone curve &lt;bbc> on DOF &lt;disp> (1 X, 2 Y, 3 Z translation) and Force Opt 4 (GMR, Error 127 otherwise). |
| **SOILREDEF** | `SOILREDEF,<soil>,<dir>` | - | P2 | no | Excavation layers from a 2D soil model (not validated: warning) -- prints 'SOILREDEF is not available in this build (tier P2)' |
| **SOLIDPILE** | `SOLIDPILE,<group>,[stiff],[soft],[stiff2]` | - | P2 | yes | Separate the SOLID pile group &lt;group> from the soil -- its interface nodes are duplicated, the pile elements reconnected to the duplicates and the two joined by springs (4 % damping) in 4 new groups: side-wall X and Y ... |
| **UNIPNL** | `UNIPNL,<group>` | - | P2 | yes | One group per element of &lt;group> (element 1 stays, the others go to new appended groups, each numbered element 1) -- one panel per shell, e.g. for curved walls. |
| **WALLFLR** | `WALLFLR` | - | P2 | yes | Separate shell walls and floors -- deletes every non-shell element of the active model, puts each set of >= 5 coplanar shells in its own group (groups 2, 3 ...; titles by orientation) and the other shells in group 1. |

## 3.4.P Binary databases

| Command | Syntax | Abbreviation / alias | Tier | Implemented | Meaning |
|---|---|---|---|---|---|
| **ACCDBANI** | `ACCDBANI,<dir>,[label]` | ACCANIDB | P2 | no | Animations from DBs -- prints 'ACCDBANI is not available in this build (tier P2)' |
| **BINFRAMEOUT** | `BINFRAMEOUT,<db>,<frame>,<TS>,[Split],<dir>` | - | P2 | no | ASCII frames from a DB (−1 = max frame) -- prints 'BINFRAMEOUT is not available in this build (tier P2)' |
| **BINOUT** | `BINOUT,[mot],[str],[reldisp]` | - | P0 | stored only | Flags (blank = unchanged); reldisp 2 = THD -- parsed and stored for WRITE; no computation in this build (tier P0) |
| **BINSTRTBL** | `BINSTRTBL,<group>,<EVar>,[step],<file>` | - | P2 | no | CSV stress table (−1 = signed abs-max) -- prints 'BINSTRTBL is not available in this build (tier P2)' |
| **COMBACCDB** | `COMBACCDB` | - | P2 | no | Combine directional DBs (algebraic sum; COMBDISPDIR assembles components) -- prints 'COMBACCDB is not available in this build (tier P2)' |
| **COMBDISPDB** | `COMBDISPDB` | COMDISPDB | P2 | no | Combine directional DBs (algebraic sum; COMBDISPDIR assembles components) -- prints 'COMBDISPDB is not available in this build (tier P2)' |
| **COMBDISPDIR** | `COMBDISPDIR` | - | P2 | no | Combine directional DBs (algebraic sum; COMBDISPDIR assembles components) -- prints 'COMBDISPDIR is not available in this build (tier P2)' |
| **COMBTHSDB** | `COMBTHSDB` | - | P2 | no | Combine directional DBs (algebraic sum; COMBDISPDIR assembles components) -- prints 'COMBTHSDB is not available in this build (tier P2)' |
| **DELDB** | `DELDB,[sel]` | - | P2 | no | ALL, ACC, THS, DISP -- prints 'DELDB is not available in this build (tier P2)' |
| **DISPDBANI** | `DISPDBANI,<dir>,[label]` | - | P2 | no | Animations from DBs -- prints 'DISPDBANI is not available in this build (tier P2)' |
| **LOADACCDB** | `LOADACCDB,<file>,[sel]` | LOADACCDBANI | P2 | no | Load one DB of each type (sel 0 SASSI, 1 ANSYS) -- prints 'LOADACCDB is not available in this build (tier P2)' |
| **LOADDISPDB** | `LOADDISPDB,<file>,[sel]` | - | P2 | no | Load one DB of each type (sel 0 SASSI, 1 ANSYS) -- prints 'LOADDISPDB is not available in this build (tier P2)' |
| **LOADTHSDB** | `LOADTHSDB,<file>,[sel]` | - | P2 | no | Load one DB of each type (sel 0 SASSI, 1 ANSYS) -- prints 'LOADTHSDB is not available in this build (tier P2)' |
| **MAXDBFRAME** | `MAXDBFRAME,<Type>,[dir]` | - | P2 | no | Single-frame maximum animation -- prints 'MAXDBFRAME is not available in this build (tier P2)' |
| **THSDBANI** | `THSDBANI,<dir>,[label]` | - | P2 | no | Animations from DBs -- prints 'THSDBANI is not available in this build (tier P2)' |

## 3.4.Q Thick shell

| Command | Syntax | Abbreviation / alias | Tier | Implemented | Meaning |
|---|---|---|---|---|---|
| **THSHLSMH** | `THSHLSMH,<passes>,[type],[workdir]` | - | P2 | yes | Smooth the TSHELL transverse shear forces of STRESS over rings of neighbouring elements (type 0 Parzen weights, 1 plain average; NON VALIDATED in this version). |
| **THSHLSTR** | `THSHLSTR,<flag>` | - | P2 | yes | TSHELL face stresses/strains in STRESS (0 the 8 basic components only, 1 also the top/bottom face values from their maxima, D-TSH-01). |

## 3.4.R Extension commands of SASSI-EDU

| Command | Syntax | Abbreviation / alias | Tier | Implemented | Meaning |
|---|---|---|---|---|---|
| **ACTIVATEPLOT** | `ACTIVATEPLOT,<id>` | - | P1 | yes | Make open plot &lt;id> the active plot (SASSI-EDU extension, rule L17). |
| **ANALYSX** | `ANALYSX,<ffm>,<delrst>` | - | P0 | yes | Free-Field Load (0) / Free-Field Motion (1) application of the incoherency factors; delete restart files after a successful restart run. |
| **BUILDFILE77** | `BUILDFILE77,<out>,<in1>,...,<inN>` | - | P1 | yes | Build_FILE77 -- combine the FILE77s of per-level HOUSE runs into one FILE77 (spec 05b section 7; nodes in the FILE4 interaction order when &lt;model>.N4 exists). |
| **CMODFORM** | `CMODFORM,<form>` | - | P0 | yes | Complex modulus 0 = 1-2b^2+2ib sqrt(1-b^2) (default), 1 = 1+2ib (benchmarks). |
| **COMBXYZSTRAIN** | `COMBXYZSTRAIN,<FILE74x>,<FILE74y>,<FILE74z>,<out>` | - | P1 | yes | COMB_XYZ_STRAIN (D-NLS-04) -- SRSS of the directional effective strains per element -> &lt;out> (default FILE74), with G and beta from FILE73; blank inputs skipped. |
| **COMBXYZTHD** | `COMBXYZTHD,<inpfile>` | - | P2 | yes | COMB_XYZ_THD (D-NON-08) -- for every line 'outfile fileX fileY fileZ' of &lt;inpfile> (default COMB_XYZ_THD.inp; line 1 = number of lines) the output is the time-step-wise sum of the X, Y and Z relative-displacement ... |
| **EDUOPT** | `EDUOPT,<key>,<value>` | - | P0 | yes | Algorithm switch of requirements section 7; a blank value resets the key. |
| **FCOPY** | `FCOPY,<src>,<dst>` | - | P0 | yes | Copy a file of the model directory (e.g. FILE1 -> FILE1X, FILE8 -> FILE81). |
| **FMOVE** | `FMOVE,<src>,<dst>` | - | P0 | yes | Rename a file of the model directory (e.g. FILE8 -> FILE82). |
| **HARMFRAME** | `HARMFRAME,<Src>,<Freq>,<OutDir>,[NFrames],[Ref]` | - | P1 | yes | Steady-state harmonic frames u = Re(H e^{iwt}) of every node over one period at the computed frequency closest to &lt;Freq> (SASSI-EDU extension; PROCFRAME and DEFORMPLOT animate them). |
| **HOUSEX** | `HOUSEX,<optimize>,<supmode>,<nsim>,<nlssi>,<ansys>` | - | P0 | yes | HOUSE dialog fields without an argument (node optimizer, Linear/Quadratic superposition, stochastic simulations, non-linear soil SSI, ANSYS model input). |
| **LGFILE** | `LGFILE,<key>,<name>` | - | P2 | yes | A file or path box of the LOADGEN dialogs (keys SSIPATH HOUSE DISP DISPROT ACC ACCROT ANSYSPATH LUMPED MASTER APDL APDLDYN GROUND); &lt;name> is the rest of the line, blank = default. |
| **LGLIST** | `LGLIST,[STATIC\|DYNAMIC]` | - | P2 | yes | List the LOADGEN (Option A) settings and the deck RUNLOADGEN would write. |
| **LGMAP** | `LGMAP,<mode>,<tol>,<file>` | - | P2 | yes | ANSYS node numbering -- IDENTITY (default), PAIRS (file of 'sassi_node ansys_node' lines) or COORD (ANSYS node at the same position in a .cdb / APDL file, tolerance &lt;tol>). |
| **LGNODE** | `LGNODE,<kind>,<nodes>` | - | P2 | yes | Add nodes to a LOADGEN node list (kind D interface nodes that receive D, M master nodes, A acceleration check nodes); LGNODE,&lt;kind>,0 clears the list. |
| **LGOPT** | `LGOPT,<digits>,<opmode>,<rest>` | - | P2 | yes | Significant digits of the APDL values (6..17, default 12), data-check mode (1 = check the inputs, no file written) and &lt;rest> (1, default: the relative displacements start at rest -- their initial value, the zero-mean ... |
| **LGTIME** | `LGTIME,<crit>,...` | - | P2 | yes | Critical times of the static loads -- V\|VX\|VY\|VZ,&lt;n>,&lt;tsep>; MX\|MY\|MZ,&lt;n>,&lt;tsep>,&lt;x0>, &lt;y0>,&lt;z0>; ACC\|DISP,&lt;n>,&lt;tsep>,&lt;node>,&lt;dof>; TIME,&lt;t1>,&lt;t2>,...; STEP,&lt;k1>,&lt;k2>,... (default V,1). |
| **LIBRARY** | `LIBRARY,[kind\|@name\|DEFAULTS]` | - | P0 | yes | List the built-in input library (files named @&lt;file>, sassi/data/library), one kind (RECORD, LOAD, SPECTRUM, PSD, DYNP) or one file, or the defaults of blank inputs the active model uses (DEFAULTS). |
| **LOADGEN** | `LOADGEN,<data>,<multi>,<rotdisp>,<rotacc>,<masstype>,<genmass>,<source>` | - | P2 | yes | ANSYS Eq. Static Load options (Option A): data 1 DISP / 2 ACC / 3 DISPACC / 4 SOILDISP, Use Multiple File List Inputs, Rotational Disp., Rotational Accel., Mass Type 1 LUMPED / 2 MASTER, Generate Mass Data, history ... |
| **LOADGENDYN** | `LOADGENDYN,<alpha>,<beta>,<method>,<refnode>,<source>,<rotdisp>,<rotacc>,<zeta>,<f1>,<f2>,<gfopt>,<gmult>` | - | P2 | yes | ANSYS Dynamic Load options (Option A): Rayleigh alpha/beta (or from zeta at f1, f2), method REL (ground ACEL + relative D, manual) / ACC (fixed base driven by &lt;refnode>), reference node (0 = control motion), source ... |
| **MOTIONX** | `MOTIONX,<f1213>,<resp>,<srss>,<savetf>,<saveacc>,<savers>,<saverot>,<rsttf>,<rstacc>,<rstrs>` | - | P0 | yes | Save FILE12/13 (0/1 FILE13/2 FILE12); vibration response 0 disp / 1 vel / 2 acc; use SRSSTF.txt; post-processing Save/Restart flags |
| **NLSSIITER** | `NLSSIITER,<var>[+<var2>...],[<maxit>],[<tolG>],[<tolB>]` | - | P1 | yes | Run the iteration commands stored in the variable(s) (one per item) until the nonlinear soil properties converge (max \|dG/G\| < tolG %, max \|d beta\| < tolB %, D-NLS-02) or &lt;maxit> passes (default 8); '#' in the body = ... |
| **NLSSIRESET** | `NLSSIRESET` | - | P1 | yes | Start a new nonlinear soil SSI analysis -- &lt;model>.liq = 0 (the next HOUSE run takes the .pin initial properties), a new NLSOIL_CONVERGENCE.TXT history, and FILE78 / FILE74 of the previous analysis renamed FILE78.prev / ... |
| **NONLINITER** | `NONLINITER,<var>[+<var2>...],[<maxit>]` | - | P2 | yes | Run the command lines stored in the variable(s) (one command per item; quote items with commas) once per pass, '#' = pass number, until NONLINEAR reports convergence (D-NON-06: \|dE/E\| < 2 %, \|d xi\| < 0.5 %) or &lt;maxit> ... |
| **NONLINRESET** | `NONLINRESET` | - | P2 | yes | Delete PANEL.NON, SPRING.NON, the \*\_EQL_Matl_Prop.txt files and the convergence history in the model directory -- the next NONLINEAR run is an elastic run (re-run AFWRITE for the elastic HOUSE deck). |
| **NONLINSAVE** | `NONLINSAVE,<tag>,[<file8>]` | - | P2 | yes | Copies of the NONLINEAR results under the Fig. 1.2 names -- PanelNNNN\_&lt;tag> .thd/.ths/.crv, Panel\_&lt;tag>.fmu, Panel_EQL_Matl_Prop\_&lt;tag>.txt (SPRING likewise), &lt;model>\_&lt;tag>.hou; &lt;file8> = 1 also copies FILE8 (FILE8X/Y/Z) ... |
| **NONLINTHD** | `NONLINTHD,<prefix>` | - | P2 | yes | Copy the .THD files NONLINEAR reads (panel corners X/Y/Z, spring ends) to &lt;prefix>\_&lt;name> (e.g. X_00011TR_X.THD after the X-input RELDISP run) for COMB_XYZ_THD. |
| **PIN** | `PIN,<esf>,[<ncurv>]` | - | P1 | yes | .pin line 1 -- effective strain factor ESF (blank: SOIL &lt;ratio>) and number of FILE73 curves NCURV (blank: number of DYNP labels); NGRP is counted by AFWRITE. The groups follow with PINGRP / PINMAT (GFAC, DFAC relative ... |
| **PINDEL** | `PINDEL,[<igrp>]` | - | P1 | yes | Delete nonlinear group &lt;igrp> (PINGRP and its PINMAT records); blank = all .pin records. |
| **PINGRP** | `PINGRP,<igrp>,<istr>,[<gfac>],[<dfac>],[<icurve>]` | - | P1 | yes | Nonlinear soil group (SOLID/PLANE structure elements, ETYPE 1, whose M-table materials hold the low-strain soil: G_max) with strain flag ISTR (0 maximum shear component, 1 octahedral / 2D maximum shear) and optional ... |
| **PINLIST** | `PINLIST` | - | P1 | yes | List the .pin records (PIN, PINGRP, PINMAT) and the .pin text AFWRITE writes. |
| **PINMAT** | `PINMAT,<igrp>,<mat>,<gfac>,<dfac>,<icurve>` | - | P1 | yes | .pin line of material &lt;mat> of nonlinear group &lt;igrp>: iteration 0 has G = GFAC G_layer and beta = DFAC beta_layer of the free-field TOPL layer at the element (1.0 = the free field), later iterations G = G_max (G/G_max) ... |
| **RELDX** | `RELDX,<saverot>,<rstframes>` | - | P1 | yes | Save Rotations for ANSYS; Restart for Frame Generation (P1). |
| **REMOVEFREQ** | `REMOVEFREQ,<infile>,<outfile>,<n1>,...` | - | P1 | yes | Remove_Frequencies_from_FILE8 (P1). |
| **RUNLOADGEN** | `RUNLOADGEN,[STATIC\|DYNAMIC],[model]` | - | P2 | yes | Write &lt;model>.lgn from the LOADGEN records and run LOADGEN in the model directory (the Ok of the ANSYS Eq. Static Load / ANSYS Dynamic Load dialogs; STATIC by default). |
| **SHOWSOIL** | `SHOWSOIL,[opt],[cut],[margin],[depth]` | - | P1 | yes | The free-field soil drawn around the foundation (SASSI-EDU, display only). |
| **SITEX** | `SITEX,<soilmode>` | - | P0 | yes | 0 Linear Soil (L table), 1 Non-Linear Soil (FILE88 from SOIL), D-SIT-06. |
| **SOILX** | `SOILX,<indir>,<mult>,<max>,[<cl>],[<file>]` | - | P0 | yes | SOIL input direction, SOIL scaling (exactly one of mult/max non-zero) and optional SOIL-only control layer and history file (D-SOL-11, D-SOL-13). |
| **STRESSX** | `STRESSX,<pzadj>,<smo>,<skip>,<savemax>,<saveth>,<rstns>,<rstsp>` | - | P0 | yes | STRESS dialog fields. |
| **VERIFY** | `VERIFY,[id\|P0\|P1\|P2\|ALL\|LIST]` | - | P0 | yes | Run verification problems and print computed, reference, error, tolerance and pass/fail (requirements section 6.1); blank = P0. |
| **VERIFYREPORT** | `VERIFYREPORT,[sel],[file],[fig]` | - | P0 | yes | Run verification problems (sel = P0 default, P1, P2, ALL, FAST or ids joined by '+') and write the Markdown verification report (default VERIFICATION_REPORT.md). |
