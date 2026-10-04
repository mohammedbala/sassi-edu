# `sassi.prep` — command interpreter and model editing

This package is the "UI" of SASSI-EDU: the `.pre` command language, the numbered models in
memory, WRITE/SAVE/RESUME, and the command handlers. The console (`sassi`), `.pre` files (INP),
macros, FOREACH and the GUI all go through one `Interpreter` (requirements UI-01, L17).

```
sassi/model/           model database (SSIModel, entities, option records, geometry, materials)
sassi/prep/
  lexer.py             L1-L9: head/argument split, blank fields, quotes, rest-of-line text,
                       legacy '(...)' token, id lists 'a-b'
  messages.py          MessageSink with the six message classes and the display filters
  registry.py          @command decorator, the catalogue of every command (abbreviation, tier,
                       class, lexical options), lookup (L10), placeholders of missing commands
  interpreter.py       Interpreter (models, CWD, variables, macros, INP stack, substitution L13),
                       Call (what a handler receives)
  writer.py            WRITE: model -> .pre (UT-03 round trip), hooks for typed writers
  commands/            handler modules, imported automatically (every *.py of the package)
    session.py         INP WRITE SAVE RESUME MDL MDLNAME TIT STATUS ACTM CPMODEL DMODEL MODELLIST
                       CD MKDIR MOPT GRAVITY GROUNDELEV FREQ LFREQ DAMP TOPL AMP NOUT EOUT RDND
    nodes.py           N NDEL NGEN FILL NMED LMOVE NMOVE NSCALE NLIST D INT INTLIST CSYS LOC LOCAL
                       GLOBAL SDEL SLIST
    elements.py        GROUP MTYPE GTIT GDEL GLIST E EGEN EDEL ECOMPR ELIST ETYPE EINT THICK KI KJ
                       MSET MACT RSET RACT
    tables.py          M L R SC MXR MXI MXM DELM DELL DELR DELSC MXDEL MLIST LLIST RLIST SCLIST MXLIST
    loads.py           F MM FDEL MMDEL FLIST MMLIST FSCALE MSCALE MT MR MTGEN MRGEN MTDEL MRDEL
                       MTSCALE MRSCALE MTLIST MUNITS
    program.py         VAR SETVAR SHOWVAR VARLIST FOREACH LOADMACRO MACRO MACROLIST LOADVAR RND
                       ADDRND RNDSEED REDUCESET
sassi/cli.py           `sassi` console, `sassi run file.pre [--quiet]`, `sassi -c LINE`, `--version`
```

## Processing of one line

1. `$k$` macro placeholders (MACRO does this before the line reaches `execute`).
2. `@NAME`, `@NAME[i]`, `@NAME++`, `@NAME+k`, `@NAME=k` and, inside FOREACH, `#` / `#+k` —
   left to right, each once (D-PAR-15/16). Skipped for commands registered with `raw=True`
   (FOREACH: its body is substituted afresh at every iteration).
3. Lexing with the command's options (`text_from`, `paren`). With `text_from=k` only the fields
   before argument k are scanned; argument k is the rest of the line **verbatim** (an unbalanced
   `"` included, L7 / D-PAR-06).
4. `max_args`: arguments after `max_args` are **dropped** (with the warning "arguments after
   argument n ignored") before the handler runs, so `c.nargs` never exceeds `max_args`. Declare
   `max_args` only when the documented argument count is a hard limit; a handler that accepts
   more than the manual's count (FREQ, DAMP) omits it and warns itself.
5. Dispatch to the handler. `CommandError` → error message, the line is skipped, processing
   continues (L12). `CommandReported` (a subclass) → the command fails *without* a further
   message, for handlers that have already reported the problem (INP of an unreadable file).
   Any other exception is reported as an internal error, never fatal; listing commands must
   therefore catch the `ValueError` of derived values (`frequency_step()`, `nodal_masses()`)
   and report them as unavailable.

`execute()` returns False for a failed command (`last_ok`, `sassi -c` exit status, and the line
is not added to the replay history). `.pre`, macro and variable files are read with
`interpreter.read_command_lines(path)` (UTF-8, a byte-order mark removed, bad bytes replaced).

While INP, a macro or a FOREACH runs, ECHO and CONFIRM messages are suppressed (L16). Each INP
file ends with an INFO summary `INP <file>: n lines, n commands, n warnings, n errors`, and the
outermost one with `INPUT FILE REACHED EOF, INPUT SWITCHED TO KEYBOARD`.

## Adding a command (for the option / CHECK / AFWRITE / RUN work package and others)

Create a module in `sassi/prep/commands/` (it is imported automatically; an import error is
reported as a warning at interpreter start-up and does not stop the other modules):

```python
from sassi.prep.registry import command, CommandError

@command("SITE")                      # abbreviation, tier, class, text options from the catalogue
def cmd_site(c):
    """SITE,<opmode>,<mode1>,...: SITE options (record setter)."""
    rec = SiteRecord("SITE", [c.raw(k) for k in range(1, 15)])   # see "Typed option records"
    if rec.nl not in (0,) and not 4 <= rec.nl <= 20:
        c.warn("<nl> must be 0 or 4..20 (CHECK Error 47)")
    c.model.options.set_record(rec)
    c.confirm("SITE options set")
```

* The **catalogue** in `registry.py` already lists every command of requirements section 3.4 with
  its abbreviation, tier, state class and lexical options. A handler normally needs only the
  name. Decorator options override the catalogue: `abbrev=("ABCD",)`, `tier`, `cls`,
  `text_from=k` (argument k takes the rest of the line, rule L7), `raw=True` (unsubstituted
  rest of line), `paren=True` (legacy final `(...)` token, L9), `max_args=n` (later arguments
  are dropped with a warning), `summary`.
* A real handler **replaces the placeholder** of a catalogued command. Registering a second real
  handler for the same name raises `ValueError` unless `replaces=True` is given — e.g. a typed
  `DAMP` handler that should replace the storage-only one in `session.py`.
* Commands that are not in the catalogue (new extension commands) can be registered freely;
  add them to the catalogue when they have an abbreviation or a tier other than P0.
* Lookup is exact: full name, then documented abbreviation/alias (L10). Never add prefix matching.

### The `Call` object (`c`)

| member | meaning |
|---|---|
| `c.model` | active `SSIModel` |
| `c.interp` | the `Interpreter` (models, `cwd`, `variables`, `macros`, `session`, `write_options`) |
| `c.int(k, default=None, required=False, what=None)` | argument **k (1-based)** as int; blank → default; `required` blank → 0 + warning (D-PAR-05); reals rounded with a warning (D-PAR-07); non-numeric → `CommandError` |
| `c.float(k, ...)`, `c.str(k, default)`, `c.word(k)` (upper-cased keyword), `c.given(k)`, `c.raw(k)` | other accessors |
| `c.ints(k1, k2)`, `c.floats(k1, k2)`, `c.id_list(k)` | ranges of arguments; `id_list` expands `1-5 7;9` (L8) |
| `c.tokens`, `c.nargs`, `c.rest`, `c.line`, `c.lexed.legacy_paren` | raw lexer output |
| `c.text(k)` | rest-of-line text argument (with `text_from=k`) |
| `c.warn(t)`, `c.error(t)`, `c.info(t)`, `c.confirm(t)` | messages (warnings/errors are prefixed with the command name) |
| `c.fail(t)` / `raise CommandError(t)` | abort the command with an error |
| `raise CommandReported(t)` | abort the command; the error was already reported (no second message) |
| `c.input_path(p)`, `c.output_path(p)` | path resolution: absolute → model path → CWD → calling file (L15) / model path or CWD |

Helpers for ranges `<n1>,[<n2>],[<inc>]` are in `sassi.prep.commands`: `select_ids(ids, a, b, inc)`
and `range_args(c, k, ids, all_default=False)` (spec 08 section 1.4). Pass a dict or set (O(1)
membership for short ranges); for "ids of several tables" use `UnionIds(m.nodes, m.tmass, ...)`
instead of building a new set per command (MUNITS), so INP of a large model stays linear.

Generators validate the numbers they create **before** creating anything: node and element
numbers are >= 1 and a generated (non-blank) element node reference is >= 1 (NGEN, LMOVE,
NMOVE, EGEN, MTGEN, MRGEN), as N and E require.

List commands: FREQ, DAMP and TOPL drop zero values after the first (spec 07 section 1.5). AMP
**keeps** them: with complex SAR (D-INC-10) the values are (Re, Im) pairs with frequent Im = 0, a
SAR of 0 is admissible (Error 118 range [0, 10]) and the count must match the frequency set
(Error 119). Only a first value of 0 clears an AMP list.

### Placeholders

Until a real handler exists, a catalogued command gets a placeholder from its class:

| class | placeholder behaviour |
|---|---|
| `record` | `model.options.records[NAME]` = generic `Record` of all tokens (record replaced); SOIL keeps a legacy `(...)` token and warns (Error 102) |
| `record_keep` | as `record` but blank fields keep the stored value (BINOUT) |
| `indexed` | `model.options.indexed[NAME][key]` = `Record` of all tokens; key = argument 1 (int or upper-case word); DYNP `(label, no)`, BBCP `(num, point)`, EDUOPT upper-case key; SYMM with node1 = 0 deletes the entry |
| `string` | `model.options.strings[NAME]` |
| `action`, `ui` | message only |

P1/P2 placeholders also print `<CMD> is not available in this build (tier Pn)` (requirements
section 0.1); P0 storage placeholders are silent, P0 action placeholders print the message with
tier P0. The list/request classes (FREQ, DAMP, TOPL, AMP, NOUT, EOUT, RDND) have real typed
handlers in `session.py`.

## Model API overview (`sassi.model`)

```python
m = interp.model                       # SSIModel
m.nodes[id] -> Node(x, y, z, csys, fix[6], flags{0..3})   # stored in the defining csys (D-MDL-04)
m.csys[id] -> CoordSys(origin, R, kind 'LOC'|'LOCAL', params / nodes, points)
m.node_global(id); m.global_coordinates(ids) -> (ids, xyz); m.node_in_system(id, s)
m.groups[g] -> Group(type, title, elements{e: Element(nodes, mat, prop, etype, eint, thick, ki, kj)})
m.elements_by_type(); m.element_node_table(type_code) -> (gids, eids, nodes)
m.materials / layers / sections / springs / matrices   # M, L, R, SC, MX (raw values)
m.forces / moments -> NodalLoad(factor[3], arrival[3]); m.tmass / rmass; m.mass_unit(n)
m.nodal_masses(); m.general_matrices(p)                 # mass units (MUNITS, MOPT <matrix>)
m.freq_sets{set: [numbers]}; m.frequency_step(); m.frequencies(set)
m.damp; m.topl; m.amp{motion: [...]}; m.nout; m.eout; m.rdnd
m.options -> OptionStore: record(NAME), entries(NAME), entry(NAME, key), string(NAME)
m.gravity; m.ground_elevation; m.mopt                   # shared variables (HOUSE args 1-2, MOPT)
m.extensions{...}  JSON-able data of other packages (SAVE keeps it; WRITE via a section hook)
m.ui_state{...}    GUI state (hide sets) saved by SAVE but not by WRITE (D-UI-11)
m.to_json() / SSIModel.from_json(); m.copy(); m.canonical(); m.same_state(other); m.model_hash()
```

Material conversions (types 1/2/3), layer constants, beam section formulas and weight → mass are
in `sassi.model.materials`; LOC/LOCAL rotation matrices in `sassi.model.geometry`.

### Typed option records

Records keep the **token text** of every argument (`None` = blank = documented default), so a
record round-trips exactly whether or not it is typed. A typed record declares its fields and is
registered once; `make_record`, `OptionStore.from_json` (RESUME) and the placeholders then create
instances of it automatically:

```python
from sassi.model import Record, Field, register_record_type

@register_record_type
class SiteRecord(Record):
    COMMAND = "SITE"
    FIELDS = (Field("opmode", int, 0), Field("mode1", int, 1), Field("fstep", float, 0.0),
              Field("nl", int, 20), ...)

rec = m.options.record("SITE")     # SiteRecord
rec.fstep; rec.get("nft"); rec.set("nft", 4096); rec.arg(13); rec.is_default()
```

Positional access is 1-based (`rec.arg(1)` = first argument). `GRAVITY` and `GROUNDELEV` set
HOUSE arguments 1 and 2 with `set_arg`, so they work with the generic and the typed record.
`SSIModel.frequency_step()` reads SITE arguments 3, 12, 13 positionally.

## WRITE

`writer.write_pre(model, mdl=False, afwr=False)` returns `(text, notes)`. The output order is
spec 07 section 9.2.47 (D-PAR-17); lists and requests are written with their clear first
(D-PAR-14); numbers use the shortest exact text, so INP of the output gives a deep-equal model
(`SSIModel.same_state`, UT-03) and WRITE of that model gives the same bytes.

* `writer.OPTION_ORDER` — the analysis-option sections and command order.
* `writer.register_writer("SITE", func)` — `func(model) -> [lines]` replaces the generic emission
  of one command (typed writers, e.g. omit defaults). X-commands with a typed record are written
  only when `is_default()` is False (requirements section 3.4.R); a handler should then *not store*
  a default X record, so that the state stays canonical.
* `writer.register_section(title, func)` — a section before the AFWRITE section (e.g. data in
  `model.extensions`).
* Records that are in the store but in no section are written under "Other stored options", so
  nothing is lost.
* Coordinate systems are re-created by the command that defined them, so the round trip is bit
  exact. A LOCAL system whose defining nodes were deleted, moved or redefined afterwards is
  written as: `CSYS,0`, `N` at the stored definition points, `LOCAL`, then the nodes restored
  (or `NDEL`). This happens in the node section, before D/INT/loads.
* Groups: data kept by MTYPE that the current type does not use (8-node elements, EINT 2,
  THICK, KI/KJ) is written under a type that accepts it and followed by `MTYPE,g,<type>`. A
  model with no active group (GDEL of the active group) ends the group section with a
  temporary `GROUP,<max+1>,SOLID` / `GDEL,<max+1>`.

## SAVE / RESUME

`<path>/<name>.sdb` is a `sassi.io.container` file of kind `SDB` holding the JSON of
`SSIModel.to_json()` (complete state: histories, extensions and GUI state included) with
`meta.sdb_version`; RESUME refuses another major version (D-MDL-02). Both need MDL (D-MDL-03).

## Tests

`tests/unit/test_prep_*.py` and `tests/unit/test_model_*.py`. Use `Interpreter(cwd=tmp_path)` and
`ui.run_text(text)` (INP semantics) or `ui.execute(line)` (keyboard semantics); read messages with
`ui.sink.texts(Kind.ERROR)`. Use command names starting with `ZZTEST` for throw-away registrations
and remove them from `registry.REGISTRY` afterwards.
