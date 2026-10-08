/* SASSI-EDU GUI -- the Properties panel (requirements 7.20, D-W6-16): the selected nodes and elements of a
 * model, their properties edited, and the input file that built the model kept in step.
 *
 * Selecting: a click on a node or an element in an element, node or cut plot selects only it; with Shift, Ctrl
 * or Cmd it adds it or removes it; with Alt, a click on an element selects its node nearest to the click
 * (plots.js: SELCLR, NODESEL, ELEMSEL -- command text, rule L17).  "Select" in the panel adds nodes or elements
 * by numbers and ranges.  The panel opens on the right when something is selected (View > Properties shows or
 * hides it at any time).
 *
 * Editing: the fields hold the values of the selection ("mixed" where the selected items differ).  Apply sends
 * the fields that were changed (POST /api/properties): the server runs them as command text -- N, D, INT, MT,
 * MR, MSET, RSET, THICK, ETYPE, EINT, KI, KJ -- and changes the input file that built the model in this session
 * so that it builds the edited model: a node's own N / MT / MR line is rewritten, everything else goes into an
 * edits section after the part of the file that last changed the model.  An open File Editor of that file
 * shows the change and marks the changed lines (with unsaved edits of its own the change goes into its buffer,
 * not to the disk).  The panel then shows the new values.
 */
"use strict";

const SASSI_PROPS = (function () {
  const S = SASSI, el = S.el;
  const R = {};
  const DOFS = ["UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ"];
  const T = {SOLID: 1, BEAMS: 2, SHELL: 3, PLANE: 4, TSHELL: 5, SPRING: 7, GENERAL: 9};
  const ETYPES = [["0", "implicit (by ground elevation)"], ["1", "structure"], ["2", "excavated soil / embedded shell"]];
  const EINTS = [["0", "default (2x2x2 / reduced)"], ["1", "skewed 3x3x3 / selective"], ["2", "distorted 4x4x4 (SOLID)"]];
  let dock = null, data = null, pinned = false, dismissed = "", timer = null, lastNote = null;

  const dockEl = () => dock || (dock = document.getElementById("propdock"));
  const keyOf = (d) => d ? JSON.stringify([d.model, d.nodes.map((n) => n.id), d.groups.map((g) => [g.id, g.elements.map((e) => e.id)])]) : "";
  const count = (d) => d ? d.nodes.length + d.groups.reduce((a, g) => a + g.elements.length, 0) : 0;
  const same = (vals) => vals.every((v) => JSON.stringify(v) === JSON.stringify(vals[0]));
  const fmt = (v) => (v === null || v === undefined ? "" : S.fmt(v));

  /** Ids as ranges: [1, 2, 3, 7] -> "1-3, 7". */
  R.ranges = function (ids) {
    const s = [...new Set(ids.map(Number))].sort((a, b) => a - b), out = [];
    for (let k = 0; k < s.length;) {
      let j = k;
      while (j + 1 < s.length && s[j + 1] === s[j] + 1) j++;
      out.push(j === k ? String(s[k]) : `${s[k]}-${s[j]}`);
      k = j + 1;
    }
    return out.join(", ");
  };
  /** "1-4, 7 9" -> [1, 2, 3, 4, 7, 9] (null when malformed). */
  R.parseIds = function (text) {
    const out = [];
    for (const tok of String(text).split(/[\s,;]+/).filter(Boolean)) {
      const m = /^(\d+)(?:-(\d+))?$/.exec(tok);
      if (!m) return null;
      const a = Number(m[1]), b = m[2] ? Number(m[2]) : a;
      if (b < a || b - a > 100000) return null;
      for (let i = a; i <= b; i++) out.push(i);
    }
    return out;
  };
  /** The command lines that add ids to a selection: ranges, at most 19 or 20 per command. */
  R.selectLines = function (head, ids, per) {
    const toks = R.ranges(ids).split(", ").filter(Boolean), out = [];
    for (let k = 0; k < toks.length; k += per) out.push(`${head},${toks.slice(k, k + per).join(",")}`);
    return out;
  };

  // ---------------------------------------------------------------- loading and showing
  R.refresh = function () { clearTimeout(timer); timer = setTimeout(load, 120); };
  async function load() {
    let d;
    try { d = await S.get("/api/selection"); } catch (e) { return; }
    const k = keyOf(d);
    const changed = k !== keyOf(data);
    data = d;
    if (count(d) && k !== dismissed) show(true);
    else if (!count(d) && !pinned) show(false);
    if (changed) lastNote = null;
    if (!dockEl().hidden) render();
  }
  function show(on) {
    const d = dockEl();
    if (!d || d.hidden === !on) return;
    d.hidden = !on;
    window.dispatchEvent(new Event("resize"));          // the plots follow their new width
    if (S.rebuildMenus) S.rebuildMenus();
  }
  /** View > Properties: show (pinned, also with nothing selected) or hide the panel. */
  R.toggle = function () {
    const on = dockEl().hidden;
    pinned = on;
    dismissed = on ? "" : keyOf(data);
    show(on);
    if (on) { render(); R.refresh(); }
  };
  R.isShown = () => !!dockEl() && !dockEl().hidden;

  // ---------------------------------------------------------------- the panel
  function render() {
    const d = data, box = dockEl();
    if (!box) return;
    const close = el("button", {class: "act-x", title: "close (opens again with the next selection)", text: "×", onclick: () => {
      pinned = false; dismissed = keyOf(data); show(false);
    }});
    const head = el("div", {class: "pp-head"}, el("strong", {text: "Properties"}), el("span", {class: "grow"}), close);
    const parts = [head];
    if (!d || d.model === null) {
      parts.push(el("div", {class: "pp-empty", text: "No model in memory."}));
      box.replaceChildren(...parts);
      return;
    }
    const n = count(d);
    parts.push(el("div", {class: "pp-sub"},
      el("span", {text: `Model ${d.model}${d.name ? " (" + d.name + ")" : ""} · ${n ? n + " selected" : "nothing selected"}`}),
      n ? el("button", {class: "btn small", text: "Clear", title: "SELCLR", onclick: () => S.command("SELCLR")}) : null));
    parts.push(selectRow(d));
    const form = {edits: []};
    if (!n) parts.push(el("div", {class: "pp-empty", text: "Click a node or an element in an element, node or cut plot (Shift, Ctrl or Cmd: add or remove; Alt on an element: its nearest node), or select by numbers above."}));
    if (d.nodes.length) parts.push(nodeSection(d, form));
    for (const g of d.groups) parts.push(groupSection(d, g, form));
    if (n) parts.push(footer(d, form));
    box.replaceChildren(...parts);
  }

  function selectRow(d) {
    const what = el("select", {class: "pp-what"}, el("option", {value: "N", text: "nodes"}), el("option", {value: "E", text: "elements of group"}));
    const grp = el("input", {class: "pp-grp", type: "text", inputmode: "numeric", value: String(d.active_group || 1), title: "group number", hidden: true});
    const ids = el("input", {class: "pp-ids", type: "text", placeholder: "e.g. 1-9, 41", title: "numbers and ranges a-b"});
    const go = el("button", {class: "btn small", text: "Select", title: "add to the selection (NODESEL / ELEMSEL)"});
    what.addEventListener("change", () => { grp.hidden = what.value !== "E"; });
    const run = () => {
      const list = R.parseIds(ids.value);
      if (!list || !list.length) { ids.classList.add("bad"); return; }
      ids.classList.remove("bad");
      let lines;
      if (what.value === "N") {
        const have = new Set(d.nodes.map((x) => x.id));
        lines = R.selectLines("NODESEL", list.filter((i) => !have.has(i)), 20);
      } else {
        const g = Number(grp.value);
        if (!(g >= 1)) { grp.classList.add("bad"); return; }
        const sel = d.groups.find((x) => x.id === g), have = new Set(sel ? sel.elements.map((e) => e.id) : []);
        lines = R.selectLines(`ELEMSEL,${g}`, list.filter((i) => !have.has(i)), 19);
      }
      if (lines.length) S.command(lines);
      ids.value = "";
    };
    go.addEventListener("click", run);
    ids.addEventListener("keydown", (ev) => { if (ev.key === "Enter") run(); });
    return el("div", {class: "pp-select"}, el("span", {class: "pp-k", text: "Select"}), what, grp, ids, go);
  }

  // ---------------------------------------------------------------- fields
  /** A number input holding the common value (blank and "mixed" when they differ); edit(): the new value or
   *  null when unchanged. */
  function numField(vals, opts) {
    const common = same(vals) ? vals[0] : undefined;
    const init = common === undefined ? "" : fmt(common);
    const inp = el("input", {type: "text", inputmode: "decimal", class: "pp-num", value: init, placeholder: common === undefined ? "mixed" : "", title: opts && opts.title});
    inp.addEventListener("input", () => {
      const v = inp.value.trim();
      inp.classList.toggle("bad", v !== "" && !Number.isFinite(Number(v)));
      inp.classList.toggle("chg", v !== init);
    });
    return {el: inp, edit: () => { const v = inp.value.trim(); return v === init || v === "" || !Number.isFinite(Number(v)) ? null : Number(v); }};
  }
  /** A tri-state check box (indeterminate when the selected items differ); edit(): 1 / 0, or null when unchanged. */
  function checkField(vals, label, title) {
    const all1 = vals.every(Boolean), all0 = vals.every((v) => !v);
    const cb = el("input", {type: "checkbox"});
    cb.checked = all1; cb.indeterminate = !all1 && !all0;
    let touched = false;
    cb.addEventListener("change", () => { touched = true; cb.parentNode.classList.toggle("chg", true); });
    return {el: el("label", {class: "pp-check", title: title || ""}, cb, label),
      edit: () => (!touched ? null : cb.checked && !all1 ? 1 : !cb.checked && !all0 ? 0 : null)};
  }
  /** A choice among numbered items (materials, layers, sections); edit(): the new number or null. */
  function choiceField(vals, items, labelOf) {
    const common = same(vals) ? String(vals[0]) : "";
    const sel = el("select", {class: "pp-sel"});
    if (!common) sel.appendChild(el("option", {value: "", text: "mixed"}));
    const known = new Set(items.map((it) => String(it.id)));
    for (const it of items) sel.appendChild(el("option", {value: String(it.id), text: labelOf(it)}));
    for (const v of new Set(vals.map(String))) if (!known.has(v)) sel.appendChild(el("option", {value: v, text: `${v} (not defined)`}));
    sel.value = common;
    sel.addEventListener("change", () => sel.classList.toggle("chg", sel.value !== common));
    return {el: sel, edit: () => (sel.value === "" || sel.value === common ? null : Number(sel.value))};
  }
  const row = (label, ...c) => el("div", {class: "pp-row"}, el("span", {class: "pp-k", text: label}), ...c);
  const chip = (text, title, onclick) => el("button", {class: "pp-chip", text, title, onclick});

  // ---------------------------------------------------------------- nodes
  function nodeSection(d, form) {
    const N = d.nodes, ids = N.map((x) => x.id);
    const sec = el("section", {class: "pp-sec"}, el("h4", {text: `Node${N.length > 1 ? "s" : ""} ${R.ranges(ids)}${N.length > 1 ? ` (${N.length})` : ""}`}));
    const xyz = [0, 1, 2].map((k) => numField(N.map((x) => x.xyz[k])));
    sec.appendChild(row("X, Y, Z", ...xyz.map((f) => f.el)));
    if (N.some((x) => x.csys !== 0)) sec.appendChild(el("div", {class: "pp-note", text: "global coordinates (the nodes were defined in a local system: N redefines them in the global one)"}));
    const fix = DOFS.map((lab, j) => checkField(N.map((x) => x.fix[j]), lab, `fixed ${lab} (D)`));
    sec.appendChild(row("Fixed", el("div", {class: "pp-checks"}, ...fix.map((f) => f.el))));
    const inter = checkField(N.map((x) => x.interaction), "interaction node", "INT code 0: the soil impedance acts here");
    sec.appendChild(row("", inter.el));
    const mt = [0, 1, 2].map((k) => numField(N.map((x) => (x.mass ? x.mass[k] : 0))));
    const mr = [0, 1, 2].map((k) => numField(N.map((x) => (x.rmass ? x.rmass[k] : 0))));
    const units = N.every((x) => x.munits === 1) ? "weight" : N.every((x) => x.munits === 0) ? "mass" : "weight or mass";
    sec.appendChild(row("Mass X, Y, Z", ...mt.map((f) => f.el)));
    sec.appendChild(row("Rot. XX, YY, ZZ", ...mr.map((f) => f.el)));
    sec.appendChild(el("div", {class: "pp-note", text: `MT / MR, as ${units} (MUNITS)`}));
    if (N.length === 1 && N[0].elements.length) {
      sec.appendChild(row("Elements", el("div", {class: "pp-chips"}, ...N[0].elements.map(([g, e]) =>
        chip(`G${g} E${e}`, `select element ${e} of group ${g}`, () => S.command(["SELCLR", `ELEMSEL,${g},${e}`]))))));
    }
    form.edits.push(() => {
      const out = [], nodes = ids;
      const c = xyz.map((f) => f.edit());
      if (c.some((v) => v !== null)) out.push({op: "xyz", nodes, x: c[0], y: c[1], z: c[2]});
      const dofs = {};
      fix.forEach((f, j) => { const v = f.edit(); if (v !== null) dofs[DOFS[j]] = v; });
      if (Object.keys(dofs).length) out.push({op: "fix", nodes, dofs});
      const iv = inter.edit();
      if (iv !== null) out.push({op: "interaction", nodes, value: iv});
      const m1 = mt.map((f) => f.edit()), m2 = mr.map((f) => f.edit());
      if (m1.some((v) => v !== null)) out.push({op: "mass", nodes, values: m1});
      if (m2.some((v) => v !== null)) out.push({op: "rmass", nodes, values: m2});
      return out;
    });
    return sec;
  }

  // ---------------------------------------------------------------- elements
  function groupSection(d, g, form) {
    const E = g.elements, ids = E.map((e) => e.id), t = g.type;
    const sec = el("section", {class: "pp-sec"}, el("h4", {text: `Group ${g.id} ${g.type_name}${g.title ? " · " + g.title : ""} · element${E.length > 1 ? "s" : ""} ${R.ranges(ids)}${E.length > 1 ? ` (${E.length})` : ""}`}));
    const fields = [];
    const add = (attr, f, label) => { fields.push([attr, f]); sec.appendChild(row(label, f.el)); };
    if ([T.SOLID, T.BEAMS, T.SHELL, T.PLANE, T.TSHELL].includes(t)) {
      const soil = (t === T.SOLID || t === T.PLANE) && E.every((e) => e.etype === 2);
      add("mat", soil ? choiceField(E.map((e) => e.mat), d.layers, (it) => `L ${it.id}: ${it.text}`)
        : choiceField(E.map((e) => e.mat), d.materials, (it) => `M ${it.id}: ${it.text}`), soil ? "Soil layer" : "Material");
    }
    if (t === T.BEAMS) add("prop", choiceField(E.map((e) => e.prop), d.sections, (it) => `R ${it.id}: ${it.text}`), "Section");
    if (t === T.SPRING) add("prop", choiceField(E.map((e) => e.prop), d.springs, (it) => `SC ${it.id}: ${it.text}`), "Spring");
    if (t === T.GENERAL) add("prop", numField(E.map((e) => e.prop), {title: "matrix property number"}), "Matrix");
    if (t === T.SHELL || t === T.TSHELL) add("thick", numField(E.map((e) => e.thick)), "Thickness");
    if ([T.SOLID, T.PLANE, T.SHELL, T.TSHELL].includes(t)) add("etype", choiceField(E.map((e) => e.etype), ETYPES.map(([id, text]) => ({id, text})), (it) => `${it.id}: ${it.text}`), "Type (ETYPE)");
    if (t === T.SOLID || t === T.TSHELL) add("eint", choiceField(E.map((e) => e.eint), EINTS.map(([id, text]) => ({id, text})), (it) => `${it.id}: ${it.text}`), "Integration");
    const rel = {};
    if (t === T.BEAMS) {
      for (const [end, key] of [["I", "ki"], ["J", "kj"]]) {
        rel[end] = ["P1", "P2", "P3", "M1", "M2", "M3"].map((lab, j) => checkField(E.map((e) => e[key][j]), lab, `release ${lab} at end ${end} (K${end})`));
        sec.appendChild(row(`Releases ${end}`, el("div", {class: "pp-checks"}, ...rel[end].map((f) => f.el))));
      }
    }
    if (E.length === 1) {
      sec.appendChild(row("Nodes", el("div", {class: "pp-chips"}, ...E[0].nodes.filter((n) => n > 0).map((n) =>
        chip(String(n), `select node ${n}`, () => S.command(["SELCLR", `NODESEL,${n}`]))))));
    }
    form.edits.push(() => {
      const out = [];
      for (const [attr, f] of fields) {
        const v = f.edit();
        if (v !== null) out.push({op: "elem", group: g.id, elements: ids, attr, value: v});
      }
      for (const end of Object.keys(rel)) {
        const vals = rel[end].map((f) => f.edit());
        if (vals.some((v) => v !== null)) {
          // a release edit sets all six codes: unchanged ones keep the common value (mixed ones: per element)
          const key = end === "I" ? "ki" : "kj";
          const groups = new Map();
          for (const e of E) {
            const codes = e[key].map((c, j) => (vals[j] === null ? c : vals[j]));
            const k = codes.join("");
            if (!groups.has(k)) groups.set(k, {codes, ids: []});
            groups.get(k).ids.push(e.id);
          }
          for (const {codes, ids: es} of groups.values()) out.push({op: "release", group: g.id, elements: es, end, codes});
        }
      }
      return out;
    });
    return sec;
  }

  // ---------------------------------------------------------------- apply
  function footer(d, form) {
    const upd = el("input", {type: "checkbox"});
    upd.checked = !!d.input;
    upd.disabled = !d.input;
    const note = el("div", {class: "pp-result"});
    if (lastNote) note.replaceChildren(...lastNote);
    const apply = el("button", {class: "btn primary small", text: "Apply"});
    const revert = el("button", {class: "btn small", text: "Revert", onclick: () => render()});
    apply.addEventListener("click", async () => {
      const edits = form.edits.flatMap((f) => f());
      if (document.querySelector("#propdock .bad")) { note.textContent = "Fix the values marked red first."; return; }
      if (!edits.length) { note.textContent = "Nothing changed."; return; }
      let input = null;
      if (d.input && upd.checked) {
        input = {path: d.input.path};
        const ed = editorOf(d.input.path);
        if (ed && ed.isDirty()) input.text = ed.text();
      }
      apply.disabled = true;
      let r;
      try { r = await S.post("/api/properties", {model: d.model, edits, input}); }
      catch (e) { note.textContent = `Not applied: ${e.message}`; note.classList.add("bad"); apply.disabled = false; return; }
      const parts = [el("div", {text: `${r.lines.length} command${r.lines.length === 1 ? "" : "s"} run${r.ok ? "" : " (with errors: see the Command History)"}.`})];
      if (r.input) {
        const ed = editorOf(r.input.path);
        if (ed) { ed.setText(r.input.text, r.input.saved); ed.flash(r.input.changed); }
        const lines = r.input.changed.map((i) => i + 1);
        parts.push(el("div", {}, `${r.input.name}: ${r.input.inplace ? `${r.input.inplace} line${r.input.inplace > 1 ? "s" : ""} rewritten, ` : ""}${r.input.added} added` +
          (lines.length ? ` (line${lines.length > 1 ? "s" : ""} ${R.ranges(lines)})` : "") + (r.input.saved ? "" : " in the File Editor, not saved yet") + " ",
          el("button", {class: "btn small", text: "Show", onclick: async () => {
            await S.openEditor(r.input.path);
            const e2 = editorOf(r.input.path);
            if (e2) e2.flash(r.input.changed);
          }})));
      }
      lastNote = parts;
      R.refresh();
    });
    const where = d.input
      ? el("label", {class: "pp-check", title: d.input.path}, upd, `also change the input file ${d.input.name}`)
      : el("div", {class: "pp-note", text: "The model was not built from an input file in this session: the changes go to the model only (and to a File Editor connected to Command Entry)."});
    return el("div", {class: "pp-foot"}, el("div", {class: "pp-btns"}, apply, revert), where, note);
  }
  function editorOf(path) {
    const E = S.editors || {};
    if (E[path]) return E[path];
    const base = path.split(/[\\/]/).pop(), same_ = Object.values(E).filter((x) => x.path.split(/[\\/]/).pop() === base);
    return same_.length === 1 ? same_[0] : null;
  }

  setTimeout(() => R.refresh(), 1500);          // a selection kept by the session (a page reload)
  return R;
})();
SASSI.Props = SASSI_PROPS;
