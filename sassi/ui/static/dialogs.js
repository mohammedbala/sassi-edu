/* SASSI-EDU GUI -- dialogs: Options (Model, Write, Check, Analysis with its 12 tabs), Load Model,
 * file picker, converters, Line Selection, Graph Plot Options, Window Options, soil plots,
 * animations, Modules Location / Extension, Check Errors, About, Help, Verification.
 *
 * Every OK that changes the model or a plot submits the equivalent command text (rule L17); the
 * Options dialogs post their values to /api/options/<NAME>, where the server validates them with the
 * CHECK rules (UI-06), builds the command text and executes it through the interpreter. */
"use strict";

(function (S) {
  const el = S.el;
  const D = {};
  S.D = D;
  const stack = [];

  // ================================================================== modal framework
  /** Open a modal dialog.  buttons: [{label, primary, action: async () => true to close}] */
  D.modal = function (opt) {
    const msg = el("span", {class: "msg"});
    const foot = el("div", {class: "modal-foot"}, msg);
    const close = () => {
      back.remove();
      const i = stack.indexOf(api);
      if (i >= 0) stack.splice(i, 1);
      if (opt.onClose) opt.onClose();
    };
    const titleBar = el("div", {class: "modal-title"}, el("span", {text: opt.title}),
      el("button", {class: "x", text: "×", title: "Close", onclick: close}));
    const body = el("div", {class: "modal-body"});
    if (opt.body) body.appendChild(opt.body);
    const box = el("div", {class: "modal", role: "dialog", "aria-label": opt.title, style: opt.width ? {width: opt.width} : null}, titleBar, body, foot);
    const back = el("div", {class: "modal-back"}, box);
    const api = {root: box, body, close, setMessage(text, note) { msg.textContent = text || ""; msg.classList.toggle("note", !!note); }, buttons: []};
    for (const b of opt.buttons || [{label: "Close", primary: true, action: () => true}]) {
      const btn = el("button", {class: "btn" + (b.primary ? " primary" : ""), text: b.label, disabled: b.disabled});
      btn.addEventListener("click", async () => {
        if (btn.disabled) return;
        btn.disabled = true;
        try {
          const r = b.action ? await b.action(api) : true;
          if (r !== false) close();
        } catch (e) {
          api.setMessage(e.message || String(e));
        } finally {
          btn.disabled = false;
        }
      });
      api.buttons.push(btn);
      foot.appendChild(btn);
    }
    // drag by the title bar
    titleBar.addEventListener("mousedown", (ev) => {
      if (ev.target.tagName === "BUTTON") return;
      const r = box.getBoundingClientRect();
      const dx = ev.clientX - r.left, dy = ev.clientY - r.top;
      back.style.alignItems = "flex-start";
      back.style.justifyContent = "flex-start";
      back.style.paddingTop = "0";
      box.style.position = "absolute";
      const move = (e) => { box.style.left = Math.max(0, e.clientX - dx) + "px"; box.style.top = Math.max(0, e.clientY - dy) + "px"; };
      const up = () => { document.removeEventListener("mousemove", move); document.removeEventListener("mouseup", up); };
      document.addEventListener("mousemove", move);
      document.addEventListener("mouseup", up);
      move(ev);
    });
    S.$("#modal-root").appendChild(back);
    stack.push(api);
    const first = body.querySelector("input:not([disabled]), select, textarea");
    if (first) setTimeout(() => first.focus(), 30);
    return api;
  };
  D.closeTop = function () { const t = stack[stack.length - 1]; if (t) t.close(); };
  D.alert = (title, text) => D.modal({title, body: el("div", {style: {whiteSpace: "pre-wrap", overflowWrap: "anywhere", maxWidth: "560px"}, text})});
  /** In-page replacement of window.prompt (which blocks the page): resolves with the text, or null
   *  when the dialog is cancelled or closed. */
  D.ask = function (title, label, value) {
    return new Promise((resolve) => {
      let answer = null;
      const inp = txt(value || "", {"aria-label": label});
      const api = D.modal({title, body: el("div", {class: "dlg", style: {width: "440px"}}, row(label, inp)), onClose: () => resolve(answer),
        buttons: [{label: "OK", primary: true, action: () => { answer = inp.value; return true; }}, {label: "Cancel", action: () => true}]});
      inp.addEventListener("keydown", (ev) => { if (ev.key === "Enter") { ev.preventDefault(); api.buttons[0].click(); } });
    });
  };
  /** In-page replacement of window.confirm: resolves true (first button) or false (second button,
   *  close box, Escape). */
  D.confirm = function (title, text, yes, no) {
    return new Promise((resolve) => {
      let answer = false;
      D.modal({title, body: el("div", {style: {whiteSpace: "pre-wrap", overflowWrap: "anywhere", width: "480px", maxWidth: "100%"}, text}), onClose: () => resolve(answer),
        buttons: [{label: yes || "Yes", primary: true, action: () => { answer = true; return true; }}, {label: no || "No", action: () => true}]});
    });
  };
  /** Accessible names for the inputs of a dialog row (the label text). */
  function nameInputs(node, label) {
    if (!label) return node;
    node.querySelectorAll("input:not([type=radio]):not([type=checkbox]), select, textarea").forEach((x) => {
      if (!x.getAttribute("aria-label")) x.setAttribute("aria-label", label);
    });
    return node;
  }
  function row(label, ...inputs) { return nameInputs(el("div", {class: "row"}, el("label", {class: "lab", text: label}), ...inputs), label); }
  function txt(value, attrs) { return el("input", Object.assign({type: "text", value: value === undefined || value === null ? "" : String(value)}, attrs || {})); }
  function chk(on, attrs) { const c = el("input", Object.assign({type: "checkbox"}, attrs || {})); c.checked = !!on; return c; }
  function fieldset(legend, ...kids) { return el("fieldset", {}, el("legend", {text: legend}), ...kids); }
  /** The value of the choice a <select> shows (numbers stay numbers, words stay words: REL, FILE8 ...). */
  function choiceValue(choices, text) {
    const c = (choices || []).find(([v]) => String(v) === String(text));
    return c ? c[0] : text;
  }
  /** A path or free-text command token: double-quoted when it holds a comma (D-PAR-06), so a folder
   *  such as "Smith, J" stays one argument and the Command History replays the same command. */
  const q = (v) => {
    const s = String(v === undefined || v === null ? "" : v).trim();
    return s.includes(",") ? `"${s}"` : s;
  };
  S.q = q;
  async function cmdOk(lines, api) {
    const r = await S.command(lines);
    const bad = r.results.filter((x) => !x.ok);
    if (bad.length) {
      const errs = r.messages.filter((m) => m.kind === "ERROR").map((m) => m.text);
      if (api) api.setMessage(errs.join("\n") || `${bad[0].line} failed`);
      return false;
    }
    return true;
  }

  // ================================================================== file picker
  /** A file dialog on the model directory / working directory (path safety on the server).
   *  opt: {title, filter: ['.pre'], create, dir, onOk(path), folder (choose a directory)} */
  D.pickFile = function (opt) {
    const list = el("div", {class: "list"});
    const name = txt(opt.initial || "", {style: {flex: "1"}, placeholder: opt.create ? "file name (a new name creates the file)" : "file name"});
    const dirLab = el("span", {class: "grow", style: {fontFamily: "var(--mono)", fontSize: "12px", overflow: "hidden", textOverflow: "ellipsis"}});
    const showAll = chk(!opt.filter, {});
    const cwd = txt(S.state ? S.state.cwd : "", {style: {flex: "1"}});
    let cur = null, selected = null;
    const load = async (dir) => {
      let d;
      try { d = await S.get("/api/files" + (dir ? `?dir=${encodeURIComponent(dir)}` : "")); } catch (e) { api.setMessage(e.message); return; }
      cur = d;
      dirLab.textContent = d.dir;
      list.innerHTML = "";
      const add = (label, size, onclick, ondbl, cls) => {
        const it = el("div", {class: cls || ""}, el("span", {text: label}), el("span", {class: "sz", text: size}));
        it.addEventListener("click", () => { list.querySelectorAll(".sel").forEach((x) => x.classList.remove("sel")); it.classList.add("sel"); onclick(); });
        if (ondbl) it.addEventListener("dblclick", ondbl);
        list.appendChild(it);
      };
      if (d.parent) add("../", "", () => {}, () => load(d.parent));
      for (const sd of d.dirs) add(sd.name + "/", "", () => { selected = sd.path; if (opt.folder) name.value = sd.path; }, () => load(sd.path));
      const flt = (f) => showAll.checked || !opt.filter || opt.filter.some((e) => f.name.toLowerCase().endsWith(e.toLowerCase()));
      if (!opt.folder) for (const f of d.files.filter(flt)) add(f.name, f.size.toLocaleString(), () => { selected = f.path; name.value = f.name; }, () => { selected = f.path; name.value = f.name; okBtn.click(); });
      api.setMessage(d.roots ? `allowed: ${d.roots.join("  |  ")}` : "", true);
    };
    showAll.addEventListener("change", () => load(cur && cur.dir));
    const body = el("div", {class: "picker"},
      el("div", {class: "row"}, el("label", {text: "Folder"}), dirLab),
      list,
      el("div", {class: "row"}, el("label", {class: "lab", text: opt.folder ? "Folder" : "File name"}), name),
      opt.filter ? el("div", {class: "row"}, el("label", {class: "chk"}, showAll, " show all files")) : null,
      el("div", {class: "row"}, el("label", {class: "lab", text: "Working directory"}), cwd,
        el("button", {class: "btn small", text: "CD", title: "CD,<dir>: change the working directory (command)", onclick: async () => {
          if (await cmdOk(`CD,${q(cwd.value)}`, api)) { await S.refreshState(); load(null); }
        }})));
    const resolve = () => {
      const v = name.value.trim();
      if (!v) return null;
      if (v.startsWith("/") || /^[A-Za-z]:[\\/]/.test(v) || v.startsWith("~")) return v;
      return (cur ? cur.dir.replace(/\/$/, "") + "/" : "") + v;
    };
    const api = D.modal({title: opt.title || "Select File", body, buttons: [
      {label: "OK", primary: true, action: async () => {
        const p = opt.folder ? (name.value.trim() || (cur && cur.dir)) : resolve();
        if (!p) { api.setMessage("choose a file"); return false; }
        const r = opt.onOk ? await opt.onOk(p) : true;
        return r !== false;
      }},
      {label: "Cancel", action: () => true}]});
    const okBtn = api.buttons[0];
    name.addEventListener("keydown", (ev) => { if (ev.key === "Enter") { ev.preventDefault(); okBtn.click(); } });
    load(opt.dir || null);
    return api;
  };

  // ================================================================== Model menu
  D.loadModel = async function () {
    let db;
    try { db = await S.get("/api/db"); } catch (e) { S.local("ERROR", e.message); return; }
    let sel = null;
    const tree = el("div", {class: "tree"});
    const info = el("div", {class: "kv", style: {marginTop: "8px", minHeight: "48px"}});
    const draw = () => {
      tree.innerHTML = "";
      if (!db.groups.length) tree.appendChild(el("div", {class: "empty", text: "empty: use Add Group, then Add Model"}));
      for (const g of db.groups) {
        const gi = el("div", {class: "g" + (sel && sel.type === "g" && sel.g === g.name ? " sel" : ""), text: "▾ " + g.name});
        gi.addEventListener("click", () => { sel = {type: "g", g: g.name}; draw(); });
        tree.appendChild(gi);
        for (const m of g.models) {
          const mi = el("div", {class: "m" + (sel && sel.type === "m" && sel.g === g.name && sel.m.name === m.name && sel.m.path === m.path ? " sel" : ""), text: m.name});
          mi.addEventListener("click", () => { sel = {type: "m", g: g.name, m}; draw(); });
          mi.addEventListener("dblclick", () => { sel = {type: "m", g: g.name, m}; api.buttons[0].click(); });
          tree.appendChild(mi);
        }
      }
      info.innerHTML = "";
      if (sel && sel.type === "g") info.append(el("span", {class: "k", text: "Group Name:"}), el("span", {text: sel.g}));
      if (sel && sel.type === "m") info.append(el("span", {class: "k", text: "Model:"}), el("span", {text: sel.m.name}),
        el("span", {class: "k", text: "Location:"}), el("span", {text: sel.m.path}), el("span", {class: "k", text: "Title:"}), el("span", {text: sel.m.title || ""}));
    };
    const post = async (body) => { db = await S.post("/api/db", body); draw(); };
    const body = el("div", {class: "dlg", style: {display: "flex", flexWrap: "wrap", gap: "12px", width: "620px"}},
      el("div", {style: {flex: "1 1 300px", minWidth: "0"}}, tree),
      el("div", {style: {display: "flex", flexDirection: "column", gap: "6px", flex: "0 1 200px", minWidth: "160px"}},
        el("button", {class: "btn", text: "Add Group", onclick: async () => {
          const n = await D.ask("Add Group", "Group name");
          if (n && n.trim()) await post({action: "add_group", group: n.trim()}).catch((e) => api.setMessage(e.message));
        }}),
        el("button", {class: "btn", text: "Remove Group", onclick: async () => {
          if (!sel || sel.type !== "g") { api.setMessage("highlight a group"); return; }
          if (!(await D.confirm("Remove Group", `Remove group ${sel.g} from the database?`))) return;
          const del = await D.confirm("Remove Group", `Also delete the files of the models of ${sel.g} from disk?\n(only files named after each model: <name>.*, <name>_*)`, "Delete files", "Keep files");
          await post({action: "remove_group", group: sel.g, delete_files: del}).catch((e) => api.setMessage(e.message));
          sel = null; draw();
        }}),
        el("button", {class: "btn", text: "Add Model", onclick: () => {
          if (!sel) { api.setMessage("highlight a group (or a model of the group)"); return; }
          D.addModel(sel.g, (b) => post(Object.assign({action: "add_model"}, b)));
        }}),
        el("button", {class: "btn", text: "Remove Model", onclick: async () => {
          if (!sel || sel.type !== "m") { api.setMessage("highlight a model"); return; }
          if (!(await D.confirm("Remove Model", `Remove model ${sel.m.name} from group ${sel.g}?`))) return;
          const del = await D.confirm("Remove Model", `Also delete the files of model ${sel.m.name} from disk (${sel.m.path}/${sel.m.name}.*)?`, "Delete files", "Keep files");
          await post({action: "remove_model", group: sel.g, name: sel.m.name, path: sel.m.path, delete_files: del}).catch((e) => api.setMessage(e.message));
          sel = null; draw();
        }}),
        info));
    const api = D.modal({title: "Load Model", body, buttons: [
      {label: "Open", primary: true, action: async () => {
        if (!sel || sel.type !== "m") { api.setMessage("highlight a model"); return false; }
        const r = await S.post("/api/db", {action: "open", name: sel.m.name, path: sel.m.path, title: sel.m.title});
        await S.refreshState();
        return r.ok !== false;
      }},
      {label: "Cancel", action: () => true}]});
    draw();
    api.setMessage(db.path, true);
  };
  D.addModel = function (group, done) {
    const name = txt(""), path = txt(S.state && S.state.cwd ? S.state.cwd : ""), title = txt("");
    const api = D.modal({title: "Add Model", body: el("div", {style: {width: "460px"}},
      row("Group", el("span", {text: group})), row("Model name", name), row("Model directory", path), row("Title", title),
      el("div", {class: "opt-note", text: "The directory is created if missing; the model files appear at the first SAVE."})), buttons: [
      {label: "OK", primary: true, action: async () => { await done({group, name: name.value, path: path.value, title: title.value}); return true; }},
      {label: "Cancel", action: () => true}]});
    return api;
  };
  /** Model > Converters (spec 09 section 6.4): "SASSI .hou to .pre Converter" -> CONVERT,SSI,<model>,
   *  <file.hou>,[<prefile>]; "ANSYS .cdb to .pre Converter" -> CONVERT,ANSYS,<model>,<file.cdb>,<gravity>,
   *  [<prefile>],[<damp>].  A blank model number converts into the active model; the dialog proposes
   *  the lowest unused number so that the active model is not replaced by accident. */
  D.converter = function (kind) {
    const ansys = kind === "ANSYS";
    const model = txt(S.lowestUnusedModel ? S.lowestUnusedModel() : 1, {class: "num", placeholder: "active"});
    const file = txt(""), pre = txt("", {placeholder: "optional: the converted model as a .pre file"});
    const grav = txt("32.2", {class: "num"}), damp = txt("", {class: "num", placeholder: "materials"});
    const ext = ansys ? [".cdb"] : [".hou"];
    const body = el("div", {class: "dlg", style: {width: "560px"}},
      row("Input File Name", file, el("button", {class: "btn small", text: "<<", title: "browse", onclick: () => D.pickFile({title: "Input file", filter: ext, onOk: (p) => { file.value = p; }})})),
      row("Output .pre File Name", pre),
      row("Save Converted Data to Model Number", model),
      ansys ? row("Enter Value for Gravity", grav, el("span", {class: "opt-note", text: "32.2 ft/s², 386.4 in/s², 9.81 m/s²"})) : null,
      ansys ? row("Damping Ratio (extension)", damp, el("span", {class: "opt-note", text: "blank: from the materials"})) : null,
      el("div", {class: "opt-note", text: (ansys ? "CONVERT,ANSYS,<model>,<file.cdb>,<gravity>,[<prefile>],[<damp>]" : "CONVERT,SSI,<model>,<file.hou>,[<prefile>]: reads the .sit / .poi decks alongside") +
        ". Blank model number = the active model."}),
      el("div", {class: "note", style: {color: "var(--warn)"}, text: "The converters have had limited testing: check every converted model (and read the warnings) before an analysis."}));
    const api = D.modal({title: ansys ? "ANSYS .cdb to .pre Converter" : "SASSI .hou to .pre Converter", body, buttons: [
      {label: "Convert", primary: true, action: (a) => {
        if (!file.value.trim()) { a.setMessage("Input File Name is required"); return false; }
        if (ansys && !(Number(grav.value) > 0)) { a.setMessage("Enter Value for Gravity: a positive number is required"); return false; }
        const args = ansys ? ["ANSYS", model.value.trim(), q(file.value), grav.value.trim(), q(pre.value), damp.value.trim()] : ["SSI", model.value.trim(), q(file.value), q(pre.value)];
        while (args.length && args[args.length - 1] === "") args.pop();
        return cmdOk(`CONVERT,${args.join(",")}`, a);
      }},
      {label: "Cancel", action: () => true}]});
    return api;
  };
  D.output = function () {
    const m = S.state.model || {};
    const file = txt(m.name ? m.name + ".pre" : "model.pre"), dir = txt(m.path || S.state.cwd);
    return D.modal({title: "Output (WRITE)", body: el("div", {style: {width: "520px"}}, row("File name", file),
      row("Path", dir, el("button", {class: "btn small", text: "<<", onclick: () => D.pickFile({title: "Folder", folder: true, onOk: (p) => { dir.value = p; }})})),
      el("div", {class: "opt-note", text: "WRITE,<file>,<path>: the model as a .pre command file (Options > Write adds MDL / AFWRITE lines)"})),
      buttons: [{label: "OK", primary: true, action: (a) => cmdOk(`WRITE,${q(file.value)},${q(dir.value)}`, a)}, {label: "Cancel", action: () => true}]});
  };
  D.exportAnsys = function () {
    const m = S.state.model || {};
    const file = txt(m.name ? m.name + ".inp" : ""), dir = txt(m.path || "");
    return D.modal({title: "Export to ANSYS", body: el("div", {style: {width: "520px"}}, row("File name", file), row("Directory", dir),
      el("div", {class: "opt-note", text: "ANSYS,[FileName],[Dir]: APDL input; run ANSYSREFORMAT first for beam end releases (manual WARNING)"})),
      buttons: [{label: "OK", primary: true, action: (a) => cmdOk(`ANSYS,${q(file.value)},${q(dir.value)}`, a)}, {label: "Cancel", action: () => true}]});
  };
  D.exit = async function () {
    let info = {unsaved: []};
    try { info = await S.get("/api/exit"); } catch (e) { /* ignore */ }
    let text = "Exit SASSI-EDU?  SASSIini.xml will be saved; models and plots are not saved.";
    if (info.unsaved.length) text += "\n\nModels changed since their last SAVE:\n" + info.unsaved.map((u) => `  Model ${u.number} ${u.name || "(no name)"}`).join("\n");
    if (!(await D.confirm("Exit", text, "Exit", "Cancel"))) return;
    try { await S.post("/api/exit"); } catch (e) { /* the server stops */ }
    document.body.innerHTML = '<div class="empty" style="margin-top:20vh;font-size:16px">SASSI-EDU has stopped. You can close this tab.</div>';
  };

  // ================================================================== File menu
  D.exportImage = function () {
    const p = S.P.activePlot();
    if (!p) { S.local("WARNING", "Export Image: no active plot"); return; }
    if (S.web) return exportImageWeb(p);
    const name = txt(`${(p.caption || "plot").replace(/[^\w.-]+/g, "_")}.bmp`);
    return D.modal({title: "Export Image", body: el("div", {style: {width: "460px"}}, row("File name", name),
      el("div", {class: "opt-note", text: "CAPTUREPLOT,<file>: PNG when the name contains .png, otherwise BMP (D-UI-10)"})),
      buttons: [{label: "OK", primary: true, action: (a) => cmdOk(`CAPTUREPLOT,${q(name.value)}`, a)}, {label: "Cancel", action: () => true}]});
  };
  /** Browser version: CAPTUREPLOT draws with matplotlib, which the page does not load -- the plot of the
   *  active tab is saved as a PNG by Plotly instead (the camera icon of the plot does the same). */
  function exportImageWeb(p) {
    const tab = S.tabs.find((t) => t.plotId === p.id);
    const gd = tab && tab.pane.querySelector(".js-plotly-plot");
    if (!gd || !window.Plotly) { S.local("WARNING", "Export Image: the active plot is not drawn yet"); return; }
    S.selectTab(tab.id);
    const name = txt(`${(p.caption || "plot").replace(/[^\w.-]+/g, "_")}.png`);
    return D.modal({title: "Export Image", body: el("div", {style: {width: "460px"}}, row("File name", name),
      el("div", {class: "opt-note", text: "Saves the plot of the active tab as a PNG file on your computer (browser download)."})),
      buttons: [{label: "OK", primary: true, action: async () => {
        const base = name.value.trim().replace(/\.png$/i, "") || "plot";
        await window.Plotly.downloadImage(gd, {format: "png", filename: base});
        S.status(`image ${base}.png saved by the browser`);
        return true;
      }}, {label: "Cancel", action: () => true}]});
  }
  // ================================================================== File > Upload / Download (browser version)
  /** Save text as a file on the visitor's computer (browser download). */
  function saveText(name, text) {
    const url = URL.createObjectURL(new Blob([text], {type: "text/plain;charset=utf-8"}));
    const a = el("a", {href: url, download: name, style: {display: "none"}});
    document.body.appendChild(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 10000);
  }
  const WEB_NOTE = "The workspace lives in the memory of this browser tab: a reload or closing the tab starts afresh " +
    "(your course progress is kept). Download the files you want to keep.";
  /** File > Upload to Workspace... (browser version): text files from the visitor's computer into a folder of the
   *  workspace (POST /api/file); a binary file is refused (the File Editor and the API carry text). */
  D.uploadFiles = function () {
    const input = el("input", {type: "file", multiple: true, "aria-label": "files to upload"});
    const dir = txt((S.state && (S.state.model.path || S.state.cwd)) || "", {style: {flex: "1"}});
    const api = D.modal({title: "Upload to Workspace", body: el("div", {style: {width: "560px"}},
      row("Files", input),
      row("Folder", dir, el("button", {class: "btn small", text: "<<", title: "choose a folder of the workspace", onclick: () => D.pickFile({title: "Folder", folder: true, dir: dir.value, onOk: (p) => { dir.value = p; }})})),
      el("div", {class: "opt-note", text: "Text files only (.pre, decks, motions, spectra ...; up to 8 MB each). A file of the same name is replaced. " + WEB_NOTE})),
      buttons: [{label: "Upload", primary: true, action: async (a) => {
        const files = Array.from(input.files || []);
        if (!files.length) { a.setMessage("choose one or more files"); return false; }
        const folder = dir.value.trim().replace(/\/+$/, "");
        const done = [], bad = [];
        for (const f of files) {
          const bytes = new Uint8Array(await f.arrayBuffer());
          if (bytes.length > 8 * 1024 * 1024) { bad.push(`${f.name}: larger than 8 MB`); continue; }
          if (bytes.subarray(0, 4096).includes(0)) { bad.push(`${f.name}: a binary file (text files only)`); continue; }
          const text = new TextDecoder("utf-8").decode(bytes);
          try {
            const r = await S.post("/api/file", {name: `${folder}/${f.name}`, text});
            done.push(r.path);
          } catch (e) { bad.push(`${f.name}: ${e.message}`); }
        }
        for (const d of done) S.local("CONFIRM", `Upload to Workspace: ${d}`);
        for (const b of bad) S.local("WARNING", `Upload to Workspace: ${b}`);
        if (done.length && S.P && S.P.refreshResults) S.P.refreshResults();
        if (bad.length) { a.setMessage(bad.join("\n")); return false; }
        return true;
      }}, {label: "Cancel", action: () => true}]});
    return api;
  };
  /** File > Download... (browser version): a text file of the workspace (GET /api/file) saved on the visitor's
   *  computer. */
  D.downloadFile = function () {
    return D.pickFile({title: "Download (a text file of the workspace)", dir: S.state && S.state.model.path || undefined, onOk: async (p) => {
      let d;
      try { d = await S.get(`/api/file?name=${encodeURIComponent(p)}`); } catch (e) { S.local("ERROR", `Download ${p}: ${e.message}`); return false; }
      saveText(d.name, d.text);
      S.status(`${d.name} saved by the browser`);
      return true;
    }});
  };

  D.exportTable = function () {
    const p = S.P.activePlot();
    if (!p) return;
    const name = txt(`plot${String(p.id).padStart(2, "0")}.csv`);
    return D.modal({title: "Export Table", body: el("div", {style: {width: "460px"}}, row("CSV file", name),
      el("div", {class: "opt-note", text: "x column, one column per line on the union grid (D-UI-17)"})),
      buttons: [{label: "OK", primary: true, action: async (a) => {
        try { const r = await S.post("/api/export_table", {name: name.value}); S.status(`table written to ${r.path}`); return true; } catch (e) { a.setMessage(e.message); return false; }
      }}, {label: "Cancel", action: () => true}]});
  };

  // ================================================================== Plot menu
  D.cutPlot = async function () {
    let cuts = {cuts: []};
    try { cuts = await S.get("/api/cuts"); } catch (e) { /* none */ }
    const cut = txt(cuts.cuts.length ? cuts.cuts[0].cut : 1, {class: "num"});
    const model = txt(S.state.active, {class: "num"});
    D.modal({title: "Select Cut to Display", body: el("div", {style: {width: "380px"}}, row("Cut Number", cut), row("Model Number", model),
      el("div", {class: "opt-note", text: cuts.cuts.length ? "defined cuts: " + cuts.cuts.map((c) => `${c.cut} (${c.elements} elements)`).join(", ") : "no cut defined yet (CUTADD / CUTVOL / SLICE)"})),
      buttons: [{label: "Ok", primary: true, action: (a) => cmdOk(`CUTPLOT,${cut.value},${model.value}`, a)}, {label: "Cancel", action: () => true}]});
  };

  /** Line Selection dialog (spec 06 section 5.2): Plot > Spectrum TFU-TFI / Time History, the toolbar,
   *  Results > Plot and SPECPLOT / THPLOT typed without arguments (D-UI-08).
   *
   *  The dialog starts from the plot it would extend: the active plot when it is of the same kind,
   *  else the most recent plot of that kind (Results > Plot of a file starts empty).  The lines of
   *  that plot are checked and its title, axis labels, log axes, minor ticks and ranges are filled
   *  in.  Lines loaded with Add Line(s) are checked.  Ok submits SPECPLOT / THPLOT with the checked
   *  lines -- a new plot, as in ACS SASSI ("Replace" closes the plot the dialog started from) -- and
   *  the settings commands that make the new plot look as shown: a new plot starts from the session
   *  defaults, so only the differences are submitted (L17). */
  D.lineSelection = async function (kind, preset) {
    kind = kind === "THPLOT" ? "THPLOT" : "SPECPLOT";
    preset = preset || {};
    const isSpec = kind === "SPECPLOT";
    // the plot state (as S.state.plots, fetched fresh): the active plot, every open plot, the 2D defaults
    let plots = (S.state && S.state.plots) || {};
    try { plots = await S.get("/api/plots"); } catch (e) { /* keep the state of the last event */ }
    const open = plots.plots || [];
    const active = open.find((p) => p.id === plots.active) || null;
    const src = preset.file ? null : (active && active.kind === kind ? active : open.filter((p) => p.kind === kind).pop() || null);
    const st = (src ? src.settings : plots.defaults2d) || {};
    const checked = new Set(src ? (src.params.lines || []) : []);
    let lines = {};
    let startTouched = false;        // a Starting Number typed by the user is kept; otherwise the next free number
    let fileAdded = false;           // the file of the Input Line File field was loaded (Ok does not load it again)
    const sval = (v) => (v === null || v === undefined ? "" : String(v));
    const listBox = el("div", {class: "checklist"}, el("div", {class: "empty", text: "loading the lines in memory ..."}));
    const title = txt(src ? src.title || "" : ""), xl = txt(st.xtitle || ""), yl = txt(st.ytitle || "");
    const num = (v, label) => txt(sval(v), {class: "num", "aria-label": label});
    const xmin = num(st.xmin, "X Min"), xmax = num(st.xmax, "X Max"), ymin = num(st.ymin, "Y Min"), ymax = num(st.ymax, "Y Max");
    const xlog = chk(st.log_x), xtick = chk(st.minor_x), ylog = chk(st.log_y), ytick = chk(st.minor_y);
    const file = txt(preset.file || ""), start = txt("1", {class: "num"}), nlines = txt("1", {class: "num"});
    const pair = el("select", {}, el("option", {value: "0", text: "0: dt first, then values (ACS SASSI)"}), el("option", {value: "1", text: "1: time / value pairs"}));
    const replace = chk(false);
    start.addEventListener("input", () => { startTouched = true; });
    const nextFree = () => { const ns = Object.keys(lines).map(Number); return (ns.length ? Math.max(...ns) : 0) + 1; };
    const refresh = async () => {
      try { lines = (await S.get("/api/lines")).lines; } catch (e) { lines = {}; }
      listBox.innerHTML = "";
      const nums = Object.keys(lines).map(Number).sort((a, b) => a - b);
      if (!nums.length) listBox.appendChild(el("div", {class: "empty", text: "no line in memory: load a file below (Add Line(s))"}));
      for (const n of nums) {
        const L = lines[n];
        const c = chk(checked.has(n), {"aria-label": `line ${n}`});
        c.addEventListener("change", () => { if (c.checked) checked.add(n); else checked.delete(n); });
        listBox.appendChild(el("label", {}, c, `${n}: ${L.name}`, el("span", {class: "chip", text: `${L.kind}, ${L.n} pts`})));
      }
      if (!startTouched) start.value = String(nextFree());
    };
    D._lineDialogRefresh = refresh;
    let resolved = null;       // absolute path found by the server (READSPEC then finds it whatever the CWD)
    const inspect = async () => {
      resolved = null;
      fileAdded = false;
      if (!file.value.trim()) return;
      try {
        const fi = await S.get(`/api/fileinfo?name=${encodeURIComponent(file.value.trim())}`);
        resolved = fi.path;
        // transfer functions hold "f amp phase": Lines in file 1 loads the amplitude; other files all columns
        if (isSpec) nlines.value = String(/\.tf[uid]$/i.test(fi.name) ? 1 : Math.max(1, fi.columns || 1));
        else pair.value = String(fi.pair || 0);
        api.setMessage(`${fi.name}: ${fi.label || fi.kind}, ${fi.columns} data column(s)` +
          (isSpec && /\.tf[uid]$/i.test(fi.name) ? " (amplitude, phase)" : "") + `, ${fi.points} points`, true);
      } catch (e) { api.setMessage(e.message); }
    };
    file.addEventListener("change", inspect);
    const addLines = async () => {
      if (!file.value.trim()) { api.setMessage("choose a file"); return false; }
      if (!resolved) await inspect();
      const fname = resolved || file.value.trim();
      // the target line numbers are taken BEFORE the command runs: READSPEC's 'lines' event refreshes
      // this list (and the next free Starting Number) while the request is still in flight
      const s0 = Number(start.value), n = isSpec ? Number(nlines.value) : 1;
      if (!Number.isInteger(s0) || s0 < 1 || !Number.isInteger(n) || n < 1) { api.setMessage("Starting Number and Lines in file must be positive integers"); return false; }
      const line = isSpec ? `READSPEC,${q(fname)},${n},${s0}` : `READTH,${q(fname)},${pair.value},${s0}`;
      const r = await S.command(line);
      if (!r.ok) {
        api.setMessage(r.messages.filter((m) => m.kind === "ERROR").map((m) => m.text).join("\n") || `${line} failed`);
        return false;
      }
      for (let k = 0; k < n; k++) checked.add(s0 + k);           // new lines are checked
      fileAdded = true;
      startTouched = false;
      await refresh();
      api.setMessage(`${line}: line(s) ${s0}${n > 1 ? `-${s0 + n - 1}` : ""} added and checked`, true);
      return true;
    };
    const axisBox = (label, mn, mx, lg, tk) => fieldset(label,
      el("div", {class: "minmax"}, el("label", {text: "Min"}), mn, el("label", {text: "Max"}), mx),
      el("div", {class: "row"}, el("label", {class: "chk"}, lg, " Logarithmic"), el("label", {class: "chk"}, tk, " Show Ticks")));
    const filter = isSpec ? [".tfu", ".tfi", ".tfd", ".rs", ".rso", ".rsi", ".psd", ".fft", ".txt", ".csv"] : [".acc", ".thd", ".ths", ".th", ".vel", ".dis", ".txt"];
    const body = el("div", {class: "dlg", style: {width: "660px"}},
      fieldset("Lines in memory (check the lines to plot)", listBox),
      el("div", {class: "dlg-cols"},
        fieldset("Title / Axis Labels", row("Title", title), row("X-Label", xl), row("Y-Label", yl)),
        el("div", {}, axisBox("X Axis Options", xmin, xmax, xlog, xtick), axisBox("Y Axis Options", ymin, ymax, ylog, ytick))),
      fieldset("Input Line File",
        row("File Name", file, el("button", {class: "btn small", text: "<<", title: "browse", onclick: () => D.pickFile({title: "Line file", filter,
          onOk: (p) => { file.value = p; inspect(); }})})),
        isSpec ? row("Starting Number", start, el("label", {text: "Lines in file"}), nlines) : row("Line Number", start, el("label", {text: "Pair"}), pair),
        el("div", {class: "row"}, el("button", {class: "btn", text: "Add Line(s)", onclick: addLines}),
          el("span", {class: "opt-note", text: isSpec ? "READSPEC: the frequency column is not counted; columns 1..n go to lines start, start+1, ..." : "READTH: one history per line"}))),
      src ? el("div", {class: "row"}, el("label", {class: "chk", title: "close the plot first (ACTIVATEPLOT + CLOSEPLOT); otherwise Ok opens a new plot"},
        replace, ` Replace plot ${src.id} (${src.caption}${src.title ? " - " + src.title : ""})`)) : null,
      el("div", {class: "opt-note", text: src ? `Checked lines, titles and axes are those of plot ${src.id}.` : "Ok opens a new plot with the checked lines."}));
    const validRange = () => {
      const v = [xmin, xmax, ymin, ymax].map((x) => x.value.trim());
      for (const [k, x] of v.entries()) if (x !== "" && !isFinite(Number(x))) return `${["X Min", "X Max", "Y Min", "Y Max"][k]}: '${x}' is not a number`;
      for (const [lo, hi, lg, ax] of [[v[0], v[1], xlog.checked, "X"], [v[2], v[3], ylog.checked, "Y"]]) {
        if (lo !== "" && hi !== "" && Number(lo) >= Number(hi)) return `${ax} Min must be smaller than ${ax} Max`;
        if (lg && lo !== "" && Number(lo) <= 0) return `${ax} axis is logarithmic: Min must be > 0`;
      }
      return null;
    };
    const api = D.modal({title: "Line Selection", body, buttons: [
      {label: "Ok", primary: true, action: async (a) => {
        const bad = validRange();
        if (bad) { a.setMessage(bad); return false; }
        if (file.value.trim() && !fileAdded) { if (!(await addLines())) return false; }
        const nums = [...checked].filter((n) => lines[n]).sort((x, y) => x - y);
        if (!nums.length) { a.setMessage("check at least one line"); return false; }
        if (nums.length > 50) { a.setMessage(`${nums.length} lines checked: a plot shows at most 50`); return false; }
        let d2 = plots.defaults2d || {};
        try { d2 = (await S.get("/api/plots")).defaults2d || d2; } catch (e) { /* last known defaults */ }
        const cmds = [];
        if (replace.checked && src) {
          if (src.id !== plots.active) cmds.push(`ACTIVATEPLOT,${src.id}`);
          cmds.push("CLOSEPLOT");
        }
        cmds.push(`${kind},${nums.join(",")}`);
        if (title.value.trim()) cmds.push(`PLOTTITLE,${title.value.trim()}`);
        if (xl.value.trim() !== (d2.xtitle || "")) cmds.push(`XTITLE,${xl.value.trim()}`);
        if (yl.value.trim() !== (d2.ytitle || "")) cmds.push(`YTITLE,${yl.value.trim()}`);
        const ax = [st.major_x === false ? 0 : 1, st.major_y === false ? 0 : 1, +xtick.checked, +ytick.checked, +xlog.checked, +ylog.checked];
        const ax0 = [d2.major_x === false ? 0 : 1, d2.major_y === false ? 0 : 1, +!!d2.minor_x, +!!d2.minor_y, +!!d2.log_x, +!!d2.log_y];
        if (ax.join() !== ax0.join()) cmds.push(`AXES,${ax.join(",")}`);
        const rg = [xmin, xmax, ymin, ymax].map((x) => x.value.trim());
        const rg0 = [d2.xmin, d2.xmax, d2.ymin, d2.ymax].map(sval);
        if (rg.map((v) => (v === "" ? "" : String(Number(v)))).join() !== rg0.map((v) => (v === "" ? "" : String(Number(v)))).join()) cmds.push(`PLOTRANGE,${rg.join(",")}`);
        return cmdOk(cmds, a);
      }},
      {label: "Cancel", action: () => true}],
      onClose: () => { if (D._lineDialogRefresh === refresh) D._lineDialogRefresh = null; }});
    D._lineDialog = {api, refresh, checked, src, fields: {title, xl, yl, xmin, xmax, ymin, ymax, xlog, ylog, xtick, ytick, start, nlines, file, replace}};
    await refresh();
    if (preset.file) inspect();
  };
  D.onLinesChanged = function () { if (D._lineDialogRefresh) D._lineDialogRefresh(); };
  D.quickPlot = function (path, plot) {
    D.lineSelection(plot === "th" ? "THPLOT" : "SPECPLOT", {file: path});
  };

  /** Select Dynamic Soil Property (spec 06 section 7.2): pick / create a DYNP property, edit its table.
   *  opt.select(label): called instead of plotting (SOIL tab "..." button). */
  D.soilProperty = async function (opt) {
    opt = opt || {};
    let props = {};
    try { props = (await S.get("/api/dynp")).properties; } catch (e) { /* none */ }
    let cur = opt.initial && props[opt.initial] ? opt.initial : Object.keys(props)[0] || null;
    const lst = el("div", {class: "checklist", style: {height: "120px"}});
    const tbody = el("tbody");
    const titleIn = txt(cur || "");
    const edits = {};          // label -> rows
    const rowsOf = (label) => edits[label] || (props[label] || []).map((r) => [r.sg, r.g, r.sd, r.d]);
    const drawList = () => {
      lst.innerHTML = "";
      const names = [...new Set(Object.keys(props).concat(Object.keys(edits)))].sort();
      if (!names.length) lst.appendChild(el("div", {class: "empty", text: "no dynamic soil property: New"}));
      for (const n of names) {
        const it = el("label", {class: n === cur ? "sel" : "", text: n});
        it.addEventListener("click", () => { cur = n; titleIn.value = n; drawList(); drawTable(); });
        lst.appendChild(it);
      }
    };
    const drawTable = () => {
      tbody.innerHTML = "";
      if (!cur) return;
      const rows = rowsOf(cur);
      for (let i = 0; i < 11; i++) {
        const r = rows[i] || ["", "", "", ""];
        const cells = r.map((v, j) => {
          const inp = txt(v === null || v === undefined ? "" : v, {style: {width: "90px"}});
          inp.addEventListener("input", () => {
            const rr = edits[cur] || (edits[cur] = rowsOf(cur).map((x) => x.slice()));
            while (rr.length <= i) rr.push(["", "", "", ""]);
            rr[i][j] = inp.value;
          });
          return el("td", {}, inp);
        });
        tbody.appendChild(el("tr", {}, el("td", {text: i + 1}), ...cells));
      }
    };
    const body = el("div", {style: {width: "560px"}},
      el("div", {style: {display: "flex", gap: "10px"}},
        el("div", {style: {flex: "1"}}, lst),
        el("div", {style: {display: "flex", flexDirection: "column", gap: "6px"}},
          el("button", {class: "btn", text: "New", onclick: async () => {
            const n = await D.ask("New Dynamic Soil Property", "Name (case-sensitive)");
            if (!n) return;
            edits[n] = edits[n] || [];
            cur = n; titleIn.value = n; drawList(); drawTable();
          }}),
          el("button", {class: "btn", text: "Edit", disabled: true, title: "renaming has no command in ACS SASSI: create a New property"}),
          el("button", {class: "btn", text: "Delete", disabled: true, title: "no command deletes a DYNP property"}))),
      row("Title", titleIn),
      el("div", {class: "tablewrap", style: {maxHeight: "300px"}}, el("table", {class: "grid"},
        el("thead", {}, el("tr", {}, el("th", {text: "#"}), el("th", {text: "Strain"}), el("th", {text: "Mod. Red."}), el("th", {text: "Strain"}), el("th", {text: "Damp"}))), tbody)),
      el("div", {class: "opt-note", text: "Strain in %, G/Gmax in [0, 1], damping in % (SHAKE convention); at most 11 points. OK submits DYNP,<no>,<sg>,<g>,<sd>,<d>,<label> for the edited points."}));
    const api = D.modal({title: "Select Dynamic Soil Property", body, buttons: [
      {label: "Ok", primary: true, action: async (a) => {
        const cmds = [];
        for (const [label, rows] of Object.entries(edits)) {
          rows.forEach((r, i) => {
            if (r.every((v) => String(v).trim() === "")) return;
            cmds.push(`DYNP,${i + 1},${r.map((v) => String(v).trim()).join(",")},${label}`);
          });
        }
        if (cmds.length && !(await cmdOk(cmds, a))) return false;
        if (!cur) { a.setMessage("select a property"); return false; }
        if (opt.select) { opt.select(cur); return true; }
        return cmdOk(`SOILPROPPLOT,${q(cur)}`, a);
      }},
      {label: "Cancel", action: () => true}]});
    drawList();
    drawTable();
    return api;
  };

  D.procFrame = function () {
    const lf = txt(""), dir = txt(""), desc = txt("");
    const types = [[0, "Bubble"], [1, "Vector"], [2, "Contour"], [3, "Time History"]];
    const radios = el("div", {class: "radios"}, ...types.map(([v, l]) => el("label", {}, el("input", {type: "radio", name: "anitype", value: v, checked: v === 0}), " " + l)),
      ...["Stress DB (Binary)", "ACC DB (Binary)", "RelDisp DB (Binary)"].map((l) => el("label", {style: {color: "#9aa3ae"}}, el("input", {type: "radio", disabled: true}), " " + l + " (P2)")));
    D.modal({title: "Parse Frame Data", body: el("div", {class: "dlg", style: {width: "600px"}},
      row("List File Name", lf, el("button", {class: "btn small", text: "<<", title: "a frame list file (*.dispani, *.tfiani ...)", onclick: () => D.pickFile({title: "Animation frame list", filter: [".dispani", ".tfiani", ".impani", ".txt", ".lst"], onOk: (p) => { lf.value = p; }})}),
        el("button", {class: "btn small", text: "Folder", title: "a frame folder written by a module restart option (TFU/, THD/, ACC/, RS/ ...)", onclick: () => D.pickFile({title: "Frame folder (MOTION / RELDISP / STRESS restart frames)", folder: true, onOk: (p) => { lf.value = p; }})})),
      row("Frame Storage Dir", dir, el("button", {class: "btn small", text: "<<", onclick: () => D.pickFile({title: "Frame storage directory", folder: true, onOk: (p) => { dir.value = p; }})})),
      row("Data Description", desc), fieldset("Plot Type", radios),
      el("div", {class: "opt-note", text: "List File Name: a frame list (line 1 ignored, then one frame file per line) or, as a SASSI-EDU extension, " +
        "a frame folder written by the module restart options (MOTION Restart for TF / ACC / RS: TFU/, ACC/, RS/; RELDISP Restart For Frame " +
        "Generation: THD/; STRESS Restart for Nodal Stress / Soil Pressure Contours). The Frame Storage Dir is created by PROCFRAME."})),
      buttons: [{label: "Ok", primary: true, action: (a) => {
        const t = (radios.querySelector("input[name=anitype]:checked") || {}).value || 0;
        if (!lf.value || !dir.value) { a.setMessage("List File Name and Frame Storage Dir are required"); return false; }
        return cmdOk(`PROCFRAME,${q(lf.value)},${q(dir.value)},${q(desc.value)},${t}`, a);
      }}, {label: "Cancel", action: () => true}]});
  };

  D.loadFrameData = async function (kind) {
    let db = {entries: []};
    try { db = await S.get("/api/animations"); } catch (e) { /* none */ }
    let sel = null;
    const tbody = el("tbody");
    const start = txt("1", {class: "num"}), end = txt("", {class: "num"}), stride = txt("1", {class: "num"});
    const scale = txt("1", {class: "num"}), cmin = txt("", {class: "num"}), cmax = txt("", {class: "num"}), col = txt("1", {class: "num"});
    const colour = kind === "BUBBLEPLOT" || kind === "CONTOURPLOT";
    const draw = () => {
      tbody.innerHTML = "";
      if (!db.entries.length) tbody.appendChild(el("tr", {}, el("td", {colspan: 4, class: "l", text: "no processed animation: Plot > Process Animation Frame List"})));
      for (const e of db.entries) {
        const tr = el("tr", {class: sel === e ? "sel" : ""}, el("td", {class: "l", text: e.description || ""}), el("td", {class: "l", text: e.directory}),
          el("td", {class: "l", text: ({0: "Bubble", 1: "Vector", 2: "Contour", 3: "Time History"})[e.type] || e.type}), el("td", {text: e.frames}));
        tr.addEventListener("click", () => {
          sel = e;
          start.value = e.start || 1; end.value = e.end || e.frames; stride.value = e.stride || 1;
          scale.value = e.scale || 1; cmin.value = e.cmin || ""; cmax.value = e.cmax || "";
          draw();
        });
        tbody.appendChild(tr);
      }
    };
    const body = el("div", {class: "dlg", style: {width: "700px"}},
      fieldset("Select From Database", el("div", {class: "tablewrap", style: {maxHeight: "200px"}}, el("table", {class: "grid", style: {width: "100%"}},
        el("thead", {}, el("tr", {}, el("th", {text: "Description"}), el("th", {text: "Animation Directory"}), el("th", {text: "Type"}), el("th", {text: "Frames"}))), tbody)),
        el("div", {class: "row"}, el("button", {class: "btn small", text: "Remove Animation", onclick: async () => {
          if (!sel) return;
          if (!(await D.confirm("Remove Animation", `Remove ${sel.directory} from the animation database?`))) return;
          const del = await D.confirm("Remove Animation", "Also delete its processed frame files (frame_*.npy, nodes*.npy, index.json)?", "Delete files", "Keep files");
          try { db = await S.post("/api/animations/remove", {directory: sel.directory, delete_files: del}); sel = null; draw(); } catch (e) { api.setMessage(e.message); }
        }}), el("span", {class: "opt-note", text: db.path || ""}))),
      fieldset("Animation Control", row("Start", start, el("label", {text: "End"}), end, el("label", {text: "Stride"}), stride),
        colour ? row("Colormap Min", cmin, el("label", {text: "Max"}), cmax, el("label", {text: "Col"}), col) : row("Scale Factor", scale)));
    const api = D.modal({title: "Load Frame Data", body, buttons: [
      {label: "Ok", primary: true, action: (a) => {
        if (!sel) { a.setMessage("click an animation of the list"); return false; }
        const line = colour ? `${kind},${q(sel.directory)},${start.value},${end.value},${stride.value},${cmin.value},${cmax.value},${col.value}`
          : `${kind},${q(sel.directory)},${start.value},${end.value},${stride.value},${scale.value}`;
        return cmdOk(line, a);
      }},
      {label: "Cancel", action: () => true}]});
    draw();
  };

  // ================================================================== settings windows of the active plot
  D.windowSettings = function () {
    const p = S.P.activePlot();
    if (!p) { S.local("WARNING", "Windows Settings: no active plot"); return; }
    if (p.kind === "SPECPLOT" || p.kind === "THPLOT") return D.graphPlotOptions(p);
    if (p.kind === "LAYERPLOT") return D.layerSettings(p);
    if (p.kind === "SOILPROPPLOT") return D.soilPropSettings(p);
    return D.windowOptions3D(p);
  };
  function axesBox(s, prefix) {
    return {
      min: txt(s[prefix + "min"] === null || s[prefix + "min"] === undefined ? "" : s[prefix + "min"], {class: "num"}),
      max: txt(s[prefix + "max"] === null || s[prefix + "max"] === undefined ? "" : s[prefix + "max"], {class: "num"}),
      log: chk(s["log_" + prefix]), tick: chk(s["minor_" + prefix]),
    };
  }
  D.graphPlotOptions = async function (p) {
    let lines = {};
    try { lines = (await S.get("/api/lines")).lines; } catch (e) { /* none */ }
    const s = p.settings;
    const shown = new Set(p.params.lines || []);
    const chosen = new Set(shown);
    let highlighted = null;
    const lst = el("div", {class: "checklist"});
    const markers = chk(false, {disabled: true});
    const drawList = () => {
      lst.innerHTML = "";
      for (const n of Object.keys(lines).map(Number).sort((a, b) => a - b)) {
        const c = chk(chosen.has(n));
        c.addEventListener("change", () => { if (c.checked) chosen.add(n); else chosen.delete(n); });
        const lab = el("label", {class: n === highlighted ? "sel" : ""}, c, `${n}: ${lines[n].name}`);
        lab.addEventListener("click", (ev) => { if (ev.target === c) return; highlighted = n; markers.disabled = false; markers.checked = !!lines[n].markers; drawList(); });
        lst.appendChild(lab);
      }
    };
    const title = txt(p.title || ""), xl = txt(s.xtitle || ""), yl = txt(s.ytitle || "");
    const X = axesBox(s, "x"), Y = axesBox(s, "y");
    const stipple = chk(s.stipple);
    const dest = txt("", {class: "num", placeholder: "next free"});
    const pk = txt("0", {class: "num"}), br = txt("0", {class: "num"});
    const c1 = txt("1", {class: "num"}), c2 = txt("1", {class: "num"}), c3 = txt("1", {class: "num"});
    const destNum = () => dest.value.trim() || String(Math.max(0, ...Object.keys(lines).map(Number)) + 1);
    const srcs = () => [...chosen].sort((a, b) => a - b);
    const calc = async (mk) => {
      const sl = srcs();
      if (!sl.length) { api.setMessage("check the source lines"); return; }
      const d = destNum();                     // fixed before the new line exists
      const r = await S.command(mk(d, sl));
      if (r.ok) { lines = (await S.get("/api/lines")).lines; chosen.add(Number(d)); drawList(); api.setMessage(`line ${d} created (no file written: use WRITESPEC / WRITETH)`, true); }
      else api.setMessage(r.messages.filter((m) => m.kind === "ERROR").map((m) => m.text).join("\n"));
    };
    const body = el("div", {class: "dlg", style: {width: "720px"}},
      el("div", {class: "dlg-cols"},
        el("div", {}, fieldset("Lines on the plot / inputs of the calculations", lst),
          fieldset("Line Options", el("label", {class: "chk"}, markers, " Data Points (highlighted line)"), el("label", {class: "chk"}, stipple, " Line Stippling (all lines)"))),
        el("div", {},
          fieldset("Title / Axis Labels", row("Title", title), row("X-Label", xl), row("Y-Label", yl)),
          fieldset("X Axis Options", row("Min", X.min, el("label", {text: "Max"}), X.max), el("div", {class: "row"}, el("label", {class: "chk"}, X.log, " Logarithmic"), el("label", {class: "chk"}, X.tick, " Show Ticks"))),
          fieldset("Y Axis Options", row("Min", Y.min, el("label", {text: "Max"}), Y.max), el("div", {class: "row"}, el("label", {class: "chk"}, Y.log, " Logarithmic"), el("label", {class: "chk"}, Y.tick, " Show Ticks"))))),
      fieldset("Spectra Analysis Postprocessing",
        el("div", {class: "row"},
          el("button", {class: "btn small", text: "Average", onclick: () => calc((d, sl) => `AVERAGE,${d},${sl.join(",")}`)}),
          el("button", {class: "btn small", text: "SRSS", onclick: () => calc((d, sl) => `SRSS,${d},${sl.join(",")}`)}),
          el("label", {text: "Peak Difference(%)"}), pk, el("label", {text: "Broaden(%)"}), br,
          el("button", {class: "btn small", text: "Broaden", onclick: () => calc((d, sl) => `BROADEN,${d},${pk.value},${br.value},${sl.join(",")}`)})),
        el("div", {class: "row"}, el("label", {text: "Spectra 1"}), c1, el("label", {text: "2"}), c2, el("label", {text: "3"}), c3,
          el("button", {class: "btn small", text: "Linear Combin", onclick: () => calc((d, sl) => `LINECOMBIN,${d},` + sl.slice(0, 3).map((n, i) => `${n},${[c1, c2, c3][i].value}`).join(","))}))),
      fieldset("Temporal Analysis Post Processing", el("div", {class: "row"},
        el("button", {class: "btn small", text: "Addition", onclick: () => calc((d, sl) => `ADDITION,${d},${sl.join(",")}`)}),
        el("button", {class: "btn small", text: "Subtraction", onclick: () => calc((d, sl) => `SUBTRACTION,${d},${sl.join(",")}`)}),
        el("label", {text: "Post Processing Results: line number"}), dest)));
    const api = D.modal({title: "Graph Plot Options", body, buttons: [
      {label: "Ok", primary: true, action: async (a) => {
        const cmds = [];
        const sl = srcs();
        const same = sl.length === shown.size && sl.every((n) => shown.has(n));
        if (!same) {            // no command changes the lines of a plot: re-plot (L17)
          if (!sl.length) { a.setMessage("check at least one line"); return false; }
          cmds.push("CLOSEPLOT", `${p.kind},${sl.join(",")}`);
        }
        if (title.value !== (p.title || "") || (!same && title.value)) cmds.push(`PLOTTITLE,${title.value}`);
        if (xl.value !== (s.xtitle || "") || (!same && xl.value)) cmds.push(`XTITLE,${xl.value}`);
        if (yl.value !== (s.ytitle || "") || (!same && yl.value)) cmds.push(`YTITLE,${yl.value}`);
        const ax = [s.major_x ? 1 : 0, s.major_y ? 1 : 0, +X.tick.checked, +Y.tick.checked, +X.log.checked, +Y.log.checked];
        const ax0 = [s.major_x ? 1 : 0, s.major_y ? 1 : 0, +s.minor_x, +s.minor_y, +s.log_x, +s.log_y];
        if (!same || ax.join() !== ax0.join()) cmds.push(`AXES,${ax.join(",")}`);
        const rg = [X.min.value, X.max.value, Y.min.value, Y.max.value].map((v) => v.trim());
        const rg0 = [s.xmin, s.xmax, s.ymin, s.ymax].map((v) => v === null || v === undefined ? "" : String(v));
        if (rg.join() !== rg0.join() || (!same && rg.some(Boolean))) cmds.push(`PLOTRANGE,${rg.join(",")}`);
        if (stipple.checked !== !!s.stipple || (!same && stipple.checked)) cmds.push(`STIPPLE,${+stipple.checked}`);
        if (highlighted !== null && !markers.disabled && markers.checked !== !!lines[highlighted].markers) cmds.push(`MARKERS,${+markers.checked},${highlighted}`);
        if (!cmds.length) return true;
        return cmdOk(cmds, a);
      }},
      {label: "Cancel", action: () => true}]});
    drawList();
  };
  D.layerSettings = function (p) {
    const pr = p.params, show = pr.show || {};
    const start = txt(pr.start || 1, {class: "num"}), end = txt(pr.end === undefined ? -1 : pr.end, {class: "num"});
    const cols = [["THICK", "Show Thickness", "thick"], ["WEIGHT", "Show Specific Weight", "weight"], ["VP", "Show P-Wave Velocity", "vp"],
      ["VS", "Show S-Wave Velocity", "vs"], ["PDAMP", "Show P-Wave Damping Ratio", "pdamp"], ["SDAMP", "Show S-Wave Damping Ratio", "sdamp"]];
    const boxes = cols.map(([k, l, key]) => [k, key, chk(show[key] !== false && (key !== "thick" || show.thick))]);
    D.modal({title: "Soil Layer Windows Setting", body: el("div", {style: {width: "360px"}}, row("Start Layer", start), row("EndLayer (-1 deepest)", end),
      ...boxes.map(([, , c], i) => el("div", {class: "row"}, el("label", {class: "chk"}, c, " " + cols[i][1])))),
      buttons: [{label: "Ok", primary: true, action: (a) => {
        const cmds = [];
        if (String(start.value) !== String(pr.start || 1)) cmds.push(`WINDOWSETTINGS,START,${start.value}`);
        if (String(end.value) !== String(pr.end === undefined ? -1 : pr.end)) cmds.push(`WINDOWSETTINGS,END,${end.value}`);
        for (const [k, key, c] of boxes) {
          const was = show[key] !== false && (key !== "thick" || !!show.thick);
          if (c.checked !== was) cmds.push(`WINDOWSETTINGS,SHOW,${k},${+c.checked}`);
        }
        return cmds.length ? cmdOk(cmds, a) : true;
      }}, {label: "Cancel", action: () => true}]});
  };
  D.soilPropSettings = function (p) {
    const s = p.settings, show = (p.params && p.params.show) || {};
    const glog = chk(s.log_y), gtick = chk(s.minor_y), slog = chk(s.log_x), stick = chk(s.minor_x);
    const mod = chk(show.modulus !== false), damp = chk(show.damping !== false);
    D.modal({title: "Soil Properties Window Settings", body: el("div", {style: {width: "340px"}},
      fieldset("G,D Axis", el("label", {class: "chk"}, glog, " Logarithmic"), el("label", {class: "chk"}, gtick, " Show Ticks")),
      fieldset("Shear Strain Axis", el("label", {class: "chk"}, slog, " Logarithmic"), el("label", {class: "chk"}, stick, " Show Ticks")),
      el("label", {class: "chk"}, mod, " Show Shear Modulus Line"), el("label", {class: "chk"}, damp, " Show Damping Line")),
      buttons: [{label: "Ok", primary: true, action: (a) => {
        const cmds = [];
        const ax = [s.major_x ? 1 : 0, s.major_y ? 1 : 0, +stick.checked, +gtick.checked, +slog.checked, +glog.checked];
        const ax0 = [s.major_x ? 1 : 0, s.major_y ? 1 : 0, +s.minor_x, +s.minor_y, +s.log_x, +s.log_y];
        if (ax.join() !== ax0.join()) cmds.push(`AXES,${ax.join(",")}`);
        if (mod.checked !== (show.modulus !== false)) cmds.push(`WINDOWSETTINGS,SHOW,MODULUS,${+mod.checked}`);
        if (damp.checked !== (show.damping !== false)) cmds.push(`WINDOWSETTINGS,SHOW,DAMPING,${+damp.checked}`);
        return cmds.length ? cmdOk(cmds, a) : true;
      }}, {label: "Cancel", action: () => true}]});
  };
  /** Window Options (spec 06 section 1.4) for 3D plots and animations. */
  D.windowOptions3D = async function (p) {
    let md = null;
    try { md = await S.get(`/api/model?number=${p.model}`); } catch (e) { S.local("ERROR", e.message); return; }
    const ui = md.ui_state || {};
    const bb = ui.display_volume || md.bbox;
    const vol = bb.map((v) => txt(S.fmt(v), {class: "num", style: {width: "80px", flex: "0 0 80px"}}));
    const hidden = new Set(ui.hide_groups || []);
    const groupBoxes = md.groups.map((g) => [g.id, chk(!hidden.has(g.id)), `Group ${g.id} ${g.type_name} - ${g.n}`]);
    const anim = ["BUBBLEPLOT", "VECTORPLOT", "CONTOURPLOT", "DEFORMPLOT"].includes(p.kind);
    const colourOn = p.kind === "BUBBLEPLOT" || p.kind === "CONTOURPLOT";
    const cmin = txt(p.params.vmin !== undefined ? p.params.vmin : "", {class: "num", disabled: !colourOn});
    const cmax = txt(p.params.vmax !== undefined ? p.params.vmax : "", {class: "num", disabled: !colourOn});
    const dirSel = el("select", {disabled: p.kind !== "VECTORPLOT"}, ...["X", "Y", "Z", "ALL"].map((d) => el("option", {value: d, text: d === "ALL" ? "All" : d, selected: p.view.direction === d})));
    const hideMode = el("select", {}, el("option", {value: "HIDE", text: "Hide"}), el("option", {value: "SHOW", text: "Show"}));
    const grp = txt("", {class: "num", disabled: p.kind === "NODEPLOT"}), ids = txt("", {placeholder: "3 7 12  or  100-250"});
    const scale = txt(p.params.scale !== undefined ? p.params.scale : 1, {class: "num", disabled: !(p.kind === "VECTORPLOT" || p.kind === "DEFORMPLOT")});
    const pause = txt(p.params.frame_pause || 33, {class: "num"});
    const undeformed = chk(p.view.show_undeformed, {disabled: p.kind !== "DEFORMPLOT"});
    const title = txt(p.title || "");
    const applyHide = async () => {
      const list = ids.value.trim();
      if (!list) { api.setMessage("enter element or node numbers"); return; }
      const toks = list.replace(/\s*-\s*/g, "-").split(/[\s,;]+/).filter(Boolean);
      if (toks.length > 1 && toks.some((x) => x.includes("-"))) { api.setMessage("one list or one range per press (spec 06 section 1.4)"); return; }
      const idsTok = toks.join(",");
      const node = p.kind === "NODEPLOT" || !grp.value.trim();
      const line = node ? `WINDOWSETTINGS,${hideMode.value}NODE,${idsTok}` : `WINDOWSETTINGS,${hideMode.value}ELEM,${grp.value.trim()},${idsTok}`;
      if (await cmdOk(line, api)) api.setMessage(`${line}: applied`, true);
    };
    const body = el("div", {class: "dlg", style: {width: "640px"}},
      fieldset("Model Display Volume", el("div", {style: {display: "flex", flexWrap: "wrap", gap: "4px 18px"}},
        ...["X", "Y", "Z"].map((ax, i) => nameInputs(el("div", {class: "row"}, el("label", {text: `${ax} Min`}), vol[2 * i], el("label", {text: "Max"}), vol[2 * i + 1]), `${ax} range`)))),
      el("div", {class: "dlg-cols"},
        fieldset("Show Element Group", el("div", {class: "checklist", style: {height: "130px"}}, ...groupBoxes.map(([, c, l]) => el("label", {}, c, l)))),
        el("div", {},
          fieldset("Colormap Value Range", row("Min", cmin, el("label", {text: "Max"}), cmax)),
          fieldset("Output Direction", dirSel),
          fieldset("Hide/Show " + (p.kind === "NODEPLOT" ? "Nodes" : "Elements"), row("Action", hideMode), row("Group", grp), row(p.kind === "NODEPLOT" ? "Node Numbers" : "Elem. Numbers", ids),
            el("button", {class: "btn small", text: `Hide/Show ${p.kind === "NODEPLOT" ? "Nodes" : "Elem."}`, onclick: applyHide}),
            el("div", {class: "note", text: "applied immediately (kept even with Cancel); blank Group = nodes"})))),
      fieldset("Animation Options", row("Scale Factor", scale, el("label", {text: "Frame Pause (ms)"}), pause)),
      el("div", {class: "row"}, el("label", {class: "chk"}, undeformed, " Show Undeformed Shape")),
      row("Title", title),
      el("div", {class: "row"}, el("button", {class: "btn small", text: "Show all (SHOWALL)", onclick: () => cmdOk("WINDOWSETTINGS,SHOWALL", api)})));
    const api = D.modal({title: "Window Options", body, buttons: [
      {label: "OK", primary: true, action: (a) => {
        const cmds = [];
        const v = vol.map((x) => x.value.trim());
        const v0 = bb.map((x) => S.fmt(x));
        if (v.join() !== v0.join()) cmds.push(`WINDOWSETTINGS,VOLUME,${v.join(",")}`);
        for (const [g, c] of groupBoxes) {
          if (c.checked && hidden.has(g)) cmds.push(`WINDOWSETTINGS,SHOWGROUP,${g}`);
          if (!c.checked && !hidden.has(g)) cmds.push(`WINDOWSETTINGS,HIDEGROUP,${g}`);
        }
        if (colourOn && (String(cmin.value) !== String(p.params.vmin) || String(cmax.value) !== String(p.params.vmax))) cmds.push(`WINDOWSETTINGS,RANGE,${cmin.value},${cmax.value}`);
        if (p.kind === "VECTORPLOT" && dirSel.value !== p.view.direction) cmds.push(`WINDOWSETTINGS,DIRECTION,${dirSel.value}`);
        if (!scale.disabled && String(scale.value) !== String(p.params.scale)) cmds.push(`WINDOWSETTINGS,SCALE,${scale.value}`);
        if (String(pause.value) !== String(p.params.frame_pause || 33)) cmds.push(`WINDOWSETTINGS,FRAMEPAUSE,${pause.value}`);
        if (!undeformed.disabled && undeformed.checked !== !!p.view.show_undeformed) cmds.push(`WINDOWSETTINGS,UNDEFORMED,${+undeformed.checked}`);
        if (title.value !== (p.title || "")) cmds.push(`PLOTTITLE,${title.value}`);
        return cmds.length ? cmdOk(cmds, a) : true;
      }},
      {label: "Cancel", action: () => true}]});
    void anim;
  };
  D.shaderOptions = function () {
    const sh = (S.state && S.state.plots && S.state.plots.shader) || {points: 10, linew: 0.02, shrink: 0.06, scale: 1};
    const f = {points: txt(sh.points, {class: "num"}), scale: txt(sh.scale, {class: "num"}), linew: txt(sh.linew, {class: "num"}), shrink: txt(sh.shrink, {class: "num"})};
    D.modal({title: "Shader Options", body: el("div", {style: {width: "420px"}}, row("Node/Bubble Node Size", f.points), row("Vector/Displacement Scale Factor", f.scale),
      row("Element Outline Thickness (% of element)", f.linew), row("Element Shrink (% of element)", f.shrink),
      el("div", {class: "opt-note", text: "fractions: 0.06 = 6 % (SHADEROPTIONS,[points],[linew],[shrink],[scale])"})),
      buttons: [{label: "OK", primary: true, action: async (a) => {
        const ok = await cmdOk(`SHADEROPTIONS,${f.points.value},${f.linew.value},${f.shrink.value},${f.scale.value}`, a);
        if (ok) S.refreshState();
        return ok;
      }}, {label: "Cancel", action: () => true}]});
  };
  D.changeView = function () {
    const p = S.P.activePlot();
    if (!p) return;
    const v = p.view;
    const f = ["rx", "ry", "rz", "px", "py", "zoom"].map((k) => txt(S.fmt(v[k]), {class: "num"}));
    D.modal({title: "Change View (CNGVIEW)", body: el("div", {style: {width: "420px"}},
      row("Rotation about X (deg)", f[0]), row("Rotation about Y (deg)", f[1]), row("Rotation about Z (deg)", f[2]),
      row("Pan X", f[3]), row("Pan Y", f[4]), row("Zoom (1 = fit)", f[5])),
      buttons: [{label: "OK", primary: true, action: (a) => cmdOk(`CNGVIEW,${f.map((x) => x.value).join(",")}`, a)}, {label: "Cancel", action: () => true}]});
  };
  D.changeCenter = function () {
    const p = S.P.activePlot();
    if (!p) return;
    const c = p.view.center || ["", "", ""];
    const f = c.map((x) => txt(x === "" ? "" : S.fmt(x), {class: "num"}));
    D.modal({title: "Change Center (CNGCENTER)", body: el("div", {style: {width: "360px"}}, row("X", f[0]), row("Y", f[1]), row("Z", f[2])),
      buttons: [{label: "OK", primary: true, action: (a) => cmdOk(`CNGCENTER,${f.map((x) => x.value).join(",")}`, a)}, {label: "Cancel", action: () => true}]});
  };
  D.showDof = function () {
    const p = S.P.activePlot();
    const cur = new Set((p && p.view.show_dof) || []);
    const labs = ["X", "Y", "Z", "XX", "YY", "ZZ"];
    const boxes = labs.map((l, i) => [l, chk(cur.has(i))]);
    D.modal({title: "Boundary Conditions", body: el("div", {style: {width: "320px"}},
      el("div", {style: {display: "grid", gridTemplateColumns: "repeat(3, auto)", gap: "4px 16px", justifyContent: "start"}},
        ...boxes.map(([l, c]) => el("label", {class: "chk"}, c, " " + l))),
      el("div", {class: "opt-note", text: "nodes with a fixed DOF among the checked ones are marked green (SHOWDOF)"})),
      buttons: [{label: "OK", primary: true, action: (a) => {
        const on = boxes.filter(([, c]) => c.checked).map(([l]) => l);
        return cmdOk(`SHOWDOF,${on.length ? on.join(",") : "NONE"}`, a);
      }}, {label: "Cancel", action: () => true}]});
  };
  D.colours = function () {
    const c = (S.state.settings && S.state.settings.colours) || {};
    const kinds = [["ECHO", "Command echo"], ["CONFIRM", "Output confirmation"], ["COMMENT", "Comments"], ["INFO", "Information"], ["WARNING", "Warnings"], ["ERROR", "Errors"]];
    const inputs = kinds.map(([k]) => el("input", {type: "color", value: c[k] || "#000000"}));
    D.modal({title: "Select Colors", body: el("div", {style: {width: "360px"}}, ...kinds.map(([, l], i) => row(l, inputs[i])),
      el("div", {class: "opt-note", text: "Command History colours (saved in SASSIini.xml). Plot palettes: COLOR,<Palette>,<Num>,<R>,<G>,<B>."})),
      buttons: [{label: "OK", primary: true, action: async () => {
        const body = {colours: {}};
        kinds.forEach(([k], i) => { body.colours[k] = inputs[i].value; });
        const r = await S.post("/api/settings/colours", body);
        S.state.settings = r.settings;
        S.applyColours();
        return true;
      }}, {label: "Cancel", action: () => true}]});
  };

  // ================================================================== Modules menu
  D.moduleLocation = function () {
    const st = S.state.settings;
    const inputs = {};
    const rows = st.location_rows.map((m) => {
      const dis = st.location_disabled.includes(m);
      const inp = txt(dis ? "" : st.locations[m] || "built-in", {disabled: dis, style: {flex: "1"}});
      inputs[m] = inp;
      return row(`${m} Module`, inp, el("button", {class: "btn small", text: "built-in", disabled: dis, onclick: () => { inp.value = "built-in"; }}));
    });
    D.modal({title: "Module Directories", body: el("div", {style: {width: "620px"}}, ...rows,
      el("div", {class: "opt-note", text: "built-in = the Python module of SASSI-EDU (default). An external executable is run in the model directory with the three-line batch protocol (model, deck, listing) on its standard input."})),
      buttons: [{label: "Ok", primary: true, action: async (a) => {
        const loc = {};
        for (const [m, inp] of Object.entries(inputs)) if (!inp.disabled) loc[m] = inp.value.trim() || "built-in";
        const r = await S.post("/api/settings/locations", {locations: loc});
        S.state.settings = r.settings;
        if (r.errors.length) { a.setMessage(r.errors.join("\n")); return false; }
        return true;
      }}, {label: "Cancel", action: () => true}]});
  };
  D.moduleExtension = function () {
    const ext = JSON.parse(JSON.stringify(S.state.settings.extensions));
    const sel = el("select", {}, ...Object.keys(ext).map((m) => el("option", {value: m, text: m})));
    const ie = txt(""), oe = txt("");
    const load = () => { ie.value = ext[sel.value][0]; oe.value = ext[sel.value][1]; };
    sel.addEventListener("change", load);
    ie.addEventListener("input", () => { ext[sel.value][0] = ie.value; });
    oe.addEventListener("input", () => { ext[sel.value][1] = oe.value; });
    D.modal({title: "File Extension Options", body: el("div", {style: {width: "420px"}}, row("Module", sel), row("Input File Extension", ie), row("Output File Extension", oe),
      el("div", {class: "opt-note", text: "Ok saves the edits of all modules (SASSIini.xml); Cancel drops them. SASSI-EDU decks use the manual's extensions; these values are passed to external executables."})),
      buttons: [{label: "Ok", primary: true, action: async () => {
        const r = await S.post("/api/settings/extensions", {extensions: ext});
        S.state.settings = r.settings;
        return true;
      }}, {label: "Cancel", action: () => true}]});
    load();
  };

  // ================================================================== View / Help
  /** View > Check Errors (requirements 5.5): the .err text of the last CHECK.  It also pops up by
   *  itself after a CHECK / AFWRITE with messages (event 'check') unless Options > Check > Suppress
   *  Error Window is set; an open window is refreshed instead of opening a second one. */
  D.checkErrors = async function (opts) {
    let d;
    try { d = await S.get("/api/check"); } catch (e) { S.local("ERROR", e.message); return; }
    if (D._checkWin && D._checkWin.open) { D._checkWin.fill(d); return D._checkWin.api; }
    const pre = el("pre", {class: "textview", style: {maxHeight: "60vh", overflow: "auto", border: "1px solid var(--border-2)", width: "720px", maxWidth: "100%"}});
    const note = el("div", {class: "opt-note", text: ""});
    let api = null;
    const fill = (x) => {
      pre.innerHTML = "";
      const lines = (x.text || "(no CHECK yet: press Run CHECK, or run AFWRITE)").split("\n");
      for (const ln of lines) {
        const cls = /^Error\b/.test(ln) ? "m-ERROR" : /^Warning\b/.test(ln) ? "m-WARNING" : /^Errors and Warnings for/.test(ln) ? "m-ECHO" : "";
        pre.appendChild(el("div", {class: cls, text: ln}));
      }
      note.textContent = x.source ? `from ${x.source}` : "";
      if (api) api.setMessage(x.summary || "", true);
    };
    const box = el("div", {class: "history"}, pre);
    const win = {open: true, fill};
    api = D.modal({title: d.title, body: el("div", {}, box, note), onClose: () => { win.open = false; if (D._checkWin === win) D._checkWin = null; },
      buttons: [
        {label: "Run CHECK", action: async () => { await S.command("CHECK"); fill(await S.get("/api/check")); return false; }},
        {label: "Close", primary: true, action: () => true}]});
    win.api = api;
    D._checkWin = win;
    fill(d);
    return api;
  };
  D.about = async function () {
    let a;
    try { a = await S.get("/api/about"); } catch (e) { S.local("ERROR", e.message); return; }
    const kv = el("div", {class: "kv"});
    for (const [k, l] of [["version", "Version"], ["build", "Build"], ["methodology", "Methodology"], ["python", "Python"], ["numpy", "NumPy"],
      ["scipy", "SciPy"], ["plotly", "Plotly"], ["matplotlib", "Matplotlib"], ["platform", "Platform"], ["settings", "Settings"]]) kv.append(el("span", {class: "k", text: l}), el("span", {text: a[k]}));
    D.modal({title: "About SASSI-EDU", body: el("div", {style: {width: "560px"}}, el("h3", {text: a.title, style: {margin: "0 0 10px"}}),
      el("p", {text: "An educational Python re-implementation of the ACS SASSI soil-structure-interaction program (module and command nomenclature of the ACS SASSI V3 manual). Not affiliated with or endorsed by the authors of ACS SASSI."}),
      el("p", {text: "A learning tool: not qualified for design or licensing work."}),
      a.web ? el("p", {text: "Browser version: Python runs in this tab (Pyodide); nothing is sent to a server. " + WEB_NOTE}) : null, kv)});
  };
  /** Help > Help (F1, requirements 5.3): the project documentation rendered to HTML by the server
   *  (sassi.ui.markdown: raw HTML never passes) -- documentation index, User Guide, GUI guide, ANSYS,
   *  Theory Manual, Verification Manual, Command Reference, tutorial examples -- with the headings of
   *  the open document, links between documents, Back / Forward, a search over all documents and
   *  the live command index of the interpreter.  opts: {doc, anchor} opens that page. */
  D.helpTab = async function (opts) {
    opts = opts || {};
    const old = S.tab("help");
    if (old) {
      S.selectTab("help");
      if (opts.doc) old.help.open(opts.doc, opts.anchor);
      return old.help;
    }
    let h;
    try { h = await S.get("/api/help"); } catch (e) { S.local("ERROR", e.message); return null; }
    const docs = h.docs || [];
    const cache = {};
    if (h.home_doc) cache[h.home_doc.id] = h.home_doc;            // rendered with the document list
    let cur = null;                // {id, title, toc} of the page shown ("@commands", "@search" for the built-in pages)
    const hist = [];
    let pos = -1;
    const content = el("article", {class: "doc"});
    const main = el("div", {class: "help-main"}, content);
    const docList = el("div", {class: "help-docs"});
    const tocBox = el("div", {class: "help-toc"});
    const tocFilter = txt("", {placeholder: "filter the headings", "aria-label": "filter the headings", class: "help-filter"});
    const search = txt("", {placeholder: "search all documents", "aria-label": "search all documents", style: {width: "220px", flex: "0 1 220px"}});
    const where = el("span", {class: "grow help-where"});
    const back = el("button", {class: "btn small", text: "◀ Back", title: "previous page", disabled: true});
    const fwd = el("button", {class: "btn small", text: "Forward ▶", title: "next page", disabled: true});
    const nav = el("nav", {class: "help-nav", "aria-label": "Help documents"},
      el("div", {class: "help-head", text: "Documents"}), docList,
      el("div", {class: "help-head", text: "Contents"}), tocFilter, tocBox);
    const groups = [];
    for (const d of docs) if (!groups.includes(d.group)) groups.push(d.group);
    const drawDocs = () => {
      docList.innerHTML = "";
      for (const g of groups) {
        docList.appendChild(el("div", {class: "help-group", text: g}));
        for (const d of docs.filter((x) => x.group === g)) {
          const it = el("a", {href: "#", class: "help-doc" + (cur && cur.id === d.id ? " sel" : ""), title: d.id, text: d.title});
          it.addEventListener("click", (ev) => { ev.preventDefault(); open(d.id); });
          docList.appendChild(it);
        }
      }
      docList.appendChild(el("div", {class: "help-group", text: "Commands"}));
      const ci = el("a", {href: "#", class: "help-doc" + (cur && cur.id === "@commands" ? " sel" : ""), text: "Command index (live)", title: "every command of the interpreter with abbreviation, tier and availability"});
      ci.addEventListener("click", (ev) => { ev.preventDefault(); open("@commands"); });
      docList.appendChild(ci);
      if (!docs.length) docList.appendChild(el("div", {class: "empty", text: `no documentation found under ${h.root || "the project"}/docs`}));
    };
    const drawToc = () => {
      tocBox.innerHTML = "";
      const f = tocFilter.value.trim().toLowerCase();
      for (const t of (cur && cur.toc) || []) {
        if (t.level > 3 || (f && !t.text.toLowerCase().includes(f))) continue;
        const it = el("a", {href: "#", class: "help-toc-item lvl" + t.level, text: t.text, title: t.text});
        it.addEventListener("click", (ev) => { ev.preventDefault(); go(t.slug); });
        tocBox.appendChild(it);
      }
    };
    tocFilter.addEventListener("input", drawToc);
    const go = (slug) => {
      if (!slug) { main.scrollTop = 0; return true; }
      const target = content.querySelector(`[data-slug="${CSS.escape(slug)}"]`);
      if (!target) return false;
      main.scrollTop += target.getBoundingClientRect().top - main.getBoundingClientRect().top - 6;
      target.classList.remove("flash");
      void target.offsetWidth;
      target.classList.add("flash");
      return true;
    };
    const commandsPage = () => {
      const filter = txt("", {placeholder: "filter commands, e.g. PLOT or RUN", "aria-label": "filter commands", style: {width: "260px"}});
      const tbody = el("tbody");
      const draw = () => {
        const f = filter.value.trim().toUpperCase();
        tbody.innerHTML = "";
        for (const c of h.commands) {
          if (f && !c.name.includes(f) && !(c.summary || "").toUpperCase().includes(f) && !c.abbrev.join(" ").includes(f)) continue;
          tbody.appendChild(el("tr", {}, el("td", {class: "l", style: {fontFamily: "var(--mono)"}, text: c.name}), el("td", {class: "l", text: c.abbrev.join(" ")}),
            el("td", {text: c.tier}), el("td", {class: "l", text: c.available ? "" : "placeholder"}), el("td", {class: "l", text: c.summary})));
        }
      };
      filter.addEventListener("input", draw);
      draw();
      return [el("h1", {text: "Command index"}),
        el("p", {text: "Every command of the interpreter with its abbreviations, tier and the first line of its description (from the running program; the Command Reference document has the full syntax)."}),
        el("div", {class: "row"}, filter),
        el("div", {class: "md-table"}, el("table", {class: "grid help-list"}, el("thead", {}, el("tr", {}, el("th", {text: "Command"}), el("th", {text: "Abbrev."}),
          el("th", {text: "Tier"}), el("th", {text: ""}), el("th", {class: "l", text: "Summary"}))), tbody))];
    };
    const searchPage = async (qtext) => {
      let r;
      try { r = await S.get(`/api/help/search?q=${encodeURIComponent(qtext)}`); } catch (e) { return [el("p", {class: "md-error", text: e.message})]; }
      const out = [el("h1", {text: `Search: ${qtext}`}),
        el("p", {text: r.results.length ? `${r.results.length}${r.truncated ? "+" : ""} section(s) contain "${qtext}".` : `No section contains "${qtext}".`})];
      const ul = el("ul", {class: "help-results"});
      for (const x of r.results) {
        const a = el("a", {href: "#", text: `${x.title}${x.heading ? " › " + x.heading : ""}`});
        a.addEventListener("click", (ev) => { ev.preventDefault(); open(x.doc, x.slug); });
        ul.appendChild(el("li", {}, a, el("div", {class: "help-snippet", text: x.snippet})));
      }
      out.push(ul);
      return out;
    };
    /** Show page id (a document id, "@commands" or "@search:<text>") and scroll to anchor. */
    const open = async (id, anchor, opt2) => {
      opt2 = opt2 || {};
      if (!opt2.noHistory && cur && hist[pos]) hist[pos].scroll = main.scrollTop;
      let page;
      if (id === "@commands") page = {id, title: "Command index", toc: [], nodes: commandsPage()};
      else if (id.startsWith("@search:")) page = {id, title: "Search", toc: [], nodes: await searchPage(id.slice(8))};
      else {
        let d = cache[id];
        if (!d) {
          try { d = await S.get(`/api/help/doc?name=${encodeURIComponent(id)}`); } catch (e) {
            S.local("ERROR", `Help: ${e.message}`);
            return false;
          }
          cache[id] = d;
        }
        page = {id: d.id, title: d.title, toc: d.toc, html: d.html};
      }
      cur = page;
      content.innerHTML = "";
      if (page.html !== undefined) content.innerHTML = page.html;     // rendered by sassi.ui.markdown (escaped, no raw HTML)
      else content.append(...page.nodes);
      where.textContent = page.id.startsWith("@") ? page.title : `${page.title}  —  ${page.id}`;
      drawDocs();
      tocFilter.value = "";
      drawToc();
      if (!opt2.noHistory) {
        hist.splice(pos + 1);
        hist.push({id: page.id, anchor: anchor || null, scroll: 0});
        pos = hist.length - 1;
      }
      back.disabled = pos <= 0;
      fwd.disabled = pos >= hist.length - 1;
      if (opt2.scroll !== undefined) main.scrollTop = opt2.scroll;
      else if (!go(anchor)) main.scrollTop = 0;
      return true;
    };
    const step = (k) => {
      const t = hist[pos + k];
      if (!t) return;
      if (hist[pos]) hist[pos].scroll = main.scrollTop;
      pos += k;
      open(t.id, t.anchor, {noHistory: true, scroll: t.scroll});
    };
    back.addEventListener("click", () => step(-1));
    fwd.addEventListener("click", () => step(1));
    content.addEventListener("click", (ev) => {
      const a = ev.target.closest ? ev.target.closest("a") : null;
      if (!a || !content.contains(a)) return;
      const doc = a.getAttribute("data-doc"), anc = a.getAttribute("data-anchor");
      if (doc) { ev.preventDefault(); open(doc, anc); return; }
      if (anc !== null && (a.getAttribute("href") || "").startsWith("#")) {
        ev.preventDefault();
        if (hist[pos]) hist[pos].scroll = main.scrollTop;
        if (go(anc)) { hist.splice(pos + 1); hist.push({id: cur.id, anchor: anc, scroll: main.scrollTop}); pos = hist.length - 1; back.disabled = false; fwd.disabled = true; }
      }
    });
    search.addEventListener("keydown", (ev) => { if (ev.key === "Enter" && search.value.trim().length >= 2) open("@search:" + search.value.trim()); });
    const pane = el("div", {class: "pane"},
      el("div", {class: "pane-bar"}, back, fwd,
        el("button", {class: "btn small", text: "Home", title: "documentation index", onclick: () => open(h.home || "@commands")}),
        where, search,
        el("button", {class: "btn small", text: "Search", onclick: () => { if (search.value.trim().length >= 2) open("@search:" + search.value.trim()); }})),
      el("div", {class: "pane-body help-layout"}, nav, main));
    const t = S.addTab({id: "help", title: "Help", kind: "help", pane});
    t.help = {open, go, docs, current: () => cur, back: () => step(-1), forward: () => step(1)};
    S.selectTab("help");
    await open(opts.doc || h.home || "@commands", opts.anchor);
    return t.help;
  };
  D.verificationTab = async function () {
    if (S.tab("verify")) { S.selectTab("verify"); return; }
    const sel = el("select", {}, ...["P0", "P1", "P2", "ALL"].map((v) => el("option", {value: v, text: v === "ALL" ? "all problems" : `tier ${v}`})));
    const one = txt("", {placeholder: "or a VP id, e.g. VP-30", style: {width: "160px"}});
    const status = el("span", {class: "grow"});
    const tableHost = el("div", {style: {padding: "8px"}});
    const listHost = el("div", {style: {padding: "8px"}});
    const run = async () => {
      const s = one.value.trim() || sel.value;
      let job;
      try { job = await S.post("/api/verify", {select: s}); } catch (e) { status.textContent = e.message; return; }
      status.textContent = `VERIFY,${s}: running (job ${job.id}) ...`;
      const jt = S.openJobTab(job);
      S.selectTab("verify");
      S.jobs[job.id].onUpdate = (d) => {
        drawTable(d.table);
        status.textContent = `${d.line}: ${d.state}  ${d.table.summary || ""}`;
      };
      void jt;
    };
    const drawTable = (t) => {
      tableHost.innerHTML = "";
      if (!t || !t.problems.length) return;
      const tb = el("table", {class: "grid", style: {width: "100%"}});
      tb.appendChild(el("thead", {}, el("tr", {}, ...["Problem", "Quantity", "Computed", "Reference", "Error", "Tolerance", "Result"].map((h) => el("th", {text: h})))));
      const body = el("tbody");
      for (const p of t.problems) {
        const st = p.status === "PASSED" ? "pass" : p.status === "RUNNING" ? "" : "fail";
        body.appendChild(el("tr", {}, el("td", {class: "l", colspan: 6}, el("strong", {text: p.id}), " " + p.title + (p.elapsed !== null && p.elapsed !== undefined ? `  (${p.elapsed.toFixed(2)} s)` : "")),
          el("td", {class: st, text: p.status})));
        for (const c of p.checks) body.appendChild(el("tr", {}, el("td"), el("td", {class: "l", text: c.quantity}), el("td", {text: S.fmt(c.computed)}), el("td", {text: S.fmt(c.reference)}),
          el("td", {text: c.error === null ? "" : `${S.fmt(c.error, 3)} ${c.kind}`}), el("td", {text: S.fmt(c.tolerance, 3)}), el("td", {class: c.passed ? "pass" : "fail", text: c.passed ? "pass" : "FAIL"})));
        for (const n of p.notes) body.appendChild(el("tr", {}, el("td"), el("td", {class: "l", colspan: 6, style: {color: "var(--muted)"}, text: "note: " + n})));
      }
      tb.appendChild(body);
      tableHost.appendChild(tb);
    };
    const listBtn = el("button", {class: "btn small", text: "List problems", onclick: async () => {
      listHost.innerHTML = "";
      try {
        const r = await S.get("/api/verify/list");
        const tb = el("table", {class: "grid"}, el("thead", {}, el("tr", {}, el("th", {text: "Id"}), el("th", {text: "Tier"}), el("th", {class: "l", text: "Title"}), el("th", {class: "l", text: "Modules"}))));
        const body = el("tbody");
        for (const p of r.listed) {
          const tr = el("tr", {}, el("td", {class: "l", text: p.id}), el("td", {text: p.tier}), el("td", {class: "l", text: p.title}), el("td", {class: "l", text: p.modules}));
          tr.addEventListener("dblclick", () => { one.value = p.id; });
          body.appendChild(tr);
        }
        tb.appendChild(body);
        listHost.appendChild(el("div", {class: "opt-note", text: `${r.listed.length} verification problems (double-click one to select it)`}), tb);
      } catch (e) { listHost.textContent = e.message; }
    }});
    const pane = el("div", {class: "pane"}, el("div", {class: "pane-bar"}, el("strong", {text: "Verification (VERIFY)"}), sel, one,
      el("button", {class: "btn small primary", text: "Run", onclick: run}), listBtn, status),
      el("div", {class: "pane-body"}, el("div", {class: "opt-note", style: {padding: "8px 8px 0"}, text: "Computed, reference, error and tolerance of every check (requirements 6.1). Runs in a worker process; Cancel in the status bar stops it."}), tableHost, listHost));
    S.addTab({id: "verify", title: "Verification", kind: "verify", pane});
    S.selectTab("verify");
  };

  /** A command given without arguments asked for its dialog (sassi.plotting 'dialog' event). */
  D.openNamedDialog = function (name, ctx) {
    const cmd = (ctx && ctx.command) || "";
    switch (name) {
      case "Line Selection": return D.lineSelection(cmd || "SPECPLOT");
      case "Select Cut to Display": return D.cutPlot();
      case "Select Dynamic Soil Property": return D.soilProperty();
      case "Parse Frame Data": return D.procFrame();
      case "Load Frame Data": return D.loadFrameData(cmd || "BUBBLEPLOT");
      case "Boundary Conditions": return D.showDof();
      case "Shader Options": return D.shaderOptions();
      case "Window Settings": return D.windowSettings();
      default: S.local("WARNING", `dialog ${name} is not available`);
    }
  };

  // ================================================================== Options dialogs (requirements 5.4)
  D.optionsDialog = async function (name, startTab) {
    let data;
    try { data = await S.get(`/api/options/${name}`); } catch (e) { S.local("ERROR", e.message); return; }
    const form = data.form;
    const F = new OptionsForm(form, data.values, name);
    const body = F.build(startTab);
    /* Ok submits the changed values; Run (Modules > ANSYS Eq. Static Load / Dynamic Load: form.run) also
     * starts the module run -- the server executes the dialog's commands first, then RUNLOADGEN. */
    const commit = async (a, run) => {
      const payload = F.payload();
      if (!F.anyDirty() && !run) return true;
      if (run) payload.run = true;
      try {
        const r = await S.post(`/api/options/${name}`, payload);
        const bad = (r.results || []).filter((x) => !x.ok);
        if (bad.length) {
          a.setMessage(r.messages.filter((m) => m.kind === "ERROR").map((m) => m.text).join("\n") || "a command failed");
          return false;
        }
        if (r.run_error) { a.setMessage(r.run_error); return false; }
        if (r.job) S.openJobTab(r.job);
        if (name === "WRITE" || name === "CHECK") S.refreshState();
        return true;
      } catch (e) {
        const probs = (e.data && e.data.problems) || [e.message];
        a.setMessage(probs.join("\n"));
        return false;
      }
    };
    const buttons = [{label: "Ok", primary: true, action: (a) => commit(a, false)}];
    if (form.run) buttons.push({label: form.run.label || "Run", action: (a) => commit(a, true)});
    buttons.push({label: "Cancel", action: () => true});
    const api = D.modal({title: form.title, body, width: name === "ANALYSIS" ? "min(1180px, 96vw)" : (form.width || null), buttons});
    if (form.run) api.buttons[1].title = `Ok, then ${form.run.line} (module run with its output tab)`;
    F.modal = api;
    return api;
  };

  function OptionsForm(form, values, name) {
    this.form = form;
    this.v = JSON.parse(JSON.stringify(values));
    for (const k of ["xrecords", "xtables", "xindexed", "xindexed_defaults"]) this.v[k] = this.v[k] || {};
    this.orig = values;
    this.name = name;
    this.dirty = {records: {}, indexed: {}, strings: {}, lists: {}, requests: {}, xrecords: {}, xtables: {}, xindexed: {}};
    this.sel = {spec: 1, layer: 1, motion: 1, wave: (values.records && values.records.SITE && values.records.SITE.wopt === 1) ? 4 : 2};
    const curves = this.familyKeys("BBC");
    this.sel.bbc = curves.length ? curves[0] : 1;
    this.rows = [];
    this.activeTab = null;
  }
  /** Keys (numbers) of the entries of a keyed command-record family that exist in the dialog (%BBC). */
  OptionsForm.prototype.familyKeys = function (fam) {
    return Object.entries(this.v.xindexed[fam] || {}).filter(([, e]) => e !== null && e !== undefined)
      .map(([k]) => Number(k)).sort((a, b) => a - b);
  };
  OptionsForm.prototype.parse = function (path) {
    // REC.field | IDX[sel].field | $NAME | #NAME | #AMP[sel] | @NAME | sel.<name>
    // %NAME.field | %NAME | %NAME[sel].field: command records outside OPTION_SPECS (sassi/ui/cmdrecords.py)
    if (path.startsWith("%")) {
      const m = /^%(\w+)(?:\[(\w+)\])?(?:\.(\w+))?$/.exec(path);
      if (!m) return null;
      if (m[2]) return {kind: "xindexed", rec: m[1], key: String(this.sel[m[2]]), field: m[3]};
      if (m[3]) return {kind: "xrecords", rec: m[1], field: m[3]};
      return {kind: "xtables", key: m[1]};
    }
    if (path.startsWith("$")) return {kind: "strings", key: path.slice(1)};
    if (path.startsWith("@")) return {kind: "requests", key: path.slice(1)};
    if (path.startsWith("#")) {
      const m = /^#(\w+)(?:\[(\w+)\])?$/.exec(path);
      return {kind: "lists", key: m[2] ? `${m[1]}[${this.sel[m[2]]}]` : m[1]};
    }
    if (path.startsWith("sel.")) return {kind: "sel", key: path.slice(4)};
    const m = /^(\w+)(?:\[(\w+)\])?\.(@?\w+)$/.exec(path);
    if (!m) return null;
    if (m[2]) return {kind: "indexed", rec: m[1], key: String(this.sel[m[2]]), field: m[3]};
    return {kind: "records", rec: m[1], field: m[3]};
  };
  OptionsForm.prototype.get = function (path) {
    const p = this.parse(path);
    const v = this.v;
    if (!p) return undefined;
    if (p.kind === "sel") return this.sel[p.key];
    if (p.kind === "records") return ((v.records || {})[p.rec] || {})[p.field];
    if (p.kind === "indexed") {
      const ent = ((v.indexed || {})[p.rec] || {})[p.key];
      if (ent) return ent[p.field];
      return (this.absentEntry(p.rec, p.key) || {})[p.field];
    }
    if (p.kind === "strings") return (v.strings || {})[p.key] || "";
    if (p.kind === "lists") return (v.lists || {})[p.key] || "";
    if (p.kind === "requests") return (v.requests || {})[p.key] || [];
    if (p.kind === "xrecords") return (v.xrecords[p.rec] || {})[p.field];
    if (p.kind === "xtables") return v.xtables[p.key] || [];
    if (p.kind === "xindexed") {
      // a curve that does not exist (or was deleted in the dialog) shows the defaults
      const ent = (v.xindexed[p.rec] || {})[p.key];
      return (ent || v.xindexed_defaults[p.rec] || {})[p.field];
    }
  };
  /** Delete the entry of a keyed family chosen by its selector (%BBC[bbc]): posted as null (DELBBC). */
  OptionsForm.prototype.deleteEntry = function (path) {
    const p = this.parse(path);
    if (!p || p.kind !== "xindexed") return;
    this.v.xindexed[p.rec] = this.v.xindexed[p.rec] || {};
    this.v.xindexed[p.rec][p.key] = null;
    (this.dirty.xindexed[p.rec] = this.dirty.xindexed[p.rec] || {})[p.key] = true;
    this.markTabDirty();
  };
  OptionsForm.prototype.set = function (path, value) {
    const p = this.parse(path);
    const v = this.v;
    if (!p) return;
    if (p.kind === "sel") { this.sel[p.key] = value; return; }
    if (p.kind === "xrecords") {
      v.xrecords[p.rec] = v.xrecords[p.rec] || {};
      v.xrecords[p.rec][p.field] = value;
      this.dirty.xrecords[p.rec] = true;
      this.markTabDirty();
      return;
    }
    if (p.kind === "xtables") {
      v.xtables[p.key] = value;
      this.dirty.xtables[p.key] = true;
      this.markTabDirty();
      return;
    }
    if (p.kind === "xindexed") {
      v.xindexed[p.rec] = v.xindexed[p.rec] || {};
      if (!v.xindexed[p.rec][p.key]) v.xindexed[p.rec][p.key] = JSON.parse(JSON.stringify(v.xindexed_defaults[p.rec] || {}));
      v.xindexed[p.rec][p.key][p.field] = value;
      (this.dirty.xindexed[p.rec] = this.dirty.xindexed[p.rec] || {})[p.key] = true;
      this.markTabDirty();
      return;
    }
    if (p.kind === "records") {
      v.records[p.rec] = v.records[p.rec] || {};
      const was = v.records[p.rec][p.field];
      v.records[p.rec][p.field] = value;
      this.dirty.records[p.rec] = true;
      // Stochastically Simulated Incoherency Input: Random Phase Angle 180 (5.4 default, D-INC-02)
      if (p.rec === "INCOH" && p.field === "@stoch" && Number(value) !== Number(was)) {
        const inc = v.records.INCOH;
        if (Number(value) === 1 && Number(inc.randphz || 0) === 0) inc.randphz = 180;
        if (Number(value) === 0) Object.assign(inc, {hseed: 0, vseed: 0, randphz: 0});   // zero phases
      }
    } else if (p.kind === "indexed") {
      v.indexed[p.rec] = v.indexed[p.rec] || {};
      if (!v.indexed[p.rec][p.key]) {
        // first edit of an entry that is not stored: a SOIL layer request starts from the 5.4 values
        // (Compute Maximum, Save RS, compute stresses + strains); other entries from what is shown
        const req = (v.request_defaults || {})[p.rec];
        const base = req || this.absentEntry(p.rec, p.key);
        v.indexed[p.rec][p.key] = JSON.parse(JSON.stringify(base || {}));
      }
      v.indexed[p.rec][p.key][p.field] = value;
      this.dirty.indexed[p.rec] = this.dirty.indexed[p.rec] || {};
      this.dirty.indexed[p.rec][p.key] = true;
    } else {
      v[p.kind] = v[p.kind] || {};
      v[p.kind][p.key] = value;
      this.dirty[p.kind][p.key] = true;
    }
    this.markTabDirty();
  };
  /** Values shown for an indexed entry that is not stored: what AFWRITE writes for it (server
   *  dialogs.absent_entry) -- WAVE pages for the <wopt> currently chosen in the dialog. */
  OptionsForm.prototype.absentEntry = function (rec, key) {
    const v = this.v;
    if (rec === "WAVE") {
      const byW = (v.wave_defaults_by_wopt || {})[String(((v.records || {}).SITE || {}).wopt)];
      return (byW || v.wave_defaults || {})[key];
    }
    return (v.indexed_defaults || {})[rec];
  };
  OptionsForm.prototype.anyDirty = function () {
    return Object.values(this.dirty).some((d) => Object.keys(d).length > 0);
  };
  OptionsForm.prototype.payload = function () {
    const out = {records: {}, indexed: {}, strings: {}, lists: {}, requests: {}, xrecords: {}, xtables: {}, xindexed: {}};
    for (const r of Object.keys(this.dirty.records)) out.records[r] = this.v.records[r];
    for (const k of ["indexed", "xindexed"]) {
      for (const [r, keys] of Object.entries(this.dirty[k])) {
        out[k][r] = {};
        for (const key of Object.keys(keys)) out[k][r][key] = this.v[k][r][key];
      }
    }
    for (const k of ["strings", "lists", "requests", "xtables"]) for (const key of Object.keys(this.dirty[k])) out[k][key] = this.v[k][key];
    for (const r of Object.keys(this.dirty.xrecords)) out.xrecords[r] = this.v.xrecords[r];
    return out;
  };
  OptionsForm.prototype.cond = function (c) {
    if (!c) return true;
    if (Array.isArray(c[0])) return c.every((x) => this.cond(x));
    const [path, op, val] = c;
    const cur = this.get(path);
    const num = (x) => (x === "" || x === null || x === undefined || isNaN(Number(x))) ? x : Number(x);
    if (op === "==") return num(cur) === num(val);
    if (op === "!=") return num(cur) !== num(val);
    if (op === "in") return val.map(num).includes(num(cur));
    if (op === "notin") return !val.map(num).includes(num(cur));
    return true;
  };
  OptionsForm.prototype.df = function () {
    const r = (this.v.records || {}).SITE || {};          // (the LOADGEN dialogs have no option records)
    const fs = Number(r.fstep), dt = Number(r.delt), n = Number(r.nft);
    if (fs > 0) return fs;
    if (dt > 0 && n > 0) return 1 / (dt * n);
    return null;
  };
  OptionsForm.prototype.markTabDirty = function () {
    if (this.activeTab && this.tabHeads[this.activeTab]) this.tabHeads[this.activeTab].classList.add("dirty");
  };
  OptionsForm.prototype.build = function (startTab) {
    const root = el("div", {});
    const tabs = this.form.tabs;
    this.tabHeads = {};
    this.pages = {};
    if (tabs.length > 1) {
      const strip = el("div", {class: "opt-tabs", role: "tablist"});
      for (const t of tabs) {
        const h = el("div", {class: "opt-tab", text: t.title, role: "tab"});
        h.addEventListener("click", () => this.showTab(t.name));
        this.tabHeads[t.name] = h;
        strip.appendChild(h);
      }
      root.appendChild(strip);
    }
    this.pageHost = el("div", {});
    root.appendChild(this.pageHost);
    this.showTab(startTab || tabs[0].name);
    return root;
  };
  OptionsForm.prototype.showTab = function (name) {
    this.activeTab = name;
    for (const [n, h] of Object.entries(this.tabHeads)) h.classList.toggle("active", n === name);
    this.renderPage();
  };
  OptionsForm.prototype.renderPage = function () {
    const t = this.form.tabs.find((x) => x.name === this.activeTab);
    this.rows = [];
    const host = this.pageHost;
    host.innerHTML = "";
    if (t.note) host.appendChild(el("div", {class: "opt-note", text: t.note}));
    const page = el("div", {class: "opt-page active"});
    const cols = [el("div", {class: "opt-col"}), el("div", {class: "opt-col"}), el("div", {class: "opt-col"})];
    // one column for the small dialogs; the Analysis window and the multi-column module dialogs
    // (form.columns: ANSYS Eq. Static Load / Dynamic Load, long labels and file boxes: wider columns,
    // two of them in a narrow window) keep the three columns of their groups
    const single = !this.form.columns && (this.form.tabs.length === 1 || !["ANALYSIS"].includes(this.name));
    if (single) page.style.gridTemplateColumns = "minmax(300px, 1fr)";
    else if (this.form.columns) page.style.gridTemplateColumns = "repeat(auto-fit, minmax(330px, 1fr))";
    cols.forEach((c) => page.appendChild(c));
    const full = [];
    for (const g of t.groups) {
      const fs = el("fieldset", {}, g.title ? el("legend", {text: g.title}) : null);
      for (const it of g.items) {
        const r = this.renderItem(it);
        if (r) fs.appendChild(r);
      }
      if (g.note) fs.appendChild(el("div", {class: "note", text: g.note}));
      if (g.full) { full.push(fs); continue; }        // below the columns, across the page (wide tables)
      const ci = single ? 0 : Math.min(g.col || 0, 2);
      cols[ci].appendChild(fs);
    }
    for (const fs of full) page.appendChild(el("div", {class: "opt-full"}, fs));
    host.appendChild(page);
    this.refresh();
  };
  OptionsForm.prototype.rerender = function () {
    const sc = this.pageHost.parentElement && this.pageHost.parentElement.parentElement;
    const top = sc ? sc.scrollTop : 0;
    this.renderPage();
    if (sc) sc.scrollTop = top;
  };
  OptionsForm.prototype.renderItem = function (it) {
    const w = it.w || this.inferWidget(it);
    const F = this;
    const label = it.label !== undefined ? it.label : it.path;
    let node, inputs = [];
    const reg = (n, extra) => { const r = Object.assign({node: n, it, inputs}, extra || {}); F.rows.push(r); return nameInputs(n, label); };
    if (w === "disabled") {
      return el("div", {class: "row", title: it.note || ""}, el("label", {class: "lab", style: {color: "#9aa3ae"}, text: it.label}),
        el("input", {type: "text", value: it.value || "", disabled: true, style: {flex: "0 0 110px"}}));
    }
    if (w === "info") {          // read-only text: a context value or the entries of a keyed family
      const span = el("span", {class: "opt-info"});
      span.textContent = F.infoText(it);
      return reg(el("div", {class: "row", title: it.tip || null}, el("label", {class: "lab", text: label}), span), {info: span});
    }
    if (w === "selector") {
      let input;
      if (it.choices_by) {
        const ch = it.choices_by.map[String(F.get(it.choices_by.path))] || [];
        if (!ch.some(([v]) => Number(v) === Number(F.sel[it.sel]))) F.sel[it.sel] = ch.length ? ch[0][0] : 1;
        input = el("select", {}, ...ch.map(([v, l]) => el("option", {value: v, text: l, selected: Number(v) === Number(F.sel[it.sel])})));
      } else {
        input = el("input", {type: "number", min: it.min || 1, max: it.max || null, value: F.sel[it.sel], style: {width: "80px"}});
      }
      input.addEventListener("change", () => { F.sel[it.sel] = Number(input.value) || 1; F.rerender(); });
      return nameInputs(el("div", {class: "row"}, el("label", {class: "lab", text: label}), input), label);
    }
    if (w === "button") {
      const b = el("button", {class: "btn small", text: it.label, title: it.tip || null});
      b.addEventListener("click", () => F.action(it.action));
      return reg(el("div", {class: "row"}, b));
    }
    if (w === "check" || w === "bit") {
      const c = chk(false);
      const cur = F.get(it.path);
      c.checked = w === "bit" ? (Number(cur) & it.bit) !== 0 : Number(cur) === 1 || cur === true;
      c.addEventListener("change", () => {
        if (w === "bit") { const v = Number(F.get(it.path)) || 0; F.set(it.path, c.checked ? (v | it.bit) : (v & ~it.bit)); }
        else F.set(it.path, c.checked ? 1 : 0);
        F.rerender();
      });
      inputs.push(c);
      const forced = el("span", {class: "opt-forced"});
      node = el("div", {class: "row", title: it.tip || null}, el("label", {class: "chk"}, c, " " + label), forced);
      return reg(node, {forced});
    }
    if (w === "radio" || w === "excl") {
      const box = el("div", {class: "radios" + ((it.choices || []).length <= 3 && w === "radio" && !label ? " inline" : "")});
      const choices = it.choices_by_key ? (it.choices_by_key[String(F.sel.wave)] || []) : it.choices;
      const cur = F.get(it.path);
      const nm = "r" + Math.random().toString(36).slice(2);
      for (const [val, lab] of choices) {
        const inp = el("input", {type: w === "radio" ? "radio" : "checkbox", name: nm, value: val});
        inp.checked = String(cur) === String(val);
        inp.dataset.val = val;
        if ((it.disabled_choices || []).map(String).includes(String(val))) inp.disabled = true;
        inp.addEventListener("change", () => {
          if (w === "excl") F.set(it.path, inp.checked ? val : 0);
          else if (inp.checked) F.set(it.path, val);
          F.rerender();
        });
        inputs.push(inp);
        box.appendChild(el("label", {}, inp, " " + lab));
      }
      node = label ? el("div", {}, el("div", {class: "row"}, el("label", {class: "lab", text: label})), box) : box;
      return reg(el("div", {title: it.tip || null}, node));
    }
    if (w === "select") {
      const cur = F.get(it.path);
      const s = el("select", {title: it.tip || null}, ...it.choices.map(([v, l]) => el("option", {value: v, text: l, selected: String(v) === String(cur)})));
      s.addEventListener("change", () => { F.set(it.path, choiceValue(it.choices, s.value)); F.rerender(); });
      inputs.push(s);
      return reg(el("div", {class: "row"}, el("label", {class: "lab", text: label}), s));
    }
    if (w === "list") {
      const ta = el("textarea", {class: "field", rows: it.rows || 2, style: {width: "100%", font: "12px var(--mono)"}, title: it.tip || null});
      ta.value = F.get(it.path);
      ta.addEventListener("input", () => F.set(it.path, ta.value));
      inputs.push(ta);
      return reg(el("div", {}, el("div", {class: "row"}, el("label", {class: "lab", text: label})), ta, it.tip ? el("div", {class: "note", text: it.tip}) : null));
    }
    if (w === "table") return reg(this.renderTable(it, inputs));
    if (w === "dynp") {
      const inp = txt(F.get(it.path) || "", {class: "path"});
      inp.addEventListener("input", () => F.set(it.path, inp.value));
      inputs.push(inp);
      const b = el("button", {class: "btn small", text: "...", title: "Select Dynamic Soil Property", onclick: () => D.soilProperty({initial: inp.value, select: (lab) => { inp.value = lab; F.set(it.path, lab); }})});
      const dl = el("datalist", {id: "dynp-" + Math.random().toString(36).slice(2)}, ...(F.v.context.dynp || []).map((n) => el("option", {value: n})));
      inp.setAttribute("list", dl.id);
      return reg(el("div", {class: "row"}, el("label", {class: "lab", text: label}), inp, b, dl));
    }
    // text-like inputs: int, real, text, path
    const cur = F.get(it.path);
    const inp = txt(cur === undefined || cur === null ? "" : cur, {class: w === "path" || w === "text" ? "path" : "", placeholder: it.placeholder || null, title: it.tip || null});
    inp.addEventListener("input", () => {
      F.set(it.path, inp.value);
      inp.classList.toggle("invalid", (w === "int" || w === "real") && inp.value.trim() !== "" && isNaN(Number(inp.value)));
      F.refresh();
    });
    inputs.push(inp);
    const extras = [];
    if (w === "path") {
      const pick = () => D.pickFile({title: label, folder: !!it.folder, onOk: (p) => { inp.value = p; F.set(it.path, p); F.refresh(); }});
      const browse = el("button", {class: "btn small", text: "<<", title: it.folder ? "browse (folder)" : "browse", onclick: pick});
      extras.push(browse);
      inputs.push(browse);                    // greyed with its field (enable rules)
      if (it.edit) extras.push(el("button", {class: "btn small", text: "Edit", title: "open in the File Editor", onclick: () => inp.value.trim() && S.openEditor(inp.value.trim())}));
    }
    const hz = it.hz ? el("span", {class: "hz"}) : null;
    return reg(el("div", {class: "row"}, el("label", {class: "lab", text: label}), inp, ...extras, hz), {hz});
  };
  /** Text of an "info" item: the defined entries of a keyed family, or a context value. */
  OptionsForm.prototype.infoText = function (it) {
    if (it.keys_of) {
      const p = this.parse(it.keys_of);
      const keys = p ? this.familyKeys(p.kind === "xtables" ? p.key : p.rec) : [];
      return keys.length ? keys.join(", ") : "none";
    }
    const c = (this.v.context || {})[it.context];
    return c === undefined || c === null ? "" : String(c);
  };
  OptionsForm.prototype.inferWidget = function (it) {
    const p = this.parse(it.path || "");
    if (p && (p.kind === "strings")) return "text";
    const v = this.get(it.path);
    if (p && p.kind === "records" && p.rec === "CHECK") return "int";
    return typeof v === "number" ? (Number.isInteger(v) ? "int" : "real") : "real";
  };
  OptionsForm.prototype.renderTable = function (it, inputs) {
    const F = this;
    const rows = JSON.parse(JSON.stringify(F.get(it.path) || []));
    const groups = (F.v.context && F.v.context.groups) || {};
    const wrap = el("div", {class: it.wide ? "wide-table" : ""});
    const draw = () => {
      wrap.innerHTML = "";
      const tb = el("table", {class: "grid", style: {width: "100%"}});
      const cols = it.columns;
      tb.appendChild(el("thead", {}, el("tr", {}, ...cols.map((c) => el("th", {text: c.label, title: c.label})), el("th", {text: ""}))));
      const body = el("tbody");
      rows.forEach((r, i) => {
        const tds = cols.map((c) => {
          if (c.w === "check") {
            const cb = chk(Number(r[c.name]) >= 1);
            cb.addEventListener("change", () => { r[c.name] = cb.checked ? 1 : 0; commit(); });
            return el("td", {}, cb);
          }
          if (c.w === "select") {
            const s = el("select", {}, ...c.choices.map(([v, l]) => el("option", {value: v, text: l, selected: String(v) === String(r[c.name])})));
            s.addEventListener("change", () => { r[c.name] = choiceValue(c.choices, s.value); commit(); });
            return el("td", {}, s);
          }
          if (c.w === "index") return el("td", {text: String(i + 1)});         // row number (BBC point)
          if (c.w === "codes") {
            const g = groups[String(r.group)] || {components: []};
            const comps = g.components.length ? g.components : Array.from({length: 12}, (_, k) => `c${k + 1}`);
            const codes = (r.codes || []).concat(Array(12).fill(0)).slice(0, 12);
            r.codes = codes;
            const cells = comps.map((nm, k) => {
              const s = el("select", {title: nm}, ...[[0, "0"], [1, "1"], [2, "2"]].map(([v, l]) => el("option", {value: v, text: l, selected: Number(codes[k]) === v})));
              s.addEventListener("change", () => { r.codes[k] = Number(s.value); commit(); });
              return el("span", {style: {display: "inline-flex", flexDirection: "column", alignItems: "center", marginRight: "3px", fontSize: "10px"}}, nm, s);
            });
            return el("td", {class: "l"}, ...cells);
          }
          const inp = txt(r[c.name] === undefined ? "" : r[c.name]);
          inp.addEventListener("input", () => { r[c.name] = inp.value; commit(true); });
          inp.addEventListener("change", () => { if (c.name === "group") draw(); });
          return el("td", {}, inp);
        });
        body.appendChild(el("tr", {}, ...tds, el("td", {}, el("button", {class: "btn small", text: "Delete", onclick: () => { rows.splice(i, 1); commit(); draw(); }}))));
      });
      if (!rows.length) body.appendChild(el("tr", {}, el("td", {colspan: cols.length + 1, class: "l", style: {color: "var(--muted)"}, text: "no entry: Add"})));
      tb.appendChild(body);
      wrap.appendChild(el("div", {class: "tablewrap"}, tb));
      const nextKey = () => Math.max(0, ...rows.map((r) => Number(r[it.key]) || 0)) + 1;
      // Fill (NLSLAYER): one row per number of a context list (the SOIL sublayers) that has none yet
      const fill = it.fill ? el("button", {class: "btn small", text: it.fill.label, onclick: () => {
        const have = new Set(rows.map((r) => Number(r[it.key])));
        for (const n of ((F.v.context || {})[it.fill.context] || [])) {
          if (!have.has(Number(n))) rows.push(Object.assign({}, it.fill.row || {}, {[it.key]: n}));
        }
        rows.sort((a, b) => (Number(a[it.key]) || 0) - (Number(b[it.key]) || 0));
        commit();
        draw();
      }}) : null;
      wrap.appendChild(el("div", {class: "row"}, el("button", {class: "btn small", text: "Add", onclick: () => {
        const nr = {};
        for (const c of it.columns) {
          nr[c.name] = c.default !== undefined ? c.default : c.w === "check" ? 0 : c.w === "select" ? (c.choices[0] || [1])[0] : c.w === "codes" ? Array(12).fill(0) : "";
        }
        if (it.key) nr[it.key] = nextKey();
        rows.push(nr);
        commit();
        draw();
      }}), fill, el("span", {class: "note", text: it.note || (it.path === "@CORR" ? "CORR rows cannot be deleted by command (stored rows stay)" : "OK replaces the list (clear + one command per row)")})));
    };
    const commit = () => { F.set(it.path, rows); F.refresh(); };
    draw();
    return el("div", {}, it.label ? el("div", {class: "row"}, el("label", {class: "lab", text: it.label})) : null, wrap);
  };
  OptionsForm.prototype.action = async function (action) {
    const ctx = this.v.context || {};
    if (action === "radius") {
      try {
        const r = await S.post("/api/options/POINT/radius");
        if (r.ok && r.radius) {
          this.set("POINT.rad", Number(r.radius.average.toPrecision(6)));
          this.rerender();
          this.modal && this.modal.setMessage(`RADIUS: min ${S.fmt(r.radius.min)}, average ${S.fmt(r.radius.average)}, max ${S.fmt(r.radius.max)}`, true);
        } else this.modal && this.modal.setMessage((r.messages || []).filter((m) => m.kind === "ERROR").map((m) => m.text).join("\n") || "RADIUS failed");
      } catch (e) { this.modal && this.modal.setMessage(e.message); }
    } else if (action.startsWith("edit:")) {
      let f = action.slice(5).replace("{model}", ctx.model || "model");
      if (ctx.path) f = ctx.path.replace(/\/$/, "") + "/" + f;
      S.openEditor(f);
    } else if (action.startsWith("newkey:")) {     // New Curve: the selector shows the next free number
      const [, sel, fam] = action.split(":");
      const p = this.parse(fam);
      const keys = p ? this.familyKeys(p.key) : [];
      this.sel[sel] = keys.length ? keys[keys.length - 1] + 1 : 1;
      this.rerender();
    } else if (action.startsWith("delentry:")) {   // Delete Curve: posted as null (DELBBC on OK)
      this.deleteEntry(action.slice(9));
      this.rerender();
    }
  };
  /** Re-evaluate enable/visible rules (D-UI-04), forced values and the Hz labels. */
  OptionsForm.prototype.refresh = function () {
    const df = this.df();
    for (const r of this.rows) {
      const it = r.it;
      if (r.info) r.info.textContent = this.infoText(it);
      let en = !it.disabled && this.cond(it.enable);
      const vis = this.cond(it.visible);
      r.node.style.display = vis ? "" : "none";
      if (it.force_on) {
        this.forcedPrev = this.forcedPrev || {};
        const why = it.force_on.filter((c) => this.cond(c));
        if (why.length) {
          if (Number(this.get(it.path)) !== 1) {
            this.forcedPrev[it.path] = this.get(it.path);       // restored when the reason goes away
            this.set(it.path, 1);
          }
          r.inputs.forEach((x) => { x.checked = true; });
          en = false;
          if (r.forced) r.forced.textContent = "(required: unlagged coherency model 2-7 or multiple excitation)";
        } else {
          if (this.forcedPrev[it.path] !== undefined) {
            this.set(it.path, this.forcedPrev[it.path]);
            r.inputs.forEach((x) => { x.checked = Number(this.forcedPrev[it.path]) === 1; });
            delete this.forcedPrev[it.path];
          }
          if (r.forced) r.forced.textContent = "";
        }
      }
      for (const inp of r.inputs) {
        let dis = !en;
        if (it.disabled_choices && (it.disabled_choices.map(String).includes(inp.dataset.val))) dis = true;
        if (it.disable_choice_if && inp.dataset.val !== undefined) {
          const c = it.disable_choice_if[inp.dataset.val];
          if (c && this.cond(c)) {
            dis = true;
            if (inp.checked) { this.set(it.path, 0); r.inputs.forEach((x) => { x.checked = x.dataset.val === "0"; }); }
          }
        }
        inp.disabled = dis;
      }
      if (r.hz) {
        const raw = this.get(it.path);
        // a blank Frequency 2 is NFFT/2 (AFWRITE), so its Hz follows the Nr. of Fourier Components
        const n = (raw === "" || raw === null || raw === undefined) && it.blank_half
          ? Math.max(Math.floor(Number(this.get(it.blank_half)) / 2), 1) : Number(raw);
        r.hz.textContent = df && !isNaN(n) ? `${(n * df).toPrecision(5)} Hz` +
          ((raw === "" || raw === null) && it.blank_half ? " (NFFT/2)" : "") : "";
      }
    }
  };
})(SASSI);
