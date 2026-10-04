"""Model JSON for the GUI (``GET /api/model``): what the 3D plots and the dialogs need.

The arrays are in **global** coordinates (local coordinate systems resolved, D-MDL-04): node ids,
coordinates, fixity codes (D: 1 = fixed, per DOF UX UY UZ ROTX ROTY ROTZ), interaction-node flags
(INT code 0 = interaction node; only code 0 affects the solution, D-MDL-07), elements by group with
their type, material (MSET) and property (RSET) numbers and ETYPE (0 implicit, 1 structure, 2
excavated soil), and the lumped masses MT / MR in the model's mass units (MUNITS).

The drawing data of an open plot (faces, colours, labels, camera) come from
:func:`sassi.plotting.state.plot_data` (``GET /api/plot/<id>``); this module is the raw model.
"""
from __future__ import annotations

from typing import Any, Dict, List

import numpy as np

from ..conventions import ELEMENT_TYPE_NAMES
from ..model.entities import INT_CODES


def _f(v: float) -> float:
    v = float(v)
    return v if np.isfinite(v) else 0.0


def model_json(interp, number: int = -1, with_elements: bool = True) -> Dict[str, Any]:
    """JSON description of model ``number`` (-1 = the active model)."""
    n = interp.active_model if number is None or number < 0 else number
    if n not in interp.models:
        raise KeyError(f"model {n} is not in memory")
    m = interp.models[n]
    ids = sorted(m.nodes)
    if ids:
        gid, xyz = m.global_coordinates(ids)
        xyz = np.asarray(xyz, dtype=float).reshape(-1, 3)
        gid = [int(v) for v in gid]
    else:
        gid, xyz = [], np.zeros((0, 3))
    fix = [list(m.nodes[i].fix) for i in gid]
    flags = {str(i): sorted(int(f) for f in m.nodes[i].flags) for i in gid if m.nodes[i].flags}
    interaction = [i for i in gid if 0 in m.nodes[i].flags]
    groups: List[Dict[str, Any]] = []
    for g in sorted(m.groups):
        grp = m.groups[g]
        d: Dict[str, Any] = {"id": grp.id, "type": grp.type, "type_name": ELEMENT_TYPE_NAMES.get(grp.type, str(grp.type)),
                             "title": grp.title, "n": len(grp.elements)}
        if with_elements:
            els = grp.sorted_elements()
            d["elements"] = {"id": [e.id for e in els], "nodes": [list(e.nodes) for e in els],
                             "mat": [e.mat for e in els], "prop": [e.prop for e in els],
                             "etype": [e.etype for e in els]}
        groups.append(d)
    mass_nodes = sorted(set(m.tmass) | set(m.rmass))
    masses = {"node": mass_nodes,
              "tmass": [[_f(v) for v in (m.tmass.get(i) or [0, 0, 0])[:3]] for i in mass_nodes],
              "rmass": [[_f(v) for v in (m.rmass.get(i) or [0, 0, 0])[:3]] for i in mass_nodes]}
    if len(xyz):
        lo, hi = xyz.min(axis=0), xyz.max(axis=0)
        bbox = [_f(lo[0]), _f(hi[0]), _f(lo[1]), _f(hi[1]), _f(lo[2]), _f(hi[2])]
    else:
        bbox = [0.0] * 6
    return {
        "number": n, "name": m.name, "path": m.path, "title": m.title, "counts": m.counts(),
        "nodes": {"id": gid, "xyz": [[_f(a), _f(b), _f(c)] for a, b, c in xyz], "fix": fix},
        "flags": flags, "flag_names": {str(k): v for k, v in INT_CODES.items()}, "interaction": interaction,
        "groups": groups, "masses": masses, "bbox": bbox,
        "materials": sorted(m.materials), "layers": sorted(m.layers), "sections": sorted(m.sections),
        "springs": sorted(m.springs), "matrices": sorted(m.matrices),
        "freq_sets": {str(k): list(v) for k, v in sorted(m.freq_sets.items())},
        "gravity": _f(m.gravity), "ground_elevation": _f(m.ground_elevation),
        "ui_state": m.ui_state,
    }


def models_summary(interp) -> List[Dict[str, Any]]:
    """The models in memory (Model menu, status bar)."""
    out = []
    for n in sorted(interp.models):
        m = interp.models[n]
        c = m.counts()
        out.append({"number": n, "name": m.name, "path": m.path, "title": m.title, "active": n == interp.active_model,
                    "nodes": c.get("nodes", 0), "elements": c.get("elements", 0), "groups": c.get("groups", 0)})
    return out
