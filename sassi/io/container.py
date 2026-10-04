"""Binary inter-module file containers (FILE1, FILE2, FILE3, FILE4, COOSK, FILE8, FILE9 ...).

Decision D-GEN-03: every binary file is a NumPy ``.npz`` archive written under the
legacy file name *without* an added extension (e.g. ``FILE8``).  The archive holds a
``__meta__`` entry (JSON text) with at least::

    {"format": "SASSI-EDU", "kind": "FILE8", "version": 1, "module": "ANALYS", ...}

Arrays are stored uncompressed (fast, inspectable).  Sparse matrices are stored as
``<name>__data``, ``<name>__indices``, ``<name>__indptr``, ``<name>__shape`` (CSR).

The array names and shapes of each kind are the normative contract documented in
docs/ARCHITECTURE.md section 4 and checked by :data:`sassi.io.files.SCHEMAS`.
"""
from __future__ import annotations

import io
import json
import os
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional, Union

import numpy as np
import scipy.sparse as sp

PathLike = Union[str, os.PathLike]
FORMAT = "SASSI-EDU"


@dataclass
class Container:
    """In-memory view of a container file."""

    kind: str
    meta: Dict[str, Any] = field(default_factory=dict)
    arrays: Dict[str, np.ndarray] = field(default_factory=dict)

    def __getitem__(self, name: str) -> np.ndarray:
        return self.arrays[name]

    def __contains__(self, name: str) -> bool:
        return name in self.arrays or (name + "__data") in self.arrays

    def get(self, name: str, default=None):
        return self.arrays.get(name, default)

    def sparse(self, name: str) -> sp.csr_matrix:
        return get_sparse(self.arrays, name)


def put_sparse(arrays: Dict[str, np.ndarray], name: str, mat) -> None:
    """Store a scipy sparse matrix in ``arrays`` under ``name`` (CSR layout)."""
    m = sp.csr_matrix(mat)
    m.sort_indices()
    arrays[name + "__data"] = m.data
    arrays[name + "__indices"] = m.indices.astype(np.int64)
    arrays[name + "__indptr"] = m.indptr.astype(np.int64)
    arrays[name + "__shape"] = np.asarray(m.shape, dtype=np.int64)


def get_sparse(arrays: Dict[str, np.ndarray], name: str) -> sp.csr_matrix:
    shape = tuple(int(v) for v in arrays[name + "__shape"])
    return sp.csr_matrix(
        (arrays[name + "__data"], arrays[name + "__indices"], arrays[name + "__indptr"]), shape=shape
    )


def write_container(path: PathLike, kind: str, arrays: Dict[str, Any], meta: Optional[Dict[str, Any]] = None,
                    module: str = "", version: int = 1) -> Path:
    """Write ``arrays`` and ``meta`` to ``path`` (exact name, no extension added).

    Values in ``arrays`` may be numpy arrays, Python scalars/lists (converted with
    ``np.asarray``) or scipy sparse matrices (stored with :func:`put_sparse`).
    The file is written atomically (temporary file + rename).
    """
    path = Path(path)
    out: Dict[str, np.ndarray] = {}
    for k, v in arrays.items():
        if sp.issparse(v):
            put_sparse(out, k, v)
        else:
            a = np.asarray(v)
            if a.dtype == object:
                raise TypeError(f"array '{k}' has object dtype; use numeric or fixed-width string arrays")
            out[k] = a
    m = {"format": FORMAT, "kind": kind, "version": version, "module": module}
    if meta:
        m.update(meta)
    buf = io.BytesIO()
    np.savez(buf, __meta__=np.asarray(json.dumps(m, default=_json_default)), **out)
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "wb") as fh:
        fh.write(buf.getvalue())
    os.replace(tmp, path)
    return path


def read_container(path: PathLike, kind: Optional[str] = None) -> Container:
    """Read a container; if ``kind`` is given it must match the stored kind."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(str(path))
    if not zipfile.is_zipfile(path):
        raise ValueError(f"{path.name} is not a SASSI-EDU binary container")
    with np.load(path, allow_pickle=False) as z:
        meta = json.loads(str(z["__meta__"]))
        arrays = {k: z[k] for k in z.files if k != "__meta__"}
    if meta.get("format") != FORMAT:
        raise ValueError(f"{path.name}: not a {FORMAT} file")
    if kind is not None and meta.get("kind") != kind:
        raise ValueError(f"{path.name}: expected kind {kind}, found {meta.get('kind')}")
    return Container(kind=meta.get("kind", ""), meta=meta, arrays=arrays)


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, complex):
        return [o.real, o.imag]
    if isinstance(o, Path):
        return str(o)
    raise TypeError(f"not JSON serialisable: {type(o)}")
