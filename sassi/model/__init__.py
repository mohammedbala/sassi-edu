"""SASSI-EDU model database (ARCHITECTURE section 8).

The interpreter (:mod:`sassi.prep`) edits :class:`SSIModel` objects; AFWRITE reads them to build
the module decks.  Modules never import this package (ARCHITECTURE principle 2).

Overview::

    from sassi.model import SSIModel
    m = SSIModel()
    m.define_node(1, (0.0, 0.0, 0.0), csys=0)
    m.node_global(1)                     # global coordinates (local systems resolved)
    m.global_coordinates()               # vectorised (ids, xyz)
    m.element_node_table(1)              # SOLID element node table
    m.nodal_masses()                     # MT/MR converted to mass units (MUNITS)
    m.options.record("SITE")             # analysis option record (generic or typed)
    m.model_hash()                       # content hash
"""
from .entities import (DOF_LABELS, DOF_NAMES, ELEMENT_DOFS, INT_CODES, INT_LETTERS, MATRIX_KINDS,
                       MAX_ELEMENT_NODES, MIN_ELEMENT_NODES, BeamSection, CoordSys, Element, ElementRequest,
                       Group, MatrixProp, Material, NodalLoad, NodalRequest, Node, RelDispRequest, SoilLayer,
                       SpringProp, type_name)
from .materials import (ElasticConstants, elastic_constants, layer_constants, section_circle,
                        section_rectangle, to_mass)
from .options import (RECORD_TYPES, Field, MoptRecord, OptionStore, Record, make_record,
                      register_record_type)
from .ssimodel import SCHEMA_VERSION, History, SSIModel
from .values import fmt_num, is_number, parse_float, parse_int

__all__ = [
    "SSIModel", "SCHEMA_VERSION", "History", "Node", "CoordSys", "Group", "Element", "Material", "SoilLayer",
    "BeamSection", "SpringProp", "MatrixProp", "NodalLoad", "NodalRequest", "ElementRequest",
    "RelDispRequest", "OptionStore", "Record", "Field", "MoptRecord", "RECORD_TYPES", "make_record",
    "register_record_type", "DOF_LABELS", "DOF_NAMES", "ELEMENT_DOFS", "INT_CODES", "INT_LETTERS",
    "MATRIX_KINDS", "MAX_ELEMENT_NODES", "MIN_ELEMENT_NODES", "type_name", "ElasticConstants",
    "elastic_constants", "layer_constants", "section_circle", "section_rectangle", "to_mass",
    "fmt_num", "is_number", "parse_float", "parse_int",
]
