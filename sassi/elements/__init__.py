"""HOUSE finite-element library (ARCHITECTURE section 6.2, requirements 1.5 and 4.1).

Importing the package registers every element type in :data:`ELEMENTS`:

====  =======  ============================================================================
code  name     formulation
====  =======  ============================================================================
1     SOLID    8-node hexahedron (prism/pyramid by repeated nodes), optional 9 incompatible
               modes, 1/2 lumped + 1/2 consistent mass
2     BEAMS    3D Timoshenko frame (I, J; K orientation), end releases, consistent mass
3     SHELL    flat Kirchhoff facet: Q4 + incompatible modes / CST membrane, DKQ / DKT
               bending, no drilling stiffness, lumped mass
4     PLANE    4-node plane strain in the X-Z plane, optional 4 incompatible modes
5     TSHELL   flat Mindlin-Reissner (thick) shell: Q4 + incompatible / CST membrane,
               MITC4 / MITC3 transverse shear, EINT 0 reduced (1-point + hourglass) / 1
               selective (2x2) bending, automatic drilling stiffness, lumped mass with
               rotary inertia (P2, D-ELM-08)
7     SPRING   six uncoupled global springs with hysteretic damping
9     GENERAL  user 12x12 complex stiffness and mass (2 nodes global / 3 nodes local)
====  =======  ============================================================================

"""
from __future__ import annotations

from .base import (ELEMENTS, ElementError, ElementSpec, MaterialProps, material_from_E_nu,  # noqa: F401
                   material_from_layer, material_from_M, nodal_mass_to_mass_units, section_circle,
                   section_rectangle)
from . import solid, beam, shell, plane, spring, general, tshell  # noqa: F401,E402  (registration)
from .assemble import (AssembledModel, DofMap, ElemRecord, assemble, build_dofmap,  # noqa: F401,E402
                       natural_frequencies, nodal_masses_from_table, unstiffened_rotations)

__all__ = ["ELEMENTS", "ElementError", "ElementSpec", "MaterialProps", "material_from_M", "material_from_layer",
           "material_from_E_nu", "nodal_mass_to_mass_units", "section_circle", "section_rectangle",
           "ElemRecord", "DofMap", "AssembledModel", "build_dofmap", "assemble", "natural_frequencies",
           "nodal_masses_from_table", "unstiffened_rotations"]
