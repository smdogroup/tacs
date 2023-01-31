__all__ = ["EgadsAim"]

import pyCAPS
from typing import TYPE_CHECKING, List

class EgadsAim:
    """
    only for structure mesh
    egadsAim.input.Edge_Point_Min = 15
    egadsAim.input.Edge_Point_Max = 20
    egadsAim.input.Mesh_Elements = "Quad"
    egadsAim.input.Tess_Params = [.25,.01,15]
    """
    def __init__(self, caps_problem:pyCAPS.Problem):
        self._aim = caps_problem.analysis.create(aim="egadsTessAIM")
        self._is_setup = False

    def set_mesh(self, edge_pt_min:int=15, edge_pt_max=20, mesh_elements:str="Quad", global_mesh_size:float=0.25, max_surf_offset:float=0.01, max_dihedral_angle:float=15):
        self._aim.input.Edge_Point_Min = edge_pt_min
        self._aim.input.Edge_Point_Max = edge_pt_max
        self._aim.input.Mesh_Elements = mesh_elements
        self._aim.input.Tess_Params = [global_mesh_size, max_surf_offset, max_dihedral_angle]
        self._is_setup = True

    @property
    def is_setup(self) -> bool:
        return self._is_setup

    @property
    def aim(self):
        return self._aim