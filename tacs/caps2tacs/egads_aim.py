__all__ = ["EgadsAim"]

from typing import TYPE_CHECKING, List


class EgadsAim:
    """
    Wrapper class for ESP/CAPS EgadsAim to build structure mesh for TacsAim
    Controls the following inputs:
    egadsAim.input.Edge_Point_Min = 15
    egadsAim.input.Edge_Point_Max = 20
    egadsAim.input.Mesh_Elements = "Quad"
    egadsAim.input.Tess_Params = [.25,.01,15]
    """

    def __init__(self, caps_problem, comm):
        self.comm = comm

        self._dictOptions = None

        if comm is None or comm.rank == 0:
            self._aim = caps_problem.analysis.create(aim="egadsTessAIM")
        self._is_setup = False

    def set_mesh(
        self,
        edge_pt_min: int = 15,
        edge_pt_max=20,
        mesh_elements: str = "Quad",
        global_mesh_size: float = 0.25,
        max_surf_offset: float = 0.01,
        max_dihedral_angle: float = 15,
    ):
        """
        cascaded method to set the mesh input settings to the egadsAim
        """
        if self.comm.rank == 0:
            self._aim.input.Edge_Point_Min = edge_pt_min
            self._aim.input.Edge_Point_Max = edge_pt_max
            self._aim.input.Mesh_Elements = mesh_elements
            self._aim.input.Tess_Params = [
                global_mesh_size,
                max_surf_offset,
                max_dihedral_angle,
            ]
        self._is_setup = True
        return self

    def save_dict_options(self, dictOptions: dict = None):
        """
        Optional method to set EGADS mesh settings using dictionaries.
        Call this before setting up the TACS model. The dictionary should take
        the form of, e.g.:

        dictOptions['egadsTessAIM']['myOption'] = myValue
        """
        self._dictOptions = dictOptions

        return self

    def _set_dict_options(self):
        """
        Set EGADS options via dictionaries.
        """
        dictOptions = self._dictOptions

        if self.root_proc:
            for ind, option in enumerate(dictOptions["egadsTessAIM"]):
                self.aim.input[option].value = dictOptions["egadsTessAIM"][option]

        return self

    @property
    def is_setup(self) -> bool:
        return self._is_setup

    @property
    def aim(self):
        return self._aim

    def register_to(self, tacs_aim):
        """
        cascade method to register the egads aim to the tacs aim wrapper class
        """
        tacs_aim.register(self)
        return self
