"""
Written by Sean Engelstad, GT SMDO Lab, 2022-2023
"""
__all__ = ["AflrAim"]


class AflrAim:
    def __init__(self, caps_problem, comm, root=0):
        """
        MPI wrapper class for AflrAIM from ESP/CAPS.
        """

        self.caps_problem = caps_problem
        self.comm = comm
        self.root = root

        # holds aflr4 aim
        self._aim = None

        self._dictOptions = None

        self._build_aim()
        return

    @property
    def root_proc(self) -> bool:
        return self.comm.rank == self.root

    @property
    def aim(self):
        """surface mesher aim aka aflr4 aim"""
        return self._aim

    @property
    def analysis_dir(self):
        _analysis_dir = None
        if self.comm.rank == self.root:
            _analysis_dir = self.aim.analysisDir
        _analysis_dir = self.comm.bcast(_analysis_dir, root=self.root)
        return _analysis_dir

    def _build_aim(self):
        if self.root_proc:
            self._aim = self.caps_problem.analysis.create(aim="aflr4AIM", name="aflr4")
        return

    def set_mesh(self, min_scale=0.05, max_scale=0.5, AFLR4_Quad=False, no_prox=True):
        """
        Set mesh properties for AFLR4 AIM. A few options are available in this routine.
        To set other options for AFLR4, use the save_dict_options routine.

        Parameters
        ----------
        min_scale: Relative scale of minimum spacing to reference length.
            The relative scale of minimum spacing to reference length (ref_len) controls
            the minimum spacing that can be set on any component/body surface.
        max_scale: Relative scale of maximum spacing to reference length.
            The relative scale of maximum spacing to reference length (ref_len) controls
            the maximum spacing that can be set on any component/body surface.
        AFLR4_Quad: Generate a mixed quad/tria-face grid.
        no_prox: Disable proximity checking.
            Proximity checking is automatically disabled if there is only one component/body defined.
        """

        if self.root_proc:
            self.aim.input.min_scale = min_scale
            self.aim.input.max_scale = max_scale
            self.aim.input.AFLR4_Quad = AFLR4_Quad
            self._aim.input.no_prox = no_prox
        return self

    def save_dict_options(self, dictOptions: dict = None):
        """
        Optional method to set AFLR4 mesh settings using dictionaries.
        Call this before setting up the TACS model. The dictionary should take
        the form of, e.g.:

        dictOptions['aflr4AIM']['myOption'] = myValue
        """
        self._dictOptions = dictOptions

        return self

    def _set_dict_options(self):
        """
        Set AFLR4 options via dictionaries.
        """
        dictOptions = self._dictOptions

        if self.root_proc:
            for ind, option in enumerate(dictOptions["aflr4AIM"]):
                self.aim.input[option].value = dictOptions["aflr4AIM"][option]

        return self

    def register_to(self, tacs_aim):
        """
        cascade method to register the egads aim to the tacs aim wrapper class
        """
        tacs_aim.register(self)
        return self
