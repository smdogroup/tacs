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

    def set_mesh(
        self, ff_growth=1.4, min_scale=0.05, max_scale=0.5, use_quad=False, no_prox=True
    ):
        # set surface mesh properties
        if self.root_proc:
            self.aim.input.ff_cdfr = ff_growth
            self.aim.input.min_scale = min_scale
            self.aim.input.max_scale = max_scale
            self.aim.input.AFLR4_Quad = use_quad
            self._aim.input.no_prox = no_prox
        return self
    
    def saveDictOptions(self, dictOptions: dict = None):
        """
        Optional method to set AFLR4 mesh settings using dictionaries.
        Call this before setting up the TACS model. The dictionary should take
        the form of, e.g.: 

        dictOptions['aflr4AIM']['myOption'] = myValue
        """
        self._dictOptions = dictOptions

        return self
    
    def _setDictOptions(self):
        """
        Set AFLR4 options via dictionaries.
        """
        dictOptions = self._dictOptions

        for ind, option in enumerate(dictOptions["aflr4AIM"]):
            self.aim.input[option].value = dictOptions["aflr4AIM"][option]

        return self

    def register_to(self, tacs_aim):
        """
        cascade method to register the egads aim to the tacs aim wrapper class
        """
        tacs_aim.register(self)
        return self
