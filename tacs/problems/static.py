"""
pyStatic_problem
"""

# =============================================================================
# Imports
# =============================================================================
import warnings


class StaticProblem:
    """
    The main purpose of this class is to represent all relevant
    information for a static analysis. This will include
    information defining the loading condition as well as various
    other pieces of information.

    Parameters
    ----------
    name : str
        Name of this tacs problem

    Examples
    --------
    >>> sp = StaticProblem('lc0')
    """

    def __init__(self, name, **kwargs):

        # Always have to have the name
        self.name = name

        # Defaults
        self.loadFile = None
        self.loadFactor = 1.0

        if "loadFactor" in kwargs:
            self.loadFactor = kwargs["loadFactor"]

        # Check for function list:
        self.evalFuncs = set()
        if "evalFuncs" in kwargs:
            self.evalFuncs = set(kwargs["evalFuncs"])
        if "funcs" in kwargs:
            warnings.warn("funcs should **not** be an argument. Use 'evalFuncs' instead.")
            if self.evalFuncs is None:
                self.evalFuncs = set(kwargs["funcs"])

        # we cast the set to a sorted list, so that each proc can loop over in the same order
        self.evalFuncs = sorted(self.evalFuncs)

        # When a solver calls its evalFunctions() it must write the
        # unique name it gives to funcNames.
        self.funcNames = {}
        self.possibleFunctions = set()


    def __getitem__(self, key):

        return self.funcNames[key]
