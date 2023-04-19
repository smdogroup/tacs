"""
pyBase_constraint
"""

# =============================================================================
# Imports
# =============================================================================
from collections import OrderedDict

from ..system import TACSSystem


class TACSConstraint(TACSSystem):
    """
    Base class for TACS constraint types. Contains methods common to all TACS constraints.
    """

    def __init__(
        self, assembler, comm=None, options=None, outputViewer=None, meshLoader=None
    ):
        # Set attributes and options
        TACSSystem.__init__(self, assembler, comm, options, outputViewer, meshLoader)

        # List of constraints
        self.constraintList = OrderedDict()

        return

    ####### Eval constraint methods ########

    def addConstraint(self, conName, lower=-1e20, upper=1e20, compIDs=None, dvIndex=0):
        """
        Generic method to adding a new constraint set for TACS.

        Parameters
        ----------
        conName : str
            The user-supplied name for the constraint set. This will
            typically be a string that is meaningful to the user

        lower: float or complex
            lower bound for constraint. Defaults to -1e20.

        upper: float or complex
            upper bound for constraint. Defaults to 1e20.

        compIDs: list or None
            List of compIDs to select. If None, all compIDs will be selected. Defaults to None.

        dvIndex : int
            Index number of element DV to be used in constraint. Defaults to 0.

        """
        raise NotImplemented(
            f"'addConstraint' method is not implemented for class '{type(self).__name__}'"
        )

    def getConstraintBounds(self, bounds, evalCons=None, ignoreMissing=False):
        """
        Get bounds for constraints. The constraints corresponding to the strings in
        `evalCons` are evaluated and updated into the provided
        dictionary.

        Parameters
        ----------
        bounds : dict
            Dictionary into which the constraint bounds are saved.
            Bounds will be saved as a tuple: (lower, upper)
        evalCons : iterable object containing strings.
            If not none, use these constraints to evaluate.
        ignoreMissing : bool
            Flag to supress checking for a valid constraint. Please use
            this option with caution.

        Examples
        --------
        >>> funcs = {}
        >>> tacsConstraint.getConstraintBounds(funcs, 'LE_SPAR')
        >>> funcs
        >>> # Result will look like (if TACSConstraint has name of 'c1'):
        >>> # {'c1_LE_SPAR': (array([-1e20]), array([1e20]))}
        """
        raise NotImplemented(
            f"'getConstraintBounds' method is not implemented for class '{type(self).__name__}'"
        )

    def getConstraintKeys(self):
        """
        Return a list of the current constraint key names
        """
        return list(self.constraintList.keys())

    def evalConstraints(self, funcs, evalCons=None, ignoreMissing=False):
        """
        Evaluate values for constraints. The constraints corresponding to the strings in
        evalCons are evaluated and updated into the provided
        dictionary.

        Parameters
        ----------
        funcs : dict
            Dictionary into which the constraints are saved.
        evalCons : iterable object containing strings.
            If not none, use these constraints to evaluate.
        ignoreMissing : bool
            Flag to supress checking for a valid constraint. Please use
            this option with caution.

        Examples
        --------
        >>> funcs = {}
        >>> tacsConstraint.evalConstraints(funcs, 'LE_SPAR')
        >>> funcs
        >>> # Result will look like (if TACSConstraint has name of 'c1'):
        >>> # {'c1_LE_SPAR': array([12354.10])}
        """
        raise NotImplemented(
            f"'evalConstraints' method is not implemented for class '{type(self).__name__}'"
        )

    def evalConstraintsSens(self, funcsSens, evalCons=None):
        """
        This is the main routine for returning useful (sensitivity)
        information from constraint. The derivatives of the constraints
        corresponding to the strings in evalCons are evaluated and
        updated into the provided dictionary. The derivitives with
        respect to all design variables and node locations are computed.

        Parameters
        ----------
        funcsSens : dict
            Dictionary into which the derivatives are saved.
        evalCons : iterable object containing strings
            The constraints the user wants returned

        Examples
        --------
        >>> funcsSens = {}
        >>> tacsConstraint.evalConstraintsSens(funcsSens, 'LE_SPAR')
        >>> funcsSens
        >>> # Result will look like (if TACSConstraint has name of 'c1'):
        >>> # {'c1_LE_SPAR':{'struct':<50x242 sparse matrix of type '<class 'numpy.float64'>' with 100 stored elements in Compressed Sparse Row format>}}
        """
        raise NotImplemented(
            f"'evalConstraintsSens' method is not implemented for class '{type(self).__name__}'"
        )
