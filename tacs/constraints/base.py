"""
pyBase_constraint
"""

# =============================================================================
# Imports
# =============================================================================
from collections import OrderedDict

import numpy as np
import scipy as sp

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

        # Setup global to local dv num map for each proc
        self._initilaizeGlobalToLocalDVDict()

        return

    def _initilaizeGlobalToLocalDVDict(self):
        size = self.comm.size
        rank = self.comm.rank
        nLocalDVs = self.getNumDesignVars()
        nLocalDVsOnProc = self.comm.allgather(nLocalDVs)
        # Figure out which DVNums belong to each processor
        ownerRange = np.zeros(size + 1, dtype=int)
        # Sum local dv ranges over each proc to get global dv ranges
        ownerRange[1:] = np.cumsum(nLocalDVsOnProc)
        self.globalToLocalDVNums = dict(
            zip(range(ownerRange[rank], ownerRange[rank + 1]), range(nLocalDVs))
        )

    ####### Eval constraint methods ########

    def addConstraint(self, conName, compIDs=None, lower=-1e20, upper=1e20, **kwargs):
        """
        Generic method to adding a new constraint set for TACS.

        Parameters
        ----------
        conName : str
            The user-supplied name for the constraint set. This will
            typically be a string that is meaningful to the user

        compIDs: list[int] or None
            List of compIDs to select. If None, all compIDs will be selected. Defaults to None.

        lower: float or complex
            lower bound for constraint. Defaults to -1e20.

        upper: float or complex
            upper bound for constraint. Defaults to 1e20.

        """
        raise NotImplemented(
            f"'addConstraint' method is not implemented for class '{type(self).__name__}'"
        )

    def getConstraintBounds(self, bounds, evalCons=None):
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

        Examples
        --------
        >>> conBounds = {}
        >>> tacsConstraint.getConstraintBounds(conBounds, 'LE_SPAR')
        >>> conBounds
        >>> # Result will look like (if TACSConstraint has name of 'c1'):
        >>> # {'c1_LE_SPAR': (array([-1e20]), array([1e20]))}
        """
        # Check if user specified which constraints to output
        # Otherwise, output them all
        evalCons = self._processEvalCons(evalCons)

        # Loop through each requested constraint set
        for conName in evalCons:
            key = f"{self.name}_{conName}"
            if hasattr(self.constraintList[conName], "getBounds"):
                bounds[key] = self.constraintList[conName].getBounds()
            else:
                bounds[key] = (None, None)

    def getConstraintSizes(self, sizes, evalCons=None):
        """
        Get number for constraint equations in each set.
        The constraints corresponding to the strings in `evalCons`
        are evaluated and updated into the provided dictionary.

        Parameters
        ----------
        sizes : dict
            Dictionary into which the constraint sizes are saved.
        evalCons : iterable object containing strings.
            If not none, use these constraints to evaluate.

        Examples
        --------
        >>> conSizes = {}
        >>> tacsConstraint.getConstraintSizes(conSizes, 'LE_SPAR')
        >>> funconSizescs
        >>> # Result will look like (if TACSConstraint has name of 'c1'):
        >>> # {'c1_LE_SPAR': 10}
        """
        # Check if user specified which constraints to output
        # Otherwise, output them all
        evalCons = self._processEvalCons(evalCons)

        # Loop through each requested constraint set
        for conName in evalCons:
            key = f"{self.name}_{conName}"
            if hasattr(self.constraintList[conName], "nCon"):
                sizes[key] = self.constraintList[conName].nCon
            else:
                sizes[key] = 0

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

    def _processEvalCons(self, evalCons, ignoreMissing=True):
        """
        Internal method for processing user-provided evalCons
        """
        # Check if user specified which constraints to output
        # Otherwise, output them all
        if evalCons is None:
            evalCons = self.constraintList
        else:
            userCons = sorted(list(evalCons))
            evalCons = {}
            for func in userCons:
                if func in self.constraintList:
                    evalCons[func] = self.constraintList[func]

        if not ignoreMissing:
            for f in evalCons:
                if f not in self.constraintList:
                    raise self._TACSError(
                        f"Supplied constraint '{f}' has not been added "
                        "using addConstraint()."
                    )

        return evalCons


# Simple class for handling sparse linear constraints in parallel
class SparseLinearConstraint(object):
    dtype = TACSConstraint.dtype

    def __init__(self, comm, rows, cols, vals, nrows, ncols, lb=-1e20, ub=1e20):
        # Sparse Jacobian for constraint
        self.A = sp.sparse.csr_matrix(
            (vals, (rows, cols)), shape=(nrows, ncols), dtype=self.dtype
        )
        # Number of constraints
        self.nCon = nrows
        # MPI comm
        self.comm = comm
        # Save bound information
        if isinstance(lb, np.ndarray) and len(lb) == self.nCon:
            self.lb = lb.astype(self.dtype)
        elif isinstance(lb, float) or isinstance(lb, complex):
            self.lb = np.array([lb] * self.nCon, dtype=self.dtype)

        if isinstance(ub, np.ndarray) and len(ub) == self.nCon:
            self.ub = ub.astype(self.dtype)
        elif isinstance(ub, float) or isinstance(ub, complex):
            self.ub = np.array([ub] * self.nCon, dtype=self.dtype)

    def evalCon(self, x):
        conVals = self.comm.allreduce(self.A.dot(x))
        return conVals

    def evalConSens(self, x=None):
        return self.A.copy()

    def getBounds(self):
        return self.lb.copy(), self.ub.copy()
