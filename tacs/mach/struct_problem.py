"""
pyBase_problem
"""

# =============================================================================
# Imports
# =============================================================================
from functools import wraps
import copy

from baseclasses import StructProblem as BaseStructProblem

import numpy as np
from mpi4py import MPI
import pyNastran.bdf as pn

import tacs.TACS


# Define decorator functions for methods that must be called before initialize
def updateDVGeo(method):
    @wraps(method)
    def wrappedMethod(self, *args, **kwargs):
        if self.DVGeo is not None:
            if not self.DVGeo.pointSetUpToDate(self.ptSetName):
                coords = self.DVGeo.update(self.ptSetName, config=self.name)
                self.staticProblem.setNodes(coords.flatten())
                for constraint in self.constraints:
                    constraint.setNodes(coords.flatten())
        return method(self, *args, **kwargs)

    return wrappedMethod


class StructProblem(BaseStructProblem):
    """
    MACH StructProblem wrapper for pyTACS StaticProblem class
    """

    def __init__(
        self,
        staticProblem,
        FEAAssembler,
        DVGeo=None,
        loadFile=None,
    ):
        """
        Initialize the StructProblem.

        Parameters
        ----------
        staticProblem : tacs.problems.StaticProblem
            StaticProblem object used for modeling and solving static cases.

        FEAAssembler : tacs.pytacs.pyTACS
            pyTACS assembler object.

        DVGeo : pygeo.DVGeometry, optional
            Object responsible for manipulating the geometry design variables.
            If None, no geometric design variables will be used.

        loadFile : str, optional
            Filename of the (static) external load file. Should be
            generated from pyAerostructure.
        """

        self.staticProblem = staticProblem
        self.FEAAssembler = FEAAssembler
        self.DVGeo = None
        self.numGeoDV = 0
        self.ptSetName = None
        self.loadFile = loadFile
        self.constraints = []

        if self.staticProblem.assembler != self.FEAAssembler.assembler:
            raise RuntimeError(
                "Provided StaticProblem does not correspond to pyTACS assembler object."
            )

        self.comm = self.staticProblem.comm
        if DVGeo:
            self.setDVGeo(DVGeo)

        self.update = self.FEAAssembler.createVec(asBVec=True)
        self.oldUpdate = self.FEAAssembler.createVec(asBVec=True)
        self.temp0 = self.FEAAssembler.createVec(asBVec=True)
        self.temp1 = self.FEAAssembler.createVec(asBVec=True)
        self._Fext = self.FEAAssembler.createVec(asBVec=True)
        self._phi = self.FEAAssembler.createVec(asBVec=True)
        self._pLoad = self.FEAAssembler.createVec(asBVec=True)
        self._dSdu = self.FEAAssembler.createVec(asBVec=True)
        self._adjRHS = self.FEAAssembler.createVec(asBVec=True)
        self._dIdu = self.FEAAssembler.createVec(asBVec=True)
        self._matVecRHS = self.FEAAssembler.createVec(asBVec=True)
        self._matVecSolve = self.FEAAssembler.createVec(asBVec=True)
        self.adjoints = {
            funcName: self.FEAAssembler.createVec(asBVec=True)
            for funcName in self.evalFuncs
        }

        self.callCounter = 0
        self.doDamp = False

        if self.loadFile:
            self.readExternalForceFile(self.loadFile)

    @property
    def name(self):
        """
        Get the name of the structural problem.

        Returns
        -------
        str
            Name of the structural problem.
        """
        return self.staticProblem.name

    @property
    def Fext(self):
        """
        Get the external aerodynamic force vector.

        Returns
        -------
        numpy.ndarray
            External aerodynamic force vector as a numpy array.
        """
        return self._Fext.getArray()

    @Fext.setter
    def Fext(self, value):
        """
        Set the external aerodynamic force vector.

        Parameters
        ----------
        value : numpy.ndarray
            External aerodynamic force vector to set.
        """
        self._Fext[:] = value

    @property
    def phi(self):
        """
        Get the adjoint vector.

        Returns
        -------
        numpy.ndarray
            Adjoint vector as a numpy array.
        """
        return self._phi.getArray()

    @phi.setter
    def phi(self, value):
        """
        Set the adjoint vector.

        Parameters
        ----------
        value : numpy.ndarray
            Adjoint vector to set.
        """
        self._phi[:] = value

    @property
    def pLoad(self):
        """
        The contribution of the aerodynamic adjoint to the structural
        adjoint RHS. This term is given by dAdu^T*psi and is subtracted
        from the structural adjoint RHS.

        Returns
        -------
        numpy.ndarray
            Aerodynamic adjoint contribution to the structural adjoint RHS.
        """
        return self._pLoad.getArray()

    @pLoad.setter
    def pLoad(self, value):
        """
        Set the aerodynamic adjoint contribution to the structural
        adjoint RHS. This term is given by dAdu^T*psi and is subtracted
        from the structural adjoint RHS.

        Parameters
        ----------
        value : numpy.ndarray
            Aerodynamic adjoint contribution to the structural adjoint RHS.
        """
        self._pLoad[:] = value

    @property
    def dSdu(self):
        """
        Get the product of the derivative of structural residual with
        respect to state variables and the structural adjoint vector.

        Returns
        -------
        numpy.ndarray
            Product of the derivative of structural residual with
            respect to state variables and the structural adjoint
            vector.
        """
        return self._dSdu.getArray()

    @dSdu.setter
    def dSdu(self, value):
        """
        Set the product of the derivative of structural residual with
        respect to state variables and the structural adjoint vector.

        Parameters
        ----------
        value : numpy.ndarray
            Product of the derivative of structural residual with
            respect to state variables and the structural adjoint
            vector.
        """
        self._dSdu[:] = value

    @property
    def adjRHS(self):
        """
        Get the adjoint right-hand side vector.

        Returns
        -------
        numpy.ndarray
            Adjoint RHS vector as a numpy array.
        """
        return self._adjRHS.getArray()

    @adjRHS.setter
    def adjRHS(self, value):
        """
        Set the adjoint right-hand side vector.

        Parameters
        ----------
        value : numpy.ndarray or tacs.TACS.Vec
            Adjoint RHS vector to set.
        """
        self._adjRHS[:] = value

    @property
    def dIdu(self):
        """
        Get the derivative of objective function with respect to state variables.

        Returns
        -------
        numpy.ndarray
            Derivative vector as a numpy array with boundary conditions applied.
        """
        self.FEAAssembler.applyBCsToVec(self._dIdu)
        return self._dIdu.getArray()

    @dIdu.setter
    def dIdu(self, value):
        """
        Set the derivative of objective function with respect to state variables.

        Parameters
        ----------
        value : numpy.ndarray or tacs.TACS.Vec
            Derivative vector to set.
        """
        self._dIdu[:] = value
        self.FEAAssembler.applyBCsToVec(self._dIdu)

    @property
    def matVecRHS(self):
        """
        Get the matrix-vector product right-hand side vector.

        Returns
        -------
        numpy.ndarray
            Matrix-vector RHS vector as a numpy array.
        """
        return self._matVecRHS.getArray()

    @matVecRHS.setter
    def matVecRHS(self, value):
        """
        Set the matrix-vector product right-hand side vector.

        Parameters
        ----------
        value : numpy.ndarray or tacs.TACS.Vec
            Matrix-vector RHS vector to set.
        """
        self._matVecRHS[:] = value

    @property
    def matVecSolve(self):
        """
        Get the matrix-vector solve vector.

        Returns
        -------
        numpy.ndarray
            Matrix-vector solve vector as a numpy array.
        """
        return self._matVecSolve.getArray()

    @matVecSolve.setter
    def matVecSolve(self, value):
        """
        Set the matrix-vector solve vector.

        Parameters
        ----------
        value : numpy.ndarray or tacs.TACS.Vec
            Matrix-vector solve vector to set.
        """
        self._matVecSolve[:] = value

    def getVarName(self):
        """
        Get name for the design variables in pyOpt. Only needed
        if more than 1 pytacs object is used in an optimization

        Returns
        ----------
        varName : str
            Name of the design variables used in setDesignVars() dict.
        """
        return self.staticProblem.getVarName()

    def setDVGeo(self, DVGeo, pointSetKwargs=None):
        """
        Set the DVGeometry object that will manipulate 'geometry' in
        this object. Note that TACS doe not **strictly** need a
        DVGeometry object, but if optimization with geometric
        changes is desired, then it is required.

        Parameters
        ----------
        dvGeo : A DVGeometry object.
            Object responsible for manipulating the constraints that
            this object is responsible for.

        pointSetKwargs : dict
            Keyword arguments to be passed to the DVGeo addPointSet call.
            Useful for DVGeometryMulti, specifying FFD projection tolerances, etc

        Examples
        --------
        >>> StructProblem.setDVGeo(DVGeo)
        """
        if self.DVGeo:
            raise RuntimeError("DVGeo has already been set.")

        if pointSetKwargs is None:
            pointSetKwargs = {}

        self.DVGeo = DVGeo
        # Get the number of geometry variables
        self.numGeoDV = self.DVGeo.getNDV()

        self.ptSetName = "tacs_%s_coords" % self.name
        coords0 = self.staticProblem.getNodes()
        self.DVGeo.addPointSet(coords0.reshape(-1, 3), self.ptSetName, **pointSetKwargs)

    def getDVGeo(self):
        """
        Get the DVGeometry object.

        Returns
        -------
        DVGeometry: pygeo.parameterization.BaseDVGeometry or None
            DVGeometry object.
        """
        return self.DVGeo

    def setDesignVars(self, x):
        """
        Set the variables in the x-dict for this object.

        Parameters
        ----------
        x : dict
            Dictionary of variables which may or may not contain the
            design variable names this object needs
        """
        if self.comm.rank != 0:
            x = np.empty(0)
        self.staticProblem.setDesignVars(x)

        for constr in self.constraints:
            constr.setDesignVars(x)

    def addConstraint(self, constr):
        """
        Add pyTACS constraint object to StructProblem.

        Parameters
        ----------
        constr : tacs.constraints.TACSConstraint
            pyTACS Constraint object
        """
        if self.staticProblem.assembler != constr.assembler:
            raise ValueError(
                f"TACSConstraint object '{constr.name}' and StaticProblem '{self.staticProblem.name}' "
                "were not created by same pyTACS assembler"
            )

        self.constraints.append(constr)

    def addVariablesPyOpt(self, optProb):
        """
        Add the current set of variables to the optProb object.

        Parameters
        ----------
        optProb : pyOpt_optimization class
            Optimization problem definition to which variables are added
        """
        ndv = self.FEAAssembler.getTotalNumDesignVars()
        dvName = self.staticProblem.getVarName()
        if ndv > 0 and dvName not in optProb.variables:
            value = self.getOrigDesignVars()
            lb, ub = self.getDesignVarRange()
            scale = self.getDesignVarScales()

            optProb.addVarGroup(
                dvName, ndv, "c", value=value, lower=lb, upper=ub, scale=scale
            )

    def addConstraintsPyOpt(self, optProb):
        """
        Add any linear constraints that were generated during setup to
        the specified pyOpt problem.

        Parameters
        ----------
        optProb : :class:`Optimization <pyoptsparse.pyOpt_optimization.Optimization>` instance
            Optimization problem object to add constraints to
        """
        fcon = {}
        fconSens = {}
        conSizes = {}
        conBounds = {}
        self.evalConstraints(fcon)
        self.evalConstraintsSens(fconSens)
        for constr in self.constraints:
            constr.getConstraintSizes(conSizes)
            constr.getConstraintBounds(conBounds)

        for conName in fcon:
            nCon = conSizes[conName]
            if nCon > 0:
                # Save the nonlinear constraint name
                lb, ub = conBounds[conName]

                # Just evaluate the constraint to get the jacobian structure
                wrt = list(fconSens[conName].keys())
                optProb.addConGroup(
                    conName, nCon, wrt=wrt, jac=fconSens[conName], lower=lb, upper=ub
                )

    @updateDVGeo
    def getNodes(self):
        """
        Return the mesh coordinates of this problem.

        Returns
        -------
        coords : numpy.ndarray
            Structural coordinate in array of size (N * 3) where N is
            the number of structural nodes on this processor.
        """
        return self.staticProblem.getNodes()

    @updateDVGeo
    def solve(self, damp=1.0, useAitkenAcceleration=False, dampLB=0.2, loadScale=1.0):
        """
        Solution of the static problem for current load set. The
        forces must already be set.

        Parameters
        ----------
        Fext : numpy.ndarray or tacs.TACS.Vec, optional
            Distributed array containing additional loads (ex. aerodynamic forces for aerostructural coupling)
            to applied to RHS of the static problem, by default None
        damp : float
            Value to use to damp the solution update.
        useAitkenAcceleration : bool
            Flag to use aitkenAcceleration. Only applicable for aerostructural problems
        dampLB : float
            Lower bound for the damping parameter used in Aitken acceleration.
        loadScale : float
            Value to scale external loads by. Only useful for load step approach on nonlinear problems.
        """
        # Save old load scale and update with user-defined one
        loadScale0 = self.staticProblem.getLoadScale()
        self.staticProblem.setLoadScale(loadScale)

        # Save initial state vector
        u0 = self.temp1
        self.staticProblem.getVariables(u0)

        # Solve static problem w/o damping
        successFlag = self.staticProblem.solve(Fext=self._Fext)

        # Compute undamped update
        self.staticProblem.getVariables(self.update)
        self.update.axpy(-1.0, u0)

        # Reset problem back to unsolved state
        self.staticProblem.setVariables(u0)

        # Apply Aitken Acceleration if necessary:
        if useAitkenAcceleration:
            if self.doDamp:
                # Compute: temp0 = update - old_update
                self.temp0.zeroEntries()
                self.temp0.axpy(1.0, self.update)
                self.temp0.axpy(-1.0, self.oldUpdate)
                dnom = self.temp0.dot(self.temp0)
                damp = damp * (1.0 - self.temp0.dot(self.update) / dnom)

                # Clip to a reasonable range
                damp = np.clip(damp, dampLB, 1.0)
            self.doDamp = True

        # Update State Variables
        u1 = self.temp0
        u1.copyValues(u0)
        u1.axpy(damp, self.update)
        self.staticProblem.setVariables(u1)

        # Set the old update
        self.oldUpdate.copyValues(self.update)

        # Set load scale back
        self.staticProblem.setLoadScale(loadScale0)

        return successFlag

    @updateDVGeo
    def evalFunctions(self, funcs, evalFuncs=None, ignoreMissing=False):
        """
        Evaluate the current functions

        Parameters
        ----------
        funcs : dict
            Dictionary into which the functions are save
        evalFuncs : iterable object containing strings
            The functions that the user wants evaluated
        """
        if evalFuncs is None:
            evalFuncs = self.evalFuncs

        return self.staticProblem.evalFunctions(funcs, evalFuncs, ignoreMissing)

    @updateDVGeo
    def evalConstraints(self, fcon, evalCons=None, ignoreMissing=False):
        """
        Evaluate values for constraints. The constraints corresponding to the strings in
        evalCons are evaluated and updated into the provided
        dictionary.

        Parameters
        ----------
        fcon : dict
            Dictionary into which the constraints are saved.
        evalCons : iterable object containing strings.
            If not none, use these constraints to evaluate.
        ignoreMissing : bool
            Flag to supress checking for a valid constraint. Please use
            this option with caution.
        """
        for constr in self.constraints:
            if evalCons is None or constr.name in evalCons:
                constr.evalConstraints(fcon, evalCons, ignoreMissing)

    @updateDVGeo
    def evalFunctionsSens(self, funcsSens, evalFuncs=None, includeXptSens=False):
        """
        Evaluate the sensitivity of the desired functions

        Parameters
        ----------
        funcsSens : dict
            Dictionary into which the function sensitivities are saved
        evalFuncs : iterable object containing strings
            The functions that the user wants evaluated
        includeXptSens : bool
            Flag to include sensitivities with respect to the node locations.
        """
        if evalFuncs is None:
            evalFuncs = self.evalFuncs

        self.staticProblem.evalFunctionsSens(funcsSens, evalFuncs)

        for funcName in evalFuncs:
            funcKey = f"{self.name}_{funcName}"
            dvKey = self.staticProblem.getVarName()
            if dvKey in funcsSens[funcKey]:
                funcsSens[funcKey][dvKey] = self.comm.bcast(
                    funcsSens[funcKey][dvKey], root=0
                )

        # Compute the DVGeo sensitivities if requested
        coordName = self.staticProblem.getCoordName()
        if self.DVGeo is not None:
            for funcName in evalFuncs:
                funcKey = f"{self.name}_{funcName}"
                if coordName in funcsSens[funcKey]:
                    dIdpt = funcsSens[funcKey].pop(coordName).reshape(-1, 3)
                    dIdx = self.DVGeo.totalSensitivity(
                        dIdpt, self.ptSetName, comm=self.comm, config=self.name
                    )
                    funcsSens[funcKey].update(dIdx)

        # Pop out the node sensitivities if requested
        elif not includeXptSens:
            for funcName in evalFuncs:
                funcKey = f"{self.name}_{funcName}"
                if coordName in funcsSens[funcKey]:
                    funcsSens[funcKey].pop(coordName).reshape(-1, 3)

    @updateDVGeo
    def evalConstraintsSens(self, fconSens, evalCons=None, includeXptSens=False):
        """
        This is the main routine for returning useful (sensitivity)
        information from constraint. The derivatives of the constraints
        corresponding to the strings in evalCons are evaluated and
        updated into the provided dictionary. The derivatives with
        respect to all design variables and node locations are computed.

        Parameters
        ----------
        fconSens : dict
            Dictionary into which the derivatives are saved.
        evalCons : iterable object containing strings
            The constraints the user wants returned
        includeXptSens : bool
            Flag to include sensitivities with respect to the node locations.
        """
        sens = {}
        for constr in self.constraints:
            if evalCons is None or constr.name in evalCons:
                constr.evalConstraintsSens(sens, evalCons)

        for conKey in sens:
            dvKey = self.staticProblem.getVarName()
            if dvKey in sens[conKey]:
                sens[conKey][dvKey] = self.comm.bcast(sens[conKey][dvKey], root=0)

        # Compute the DVGeo sensitivities if requested
        coordName = self.staticProblem.getCoordName()
        if self.DVGeo is not None:
            self.DVGeo.computeTotalJacobian(self.ptSetName, config=self.name)
            for conKey in sens:
                if coordName in sens[conKey]:
                    # Pop out the constraint sensitivities wrt TACS coords
                    dIdpt = sens[conKey].pop(coordName)
                    # Check if the constraint sensitivities wrt TACS coords are zero on all procs
                    total_nnz = self.comm.allreduce(dIdpt.nnz, op=MPI.SUM)
                    if total_nnz == 0:
                        # if so, skip DVGeo sensitivities
                        continue
                    # Get the Jacobian
                    Jacobian = self.DVGeo.JT[self.ptSetName]
                    # Compute the local Jacobian product
                    dIdx_local = dIdpt.dot(Jacobian.T)
                    # Add dvgeo contribution across all procs
                    dIdx = self.comm.allreduce(dIdx_local.toarray(), op=MPI.SUM)
                    # Convert to dict
                    dIdx_dict = self.DVGeo.convertSensitivityToDict(np.atleast_2d(dIdx))
                    # Update sensitivity dict with new DVGeo sensitivities
                    sens[conKey].update(dIdx_dict)

        # Pop out the node sensitivities if requested
        elif not includeXptSens:
            for conKey in sens:
                if coordName in sens[conKey]:
                    sens[conKey].pop(coordName).reshape(-1, 3)

        fconSens.update(sens)

    def getVarsPerNode(self):
        """
        The number of degrees of freedom used at each output location.

        Returns
        -------
        int
            Number of degrees of freedom of each node in the domain.
        """
        return self.staticProblem.getVarsPerNode()

    def getNumNodes(self):
        """
        Get the number of nodes on this processor.

        Returns
        -------
        int
            Number of nodes on this processor.
        """
        return self.staticProblem.getNumOwnedNodes()

    def getOrigDesignVars(self):
        """
        Get an array holding the original values for all design
        variables added to TACS at time of initialization.

        Returns
        -------
        numpy.ndarray
            Array containing all design variable values across all processors.
        """
        local_dvs = self.FEAAssembler.getOrigDesignVars()
        all_local_dvs = self.comm.allgather(local_dvs)
        global_dvs = np.concatenate(all_local_dvs)
        return global_dvs.astype(float)

    def getDesignVarRange(self):
        """
        Get arrays containing the lower and upper bounds for the design variables,
        in the form needed by pyoptsparse.

        Returns
        -------
        tuple of numpy.ndarray
            Lower and upper bounds for the design variables.
        """
        local_lb, local_ub = self.FEAAssembler.getDesignVarRange()
        all_lb = self.comm.allgather(local_lb)
        global_lbs = np.concatenate(all_lb)
        all_ub = self.comm.allgather(local_ub)
        global_ubs = np.concatenate(all_ub)
        return global_lbs.astype(float), global_ubs.astype(float)

    def getDesignVarScales(self):
        """
        Get an array containing the scaling factors for the design
        variables, in the form needed by pyoptsparse.
        method.

        Returns
        -------
        numpy.ndarray
            Scaling values for design variables.
        """
        return np.array(self.FEAAssembler.scaleList)

    def getNumDesignVars(self):
        """
        Get total number of structural design variables across all processors.

        Returns
        -------
        int
            Total number of structural design variables.
        """
        return self.FEAAssembler.getTotalNumDesignVars()

    def getStates(self, states=None):
        """
        Return the current state values for the problem.

        Parameters
        ----------
        states : tacs.TACS.Vec or numpy.ndarray, optional
            Vector to place current state variables into.

        Returns
        -------
        numpy.ndarray
            Current state vector.
        """
        return self.staticProblem.getVariables(states)

    def setStates(self, states=None):
        """
        Set the current state values for the problem.

        Parameters
        ----------
        states : tacs.TACS.Vec or numpy.ndarray
            Vector to replace current state variables with.
        """
        self.staticProblem.setVariables(states)

    def getExternalForce(self, Fext=None):
        """
        Return the current external coupling force vector for the problem.

        Parameters
        ----------
        Fext : tacs.TACS.Vec or numpy.ndarray, optional
            Vector to place current external coupling force vector into.

        Returns
        -------
        numpy.ndarray
            Current external coupling force vector.
        """
        if isinstance(Fext, tacs.TACS.Vec):
            Fext.copyValues(self._Fext)
        elif isinstance(Fext, np.ndarray):
            Fext[:] = self._Fext[:]

        return self._Fext.getArray().copy()

    def setExternalForce(self, Fext):
        """
        Set the external coupling force vector for the problem.

        Parameters
        ----------
        Fext : tacs.TACS.Vec or numpy.ndarray
            External coupling force vector to set.
        """
        if isinstance(Fext, tacs.TACS.Vec):
            self._Fext.copyValues(Fext)
        elif isinstance(Fext, np.ndarray):
            self._Fext.getArray()[:] = Fext[:]

    def getStateNorm(self):
        """
        Get the norm of the current state vector.

        Returns
        -------
        float
            Norm of the current state vector.
        """
        u = self.temp0
        self.staticProblem.getVariables(u)
        return u.norm()

    def getResidual(self):
        """
        This routine is used to evaluate directly the structural
        residual. Only typically used with aerostructural analysis.

        Returns
        -------
        numpy.ndarray
            Computed residuals for problem.
        """
        res = self.temp0
        self.staticProblem.getResidual(res, Fext=self._Fext)
        resArray = res.getArray().copy()
        self.temp0.zeroEntries()
        return resArray

    def getResNorms(self):
        """
        Return the initial, starting and final residual norms. Note that
        the same norms are used for both solution and adjoint
        computations.

        Returns
        -------
        tuple of float
            Initial, starting, and final residual norms.
        """
        initNorm = np.real(self.staticProblem.initNorm)
        startNorm = np.real(self.staticProblem.startNorm)
        finalNorm = np.real(self.staticProblem.finalNorm)

        return (initNorm, startNorm, finalNorm)

    def setResNorms(self, initNorm=None, startNorm=None, finalNorm=None):
        """
        Set one of these norms if not None.
        This is typically only used with aerostructural analysis.

        Parameters
        ----------
        initNorm : float, optional
            Initial residual norm to set.
        startNorm : float, optional
            Starting residual norm to set.
        finalNorm : float, optional
            Final residual norm to set.
        """
        if initNorm is not None:
            self.staticProblem.initNorm = initNorm
        if startNorm is not None:
            self.staticProblem.startNorm = startNorm
        if finalNorm is not None:
            self.staticProblem.finalNorm = finalNorm

    def zeroVectors(self):
        """
        Zero all the TACS b-vectors.
        """
        self.staticProblem.zeroVariables()
        self.update.zeroEntries()

    @property
    def evalFuncs(self):
        """
        Get strings for added evaluation functions.

        Returns
        -------
        list of str
            List of function names that can be evaluated.
        """
        return list(self.staticProblem.getFunctionKeys())

    def getAdjoint(self, objective):
        """
        Return the adjoint values for objective if they exist. Otherwise just return zeros.

        Parameters
        ----------
        objective : str
            Name of the objective function.

        Returns
        -------
        numpy.ndarray
            Adjoint values for the specified objective.
        """
        if objective in self.adjoints:
            vals = self.adjoints[objective].getArray().copy()
        else:
            vals = self.FEAAssembler.createVec()
        return vals

    def setAdjoint(self, adjoint, objective=None):
        """
        Set externally supplied adjoint values.

        Parameters
        ----------
        adjoint : float, numpy array
            An array of size getStateSize() to
            be copied to the structural adjoint variables
        objective : str
            Name of objective to set adjoint for
        """
        self.phi[:] = adjoint
        if objective is not None:
            if objective not in self.adjoints:
                self.adjoints[objective] = self.FEAAssembler.createVec(asBVec=True)
            self.adjoints[objective].getArray()[:] = adjoint[:]

    def setAdjointRHS(self, func):
        """
        Set the adjoint right-hand side for given function.

        Parameters
        ----------
        func : str
            Name of the function to set adjoint RHS for.
        """
        self._dIdu.zeroEntries()
        if func in self.staticProblem.functionList:
            self.staticProblem.addSVSens([func], [self._dIdu])

    def getdSduVec(self, inVec):
        """
        Evaluate the direct residual.

        Parameters
        ----------
        inVec : numpy.ndarray
            Input vector.

        Returns
        -------
        numpy.ndarray
            Direct residual vector.
        """
        self.phi[:] = inVec
        self.staticProblem.K.mult(self._phi, self.temp0)
        outVec = self.temp0.getArray().copy()
        self.temp0.zeroEntries()

        return outVec

    def getdSduTVec(self, inVec):
        """
        Evaluate the adjoint residual.

        Parameters
        ----------
        inVec : numpy.ndarray
            Input vector.

        Returns
        -------
        numpy.ndarray
            Adjoint residual vector.
        """
        self.phi[:] = inVec
        self.temp0.zeroEntries()
        self.staticProblem.addTransposeJacVecProduct(self._phi, self.temp0, scale=1.0)
        self.FEAAssembler.applyBCsToVec(self.temp0)
        outVec = self.temp0.getArray().copy()
        self.temp0.zeroEntries()

        return outVec

    def getdIdXpt(self, func):
        """
        Get the (partial) derivative of the structural objective
        with respect to all structural nodes for the current
        structProblem.

        Parameters
        ----------
        func : str
            Name of the function to compute derivatives for.

        Returns
        -------
        tacs.TACS.Vec
            Derivative vector with respect to structural nodes.
        """
        dIdXpt = self.FEAAssembler.createNodeVec()
        # Check to see if func is actually a valid TACS function
        if func in self.evalFuncs:
            self.staticProblem.addXptSens([func], [dIdXpt], scale=1.0)

        return dIdXpt

    def getdIdXdv(self, func):
        """
        Get the (partial) derivative of the structural objective
        with respect to all structural design variables.

        Parameters
        ----------
        func : str
            Name of the function to compute derivatives for.

        Returns
        -------
        tacs.TACS.Vec
            Derivative vector with respect to structural design variables.
        """
        dIdXdv = self.FEAAssembler.createDesignVec()

        # Check to see if func is actually a valid TACS function
        if func in self.evalFuncs:
            self.staticProblem.addDVSens([func], [dIdXdv], scale=1.0)

        dvVal = self.comm.bcast(dIdXdv, root=0)

        return dvVal

    def getdIduTransProd(self, vecT, evalFuncs=None):
        """
        Perform the transpose matvec of the getdIduProd function and
        return a TACS vector. The result is stored locally.

        We assume that structProblem is set and that evalFuncs is valid.

        Parameters
        ----------
        vecT : dict
            Dictionary containing transpose vector values.
        evalFuncs : list of str, optional
            List of function names to evaluate.

        Returns
        -------
        numpy.ndarray
            Transpose matrix-vector product result.
        """
        # The old way (for testing purposes)
        self._matVecRHS.zeroEntries()
        svList = [self.FEAAssembler.createVec() for f in evalFuncs]
        self.staticProblem.addSVSens(evalFuncs, svList)
        for i, f in enumerate(evalFuncs):
            f_mangled = self.name + "_%s" % f
            self._matVecRHS.axpy(vecT[f_mangled][0], svList[i])

        return self.matVecRHS.copy()

    def getdIdXdvTransProd(self, vecT, evalFuncs=None):
        """
        Perform the transpose matvec of the getdIdXdvProd function and
        return a vector.

        We assume that structProblem is set and that handles is valid.

        Note that handles contains the list of function handles, whereas the
        names are contained in evalFuncs.

        Parameters
        ----------
        vecT : dict
            Dictionary containing transpose vector values.
        evalFuncs : list of str, optional
            List of function names to evaluate.

        Returns
        -------
        dict
            Dictionary containing transpose matrix-vector product results.
        """
        # # The old way (for testing purposes)
        prodDV = self.FEAAssembler.createDesignVec(asBVec=True)
        for f in evalFuncs:
            f_mangled = self.name + "_%s" % f
            self.staticProblem.addDVSens(
                [evalFuncs], [prodDV], scale=vecT[f_mangled][0]
            )

        # Convert result back into a dictionary
        prodDict = self.convertDesignVecToDict(prodDV.getArray())

        if self.DVGeo is not None:
            prodXpt = self.FEAAssembler.createNodeVec(asBVec=True)
            for f in evalFuncs:
                f_mangled = self.name + "_%s" % f
                self.staticProblem.addXptSens(
                    [evalFuncs], [prodXpt], scale=vecT[f_mangled][0]
                )
            xArray = prodXpt.getArray()
            xdot = self.DVGeo.totalSensitivity(
                xArray.reshape(-1, 3), self.ptSetName, comm=self.comm, config=self.name
            )
            prodDict.update(xdot)

        return prodDict

    def globalNKPreCon(self, inVec):
        """
        This function is ONLY used as a preconditioner to the
        global Aero-Structural system.

        Parameters
        ----------
        inVec : numpy.ndarray
            Input vector for preconditioning.

        Returns
        -------
        numpy.ndarray
            Preconditioned output vector.
        """
        # Place the in_vec into the residual vector
        res = self.temp0
        update = self.temp1
        res.getArray()[:] = inVec

        # Solve
        self.staticProblem.linearSolver.solve(res, update)
        outVec = update.getArray()

        # Restore Pointers
        res.zeroEntries()
        update.zeroEntries()

        return outVec

    def globalAdjointPreCon(self, inVec):
        """
        This function is ONLY used as a preconditioner to the
        global Aero-Structural adjoint system. For the linear
        symmetric case this is actually identical to the NKPreCon
        function above.

        Parameters
        ----------
        inVec : numpy.ndarray
            Input vector for preconditioning.

        Returns
        -------
        numpy.ndarray
            Preconditioned output vector.
        """
        # Place the in_vec into the residual vector
        self.temp0.getArray()[:] = inVec[:]
        self.FEAAssembler.applyBCsToVec(self.temp0)

        # Apply the preconditioner ONLY.
        self.staticProblem.PC.applyFactor(self.temp0, self.temp1)
        self.FEAAssembler.applyBCsToVec(self.temp1)
        outVec = self.temp1.getArray().copy()

        # Zero values
        self.temp0.zeroEntries()
        self.temp1.zeroEntries()

        return outVec

    def globalDirectPreCon(self, inVec):
        """
        This function is ONLY used as a preconditioner to the
        global Aero-Structural direct system. For the linear
        symmetric case this is actually identical to the NKPreCon
        function above.

        Parameters
        ----------
        inVec : numpy.ndarray
            Input vector for preconditioning.

        Returns
        -------
        numpy.ndarray
            Preconditioned output vector.
        """
        return self.globalAdjointPreCon(inVec)

    def assembleAdjointRHS(self):
        """
        Compute the final adjoint RHS:
        Adjoint RHS should be: dIdu - dAdu^T*psi - dSdu^T*phi
        """
        # Set RHS to dIdu
        self._adjRHS.copyValues(self._dIdu)

        # Subtract dSdu if non-zero
        self._adjRHS.axpy(-1.0, self._dSdu)

        # Subtract pLoad
        self._adjRHS.axpy(-1.0, self._pLoad)

        self.staticProblem.initNorm = np.real(self._adjRHS.norm())

        self.oldUpdate.zeroEntries()

    def solveAdjoint(self, damp=1.0):
        """
        Solve the structural adjoint.

        Parameters
        ----------
        damp : float
            A damping variable for adjoint update. Typically only used
            in multidisciplinary analysis
        """
        res = self.temp0
        res.zeroEntries()

        self.FEAAssembler.applyBCsToVec(self._adjRHS)

        # First compute the residual
        self.staticProblem.addTransposeJacVecProduct(self._phi, res)
        res.axpy(-1.0, self._adjRHS)  # Add the -RHS

        # Starting Norm for this computation
        self.staticProblem.startNorm = np.real(res.norm())

        # Solve Linear System
        self.update.zeroEntries()
        self.staticProblem.solveAdjoint(res, self.update)

        self.FEAAssembler.applyBCsToVec(self.update)

        # Update the adjoint vector with the (damped) update
        self._phi.axpy(-damp, self.update)

        # Compute actual final FEA Norm
        res.zeroEntries()
        self.staticProblem.addTransposeJacVecProduct(self._phi, res)
        res.axpy(-1.0, self._adjRHS)  # Add the RHS
        self.staticProblem.finalNorm = np.real(res.norm())

    def getdRdXptPhi(self, objectives):
        """
        Get the result of :math:`[dR/dX_{nodes}]^T \phi` for each of the objectives in the objective list.

        This is the total sensitivity calculation.

        Parameters
        ----------
        objectives : list of str
            List of objective function names.

        Returns
        -------
        list of tacs.TACS.Vec
            List of sensitivity products for each objective.
        """
        products = [self.FEAAssembler.createNodeVec() for obj in objectives]
        phiList = [self.getAdjoint(obj) for obj in objectives]
        self.staticProblem.addAdjointResXptSensProducts(phiList, products, scale=1.0)

        return products

    def getdRdXdvPhi(self, objectives):
        """
        Get the result of :math:`[dR/dX_{dv}]^T \phi` for the total sensitivity calculation.

        Parameters
        ----------
        objectives : list of str
            List of objective function names.

        Returns
        -------
        list of tacs.TACS.Vec
            List of sensitivity products for each objective.
        """
        products = [self.FEAAssembler.createDesignVec() for obj in objectives]
        phiList = [self.getAdjoint(obj) for obj in objectives]
        self.staticProblem.addAdjointResProducts(phiList, products, scale=1.0)
        dvVals = [self.comm.bcast(dvVec, root=0) for dvVec in products]

        return dvVals

    def getdRdXdvTransProd(self):
        """
        Perform the transpose matvec of getdRdXdvProd. The result is
        returned as a dict of numpy arrays.

        Returns
        -------
        dict
            Dictionary containing transpose matrix-vector product results.
        """
        dvProd = self.FEAAssembler.createDesignVec(asBVec=True)
        self.staticProblem.addAdjointResProducts(
            [self._matVecSolve], [dvProd], scale=1.0
        )
        # Convert result back into a dictionary
        prodDict = self.convertDesignVecToDict(dvProd.getArray())

        if self.DVGeo is not None:
            xptProd = self.FEAAssembler.createNodeVec(asBVec=True)
            self.staticProblem.addAdjointResXptSensProducts(
                [self._matVecSolve], [xptProd], scale=1.0
            )
            xArray = xptProd.getArray()
            xdot = self.DVGeo.totalSensitivity(
                xArray.reshape(-1, 3), self.ptSetName, comm=self.comm, config=self.name
            )
            prodDict.update(xdot)

        return prodDict

    def convertDesignVecToDict(self, dvVec):
        """
        Convert a design vector to a dictionary format.

        Parameters
        ----------
        dvVec : tacs.TACS.Vec or numpy.ndarray
            Design vector to convert.

        Returns
        -------
        dict
            Dictionary containing the design vector with variable name as key.
        """
        if isinstance(dvVec, tacs.TACS.Vec):
            dvVec = dvVec.getArray()

        dvVals = self.comm.bcast(dvVec, root=0)
        dvDict = {self.staticProblem.getVarName(): dvVals}

        return dvDict

    def writeSolution(self, outputDir=None, baseName=None, number=None):
        """
        This is a generic shell function that writes the output
        file(s).  The intent is that the user or calling program can
        call this function and pyTACS writes all the files that the
        user has defined. It is recommended that this function is used
        along with the associated logical flags in the options to
        determine the desired writing procedure

        Parameters
        ----------
        outputDir : str or None
            Use the supplied output directory
        baseName : str or None
            Use this supplied string for the base filename. Typically
            only used from an external solver.
        number : int or None
            Use the user supplied number to index solution. Again, only
            typically used from an external solver
        """
        self.staticProblem.writeSolution(outputDir, baseName, number)

    def writeExternalForceFile(self, outputDir=None, baseName=None, number=None):
        """
        This function writes the external loads to a file.
        This is typically used to save loads from an aerostructural run.

        Parameters
        ----------
        outputDir : str, optional
            Output directory for the force file.
        baseName : str, optional
            Base name for the force file. If None, uses problem name + "_external_forces".
        number : int, optional
            Number to append to the filename.
        """
        if baseName is None:
            baseName = self.name + "_external_forces"
        # Figure out the output file base name
        fileName = (
            self.staticProblem.getOutputFileName(outputDir, baseName, number) + ".dat"
        )
        # We want to isolate only the external loads in the rhs before writing the loads out
        rhs = self.staticProblem.rhs
        # Save a copy of the rhs vector holding the full loads
        self.temp0.copyValues(rhs)
        # Replace vector with only external loads
        rhs.copyValues(self._Fext)
        # Write external loads to bdf
        self.staticProblem.writeLoadToBDF(fileName, loadCaseID=0)
        # Reset rhs back to full loads
        rhs.copyValues(self.temp0)

    def readExternalForceFile(self, fileName):
        """
        Reads in a force file and sets the external forces in the structural problem.
        This is typically used to read in loads saved from an aerostructural run.

        Parameters
        ----------
        fileName : str
            Filename for force file. Should have .dat or .bdf extension.
        """
        forceInfo = pn.bdf.read_bdf(
            fileName, validate=False, xref=False, debug=False, punch=True
        )
        bdfInfo = self.staticProblem.bdfInfo

        # Step 1: Store original loads from bdfInfo
        originalLoads = copy.deepcopy(bdfInfo.loads)
        originalLoadCombinations = copy.deepcopy(bdfInfo.load_combinations)

        # Step 2: Overwrite bdfInfo loads with forceInfo loads
        bdfInfo.loads = copy.deepcopy(forceInfo.loads)
        bdfInfo.load_combinations = copy.deepcopy(forceInfo.load_combinations)

        # Create a copy of the internal loads already added to model
        F = self.staticProblem.F
        self.temp0.copyValues(F)
        # Zero out the loads
        F.zeroEntries()
        # Read the external loads
        self.staticProblem.addLoadFromBDF(0)
        # Copy new loads to ext load vector
        self._Fext.copyValues(F)
        # Set internal loads back to previous values
        F.copyValues(self.temp0)

        # Step 3: Restore original loads back into bdfInfo
        bdfInfo.loads = originalLoads
        bdfInfo.load_combinations = originalLoadCombinations
