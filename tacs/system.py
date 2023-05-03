"""
pySystem
"""

# =============================================================================
# Imports
# =============================================================================
import numpy as np

import tacs.TACS
from .utilities import BaseUI


class TACSSystem(BaseUI):
    """
    Base class for TACS problem/constraint types. Contains methods common to all TACS systems dealing with design variables.
    """

    def __init__(
        self, assembler, comm=None, options=None, outputViewer=None, meshLoader=None
    ):
        # TACS assembler object
        self.assembler = assembler
        # TACS F5 output writer
        self.outputViewer = outputViewer
        # TACS pyMeshLoader object
        self.meshLoader = meshLoader
        # pyNastran BDF object
        if self.meshLoader:
            self.bdfInfo = self.meshLoader.getBDFInfo()

        # Create Design variable vector
        self.x = self.assembler.createDesignVec()
        self.assembler.getDesignVars(self.x)
        self.varName = "struct"
        # Create Nodal coordinate vector
        self.Xpts = self.assembler.createNodeVec()
        self.assembler.getNodes(self.Xpts)
        self.coordName = "Xpts"

        # Setup comm and options
        BaseUI.__init__(self, options=options, comm=comm)

        return

    ####### Design variable methods ########

    def setVarName(self, varName):
        """
        Set a name for the structural variables in pyOpt. Only needs
        to be changed if more than 1 pytacs object is used in an
        optimization

        Parameters
        ----------
        varName : str
            Name of the structural variable used in addVarGroup().
        """
        self.varName = varName

    def getDesignVars(self):
        """
        Get the current set of  design variables for this problem.

        Returns
        ----------
        x : numpy.ndarray
            The current design variable vector set in tacs.

        """
        return self.x.getArray().copy()

    def setDesignVars(self, x):
        """
        Update the design variables used by tacs.

        Parameters
        ----------
        x : numpy.ndarray or dict or TACS.Vec
            The variables (typically from the optimizer) to set. It
            looks for variable in the ``self.varName`` attribute if in dict.

        """
        # Check if the design variables are being handed in a dict
        if isinstance(x, dict):
            if self.varName in x:
                self.x.getArray()[:] = x[self.varName]
        # or array
        elif isinstance(x, np.ndarray):
            self.x.getArray()[:] = x
        # Or TACS BVec
        elif isinstance(x, tacs.TACS.Vec):
            self.x.copyValues(x)
        else:
            raise ValueError(
                "setDesignVars must be called with either a numpy array, dict, or TACS Vec as input."
            )

        # Set the variables in tacs
        self.assembler.setDesignVars(self.x)

    def getDesignVarRange(self):
        """
        get the lower/upper bounds for the design variables.

        Returns
        ----------
        xlb : numpy.ndarray
            The design variable lower bound.
        xub : numpy.ndarray
            The design variable upper bound.

        """
        xlb = self.assembler.createDesignVec()
        xub = self.assembler.createDesignVec()
        self.assembler.getDesignVarRange(xlb, xub)
        return xlb.getArray(), xub.getArray()

    def _arrayToDesignVec(self, dvArray):
        """
        Converts a distributed numpy array into a TACS design variable BVec.

        Parameters
        ----------
        dvArray : numpy.ndarray
                  Numpy array for which to convert to TACS designVec.

        Returns
        -------
        xVec : TACS.Vec
               Converted TACS designVec.

        Notes
        -----
        dvArray must have correct size on each processor.
        """
        xVec = self.assembler.createDesignVec()

        # Set values
        xVec.getArray()[:] = dvArray

        # Return as tacs bvec object
        return xVec

    def getNumDesignVars(self):
        """
        Return the number of design variables on this processor.
        """
        return self.x.getSize()

    def getNodes(self):
        """
        Return the mesh coordinates of this problem.

        Returns
        -------
        coords : array
            Structural coordinate in array of size (N * 3) where N is
            the number of structural nodes on this processor.
        """
        return self.Xpts.getArray().copy()

    def setNodes(self, Xpts):
        """
        Set the mesh coordinates of the structure.

        Parameters
        ----------
        coords : numpy.ndarray
            Structural coordinate in array of size (N * 3) where N is
            the number of structural nodes on this processor.
        """
        # Check if the design variables are being handed in a dict
        if isinstance(Xpts, dict):
            if self.coordName in Xpts:
                self.Xpts.getArray()[:] = Xpts[self.coordName]
        # or array
        elif isinstance(Xpts, np.ndarray):
            self.Xpts.getArray()[:] = Xpts
        # Or TACS BVec
        elif isinstance(Xpts, tacs.TACS.Vec):
            self.Xpts.copyValues(Xpts)
        else:
            raise ValueError(
                "setNodes must be called with either a numpy array, dict, or TACS Vec as input."
            )
        self.assembler.setNodes(self.Xpts)

    def _arrayToNodeVec(self, xptsArray):
        """
        Converts a distributed numpy array into a TACS node BVec.

        Parameters
        ----------
        xptsArray : numpy.ndarray
                    Numpy array for which to convert to TACS nodeVec.

        Returns
        -------
        Xptsvec : TACS.Vec
                  Converted TACS nodeVec.

        Notes
        -----
        xptsArray must have correct size on each processor.
        """
        Xptsvec = self.assembler.createNodeVec()

        # Set values
        Xptsvec.getArray()[:] = xptsArray

        # Return as tacs bvec object
        return Xptsvec

    def getNumCoordinates(self):
        """
        Return the number of mesh coordinates on this processor.
        """
        return self.Xpts.getSize()

    ####### Variable methods ########

    def getVarsPerNode(self):
        """
        Get the number of variables per node for the model.
        """
        return self.assembler.getVarsPerNode()

    def getNumOwnedNodes(self):
        """
        Get the number of nodes owned by this processor.
        """
        return self.assembler.getNumOwnedNodes()

    def _arrayToVec(self, varArray):
        """
        Converts a distributed numpy array into a TACS state variable BVec.

        Parameters
        ----------
        varArray : numpy.ndarray
                   Numpy array for which to convert to TACS Vec.

        Returns
        -------
        varVec : TACS.Vec
                 Converted TACS Vec.

        Notes
        -----
        varArray must have correct size on each processor.
        """
        varVec = self.assembler.createVec()

        # Set values
        varVec.getArray()[:] = varArray

        # Return as tacs bvec object
        return varVec

    def getNumVariables(self):
        """
        Return the number of degrees of freedom (states) that are
        on this processor

        Returns
        -------
        nstate : int
            number of states.
        """
        vpn = self.getVarsPerNode()
        nnodes = self.getNumOwnedNodes()
        return vpn * nnodes
