"""
Author:
    - Sean P Engelstad, Alasdair Christison Gray

This class implements a constraint which enforces the panel
width design variable values passed to elements using the GPBladeStiffenedShell
constitutive model to be consistent with the true width of the panel they are
a part of.

.. note:: This class should be created using the
    :meth:`pyTACS.createPanelWidthConstraint <tacs.pytacs.pyTACS.createPanelWidthConstraint>` method.
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================

# ==============================================================================
# External Python modules
# ==============================================================================
import numpy as np
import scipy as sp

# ==============================================================================
# Extension modules
# ==============================================================================
from .panel_length import *


class PanelWidthConstraint(PanelLengthConstraint):
    def __init__(
        self,
        name,
        assembler,
        comm,
        outputViewer=None,
        meshLoader=None,
        options=None,
    ):
        """
        NOTE: This class should not be initialized directly by the user.
        Use pyTACS.createPanelWidthConstraint instead.

        Parameters
        ----------
        name : str
            Name of this tacs problem

        assembler : TACS.Assembler
            Cython object responsible for creating and setting tacs objects used to solve problem

        comm : mpi4py.MPI.Intracomm
            The comm object on which to create the pyTACS object.

        outputViewer : TACS.TACSToFH5
            Cython object used to write out f5 files that can be converted and used for postprocessing.

        meshLoader : pymeshloader.pyMeshLoader
            pyMeshLoader object used to create the assembler.

        options : dict
            Dictionary holding problem-specific option parameters (case-insensitive).
        """

        super(PanelWidthConstraint, self).__init__(
            name, assembler, comm, outputViewer, meshLoader, options
        )

    def computeRefAxis(self, refAxis: np.ndarray, comp_bndry_node_coords: np.ndarray):
        """
        remove the panelNormal from the refAxis

        Parameters
        ----------
        refAxis : numpy.ndarray
            an array of size (3,) for the xyz coordinates of the original refAxis from transform object
        comp_bndry_node_coords : numpy.ndarray
            an array of size (3*N,) for the boundary nodal coordinates on the current panel / current TACS component

        Returns
        -------
        numpy.ndarray
            an array of size (3,) for xyz coords of the panel width axis in the panel
        """

        # For a more accurate length calculation, roject the ref axis
        # onto the "average" plane of the baseline panel geometry by
        # using an SVD to compute a normal vector
        centroid = np.mean(comp_bndry_node_coords, axis=0, keepdims=True)
        centredPoints = comp_bndry_node_coords - centroid
        _, _, VT = np.linalg.svd(centredPoints, full_matrices=False)
        panelNormal = VT[-1]
        refAxis -= np.dot(refAxis, panelNormal) * panelNormal
        refAxis /= np.linalg.norm(refAxis)

        # now this refAxis is the panel length axis, compute the width axis which is normal to the length + panel normal axis
        widthAxis = np.cross(refAxis, panelNormal)
        widthAxis /= np.linalg.norm(widthAxis)
        return widthAxis
