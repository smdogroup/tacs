import numpy as np
import os
from tacs import pyTACS, TACS, elements, constitutive, functions, problems
from pytacs_analysis_base_test import PyTACSTestCase

"""
Test phase change material thermal model using a single-element mesh to verify phase transition
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/single_element.bdf")

#FUNC_REFS = {'Transient_ks_temp': 7.976969029014029}
FUNC_REFS = {'Transient_ks_temp': 20.70202984802095}

class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    def setup_pytacs(self, comm, dtype):
        """
        Setup mesh and pytacs object for problem we will be testing.
        """

        # Instantiate FEA Assembler
        struct_options = {'outputElement': TACS.PCM_ELEMENT}

        # Instantiate the pyTACS object
        fea_assembler = pyTACS(bdf_file, comm, options=struct_options)

        # Specify the plate thickness
        tplate = 1.0

        # Define material properties
        rho = 1.0  # Density kg/m^3
        kappa = 1.0  # Thermal conductivity W/(m⋅K)
        cp = 1.0  # Specific heat J/(kg⋅K)
        lh = 10.0  # Latent heat J/kg
        Tm = 10.0  # Melting temperature (relative) K

        # The callback function to define the element properties
        def elemCallBack(dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs):

            # Setup property and constitutive objects
            prop = constitutive.MaterialProperties(rho=rho, kappa=kappa, specific_heat=cp)

            # Set one thickness value for every component
            con = constitutive.PhaseChangeMaterialConstitutive(prop, prop, lh=lh, Tm=Tm, t=tplate, tNum=-1)

            # For each element type in this component,
            # pass back the appropriate tacs element object
            elemList = []
            model = elements.PCMHeatConduction2D(con)
            for elemDescript in elemDescripts:
                basis = elements.LinearQuadBasis()
                elem = elements.Element2D(model, basis)
                elemList.append(elem)

            return elemList

        # Set up constitutive objects and elements
        fea_assembler.initialize(elemCallBack)

        return fea_assembler

    def setup_tacs_vecs(self, fea_assembler, dv_pert_vec, xpts_pert_vec):
        """
        Setup user-defined vectors for analysis and fd/cs sensitivity verification
        """
        # Create temporary dv vec for doing fd/cs
        dv_pert_vec[:] = 1.0

        # Define perturbation array that moves all nodes on plate
        xpts = fea_assembler.getOrigNodes()
        xpts_pert_vec[:] = xpts

        return

    def setup_funcs(self, fea_assembler, problems):
        """
        Create a list of functions to be tested and their reference values for the problem
        """
        # Add Functions
        for problem in problems:
            problem.addFunction('ks_temp', functions.KSTemperature,
                                ksWeight=100.0)
        func_list = ['ks_temp']
        return func_list, FUNC_REFS

    def setup_tacs_problems(self, fea_assembler):
        """
        Setup pytacs object for problems we will be testing.
        """
        tacs_probs = []

        # Create transient problem, loads are already applied through BCs
        tp = fea_assembler.createTransientProblem('Transient', tInit=0.0, tFinal=30.0, numSteps=40, options={'L2Convergence':1e-16})
        tacs_probs.append(tp)

        # Get the time steps and define the loads
        timeSteps = tp.getTimeSteps()
        for i, t in enumerate(timeSteps):
            # select the component of the battery undergoing thermal runaway
            compIDs = fea_assembler.selectCompIDs(include=["component"])
            tp.addLoadToComponents(i, compIDs, [1.0])

        return tacs_probs