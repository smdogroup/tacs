"""
Mass minimization of uCRM wingbox subject to a constant vertical force
"""
from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt
import os

import openmdao.api as om

from mphys import Multipoint
from mphys.scenario_structural import ScenarioStructural
from tacs.mphys import TacsBuilder
from tacs import elements, constitutive, functions

bdf_file = os.path.join(os.path.dirname(__file__), 'beam_opt.bdf')

# Beam thickness (initial)
t = 0.01            # m
# Beam width
w = 0.05        # m
# Length of beam
L = 1.0

# Material properties
rho = 2500.0        # density kg/m^3
E = 70.0e9          # Young's modulus (Pa)
nu = 0.0           # Poisson's ratio
ys = 350e6        # yield stress

# Shear force applied at tip
V = 1e3

# Callback function used to setup TACS element objects and DVs
def element_callback(dvNum, compID, compDescript, elemDescripts, specialDVs, **kwargs):
    # Setup (isotropic) property and constitutive objects
    prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
    # Set one thickness dv for every property group
    con = constitutive.IsoRectangleBeamConstitutive(prop, t=t, w=w, tNum=dvNum)

    refAxis = np.array([0.0, 1.0, 0.0])
    transform = elements.BeamRefAxisTransform(refAxis)

    # Pass back the appropriate tacs element object
    elem = elements.LinearBeam(transform, con)
    return elem

def problem_setup(scenario_name, fea_assembler, problem):
    """
    Helper function to add fixed forces and eval functions
    to structural problems used in tacs builder
    """

    # Add TACS Functions
    problem.addFunction('mass', functions.StructuralMass)
    problem.addFunction('ks_vmfailure', functions.KSFailure, safetyFactor=1.0,
                        ksWeight=100.0)

    # Add forces to static problem
    problem.addLoadToNodes(101, [0.0, V, 0.0, 0.0, 0.0, 0.0], nastranOrdering=True)


class Top(Multipoint):
    def setup(self):
        tacs_options = {'element_callback': element_callback,
                        'problem_setup': problem_setup,
                        'mesh_file': bdf_file}

        struct_builder = TacsBuilder(tacs_options, coupled=False, write_solution=False)
        struct_builder.initialize(self.comm)
        dv_array = struct_builder.get_initial_dvs()

        dvs = self.add_subsystem('dvs', om.IndepVarComp(), promotes=['*'])
        dvs.add_output('dv_struct', dv_array)

        self.add_subsystem('mesh', struct_builder.get_mesh_coordinate_subsystem())
        self.mphys_add_scenario('tip_shear', ScenarioStructural(struct_builder=struct_builder))
        self.mphys_connect_scenario_coordinate_source('mesh', 'tip_shear', 'struct')

        self.connect('dv_struct', 'tip_shear.dv_struct')


################################################################################
# OpenMDAO setup
################################################################################

prob = om.Problem()
prob.model = Top()
model = prob.model

model.add_design_var('dv_struct', lower=0.001, upper=0.1, scaler=100.0)
model.add_objective('tip_shear.mass', scaler=1.0)
model.add_constraint('tip_shear.ks_vmfailure', lower=0.0, upper=1.0, scaler=1.0)

prob.driver = om.ScipyOptimizeDriver(debug_print=['objs', 'nl_cons'], maxiter=1000)
prob.driver.options['optimizer'] = 'SLSQP'

prob.setup()

prob.run_driver()

x = prob.get_val('mesh.x_struct0', get_remote=True)[:-3:3]
t_opt = prob['dv_struct']
m_opt = prob['tip_shear.mass']
t_exact = np.sqrt(6*(L-x)*V/w/ys)

if __name__ == "__main__" and prob.comm.size == 1:
    # Output N2 representation of OpenMDAO model
    om.n2(prob, show_browser=False, outfile='beam_opt_n2.html')

    # Plot results for solution
    plt.plot(x, t_opt, 'o', x, t_exact)
    plt.legend(['opt', 'exact'])
    plt.ylabel('t(x)')
    plt.xlabel('x')
    plt.title('Optimal beam thickness profile')
    plt.show()
