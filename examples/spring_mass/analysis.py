"""
A simple point mass attached to a spring. A uniform load of 1 N is
applied in every direction on the mass for 0.5 seconds. The response of
the system is then observed for 2*pi seconds. The results are then plotted.

The stiffness values of the spring are:
kx = 1 N/m
ky = 4 N/m
kz = 9 N/m
krx = 16 N*m
kry = 25 N*m
krz = 36 N*m

The inertial values of the mass are:
mass = 1.0 kg
Ixx = Iyy = Izz = 1.0 kg * m^2

The force applied to the mass is:
fx = fy = fz = 1.0 N
mx = my = mz = 1.0 N * m
"""
# ==============================================================================
# Standard Python modules
# ==============================================================================
from __future__ import print_function
import os
import matplotlib.pyplot as plt

# ==============================================================================
# External Python modules
# ==============================================================================
from pprint import pprint
from mpi4py import MPI
import numpy as np

# ==============================================================================
# Extension modules
# ==============================================================================
from tacs import pyTACS, functions, constitutive, elements

# Force vector to apply to mass
f = np.ones(6)

comm = MPI.COMM_WORLD

# Instantiate FEAAssembler
structOptions = {
    "printtiming": True,
}

bdfFile = os.path.join(os.path.dirname(__file__), "mass_spring.bdf")
FEAAssembler = pyTACS(bdfFile, comm, options=structOptions)


def elemCallBack(dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs):
    # Setup spring element
    # NOTE: Mass elements are always setup automatically for users outside of elemCallBack
    k = [1, 4, 9, 16, 25, 36]
    con = constitutive.DOFSpringConstitutive(k=k)
    # Stiffness components are in global coordinate system
    # We don't need a transform
    transform = None
    elem = elements.SpringElement(transform, con)
    return elem


# Set up TACS Assembler
FEAAssembler.initialize(elemCallBack)

# create tacs transient problems
problem = FEAAssembler.createTransientProblem("step_force", 0.0, 2 * np.pi, 100)

# Add functions
problem.addFunction("mass", functions.StructuralMass)
problem.addFunction(
    "max_x_disp", functions.KSDisplacement, ksWeight=100.0, direction=[1.0, 0.0, 0.0]
)

timeSteps = problem.getTimeSteps()
for step_i, time in enumerate(timeSteps):
    # Apply step load for 0.5 seconds
    if time < 0.5:
        problem.addLoadToNodes(step_i, 2, f, nastranOrdering=True)

# Solve problems
problem.solve()
funcs = {}
problem.evalFunctions(funcs)
print(funcs)

# Pick off state variable values at each time step
nnodes = FEAAssembler.getNumOwnedNodes()
vpn = FEAAssembler.getVarsPerNode()
stateHistory = np.zeros([len(timeSteps), nnodes, vpn])
for step_i in range(len(timeSteps)):
    problem.getVariables(step_i, states=stateHistory[step_i, :, :].reshape(-1))

# Plot results for first 3 dofs
plt.plot(
    timeSteps,
    stateHistory[:, 1, 0],
    timeSteps,
    stateHistory[:, 1, 1],
    timeSteps,
    stateHistory[:, 1, 2],
)
plt.legend(["dof 1", "dof 2", "dof 3"])
plt.ylabel("displacement (m)")
plt.xlabel("time (s)")
if __name__ == "__main__":
    plt.show()
